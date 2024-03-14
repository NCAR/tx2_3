!***********************************************************************

program create_lmask

  !-----------------------------------------------------------------------
  !
  !     This program creates a LMASK field for a MOM6 ocean grid by 
  !     interpolating/averaging from binary lat-lon land/ocean mask
  !
  !     This version uses a Monte-Carlo-like method 
  !     (based on a code from Mat Maltrud used for POP)
  !
  !-----------------------------------------------------------------------

  use kinds
  use constants
  use ncdf_wrapper
  use omp_lib

  implicit none

  integer, parameter :: int1=selected_int_kind(1),&
       &                int2=selected_int_kind(4),&
       &                int4=selected_int_kind(6)

  !-----------------------------------------------------------------------

  ! ETOPO 2022 / SRTM /MODIS 15-second grid
  integer (kind=int_kind) :: nx_lmask, ny_lmask
  real(kind=dbl_kind) :: dlon_lmask, dlat_lmask
  real (kind=dbl_kind), dimension(:), allocatable :: LMASK_IN_LON
  real (kind=dbl_kind), dimension(:), allocatable :: LMASK_IN_LAT
  integer(kind=int_kind), dimension(:,:), allocatable :: INDAT, LMASK_IN ! input data

  integer (kind=int_kind) :: &
       &  nx_sg, ny_sg &         ! size of MOM6 super-grid
       &, nx_o, ny_o   &         ! size of MOM6 model grid
       &, nx_sub, ny_sub         ! number of sub-points to use in x,y dir

  integer (kind=int_kind), dimension(:,:), allocatable :: &
       &  MASK_OUT                    ! LMASK field for ocean grid

  real (kind=real_kind), dimension(:,:), allocatable :: &
       &  OCEAN_FRAC                    ! LMASK field for ocean grid

  real (kind=real_kind) :: lf_threshold

  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       &  X_SG  &    ! lon of super-grid points
       &, Y_SG        ! lat of super-grid points

  real (kind=dbl_kind), dimension(:,:,:), allocatable :: &
       &  OCN_CORNER_LAT  &    ! lat of ocean T-cell corner (Q) points
       &, OCN_CORNER_LON       ! lon of ocean T-cell corner (Q) points

  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       &  OCN_CENTER_LAT  &    ! lat of ocean T-cell corner (Q) points
       &, OCN_CENTER_LON       ! lon of ocean T-cell corner (Q) points

  real (kind=dbl_kind), dimension(4) :: &
       &  corner_lat, corner_lon  ! temps to hold local corner values

  !---- netCDF stuff
  integer :: ierr, ncid, dimid, vid, ndim
  integer :: dimid_lonh, dimid_lath, dimid_lonq, dimid_latq 
  integer :: vid_geolon, vid_geolonb, &
       &      vid_geolat, vid_geolatb, &
       &      vid_mask, vid_ofrac
  integer, dimension(2) :: indcnt, indstr, dims
  character (char_len) :: cbuffer
  character (char_len) :: vname_lon, vname_lat, vname_mask
  character(char_len) :: dim_name, var_name, attr_name, attr_valc

  character (char_len) :: &
       &  grid_version     &
       &, horiz_grid_file   &    ! file with horizontal ocean grid
       &, lmask_in_file    &     ! file with input binary mask data
       &, lmask_model_file       ! file for resulting KMT field

  logical :: is_tripole=.true.

  integer (kind=int_kind) :: &
       &  i,j,k,im1,jm1, jrev,jprnt=1000   &
       &, ocn_i, ocn_j    &
       &, lmask_i, lmask_j  &
       &, isub, jsub  &
       &, nsum, nsumo &
       &, np

  real (kind=dbl_kind) :: &
       &  ifrac, jfrac         &
       &, rlat, rlon           &
       &, dlon, dlat     &
       &, frac

  namelist /model_grid_nml/nx_sub, ny_sub, grid_version, horiz_grid_file, is_tripole

  namelist /lmask_in_nml/ lmask_in_file, vname_lon, vname_lat, vname_mask

  namelist /output_nml/ lf_threshold, jprnt, lmask_model_file

  !-----------------------------------------------------------------------
  !     Initialize and get namelist input
  !-----------------------------------------------------------------------

  call init_constants

  nx_sub = 1
  ny_sub = 1
  horiz_grid_file = 'unknown-horiz-grid'
  lmask_in_file    = 'unknown-lmask_data'
  lmask_model_file = 'unknown-lmask-file'
  lf_threshold = 0.5

  read(*,nml=model_grid_nml)
  write(*,nml=model_grid_nml)

  read(*,nml=lmask_in_nml)
  write(*,nml=lmask_in_nml)

  read(*,nml=output_nml)
  write(*,nml=output_nml)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Get model super grid and allocate model grid arrays
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  print *,'Opening ',horiz_grid_file
  ierr = nf_open(horiz_grid_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr=',ierr
     print *,'could not open ',horiz_grid_file
     stop
  else
     print *,' grid file ncid = ',ncid
  endif

  ierr = nf_inq_dimid(ncid,'nx',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension nx. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_sg)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  else
     print *,' super grid nx=',nx_sg
  endif

  ierr = nf_inq_dimid(ncid,'ny',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ny. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_sg)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  else
     print *,' super grid ny=',ny_sg
  endif

  allocate(X_SG(nx_sg+1,ny_sg+1),Y_SG(nx_sg+1,ny_sg+1))
  
  ierr = nf_inq_varid(ncid,'x',vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable x. ierr= ',ierr
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_sg+1,ny_sg+1/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,X_SG)
  if ( ierr /= 0 ) then
     print *,'could not get variable x ierr= ',ierr
     stop
  endif

  ierr = nf_inq_varid(ncid,'y',vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable y ierr= ',ierr
     stop
  endif

  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,Y_SG)
  if ( ierr /= 0 ) then
     print *,'could not get variable y ierr= ',ierr
     stop
  endif

  !---- The output array
  nx_o = nx_sg/2
  ny_o = ny_sg/2
  allocate(OCEAN_FRAC(nx_o,ny_o),MASK_OUT(nx_o,ny_o),stat=ierr)
  if (ierr /=0 ) then
     print *,'Could not allocate model grid mask ierr=',ierr
  else
     print *,'Model grid dimension = ',nx_o,ny_o
  endif

  !---- data structure to hold corners of each model cell
  allocate(OCN_CORNER_LAT(4,nx_o,ny_o),OCN_CORNER_LON(4,nx_o,ny_o),&
       &   OCN_CENTER_LON(nx_o,ny_o),OCN_CENTER_LAT(nx_o,ny_o))

  !-----------------------------------------------------------------------
  !     read binary land/ocean mask
  !-----------------------------------------------------------------------

  print *,'Opening ',lmask_in_file
  ierr = nf_open(lmask_in_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr = ',ierr,' ncid = ',ncid
     print *,'could not open ',lmask_in_file
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_lon,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',vname_lon(1:len_trim(vname_lon)),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_lmask)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',vname_lon(1:len_trim(vname_lon)),' ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_lat,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',vname_lat(1:len_trim(vname_lat)),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_lmask)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',vname_lat(1:len_trim(vname_lat)),' ierr= ',ierr
     stop
  endif
  print *,'Input mask dimension = ',nx_lmask,ny_lmask

  allocate(LMASK_IN_LON(0:nx_lmask+1))
  
  ierr = nf_inq_varid(ncid,vname_lon,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',vname_lon(1:len_trim(vname_lon)),' ierr= ',ierr
     stop
  endif

  indstr = (/1,0/)
  indcnt = (/nx_lmask,0/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,LMASK_IN_LON(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',vname_lon(1:len_trim(vname_lon)),' ierr= ',ierr
     stop
  endif

  LMASK_IN_LON(0) = LMASK_IN_LON(nx_lmask)-360.
  LMASK_IN_LON(nx_lmask+1) = LMASK_IN_LON(1)+360.
  dlon_lmask = LMASK_IN_LON(2)-LMASK_IN_LON(1)

  write(*,"(a,2f16.7)") 'Extended LMASK LON min/max=',minval(LMASK_IN_LON),maxval(LMASK_IN_LON)
  write(*,"(a,f20.10)")' dlon_lmask = ',dlon_lmask,' degrees'
  write(*,"(a,4f9.4,a,4f9.4)") &
       & 'LON = ',LMASK_IN_LON(0:3),' ... ',LMASK_IN_LON(nx_lmask-2:nx_lmask+1)
  print *

  allocate(LMASK_IN_LAT(0:ny_lmask+1))

  ierr = nf_inq_varid(ncid,vname_lat,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',vname_lat(1:len_trim(vname_lat)),' ierr= ',ierr
     stop
  endif

  indstr = (/1,0/)
  indcnt = (/ny_lmask,0/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,LMASK_IN_LAT(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',vname_lat(1:len_trim(vname_lat)),' ierr= ',ierr
     stop
  endif

  ierr = nf_inq_varid(ncid,vname_lat,vid)
  indstr = (/1,0/)
  indcnt = (/ny_lmask,0/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,LMASK_IN_LAT(1))
  dlat_lmask = LMASK_IN_LAT(2)-LMASK_IN_LAT(1)
  LMASK_IN_LAT(0) = LMASK_IN_LAT(1) - dlat_lmask
  LMASK_IN_LAT(ny_lmask+1) = LMASK_IN_LAT(ny_lmask) + dlat_lmask

  write(*,"(a,2f16.7)") 'Extended LMASK LAT min/max=',minval(LMASK_IN_LAT),maxval(LMASK_IN_LAT)
  write(*,"(a,f20.10)")' dlat = ',dlat_lmask,' degrees'
  write(*,"(a,4f9.4,a,4f9.4)") &
       & 'LAT = ',LMASK_IN_LAT(0:3),' ... ',LMASK_IN_LAT(ny_lmask-2:ny_lmask+1)
  print *


  allocate(INDAT(nx_lmask,ny_lmask),LMASK_IN(0:nx_lmask+1,0:ny_lmask+1),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'allocate stat=',ierr
     stop
  end if

  ierr = nf_inq_varid(ncid,vname_mask,vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find mask variable ',vname_mask
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_lmask,ny_lmask/)
  ierr = nf_get_vara_int(ncid,vid,indstr,indcnt,INDAT)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get mask variable ',vname_mask
     stop
  endif
  LMASK_IN(1:nx_lmask,1:ny_lmask) = INDAT

  !---- Extend the intput grid one point in each direction
  !     so all possible lat-lons are included
  LMASK_IN(0,1:ny_lmask) = INDAT(nx_lmask,:)
  LMASK_IN(nx_lmask,1:ny_lmask) = INDAT(1,:)
  LMASK_IN(:,0) = LMASK_IN(:,1)
  LMASK_IN(:,ny_lmask+1) = LMASK_IN(:,ny_lmask)
  
  print *,' LMASK_IN min=',minval(LMASK_IN),' max=',maxval(LMASK_IN)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Create data structure with corner points for each cell
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !---- corners are numbered 1-4 counterclockwise from bottom left
  OCN_CORNER_LON(1,:,:) = x_sg(1::2,1::2)
  OCN_CORNER_LAT(1,:,:) = y_sg(1::2,1::2)
  OCN_CORNER_LON(2,:,:) = x_sg(3::2,1::2)
  OCN_CORNER_LAT(2,:,:) = y_sg(3::2,1::2)
  OCN_CORNER_LON(3,:,:) = x_sg(3::2,3::2)
  OCN_CORNER_LAT(3,:,:) = y_sg(3::2,3::2)
  OCN_CORNER_LON(4,:,:) = x_sg(1::2,3::2)
  OCN_CORNER_LAT(4,:,:) = y_sg(1::2,3::2)

  OCN_CENTER_LON = x_sg(2::2,2::2)
  OCN_CENTER_LAT = y_sg(2::2,2::2)
  
  !---- Convert corner point lat/lon to radians
  OCN_CORNER_LON = OCN_CORNER_LON/radian
  OCN_CORNER_LAT = OCN_CORNER_LAT/radian

  !---- Make sure all corners are in the same period of 2 pi as corner 3
  !     This logic may not work near the tripole grid poles?
  !     Handle tripole pole points as special case below. They should always be land
  do j=1,ny_o
     do i=1,nx_o

        dlon = OCN_CORNER_LON(4,i,j) - OCN_CORNER_LON(3,i,j)
        if (dlon > c3*pih) OCN_CORNER_LON(4,i,j) = OCN_CORNER_LON(4,i,j) - pi2
        if (dlon <-c3*pih) OCN_CORNER_LON(4,i,j) = OCN_CORNER_LON(4,i,j) + pi2

        dlon = OCN_CORNER_LON(2,i,j) - OCN_CORNER_LON(3,i,j)
        if (dlon > c3*pih) OCN_CORNER_LON(2,i,j) = OCN_CORNER_LON(2,i,j) - pi2
        if (dlon <-c3*pih) OCN_CORNER_LON(2,i,j) = OCN_CORNER_LON(2,i,j) + pi2

        dlon = OCN_CORNER_LON(1,i,j) - OCN_CORNER_LON(3,i,j)
        if (dlon > c3*pih) OCN_CORNER_LON(1,i,j) = OCN_CORNER_LON(1,i,j) - pi2
        if (dlon <-c3*pih) OCN_CORNER_LON(1,i,j) = OCN_CORNER_LON(1,i,j) + pi2

     enddo
  enddo

  print *,' after cell consistency check :'
  print *,' OCN_CORNER LON min/max = ',minval(OCN_CORNER_LON),&
       &                               maxval(OCN_CORNER_LON)
  print *,' OCN_CORNER LAT min/max = ',minval(OCN_CORNER_LAT),&
       &                               maxval(OCN_CORNER_LAT)


  !-----------------------------------------------------------------------
  !
  !     loop through each ocean point to determine LMASK at each ocean
  !     grid cell
  !
  !-----------------------------------------------------------------------

  !---- Convert input data coordinates to radians
  LMASK_IN_LON = LMASK_IN_LON/radian
  dlon_lmask = dlon_lmask/radian
  LMASK_IN_LAT = LMASK_IN_LAT/radian
  dlat_lmask = dlat_lmask/radian

  MASK_OUT = 0
  OCEAN_FRAC = 0.

  print *,' Before threaded loop. max_threads=',omp_get_max_threads()

  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(DYNAMIC,8) &
  !$OMP& PRIVATE(ocn_j,ocn_i,isub,jsub,ifrac,jfrac,rlat,rlon, &
  !$OMP&         corner_lat,corner_lon,lmask_i,lmask_j,nsum,nsumo,frac)
  do ocn_j=1,ny_o

     if ( ocn_j == 1 ) then
        print *,'Beginning threaded loop. num_threads=',omp_get_num_threads()
     endif
     if ( mod(ocn_j-1,jprnt) == 0 ) then
        print *,'Doing model j = ',ocn_j,' thread=',OMP_get_thread_num(),&
             &    'OCN_CORNER_LON(:,1,',ocn_j,')=',OCN_CORNER_LON(:,1,ocn_j),&
             &    'OCN_CORNER_LAT(:,1,',ocn_j,')=',OCN_CORNER_LAT(:,1,ocn_j)
     endif

     do ocn_i=1,nx_o

        corner_lat(:) = OCN_CORNER_LAT(:,ocn_i,ocn_j)
        corner_lon(:) = OCN_CORNER_LON(:,ocn_i,ocn_j)

        !*** distribute sub points
        nsum = 0
        nsumo = 0
        do jsub=1,ny_sub
           jfrac = real(jsub)/real(ny_sub+1)
           do isub=1,nx_sub
              ifrac = real(isub)/real(nx_sub+1)

              rlat = (c1-ifrac)*(c1-jfrac)*corner_lat(1)         &
                   &           +     ifrac *(c1-jfrac)*corner_lat(2)  &
                   &           +     ifrac *    jfrac *corner_lat(3)  &
                   &           + (c1-ifrac)*    jfrac *corner_lat(4)

              rlon = (c1-ifrac)*(c1-jfrac)*corner_lon(1)         &
                   &           +     ifrac *(c1-jfrac)*corner_lon(2)  &
                   &           +     ifrac *    jfrac *corner_lon(3)  &
                   &           + (c1-ifrac)*    jfrac *corner_lon(4)

              ! Remap longitudes into -pi to pi
              do while ( rlon < LMASK_IN_LON(0) )
                 rlon = rlon + pi2
              enddo
              do while ( rlon > LMASK_IN_LON(nx_lmask+1)) 
                 rlon = rlon - pi2
              enddo

              if ( rlat < LMASK_IN_LAT(0) .OR. rlat > LMASK_IN_LAT(ny_lmask+1)) then
                 print *,' Invalid search latitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lon = ',corner_lon
                 print *,' corner lat = ',corner_lat
                 print *,' rlon = ',rlon
                 print *,' rlat = ',rlat
                 stop
              endif

              !*** find location of this point on hi-res lmask grid
              lmask_i = floor( (rlon - LMASK_IN_LON(0))/dlon_lmask )+1
              lmask_j = floor( (rlat - LMASK_IN_LAT(0))/dlat_lmask )+1

              if ( lmask_i < 0 .OR. lmask_i > nx_lmask+1) then
                 print *,' Invalid search longitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lon = ',corner_lon
                 print *,' corner lat = ',corner_lat
                 print *,' rlon = ',rlon
                 print *,' lmask_i= ',lmask_i
                 stop
              endif

              if ( lmask_j < 0 .OR. lmask_j > ny_lmask+1) then
                 print *,' Invalid search latitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lat = ',corner_lat
                 print *,' corner lat = ',corner_lat
                 print *,' rlat = ',rlat
                 print *,' lmask_j = ',lmask_j
                 stop
              endif
!!$              print *,' lmask_i = ',lmask_i,' lmask_j = ',lmask_j

              !*** accumulate a wet point sum if search was successful
              nsum = nsum + 1
              if (LMASK_IN(lmask_i,lmask_j) > 0) then
                 nsumo = nsumo + 1
              endif

           enddo
        enddo

        !if more than half of the sub-grid points are ocean, call this model point ocean
        if (nsum /= 0) then
           OCEAN_FRAC(ocn_i,ocn_j) = float(nsumo)/float(nsum)
           if (OCEAN_FRAC(ocn_i,ocn_j) >= lf_threshold) MASK_OUT(ocn_i,ocn_j) = 1
        else
           OCEAN_FRAC(ocn_i,ocn_j) = 0.
           MASK_OUT(ocn_i,ocn_j) = 0
        endif

        if(mod(ocn_j-1,jprnt) == 0 .AND. ocn_i == 1) then
           print *,' MASK_OUT(1,j)=',MASK_OUT(1,ocn_j),' # points =',nsum,' # ocn points =',nsumo
        endif

     enddo
  enddo
!$OMP END PARALLEL DO

  !---- Fix possible problems at grid poles
  if ( is_tripole ) then
     np = nx_o/2
     print *,'mid pole i=',np

     print *,' Before tripole fix LMASK(1:3)=',MASK_OUT(1:3,ny_o)
     print *,' Before tripole fix LMASK(np-1:np+1)=',MASK_OUT(np-1:np+1,ny_o)
     print *,' Before tripole fix LMASK(nx_o-2:nx_o)=',MASK_OUT(nx_o-2:nx_o,ny_o)

     OCEAN_FRAC(1,ny_o)=0
     OCEAN_FRAC(np,ny_o)=0
     OCEAN_FRAC(nx_o,ny_o) = 0

     MASK_OUT(1,ny_o)=0
     MASK_OUT(np,ny_o)=0
     MASK_OUT(nx_o,ny_o) = 0

     print *,' After tripole fix LMASK(1:3)=',MASK_OUT(1:3,ny_o)
     print *,' After tripole fix LMASK(np-1:np+1)=',MASK_OUT(np-1:np+1,ny_o)
     print *,' After tripole fix LMASK(nx_o-2:nx_o)=',MASK_OUT(nx_o-2:nx_o,ny_o)
  endif

  print *,' done finding MASK_OUT'


  !-----------------------------------------------------------------------
  !
  !     write LMASK field to file
  !
  !-----------------------------------------------------------------------

  print *,'OCEAN_FRAC min,max,mean ',minval(OCEAN_FRAC),maxval(OCEAN_FRAC),&
       &   sum(OCEAN_FRAC)/real(nx_o*ny_o)

  print *,'MASK_OUT min,max,mean ',minval(MASK_OUT),maxval(MASK_OUT),&
       &   sum(MASK_OUT)/real(nx_o*ny_o)

  ierr = nf_create(trim(lmask_model_file),NF_64bit_OFFSET,ncid)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_create')
     stop
  endif

  !---- Global attributes
  attr_name = 'title'
  attr_valc = 'Ocean Landmask'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'Generating Code'
  attr_valc = 'create_model_lmask.f90'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'Grid Version'
  attr_valc = grid_version
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'Source Mask'
  attr_valc = lmask_in_file
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'Revision History'
  call date_and_time(cbuffer)
  attr_valc = 'Created ' // trim(cbuffer)
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  !---- Define dimesnions
  ierr = nf_def_dim(ncid,'lonh',nx_o,dimid_lonh)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_dim lon')
     stop
  endif
  ierr = nf_def_dim(ncid,'lath',ny_o,dimid_lath)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_dim lat')
     stop
  endif
  ierr = nf_def_dim(ncid,'lonq',nx_o+1,dimid_lonq)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_dim lon')
     stop
  endif
  ierr = nf_def_dim(ncid,'latq',ny_o+1,dimid_latq)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_dim lat')
     stop
  endif

  !---- Define cell center variables
  ndim = 2
  dims = (/dimid_lonh,dimid_lath/)
  ierr = nf_def_var(ncid,'geolon',NF_DOUBLE,ndim,dims,vid_geolon)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolon')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Longitude of tracer (T) points'
  ierr = nf_put_att_text(ncid,vid_geolon,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_geolon,attr_name,len_trim(attr_valc),trim(attr_valc))


  ierr = nf_def_var(ncid,'geolat',NF_DOUBLE,ndim,dims,vid_geolat)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolat')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Latitude of tracer (T) points'
  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_north'
  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'ocean_frac',NF_REAL,ndim,dims,vid_ofrac)
  attr_name = 'longname'
  attr_valc = 'fraction of cell area covered by ocean'
  ierr = nf_put_att_text(ncid,vid_ofrac,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'dimensionless'
  ierr = nf_put_att_text(ncid,vid_ofrac,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'mask',NF_INT,ndim,dims,vid_mask)
  attr_name = 'longname'
  attr_valc = 'Land/Ocean Mask: 0 if land, 1 if ocean at tracer points'
  ierr = nf_put_att_text(ncid,vid_mask,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'Ocean Area Threshold'
  ierr = nf_put_att_real(ncid,NF_GLOBAL,len_trim(attr_name),1,lf_threshold)
  attr_name = 'units'
  attr_valc = 'dimensionless'
  ierr = nf_put_att_text(ncid,vid_mask,attr_name,len_trim(attr_valc),trim(attr_valc))

  !---- Define cell corner variables
  ndim = 2
  dims = (/dimid_lonq,dimid_latq/)
  ierr = nf_def_var(ncid,'geolonb',NF_DOUBLE,ndim,dims,vid_geolonb)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolon')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Longitude of corner (Q) points'
  ierr = nf_put_att_text(ncid,vid_geolonb,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_geolonb,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'geolatb',NF_DOUBLE,ndim,dims,vid_geolatb)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolatb')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Latitude of corner (Q) points'
  ierr = nf_put_att_text(ncid,vid_geolatb,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_north'
  ierr = nf_put_att_text(ncid,vid_geolatb,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_enddef(ncid)

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_double(ncid,vid_geolon,indstr,indcnt,OCN_CENTER_LON)
  ierr = nf_put_vara_double(ncid,vid_geolat,indstr,indcnt,OCN_CENTER_LAT)

  ierr = nf_put_vara_int(ncid,vid_mask,indstr,indcnt,MASK_OUT)
  ierr = nf_put_vara_real(ncid,vid_ofrac,indstr,indcnt,OCEAN_FRAC)

  indstr = (/1,1/)
  indcnt = (/nx_o+1,ny_o+1/)
  ierr = nf_put_vara_double(ncid,vid_geolonb,indstr,indcnt,X_SG(1::2,1::2))
  ierr = nf_put_vara_double(ncid,vid_geolatb,indstr,indcnt,Y_SG(1::2,1::2))

  ierr = nf_close(ncid)

  !-----------------------------------------------------------------------

end program create_lmask

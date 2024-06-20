!***********************************************************************

program create_model_topo

  !-----------------------------------------------------------------------
  !
  !     This program creates MASK and TOPO fields for a MOM6 ocean grid by 
  !     interpolating/averaging from binary lat-lon land/ocean mask
  !
  !     This version uses a Monte-Carlo-like method 
  !     (based on a code from Mat Maltrud used for POP)
  ! 
  !     The input topography data set is assumed to be on a uniformly
  !     spaced lat-long grid.
  !     It is assumed that positive input elevations are above sea-level and
  !     negative elevations are below sea-level
  !
  !-----------------------------------------------------------------------

  use kinds
  use constants
  use ncdf_wrapper
  use mom6_grid
  use omp_lib

  implicit none

  interface
     subroutine sort(x,N)
       use kinds
       real(kind=real_kind), dimension(1:), intent(inout) :: x
       integer, intent(in)                   :: N
     end subroutine sort
     recursive subroutine quicksort(a)
       use kinds
       implicit none
       real(kind=dbl_kind) :: a(:)
     end subroutine quicksort
  end interface

  !-----------------------------------------------------------------------

  !---- Input topography data set
  integer (kind=int_kind) :: nx_topo, ny_topo
  real(kind=dbl_kind) :: dlon_topo, dlat_topo
  real (kind=dbl_kind), dimension(:), allocatable :: TOPO_LON
  real (kind=dbl_kind), dimension(:), allocatable :: TOPO_LAT
  real(kind=real_kind), dimension(:,:), allocatable :: INDAT, TOPO_Z ! input data

  !---- Results on model grid
  integer (kind=int_kind), dimension(:,:), allocatable :: &
       &  OCN_MASK                    ! Binary mask field for ocean (0=land, 1=ocean)

  real (kind=real_kind), dimension(:,:), allocatable :: &
       &  OCN_FRAC, &                     ! Fraction of cell area covered by ocean
       &  AVG_TOPO,  &                      ! Avg elevation (land and ocean) in cell
       &  OCN_D_MEAN, OCN_D_MEDIAN, &        ! Mean and median ocean depth in cell
       &  OCN_D_MIN,  OCN_D_MAX,&      ! Min and Max depth in cell
       &  OCN_D2_MEAN                      ! Mean of depth squared in cell

  !---- parameters
  real (kind=real_kind) :: &
       & mask_hmin, &
       & mask_threshold, &
       & depth_clip_min, depth_clip_max

  integer(kind=int_kind) ::&
       & nx_sub, ny_sub, &
       & ichk=1, jprnt=1000

  logical :: verbose=.false., do_median=.false.

  character (char_len) :: &
       &  mom6_grid_version     &
       &, mom6_horiz_grid_file   &    ! file with horizontal ocean grid
       &, topo_in_file    &     ! file with input high-res topo data
       &, output_file    &    ! file for results
       &, user_name
  character (char_len) :: vname_topo_lon, vname_topo_lat, vname_topo_z

  logical :: grid_is_tripole=.true.

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !    Work variables
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  real (kind=dbl_kind), dimension(:,:,:), allocatable :: OCN_CORNER_LAT, OCN_CORNER_LON
  real(kind=dbl_kind), dimension(:), allocatable :: rbuffer
  character (char_len) :: cbuffer

  !---- netCDF stuff
  integer :: ierr, ncid, dimid, vid, ndim
  integer :: dimid_lonh, dimid_lath, dimid_lonq, dimid_latq 
  integer, dimension(2) :: indcnt, indstr, dims
  character(char_len) :: dim_name, var_name, attr_name, attr_valc
  integer :: vid_topo, vid_ofrac, vid_mask, &
       &     vid_lonh, vid_lath, vid_lonq, vid_latq, &
       &     vid_geolon, vid_geolat, vid_geolonb, vid_geolatb, &
       &     vid_z, vid_d_mean, vid_d_median, vid_d2_mean, vid_d_min, vid_d_max

  integer (kind=int_kind) :: &
       &  i,j,k,im1,jm1 &
       &, ocn_i, ocn_j    &
       &, topo_i, topo_j  &
       &, isub, jsub  &
       &, nsum, nsumo &
       &, np

  real (kind=dbl_kind) :: &
       &  ifrac, jfrac         &
       &, rlat, rlon           &
       &, dlon, dlat     &
       &, frac, rsum

  namelist /model_grid_nml/nx_sub, ny_sub, mom6_grid_version, mom6_horiz_grid_file, grid_is_tripole

  namelist /topo_in_nml/ topo_in_file, vname_topo_lon, vname_topo_lat, vname_topo_z

  namelist /output_nml/ verbose, jprnt, ichk, do_median, &
       &                mask_hmin, mask_threshold, depth_clip_min, depth_clip_max, &
       &                user_name, output_file

  !-----------------------------------------------------------------------
  !     Initialize and get namelist input
  !-----------------------------------------------------------------------

  call init_constants

  nx_sub = 1
  ny_sub = 1
  mom6_grid_version = 'unknown-grid-version'
  mom6_horiz_grid_file = 'unknown-horiz-grid'
  topo_in_file    = 'unknown-topo-data-file'
  vname_topo_lon = 'lon'
  vname_topo_lat = 'lat'
  vname_topo_z = 'z'
  user_name = 'unknown-user'
  output_file = 'unknown-output-file'

  mask_threshold = 0.5
  mask_hmin = 0.0
  depth_clip_min = 0.0
  depth_clip_max = 1.0e5

  read(*,nml=model_grid_nml)
  write(*,nml=model_grid_nml)

  read(*,nml=topo_in_nml)
  write(*,nml=topo_in_nml)

  read(*,nml=output_nml)
  write(*,nml=output_nml)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Get model super grid and allocate model grid arrays
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call read_super_grid(mom6_horiz_grid_file)
  call compute_mom6_grid()

  !-----------------------------------------------------------------------
  !     read hi-res topography data
  !-----------------------------------------------------------------------

  print *,'Opening ',topo_in_file
  ierr = nf_open(topo_in_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr = ',ierr,'could not open ',topo_in_file
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_topo_lon,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_topo_lon),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_topo)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_topo_lon),' ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_topo_lat,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_topo_lat),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_topo)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_topo_lat),' ierr= ',ierr
     stop
  endif
  print *,'Input topo dimension = ',nx_topo,ny_topo

  allocate(TOPO_LON(0:nx_topo+1))

  ierr = nf_inq_varid(ncid,vname_topo_lon,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',trim(vname_topo_lon),' ierr= ',ierr
     stop
  endif

  indstr = (/1,0/)
  indcnt = (/nx_topo,0/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,TOPO_LON(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',vname_topo_lon(1:len_trim(vname_topo_lon)),' ierr= ',ierr
     stop
  endif

  TOPO_LON(0) = TOPO_LON(nx_topo)-360.
  TOPO_LON(nx_topo+1) = TOPO_LON(1)+360.
  dlon_topo = TOPO_LON(2)-TOPO_LON(1)

  if ( verbose ) then
     write(*,"(a,2f16.7)") 'Extended TOPO LON min/max=',minval(TOPO_LON),maxval(TOPO_LON)
     write(*,"(a,f20.10)") 'dlon_topo = ',dlon_topo,' degrees'
     write(*,"(a,4f9.4,a,4f9.4)") &
          & 'LON = ',TOPO_LON(0:3),' ... ',TOPO_LON(nx_topo-2:nx_topo+1)
     print *
  endif

  allocate(TOPO_LAT(0:ny_topo+1))

  ierr = nf_inq_varid(ncid,vname_topo_lat,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',vname_topo_lat(1:len_trim(vname_topo_lat)),' ierr= ',ierr
     stop
  endif

  indstr = (/1,0/)
  indcnt = (/ny_topo,0/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,TOPO_LAT(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',vname_topo_lat(1:len_trim(vname_topo_lat)),' ierr= ',ierr
     stop
  endif

  dlat_topo = TOPO_LAT(2)-TOPO_LAT(1)
  TOPO_LAT(0) = TOPO_LAT(1) - dlat_topo
  TOPO_LAT(ny_topo+1) = TOPO_LAT(ny_topo) + dlat_topo

  if ( verbose ) then
     write(*,"(a,2f16.7)") 'Extended TOPO LAT min/max=',minval(TOPO_LAT),maxval(TOPO_LAT)
     write(*,"(a,f20.10)")' dlat = ',dlat_topo,' degrees'
     write(*,"(a,4f9.4,a,4f9.4)") &
          & 'LAT = ',TOPO_LAT(0:3),' ... ',TOPO_LAT(ny_topo-2:ny_topo+1)
     print *
  endif


  allocate(INDAT(nx_topo,ny_topo),TOPO_Z(0:nx_topo+1,0:ny_topo+1),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'allocate stat=',ierr
     stop
  end if

  ierr = nf_inq_varid(ncid,vname_topo_z,vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find topo height variable ',vname_topo_z
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_topo,ny_topo/)
  ierr = nf_get_vara_real(ncid,vid,indstr,indcnt,INDAT)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get topo height variable ',vname_topo_z
     stop
  endif
  TOPO_Z(1:nx_topo,1:ny_topo) = INDAT

  !---- Extend the intput grid one point in each direction
  !     so all possible lat-lons are included
  TOPO_Z(0,1:ny_topo) = INDAT(nx_topo,:)
  TOPO_Z(nx_topo,1:ny_topo) = INDAT(1,:)
  TOPO_Z(:,0) = TOPO_Z(:,1)
  TOPO_Z(:,ny_topo+1) = TOPO_Z(:,ny_topo)

  print *,' TOPO_LON min=',minval(TOPO_LON),' max=',maxval(TOPO_LON)
  print *,' TOPO_LAT min=',minval(TOPO_LAT),' max=',maxval(TOPO_LAT)
  print *,' TOPO_Z min=',minval(TOPO_Z),' max=',maxval(TOPO_Z)

  !---- Convert input data coordinates to radians
  TOPO_LON = TOPO_LON/radian
  dlon_topo = dlon_topo/radian
  TOPO_LAT = TOPO_LAT/radian
  dlat_topo = dlat_topo/radian

  print *,'Avg. number of data points / grid cell x dir = ',float(nx_topo)/float(nx_o)
  print *,'Avg. number of data points / grid cell y dir = ',float(ny_topo)/float(ny_o)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !    Create data structure with corner points for each cell
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !---- data structure to hold corners of each model cell
  allocate(OCN_CORNER_LAT(4,nx_o,ny_o),OCN_CORNER_LON(4,nx_o,ny_o),stat=ierr)

  !---- corners are numbered 1-4 counterclockwise from bottom left
  OCN_CORNER_LON(1,:,:) = X_SG(1::2,1::2)
  OCN_CORNER_LAT(1,:,:) = Y_SG(1::2,1::2)
  OCN_CORNER_LON(2,:,:) = X_SG(3::2,1::2)
  OCN_CORNER_LAT(2,:,:) = Y_SG(3::2,1::2)
  OCN_CORNER_LON(3,:,:) = X_SG(3::2,3::2)
  OCN_CORNER_LAT(3,:,:) = Y_SG(3::2,3::2)
  OCN_CORNER_LON(4,:,:) = X_SG(1::2,3::2)
  OCN_CORNER_LAT(4,:,:) = Y_SG(1::2,3::2)

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

  if  (verbose ) then
     print *,' after cell consistency check :'
     print *,' OCN_CORNER LON min/max = ',minval(OCN_CORNER_LON),&
          &                               maxval(OCN_CORNER_LON)
     print *,' OCN_CORNER LAT min/max = ',minval(OCN_CORNER_LAT),&
          &                               maxval(OCN_CORNER_LAT)
  endif

  !-----------------------------------------------------------------------
  !
  !     loop through each ocean point to compute topo statistics at each
  !     ocean grid cell
  !
  !-----------------------------------------------------------------------

  !---- allocate the output arrays
  allocate(OCN_MASK(nx_o,ny_o),OCN_FRAC(nx_o,ny_o),&
       &   AVG_TOPO(nx_o,ny_o),&
       &   OCN_D_MEAN(nx_o,ny_o), OCN_D2_MEAN(nx_o,ny_o), &
       &   OCN_D_MIN(nx_o,ny_o),  OCN_D_MAX(nx_o,ny_o),&
       &    stat=ierr)
  if ( ierr /= 0 ) then
     print *,' Could not allocate output arrays. ierr=',ierr
     stop
  end if
  if ( do_median ) then
     allocate(OCN_D_MEDIAN(nx_o,ny_o),stat=ierr)
     if ( ierr /= 0 ) then
        print *,' Could not allocate output arrays. ierr=',ierr
        stop
     end if
  endif

  !---- allocate a buffer to hold sub-sample values
  allocate(rbuffer(nx_sub*ny_sub),stat=ierr)
  if ( ierr /= 0 ) then
     print *,' Could not allocate rbuffer work array. ierr=',ierr
     stop
  end if

  !---- Make sure everything is zeroed
  OCN_MASK = 0
  OCN_FRAC = 0.
  AVG_TOPO = 0.
  OCN_D_MEAN = 0.
  if ( do_median ) OCN_D_MEDIAN = 0.
  OCN_D_MIN = 0.
  OCN_D_MAX = 0.
  OCN_D2_MEAN = 0.

  if ( verbose .AND. verbose ) &
       & print *,' Before threaded loop. max_threads=',omp_get_max_threads()

  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(DYNAMIC,8) &
  !$OMP& PRIVATE(ocn_j,ocn_i,isub,jsub,ifrac,jfrac,rlat,rlon, &
  !$OMP&         topo_i,topo_j,nsum,nsumo,frac,rbuffer)
  do ocn_j=1,ny_o

     if ( ocn_j == 1 .AND. verbose ) &
          &  print *,'Beginning threaded loop. num_threads=',omp_get_num_threads()

     if ( mod(ocn_j-1,jprnt) == 0 .AND. verbose ) &
          &   print *,'Doing model j = ',ocn_j,' thread=',OMP_get_thread_num(),&
          &    'OCN_CORNER_LON(:,1,',ocn_j,')=',OCN_CORNER_LON(:,1,ocn_j),&
          &    'OCN_CORNER_LAT(:,1,',ocn_j,')=',OCN_CORNER_LAT(:,1,ocn_j)

     do ocn_i=1,nx_o

        !*** distribute sub points
        nsum = 0
        nsumo = 0
        rbuffer = 0.d0
        do jsub=1,ny_sub
           jfrac = real(jsub)/real(ny_sub+1)
           do isub=1,nx_sub
              ifrac = real(isub)/real(nx_sub+1)

              rlat = (c1-ifrac)*(c1-jfrac) *OCN_CORNER_LAT(1,ocn_i,ocn_j)         &
                   &           +     ifrac *(c1-jfrac)*OCN_CORNER_LAT(2,ocn_i,ocn_j)  &
                   &           +     ifrac *    jfrac *OCN_CORNER_LAT(3,ocn_i,ocn_j)  &
                   &           + (c1-ifrac)*    jfrac *OCN_CORNER_LAT(4,ocn_i,ocn_j)

              rlon = (c1-ifrac)*(c1-jfrac) *OCN_CORNER_LON(1,ocn_i,ocn_j)         &
                   &           +     ifrac *(c1-jfrac)*OCN_CORNER_LON(2,ocn_i,ocn_j)  &
                   &           +     ifrac *    jfrac *OCN_CORNER_LON(3,ocn_i,ocn_j)  &
                   &           + (c1-ifrac)*    jfrac *OCN_CORNER_LON(4,ocn_i,ocn_j)

              ! Remap longitudes into -pi to pi
              do while ( rlon < TOPO_LON(0) )
                 rlon = rlon + pi2
              enddo
              do while ( rlon > TOPO_LON(nx_topo+1)) 
                 rlon = rlon - pi2
              enddo

              if ( rlat < TOPO_LAT(0) .OR. rlat > TOPO_LAT(ny_topo+1)) then
                 print *,' Invalid search lat/lon'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lon = ',OCN_CORNER_LON(:,ocn_i,ocn_j)
                 print *,' corner lat = ',OCN_CORNER_LAT(:,ocn_i,ocn_j)
                 print *,' rlon = ',rlon
                 print *,' rlat = ',rlat
                 stop
              endif

              !*** find location (nearest neighbor) of this point on hi-res topo grid
!!$              topo_i = floor( (rlon - TOPO_LON(0))/dlon_topo )+1
!!$              topo_j = floor( (rlat - TOPO_LAT(0))/dlat_topo )+1
              topo_i = nint( (rlon - TOPO_LON(0))/dlon_topo )
              topo_j = nint( (rlat - TOPO_LAT(0))/dlat_topo )

              if ( topo_i < 0 .OR. topo_i > nx_topo+1) then
                 print *,' Invalid search longitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lon = ',OCN_CORNER_LON(:,ocn_i,ocn_j)
                 print *,' corner lat = ',OCN_CORNER_LAT(:,ocn_i,ocn_j)
                 print *,' rlon = ',rlon
                 print *,' topo_i= ',topo_i
                 stop
              endif

              if ( topo_j < 0 .OR. topo_j > ny_topo+1) then
                 print *,' Invalid search latitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lat = ',OCN_CORNER_LON(:,ocn_i,ocn_j)
                 print *,' corner lat = ',OCN_CORNER_LAT(:,ocn_i,ocn_j)
                 print *,' rlat = ',rlat
                 print *,' topo_j = ',topo_j
                 stop
              endif

              !---- Compute average elevation (land + ocean) for this cell
              nsum = nsum + 1
              AVG_TOPO(ocn_i,ocn_j) = AVG_TOPO(ocn_i,ocn_j) + TOPO_Z(topo_i,topo_j)

              !---- gather subcell ocean depths for statistics
              if (-TOPO_Z(topo_i,topo_j) > mask_hmin) then
                 nsumo = nsumo + 1
                 rbuffer(nsumo) = -TOPO_Z(topo_i,topo_j)
              endif

           enddo   ! isub
        enddo   ! jsub

        if (nsum /= 0) then

           rsum = 1.0d0/dble(nsum)
           AVG_TOPO(ocn_i,ocn_j) = AVG_TOPO(ocn_i,ocn_j)*rsum

           !--- if more than given fraction of the sub-grid points are ocean
           !    call this model point ocean
           OCN_FRAC(ocn_i,ocn_j) = dble(nsumo)*rsum
           if (OCN_FRAC(ocn_i,ocn_j) >= mask_threshold) OCN_MASK(ocn_i,ocn_j) = 1

           !---- Compute depth statistics
           if ( nsumo > 0 ) then
              rsum = 1.0d0/dble(nsumo)
              OCN_D_MEAN(ocn_i,ocn_j) = sum(rbuffer(1:nsumo))*rsum
              OCN_D2_MEAN(ocn_i,ocn_j) = sum(rbuffer(1:nsumo)**2)*rsum
              OCN_D_MIN(ocn_i,ocn_j) = minval(rbuffer(1:nsumo))
              OCN_D_MAX(ocn_i,ocn_j) = maxval(rbuffer(1:nsumo))

              if ( do_median ) then
                 call quicksort(rbuffer(1:nsumo))
                 if ( mod(nsumo,2) == 0 ) then
                    OCN_D_MEDIAN(ocn_i,ocn_j) = 0.5*(rbuffer(nsumo/2) + rbuffer(nsumo/2+1))
                 else
                    OCN_D_MEDIAN(ocn_i,ocn_j) = rbuffer(nsumo/2)
                 endif
              endif
           endif


        endif

        if(mod(ocn_j-1,jprnt) == 0 .AND. ocn_i == ichk .AND. verbose ) then
           print *,'i,j=',ichk,ocn_j,' # points =',nsum,' # ocn points =',nsumo
           print *,'i,j=',ichk,ocn_j,' OCN_FRAC=',OCN_FRAC(ichk,ocn_j)
           print *,'i,j=',ichk,ocn_j,' OCN_D_MEAN=',OCN_D_MEAN(ichk,ocn_j)
           if ( do_median ) print *,'i,j=',ichk,ocn_j,' OCN_D_MEDIAN=',OCN_D_MEDIAN(ichk,ocn_j)
           print *,'i,j=',ichk,ocn_j,' OCN_D_MIN=',OCN_D_MIN(ichk,ocn_j)
           print *,'i,j=',ichk,ocn_j,' OCN_D_MAX=',OCN_D_MAX(ichk,ocn_j)
        endif

     enddo ! ocn_i


  enddo   ! ocn_j
  !$OMP END PARALLEL DO


  print *,' done finding OCN_MASK'

  nsum = nx_o*ny_o

  print *,'AVG_TOPO min,max,mean ',minval(AVG_TOPO),maxval(AVG_TOPO),&
       &   sum(AVG_TOPO)/real(nsum)
  print *,'OCN_MASK min,max,mean ',minval(OCN_MASK),maxval(OCN_MASK),&
       &   sum(OCN_MASK)/real(nsum)
  print *,'OCN_FRAC min,max,mean ',minval(OCN_FRAC),maxval(OCN_FRAC),&
       &   sum(OCN_FRAC)/real(nsum)
  print *,'OCN_D_MEAN min,max,mean ',minval(OCN_D_MEAN),maxval(OCN_D_MEAN),&
       &   sum(OCN_D_MEAN)/real(nsum)
  if ( do_median ) print *,'OCN_D_MEDIAN min,max,mean ',minval(OCN_D_MEDIAN),maxval(OCN_D_MEDIAN),&
       &   sum(OCN_D_MEDIAN)/real(nsum)
  print *,'OCN_D2_MEAN min,max,mean ',minval(OCN_D2_MEAN),maxval(OCN_D2_MEAN),&
       &   sum(OCN_D2_MEAN)/real(nsum)
  print *,'OCN_D_MIN min,max,mean ',minval(OCN_D_MIN),maxval(OCN_D_MIN),&
       &   sum(OCN_D_MIN)/real(nsum)
  print *,'OCN_D_MAX min,max,mean ',minval(OCN_D_MAX),maxval(OCN_D_MAX),&
       &   sum(OCN_D_MAX)/real(nsum)


  !-----------------------------------------------------------------------
  !
  !     write Output fields to file
  !
  !-----------------------------------------------------------------------

  ierr = nf_create(trim(output_file),NF_64bit_OFFSET,ncid)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_create')
     stop
  endif

  !---- Global attributes
  attr_name = 'Description'
  attr_valc = 'Ocean Topography Statistics on MOM6 Grid'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Creator'
  attr_valc = user_name
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  call date_and_time(cbuffer)
  attr_name = 'Created'
  attr_valc = trim(cbuffer)
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Generating Code'
  attr_valc = 'create_model_topo.f90'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Model Grid Version'
  attr_valc = mom6_grid_version
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Source Topography Data'
  attr_valc = topo_in_file
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  !---- Define dimensions
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

  !---- Define 1D nominal coordinate variables
  ndim = 1
  dims = (/dimid_lonh,1/)
  ierr = nf_def_var(ncid,'lonh',NF_DOUBLE,ndim,dims,vid_lonh)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var lonh')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal longitude of cell center points'
  ierr = nf_put_att_text(ncid,vid_lonh,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_lonh,attr_name,len_trim(attr_valc),trim(attr_valc))

  dims = (/dimid_lath,1/)
  ierr = nf_def_var(ncid,'lath',NF_DOUBLE,ndim,dims,vid_lath)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var lath')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal latitude of cell center points'
  ierr = nf_put_att_text(ncid,vid_lath,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_lath,attr_name,len_trim(attr_valc),trim(attr_valc))

  dims = (/dimid_lonq,1/)
  ierr = nf_def_var(ncid,'lonq',NF_DOUBLE,ndim,dims,vid_lonq)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var lonq')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal longitude of cell corner points'
  ierr = nf_put_att_text(ncid,vid_lonq,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_lonq,attr_name,len_trim(attr_valc),trim(attr_valc))

  dims = (/dimid_latq,1/)
  ierr = nf_def_var(ncid,'latq',NF_DOUBLE,ndim,dims,vid_latq)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var latq')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal latitude of cell corner points'
  ierr = nf_put_att_text(ncid,vid_latq,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_latq,attr_name,len_trim(attr_valc),trim(attr_valc))


  !---- Define 2D coordinate variables
  ndim = 2
  dims = (/dimid_lonh,dimid_lath/)
  ierr = nf_def_var(ncid,'geolon',NF_DOUBLE,ndim,dims,vid_geolon)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolon')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Longitude of cell center points'
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
  attr_valc = 'Latitude of cell center points'
  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_north'
  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))

  ndim = 2
  dims = (/dimid_lonq,dimid_latq/)
  ierr = nf_def_var(ncid,'geolonb',NF_DOUBLE,ndim,dims,vid_geolonb)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var geolon')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Longitude of cell corner points'
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
  attr_valc = 'Latitude of cell corner points'
  ierr = nf_put_att_text(ncid,vid_geolatb,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_north'
  ierr = nf_put_att_text(ncid,vid_geolatb,attr_name,len_trim(attr_valc),trim(attr_valc))


  !---- Define field variables
  ndim = 2
  dims = (/dimid_lonh,dimid_lath/)

  ierr = nf_def_var(ncid,'z',NF_REAL,ndim,dims,vid_z)
  attr_name = 'longname'
  attr_valc = 'Average Elevation (land and ocean) in cell'
  ierr = nf_put_att_text(ncid,vid_z,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'm'
  ierr = nf_put_att_text(ncid,vid_z,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'ocn_frac',NF_REAL,ndim,dims,vid_ofrac)
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
  attr_name = 'Ocean area threshold for mask=ocean'
  ierr = nf_put_att_real(ncid,vid_mask,len_trim(attr_name),1,mask_threshold)
  attr_name = 'Depth min to contribute to ocean (m)'
  ierr = nf_put_att_real(ncid,vid_mask,len_trim(attr_name),1,mask_hmin)

  ierr = nf_def_var(ncid,'D_mean',NF_REAL,ndim,dims,vid_d_mean)
  attr_name = 'longname'
  attr_valc = 'mean ocean depth in cell'
  ierr = nf_put_att_text(ncid,vid_d_mean,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'meters'
  ierr = nf_put_att_text(ncid,vid_d_mean,attr_name,len_trim(attr_valc),trim(attr_valc))

  if ( do_median ) then
     ierr = nf_def_var(ncid,'D_median',NF_REAL,ndim,dims,vid_d_median)
     attr_name = 'longname'
     attr_valc = 'median ocean depth in cell'
     ierr = nf_put_att_text(ncid,vid_d_median,attr_name,len_trim(attr_valc),trim(attr_valc))
     attr_name = 'units'
     attr_valc = 'm'
     ierr = nf_put_att_text(ncid,vid_d_median,attr_name,len_trim(attr_valc),trim(attr_valc))
  endif

  ierr = nf_def_var(ncid,'D2_mean',NF_REAL,ndim,dims,vid_d2_mean)
  attr_name = 'longname'
  attr_valc = 'mean squared ocean depth in cell'
  ierr = nf_put_att_text(ncid,vid_d2_mean,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'meters^2'
  ierr = nf_put_att_text(ncid,vid_d2_mean,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'D_min',NF_REAL,ndim,dims,vid_d_min)
  attr_name = 'longname'
  attr_valc = 'minimum ocean depth in cell'
  ierr = nf_put_att_text(ncid,vid_d_min,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'm'
  ierr = nf_put_att_text(ncid,vid_d_min,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'D_max',NF_REAL,ndim,dims,vid_d_max)
  attr_name = 'longname'
  attr_valc = 'maximum ocean depth in cell'
  ierr = nf_put_att_text(ncid,vid_d_max,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'm'
  ierr = nf_put_att_text(ncid,vid_d_max,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_enddef(ncid)

  ierr = nf_put_vara_double(ncid,vid_lonh,1,nx_o,lonh)
  ierr = nf_put_vara_double(ncid,vid_lonq,1,nx_o+1,lonq)
  ierr = nf_put_vara_double(ncid,vid_lath,1,ny_o,lath)
  ierr = nf_put_vara_double(ncid,vid_latq,1,ny_o+1,latq)

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_double(ncid,vid_geolon,indstr,indcnt,geolon)
  ierr = nf_put_vara_double(ncid,vid_geolat,indstr,indcnt,geolat)

  indstr = (/1,1/)
  indcnt = (/nx_o+1,ny_o+1/)
  ierr = nf_put_vara_double(ncid,vid_geolonb,indstr,indcnt,geolonb)
  ierr = nf_put_vara_double(ncid,vid_geolatb,indstr,indcnt,geolatb)

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_real(ncid,vid_z,indstr,indcnt,AVG_TOPO)
  ierr = nf_put_vara_int(ncid,vid_mask,indstr,indcnt,OCN_MASK)
  ierr = nf_put_vara_real(ncid,vid_ofrac,indstr,indcnt,OCN_FRAC)
  ierr = nf_put_vara_real(ncid,vid_d_mean,indstr,indcnt,OCN_D_mean)
  if ( do_median )   ierr = nf_put_vara_real(ncid,vid_d_median,indstr,indcnt,OCN_D_median)
  ierr = nf_put_vara_real(ncid,vid_d2_mean,indstr,indcnt,OCN_D2_mean)
  ierr = nf_put_vara_real(ncid,vid_d_min,indstr,indcnt,OCN_D_min)
  ierr = nf_put_vara_real(ncid,vid_d_max,indstr,indcnt,OCN_D_max)


  ierr = nf_close(ncid)

  !-----------------------------------------------------------------------

end program create_model_topo

SUBROUTINE  Sort(x, N)
  use kinds

  IMPLICIT  NONE
  real(kind=real_kind), DIMENSION(1:), INTENT(INOUT) :: x
  INTEGER, INTENT(IN)                   :: N
  INTEGER                               :: i
  INTEGER,dimension(1)                  :: imin
  real(kind=real_kind) :: temp

  DO i = 1, N-1             ! except for the last
     imin = minloc(x(i:N))  ! find min from this to last

     temp = x(i)             ! swap this and min
     x(i) = x(imin(1)+i-1)
     x(imin(1)+i-1) = temp
  END DO
END SUBROUTINE  Sort
! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a)
  use kinds
  implicit none
  real(kind=dbl_kind) :: a(:)
  real x, t
  integer :: first = 1, last
  integer i, j

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last))
end subroutine quicksort

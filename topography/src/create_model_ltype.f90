
!***********************************************************************

program create_model_ltype

  !-----------------------------------------------------------------------
  !
  !     This program creates MASK field for surface type for a MOM6 ocean grid by 
  !     interpolating/averaging from a MODIS derived dataset
  !
  !     This version uses a Monte-Carlo-like method 
  !     (based on a code from Mat Maltrud used for POP)
  ! 
  !     Compile with gnu
  !-----------------------------------------------------------------------

  use kinds
  use constants
  use ncdf_wrapper
  use mom6_grid
  use omp_lib

  implicit none

  !-----------------------------------------------------------------------

  !---- Input topography data set
  integer (kind=int_kind) :: nx_mask, ny_mask
  real(kind=dbl_kind) :: dlon_mask, dlat_mask
  real (kind=real_kind), dimension(:), allocatable :: MASK_LON, MASK_LAT
  integer(kind=int_kind), dimension(:,:), allocatable :: INDAT, MASK_Z ! input data

  !---- Results on model grid
  integer, parameter :: ntype=8
  integer, dimension(ntype) :: cnt_type
  real (kind=real_kind), dimension(:,:,:), allocatable :: MASK_FRAC
  real (kind=real_kind), dimension(:,:), allocatable :: OCN_FRAC

  !---- parameters
  integer(kind=int_kind) ::&
       & nx_sub, ny_sub, &
       & ichk=1, jprnt=1000

  logical :: verbose=.false., do_median=.false.

  character (char_len) :: &
       &  mom6_grid_version     &
       &, mom6_horiz_grid_file   &    ! file with horizontal ocean grid
       &, mask_in_file    &     ! file with input high-res topo data
       &, output_file    &    ! file for results
       &, user_name
  character (char_len) :: vname_mask_lon, vname_mask_lat, vname_mask_z

  logical :: grid_is_tripole=.true.

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !    Work variables
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  real (kind=dbl_kind), dimension(:,:,:), allocatable :: OCN_CORNER_LAT, OCN_CORNER_LON
  real(kind=real_kind), dimension(:), allocatable :: rbuffer
  character (char_len) :: cbuffer

  !---- netCDF stuff
  integer :: ierr, ncid, dimid, vid, ndim
  integer :: dimid_lonh, dimid_lath, dimid_lonq, dimid_latq 
  integer, dimension(3) :: indcnt, indstr, dims
  character(char_len) :: dim_name, var_name, attr_name, attr_valc
  integer :: vid_topo, vid_ofrac, vid_mask, &
       &     vid_lonh, vid_lath, vid_lonq, vid_latq, &
       &     vid_geolon, vid_geolat, vid_geolonb, vid_geolatb, &
       &     vid_z

  integer (kind=int_kind) :: &
       &  i,j,k,im1,jm1 &
       &, ocn_i, ocn_j    &
       &, topo_i, topo_j  &
       &, isub, jsub  &
       &, nsum, nsumo &
       &, n

  real (kind=dbl_kind) :: &
       &  ifrac, jfrac         &
       &, rlat, rlon           &
       &, dlon, dlat     &
       &, frac

  namelist /model_grid_nml/nx_sub, ny_sub, mom6_grid_version, mom6_horiz_grid_file, grid_is_tripole

  namelist /ltype_in_nml/ mask_in_file, vname_mask_lon, vname_mask_lat, vname_mask_z

  namelist /output_nml/ verbose, jprnt, ichk, &
       &                user_name, output_file

  !-----------------------------------------------------------------------
  !     Initialize and get namelist input
  !-----------------------------------------------------------------------

  call init_constants

  nx_sub = 1
  ny_sub = 1
  mom6_grid_version = 'unknown-grid-version'
  mom6_horiz_grid_file = 'unknown-horiz-grid'
  mask_in_file    = 'unknown-topo-data-file'
  vname_mask_lon = 'lon'
  vname_mask_lat = 'lat'
  vname_mask_z = 'z'
  user_name = 'unknown-user'
  output_file = 'unknown-output-file'

  read(*,nml=model_grid_nml)
  write(*,nml=model_grid_nml)

  read(*,nml=ltype_in_nml)
  write(*,nml=ltype_in_nml)

  read(*,nml=output_nml)
  write(*,nml=output_nml)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Get model super grid and allocate model grid arrays
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  call read_super_grid(mom6_horiz_grid_file)
  call compute_mom6_grid()

  !-----------------------------------------------------------------------
  !     read hi-res mask data
  !-----------------------------------------------------------------------

  print *,'Opening ',mask_in_file
  ierr = nf_open(mask_in_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr = ',ierr,'could not open ',mask_in_file
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_mask_lon,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_mask_lon),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_mask)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_mask_lon),' ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dimid(ncid,vname_mask_lat,dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_mask_lat),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_mask)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_mask_lat),' ierr= ',ierr
     stop
  endif
  print *,'Input topo dimension = ',nx_mask,ny_mask

  allocate(MASK_LON(0:nx_mask+1))

  ierr = nf_inq_varid(ncid,vname_mask_lon,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',trim(vname_mask_lon),' ierr= ',ierr
     stop
  endif

  indstr(1) = 1
  indcnt(1) = nx_mask
  ierr = nf_get_vara_real(ncid,vid,indstr,indcnt,MASK_LON(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',trim(vname_mask_lon),' ierr= ',ierr
     stop
  endif

  MASK_LON(0) = MASK_LON(nx_mask)-360.
  MASK_LON(nx_mask+1) = MASK_LON(1)+360.
  dlon_mask = MASK_LON(2)-MASK_LON(1)

  if ( verbose ) then
     write(*,"(a,2f16.7)") 'Extended TOPO LON min/max=',minval(MASK_LON),maxval(MASK_LON)
     write(*,"(a,f20.10)") 'dlon_mask = ',dlon_mask,' degrees'
     write(*,"(a,4f9.4,a,4f9.4)") &
          & 'LON = ',MASK_LON(0:3),' ... ',MASK_LON(nx_mask-2:nx_mask+1)
     print *
  endif

  allocate(MASK_LAT(0:ny_mask+1))

  ierr = nf_inq_varid(ncid,vname_mask_lat,vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',vname_mask_lat(1:len_trim(vname_mask_lat)),' ierr= ',ierr
     stop
  endif

  indstr(1) = 1
  indcnt(1) = ny_mask
  ierr = nf_get_vara_real(ncid,vid,indstr,indcnt,MASK_LAT(1))
  if ( ierr /= 0 ) then
     print *,'could not get variable ',trim(vname_mask_lat),' ierr= ',ierr
     stop
  endif

  dlat_mask = MASK_LAT(2)-MASK_LAT(1)
  MASK_LAT(0) = MASK_LAT(1) - dlat_mask
  MASK_LAT(ny_mask+1) = MASK_LAT(ny_mask) + dlat_mask

  if ( verbose ) then
     write(*,"(a,2f16.7)") 'Extended TOPO LAT min/max=',minval(MASK_LAT),maxval(MASK_LAT)
     write(*,"(a,f20.10)")' dlat = ',dlat_mask,' degrees'
     write(*,"(a,4f9.4,a,4f9.4)") &
          & 'LAT = ',MASK_LAT(0:3),' ... ',MASK_LAT(ny_mask-2:ny_mask+1)
     print *
  endif

  allocate(INDAT(nx_mask,ny_mask),MASK_Z(0:nx_mask+1,0:ny_mask+1),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'allocate stat=',ierr
     stop
  end if

  ierr = nf_inq_varid(ncid,vname_mask_z,vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find topo height variable ',vname_mask_z
     stop
  endif
  call flush()

  indstr(1:2) = (/1,1/)
  indcnt(1:2) = (/nx_mask,ny_mask/)
  ierr = nf_get_vara_int(ncid,vid,indstr,indcnt,INDAT)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get mask variable ',vname_mask_z
     stop
  endif
  print *,'Read INDAT'
  call flush()

  MASK_Z(1:nx_mask,1:ny_mask) = INDAT

  !---- Extend the intput grid one point in each direction
  !     so all possible lat-lons are included
  MASK_Z(0,1:ny_mask) = INDAT(nx_mask,:)
  MASK_Z(nx_mask,1:ny_mask) = INDAT(1,:)
  MASK_Z(:,0) = MASK_Z(:,1)
  MASK_Z(:,ny_mask+1) = MASK_Z(:,ny_mask)

  print *,' MASK_LON min=',minval(MASK_LON),' max=',maxval(MASK_LON)
  print *,' MASK_LAT min=',minval(MASK_LAT),' max=',maxval(MASK_LAT)
  print *,' MASK_Z min=',minval(MASK_Z),' max=',maxval(MASK_Z)

  !---- Convert input data coordinates to radians
  MASK_LON = MASK_LON/radian
  dlon_mask = dlon_mask/radian
  MASK_LAT = MASK_LAT/radian
  dlat_mask = dlat_mask/radian

  print *,'Number of data points / grid cell x dir = ',nx_mask/nx_o
  print *,'Number of data points / grid cell y dir = ',ny_mask/ny_o
  call flush()

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
  allocate(MASK_FRAC(nx_o,ny_o,ntype),OCN_FRAC(nx_o,ny_o),stat=ierr)
  if ( ierr /= 0 ) then
     print *,' Could not allocate output arrays. ierr=',ierr
     stop
  end if


  !---- Make sure everything is zeroed
  MASK_FRAC = 0

  if ( verbose ) print *,' Before threaded loop. max_threads=',omp_get_max_threads()

  !$OMP PARALLEL DO &
  !$OMP SCHEDULE(DYNAMIC,8) &
  !$OMP& PRIVATE(ocn_j,ocn_i,isub,jsub,ifrac,jfrac,rlat,rlon, &
  !$OMP&         topo_i,topo_j,nsum,nsumo,frac,rbuffer)
  do ocn_j=1,ny_o

     if ( ocn_j == 1 ) then
        print *,'Beginning threaded loop. num_threads=',omp_get_num_threads()
     endif

     if ( mod(ocn_j-1,jprnt) == 0 .AND. verbose ) then
        print *,'Doing model j = ',ocn_j,' thread=',OMP_get_thread_num(),&
             &    'OCN_CORNER_LON(:,1,',ocn_j,')=',OCN_CORNER_LON(:,1,ocn_j),&
             &    'OCN_CORNER_LAT(:,1,',ocn_j,')=',OCN_CORNER_LAT(:,1,ocn_j)
     endif

     do ocn_i=1,nx_o

        !*** distribute sub points
        nsum = 0
        nsumo = 0
        cnt_type = 0
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
              do while ( rlon < MASK_LON(0) )
                 rlon = rlon + pi2
              enddo
              do while ( rlon > MASK_LON(nx_mask+1)) 
                 rlon = rlon - pi2
              enddo

              if ( rlat < MASK_LAT(0) .OR. rlat > MASK_LAT(ny_mask+1)) then
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
              topo_i = floor( (rlon - MASK_LON(0))/dlon_mask )+1
              topo_j = floor( (rlat - MASK_LAT(0))/dlat_mask )+1

              if ( topo_i < 0 .OR. topo_i > nx_mask+1) then
                 print *,' Invalid search longitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lon = ',OCN_CORNER_LON(:,ocn_i,ocn_j)
                 print *,' corner lat = ',OCN_CORNER_LAT(:,ocn_i,ocn_j)
                 print *,' rlon = ',rlon
                 print *,' topo_i= ',topo_i
                 stop
              endif

              if ( topo_j < 0 .OR. topo_j > ny_mask+1) then
                 print *,' Invalid search latitude'
                 print *,' model i,j, = ',ocn_i, ocn_j
                 print *,' subpoint i,j, = ',isub, jsub
                 print *,' corner lat = ',OCN_CORNER_LON(:,ocn_i,ocn_j)
                 print *,' corner lat = ',OCN_CORNER_LAT(:,ocn_i,ocn_j)
                 print *,' rlat = ',rlat
                 print *,' topo_j = ',topo_j
                 stop
              endif

              !---- Compute fraction of each type for this cell
              nsum = nsum + 1
              n = MASK_Z(topo_i,topo_j)+1  !! zero base
              cnt_type(n) = cnt_type(n) + 1

           enddo   ! isub
        enddo   ! jsub

        if (nsum /= 0) then

           do n=1,ntype
              MASK_FRAC(ocn_i,ocn_j,n) = cnt_type(n)/float(nsumo)
           enddo
           OCN_FRAC(ocn_i,ocn_j) = MASK_FRAC(ocn_i,ocn_j,1) + MASK_FRAC(ocn_i,ocn_j,7) + &
                &                  MASK_FRAC(ocn_i,ocn_j,8)
        endif

        if(mod(ocn_j-1,jprnt) == 0 .AND. ocn_i == ichk .AND. verbose ) then
           print *,'j=',j,' # points =',nsum
           print *,'j=',j,' MASK_FRAC(ichk,j)=',MASK_FRAC(ichk,ocn_j,:)
           print *,'j=',j,' OCN_FRAC(ichk,j)=',OCN_FRAC(ichk,ocn_j)
        endif

     enddo ! ocn_i


  enddo   ! ocn_j
  !$OMP END PARALLEL DO


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
  attr_valc = 'Land Water Mask on MOM6 Grid'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Creator'
  attr_valc = user_name
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  call date_and_time(cbuffer)
  attr_name = 'Created'
  attr_valc = trim(cbuffer)
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Generating Code'
  attr_valc = 'create_model_ltype.f90'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Model Grid Version'
  attr_valc = mom6_grid_version
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Source Topography Data'
  attr_valc = mask_in_file
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
  dims(1) = dimid_lonh
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

  dims(1) = dimid_lath
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

  dims(1) = dimid_lonq
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

  dims(1) = dimid_latq
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


  dims(1) = ntype
  ierr = nf_def_var(ncid,'stype',NF_INT,ndim,dims,vid_z)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var type')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Surface Type'
  ierr = nf_put_att_text(ncid,vid_z,attr_name,len_trim(attr_valc),trim(attr_valc))


  !---- Define 2D coordinate variables
  ndim = 2
  dims(1:2) = (/dimid_lonh,dimid_lath/)
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
  dims(1:2) = (/dimid_lonq,dimid_latq/)
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
  dims(1:2) = (/dimid_lonh,dimid_lath/)
  ierr = nf_def_var(ncid,'ocn_frac',NF_REAL,ndim,dims,vid_ofrac)
  attr_name = 'longname'
  attr_valc = 'fraction of cell area covered by ocean'
  ierr = nf_put_att_text(ncid,vid_ofrac,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'dimensionless'
  ierr = nf_put_att_text(ncid,vid_ofrac,attr_name,len_trim(attr_valc),trim(attr_valc))

  ndim = 3
  dims = (/dimid_lonh,dimid_lath,ntype/)
  ierr = nf_def_var(ncid,'LandWaterFraction',NF_REAL,ndim,dims,vid_z)
  attr_name = 'longname'
  attr_valc = 'fraction of Each Land/Water Type in cell'
  ierr = nf_put_att_text(ncid,vid_z,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'dimensionless'
  ierr = nf_put_att_text(ncid,vid_z,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_enddef(ncid)

  ierr = nf_put_vara_double(ncid,vid_lonh,1,nx_o,lonh)
  ierr = nf_put_vara_double(ncid,vid_lonq,1,nx_o+1,lonq)
  ierr = nf_put_vara_double(ncid,vid_lath,1,ny_o,lath)
  ierr = nf_put_vara_double(ncid,vid_latq,1,ny_o+1,latq)

  indstr(1:2) = (/1,1/)
  indcnt(1:2) = (/nx_o,ny_o/)
  ierr = nf_put_vara_double(ncid,vid_geolon,indstr,indcnt,geolon)
  ierr = nf_put_vara_double(ncid,vid_geolat,indstr,indcnt,geolat)

  indstr(1:2) = (/1,1/)
  indcnt(1:2) = (/nx_o+1,ny_o+1/)
  ierr = nf_put_vara_double(ncid,vid_geolonb,indstr,indcnt,geolonb)
  ierr = nf_put_vara_double(ncid,vid_geolatb,indstr,indcnt,geolatb)

  indstr(1:2) = (/1,1/)
  indcnt(1:2) = (/nx_o,ny_o/)
  ierr = nf_put_vara_real(ncid,vid_ofrac,indstr,indcnt,OCN_FRAC)

  indstr = (/1,1,1/)
  indcnt = (/nx_o,ny_o,ntype/)
  ierr = nf_put_vara_real(ncid,vid_z,indstr,indcnt,MASK_FRAC)

  ierr = nf_close(ncid)

  !-----------------------------------------------------------------------

end program create_model_ltype

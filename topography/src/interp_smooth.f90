!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Find the depth on MOM6 grid using a Cressman type distance weighted
!    smoothing interpolation. Expects a pre-computed land/ocean mask as
!    input.
!    Observed topography data is expected on a uniform lat-lon grid.
!    Typically used with ETOPO or SRTM data.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

program interp_smooth

  use kinds
  use constants
  use ncdf_wrapper
  use mom6_grid
  use omp_lib

  implicit none

  !-----------------------------------------------------------------------

  real(kind=dbl_kind), external :: SphericalDistance

  integer(kind=int_kind)  :: &
       & nx_lmask, ny_lmask

  integer (kind=int_kind)  :: &
       & nx_topo,               & ! zonal dimension of input topo data
       & ny_topo                   ! merid.dimension of ET topo

  integer (kind=int_kind), dimension(:,:), allocatable :: &
       &  MASK              ! land-ocean mask for ocean grid

  real (kind=real_kind), dimension(:,:), allocatable :: &
       &  MODEL_DEPTH              ! DEPTH field for ocean grid
  integer (kind=int_kind), dimension(:,:), allocatable :: &
       &  MODEL_COUNT
  logical (kind=log_kind), dimension(:,:), allocatable :: &
       &  MODEL_FILLED              ! land-ocean mask for ocean grid

  real (kind=real_kind), dimension(:,:), allocatable :: &
       &  TOPO_DEPTH                             ! OBS. TOPO depth data
  real (kind=dbl_kind), dimension(:), allocatable :: &
       &  TOPO_LON                               ! OBS. TOPO longitudes
  real (kind=dbl_kind), dimension(:), allocatable :: &
       &  TOPO_LAT                                    ! OBS. TOPO latitudes

  integer, dimension(1) :: jmid, imid

  integer (kind=int_kind) :: &
       & i,j, iwid, ibeg, iend, jwid, jbeg, jend, &
       & ie, je, &
       & nfilled, nempty, max_iter=100, iter,&
       & nprint

  real (kind=dbl_kind) :: &
       & w, w_sum, h_sum, d, d2, L2, lon_mod, lat_mod, dlon_topo, dlat_topo

  integer(kind=int_kind) :: cnt

  real (kind=real_kind) :: hmin=0.,&
       &                   smooth_scl=1.0,&
       &                   cressman_exp = 2.0

  !---- netCDF stuff
  integer :: ierr, ncid, dimid, vid, ndim
  integer :: dimid_lonh, dimid_lath
  integer :: dimid_lonq, dimid_latq
  integer :: did_lon, did_lat
  integer :: vid_lonh, vid_lath, vid_geolon, vid_geolat
  integer :: vid_lonq, vid_latq, vid_geolonb, vid_geolatb
  integer :: vid_lmask, vid_depth
  integer, dimension(2) :: indcnt, indstr, dims
  character (char_len) :: cbuffer
  character (char_len) :: vname_lon, vname_lat, vname_lmask, vname_depth, vname_z
  character(char_len) :: dim_name, var_name, attr_name, attr_valc

  character (char_len) :: &
       &  mom6_grid_version,     &
       &  mom6_horiz_grid_file, &       ! file with horizontal ocean grid
       & topo_in_file,    &     ! file with input high-res topo data
       &  topo_data_name, &         ! label for topography data
       &  lmask_file, &           ! file with land-ocean mask
       &  user_name
  character (char_len) :: vname_topo_lon, vname_topo_lat, vname_topo_z

  character (char_len) ::  &
       &  topo_data_file, &         ! file with OBS. TOPO data
       &  model_topo_file 


  namelist /model_input_nml/ mom6_grid_version, mom6_horiz_grid_file, lmask_file, vname_lmask

  namelist /topo_data_nml/ topo_data_name, topo_data_file, vname_lon, vname_lat, vname_z

  namelist /model_output_nml/ hmin, smooth_scl, cressman_exp,&
       &    model_topo_file, vname_depth, nprint, user_name

  !-----------------------------------------------------------------------
  !
  !     read namelist input to obtain ocean grid info and
  !     allocate ocean arrays
  !
  !-----------------------------------------------------------------------

  print *,' Before threaded loop. max_threads=',omp_get_max_threads()

  hmin = 3.
  smooth_scl = 2.
  lmask_file = 'unknown-lmask-file'
  mom6_horiz_grid_file = 'unknown-hgrid-file'

  TOPO_data_file    = 'unknown-TOPO_data'
  model_topo_file      = 'unknown-depth-file'

  read (*,nml=model_input_nml)
  write(*,nml=model_input_nml)
  read (*,nml=topo_data_nml)
  write(*,nml=topo_data_nml)
  read (*,nml=model_output_nml)
  write(*,nml=model_output_nml)

  call init_constants
  call read_super_grid(mom6_horiz_grid_file)
  call compute_mom6_grid()
  call compute_mom6_grid_metrics

  allocate(MODEL_DEPTH(nx_o,ny_o),MASK(nx_o,ny_o),&
       & MODEL_FILLED(nx_o,ny_o),stat=ierr)
  if (ierr /=0 ) then
     print *,'Could not allocate model output arrays=',ierr
  endif

  print *

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Get pre-computed model land/ocean mask
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  print *,'Opening ',trim(lmask_file)
  ierr = nf_open(lmask_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr=',ierr
     print *,'could not open ',lmask_file
     stop
  endif

  ierr = nf_inq_dimid(ncid,'lonh',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension lonh. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_lmask)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  endif
  if ( nx_lmask /= nx_o) then
     print *,'ERROR x dimension of grid and lmask inconsistent:',nx_o,nx_lmask
     stop
  endif

  ierr = nf_inq_dimid(ncid,'lath',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension lath. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_lmask)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  endif
  if ( ny_lmask /= ny_o) then
     print *,'ERROR y dimension of grid and lmask inconsistent:',ny_o,ny_lmask
     stop
  endif

  ierr = nf_inq_varid(ncid,trim(vname_lmask),vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find mask variable ',trim(vname_lmask)
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_get_vara_int(ncid,vid,indstr,indcnt,MASK)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get mask variable ',vname_lmask
     stop
  endif

  ierr = nf_close(ncid)

  print*,' LMASK mix/max=',minval(MASK),maxval(MASK)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !   Get Observational depths
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  print *,'Opening ',topo_data_file
  ierr = nf_open(topo_data_file,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr = ',ierr,' ncid = ',ncid
     print *,'could not open ',topo_data_file
     stop
  endif

  ierr = nf_inq_dimid(ncid,trim(vname_lon),dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_lon),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_topo)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_lon),' ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dimid(ncid,trim(vname_lat),dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ',trim(vname_lat),' ierr= ',ierr
     stop
  endif
  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_topo)
  if ( ierr /= 0 ) then
     print *,'could not get dimension ',trim(vname_lat),' ierr= ',ierr
     stop
  endif
  print *,'Input topo dimension = ',nx_topo,ny_topo

  allocate(TOPO_LON(nx_topo),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'ERROR Could not allocate TOPO_LON. ierr=',ierr
  endif
  ierr = nf_inq_varid(ncid,trim(vname_lon),vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',trim(vname_lon),' ierr= ',ierr
     stop
  endif
  indstr(1) = 1
  indcnt(1) = nx_topo
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,TOPO_LON)
  if ( ierr /= 0 ) then
     print *,'could not get variable ',trim(vname_lon),' ierr= ',ierr
     stop
  endif
  dlon_topo = TOPO_LON(2)-TOPO_LON(1)

  print *,' TOPO LON min/max=',minval(TOPO_LON),maxval(TOPO_LON)
  print *,' LON = ',TOPO_LON(1:3),' ... ',TOPO_LON(nx_topo-2:nx_topo)
  print *,' dlon_topo = ',dlon_topo,' degrees'
  print *

  allocate(TOPO_LAT(ny_topo),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'ERROR Could not allocate TOPO_LON. ierr=',ierr
  endif
  ierr = nf_inq_varid(ncid,trim(vname_lat),vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable ',trim(vname_lat),' ierr= ',ierr
     stop
  endif
  indstr(1) = 1
  indcnt(1) = ny_topo
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,TOPO_LAT)
  if ( ierr /= 0 ) then
     print *,'could not get variable ',trim(vname_lat),' ierr= ',ierr
     stop
  endif
  dlat_topo = TOPO_LAT(2)-TOPO_LAT(1)

  print *,' TOPO LAT min/max=',minval(TOPO_LAT),maxval(TOPO_LAT)
  print *,' LAT = ',TOPO_LAT(1:3),' ... ',TOPO_LAT(ny_topo-2:ny_topo)
  print *,' dlat = ',dlat_topo,' degrees'
  print *

  allocate(TOPO_DEPTH(nx_topo,ny_topo),stat=ierr)
  if ( ierr /= 0 ) then
     print *,'ERROR Could not allocate TOPO_DEPTH ierr=',ierr
     stop
  else
     print *,'Allocated TOPO_DEPTH!!!!'
  end if

  ierr = nf_inq_varid(ncid,trim(vname_z),vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find mask variable ',vname_z
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_topo,ny_topo/)
  ierr = nf_get_vara_real(ncid,vid,indstr,indcnt,TOPO_DEPTH)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get mask variable ',vname_z
     stop
  endif

  ierr = nf_close(ncid)

  print *,' Obs TOPO min=',minval(TOPO_DEPTH),' max=',maxval(TOPO_DEPTH)
  print *
  call flush()

  !---- Flip sign of topography so depths > 0
  TOPO_DEPTH = -TOPO_DEPTH

  !-----------------------------------------------------------------------
  !
  ! Interpolate to each ocean point
  !
  !-----------------------------------------------------------------------

  print *,' Starting depth interpolation'

  !---- Find the max possible width of the patch to apply averaging kernel
  !---- Use equatorial values of model dlon
  jmid = minloc(abs(GEOLAT(1,:)))
  print *,' j equator = ',jmid,' lat=',GEOLAT(1,jmid)

  iwid = int(maxval(DLONT(:,jmid))/dlon_topo)
  jwid = int(maxval(DLATT)/dlat_topo)
  print *,'Maximum ratio model/topo grid spacing : x',iwid,' y=',jwid

  iwid = int(iwid*smooth_scl*1.25)
  jwid = int(jwid*smooth_scl*1.25)
  print *,'Search patch for averaging topo iwid=',iwid,' jwid=',jwid

  !---- Convert Obs data LAT/LON to radians
  TOPO_LAT = TOPO_LAT/radian
  TOPO_LON = TOPO_LON/radian
  dlon_topo = dlon_topo/radian
  dlat_topo = dlat_topo/radian

  MODEL_DEPTH = 0.
  MODEL_FILLED = .false.

  nempty = 0
  nfilled = 0

  !---- Loop over model grid
  !$OMP PARALLEL DO &
  !$OMP& SCHEDULE(DYNAMIC,8) &
  !$OMP& REDUCTION(+:nfilled,nempty) &
  !$OMP& PRIVATE(i,j,imid,ibeg,iend,jmid,jbeg,jend, &
  !$OMP&         lon_mod,lat_mod,L2,d,d2,w,h_sum,w_sum,cnt)
  do j=1,ny_o

     if ( j == 1 ) then
        print *,'Beginning threaded loop. num_threads=',omp_get_num_threads()
        call flush()
     endif

     if (mod(j-1,nprint) == 0 ) &
          & print *,'Doing model j = ',j,' thread=',OMP_get_thread_num()

     do i=1,nx_o

        if ( MASK(i,j) == 0 ) cycle  ! this is a land point depth=0.

        lon_mod = GEOLON(i,j)/radian
        lat_mod = GEOLAT(i,j)/radian

        !---- Map the model point into +/- pi
        do while ( lon_mod > pi) 
           lon_mod = lon_mod - pi2
        enddo
        do while ( lon_mod < -pi) 
           lon_mod = lon_mod + pi2
        enddo

        !---- Find a search patch
        jmid = minloc(abs(TOPO_LAT-lat_mod))   ! closest point, center of patch
        jbeg = max(1,jmid(1) - jwid)
        jend = min(jmid(1) + jwid,ny_topo)

        imid = minloc(abs(TOPO_LON-lon_mod))
        if ( imid(1) <= iwid .OR. &
             & imid(1) >= nx_topo - iwid .OR. &
             & jend == ny_topo) then
           !---- If close to a cyclic boundary use whole longitude
           ibeg = 1
           iend = nx_topo
        else
           ibeg = imid(1) - iwid
           iend = imid(1) + iwid
        endif

        !---- Find local smoothing radius = smooth_scl * grid spacing
        !     (only need square of it)
        L2 = smooth_scl**2*DAREA(i,j)

        !---- Integrate over the patch
        h_sum = 0.
        w_sum = 0.
        cnt = 0
        do je=jbeg,jend
           do ie=ibeg,iend

              !---- Only use points that are in the ocean (depth > 0)
              if ( TOPO_DEPTH(ie,je) <= 0. ) cycle  

              !---- Compute great circle dist model grid pt to TOPO pt
              d = SphericalDistance(TOPO_LON(ie),TOPO_LAT(je),lon_mod,lat_mod)
              d2 = (d*radius)**2

              !---- Compute interpolation weights
              if ( d2 <= L2 ) then
                 w = ((L2 - d2)/(L2 + d2) )
                 if ( cressman_exp /= 1.0 ) w = w**cressman_exp

                 !---- Sum weighted depths
                 h_sum = h_sum + w*TOPO_DEPTH(ie,je)
                 w_sum = w_sum + w
                 cnt = cnt + 1
              endif

           enddo
        enddo

        !---- Assign to model depth field
        if (h_sum == 0. .OR. w_sum == 0.) then
           print *,'No estimate for i,j=',i,j,&
                &  ' lon_mod,lat_mod=',lon_mod*radian,lat_mod*radian,&
                &  ' ibeg,jbeg=',ibeg,jbeg,&
                &  ' imid,jmid=',imid,jmid,&
                &  ' iend,jend=',iend,jend,&
                &  ' h_sum,w_sum=',h_sum,w_sum
           nempty = nempty + 1
        else
           MODEL_DEPTH(i,j) = h_sum/w_sum
           MODEL_FILLED(i,j) = .true.
           nfilled = nfilled + 1
        endif

     enddo
      if ( mod(j-1,nprint) == 0 ) then
         print *,' Thread=',OMP_get_thread_num(), &
              & ' finished j=',j,' max depth = ',maxval(MODEL_DEPTH(:,j))
          call flush()
       endif
  enddo
  !$OMP END PARALLEL DO

  print *
  if ( nempty /= 0 ) print *,'Missed ',nempty,' ocean points'
  print *

  !---- Try to fill points where interpolation failed
  nempty = count(MASK == 1 .AND. (.NOT. MODEL_FILLED))
  print *,'Need to fill in ',nempty,' points by grid averaging.'
  print *

  iter = 0
  do while (nempty > 0 .AND. iter < max_iter )
     do j=2,ny_o-1
        do i=2,ny_o-1
           if (   MASK(i,j) == 1 .AND.  (.NOT. MODEL_FILLED(i,j))) then

              h_sum = h_sum + &
                   & MODEL_DEPTH(i,j-1)*MASK(i,j-1) + &
                   & MODEL_DEPTH(i,j+1)*MASK(i,j+1) + &
                   & MODEL_DEPTH(i-1,j)*MASK(i-1,j) + &
                   & MODEL_DEPTH(i+1,j)*MASK(i+1,j)
              w_sum = w_sum + &
                   & MASK(i,j-1) + &
                   & MASK(i,j+1) + &
                   & MASK(i-1,j) + &
                   & MASK(i+1,j)

              if ( w_sum >= 1. ) then
                 MODEL_DEPTH(i,j) = h_sum/w_sum
                 nempty = nempty - 1
                 MODEL_FILLED(i,j) = .true.
              endif

           endif
        enddo
     enddo
     print *,'Finished iteration ',iter,' nempty=',nempty
     iter = iter + 1
  enddo

  nempty = count(MASK == 1 .AND. (.NOT. MODEL_FILLED))
  print *,'Setting ',nempty,' unfilled points to depth=',hmin
  where(MASK == 1 .AND. (.NOT. MODEL_FILLED)) MODEL_DEPTH = hmin

  !---- Enforce depth >= hmin
  nempty = count(MASK == 1 .AND. MODEL_DEPTH < hmin)
  print *,'Resetting ',nempty,' points to min. depth=',hmin
  where(MASK == 1 .AND. MODEL_DEPTH < hmin) MODEL_DEPTH=hmin

  !-----------------------------------------------------------------------
  !
  !  Write out result
  !
  !-----------------------------------------------------------------------

  print *,'MODEL_DEPTH min,max,mean ',minval(MODEL_DEPTH),maxval(MODEL_DEPTH),&
       &   sum(MODEL_DEPTH)/real(nx_o*ny_o)

  !---- Create the file
  ierr = nf_create(trim(model_topo_file),NF_64bit_OFFSET,ncid)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_create')
     stop
  endif

  !---- Global attributes
  attr_name = 'Description'
  attr_valc = 'Ocean Topography on MOM6 Grid'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Creator'
  attr_valc = user_name
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  call date_and_time(cbuffer)
  attr_name = 'Created'
  attr_valc = trim(cbuffer)
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Generating Code'
  attr_valc = 'interp_smooth.f90'
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Model Grid Version'
  attr_valc = mom6_grid_version
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Source Topography Data'
  attr_valc = topo_data_file
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Source Mask Data'
  attr_valc = lmask_file
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Smoothing Scale L (grid lengths)'
  ierr = nf_put_att_real(ncid,NF_GLOBAL,attr_name,NF_FLOAT,1,smooth_scl)

  attr_name = 'Smoothing Function'
  if ( cressman_exp == 1.0 ) then
     attr_valc = '(L^2 - r^2)/(L^2 + r^2)'
  else
     write(attr_valc,"(a,f3.1)") '(L^2 - r^2)/(L^2 + r^2)^', cressman_exp
  endif
  attr_valc = trim(attr_valc)
  ierr = nf_put_att_text(ncid,NF_GLOBAL,attr_name,len_trim(attr_valc),trim(attr_valc))

  attr_name = 'Minimum Depth (m)'
  ierr = nf_put_att_real(ncid,NF_GLOBAL,attr_name,NF_FLOAT,1,hmin)

  !---- Define the dimensions
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

  !---- Define nominal (1D) coordinate variables
  ierr = nf_def_var(ncid,'lonh',NF_DOUBLE,1,dimid_lonh,vid_lonh)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var lonh')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal (1D) longitude of tracer (T) points'
  ierr = nf_put_att_text(ncid,vid_lonh,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_lonh,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'lath',NF_DOUBLE,1,dimid_lath,vid_lath)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_def_var lath')
     stop
  endif
  attr_name = 'longname'
  attr_valc = 'Nominal (1D) latitude of tracer (T) points'
  ierr = nf_put_att_text(ncid,vid_lath,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'degrees_east'
  ierr = nf_put_att_text(ncid,vid_lath,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'lonq',NF_DOUBLE,1,dimid_lonq,vid_lonq)
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

  ierr = nf_def_var(ncid,'latq',NF_DOUBLE,1,dimid_latq,vid_latq)
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


  !---- Define 2D geographical coordinate variables
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
  ierr = nf_def_var(ncid,'mask',NF_INT,ndim,dims,vid_lmask)
  attr_name = 'longname'
  attr_valc = '0 if land, 1 if ocean at tracer points'
  ierr = nf_put_att_text(ncid,vid_lmask,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'nondimensional'
  ierr = nf_put_att_text(ncid,vid_lmask,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_def_var(ncid,'depth',NF_REAL,ndim,dims,vid_depth)
  attr_name = 'longname'
  attr_valc = 'grid cell depth'
  ierr = nf_put_att_text(ncid,vid_depth,attr_name,len_trim(attr_valc),trim(attr_valc))
  attr_name = 'units'
  attr_valc = 'm'
  ierr = nf_put_att_text(ncid,vid_depth,attr_name,len_trim(attr_valc),trim(attr_valc))

  ierr = nf_enddef(ncid)

  !---- Write the variables to file
  ierr = nf_put_vara_double(ncid,vid_lonh,1,nx_o,LONH)
  ierr = nf_put_vara_double(ncid,vid_lonq,1,nx_o+1,LONQ)
  ierr = nf_put_vara_double(ncid,vid_lath,1,ny_o,LATH)
  ierr = nf_put_vara_double(ncid,vid_latq,1,ny_o+1,LATQ)

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_double(ncid,vid_geolon,indstr,indcnt,GEOLON)
  ierr = nf_put_vara_double(ncid,vid_geolat,indstr,indcnt,GEOLAT)

  indstr = (/1,1/)
  indcnt = (/nx_o+1,ny_o+1/)
  ierr = nf_put_vara_double(ncid,vid_geolonb,indstr,indcnt,GEOLONB)
  ierr = nf_put_vara_double(ncid,vid_geolatb,indstr,indcnt,GEOLATB)

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_int(ncid,vid_lmask,indstr,indcnt,MASK)
  ierr = nf_put_vara_real(ncid,vid_depth,indstr,indcnt,MODEL_DEPTH)

  ierr = nf_close(ncid)



  stop
end program interp_smooth

function SphericalDistance(lonA, latA, lonB, latB)

  use kinds
  use constants

!!!! Haversine formula (sin squared form)
  ! Calculate distance in meters between two points on idealized spherical earth
  ! given longitudes and lattudes in radians

  implicit none

  real(kind=dbl_kind) SphericalDistance
  real(kind=dbl_kind) :: lonA, latA, lonB, latB


  real(kind=dbl_kind) :: delta_lon, delta_lat
  real(kind=dbl_kind) :: term1, term2, temp, gamma

  !---- Angular separation
  delta_lon = abs(lonA - lonB)
  if ( delta_lon > pi2 ) delta_lon = delta_lon - pi2
  delta_lat = latA - latB

  term1 = sin(p5*delta_lat)**2
  term2 = cos(latA)*cos(latB)*sin(p5*delta_lon)**2
  temp = term1 + term2
  gamma = c2*atan2(sqrt(temp),sqrt(c1-temp))

  SphericalDistance = gamma

end function SphericalDistance

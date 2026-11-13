program lake_fill

  use kinds
  use ncdf_wrapper

  implicit none

  integer :: nx_o, ny_o

  integer (kind=int_kind), dimension(:,:), allocatable :: LMASK

  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       &  GEOLAT  &    ! lat of ocean T-cell corner (Q) points
       &, GEOLON       ! lon of ocean T-cell corner (Q) points

  logical, dimension(:,:), allocatable :: in_ocean, work
  integer:: n_in_ocean, i, j, ibegin, jbegin, iter, n_last_iter, max_iter

  !---- netCDF stuff
  integer :: ierr, ncid, dimid, vid, attid, ndim
  integer :: dimid_lon, dimid_lat
  integer :: did_lon, did_lat
  integer :: vid_geolon, vid_geolat, vid_lmask
  integer, dimension(2) :: indcnt, indstr, dims
  character (char_len) :: cbuffer
  character (char_len) :: vname_lon, vname_lat, vname_mask
  character(char_len) :: dim_name, var_name, attr_name, attr_valc


  character(char_len) :: file_in, file_out

  namelist /file_info/ file_in, vname_mask, file_out, ibegin, jbegin, max_iter

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !---- Intialization

  read(*,file_info)
  write(*,file_info)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Get starting land mask
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  print *,'Opening ',file_in
  ierr = nf_open(file_in,NF_NOWRITE,ncid)
  if ( ierr /= 0 ) then
     print *,'ierr=',ierr,'could not open ',file_in
     stop
  endif

  ierr = nf_inq_dimid(ncid,'lonh',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension nx_o. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_o)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  endif

  ierr = nf_inq_dimid(ncid,'lath',dimid)
  if ( ierr /= 0 ) then
     print *,'could not find dimension ny_o. ierr= ',ierr
     stop
  endif

  ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_o)
  if ( ierr /= 0 ) then
     print *,ierr
     print *,'could not get dimension id= ',dimid
     stop
  endif

  allocate(GEOLON(nx_o,ny_o),GEOLAT(nx_o,ny_o))
  
  ierr = nf_inq_varid(ncid,'geolon',vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable geolon. ierr= ',ierr
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,GEOLON)
  if ( ierr /= 0 ) then
     print *,'could not get variable geolon ierr= ',ierr
     stop
  endif

  ierr = nf_inq_varid(ncid,'geolat',vid)
  if ( ierr /= 0 ) then
     print *,'could not find variable geolat ierr= ',ierr
     stop
  endif

  ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,GEOLAT)
  if ( ierr /= 0 ) then
     print *,'could not get variable geolat ierr= ',ierr
     stop
  endif

  !---- The mask array
  allocate(LMASK(nx_o,ny_o),stat=ierr)
  if (ierr /=0 ) then
     print *,'Could not allocate model grid mask ierr=',ierr
  else
     print *,'Model grid dimension = ',nx_o,ny_o
  endif

  ierr = nf_inq_varid(ncid,trim(vname_mask),vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find mask variable ',trim(vname_mask)
     stop
  endif

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_get_vara_int(ncid,vid,indstr,indcnt,LMASK)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not get mask variable ',vname_mask
     stop
  endif

  ierr = nf_close(ncid)

  write(*,"(a,2f16.7)") ' LMASK LON mix/max=',minval(GEOLON),maxval(GEOLON)
  write(*,"(a,2f16.7)") ' LMASK LAT mix/max=',minval(GEOLAT),maxval(GEOLAT)
  write(*,"(a,2i10)") ' LMASK mix/max=',minval(LMASK),maxval(LMASK)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if ( LMASK(ibegin,jbegin) /= 1 ) then
     print *,' Starting point not in ocean'
     stop 901
  endif

  print *,'Initial number of wet points = ',count(LMASK==1)

  allocate(in_ocean(nx_o,ny_o))

  in_ocean = .false.
  !---- Move north and south from starting point
  do j=jbegin,ny_o
     if ( lmask(ibegin,j) == 0 ) exit
     in_ocean(ibegin,j) = .true.

     do i=ibegin+1,nx_o
        if ( lmask(i,j) == 0 ) exit
        in_ocean(i,j) = .true.
     enddo
     do i=ibegin-1,1,-1
        if ( lmask(i,j) == 0 ) exit
        in_ocean(i,j) = .true.
     enddo
  enddo

  do j=jbegin-1,1,-1
     if ( lmask(ibegin,j) == 0 ) exit
     in_ocean(ibegin,j) = .true.

     do i=ibegin+1,nx_o
        if ( lmask(i,j) == 0 ) exit
        in_ocean(i,j) = .true.
     enddo
     do i=ibegin-1,1,-1
        if ( lmask(i,j) == 0 ) exit
        in_ocean(i,j) = .true.
     enddo
  enddo

  n_last_iter = count(in_ocean)
  print *,'After initial sweep # of points in ocean = ',n_last_iter

  do iter=1,max_iter
     
     do j=2,ny_o-1
        do i=2,nx_o-1
           if (in_ocean(i,j) .AND. lmask(i+1,j)==1) in_ocean(i+1,j) = .true.
           if (in_ocean(i,j) .AND. lmask(i-1,j)==1) in_ocean(i-1,j) = .true.
           if (in_ocean(i,j) .AND. lmask(i,j+1)==1) in_ocean(i,j+1) = .true.
           if (in_ocean(i,j) .AND. lmask(i,j-1)==1) in_ocean(i,j-1) = .true.
        enddo
     enddo

     n_in_ocean = count(in_ocean)
     if ( n_last_iter == n_in_ocean ) exit

     if ( MOD(iter,10) == 0 ) print *,' iter=',iter,' # in ocean = ',n_in_ocean
     n_last_iter = n_in_ocean

  enddo
  print *,'Finished in ',iter,' iterations'

  LMASK = 0
  where(in_ocean) LMASK = 1

  print *,'Final number of wet points = ',count(LMASK==1)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Output result
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ierr = nf_open(trim(file_in),NF_WRITE,ncid)
  if ( ierr /= NF_NOERR ) then
     call handle_err(ierr,'nf_open')
     stop
  endif

  ierr = nf_inq_varid(ncid,trim(vname_mask),vid)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not find mask variable ',trim(vname_mask)
     stop
  endif

  attr_name = 'Revision History'
  ierr = nf_inq_attid(ncid,NF_GLOBAL,trim(attr_name),attid)
  print *,'rev att id = ',attid
  ierr = nf_get_att_text(ncid,NF_GLOBAL,trim(attr_name),attr_valc)
  print *,'attr val=',trim(attr_valc)


  call date_and_time(cbuffer)
  attr_valc = trim(attr_valc) // ': Modified by lake fill ' // trim(cbuffer)

  ierr = nf_put_att_text(ncid,vid_geolon,attr_name,len_trim(attr_valc),trim(attr_valc))

  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_int(ncid,vid,indstr,indcnt,LMASK)
  if ( ierr /= 0 ) then
     print *,' ierr=',ierr
     print *,' Could not put mask variable ',vname_mask
     stop
  endif

  ierr = nf_close(ncid)


!!$  ierr = nf_def_dim(ncid,'xh',nx_o,dimid_lon)
!!$  if ( ierr /= NF_NOERR ) then
!!$     call handle_err(ierr,'nf_def_dim lon')
!!$     stop
!!$  endif
!!$
!!$  ierr = nf_def_dim(ncid,'yh',ny_o,dimid_lat)
!!$  if ( ierr /= NF_NOERR ) then
!!$     call handle_err(ierr,'nf_def_dim lat')
!!$     stop
!!$  endif
!!$
!!$  ndim = 2
!!$  dims = (/dimid_lon,dimid_lat/)
!!$  ierr = nf_def_var(ncid,'geolon',NF_DOUBLE,ndim,dims,vid_geolon)
!!$  if ( ierr /= NF_NOERR ) then
!!$     call handle_err(ierr,'nf_def_var geolon')
!!$     stop
!!$  endif
!!$  attr_name = 'longname'
!!$  attr_valc = 'Longitude of tracer (T) points'
!!$  ierr = nf_put_att_text(ncid,vid_geolon,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$  attr_name = 'units'
!!$  attr_valc = 'degrees_east'
!!$  ierr = nf_put_att_text(ncid,vid_geolon,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$
!!$  ierr = nf_def_var(ncid,'geolat',NF_DOUBLE,ndim,dims,vid_geolat)
!!$  if ( ierr /= NF_NOERR ) then
!!$     call handle_err(ierr,'nf_def_var geolat')
!!$     stop
!!$  endif
!!$  attr_name = 'longname'
!!$  attr_valc = 'Latitude of tracer (T) points'
!!$  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$  attr_name = 'units'
!!$  attr_valc = 'degrees_north'
!!$  ierr = nf_put_att_text(ncid,vid_geolat,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$
!!$  ierr = nf_def_var(ncid,'wet',NF_INT,ndim,dims,vid_lmask)
!!$  attr_name = 'longname'
!!$  attr_valc = '0 if land, 1 if ocean at tracer points'
!!$  ierr = nf_put_att_text(ncid,vid_lmask,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$  attr_name = 'units'
!!$  attr_valc = 'dimensionless'
!!$  ierr = nf_put_att_text(ncid,vid_lmask,attr_name,len_trim(attr_valc),trim(attr_valc))
!!$
!!$  ierr = nf_enddef(ncid)


  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
!!$  ierr = nf_put_vara_double(ncid,vid_geolon,indstr,indcnt,GEOLON)
!!$  ierr = nf_put_vara_double(ncid,vid_geolat,indstr,indcnt,GEOLAT)
!!$
  indstr = (/1,1/)
  indcnt = (/nx_o,ny_o/)
  ierr = nf_put_vara_int(ncid,vid_lmask,indstr,indcnt,LMASK)

  ierr = nf_close(ncid)



stop
end program lake_fill

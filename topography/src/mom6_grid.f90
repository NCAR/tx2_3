module mom6_grid

  use kinds
  use constants

  implicit none
  save

  integer(kind=int_kind)  :: &
       & nx_sg, ny_sg, &         ! super grid (mosaic) dimensions
       & nx_o, ny_o              ! model grid dimensions
       
  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       & X_SG, Y_SG, &                            ! MODEL Super grid
       & DX_SG, DY_SG  ,&
       & DAREA_SG, ANGLE_SG


  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       &  GEOLAT, GEOLON  &     ! lat/lon of ocean T-cell center points
       &, GEOLATB, GEOLONB      ! lat/lon of ocean T-cell corner points

  real (kind=dbl_kind), dimension(:), allocatable :: &
       &  LATH, LONH  &    ! nominal (1D) lat/lon of ocean T-cell center points
       &, LATQ, LONQ       ! nominal (1D) lat/lon of ocean T-cell corner points

  real (kind=dbl_kind), dimension(:,:), allocatable :: &
       &  DLATT, DYT  &
       &, DLONT, DXT  &
       &, DAREA

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine read_super_grid(file_super_grid)

    use ncdf_wrapper

    character(char_len) :: file_super_grid

    character (char_len) :: cbuffer
    integer :: ierr, ncid, dimid, vid
    integer, dimension(2) :: indstr, indcnt


    !---- Open the file
    print *,'Opening ',trim(file_super_grid)
    ierr = nf_open(file_super_grid,NF_NOWRITE,ncid)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,'could not open ',file_super_grid
       stop
    endif

    !---- Get the grid size
    ierr = nf_inq_dimid(ncid,'nx',dimid)
    if ( ierr /= 0 ) then
       print *,'could not find dimension nx. ierr= ',ierr
       stop
    endif

    ierr = nf_inq_dim(ncid,dimid,cbuffer,nx_sg)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not get dimension id= ',dimid
       stop
    endif

    ierr = nf_inq_dimid(ncid,'ny',dimid)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not find dimension ny.'
       stop
    endif

    ierr = nf_inq_dim(ncid,dimid,cbuffer,ny_sg)
    if ( ierr /= 0 ) then
       print *,'ierr =',ierr,'could not get dimension id= ',dimid
       stop
    endif

    !---- Allocate space for super grid variables
    allocate(X_SG(nx_sg+1,ny_sg+1),Y_SG(nx_sg+1,ny_sg+1),&
         &   DX_SG(nx_sg,ny_sg+1), DY_SG(nx_sg+1,ny_sg),&
         &   DAREA_SG(nx_sg,ny_sg), ANGLE_SG(nx_sg+1,ny_sg+1),&
         &   stat=ierr)
    if ( ierr /= 0 ) then
       print *,'ERROR Could not allocate super grid ierr=',ierr
    else
       print *,'Super grid dimensions = ',nx_sg,ny_sg
    endif

    !---- Read the super grid variables
    indstr = (/1,1/)
    indcnt = (/nx_sg+1,ny_sg+1/)

    ierr = nf_inq_varid(ncid,'x',vid)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not find variable x'
       stop
    endif
    ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,X_SG)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not get variable x'
       stop
    endif
    print *,'X_SG min/max = ',minval(X_SG),maxval(X_SG)

    ierr = nf_inq_varid(ncid,'y',vid)
    if ( ierr /= 0 ) then
       print *,'ierr = ',ierr,' could not find variable y'
       stop
    endif
    ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,Y_SG)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not get variable y'
       stop
    endif
    print *,'Y_SG min/max = ',minval(Y_SG),maxval(Y_SG)

    indstr = (/1,1/)
    indcnt = (/nx_sg,ny_sg+1/)
    ierr = nf_inq_varid(ncid,'dx',vid)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not find variable dx'
       stop
    endif
    ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,DX_SG)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not get variable dx'
       stop
    endif
    print *,'DX_SG min/max = ',minval(DX_SG),maxval(DX_SG)

    indstr = (/1,1/)
    indcnt = (/nx_sg+1,ny_sg/)
    ierr = nf_inq_varid(ncid,'dy',vid)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not find variable dy'
       stop
    endif
    ierr = nf_get_vara_double(ncid,vid,indstr,indcnt,DY_SG)
    if ( ierr /= 0 ) then
       print *,'ierr=',ierr,' could not get variable dy'
       stop
    endif
    print *,'DY_SG min/max = ',minval(DY_SG),maxval(DY_SG)

    ierr = nf_close(ncid)

    return
  end subroutine read_super_grid

  subroutine compute_mom6_grid()

    integer(kind=int_kind) :: ierr
    integer(kind=int_kind), dimension(1) :: imax

    !---- Allocate space for the model grid arrays
    nx_o = nx_sg/2
    ny_o = ny_sg/2

    allocate(GEOLON(nx_o,ny_o),GEOLAT(nx_o,ny_o), &
         &   GEOLONB(nx_o+1,ny_o+1),GEOLATB(nx_o+1,ny_o+1), &
         &   LONH(nx_o),LATH(ny_o), &
         &   LONQ(nx_o+1),LATQ(ny_o+1), &
         &  stat=ierr)
    if (ierr /=0 ) then
       print *,'Could not allocate model grid mask ierr=',ierr
    else
       print *,'Model grid dimension = ',nx_o,ny_o
    endif

    GEOLON = X_SG(2::2,2::2)
    GEOLAT = Y_SG(2::2,2::2)

    GEOLONB = X_SG(1::2,1::2)
    GEOLATB = Y_SG(1::2,1::2)

    print *,'T-point LON min/max = ',minval(GEOLON),maxval(GEOLON)
    print *,'T-point LAT min/max = ',minval(GEOLAT),maxval(GEOLAT)
    print *,'Q-point LON min/max = ',minval(GEOLONB),maxval(GEOLONB)
    print *,'Q-point LAT min/max = ',minval(GEOLATB),maxval(GEOLATB)

    LONH = GEOLON(:,1)
    LONQ = GEOLONB(:,1)

    imax = maxloc(GEOLAT(:,ny_o))
    LATH = GEOLAT(imax(1),:)
    imax = maxloc(GEOLATB(:,ny_o+1))
    LATQ = GEOLATB(imax(1),:)

    return
  end subroutine compute_mom6_grid

  subroutine compute_mom6_grid_metrics()

    integer(kind=int_kind) :: ierr

    allocate(DLONT(nx_o,ny_o),DLATT(nx_o,ny_o), &
         &   DXT(nx_o,ny_o),DYT(nx_o,ny_o),DAREA(nx_o,ny_o), &
         &  stat=ierr)
    if (ierr /=0 ) then
       print *,'Could not allocate model grid mask ierr=',ierr
    else
       print *,'Model grid dimension = ',nx_o,ny_o
    endif

    DLONT = (X_SG(3::2,2::2) - X_SG(1::2,2::2))
    DLATT = (Y_SG(2::2,3::2) - Y_SG(2::2,1::2))

    DXT = DX_SG(1::2,2::2) + DX_SG(2::2,2::2)
    DYT = DY_SG(2::2,1::2) + DX_SG(2::2,2::2)

    DAREA = DXT*DYT

    print *,'DLONT min/max = ',minval(DLONT),maxval(DLONT)
    print *,'DLATT min/max = ',minval(DLATT),maxval(DLATT)
    print *,'DX min/max = ',minval(DXT),maxval(DXT)
    print *,'DY min/max = ',minval(DYT),maxval(DYT)
    print *,'DAREA min/max = ',minval(DAREA),maxval(DAREA)
    print *

    return
  end subroutine compute_mom6_grid_metrics

end module mom6_grid

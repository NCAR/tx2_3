MODULE gridio
  USE param
  USE common
  USE nf90util
  USE netcdf, ONLY: nf90_put_att
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: hgrid_file = 'ocean_hgrid.nc'
  CHARACTER(LEN=*), PARAMETER :: global_file = 'coordinates.nc'
  CHARACTER(LEN=*), PARAMETER :: north_file = 'coordinates_north.nc'

CONTAINS


  SUBROUTINE write_hgrid

    INTEGER :: ncid, stat, nx, ny, nxs, nys, nxsp, nysp, string, i, j
    TYPE(ncdimtype) :: dim_nxp, dim_nyp, dim_nx, dim_ny, dim_str
    TYPE(ncvartype) :: var_str, var_x, var_y, var_dx, var_dy, var_area, var_angle
    REAL(wp) :: s
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: x, y, dx, dy, area, angle
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: lon, lat
    CHARACTER(LEN=8) :: date
    CHARACTER(LEN=10) :: time
    CHARACTER(LEN=5) :: zone
    CHARACTER(LEN=50) :: date_att

    nx = SIZE(glam,1)
    ny = SIZE(glam,2)

    nxs = 2 * (nx-1) - 2
    nys = 2 * (ny-1) - 2
    nxsp = nxs+1
    nysp = nys+1
    string = 255

    ALLOCATE ( x(nxsp, nysp), y(nxsp, nysp), dx(nxs, nysp), dy(nxsp, nys), area(nxs,nys), angle(nxsp, nysp) )
    ALLOCATE ( lon(nxsp+1, nysp+1), lat(nxsp+1, nysp+1) )


    ! MOM6 grid coordinates
    x(1:nxsp:2, 1:nysp:2) = glam(1:nx-1, 1:ny-1, 4)
    x(2:nxsp:2, 1:nysp:2) = glam(2:nx-1, 1:ny-1, 3)
    x(1:nxsp:2, 2:nysp:2) = glam(1:nx-1, 2:ny-1, 2)
    x(2:nxsp:2, 2:nysp:2) = glam(2:nx-1, 2:ny-1, 1)

    y(1:nxsp:2, 1:nysp:2) = gphi(1:nx-1, 1:ny-1, 4)
    y(2:nxsp:2, 1:nysp:2) = gphi(2:nx-1, 1:ny-1, 3)
    y(1:nxsp:2, 2:nysp:2) = gphi(1:nx-1, 2:ny-1, 2)
    y(2:nxsp:2, 2:nysp:2) = gphi(2:nx-1, 2:ny-1, 1)

    WHERE(x > MINVAL(x(1,:))) x = x - 360_wp !longitude monotonically increasing
    x(1,:) = x(1,:) - 360_wp
    WHERE(x(1:nxsp/2,nysp) == 73.) x(1:nxsp/2,nysp) = x(1:nxsp/2,nysp) - 360_wp


    ! MOM6 grid coordinates
    dx(1:nxs:2, 1:nysp:2) = 0.25 * ( e1(1:nx-2, 1:ny-1, 4) + e1(2:nx-1, 1:ny-1, 3) )
    dx(2:nxs:2, 1:nysp:2) = 0.25 * ( e1(2:nx-1, 1:ny-1, 3) + e1(2:nx-1, 1:ny-1, 4) )
    dx(1:nxs:2, 2:nysp:2) = 0.25 * ( e1(1:nx-2, 2:ny-1, 2) + e1(2:nx-1, 2:ny-1, 1) )
    dx(2:nxs:2, 2:nysp:2) = 0.25 * ( e1(2:nx-1, 2:ny-1, 1) + e1(2:nx-1, 2:ny-1, 2) )

    dy(1:nxsp:2, 1:nys:2) = 0.25 * ( e2(1:nx-1, 1:ny-2, 4) + e2(1:nx-1, 2:ny-1, 2) )
    dy(1:nxsp:2, 2:nys:2) = 0.25 * ( e2(1:nx-1, 2:ny-1, 2) + e2(1:nx-1, 2:ny-1, 4) )
    dy(2:nxsp:2, 1:nys:2) = 0.25 * ( e2(2:nx-1, 1:ny-2, 3) + e2(2:nx-1, 2:ny-1, 1) )
    dy(2:nxsp:2, 2:nys:2) = 0.25 * ( e2(2:nx-1, 2:ny-1, 1) + e2(2:nx-1, 2:ny-1, 3) )
   
 
    !area assuming cell are cyclic quadrilateral
    !area = sqrt((s-a)*(s-b)*(s-c)*(s-d))
    !with s = 0.5 * (a+b+c+d)
    DO i = 1, nxs 
      DO j = 1, nys
        s = 0.5 * (dx(i,j) + dx(i,j+1) + dy(i,j) + dy(i+1,j))
        area(i,j) = SQRT((s-dx(i,j))*(s-dx(i,j+1))*(s-dy(i,j))*(s-dy(i+1,j)))
        IF (ISNAN(area(i,j))) area(i,j) = (dx(i,j) + dx(i,j+1)) * (dy(i,j) + dy(i+1,j)) / 4_wp
      END DO
    END DO


    !forward azimuth
    !theta = atan2( sin(lamda2-lambda1) * cos(phi2), 
    !               cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(lamda2-lambda1) )
    ! where lambda1/phi1 is the start point,lamda2/phi2 the enpoint
    lon(1:nxsp+1:2,   1:nysp+1:2)   = glam(1:nx-1, 1:ny-1, 4)
    lon(2:nxsp+1:2,   1:nysp+1:2)   = glam(2:nx, 1:ny-1, 3)
    lon(1:nxsp+1:2,   2:nysp+1:2)   = glam(1:nx-1, 2:ny, 2)
    lon(2:nxsp+1:2,   2:nysp+1:2)   = glam(2:nx, 2:ny, 1)

    lat(1:nxsp+1:2,   1:nysp+1:2)   = gphi(1:nx-1, 1:ny-1, 4)
    lat(2:nxsp+1:2,   1:nysp+1:2)   = gphi(2:nx, 1:ny-1, 3)
    lat(1:nxsp+1:2,   2:nysp+1:2)   = gphi(1:nx-1, 2:ny, 2)
    lat(2:nxsp+1:2,   2:nysp+1:2)   = gphi(2:nx, 2:ny, 1)

    lon = lon * rad
    lat = lat * rad

    DO i = 1, nxsp 
      DO j = 1, nysp
        angle(i,j) = ATAN2(SIN(lon(i+1,j)-lon(i,j)) * COS(lat(i+1,j)), &
                           COS(lat(i,j)) * SIN(lat(i+1,j)) - SIN(lat(i,j)) * &
                           COS(lat(i+1,j)) * COS(lon(i+1,j)-lon(i,j))) * 1_wp / rad 
      END DO
    END DO
    angle(:,:) = 90_wp - angle(:,:) !angle between X-axis and true EAST


    CALL create_file( hgrid_file, nf90_netcdf4, ncid )

    !add global attribut
    stat = nf90_put_att ( ncid, nf90_global, "Description", &
            "MOM6 2/3 degree tripolar grid (ORCA type)" )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "write_hgrid", &
         "putting global attribute: Description=MOM6 2/3 degree tripolar grid &
         (ORCA type)", ncid=ncid )
    stat = nf90_put_att ( ncid, nf90_global, "Author", "Fred Castruccio (fredc@ucar.edu)" )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "write_hgrid", &
         "putting global attribute: Author=Fred Castruccio (fredc@ucar.edu)", ncid=ncid )
    CALL DATE_AND_TIME ( date=date, time=time, zone=zone )                                         
    date_att = &
           date(7:8) // "/" // date(5:6) // "/" // date(1:4) // " " // &                              
           time(1:2) // ":" // time(3:4) // ":" // time(5:6) // " " // &                              
           zone
    stat = nf90_put_att ( ncid, nf90_global, "Created", date_att )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "write_hgrid", &
         "putting global attribute: Created=" // date_att, ncid=ncid )
    stat = nf90_put_att ( ncid, nf90_global, "type", "MOM6 supergrid file" )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "write_hgrid", &
         "putting global attribute: type=MOM6 supergrid file", ncid=ncid )

    dim_ny%ncid = ncid
    dim_nyp%ncid = ncid
    dim_nyp%name = 'nyp'
    dim_nyp%len  = nysp
    CALL def_dim( dim_nyp )

    dim_nxp%ncid = ncid
    dim_nxp%name = 'nxp'
    dim_nxp%len  = nxsp
    CALL def_dim( dim_nxp )

    dim_ny%ncid = ncid
    dim_ny%name = 'ny'
    dim_ny%len  = nys
    CALL def_dim( dim_ny )

    dim_nx%ncid = ncid
    dim_nx%name = 'nx'
    dim_nx%len  = nxs
    CALL def_dim( dim_nx )

    dim_str%ncid = ncid
    dim_str%name = 'string'
    dim_str%len  = string
    CALL def_dim( dim_str )

    var_str%ncid = ncid
    var_str%name = 'tile'
    var_str%xtype = nf90_char
    CALL def_var( var_str, dims=(/dim_str/) )

    var_dy%ncid = ncid
    var_y%ncid = ncid
    var_y%name = 'y'
    var_y%xtype = nf90_double
    CALL def_var( var_y, dims=(/dim_nxp,dim_nyp/), units="degrees" )

    var_x%ncid = ncid
    var_x%name = 'x'
    var_x%xtype = nf90_double
    CALL def_var( var_x, dims=(/dim_nxp,dim_nyp/), units="degrees" )

    var_dy%ncid = ncid
    var_dy%name = 'dy'
    var_dy%xtype = nf90_double
    CALL def_var( var_dy, dims=(/dim_nxp,dim_ny/), units="meters" )

    var_dx%ncid = ncid
    var_dx%name = 'dx'
    var_dx%xtype = nf90_double
    CALL def_var( var_dx, dims=(/dim_nx,dim_nyp/), units="meters" )

    var_area%ncid = ncid
    var_area%name = 'area'
    var_area%xtype = nf90_double
    CALL def_var( var_area, dims=(/dim_nx,dim_ny/), units="m2" )

    var_angle%ncid = ncid
    var_angle%name = 'angle_dx'
    var_angle%xtype = nf90_double
    CALL def_var( var_angle, dims=(/dim_nxp,dim_nyp/), units="meters" )

    CALL put_var( var_str, 'tile1' )
    CALL put_var( var_y, y(:,:) )
    CALL put_var( var_x, x(:,:) )
    CALL put_var( var_dy, dy(:,:) )
    CALL put_var( var_dx, dx(:,:) )
    CALL put_var( var_area, area(:,:) )
    CALL put_var( var_angle, angle(:,:) )

    CALL close_file( ncid )

    DEALLOCATE (x, y, dx, dy, area, angle)

  END SUBROUTINE write_hgrid


  SUBROUTINE write_glo

    INTEGER :: ncid, stat, nx, ny, ng, jv, jg
    TYPE(ncdimtype) :: dimx, dimy, dimg
    TYPE(ncvartype), DIMENSION(2) :: navs
    TYPE(ncvartype), DIMENSION(16) :: vars
    TYPE(ncvartype) :: aniso

    ! To do, add some scalars defining the grid

    nx = SIZE(glam,1)
    ny = SIZE(glam,2)
    ng = SIZE(glam,3)

    CALL create_file( global_file, nf90_netcdf4, ncid )

    dimx%ncid = ncid
    dimx%name = 'x'
    dimx%len  = nx
    CALL def_dim( dimx )

    dimy%ncid = ncid
    dimy%name = 'y'
    dimy%len  = ny
    CALL def_dim( dimy )

    dimg%ncid = ncid
    dimg%name = 'grid'
    dimg%len  = 4
    CALL def_dim( dimg )

    navs(:)%ncid = ncid
    navs(:)%name = (/ 'nav_lon', 'nav_lat' /)
    navs(:)%xtype = nf90_float
    CALL def_var( navs(1), dims=(/dimx,dimy/), units="degrees_east" )
    CALL def_var( navs(2), dims=(/dimx,dimy/), units="degrees_north" )
    
    vars(:)%ncid = ncid
    vars(1:8)%name = (/ 'glamt', 'glamu', 'glamv', 'glamf', &
                        'gphit', 'gphiu', 'gphiv', 'gphif'   /)
    vars(9:16)%name = (/ 'e1t',   'e1u',   'e1v',   'e1f',  &
                         'e2t',   'e2u',   'e2v',   'e2f'    /)
    vars(:)%xtype = nf90_double
    DO jv = 1, 16
       CALL def_var( vars(jv), dims=(/dimx,dimy/) )
    END DO

    aniso%ncid = ncid
    aniso%name = 'anisotropy'
    aniso%xtype = nf90_double
    CALL def_var( aniso, dims=(/dimx,dimy,dimg/) )

    CALL put_var( navs(1), glam(:,:,1) )
    CALL put_var( navs(2), gphi(:,:,1) )
    DO jg = 1, 4
       CALL put_var( vars(jg),    glam(:,:,jg) )
       CALL put_var( vars(jg+4),  gphi(:,:,jg) )
       CALL put_var( vars(jg+8),    e1(:,:,jg) )
       CALL put_var( vars(jg+12),   e2(:,:,jg) )
    END DO
    CALL put_var( aniso, e1(:,:,:)/e2(:,:,:) )

    CALL close_file( ncid )


  END SUBROUTINE write_glo



  SUBROUTINE write_nth

     INTEGER :: ncid, jv
     TYPE(ncdimtype) :: dimx, dimy
     TYPE(ncvartype), DIMENSION(7) :: paras
     TYPE(ncvartype), DIMENSION(4) :: vars


     IF ( l_write_binary ) CALL write_nth_binary
    
     CALL create_file( north_file, nf90_netcdf4, ncid )
     
     dimx%ncid = ncid
     dimx%name = 'x'
     dimx%len  = jpin
     CALL def_dim( dimx )

     dimy%ncid = ncid
     dimy%name = 'y'
     dimy%len  = jpjn
     CALL def_dim( dimy )

     paras(:)%ncid = ncid
     paras(:)%name = (/ "ra    ", "jpi   ", "jpj   ", "jpeq  ", &
                     "jpnord", "jpin  ", "jpjn  " /)
     paras(1  )%xtype = nf90_double
     paras(2:7)%xtype = nf90_int
     paras(:)%ndims = 0
     DO jv = 1, 7
        CALL def_var( paras(jv) )
     END DO
     
     vars(:)%ncid = ncid
     vars(:)%xtype = nf90_double
     vars(:)%ndims = 2
     vars(:)%name = (/ 'glamnth', 'gphinth', 'e1nth  ', 'e2nth  ' /)
     DO jv = 1, 4
        CALL def_var( vars(jv), dims=(/dimx,dimy/) )
     END DO
     
     CALL put_var( paras(1), ra )
     CALL put_var( paras(2), jpi )
     CALL put_var( paras(3), jpj )
     CALL put_var( paras(4), jpeq )
     CALL put_var( paras(5), jpnord )
     CALL put_var( paras(6), jpin )
     CALL put_var( paras(7), jpjn )

     CALL put_var( vars(1), glamnth )
     CALL put_var( vars(2), gphinth )
     CALL put_var( vars(3), e1nth )
     CALL put_var( vars(4), e2nth )

     CALL close_file( ncid )


  END SUBROUTINE write_nth




  SUBROUTINE read_nth

     INTEGER :: ncid, jv
     TYPE(ncdimtype) :: dimx, dimy
     TYPE(ncvartype), DIMENSION(7) :: paras
     TYPE(ncvartype), DIMENSION(4) :: vars
     REAL(wp) :: zra
     INTEGER :: ijpi, ijpj, ijpeq, ijpnord, ijpin, ijpjn

     write(0,*)
     write(0,*)
     write(0,*) 'appel de hgrlec'
     write(0,*) '==============='
     write(0,*)
     !
     !
     CALL open_file( north_file, nf90_nowrite, ncid )

     paras(:)%ncid = ncid
     paras(:)%name =  (/ "ra    ", "jpi   ", "jpj   ", "jpeq  ", &
                         "jpnord", "jpin  ", "jpjn  " /)
     DO jv = 1, 7
        CALL fill_var_info( paras(jv) )
     END DO

     CALL get_var( paras(1), zra )
     CALL get_var( paras(2), ijpi )
     CALL get_var( paras(3), ijpj )
     CALL get_var( paras(4), ijpeq )
     CALL get_var( paras(5), ijpnord )
     CALL get_var( paras(6), ijpin )
     CALL get_var( paras(7), ijpjn )

! ... controle
! ... check (it would be good to allow chages to south part only)
!           (cant change jpi, jpnord or jpeqn=jpeq-jpeqt)
     IF (  (     zra .NE.     ra ) .OR. (  ijpi .NE.  jpi ) .OR. &
           (    ijpj .NE.    jpj ) .OR. ( ijpeq .NE. jpeq ) .OR. &
           ( ijpnord .NE. jpnord ) .OR. ( ijpin .NE. jpin ) .OR. &
           (    jpjn .NE.   jpjn )   ) THEN
        WRITE(0,*) 
        WRITE(0,*) 'inconsitency, input mesh and parameter:'
        WRITE(0,*) '======================================='
        WRITE(0,*) 
        WRITE(0,*) '    ra    = ', ra    ,'   ra     read= ', zra
        WRITE(0,*) '    jpi   = ', jpi   ,'   jpi    read= ', ijpi
        WRITE(0,*) '    jpj   = ', jpj   ,'   jpj    read= ', ijpj
        WRITE(0,*) '    jpeq  = ', jpeq  ,'   jpeq   read= ', ijpeq
        WRITE(0,*) '    jpnord= ', jpnord,'   jpnord read= ', ijpnord
        WRITE(0,*) '    jpin  = ', jpin  ,'   jpin   read= ', ijpin
        WRITE(0,*) '    jpjn  = ', jpjn  ,'   jpjn   read= ', ijpjn

!!$ To do:
!!$        WRITE(0,*) 'Run again with the values in "consistentNML"'
!!$        OPEN( FILE="consistentNML", UNIT=10, FORM="FORMATTED" )
!!$        WRITE( 10, NML=consistentNML )
!!$        CLOSE( 10 )
        STOP 'Need to change parameters'
     ENDIF

     vars(:)%ncid = ncid
     vars(:)%name = (/ 'glamnth', 'gphinth', 'e1nth  ', 'e2nth  ' /)
     DO jv = 1, 4
        CALL fill_var_info( vars(jv) )
     END DO

     CALL get_var( vars(1), glamnth )
     CALL get_var( vars(2), gphinth )
     CALL get_var( vars(3), e1nth )
     CALL get_var( vars(4), e2nth )

     CALL close_file( ncid )


  END SUBROUTINE read_nth


  SUBROUTINE write_nth_binary
    
    INTEGER :: ipart, ideb, ifin, numwri
    CHARACTER(LEN=nf90_max_name) :: clname

    !
    !
    ! II. Output de la demi-grille nord
    ! -----------------------------
    !
    WRITE(0,*)
    WRITE(0,*) ' hgr_nth: ecriture de la demi-grille nord dans numwri'
    WRITE(0,*) '         partie: ', ipart,' de ji= ',ideb,' a ', ifin
    WRITE(0,*)
    !
    ! open file numwri
    numwri=30
    clname='hgr_nth'
!     OPEN (numwri,FILE='hgr_nth.output',
    OPEN (numwri,FILE=clname, &
                  FORM='UNFORMATTED',STATUS='UNKNOWN')
    WRITE(0,*)
    WRITE(0,*) ' open file '//clname//'.output OK, unit ',numwri
    WRITE(0,*)
    !
    REWIND(numwri)
    WRITE (numwri) ra, jpi, jpj, jpeq, jpnord, jpin, jpjn
!      DO ji = 1, 10
!        iji = 10*(ipart-1)+ji
!        DO jj = 1, jpjn
!          zglam(ji,jj) = glam(iji,jj)
!          zgphi(ji,jj) = gphi(iji,jj)
!          ze1  (ji,jj) = e1  (iji,jj)
!          ze2  (ji,jj) = e2  (iji,jj)
!        END DO
!      END DO
    WRITE (numwri) glamnth
    WRITE (numwri) gphinth
    WRITE (numwri) e1nth
    WRITE (numwri) e2nth
    CLOSE(numwri)


  END SUBROUTINE write_nth_binary


END MODULE gridio

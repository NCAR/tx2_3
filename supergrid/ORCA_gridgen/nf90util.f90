MODULE nf90util
  ! Utility module to make the use of netcdf a bit more high-level
  ! Provides structures for storing object information and routines for 
  ! using these structures.
  ! Takes care of the error handling, which makes the calling code far more
  ! readable (which is my only dislike of the netcdf interface)

  USE kinds
  USE netcdf
  !!----------------------------------------------------------------------
  !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nc_files = 1000  ! Max concurrently open files
  INTEGER, PARAMETER :: max_var_dims = 4     ! Max dimensions per variable
  INTEGER, PARAMETER :: eNoEnt = 2  ! System err val for "No such file or dir"
  REAL(wp), PARAMETER :: mdi_unset = 21.202  ! A value unlikely to be used as mdi

  TYPE ncdimtype
     INTEGER :: ncid=-1, id=-1, len
     CHARACTER(LEN=nf90_max_name) :: name=""
     LOGICAL :: lunlim=.FALSE.
  END TYPE ncdimtype

  TYPE ncvartype
     ! contains all the info to allocate, read, write a netcdf variable
     INTEGER :: ncid, id, ndims, xtype=nf90_float
     CHARACTER(LEN=nf90_max_name) :: name="", grid="", levs=""
     REAL(wp) :: mdi = mdi_unset
     INTEGER, DIMENSION(max_var_dims) :: dimids, dimlens
     ! warning: duplicating dim info - be sure to keep consistent
  END TYPE ncvartype

  INTERFACE get_var
     MODULE PROCEDURE get_var_0d_real, get_var_0d_double, get_var_0d_int, &
                      get_var_1d_real, get_var_1d_double, get_var_1d_int, &
                      get_var_2d_real, get_var_2d_double, get_var_2d_int, &
                      get_var_3d_real, get_var_3d_double, get_var_3d_int, &
                      get_var_4d_real, get_var_4d_double, get_var_4d_int, &
                      get_var_1d_char, get_var_2d_char
  END INTERFACE
  
  INTERFACE put_var
     MODULE PROCEDURE put_var_0d_real, put_var_0d_double, put_var_0d_int, &
                      put_var_1d_real, put_var_1d_double, put_var_1d_int, &
                      put_var_2d_real, put_var_2d_double, put_var_2d_int, &
                      put_var_3d_real, put_var_3d_double, put_var_3d_int, &
                      put_var_4d_real, put_var_4d_double, put_var_4d_int, &
                      put_var_1d_char, put_var_2d_char
  END INTERFACE
  

  ! Names of currently open datasets
  CHARACTER(LEN=nf90_max_name), DIMENSION(0:max_nc_files-1) :: dataset_filename=""

  INTEGER :: nverbose = 2   ! Level of verbosity - default minimal
  INTEGER :: ncid_min=1000, ncid_max=-1 ! Extrema of previously used NC ids

CONTAINS

  SUBROUTINE handle_error( stat, routine, action1, action2, action3, &
       dim, var, var2, ncid, nostop )
    INTEGER, INTENT(IN) :: stat
    CHARACTER(LEN=*), INTENT(IN) :: routine
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) ::action1, action2, action3
    TYPE(ncdimtype), OPTIONAL, INTENT(IN) :: dim
    TYPE(ncvartype), OPTIONAL, INTENT(IN) :: var, var2
    INTEGER, OPTIONAL, INTENT(IN) :: ncid
    LOGICAL, OPTIONAL, INTENT(IN) :: nostop
    LOGICAL :: stop_here

    IF (PRESENT(nostop)) THEN
       stop_here = .NOT.nostop
    ELSE
       stop_here = .TRUE.
    END IF

    WRITE(*,*)
    WRITE(*,*) "NetCDF error: ",TRIM(nf90_strerror(stat))
    WRITE(*,*) "in routine ",routine
    IF (PRESENT(action1)) WRITE(*,*) "while ", action1
    IF (PRESENT(action2)) WRITE(*,*) action2
    IF (PRESENT(action3)) WRITE(*,*) action3
    WRITE(*,*) ""

    IF (PRESENT(ncid)) CALL print_ncid( ncid )

    IF (PRESENT(dim)) CALL print_dim( dim )

    IF (PRESENT(var)) CALL print_var( var )
    IF (PRESENT(var2)) CALL print_var( var2 )

    IF (stop_here) THEN
       STOP 9
    ELSE
       WRITE(*,*) "Continuing despite the error"
    END IF

    WRITE(*,*) ""

  END SUBROUTINE handle_error


  SUBROUTINE print_dimlens( var, array_shape, string_length )
    TYPE(ncvartype) :: var
    INTEGER, DIMENSION(:), INTENT(IN) :: array_shape
    INTEGER, OPTIONAL, INTENT(IN) :: string_length
    INTEGER :: rank, ndim

    rank = SIZE( array_shape )
    ndim = var%ndims

    WRITE(*,*) "Rank of array:   ", rank
    WRITE(*,*) "NetCDF var ndim: ", ndim

    IF (PRESENT(string_length)) WRITE(*,*) "String length:   ", string_length
    IF (rank > 0) WRITE(*,*) "Array dimensions:", array_shape
    IF (ndim > 0) WRITE(*,*) "NetCDF dimlenths:", var%dimlens(:ndim)

    CALL print_var( var )

  END SUBROUTINE print_dimlens

  
  SUBROUTINE print_ncid( ncid )
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: stat 

    WRITE(*,*) "Dataset details"
    WRITE(*,*) "-----------------"
    WRITE(*,*) "ID:        ", ncid
    !IF ( ncid > 0 ) WRITE(*,*) "Filename:  ", TRIM(dataset_filename(ncid)) !FC
    WRITE(*,*) ""

  END SUBROUTINE print_ncid



  SUBROUTINE print_dim( dim )
    TYPE(ncdimtype), INTENT(IN) :: dim
    CHARACTER(LEN=3) :: yesno

    CALL print_ncid( dim%ncid )

    WRITE(*,*) "Dimension details"
    WRITE(*,*) "-----------------"
    WRITE(*,*) "name:     ", TRIM(dim%name)
    WRITE(*,*) "ID:       ", dim%id
    WRITE(*,*) "length:   ", dim%len
    IF (dim%lunlim) THEN
       yesno = "yes"
    ELSE
       yesno = "no"
    END IF
    WRITE(*,*) "unlimited:     ", yesno
    WRITE(*,*) ""

  END SUBROUTINE print_dim


  SUBROUTINE print_var( var )
    TYPE(ncvartype), INTENT(IN) :: var
    CHARACTER(LEN=6) :: ctype
    CHARACTER(LEN=*), PARAMETER :: cfmt='(10x,"(dimension",i3," is the unlimited dim)")'
    INTEGER :: stat, jd, idunlim

    CALL print_ncid( var%ncid )

    WRITE(*,*) "Variable details"
    WRITE(*,*) "----------------"
    WRITE(*,*) "name:     ", TRIM(var%name)
    WRITE(*,*) "ID:       ", var%id
    WRITE(*,*) "ndims:    ", var%ndims
    WRITE(*,FMT='(" dimids:  ",4i8)') var%dimids(:)
    WRITE(*,FMT='(" dimlens: ",4i8)') var%dimlens(:)

    stat = nf90_inquire( var%ncid, unlimitedDimId=idunlim )
    SELECT CASE( stat )
    CASE( nf90_noerr )
       DO jd = 1, var%ndims
          IF ( var%dimids(jd) == idunlim )  WRITE(*,cfmt) jd
       END DO
    CASE( nf90_ebadid )
       ! Do nothing
    CASE DEFAULT
       IF (stat /= nf90_noerr)   CALL handle_error( stat, "print_var", &
            "inquiring file for unlim ID", ncid=var%ncid )
    END SELECT

    SELECT CASE( var%xtype )
    CASE( nf90_byte )
       ctype = "byte"
    CASE( nf90_char )
       ctype = "char"
    CASE( nf90_short )
       ctype = "short"
    CASE( nf90_int )
       ctype = "int"
    CASE( nf90_float )
       ctype = "float"
    CASE( nf90_double )
       ctype = "double"
    CASE DEFAULT
       ctype = "stupid"
    END SELECT
    WRITE(*,*) "xtype:   ", ctype
    WRITE(*,*) "grid:    ", TRIM(var%grid)

    IF ( var%mdi /= mdi_unset ) WRITE(*,*) "mdi:     ", var%mdi

    WRITE(*,*) ""


  END SUBROUTINE print_var


  SUBROUTINE open_file( filename, mode, ncid, exist, define )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: mode
    INTEGER, INTENT(OUT) :: ncid
    LOGICAL, INTENT(OUT), OPTIONAL :: exist
    LOGICAL, OPTIONAL, INTENT(IN) :: define

    INTEGER :: stat, jid
    CHARACTER(LEN=7) :: cmode
    LOGICAL :: ldefine_mode, ll_old_file

    IF (PRESENT(exist)) exist = .TRUE.

    SELECT CASE (mode)
    CASE( nf90_write )
       cmode = "write"
    CASE( nf90_nowrite )
       cmode = "nowrite"
    CASE( nf90_share )
       cmode = "share"
    CASE DEFAULT
       cmode = "stupid"  ! has to be one of the above
    END SELECT

    IF (PRESENT(define)) THEN
       ldefine_mode = define
    ELSE
       ldefine_mode = .FALSE.
    END IF

    ! Check if this file is already open, and just return that id
    ll_old_file = .FALSE.
    DO jid = ncid_min, ncid_max
       IF ( filename == dataset_filename(jid) ) THEN
          ll_old_file = .TRUE.
          ncid = jid
          IF ( mode == nf90_write ) THEN
             WRITE(*,*) "Warning: trying to open a file in write mode:",TRIM(filename)
             WRITE(*,*) "But file is already opened - I think this could be bad"
             ! Perhaps only bad if multiple opens in write mode (or plus create)
             ! Maybe should store mode in a structure with the filename?
          END IF
          EXIT
       END IF
    END DO

    IF ( .NOT. ll_old_file ) THEN
       stat = nf90_open( filename, mode, ncid )
       IF (stat /= nf90_noerr) THEN
          ! Problem opening file
          IF (PRESENT(exist) .AND. (stat == eNoEnt .OR. stat == nf90_eNotNc)) THEN
             exist = .FALSE.
          ELSE
             CALL handle_error( stat, "open_file", "opening file: ", &
                  TRIM(filename), "in "//TRIM(cmode)//" mode" )
          END IF

       ELSE
          ! File opened OK
          !dataset_filename(ncid) = filename !FC
          ncid_min = MIN( ncid, ncid_min )
          ncid_max = MAX( ncid, ncid_max )

          IF (nverbose >=2 ) THEN
             WRITE(*,FMT= '( 10x,"Dataset",i3," has been opened in ",a7," mode &
                  &for file:" )' ) ncid, cmode
             WRITE(*,FMT='(10x,a)') TRIM(filename)
             WRITE(*,*) ""
          END IF

       END IF
    END IF

    IF (ldefine_mode) THEN
       stat = nf90_redef( ncid )
       IF (stat /= nf90_noerr)   CALL handle_error( stat, "open_file", &
            "putting dataset into define mode", ncid=ncid )
    END IF

  END SUBROUTINE open_file


  SUBROUTINE create_file( filename, mode, ncid, title )
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: mode
    INTEGER, INTENT(OUT) :: ncid
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: title
    INTEGER :: stat
    CHARACTER(LEN=12) :: cmode

    SELECT CASE (mode)
    CASE( nf90_clobber )
       cmode = "clobber"
    CASE( nf90_noclobber )
       cmode = "noclobber"
    CASE( nf90_share )
       cmode = "share"
    CASE( nf90_64bit_offset )
       cmode = "64bit_offset"
    CASE DEFAULT
       cmode = "stupid"  ! has to be one of the above
    END SELECT

    stat = nf90_create( filename, mode, ncid )
    IF (stat /= nf90_noerr)   CALL handle_error( stat, "create_file", &
            "creating file: ", TRIM(filename), " in "//TRIM(cmode)//" mode" )

    !dataset_filename(ncid) = filename  !FC
    IF (nverbose >= 2) THEN
       WRITE(*,FMT= '( 10x,"Dataset",i3," has been created in ",a9," mode &
            &for file:" )' ) ncid, cmode
       WRITE(*,FMT='(10x,a)') TRIM(filename)
       WRITE(*,*) ""
    END IF

    ! created file ok - add standard global attributes
    stat = nf90_put_att ( ncid, nf90_global, "file_name", filename )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "create_file", &
         "putting global attribute: file_name="//filename, ncid=ncid )

    IF (PRESENT(title)) THEN
       stat = nf90_put_att( ncid, nf90_global, "title", title )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "create_file", &
            "putting global attribute: title="//title, ncid=ncid )
    END IF

  END SUBROUTINE create_file

  SUBROUTINE close_file( ncid )
    INTEGER, INTENT(INOUT) :: ncid
    INTEGER :: stat

    stat = nf90_close( ncid )
    SELECT CASE( stat )
    CASE( nf90_noerr )
       ! Closed with no errors
       IF (nverbose >= 2) WRITE(*,FMT='(10x,"Closed dataset",i3,/)') ncid
       !dataset_filename(ncid) = "" !FC
       ncid = -9 ! prevent ncid matching a file that is assigned this ID later
    CASE( nf90_ebadid )
       ! Not an open dataset - just be quiet
    CASE DEFAULT
       ! Some other error - call the handler
       CALL handle_error( stat, "close_file", "closing file", &
            ncid=ncid )
    END SELECT

  END SUBROUTINE close_file


  SUBROUTINE fill_dim_info ( dim, ncid, name, id, exist )
    ! fills in the missing details in an ncvartype structure
    ! needs the ncid and either the dimension id or name
    ! Behavior: if id is given in argument list this is used to find dim,
    !           else name is used (whether passed in args or structure)
    TYPE(ncdimtype), INTENT(INOUT) :: dim
    INTEGER, INTENT(IN), OPTIONAL :: ncid, id
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name
    LOGICAL, INTENT(OUT), OPTIONAL :: exist

    INTEGER :: stat, ndims, unlimid

    IF (PRESENT( exist )) exist = .TRUE.

    IF (PRESENT( ncid )) dim%ncid = ncid

    stat = nf90_inquire( dim%ncid, nDimensions=ndims, unlimitedDimID=unlimid )
    IF (stat /= nf90_noerr)  &
   &   CALL handle_error( stat, "fill_dim_info", "inquiring dataset", dim=dim )

    IF (PRESENT( id )) THEN
       dim%id = id
    ELSE
       ! get the dim id from the name
       IF (PRESENT(name)) dim%name = name
       stat = nf90_inq_dimid( dim%ncid, dim%name, dim%id )
       IF ( stat /= nf90_noerr ) dim%id = -1 ! ensure inq_dim fails
    END IF

    stat=nf90_inquire_dimension( dim%ncid, dim%id, name=dim%name, len=dim%len )

    IF ( stat == nf90_ebaddim .AND. PRESENT(exist) ) THEN
       exist = .FALSE.
    ELSE IF ( stat /= nf90_noerr ) THEN
       CALL handle_error( stat, "fill_dim_info", &
            "inquiring dim id or name and length", dim=dim )
    ELSE
       ! Dimension exists:
       ! is this the unlimitted dimension?
       dim%lunlim = (dim%id == unlimid)
    END IF

  END SUBROUTINE fill_dim_info


  SUBROUTINE fill_var_info ( var, ncid, name, id, exist )
    ! fills in the missing details in an ncvartype structure
    ! needs the ncid and either the dimension id or name
    ! Behavior: if id is given in argument list this is used to find var,
    !           else name is used (whether passed in args or structure)
    TYPE(ncvartype), INTENT(INOUT) :: var
    INTEGER, INTENT(IN), OPTIONAL :: ncid, id
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: name
    LOGICAL, INTENT(OUT), OPTIONAL :: exist

    TYPE(ncdimtype) :: tmpdim
    REAL :: zmdi
    INTEGER :: stat, nvars, jd
    LOGICAL :: ll_notvar

    IF (PRESENT( exist )) exist = .TRUE.

    IF (PRESENT( ncid )) var%ncid = ncid
    stat = nf90_inquire( var%ncid, nVariables=nvars )
    IF (stat /= nf90_noerr) CALL handle_error( stat, "fill_var_info", &
         "inquiring dataset", var=var )

    IF (PRESENT( id )) THEN
       var%id = id
    ELSE
       ! getting the var id from the name
       IF (PRESENT( name )) var%name = name
       stat = nf90_inq_varid( var%ncid, var%name, var%id )
       IF ( stat /= nf90_noerr ) var%id = -1 ! ensure inq_var fails
    END IF

    stat=nf90_inquire_variable( var%ncid, var%id, name=var%name, &
         xtype=var%xtype, ndims=var%ndims, dimids=var%dimids(:) )

    IF ( stat == nf90_enotvar .AND. PRESENT(exist) ) THEN
       exist = .FALSE.
    ELSE IF ( stat /= nf90_noerr ) THEN
       CALL handle_error( stat, "fill_var_info", &
            "inquiring var name, ndims and dimids", var=var )
    ELSE
       ! Variable exists:
       ! put in the dimension lengths for aiding memory allocation
       ! (put length 1 if dimension not used by variable)
       var%dimlens(:) = 1
       DO jd = 1, var%ndims
          CALL fill_dim_info( tmpdim, var%ncid, id=var%dimids(jd) )
          var%dimlens(jd) = tmpdim%len
       END DO

       ! Pick up the value of the missing data indicator
       stat = nf90_get_att( var%ncid, var%id, "_FillValue", zmdi )
       SELECT CASE ( stat )
       CASE ( nf90_noerr )
          var%mdi = zmdi
       CASE ( nf90_enotatt )
          ! This variable has no _FillValue attribute
          ! Perhaps issue warning, 
          ! or try missing_value (with appropriate change to copy_atts)
          var%mdi = mdi_unset
       CASE DEFAULT
          CALL handle_error( stat, "fill_var_info", &
               "getting var attribute: _FillValue", var=var )
       END SELECT
    END IF

    ! Allow a second method of checking if var exists: invalid id
    IF (PRESENT(exist)) THEN
       IF (.NOT.exist) var%id = -2
    END IF
    
  END SUBROUTINE fill_var_info



  SUBROUTINE def_dim ( dim, ncid, name, len )
    ! Defines the dimension if required and returns the dimid
    ! Required inputs are dim%ncid and dim%name
    TYPE (ncdimtype), INTENT(INOUT) :: dim
    INTEGER, OPTIONAL, INTENT(IN) :: ncid   !  ncid of dim
    INTEGER, OPTIONAL, INTENT(IN) :: len    !  length of dim
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: name   ! name of dim

    INTEGER :: dimid, stat, in_length

    IF (PRESENT(ncid)) dim%ncid = ncid
    IF (PRESENT(name)) dim%name = name
    IF (PRESENT(len )) dim%len  = len
    dim%id = -1

    stat = nf90_inq_dimid( dim%ncid, dim%name, dim%id )
    SELECT CASE (stat)
    CASE(nf90_noerr)
       ! Dimension of this name already exists:
       ! fill in dim info and check length agrees
       in_length = dim%len
       CALL fill_dim_info( dim )
       
       IF ( in_length /= dim%len ) THEN
          WRITE(*,*) "Error (def_dim): Dimension of same name exists already"
          WRITE(*,*) " but length requested is different:", in_length
          CALL print_dim( dim )
          STOP
       END IF

    CASE (nf90_ebaddim)
       ! dimension needs defined

       IF (dim%lunlim) THEN
          stat = nf90_def_dim( dim%ncid, dim%name, nf90_unlimited, dim%id )
       ELSE
          stat = nf90_def_dim( dim%ncid, dim%name, dim%len, dim%id )
       END IF
       IF (stat /= nf90_noerr) CALL handle_error( stat, "def_dim", &
            "defining dimension", dim=dim )

    CASE DEFAULT 
       ! There was an error inquiring the dim id
       CALL handle_error( stat, "def_dim", "inquiring dim id", dim=dim )
    END SELECT


  END SUBROUTINE def_dim


  SUBROUTINE def_var ( var, ncid, name, dims, units, mdi, type )
    ! Defines the variable and returns the varid
    ! Note that ndims is taken from the number of +ve dimids
    TYPE(ncvartype), INTENT(INOUT) :: var
    INTEGER, OPTIONAL, INTENT(IN) :: ncid            ! ncid for var
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: name   ! name of var
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: units
    TYPE(ncdimtype), DIMENSION(:), OPTIONAL, INTENT(IN) :: dims
    REAL(wp), OPTIONAL, INTENT(IN) :: mdi
    INTEGER, OPTIONAL, INTENT(IN) :: type

    INTEGER :: jd, nd, stat
    TYPE(ncdimtype) :: dim

    IF (PRESENT(ncid)) var%ncid  = ncid
    IF (PRESENT(name)) var%name  = name
    IF (PRESENT(mdi )) var%mdi   = mdi
    IF (PRESENT(type)) var%xtype = type
    IF (PRESENT(dims)) THEN
       var%dimids(:) = -1
       var%dimids(1:SIZE(dims)) = dims(:)%id
    END IF

    var%ndims = 0
    DO jd = 1, max_var_dims
       IF ( var%dimids(jd) > 0 ) THEN
          var%ndims = jd
       ELSE
          EXIT
       END IF
    END DO

    ! put into define mode (ignore error if already in define mode)
    stat = nf90_redef ( var%ncid )
    IF ((stat /= nf90_noerr).AND.(stat /= nf90_eindefine)) THEN
       CALL handle_error( stat, "def_var", &
            "entering define mode to define a variable", var=var )
    END IF

    nd = var%ndims
    IF (nd == 0) THEN
       ! define scalar variable
       stat = nf90_def_var( var%ncid, var%name, var%xtype, var%id )
    ELSE
       ! define array variable
       stat = nf90_def_var( var%ncid, var%name, var%xtype,  &
                            var%dimids(1:nd), var%id  )
    END IF
    IF (stat /= nf90_noerr) CALL handle_error( stat, "def_var", &
         "defining variable", var=var )

    ! add standard attributes
    IF (PRESENT(units)) THEN
       stat = nf90_put_att( var%ncid, var%id, "units", units )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "def_var", &
            "putting var attribute: units="//units, var=var )
    END IF

    IF ( var%mdi /= mdi_unset ) THEN
       stat = nf90_put_att( var%ncid, var%id, "_FillValue", var%mdi )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "def_var", &
            "putting var attribute: mdi=", var=var )
    END IF

    ! fill in dim lengths info (useful for allocating)
    var%dimlens(:) = 1
    DO jd = 1, var%ndims
       CALL fill_dim_info( dim, ncid=var%ncid, id=var%dimids(jd) )
       var%dimlens(jd) = dim%len
    END DO

    
  END SUBROUTINE def_var




  ! To do: alter copy_atts so that it automatically copies global attributes
  ! if they don't already exist

  SUBROUTINE copy_atts ( src, dst, nogrid )
    ! copy attributes from src to dst
    ! initially copy only text, as numerical may no longer be valid
    ! later consider doing something with the numerical ones


    !could do global atts by passing in vars with ids = nf90_global
    ! would reuire a few extra if tests (copy all atts for global, history)
    TYPE(ncvartype), INTENT(IN) :: src, dst
    LOGICAL, OPTIONAL, INTENT(IN) :: nogrid
    INTEGER, PARAMETER :: max_hist = 1000
    CHARACTER(LEN=1), PARAMETER :: newline = ACHAR(10)
    INTEGER :: natts, attid, atttype, lh, i, stat
    CHARACTER(LEN=nf90_max_name) :: attname, history_append
    CHARACTER(LEN=max_hist) :: history
    CHARACTER(LEN=8) :: date
    CHARACTER(LEN=10) :: time
    CHARACTER(LEN=5) :: zone
    LOGICAL :: add_grid

    IF (PRESENT(nogrid)) THEN
       add_grid = .NOT.nogrid
    ELSE
       add_grid = .TRUE.
    END IF
    

    IF (src%id == nf90_global) THEN
       stat = nf90_inquire( src%ncid, nAttributes=natts )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
            "inquiring dataset", ncid=src%ncid )
    ELSE
       stat = nf90_inquire_variable( src%ncid, src%id, nAtts=natts )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
            "inquiring variable", var=src )
    END IF
    

    DO attid = 1, natts
       ! get attribute type (atts accessed by name, so need this first)
       stat = nf90_inq_attname ( src%ncid, src%id, attid, attname )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
            "inquiring attribute name", var=src )

       
       SELECT CASE( attname )
       CASE ( "grid", "history", "file_name", "coordinates", "associate" )
          ! ignore these attrubutes
          ! we may want to append to "history" in future
          
       CASE DEFAULT
          ! copy this attributes
          stat = nf90_inquire_attribute ( src%ncid, src%id, attname, &
               xtype=atttype )
          IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
               "inquiring attribute", var=src )

          IF ( (atttype == nf90_char)         .OR. &
               (attname == "_FillValue")        ) THEN
             stat = nf90_copy_att ( src%ncid, src%id, attname, dst%ncid, dst%id )
             IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
                  "copying attributes from var to var2", var=src, var2=dst )
          END IF

       END SELECT

    ENDDO

    ! add new attributes (global and variable)
    ! file_name attribute is added when file is created
    IF (src%id == nf90_global) THEN
       ! history:
       CALL DATE_AND_TIME ( date=date, time=time, zone=zone )
       history_append = &
           date(7:8) // "/" // date(5:6) // "/" // date(1:4) // " " // &
           time(1:2) // ":" // time(3:4) // ":" // time(5:6) // " " // &
           zone // ": scrip_remap (see variables for grid info)"
       
       history = ""
       stat = nf90_get_att ( src%ncid, nf90_global, "history", history )
       SELECT CASE( stat )
          CASE( nf90_noerr, nf90_enotatt )
             ! either found attribute or attribute not found
             lh = LEN_TRIM(history)
             IF (lh == 0) THEN
                history = history_append
             ELSE
                history = history(1:lh) // newline // history_append
             END IF
             
             stat = nf90_put_att ( dst%ncid, nf90_global, "history", history )
             IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
                  "putting history string global attribute", ncid=dst%ncid )
             
          CASE DEFAULT
             ! There was an error looking for the attribute
             CALL handle_error( stat, "copy_atts", &
                  "getting history string global attribute", ncid=src%ncid )
          END SELECT

    ELSE ! attributes are non-global

       IF (add_grid .AND. dst%grid /= "") THEN
          stat = nf90_put_att( dst%ncid, dst%id, "grid", dst%grid )
          IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_atts", &
               "putting grid attribute", var=dst )
       END IF

    END IF

  END SUBROUTINE copy_atts


  RECURSIVE SUBROUTINE copy_var( src, dims, dst_name, dst_ncid )
    ! Define a variable with dims and copy atts and data from a source var.
    ! The destination dataset is taken to be that of the first dimension
    ! unless dst_ncid is provided, in which case dims are ignored and 
    ! dim names and lengths taken from src file.
    ! Dest var name taken from src unless optional dst_name is provided.
    ! Associated, coordinate variables also copied, hence routine recursive
    TYPE(ncvartype), INTENT(IN) :: src
    TYPE(ncdimtype), DIMENSION(:), OPTIONAL, INTENT(IN) :: dims
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: dst_name
    INTEGER, OPTIONAL, INTENT(IN) :: dst_ncid

    INTEGER, PARAMETER :: max_string_data = 256

    TYPE(ncvartype) :: dst, srccoord
    TYPE(ncdimtype), ALLOCATABLE, DIMENSION(:) :: dst_dims
    TYPE(ncdimtype) :: srcdim, dstdim
    REAL(dp), ALLOCATABLE, DIMENSION(:)       :: vardta1
    REAL(dp), ALLOCATABLE, DIMENSION(:,:)     :: vardta2
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:)   :: vardta3
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:) :: vardta4
    CHARACTER(LEN=max_string_data) :: varstr1
    CHARACTER(LEN=max_string_data), ALLOCATABLE, DIMENSION(:) :: varstr2
    INTEGER :: ncid, stat, jd, nd, n1, n2, n3, n4
    CHARACTER(LEN=nf90_max_name) :: name
    LOGICAL :: ll_var_exists, ll_bad, ll_exist


    nd = src%ndims
    ALLOCATE( dst_dims( nd ) )

    IF ( PRESENT(dst_ncid) ) THEN

       ncid = dst_ncid
       ! get source dim info and use name and length to define dest dim
       DO jd = 1, nd
          CALL fill_dim_info( dst_dims(jd), src%ncid, id=src%dimids(jd) )
          CALL def_dim( dst_dims(jd), ncid )
       END DO

    ELSE IF ( PRESENT(dims) ) THEN

       ncid = dims(1)%ncid

       ! Check that dst_dims consistent with src.
       ll_bad = .FALSE.
       IF ( SIZE(dims) == nd ) THEN
          dst_dims = dims
       ELSE IF ( src%xtype == nf90_char .AND. nd == SIZE(dims) + 1 ) THEN
          ! String length dim not included in dims(:) - add it in as 1st dim
          ! Get dim name and length from source and define in dest if not exist

          CALL fill_dim_info( srcdim, src%ncid, id=src%dimids(1) )
          CALL fill_dim_info( dst_dims(1), ncid, name=srcdim%name, exist=ll_exist )

          IF (.NOT. ll_exist) CALL def_dim( dst_dims(1), ncid, name=srcdim%name, &
                                            len=srcdim%len )

          dst_dims(2:) = dims(:)
       ELSE
          ll_bad = .TRUE.
       END IF

       IF (.NOT. ll_bad) THEN
          DO jd = 1, nd
             IF ( dst_dims(jd)%len /= src%dimlens(jd) ) ll_bad = .TRUE.
             IF ( dst_dims(jd)%id < 1 ) ll_bad = .TRUE.
          END DO
       END IF

       IF ( ll_bad ) THEN
          WRITE(*,*) "Error (copy var): dims bad or don't match source var shape"
          WRITE(*,*) "Source var dimlens array:  ", src%dimlens(:nd)
          WRITE(*,*) "Lengths of dest dims given:", dst_dims(:)%len
          WRITE(*,fmt='("Source variable...")')
          CALL print_var( src )
          DO jd = 1, nd
             WRITE(*,fmt='("Dimension ",i2,"...")') jd
             CALL print_dim( dst_dims(jd) )
          END DO
          STOP
       END IF


    ELSE
       WRITE(*,*) "Error (copy_var): either dims or dst_ncid must be supplied"
       STOP  
    END IF

    dst%ncid = ncid
    IF (PRESENT(dst_name)) THEN
       dst%name = dst_name
    ELSE
       dst%name = src%name
    END IF
    ! Check if var exists before defining it
    CALL fill_var_info( dst, exist=ll_var_exists )

    IF ( ll_var_exists ) THEN
       ! Variable already exists - do nothing
       ! In future, could add an option to fail at this point
       !  or check that this is the same var (from dims & data)
    ELSE

       ! Define dest variable and copy atts and data from source
       CALL def_var( dst, dims=dst_dims, type=src%xtype )
       CALL copy_atts( src, dst, nogrid=.TRUE. )
       DO jd = 1, nd
          dstdim = dst_dims(jd)
          CALL fill_dim_info( srcdim, src%ncid, id=src%dimids(jd) )
          srccoord = find_coordinate( srcdim )
          IF (srccoord%id > 0) CALL copy_var( srccoord, (/dstdim/), dstdim%name )
       END DO

       n1 = src%dimlens(1)
       n2 = src%dimlens(2)
       n3 = src%dimlens(3)
       n4 = src%dimlens(4)
       SELECT CASE ( src%xtype )
       CASE( NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE )
          ! All these can be represented accurately by a REAL(8) variable
          SELECT CASE ( src%ndims )
          CASE( 1 )
             ALLOCATE( vardta1( n1 ) )
             CALL get_var( src, vardta1 )
             CALL put_var( dst, vardta1 )
             DEALLOCATE( vardta1 )
          CASE( 2 )
             ALLOCATE( vardta2( n1, n2 ) )
             CALL get_var( src, vardta2 )
             CALL put_var( dst, vardta2 )
             DEALLOCATE( vardta2 )
          CASE( 3 )
             ALLOCATE( vardta3( n1, n2, n3 ) )
             CALL get_var( src, vardta3 )
             CALL put_var( dst, vardta3 )
             DEALLOCATE( vardta3 )
          CASE( 4 )
             ALLOCATE( vardta4( n1, n2, n3, n4 ) )
             CALL get_var( src, vardta4 )
             CALL put_var( dst, vardta4 )
             DEALLOCATE( vardta4 )
          END SELECT
       CASE( NF90_CHAR )
          SELECT CASE( src%ndims ) ! 1st dim in char var is string length
          CASE( 1 )
             CALL get_var( src, varstr1(1:n1) )
             CALL put_var( src, varstr1(1:n1) )
          CASE( 2 )
             ALLOCATE( varstr2( n2 ) )
             CALL get_var( src, varstr2(:)(1:n1) )
             CALL put_var( dst, varstr2(:)(1:n1) )
             DEALLOCATE( varstr2 )
          END SELECT
       END SELECT
             

       ! put_var switches to data mode - switch back to define mode
       stat = nf90_redef( dst%ncid )
       IF (stat /= nf90_noerr) CALL handle_error( stat, "copy_var",&
            "Putting back into define mode", ncid=dst%ncid )

    END IF ! var exists

  END SUBROUTINE copy_var


  FUNCTION find_coordinate ( dim )
    ! Look for a coordinate variable for dim
    ! It must be 1D and have same name as dim
    ! We assume that dim%ncid, dim%id and dim%name are correctly set
    ! If no coordinate is found, return coord%id = -3

    ! In future, could find coords using CF conventions ("coordinates" att?)
    ! They might have more than 1 dimension, e.g. nav_lon/lat

    TYPE(ncvartype) :: find_coordinate
    TYPE(ncdimtype), INTENT(IN) :: dim

    TYPE(ncvartype) :: coord
    INTEGER :: varid, dimid, vartype, ndims, stat
    LOGICAL :: ll_exist, ll_coord
    
    ! Firstly check if a variable of the same name exists
    ! If it does, then see if it is 1D and is dimensioned on the 
    !   dimension of the same name
    ! If we find it, define a new 1D variable with same name as dest dim
    !   and copy attributes and data from source variable
    ! Warning: putting into data mode and back again repeatedly may
    ! result in a needless file copy.  We haven't put the main data in
    ! yet, so this should be cheap.
    
    CALL fill_var_info( coord, dim%ncid, dim%name, exist=ll_exist )

    ll_coord = .FALSE.
    IF( ll_exist ) THEN
       IF ( coord%xtype /= nf90_char ) THEN
          IF ( coord%ndims == 1 .AND. coord%dimids(1) == dim%id ) ll_coord = .TRUE.
       ELSE
          ! string variable - 1st dim is string len, 2nd is dimension
          IF ( coord%ndims == 2 .AND. coord%dimids(2) == dim%id ) ll_coord = .TRUE.
       END IF

    END IF
       
    IF (.NOT. ll_coord) coord%id = -3

    find_coordinate = coord

  END FUNCTION find_coordinate


  FUNCTION dims_equal( dim1, dim2 )
    ! For the dimensions to be equal they must have the same length
    ! If they both have coordinates, then in addition, the coords be be equal
    LOGICAL :: dims_equal
    TYPE(ncdimtype), INTENT(IN) :: dim1, dim2

    REAL(wp), PARAMETER :: pptol = 1.0e-3  ! relative difference allowed
                                           ! between coordinate values
    TYPE(ncvartype) :: coord1, coord2
    REAL(wp), DIMENSION(dim1%len) :: zdata1, zdata2, zdiff
    REAL(wp) :: zdiffmax

    dims_equal = .TRUE.
    IF ( dim1%len /= dim2%len ) THEN
       dims_equal = .FALSE.
    ELSE
       coord1 = find_coordinate( dim1 )
       coord2 = find_coordinate( dim2 )
       IF ( coord1%id > 0 .AND. coord2%id > 0 ) THEN
          CALL get_var( coord1, zdata1 )
          CALL get_var( coord2, zdata2 )
          
          zdiff(:) = ( zdata1(:) - zdata2(:) ) / zdata1(:)
          zdiffmax = MAXVAL( ABS( zdiff(:) )  )
          IF ( zdiffmax > pptol ) dims_equal = .FALSE.
       END IF

    END IF


  END FUNCTION dims_equal


  SUBROUTINE get_var_any( var, pvr0, pvr1, pvr2, pvr3, pvr4, &
                               pvd0, pvd1, pvd2, pvd3, pvd4, &
                               pvi0, pvi1, pvi2, pvi3, pvi4, &
                               pvc1, pvc2, &
                               start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp),                     OPTIONAL, INTENT(OUT) :: pvr0
    REAL(sp), DIMENSION(:),       OPTIONAL, INTENT(OUT) :: pvr1
    REAL(sp), DIMENSION(:,:),     OPTIONAL, INTENT(OUT) :: pvr2
    REAL(sp), DIMENSION(:,:,:),   OPTIONAL, INTENT(OUT) :: pvr3
    REAL(sp), DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT) :: pvr4
    REAL(dp),                     OPTIONAL, INTENT(OUT) :: pvd0
    REAL(dp), DIMENSION(:),       OPTIONAL, INTENT(OUT) :: pvd1
    REAL(dp), DIMENSION(:,:),     OPTIONAL, INTENT(OUT) :: pvd2
    REAL(dp), DIMENSION(:,:,:),   OPTIONAL, INTENT(OUT) :: pvd3
    REAL(dp), DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT) :: pvd4
    INTEGER,                      OPTIONAL, INTENT(OUT) :: pvi0
    INTEGER,  DIMENSION(:),       OPTIONAL, INTENT(OUT) :: pvi1
    INTEGER,  DIMENSION(:,:),     OPTIONAL, INTENT(OUT) :: pvi2
    INTEGER,  DIMENSION(:,:,:),   OPTIONAL, INTENT(OUT) :: pvi3
    INTEGER,  DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT) :: pvi4
    CHARACTER(LEN=*),             OPTIONAL, INTENT(OUT) :: pvc1
    CHARACTER(LEN=*),DIMENSION(:),OPTIONAL, INTENT(OUT) :: pvc2
    INTEGER,  DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map

    INTEGER :: stat, ncid, vid
    CHARACTER(LEN=*), PARAMETER :: cl_routine="get_var_any"

    ncid = var%ncid
    vid  = var%id

    ! put into data mode (ignore error if not in define mode)
    stat = nf90_enddef ( ncid )
    IF ((stat /= nf90_noerr).AND.(stat /= nf90_enotindefine)) THEN
       CALL handle_error( stat, cl_routine, &
            "entering data mode to put a variable", var=var )
    END IF
    
    IF (PRESENT(pvr0)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvr0 )
    ELSE IF (PRESENT(pvr1)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvr1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr2)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvr2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr3)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvr3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr4)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvr4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd0)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvd0 )
    ELSE IF (PRESENT(pvd1)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvd1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd2)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvd2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd3)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvd3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd4)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvd4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi0)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvi0 )
    ELSE IF (PRESENT(pvi1)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvi1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi2)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvi2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi3)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvi3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi4)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvi4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvc1)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvc1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvc2)) THEN
       stat = nf90_get_var ( var%ncid, var%id, pvc2, &
                             start=start, count=count, stride=stride, map=map )
    END IF

    SELECT CASE ( stat )
    CASE ( nf90_noerr )
       ! Success
    CASE ( nf90_einvalcoords, nf90_eedge )
       ! Something wrong with the dimension lengths
       WRITE(*,*)
       WRITE(*,*) "NetCDF error: ",TRIM(nf90_strerror(stat))
       WRITE(*,*) "in routine ",cl_routine

       IF      (PRESENT(pvr0)) THEN ; CALL print_dimlens( var, SHAPE(pvr0) )
       ELSE IF (PRESENT(pvr1)) THEN ; CALL print_dimlens( var, SHAPE(pvr1) )
       ELSE IF (PRESENT(pvr2)) THEN ; CALL print_dimlens( var, SHAPE(pvr2) )
       ELSE IF (PRESENT(pvr3)) THEN ; CALL print_dimlens( var, SHAPE(pvr3) )
       ELSE IF (PRESENT(pvr4)) THEN ; CALL print_dimlens( var, SHAPE(pvr4) )
       ELSE IF (PRESENT(pvd0)) THEN ; CALL print_dimlens( var, SHAPE(pvd0) )
       ELSE IF (PRESENT(pvd1)) THEN ; CALL print_dimlens( var, SHAPE(pvd1) )
       ELSE IF (PRESENT(pvd2)) THEN ; CALL print_dimlens( var, SHAPE(pvd2) )
       ELSE IF (PRESENT(pvd3)) THEN ; CALL print_dimlens( var, SHAPE(pvd3) )
       ELSE IF (PRESENT(pvd4)) THEN ; CALL print_dimlens( var, SHAPE(pvd4) )
       ELSE IF (PRESENT(pvi0)) THEN ; CALL print_dimlens( var, SHAPE(pvi0) )
       ELSE IF (PRESENT(pvi1)) THEN ; CALL print_dimlens( var, SHAPE(pvi1) )
       ELSE IF (PRESENT(pvi2)) THEN ; CALL print_dimlens( var, SHAPE(pvi2) )
       ELSE IF (PRESENT(pvi3)) THEN ; CALL print_dimlens( var, SHAPE(pvi3) )
       ELSE IF (PRESENT(pvi4)) THEN ; CALL print_dimlens( var, SHAPE(pvi4) )
       ELSE IF (PRESENT(pvc1)) THEN ; CALL print_dimlens( var, SHAPE(pvc1), LEN(pvc1) )
       ELSE IF (PRESENT(pvc2)) THEN ; CALL print_dimlens( var, SHAPE(pvc2), LEN(pvc2) )
       END IF

       STOP 

    CASE DEFAULT
       CALL handle_error( stat, cl_routine, "getting variable", var=var )
    END SELECT

  END SUBROUTINE get_var_any


  SUBROUTINE put_var_any( var, pvr0, pvr1, pvr2, pvr3, pvr4, &
                               pvd0, pvd1, pvd2, pvd3, pvd4, &
                               pvi0, pvi1, pvi2, pvi3, pvi4, &
                               pvc1, pvc2, &
                               start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp),                     OPTIONAL, INTENT(IN) :: pvr0
    REAL(sp), DIMENSION(:),       OPTIONAL, INTENT(IN) :: pvr1
    REAL(sp), DIMENSION(:,:),     OPTIONAL, INTENT(IN) :: pvr2
    REAL(sp), DIMENSION(:,:,:),   OPTIONAL, INTENT(IN) :: pvr3
    REAL(sp), DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: pvr4
    REAL(dp),                     OPTIONAL, INTENT(IN) :: pvd0
    REAL(dp), DIMENSION(:),       OPTIONAL, INTENT(IN) :: pvd1
    REAL(dp), DIMENSION(:,:),     OPTIONAL, INTENT(IN) :: pvd2
    REAL(dp), DIMENSION(:,:,:),   OPTIONAL, INTENT(IN) :: pvd3
    REAL(dp), DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: pvd4
    INTEGER,                      OPTIONAL, INTENT(IN) :: pvi0
    INTEGER,  DIMENSION(:),       OPTIONAL, INTENT(IN) :: pvi1
    INTEGER,  DIMENSION(:,:),     OPTIONAL, INTENT(IN) :: pvi2
    INTEGER,  DIMENSION(:,:,:),   OPTIONAL, INTENT(IN) :: pvi3
    INTEGER,  DIMENSION(:,:,:,:), OPTIONAL, INTENT(IN) :: pvi4
    CHARACTER(LEN=*),             OPTIONAL, INTENT(IN) :: pvc1
    CHARACTER(LEN=*),DIMENSION(:),OPTIONAL, INTENT(IN) :: pvc2
    INTEGER,  DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map

    INTEGER :: stat, ncid, vid
    CHARACTER(LEN=*), PARAMETER :: cl_routine="put_var_any"

    ncid = var%ncid
    vid  = var%id

    ! put into data mode (ignore error if not in define mode)
    stat = nf90_enddef ( ncid )
    IF ((stat /= nf90_noerr).AND.(stat /= nf90_enotindefine)) THEN
       CALL handle_error( stat, cl_routine, &
            "entering data mode to put a variable", var=var )
    END IF

    IF (PRESENT(pvr0)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvr0 )
    ELSE IF (PRESENT(pvr1)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvr1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr2)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvr2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr3)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvr3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvr4)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvr4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd0)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvd0 )
    ELSE IF (PRESENT(pvd1)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvd1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd2)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvd2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd3)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvd3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvd4)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvd4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi0)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvi0 )
    ELSE IF (PRESENT(pvi1)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvi1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi2)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvi2, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi3)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvi3, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvi4)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvi4, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvc1)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvc1, &
                             start=start, count=count, stride=stride, map=map )
    ELSE IF (PRESENT(pvc2)) THEN
       stat = nf90_put_var ( var%ncid, var%id, pvc2, &
                             start=start, count=count, stride=stride, map=map )
    END IF

    SELECT CASE ( stat )
    CASE ( nf90_noerr )
       ! Success
    CASE ( nf90_einvalcoords, nf90_eedge )
       ! Something wrong with the dimension lengths
       WRITE(*,*)
       WRITE(*,*) "NetCDF error: ",TRIM(nf90_strerror(stat))
       WRITE(*,*) "in routine ",cl_routine

       IF      (PRESENT(pvr0)) THEN ; CALL print_dimlens( var, SHAPE(pvr0) )
       ELSE IF (PRESENT(pvr1)) THEN ; CALL print_dimlens( var, SHAPE(pvr1) )
       ELSE IF (PRESENT(pvr2)) THEN ; CALL print_dimlens( var, SHAPE(pvr2) )
       ELSE IF (PRESENT(pvr3)) THEN ; CALL print_dimlens( var, SHAPE(pvr3) )
       ELSE IF (PRESENT(pvr4)) THEN ; CALL print_dimlens( var, SHAPE(pvr4) )
       ELSE IF (PRESENT(pvd0)) THEN ; CALL print_dimlens( var, SHAPE(pvd0) )
       ELSE IF (PRESENT(pvd1)) THEN ; CALL print_dimlens( var, SHAPE(pvd1) )
       ELSE IF (PRESENT(pvd2)) THEN ; CALL print_dimlens( var, SHAPE(pvd2) )
       ELSE IF (PRESENT(pvd3)) THEN ; CALL print_dimlens( var, SHAPE(pvd3) )
       ELSE IF (PRESENT(pvd4)) THEN ; CALL print_dimlens( var, SHAPE(pvd4) )
       ELSE IF (PRESENT(pvi0)) THEN ; CALL print_dimlens( var, SHAPE(pvi0) )
       ELSE IF (PRESENT(pvi1)) THEN ; CALL print_dimlens( var, SHAPE(pvi1) )
       ELSE IF (PRESENT(pvi2)) THEN ; CALL print_dimlens( var, SHAPE(pvi2) )
       ELSE IF (PRESENT(pvi3)) THEN ; CALL print_dimlens( var, SHAPE(pvi3) )
       ELSE IF (PRESENT(pvc1)) THEN ; CALL print_dimlens( var, SHAPE(pvc1), LEN(pvc1) )
       ELSE IF (PRESENT(pvc2)) THEN ; CALL print_dimlens( var, SHAPE(pvc2), LEN(pvc2) )
       END IF

       STOP 

    CASE DEFAULT
       CALL handle_error( stat, cl_routine, "getting variable", var=var )
    END SELECT

  END SUBROUTINE put_var_any



  SUBROUTINE quickdump( name, data )
    ! For debugging purposes, puts data into a file
    ! just 2d for now
    CHARACTER(LEN=*), INTENT(IN) :: name ! varname and filename
    REAL(wp), DIMENSION(:,:) :: data

    REAL(wp), DIMENSION( SIZE(data,1), SIZE(data,2) ) :: zdata
    INTEGER :: rank, i
    TYPE(ncvartype) :: var
    TYPE(ncdimtype) :: dim

    CALL create_file( name//".nc", nf90_clobber, var%ncid )
    dim%ncid = var%ncid

    rank = SIZE(SHAPE(data))
    DO i = 1, rank
       WRITE( dim%name, FMT='("dim",i1)' ) i
       dim%len = SIZE(data,i)
       CALL def_dim( dim )
       var%dimids(i) = dim%id
       var%dimlens(i) = dim%len
    END DO

    var%name = name
    var%ndims = rank
    CALL def_var( var )

    zdata = data  ! this allows 'data' to be an array-subsection
    CALL put_var( var, zdata )

    CALL close_file( var%ncid )
    
  END SUBROUTINE quickdump


  ! Overloads for get_var.  All of them call get_var_any

  SUBROUTINE get_var_0d_real( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), INTENT(OUT) :: var_dta
    CALL get_var_any( var, pvr0=var_dta )
  END SUBROUTINE get_var_0d_real

  SUBROUTINE get_var_1d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvr1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_1d_real

  SUBROUTINE get_var_2d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvr2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_2d_real

  SUBROUTINE get_var_3d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvr3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_3d_real

  SUBROUTINE get_var_4d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvr4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_4d_real

  SUBROUTINE get_var_0d_double( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), INTENT(OUT) :: var_dta
    CALL get_var_any( var, pvd0=var_dta )
  END SUBROUTINE get_var_0d_double

  SUBROUTINE get_var_1d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvd1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_1d_double

  SUBROUTINE get_var_2d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvd2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_2d_double

  SUBROUTINE get_var_3d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvd3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_3d_double

  SUBROUTINE get_var_4d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvd4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_4d_double

  SUBROUTINE get_var_0d_int( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, INTENT(OUT) :: var_dta
    CALL get_var_any( var, pvi0=var_dta )
  END SUBROUTINE get_var_0d_int

  SUBROUTINE get_var_1d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvi1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_1d_int

  SUBROUTINE get_var_2d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvi2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_2d_int

  SUBROUTINE get_var_3d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvi3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_3d_int

  SUBROUTINE get_var_4d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:,:,:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvi4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_4d_int

  ! Character arrays:
  !  0d = single character = string length 1 => handled by 1d case
  !  1d char array = a string
  !  2d char array = vector of strings
  !  etc...
  !  For 2d and above, if calling arg is a vector of substrings, the arg is
  ! returned as if the full strings had been passed.  This differs from the
  ! behaviour for numerical arrays.  My solution is to copy the data into a
  ! local array of the correct shape.  This will be wasteful of memory for very
  ! large arrays.  The alternative is to use the map argument, but I can't be
  ! bothered.
 
  SUBROUTINE get_var_1d_char( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    CHARACTER(LEN=*), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL get_var_any( var,pvc1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE get_var_1d_char

  SUBROUTINE get_var_2d_char( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map

    CHARACTER(LEN=LEN(var_dta)), DIMENSION(SIZE(var_dta)) :: tmp_dta

    CALL get_var_any( var,pvc2=tmp_dta,start=start,count=count,stride=stride,map=map )
    var_dta = tmp_dta

  END SUBROUTINE get_var_2d_char



  ! Overloads for put_var.  All of them call put_var_any

  SUBROUTINE put_var_0d_real( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), INTENT(IN) :: var_dta
    CALL put_var_any( var, pvr0=var_dta )
  END SUBROUTINE put_var_0d_real

  SUBROUTINE put_var_1d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvr1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_1d_real

  SUBROUTINE put_var_2d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvr2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_2d_real

  SUBROUTINE put_var_3d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvr3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_3d_real

  SUBROUTINE put_var_4d_real( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(sp), DIMENSION(:,:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvr4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_4d_real

  SUBROUTINE put_var_0d_double( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), INTENT(IN) :: var_dta
    CALL put_var_any( var, pvd0=var_dta )
  END SUBROUTINE put_var_0d_double

  SUBROUTINE put_var_1d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvd1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_1d_double

  SUBROUTINE put_var_2d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvd2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_2d_double

  SUBROUTINE put_var_3d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvd3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_3d_double

  SUBROUTINE put_var_4d_double( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    REAL(dp), DIMENSION(:,:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvd4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_4d_double

  SUBROUTINE put_var_0d_int( var, var_dta )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, INTENT(IN) :: var_dta
    CALL put_var_any( var, pvi0=var_dta )
  END SUBROUTINE put_var_0d_int

  SUBROUTINE put_var_1d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvi1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_1d_int

  SUBROUTINE put_var_2d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvi2=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_2d_int

  SUBROUTINE put_var_3d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvi3=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_3d_int

  SUBROUTINE put_var_4d_int( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    INTEGER, DIMENSION(:,:,:,:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvi4=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_4d_int

  ! Character arrays: see notes for get_var_1d_char etc

  SUBROUTINE put_var_1d_char( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    CHARACTER(LEN=*), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map
    CALL put_var_any( var,pvc1=var_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_1d_char

  SUBROUTINE put_var_2d_char( var, var_dta, start, count, stride, map )
    TYPE (ncvartype), INTENT(IN) :: var
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: var_dta
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count, stride, map

    CHARACTER(LEN=LEN(var_dta)), DIMENSION(SIZE(var_dta)) :: tmp_dta

    tmp_dta = var_dta
    CALL put_var_any( var,pvc2=tmp_dta,start=start,count=count,stride=stride,map=map )
  END SUBROUTINE put_var_2d_char

END MODULE nf90util

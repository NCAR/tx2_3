!=======================================================================
!     calculate grid in the tropics using a stretch function for
!     increased resolution.
!
!     Either call this routine from mshglo after the grid is complete,
!     or call from a seperate program if memory is tight.
! 
!     Everything done in degrees unless otherwise stated
!---------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------

MODULE mshtrp
  USE kinds
  USE nf90util
  USE shoots, ONLY: newton
  USE trop
  USE common, ONLY: glam, gphi, e1,e2
  USE param, ONLY: jpeqt, jpim, jpjm, rad
  USE functions, ONLY:fsdila
  
  IMPLICIT NONE

CONTAINS


  SUBROUTINE msh_trp

    INTEGER :: jj, ijtarg, allostat, ijntv, ijm, iphic
    INTEGER :: ijgphic, ijtrop, ijgtrop, ijgcanc, ijgcapr, ijgnth
    INTEGER :: iig, iis, ijsnth, ijscapr, ijseqt, ijscanc
    INTEGER :: con, iters, jg, jd, jv

    REAL(wp) :: zphic, zc, zdj, zerr, ztarg, wgt, zlat, zlatp

    ! The tropics grid size (and hence whole grid size) is unknown
    ! until this routine, so it is simplest to write the netcdf file here
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: zphitrop, ze1trop, ze2trop
    REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwlam, zwphi, zwe1, zwe2
    LOGICAL :: llfound

    WRITE(*,*) ""
    WRITE(*,*) "***************************"
    WRITE(*,*) " msh_trp: calcuting tropical grid"

    ! glam simply gets extra rows
    ! gphi, e1 and e2 must be recalculated

    ! Identify the lat range requested and check that phi,e2 sufficently
    ! symmetrical and that lon and e1 are suffiently constant.
    ! Store closest phi to phi_join (=phic) and its index ijgcanc.
    !     
    IF (MAXVAL(ABS(gphi(:,jpeqt,1:2))) > rzero) THEN
       WRITE(*,*) "Error: Equator latitude not zero"
       WRITE(*,*) "max lat, zero:",MAXVAL(ABS(gphi(:,jpeqt,1:2))), rzero
       STOP
    END IF

    ijm = SIZE(glam,2)
    zlat = 0.0
    llfound=.FALSE.
    DO jj = jpeqt, ijm
       zlatp = zlat
       zlat = gphi(1,jj,1)
       IF (zlat >= pphi_join) THEN
          llfound=.TRUE.
          IF ( zlat-pphi_join < pphi_join-zlatp ) THEN
             zphic = zlat
             ijgcanc = jj
          ELSE
             zphic = zlatp
             ijgcanc = jj-1
          END IF
          EXIT
       END IF
    END DO
    IF (llfound) THEN
       WRITE(*,*) "Requested latitude to join tropical grid:",pphi_join
       WRITE(*,*) "Closest latitude in input grid, phic=",zphic
    ELSE
       STOP "Didn't find phi_join" 
    END IF

    !     Use phi as the independent variable and j the dependent variable
    !     Generate high res grid of phi, with phic as final point.
    !     Take recip of input function dphi/dj(phi) to generate dj/dphi(phi).
    !     Integrate to get j(phi).  
    !     Store j(phic) and the nearest integer (=jc).

    nhp = NINT( REAL(jptint,KIND=wp) * zphic / ppresmin )

    ALLOCATE( hphi(nhp), hj(nhp) ) !, zdjdphi(nhp)
    hphi = (/ ( zphic * REAL(jj,KIND=wp)/REAL(nhp-1), jj=0,nhp-1 ) /)


    zc = 1.0

    CALL integ ( zc, zerr, ld_set_targ=.TRUE.)

    WRITE(*,*) ""
    WRITE(*,&
         FMT='(" With c=",f7.4,", nearest integer is",i4,", mismatch is",e14.6)' &
         ) zc, njc, zerr
    !     Alter multiplier of phi in dphi/dj and iterate until j(phic)=jc.
    !     
    CALL newton ( zc, zerr, con, iters )
    IF (con > 0) THEN
       write(*,*) ""
       WRITE(*,FMT='(" Found root in ",i4," iters (con=",i1,")")') iters, con
       WRITE(*,FMT='(" Final mismatch in j:",e12.4)') zerr
       WRITE(*,FMT='(" Value of c used in function:",f14.8)') zc
    ELSE
       WRITE(*,*) "Error, failed to find root, con = ",con
       STOP
    END IF

    ! Now that correct function has been found and j(phi) calculated on
    ! high res phi grid, get values of phi at each half-integer value
    ! of j (for T-U and V-F grids) and evaluate e2=ra*dphi/dj at each
    ! of these lats.  Tricky bit done. 
    !
    ijntv = 2*(njc) - 1   ! npoints:  njc T points + njc-1 V points
    ALLOCATE( zphitrop(ijntv), ze1trop(ijntv), ze2trop(ijntv), STAT=allostat )
    IF (allostat /= 0) STOP "Err allocating zphitrop, ze1trop, ze2trop"

    zphitrop(1) = 0.0
    ijtarg = 2
    ztarg = 1.5 
    jj = 2
    DO WHILE ((jj <= nhp-1) .AND. (INT(ztarg) < njc))  !!!! not
!!!! quite right if ztarg is 49.999, but NINT wrong if 49.5 or 49.50001
!!!! - could multiply by 2 and do NINT then integer-divide by 2 to
!!!! force round-down

       IF (hj(jj) > ztarg) THEN
          wgt = ( ztarg - hj(jj-1) ) / ( hj(jj) - hj(jj-1) )
          zphitrop(ijtarg) = (1.0_wp - wgt) * hphi(jj-1) + wgt * hphi(jj)
          ijtarg = ijtarg + 1
          ztarg = 1.0_wp + rhalf*REAL(ijtarg-1, KIND=wp)
       END IF

       jj = jj + 1
    END DO

    IF ((ijtarg /= ijntv) .OR. (NINT(ztarg) /= njc)) THEN ! tests are equivalent
       WRITE(*,*) "Error in msh_trp, next targ point not right"
       WRITE(*,*) "Array index: ",ijtarg, " should be", ijntv
       WRITE(*,*) "j value: ",ztarg, " should be", njc
       STOP
    END IF

    zphitrop(ijntv) = zphic

    ze1trop(:) = ra * rad * COS(zphitrop(:)*rad) * ppres !FC
    ze2trop(:) = ra * rad * dphidj(zphitrop(:), zc)

    ! Create new global grids of the correct size and paste in
    ! the correct sections of the input grid and the tropical grid.
    ! Output to netcdf and deallocate.
    !
    ! Note: we do not replace the row corresponding to phic, we just
    ! check that e2 agrees well there
    ! 
    ! Calculating new sizes:
    ! 1. in each hemisphere, replacing ijgtrop-1
    !    rows with njc-2 rows (not including equator or phic lines)
    ! 2. capr[icorn] and canc[er] denote the southern and northern limits
    !    of the tropical grid, respectively (lat=+/-phic).
    ! 3. Prefix iig, ijg or jp refers to unstreched grid, ij refers to
    !    tropical grid, and ijs refers to the full stretched grid

    ijgtrop = ijgcanc - jpeqt + 1  ! from equator to phic, inclusive
    ijgcapr = jpeqt - ijgtrop + 1
    ijgnth = SIZE(glam,2)
    iig = SIZE(glam,1)

    iis = iig
    ijsnth = ijgnth - (ijgtrop-2)*2 + (njc-2)*2
    ijscapr = ijgcapr
    ijseqt  = ijscapr + njc - 1
    ijscanc = ijseqt  + njc - 1

    ! Allocate temporary workspace to save existing global arrays
    ! Then reallocate global arrays with new larger size.
    ALLOCATE ( zwlam(jpim,jpjm,4), zwphi(jpim,jpjm,4), &
                zwe1(jpim,jpjm,4),  zwe2(jpim,jpjm,4)  )

    zwlam = glam
    zwphi = gphi
    zwe1  = e1
    zwe2  = e2

    DEALLOCATE( glam, gphi, e1, e2 )

    ALLOCATE (glam(iis,ijsnth,4), gphi(iis,ijsnth,4), &
                e1(iis,ijsnth,4),   e2(iis,ijsnth,4)  )

    ! south
    glam(:, 1:ijscapr, :) = zwlam(:, 1:ijgcapr, :)
    gphi(:, 1:ijscapr, :) = zwphi(:, 1:ijgcapr, :)
    e1(:, 1:ijscapr, :) =   zwe1(:, 1:ijgcapr, :)
    e2(:, 1:ijscapr, :) =   zwe2(:, 1:ijgcapr, :)

    !north
    glam(:, ijscanc:ijsnth, :) = zwlam(:, ijgcanc:ijgnth, :)
    gphi(:, ijscanc:ijsnth, :) = zwphi(:, ijgcanc:ijgnth, :)
    e1(:, ijscanc:ijsnth, :) =   zwe1 (:, ijgcanc:ijgnth, :)
    e2(:, ijscanc:ijsnth, :) =   zwe2 (:, ijgcanc:ijgnth, :)

    ! tropics:
    ! lam just copied into the extra rows
    ! phi, e1, e2 copied from tropical grid, using symmetry about equator

    glam(:, ijscapr+1:ijscanc-1, :) = &
         SPREAD( glam(:, ijscapr, :), 2, 2*njc-1 )

    CALL filltrop ( gphi, zphitrop, ijscapr, ijseqt, ijscanc, 2 )
    CALL filltrop ( e1,   ze1trop,  ijscapr, ijseqt, ijscanc, 3 )
    CALL filltrop ( e2,   ze2trop,  ijscapr, ijseqt, ijscanc, 4 )



!!!! To do:  finally check that e2(:, ijscanc, 1) ~= ze2trop(last)

    DEALLOCATE( hphi, hj, zphitrop, ze2trop, STAT=allostat )
    IF (allostat /= 0) STOP "Error deallocating batch 1"
    DEALLOCATE( zwlam, zwphi, zwe1, zwe2, STAT=allostat )
    IF (allostat /= 0) STOP "Error deallocating batch 2"



  END SUBROUTINE msh_trp



  SUBROUTINE  filltrop ( pgrid, ptrop, kjscapr, kjseqt, kjscanc, ktype )
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT) :: pgrid
    REAL(wp), DIMENSION(:), INTENT(IN) :: ptrop
    INTEGER, INTENT(IN) :: kjscapr, kjseqt, kjscanc, ktype

    REAL(wp) :: zsym
    INTEGER :: ijntv, iis

    ijntv = SIZE(ptrop)
    iis = SIZE(pgrid,1)

    ! The procedure for filling in phi, e1, e2 is as follows:
    ! 1. get north half of 1st column of T-grid from ptrop, incl. equator
    ! 2. get south half of 1st column of T-grid from ptrop, using symmetry
    ! 3. copy 1st column to all columns of T-grid with SPREAD
    ! 4. copy T-grid to U-grid
    ! 5. get north half of 1st column of V-grid from ptrop, no equator
    ! 6. get south half of 1st column of V-grid from ptrop, using symmetry
    ! 7. copy 1st column to all columns of V-grid with SPREAD
    ! 8. copy V-grid to F-grid

    IF (ktype == 2) THEN 
       zsym=-1.0 
    ELSE 
       zsym=1.0
    END IF

    pgrid(1,    kjseqt:kjscanc-1, 1) =  ptrop(1:ijntv-2: 2) 
    pgrid(1, kjscapr+1:kjseqt -1, 1) =  ptrop(ijntv-2:3:-2) * zsym
    pgrid(:, kjscapr+1:kjscanc-1, 1) = &
         SPREAD( pgrid(1, kjscapr+1:kjscanc-1, 1), 1, iis )    
    pgrid(:, kjscapr+1:kjscanc-1, 2) = pgrid(:, kjscapr+1:kjscanc-1, 1)


    pgrid(1,    kjseqt:kjscanc-1, 3) =  ptrop(2:ijntv-1: 2)
    pgrid(1,   kjscapr:kjseqt -1, 3) =  ptrop(ijntv-1:2:-2) * zsym
    pgrid(:,   kjscapr:kjscanc-1, 3) = &
         SPREAD( pgrid(1, kjscapr:kjscanc-1, 3), 1, iis )
    pgrid(:, kjscapr:kjscanc-1, 4) = pgrid(:, kjscapr:kjscanc-1, 3)

  END SUBROUTINE filltrop





END MODULE mshtrp

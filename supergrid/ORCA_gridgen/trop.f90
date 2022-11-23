MODULE trop
  USE kinds
  USE param, ONLY: rad, ra
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
IMPLICIT NONE
! This module contains the routines specific to the tropical grid we want

REAL(wp), PARAMETER :: &  !: grid parameters:
     ppres=2.0_wp/3.0_wp, &     !: resolution of baseline grid at equator
     ppresmin=1.0_wp/4.0_wp, & !: smallest grid resolution
     pphi_join=20.0   !: requested northern limit of tropics grid (degrees)


REAL(wp), PARAMETER ::           & !: numerical parameters:
     rzero=1.0e-14_wp,           & !: very small number
     rhalf=0.5_wp

INTEGER, PARAMETER :: jptint = 100  ! high res points per model grid box

!:  Arrays allocated in msh_trp, but global to allow use here:
REAL(wp), ALLOCATABLE, DIMENSION(:) :: & !:
     hphi, & !: high res phi grid
     hj      !: values of j on the phi grid

!: Scalars initialised in msh_trp, but global to allow use here:
INTEGER :: & !:
     njc, &  !: target value of j at which phi(j)=phic
     nhp     !: size of the high res phi grid

CONTAINS


  SUBROUTINE integ (pc, perr, ld_set_targ)
    REAL(wp), INTENT(IN) :: pc
    REAL(wp), INTENT(OUT) :: perr
    LOGICAL, INTENT(IN), OPTIONAL :: ld_set_targ

    REAL(wp), DIMENSION(size(hphi)) :: zdjdphi
    REAL(wp) :: zdj, zstep
    INTEGER :: jj
    LOGICAL :: ll_set_targ
    LOGICAL, PARAMETER :: verbose=.FALSE.

    IF (PRESENT(ld_set_targ)) THEN
       ll_set_targ = ld_set_targ
    ELSE
       ll_set_targ = .FALSE.
    END IF

    zstep = hphi(nhp) / REAL(nhp-1, KIND=wp)

    zdjdphi(:) = 1.0_wp / dphidj(hphi(:), pc)


    hj(1) = 1.0  ! index of equator point is 1
    DO jj = 2, nhp  ! Trapezium integration 0 -> phi
       zstep = hphi(jj) - hphi(jj-1)
       zdj = zstep  *  rhalf*( zdjdphi(jj-1) + zdjdphi(jj) )
       hj(jj) = hj(jj-1) + zdj
    END DO

    IF (ll_set_targ)  njc = NINT( hj(nhp) )

    perr = REAL(njc, KIND=wp) - hj(nhp)

    IF (verbose) WRITE(*,&
         FMT='(" Integration: c=",g15.8,"j(phic)=",g15.8," miss=",e16.8)') &
         pc, hj(nhp), perr


  END SUBROUTINE integ


  FUNCTION dphidj (pphi, pc)
    ! dphidj defines the desired grid sizes as a function of latitude

    !! WARNING, spacings actually drop below minres because of decrease
    !! in baseline spacings.  Drop is only in 4th dp - we could force
    !! up to minres if this is important.

    ! dummy arguments
    REAL(wp), DIMENSION(:), INTENT(IN) :: pphi   !: latitude
    REAL(wp), INTENT(IN) :: pc                  !: scale factor in function
    REAL(wp), DIMENSION(SIZE(pphi)) ::  dphidj   !: result of function

    ! parameters
    !INTEGER, PARAMETER :: jpexp=4 !: exponent used in "gaussian"
    !REAL(wp), PARAMETER :: ppa=12.0, ppb=8.0 !: quotients in function
    INTEGER, PARAMETER :: jpexp=8 !: exponent used in "gaussian"
    REAL(wp), PARAMETER :: ppa=12.0, ppb=16.0 !: quotients in function

    ! derived parameters
    REAL(wp), PARAMETER :: ppaa = ppa**jpexp, ppbb = ppb**jpexp

    ! local variables
    REAL(wp) :: zdpdj, zfac
    REAL(wp), DIMENSION(SIZE(pphi)) :: & !:
             zdpbase,                  & !: baseline dphidj values
             zcphi                       !:

    zcphi(:) = (pc * pphi(:))**jpexp
    zdpbase(:) = ppres * COS(pphi(:)*rad)  !  Mercator meridional grid spacing
    zfac = (ppres - ppresmin)

    dphidj = zdpbase - zfac * rhalf*( EXP(-zcphi/ppaa) + EXP(-zcphi/ppbb) )

  END FUNCTION dphidj



END MODULE trop

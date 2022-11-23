MODULE mshsth

  USE param
  USE common
  USE functions
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

CONTAINS

  SUBROUTINE msh_sth
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE mshsth
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     calcul de la grille globale et ses facteurs d echelle dans la
    !!	partie sud a partir des fonctions analytiques.
    !!
    !!   Method :
    !!   --------
    !!
    !!
    !!   Input :
    !!   -------
    !!      argument
    !!		ktype           : = 1  longitude
    !!				  = 2  latitude
    !!				  = 3  e1, facteur d echelle zonal
    !!				  = 4  e2, facteur d echelle meridien
    !!
    !!   Output :
    !!   -------
    !!      argument
    !!		ptab(jpim,jpjm,4): element de la grille defini dans la
    !!				  partie sud aux points T, U, V et F
    !!
    !!
    !!   Modifications :
    !!   ---------------
    !!       original  : 99 (G. Madec)
!!!---------------------------------------------------------------------
    !! parameters and commons
    !! ======================
    !!
    !!----------------------------------------------------------------------
    !! local declarations
    !! ==================
    !
    INTEGER :: ji, jj

!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1999)
!!!---------------------------------------------------------------------
    !
    !
    ! Calcul au sud de l equateur de la grille nord
    ! ---------------------------------------------
    ! (calcul a partir des fonctions analytiques)
    ! attention, decalage de -0.5 en ji
    !
    ! ... T-U points (jusqu a l equateur de la grille nord)
    !
    !   ... longitude (ktype=1)      
    DO jj = 1, jpeq
       DO ji = 1, jpim
          ! longitude
          glam(ji,jj,1) = fslam( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)     )
          glam(ji,jj,2) = fslam( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)     )
          glam(ji,jj,3) = fslam( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)+0.5 )
          glam(ji,jj,4) = fslam( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)+0.5 )
          ! latitude
          gphi(ji,jj,1) = fsphi( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)     )
          gphi(ji,jj,2) = fsphi( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)     )
          gphi(ji,jj,4) = fsphi( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)+0.5 )
          gphi(ji,jj,3) = fsphi( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)+0.5 )
          ! e1
          e1(ji,jj,1) = ra * rad * COS( rad*gphi(ji,jj,1) ) &
               * fsdila( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)     )
          e1(ji,jj,2) = ra * rad * COS( rad*gphi(ji,jj,2) ) &
               * fsdila( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)     )
          e1(ji,jj,3) = ra * rad * COS( rad*gphi(ji,jj,3) ) &
               * fsdila( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)+0.5 )
          e1(ji,jj,4) = ra * rad * COS( rad*gphi(ji,jj,4) ) &
               * fsdila( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)+0.5 )
          ! e2
          e2(ji,jj,1) = ra * rad &
               * fsdjph( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)     )
          e2(ji,jj,2) = ra * rad &
               * fsdjph( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)     )
          e2(ji,jj,3) = ra * rad &
               * fsdjph( REAL(ji,KIND=wp)-0.5 , REAL(jj,KIND=wp)+0.5 )
          e2(ji,jj,4) = ra * rad &
               * fsdjph( REAL(ji,KIND=wp)     , REAL(jj,KIND=wp)+0.5 )
       END DO
    END DO

  END SUBROUTINE msh_sth
END MODULE mshsth

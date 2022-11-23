MODULE mshnth
  USE kinds
  USE param
  USE common
  USE functions
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

CONTAINS

  SUBROUTINE msh_nth( ktype, pin, pout )
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE msh_nth
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     calcul de la grille globale et ses facteur d echelle dans la 
!!!	partie north (double pole) par extraction et symetrie de la
!!!	demi-grille nord et de ses facteur d echelle
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
    !!		pin(jpin,jpjn)  : element de la demi-grille nord a
    !!				  etendre sur la grille globale
    !!
    !!   Output :
    !!   -------
    !!      argument
    !!		ptab(jpim,jpjm,4): element de la grille defini sur
    !!				  tout le domaine, symetrique east-west
    !!				  et rempliement au nord au points T, 
    !!				  U, V et F.
    !!
    !!   External :
    !!   ----------
    !!              prihre          : 
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
    INTEGER ktype, jj, ji, ijit, ijjt, iis, ijist, ijisu, jpt
    !
    REAL(wp) :: pin (jpin,jpjn), pout(jpim,jpjm,4)

    !!----------------------------------------------------------------------
    !! statement functions
    !! ===================
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1999)
!!!---------------------------------------------------------------------
    !
    !
    ! 1. remplissage de la demi-grille nord pin dans pout
    ! ---------------------------------------------------
    !
    WRITE(0,*) 'mshnth: remplissage 1/2 grille nord ktype = ', ktype
    WRITE(0,*) '     i: ', 2/2 + 1, ' a ',(jpin-1-2 + 1 - 1)/2 + 1
    WRITE(0,*) '     j: ', (1 - 1)/2 + 1, ' a ',(jpjn-1-1 + 1)/2 + 1
    !
    ! ... points interieurs
    !     (ijit:    2 -> (jpi-1)/2
    !     (ijjt: jpeq ->  jpj
    !
    DO jj = 1, jpjn-1, 2
       DO ji = 2, jpin-1, 2
          !   ... indices en T
          ijit = ji/2 + 1
          ijjt = (jj-1)/2 + jpeq
          !   ... T-point
          pout(ijit,ijjt,1) = pin(ji  ,jj  )
          !   ... U-point
          pout(ijit,ijjt,2) = pin(ji+1,jj  )
          !   ... V-point
          pout(ijit,ijjt,3) = pin(ji  ,jj+1)
          !   ... F-point
          pout(ijit,ijjt,4) = pin(ji+1,jj+1)
       END DO
    END DO
    !
    ! ... ligne ji = 1
    IF ( ktype .EQ. 1 ) THEN
       !   ... longitude (ktype=1) changement de signe
       DO jj = 1, jpjn-1, 2
          ijjt = (jj-1)/2 + jpeq
          pout(1,ijjt,1) = -pout(2,ijjt,1)
          pout(1,ijjt,2) =  pin (1,jj)
          pout(1,ijjt,3) = -pout(2,ijjt,3)
          pout(1,ijjt,4) =  pin (1,jj+1)
       END DO
    ELSE
       !   ... latitude et e1, e2  pas de changement de signe)
       DO jj = 1, jpjn-1, 2
          ijjt = (jj-1)/2 + jpeq
          pout(1,ijjt,1) =  pout(2,ijjt,1)
          pout(1,ijjt,2) =  pin (1,jj)
          pout(1,ijjt,3) =  pout(2,ijjt,3)
          pout(1,ijjt,4) =  pin (1,jj+1)
       END DO
    ENDIF
    !
    ! 2. grille nord
    ! --------------
    ! symetrie grille nord par rapport a la ligne de U-F (indice iis)
    ! (cette symetrie assure la periodicite east-west)
    !
    WRITE(0,*) 'mshnth: extention a la grille nord ktype = ', ktype
    WRITE(0,*) '        symetrie par rapport a ji= ', (jpim-1)/2+1
    !
    IF ( ktype .EQ. 1 ) THEN
       !
       ! ... longitude (ktype=1)
       !      print*, "iis, jj, ji, ijist, pout(iis), pout(ijist):"
       iis = (jpim-1)/2+1
       DO jj = jpeq, jpjm-1
          DO ji = iis+1, jpim-1
             !   ... indice de symetrie
             ijist = 2*iis-ji+1
             ijisu = 2*iis-ji
             !   ... T-V points
             pout(ji,jj,1) = 2.*pout(iis,jj,2) - pout(ijist,jj,1)
             pout(ji,jj,3) = 2.*pout(iis,jj,4) - pout(ijist,jj,3)
             !   ... U-F points
             pout(ji,jj,2) = 2.*pout(iis,jj,2) - pout(ijisu,jj,2)
             pout(ji,jj,4) = 2.*pout(iis,jj,4) - pout(ijisu,jj,4)
          END DO
       END DO
       !

    ELSE
       !
       ! ... latitude et e1, e2
       iis = (jpim-1)/2+1
       DO jj = jpeq, jpjm-1
          DO ji = iis+1, jpim-1
             !   ... indice de symetrie
             ijist = 2*iis-ji+1
             ijisu = 2*iis-ji
             !   ... T-V points
             pout(ji,jj,1) = pout(ijist,jj,1)
             pout(ji,jj,3) = pout(ijist,jj,3)
             !   ... U-F points
             pout(ji,jj,2) = pout(ijisu,jj,2)
             pout(ji,jj,4) = pout(ijisu,jj,4)
          END DO
       END DO
       !
    ENDIF
    !
    ! ... Wrap rows for the north fold 
    !     (would be cleaner to do this before extracting the T,U,V,F-grids?)
    !
    ! ... repliement au nord ligne definir la valeur en jj=jpjm=jpj+1
    !     equivalent a une symetrie autour de la ligne F jj=jpj
    iis = (jpim-1)/2+1
    IF (ktype == 1) THEN
       !   ... longitude (needs mirrored about line of symmetry iis) 
       DO ji = 1, jpim
          !   ... indice de symetrie
          ijist = 2*iis-ji+1
          ijisu = 2*iis-ji
          !   ... T-U points
          pout(ji,jpjm,1) = pout(ijist,jpjm-1,1)
          pout(ji,jpjm,2) = pout(ijisu,jpjm-1,2)
          !   ... V-F points
          pout(ji,jpjm,3) = pout(ijist,jpjm-2,3)
          pout(ji,jpjm,4) = pout(ijisu,jpjm-2,4)
       END DO
    ELSE
       !   ... lat, e1, e2 (simply copy the appropriate rows)
       DO ji = 1, jpim
          !   ... T-U points
          pout(ji,jpjm,1) = pout(ji,jpjm-1,1)
          pout(ji,jpjm,2) = pout(ji,jpjm-1,2)
          !   ... V-F points
          pout(ji,jpjm,3) = pout(ji,jpjm-2,3)
          pout(ji,jpjm,4) = pout(ji,jpjm-2,4)
       END DO
    END IF
    !
    ! ... symetrie east-west pour obtenir des long. croissantes
    WRITE(0,*) 'mshnth: longitude croissante x -1  '
    IF ( ktype .EQ. 1 ) THEN
       DO jpt = 1, 4
          DO jj = jpeq, jpjm
             DO ji = 1, jpim
                pout(ji,jj,jpt) = - pout(ji,jj,jpt)
             END DO
          END DO
       END DO
       !
    ENDIF
    !
    !
    ! ... cyclic east-west sur toute la grille (nord & sud)
    WRITE(0,*) 'mshnth: cyclic east-west sur mesh nord & sud'
    IF ( ktype .EQ. 1 ) THEN
       DO jpt = 1, 4
          DO jj = 1, jpjm
             pout(jpim-1,jj,jpt) = pout( 1 ,jj,jpt) + 360.
             pout(jpim  ,jj,jpt) = pout( 2 ,jj,jpt) + 360.
          END DO
       END DO
    ELSE
       DO jpt = 1, 4
          DO jj = 1, jpjm
             pout(jpim-1,jj,jpt) = pout( 1 ,jj,jpt)
             pout(jpim  ,jj,jpt) = pout( 2 ,jj,jpt)
          END DO
       END DO
    ENDIF
    !
    ! ... and 'rotate' so that central longitude is LRAMPO
    !     and put longs in range [-180, 180] (NB: because MOD leaves -ve
    !     values, add 720 before applying)
    IF ( ktype == 1 ) THEN
       pout(:,:,:) = pout(:,:,:) + RLAMPO 
       pout(:,:,:) = MOD(pout(:,:,:)+720.0, 360.0) - 180.0
    ENDIF
    !

  END SUBROUTINE msh_nth
END MODULE mshnth

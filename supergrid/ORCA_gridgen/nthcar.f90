MODULE nthcar
  USE kinds
  USE param
  USE common
  USE functions
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

CONTAINS

  SUBROUTINE nth_car (kji)
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE nth_car
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!        Construction d une ligne meridienne et de ses facteurs 
!!!	d echelle partant du parallele geographique le plus nord
!!!	(i.e. l"equateur" de la grille).
!!!
    !!   Method :
    !!   --------
    !!	integration de l eq. diff. fait sur jph = 20 *(jpjf-1)+1 points
    !!	point lat/long sur la grille fine jpjf point et sur les deux 
    !!	meridiens avant et apres pour evaluation des e1, e2nth (grille fine
    !!	au moins 2 fois plus fine que la grille a calculer):
    !!	jpjf=jpfacj*(jpnord-1)+1+1
    !!      stockage des lat/lon et e1/e2nth sur la grille nord (jpin, jpjn)
    !!	jpin= jpi  jpjn= 2*jpnord
    !!      
    !!	Approximation du premier ou du second ordre
    !!      On se place dans le plan stereographique polaire, en coordonnees
    !!      cartesiennes avec le cercle equateur de la grille de rayon REQ
    !!
    !!	ETAPE I : on integre depuis le cercle equatorial de rayon REQ
    !!		  l equation differentielle:
    !!			Y = FSDER(J(X,Y),X,Y)
    !!		  integration avec un pas d espace jph fois plus petit
    !!		  que la demi grille nord fine, avec un schema leap-frog
    !!		  (2eme ordre).
    !!	ETAPE II: on interpole la courbe resultante sur la demi grille
    !!		  nord fine avec une approximation du premier ou second
    !!		  ordre. (attention interpolation a partir de pts non
    !!		  equirepartis)
    !!
    !!
    !!   Input :
    !!   -------
    !!      argument		:
    !!		kji		: indice ji de la ligne meridienne a
    !!				  construire.
    !!      common                  : no
    !!
    !!   Output :
    !!   --------
    !!      argument                : NO
    !!      common
    !!              /nthdgf/ xf, yf  : abscisses et ordonnees des pts de la
    !!				  demi grille nord fine (en coord. car-
    !!				  TESIENNES DANS LE PLAN STEREOGRAPHIQUE
    !!				  POLAIRE AVEC LE CERCLE EQUATEUR = LE
    !!                		  CERCLE UNITE)
    !!
    !!   Modifications :
    !!   ---------------
    !!       original  : 92 (G. Madec)
    !!       update    : 99 (G. Madec) new computation (2 poles)
!!!---------------------------------------------------------------------
    !! parameters and commons
    !! ======================
    !!----------------------------------------------------------------------
    !! dummy variables
    !! ===============
    INTEGER, INTENT(IN) :: kji     ! index of meridien to compute

    !! local declarations
    !! ==================
    INTEGER itab(jpjf)
    INTEGER jcmin, jcmax
    !
    REAL(wp) :: zwj(jph), zwx(jph), zwy(jph), zwp(jph)
    REAL(wp) :: zlamf(jpjf,-2:2), zphif(jpjf,-2:2)
    REAL(wp) :: zster(-2:2)

    INTEGER :: jc, jh, jn, jj, ijh, ij, ii

    REAL(wp) :: zi, zlam, zdx, zj1, zj2, zj3, zfs, zj, zpente, zcen, za, zb
    REAL(wp) :: ztheta, zyp1, zyp2, zypmin, zypmax, zy, zyf, zxf, zcosp, zdldi
    REAL(wp) :: zdpdi, zdldj, zdpdj, gamma

!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1999)
!!!---------------------------------------------------------------------
    !
    !
    ndebug = 0
    !
    ! =======================================
    ! I. Contruction d une courbe meridienne
    !    et de 4 autres de part et d autre
    ! =======================================
    !
    ! Boucle sur les 5 courbes
    ! ========================
    !
    IF (l_e1_calc) THEN  ! Need 5 curves to get zonal scale factors
       jcmin = -2  
       jcmax = 2
    ELSE  ! Only need 1 curve if just doing lats, lons + zonal scale factors
       jcmin = 0  
       jcmax = 0
    ENDIF

    DO jc = jcmin, jcmax
       !
       !
       !    1. ABSCISSES DES POINTS D UNE COURBE
       !    ====================================
       !    
       ! ... Indice j du cercle equateur egal a 1
       zwj(1) = 1.
       ! ... indice sur la ligne pole: jpnord+0.5 par construction
       !     (i.e. la ligne nord est une ligne passant par des points V et F)
       zwj(jph) = FLOAT(jpnord) + 0.5
       !
       ! ... Coordonnees initiales: equateur de la grille de rayon req
       !     position en indice kji varie avec un pas de 2*jpfac autour
       !     de kji
       zi= 1. + FLOAT(kji-1) / 2. + FLOAT(jc) / FLOAT(2*jpfac)
       zlam = rad * fsx0p( zi )
       zwx(1) = req * COS( zlam )
       zwy(1) = req * SIN( zlam )
       !
       !     en radian, les coordonnees du cercle eq de la grille sont:
       !     (longitude de 0 a -180 pas de rdx/2 degre)
       zlamf(1,jc) = rad * ( -rdx*(zi-1)/2. )
       zphif(1,jc) = rad * rphieq
       !
       IF ( ndebug .EQ. 2 ) THEN
          write(0,*) 'kji= ', kji,' jc= ', jc, &
               ' lam/phi= ', zlamf(1,jc), zphif(1,jc), &
               ' lam/phi= ', zlamf(1,jc)/rad, zphif(1,jc)/rad
       ENDIF
       !
       ! ... Pas d espace: (toujours le meme nb de pts par courbe)
       zdx  = zwx(1) / FLOAT(jph-1)
       !
       ! ... abscisse des points
       DO jh = 2, jph
          zwx(jh) = zwx(1) - zdx * FLOAT(jh-1)
       END DO
       !
       !
       !    2. ORDONNEE ET INDICE DES POINTS D UNE COURBE
       !    =============================================
       ! (integration de l equation differentielle en leap-frog)
       !
       ! 1) Demarage Matsuno
       ! -------------------
       !
       zwy(2) = zwy(1) - zdx * fsder( zwj(1), zwx(1), zwy(1) )
       !
       ! 2) Integration Leapfrog
       ! -----------------------
       !
       ! ... Recherche de zwj(jh) tel que fsfun(zwj,zwx,zwy)=0.
       !     (zero atteint a la precision de la machine pres)
       !
       DO jh = 2, jph-1
          !
          !   ... Indice varie au maximum de 1 a jpjnord+0.5
          zj1 = 1.
          zj2 = FLOAT(jpnord) + 0.5
          !   ... methode de dicotomie
          DO jn = 1, ncompt
             zj3     = ( zj1 + zj2 ) * 0.5
             zfs     = fsfun( zj3, zwx(jh), zwy(jh) )
             !     ... Si ABS(zfs)=< a rzero le zero est atteint on ne change
             !         plus d indice zj1=zj2=zj3
             !ccc        zj2 = CVMGP( zj3, zj2, (rzero-ABS(zfs)) )
             !ccc        zj1 = CVMGP( zj3, zj1, (rzero-ABS(zfs)) )
             IF( ABS(zfs) .LT. rzero ) zj2 = zj3
             IF( ABS(zfs) .LT. rzero ) zj1 = zj3
             !     ... sinon, on change zj1 ou zj2 suivant le signe de zfs
             !ccc        zj1 = CVMGM( zj3, zj1, zfs )
             !ccc        zj2 = CVMGP( zj3, zj2, zfs )
             IF( zfs .LT. 0. ) zj1 = zj3
             IF( zfs .GT. 0. ) zj2 = zj3
             !
          END DO
          !
          !   ... stockage de l indice du zero
          zwj(jh) = ( zj1 + zj2 ) * 0.5
          !
          !   ... Integration spatiale en jh+1
          zwy(jh+1) = zwy(jh-1) &
                         - 2.* zdx * fsder( zwj(jh), zwx(jh), zwy(jh) )
          !
          !   ... filtre d Asselin sur zwy(jh)
          gamma=1.e-1
          zwy(jh) = zwy(jh) + gamma*( zwy(jh-1) -2*zwy(jh) + zwy(jh+1) )
          !       
       END DO
       !
       IF ( ndebug .EQ. 3 ) THEN
          write(70,*)
          write(70,*) 'kji= ', kji,' jc= ', jc
          DO jh= jph, 1, -1
             write(70,*) jh, zwj(jh),zwx(jh),zwy(jh)
          END DO
       ENDIF
       !
       !
       ! ==================================
       ! II. CONSTRUCTION DE LA GRILLE FINE
       ! ==================================
       !     
       ! Grille fine en jpjf pour le calcul des facteurs d echelle
       !
       !
       ! 1) Recherche de l indice du point juste avant la "latitude" JJ
       ! --------------------------------------------------------------
       ! (le point jh=1+(jj-1)/jpfacj correspondant au cercle d indice jj)
       !     
       DO jj = 1, jpjf-1
          itab(jj) = 0
          zj =  1. + FLOAT(jj-1)/FLOAT(2*jpfac)
          DO jh = 1, jph-1
             IF ( zj.GE.zwj(jh) .AND. zj.LT.zwj(jh+1) ) itab(jj) = jh
          END DO
       END DO
       itab(jpjf) = jph
       !
       ! ... controle
       ijh = 0
       DO jj = 1, jpjf-1
          IF ( itab(jj) .EQ. 0 ) ijh = jj
       END DO
       IF ( ijh .NE. 0 ) THEN
          WRITE(0,*) 'nth_car: error ligne ji=', kji, ' jj=', ijh
       ENDIF
       !
       !
       !  2) calcul de la pente en tout point de la courbe
       !  ------------------------------------------------
!!!! pente = slope
       !
       ! ... pente sur la ligne pole: la pente est nulle
       zwp(jph) = 0.e0
       !
       ! ... pente en jh: donnee par fsder   
       DO jh = 1, jph-1
          zwp(jh) = fsder( zwj(jh), zwx(jh), zwy(jh) )
       END DO
       !
       !  3) boucle sur les points de grille
       !  ----------------------------------
       ! on ne prend pas en compte le point d arrivee: il est sur la ligne pole
       !
       DO jj = 1, jpjf-1
          !
          ij = itab(jj)
          zj = 1. + FLOAT(jj-1) / FLOAT(2*jpfac)
          !
          ! ... calcul de la pente de la courbe au point du maillage (i.e. a
          !     l intersection du cercle indice ZJ et de la courbe kji ) par
          !     interpolation lineaire
          zpente = zwp(ij) &
                      + ( zj - zwj(ij) ) * ( zwp(ij+1) - zwp(ij) ) &
                                         / ( zwj(ij+1) - zwj(ij) )

          !
          ! ... Abscisses et ordonnees du point: il est sur l ellipse d indice
          !     zj et peut etre reperi par la valeur de sa normale zpente*za/zb
          !   ... centre de l ellipse:
          zcen = fsy0(zj)
          !   ... grand et petit axes:
          za = fsa(zj)
          zb = fsb(zj)
          !   ... position des points:
          ztheta = ATAN( zpente*za/zb )
          zlamf(jj,jc) = zb * COS( ztheta ) 
          zphif(jj,jc) = za * SIN( ztheta ) + zcen 
          !
          !         zlamf, zphif here are x, y in polar st. plane
       END DO
       !
       ! ... point d arrivee
       !     le dernier point est sur la ligne pole:
       !     zwx=0, zwy(jph), zwp=0, et zwj=jpnord+0.5
       !     On controle que se point trouve sur le segment final, si il en
       !     sort, on le place au pole le plus proche de sq position
       !     
       !   ... abscisse: sur la ligne pole, XF(jph)=0
       zlamf(jpjf,jc) = 0.
       !
       !   ... ordonnees des 2 poles (classes en ordre croissant)
       zyp1 = fsy0( rjpnth ) - fsa( rjpnth )
       zyp2 = fsy0( rjpnth ) + fsa( rjpnth )
       zypmin = min( zyp1, zyp2 )
       zypmax = max( zyp1, zyp2 )
       !
       !   ... verification: zypmin <= yf(jph) <= zypmax
       zy = min ( zwy(jph), zypmax )
       zphif(jpjf,jc) = max ( zy      , zypmin )
       !
       !
       !
       !  4) passage en coordonnees geographiques (pour la grille fine)
       !  ---------------------------------------
       !  attention, sauf pour la ligne eq. de la grille deja initialisee en
       !             coord geographique (en radian)
       ! ... meridienne maitresse (jc=0)
       !     (tous les points jpjf sont calcules)
       !
       !     IF ( zlamf(1,jc) .EQ. 0. ) THEN

       IF ( ( kji.EQ.1 .OR. kji.EQ.jpin ) .AND. ( jc.EQ.0 ) ) THEN
          WRITE(0,*) ' meridienne maitresse kji= ', kji
          !
          ! ... meridien demarrant sur l axe des y 
          ! Warning: if x happens to be -ve, then the sign of y (=zlamf)
          ! is reversed and this IF test will not work
          DO jj = 1, jpjf
             IF ( zphif(jj,jc).GT.0. ) THEN
                zlamf(jj,jc) =    0.
             ELSE
                zlamf(jj,jc) =-rpi
             ENDIF
          END DO
          DO jj = 1, jpjf
             zyf = zphif(jj,jc)
             zphif(jj,jc) = rpi/2. - 2.*ATAN( ABS( zyf ) )
          END DO
          !
       ELSE
          !
          ! ... segment interpole (jj=jpjf)
          IF ( zphif(jpjf,jc) .GE. 0. ) THEN
             zlamf(jpjf,jc) =  0.
          ELSE
             zlamf(jpjf,jc) =-rpi
          ENDIF
          zyf = zphif(jpjf,jc)
          zphif(jpjf,jc) = rpi/2. - 2.*ATAN( ABS( zyf ) )
          !
          ! ... hors ligne pole (kji=1 ou jpi et jj=jpjf
          DO jj = 1, jpjf-1
             zxf = zlamf(jj,jc)
             zyf = zphif(jj,jc)
             zlamf(jj,jc) =-rpi/2. + ATAN( zyf / zxf )
             zphif(jj,jc) = rpi/2. - 2.*ATAN( SQRT( zxf*zxf+zyf*zyf ) )
          END DO
          !
       ENDIF
       !

       IF ( ndebug .EQ. 3 ) THEN
          write(70,*)
          write(70,*) 'kji= ', kji,' jc= ', jc
          DO jj= jpjf, 1, -1
             write(70,*) jj, zlamf(jj,jc)/rad,zphif(jj,jc)/rad
          END DO
       ENDIF
       !
       !
       !
       ! Fin de boucle sur les 5 lignes meridiennes
       ! ==========================================
       !
    END DO
    !
    !
    ! ... ligne autour des meridiennes maitresses (kji=1 et kji=jpj)
    IF( kji .EQ. 1 ) THEN
       DO jj = 1, jpjf
          zlamf(jj,-2) = -zlamf(jj, 2)
          zlamf(jj,-1) = -zlamf(jj, 1)
          zphif(jj,-2) =  zphif(jj, 2)
          zphif(jj,-1) =  zphif(jj, 1)
       END DO
    ENDIF
    IF( kji .EQ. jpin ) THEN
       DO jj = 1, jpjf
          zlamf(jj, 2) = 2.*zlamf(jj, 0) - zlamf(jj,-2)
          zlamf(jj, 1) = 2.*zlamf(jj, 0) - zlamf(jj,-1)
          zphif(jj, 2) =  zphif(jj,-2)
          zphif(jj, 1) =  zphif(jj,-1)
       END DO
    ENDIF
    !
    ! affichage des zlamf, zphif calcules
    IF ( ndebug .GE. 2 ) THEN
       write(70,*)
       write(70,*) 'kji= ', kji,' zlamf '
       DO jj= jpjf, 1, -1
          write(70,*)     jj, zlamf(jj,-2)/rad, zlamf(jj,-1)/rad, &
                   zlamf(jj, 0)/rad, zlamf(jj, 1)/rad, zlamf(jj, 2)/rad
       END DO
       write(70,*)
       write(70,*) 'kji= ', kji,' zphif '
       DO jj= jpjf, 1, -1
          write(70,*)     jj, zphif(jj,-2)/rad, zphif(jj,-1)/rad, &
                   zphif(jj, 0)/rad, zphif(jj, 1)/rad, zphif(jj, 2)/rad
       END DO
    ENDIF
    !
    !
    ! ==============================================================
    ! III. Sauvegarde de la grille nord et de ses facteurs d echelle
    ! ==============================================================
    !
    ! 1) Sauvegarde des lat/lon en degre
    ! ----------------------------------
    !
    print*,'min/max longitudes: ',zlamf(1,0)/rad,zlamf(jpjf-1,0)/rad
    print*,'min/max latitude:   ',zphif(1,0)/rad,zphif(jpjf,0)/rad

    DO jj = 1, jpjn
       ij = jpfac*(jj-1) + 1
       glamnth(kji,jj) =  zlamf(ij,0) / rad
       gphinth(kji,jj) =  zphif(ij,0) / rad
    END DO
    !
    IF (l_e1_calc) THEN
       ! 2) Calul du facteur d echelle zonal
       ! -----------------------------------
       ! (attention: zlamf, zphif sont en radian pas de mult. par rad)
       !
       ! ... en jj=1 la valeur de la derivee est connue analytiquement
       !!      e1(kji,1) = ra * rad * COS( rad*gphi(kji,1) )
       !!     $              * fsdila( FLOAT(kji), FLOAT(JPEQ) )
       !     on la calcule qaund meme pour avoir une idee de l erreur du
       !     calcul des derivee
       !
       ! ... Calcul de derivee par rapport a I, a J constant
       !     (derivee suivant I d ordre 4, pas d integration 1/(2jpfac) )
       DO jj = 1, jpjn-1
          ij = jpfac*(jj-1) + 1
          zcosp = COS( zphif(ij, 0) )
          zdldi =    zlamf(ij,-2) - 8.*zlamf(ij,-1) &
                   + 8.*zlamf(ij,+1) -    zlamf(ij,+2)
          zdpdi =    zphif(ij,-2) - 8.*zphif(ij,-1) &
                   + 8.*zphif(ij,+1) -    zphif(ij,+2)
          e1nth(kji,jj) = ra * FLOAT(2*jpfac)/12. &
                            * SQRT( zcosp*zcosp*zdldi*zdldi + zdpdi*zdpdi ) 
       END DO
       !
       ! ... en jpjn, ligne nord, zlamf=cste mais saut de 180 degres au passage
       !     du pole nord geographique. On impose zdldi = 0
       ! ... However, for the point at the pole, lambda is discontinuous so
       !     does not exist, and zdpdi is zero(ish).  We need to transform to 
       !     coordinates which are continuous across the pole.  Here we use
       !     stereographic coords.  The coord normal to the fold is constant.

       ij = jpjf
       DO ii = -2, 2
          IF (ABS(zlamf(ij,ii)) < rpi/2.) THEN
             zster(ii) = zphif(ij,ii) - rpi/2.0
          ELSE
             zster(ii) = rpi/2.0 - zphif(ij,ii)
          END IF
       END DO

       zdpdi =    zster(-2) - 8.*zster(-1) &
            + 8.*zster(+1) -    zster(+2)
       e1nth(kji,jpjn) = ra * FLOAT(2*jpfac)/12. * ABS( zdpdi ) 

    ENDIF     ! l_e1_calc

    !
    ! 3) Calul du facteur d echelle meridien
    ! --------------------------------------
    ! (attention: zlamf, zphif sont en radian pas de mult. par rad)
    !
    ! ... en jj=1 la valeur de la derivee est connue analytiquement
    e2nth(kji,1) = ra * rad &
         * fsdjph( REAL(kji,KIND=wp), REAL(JPEQ,KIND=wp) )
    ! 
    ! ... Calcul de derivee par rapport a J, a I constant
    !     (derivee suivant J d ordre 4, pas d integration 1/(2*jpfac) )
    DO jj = 2, jpjn-1
       ij = jpfac*(jj-1) + 1
       zcosp = COS( zphif(ij, 0) )
       zdldj =    zlamf(ij-2,0) - 8.*zlamf(ij-1,0) &
            + 8.*zlamf(ij+1,0) -    zlamf(ij+2,0)
       zdpdj =    zphif(ij-2,0) - 8.*zphif(ij-1,0) &
            + 8.*zphif(ij+1,0) -    zphif(ij+2,0)
       e2nth(kji,jj) = ra * FLOAT(2*jpfac)/12. &
            * SQRT( zcosp*zcosp*zdldj*zdldj + zdpdj*zdpdj )
    END DO
    !
    ! ... en jpjn, ligne nord  symetrie
    !       ==> meme latitude ==> zdpdj=0
    !       ==> en longitude:  zlamf(ij-2,0)=2*zlamf(ij,0)-zlamf(ij+2,0)
    !                          zlamf(ij-1,0)=2*zlamf(ij,0)-zlamf(ij+1,0)
    !
    ! ... this is true provided the point is not at the pole.  If it is
    !     we need coordinates which are continuous at the pole.  We use
    !     stereographic coordinates.  The coord along the fold is constant
    !     and the meridien is a geographic meridien so we can derive the
    !     other coordinate from latitude alone.
    !     ( I'm being a bit lazy here - to reduce errors, we should
    !       transform to full stereographic (or rotated-pole) for all points
    !       above a certain latitude. )
    ij = jpjf
    zcosp = COS( zphif(ij, 0) )
    IF ( ABS(zphif(ij,0) - rpi/2.0) > 1.0e-6 ) THEN
       ! not at the pole
       zdldj = 2.*zlamf(ij-2,0) - 16.*zlamf(ij-1,0) + 14.*zlamf(ij,0)
       e2nth(kji,jpjn) = ra * FLOAT(2*jpfac)/12. * ABS( zcosp*zdldj ) 
    ELSE
       ! at pole
       zster(-2) = zphif(ij-2,0) - rpi/2.0
       zster(-1) = zphif(ij-1,0) - rpi/2.0
       zster(0)  = zphif(ij  ,0) - rpi/2.0
       zster(1)  = rpi/2.0 - zphif(ij-1,0)
       zster(2)  = rpi/2.0 - zphif(ij-2,0)
       zdpdj =    zster(-2) - 8.*zster(-1) &
            + 8.*zster(+1) -    zster(+2)
       e2nth(kji,jpjn) = ra * FLOAT(2*jpfac)/12. * ABS( zdpdj ) 
    END IF

    !
    !     
    ! There are problems if scale factors are zero, which happened for ORCA1
    e1nth(kji,jpjn) = MAX( e1nth(kji,jj), 1.0 )
    e2nth(kji,jpjn) = MAX( e2nth(kji,jj), 1.0 )
    !
    IF ( ndebug .EQ. 2 ) THEN
       write(70,*)
       write(70,*) 'kji= ', kji,' e1nth  e2nth'
       DO jj= jpjn, 1, -1
          write(70,*) jj, e1nth(kji,jj),e2nth(kji,jj)
       END DO
    ENDIF
       !
  END SUBROUTINE nth_car

END MODULE nthcar

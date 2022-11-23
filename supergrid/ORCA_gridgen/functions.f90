MODULE functions
  USE kinds
  USE param
  USE common
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

  !! COMMON comcff: coefficients des functions
  !! ----------------------------------------------------
  !!      rjpeqt  : indice de la position de l equateur terrestre et sa
  !!                valeur en reel
  !!      rrotlam : Rotation de la grille autour de l axe nord-sud 
  !!      rjpnthm1: jpnth -1 
  !!      nrota	  : Angle de rotation (partie entiere de RLAMPO)
  !!	  rjp..	  : Conversion reelle des parametres
  !!      ra.,rb.,: Coefficients des fonctions de base
  !!      ry0
  !!      rphif.  : Valeurs finales des fonctions de base
  !!
  REAL(wp) rjpeqt, nrota, rjpnthm1, rrotlam, rjpeq, rjpnp, rjpj,&
       ramin, rajust, rphifa, ra1, ra2,   rphifb, rbmax,   &
       rb1, rb2, ry0
  
  REAL(wp) :: rphi1, rphi2, ry1, ry2, rphify0

CONTAINS

  FUNCTION atanh (px)
    REAL(wp) :: atanh
    REAL(wp), INTENT(IN) :: px

!!$C
!!$C  FUNCTION : POSITION GEOGRAPHIQUE
!!$C ******************  (PSEUDO LATITUDE ET LONGITUDE)
!!$C
!!$C ==========================
!!$C TRIGONOMETRIE HYPERBOLIQUE
!!$C ==========================
!!$C
!!$C Argument de tangente hyperbolique
!!$C

      atanh = 0.5 * LOG( (1.+px) / (1.-px) )

    END FUNCTION atanh
!!$C
!!$C ===========
!!$C CONVERSIONS
!!$C ===========
!!$C
!!$C I. Plan stereographique -> geographiques
!!$C ========================================
!!$C
!!$C Passage des coordonees du plan stereographique aux coordonnees
!!$C geographiques
!!$C
!!$C   Fonction
!!$C   --------
!!$C
    FUNCTION fsyphi(py)
      REAL(wp) :: fsyphi
      REAL(wp), INTENT(IN) :: py

      fsyphi = 90. - 2 * ATAN(py) / rad

    END FUNCTION fsyphi
!!$C
!!$C   Derivee premiere
!!$C   ----------------
!!$C     
    FUNCTION fsdyphi(py)
      REAL(wp) :: fsdyphi
      REAL(wp), INTENT(IN) :: py

      fsdyphi = -2 / ( rad * (1+py**2) )

    END FUNCTION fsdyphi
!!$C
!!$C II. Geographiques -> stereographique
!!$C ====================================
!!$C
!!$C Passage des coordonnees geographiques en degre au plan stereographique
!!$C polaire nord
!!$C
    FUNCTION fsphiy(pphi)
      REAL(wp) :: fsphiy
      REAL(wp), INTENT(IN) :: pphi

      fsphiy=TAN( rpi/4. - rad * pphi/2. )

    END FUNCTION fsphiy
!C
    FUNCTION fsdphiy(pphi)
      REAL(wp) :: fsdphiy
      REAL(wp), INTENT(IN) :: pphi

      FSDPHIY=-0.5*RAD*(1+(FSPHIY(PPHI))**2)

    END FUNCTION fsdphiy
!!$C
!!$C ===================================
!!$C FONCTIONS DE LA PARTIE SUD DU GLOBE
!!$C ===================================
!!$C
!!$C Pour ces fonctions , PJ varie de 1 (au pole sud ) 
!!$C a JPEQ ( a lequateur).
!!$C
!!$C
!!$C Fonction longitude (hemisphere sud)
!!$C ===================================
!!$C
!!$C Fonction
!!$C
    FUNCTION fsx0p(pi)
      REAL(wp) :: fsx0p
      REAL(wp), INTENT(IN) :: pi

      FSX0P = 90.-RDX*(PI-1.)

    END FUNCTION fsx0p
!!$C 
!!$C Derivee premiere
!!$C
    FUNCTION fsdx0p(pi)
      REAL(wp) :: fsdx0p
      REAL(wp), INTENT(IN) :: pi

      FSDX0P = -RDX

    END FUNCTION fsdx0p
!!$C
!!$C Fonction latitude (hemisphere sud)
!!$C ==================================
!!$C     
    FUNCTION theta(pj)
      REAL(wp) :: theta
      REAL(wp), INTENT(IN) :: pj

      THETA = RDX*RAD*(PJ-RJPEQT)

    END FUNCTION theta
!!$C
!!$C Fonction
!!$C
    FUNCTION fsy0p(pj)
      REAL(wp) :: fsy0p
      REAL(wp), INTENT(IN) :: pj

      FSY0P = asin(tanh(THETA(PJ)))/RAD

    END FUNCTION fsy0p
!!$C     
!!$C Derivee premiere
!!$C   
    FUNCTION fsdy0p(pj)
      REAL(wp) :: fsdy0p
      REAL(wp), INTENT(IN) :: pj

      FSDY0P = RDX/cosh(THETA(PJ))

    END FUNCTION fsdy0p
!!$C
!!$C Derivee seconde
!!$C
    FUNCTION fsd2y0p(pj)
      REAL(wp) :: fsd2y0p
      REAL(wp), INTENT(IN) :: pj

      FSD2Y0P = -(RDX**2)*RAD* &
          sinh(THETA(PJ))/(cosh(THETA(PJ)**2)) 

   END FUNCTION fsd2y0p
!!$C     
!!$C Derivee troisieme
!!$C  
    FUNCTION fsd3y0p(pj)
      REAL(wp) :: fsd3y0p
      REAL(wp), INTENT(IN) :: pj
   
      FSD3Y0P = (RDX**3)*(RAD**2)* &
          ((sinh(THETA(PJ)))**2-1.)/((cosh(THETA(PJ)))**3)

   END FUNCTION fsd3y0p
!!$C     
!!$C Derivee quatrieme     
!!$C  
    FUNCTION fsd4y0p(pj)
      REAL(wp) :: fsd4y0p
      REAL(wp), INTENT(IN) :: pj
 
      FSD4Y0P = (RDX**4)*(RAD**3)*sinh(THETA(PJ))* &
          (5.-(sinh(THETA(PJ)))**2)/ &
          (cosh(THETA(PJ)))**4

   END FUNCTION fsd4y0p
!!$C
!!$C =====================    
!!$C CHANGEMENTS D ORIGINE
!!$C =====================
!!$C     
!!$C On se ramene a JPEQ=0
!!$C  
    FUNCTION fsj(pj)
      REAL(wp) :: fsj
      REAL(wp), INTENT(IN) :: pj
  
      FSJ=PJ-RJPEQ+1

    END FUNCTION fsj
!!$C
!!$C ========================
!!$C FONCTIONS INTERMEDIAIRES
!!$C ========================
!!$C     
!!$C Fonctions associees a FSPHIA
!!$C    
    FUNCTION fsfa(pj)
      REAL(wp) :: fsfa
      REAL(wp), INTENT(IN) :: pj

      FSFA=(ATAN(RA2*(PJ-1.)/rjpnthm1)+(RA2*(PJ-1.)/rjpnthm1)**3/ &
          (3*((RA2*(PJ-1.)/rjpnthm1)**3/(RA1*RA2)+1.)))*rjpnthm1/RA2
!!$C      FSFA=(ATAN(RA2*(PJ-1.))+(RA2*(PJ-1.))**3/
!!$C     $     (3*((PJ-1.)**3/(RA1*RA2)+1.)))/RA2
!!$C

   END FUNCTION fsfa

    FUNCTION fsdfa(pj)
      REAL(wp) :: fsdfa
      REAL(wp), INTENT(IN) :: pj

      FSDFA=1./(1.+(RA2*(PJ-1.)/rjpnthm1)**2)+ &
          (RA2*(PJ-1.)/rjpnthm1)**2/(((PJ-1.)/rjpnthm1)**3/(RA1*RA2)+ &
          1.)**2

   END FUNCTION fsdfa
!!$C     
!!$C Fonctions associees a FSPHIB     
!!$C
    FUNCTION fsfb(pj)
      REAL(wp) :: fsfb
      REAL(wp), INTENT(IN) :: pj

      FSFB=1.+RB1*tanh(RB2*((PJ-1.)/rjpnthm1)**3)

    END FUNCTION fsfb
!!$C   
    FUNCTION fsdfb(pj)
      REAL(wp) :: fsdfb
      REAL(wp), INTENT(IN) :: pj
  
      FSDFB=RB1*RB2/rjpnthm1*((PJ-1.)/rjpnthm1)**2/ &
          (cosh(RB2*((PJ-1.)/rjpnthm1)**3))**2

   END FUNCTION fsdfb
!!$C
!!$C =================
!!$C FONCTIONS DE BASE
!!$C =================
!!$C
!!$C Pour ces fonctions, l indice PJ varie de 1 (a l equateur) a JPMAX+1
!!$C (au pole)     
!!$C     
!!$C GRAND AXE DES ELLIPSES
!!$C ======================    
!!$C     
!!$C   Fonction
!!$C   --------  
!!$C
    FUNCTION fsphia(pj)
      REAL(wp) :: fsphia
      REAL(wp), INTENT(IN) :: pj

      FSPHIA = FSY0P(FSFA(PJ)+RJPEQ)

    END FUNCTION fsphia

    FUNCTION fsa(pj)
      REAL(wp) :: fsa
      REAL(wp), INTENT(IN) :: pj

      FSA = FSPHIY(FSPHIA(PJ))

    END FUNCTION fsa
!!$C     
!!$C   Derivee premiere
!!$C   ----------------
!!$C     
    FUNCTION fsdphia(pj)
      REAL(wp) :: fsdphia
      REAL(wp), INTENT(IN) :: pj

      FSDPHIA = FSDFA(PJ)*FSDY0P(FSFA(PJ)+RJPEQ)

    END FUNCTION fsdphia

    FUNCTION fsda(pj)
      REAL(wp) :: fsda
      REAL(wp), INTENT(IN) :: pj

      FSDA = FSDPHIA(PJ)*FSDPHIY(FSPHIA(PJ))

    END FUNCTION fsda
!!$C     
!!$C PETIT AXE DES ELLIPSES
!!$C ======================    
!!$C     
!!$C   Fonction
!!$C   --------  
!!$C   
    FUNCTION fsphib(pj)
      REAL(wp) :: fsphib
      REAL(wp), INTENT(IN) :: pj
  
      FSPHIB = FSFB(PJ)*FSY0P(PJ+RJPEQ-1.)

    END FUNCTION fsphib

    FUNCTION fsb(pj)
      REAL(wp) :: fsb
      REAL(wp), INTENT(IN) :: pj

      FSB = FSPHIY(FSPHIB(PJ))

    END FUNCTION fsb
!!$C     
!!$C   Derivee premiere
!!$C   ----------------
!!$C   
    FUNCTION fsdphib(pj)
      REAL(wp) :: fsdphib
      REAL(wp), INTENT(IN) :: pj

      FSDPHIB = FSDFB(PJ)*FSY0P(PJ+RJPEQ-1.)+ &
          FSFB(PJ)*FSDY0P(PJ+RJPEQ-1.)

   END FUNCTION fsdphib

    FUNCTION fsdb(pj)
      REAL(wp) :: fsdb
      REAL(wp), INTENT(IN) :: pj

      FSDB = FSDPHIB(PJ)*FSDPHIY(FSPHIB(PJ))

    END FUNCTION fsdb
!!$C
!!$C
!!$C CENTRE DES ELLIPSES
!!$C ===================
!!$C     
!!$C   Fonction
!!$C   --------
!!$C
    FUNCTION fsy0(pj)
      REAL(wp) :: fsy0
      REAL(wp), INTENT(IN) :: pj

      FSY0 = (FSPHIY(FSY0P(PJ+RJPEQ-1.))-FSA(PJ))*RY0

    END FUNCTION fsy0

    FUNCTION fsphiy0(pj)
      REAL(wp) :: fsphiy0
      REAL(wp), INTENT(IN) :: pj

      FSPHIY0 = FSYPHI(FSY0(PJ))

    END FUNCTION fsphiy0
!!$C     
!!$C   Derivee premiere
!!$C   ----------------
!!$C   
    FUNCTION fsdy0(pj)
      REAL(wp) :: fsdy0
      REAL(wp), INTENT(IN) :: pj
  
      FSDY0 = RY0*(FSDY0P(PJ+RJPEQ-1.)* &
          FSDPHIY(FSY0P(PJ+RJPEQ-1.))-FSDA(PJ))

   END FUNCTION fsdy0

    FUNCTION fsdphiy0(pj)
      REAL(wp) :: fsdphiy0
      REAL(wp), INTENT(IN) :: pj

      FSDPHIY0 = FSDY0(PJ)*FSDYPHI(FSY0(PJ))

    END FUNCTION fsdphiy0
!!$C
!!$C ===================
!!$C FONCTIONS ASSOCIEES
!!$C ===================
!!$C
!!$C    Premiere fonction
!!$C    -----------------
!!$C Intersection des ellipses avec la premiere meridienne
!!$C
    FUNCTION fsy1p(pj)
      REAL(wp) :: fsy1p
      REAL(wp), INTENT(IN) :: pj

      FSY1P = 180.-FSYPHI(FSY0(FSJ(PJ))-FSA(FSJ(PJ)))

    END FUNCTION fsy1p
!!$C
!!$C
!!$C Derivee premiere
!!$C
    FUNCTION fsdy1(pj)
      REAL(wp) :: fsdy1
      REAL(wp), INTENT(IN) :: pj

      FSDY1 = -(FSDY0(FSJ(PJ))-FSDA(FSJ(PJ)))* &
                 FSDYPHI(FSY0(FSJ(PJ))-FSA(FSJ(PJ)))

   END FUNCTION fsdy1
!!$C
!!$C
!!$C Derivee seconde
!!$C
   FUNCTION fsddy1(pj)
      REAL(wp) :: fsddy1
      REAL(wp), INTENT(IN) :: pj

      FSDDY1 = 0.

    END FUNCTION fsddy1
!!$C
!!$C
!!$C    Deuxieme fonction
!!$C    -----------------
!!$C Intersection des ellipses avec la seconde meridienne
!!$C
    FUNCTION fsy2p(pj)
      REAL(wp) :: fsy2p
      REAL(wp), INTENT(IN) :: pj

      FSY2P = FSYPHI(FSY0(FSJ(PJ))+FSA(FSJ(PJ)))

    END FUNCTION fsy2p
!!$C
!!$C
!!$C Derivee premiere
!!$C
    FUNCTION fsdy2(pj)
      REAL(wp) :: fsdy2
      REAL(wp), INTENT(IN) :: pj
    
      FSDY2 = (FSDY0(FSJ(PJ))+FSDA(FSJ(PJ)))* &
                 FSDYPHI(FSY0(FSJ(PJ))+FSA(FSJ(PJ)))

   END FUNCTION fsdy2
!!$C
!!$C     
!!$C Derivee seconde
!!$C  
   FUNCTION fsddy2(pj)
      REAL(wp) :: fsddy2
      REAL(wp), INTENT(IN) :: pj

      FSDDY2 = 0.

    END FUNCTION fsddy2
!!$C
!!$C
!!$C =========================
!!$C UTILITAIRES D INTEGRATION
!!$C =========================
!!$C
!!$C
!!$C fonctionnelle definissant les ellipsess de "latitude"=cste
!!$C dans le repere Ox,Oy du plan stereographique
!!$C
    FUNCTION FSFUN(PJ,PX,PY)
      REAL(wp) :: fsfun
      REAL(wp), INTENT(IN) :: pj, px, py

      FSFUN = ( FSA(PJ)*PX )**2 + &
                       ( FSB(PJ)*(PY-FSY0(PJ)) )**2 &
                       - (FSA(PJ)*FSB(PJ))**2

   END FUNCTION FSFUN
!!$C     
!!$C fonctionnelle de derivee (repere Ox,Oy)
!!$C
    FUNCTION FSDER(PJ,PX,PY)
      REAL(wp) :: fsder
      REAL(wp), INTENT(IN) :: pj, px, py

      FSDER = (( FSB(PJ)/FSA(PJ) )**2)*( PY-FSY0(PJ)) / PX

    END FUNCTION FSDER
!!$C
!!$C
!!$C-----------------------------------------------------------------------
  !
  ! SOUTHERN HEMISPHERE FUNCTIONS :
  ! ---------------------
  !
  !	Fonctions de maillage dans l hemisphere sud
  !	Attention au changement de signe sur FSDILA et au decalage de 
  !	90 degre sur FSLAM lie a la construction particuliere de la
  !	grille.
  !
  !  A.   FSLAM : longitude en degres
  !       =====
  FUNCTION fslam(pi,pj)
    REAL(wp) fslam
    REAL(wp), INTENT(IN) :: pi, pj

    fslam = 90. - fsx0p(pi)

  END FUNCTION fslam
  !
  !
  !  B.   FSPHI : latitude en degres
  !       =====
  FUNCTION fsphi(pi,pj)
    REAL(wp) fsphi
    REAL(wp), INTENT(IN) :: pi, pj

    fsphi = fsy0p(pj)

  END FUNCTION fsphi
  !
  !
  !  C.   FSDILA, FSDJLA : Derivees premieres de FSLAM par rapport a I et J
  !       ==============
  !
  FUNCTION fsdila(pi,pj)
    REAL(wp) fsdila
    REAL(wp), INTENT(IN) :: pi, pj

    fsdila = - fsdx0p(pi)

  END FUNCTION fsdila

  FUNCTION fsdjla(pi,pj)
    REAL(wp) fsdjla
    REAL(wp), INTENT(IN) :: pi, pj

    fsdjla = 0.e0

  END FUNCTION fsdjla
  !
  !
  !  D.   FSDIPH, FSDJPH : Derivee premieres de FSPHI par rapport a I et J
  !       ==============
  !
  FUNCTION fsdiph(pi,pj)
    REAL(wp) fsdiph
    REAL(wp), INTENT(IN) :: pi, pj

    fsdiph = 0.e0

  END FUNCTION fsdiph

  FUNCTION fsdjph(pi,pj)
    REAL(wp) fsdjph
    REAL(wp), INTENT(IN) :: pi, pj

    fsdjph = fsdy0p(pj)

  END FUNCTION fsdjph
  !


SUBROUTINE set_coefficients


  REAL(wp) :: zra0, zra1, zf0, zf1, zf2, zval, zy0pjn

  INTEGER :: icompt


  !-----------------------------------------------------------------------
  !      
  ! VALEUR DES COEFFICIENTS DES FUNCTIONS
  ! 
  !-----------------------------------------------------------------------
  !
  ! Position du pole de la grille en degres
  !     
  RLAMPO  = 73.+180.  !  grid starts at RLAMP0 - 180
  RPHIPO  = 90.
  NROTA   = INT(RLAMPO)
  !
  !
  ! Conversion reelle des parametres de base
  !
  RJPEQT  = float(JPEQT )
  RJPEQ   = float(JPEQ  )
  RJPNP   = float(JPNP  )
  ! ... indice de la ligne des poles (le +0.5 correspond a une ligne V-F)
  rjpnth  = FLOAT(jpnord) + 0.5
  rjpnthm1= rjpnth - 1.
  ! ... indice de la ligne des poles pour la grille globale
  rjpj    = FLOAT(jpj) + 0.5
  !
  ! Pas en latitude de la grille 
  !     
  RDX     = 360./float(JPI-1)
  !
  ! Position du pole sud de la grille
  !
  RPHISUD = FSY0P(1.0_wp)
  !
  ! Position reelle de l equateur de la grille en degres
  !
  RPHIEQ  = FSY0P(RJPEQ)
  ! Rayon de l equateur de la grille dans le plan stereographique
  ! (on se place dans un repere de telle sorte que l equateur soit
  ! de rayon 1)
  !
  REQ     = tan(RPI/4.-RAD*RPHIEQ/2.)
  !
  ! Coordonnee du pole nord de la grille
  !
  RPOL    = tan(RPI/4.-RAD*RPHIPO/2.)
  !
  !-----------------------------------------------------------------------
  !
  ! Coefficients des fonctions de bases/
  !
  !-----------------------------------------------------------------------
  !
  ! Positions finales des deux poles
  ! rphi1 sur l Asie ; rphi2 sur le Canada
  !
  rphi1 = 50.       ! 50N
  !      rphi2 = 108
  rphi2 = 180.-66.  ! 66N
  !
  ry1 = fsphiy(rphi1)
  ry2 = fsphiy(rphi2)
  !
  !    Premiere fonction de base
  !    -------------------------
  rphifa  = fsyphi( 0.5*ABS( ry1-ry2 ) )
  ra1     = 1e-5
  !
  ! Dichotonie pour trouver RA2 tel que FSPHIA(JPNORD)=RPHIFA
  zra0    = 1e-1
  zra1    = 10.
  icompt=0
  WRITE (0,*) 'rzero=',rzero
  !
  DO WHILE (ABS(zf2-rphifa).GE.rzero)
     ra2=zra0
     zf0=fsphia(rjpnth)
     ra2=zra1
     zf1=fsphia(rjpnth)
     ra2=(zra0+zra1)/2.
     zf2=fsphia(rjpnth)
     zval=(zf0-rphifa)*(zf2-rphifa)
     icompt=icompt+1
     !      WRITE (0,*) icompt,'  zra0=',zra0,'   zf0=',zf0,
     !     $     '   zra1=',zra1,'   zf1=',zf1
     IF ((ABS(zra0-zra1).LT.rzero).OR.(icompt.GE.2000)) THEN
        WRITE (0,*) 'Fin anormale de convergence'
        WRITE (0,*) 'zra0=',zra0,'   zra1=',zra1
        WRITE (0,*) 'ABS(zra1-zra2)=',ABS(zra1-zra1)
        WRITE (0,*) 'fsphiamax=',fsphia(rjpnth)
        WRITE (0,*) 'fsphiamax-fsphia=',ABS(fsphia(rjpnth)-rphifa)
        EXIT
     END IF

     IF ((zf0-rphifa)*(zf2-rphifa).GT.0) THEN
        zra0=ra2
     ELSE
        zra1=ra2
     ENDIF
  END DO
  WRITE (0,*) 'Calcul de ra2 effectue : ra2=',ra2
  WRITE (0,*) 'fsa = ',fsa(rjpnth)
  !
  !    Deuxieme fonction de base
  !    -------------------------
  zy0pjn  = fsy0p(rjpj)
  !
  rphifb  = 90.
  rbmax   = ( rphifb - zy0pjn ) / zy0pjn
  rb2     = 1.
  rb1     = rbmax / TANH(rb2)
  !
  !    Troisieme fonction de base
  !    --------------------------
  rphify0 = fsyphi( 0.5*( ry1+ry2 ) )
  ry0     = fsphiy(rphify0) / ( fsphiy(zy0pjn) - fsa(rjpnth) )
  !
  !=======================================================================
END SUBROUTINE set_coefficients


END MODULE functions

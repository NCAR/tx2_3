PROGRAM tripole
  USE kinds
  USE nf90util
  USE param
  USE common
  USE functions
  USE hgrnth
  USE hgrctl
  USE mshtrp , ONLY: msh_trp
  USE mshnth
  USE mshsth
  USE gridio

IMPLICIT NONE

!!!---------------------------------------------------------------------
!!!
!!!                       PROGRAM GLOHGR
!!!                     ******************
!!!
!!!  PURPOSE : 
!!!  --------- 
!!!	CONSTRUCTION DE LA GRILLE DU MODEL GLOBAL: LE POLE NORD
!!!	DE LA GRILLE DU MODELE EST DEPLACE A LA LATITUDE RPHIPO 
!!!	       ET A LA LONGITUDE RLAMPO.
!!!
!!
!!   MODIFICATIONS:
!!   --------------
!!       ORIGINAL  : 93-13 (G. MADEC)
!!!--------------------------------------------------------------------- 
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1/07/99) 
!!!---------------------------------------------------------------------
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
!! mshglo local declarations
!! ==================
      INTEGER :: ijpi, ijpj, ijpeq, ijpnord, ijpin, ijpjn
      INTEGER :: ilu, ilseq, numin
      INTEGER :: iused(1,100)

      CHARACTER(LEN=5) :: clfield
      CHARACTER(LEN=21) :: clold, clfor, clseq, clnew
      CHARACTER(LEN=32) :: clname
      CHARACTER(LEN=21) :: cldir, clunf, clunk

      REAL(wp) :: zra
      REAL(wp) :: zlu (jpin,jpjn), zani(jpim,jpjm)
      REAL(wp) :: zout(jpim,jpjm)

      LOGICAL :: lwp

      INTEGER :: ncid, jd, jv
      TYPE(ncdimtype), DIMENSION(2) :: dims
      TYPE(ncvartype), DIMENSION(4) :: vars
      TYPE(ncvartype), DIMENSION(7) :: paras
!!----------------------------------------------------------------------
  !
  !
  ! 0. Initialisation
  ! =================
  !
  !
  ! 0.1 Parametres de controle (set in MODULE param)
  ! --------------------------
  !     
  write(0,*)
  write(0,*)
  write(0,*) '    Parametres de controle:'
  write(0,*) '    -----------------------'
  write(0,*) 
  write(0,*) 'Parametre de debugage                      ', ndebug
  write(0,*)
  write(0,*)
  !
  !
  ! 0.2 Constantes physiques
  ! ------------------------
  !
  write(0,*) '    Constantes mathematiques'
  write(0,*) '    ------------------------'
  write(0,*)
  write(0,*) 'Rayon de la terre                   ra    = ',ra
  write(0,*) 'PI                                  pi    = ',rpi
  write(0,*) 'Conversion degre radian             rad   = ',rad
  write(0,*) 'Zero machine                        rzero = ',rzero
  write(0,*)
  write(0,*)
  !
  ! 0.3 Parametres des fonctions
  ! ----------------------------
  !
  !
  CALL set_coefficients
  !
  !
  ! 0.4 Parametres de la grille
  ! ---------------------------
  !
  write(0,*) '    Parametres de la grille'
  write(0,*) '    -----------------------'
  write(0,*)
  write(0,*) 'Taille de la grille globale:'
  write(0,*) '    Nb de pts de grille -zonal-             ', JPI
  write(0,*) '    Nb de pts de grille -meridien-          ', JPJ
  write(0,*) '    Nb de pts de grille partie sud          ', JPEQ
  write(0,*) '    Nb de pts de grille partie nord         ', JPNORD
  write(0,*)
  write(0,*) 'Latitude et longitude du pole nord modele:'
  write(0,*) '    RLAMPO =                                ', RLAMPO
  write(0,*) '    RPHIPO =                                ', RPHIPO
  write(0,*) 'Latitude du premier pole                    ', RPHI1
  write(0,*) 'Latitude du second pole                     ', RPHI2
  write(0,*) 'Position du pole nord modele sur l axe oy   ', RPOL
  write(0,*)
  write(0,*) 'Position de l equateur de la grille         ', RPHIEQ
  write(0,*) 'Rayon de l equateur (plan stereographique)  ', REQ
  write(0,*)
  write(0,*) 'Limite sud de la grille                     ', RPHISUD
  !
  !
  ! I. Compute or read the north half-grid
  ! ======================================
  !
  IF ( l_nth_calc ) THEN
     !
     ! calculate the north part and save to file
     ! I.1 Calcul de la grille globale et ses facteurs d echelle
     ! ---------------------------------------------------------
     !
     !     Controle des caracteristiques de la grille a construire
     !     - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     CALL hgr_ctl


     !     Construction de la grille de l hemisphere nord
     !     - - - - - - - - - - - - - - - - - - - - - - - -

     CALL hgr_nth


     !     Save north grid to file
     !     - - - - - - - - - - - -
     CALL write_nth

  ELSE

     ! Read the north part from file
     ! I.2 Lecture de la grille globale et ses facteurs d echelle
     ! ----------------------------------------------------------

     CALL read_nth

  ENDIF


  !
  ! II. Construct the global grid
  ! =============================

  ! Allocate the global arrays under the assumption that there will be 
  ! no tropical stretching
  ALLOCATE ( glam(jpim,jpjm,4), gphi(jpim,jpjm,4), &
               e1(jpim,jpjm,4),   e2(jpim,jpjm,4)  )

  ! Compute the southern part of the grid (including unstretched tropics)
  CALL msh_sth

  ! Write the northern part to the global arrays
  CALL msh_nth( 1, glamnth, glam )
  CALL msh_nth( 2, gphinth, gphi )
  CALL msh_nth( 3, e1nth, e1 )
  CALL msh_nth( 4, e2nth, e2 )

  ! Add streched tropics
  IF ( l_trop_stretch ) CALL msh_trp

  !  Write the global arrays to file
  CALL write_glo
  CALL write_hgrid

  DEALLOCATE( glam, gphi, e1, e2 )

END PROGRAM tripole

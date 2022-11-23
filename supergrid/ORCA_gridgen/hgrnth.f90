MODULE hgrnth
  USE kinds
  USE param
  USE common
  USE functions
  USE nthcar
  USE OMP_LIB
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------


  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE hgr_nth
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE hgr_nth
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!	Construction de la grille du model global dans l hemisphere
!!!	nord : le pole nord de la grille du modele est deplace a la
!!!	latitude rphipo et a la longitude rlampo.
!!!
    !!
    !!   Method :
    !!   --------
    !!
    !!
    !!   Input :
    !!   -------
    !!      argument                : NO
    !!      common                  : NO
    !!
    !!   Output :
    !!   -------
    !!      file
    !!              numgri=30       : horizontal grid and scale factors
    !!                                (file= glohgr.output )
    !!
    !!   External :
    !!   ----------
    !!              nth_car		: grille fine en coord. cartesiennes
    !!		nthgeo		: grille nord en coord. geographiques
    !!		nthfac		: facteurs d echelle
    !!
    !!   Modifications :
    !!   ---------------
    !!       original  : 93 (G. Madec)
    !!       update    : 99 (G. Madec) new computation (2 poles)
    !!
    !!----------------------------------------------------------------------

    !! local declarations
    !! ==================
    INTEGER :: jn, ji, ipart, ideb, ifin, numwri
    INTEGER :: indic
    !
    REAL(wp) :: zglam(10,jpjn), zgphi(10,jpjn), ze1(10,jpjn), ze2(10,jpjn)
    !
    CHARACTER(LEN=32) :: clname
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1999)
!!!---------------------------------------------------------------------
    !
    !
    ! 0. VALEURS DE DEPART SUR LE DEMI CERCLE EQUATEUR
    ! ------------------------------------------------
    !
    ! Recherche du nombre maximal de dicotonie tel que la precision atteinte
    ! soit de RZERO (maximum d iteration imposi a 100)
    !
    ncompt = 100
    indic  = 0
    DO jn = 1, 100
       IF( (indic .EQ. 0) .AND. ((0.5**jn)*(jpjn-1)) .LT. rzero ) then
          ncompt=jn
          indic=1
       ENDIF
    END DO
    WRITE(0,*) 'hgr_nth: nb max de dicotonie    ncompt= ',ncompt
    WRITE(0,*) '        precision requise      rzero = ',rzero
    !
    !
    ! I. Construction de la demi-grille nord et ses facteurs d echelle
    ! -----------------------------------------------------------------
    !
    WRITE(0,*)
    WRITE(0,*) '    Appel de nth_car'
    WRITE(0,*) '    ---------------'
    WRITE(0,*) 
    !
    ! ... boucle sur les courbes de pseudo-longitude=cstes
!$OMP PARALLEL DO PRIVATE(ji)
    DO ji = 1, jpin
       WRITE(0,*) 'ligne no = ', ji
       CALL nth_car(ji)
    END DO
!$OMP END PARALLEL DO
    !!


    ! The following code is for computing sections of the grid
    ! as promted at run time (include as an option?)
    !      clname = 'hgr_nth_p'
    !      READ(5,*) ipart
    !      WRITE(0,*) ' partie ', ipart
    !      ideb = 10*(ipart-1)+1
    !      ifin = 10*(ipart-1)+10
    !!   ... derniere partie, que premier ligne!
    !      ifin = MIN(ifin,jpin)
    !
    !      DO ji = ideb, ifin
    !        WRITE(0,*) 'ligne no = ', ji
    !        CALL nth_car(ji)
    !      END DO


!!$ To do:
!!$ the following code used to be executed by msh_glo after the north
!!$ grid was read in.  Work out why it's needed and perhaps
!!$ reinstate it here (translation to new vars is below).
!!$
!!$        write(0,*) 'bouge pole nord amerique...' ! move the pole
!!$        IF ( jtyp .EQ. 1 ) zlu(jpin,jpjn) = zlu(jpin,jpjn-1)
!!$        IF ( jtyp .EQ. 2 ) zlu(jpin,jpjn) = 66.   !!!???!!!
!!$        IF ( jtyp .EQ. 3 ) zlu(jpin,jpjn) = zlu(jpin,jpjn-1)
!!$        IF ( jtyp .EQ. 4 ) zlu(jpin,jpjn) = zlu(jpin,jpjn-1)
!!$
!!$ Translates to:
!!$   glamnth(jpin,jpjn) = glamnth(jpin,jpjn-1)
!!$   gphinth(jpin,jpjn) = 66. ! hardwired
!!$   e1nth(jpin,jpjn)   = e1nth(jpin,jpjn-1)
!!$   e2nth(jpin,jpjn)   = e2nth(jpin,jpjn-1)



    !!
    !!  Binary output code has been moved to write_nth in module gridio


  END SUBROUTINE hgr_nth
END MODULE hgrnth

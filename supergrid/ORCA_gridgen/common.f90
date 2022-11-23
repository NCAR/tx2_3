MODULE common
  USE kinds
  USE param
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------

  
IMPLICIT NONE

  !!
  !!----------------------------------------------------------------------
  !!
  !!      rlampo  : Longitude du centre de la partie nord de la grille 
  !!      rphipo  : Latitude du centre de la partie nord de la grille
  !!      rphieq  : Latitude de l equateur de la grille
  !!      rphisud : Latitude la plus au sud de la grille 
  !!      req     : Rayon de l equateur dans le plan stereographique
  !!      rpol    : Ordonee de l equateur dans le plan stereographique
  !!      rjpnth  : Valeur reelle de l indice dans la partie nord
  !!      rdx     : Pas de la grille en latitude 
  !!     
  REAL(wp) :: rlampo, rphipo, rphieq, rphisud, req, rpol, rjpnth, rdx

  !!
  !!----------------------------------------------------------------------
  !!
  !! COMMON comdgn: demi-grille nord (points T, U, V et F)
  !! --------------------------------
  !!      glam,gphi : lat/lon en degre
  !!      e1, e2    : facteur d echelle (metres)
  !!
  REAL(wp) ::  glamnth(jpin,jpjn), gphinth(jpin,jpjn)
  REAL(wp) ::  e1nth  (jpin,jpjn), e2nth  (jpin,jpjn)

  !! COMMON commsh: global mesh (points T, U, V et F)
  !! ------------------------------------------------
  !!      glam,gphi : lat/lon en degre pour les 4 points T, U, V et F
  !!      e1, e2    : facteur d echelle (m) pour les 4 points T, U, V et F
  !!
  REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: glam, gphi, e1, e2

  !!
  !!----------------------------------------------------------------------
END MODULE common

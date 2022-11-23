MODULE kinds
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------

  IMPLICIT NONE
  
  INTEGER, PUBLIC, PARAMETER ::          &     !: Floating point section
       sp = SELECTED_REAL_KIND( 6, 37),  &     !: single precision (real 4)
       dp = SELECTED_REAL_KIND(12,307)         !: double precision (real 8)

#if defined key_64bit
  INTEGER, PUBLIC, PARAMETER ::  wp =  dp      !: working precision
#else
  INTEGER, PUBLIC, PARAMETER ::  wp =  sp      !: working precision
#endif

END MODULE kinds

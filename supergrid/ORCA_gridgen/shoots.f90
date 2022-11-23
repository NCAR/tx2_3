MODULE shoots
  USE kinds
  USE trop, ONLY: integ
IMPLICIT NONE
!!----------------------------------------------------------------------
!! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
!!----------------------------------------------------------------------
! instructions to use as general utility:
! (to do: add matrix solver for >1d)
!
! Change above USE statement to include the module+routine which
! defines the function y(x).  Change subroutine XtoY below to call the
! routine which defines the function y(x).


INTERFACE XtoY
   MODULE PROCEDURE XtoY, XtoY_scalar
END INTERFACE

INTERFACE newton
   MODULE PROCEDURE newton, newton_scalar
END INTERFACE

REAL(wp), PARAMETER :: Ytol = 1.0e-6, Xcorrtol = 1.0e-6


CONTAINS

! would be cleaner to convert XtoY into a function Y(X)
SUBROUTINE XtoY (X, Y)
  ! subroutine which supplies vector Y given vector X
  ! X, Y same size
  REAL(wp), DIMENSION(:), INTENT(IN) :: X
  REAL(wp), DIMENSION(SIZE(X)), INTENT(OUT) :: Y

  CALL integ (X(1), Y(1))

END SUBROUTINE XtoY


SUBROUTINE newton (X, Y, con, iters)
    ! newton-raphson rootfinder for a system of functions Y(X) 
    ! X and Y are equal sized vectors
  REAL(wp), DIMENSION(:), INTENT(INOUT) :: X
  REAL(wp), DIMENSION(:), INTENT(OUT) :: Y
  INTEGER, INTENT(OUT) :: con, iters
  REAL(wp), DIMENSION(SIZE(X)) :: Xcorr
  REAL(wp), DIMENSION(SIZE(X),SIZE(X)) :: dY
  INTEGER :: i, imax=100, j, jmax = 4, bad
  LOGICAL, PARAMETER :: testbad = .FALSE.
  con = 0

  DO i  = 1, imax
     CALL XtoYdY (X, Y, dY)

!     CALL solve (dY, Xcorr, -Y)   ! Solve matrix equation for Xcorr 
! In this case the unknown vector is 1D so solve is simply:
     Xcorr(1) = -Y(1) / dY(1,1)

     X(:) = X(:) + Xcorr(:)

     ! test new iterate for system badness
     IF (testbad) THEN
        DO j = 1, jmax  
          ! CALL baddie (0.0, load(X), bad, "newton")
           IF (bad == 1) THEN
              WRITE(13,*) 'warning newton: bad new iterate. X = ',X,&
                   & 'Trying smaller correction.'
              X = X -  Xcorr / 2**j
           ELSE
              EXIT
           END IF
        END DO
        IF (j == jmax+1) STOP 'newton: bad new iterate, jmax too small'
     END IF ! testbad


     ! possible exits from DO loop

     IF (SUM(ABS(Y)) < Ytol) THEN
        con = 1    ! Y convergence
        iters = i
        EXIT
     ELSE IF (SUM(ABS(Xcorr/X)) < Xcorrtol) THEN
        con = 2    ! X convergence
        iters = i
        EXIT
     END IF

  END DO


  IF (con == 0) THEN
     WRITE(*,*) 'Warning newton: N-R not converged. iters =',imax
     con = -1
  END IF


END SUBROUTINE newton
SUBROUTINE newton_scalar (X, Y, con, iters)
  ! scalar wrapper for subroutine newton
  REAL(wp), INTENT(INOUT) :: X
  REAL(wp), INTENT(OUT) :: Y
  INTEGER, INTENT(OUT) :: con, iters
  
  REAL(wp), DIMENSION(1) :: Xarr, Yarr
  
  Xarr(1) = X
  
  CALL newton (Xarr, Yarr, con, iters)
  
  X = Xarr(1)
  Y = Yarr(1)

END SUBROUTINE newton_scalar



SUBROUTINE XtoYdY (X, Y, dY)
  ! Uses XtoY to get Y and dY/dX for N-R
  REAL(wp), DIMENSION(:), INTENT(IN) :: X
  REAL(wp), DIMENSION(:), INTENT(OUT) :: Y
  REAL(wp), DIMENSION(SIZE(X),SIZE(Y)), INTENT(OUT) :: dY
  REAL(wp), DIMENSION(SIZE(X)) :: d, X1, Y1
  REAL(wp), PARAMETER :: incr=1.0e-6  ! increment for derivs (tunable?)
  INTEGER :: i

  CALL XtoY (X, Y)

  d = incr * X
  DO i = 1, SIZE(X)  ! find dY / dX(i)
     
     X1 = X
     X1(i) = X(i) + d(i)

     CALL XtoY (X1, Y1)

     dY(i,:) = (Y1 - Y) / d(i)
     
  END DO


END SUBROUTINE XtoYdY


SUBROUTINE XtoY_scalar (x, y)
  ! scalar front end for XtoY, overloaded
  REAL(wp), INTENT(IN) :: x
  REAL(wp), INTENT(OUT) :: y
  REAL(wp), DIMENSION(1) :: Xarr, Yarr

  Xarr(1) = x
  CALL XtoY (Xarr, Yarr)
  y = Yarr(1)

END SUBROUTINE XtoY_scalar




SUBROUTINE newtsurf (X, Y, con, iters)
  !newton-raphson rootfinder for a single function Y of several
  ! variables X.
  ! For equal-sized X and Y, use routine newton below.
  REAL(wp), DIMENSION(:), INTENT(INOUT) :: X
  REAL(wp), DIMENSION(1), INTENT(OUT) :: Y
  INTEGER, INTENT(OUT) :: con, iters
  REAL(wp), DIMENSION(SIZE(X)) :: Xcorr
  REAL(wp), DIMENSION(SIZE(X),SIZE(Y)) :: dY
  REAL(wp) :: magdY, Yerr
  INTEGER :: itmax=100, j, jmax=4, bad
  LOGICAL, PARAMETER :: testbad = .FALSE.

  con = 0

  CALL XtoYdY (X, Y, dY)
  IF (SUM(ABS(Y)) < Ytol) THEN
     iters = 0
     con = 1
  END IF


  DO iters = 1, itmax

     ! use dY to head toward the nearest root
     magdY = SQRT(SUM (dY**2))  
     Xcorr = - Y/magdY**2 * dY(:,1)

     X = X + Xcorr


     IF (testbad) THEN
        ! for zprof: test new iterate for system badness
        DO j = 1, jmax  
          ! CALL baddie (0.0, load(X), bad, "newton")
           IF (bad == 1) THEN
              WRITE(13,*) 'warning newton: bad new iterate. X = ',X,&
                   & 'Trying smaller correction.'
              X = X -  Xcorr / 2**j
           ELSE
              EXIT
           END IF
        END DO
        IF (j == jmax+1) STOP 'newton: bad new iterate, jmax too small'
     END IF ! testbad

     CALL XtoYdY (X, Y, dY)

     ! possible exits from DO loop
     IF (SUM(ABS(Y)) < Ytol) THEN
        con = 1    ! Y convergence
        EXIT
     ELSE IF (SUM(ABS(Xcorr/X)) < Xcorrtol) THEN
        con = 2    ! X convergence
        EXIT
     END IF

  END DO


  IF (con == 0) THEN
     WRITE(*,*) 'Warning newton: N-R not converged. iters =',itmax
     con = -1
  END IF


END SUBROUTINE newtsurf




SUBROUTINE bisection (X, Y, con, iters)
  ! bisection method to find root of Y(X). X, Y 1D vectors
  ! Global varible used: Ytol, Xcorrtol

  REAL(wp), DIMENSION(1), INTENT(INOUT) :: X
  REAL(wp), DIMENSION(1), INTENT(OUT) :: Y
  INTEGER, INTENT(OUT) :: con, iters
  REAL(wp) :: right, left, mid, Yright, Yleft, Ymid, range
  REAL(wp), PARAMETER :: jump=1.2
  INTEGER :: i, findmax = 100, itmax = 100
  REAL(wp), DIMENSION(SIZE(X)) :: miss
  con = 0

  right = X(1)
  left = X(1)

  ! expand interval until root is enclosed
  DO i = 1, findmax  
     left = left / jump
     right = right * jump

     CALL XtoY (left, Yleft)
     CALL XtoY (right, Yright)
     IF (Yleft*Yright < 0) EXIT
  END DO
  IF (i == findmax+1)  THEN
     con = -1
     WRITE(*,*) "bisection: root not enclosed by interval:",left,right 
     RETURN
  END IF

  ! Perform the bisection
  DO i = 1, itmax
     ! find midpoint and evaluate Ymid
     mid = (left + right)/2.0 
     CALL XtoY (mid, Ymid)


     ! form new interval containing root
     IF (Ymid*Yright < 0) THEN
        left = mid
        Yleft = Ymid
     ELSE
        right = mid
        Yright = Ymid
     END IF
     range = right - left


     ! possible exits
     IF (ABS(Ymid/mid) < Ytol) THEN
        con = 1    ! Y convergence
        iters = i
        X = mid
        Y = Ymid
        EXIT
     ELSE IF (ABS(range/mid) < Xcorrtol) THEN
        con = 2    ! X convergence
        iters = i
        X = mid
        Y = Ymid
        EXIT
     END IF

  END DO


  IF (con == 0) THEN
     WRITE(*,*) "Warning bisection: not converged after",itmax,"iterations"
     con = -1
     X = mid
     Y = Ymid
  END IF

END SUBROUTINE bisection



END MODULE shoots



!!$SUBROUTINE finder (a, b, h)
!!$  REAL(wp), INTENT(IN) :: a, b, h
!!$  REAL(wp), DIMENSION(3) :: S
!!$  REAL(wp), DIMENSION(1) :: miss, lmiss, v
!!$  INTEGER :: n
!!$
!!$
!!$
!!$  DO n = 1, 1000
!!$
!!$     v = 0.001* REAL(n,KIND=wp)
!!$   
!!$     CALL load (v, S)
!!$     CALL odeint (S, a, b, h) 
!!$     CALL missed (S, miss)
!!$
!!$     IF (miss(1)*lmiss(1) <= 0) THEN
!!$        WRITE(*,*) v, miss
!!$!        WRITE(UNIT=10,FMT = '(2e12.4)') v, miss
!!$     END IF
!!$
!!$     WRITE(UNIT=11,FMT = '(2e12.4)') v, miss
!!$     lmiss = miss
!!$
!!$  END DO
!!$END SUBROUTINE finder





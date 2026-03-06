program test_grid
  use kinds
  use constants
  use ncdf_wrapper
  use mom6_grid
  use sorting

  implicit none

  interface
     subroutine sort(x,N)
       use kinds
       real(kind=real_kind), dimension(1:), intent(inout) :: x
       integer, intent(in)                   :: N
     end subroutine sort
     recursive subroutine quicksort(a)
       use kinds
       implicit none
       real(kind=real_kind) :: a(:)
     end subroutine quicksort
  end interface
  !-----------------------------------------------------------------------

  character (char_len) ::   mom6_grid_version, mom6_horiz_grid_file

  logical :: grid_is_tripole=.true.


  real(kind=real_kind), dimension(10) :: x = (/ 9., 8., 7., 6., 5., 4., 3., 2., -10., -20./)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !    Work variables
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  integer (kind=int_kind) :: i,j,k,im1,jm1

  namelist /model_grid_nml/ mom6_grid_version, mom6_horiz_grid_file, grid_is_tripole

  !-----------------------------------------------------------------------
  !     Initialize and get namelist input
  !-----------------------------------------------------------------------

!!$  call init_constants
!!$
!!$  mom6_grid_version = 'unknown-grid-version'
!!$  mom6_horiz_grid_file = 'unknown-horiz-grid'
!!$
!!$  read(*,nml=model_grid_nml)
!!$  write(*,nml=model_grid_nml)
!!$
!!$  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$  ! Get model super grid and allocate model grid arrays
!!$  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!$
!!$  call read_super_grid(mom6_horiz_grid_file)
!!$  call compute_mom6_grid()
!!$
!!$  print *,' nx_o=',nx_o,' nx_o/4=',nx_o/4,'  nx_o/2=',nx_o/2
!!$
!!$  print *,' LATH(ny_o-3:) = ',LATH(ny_o-3:)
!!$  print *,' LATQ(ny_o-3:) = ',LATH(ny_o-3:)

  print *,'before sort:',x
  call sort(x,10)
  print *,'after sort:',x

  x = (/ 9., 8., 7., 6., 5., 4., 3., 2., -10., -20./)
  print *,'before quicksort:',x
  call quicksort(x(1:8))
  print *,'after quicksort:',x


end program test_grid

SUBROUTINE  Sort(x, N)
  use kinds

  IMPLICIT  NONE
  real(kind=real_kind), DIMENSION(1:), INTENT(INOUT) :: x
  INTEGER, INTENT(IN)                   :: N
  INTEGER                               :: i
  INTEGER,dimension(1)                  :: imin
  real(kind=real_kind) :: temp

  DO i = 1, N-1             ! except for the last
     imin = minloc(x(i:N))  ! find min from this to last

     temp = x(i)             ! swap this and min
     x(i) = x(imin(1)+i-1)
     x(imin(1)+i-1) = temp
  END DO
END SUBROUTINE  Sort
! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a)
  use kinds
  implicit none
  real :: a(:)
  real x, t
  integer :: first = 1, last
  integer i, j

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last))
end subroutine quicksort

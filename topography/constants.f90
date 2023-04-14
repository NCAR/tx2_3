      module constants

      use kinds

      implicit none
      save

!-----------------------------------------------------------------------
!
!     physical constants (all in cgs units except for those in MKS)
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), parameter ::     &
     &  grav      = 9.806_dbl_kind,          & ! gravit. accel. (cm/s**2)
     &  omega     = 7.292123625e-5_dbl_kind, & ! angular vel. of Earth 1/s
     &  radius    = 6370.0e3_dbl_kind          ! radius of Earth (m)

      real(kind=dbl_kind), parameter :: &
           & sec_per_minute = 60.0d0,&
           & minute_per_hour = 60.0d0,&
           & hour_per_day = 24.0d0,&
           & day_per_year = 365.0d0,&
           & sec_per_hour = sec_per_minute*minute_per_hour,&
           & sec_per_day = sec_per_hour*hour_per_day,&
           & sec_per_year = sec_per_day*day_per_year

!-----------------------------------------------------------------------
!
!     numbers
!
!-----------------------------------------------------------------------

      character (char_len) :: char_blank          ! empty character string

      real (kind=real_kind), allocatable, dimension(:,:) :: ONE_R
      real (kind=dbl_kind), allocatable, dimension(:,:) :: ONE_D

      real (kind=dbl_kind), parameter ::  &
     &  c0   = 0.0_dbl_kind, &
     &  c1   = 1.0_dbl_kind, &
     &  c1p5 = 1.5_dbl_kind, &
     &  c2   = 2.0_dbl_kind, &
     &  c3   = 3.0_dbl_kind, &
     &  c4   = 4.0_dbl_kind, &
     &  c5   = 5.0_dbl_kind, &
     &  c8   = 8.0_dbl_kind, &
     &  c10  = 10.0_dbl_kind, &
     &  c16  = 16.0_dbl_kind, &
     &  c1000= 1000.0_dbl_kind, &
     &  p33  = c1/c3,           &
     &  p5   = 0.5_dbl_kind, &
     &  p25  = 0.25_dbl_kind, &
     &  p125 = 0.125_dbl_kind, &
     &  p001 = 0.001_dbl_kind, &
     &  eps  = 1.0e-10_dbl_kind, &
     &  eps2  = 1.0e-20_dbl_kind

      real (kind=dbl_kind) :: pi, pih, pi2         ! pi, pi/2 and 2pi

!-----------------------------------------------------------------------
!
!     conversion factors
!
!-----------------------------------------------------------------------

      real (kind=dbl_kind), parameter ::     &
     &  T0_Kelvin     = 273.16_dbl_kind      & ! zero point for Celcius
     &, mpercm        = .01_dbl_kind         & ! meters per cm
     &, cmperm        = 100._dbl_kind        & ! cm per meter
     &, salt_to_ppt   = 1000._dbl_kind       & ! salt (g/g) to ppt
     &, ppt_to_salt   = 1.e-3_dbl_kind       & ! salt ppt to g/g
     &, mass_to_Sv    = 1.0e-12_dbl_kind     & ! mass flux to Sverdrups   
     &, heat_to_PW    = 4.186e-15_dbl_kind   & ! heat flux to Petawatts 
     &, salt_to_Svppt = 1.0e-9_dbl_kind      & ! salt flux to Sv*ppt
     &, salt_to_mmday = 3.1536e+5_dbl_kind     ! salt to water (mm/day)

      real (kind=dbl_kind) :: radian           ! degree-radian conversion


!***********************************************************************

      contains

!***********************************************************************

      subroutine init_constants

      implicit none

!-----------------------------------------------------------------------
!
!     This subroutine initializes constants that are best defined
!     at run time (e.g. pi).
!
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: n, ierr

!-----------------------------------------------------------------------

      pi  = c4*atan(c1)
      pi2 = c2*pi
      pih = p5*pi

      radian = 180.0_dbl_kind/pi

      do n=1,char_len
        char_blank(n:n) = ' '
      end do

!-----------------------------------------------------------------------

      end subroutine init_constants

!***********************************************************************

      end module constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

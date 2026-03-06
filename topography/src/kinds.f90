!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      module kinds

!***********************************************************************
!
!     This module defines variable precision for all common data
!     types.
!
!     CVS:$Id: kinds_mod.F,v 1.2 2000/05/05 18:57:27 pwjones Exp $
!     CVS:$Name: POP_1_2 $
!
!-----------------------------------------------------------------------

      implicit none
      save

!-----------------------------------------------------------------------

      integer, parameter :: char_len  = 256,                     &
     &                      int_kind  = kind(1),                &
     &                      log_kind  = kind(.true.),           &
     &                      real_kind = selected_real_kind(6),  &
     &                      dbl_kind  = selected_real_kind(13)

      integer, parameter :: r8 = dbl_kind,&
           &                r4 = real_kind

!***********************************************************************

      end module kinds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

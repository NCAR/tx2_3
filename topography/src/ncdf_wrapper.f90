      module ncdf_wrapper

      include 'netcdf.inc'
      contains

      subroutine handle_err (status, string)

      implicit none

      integer status
      character*(*) string

      write (*,*) nf_strerror(status),': ',string

      return

      end subroutine handle_err

      end module ncdf_wrapper

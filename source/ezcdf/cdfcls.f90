! automatic conversion to free f90 compatible form
! free.pl cdfcls.for
! linewidth: 72
! file names: cdfcls.for
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
      subroutine cdfcls(ncid)
      include "netcdf.inc"
      INTEGER ncid,status
      status = nf_close(ncid)
      call handle_err(status,' ','cdfcls','nf_close')
      return
      end
!     _____ error handling.  If there is any error the execution of
!     ----- the program is stopped

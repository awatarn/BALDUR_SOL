! automatic conversion to free f90 compatible form
! free.pl cdfopn.for
! linewidth: 72
! file names: cdfopn.for
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
      subroutine cdfopn(ncid,filename,opt)
! Create/open cdf file
! 03/09/99 C.Ludescher
!
      include "netcdf.inc"
      INTEGER,       intent(in) :: ncid
      character*(*), intent(in) :: filename
      character*1,   intent(in) :: opt
      integer status
 
!     --- changed "noclobber" to "clobber" - J. Menard 10/12/98
!      status = nf_create(filename,nf_noclobber,ncid)
 
      if (opt .eq. 'w' .or. opt .eq. 'W') then
         status = nf_create(filename,nf_clobber,ncid)
         call handle_err(status,filename,'cdfcrt','nf_create')
      else
! Open for read
         status = nf_open(filename, nf_nowrite, ncid)
         call handle_err(status,'filename','cdfopn','nf_open')
      end if
      return
      end

      subroutine ezcdf_open(ncid,filename,opt,ier)
! Create/open cdf file
! 03/09/99 C.Ludescher
!
      include "netcdf.inc"
      INTEGER,       intent(out) :: ncid
      character*(*), intent(in) :: filename
      character*1,   intent(in) :: opt
      integer, optional, intent(out) :: ier
      integer :: status
 
      if (opt .eq. 'w' .or. opt .eq. 'W') then
! New file...
         status = nf_create(filename,nf_clobber,ncid)
         call handle_err(status,filename,'cdfcrt','nf_create')
      else if (opt .eq. 'm' .or. opt .eq. 'M') then
! Open for read/write
         status = nf_open(filename, nf_write, ncid)
         call handle_err(status,'filename','cdfopn','nf_open')
      else
! Open for readonly
         status = nf_open(filename, nf_nowrite, ncid)
         call handle_err(status,'filename','cdfopn','nf_open')
      end if
      if (PRESENT (ier)) then
         if (status .ne. NF_NOERR) then
            ier = 1
         else
            ier = 0
         endif
      endif
      return
      end

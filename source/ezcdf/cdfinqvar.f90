subroutine cdfInqVar(ncid,varnam,dimlens,eztype,ier)
  ! Inquire a Variable and its dimensions
  ! 03/08/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex type (ap) Wed May 16 15:18:05 EDT 2001
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in)          :: ncid
  character*(*), intent(in)          :: varnam
  ! Output
  integer, dimension(:), intent(out) :: dimlens
  character*4,           intent(out) :: eztype
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: ndims, varid, natts, xtype
  integer                 :: status, i
  integer, dimension(3)   :: dimids
  character*(nf_max_name) :: name, varnam_long
  integer, parameter      :: cmplx_len = 13
  character(cmplx_len), parameter :: cmplx_name = '__CmPlx_Re_Im'
  logical :: is_complex
 
  if (PRESENT (ier)) ier = 1
  is_complex = .false.
 
  status = nf_inq_varid(ncid,varnam,varid)
  if (status .ne. 0) then
     ! perhaps varnam is complex, try...
     status = nf_inq_varid(ncid,varnam//cmplx_name,varid)
     if(status .eq. 0) is_complex = .true.
  endif
  varnam_long = varnam
  if(is_complex)  varnam_long = varnam//cmplx_name
  call handle_err(status,varnam_long,'cdfInqVar','nf_inq_varid')
  if(status .ne. 0) return
 
  status = nf_inq_var(ncid,varid,name,xtype,ndims,dimids,natts)
  call handle_err(status,varnam_long,'cdfInqVar','nf_inq_var')
  if (status .ne. 0) return
 
  if (size(dimlens) .lt. ndims) return
  dimlens = 0
  select case (xtype)
  case (nf_double)
     eztype = 'R8'
     if(is_complex) eztype = 'C16'
  case (nf_int)
     eztype = 'INT'
  case (nf_float)
     eztype = 'R4'
     if(is_complex) eztype = 'C8'
  case (nf_char)
     eztype = 'CHAR'
  end select
  do i=1,ndims
     status = nf_inq_dim(ncid,dimids(i),name,dimlens(i))
     call handle_err(status,varnam_long,'cdfInqVar','nf_inq_dim')
  end do
 
  if(is_complex) then
     dimlens(1) = dimlens(1)/2
  endif
  if (PRESENT (ier)) ier = status
 
end subroutine cdfInqVar

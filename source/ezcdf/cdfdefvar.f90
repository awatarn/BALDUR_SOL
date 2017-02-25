subroutine cdfDefVar(ncid,varnam,dimlens,type,ier)
  ! Define a Variable and its dimensions
  ! 03/08/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! added support for complex types (ap) Wed May 16 15:06:34 EDT 2001
 
  implicit none
  include "netcdf.inc"
  ! Input
  integer,              intent(in)  :: ncid
  character*(*),        intent(in)  :: varnam
  integer, dimension(:),intent(in)  :: dimlens
  character*(*),        intent(in)  :: type
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer                 :: ndims, varid, xtype, maxdims
  integer                 :: i,     n,     status
  integer, dimension(3)   :: dimids
  character*(nf_max_name) :: dimnams(3)
  character*(5)           :: zdim
  character*(9)           :: zdim1
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  integer dims(3)
  integer flag
 
  if (PRESENT (ier)) ier = 1
  status = 0
 
  dims = dimlens
  if (type .eq. 'R8' .or. type .eq. 'r8') then
     xtype=nf_double
  else if (type .eq. 'C16' .or. type .eq. 'c16') then
     ! complex array = real array with double leading length
     xtype=nf_double
     dims(1) = max(2*dimlens(1), 2)
  else if (type .eq. 'INT' .or. type .eq. 'int') then
     xtype=nf_int
  else if (type .eq. 'R4' .or. type .eq. 'r4') then
     xtype=nf_float
  else if (type .eq. 'C8' .or. type .eq. 'c8') then
     ! complex array = real array with double leading length
     xtype=nf_float
     dims(1) = max(2*dimlens(1), 2)
  else if (type .eq. 'CHAR' .or. type .eq. 'char') then
     xtype=nf_char
  end if
 
  ! The normal behavior is to decrement the rank if all
  ! sizes to the right are 1. However, in some cases it
  ! may be desirable to keep this dimension, in which the
  ! user needs to pass -1 as dimension.
 
  ndims = size(dims)
  n=ndims
  flag = 1
  do i=n,1,-1
     if ((dims(i) == 1 .or. dims(i) == 0) .and. flag==1 ) then
        ndims = i-1
     else
        flag=0
     endif
     dims(i) = abs(dims(i))
  end do
 
  ! Check ndims <= maxdims
  if (type .eq. 'CHAR' .or. type .eq. 'char') then
     maxdims = 2
  else
     maxdims = 3
  endif
  if (ndims .gt. maxdims) then
     WRITE(*,10)ndims,type
10   format('% cdfDefVar --E--  Rank',i3,' not supported for ',a)
     return
  endif
 
  do i=1,ndims
     if (dims(i).gt. 999999999) then
        print *,'% cdfDefVar --E-- dimension >= 1.e9 not supported'
        return
     else if (dims(i) .lt. 100000) then
        write(zdim,'(i5.5)') dims(i)
        dimnams(i) = 'dim_'//zdim
     else
        write(zdim1,'(i9.9)') dims(i)
        dimnams(i) = 'dim_'//zdim1
     endif
     status = nf_inq_dimid(ncid,dimnams(i),dimids(i))
     if (status .ne. NF_NOERR) then
        status = nf_def_dim(ncid,dimnams(i),dims(i),dimids(i))
        call handle_err(status,dimnams(i),'cdfdef','nf_def_dim')
     endif
  end do
 
  if (status .ne. 0) return
 
 
  if(type/='C16' .and. type/='c16' .and. &
       & type/='C8' .and. type/='c8') then
 
     ! general case
 
     status = nf_def_var(ncid,varnam,xtype,ndims,dimids,varid)
 
  else
 
     ! do some name mangling to keep track of complex types, since
     ! and store them as reals.
 
     status = nf_def_var(ncid,varnam//cmplx_name, &
          xtype,ndims,dimids,varid)
  endif
 
  call handle_err(status,'varnam','cdfDefVar','nf_def_var')
  if (PRESENT (ier)) ier = status
end subroutine cdfDefVar

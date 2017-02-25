! automatic conversion to free f90 compatible form
! free.pl cdfinqv.for
! linewidth: 72
! file names: cdfinqv.for
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
      subroutine cdfInqV(ncid,varnam,varid,dimlens,ndims,status)
! Inquire variable-id and dimlens
! 03/09/99 C.Ludescher
!
      implicit none
!
      include "netcdf.inc"
! Input
      integer,       intent(in)  :: ncid
      character*(*), intent(in)  :: varnam
! Returns
!      integer, dimension(3)  ::  dimlens
!      integer ::  ndims, varid
!      integer ::  status
      integer, dimension(:), intent(out) ::  dimlens
      integer,               intent(out) ::  ndims, varid
      integer,               intent(out) ::  status
! Local
      integer                 :: natts, xtype, i
      integer, dimension(3)   :: dimids
      character*(nf_max_name) :: name
!---------------------------------------------------------------------------
 
      status = nf_inq_varid(ncid,varnam,varid)
      call handle_err(status,varnam,'cdfInqV','nf_inq_varid')
      if (status .ne. 0) return
      status = nf_inq_var(ncid,varid,name,xtype,ndims,dimids,natts)
      call handle_err(status,varnam,'cdfInqV','nf_inq_var')
      do i=1,ndims
         status = nf_inq_dimlen(ncid,dimids(i),dimlens(i))
         call handle_err(status,varnam,'cdfInqV','nf_inq_dimlen')
      end do
      END SUBROUTINE cdfInqV

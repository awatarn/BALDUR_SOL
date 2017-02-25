! automatic conversion to free f90 compatible form
! free.pl cdfinqint.for
! linewidth: 72
! file names: cdfinqint.for
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
      MODULE ezcdf_InqInt
! 03/09/99 C.Ludescher
!
      INTERFACE cdfInqVarDim
      SUBROUTINE cdfInqV(ncid,varnam,varid,dimlens,ndims,status)
      implicit none
!
! Input
      integer,       intent(in)  :: ncid
      character*(*), intent(in)  :: varnam
! Returns
      integer, dimension(:), intent(out) ::  dimlens
      integer,               intent(out) ::  ndims, varid
      integer,               intent(out) ::  status
      END SUBROUTINE cdfInqV
      SUBROUTINE cdfgv(ncid,varnam,varid,dimlens,sizes,type,status)
      implicit none
! Input
      integer                 :: ncid
      character*(*)           :: varnam
      character*1             :: type
      integer, dimension(:)   :: sizes
! Output
      integer, dimension(:)   :: dimlens
      integer                 :: varid, status
      END SUBROUTINE cdfgv
      END INTERFACE
      END MODULE ezcdf_InqInt
 
 

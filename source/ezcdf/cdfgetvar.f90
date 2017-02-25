! Read netcdf data variable
! 02/15/99 C. Ludescher
! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
 
SUBROUTINE cdfr_3i(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  ! Read 2 dimensional Integer array
  !
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:,:,:), intent(out) :: varval
  integer, optional,         intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  cnt(3) = 1
  do k = 1,min(dimlens(2),ldim(2))
     st(2) = k
     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab in varval
        st(3) = j           ! Start of slab
        status = nf_get_vara_int(ncid,varid,st,cnt,varval(1,k,j))
        if (status .ne. NF_NOERR) then
           call handle_err(status,varnam,'cdfr_3i',                 &
                &              'nf_get_vara_int')
           return
        end if
     end do
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3i
SUBROUTINE cdfr_3d(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:,:,:), intent(out) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  cnt(3) = 1
  do k = 1,min(dimlens(2),ldim(2))
     st(2) = k
     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
        st(3) = j                      ! Start of slab
        status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,k,j))
        if (status .ne. NF_NOERR) then
           call handle_err(status,varnam,'cdfr_3d',                    &
                &           'nf_get_vara_double')
           return
        end if
     end do
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3d
SUBROUTINE cdfr_3c16(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:,:,:), intent(out) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  cnt(3) = 1
  do k = 1,min(dimlens(2),ldim(2))
     st(2) = k
     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
        st(3) = j                      ! Start of slab
        status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,k,j))
        if (status .ne. NF_NOERR) then
           call handle_err(status,varnam,'cdfr_3c16',                    &
                &           'nf_get_vara_double')
           return
        end if
     end do
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3c16
SUBROUTINE cdfr_3f(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:,:,:), intent(out) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  cnt(3) = 1
  do k = 1,min(dimlens(2),ldim(2))
     st(2) = k
     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
        st(3) = j                      ! Start of slab
        status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,k,j))
        if (status .ne. NF_NOERR) then
           call handle_err(status,varnam,'cdfr_3f',                    &
                &           'nf_get_vara_real')
           return
        end if
     end do
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3f
SUBROUTINE cdfr_3c8(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:,:,:), intent(out) ::  varval
  integer, optional,               intent(out) :: ier
  ! Local
  integer, dimension(3)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j,  k
  integer, dimension(3)   :: dimlens
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  cnt(3) = 1
  do k = 1,min(dimlens(2),ldim(2))
     st(2) = k
     do j = 1,min(dimlens(3),ldim(3)) ! For each Z : read slab into varval
        st(3) = j                      ! Start of slab
        status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,k,j))
        if (status .ne. NF_NOERR) then
           call handle_err(status,varnam,'cdfr_3c8',                    &
                &           'nf_get_vara_real')
           return
        end if
     end do
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_3c8
 
SUBROUTINE cdfr_2i(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  ! Read 2 dimensional Integer array
  !
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:,:), intent(out) :: varval
  integer, optional,       intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_int(ncid,varid,st,cnt,varval(1,j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2i',                    &
             &           'nf_get_vara_int')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2i
SUBROUTINE cdfr_2d(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:,:), intent(out) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2d',                    &
             &           'nf_get_vara_double')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2d
SUBROUTINE cdfr_2c16(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:,:), intent(out) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_double(ncid,varid,st,cnt,varval(1,j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c16',                    &
             &           'nf_get_vara_double')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c16
SUBROUTINE cdfr_2f(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:,:), intent(out) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2f',                    &
             &           'nf_get_vara_real')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2f
SUBROUTINE cdfr_2c8(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:,:), intent(out) ::  varval
  integer, optional,             intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, j
  integer, dimension(2)   :: dimlens
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1) ! Re/Pairs
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_real(ncid,varid,st,cnt,varval(1,j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c8',                    &
             &           'nf_get_vara_real')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c8
SUBROUTINE cdfr_2c(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  character*(*), dimension(:),intent(out) :: varval
  integer, optional,          intent(out) :: ier
  ! Local
  integer, dimension(2)   :: st,     cnt,    ldim
  integer                 :: varid,  status, charlen, j
  integer, dimension(2)   :: dimlens
  if (PRESENT (ier)) ier = 1
  ldim(1) = len(varval)
  ldim(2)=size(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'c',status)
  if (status .ne. 0) return
  !
  st(1) = 1
  cnt(1) = min(dimlens(1),ldim(1)) ! x count
  cnt(2) = 1
  do j = 1,min(dimlens(2),ldim(2)) ! For each Y : read slab into varval
     st(2) = j              ! Start of slab
     status = nf_get_vara_text(ncid,varid,st,cnt,varval(j))
     if (status .ne. NF_NOERR) then
        call handle_err(status,varnam,'cdfr_2c',                       &
             &        'nf_get_var_text')
        return
     end if
  end do
  if (PRESENT (ier)) ier = 0
END SUBROUTINE cdfr_2c
 
SUBROUTINE cdfr_1i(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer, dimension(:), intent(out) :: varval
  integer, optional,     intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_int(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1i','nf_get_var_int')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1i
SUBROUTINE cdfr_1d(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8), dimension(:), intent(out) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1d
SUBROUTINE cdfr_1c16(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r8), dimension(:), intent(out) :: varval
  integer, optional,           intent(out) :: ier
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1) ! Re/Im pairs
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c16','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c16
SUBROUTINE cdfr_1f(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4), dimension(:), intent(out) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1f
SUBROUTINE cdfr_1c8(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4), dimension(:), intent(out) :: varval
  integer, optional,           intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim=ubound(varval)
  ldim(1) = 2*ldim(1)
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c8','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c8
SUBROUTINE cdfr_1c(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in)  :: ncid
  character*(*), intent(in)  :: varnam
  ! Output
  character*(*),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = len(varval)
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'c',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_text(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_1c','nf_get_var_text')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_1c
 
SUBROUTINE cdfr_0i(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  integer,           intent(out) ::  varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'i',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_int(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0i','nf_get_var_int')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0i
SUBROUTINE cdfr_0d(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0d
SUBROUTINE cdfr_0c16(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r8),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim(1) = 2 ! Re/Im pair
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'d',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_double(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0d','nf_get_var_double')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0c16
SUBROUTINE cdfr_0f(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  REAL(KIND=r4),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  if (PRESENT (ier)) ier = 1
  ldim(1) = 0
  call cdfgv(ncid,varnam,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0f
SUBROUTINE cdfr_0c8(ncid,varnam,varval,ier)
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  ! Output
  COMPLEX(KIND=r4),     intent(out) :: varval
  integer, optional, intent(out) :: ier
  ! Local
  integer                 :: varid,  status
  integer, dimension(1)   :: dimlens,ldim
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  ldim(1) = 2 ! Re/Im pair
  call cdfgv(ncid,varnam//cmplx_name,varid,dimlens,ldim,'r',status)
  if (status .ne. 0) return
  !
  status = nf_get_var_real(ncid,varid,varval)
  call handle_err(status,varnam,'cdfr_0f','nf_get_var_real')
  if (PRESENT (ier)) ier = status
END SUBROUTINE cdfr_0c8

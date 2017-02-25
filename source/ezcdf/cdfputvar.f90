!
! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
!
! added support (ap) Wed May 16 15:06:34 EDT 2001
SUBROUTINE cdfw_3i(ncid,varnam,varval,ier)
  !     Write 3 dimensional Integer data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,                   intent(in) :: ncid
  character*(*),             intent(in) :: varnam
  integer, dimension(:,:,:), intent(in) :: varval
  ! Output
  integer, optional,        intent(out) :: ier
  ! Local
  integer, dimension(3) :: st = (/1,1,1/)
  integer, dimension(3) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 3) then
     print "('% cdfPutVar_3i: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_int(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdfPutVar_3i','nf_put_var_int')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_3i
 
SUBROUTINE cdfw_3d(ncid,varnam,varval,ier)
  !     Write 3 dimensional 64-bit Real data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                         intent(in) :: ncid
  character*(*),                   intent(in) :: varnam
  REAL(KIND=r8), dimension(:,:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(3) :: st = (/1,1,1/)
  integer, dimension(3) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 3) then
     print "('% cdfPutVar_3d: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_3d','nf_put_vara_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_3d
 
SUBROUTINE cdfw_3c16(ncid,varnam,varval,ier)
  !     Write 3 dimensional 128-bit complex data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                         intent(in) :: ncid
  character*(*),                   intent(in) :: varnam
  COMPLEX(KIND=r8), dimension(:,:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(3) :: st = (/1,1,1/)
  integer, dimension(3) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 3) then
     print "('% cdfPutVar_3c16: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_3c16','nf_put_vara_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_3c16
 
SUBROUTINE cdfw_3f(ncid,varnam,varval,ier)
  !     Write 3 dimensional default Real data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                         intent(in) :: ncid
  character*(*),                   intent(in) :: varnam
  REAL(KIND=r4), dimension(:,:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(3) :: st = (/1,1,1/)
  integer, dimension(3) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 3) then
     print "('% cdfPutVar_3f: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_3f','nf_put_vara_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_3f
 
SUBROUTINE cdfw_3c8(ncid,varnam,varval,ier)
  !     Write 3 dimensional complex*8 data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                         intent(in) :: ncid
  character*(*),                   intent(in) :: varnam
  COMPLEX(KIND=r4), dimension(:,:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(3) :: st = (/1,1,1/)
  integer, dimension(3) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 3) then
     print "('% cdfPutVar_3c8: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_3c8','nf_put_vara_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_3c8
 
SUBROUTINE cdfw_2i(ncid,varnam,varval,ier)
  !     Write 2 dimensional Integer data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,                 intent(in) :: ncid
  character*(*),           intent(in) :: varnam
  integer, dimension(:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2i: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_int(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdfPutVar_2i','nf_put_var_int')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2i
 
SUBROUTINE cdfw_2d(ncid,varnam,varval,ier)
  !     Write 2 dimensional 64-bit Real data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                       intent(in) :: ncid
  character*(*),                 intent(in) :: varnam
  REAL(KIND=r8), dimension(:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2d: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_2d','nf_put_vara_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2d
 
SUBROUTINE cdfw_2c16(ncid,varnam,varval,ier)
  !     Write 2 dimensional 128-bit complex data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                       intent(in) :: ncid
  character*(*),                 intent(in) :: varnam
  COMPLEX(KIND=r8), dimension(:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2c16: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_double(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_2c16','nf_put_vara_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2c16
 
SUBROUTINE cdfw_2f(ncid,varnam,varval,ier)
  !     Write 2 dimensional default Real data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                       intent(in) :: ncid
  character*(*),                 intent(in) :: varnam
  REAL(KIND=r4), dimension(:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2f: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_2f','nf_put_vara_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2f
 
SUBROUTINE cdfw_2c8(ncid,varnam,varval,ier)
  !     Write 2 dimensional complex*8 data array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                       intent(in) :: ncid
  character*(*),                 intent(in) :: varnam
  COMPLEX(KIND=r4), dimension(:,:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts= nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0)  return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2c8: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_real(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdf_2c8','nf_put_vara_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2c8
 
SUBROUTINE cdfw_2c(ncid,varnam,varval,ier)
  !     Write 2 dimensional character array
  !     Use cdfDefVar to define the Variable
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,                     intent(in) :: ncid
  character*(*),               intent(in) :: varnam
  character*(*), dimension(*), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(2) :: st = (/1,1/)
  integer, dimension(2) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 2) then
     print "('% cdfPutVar_2c: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_vara_text(ncid,varid,st,dimlens,varval)
  call handle_err(sts,varnam,'cdfPutVar_2c','nf_put_var_text')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_2c
 
SUBROUTINE cdfw_1i(ncid,varnam,varval,ier)
  !     write 1 dimensional Integer array
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,                 intent(in) :: ncid
  character*(*),           intent(in) :: varnam
  integer, dimension(:),   intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1i: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_int(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_1i','nf_put_var_int')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1i
 
SUBROUTINE cdfw_1d(ncid,varnam,varval,ier)
  !     Write 1 dimensional 64-bit data array
  !     Use cdfDefVar to define the Variable
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                     intent(in) :: ncid
  character*(*),               intent(in) :: varnam
  REAL(KIND=r8), dimension(:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1d: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_double(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPurVar_1d','nf_put_var_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1d
 
SUBROUTINE cdfw_1c16(ncid,varnam,varval,ier)
  !     Write 1 dimensional 128-bit complex data array
  !     Use cdfDefVar to define the Variable
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,                     intent(in) :: ncid
  character*(*),               intent(in) :: varnam
  COMPLEX(KIND=r8), dimension(:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1c16: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_double(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPurVar_1c16','nf_put_var_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1c16
 
SUBROUTINE cdfw_1f(ncid,varnam,varval,ier)
  !     Write 1 dimensional default real array
  !     Use cdfDefVar to define the Variable
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                     intent(in) :: ncid
  character*(*),               intent(in) :: varnam
  REAL(KIND=r4), dimension(:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1f: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_real(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPurVar_1f','nf_put_var_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1f
 
SUBROUTINE cdfw_1c8(ncid,varnam,varval,ier)
  !     Write 1 dimensional complex*8 array
  !     Use cdfDefVar to define the Variable
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,                     intent(in) :: ncid
  character*(*),               intent(in) :: varnam
  COMPLEX(KIND=r4), dimension(:), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1c8: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_real(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPurVar_1c8','nf_put_var_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1c8
 
SUBROUTINE cdfw_1c(ncid,varnam,varval,ier)
  !     Write 1 dimensional charcter data array
  !     Use cdfDefVar to define the Variable
  !
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  character*(*), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st = (/1/)
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_1c: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_text(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_1c','nf_put_var_text')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_1c
 
SUBROUTINE cdfw_0i(ncid,varnam,varval,ier)
  !     write integer scalar
  USE ezcdf_InqInt
  implicit none
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  integer,       intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 0) then
     print "('% cdfPutVar_0i: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_int(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_0i','nf_put_var_int')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_0i
 
SUBROUTINE cdfw_0d(ncid,varnam,varval,ier)
  !     Write Real*8 Scalar
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  REAL(KIND=r8), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 0) then
     print "('% cdfPutVar_0d: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_double(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_0d','nf_put_var_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_0d
 
SUBROUTINE cdfw_0c16(ncid,varnam,varval,ier)
  !     Write Complex*16 Scalar
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  COMPLEX(KIND=r8), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then ! scalar complex stored as 2 element real array
     print "('% cdfPutVar_0c16: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_double(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_0c16','nf_put_var_double')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_0c16
 
SUBROUTINE cdfw_0f(ncid,varnam,varval,ier)
  !     Write default real
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  REAL(KIND=r4), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 0) then
     print "('% cdfPutVar_0f: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_real(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_0f','nf_put_var_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_0f
 
SUBROUTINE cdfw_0c8(ncid,varnam,varval,ier)
  !     Write complex*8
  !
  USE ezcdf_InqInt
  implicit none
  integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
  include "netcdf.inc"
  ! Input
  integer,       intent(in) :: ncid
  character*(*), intent(in) :: varnam
  complex(KIND=r4), intent(in) :: varval
  ! Output
  integer, optional,      intent(out) :: ier
  ! Local
  integer, dimension(1) :: st
  integer, dimension(1) :: dimlens
  integer               :: varid, ndims, sts
  character*13, parameter :: cmplx_name = '__CmPlx_Re_Im'
 
  if (PRESENT (ier)) ier = 1
  sts = nf_enddef(ncid)
  call cdfInqV(ncid,varnam//cmplx_name,varid,dimlens,ndims,sts)
  if (sts .ne. 0) return
  if (ndims .ne. 1) then
     print "('% cdfPutVar_0c8: --E-- The variable ',a,               &
          &         ' was defined as',i2,' dimensional')",varnam,ndims
     return
  end if
  sts = nf_put_var_real(ncid,varid,varval)
  call handle_err(sts,varnam,'cdfPutVar_0c8','nf_put_var_real')
  if (PRESENT (ier)) ier = sts
END SUBROUTINE cdfw_0c8

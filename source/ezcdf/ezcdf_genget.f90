MODULE ezcdf_GenGet
  ! Generic Interface to Read netcdf data Variables
  ! 03/10/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex types (ap) Wed May 16 15:18:05 EDT 2001
  INTERFACE  cdfinquire
     SUBROUTINE cdfInqVar(ncid,varnam,dimlens,type,ier)
       implicit none
       ! Input
       integer,       intent(in)          :: ncid
       character*(*), intent(in)          :: varnam
       ! Output
       integer, dimension(:), intent(out) :: dimlens
       character*(*),         intent(out) :: type
       integer, optional,     intent(out) :: ier
     END SUBROUTINE cdfInqVar
  END INTERFACE
  !====================================================================
  ! Generic Read Routines: cdfGetVar
  INTERFACE cdfGetVar
 
     SUBROUTINE cdfr_3i(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       integer, dimension(:,:,:), intent(out) :: varval
       integer, optional,         intent(out) :: ier
     END SUBROUTINE cdfr_3i
     SUBROUTINE cdfr_3d(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r8), dimension(:,:,:), intent(out) ::  varval
       integer, optional,               intent(out) :: ier
     END  SUBROUTINE cdfr_3d
     SUBROUTINE cdfr_3c16(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r8), dimension(:,:,:), intent(out) ::  varval
       integer, optional,               intent(out) :: ier
     END  SUBROUTINE cdfr_3c16
     SUBROUTINE cdfr_3f(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r4), dimension(:,:,:), intent(out) ::  varval
       integer, optional,               intent(out) :: ier
     END  SUBROUTINE cdfr_3f
     SUBROUTINE cdfr_3c8(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r4), dimension(:,:,:), intent(out) ::  varval
       integer, optional,               intent(out) :: ier
     END  SUBROUTINE cdfr_3c8
 
     SUBROUTINE cdfr_2i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       integer, dimension(:,:), intent(out) :: varval
       integer, optional,       intent(out) :: ier
     END SUBROUTINE cdfr_2i
     SUBROUTINE cdfr_2d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r8), dimension(:,:), intent(out) ::  varval
       integer, optional,             intent(out) :: ier
     END SUBROUTINE cdfr_2d
     SUBROUTINE cdfr_2c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r8), dimension(:,:), intent(out) ::  varval
       integer, optional,             intent(out) :: ier
     END SUBROUTINE cdfr_2c16
     SUBROUTINE cdfr_2f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r4), dimension(:,:), intent(out) ::  varval
       integer, optional,             intent(out) :: ier
     END SUBROUTINE cdfr_2f
     SUBROUTINE cdfr_2c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r4), dimension(:,:), intent(out) ::  varval
       integer, optional,             intent(out) :: ier
     END SUBROUTINE cdfr_2c8
     SUBROUTINE cdfr_2c(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       character*(*), dimension(:),intent(out) :: varval
       integer, optional,          intent(out) :: ier
     END SUBROUTINE cdfr_2c
 
     SUBROUTINE cdfr_1i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       integer, dimension(:), intent(out) :: varval
       integer, optional,     intent(out) :: ier
     END SUBROUTINE cdfr_1i
     SUBROUTINE cdfr_1d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r8), dimension(:), intent(out) :: varval
       integer, optional,           intent(out) :: ier
     END SUBROUTINE cdfr_1d
     SUBROUTINE cdfr_1c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r8), dimension(:), intent(out) :: varval
       integer, optional,           intent(out) :: ier
     END SUBROUTINE cdfr_1c16
     SUBROUTINE cdfr_1f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r4), dimension(:), intent(out) :: varval
       integer, optional,           intent(out) :: ier
     END SUBROUTINE cdfr_1f
     SUBROUTINE cdfr_1c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r4), dimension(:), intent(out) :: varval
       integer, optional,           intent(out) :: ier
     END SUBROUTINE cdfr_1c8
     SUBROUTINE cdfr_1c(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in)  :: ncid
       character*(*), intent(in)  :: varnam
       ! Output
       character*(*),     intent(out) :: varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_1c
 
     SUBROUTINE cdfr_0i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       integer,           intent(out) ::  varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_0i
     SUBROUTINE cdfr_0d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       include "netcdf.inc"
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r8),     intent(out) :: varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_0d
     SUBROUTINE cdfr_0c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       include "netcdf.inc"
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r8),     intent(out) :: varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_0c16
     SUBROUTINE cdfr_0f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       include "netcdf.inc"
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       REAL(KIND=r4),     intent(out) :: varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_0f
     SUBROUTINE cdfr_0c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       include "netcdf.inc"
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       ! Output
       COMPLEX(KIND=r4),     intent(out) :: varval
       integer, optional, intent(out) :: ier
     END SUBROUTINE cdfr_0c8
  END INTERFACE
END MODULE ezcdf_GenGet

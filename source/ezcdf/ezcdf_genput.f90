MODULE ezcdf_GenPut
  ! Generic Interface to Write netcdf data Variables
  ! 03/10/99 C. Ludescher
  ! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
  ! + support for complex numbers (ap) Wed May 16 15:18:05 EDT 2001
  implicit none
  INTERFACE cdfdefine
     SUBROUTINE cdfDefVar(ncid,varnam,dimlens,type,ier)
       implicit none
       ! Input
       integer,              intent(in)  :: ncid
       character*(*),        intent(in)  :: varnam
       integer, dimension(:),intent(in)  :: dimlens
       character*(*),        intent(in)  :: type
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfDefVar
  END INTERFACE
 
  INTERFACE cdfPutVar
 
     SUBROUTINE cdfw_3i(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       ! Input
       integer,                   intent(in) :: ncid
       character*(*),             intent(in) :: varnam
       integer, dimension(:,:,:), intent(in) :: varval
       ! Output
       integer, optional,        intent(out) :: ier
     END SUBROUTINE cdfw_3i
     SUBROUTINE cdfw_3d(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                         intent(in) :: ncid
       character*(*),                   intent(in) :: varnam
       REAL(KIND=r8), dimension(:,:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_3d
     SUBROUTINE cdfw_3c16(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                         intent(in) :: ncid
       character*(*),                   intent(in) :: varnam
       COMPLEX(KIND=r8), dimension(:,:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_3c16
     SUBROUTINE cdfw_3f(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                         intent(in) :: ncid
       character*(*),                   intent(in) :: varnam
       REAL(KIND=r4), dimension(:,:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_3f
     SUBROUTINE cdfw_3c8(ncid,varnam,varval,ier)
       USE ezcdf_InqInt
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                         intent(in) :: ncid
       character*(*),                   intent(in) :: varnam
       COMPLEX(KIND=r4), dimension(:,:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_3c8
 
     SUBROUTINE cdfw_2i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,                 intent(in) :: ncid
       character*(*),           intent(in) :: varnam
       integer, dimension(:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2i
     SUBROUTINE cdfw_2d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                       intent(in) :: ncid
       character*(*),                 intent(in) :: varnam
       REAL(KIND=r8), dimension(:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2d
     SUBROUTINE cdfw_2c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                       intent(in) :: ncid
       character*(*),                 intent(in) :: varnam
       COMPLEX(KIND=r8), dimension(:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2c16
     SUBROUTINE cdfw_2f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                       intent(in) :: ncid
       character*(*),                 intent(in) :: varnam
       REAL(KIND=r4), dimension(:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2f
     SUBROUTINE cdfw_2c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                       intent(in) :: ncid
       character*(*),                 intent(in) :: varnam
       COMPLEX(KIND=r4), dimension(:,:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2c8
     SUBROUTINE cdfw_2c(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,                     intent(in) :: ncid
       character*(*),               intent(in) :: varnam
       character*(*), dimension(*), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_2c
 
     SUBROUTINE cdfw_1i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,                 intent(in) :: ncid
       character*(*),           intent(in) :: varnam
       integer, dimension(:),   intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1i
     SUBROUTINE cdfw_1d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                     intent(in) :: ncid
       character*(*),               intent(in) :: varnam
       REAL(KIND=r8), dimension(:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1d
     SUBROUTINE cdfw_1c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,                     intent(in) :: ncid
       character*(*),               intent(in) :: varnam
       COMPLEX(KIND=r8), dimension(:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1c16
     SUBROUTINE cdfw_1f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                     intent(in) :: ncid
       character*(*),               intent(in) :: varnam
       REAL(KIND=r4), dimension(:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1f
     SUBROUTINE cdfw_1c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,                     intent(in) :: ncid
       character*(*),               intent(in) :: varnam
       COMPLEX(KIND=r4), dimension(:), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1c8
     SUBROUTINE cdfw_1c(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       character*(*), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_1c
 
     SUBROUTINE cdfw_0i(ncid,varnam,varval,ier)
       implicit none
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       integer,       intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_0i
     SUBROUTINE cdfw_0d(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       REAL(KIND=r8), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_0d
     SUBROUTINE cdfw_0c16(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       COMPLEX(KIND=r8), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_0c16
     SUBROUTINE cdfw_0f(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       REAL(KIND=r4), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_0f
     SUBROUTINE cdfw_0c8(ncid,varnam,varval,ier)
       implicit none
       integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
       ! Input
       integer,       intent(in) :: ncid
       character*(*), intent(in) :: varnam
       COMPLEX(KIND=r4), intent(in) :: varval
       ! Output
       integer, optional,      intent(out) :: ier
     END SUBROUTINE cdfw_0c8
  END INTERFACE
END MODULE ezcdf_GenPut

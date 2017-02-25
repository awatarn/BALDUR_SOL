! Interface for cdfopn to handle optional argument ier
! 04/28/00 C.Ludescher
! + ezcdf_close for symmetry (ap)
MODULE ezcdf
  USE ezcdf_GenPut
  USE ezcdf_GenGet
  INTERFACE ezcdf_open
     SUBROUTINE ezcdf_open(ncid,filename,opt,ier)
       include "netcdf.inc"
       INTEGER,       intent(out) :: ncid
       character*(*), intent(in) :: filename
       character*1,   intent(in) :: opt
       integer, optional,         intent(out) :: ier
     END SUBROUTINE ezcdf_open
  END INTERFACE
  INTERFACE ezcdf_close
     SUBROUTINE ezcdf_close(ncid,ier)
       include "netcdf.inc"
       INTEGER,       intent(in) :: ncid
       integer, optional,         intent(out) :: ier
     END SUBROUTINE ezcdf_close
  END INTERFACE
END MODULE ezcdf

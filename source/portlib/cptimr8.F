#include "fpreproc/fp_const.h"
C
      subroutine cptimr8(ztime)
C
      IMPLICIT NONE
C
      REAL*8 ztime
      REAL zctime
C
C  this routine returns the elapsed CPU time.  Valid for UNIX or VMS
C  value returned is REAL*8.
C
      call cptimi(zctime)
      ztime=CONST(0.01,0)*zctime
C
      return
      end
 

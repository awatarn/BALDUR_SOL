      subroutine c9date(zdate)
      character*(*) zdate
C
C  get a 9 character date -- portably -- by calling C code
C
#include "fpreproc/byte_declare.h"
C
      BYTE_DECLARE z9date(9)
C
      ilen=9
      call c9dat0(z9date,ilen)
C
      call cstring(zdate,z9date,'2F')
C
      return
      end

c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bdhden (
     &       lbound,     cbound,     znebar,     zcurrent,
     &       zbtor,      zneped)
c
      IMPLICIT NONE
c 
      integer    lbound(32), nout
c
      real       cbound(32),     znebar,     zneped,      zcurrent,
     &           zbtor,          Cn
c
c..calculate the pedestal density for simple empirical model
c  constructed from fitting 533 Type I ELMy H-mode plasma
c
      if ((cbound(6) .gt. 0.01) .and. (cbound(6) .lt. 2.00)) then
      	 Cn = cbound(6)
      else
        Cn = 0.71
      endif      
c
      if (lbound(3) .le. 1) then
c
        zneped = Cn * znebar
c
        write(6,*)'nped model 1'
c
c..calculate the pedestal density for simple empirical model
c  constructed from line average density (flat profile)
c

      elseif (lbound(3) .eq. 2) then
c
        zneped = 0.74 * ((znebar*1.e-20)**0.99) * 
     &           (zcurrent**0.15) * (zbtor**(-0.12))
c
        write(6,*)'nped model 2'
c  
       endif 
c
      return
      end

c@chisqr.f  22:00 27-Feb-96  Glenn Bateman
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    chisqr.f computes the normalized chi square statistic
c
      subroutine chisqr ( kdata, preldev, prelerr
     &  , epslon, pvarian, pchisqr, kerr )
c
c  kdata    = input number of data points
c  preldev  = input array of normalized deviations
c  prelerr  = input array of normalized error bar values
c  epslon   = input sqrt ( machine epsilon )
c  pvarian  = input value of normalized variance
c  pchisqr  = output value of normalized chi squared value
c  kerr     = 0 for normal output
c           = 1 if not converged
c
c    Compute:
c
c  pchisqr = 
c    sum ( preldev(j) - poffset )**2 / ( pvarian + prelerr(j)**2 )
c    / max ( kdata - 1, 1 )
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      integer kdata, kerr, j
c
      real preldev(*), prelerr(*), epslon, pvarian, pchisqr
     &  , zoffset, ztemp0, ztemp1, ztemp2, ztemp3
c
c
c
c..evaluate equations
c
      ztemp1 = 0.0
      ztemp2 = 0.0
      ztemp3 = 0.0
      do j=1,kdata
        ztemp0 = 1.0 / ( max(pvarian, epslon) + prelerr(j)**2 )
        ztemp1 = ztemp1 + ztemp0
        ztemp2 = ztemp2 + preldev(j) * ztemp0
      enddo
c
      zoffset = ztemp2 / ztemp1
c
      do j=1,kdata
        ztemp0 = 1.0 / ( max(pvarian, epslon) + prelerr(j)**2 )
        ztemp3 = ztemp3 + ( preldev(j) - zoffset )**2 * ztemp0
      enddo
c
      pchisqr = ztemp3 / max ( kdata - 1, 1 )
c
c
      return
      end


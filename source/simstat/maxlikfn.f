c@maxlikfn.f  19:00 27-Feb-96  Glenn Bateman
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine maxlikfn ( kdata, preldev, prelerr
     &  , epslon, pvarian, poffset, pmlfunc, kerr )
c
c  kdata    = number of data points
c  preldev  = array of deviations normalized by pmax
c  prelerr  = array of error bar values normalized by pmax
c  epslon   = sqrt ( machine epsilon )
c  pvarian  = output value of normalized variance
c  poffset  = output value of offset normalized by pmax
c  pmlfunc  = output value of maximum likelihood function (see below)
c  kerr     = 0 for normal output
c           = 1 if not converged
c
c    Evaluate:
c
c  sum [ pvarian + prelerr(j)**2 - ( preldev(j) - poffset )**2 )
c     / ( pvarian + prelerr(j)**2 )**2 ]
c
c  using
c
c  poffset =   sum [ preldev(j) / ( pvarian + prelerr(j)**2 ) ]
c            / sum [ 1.0 / ( pvarian + prelerr(j)**2 ) ]
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      integer kdata, kerr, j, ji
c
      real preldev(*), prelerr(*), epslon, poffset, pvarian, pmlfunc
     &  , ztemp0, ztemp1, ztemp2, ztemp3
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
      poffset = ztemp2 / ztemp1
c
      do j=1,kdata
        ztemp0 = 1.0 / ( max(pvarian, epslon) + prelerr(j)**2 )
        ztemp1 = ztemp1 + ztemp0
        ztemp3 = ztemp3 + ( preldev(j) - poffset )**2 * ztemp0**2
      enddo
c
      pmlfunc = ztemp1 - ztemp3
c
c
      return
      end

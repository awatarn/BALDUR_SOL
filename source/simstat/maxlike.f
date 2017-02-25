c@maxlike.f  19:00 27-Feb-96  Glenn Bateman
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    maxlike.f computes the maximum likelihood offset and RMS deviations
c
      subroutine maxlike ( kdata, preldev, prelerr
     &  , epslon, poffset, prmsdev, pchisqr, kerr )
c
c  kdata    = input number of data points
c  preldev  = input array of normalized deviations
c  prelerr  = input array of normalized error bar values
c  epslon   = input sqrt ( machine epsilon )
c  poffset  = output value of normalized offset
c  prmsdev  = output vlaue of normalized RMS deviation
c  pchisqr  = normalized chi squared value
c  kerr     = 0 for normal output
c           = 1 if not converged
c
c    Solve the equations
c
c  sum ( preldev(j) - poffset ) / ( prmsdev**2 + prelerr(j)**2 ) = 0
c
c  sum [ prmsdev**2 + prelerr(j)**2 - ( preldev(j) - poffset )**2 )
c     / ( prmsdev**2 + prelerr(j)**2 )**2 ]  = 0
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      integer kdata, kerr, j, ji
c
      real preldev(*), prelerr(*), epslon, poffset, prmsdev, pchisqr
     &  , zoffnew, zoffold, zvarnew, zvarold, zmlfnew, zerrmax
     &  , zvarlow, zofflow, zmlflow, zvarhigh, zoffhigh, zmlfhigh
     &  , ztemp0, ztemp1, ztemp2, ztemp3
c
c  zoffnew, zoffold = offset values
c  zvarnew, zvarold = variance values
c
c
c..check input
c
      if ( kdata .lt. 1 ) then
        kerr = 2
        write (*,*) ' kdata .lt. 1 in sbrtn maxlike'
        poffset = 0.0
        prmsdev = 0.0
        pchisqr = 0.0
        return
      endif
c
c..initial values
c
      poffset = 0.0
      prmsdev = 0.0
c
      zerrmax = 0.0
      zoffnew = 0.0
      do j=1,kdata
        zerrmax = max ( zerrmax, prelerr(j) )
        zoffnew = zoffnew + preldev(j)
      enddo
      zoffnew = zoffnew / real(kdata)
c
      zvarnew = 0.0
      do j=1,kdata
        zvarnew = zvarnew + ( preldev(j) - zoffnew )**2
      enddo
      zvarnew = zvarnew / real(kdata)
c
c..intermediate printout
c
cbate      write (*,*)
cbate      write (*,100)
cbate  100 format ('  ji',t8,'zoffnew',t23,'zvarnew'
cbate     &  ,t38,'zmlfunc')
c
      ji = 0
cbate      write (*,110) ji, zoffnew, zvarnew
c
c
c..if max prelerr(j) .lt. epslon, compute output and return
c
      if ( zerrmax .lt. epslon ) then
        poffset = zoffnew
        prmsdev = 0.0
        if ( zvarnew .gt. epslon ) prmsdev = sqrt ( zvarnew )
        call chisqr ( kdata, preldev, prelerr
     &  , epslon, zvarnew, pchisqr, kerr )
        return
      endif
c
c
c..solve nonlinear equation for variance using the bisection method
c  At sufficiently small values of the variance,
c  the functin should be negative.
c  At sufficiently large values of the variance,
c  the function should become positive.
c
      zvarlow = 0.0
      call maxlikfn ( kdata, preldev, prelerr
     &  , epslon, zvarlow, zofflow, zmlflow, kerr )
c
cbate        write (*,110) ji, zofflow, zvarlow, zmlflow
c
      if ( zmlflow .gt. 0.0 ) then
        poffset = zofflow
        prmsdev = 0.0
        call chisqr ( kdata, preldev, prelerr
     &  , epslon, zvarlow, pchisqr, kerr )
        return
      endif
c
c..search for the upper bound of variance
c
      do ji=1,5
        zvarhigh = zvarnew * real(ji)**2
        call maxlikfn ( kdata, preldev, prelerr
     &    , epslon, zvarhigh, zoffhigh, zmlfhigh, kerr )
c
cbate        write (*,110) ji, zoffhigh, zvarhigh, zmlfhigh
c
        if ( zmlfhigh .gt. 0.0 ) go to 20
        zvarlow = zvarhigh
        zofflow = zoffhigh
        zmlflow = zmlfhigh
      enddo
c
      kerr = 1
      write (*,*)
      write (*,*) ' maxlike did not converge '
      poffset = zoffhigh
      prmsdev = sqrt ( zvarhigh )
      call chisqr ( kdata, preldev, prelerr
     &  , epslon, zvarhigh, pchisqr, kerr )
      return
c
c..now do bisection
c
  20  continue
c
      do ji=1,20
c
        zvarnew = 0.5 * ( zvarlow + zvarhigh )
c
        call maxlikfn ( kdata, preldev, prelerr
     &    , epslon, zvarnew, zoffnew, zmlfnew, kerr )
c
        if ( zmlfnew .gt. 0.0 ) then
          zvarhigh = zvarnew
          zoffhigh = zoffnew
          zmlfhigh = zmlfnew
        else
          zvarlow  = zvarnew
          zofflow  = zoffnew
          zmlflow  = zmlfnew
        endif
c
c
cbate        write (*,110) ji, zoffnew, zvarnew, zmlfnew
cbate  110   format (i4,1p5e15.6)
c
      enddo
c
      kerr = 0
      poffset = zoffnew
      prmsdev = 0.0
      if ( zvarnew .gt. epslon ) prmsdev = sqrt ( zvarnew )
c
      call chisqr ( kdata, preldev, prelerr
     &  , epslon, zvarnew, pchisqr, kerr )
c
      return
      end


c@spntgl   .../baldur/code/util/spntgl.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      real function spntgl( n, ulower, upper, x, y, b, c, d )
c
      implicit none
c
      integer n
      real  ulower, upper, x(n), y(n), b(n), c(n), d(n)
c
c  this subroutine evaluates the integral of a cubic spline function
c
c  spntgl = integral_{ulower}^{upper} f(u) du
c
c  f(u) = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
      integer i, j, k, iend
      real  zsign, zlower, zupper, zlow, zupp, zint
      data i/1/
      save i
c
c..set lower and upper bound and sense of integration
c
      zsign  = 1.0
      if ( ulower .gt . upper ) then
        zsign = - 1.0
      else if ( ulower .eq. upper ) then
        spntgl = 0.0
        return
      endif
      zlower = min ( ulower, upper )
      zupper = max ( ulower, upper )
c
c  binary search
c
      if ( i .ge. n ) i = 1
      if ( zlower .lt. x(i)  .or.  zlower .gt. x(i+1) ) then
        i = 1
        j = n+1
   20   k = (i+j)/2
        if ( zlower .lt. x(k) ) j = k
        if ( zlower .ge. x(k) ) i = k
        if ( j .gt. i+1 ) go to 20
      endif
c
c  evaluate integral
c
      zint = 0.0
      zupp = zlower
c
  30  zlow = zupp
c
        zupp = zupper
        iend = 1
      if ( i .lt. n  .and.  zupper .gt. x(i+1) ) then
        zupp = x(i+1)
        iend = 0
      endif
c
      zint = zint + ( zupp - zlow ) * y(i)
     &   + ( (zupp-x(i))**2 - (zlow-x(i))**2 ) * b(i)/2.0
     &   + ( (zupp-x(i))**3 - (zlow-x(i))**3 ) * c(i)/3.0
     &   + ( (zupp-x(i))**4 - (zlow-x(i))**4 ) * d(i)/4.0 
c
      if ( iend .eq. 0 ) then
        i = i + 1
        go to 30
      endif
c
      spntgl = zsign * zint
c
      return
      end

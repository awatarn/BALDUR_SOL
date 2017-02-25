c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@sevalb   .../baldur/code/util/sevalb.f
c
c   Evaluate piecewise cubic polynomial or its derivative
c
      subroutine sevalb ( lderiv, n, x, y, b, c, d, m, u, s )
c
      integer lderiv, n, m
      real    x(n), y(n), b(n), c(n), d(n), u(*), s(*)
c
ccccccccccccccc
ccccccccccccccc
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    lderiv = 0 to evaluate polynomial
c           = 1 for first derivative
c           = 2 for second derivative
c    n = the number of data points
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c    m = number of points to be evaluated
c    u = array of abscissa at which the spline is to be evaluated
c    s = array of output values (output)
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
c
      integer i, ip, j, k
      real zdx
      data i/1/
      save i
c
      do 40 j=1,m
c
      if ( i .ge. n ) i = 1
      if ( u(j) .lt. x(i) ) go to 10
      if ( u(j) .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      ip = n+1
   20 k = (i+ip)/2
      if ( u(j) .lt. x(k) ) ip = k
      if ( u(j) .ge. x(k) ) i = k
      if ( ip .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 zdx = u(j) - x(i)
c
      if ( lderiv .lt. 1 ) then
        s(j) = y(i) + zdx*(b(i) + zdx*(c(i) + zdx*d(i)))
      else if ( lderiv .eq. 1 ) then
        s(j) = b(i) + zdx * ( 2.0 * c(i) + zdx * 3.0 * d(i) )
      else if ( lderiv .eq. 2 ) then
        s(j) = 2.0 * c(i) + zdx * 6.0 * d(i)
      else
        s(j) = 6.0 * d(i)
      endif
c
 40   continue
c
      return
      end

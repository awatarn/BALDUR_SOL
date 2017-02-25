c/ 15:00 06-oct-91 /11040/bald91/wbaldn1(LIB) wsutil(LIB) DINTLIB, Bateman
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c**********************************************************************c
c /11040/bald91/wbaldn1(LIB) wsutil(LIB) DINTLIB, Bateman
c
c    To obtain this file, type:
c cfs get /11040/bald91/wbaldn1 ^ end                   ( ^ = linefeed )
c lib wbaldn1 ^ x wsutil ^ n wsutil ^ x dintlib ^ end
c
c    Then, typing COSMOS DINTLIB compiles these subroutines and adds
c them to the library YUTILIB which can then be loaded together with
c the rest of the BALDUR code or other programs.
c
c Library of interpolation routines by Glenn Bateman
c
c             From Forsythe, Malcolm and Moler, "Computer Methods for
c                  Mathematical Computations", Prentice-Hall, 1977:
c  SPLINE ... to determine cubic spline coefficients
c  SEVAL  ... to evaluate cubic splines at array of points
c
c             Routines adapted from SPLINE and SEVAL by R. M. Wieland:
c  SPEVAL ... to evaluate derivative of cubic splines
c  SPLAAN ... to determine cubic spline coefficients when s'(x1)=0
c  SPLEEN ... to determine cubic spline coefficients when s'(xn)=0
c
c             Routines by Glenn Bateman:
c  QDINT1 ... quadratic interpolation given x(j), f(j), and f'(1).
c  QDEVAL ... Evalutation of a piecewise quadratic interpolation,
c                 given a close estimate of the interval
c                 in which to evaluate.
c  QDSOLV ... Solution of the quadratic equation in a given range.
c
c  CUBINT ... Perform and evaluate piecewise cubic interpolations
c                 using the IMSL routines IQHSCU, ICSEVU,DCSEVU,DCSQDU
c               moved to file cubint.f
c
c  TIMINT ... Piecewise linear interpolation with flat extrapolation
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@spline
c rgb 06-oct-91 replaced call tksplx with call abort when abscissa
c   is not in ascending order
c******************** start file spline.for ; group trkrlib ******************
cc
cccccccccccccccc
c       the codes (spline & seval) are taken from:
c       forsythe,malcolm and moler, "computer methods for
c       mathematical computations",prentice-hall, 1977.
c
c       the codes (spleen,splaan & speval) are adaptations
c       by r.m. wieland for special cases ... see comments


      subroutine spline (n, x, y, b, c, d)
      integer n
      real x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      real t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
c
c dmc - check ordinates
        do 5 i=1,nm1
          if ( x(i+1) .le. x(i) ) then
           call abortb (6,
     & 'abscissa not in ascending order in sbrtn spline file dintlib')
          endif
 5      continue
c
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end
c******************** end file spline.for ; group trkrlib ******************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@seval
c******************** start file seval.for ; group trkrlib ******************
c  mod dmc summer 1987 jet/garching
c   routine may evaluate linear instead of spline interpolation, if
c   common switch ilin is set ***
c....................................................
      real function seval(n, u, x, y, b, c, d)
      integer n
      real  u, x(n), y(n), b(n), c(n), d(n)
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
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c
c  dmc garching oct 1985 kluge
        common/zseval/ iput,dx,ilin
c  return zone index to caller via this common block ***
c
      integer i, j, k
      real dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
      if(ilin.eq.0) then 
        seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      else
        if(i.eq.n) then
          zslop=(y(n)-y(n-1))/(x(n)-x(n-1))
        else
          zslop=(y(i+1)-y(i))/(x(i+1)-x(i))
        endif
        seval=y(i)+dx*zslop
      endif
c
      iput=i
c
      return
      end
c******************** end file seval.for ; group trkrlib ******************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@speval
c******************** start file speval.for ; group trkrlib ******************
c....................................................
ccccccccccccccc
      real function speval(n, u, x, y, b, c, d)
      integer n
      real  u, x(n), y(n), b(n), c(n), d(n)
c
ccccccccccccccc
ccccccccccccccc
c  this subroutine evaluates the derivative of the cubic spline function
c
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
      integer i, j, k
      real dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
c
c  binary search
c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
c
c  evaluate spline
c
   30 dx = u - x(i)
        speval = b(i) + dx*(2.*c(i) + 3.*dx*d(i))
      return
      end
c******************** end file speval.for ; group trkrlib ******************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@splaan
c******************** start file splaan.for ; group trkrlib ******************
c-------------------------------------------
ccccccccccccccc
      subroutine splaan (n, x, y, b, c, d)
      integer n
      real x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline for which s'(x1)=0.
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      real t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1)=2.*d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = (y(2)-y(1))/d(1)
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end
c******************** end file splaan.for ; group trkrlib ******************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@spleen
      subroutine spleen (n, x, y, b, c, d)                                      
      integer n
      real x(n), y(n), b(n), c(n), d(n)
c
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline for which s'(xn)=0.
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
ccccccccccccccc
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c
c
      integer nm1, ib, i
      real t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = +2.*d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
        c(n)= -(y(n)-y(n-1))/d(n-1)
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@qdint1
c**********************************************************************c
c
c     *********************
c     * SUBROUTINE QDINT1 *
c     *********************
c
c       Set up quadratic interpolation given x(j), f(j), and f'(1):
c
c     usage:
c           CALL QDINT1 (x,f,nx,fp1,cf)
c
c  x   = array of abscissa points, j=1,nx
c  f   = array or ordinate points, j=1,nx
c  nx  = number of points used
c  fp1 = d f(x) / dx at the first point x(1)
c  cf  = output array
c       cf(1,j) = x(j)
c       cf(2,j) = f(j)
c       cf(3,j) = f'(x)  at x(j)+
c       cf(4,j) = 0.5 * f''(x)  at x(j)+
c
c**********************************************************************c
c
c       Assuming the second derivative is uniform in each interval
c
c  d2f/dx2 = cf(4,j) = constant in x(j) =< x < x(j+1)
c
c  then df/dx = cf(3,j) + cf(4,j) * ( x - x(j) ) in x(j) =< x < x(j+1)
c
c  and  f(x) = cf(2,j) + (x-x(j)) * ( cf(3,j) + (x-x(j)) * cf(4,j) )
c
c  in x(j) =< x < x(j+1).  This routine finds the coefficients
c
c  cf(i,j), j=1,nx, i=1,4, given x(j) and f(j) at j=1,nx
c
c  as well as fp1 = df/dx at x(1).
c
c  Continuity of f(x) and f'(x) are imposed throughout the interval.
c
c**********************************************************************c
c
        subroutine QDINT1 (x,f,nx,fp1,cf)
c
        dimension x(*), f(*), cf(4,*)
c
        if ( nx .lt. 2 ) call abortb (6
     & ,'nx .lt. 2 in sbrtn QDINT1')
c
c..the following is an unavoidably recursive loop
c
        cf(3,1) = fp1           ! df/dx at x=x(1)
c
      do 10 j=1,nx-1
        cf(3,j+1) = 2 * ( f(j+1) - f(j) ) / ( x(j+1) - x(j) ) - cf(3,j)
  10  continue
c
      do 20 j=1,nx-1
        cf(1,j) = x(j)
        cf(2,j) = f(j)
        cf(4,j) = ( cf(3,j+1) - cf(3,j) ) / ( 2. * (x(j+1)-x(j)) )
  20  continue
c
      return
      end
c
c@qdeval
c**********************************************************************c
c
c     *********************
c     * SUBROUTINE QDEVAL *
c     *********************
c
c  Evalutation of a piecewise quadratic interpolation,
c  given a close estimate of the interval in which to evaluate.
c
c     usage:
c           CALL QDEVAL (cf,ibound,ib,x,f)
c
c  cf = input:  quadratic interpolation coefficients
c       cf(1,j) = x(j)  , j=1,ibound
c       cf(2,j) = f(j)
c       cf(3,j) = f'(x)  at x(j)+
c       cf(4,j) = 0.5 * f''(x)  at x(j)+
c  ibound = input:  maximum number of coeffients cf(*,j) used
c  ib   = input: estimate of nearest boundary just below x
c       = output: index of cf(1,j) just below x
c  x    = input:   ordinate at which interpolation is to be evaluated
c  f    = output, evaluation of quadratic interpolation
c       = c(2,j) + (x-c(1,j))*(c(3,j) + (x-c(1,j))*c(4,j))
c               for c(1,j) =< x < c(1,j+1)
c
c**********************************************************************c
c
      subroutine QDEVAL (cf,ibound,ib,x,f)
c
      dimension  cf(4,*)
c
c
        ib = max ( ib, 1 )
        ib = min ( ib, ibound-1 )

  10  continue
c
        if ( x .lt. cf(1,ib) ) then
          ib = ib - 1
         if ( ib .lt. 1) call abortb (6
     &      ,'x .lt. cf(1,1) in sbrtn QDEVAL')
          go to 10
        endif
c
  20  continue
c
        if (x .gt. cf(1,ib+1) ) then
          ib = ib + 1
          if (ib .gt. ibound-1 )
     &  call abortb (6,'x .gt. cf(1,ibound-1) in sbrtn QDEVAL')
          go to 20
        endif
c
c  Evaluate quadratic interpolation
c
        f = cf(2,ib) + (x-cf(1,ib)) * ( cf(3,ib) 
     &       + (x-cf(1,ib)) * cf(4,ib) )
c
      return
      end
c
c@qdsolv
c**********************************************************************c
c
c     *********************
c     * SUBROUTINE QDSOLV *
c     *********************
c
c     usage:
c           CALL QDSOLV (cf,x1,x2,f,x)
c
c  x   = solution of the quadratic equation
c
c               f = cf(4)*(x-cf(1))**2 + cf(3)*(x-cf(1)) + cf(2)
c
c  in the range  x1 =< x =< x2 , otherwise, abort.
c
c**********************************************************************c
c
      subroutine QDSOLV (cf,x1,x2,f,x)
c
      dimension  cf(*)
c
c      write (6,110) cf(4),cf(3),cf(2),cf(1),x1,x2,f
c 110  format (' QDSOLV: CF =',1p4e13.5/
c     &        '         x1 =',   e13.5,' x2 =',e13.5,' f =',e13.5)
c
        if ( cf(4) .eq. 0.) then
c
          if (cf(3) .eq. 0. ) 
     &  call abortb (6,'both cf(4) and cf(3) = 0. in sbrtn QDSOLV')
c
          x = cf(1) - ( cf(2) - f ) / cf(3)
c
        else
c
          zsurd = cf(3) * cf(3) - 4. * cf(4) * ( cf(2) - f )
c
          if ( zsurd .lt. 0.)
     &  call abortb (6,'zsurd .lt. 0. in sbrtn QDSOLV')
c
          x = cf(1) + ( - cf(3) + sqrt (zsurd ) ) / ( 2. * cf(4) )
c
          if ( x .lt. min(x1,x2) .or. x .gt. max(x1,x2) )
     &  x = cf(1) + ( cf(2) - f ) / ( cf(4) * ( x - cf(1) ) )
c
        endif
c
        if ( x .lt. min(x1,x2) .or. x .gt. max(x1,x2) )
     &  call abortb (6,'solution out of range in sbrtn QDSOLV')
c
c      write (6,112) x
c 112  format ('          x =',1pe13.5/)
c
      return
      end
c@timint  .../baldur/bald/dintlib.f
c rgb 21 feb 2001 set py=pybkp(iy,1) and return if itimax .lt. 2
c   and cleaned up the routine
c--------1---------2---------3---------4---------5---------6---------7-c
cdoc
c       ------------
c       sbrtn TIMINT  in file DINTLIB
c       ------------
c
c       Sbrtn TIMINT performs a linear interpolation in time,
c  with flat extrapolation outside the range of the
c  ascending portion of the sequence of time breakpoints.
c       The argument list is:
c       call timint (ptime, py, ptbkp, itimax, pybkp, iy, iydim)
c
c  PTIME        = current time.                         (input)
c  PY           = interpolated value                    (output).
c  PTBKP(JT)    = 1-D array of time breakpoints.        (input)
c  ITIMAX       = maximum number of time breakpoints.   (input)
c  PYBKP(JY,JT) = Y-values at the breakpoints.          (input)
c  IY           = value of JY index                     (input)
c  IYDIM        = dimension of 1st index in PYBKP       (input)
c
c       NOTE:  IF ( PTBKP(2) .LE. PTBKP(1) )  THEN PY = PYBKP(IY,1)
c              IF ( PTIME .LE. PTBKP(1) )  THEN PY = PYBKP(IY,1)
c              IF J IS THE INDEX OF THE LAST ELEMENT IN AN ASCENDING
c                       SEQUENCE OF TIME BREAKPOINTS
c                       AND IF ( PTIME .GE. PTBKP(J) )
c                       THEN PY = PYBKP(IY,J)
c
cend
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine timint (ptime,py,ptbkp,itimax,pybkp,iy,iydim)
c
        dimension  ptbkp(*), pybkp(iydim,*)
c
        if ( itimax .lt. 2 .or.
     &        ptbkp(2) .le. ptbkp(1) .or. ptime .le. ptbkp(1) ) then
          py = pybkp(iy,1)
          return
        endif
c
c..Find the appropriate interval
c
c        py = pybkp(iy,itimax)
c
c        if ( ptime .ge. ptbkp(itimax) ) return
c
c..time array must be monotonically increasing
c
        do jt=2,itimax
          if ( ptbkp(jt) .le. ptbkp(jt-1) ) then
            py = pybkp(iy,jt-1)
            go to 12
          endif
c
c..Interpolate when interval is found
c
          if ( ptbkp(jt) .ge. ptime ) then
            zint = ( ptbkp(jt) - ptime ) / ( ptbkp(jt) - ptbkp(jt-1) )
            py = pybkp(iy,jt-1) * zint + pybkp(iy,jt) * ( 1.-zint )
            go to 12
          endif
c
        enddo
c
  12  continue
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c

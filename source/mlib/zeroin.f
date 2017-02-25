c@zeroin   .../baldur/code/util/zeroin.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      real function zeroin ( ax, bx, func, tol )
c
      implicit none
c
      real ax, bx, func, tol
c
      external func
c
c     zeroin computes the zero of function func(x) for ax < x < bx
c
c  ax    left endpoint of the initial interval  (input)
c  bx    right endpoint of the initial interval (input)
c  func  name of user supplied function subprogram
c        for any ax < x < bx
c  tol   desired length of the interval of uncertainty
c        of the final result (.ge. 0.0)
c
c        It is assumed that func(ax) and func(bx) have opposite signs
c  without a check.  zeroin returns a zero x in the given interval
c  ax, bx  to within a tolerance 4 * macheps * abs(x) + tol, where
c  macheps is the relative machine precision.
c        This function subprogram is a slightly modified translation
c  of the ALGOL 60 procedure zero given in Richard Brent, Algorithms
c  For Minimization Without Derivatives, Prentice Hall, inc. (1973).
c        From G. E. Forsythe, M. A. Malcolm, and C. B. Moler,
c  "Computer Methods for Mathematical Computations,"
c  Prentice-Hall (1977) pp. 163-166.
c
      real a, b, c, d, e, eps, fa, fb, fc, tol1, xm, p, q, r, s
c
      integer inital
      data inital / 0 /
      save inital, eps
c
c..compute eps, the relative machine precision
c
      if ( inital .lt. 1 ) then
        eps = 1.0
  10    eps = eps / 2.0
        tol1 = 1.0 + eps
        if ( tol1 .gt. 1.0 ) go to 10
        inital = 1
      endif
c
c..initialization
c
      a  = ax
      b  = bx
      fa = func(a)
      fb = func(b)
c
c..begin step
c
  20  c  = a
      fc = fa
      d  = b - a
      e  = d
  30  continue
      if ( abs(fc) .lt. abs(fb) ) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      endif
c
c..convergence test
c
      tol1 = 2.0 * eps * abs(b) + 0.5 * tol
      xm   = 0.5 * ( c - b )
      if ( abs(xm) .le. tol1 ) go to 90
      if ( fb .eq. 0.0 ) go to 90
c
c..is bisection necessary ?
c
      if ( abs(e) .lt. tol1 ) go to 70
      if ( abs(fa) .le. abs(fb) ) go to 70
c
c..is quadratic interpolation possible ?
c
      if ( a .eq. c ) then
c
c..use linear interpolation
c
        s = fb / fa
        p = 2.0 * xm * s
        q = 1.0 - s
      else
c
c..use inverse quadratic interpolation
c
        q = fa / fc
        r = fb / fc
        s = fb / fa
        p = s * ( 2.0 * xm * q * ( q - r ) - ( b - a )*( r - 1.0 ) )
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 )
      endif
c
c..adjust signs
c
      if ( p .gt. 0.0 ) q = -q
      p = abs(p)
c
c..is interpolation acceptable ?
c
      if ( ( (2.0*p) .ge. ( 3.0*xm*q - abs(tol1*q) ) )
     &   .or.  ( p .ge. abs(0.5*e*q) ) ) go to 70
c
c..accept interpolation
c
        e = d
        d = p / q
      go to 80
c
c..bisection
c
  70  continue
        d = xm
        e = d
c
c..complete the step
c
  80  continue
      a  = b
      fa = fb
      if ( abs(d) .gt. tol1 ) b = b + d
      if ( abs(d) .le. tol1 ) b = b + sign(tol1,xm)
      fb = func(b)
      if ( ( fb*(fc/abs(fc)) ) .gt. 0.0 ) go to 20
      go to 30
c
c..done
c
  90  continue
      zeroin = b
c
      return
      end

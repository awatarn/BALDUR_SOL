c@cubint
c rgb 03-jun-96 removed if (ier .gt. 32) call abortb (6,text)
c rgb 30-apr-96 implemented spline, sevalb, and spntgl routines
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE cubint *
c     *********************
c
c     Purpose:
c           perform and evaluate a piecewise cubic interpolation
c
c     usage:
c           CALL cubint (X,Y,nxin, lintrp, C,ic,
c
c                                  lintrp< 0  : C(ic,3)  input
c                                        = 0  : C(ic,3)  calc. Cub.Splns
c                                        > 0  : C(ic,3)  calc. QuasHerms
c
c     >                  U,S,m   , levalu,
c
c                                  levalu= 0 ...S= s(u)
c                                  levalu= 1 ...S= ds/du
c                                  levalu= 2 ...S= d2s/du2
c                                  levalu= 3 ...S= integr(s.du)
c
c     >                  xsym,lsym,text)
c
c                             lsym= 0  ...no symmetry on input fn y(x)
c                                 = +1 ...EVEN symmetry at xsym
c                                 = -1 ...ODD  symmetry at xsym
c
c
c     ..................................................................
c
c     Purpose: perform and evaluate a piecewise cubic interpolation
c
c     Sbrtn cubint calls IMSL sbrtns IQHSCU, ICSEVU, DCSEVU, DCSQDU
c     in order to both set up and evaluate a cubic spline
c     piecewise-cubic interpolation through a given set of data points.
c
c     There are 13 arguments; 11 input; 1 output; 1 input/output
c
c x   - vector of length nx containing the abscissae of nx data points
c       must be in ascending order
c y   - input vector of length nx containing ordinates of data points
c nxin - input number of elements of x and y.  must be > 2.
c lintrp - controls the method of interpolation
c       if lintrp is negative, spline coefficients c are not computed
c                           they must then be provieded by the user
c       if lintrp = 0, cubic splines are used to compute coefficients c
c                      third derivative is made to be continuous
c                      at the second and penultimate knots
c       if lintrp > 0, quaisi-Hermite piecewise cubic interpolation
c                      is used to compute coeffieients c
c c   - output nx by 3 matrix of spline coefficients
c        input if lintrp < 0
c ic  - input row dimension of matrix c. must be .ge. nx-1
c u   - input vector of abscissae where interpolation is desired
c s   - output vector of interpolated values at x=u
c m   - input number of interpolated values u and s
c       evaluation is skipped if m .le. 0.
c levalu - controls evaluation of the interpolation
c     = 0 for s(j) to be the interpolated values at u(j), j=1,...,m
c     = 1 for s(j) to be the interpolated first diervatives at u(j)
c     = 2 for s(j) to be the interpolated second derivatives at u(j)
c     = 3 for s(j) to represent the integral from u(1) to u(j)
c xsym - input abscissa around which symmetry is imposed
c lsym - input integer to control evaluation near the endpoints
c     = 0 for no symmetry condition on input function y(x)
c     = 1 for even symmetry around xsym  y(xsym-x) = y(x-xsym)
c     = -1 for odd symmetry around xsym  y(xsym-x) = - y(x-xsym)
c text - input character string which will appear on output unit 6
c       in the event of a fatal or warining error.
c xsrc, yscr, uscr - internal arrays, must keep nx .le. 101
c
c     ------------------------------------------------------------------
c  RGB 02-FEB-87 16:00 Moved sbrtn to 11040 .wbaldn1 wsutil dintlib
c         changed dimension from 101 to 301
c     ------------------------------------------------------------------
c***********************************************************************
c
c
      subroutine cubint (x,y,nxin,lintrp,c,ic,
     >                   u,s,m   ,levalu,xsym,lsym, text)
c
      dimension  x(*), y(*), c(*), u(*), s(*)
      dimension  xscr(301),yscr(301),  uscr(301)
      character  *(*)  text
c
      save  i1, i2, i3
c
c     ------------------------------------------------------------------
c
c
      if (m .gt. 301)
     & call abortb (6,'abort in sbrtn cubint, m .gt. 301')
c
      nx = abs ( nxin )
c
CLL   1.   set up reflection symmetry
c
c  extend the input arrays (x,y)
c  to take advantage of symmetry, if possible
c
      do 10 j=1,nx
        xscr(j) = x(j)
        yscr(j) = y(j)
  10  continue
c
c
      if (lsym .ne. 0) then
         xref = xsym   ! reflection point for symmetry conditions
         ssym = sign (1,lsym)
        dxmin = 0.2 * (x(2) - x(1))
c
c  dxmin is the minimum distance allowed between extended points
c
      if (xsym .gt. 2.*x(1) - x(3) + dxmin .and. xsym .lt. x(3)) then
c
c..reflect up to the first three points across the symmetry line
c
      i = 0
      do 12 j=3,1,-1
      if (xsym .le. x(j)) then
        i = i + 1
        xscr(i) = 2. * xsym - x(j)
        yscr(i) = ssym * y(j)
      endif
  12  continue
c
c..fill in the rest of the scratch array after reflected points
c
      do 14 j=1,nx
      if (x(j) .gt. xsym) then
        i = i + 1
        xscr(i) = x(j)
        yscr(i) = y(j)
      endif
  14  continue
      nx = i   ! last point in the extended array
c
      endif   ! end of reflecting first three points
c
         else   ! no symmetry condition (lsym = 0)
c
         xref = x(1)   ! this will not interfere with extrap beyond ends
c
      endif   ! end of setting up reflection symmetry
c
      if ( nx .gt. ic ) call abortb (6,' nx .gt. ic in sbrtn cubint')
c     ==================================================================
c
CLL   2.   compute interpolating coefficients c(i,j) i=1,ic j=1,3
c
      i1 = 1
      i2 = 1 + nx
      i3 = 1 + 2*nx
c
c  for now, call spline to set up cubic spline coefficients
c  implement quasi-Hermite later
c
cbate      if ( lintrp .le. 0 )
      call spline (nx,xscr,yscr,c(i1),c(i2),c(i3))
c
c      if (lintrp .le. 0) call csint (nx,xscr,yscr,c,c(nx+1))
c      if (lintrp .ge. 1) call csakm (nx,xscr,yscr,c,c(nx+1))
c
c      if (lintrp .eq. 0) call icsccu(xscr,yscr,nx,c,ic,ier) ! cubic splines
c      if (lintrp .gt. 0) call iqhscu(xscr,yscr,nx,c,ic,ier) ! quasi-Hermite
c      if (ier .gt. 128) call abortb (6,text)
c     ------------------------------------------------------------------
c
c..should the interpolation be evaluated?
c
      if (m .lt. 1) return
c
      if (levalu .lt. 0 .or. levalu .gt. 3) then
        write (6,*) ' levalu out of range in sbrtn cubint'
        write (6,*) text
        return
      endif
c     ==================================================================
c
c
CLL   3.   evaluate interpolation
c
      if (levalu .eq. 0) then
c
c  force abscissae to be within bounds
c
      do j=1,m
        uscr(j) = max (u(j), xscr(1), 2*xref-u(j))
        uscr(j) = min (uscr(j), xscr(nx))
      enddo
c
c  evaluate piecewise-cubic interpolation
c
      call sevalb ( 0, nx, xscr, yscr, c(i1), c(i2), c(i3), m, uscr, s )
c
c      do j=1,m
c        zu = uscr(j)
c            s(j) = csval (uscr(j),nx-1,c,c(nx+1))
c      enddo
c
c      call icsevu (xscr,yscr,nx,c,ic,uscr,s,m,ier)
c      if (ier .gt. 32) call abortb (6,text)
c
c  apply symmetry condition around xref and
c  linearly extrapolate beyond end points
c  later, linear extrapolation may be replaced with quadratic extrapolation
c
      do 22 j=1,m
c
c  extrapolate off high end
c
      x2 = max ( u(j), 2*xref - u(j) ) ! abscissa or reflected abscissa
      dydx2 =
     & (yscr(nx) - yscr(nx-1)) / (xscr(nx) - xscr(nx-1))   ! slope at xscr(nx)
      if (x2 .gt. xscr(nx)) s(j) = yscr(nx) + dydx2 * (x2 - xscr(nx))
c
c  extrapolate and apply symmetry off low end
c
      if (u(j) .lt. xref) then
         if (lsym .eq. 0) then  ! extrapolate, note xref=xscr(1)
         dydx1 =
     &     (yscr(2) - yscr(1)) / (xscr(2) - xscr(1))   ! slope at xscr(1)
         s(j) = yscr(1) + dydx1 * (u(j) - xscr(1))
         else                !apply symmetry condition across xref=xsym
         s(j) = ssym * s(j)
         endif
      endif
  22  continue
c
      return
      endif
c
c     ==================================================================
c
CLL   4.   evaluate first derivatives
c
      if (levalu .eq. 1) then
c
      do j=1,m   ! force abscissa to be within bounds
        uscr(j) = max (u(j), xscr(1), 2*xref-u(j))
        uscr(j) = min (uscr(j), xscr(nx))
      enddo
c
      call sevalb ( 1, nx, xscr, yscr, c(i1), c(i2), c(i3), m, uscr, s )
c
c      do j=1,m
c        s(j) = speval ( nx, uscr, xsrc, yscr, c(i1), c(i2), c(i3) )
c            s(j) = csder (1,uscr(j),nx-1,c,c(nx+1))
c      enddo
c
c      m2 = 0
c      call dcsevu (xscr,yscr,nx,c,ic,uscr,s,m,sd,m2,ier)
c      if (ier .gt. 32) call abortb (6,text)
c
c  apply symmetry condition around xref and
c  linearly extrapolate beyond end points
c  later, linear extrapolation may be replaced with quadratic extrapolation
c
      do 32 j=1,m
      x2 = 2*xref - u(j)   ! reflected abscissa
c
c  extrapolate off high end
c
      if (x2 .gt. xscr(nx)) s(j) =
     @ (yscr(nx) - yscr(nx-1)) / (xscr(nx) - xscr(nx-1))  ! slope at xscr(nx)
c
c  extrapolate and apply symmetry off low end
c
      if (u(j) .lt. xref) then
         if (lsym .eq. 0) then  ! extrapolate, note xref=xscr(1)
         s(j) = (yscr(2) - yscr(1)) / (xscr(2) - xscr(1))  ! slope at xscr(1)
         else                !apply symmetry condition across xref=xsym
         s(j) = - ssym * s(j)
         endif
      endif
c
  32  continue
c
      return
      endif
c
c     ==================================================================
c
CLL   5.   evaluate second derivatives
c
      if (levalu .eq. 2) then
c
      do j=1,m   ! force abscissa to be within bounds
        uscr(j) = max (u(j), xscr(1), 2*xref-u(j))
        uscr(j) = min (uscr(j), xscr(nx))
      enddo
c
      call sevalb ( 2, nx, xscr, yscr, c(i1), c(i2), c(i3), m, uscr, s )
c
c      do 41 j=1,m
c            s(j) = csder (2,uscr(j),nx-1,c,c(nx+1))
c  41  continue
c
c      m1 = 0
c      call dcsevu (xscr,yscr,nx,c,ic,uscr,sd,m1,s,m,ier)
c
crgb960603      if (ier .gt. 32) call abortb (6,text)
c
c  apply symmetry condition around xref and
c  linearly extrapolate beyond end points
c  later, linear extrapolation may be replaced with quadratic extrapolation
c
      do 42 j=1,m
        x2 = 2*xref - u(j)   ! reflected abscissa
c
c  extrapolate off high end
c
      if (x2 .gt. xscr(nx)) s(j) = 0.
c
c  extrapolate and apply symmetry off low end
c
      if (u(j) .lt. xref) then
         if (lsym .eq. 0) then  ! extrapolate, note xref=xscr(1)
         s(j) = 0.
         else                   ! apply symmetry condition across xref=xsym
         s(j) = ssym * s(j)
         endif
      endif
c
  42  continue
c
      return
      endif
c
c     ==================================================================
c
CLL   6.   integrate from u(1),s(1) to u(j),s(j), j=2,...,m
c
      if (levalu .eq. 3 .and. m .gt. 1) then
c
      s(1) = 0.0
c
      do j=2,m
c
      ds = spntgl ( nx, u(j-1), u(j), xscr, yscr, c(i1), c(i2), c(i3) )
c
cbate        ds = csitg (u(j-1),u(j),nx-1,c,c(nx+1))
c
c      call dcsqdu (xscr,yscr,nx,c,ic,u(j-1),u(j),ds,ier)
c      if (ier .gt. 32) call abortb (6,text)
c
        s(j) = s(j-1) + ds
c
      enddo
c
      return
      endif
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c

c@int1d.f
c  rgb 12-dec-95 fixed extrapolation of derivatives
c  rgb 06-apr-95 generalized extrapolation of Akima's algorithm
c  rgb 05-apr-95 added option to compute derivatives (kderiv)
c  rgb 31-mar-95 added kestrap extrapolation options (kextrap)
c--------1---------2---------3---------4---------5---------6---------7-c
c  Interpolation of data over a 1-D array
c
      subroutine int1d (kintrp, kextrap, kderiv, kx, kskip, px, pf
     &  , ku, pu, kxlast, pw)
c
c  kintrp  controls the type of interpolation used (see below)   (input)
c  kexptrap controls the type of extrapolation used at the ends  (input)
c  kderiv  compute derivative of order kderiv (0, 1, 2, 3)
c  kx     = number of elements in the 1-D array px   (input)
c  kskip  = first dimension of 2-D array pf          (input)
c  px     = 1-D array of points on 1-axis of f(x)    (input)
c  pf     = 2-D array of values of f(y,x)            (input)
c  ku     = number of points in 1-D array pu         (input)
c         = number of values given for pf(kx,ku)
c  pu     = 1-D array of x values for interpolation  (input)
c  kxlast  = last value for which                     (input and output)
c             px(kxlast) .le. pu(kx) .lt. px(kxlast+1)
c  pw     = 1-D array of interpolated values         (output)
c
c  Note:  if ku .le. kskip, then ku is taken to be the number of values
c  in the second dimension of pf(kx,ku)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  if kintrp .ne. 2, linear interpolation is used
c
c  if kintrp = 2, then Hiroshi Akima's algorithm ACM 433 is used
c from Comm. ACM vol 15, no 10, page 914--918 (Oct. 1972).
c
c  kextrap = 0 for flat extrapolation off the ends
c          = 1 for linear extrapolation
c          = 2 for quadratic extrapolation
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      real px(*), pf(kskip,*), pu(*), pw(*)
c
      integer kintrp, kderiv, kextrap, kx, kskip, ku, kxlast
     &  , iskip, iu
     &  , ixpu, ix, imn, imx, jx, ju
c
      real zm0, zm1, zm2, zm3, zm4, zm5, zm6
     &   , zw1, zw2, zw3, zw4, zw5, zw6
     &   , zt3, zt4, za3
     &   , zx3, zx4, zy3, zy4
     &   , zq0, zq1, zq2, zq3
     &   , zdx, zdref, zepslon
c
c
c
c..make sure input values make sense
c
      if ( kxlast .gt. kx ) kxlast = 2
      kxlast = max ( 2, kxlast )
c
      iskip = max ( kskip, 1 )
      iu    = max ( ku, 1 )
c
      ixpu = 0
c
      zepslon = 1.e-9
c
c
c..case with one point given
c
      if ( kx .lt. 2 ) then
        do ju=1,ku
          if ( kderiv .lt. 1 ) then
            pw(ju) = pf(1,1)
          else
            pw(ju) = 0.0
          endif
        enddo
        kxlast = 1
        return
      endif
c
c
c..main do loop over values of pu
c
      do ju=1,ku
c
c..compute ix value for which
c    pu(ju) .ge. px(ix-1) .and. pu(ju) .lt. px(ix)
c
        if (pu(ju) .ge. px(kxlast-1) .and. pu(ju) .lt. px(kxlast)) then
          ix = kxlast
        elseif ( kx .eq. 2 ) then
          ix = 2
        elseif ( pu(ju) .ge. px(kx) ) then
          ix = kx
        elseif ( pu(ju) .lt. px(2) ) then
          ix = 2
        else
          imn = 2
          imx = kx
  10      ix = ( imn + imx ) / 2
          if (pu(ju) .ge. px(ix) ) then
            imn = ix + 1
          else
            imx = ix
          endif
          if ( imx .gt. imn ) go to 10
          ix = imx
        endif
c
c
c..flat extrapolation outside interval
c
      if ( kextrap .lt. 1 .and. pu(ju) .lt. px(1) ) then
c
        if ( kderiv .lt. 1 ) then
          pw(ju) = pf(1,1)
        else
          pw(ju) = 0.0
        endif
c
      elseif ( kextrap .lt. 1 .and. pu(ju) .gt. px(kx) ) then
c
        if ( kderiv .lt. 1 ) then
          pw(ju) = pf(1,kx)
        else
          pw(ju) = 0.0
        endif
c
c
c..linear interpolation
c
      elseif ( kintrp .lt. 2 ) then
c
        if ( kderiv .lt. 1 ) then
          pw(ju) = ( pf(1,ix-1) * ( px(ix) - pu(ju) )
     &             + pf(1,ix) * ( pu(ju) - px(ix-1) ) )
     &             / ( px(ix) - px(ix-1) )
c
          if ( pu(ju) .lt. px(1) .and. kextrap .lt. 1 )
     &      pw(ju) = pf(1,1)
c
          if ( pu(ju) .gt. px(kx) .and. kextrap .lt. 1 )
     &      pw(ju) = pf(1,kx)
c
        elseif ( kderiv .eq. 1 ) then
          pw(ju) = (pf(1,ix)-pf(1,ix-1))/(px(ix)-px(ix-1))
          if ( pu(ju) .lt. px(1) .and. kextrap .lt. 1 ) pw(ju) = 0.0
          if ( pu(ju) .gt. px(kx) .and. kextrap .lt. 1 ) pw(ju) = 0.0
c
        else
          pw(ju) = 0.0
c
        endif
c
c
c..interpolation by Akima's method
c
      else
c
c    This algorithm computes a cubic polynomial in each interval
c  f(x) = a_0 + a_1 dx + a_2 dx**2 + a_3 dx**3
c  where dx = x - x_j in the interval x_{j-1} .le. x .lt. x_j.
c  Let zm3 = zm(j-1/2) = ( f_j - f_{j-1} ) / ( x_j - x_{j-1} )
c  (shift left for zm2, shift right for zm4, ...)
c  Let w_j = | zm(j+1/2) + zm(j-1/2) |
c  Note, if w_{j-1} = w_{j+1} = 0.0, then set w_{j-1} = w_{j+1} = 0.5
c  At each end of the interval, approximate the derivative by
c  df/dx_j = (w_{j-1} zm(j+1/2) + w_{j+1} zm(j-1/2))/(w_{j-1}+w_{j+1})
c  Then
c  a_0 = f_{j-1}
c  a_1 = df/dx_{j-1}
c  a_3 = 3 zm(j-1/2) - 2 df/dx_{j-1} - df/dx_j
c  a_4 = df/dx_{j-1} + df/dx_j - 2 zm(j-1/2)
c
c    Notation:
c  zm3 = slope within the interval
c  zt3 = d f / d x at jx-1
c
c..check if the desired point is in the same interval as previous point
c  If yes, then skip to the computation of the polynomial
c
        if ( ix .eq. ixpu ) go to 20
        ixpu = ix
c
c..compute finite differences
c
        jx = max ( min ( ix, kx) , 2 )
c
        zx3 = px(jx-1)
        zy3 = pf(1,jx-1)
        zx4 = px(jx)
        zy4 = pf(1,jx)
        za3 = px(jx) - px(jx-1)
        zm3 = ( pf(1,jx) - pf(1,jx-1) ) / ( px(jx) - px(jx-1) )
c
        if ( kextrap .lt. 1 ) then
          zm0 = 0.0
          zm1 = 0.0
          zm2 = 0.0
          zm4 = 0.0
          zm5 = 0.0
          zm6 = 0.0
        else
          zm0 = zm3
          zm1 = zm3
          zm2 = zm3
          zm4 = zm3
          zm5 = zm3
          zm6 = zm3
        endif
c
        if ( jx-4 .ge. 1 )
     &    zm0 = ( pf(1,jx-3) - pf(1,jx-4) ) / ( px(jx-3) - px(jx-4) )
c
        if ( jx-3 .ge. 1 )
     &    zm1 = ( pf(1,jx-2) - pf(1,jx-3) ) / ( px(jx-2) - px(jx-3) )
c
        if ( jx-2 .ge. 1 )
     &    zm2 = ( pf(1,jx-1) - pf(1,jx-2) ) / ( px(jx-1) - px(jx-2) )
c
        if ( jx+1 .le. kx )
     &    zm4 = ( pf(1,jx+1) - pf(1,jx) ) / ( px(jx+1) - px(jx) )
c
        if ( jx+2 .le. kx )
     &    zm5 = ( pf(1,jx+2) - pf(1,jx+1) ) / ( px(jx+2) - px(jx+1) )
c
        if ( jx+3 .le. kx )
     &    zm6 = ( pf(1,jx+3) - pf(1,jx+2) ) / ( px(jx+3) - px(jx+2) )
c
c..extrapolate derivatives if necessary
c
        if ( kextrap .eq. 1 ) then
          if ( jx-2 .lt. 1  ) zm2 = zm3
          if ( jx+1 .gt. kx ) zm4 = zm3
          if ( jx-3 .lt. 1  ) zm1 = zm2
          if ( jx+2 .gt. kx ) zm5 = zm4
          if ( jx-4 .lt. 1  ) zm0 = zm1
          if ( jx+3 .gt. kx ) zm6 = zm5
        elseif ( kextrap .gt. 1 ) then
          if ( jx-2 .lt. 1  ) zm2 = 2.0 * zm3 - zm4
          if ( jx+1 .gt. kx ) zm4 = 2.0 * zm3 - zm2
          if ( jx-3 .lt. 1  ) zm1 = 2.0 * zm2 - zm3
          if ( jx+2 .gt. kx ) zm5 = 2.0 * zm4 - zm3
          if ( jx-4 .lt. 1  ) zm0 = 2.0 * zm1 - zm2
          if ( jx+3 .gt. kx ) zm6 = 2.0 * zm5 - zm4
        endif
c
c..now estimate the derivatives at the ends of the interval zt3 and zt4
c
        zw1 = abs ( zm1 - zm0 )
        zw2 = abs ( zm2 - zm1 )
        zw3 = abs ( zm3 - zm2 )
        zw4 = abs ( zm4 - zm3 )
        zw5 = abs ( zm5 - zm4 )
        zw6 = abs ( zm6 - zm5 )
c
        zdref = abs ( (abs(zy3) + abs(zy4))/abs(zx4-zx3) )
c
        zt3 = 0.5 * ( zm2 + zm3 )
        if ( zw4 + zw2 .gt. zepslon*zdref )
     &    zt3 = ( zw4 * zm2 + zw2 * zm3 ) / ( zw4 + zw2 )
        if ( jx-1 .le. 1 .and. kextrap .lt. 1 ) zt3 = 0.0
c
        zt4 = 0.5 * ( zm3 + zm4 )
        if ( zw5 + zw3 .gt. zepslon*zdref )
     &    zt4 = ( zw5 * zm3 + zw3 * zm4 ) / ( zw5 + zw3 )
        if ( jx .ge. kx .and. kextrap .lt. 1 ) zt4 = 0.0
c
c..extrapolate the derivatives if necessary
c
c        if ( jx-1 .le. 1  .and.  kextrap .gt. 0 ) then
c
c          zt5 = 0.5 * ( zm4 + zm5 )
c          if ( zw6 + zw4 .gt. zepslon*zdref )
c     &      zt5 = ( zw6 * zm4 + zw4 * zm5 ) / ( zw6 + zw4 )
c          if ( jx+1 .ge. kx .and. kextrap .lt. 1 ) zt5 = 0.0
c
c          zt3 = 2.0 * zt4 - zt5
c
c        endif
c
c        if ( jx .ge. kx  .and.  kextrap .gt. 0 ) then
c
c          zt2 = 0.5 * ( zm1 + zm2 )
c          if ( zw3 + zw1 .gt. zepslon*zdref )
c     &      zt2 = ( zw3 * zm1 + zw1 * zm2 ) / ( zw3 + zw1 )
c          if ( jx-2 .le. 1 .and. kextrap .lt. 1 ) zt2 = 0.0
c
c          zt4 = 2.0 * zt3 - zt2
c
c        endif
c
c
c
c..determine coefficients
c
        zq0 = zy3
        zq1 = zt3
        zq2 = ( 2.0 * ( zm3 - zt3 ) + zm3 - zt4 ) / za3
        zq3 = ( -2.0 * zm3 + zt3 + zt4 ) / ( za3**2 )
c
  20    zdx = pu(ju) - zx3
c
        if ( kderiv .lt. 1 ) then
c
          pw(ju) = zq0 + zdx * ( zq1 + zdx * ( zq2 + zdx * zq3 ) )
c
          if ( pu(ju) .lt. px(1) ) then
            if ( kextrap .lt. 1 ) then
              pw(ju) = pf(1,1)
            elseif ( kextrap .eq. 1 ) then
              pw(ju) = pf(1,1) + zt3 * ( pu(ju) - zx3 )
            endif
          endif
c
          if ( pu(ju) .gt. px(kx) ) then
            if ( kextrap .lt. 1 ) then
              pw(ju) = pf(1,kx)
            elseif ( kextrap .eq. 1 ) then
              pw(ju) = pf(1,kx) + zt4 * ( pu(ju) - zx4 )
            endif
          endif
c
        elseif ( kderiv .eq. 1 ) then
c
          pw(ju) = zq1 + 2.0 * zq2 * zdx + 3.0 * zq3 * zdx**2
c
          if ( pu(ju) .lt. px(1) ) then
            if ( kextrap .lt. 1 ) then
              pw(ju) = 0.0
            elseif ( kextrap .eq. 1 ) then
              pw(ju) = zt3
            endif
          endif
c
          if ( pu(ju) .gt. px(kx) ) then
            if ( kextrap .lt. 1 ) then
              pw(ju) = 0.0
            elseif ( kextrap .eq. 1 ) then
              pw(ju) = zt4
            endif
          endif
c
        elseif ( kderiv .eq. 2 ) then
c
          pw(ju) = 2.0 * zq2 + 6.0 * zq3 * zdx
          if ( pu(ju) .lt. px(1) .and. kextrap .lt. 2 ) pw(ju) = 0.0
          if ( pu(ju) .gt. px(kx) .and. kextrap .lt. 2 ) pw(ju) = 0.0
c
        elseif ( kderiv .eq. 3 ) then
c
          pw(ju) = 6.0 * zq3
          if ( pu(ju) .lt. px(1) .and. kextrap .lt. 2 ) pw(ju) = 0.0
          if ( pu(ju) .gt. px(kx) .and. kextrap .lt. 2 ) pw(ju) = 0.0
c
        else
c
          pw(ju) = 0.0
c
        endif
c
      endif
c
c..end of main do loop
c
      enddo
c
c..set kxlast for next call
c
      kxlast = ix
c
c
      return
      end


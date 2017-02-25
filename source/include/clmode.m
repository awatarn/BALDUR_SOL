c@clmode for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cmodes/
     &      nmodes,           modem(kpmode),    moden(kpmode)
     & ,    mdmaxn,           mdmaxm(kpmode),   mdtrix(kpmode,kpmode)
     & ,    qmode(kpmode),    xhalf(kpmode),    xdelta(kpmode)
     & ,    dprime(kpmode),   lxmode(kpmode)
     & ,    xdelti(kpmode),   resid(kpmode)
     & ,    qhalf1(kpmode),   resid1(kpmode),   dqhalf(kpmode)
     & ,    qhalf2(kpmode),   resid2(kpmode),   dqfact
     & ,    b1amp(kpmode),    b1scal(kpmode),   ncount
     & ,    deltax,           xhfmin,           xhfmax
     & ,    xwall,            rxwall(kpmode),   xistop
     & ,    xbkpnt(kpmode),   ibkpnt(kpmode),   nbkpnt
     & ,    xmxisl
     & ,    xmode(kpmode),    xinner(kpmode),   xouter(kpmode)
     & ,    xinb(kpmode),     rtinb(kpmode),    djinb(kpmode)
     & ,    xoutb(kpmode),    rtoutb(kpmode),   djoutb(kpmode)
     & ,    b1xqhf(kpmode),   b1xisl(kpmode),   sgnisl(kpmode)
     & ,    b1xedg(kpmode),   b1tedg(kpmode),   b1zedg(kpmode)
c
c   nmodes        number of perturbed helical modes
c
c   modem(j)      mode number m (poloidal harmonics) for j=1,nmodes
c
c   moden(j)      mode number n (toroidal harmonics) for j=1,nmodes
c
c                 Note:  These must be in order
c                        modem(j)/moden(j) .lt. modem(j+1)/moden(j+1)
c
c   mdmaxn        maximum of the mode numbers n (toroid harmonics)
c
c   mdmaxm(n)     maximum of the mode numbers m (poloidal harmonics)
c                        for each toroidal harmonic n
c
c   mdtrix(m,n)   = index j associated with m=modem(j) and n=moden(j)
c
c                 These last three variables are set in sbrtn modset
c
c   qmode(j)      = m / n  for each mode
c
c   xhalf(j)      the width of the jth magnetic island in
c                     x space (input/output)
c
c   xdelta(j)     distance between the x-point of each island from
c                     the mode rational surface divided by xhalf(j)
c
c   dprime(j)     delta-prime for tearing mode
c
c   lxmode(j)     = -1 if there is no mode rational surface at all
c                      (ie xmode(j) .lt. 0.)
c                 = 0  if the mode rational surface is within the plasma
c                 = 1  if the mode rational surface
c                         is between the plasma and the wall
c                 = 2  if the mode rational surface is outside the wall
c
c   xdelti(j)     = temporary value of xdelta(j) computed in sbrtn ibcomp
c                   and transferred to xdelta(j) in sbrtn iwidth
c
c   resid(j)      = residual computed in sbrtn ibcomp
c                 = B^{1 x}_{mn} at the wall
c
c   qhalf1(j)     island half-width during first iteration
c
c   resid1(j)     residual resulting from iteration with qhalf1(j)
c
c   dqhalf(j)     increment to be used for next iteration
c                       of island half-width
c
c   qhalf2(j)     island half-width during second iteration
c
c   dqfact        factor by which dqhalf(j) is amplified
c                       when residuals remain the same sign
c                       defaulted to 1.5
c                       dqfact = cislnd(20) if cislnd(20) .gt. 1.
c
c   resid2(j)     residual resulting from iteration with qhalf2(j)
c
c   b1amp(j)      the amplitute of the asymptotic expansion
c                     at xi = deltaxi.  bv(j)=b1amp(j)*xi**(modem(j)-1)
c
c   b1scal(j)     scale factor applied to perturbed magnetic field
c                 determined from magnetic island width
c
c   ncount        number of calls to sbrtn iwidth from zero finder
c                       reset to 0 before each new call to zz#nbf
c
c   deltax        small value of xi at which b1amp(j) is used (input)
c
c   xhfmin        minimum island width as a function of q-value (input)
c
c   xhfmax        maximum allowable xhalf(j)  (set in sbrtn island)
c
c   xwall         r_wall / r_plasma  (input)
c
c   rxwall(j)     ( 1. - xwall ** (-2*m) ) / ( 1. + xwall ** (-2*m) )
c                     computed in sbrtn island and used in sbrtn iwidth
c
c   xistop        stop integration in sbrtn idiftr when pxi .gt. xistop
c                   for diagnostic purposes
c
c            Variables computed in sbrtn ixmode:
c
c   xmode(j)      xi values at mode rational surfaces (output)
c
c   xbkpnt(j)     xi breakpoints for use by the ode solver in ibcomp
c
c   ibkpnt(j)     index with each breakpoint keyed to mode number
c                 = 0 if breakpoint does not correspond to a mode
c
c   nbkpnt        = number of breakpoints xbkpnt(j)
c
c   xmxisl        maximum value of xmode(j) for which island width computed
c
c            Variables computed in sbrtn islbnd:
c
c   xinner(j)     xi values at inner edges of magnetic islands
c
c   xouter(j)     xi values at outer edges of magnetic islands
c
c   xinb(j)       interpolation breakpoint at inboard edge of island
c   xoutb(j)      interpolation breakpoint at outboard edge of island
c
c   rtinb(j)      value of (xi-xmode(j)) / ( xhalf(j) * (nq-m) ) at
c   rtoutb(j)     inboard and outboard breakpoints xinb(j) and xoutb(j)
c
c   djinb(j)      value of d J^{0\zeta} / d xi at
c   djoutb(j)     inboard and outboard breakpoints xinb(j) and xoutb(j)
c
c            Variables computed in sbrtn ibcomp:
c
c   b1xqhf(j)     b1uxi(j) computed from island half-width
c                   at each mode rational surface
c
c   b1xisl(j)     b1uxi(j) computed by integrating tearing mode equations
c                   up to each mode rational surface
c
c   sgnisl(jm)    sign of b1xqhf(jm).  This is passed to sbrtn island
c                 through the sign of b1scal(jm) (ie through pb1scl(jm))
c
c   b1xedg(j)     b1uxi(j) at the edge of the plasma  xi=1.
c
c   b1tedg(j)     b1ltr(j) at the edge of the plasma  xi=1.
c
c   b1zedg(j)     b1lzr(j) at the edge of the plasma  xi=1.

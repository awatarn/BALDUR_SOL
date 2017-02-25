c@clislb.m
c
      common /comisl/
     &      timisl,           nmodes
     & ,    modem(kpmode),    moden(kpmode),    qmode(kpmode)
     & ,    xmode(kpmode),    xinner(kpmode),   xouter(kpmode)
     & ,    xhalf(kpmode),    xdelta(kpmode),   dprime(kpmode)
     & ,    xhfmin,           xhfmax
     & ,    qb1amp(kpmode),   qb1scl(kpmode)
     & ,    qdeltx,           qxwall
     & ,    qpeakc(kpgrid),   qshldc(kpgrid),   qwexc
     & ,    nxiham,           nhamad,           nithta
     & ,    cislnd(32),       lislnd(32)
     & ,    equixi(kpgrid),   equir0(kpgrid)
     & ,    eqrmom(kpgrid,kpharm),   eqymom(kpgrid,kpharm)
c
c
c..This common block is used to interface
c     the tearing mode package with BALDUR:
c
c
c   timisl        time at which tearing mode island is computed [sec]
c                       initialize to less than tinit
c
c   nmodes        number of perturbed helical modes
c
c   modem(j)      mode number m for j=1,nmodes
c
c   moden(j)      mode number n for j=1,nmodes
c
c                 Note:  These must be in order
c                        modem(j)/moden(j) .lt. modem(j+1)/moden(j+1)
c
c   qmode(j)      = m / n  for each mode
c
c   xmode(j)      xi values at mode rational surfaces (output)
c
c   xinner(j)     xi values at inner edges of magnetic islands
c
c   xouter(j)     xi values at outer edges of magnetic islands
c
c   xhalf(j)      the width of the jth magnetic island in xi space
c
c   xdelta(j)     shift of x-point out from mode rational surface
c
c   dprime(j)     delta-prime for tearing mode
c
c   xhfmin        minimum island half-width allowed for computation
c
c   xhfmax        maximum island half-width allowed for computation
c
c   qb1amp(j)     the amplitute of the asymptotic expansion
c                     at xi = deltaxi.  bv(j)=b1amp(j)*xi**(modem(j)-1)
c
c   qb1scl(j)     scale factor applied to perturbed magnetic field
c                 determined from magnetic island width
c
c   qdeltx        small value of xi at which b1amp(j) is used (input)
c
c   qxwall        r_wall / r_plasma
c
c   qpeakc(j)     current peaking factor inside magnetic island
c                       same as \gamma in pf 29 (1987) 756
c
c   qshldc(j)     controls slope of current off island shoulders
c                       see sbrtn idiftr
c
c   qwexc         controls slope of current off island shoulders
c                       see sbrtn idiftr
c
c   nxiham        number of flux surfaces on which
c                       Hamada coordinate metric elements are computed
c                       = nxiequ at present
c
c   nhamad        number of hamada harmonics desired
c
c   nithta        number of angles or Fourier haromonics
c                       around each flux surface
c
c..the following should be placed in commhd
c
c   equixi(j)     xi on equilibrium moments grid  j=1,mflxs
c
c   equir0(j)     R_0 harmonics
c
c   eqrmom(j,m)   R_{j,m} harmonics               m=1,mhrms
c
c   eqymom(j,m)   Y_{j,m} harmonics
c
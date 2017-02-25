c@clintf.m
c
c     ------------------------------------------------------------------
cl    fluxsurface quantities interpolated from BALDUR-grid .... 86-05-23
c     ------------------------------------------------------------------
      common /comBFX/
     -        fxsib  , fxsiz  , fxp    , fxpp   , fxffp  , fxpsip ,
     -                          fxq    , fxiota ,
     &        eqprzi , eqiotb , eqxibi , eqxizi ,
     <        finbfx
      dimension  fxsib (Kjflx),fxsiz (Kjflx)
      dimension  fxp   (Kjflx),fxpp  (Kjflx),fxffp (Kjflx),fxpsip(Kjflx)
      dimension  fxq   (Kjflx),fxiota(Kjflx)
c     ------------------------------------------------------------------
      dimension  eqprzi(Kjflx),eqiotb(Kjflx),eqxibi(Kjflx),eqxizi(Kjflx)
c     ------------------------------------------------------------------
c
c     ------------------------------------------------------------------
cl    equilibrium and flux surface average -interface -- RGB -- 86-05-22
c     ------------------------------------------------------------------
      parameter  (Kpbt=20)
cl    parameter  (Kpbt=20)
c     ------------------------------------------------------------------
      common /comRGB/
     &        neqdim , mombnd , leqbnd , ntbk   , ntbkmx ,
     &        r0b    , rmb    , ymb    , r0bt   , rmbt   , ymbt   ,
     &        rmajbt , rminbt , elngbt , trngbt ,
     &        rmajb  , rminb  , elongb , tringb , eqtime ,
     &        tbk    , neqdt  , curmat , btort  , rtort  , rbtort ,
     <        finrgb
c
      dimension  rmb(Kmhrm),ymb(Kmhrm)
      dimension  r0bt(Kpbt),rmbt(Kpbt,Kmhrm),ymbt(Kpbt,Kmhrm)
      dimension  tbk(Kpbt),neqdt(Kpbt),curmat(Kpbt),btort(Kpbt),
     &           rtort(Kpbt),rbtort(Kpbt)
     & , rmajbt(kpbt), rminbt(kpbt), elngbt(kpbt), trngbt(kpbt)
c     ------------------------------------------------------------------
      common /comb08/  elonga, vminit(6)
c
c neqdim = dimension of eqtmp, ceqtmp, deqtmp, eeqtmp
c mombnd = number of harmonics specified for boundary
c leqbnd = 0 (default) for circular cylinder equilibrium (1-D BALDUR)
c        = 1 for geometric representation of boundary shape
c            (in terms of rmajbt, rminbt, elngbt, trngbt)
c        = 2 for harmonic represention of boundary shape
c            (in terms of r0bt, rmbt, ymbt)
c ntbk   = current value of tbk(ntbk) .le. eqtime .lt. tbk(ntbk+1)
c ntbkmx = maximum number of breakpoints prescribed
c    Input boundary values at breakpoint times tbk(j):
c tbk(j)      = breakpoint times
c neqdt(j)    = number of equilibria computed before the next tbk(j+1)
c curmat(j)   = toroidal plasma current in megamperes
c btort(j)    = vacuum toroidal field in tesla measured at rtort(j)
c rtort(j)    = radius at which btort(j) is defined
c rbtort(j)   = value of R Btor at time tbk(j) ( no longer used)
c    Harmonic representation of boundary shape:
c r0bt(j)     = 0th harmonic
c rmbt(j,m)   = mth harmonic of R
c ymbt(j,m)   = mth harmonic of Y
c    Geometric representation of boundary shape:
c rmajbt(j)   = major radius
c rminbt(j)   = half-width
c elngbt(j)   = elongation
c trngbt(j)   = triangularity
c    Boundary values at the current time:
c eqtime = time at which the following values are computed (in seconds)
c r0b    = boundary value of the 0th harmonic of R(theta)
c rmb    = boundary value of the mth harmonic set in sbrtn mhd
c ymb    = boundary value of the mth harmonic of Y(theta)
c rmajb  = interpolated value of rmajbt(jt)
c rminb  = interpolated value of rminbt(jt)
c elongb = interpolated value of elngbt(jt)
c tringb = interpolated value of trngbt(jt)
c    Variables needed for the VMOMS package:
c elonga = initial elongation at the magnetic axis
c vminit(j), j=1,6 values used to initialize param(j), use defaults
c
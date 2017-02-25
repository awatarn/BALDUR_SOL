c@cmonte.m
cl                  c2.4     user interface
c     version  aes 03-mar-82  add nlpion
       common/comusr/
     r  den0  , den0in, denein, denion, e0in  , eneut , fluxin, fluxn ,
     r  fluxr , fract , outflx, rmass , rnu   , rplsma, rsplit, rsurf ,
     r  sflux , sigma , sigvbx, sne   , sni   , snvol , spe   , spi   ,
     r  spii  , tein  , tiin  , veff  , wmin  ,
     i  maxgas, maxrad, maxsrc, maxsur, ngases, ngasi , npts  , nrandm,
     i  nsrces, nsrci , nsurf , nvsrc ,
     l  nlerr , nlmod , nlmono, nlpion, nlscat, nlsput
       logical
     l  nlerr , nlmod , nlmono, nlpion, nlscat, nlsput
       dimension
     r   den0(2,25,6),       den0in(2,6),        denein(25),
     r   denion(2,25),       e0in(6),  eneut(2,25,6),
     r   fluxin(6),          fluxn(6), fluxr(2,6),
     r   fract(2), outflx(2,6),        rmass(2), rnu(10),
     r   rsplit(10),         rsurf(25),          sflux(6),
     r   sigma(2,25,6),      sigvbx(25),         sne(2,25,6),
     r   sni(2,25,6),        snvol(2,25),        spe(25,6),
     r   spi(25,6),          spii(25,6),         tein(25),
     r   tiin(25), veff(2,6),
     i   ngasi(6), npts(6),  nsrci(3,2),         nvsrc(6),
     l   nlmono(6)
c
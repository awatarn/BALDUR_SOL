c@cltprf for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cprof1/
     &            exzoni(kpgrid),   ej7r(kpgrid),     epres(kpgrid)
     & ,          erbphi(kpgrid),   e17r2(kpgrid),    exbndi(kpgrid)
     & ,          eqprof(kpgrid),   eqxisl(kpgrid)
     & , eqrci0(kpgrid), eqrcim(kpgrid,kpharm), eqrsim(kpgrid,kpharm)
     & , eqyci0(kpgrid), eqycim(kpgrid,kpharm), eqysim(kpgrid,kpharm)
c
c  exzoni(j)      xi on BALDUR zone centers ( = xzoni(j) in BALDUR )
c  ej7r(j)        <J_\phi / R^2> on exzoni(j) grid
c  epres(j)       total plasma pressure on exzoni(j) grid
c  erbphi(j)      R B_\phi on exzoni(j) grid
c  e17r2(j)       < 1 / R^2 > on exzoni(j) grid
c
c  exbndi(j)      xi on BALDUR zone boundaries ( = xbouni(j) in BALDUR )
c  eqprof(j)      magnetic q value on exbndi(j) grid
c
c  eqxisl(jx)     xi on equilibrium moments grid from jx=1 (at mag axis)
c                   to jx = nxiequ (at edge of plasma)
c
c    The equilibrium flux surface shapes are given in harmonic form:
c  R(eqxisl(jx),theta) = eqrci0(jx) + sum_{jm=1}^{nmomeq}
c    [ eqrcim(jx,jm) * cos(jm*theta) + eqrsim(jx,jm) * sin(jm*theta) ]
c  Y(eqxisl(jx),theta) = eqyci0(jx) + sum_{jm=1}^{nmomeq}
c    [ eqycim(jx,jm) * cos(jm*theta) + eqysim(jx,jm) * sin(jm*theta) ]
c
c    This common block is intended only for use in the driver program
c  and its sbrtns.

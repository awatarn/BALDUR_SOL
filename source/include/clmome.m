c@clmome for the island package DISLAND2, Bateman, PPPL
c rgb 23-nov-91 changed ceqrcm((kpgrid,3,kpharm) to ceqrcm(kpgrid*5,kpharm)
c
      common /cmomeq/
     &   nxiequ,         nmomeq,               eqxibi(kpgrid)
     & , rc0bc, rcmbc(kpharm), rsmbc(kpharm), shiftr
     & , yc0bc, ycmbc(kpharm), ysmbc(kpharm), shifty
     & , eqrcj0(kpgrid), eqrcjm(kpgrid,kpharm), eqrsjm(kpgrid,kpharm)
     & , eqycj0(kpgrid), eqycjm(kpgrid,kpharm), eqysjm(kpgrid,kpharm)
     & ,ceqrc0(kpgrid*5),ceqrcm(kpgrid*5,kpharm),ceqysm(kpgrid*5,kpharm)
     & , drc0dx(kpgrid), drcmdx(kpgrid,kpharm), drsmdx(kpgrid,kpharm)
     & , dyc0dx(kpgrid), dycmdx(kpgrid,kpharm), dysmdx(kpgrid,kpharm)
c
c     The equilibrium flux surfaces are represented by
c
c  R(xi,theta) = eqrcj0(j) 
c                + sum over m { eqrcjm(j,m) * cos ( m * theta ) }
c                + sum over m { eqrsjm(j,m) * sin ( m * theta ) }
c  Y(xi,theta) = eqycj0(j)
c                + sum over m { eqycjm(j,m) * cos ( m * theta ) }
c                + sum over m { eqysjm(j,m) * sin ( m * theta ) }
c
c     for xi = eqxibi(j)
c
c     The independent variable xi is taken to be the same as
c  xbouni(j) on BALDUR zone bndaries and xzoni(j) on BALDUR zone centers.
c
c  nxiequ         = number of flux surfaces                 (input)
c                       on which the equilibrium moments are prescribed
c
c  nmomeq         = number of equilibrium moments
c
c  eqxibi(j)      = xi values for those flux surfaces
c                       on which the equilibrium moments are prescribed
c
c  rc0bc, rcmbc(m), and ysmbc(m) are the outer boundary values
c
c    ceqr0, ceqrcm, and ceqysm are the cubic spline coeffs 
c  corresponding to eqrcj0, eqrcjm, and eqysjm.
c    Note:  These are computed using sbrtn cubint, so they may by offset
c  by four elements to ensure proper symmetry at r=0.
c
c  drc0dx(j) = d eqrcj0(j  ) / d xi    for xi = eqxibi(j)
c  drcmdx(j) = d eqrcjm(j,m) / d xi
c  drsmdx(j) = d eqrsjm(j,m) / d xi
c  dyc0dx(j) = d eqycj0(j  ) / d xi
c  dycmdx(j) = d eqycjm(j,m) / d xi
c  dysmdx(j) = d eqysjm(j,m) / d xi

c@cd3he 15:00 12-feb-92 from Sugiyama
c
c     d-3he fusion
c
      parameter(nsv=12,n9=55)
c     common/sigval/sx(nsv),syddt(nsv),sbddt(nsv),
c    a  scddt(nsv),sdddt(nsv),
c    b  sydd3(nsv),sbdd3(nsv),scdd3(nsv),sddd3(nsv),
c    c  sydt(nsv),sbdt(nsv),scdt(nsv),sddt(nsv),
c    d  syd3(nsv),sbd3(nsv),scd3(nsv),sdd3(nsv),
c    e  sytt4(nsv),sbtt4(nsv),sctt4(nsv),sdtt4(nsv),
c    f  sxb(7),sydtb(6),sbdtb(6),scdtb(6),sddtb(6),syttb(7),sbttb(7),
c    g  scttb(7),sdttb(7),syd3b(6),sbd3b(6),scd3b(6),sdd3b(6)
c
      common/testf/tppd(n9),tpp3(n9),tpt(n9),tp3(n9),tp4t(n9),
     1 tp43(n9),tp42t(n9),tepd(n9),tep3(n9),tet(n9),te3(n9),te4t(n9),
     2 te43(n9),te42t(n9),fpde(n9),fp3e(n9),fte(n9),f3e(n9),f4te(n9),
     3 f43e(n9),f42te(n9),fddt(n9),fdd3(n9),fdt(n9),fdtb(n9),fd3(n9),
     4 fd3b(n9),ftt4(n9),fttb(n9),
     r   epde(n9),ep3e(n9),etre(n9),e3e(n9),e4te(n9),e43e(n9),
     r   e42te(n9),epdi(n9),ep3i(n9),etri(n9),e3i(n9),e4ti(n9),
     r   e43i(n9),e42ti(n9),
     r   drnp(n9),drnd(n9),drnt(n9),drn3(n9),drn4(n9)
        common/d3fus/q1fus(n9),q2fus(n9),q1lfus(n9),q2lfus(n9),
     r   qefus(n9),wntn(n9),d3fast(n9)
       common/d3fstp/fnpd0(n9),fnp30(n9),fnt0(n9),fn30(n9),fn4t0(n9),
     1 fn430(n9),fn42t0(n9),fnpd(n9),fnp3(n9),fnt(n9),fn3(n9),
     2 fn4t(n9),fn43(n9),fn42t(n9),rh1fst(2,n9),rh2fst(2,n9),
     3 epd0(n9),ep30(n9),et0(n9),e30(n9),e4t0(n9),e430(n9),e42t0(n9),
     4 epd(n9),ep3(n9),et(n9),e3(n9),e4t(n9),e43(n9),e42t(n9)

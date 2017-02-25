c@.../baldur/code/com/cmonspl.m  Common block for Monte Carlo package.
cl                  c2.6     splitting variables
c     version 2p  mhh  25/8/76  princeton
       common/comspl/
     r  rnumb , vl    , vxl   , vyl   , vzl   , wl    , xl    , yl    ,
     r  zl    ,
     i  maxlev, mcell , msurfl, ngasl , nlevel, nodes , nu    ,
     l  nlsplt, nlsurf
       logical
     l  nlsplt, nlsurf
       dimension
     r   rnumb(25),          vl(15),   vxl(15),  vyl(15),
     r   vzl(15),  wl(15),   xl(15),   yl(15),   zl(15),
     i   mcell(15),          msurfl(15),         ngasl(15),
     i   nodes(15),
     l   nlsplt(25)

c@cfokkr.m
cl                  c7.1     comfok--fokker-plank neutral beam
c
       common/comfok/
     r  bbdd,   bbddrs, bbddsv, bbdt,   bbdtrs, bbdtsv,
     r  chscat, czgamm, cznue , halsps, hchexs, hdei  , hdeps , hdmu3 ,
     r  hdmub , hdmuc , hecmps, heelec, hei   , heinjs, heion , heloss,
     r  heplas, hesrci, hfi   , hlogbs, hmub  , hmuc  , hncmps, hninjs,
     r  hnloss, hnplas, hnsrcs, hscats, hsigv , hsin  , hslows, htbi  ,
     r  htlast, htspon, websps, wibsps,
     i  lhemax, lhemin, libeam, mhsp  , mxhe  , mxhmu , mxhsp , mxhsrc,
     i  nbbsp1, nbbsp2, nbbspd,
     i  nhbem1, nilast
       real ::
     r   bbdtsv(10,10,20,20),bbddsv(10,10,20,20),
     r   halsps(55,2),       hchexs(20,55,2),    hdei(20,2),
     r   hdeps(10,10,2),     hdmu3(10),          hdmub(11),
     r   hdmuc(10),          hecmps(55),         heelec(20,55,2),
     r   hei(20,2),          heinjs(55),         heion(20,55,2),
     r   heloss(55),         heplas(55),         hesrci(10,2),
     r   hfi(10,20,55,2),    hlogbs(20,55),      hmub(11),
     r   hmuc(10),           hncmps(55),         hninjs(55),
     r   hnloss(55),         hnplas(55),         hnsrcs(10,2),
     r   hscats(20,55,2),    hsigv(20,2),        hsin(11),
     r   hslows(20,55,2),    htspon(2),          websps(55,2),
     r   wibsps(55,2)
       integer ::
     i   lhemin(55,2),       libeam(12),         nhbem1(2)
c
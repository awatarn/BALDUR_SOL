c@come04.m
c
      parameter (nxi=64,nyi=64,md=nxi+1,nd=nyi+1)
      parameter (im=24,jm=43,imxp=2*im+1,mmxp=7,mmmp=2*mmxp+1)
c
      common  /come04/                                     morder ,
     1        j      , jmax   , imax   , mmax   , ndim   , m , n  ,
     2        dx     , dy     , pds(6) ,  xma   , yma    , psi    ,
     3        x      , y      , xtp    , ytp    , ch     , cff    ,
     a                          dxtp   , dytp   ,
     4        xftp   , yftp   , cfcs   , cfsn   , hcx    , W2xy   ,
     5        l1     , l2     , l3     , l4     ,          M2xy   ,
     a        cfl    , cfc    , cfs    ,
     b        xbc    , ybc    , cbc    , wk     , lflgbc ,
     6        eps    , lp     , lq     , zmm    , zm2    , zdm2
c     .............................................................
      REAL         M2xy
      dimension    x(md),y(nd)
      dimension    xtp (imxp,jm),ytp (imxp,jm)
      dimension    dxtp(imxp,jm),dytp(imxp,jm)
      dimension    ch(jm),cff(mmxp)
      dimension    xftp(imxp),yftp(imxp)
      dimension    cfcs(mmxp,imxp),cfsn(mmxp,imxp),hcx(mmmp)
      dimension    w2xy(mmmp,mmmp),zdm2(mmmp),M2xy(mmmp,mmmp)
      dimension    cfl(mmmp),cfc(mmmp),cfs(mmmp)
      dimension    xbc(md),ybc(nd)
      dimension    cbc(2,md,2,md),wk(2*md*nd+2*md)
      dimension    psi(md,nd)
c
c@come14.m
c
      parameter (na=25,nb=21)
c
      common  /come14/
     1        ndima  , ncl    , nin    , nout   , ndsk   , iflux  ,
     2        cutoff , cutfak , psmxx  , psimin , psimax , psilim ,
     3        x0     , y0     , xa     , yb     , yt     ,
     4        pfx    , pfy    , pfxy   , theta  ,
     6        cxorg  , cx     , cy     , cxyst  , suavrg
c     .............................................................
      dimension    pf(na),pfx(na),pfy(na),pfxy(na,na)
      dimension    theta(imxp)
      dimension    cxorg(jm),cx(mmxp,jm),cy(mmxp,jm),cxyst(mmmp,jm)
      dimension    suavan(jm,16),suavrg(jm,16)
c     -------------------------------------------------------------

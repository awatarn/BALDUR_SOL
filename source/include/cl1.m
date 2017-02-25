c@cl1.m
c     ------------------------------------------------------------------
cl    local commons for control of equilibr. & averages........ 86-04-22
c     ------------------------------------------------------------------
      parameter  (Knf=32,Knfh=1+Knf/2,Kpw=(3*Knf/2)+2)
cl    parameter  (Knf=32,Knfh=1+Knf/2,Kpw=(3*Knf/2)+2)
c     ------------------------------------------------------------------
c
      common /c1/   lsw(8), neqout,ntty, jr, mom
     &      , mdim, jdim, lprt, nout
     &      , ftol, relerr, gam, gamma, itmom, nskip
     &      , initin, raxin, deltin, shift
c
      common /c2/
     &        xi(Kjflx), xiz(Kjflx)
     &      ,  r0(Kjflx),    rm(Kmhrm,Kjflx),    ym(Kmhrm,Kjflx)
     &      ,                rj(Kjflx,Kmhrm),    yj(Kjflx,Kmhrm)
     &      , cr0(5*Kjflx), crj(5*Kjflx,Kmhrm), cyj(5*Kjflx,Kmhrm)
c
      common /cave/
     &        dr0dxi(Kjflx), drmdxi(Kjflx,Kmhrm), dymdxi(Kjflx,Kmhrm)
     &      , r (Knf,Kjflx), y (Knf,Kjflx)
     &      , drdth(Knf), dydth(Knf)
     &      , drdxi(Knf), dydxi(Knf), aj(Knf)
     &      , rgrmin, rgrmax, zgrmin, zgrmax, pgrmin, pgrmax
     &      , tgrmin, tgrmax, cgrmin, cgrmax
c
      common /compf/ rf,       yf,       work
      complex        rf(Knfh), yf(Knfh), work(Kpw)
c
      common /cgeom/
     &               rin (Kjflx), yin (Kjflx)
     &             , rind(Kjflx), yind(Kjflx)
     &             , rtop(Kjflx), ytop(Kjflx)
     &             , rbot(Kjflx), ybot(Kjflx)
     &             , rout(Kjflx), yout(Kjflx)
c     ------------------------------------------------------------------
c
c     ------------------------------------------------------------------
cl    fluxsurface average quantities def. on flux-grid ........ 86-04-24
c     ------------------------------------------------------------------
      common /comFXA/
     -        rho    , dvdxi  , avdxi  , avdxi2 , adhir2 ,
     -        rbtorz , avir   , avir2  , vol    , area   , avjac  ,
     -        rbtorb , df2dxi ,
     -        pb     , pz     , dpdxi  , eiota  , cpz    ,
     -        xipsi  , xierr  , dsidxi ,
     -        rhoz   , dvdrhz , drhdxi , xia    ,
     -        ctemp  , temp   , dtemp  ,
     -        beta   , torcur , psitor ,
     >        ngeom  , ngskip ,
     -        xig    , qrmid  , qahalf , qelong , qtrian , qdent  ,
     -        rminar , rmajgc ,
     <        finfxa
      dimension  rho   (Kjflx),dvdxi (Kjflx),avdxi (Kjflx),avdxi2(Kjflx)
      dimension  adhir2(Kjflx),rbtorz(Kjflx),avir  (kjflx),avir2 (Kjflx)
      dimension  vol   (Kjflx),area  (Kjflx),avjac (Kjflx)
      dimension  rbtorb(Kjflx),df2dxi(Kjflx)
      dimension  pb    (Kjflx),pz    (Kjflx),dpdxi (Kjflx),eiota (Kjflx)
      dimension  cpz   ((Kjflx+4)*5)
      dimension  xipsi (Kjflx),xierr (Kjflx),dsidxi(Kjflx)
      dimension  ctemp ((Kjflx+4)*5)        ,temp  (Kjflx),dtemp (Kjflx)
      dimension  beta  (Kjflx),torcur(Kjflx),psitor(Kjflx)
      dimension  xig   (Kjflx),qrmid (Kjflx),qahalf(Kjflx),qelong(Kjflx)
      dimension  rhoz  (Kjflx),dvdrhz(Kjflx),drhdxi(Kjflx),xia   (Kjflx)
      dimension  qtrian(Kjflx),qdent (Kjflx),rminar(Kjflx),rmajgc(Kjflx)
c


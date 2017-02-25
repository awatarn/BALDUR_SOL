c@csynrd
c  cliche csynrd: tamor's synchrotron radiation internal block
c     and input (comsync)
c rgb 2-jun-92 added teqp, deny, voq, ploss for the original routines
c   these have units different from stemp, sdens, svol, wsyn which
c   are used in Sugiyama's routines
c
      common /comsyn/ree,roo,reo,roe,bmin,bmax,pwall,ptot,areaw
     &  ,taucrt
     &  ,stemp(55),sdens(55),bavg(55),svol(55),surf(55),wsyn(55)
     &  ,teqp(55), deny(55), voq(55), ploss(55)
     &  ,nzons,nw1,ndbug,iedit

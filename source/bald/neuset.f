c@neuset  /11040/baldur/code/bald/dimprad.f
c  rgb 31-jan-95 dimension pfutz(1) -> dimension pfutz(*)
c  dps 15-aug-88 Subroutine formed using statements formerly in IMPRAD.
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
      subroutine neuset(pfutz,peps,pepsin,k)
c
        dimension pfutz(*)
c
        data    iheflx,icfout,idvout,inpath /70,71,72,73/,
     1          ionlos,islant,itmhe0,istart /74,75,76,77/
        data  impneu,nergy1,nergy2,izlos1,izlos2 /200,201,202,203,204/,
     1        ixlos1,ixlos2,istar1,istar2,inpth1 /205,206,207,208,209/,
     2        inpth2,islnt1,islnt2,iajust,iprint /210,211,212,213,214/,
     3        ifact1,ifact2 /215,216/
        data  ileve1,ileve2,ixflx1,ixflx2,minti1 /220,221,222,223,224/,
     1        ihnce1,ihnce2,maxti1,ilower,iraise /225,226,227,228,229/,
     2        itime2,minti2,maxti2,itime3,ichng1 /260,261,262,270,271/
c
        if (k.eq.2) go to 10
c
c  Parameters for neutral impurity influx
c
        if(pfutz(impneu).le.peps) pfutz(impneu)=0.0
        if(pfutz(nergy1).le.peps) pfutz(nergy1)=.01
        if(pfutz(nergy2).le.peps) pfutz(nergy2)=.01
        if(pfutz(izlos1).le.peps) pfutz(izlos1)=.0136
        if(pfutz(izlos2).le.peps) pfutz(izlos2)=.0136
        if(pfutz(istar1).le.peps) pfutz(istar1)=1.0
        if(pfutz(istar2).le.peps) pfutz(istar2)=1.0
        if(pfutz(inpth1).le.peps) pfutz(inpth1)=0.9
        if(pfutz(inpth2).le.peps) pfutz(inpth2)=0.9
        if(pfutz(islnt1).le.peps) pfutz(islnt1)=1.0
        if(pfutz(islnt1).lt.0.5   ) pfutz(islnt1)=0.5
        if(pfutz(islnt2).le.peps) pfutz(islnt2)=1.0
        if(pfutz(islnt2).lt.0.5   ) pfutz(islnt2)=0.5
        if(pfutz(iajust).le.peps) pfutz(iajust)=0.0
        if(pfutz(ifact1).le.peps) pfutz(ifact1)=1.0
        if(pfutz(ifact2).le.peps) pfutz(ifact2)=1.0
        if(pfutz(iprint).le.peps) pfutz(iprint)=0.0
        if(pfutz(ileve1).le.peps) pfutz(ileve1)=1.0
        if(pfutz(ileve2).le.peps) pfutz(ileve2)=1.0
        if(pfutz(ixflx1).le.peps) pfutz(ixflx1)=1.e17
        if(pfutz(ixflx2).le.peps) pfutz(ixflx2)=1.e17
        if(pfutz(ihnce1).le.peps) pfutz(ihnce1)=10.
        if(pfutz(ihnce2).le.peps) pfutz(ihnce2)=10.
        if(pfutz(minti1).le.peps) pfutz(minti1)=0.0
        if(pfutz(maxti1).le.peps) pfutz(maxti1)=pepsin
        pfutz(maxti1)=max(pfutz(maxti1),(1.5*pfutz(minti1)))
        if(pfutz(ilower).le.peps) pfutz(ilower)=0.90
        pfutz(ilower)=min(pfutz(ilower),1.0)
        pfutz(ilower)=max(pfutz(ilower),0.2)
        if(pfutz(iraise).le.peps) pfutz(iraise)=1.10
        pfutz(iraise)=max(pfutz(iraise),1.0)
        pfutz(iraise)=min(pfutz(iraise),5.0)
        if(pfutz(itime2).le.peps) pfutz(itime2)=pepsin
        if(pfutz(minti2).le.peps) pfutz(minti2)=pfutz(minti1)
        if(pfutz(maxti2).le.peps) pfutz(maxti2)=pfutz(maxti1)
c
      return
c
   10 continue
c
c  Neutral impurity recycling parameters
c
        if(pfutz(icfout).le.peps) pfutz(icfout)=.9
        if(pfutz(idvout).le.peps) pfutz(idvout)=.9
        if(pfutz(inpath).le.peps) pfutz(inpath)=.9
        if(pfutz(ionlos).le.peps) pfutz(ionlos)=0.1
        if(pfutz(islant).le.peps) pfutz(islant)=1.0
        if(pfutz(itmhe0).le.peps) pfutz(itmhe0)=0.04
        if(pfutz(istart).le.peps) pfutz(istart)=1.0
c
      return
      end

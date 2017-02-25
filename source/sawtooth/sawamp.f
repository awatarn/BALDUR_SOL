c@sawamp  .../baldur/code/bald/sawamp.f
c  rgb 09-jun-87  skipped over thermal d + d reactions if ldeut .le. 0
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
c               %%%%%%%%%%%%%
c               %
c               awamp          %
c               %
c               %%%%%%%%%%%%%
c
c
        subroutine sawamp(pte,pne,pdd,ppalf,palf,puthrm)
c
c
c       set input variables to the current values
c
c  pte    = peak electron temperature [keV]
c  pne    = line average density [cm**(-3)]
c  ppalf  = alpha power produced by fusion
c  puthrm = electron and ion thermal energy
c
c
      include 'cparm.m'
      include 'cbparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c       peak electron temperature and line average density
c
!cap
        zdtfus = 0.
!
        ztemax=0.
        zsum=0.
        do 150 jz=lcentr,ledge
        if(tes(2,jz).gt.ztemax) ztemax=tes(2,jz)
        zsum=zsum+rhoels(2,jz)*dxzoni(jz)
150     continue
        pte=useh*ztemax
        pne=zsum
c
c       neutron rates
c
        zdts=0.
        zdds = 0.0
c
c
        do 358 jz = lcentr, ledge
c
        if (ltrit.le.0) go to 320
c
c               (thermal d-t sigma-v, in cm**3/sec)
c
        zti = tis(2,jz) * evsinv * 0.001  ! ion temp [keV] zone center
        if(zti.lt.0.5.or.zti.gt.240.) go to 305
        if(zti.gt.80.0) go to 310
c
c       hively, nuclear fusion 17, 873 (1977)
c       (misprint in hively's table i; coeff of ti**3 must be negative)
c
        zsigv = exp(-22.711874 * zti**(-.27472414) +
     2  (-23.836412) + (-9.3925317e-02)*zti + 7.9944173e-04*zti**2 -
     3  3.1438516e-06*zti**3)
        go to 315
 305    continue
        zsigv=0.0
        go to 315
 310    continue
        zlnzti=log(zti)
        zsigv=-.187196*zlnzti*zlnzti+1.44664*zlnzti-37.4296
        zsigv=1.027*exp(zsigv)
 315    continue
c
        zdtfus = rhohs(ldeut,2,jz) * rhohs(ltrit,2,jz) * zsigv
c
  320   continue
c
c
c               thermal d + d  to  he3 + neutron  sigma-v bar,
c               in cm**3/sec
c
      if ( ldeut .le. 0 ) go to 358
c
        zti = tis(2,jz) * evsinv * 0.001
        if(zti.lt.0.5.or.zti.gt.240.) go to 330
        if(zti.gt.80.0) go to 335
        if(zti.gt.30.) go to 325
c
c       mikkelsen fit to the <sigma-v> from greene
c
        zsigv=exp(-35.509-15.184/(zti**0.37))
        go to 340
325     continue
c
c       hively, nuclear fusion 17, 873 (1977)
c
        zsigv = exp(-15.314707 * zti**(-.40568858) +
     2          (-35.903973) + .19465321*zti**(-1.5743113))
        go to 340
 330    continue
        zsigv=0.0
        go to 340
 335    continue
        zlnzti=log(zti)
        zsigv=-.160990*zlnzti*zlnzti+2.55872*zlnzti-46.5088
        zsigv=.90078*exp(zsigv)
 340    continue
c
        zddfus = 0.5*rhohs(ldeut,2,jz)**2 * zsigv
c
c
        zdts = zdts + dvoli(jz)*zdtfus
        zdds = zdds + dvoli(jz)*zddfus
c
c
  358   continue
c
        pdt = zdts * usid
        pdd = zdds * usid
c
        ppalf=pdt*efusei*uise*usep
c
c               next get total beta
c
        if(nadump(1).le.lcentr) then
        isep = mzones
        isepm1 = ledge
        else
        isep = nadump(1)
        isepm1 = isep - 1
        endif
c
        zelec=0.
        zion=0.
        zbeam=0.
        zalpha=0.
        do 576 jzz=lcentr,ledge
        zelec=zelec+rhoels(2,jzz)*tes(2,jzz)*dvoli(jzz)
        zion=zion+rhoins(2,jzz)*tis(2,jzz)*dvoli(jzz)
        zbeam=zbeam+rhobis(2,jzz)*hebems(jzz)*dvoli(jzz)
        zalpha=zalpha+alphai(jzz)*ealfai(jzz)*dvoli(jzz)
576     continue
c
        puthrm=usee*1.5*(zelec+zion)*usid
c
        z00 = usid / vols(mzones,1)
c
        zelec=z00*zelec
        zion=z00*zion
        zbeam=z00*zbeam
        zalpha=z00*zalpha
c
        znorm=2.*(8.*fcpi/bzs**2)/xbouni(isep)**2
        zelec=znorm*zelec
        zion=znorm*zion
        zbeam=znorm*zbeam*(2./3.)
        zalpha=zalpha*znorm*(2./3.)*uisd*uise
        ztbeta=zelec+zion+zbeam+zalpha
c
        palf=(ppalf/1.e6)/((zelec+zion)**2*(bzs/1.e4)**4
     &  * vols(isep,1) * 1.e-6)
c
        return
        end

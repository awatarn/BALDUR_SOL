c/dfusion.f 28-Jan-97 maintained by Glenn Bateman, Lehigh 
c [plasma.physics.lehigh.edu] /u2/baldur/bald/code1/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c..This file contains the following subroutines:
c
c  alfini - initialize common variables used by the alpha sbrtns
c  alphas - compute alpha particle sources and sinks
c  arlput - sets array arlong = time per zone for each sample alpha particle
c  squeze - repacks arlong array
c  alfend - messages afrom sbrtns aflini, alphas, arlput, squeze
c  alfset - initialies geometric and plasma profile arrays 
c              used by sbrtn bounce
c  xinter - 
c  bounce - alpha trajectory and collisions with background plasma
c  aphif  - function
c  endzon - function
c
c**********************************************************************c
c@aclear  .../bald/code/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine aclear
        return
        end
c
c**********************************************************************c
c@alfini  .../bald/code/bald/dfusion.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 18.86 11-feb-91  initialize atbi to tai.  cleaned up if statements
c  drm 4-dec-84 changed nilast to nalast, add fit to greene's d-d
c       neutron cross section
c  mhr 16-may-84 protected against dividing by zero by limiting
c               ipne to be .lt. nsampl
c  aes 11-feb-82 add afslow: multiplies slowdown rate. (default=1.)
c  aes 6-oct-81 moved adds and addtot lines to after label 358
c  aes 23-nov-80 add zecrh at do 1025
c  dhfz 23-jan-80 add d-t & d-d rate fits for 80<t<240.
c  dhfz 23-jan-80 correct sign to - at end of $305-2
c  dhfz 23-jan-80 delete all 'c7600' lines
c  dhfz 19-dec-79 date of present version
c**********************************************************************
        subroutine alfini
c
c  initialize common variables used by alpha sbrtns
c
c
       include 'cparm.m'
       include 'cbaldr.m'
       include 'calpht.m'
c
c               initialization
c
c
        aalpha = 4.00260
        azalfa = 2.0
        efusei = 3.52e+06 * evs * usie
        nalast = 0
        atbi   = tai
      do 19 jz=lcentr,ledge
        ealfai(jz)=efusei
19    continue
        efini=3520000.
        afast=4.
        qfast=2.
c
        fnoutr=0.01
        ftoutr=0.01
        lleror=.false.
        llslow=.true.
        llptch=.false.
        nlost=0
        nstop=0
        nstart=0
        narrow=0
        ndzerr=0
        mskip=0
        nsqeze=0
        nslost=0
        mxlong=17000
        mxshrt=500
        mxznes=75
        mxmnzn=1500
        npznes=20
        nsubzn=24
        ndrad=7
        ndang=4
        ndpit=10
        frcang=0.6
        frcpit=1.5
        nxsran=1
        nxcorl=1
        nxwght=1
        nxsmpl=1
        nlsmpl=1000
        espred=0.
        frcnew=0.02
        nlnewo=.true.
c
      if ( (rcwals-rdwals.le.rmajs-rmins) .and.
     &       (rcwals+rdwals.ge.rmajs+rmins) ) go to 15
        rcwals=rmajs
        rdwals=rmins
  15  continue
c
      if ( (ntype .lt. 1 ) .or. (ntype .gt. 3 ) ) ntype=2
      if (ntype .eq. 1 ) then
        frslow=0.80
        elecut=efini
        if ( ledge .gt. mxznes ) call alfend(1)
      else
        if ( frslow .le. elecut/efini) ntype=2
        if ( ntype .eq. 2 ) then
          nlnewo = .false.
          frslow = 0.80
          elecut = 700000.
        else if ( ntype .eq. 3 ) then
          frslow = 0.80
          elecut = 700000.
        endif
      endif
c
        actvfe=0.
        actvfw=0.
        actvse=0.
        actvsw=0.
        llset=.false.
        nsampl=0
        nbounc=0
        npactv=0
        nldead=0
        ntperm=0
        nxper=mxshrt+1
        nxorb=0
        npdead=0
        nulong=1
c
      do 110 j=1,mxshrt
        nblong(j)=0
        nlbrth(j)=.false.
 110  continue
c
      do 115 jz=1,mxznes
        aeslow(jz)=0.
        awslow(jz)=0.
 115  continue
c
        return
        end
c
c**********************************************************************c
c@alphas  .../bald/code/bald/dfusion.f
c  rgb 18.86 11-feb-91 removed all calls to expert
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine alphas(k)
c
c
cl      2.19    compute alpha particle sources and sinks
c
c
        logical ilpass
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'calpht.m'
      include 'commhd.m'
c
c
        logical ilstep,ilongt
        dimension zcoefe(75),zcoefi(75),zcoefx(75),zawe(75),zawi(75),
     x  zbcofe(mj),zbcofi(mj),zbcofx(mj),zbifrc(mj),zbtaus(mj),
     x  zbwe(mj),zbwi(mj),zarns(75),zares(75),
     x  zasalf(75),ikmax(75),ikill(75)
 
c
c
        data
     i  incons,intot,inflxi,inflxo,inloss,insorc,incomp /1,2,3,4,5,6,7/,
     i  iecons,ietot,ieflxi,ieflxo,ieloss,iesorc,iecomp /1,2,3,4,5,6,7/
c
c
c------------------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /19/
c
c
        if (.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.19 subroutine alphas bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       comalf, comorb, wealfs, wialfs, salfs, shfus (comdif)
c
c
c------------------------------------------------------------------------------
c
c
c       ie----          indices in econsi
c       ilstep          .true. if equations need to be solved
c       in----          indices in acons
c       iz              backwards zone index
c
c       z2voli          volume of torus * 2
c       zbe             slowing down rate due to electrons
c       zbi             slowing down rate due to ions * v**3
c       zzbe            constant part of zbe
c       zzbi            constant part of zbi
c       zdti1           part of timestep for which slowing down, etc.,
c                       has already been computed. (0 unless there
c                       has been a t.s. repetition)
c       zdti2           part of timestep for which slowing down, etc.,
c                       is to be computed. (dti unless there has been
c                       a t.s. repetition)
c       zdt10           zdti1 / dti
c       zdt20           zdti2 / dti
c       zslmax          maximum val. of aslows
c       zddfus          d+d to he3+n reaction rate, reac./cc/sec
c       zdds            accumulator for adds integr. over torus
c
c
c
c------------------------------------------------------------------------------
c
c
c
c
        if (ldeut.le.0) return
        if (k.eq.2) go to 200
        if (k.gt.2) go to 7000
c
  100   continue
        call alfend(2)
        return
c
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
cl      2)      adjust for timestep repetition
c
c
c
  200   continue
c
        if (tbi.le.atbi) return
c
        zdti2 = max(0.0, tbi-atbi )
        zdti1 = max(0.0, dti-zdti2 )
        zdt20 = zdti2 / dti
        zdt10 = max(0.0, 1.0-zdt20 )
c
c
c               scale source / sink terms
c
c
        call scaler( wealfs ,mxzone,zdt10)
        call scaler( wialfs ,mxzone,zdt10)
        call scaler( alphai ,mxzone,zdt10)
        call scaler( ealfai ,mxzone,zdt10)
        call scaler( salfs  ,mxzone,zdt10)
        call scaler( shfus  ,mxzone*mxhyd,zdt10)
c
c
        atbi = tbi
c
c
c
cl      3)      compute diffusion coefficient, slowing down rate,
c               thermal d-t reaction rate.  return if no alphas and no
c               fusion
c
c
c
  300   continue
c
c
        zdds = 0.0
c
c
c
        do 358 jz = lcentr, mzones
c
        if (jz.gt.ledge) go to 357
c
c
c
c
c               thermal d-t reaction rate
c
c
        if (ltrit.le.0) go to 320
c
c
c
c
c               (thermal d-t sigma-v, in cm**3/sec)
c
c
        zti = tis(2,jz) * evsinv * 0.001
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
        afuses(jz) = rhohs(ldeut,2,jz) * rhohs(ltrit,2,jz) * zsigv
        shfus(ldeut,jz) = shfus(ldeut,jz)  + zdt20*afuses(jz)
        shfus(ltrit,jz) = shfus(ltrit,jz)  + zdt20*afuses(jz)
c
  320   continue
c
c
c               thermal d + d  to  he3 + neutron  sigma-v bar,
c               in cm**3/sec
c
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
        shfus(ldeut,jz) = shfus(ldeut,jz) + zddfus*2.0
        zdds = zdds + dx2i(jz)*zddfus
c
c
c
c
 357    continue
  358   continue
        adds = zdds * 2.0 * vols(mzones,1)
        addtot = addtot + adds * zdti2*uist
c
c
c
c       check for alpha production
c
c
        do 370 jz=lcentr,ledge
        if(halfas(jz)+afuses(jz).gt.epslon) go to 375
370     continue
        go to 5000
c
c
 375    continue
c
        efast=efini+abs(espred)
        call alfset
        acons(2)=1.1
c
c
c       rgb 14-nov-85  the following statements may need to be
c        generalized to noncircular geometry
c
c
        zgamma=fcc*afast*fcmp*(10.**(fxc+fxnucl-fxes))/(fces*
     x  qfast*bpols(1,ledge+1)*rmins)
        zslrm=fcc*afast*fcmp*(10.**(fxc+fxnucl-fxes))/(fces*
     x  qfast*bzs*rmins)
        zv0drf=fcc*afast*fcmp*(10.**(fxc+fxnucl-fxes))/(fces*
     x  qfast*bzs*rmajs*2.)
        zsped=sqrt(2.*evs/(fcmp*(10.**fxnucl)*afast))
c
c
        if(nalast.eq.0) go to 380
        if (nalast+niprof.gt.nstep) go to 470
 380    continue
c
c%%%%%%%%%
c                                       begin aoloss section
c
 400    nalast = nstep
c
c       set up for bounce calls
c
        efast=efini
        speed=zsped*sqrt(efast)
        gamma=zgamma*speed
        slarmp=zslrm*speed
        v0drft=zv0drf*speed**2
        iquick=1
        ibounc=1
c
c       begin zone loop
c
        ie=ledge+1-lcentr
        do 430 jsudo=1,ie
        jz=ledge+1-jsudo
        zbsmn=xzoni(jz)
        iqang=8
        if(zbsmn.lt.0.8) iqang=6
        if(zbsmn.lt.0.5) iqang=4
        if(zbsmn.lt.0.2) iqang=4
        iqpit=96/iqang
        zdang=3.14159/float(iqang)
        zdcpit=2./float(iqpit)
        ilost=0
c
c       begin poloidal angle loop
c
        do 420 jang=1,iqang
        zcosa=cos((float(jang)-0.5)*zdang)
        zsina=sin((float(jang)-0.5)*zdang)
c
c       begin pitch angle loop
c
        do 410 jpit=1,iqpit
        sxion=zbsmn*zcosa
        szion=zbsmn*zsina
        cospit=(float(jpit)-0.5)*zdcpit-1.
        if(espred.eq.0.) go to 405
        efast=efini+espred*(2.*ranz()-1.)
        speed=zsped*sqrt(efast)
        gamma=zgamma*speed
        slarmp=zslrm*speed
        v0drft=zv0drf*speed**2
405     continue
        call bounce(iquick,ibounc)
        if(mlost) ilost=ilost+1
 410    continue
 420    continue
        aoloss(jz)=float(ilost)/float(iqang*iqpit)
        if(ilost.eq.0) go to 440
 430    continue
        go to 460
c
c       set aoloss to zero in inner zones
c
 440    continue
        do 450 j=lcentr,jz
        aoloss(j)=0.
 450    continue
 460    continue
  470   continue
c
c
        do 480 jz = 1, ledge
        iz = 1 + ledge - jz
c
c
c
        if (iz.le.1) go to 475
c
c
c
c
c
        acons(inloss) = acons(inloss)
        acons(insorc) = acons(insorc)
        acons(intot) = acons(intot)
c
  475   continue
  480   continue
c
c
c               compute influxes at outer and inner boundaries
c
c
        acons(inflxo) = acons(inflxo)
c
c
        acons(inflxi) = acons(inflxi)
c
c
c
  500   continue
        do 578 jz = 1, ledge
        iz = 1 + ledge - jz
        if (iz.le.1) go to 568
c
c
        econsi(ieloss) = econsi(ieloss)
        econsi(iesorc) = econsi(iesorc)
        econsi(ietot) = econsi(ietot)
c
        if (iz.ne.ledge) go to 577
c
  568   continue
c
c
c
c
c
        z1=1.
c
c
        if (iz.eq.ledge) econsi(ieflxo) = econsi(ieflxo) - z1
        if (iz.eq.1) econsi(ieflxi) = econsi(ieflxi) + z1
 577    continue
  578   continue
c
c
c       6)      set cumulative totals
c
c
c
  600   continue
c
c
        z0 = 2.0 * vols(mzones,1) * uist*zdti2
c
        do 608 jz = lcentr, ledge
        atfuse(jz) = atfuse(jz) + afuses(jz)*z0*dx2i(jz)
        abfuse(jz) = abfuse(jz) + halfas(jz)*z0*dx2i(jz)
  608   continue
c
c
c%%%%%%%%%
c                                       begin orbit following
c
        if(ntype.eq.1) go to 2000
c
c
c       integrate the alpha birth rates in time and space
c
        zzsrv=0.
        zzbir=0.
        zvol = 2.0 * vols(mzones,1)
        do 1010 jz=lcentr,ledge
        zbir=(halfas(jz)+afuses(jz))*dx2i(jz)*zvol*zdti2*uist
        biralf(jz)=biralf(jz)+zbir
        zzbir=zzbir+zbir
        zzsrv=zzsrv+zbir*(1.-aoloss(jz))
 1010   continue
        totbir=totbir+zzbir
        totsrv=totsrv+zzsrv
c
c       set ntsmpl
c
        mxsmpl=mxshrt
        if(npactv.eq.0) go to 1020
        ii=0
        ipactv=0
        do 1015 jp=1,nsampl
        if(nblong(jp).eq.0) go to 1014
        ii=ii+2*nzlong(jp)+2
        ipactv=ipactv+1
1014    continue
1015    continue
        mxsmpl=min(mxshrt,int(float(ipactv)*float(mxlong)
     x  *0.80/float(ii)))
1020    continue
        zohm=0.
        zbem=0.
        zaux=0.
        zecrh=0.
        zalf=0.
        do 1025 jz=lcentr,ledge
        zohm=zohm+weohms(jz)*dx2i(jz)
        zbem=zbem+(webems(jz)+wibems(jz))*dx2i(jz)
        zaux=zaux+(weauxs(jz)+wiauxs(jz))*dx2i(jz)
        zecrh=zecrh+(weecrh(jz)+wiecrh(jz))*dx2i(jz)
        zalf=zalf+(halfas(jz)+afuses(jz))*dx2i(jz)
1025    continue
        zalf=zalf*efini*evs
        zz=zalf/max(zohm,zbem,zaux,zecrh,zalf)
        ntsmpl=min(mxsmpl,int(float(mxsmpl)*4.*zz))
        ntsmpl=max0(50,ntsmpl)
        if(ntsmpl.le.nsampl) go to 1029
        ie=ntsmpl-nsampl
        do 1027 j=1,ie
        npdead=npdead+1
        nadead(npdead)=nsampl+j
1027    continue
        if(npactv.eq.0) nsampl=ntsmpl
1029    continue
        nsampl=max0(ntsmpl,nsampl)
        nsampl=min(mxshrt,nsampl)
c
c       how many particles should be started ?
c
        if(actvfw.ne.0.) go to 1030
        ipnew=nsampl
        go to 1035
1030    continue
        zz=(totsrv*efini)/(totsrv*efini+actvfe)
        if(nxwght.eq.1) zz=totsrv/(totsrv+actvfw)
        ii=ntsmpl-3*int(sqrt(float(nldead)))
        ipnew=min(ntsmpl,int(float(ii)*zz+ranz()))
        ipnew=min(nsampl-1,ipnew)
1035    continue
        igtnew=int(frcnew*float(ntsmpl))
        if(ipnew.lt.igtnew) go to 1600
c
c
c
c       initialize for starting new sample particles
c
        zerlte=0.
        iptcle=0
        itdead=npdead
        ibirth=0
        inhalf=0
        infull=0
c
c               initialize for bounce calls
c
        speed=zsped*sqrt(efini)
        gamma=zgamma*speed
        slarmp=zslrm*speed
        v0drft=zv0drf*speed**2
        iquick=0
        ibounc=1
c
c
c               initialize for npermu(ntperm)
c
c
        if(ipnew.gt.ndrad*ndang*ndpit) go to 1045
1040    continue
        idrad=ndrad
        idang=ndang
        idpit=ndpit
        if(ntperm.eq.ndrad*ndang*ndpit) go to 1070
        ntperm=ndrad*ndang*ndpit
        nxper=ntperm+1
        go to 1070
1045    continue
        zz=(float(ipnew)/(frcang*frcpit))**0.3333
        idrad=int(zz+0.5)
        idang=int(zz*frcang+0.5)
        idpit=int(float(ipnew)/float(idrad*idang)+0.5)
1050    if(idrad*idang*idpit.le.ndrad*ndang*ndpit) go to 1040
        ntperm=idrad*idang*idpit
        if(ntperm.le.mxshrt) go to 1060
        idpit=idpit-1
        go to 1050
1060    continue
        nxper=ntperm+1
1070    continue
c
c       set up ikmax
c
        if(ipnew.le.npdead) go to 1095
        if(npactv.eq.0) go to 1095
        do 1075 jz=1,ntznes
        ikmax(jz)=0
        ikill(jz)=0
1075    continue
        ipactv=0
        do 1080 jp=1,nsampl
        if(nblong(jp).eq.0) go to 1079
        ii=nzbgin(jp)
        ikmax(ii)=ikmax(ii)+1
        ipactv=ipactv+1
1079    continue
1080    continue
        zz=float(ipnew-npdead)/float(ipactv)
        do 1090 jz=1,ntznes
        ikmax(jz)=int(zz*float(ikmax(jz)))+1
1090    continue
1095    continue
c
c       start loop for new particles
c
        do 1500 jstart=1,ipnew
c
c
        if(npdead.lt.1) go to 1100
c
c       use a dead particle number
c
        iptcle=nadead(npdead)
        npdead=npdead-1
        go to 1130
c
c       kill off an active particle
c
 1100   continue
        npdead=0
c
        if(npactv.ne.0) go to 1110
        iptcle=iptcle+1
        go to 1140
1110    iptcle=min(nsampl,int(float(nsampl)*ranz())+1)
        if(nlbrth(iptcle)) go to 1110
        ii=nzbgin(iptcle)
        if(ikill(ii).ge.ikmax(ii)) go to 1110
        ikill(ii)=ikill(ii)+1
        if(nblong(iptcle).ne.0) go to 1120
        call mesage('  alphas 4130:  nblong(iptcle)=0 ')
        go to 1130
1120    continue
        zerlte=zerlte+arwght(iptcle)*(areion(iptcle)-adeion(iptcle))
c
1130    continue
        ibirth=ibirth+1
        nlbrth(iptcle)=.true.
c
1140    continue
c
c
c
 1200   if(nxper.le.ntperm) go to 1300
c
c               load up npermu(j)
c
        nxper=1
        do 1210 j=1,ntperm
        npermu(j)=j
 1210   continue
        do 1220 j=1,ntperm
        iflip=j+int(ranz()*float(ntperm+1-j))
        iflip=min(iflip,ntperm)
        ii=npermu(j)
        npermu(j)=npermu(iflip)
        npermu(iflip)=ii
1220    continue
c
c       set birth position
c
1300    continue
        istart=npermu(nxper)
        nxper=nxper+1
        ii=idrad*idang*idpit
        if(nxcorl.eq.0) istart=min(ii,int(float(ii)*ranz()+1.))
        israd=(istart-1)/(idang*idpit)+1
        ii=istart-(israd-1)*(idang*idpit)
        isang=(ii-1)/idpit+1
        ispit=ii-(isang-1)*idpit
        jz=lcentr-1
        zbir=0.
        zz=0.5
        if(nxsran.eq.1) zz=ranz()
        ztest=totbir*(float(israd)-zz)/float(idrad)
1310    if(zbir.gt.ztest) go to 1330
        jz=jz+1
        zbir=zbir+biralf(jz)
        if(jz.gt.ledge) call alfend(3)
        go to 1310
1330    continue
        zbrad=xbouni(jz+1)-(zbir-ztest)*dxzoni(jz)/biralf(jz)
        if(jz.le.lcentr+1) zbrad=xbouni(jz+1)*sqrt(ztest/zbir)
        zz=0.5
        if(nxsran.eq.1) zz=ranz()
        zang=fcpi*(float(isang)-zz)/float(idang)
        zz=0.5
        if(nxsran.eq.1) zz=ranz()
        cospit=2.*(float(ispit)-zz)/float(idpit)-1.
        sxion=zbrad*cos(zang)
        szion=zbrad*sin(zang)
        if(espred.eq.0.) go to 1400
        efast=efini+espred*(2.*ranz()-1.)
        speed=zsped*sqrt(efast)
        gamma=zgamma*speed
        slarmp=zslrm*speed
        v0drft=zv0drf*speed**2
c
c
c
c
c       call bounce and set particle arrays
c
1400    continue
        call bounce(iquick,ibounc)
        if(mlost) nslost=nslost+1
        if(llskip) go to 1200
        arwght(iptcle)=totsrv/float(ipnew)
        arcnap(iptcle)=canap
        arstar(iptcle)=smj0
        areion(iptcle)=efast
        arsped(iptcle)=speed
        adcnap(iptcle)=0.
        adstar(iptcle)=0.
        adeion(iptcle)=0.
        call arlput(iptcle)
        if(lzinnr.ne.1) infull=infull+1
        if(lzinnr.eq.1) inhalf=inhalf+1
 1500   continue
c
c
c
c       put roulette weight back in old particles
c
        if(zerlte.eq.0) go to 1511
        zz=actvfe/(actvfe-zerlte)
        do 1510 jptcle=1,nsampl
        if(.not.nlbrth(jptcle)) arwght(jptcle)=zz*arwght(jptcle)
1510    continue
1511    continue
c
c       reset to accumulate new alpha production
c
        totbir=0.
        totsrv=0.
        do 1520 jz=lcentr,ledge
        biralf(jz)=0.
1520    continue
        llset=.true.
1600    continue
c
c%%%%%%%%%
c                               begin heating and orbit updates
c
c       set up for energy loss coeficients
c
2000    continue
        zcfecr=3.*sqrt(fcpi*1836)/4.
        zcflsi=4.*fcpi*(qfast**2)*((fces**4)/fcmp)*(10.**
     x  (fxes*4.-fxnucl))*afslow
        zcflse=(16./3.)*sqrt(fcpi)*(qfast)**2*((fces)**4/fcme)
     x  *(10.**(fxes*4.-fxme))*afslow
        iz=2
c
c       set up loss rate arrays in baldur zones
c
        do 2030 jz=lcentr,ledge
        zloge=(51.2-log(sqrt(rhoels(iz,jz))/tes(iz,jz)))
        zbcofe(jz)=rhoels(iz,jz)*zcflse*zloge
        zbcofx(jz)=1./sqrt(2.*tes(iz,jz)/(fcme*10.**fxme))
        zions=0.
        do 2010 jion=lhyd1,lhydn
        zions=zions+nzspec(jion)**2*rhohs(jion,2,jz)/aspec(jion)
 2010   continue
        if(mimp.eq.0) go to 2025
        do 2020 jion=1,mimp
        ii=limp1-1+jion
        zions=zions+c2mean(jion,2,jz)*rhois(jion,2,jz)/aspec(ii)
 2020   continue
2025    continue
        zlogi=(49.7-0.5*log(rhoels(iz,jz)/tes(iz,jz)))
        zbcofi(jz)=zcflsi*zions*zlogi
        zecrit=afast*(zcfecr*zions*zlogi/(zloge*rhoels(iz,jz)))**
     x  (2./3.)*tes(iz,jz)*evsinv
        zerat=elecut/zecrit
        zz=sqrt(zerat)
        zsq3=1./sqrt(3.)
        zbifrc(jz)=(2./zerat)*((atan((2.*zz-1.)*zsq3)+fcpi/6.)
     x  *zsq3-log(((1.+zz)**2)/(1.-zz+zz**2))/6.)
        zbtaus(jz)=log(1.+zerat**1.5)*afast*fcmp*(10.**fxnucl)
     x  /(3.*(zbcofx(jz)**3)*zbcofe(jz))
        if(ntype.eq.1) zbtaus(jz)=zbtaus(jz)*1.5/log(1.+zerat**1.5)
        zerat=efini/zecrit
        aslows(jz)=1./(log(1.+zerat**1.5)*afast*fcmp*(10.**fxnucl)
     x  /(3.*(zbcofx(jz)**3)*zbcofe(jz)))
 2030   continue
        if(ntype.ne.1) go to 2035
        do 3050 jz=lcentr,ledge
        ealfai(jz)=efusei
        ztaus=zbtaus(jz)
        zitaus=1./ztaus
        znstri=ztaus*(1.-aoloss(jz))*(afuses(jz)+halfas(jz))*usid
        zdni=awslow(jz)-znstri
        zaxexp=ztaus/(ztaus+zdti2*uist)
        zavexp=2.*ztaus/(2.*ztaus+zdti2*uist)
        zavgni=zdt20*(znstri+zdni*zavexp)
        alphai(jz)=alphai(jz)+zavgni
        salfs(jz)=salfs(jz)+zavgni*uisd*zitaus
        zwalfs=zavgni*uisd*ealfai(jz)*uise*zitaus
        wealfs(jz)=wealfs(jz)+(1.-zbifrc(jz))*zwalfs
        wialfs(jz)=wialfs(jz)+zbifrc(jz)*zwalfs
        awslow(jz)=znstri+zdni*zaxexp
3050    continue
        go to 5000
c
c       set up loss rate arrays in alphas zones
c
2035    continue
        ii=lcentr
        zdr=1./float(npznes)
        do 2080 jz=1,ntznes
        zsmnc=zdr*(float(jz)-0.5)
2040    if(zsmnc.le.xzoni(ii)) go to 2045
        if(ii.eq.ledge) go to 2060
        ii=min(ii+1,ledge)
        go to 2040
2045    continue
        if(ii.ne.lcentr) go to 2050
        zcoefe(jz)=zbcofe(lcentr)
        zcoefi(jz)=zbcofi(lcentr)
        zcoefx(jz)=zbcofx(lcentr)
        wifrac(jz)=zbifrc(lcentr)
        sltaus(jz)=zbtaus(lcentr)
        go to 2070
 2050   continue
        zz=(zsmnc-xzoni(ii-1))/(xzoni(ii)-xzoni(ii-1))
        zcoefe(jz)=(1.-zz)*zbcofe(ii-1)+zz*zbcofe(ii)
        zcoefi(jz)=(1.-zz)*zbcofi(ii-1)+zz*zbcofi(ii)
        zcoefx(jz)=(1.-zz)*zbcofx(ii-1)+zz*zbcofx(ii)
        wifrac(jz)=(1.-zz)*zbifrc(ii-1)+zz*zbifrc(ii)
        sltaus(jz)=(1.-zz)*zbtaus(ii-1)+zz*zbtaus(ii)
        go to 2070
2060    continue
        zcoefe(jz)=fnoutr*zbcofe(ledge)
        zcoefi(jz)=fnoutr*zbcofi(ledge)
        zcoefx(jz)=zbcofx(ledge)
        wifrac(jz)=zbifrc(ledge)
        sltaus(jz)=zbtaus(ledge)
2070    continue
2080    continue
c
c       zero zawe and zawi
c
        do 2110 jz=1,ntznes
        zawe(jz)=0.
        zawi(jz)=0.
        zarns(jz)=0.
        zares(jz)=0.
        zasalf(jz)=0.
2110    continue
c
c
c       increment zawi,zawe,zares,zarns,and zasalf
c               and advance aeslow and awslow in time
c
        actvse=0.
        actvsw=0.
        do 2120 jz=1,ntznes
        actvse=actvse+aeslow(jz)/evs
        actvsw=actvsw+awslow(jz)
        zares(jz)=zares(jz)+aeslow(jz)
        zarns(jz)=zarns(jz)+awslow(jz)
        zeslow=(1.-1./(1.+zdti2*uist/sltaus(jz)))*aeslow(jz)
        zawi(jz)=zawi(jz)+wifrac(jz)*zeslow/(zdti2*uist)
        zawe(jz)=zawe(jz)+(1.-wifrac(jz))*zeslow/(zdti2*uist)
        aeslow(jz)=aeslow(jz)-zeslow
        zwslow=(1.-1./(1.+zdti2*uist/sltaus(jz)))*awslow(jz)
        zasalf(jz)=zasalf(jz)+zwslow/(zdti2*uist)
        awslow(jz)=awslow(jz)-zwslow
 2120   continue
c
c
c       increment adeion,adcnap and update particle orbits
c
c               initialize
c
        iquick=0
        ibounc=2
        itjump=2
        actvfe=0.
        actvfw=0.
        npactv=0
        nldead=0
c
c               start particle loop for heating and orbit update
c
        do 2500 jptcle=1,nsampl
        zdt=zdti2
        if(nlbrth(jptcle)) zdt=0.5*zdti2
        ztleft=zdt
2210    continue
        if(zdt.le.0.) go to 2270
        iblong=nblong(jptcle)
        if(iblong.le.0) go to 2490
        if(arlong(iblong).gt.0.) go to 2220
2215    call mesage('  alphas:2215 arlong(iblong).le.0. ')
        go to 2490
2220    continue
        izbgin=nzbgin(jptcle)
        izlong=nzlong(jptcle)
        izstop=izbgin+izlong-1
        speed=arsped(jptcle)
        zwght=arwght(jptcle)
        ii=itjump+iblong-izbgin
c               first estimate of zelosr
        zelosr=0.
        do 2230 jz=izbgin,izstop
        zxe3=(speed*zcoefx(jz))**3
        zelosr=zelosr+arlong(ii+jz)*(zcoefi(jz)+zcoefe(jz)*zxe3
     x  /(1.+0.75*zxe3))/speed
2230    continue
        zelosr=zelosr*evsinv
c       ztest=min(areion(jptcle)-elecut*0.99,
c     x  areion(jptcle)*(1.-frslow))
        ztest=areion(jptcle)*(1.-frslow)
        zde=min(adeion(jptcle)+zelosr*zdt*uist,ztest)
        zasped=0.5*(speed+zsped*sqrt(areion(jptcle)-zde))
c               second estimate of zelosr
        zelosr=0.
        zdcnpr=0.
        do 2240 jz=izbgin,izstop
        zxe3=(zasped*zcoefx(jz))**3
        zde=(zcoefi(jz)+zcoefe(jz)*zxe3/(1.+0.75*zxe3))/zasped
        zelosr=zelosr+zde*arlong(ii+jz)
        zdcnpr=zdcnpr+zde*arlong(ii+jz+izlong)
2240    continue
        zdcnpr=zdcnpr*zgamma/zelosr
        zelosr=zelosr*evsinv
c               does orbit need updating ?
        ilpass=.true.
        zde=adeion(jptcle)+zelosr*zdt*uist
        if(zde.lt.ztest) go to 2250
        zdt=min(ztleft,(ztest-adeion(jptcle))/(zelosr*uist))
        if(areion(jptcle)-adeion(jptcle).lt.elecut) go to 2250
        ilpass=.false.
2250    continue
c               put heat into zawe and zawi
        zaengy=areion(jptcle)-adeion(jptcle)-0.5*zelosr*zdt*uist
        ztfrac=zdt/zdti2
        do 2260 jz=izbgin,izstop
        zt=zwght*ztfrac*arlong(ii+jz)
        zawi(jz)=zawi(jz)+zt*zcoefi(jz)/zasped
        zxe3=(zasped*zcoefx(jz))**3
        zawe(jz)=zawe(jz)+zt*zcoefe(jz)*zxe3/(zasped*(1.+0.75*zxe3))
        zarns(jz)=zarns(jz)+zt
        zares(jz)=zares(jz)+zaengy*evs*zt
2260    continue
c               decide to increment adeion, update, or dump into aeslow
        zde=zelosr*zdt*uist
        efast=areion(jptcle)-adeion(jptcle)-zde
        if(.not.ilpass) go to 2300
        if(efast.le.elecut) go to 2400
        speed=zsped*sqrt(efast)
        zdc=zdcnpr*(arsped(jptcle)-speed)
        adeion(jptcle)=adeion(jptcle)+zde
        adcnap(jptcle)=adcnap(jptcle)+zdc
        arsped(jptcle)=speed
        if(zdt.eq.ztleft) go to 2270
        ztleft=ztleft-zdt
        zdt=ztleft
        go to 2210
 2270   continue
        actvfe=actvfe+arwght(jptcle)*(areion(jptcle)-adeion(jptcle))
        actvfw=actvfw+arwght(jptcle)
        npactv=npactv+1
        go to 2490
c
c               update particle orbit
c
 2300   continue
        speed=zsped*sqrt(efast)
        zdc=zdcnpr*(arsped(jptcle)-speed)
        gamma=zgamma*speed
        slarmp=zslrm*speed
        v0drft=zv0drf*speed**2
        lzinnr=izbgin
        inrold=lzinnr
        lzoutr=izstop
        smj0=arstar(jptcle)
        canap=arcnap(jptcle)-(adcnap(jptcle)+zdc)
        llskip=.false.
        if(nlnewo) call bounce(iquick,ibounc)
        if(llskip) go to 2400
        areion(jptcle)=efast
        arcnap(jptcle)=canap
        arsped(jptcle)=speed
        adcnap(jptcle)=0.
        adeion(jptcle)=0.
        iptcle=jptcle
        if(nlnewo) call arlput(iptcle)
        ztleft=ztleft-zdt
        zdt=ztleft
        go to 2210
 
c
c               dump particle into aeslow
c
 2400   if(llskip) mskip=mskip+1
        if(.not.llstop) go to 2410
2410    continue
        nldead=nldead+1
        npdead=npdead+1
        nadead(npdead)=jptcle
        arlong(iblong)=-abs(arlong(iblong))
        nblong(jptcle)=0
        if(mlost) go to 2490
        zz=1.
        if(izstop.lt.npznes) go to 2417
        if(fnoutr.eq.1.) go to 2417
        zz=0.
        ii=iblong+itjump-izbgin
        do 2415 jz=izbgin,izstop
        zz=zz+arlong(ii+jz)
2415    continue
        zz=1./zz
2417    continue
        zeslow=zz*arwght(jptcle)*efast*evs
        zwslow=zz*arwght(jptcle)
        ii=iblong+itjump-izbgin
        do 2420 jz=izbgin,izstop
        aeslow(jz)=aeslow(jz)+zeslow*arlong(ii+jz)
        awslow(jz)=awslow(jz)+zwslow*arlong(ii+jz)
 2420   continue
 2490   continue
2500    continue
c
c
c
c
c
c
c       put blanket energy in outermost plasma zone
c
        if(fnoutr.ne.1.) go to 2595
        is=npznes+1
        if(is.gt.ntznes) go to 2595
        do 2590 jz=is,ntznes
        zawe(npznes)=zawe(npznes)+zawe(jz)
        zawi(npznes)=zawi(npznes)+zawi(jz)
2590    continue
2595    continue
c
c       convert zawe,zawi,zarns,and zasalf to densities
c       convert zares to average energy per particle
c
c
        zz = 2.0 * vols(mzones,1)
        zsum=0.
        zdr=1./float(npznes)
        do 2600 jz=1,ntznes
        zsum=zsum+zawe(jz)+zawi(jz)
        zvol=zz*(float(jz)-0.5)*zdr**2
        zawe(jz)=zawe(jz)/zvol
        zawi(jz)=zawi(jz)/zvol
        if(zarns(jz).ne.0.) zares(jz)=zares(jz)/zarns(jz)
        zarns(jz)=zarns(jz)/zvol
        zasalf(jz)=zasalf(jz)/zvol
 2600   continue
c
c       fill zbwe,zbwi,zbrns,and zbres to convert to baldur zones
c
        ii=1
        iuznes=npznes
        zdr=1./float(npznes)
        zvol = 2.0 * vols(mzones,1)
        zbsum=0.
        do 2660 jz=lcentr,ledge
 2610   zsmnc=zdr*(float(ii)-0.5)
        if(xzoni(jz).le.zsmnc) go to 2620
        if(ii.eq.iuznes) go to 2640
        ii=min(ii+1,iuznes)
        go to 2610
 2620   continue
        if(ii.gt.1) go to 2630
        zbwi(jz)=zawi(1)
        zbwe(jz)=zawe(1)
        zbrns=zarns(1)
        zbres=zares(1)
        zbsalf=zasalf(1)
        go to 2650
 2630   continue
        zz=(xzoni(jz)-(zsmnc-zdr))/zdr
        zbwi(jz)=(1.-zz)*zawi(ii-1)+zz*zawi(ii)
        zbwe(jz)=(1.-zz)*zawe(ii-1)+zz*zawe(ii)
        zbrns=(1.-zz)*zarns(ii-1)+zz*zarns(ii)
        zbres=(1.-zz)*zares(ii-1)+zz*zares(ii)
        zbsalf=(1.-zz)*zasalf(ii-1)+zz*zasalf(ii)
        go to 2650
2640    continue
        zbwi(jz)=zawi(iuznes)
        zbwe(jz)=zawe(iuznes)
        zbrns=zarns(iuznes)
        zbres=zares(iuznes)
        zbsalf=zasalf(iuznes)
 2650   continue
        alphai(jz)=alphai(jz)+zdt20*usid*zbrns
        ealfai(jz)=ealfai(jz)+zdt20*usie*zbres
        salfs(jz)=salfs(jz)+zdt20*zbsalf
        zbsum=zbsum+(zbwe(jz)+zbwi(jz))*dx2i(jz)*zvol
 2660   continue
c
c       increment wealfs,wialf,alphai,and ealfai
c
        zz=zsum/zbsum
        do 2670 jz=lcentr,ledge
        wealfs(jz)=wealfs(jz)+zbwe(jz)*zz*zdt20
        wialfs(jz)=wialfs(jz)+zbwi(jz)*zz*zdt20
2670    continue
c
c       reset nlbrth
c
        do 4950 j=1,mxshrt
        nlbrth(j)=.false.
4950    continue
5000    continue
        return
c
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
c
cl      10)     compression
c
c
c
 7000   continue
c
        if (.not.nlcomp) return
c
c               compression of density and energy
c
        call mesage(' compression of alphas not implemented ')
        return
c
c
c
        end
c
c**********************************************************************c
c@arlput  .../bald/code/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine arlput(kptcle)
c
c
c  stores information in arlong array
c  arlong(17000) = time per zone for each sample alpha particle
c
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
        iptcle=kptcle
        if((iptcle.lt.1).or.(iptcle.gt.nsampl)) call alfend(4)
        if((lzinnr.lt.1).or.(lzinnr.gt.ntznes)) call alfend(5)
        if((lzoutr.lt.1).or.(lzoutr.gt.ntznes)) call alfend(6)
        nzbgin(iptcle)=lzinnr
        izlong=lzoutr+1-lzinnr
        nzlong(iptcle)=izlong
        ilngth=2+2*izlong
        iblong=nblong(iptcle)
        if(iblong.eq.0) go to 200
        if(arlong(iblong).le.0.) call alfend(7)
        if(int(arlong(iblong+1)).ne.iptcle) call alfend(8)
        if(ilngth.le.int(abs(arlong(iblong)))) go to 310
 100    arlong(iblong)=-abs(arlong(iblong))
200     continue
        if(nulong+ilngth.le.mxlong) go to 300
        nsqeze=nsqeze+1
        call squeze
        if(nulong+ilngth.gt.mxlong) call alfend(9)
c
c
c
300     continue
        iblong=nulong
        nulong=nulong+ilngth
        arlong(iblong)=float(ilngth)
        nblong(iptcle)=iblong
        arlong(iblong+1)=float(iptcle)
 310    continue
        if(iblong+2*izlong+2.gt.mxlong) call alfend(10)
        arlong(iblong)=abs(arlong(iblong))
        do 320 j=1,izlong
        ii=iblong+1+j
        iz=lzinnr-1+j
        arlong(ii)=timlar(iz)
        arlong(ii+izlong)=timdcp(iz)
320     continue
        return
        end
c
c**********************************************************************c
c@squeze  .../bald/code/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine squeze
c
c  Repacks arlong array
c
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
        if(nulong.gt.mxlong) call alfend(11)
        iulold=nulong
        iput=1
        ilook=1
 100    if(ilook.ge.iulold) return
        ilngth=int(arlong(ilook))
        if(ilngth.gt.0) go to 200
        if(ilngth.eq.0) call alfend(12)
        ilook=ilook-ilngth
        go to 100
 200    ixptcl=int(arlong(ilook+1))
        if(ilook.ne.nblong(ixptcl)) call alfend(13)
        if(ilook.eq.iput) go to 500
        imove=1+2*nzlong(ixptcl)
        if(ilngth.lt.imove+1) call alfend(14)
        if(iput+imove.gt.mxlong) call alfend(15)
        do 250 j=1,imove
        arlong(iput+j)=arlong(ilook+j)
 250    continue
        arlong(iput)=imove+1
        nblong(ixptcl)=iput
        iput=iput+imove+1
        nulong=iput
        ilook=ilook+ilngth
        if(iput.gt.ilook) call alfend(16)
        go to 100
 500    continue
        if(ilngth.lt.2+2*nzlong(ixptcl)) call alfend(14)
        imove=min(ilngth,2*nzlong(ixptcl)+4)
        iput=iput+imove
        nulong=iput
        arlong(ilook)=imove
        ilook=ilook+ilngth
        go to 100
        end
c
c**********************************************************************c
c@alfend  .../bald/code/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine alfend(k)
       include 'cparm.m'
      include 'cbaldr.m'
        ierror=k
        if((ierror.lt.1).or.(ierror.gt.19)) go to 5000
        go to (100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,
     x  1400,1500,1600,1700,1800,1900,2000),k
100     continue
        call mesage(' alfini:  ledge .gt. mxznes ')
        call error_olymp(1,2,19,ierror, 
     x	' alfini:  ledge .gt. mxznes ')
        go to 6000
200     continue
        call mesage(' using old start. call alfini, not alphas(1) ')
        call error_olymp(1,2,19,ierror,
     x  ' using old start. call alfini, not alphas(1) ')
        go to 6000
300     continue
        call mesage(' alphas: can''t get zbrad ')
        call error_olymp(1,2,19,ierror, ' alphas: can''t get zbrad ')
        go to 6000
400     continue
        call mesage(' arlput:  iptcle out of legal range ')
        call error_olymp(1,2,19,ierror, 
     >	' arlput:  iptcle out of legal range ')
        go to 6000
 500    continue
        call mesage('  arlput:  lzinnr out of bounds ')
        call error_olymp(1,2,19,ierror, 
     x	'  arlput:  lzinnr out of bounds ')
        go to 6000
 600    continue
        call mesage('  arlput:  lzoutr out of bounds ')
        call error_olymp(1,2,19,ierror, 
     >  ' arlput:  lzoutr out of bounds ')
        go to 6000
700     continue
        call mesage(
     >    '  arlput:arlong(iblong).le.0. ')
        call error_olymp(1,2,19,ierror, 
     >    '  arlput:arlong(iblong).le.0. ')
        go to 6000
800     continue
        call mesage(' arlput:  arlong(iblong+1).ne.iptcle ')
        call error_olymp(1,2,19,ierror,
     >           '  arlput: arlong(iblong+1).ne.iptcle')
        go to 6000
900     continue
        call mesage('  arlput:  no room even after squeeze ')
        call error_olymp(1,2,19,ierror,
     >     '  arlput: no room even after squeeze')
        go to 6000
1000    continue
        call mesage(' arlput:  iblong+2*izlong+2.gt.mxlong ')
        call error_olymp(1,2,19,ierror,
     x	' arlput: iblong+2*izlong+2.gt.mxlong')
        go to 6000
1100    continue
        call mesage(' squeze disaster:nulong.gt.mxlong ')
        call error_olymp(1,2,19,ierror,
     x	' squeze disaster:nulong.gt.mxlong ') 
        go to 6000
1200    continue
        call mesage(' squeze:  ilngth.eq.0 ')
        call error_olymp(1,2,19,ierror, ' squeze:  ilngth.eq.0 ')
        go to 6000
1300    continue
        call mesage('  squeze disaster:  ilook.ne.nblong(ixptcl) ')
        call error_olymp(1,2,19,ierror,
     >         ' squeze disaster: ilook.ne.nblong(ixptcl)')
        go to 6000
1400    continue
        call mesage(' squeze disaster:  ilngth.lt.imove+1 ')
        call error_olymp(1,2,19,ierror,
     >          ' squeze disaster:  ilngth.lt.imove+1 ')
        go to 6000
1500    continue
        call mesage(' squeze disaster:mxlong is too small ')
        call error_olymp(1,2,19,ierror,
     >          ' squeze disaster:mxlong is too small ')
        go to 6000
1600    continue
        call mesage(' squeze disaster:iput.gt.ilook ')
        call error_olymp(1,2,19,ierror, 
     x	' squeze disaster:iput.gt.ilook ')
        go to 6000
1700    continue
        call mesage(' alfset:121   mxznes is too small ')
        call error_olymp(1,2,19,ierror, 
     x	   ' alfset:121   mxznes is too small ')																					
        go to 6000
1800    continue
        call mesage(' xinter error:no solution ')
        call error_olymp(1,2,19,ierror, ' xinter error:no solution ')
        go to 6000
1900    continue
        call mesage(' endzon disaster: no solution possible ')
        call error_olymp(1,2,19,ierror, 
     >         ' endzon disaster: no solution possible ')
        go to 6000
2000    continue
5000    continue
        call mesage(' alfend:  called with bad ierror ')
        call error_olymp(1,2,19,ierror, 
     >         ' alfend:  called with bad ierror ')
6000    continue
        nlend=.false.
        call endrun
        stop
        end
c
c**********************************************************************c
c@alfset  .../bald/code/bald/dfusion.f
c rgb 28-Jan-97 rmins*rndup before 123, cleaned up
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine alfset
c
c
c  Sbrtn alfset initializes geometric and plasma profile arrays
c        used by sbrtn bounce
c
c
        logical ilquit
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
      include 'commhd.m'
c
c
        iz=1
        zeion=efast
        zqb=qfast
        zab=afast
        zslrm=fcc*afast*fcmp*(10.**(fxc+fxnucl-fxes))/(fces*
     x  qfast*bzs*rmins)
        zsped=sqrt(2.*evs/(fcmp*(10.**fxnucl)*afast))
        speed=zsped*sqrt(zeion)
        vifast=1./speed
        aspect=rmajs/rmins
        aspeci=1/aspect
        zfwgt=sqrt(1.6e-12*zeion/(zab*1836.))
c
c
c       set up constants related to first orbit loss
c
        zrdif = rcwals - rmajs
        if ( rmins+abs(zrdif) .gt. rdwals*rndup ) then
123       call mesage(' alfset error:123  rdwals is inside plasma')
          write (6,*) 'rmins + abs(zrdif) = ',
     &      rmins + abs(zrdif)
          write (6,*) ' rmins = ', rmins
          write (6,*) ' rmajs = ', rmajs
          write (6,*) ' rcwals = ', rcwals
          write (6,*) ' rdwals*rndup = ',rdwals * rndup
          rdwals = ( rmins + abs(zrdif) ) / rndup
        endif
        soffst=zrdif/rmins
        srwall=rdwals/rmins
100     npmnzn=nsubzn*npznes
        zdmnzn=1./float(npmnzn)
        zz=srwall+abs(soffst)
        ntmnzn=int(zz/zdmnzn)
        if(ntmnzn.gt.npmnzn) go to 105
        ntmnzn=npmnzn
        go to 110
105     continue
        if(zdmnzn*float(ntmnzn+1).ge.zz) go to 110
        ntmnzn=ntmnzn+1
        go to 105
110     continue
        if(ntmnzn+1.le.mxmnzn) go to 115
        nsubzn=nsubzn-2
        go to 100
115     continue
        ntznes=ntmnzn/nsubzn
116     continue
        if(ntznes*nsubzn.ge.ntmnzn) go to 120
        ntznes=ntznes+1
        go to 116
120     continue
        if(ntznes.gt.mxznes) call alfend(17)
        slarmp=zslrm*speed
        ztaucr=soffst/(aspeci*slarmp)
        zls1=soffst*(soffst-2*slarmp*aspeci*(srwall-slarmp))+
     x  (soffst*aspeci*slarmp)**2
        zls2=(aspeci*slarmp)**2
        zls3=(srwall-slarmp)**2-soffst**2
        zls4=soffst-slarmp*aspeci*(srwall-slarmp)
c
c
c       set up argrtr(j),arlssr(j) arrays for j = 1,ntmnzn+1,iniskp
c
        iniskp=nsubzn/2
        ie=ntmnzn/iniskp+1
        do 150 jq=1,ie
        zsmn=zdmnzn*float(ntmnzn-(jq-1)*iniskp)
        if(ztaucr.lt.abs(zsmn-abs(soffst))) go to 130
        if(ztaucr.gt.abs(zsmn+abs(soffst))) go to 140
125     ilquit=.false.
        zxlssr=+100.
        zxgrtr=-100.
        go to 145
130     ztaumx=abs(zsmn-soffst)+slarmp*(1.+aspeci*zsmn)
        zxlssr=-zsmn
        zxgrtr=zsmn
        ilquit=.true.
        if(ztaumx.le.srwall) go to 145
        ilquit=.false.
        ztaumn=abs(zsmn+soffst)+slarmp*(1.-aspeci*zsmn)
        if(ztaumn.ge.srwall) go to 125
        zxlgit=1.
        call xinter(zsmn,zls1,zls2,zls3,zls4,zxlgit)
        zxgrtr=zxlgit
        go to 145
140     ztaumx=abs(zsmn+soffst)+slarmp*(1.-aspeci*zsmn)
        zxgrtr=zsmn
        zxlssr=-zsmn
        ilquit=.true.
        if(ztaumx.lt.srwall) go to 145
        ilquit=.false.
        ztaumn=abs(zsmn-soffst)+slarmp*(1.+aspeci*zsmn)
        if(ztaumn.ge.srwall) go to 125
        call xinter(zsmn,zls1,zls2,zls3,zls4,zxlgit)
        zxlssr=zxlgit
        go to 145
145     argrtr(jq)=1.+aspeci*zxgrtr
        arlssr(jq)=1.+aspeci*zxlssr
        llquit(jq)=ilquit
150     continue
c
c
c
c       set up arfphi(j) for j = 1,npmnzn+1
        iz=1
        zaphi=0.
        znorm=bpols(iz,ledge+1)
        ii=lcentr
        ie=npmnzn+1
        arfphi(1)=0.
        do 250 j=2,ie
        zsmnc=zdmnzn*(float(j)-1.5)
210     if(zsmnc.lt.xbouni(ii+1)) go to 220
        if(ii.ge.ledge) go to 220
        ii=ii+1
        go to 210
220     continue
        zsbpol=bpols(iz,ii)+(bpols(iz,ii+1)-bpols(iz,ii))
     x  *(zsmnc-xbouni(ii))/dxzoni(ii)
        zaphi=zaphi-zdmnzn*zsbpol
        arfphi(j)=zaphi/znorm
250     continue
        ie=ntmnzn+1
        is=npmnzn+2
        if(ie.lt.is) go to 300
        do 275 j=is,ie
        zsmn=zdmnzn*float(j-1)
        arfphi(j)=aphif(zsmn)
275     continue
c
c
c
c       set up raduis(j1), brfact(j1), j1=1,ntmnzn+1
300     continue
        ii=lcentr
        ie=ntmnzn+1
        do 350 j=1,ie
        zsmn=zdmnzn*float(j-1)
        raduis(j)=zsmn
310     if(zsmn.le.xbouni(ii+1)) go to 320
        if(ii.ge.ledge) go to 330
        ii=min(ii+1,ledge)
        go to 310
320     continue
        zsbpol=bpols(iz,ii)+(bpols(iz,ii+1)
     x  -bpols(iz,ii))*(zsmn-xbouni(ii))/dxzoni(ii)
        go to 340
330     continue
        zsbpol=bpols(iz,ledge+1)/zsmn
340     continue
        brfact(j)=zsbpol/bzs
350     continue
c
        return
        end
c
c**********************************************************************c
c@xinter  .../bald/code/bald/dfusion.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine xinter(parg1,pls1,pls2,pls3,pls4,plegit)
        zc=pls3-parg1**2
        zd=sqrt(pls1+pls2*parg1**2)
        if(pls4) 100,100,200
 100    zx1=(-pls4+zd)/pls2
        zx2=-zc/(pls4-zd)
        go to 300
 200    zx1=-zc/(pls4+zd)
        zx2=-(pls4+zd)/pls2
 300    zmin=min(abs(zx1),abs(zx2))
        if(zmin.gt.parg1) go to 500
        zmax=max(abs(zx1),abs(zx2))
        if(zmax.lt.parg1) call mesage(' xinter error: double solution ')
        if(abs(zx2).lt.abs(zx1)) go to 400
        plegit=zx1
        return
 400    plegit=zx2
        return
500     call alfend(18)
        end
c**********************************************************************c
c@bounce  .../bald/code/bald/dfusion.f
c  rgb 30-apr-87 removed all references to variable IXLARM
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
        subroutine bounce(k1,k2)
c
c
c  Sbrtn bounce calculates alpha trajectory
c        and collisions with background plasma
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
      include 'commhd.m'
c
        dimension zrvtot(155),zrvrad(155),zrx(75),zrcos(75),ztmdcp(75)
        dimension times(75),errors(75),zrsx(155),zrsz(155),
     x  zrcosp(155),
     x  ztime(75),zvtim(75),zerad(75),zertot(75),toterr(75)
        asinpi(z)=z*(0.31513+0.0754*z**2)
c
c
        iponts=10
        izntot=1
cbate        ixlarm=1
        iquick=k1
        ibounc=k2
        if(iquick.ne.1) nbounc=nbounc+1
        iutty=5
        mlost=.false.
        llstop=.false.
        llthin=.false.
        llstuk=.false.
        if(ibounc.ne.1) go to 100
c
c
c       initialize for ibounc equal to one
c
c
        zismn=sqrt(sxion**2+szion**2)
        zismj=(1.+aspeci*sxion)
        z0smj=zismj*(1.-cospit**2)
        canap=aphif(zismn)+sign(1.,cospit)*gamma
     x  *sqrt(zismj*(zismj-z0smj))
        zdsmn=raduis(2)-raduis(1)
        iion=int(zismn/zdsmn)+1
        if(zismn/zdsmn-iion+1.gt.0.5) iion=iion+1
        iion=min(ntmnzn,iion)
        iion=max0(2,iion)
        do 50 j=1,3
        iq=iion-2+j
        zc=(canap-arfphi(iq))/gamma
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        zx=aspect*(zsmj-1.)
        if(abs(zx).gt.raduis(iq)) go to 49
        zcospt=zc/zsmj
        zlarm=slarmp*speed*vifast*(1.+aspeci*zx)*sqrt(1.-zcospt**2)
        zsqtau=raduis(iq)**2+soffst*(soffst-2.*zx)
        zsqtst=(srwall-zlarm)**2
        if(zsqtau.gt.zsqtst) mlost=.true.
        iion=iq
        if(.not.mlost) go to 320
        if(nxorb.eq.1) go to 320
        go to 850
 49     continue
  50    continue
        llstuk=.true.
  55    nstart=nstart+1
        go to 850
c
c
c       initialize for ibounc not equal to one
c
c
 100    z0smj=smj0
        llstop=.false.
        linnr=1+nsubzn*(lzinnr-1)
        loutr=nsubzn*lzoutr
        do 125 j=1,2
        if(j.eq.2) go to 110
        js=(loutr+linnr)/2
        je=loutr
        go to 112
 110    je=js
        js=linnr
 112    continue
        do 120 jini=js,je
        iion=jini
        zc=(canap-arfphi(jini))/gamma
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        sxion=aspect*(zsmj-1.)
        if(abs(sxion).gt.raduis(jini)) go to 119
        szion=sqrt(raduis(jini)**2-sxion**2)
        zismj=1.+aspeci*sxion
        cospit=zc/zsmj
        z0smj=zismj*(1.-cospit**2)
        zismn=raduis(iion)
        go to 320
119     continue
 120    continue
 125    continue
        nstop=nstop+1
        llstop=.true.
        go to 850
c
c       do first pass to get orbit width
c
 320    iniskp=nsubzn/2
 325    continue
        do 290 jtwice=1,2
        if(jtwice.eq.2) go to 225
c
c       fast loop over outer zones
c
        isfast=ntmnzn+1-iniskp*((ntmnzn+1-iion)/iniskp+1)
        js=min(isfast+iniskp,ntmnzn+1)
        je=ntmnzn+1
        ioutr1=isfast
        mlost=.false.
        go to 230
c
c       fast loop over inner zones
c
 225    je=max0(2,isfast)
        js=2
        iinnr1=isfast+iniskp
 230    continue
        do 340 jloop=js,je,iniskp
        j=jloop
        if(jtwice.eq.2) j=js+je-jloop
        zc=(canap-arfphi(j))/gamma
        zt=zc**2
        iq=1+(ntmnzn+1-j)/iniskp
        if(zt.lt.arlssr(iq)*(arlssr(iq)-z0smj)) go to 332
        if(zt.gt.argrtr(iq)*(argrtr(iq)-z0smj)) go to 332
        go to 338
 332    if(llquit(iq)) go to 345
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        zx=aspect*(zsmj-1.)
        if(abs(zx).gt.raduis(j)) go to 345
        zcospt=zc/zsmj
        zlarm=slarmp*speed*vifast*(1.+aspeci*zx)*sqrt(1.-zcospt**2)
        zsqtau=raduis(j)**2+soffst*(soffst-2.*zx)
        zsqtst=(srwall-zlarm)**2
        if(zsqtau.gt.zsqtst) go to 341
 338    if(jtwice.eq.1)ioutr1=j
        if(jtwice.eq.2) iinnr1=j
 340    continue
        if(jtwice.eq.2) go to 345
 341    mlost=.true.
        if(nxorb.eq.1) go to 345
        go to 850
 345    if(iquick.ne.1) go to 290
        loutr=ioutr1
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        cosmax=zc/zsmj
        go to 850
 290    continue
c
c       is the orbit only one zone wide?
c
        if(ioutr1.gt.iinnr1) go to 360
        iinnr1=iion-1
        ioutr1=iion+1
        ibgskp=1
        if(lleror) ibgskp=2
        ityskp=1
        go to 400
c
c
c       set ibgskp so that there are at least iponts grid crossings
c
 360    continue
        idiff=ioutr1-iinnr1
        do 370 jsudo=1,nsubzn
        ibgskp=nsubzn+1-jsudo
        if(lleror) ibgskp=2*((ibgskp+1)/2)
        if(idiff.gt.ibgskp*iponts) go to 380
370     continue
 380    continue
        ityskp=ibgskp
        if(lleror) ityskp=ibgskp/2
c
c       set up zrcosp(j) and zrsx(j) and set iinnr and ioutr
c
 400    continue
        imidle=(iinnr1+ioutr1)/2
        imidle=ntmnzn+1-ityskp*((ntmnzn+1-imidle)/ityskp)
        ims=mxznes
        do 490 jtwice=1,2
        if(jtwice.eq.2) go to 435
c
c       set up for outer zones
c
        js=imidle+ityskp
        je=ntmnzn+1
        ioutr=imidle
        go to 440
c
c       set up for inner zones
c
 435    js=2
        je=imidle
        iinnr=imidle+ityskp
 440    continue
        do 420 jloop=js,je,ityskp
        j=jloop
        if(jtwice.eq.2) j=je+js-jloop
        ishort=ims+(j-imidle)/ityskp
        zc=(canap-arfphi(j))/gamma
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        zrcosp(ishort)=zc/zsmj
        zsx=(zsmj-1.)*aspect
        if(abs(zsx).gt.raduis(j)) go to 489
        zrsx(ishort)=zsx
        if(jtwice.eq.1) ioutr=j
        if(jtwice.eq.2) iinnr=j
 420    continue
        if(jtwice.eq.2) go to 489
        mlost=.true.
        if(nxorb.eq.1) go to 489
        go to 850
 489    continue
 490    continue
        if(ioutr.ge.iinnr+(izntot+1)*ibgskp) go to 500
        if((imidle.eq.iion).and.(ityskp.eq.1)) go to 494
        iinnr1=iion-1
        ioutr1=iion+1
        ibgskp=1
        if(lleror) ibgskp=2
        ityskp=1
        go to 400
494     continue
        llthin=.true.
495     narrow=narrow+1
        go to 850
c
c
c       set up end zones
c
c
 500    js=iinnr
        je=ioutr
        jd=je-js
        do 510 j=js,je,jd
        if(j.eq.js) go to 505
        if(mlost) go to 509
        je1=je
        je2=min(ntmnzn+1,je+ityskp)
        zc=raduis(je1)
        zx=raduis(je2)
        zxe=endzon(zc,zx,z0smj)
502     if(abs(zxe).gt.raduis(ioutr)) go to 509
        ioutr=ioutr-ityskp
        go to 502
 505    js1=js
        js2=max0(1,js-ityskp)
        zc=raduis(js1)
        zx=raduis(js2)
        zxs=endzon(zc,zx,z0smj)
507     if(abs(zxs).lt.raduis(iinnr)) go to 509
        iinnr=iinnr+ityskp
        go to 507
 509    continue
 510    continue
c
c
c
c
c       insure that iinnr and ioutr are not too close to 1 and ntmnzn
c
c
        itest=ntmnzn+1-ibgskp
 515    if(ioutr.le.itest) go to 516
        ioutr=ioutr-ityskp
        go to 515
 516    if(iinnr.ge.1+ibgskp) go to 520
        iinnr=iinnr+ityskp
        go to 516
 520    continue
c
c
c       insure that there are an even number of ityskp's
c
c
        if(.not.lleror) go to 525
        iq=(ioutr-iinnr)/ityskp
        if(iq.eq.2*(iq/2)) go to 525
        if(abs(zrsx(iinnr)/raduis(iinnr)).gt.abs(zrsx(ioutr)/
     x  raduis(ioutr))) go to 521
        ioutr=ioutr-ityskp
        go to 525
 521    iinnr=iinnr+ityskp
 525    continue
c
c
c       are there izntot ibgskp's?
c
c
        if(ioutr.ge.iinnr+(1+izntot)*ibgskp) go to 526
        llthin=.true.
527     narrow=narrow+1
        go to 850
 526    continue
        js=iinnr
        je=ioutr
        jd=je-js
 
 
        do 550 j=js,je,jd
        if(j.eq.js) go to 530
        if(mlost) go to 551
        zx=zxe
        j0=ioutr
        j1=ioutr+ityskp
        is1=ims+(j1-imidle)/ityskp
        j2=ioutr+ibgskp
        is2=ims+(j2-imidle)/ityskp
        ioutr=ioutr+ibgskp
        go to 540
 530    zx=zxs
        j0=iinnr
        j1=iinnr-ityskp
        is1=ims+(j1-imidle)/ityskp
        j2=iinnr-ibgskp
        is2=ims+(j2-imidle)/ityskp
        iinnr=iinnr-ibgskp
 540    continue
        do 549 ji=1,2
        if(ji.eq.2) go to 545
        if(.not.lleror) go to 548
        ja=is1
        zsmn=raduis(j0)+0.75*(abs(zx)-raduis(j0))
        go to 546
545     ja=is2
        zsmn=abs(zx)
 546    zq=npmnzn*zsmn+1.
        iq=min(ntmnzn,int(zq))
        zbfact=brfact(iq)+(zq-iq)*(brfact(iq+1)-brfact(iq))
        zc=(canap-aphif(zsmn))/gamma
        zsmj=0.5*(z0smj+sqrt(z0smj**2+4.*zc**2))
        zsx=(zsmj-1.)*aspect
        zrsx(ja)=zsx
        zrsz(ja)=0.
        if(ji.eq.1) zrsz(ja)=sqrt(zsmn**2-zsx**2)
        zcospt=zc/zsmj
        zrcosp(ja)=zcospt
        vdrift=v0drft*(1.+zcospt**2)
        vqpar=-zcospt*zbfact*speed
        zrvtot(ja)=sqrt(vqpar**2+vdrift**2+2.*vqpar*vdrift*zsx/zsmn)
 548    continue
 549    continue
 550    continue
 551    continue
c
c       set up zrvtot in the appropriate zones
c
        do 620 jl=1,2
        if(jl.eq.2) go to 605
        js=iinnr+ibgskp
        je=iinnr+ibgskp*(1+izntot)
        go to 607
 605    js=ioutr-ibgskp*(1+izntot)
        je=ioutr-ibgskp
        if(mlost) je=ioutr
 607    continue
        do 610 j=js,je,ityskp
        ishort=ims+(j-imidle)/ityskp
        zcospt=zrcosp(ishort)
        zsx=zrsx(ishort)
        zsmn=raduis(j)
        zrsz(ishort)=sqrt(zsmn**2-zsx**2)
        vdrift=v0drft*(1.+zcospt**2)
        vqpar=-zcospt*brfact(j)*speed
        zrvtot(ishort)=sqrt(vqpar**2+vdrift**2+2.*vqpar*vdrift*zsx/zsmn)
 610    continue
 620    continue
c
c       set up zrvrad in the appropriate zones
c
        js=iinnr+ibgskp*(1+izntot)
        je=ioutr-ibgskp*(1+izntot)
        do 630 j=js,je,ityskp
        ishort=ims+(j-imidle)/ityskp
        zsmn=raduis(j)
        zsz=sqrt(zsmn**2-zrsx(ishort)**2)
        zrsz(ishort)=zsz
        zrvrad(ishort)=v0drft*(1.+zrcosp(ishort)**2)*zsz/zsmn
 630    continue
c
c       time distribution
c
        if(mlost) go to 850
        do 700 j=1,ntznes
        ztmdcp(j)=0.
        times(j)=0.
        errors(j)=0.
 700    continue
        zvmx=0.
        zmax=0.
        itim=0
        zdsmn=raduis(1+nsubzn)-raduis(1)
        zdr=raduis(ibgskp+1)-raduis(1)
        ijle=iinnr+ibgskp*izntot
        ijge=ioutr-ibgskp*(1+izntot)
        ie=ioutr-ibgskp
        do 770 j=iinnr,ie,ibgskp
        il=ims+(j-imidle)/ityskp
        im=il+1
        iu=il+ibgskp/ityskp
        if(j.le.ijle) go to 730
        if(j.ge.ijge) go to 730
        ztm=zdr*0.5*(1./zrvrad(il)+1./zrvrad(iu))
        zd1=zrcosp(il)*(1.+aspeci*zrsx(il))/zrvrad(il)
        zd2=zrcosp(iu)*(1.+aspeci*zrsx(iu))/zrvrad(iu)
        zsmjcs=0.5*zdr*(zd1+zd2)
        go to 750
 730    continue
        zdist=sqrt((zrsx(iu)-zrsx(il))**2+(zrsz(iu)-zrsz(il))**2)
        ztm=zdist*0.5*(1./zrvtot(il)+1./zrvtot(iu))
        zd1=zrcosp(il)*(1.+aspeci*zrsx(il))/zrvtot(il)
        zd2=zrcosp(iu)*(1.+aspeci*zrsx(iu))/zrvtot(iu)
        zsmjcs=0.5*zdist*(zd1+zd2)
 750    continue
        if(j.ne.iinnr) go to 751
c
c       resetting zxs
c
        ii=ims+(iinnr-imidle)/ityskp
        zxs=abs(zrsx(ii))
        izone=int(abs(zxs)/zdsmn)+1
        zwgt=sqrt(abs(raduis(1+nsubzn*izone)**2-zxs**2))/
     x  sqrt(abs(raduis(iinnr+ibgskp)**2-zxs**2))
        zwgt=min(1.,zwgt)
        go to 756
751     if(j.eq.ie) go to 752
        zd2=raduis(j+ibgskp)
        go to 753
752     continue
        zd2=abs(zxe)
753     continue
        izone=(j-1)/nsubzn+1
        zd1=raduis(j)
 755    continue
        zwgt=min(1.,(raduis(izone*nsubzn+1)-zd1)/(zd2-zd1))
756     continue
        times(izone)=times(izone)+zwgt*ztm
        ztmdcp(izone)=ztmdcp(izone)+zsmjcs*zwgt
        if(zwgt.eq.1) go to 760
        times(izone+1)=times(izone+1)+(1.-zwgt)*ztm
        ztmdcp(izone+1)=ztmdcp(izone+1)+(1.-zwgt)*zsmjcs
 760    continue
 765    continue
        zrx(izone)=zrsx(il)
        zrcos(izone)=zrcosp(il)
        if(j.lt.ioutr-ibgskp) go to 769
        zrx(izone+1)=zrsx(il)
        zrcos(izone+1)=zrcosp(il)
 769    continue
 770    continue
c
c
c
c       larmor spread
c
 1830   continue
        do 1828 j=1,ntznes
        timlar(j)=times(j)
        timdcp(j)=ztmdcp(j)
 1828   continue
cbate        if(ixlarm.eq.0) go to 839
        zdsmn=raduis(nsubzn/2+1)-raduis(1)
        do 1835 j=1,ntznes
        timlar(j)=0.
        timdcp(j)=0.
 1835   continue
        izn=1-nsubzn/2
        do 840 jzn=1,ntznes
        if(times(jzn).eq.0.) go to 839
        zlarm=slarmp*speed*vifast*(1.+aspeci*zrx(jzn))*
     x  sqrt(1.-zrcos(jzn)**2)
        icentr=1
        if(jzn.ne.1) go to 1836
        if(zlarm.gt.4.*zdsmn) go to 1836
        icentr=5
        zlz=0.
        zlc=abs(zxs)
        zdc=sqrt((2.*zdsmn)**2-zxs**2)/float(icentr)
        times(jzn)=times(jzn)/float(icentr)
        ztmdcp(jzn)=ztmdcp(jzn)/float(icentr)
        go to 1840
1836    continue
        if(jzn.ge.npznes) go to 1840
        if(zlarm.gt.2.*zdsmn) go to 1840
        icentr=5
        zlc=abs(zxs)
        zsmnl=2.*zdsmn*float(jzn-1)
        zsmnu=2.*zdsmn*float(jzn)
        if(times(jzn+1).eq.0.) zsmnu=abs(zxe)
        zlz=0.
        if(zlc.lt.zsmnl) zlz=sqrt((zsmnl)**2-zlc**2)
        zdc=(sqrt((zsmnu)**2-zlc**2)-zlz)/float(icentr)
        times(jzn)=times(jzn)/float(icentr)
        ztmdcp(jzn)=ztmdcp(jzn)/float(icentr)
1840    continue
        do 838 jc=1,icentr
        zcentr=raduis(nsubzn*jzn+izn)
        if(icentr.eq.1) go to 1849
        zcentr=sqrt(zlc**2+(zlz+zdc*(float(jc)-0.5))**2)
1849    continue
        iznmin=int(float(npznes)*abs(zcentr-zlarm))+1
        iznmax=int(float(npznes)*abs(zcentr+zlarm))+1
        iznmax=min(iznmax,ntznes)
        if(iznmax.le.iznmin) go to 835
        zx1=1./(2.*zcentr*zlarm)
        zx2=(zcentr**2+zlarm**2)*zx1
        zalast=1.
        ie=iznmax-1
        do 830 jl=iznmin,ie
        zzx=zx1*(raduis(jl*nsubzn+1)**2)-zx2
        if(int(1.414*zzx)+0) 910,920,930
 910    zangle=1.-asinpi(sqrt(1.+1.e-6-(zzx)**2))
        go to 940
 920    zangle=0.5-asinpi(zzx)
        go to 940
 930    zangle=asinpi(sqrt(1.+1.e-6-(zzx)**2))
 940    continue
        timlar(jl)=timlar(jl)+times(jzn)*(zalast-zangle)
        timdcp(jl)=timdcp(jl)+ztmdcp(jzn)*(zalast-zangle)
        zalast=zangle
 830    continue
        timlar(iznmax)=timlar(iznmax)+times(jzn)*zalast
        timdcp(iznmax)=timdcp(iznmax)+ztmdcp(jzn)*zalast
        go to 836
 835    timlar(iznmin)=timlar(iznmin)+times(jzn)
        timdcp(iznmin)=timdcp(iznmin)+ztmdcp(jzn)
 836    continue
838     continue
839     continue
840     continue
        smj0=z0smj
        lzinnr=ntznes
        lzoutr=1
        zz=0.
        do 1800 jz=1,ntznes
        if(timlar(jz).eq.0.) go to 1799
        zz=zz+timlar(jz)
        if(jz.lt.lzinnr) lzinnr=jz
        if(jz.gt.lzoutr) lzoutr=jz
 1799   continue
 1800   continue
        do 1850 jz=lzinnr,lzoutr
        timlar(jz)=timlar(jz)/zz
        timdcp(jz)=timdcp(jz)/zz
 1850   continue
        if(fnoutr.ne.1.) lzoutr=min(npznes,lzoutr)
 850    llskip=.false.
        llskip=(mlost.or.llstuk).or.(llthin.or.llstop)
        return
        end
c
c
c%%%%%%%%%
c
c
        function aphif(parg1)
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
        zsmn=parg1
        zj=float(npmnzn)*zsmn+1.
        ij=int(zj)
        if(ij.gt.1) go to 200
        if(ij.lt.1) call mesage(' aphif error: negative zsmn ')
        aphif=(arfphi(2)-arfphi(1))*(zj-1.)**2
        return
 200    if(ij.gt.npmnzn) go to 500
        aphif=arfphi(ij)+(arfphi(ij+1)-arfphi(ij))
     x  *(zj-float(ij))
        return
 500    aphif=arfphi(npmnzn+1)-log(zsmn)
        return
        end
c
c
c%%%%%%%%%
c
c
        function endzon(z1,z2,z3)
       include 'cparm.m'
      include 'cbaldr.m'
       include 'calpht.m'
        fcp(z)=z-aspect*abs(0.5*(z0smj+sqrt(z0smj**2+4.*((canap-
     x  aphif(z))/gamma)**2))-1.)
        smnin=z1
        smnout=z2
        z0smj=z3
        iutty=5
        sgrtr=max(smnin,smnout)
        slssr=min(smnin,smnout)
        smnlst=smnout
        smnow=smnin
        fcplst=fcp(smnlst)
        fcpnow=fcp(smnow)
        if(fcplst.gt.0.) go to 900
        if(fcpnow.lt.0.) go to 910
        itrial=0
 400    smntrl=smnow-fcpnow*(smnlst-smnow)/(fcplst-fcpnow)
        if(smntrl.gt.sgrtr) go to 500
        if(smntrl.lt.slssr) go to 500
        smnlst=smnow
        smnow=smntrl
        fcplst=fcpnow
        fcpnow=fcp(smnow)
        if(abs(fcpnow).lt.1.e-6) go to 800
        itrial=itrial+1
        if(itrial.gt.5) go to 520
        if(fcpnow.lt.0.) go to 400
        if(abs(smnow-smnlst).lt.1.e-4) go to 800
 450    zzz=1.
        go to 400
 500    continue
        ndzerr=ndzerr+1
 520    zdlsmn=smnout-smnin
        smntrl=smnin
 600    zdlsmn=0.5*zdlsmn
        if(abs(zdlsmn).lt.1.e-4) go to 800
        if (fcp(smntrl+zdlsmn).gt.0.) smntrl=smntrl+zdlsmn
        go to 600
 800    endzon=smntrl
        return
c
c
c
 900    continue
        if(fcplst.gt.3.e-6) go to 920
        endzon=smnlst
        return
 910    if(fcpnow.lt.-3.e-6) go to 920
        endzon=smnow
        return
 920    continue
        call alfend(19)
        end
c--------1---------2---------3---------4---------5---------6---------7-
c---2---------3---------4---------5---------6---------7-c

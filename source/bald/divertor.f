c  11:00 25-nov-91 divertor
c/ 16.10 02-jun-89 /11040/bald89/wbaldn1 DIVERTOR, Bateman Stotler, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type
c cfs get /11040/bald89/wbaldn1
c end
c lib wbaldn1 ^ x divertor ^ end
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  File DIVERTOR consists of the following subprograms:
c
c pdx    : divertor scrape-off model
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
         subroutine pdx
c
c
c        2.20  divertor region sink and source terms
c
c**********************************************************************c
c@pdx
c
c  rap  08-mar-03 rid of equivalence statement
c  rap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c       dps 26-sep-88 15.03 Added option for constant sound speed /
c                     connection length with cimprd(1). It is implemented
c                     if it is > 0.
c       fgps 18-feb-83 made subrtn pdx compatible with having its
c                      main output array, scroff(mxchi,55), in common
c                      /comdf2/.
c       fgps 7-apr-80 introduced factors to adjust loss-rates for
c                     hydrogenic ions and impurity ions
c       fgps 12-jul-79 (1) included impurities in ave ion mass
c                      used for sound spd, zvs; (2) improved ch-ex
c                      energy dependence; and, (3) extended scrapeoff
c                      losses to impurity equations
c       fgps 16-apr-79 introduced the average hydrogen mass
c       jmo     feb-23-78       include density sink for more than
c               one hydrogen species,i.e. define pdiv(ihyd,j),
c               ihyd=lhyd1,lyhdn
c       jmo     jan-8-78 add charge exchange friction effects to
c               divertor-limiter model
c**********************************************************************c
c
c  output:    
c  scroff(mchi,55) = loss rate in scrapeoff due to divertor (sec**(-1))
c                    standard units
c
c**********************************************************************c
c
       include 'cparm.m'
       include 'cbaldr.m'
c
c
        real :: zarray(9)
        data iclass /2/,  isub /20/
!ap
        zarray(1:9) = cfutz(380:388)
!
c       if(.not.nlomt2(isub)) go to 10
        if(nadump(1).gt.epslon) go to  10
        call mesage(' *** 2.20 pdx subroutine bypassed ')
        return
10      continue
c
        if(nstep.gt.1) go to 12
        jmin=nadump(1)
        zpmass=fcau*(10.0**fxau)
   12   continue
c
c  start of loop over scrapeoff region
c
        do 20 j=jmin,ledge
        if(j.lt.nadump(3)) zl=cfutz(127)
        if(j.ge.nadump(3)) zl=cfutz(128)
        zhydm=ahmean(2,j)*zpmass
        zihydm=1./zhydm
c
c  inverse average ion mass
c
        z0=0.0
        z1=0.0
        do 110 jh=1,mhyd
        z0=z0+rhohs(jh,2,j)
        z1=z1+rhohs(jh,2,j)*aspec(jh)
  110   continue
        if(mimp.le.0) go to 122
        do 120 ji=1,mimp
        ii=ji+lhydn
        z0=z0+rhois(ji,2,j)
        z1=z1+rhois(ji,2,j)*aspec(ii)
  120   continue
  122   continue
        ziionm=z0/(z1*zpmass)
c
c
c  if nadump(4)=0, use sound speed model for plasma
c  flow into limiter/divertor
c  if nadump(4)=1, use model with neutral friction included
c
c  compute parallel flow velocity zvs of plasma into limiter/div.
c
c  sound speed model
c
        zvs=sqrt((tis(2,j)+tes(2,j))*ziionm)
        if(nadump(4).lt.epslon) go to 15
c  include charge exchange "friction" on plasma.
c
c  total density, average mass and temperature of the neutrals
c
        zsneu=0.0
        zmneu=0.0
        ztneu=0.0
        do 14 jh=1,mhyd
        zsneu=zsneu+rhons(jh,j)
        zmneu=zmneu+rhons(jh,j)*aspec(jh)
        ztneu=ztneu+rhons(jh,j)*tns(jh,j)
   14   continue
        z0=1./(zsneu+epslon)
        zmneu=z0*zmneu*zpmass
        ztneu=z0*ztneu
c
c  find average neutral velocity zv0
c
        zv0=0.0
        if(zsneu.gt.epslon) zv0=sqrt(3.*ztneu/zmneu)
c
c  compute average hydrogen ion velocity zvi
c
        zvi=sqrt(2.*tis(2,j)*zihydm)
c
c  compute average collision energy ze in terms of h-1
c
        zrelvv=zv0*zv0+zvi*zvi
        ze=0.5*zpmass*zrelvv
c
c  convert ze in units of ergs to units of evs
c
        ze=ze*evsinv
c
c  compute c-x cross-section zsig in cm**2
c
        zsig=0.6937e-14*(1.-0.155*log10(ze))**2/
     x  (1.+0.1112e-14*ze**3.3)
c
c  compute charge exchange frequency zfreq
c
        zfreq=zsneu*zsig*sqrt(zrelvv)
c
c  compute flow velocity of plasma along field zvs
c
        z00=sqrt(2.)*zvs
        zvs=z00*zvs/(z00+zl*zfreq)
c
15      continue
c
c       calculate transport coefficients= vs/l
c
        ztx=0.0
        if(zl.gt.epslon) ztx=zvs/zl
c
c  Constant value for vs/l if cimprd(1) > 0     
c
        if (cimprd(1).gt.epslon) ztx = cimprd(1)
        scroff(lelec,j)=-5.8*ztx
        scroff(lion,j)=-2.*ztx
        amach1=1.
        if(nstep.lt.1) go to 230
        if(cfutz(383).eq.0.0) go to 230

        zgas=float(ngas(1))
        call diver(j, tes(2,j),t2,rhoels(2,j),d2,u1,amach1,ztau2,rmajs,
     x  xzeff(1,j),q(j),zgas,zarray)
 230    continue
        do 17 jh=1,mhyd
        scroff(jh,j)=-cfutz(124)*ztx*amach1
 17    continue
c
        if(mimp.le.0) go to 19
c
        do 18 ji=1,mimp
        ii=ji+lhydn
        scroff(ii,j)=-cfutz(125)*ztx*amach1
   18   continue
   19   continue
   20   continue
c
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
        subroutine diver(pj,ptemp1,ptemp2,pdens1,pdens2,pu1,pamach1,
     x  ptau2,prdius,pzzeff,pq,pai,par)
        external case15
c**********************************************************************c
c@diver
c  rgb 02-jun-89 spread out data statements to avoid civic error
c
c       diver is a program to solve for the plasma transport
c       in a 1/2 dim model of a divertor (ie 2 points)
c       including sources.
c       the particle, momentum, and heat transport are solved
c       for the throat and plate, assuming all derivatives to
c       be between these points.
c       defining parameters, constants, and variables.
c       al2= neutral channel length
c       al1= field line length
c       temp1, temp2 = temperature entrance, plate
c       1=  entrance ,    2= plate
c       dens1, dens2= density
c       fpump= fraction of neutrals pumped
c       dltae= energy loss per ionization
c       r= outer rdius
c       amion= ion mass, ameltn= electron mass, am0= neutral mass  in ev
c       u1, u2  fluid velocity at 1 and 2  in cm/s
c       vthe= thermal velocity of electron  in cm/s
c       alam2= attenuation length of neutrals due to ionization
c       gamma1, gamma2 = fluxes of plasma
c       source = source function due to ionization
c       q1a= electron heat transfer due to advection
c       q1b= electron heat transfer due to conductivity
c       q1= minimum of q1a or q1b
c
c
c       evtmw converts ev/sec into megawatts
c
c
c       ai=1 for hydrogen plasma, =2 for deuterium
c       ispcie=1 for atomic neutrals, =2 for molecular neutrals
c
c
c
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxx
cxxxx      input to program
cxxxx       density at throat = dens1
cxxxx       temperature at throat = temp1
cxxxx       pumping in diverter (fractional) = fpump
cxxxx
cxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c
        common/cnstnt/qval,rdius,al1,al2,dltae,ai,amp,amion,am0,
     1  c,rion,e0,gi,ge,iunit,ameltn,bk,ispcie,iprint,secelt,evtmw
        common/pramtr/temp1,dens1,temp2,dens2,gamma1,gamma2,vthe,v0,u1,
     1  u2,q1a,q1b,q1,q2,fpump,alpha,zzeff,thetl,thetk
        common/outpt/alam2,tau2,source,akappa,amach1,amach2,deltaq,
     1  dflow,dsorce,alamda,usound,u1tst,gm1tst,qtest
c       convergence parameters
        data dlmn/1.e-3/, dtest/1.e-4/, x0/2./, delx/10./,
     &     xl/1.1/, xu/600./
c       bk= boltzmanns constant
        dimension par(9)
        rdius=prdius
        zzeff=pzzeff
        qval=pq
        ai=abs(pai)
        al2=par(1)
        fpump=par(2)
        dltae=par(3)
        ispcie=par(4)
        secelt=par(5)
        e0=par(6)
        alpha=par(7)
        thetl=par(8)
        thetk=par(9)
        temp1=ptemp1/bk
        dens1=pdens1
        al1=3.14159*qval*rdius
        aspcie=ispcie
c       masses in ev
        amion=ai*amp
        am0=aspcie*ai*amp
c       secelt= secondary electron emission
c       e0= desorbed molecule
c
c
c
c
        zdens1=dens1
        xu=temp1
        delx=(xu-xl)/10.
        delmn=dlmn*delx
        u1tst=(2.*bk*temp1/amion)**0.5
        gm1tst=dens1*u1tst
c       qtest=5.*gamma1*temp1 assuming u1=vth(i0n)
        qtest=5.*gm1tst*temp1*evtmw
        ztest=dtest*qtest
c       ztest is printed out as ptest
c
c
        call
     x  solve1(x0,delx,xl,xu,ztest,delmn,xzero,kfail1,kfail2)
c
c
        pu1=u1
        pamach1=amach1
        ptemp2=temp2*bk
        pdens2=dens2
        ptau2=tau2
        return
        end
c******************
      BLOCK DATA case15
        common/cnstnt/qval,rdius,al1,al2,dltae,ai,amp,amion,am0,
     1  c,rion,e0,gi,ge,iunit,ameltn,bk,ispcie,iprint,secelt,evtmw
c       physical parameters
        data bk/1.6e-12/, c/3.e10/, amp/1.67e-24/, ameltn/9.11e-27/,
     1  evtmw/1.6e-25/
        end

        function funqq(kfail1,kfail2,parg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c@funqq
c       solve 2 point heat transfer equation: substituting particle
c       and momentum balance and using definition of source function
        common/cnstnt/qval,rdius,al1,al2,dltae,ai,amp,amion,am0,
     1  c,rion,e0,gi,ge,iunit,ameltn,bk,ispcie,iprint,secelt,evtmw
        common/pramtr/temp1,dens1,temp2,dens2,gamma1,gamma2,vthe,v0,u1,
     1  u2,q1a,q1b,q1,q2,fpump,alpha,zzeff,thetl,thetk
        common/outpt/alam2,tau2,source,akappa,amach1,amach2,deltaq,
     1  dflow,dsorce,alamda,usound,u1tst,gm1tst,qtest
c       evtmw converts ev/sec to megawatts
c       deltaq is function to be solved funqq=0.0
        amach2=1.
        gi=1.
        ge=2.9*secelt
        temp2=parg
        rion=frate(ispcie,temp2)
        dens2=dens1*temp1/temp2/2.
        u2=(2.*bk*temp2/amion)**0.5
        gamma2=dens2*u2
        v0=(2.*bk*e0/am0)**0.5
        alam2=dens2*rion/v0
        tau2=alam2*al2
        source=gamma2*(1.-fpump)*(1.-exp(-tau2))
        gamma1=gamma2-source
        u1=gamma1/dens1
        usound=(2.*bk*temp1/amion)**0.5
        amach1=u1/usound
c       q (heat flux) evaluated in ev/cm**2/sec
c       flux limited theories
c       q1a=flux limited electron heat flow
c       q1b=electron heat flow by conductivity
c
        vthe=((2.*bk/ameltn)**0.5)*((1-thetl)*temp1+thetl*temp2)**0.5
        q1a=alpha*dens1*vthe*(temp1-temp2)*evtmw
        alamda=25.3-1.15*log10(dens1)+2.3*log10(temp1)
        if(temp1.gt.50.) alamda=23.4-1.15*log10(dens1)
     1  + 3.45*log10(temp1)
        zfzeff=(12.5-9.34/sqrt(zzeff))/(zzeff*3.161)
        akappa=(1.8e20)/(alamda/10.)*zfzeff*((1-thetk)*temp1+
     1  thetk*temp2)**2.5
        q1b=akappa*(temp1-temp2)/al1*evtmw
        q1=min(q1a,q1b)
        q2=2.*(gi+ge)*gamma2*temp2*evtmw
        dsorce=dltae*source*evtmw
        dflow=5.*gamma1*temp1*evtmw
        deltaq=q2-(dflow+q1) -dsorce
        funqq=deltaq
        return
        end
        function frate(kspcie,ptemp2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c@frate
c        frate calculates the rate coefficients for plasma-neutral
c        processes (ispcie=1 for atoms,  =2 for molecules)
c
c
c
c
        dimension zr1(7),zr2(7),zrasv(2),zrsv(2)
c       zr1  e + h ionization,  zr2  e +  h2 ionization
c       references: freeman and jones thermal average maxwellian
        data zr1/-.317385e2,.1143818e2,-.3833998e1,.7046692,
     1  -.7431486e-1,.4153749e-2,-.9486967e-4/
        data zr2/-.3265772e2,.1233251e2,-.4059211e1,.7328010,
     1  -.7596204e-1,.4181685e-2,-.9427185e-4/
        zye=log(ptemp2)
        do 2 k=1,2
 2      zrasv(k)=0.
        do 3 ka=1,7
        kb=ka-1
        zrasv(1)=zrasv(1)+zr1(ka)*zye**kb
        zrasv(2)=zrasv(2)+zr2(ka)*zye**kb
 3      continue
        do 4 ki=1,2
 4      zrsv(ki)=exp(zrasv(ki))
        frate=zrsv(kspcie)
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
        subroutine
     x  solve1(px0,pdelx,pxl,pxu,ptest,pdelmn,pxzero,kfail1,kfail2)
c
c**********************************************************************c
c@solve1
c       solver.for[477,4121] is a file prepared by c. singer
c       for w. langer on 17 may, 1983.  this file contains a
c       function solving routines used to solve either one
c       or two algebraic equations, the subroutines
c       solve1 and solve2 used to find the roots, and
c       of subroutines fsolv1 and fsolv2 used to call the
c       function(s) whose root(s) are w
c
c
c
c       calling sequence (neglecting calls to output, error, and page):
c             solve1
c               fsolv1
c                 funtec
c               solve2
c                 fsolv2
c                       solve1
c                         fsolv1
c
c  solve1--solve1 is virtually identical to mikkelsen's function
c       solver [1].  here it is used to search downward from a
c       specified maximum temperature to find solutions to
c       power balance equations.  it is also used to search upward
c       from zero ripple to find the ripple which balances ion
c       power in the coast phase.
c  fsolv1--the function passed to solve1 by the version of
c       fsolve used here is controlled by the common variable lfun.
c
c
c
c
c
c
c
c       general purpose algebraic equation solver
c       ces--solve1 solves the inner equation in a set of 1 or 2.
c       ces--divided by abs(...) to normalize sign test.
c       ces--additional test after go to 902
c       ces--overflow protection on zffr
c
c       px0     initial value of x for beginning search for root
c       pdelx   initial value of step size for x
c       pxl     lower bound of search region
c       pxu     upper bound of search region
c       ptest   if abs(f(x)) .lt. ptest then the root is found
c       pdelmn  lower limit on step size
c       pxzero  return variable for the root
c       kfail1  return flag for failure of outer equation
c       kfail2  return flag for failure of inner equation
c
c       set internal variables
c
        kfail1=0
        iscan=1
        icount=0
        zx0=px0
        zdelx=pdelx
        zxl=pxl
        zxu=pxu
        ztest=abs(ptest)
        zdelmn=pdelmn
c
c       test for consistent input values
c
        if(abs(zdelx).lt.zdelmn) go to 901
        if((zx0-zxl)*(zx0-zxu).gt.0.) go to 902
        if((zx0-zxl)*(zx0-zxu).eq.0..and.zx0.lt.zxl) go to 902
        if((zx0-zxl)*(zx0-zxu).eq.0..and.zx0.gt.zxu) go to 902
c
c       take initial step and set direction to go toward the root
c
100     continue
        icount=icount+1
        zxold=zx0
        zfold=fsolv1(kfail1,kfail2,zxold)
        if(zfold.eq.0.) go to 801
        zxnew=zx0+zdelx
        if(zxnew.lt.zxl) zxnew=zxl
        if(zxnew.gt.zxu) zxnew=zxu
        zfnew=fsolv1(kfail1,kfail2,zxnew)
        if(zfnew.eq.0.) go to 802
        if((zfold/abs(zfold))*(zfnew/abs(zfnew)).lt.0.) go to 500
        if(abs(zfnew).gt.abs(zfold))zdelx=-zdelx
        zstart=zxnew
c
c       begin regular stepping to find two points on
c       either side of the root
c
200     continue
        zxold=zxnew
        zfold=zfnew
        zxnew=zxnew+zdelx
        istop=0
        if((zxnew-zxl)*(zxnew-zxu).gt.0.) istop=1
        if(zxnew.lt.zxl) zxnew=zxl
        if(zxnew.gt.zxu) zxnew=zxu
        zfnew=fsolv1(kfail1,kfail2,zxnew)
        if(zfnew.eq.0.) go to 802
        if((zfold/abs(zfold))*(zfnew/abs(zfnew)).lt.0.) go to 500
        if(istop.eq.1) go to 903
        go to 200
c
c       begin using regula falsi method to rapidly approach the root
c
500     continue
c
c       set up in case we return to 100
c
        if(abs(zdelx).le.zdelmn) go to 906
        iscan=1
        zx0=zxold
        zxl=min(zxold,zxnew)-zdelx
        zxl=max(zxl,pxl)
        zxu=max(zxold,zxnew)+zdelx
        zxu=min(zxu,pxu)
        zdelx=zdelx/10.
        if(abs(zdelx).lt.zdelmn) zdelx=zdelmn
c
c       begin regula falsi
c
530     continue
        icount=icount+1
        if(icount.ge.50) go to 907
        zndiff=(zfnew-zfold)/zfold
        if(abs(zndiff).lt.1.e-10) zndiff=1.e-10
        zxfr=zxold-(zxnew-zxold)/zndiff
c
c       if outside the new limited search region go back to 100
c
        if((zxfr-zxl)*(zxfr-zxu).gt.0.) go to 100
        zffr=fsolv1(kfail1,kfail2,zxfr)
        if(abs(zffr).lt.ztest) go to 803
        zxold=zxnew
        zfold=zfnew
        zxnew=zxfr
        zfnew=zffr
        go to 530
c
c       succesful returns
c
801     continue
        pxzero=zxold
        return
802     continue
        pxzero=zxnew
        return
803     continue
        pxzero=zxfr
        return
c
c       unsuccesful returns
c
901     continue
c
c       initial step size is too small
c
        kfail1=1
        pxzero=0.
        return
c
c       starting point is out of bounds
c
902     continue
        kfail1=2
        pxzero=0.
        return
903     continue
c
c       if we have scanned only to one side of zx0, try the other
c
        if(iscan.eq.2) go to 905
        iscan=2
        zxnew=zstart
        zdelx=-zdelx
        go to 200
c
c       both sides scanned; failure
c
905     continue
        kfail1=3
        pxzero=0.
        return
c
c       step size has become too small
c
906     continue
        kfail1=4
        pxzero=zxnew
        return
c
c       it should never take this long
c
907     continue
        kfail1=5
        pxzero=0.
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
        subroutine
     x  solve2(px0,pdelx,pxl,pxu,ptest,pdelmn,pxzero,kfail1,kfail2)
c
c**********************************************************************c
c@solve2
c       general purpose algebraic equation solver
c       ces--solve2 solves the outermost equation in a set of 1 or 2.
c       ces--divided by abs(...) to normalize sign test.
c       ces--additional test after go to 902
c       ces--overflow protection on zffr
c
c       px0     initial value of x for beginning search for root
c       pdelx   initial value of step size for x
c       pxl     lower bound of search region
c       pxu     upper bound of search region
c       ptest   if abs(f(x)) .lt. ptest then the root is found
c       pdelmn  lower limit on step size
c       pxzero  return variable for the root
c       kfail1  return flag for failure of outer equation
c       kfail2  return flag for failure of inner equation
c
c       set internal variables
c
        kfail1=0
        kfail2=0
        iscan=1
        icount=0
        zx0=px0
        zdelx=pdelx
        zxl=pxl
        zxu=pxu
        ztest=abs(ptest)
        zdelmn=pdelmn
c
c       test for consistent input values
c
        if(abs(zdelx).lt.zdelmn) go to 901
        if((zx0-zxl)*(zx0-zxu).gt.0.) go to 902
        if((zx0-zxl)*(zx0-zxu).eq.0..and.zx0.lt.zxl) go to 902
        if((zx0-zxl)*(zx0-zxu).eq.0..and.zx0.gt.zxu) go to 902
c
c       take initial step and set direction to go toward the root
c
100     continue
        icount=icount+1
        zxold=zx0
        zfold=fsolv2(kfail1,kfail2,zxold)
        if1old=kfail1
        if(zfold.eq.0..and.if1old.eq.0) go to 801
        zxnew=zx0+zdelx
        if(zxnew.lt.zxl) zxnew=zxl
        if(zxnew.gt.zxu) zxnew=zxu
        zfnew=fsolv2(kfail1,kfail2,zxnew)
        if1new=kfail1
        if(zfnew.eq.0..and.if1new.eq.0) go to 802
        if((zfold/abs(zfold))*(zfnew/abs(zfnew)).lt.0.
     x  .and.if1old.eq.0.and.if1new.eq.0) go to 500
        if(abs(zfnew).gt.abs(zfold)
     x  .and.if1old.eq.0.and.if1new.eq.0)zdelx=-zdelx
        zstart=zxnew
c
c       begin regular stepping to find two points on
c       either side of the root
c
200     continue
        zxold=zxnew
        zfold=zfnew
        if1old=if1new
        zxnew=zxnew+zdelx
        istop=0
        if((zxnew-zxl)*(zxnew-zxu).gt.0.) istop=1
        if(zxnew.lt.zxl) zxnew=zxl
        if(zxnew.gt.zxu) zxnew=zxu
        zfnew=fsolv2(kfail1,kfail2,zxnew)
        if1new=kfail1
        if(zfnew.eq.0..and.if1new.eq.0) go to 802
        if((zfold/abs(zfold))*(zfnew/abs(zfnew)).lt.0.
     x  .and.if1old.eq.0.and.if1new.eq.0) go to 500
        if(istop.eq.1) go to 903
        go to 200
c
c       begin using regula falsi method to rapidly approach the root
c
500     continue
c
c       set up in case we return to 100
c
        if(abs(zdelx).le.zdelmn) go to 906
        iscan=1
        zx0=zxold
        zxl=min(zxold,zxnew)-zdelx
        zxl=max(zxl,pxl)
        zxu=max(zxold,zxnew)+zdelx
        zxu=min(zxu,pxu)
        zdelx=zdelx/10.
        if(abs(zdelx).lt.zdelmn) zdelx=zdelmn
c
c       begin regula falsi
c
530     continue
        icount=icount+1
        if(icount.ge.50) go to 907
        zndiff=(zfnew-zfold)/zfold
        if(abs(zndiff).lt.1.e-10) zndiff=1.e-10
        zxfr=zxold-(zxnew-zxold)/zndiff
c
c       if outside the new limited search region go back to 100
c
        if((zxfr-zxl)*(zxfr-zxu).gt.0.) go to 100
        zffr=fsolv2(kfail1,kfail2,zxfr)
        if(kfail1.ne.0) go to 908
        if(abs(zffr).lt.ztest) go to 803
        zxold=zxnew
        zfold=zfnew
        zxnew=zxfr
        zfnew=zffr
        go to 530
c
c
c       if(kfail1.ne.0) go to 908
c
c       succesful returns
c
801     continue
        pxzero=zxold
        return
802     continue
        pxzero=zxnew
        return
803     continue
        pxzero=zxfr
        return
c
c       unsuccesful returns
c
901     continue
c
c       initial step size is too small
c
        kfail2=1
        pxzero=0.
        return
c
c       starting point is out of bounds
c
902     continue
        kfail2=2
        pxzero=0.
        return
903     continue
c
c       if we have scanned only to one side of zx0, try the other
c
        if(iscan.eq.2) go to 905
        iscan=2
        zxnew=zstart
        zdelx=-zdelx
        go to 200
c
c       both sides scanned; failure
c
905     continue
        kfail2=3
        pxzero=0.
        return
c
c       step size has become too small
c
906     continue
        kfail2=4
        pxzero=zxnew
        return
c
c       it should never take this long
c
907     continue
        kfail2=5
        pxzero=0.
        return
c
c       no root of fsolv1 at nominal root of fsolv2
c
908     continue
        kfail2=6
        pxzero=0.
        return
        end
        function fsolv1(kfail1,kfail2,parg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        common/cnstnt/qval,rdius,al1,al2,dltae,ai,amp,amion,am0,
     1  c,rion,e0,gi,ge,iunit,ameltn,bk,ispcie,iprint,secelt,evtmw
        common/pramtr/temp1,dens1,temp2,dens2,gamma1,gamma2,vthe,v0,u1,
     1  u2,q1a,q1b,q1,q2,fpump,alpha,zzeff,thetl,thetk
        common/outpt/alam2,tau2,source,akappa,amach1,amach2,deltaq,
     1  dflow,dsorce,alamda,usound,u1tst,gm1tst,qtest
c       fsolv1 is called by solve1 to solve funtec, fundel,
c       or funted.
c
        fsolv1=funqq(kfail1,kfail2,parg)
        return
        end
        function fsolv2(kfail1,kfail2,parg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        common/cnstnt/qval,rdius,al1,al2,dltae,ai,amp,amion,am0,
     1  c,rion,e0,gi,ge,iunit,ameltn,bk,ispcie,iprint,secelt,evtmw
        common/pramtr/temp1,dens1,temp2,dens2,gamma1,gamma2,vthe,v0,u1,
     1  u2,q1a,q1b,q1,q2,fpump,alpha,zzeff,thetl,thetk
        common/outpt/alam2,tau2,source,akappa,amach1,amach2,deltaq,
     1  dflow,dsorce,alamda,usound,u1tst,gm1tst,qtest
c       fsolv2 is called by solve2 to solve funtid
c
c
        fsolv2=1.
        return
        end

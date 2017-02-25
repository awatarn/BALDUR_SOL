c  19:00 28-nov-91 dauxheat.f
c/ 16.00 17-dec-88 11040/bald89/wbaldn1 DAUXHEAT, Bateman, PPPL
c  BALDUR  file DAUXHEAT by Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type
c cfs get /11040/bald89/wbaldn1
c end
c lib wbaldn1 ^ x dauxheat ^ end
c
c**********************************************************************c
c
c     this file contains the following packages:
c  heat   - auxilliary heating model
c
c  ecrh   - electron cyclotron resonance heating 
c              to compute weecrh and wiecrh
c  rprint - print r-f heating summary
c
c  cmpres - adiabatic compression - not used in 1-1/2-D BALDUR
c              this sbrtn adds or subtracts zones as the plasma
c              expands or compresses
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c@heat in ../baldur/code/bald/dauxheat
c  rgb 10-feb-01 interpolate rlepwr over nrlepwr points if nrlepwr > 0
c    and rlipwr over nrlipwr points if nrlipwr > 0
c  rgb 12-jun-00 interpolate rlepwr and rlipwr as a function of bdtime
c    to rlepwrt and rlipwrt scalars at time tai
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 07-feb-00 removed equivalence (nlzzzz(1),zdn(1)), ...
c  rgb changed hollerith to '...' and revised printout of label
c     rgb 17-nov-85    1 1/2 D upgrade using common commhd
c        plasma surface area = <|del xi|> V'(xi)
c        plasma surface area = <|del xi|> V'(xi)
c       drm 4-dec-84 put mikkelsen's heat in standard area
c       drm 30-may-84 compare tai (not tbi) to t-on and t-off
c       drm 8-may-84 removed the line which jumped over the initialization
c               so that the profiles can be updated in mid-run
c       aes 19-apr-82 fix handling of zkprim in exceptional cases --
c               fix some comments -- use evsinv, not 6.2418e08
c       aes 12-mar-82 make times printed in hetprt same as in mprint
c               -- also fix turning off of hetprt before, after heating
c       aes 20-jan-82 added label5 to page headers
c       aes 01-dec-81 add extra linefeed in format 9010
c       aes 19-nov-81 turn off hetprt before and after heating
c       aes 19-nov-81 move error message 9000 to before entry hetprt
c               -- change label to 8100
c       aes 18-nov-81 commented out part of equivalence statement
c               to make entry hetprt work correctly
c       aes 17-nov-81 made message length = 48h in error call 9000
c       aes 17-nov-81 edit printout --> entry hetprt
c       aes 28-oct-81 put all data statements before executable code
c       aes 28-aug-81 change "sub" to "isub" in call to error
c       aes 20-nov-80 fixed normalization again; k-switch
c         for goto 20; renamed rleaux --> rleprf,
c         zweaux --> rleaux, etc.
c       aes 18-nov-80 reinstated zweaux,zwiaux
c       aes 18-nov-80 changed all 52's to mzones
c       aes 18-nov-80 fixed normalization near lines 210,310
c       aes 18-nov-80 reinstalled i2,i3 at initialization
c       aes 4-nov-80 removed lines around line 100 -- substituted
c               appropriate lines in preset, clear
c       aes 4-nov-80 changed variable names to rl,nrl prefixes
c       fgps 9-may-80 reinitialized zdn(*) to 0.0 each time step
c       fgps 21-feb-80 revised the print-out
c       fgps 12-feb-80 set up a calculation of lower-hybrid heating
c                      using jom's subroutine auxhet of 11-oct-78
c******************************************************************************
c
c
        subroutine heat(k)
c
c       2.5  computes auxiliary heating effects
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
        dimension
     a   zvols(mj)    , zwefit(mj)   ,
     a   zwifit(mj)   , zdn(mj)      , zleft(mj)    , zabs(mj)     ,
     a   zatten(mj)   , zepowr(mj)   ,
     a   zipowr(mj)   , zzzz1(mj)    , zzzz2(mj)    , zzzz3(mj)    ,
     i   ilmc(50)     , kzzz1(mj)
c
cgb        equivalence
cgb     x   (rleaux(1),zwefit(1))       , (rliaux(1),zwifit(1))       ,
cgb     x   (nlzzzz(1),zdn(1))          , (nlzzzz(56),zleft(1))       ,
cgb     x   (nlzzzz(111),zabs(1))       , (nlzzzz(166),zatten(1))
c*****   (nlzzzz(221),zepowr(1))     , (nlzzzz(276),zipowr(1))     ,
c*****   (nlzzzz(331),zzzz1(1))      , (nlzzzz(386),zzzz2(1))      ,
c*****   (nlzzzz(441),zzzz3(1))      , (nlzzzz(496),kzzz1(1))
c
c
c------------------------------------------------------------------------------
c
        data    iclass /2/,     isub /5/
        data  iauxht,ihton,ihtoff,ihtout/60,61,62,63/
c
        if(.not.nlomt2(isub))   go to 10
        call mesage(' *** 2.5 subroutine heat bypassed')
        return
   10   continue
c
c------------------------------------------------------------------------------
cinput
c       Input variables for the auxiliary heating package, sbrtr HEAT
c       -------------------------------------------------
c
c       input parameters for auxiliary heating are read in from two
c       sources:  (1) cfutz(60) through cfutz(62); and, (2) namelist
c       /nurun1/.  the cfutz-factors activate the heating effects,
c       while the namelist/nurun1/ introduces parameters which define
c       the heating process.
c
c        cfutz(iauxht) = switch which activates the auxiliary heating;
c              60        1.="on".
c
c        cfutz(ihton)  = time (secs) when the auxiliary heating is
c              61        turned on; defaulted to 0.
c
c        cfutz(ihtoff) = time (secs) when the auxiliary heating is
c              62        turned off.
c
c        cfutz(ihtout) = selects time step after which extensive print-
c              63        outs are continually made of quantities used
c                        in model 3; defaulted to infinity.
c
c        data  iauxht,ihton,ihtoff,ihtout/60,61,62,63/
c
c       the parameters included in namelist/nurun1/ are set out below:
c
c        nrlmod   = switch to select the auxiliary heating mechanism;
c                   defaulted to 3.
c
c       the following quantities apply when nrlmod=1 (tabular input)
c        rlipwr(j) = total power (watts) deposited in the ions;
c                   as a function of time at bdtime(j) breakpoints
c                   defaulted to 1.0e+06.
c        rlepwr(j) = total power (watts) deposited in the electrons;
c                   as a function of time at bdtime(j) breakpoints
c                   defaulted to 1.0e+06.
c        rliprf(j)= input profile determining rliaux(j)
c        rleprf(j)= input profile determining rleaux(j)
c        rliaux(j)= relative magnitude of power given to the ions in
c                   zone j.  here j=1 corresponds to lcentr and j goes
c                   from 1 to mzones.  the default value is 0.
c        rleaux(j)= relative magnitude of power given to the electrons
c                   in zone j.  again j=1 corresponds to lcentr, j goes
c                   from 1 to mzones, and the default value is 0.
c
c       the following quantites are used for nrlmod=2 (polynomial fit)
c        rlipwr(j) = total power (watts) deposited in the ions;
c                   as a function of time at bdtime(j) breakpoints
c                   defaulted to 1.0e+06.
c        rlepwr(j) = total power (watts) deposited in the electrons;
c                   as a function of time at bdtime(j) breakpoints
c                   defaulted to 1.0e+06.
c        rle0      = relative power given to electrons at center of the
c                   plasma.
c        rle1      = relative power given to electrons at the outer
c                   edge of the plasma.
c        nrleex    = exponent in the radial variation of relative power
c                   given to the electrons.
c        rli0      = relative power given to the ions at the center of
c                   the plasma.
c        rli1      = relative power given to the ions at the outer edge
c                   of the plasma.
c        nrliex    = exponent in the radial variation of relative power
c                   given to the ions.
c
c       the following quantities are used to set up nrlmod=3
c       ***************(lower-hybrid heating)***************
c        rlfreq   = frequency (hz) of rf source; defaulted to 5.0e
c                   +08.
c        rlipwr   = total power (watts) given in nrlmod=3 to both ions
c                   and electrons; defaulted to 1.0e+06.
c        nrlpar    = number of modes in the wave-number spectrum, the
c                   prescriptible number being limited to 50.
c        rlpara(i) = index of refraction for mode i, where the maximum
c                   value of i is nrlpar.
c        rlpowr(i)= relative input power prescribed for mode i; the de-
c                   fault value is 0.0 for each i.
c        nrldmp    = switch to select a model for the energy deposition:
c                   nrldmp=1 allows energy loss only after linear mode
c                   conversion; nrldmp=2 allows only electron landau
c                   damping; nrldmp=3 allows both mechanisms and is the
c                   default value.
c
cend
        if(k.eq.2) go to 20
        i1 = mxzone
        i2 = mzones
        i3 = 50
        call resetr(weauxs,i1,0.0)
        call resetr(wiauxs,i1,0.0)
   20   continue
        if(nstep.eq.0) return
        if(cfutz(iauxht).le.epslon) return
c
  100   continue
c
c       initialization
c
        if(cfutz(ihtoff).le.epslon) go to 8100
        zhton=cfutz(ihton)*usit
        zhtoff=cfutz(ihtoff)*usit
        ncheck=5000
        if(cfutz(ihtout).gt.epslon) ncheck=int(cfutz(ihtout)+0.1)
c
c       compute twice the total volume of the cylinder for use below
        ztovol = 2.0 * vols(mzones,1)
c
c       compute reciprocal zone volume as a fcn. of radius
c
        do 50 j=lcentr,mzones
          zvols(j)=dx2i(j)*ztovol
          zvols(j) = 1.0 / ( dx2i(j) * ztovol )
   50   continue
c
c       get reciprocal of ztovol
        ztovol = 1.0/ztovol
c
c       this subrountine calculates auxiliary heating terms
c       three models are used:
c
c               ** the power/zone is input at the beginning
c               of the run and held constant throughout
c               in models 1 and 2
c
c               ** the power/mode is input as a function of
c               parallel index of refraction rlpara=kpara*c/freq
c               the linear mode conversion point is calculated
c               for each wave no. and the amount of energy in that
c               mode is input to the radial zone which
c               contains the l.m.c. point. (model 3)
c               in this way the power deposition profile
c               changes self-consistently with the plasma profiles
c               the energy lost in each zone via electron-landau
c               damping is also calcaluted.
c                       if nrldmp=1, only lmc is used
c                       if nrldmp=2, only e-landau damping
c                       if nrldmp=3, both are used(default)
c
c
c       select model.
c       if  nrlmod = 1 use tabular input for power/zone
c                       in watts. rleprf,rliprf are
c                       set in the input namelist.
c
c                    2  use polynomial fit for power/zone
c                       in watts. pwr/zone=
c                       rle1+(rle0-rle1)*(1-(r/a)**nrleex)
c                       rli1+(rli0=rli1)*(1-(r/a)**nrliex),
c                       for electrons and ions
c
c                    3  input power spectrum, giving
c                       power/mode in watts. nrlpar=# of modes
c                       rlpara(j)=index of refraction of mode j
c                       rlpowr(j)=power in mode j
c                       j=1,nrlpar
c                       the radial zone which contains the linear mode
c                       conversion point is calculated self-
c                       consistently from the density,temp,b profiles
c                       for each mode at each timestep
c                       similarly the amount of power absorbed by
c                       the electrons via landau damping is
c                       calculated for each mode in each radial zone.
c                       if  nrldmp = 1, only the lmc model is used,
c                                       and energy goes only to the
c                                       ions at the zone coresponding
c                                       to the lmc points for each mode.
c                                 = 2, only electron landau damping is
c                                       considered, and all the energy
c                                       is absorbed by the electrons.
c                                 = 3, both loss mechanisms used.
c                                       (default)
c                                       the wave loses energy to the
c                                       electrons via landau damping
c                                       until it gets to the lmc point,
c                                       where it loses the rest
c                                       to the ions.
c
c
c
c
        go to (200,300,400), nrlmod
c
  200   continue
c
c       tabular input for power/zone
c               rleaux(j)=power to electrons in zone j
c               rliaux(j)=power to ions      in zone j
c
        zetot=0.0
        zitot=0.0
        zlhetot=0.0
        do 210 j=lcentr,mzones
          rleaux(j)  = rleprf(j-lcentr+1)
          rliaux(j)  = rliprf(j-lcentr+1)
          rlheaux(j) = rlheprf(j-lcentr+1)
          zetot=zetot+rleaux(j)*dx2i(j)
          zitot=zitot+rliaux(j)*dx2i(j)
          zlhetot=zlhetot+rlheaux(j)*dx2i(j)
  210   continue
        if ( zlhetot .eq. 0 ) zlhetot=1.e00
c
c..interpolate rlepwr and rlipwr as a function of bdtime
c  to rlepwrt and rlipwrt scalars at time tai
c
        zt = tai * uist
c
        if ( nrlepwr .gt. 1 ) then
          call timint (zt, rlepwrt, bdtime, nrlepwr, rlepwr, 1, 1)
        elseif ( nrlepwr .eq. 1 ) then
          rlepwrt = rlepwr(1)
        else
          rlepwrt = 0.0
        endif
c
        if ( nrlipwr .gt. 1 ) then
          call timint (zt, rlipwrt, bdtime, nrlipwr, rlipwr, 1, 1)
        elseif ( nrlipwr .eq. 1 ) then
          rlipwrt = rlipwr(1)
        else
          rlipwrt = 0.0
        endif
c
        zenorm = 1.e+7 * rlepwrt / zetot
        zinorm = 1.e+7 * rlipwrt / zitot
        zlhenorm = 1.e+7 * rlhepwr / zlhetot
c
c       normalize rleaux,rliaux in terms of ergs/(cm**3*sec)
c       (adjust indexing)
c
        do 220 j=mzones,lcentr,-1
          zlhe=zlhenorm*rlheaux(j)*ztovol
          rleaux(j)=zenorm*rleaux(j)*ztovol + zlhe
          rliaux(j)=zinorm*rliaux(j)*ztovol
  220   continue
        go to 400
c
  300   continue
c
c               polynomial fit
c
        zetot=0.0
        zitot=0.0
        do 310 j=lcentr,mzones
          zwefit(j)=rle1+(rle0-rle1)*(1.-xzoni(j)**nrleex)
          zwifit(j)=rli1+(rli0-rli1)*(1.-xzoni(j)**nrliex)
          zetot=zetot+zwefit(j)*dx2i(j)
          zitot=zitot+zwifit(j)*dx2i(j)
  310   continue
c
c..interpolate rlepwr and rlipwr as a function of bdtime
c  to rlepwrt and rlipwrt scalars at time tai
c
        zt = tai * uist
        call timint (zt, rlepwrt, bdtime, 20, rlepwr, 1, 1)
        call timint (zt, rlipwrt, bdtime, 20, rlipwr, 1, 1)
c
        zenorm = 1.e+7 * rlepwrt / zetot
        zinorm = 1.e+7 * rlipwrt / zitot
c
c       load rleaux,rliaux in terms of ergs/(cm**3*sec)
c
        do 320 j=lcentr,mzones
          rleaux(j)=zenorm*zwefit(j)*ztovol
          rliaux(j)=zinorm*zwifit(j)*ztovol
  320   continue
c
  400   continue
c
c
 1000   continue
c
c       on-off logic
c
        if(tai.lt.zhton)                   go to 3000
        if(tai.ge.zhtoff)                  go to 3000
c
c
 2000   continue
c
c       load heat sources per radial zone
c
        if(nrlmod.eq.3) go to 2020
c
        do 2010 jz=lcentr,ledge
          weauxs(jz)=rleaux(jz)
          wiauxs(jz)=rliaux(jz)
 2010   continue
c
        go to 8000
c
 2020   continue
c
c
c
c               self-consistent model
c
        zabs(mzones)=0.0
        zatten(1)=1.0
        zleft(mzones)=1.0
        ztot=0.0
        call resetr(zdn,i2,0.0)
        call resetr(weauxs,i1,0.0)
        call resetr(wiauxs,i1,0.0)
c
c       find linear mode conversion point for each
c       wavenumber rlpara in the input spectrum
c
        do 2100 j=1,nrlpar
        do 2030 jz=lcentr,mzones
c
c               compute dimensionless quantities
c
        zghz=rlfreq*1.e-9
        tikev=tis(2,jz)*evsinv*1.e-3
        teev=tes(2,jz)*evsinv
        tekev=teev*1.e-3
        bgauss=sqrt(bpols(2,jz)*bpols(2,jz)+bzs*bzs)
        btesla=bgauss*1.0e-04
c
c               compute linear mode conversion density znelmc
c
        zfct=2.35*(zghz/btesla)**2
        znelmc=2.27e+13*(zghz)**2/
     1  (0.5-zfct+7.75e-2*rlpara(j)*sqrt(tikev+tekev*zfct**2))
c               define zdn=ne-nelmc
c       zdn changes sign at the zone containing the lmc point
c
        zdn(jz)=rhoels(2,jz)-znelmc
        if(zdn(jz-1)*zdn(jz).lt.0.)
     1  ilmc(j)=jz
c
c       compute electron landau damping losses
c
c               use ful quantities
c
c       z1 = c/(nz*ve*sqrt(2.))
c
        z1=5.04e2/(rlpara(j)*sqrt(teev))
c
c       z2=rlfreq/wlh(in hz)
c
        zwpihz=210.*sqrt(rhoels(2,jz)/float(iabs(ngas(1))))
        zwrat=3.21e-3*sqrt(rhoels(2,jz))/bgauss
        zwlhhz=zwpihz/sqrt(1.+zwrat**2)
        z2=rlfreq/zwlhhz
c
c       calcaulate zkprim, attenuation "k" for e-landau power absorb.
c
c       if rlfreq/zwlh<1.00001, set zprim=20., which sets zatten<<1.
c               (i.e. as w->wlh, wave is almost totally attenuated)
c       if rlfreq/zwlh> 1.00001, calculate zkprim using formula
c
        if(z2.gt.1.00001.and.z1.lt.10.)
     x  zkprim=sqrt(fcpi*1836.*float(iabs(ngas(1))))*
     1  (rlfreq*rlpara(j)/3.e10)*(z1**3/sqrt(z2**2-1.))*
     2  exp(-z1**2)
        if(z2.gt.1.00001.and.z1.ge.10.) zkprim=0.
        if(z2.le.1.00001) zkprim=20.
c
c       calculate zatten, attenuation factor for each radial zone
c
        zatten(jz)=exp(-2.*zkprim * ( ahalfs(jz+1,1) - ahalfs(jz,1) ) )
c
 2030   continue
c
c       calcaulate the amount of energy zleft(jz), entering zone
c       jz. the amount of power absorbed is
c       zabs(jz). note that the wave propagates inward in minor radius.
c
        do 2040 jb=mzones,lcentr,-1
        zleft(jb-1)=zleft(jb)*zatten(jb-1)
        zabs(jb-1)=zleft(jb)-zleft(jb-1)
 2040   continue
c
c       select power deposition model
c        nrldmp = 1, lmc only
c                2, e-landau damping only
c                3, both
c
        go to (2048,2050,2070), nrldmp
c
 2048   continue
c
c       lmc only
c
c               deposit the energy in the zone
c               containing lmc point
c
        iilmc=ilmc(j)
        wiauxs(iilmc)=rlpowr(j)+wiauxs(iilmc)
        go to 2090
c
 2050   continue
c
c       electron landau damping only
c
        do 2060 jz=lcentr,mzones
        weauxs(jz)=zabs(jz)*rlpowr(j)+weauxs(jz)
 2060   continue
        go to 2090
c
 2070   continue
c
c       both lmc and e-landau damping losses
c
        do 2080 jz=lcentr,mzones
        if(jz.eq.ilmc(j)) wiauxs(jz)=zleft(jz)*rlpowr(j)+wiauxs(jz)
        if(jz.gt.ilmc(j)) weauxs(jz)=zabs(jz)*rlpowr(j)+weauxs(jz)
 2080   continue
 2090   continue
        ztot=rlpowr(j)+ztot
        if(nstep.le.ncheck) go to 2098
        ii=j
        write(nprint,9100) nstep, ii, ilmc(ii), ztot, rlipwrt
        call rarray(' zdn(jz)',zdn,mzones)
        call rarray(' zatten ',zatten,mzones)
        call rarray(' zleft  ',zleft,mzones)
        call rarray('zabs(jz)',zabs,mzones)
 2098   continue
 2100   continue
        znorm=1.e+7 * rlipwrt / ztot
c
c       normalize the power input per zone in terms of ergs/(cm**3*sec)
c
        do 2110 jz=lcentr,mzones
        weauxs(jz)=weauxs(jz)*znorm*zvols(jz)
        wiauxs(jz)=wiauxs(jz)*znorm*zvols(jz)
 2110   continue
c
        go to 8000
c
 3000   continue
c
c       set heat sources to zero in each radial zone
c
        do 3010 jz=lcentr,mzones
        weauxs(jz)=0.0
        wiauxs(jz)=0.0
 3010   continue
c
 8000   continue
        return
c
 8100   continue
        call error_olymp(1,iclass,isub,2,
     .  '** error ** off-time not defined for aux heating')
c
        return
c**********************************************************************
c
        entry hetprt
c
c       edit print-out of auxiliary heating terms
c
        if(cfutz(iauxht).le.epslon) return
        if(tai.lt.zhton .or. tai.gt.zhtoff) return
c
        lpage = lpage + 1
c
        z0=uist*1.e+03
        zt=tai*z0
        zdt=dtoldi*z0
c
        write (nprint,9010) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
        write(nprint,9020)
c
        zextot=0.0
        zixtot=0.0
        do 8400 i=lcentr,mzones
          j=i
          ik=min0(i,i3)
          zrad = ahalfs(i,2)
          z0=1.0e-07/zvols(i)
          zepowr(i)=z0*weauxs(i)
          zipowr(i)=z0*wiauxs(i)
          zextot=zepowr(i)+zextot
          zixtot=zipowr(i)+zixtot
          zzzz1(i)=0.0
          zzzz2(i)=0.0
          zzzz3(i)=0.0
          kzzz1(i)=0
          if (i .le. nrlpar) then
            zzzz1(i)=rlpara(ik)
            zzzz2(i)=rlpowr(ik)
            klmc=ilmc(ik)
            kzzz1(i)=klmc
            zzzz3(i) = ahalfs(klmc,2)
            endif
c
        write(nprint,9030) j,zrad,weauxs(i),wiauxs(i),zepowr(i),
     &                     zipowr(i),ik,zzzz1(i),zzzz2(i),kzzz1(i),
     &                     zzzz3(i)
 8400   continue
c
        write(nprint,9040)
        write(nprint,9050) zextot,zixtot
c
        return
c
 9010   format(1h1,2x,a48,10x,a72//
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.'/
     3  /35x,' quantities having to do with auxiliary heating')
 9020 format(1x,'j',2x,' radius',2x,'weauxs(j)',1x,'wiauxs(j)',1x,
     & 'elec heat',1x,' ion heat',1x,'*******',1x,'mode index',1x,
     & 'wave numr',1x,' power in',2x,'jz at lmc',2x,'rad at lmc'/)
 9030 format(i3,1x,f7.2,1x,4(1x,1pe9.2),14x,i3,2x,2(1x,1pe9.2),6x,i3,
     & 6x,0pf7.2)
 9040 format(/1x,' "elec heat" and "ion heat" are the watts per radial'
     & ,' zone, while "weauxs" and "wiauxs" are ergs/(cm**3*sec)'/)
 9050 format(1x,' total power being given the electrons =',1pe9.2,
     & 'watts;',2x,' total power being given to the ions =',1pe9.2,
     & 'watts')
 9100 format(/1x,' nstep=',i4,3x,' mode index=',i3,3x,' ilmc(i)=',
     & i3,3x,' ztot=',1pe9.2,3x,' rlipwrt=',1pe9.2)
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@he3
c       aes 15-apr-82 '240.', not '240', in if statement below do 260
c       dhfz 26-apr-79 fix ti>80. and 240. bugs
c       dhfz 4-4-79 added fits for ti>80  kev
c       dhfz 20-mar-79 first writeup
c-----------------------------------------------------------
c
        subroutine he3(k)
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
        dimension zene(8),zamu(8),zdens(8)
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        common/comtst/zll(8,55),zione(8,55),zelece(8,55)
c+++++++++++++++++++++++++++++++++++++++++++++++++++
        data   iclass/2/,    isub/19/
c
        data zene/3520.,3030.,820.,3670.,14670.,1010.,0.,0./,
     1          zamu/4.,1.,3.,4.,1.,3.,2.,3./
c
        if(.not.nlomt2(isub)) go to 10
        return
c
 10     continue
c
        go to (100,200,300),k
c
                                call expert(iclass,isub,1)
c
 100    continue
c
        return
c
 200    continue
c
        ihyd = 0
        ideut = 0
        ihe3 = 0
        ihe4 = 0
c
c
        do 205 ih = 1,mhyd
        if(ngas(ih).eq.1) ihyd = ih
        if(ngas(ih).eq.-2) ideut = ih
 205    continue
c
c
        if(mimp.eq.0) go to 215
        do 210 ii = 1,mimp
        if(nimp(ii).eq.-4) ihe3 = ii
        if(nimp(ii).eq.2) ihe4 = ii
 210    continue
c
 215    continue
c
        if(ideut*ihe3.eq.0) go to 299
c
        zones = float(nzones)
        adds = 0.
        ztest=tis(2,2)*useh
c
        do 260 jz = lcentr,mzones
c
        zt = tis(2,jz)*useh
        if(ztest.gt.240.) go to 9000
        if(ztest.gt.80.) go to 220
c
c
c               reaction rate fits from "convenient computational
c               forms for maxwellian reactivities",l.hively,
c               nuc. fus. 1 7 4 (1977)
c
        zt3 = zt**.3333333
        zt32 = zt3*zt3
c
c               zrho1s=sigma-v-bar for
c       d + t -> he4(3.52mev) + n(14.06mev)
c
        za1 = -21.377692
        za2 = -25.204054
        za3 = -7.1013427e-2
        za4 = 1.9375451e-4
        za5 = 4.9246592e-6
        za6 = -3.9836572e-8
        zr = .2935
c
        zrho1s =exp(za1/zt**zr+za2+zt*(za3+zt*(za4+zt*(za5+zt*za6))))
c
c               zrho2s=sigma-v-bar for
c       d + he3 -> he4(3.67mev) + p(14.67mev)
c
        za1 = -27.764468
        za2 = -31.023898
        za3 = 2.7889999e-2
        za4 = -5.5321633e-4
        za5 = 3.0293927e-6
        za6 = -2.5233325e-9
        zr = .3597
c
        zrho2s =exp(za1/zt**zr+za2+zt*(za3+zt*(za4+zt*(za5+zt*za6))))
c
c               zrho3s=sigma-v-bar for
c       d + d -> t(1.01mev) + p(3.03mev)
c
        za1 = 2.0018602e-14
        za2 = 19.307336
        za3 = 5.7756259e-3
        zr = .94955669
c
        zrho3s = za1*(1.+za3*(zt**zr))*exp(-za2/zt3)/zt32
c
c               zrho4s=sigma-v-bar for
c       d + d -> he3(.82mev) + n(2.45mev)
c
        za1 = 2.721942e-14
        za2 = 19.795543
        za3 = 5.3891144e-3
        zr = .91723383
        zrho4s = za1*(1.+za3*(zt**zr))*exp(-za2/zt3)/zt32
c
        go to 240
 220    continue
c
c               fits for t>80. kev
c
        zlog=log(zt)
        zlog2=zlog*zlog
c
c               zrho1s
c
        za=-.187196
        zb=1.44664
        zc=-37.4296
c
        zrho1s=exp(za*zlog2+zb*zlog+zc)
c
c               zrho2s
c
        za=-.525014
        zb=5.83727
        zc=-52.1
c
        zrho2s=exp(za*zlog2+zb*zlog+zc)
c
c               zrho3s
c
        za=-.120345
        zb=2.15468
        zc=-45.7268
c
        zrho3s=exp(za*zlog2+zb*zlog+zc)
c
c               zrho4s
c
        za=-.16099
        zb=2.55872
        zc=-46.5088
c
        zrho4s=exp(za*zlog2+zb*zlog+zc)
c
 240    continue
c
c               auxiliary terms for zlamda
c
        zden1s = 0.
        zden2s = 0.
        zmean = 0.
c
        if(ihyd.ne.0) zden1s = rhohs(ihyd,2,jz)
        if(ihe4.ne.0) zden2s = rhois(ihe4,2,jz)
        if(ihe4.ne.0) zmean = cmean(ihe4,2,jz)
c
        ztee = tes(2,jz)*useh
c
        z7 = .5*rhohs(ideut,2,jz)+1.3333*rhois(ihe3,2,jz)
     1          +zden1s + .25*zden2s*zmean**2
        z8 = (rhoels(2,jz)/z7)**.666667*.06757/ztee
c
c               shfus
c
        zene(7) = tis(2,jz)*useh
        zene(8) = zene(7)
c
        z11 = .5*rhohs(ideut,2,jz)*rhohs(ideut,2,jz)*zrho3s
        z12 = .5*rhohs(ideut,2,jz)*rhohs(ideut,2,jz)*zrho4s
        z13 = rhohs(ideut,2,jz)*rhois(ihe3,2,jz)*zrho2s
c
        zdens(1) = z11
        zdens(2) = z11
        zdens(3) = z12
        zdens(4) = z13
        zdens(5) = z13
        zdens(6) = z11
        zdens(7) = z13+3.*z11+2.*z12
        zdens(8) = z13
c
        shfus(ihyd,jz) = -z13-z11
        shfus(ideut,jz) = z11-z13+z12
        salfs(jz) = z13 + z11
        aoloss(jz) = z12 - z13
c
        wealfs(jz) = 0.
        wialfs(jz) = 0.
c
        do 250 ind = 1,8
c
c               compute zlamda
c
        zx = zene(ind)*z8/zamu(ind)**.666667
        z = sqrt(zx)
c
        z9 = log((1.+2.*z+z*z)/(1.-z+z*z))*.166667
        z10 = (atan((2.*z - 1.)/1.7321)+.5236)*.5774
        zl = (2./zx)*(z10-z9)
        zlamda = min(zl,1.)
c
c               compute wialfs,wealfs
c
        zions = zdens(ind)*zene(ind)*uesh*zlamda
c
        wialfs(jz) = wialfs(jz) + zions
c
        zelecs = zdens(ind)*zene(ind)*uesh*(1.-zlamda)
c
        wealfs(jz) = zelecs + wealfs(jz)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++
        zll(ind,jz) = zlamda
c
        zione(ind,jz) = zions*usee
c
        zelece(ind,jz) = zelecs*usee
c
        afuses(jz) = zrho2s
        aslows(jz) = zrho3s
        atfuse(jz) = zrho4s
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
 250    continue
c
c               adds :
c
        zterm = float(jz-2)+.5
        zvol = 2.0 * vols(mzones,1) * zterm/(zones*zones)
        zadds=(z11+z12)*zvol*dti*uist
        adds=adds+zadds
c
 260    continue
c
c               addtot:
c
        addtot  = addtot+adds
c
 299    continue
c
                                call expert(iclass,isub,2)
        return
c
 300    continue
c
                                call expert(iclass,isub,3)
c
        return
 9000   continue
c
        write(nprint,9010)
 9010   format(1x,36h*** ti > 240. kev - run is ended ***)
c
        nlend=.true.
        call endrun
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c
c       aes 29-oct-81
c**********************************************************************
c
        subroutine icrf(k)
c
      include 'cparm.m'
      include 'cbaldr.m'
c
c
c       this is a dummy version --does nothing but print warning.
c
        if (k.ne.1) return
        write(nprint,5)
        write(ntychl,5)
    5   format(1h1//10x,42h*** warning: dummy icrf routine loaded ***)
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@ecrh
c       aes 15-sep-81 allow 20 timezones (just after 240)
c       aes 20-nov-80 mzones+1 -->mzones, etc. after 110
c       aes 20-nov-80 changed nzones to mzones in initialization,
c               and reinitialize when ecrh is called with k=3
c       aes 4-nov-80 changed variable names to re,nre prefixes;
c               also weauxs-->weecrh, wiauxs-->wiecrh
c       aes 4-nov-80 changed subroutine name to ecrh
c      dhfz 3-oct-79 refine time control
c       dhfz 24-sept-79 second write-up
c       dhfz 26-jun-79 first writeup
c**********************************************************************c
c
        subroutine ecrh(k)
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
        dimension zb2(64),ip(55)
c----------------------------------------------------------------------
c
c    variables used:
c
c    nrenp          --- no. monte-carlo tests for area
c    rehe(it)    ---*total r-f power to electrons between time=reon(it)
c                    and time=reoff(it) (watts)
c    rehi(it)    ---*total r-f power to ions between time=reon(it)
c                    and time=reoff(it) (watts)
c    reoff(it)   ---*turn heating off at time=reoff(it) (sec)
c    reon(it)    ---*turn heating on at time=reon(it) (sec)
c    revols(jz)   --- area(window intersect zone jz)/area(zone jz)
c    rex1,rex2,rey1,rey2 ---*coords of rf window vertices are
c                     (rex1,rey1),(rex1,rey2),(rex2,rey1),(rex2,rey2)
c                     with rex1 < rex2, rey1 < rey2.
c                     the origin (0.,0.) corresponds to the center of
c                     the radial zones. (cm)
c
c    zelecs      --- total no. electrons in window volume
c    zint        --- time interpolation factor
c    zions       --- total no. ions in window volume
c    zrehe       --- average total rf heating to electrons
c    zrehi       --- average total rf heating to ions
c
c    output variables:
c
c    weecrh(jz)  --- r-f electron heating rate density in zone jz
c                    enters as +weecrh(jz) in diffusion eq'n.
c                    (ergs/cc*sec)
c    wiecrh(jz)  --- r-f ion r-f heating rate density in zone jz
c                    enters as +wiecrh in diffusion eq'n.
c                    (ergs/cc*sec)
c
c    * input in nurun1
c
c    this subroutine simulates r-f heating (in particular electron
c    synchotron heating) to compute the diffusion eq'n. terms +weecrh
c    and +wiecrh. we assume that rehe (resp. rehi) watts
c    are deposited in a rectangular window with the power deposited
c    in zone jz to the electrons (resp. ions) being proportional
c    to the number of electrons (resp. ions) in the intersection
c    of the window and zone jz divided by total number of electrons
c    (resp. ions) in the window.
c
c----------------------------------------------------------------------
 10     if(rex1.eq.rex2)return
        go to (100,200,100)k
c----------------------------------------------------------------------
c                                                      compute revols
 100    if(-rex1.ge.rmins.and.rex2.ge.rmins.and.-rey1.ge.rmins
     &   .and.rey2.ge.rmins)go to 180
        za=fcpi*rmins*rmins
        za1=(rex2-rex1)*(rey2-rey1)
        zq=za1/za
        nrenp=int(zq*3.e6)
        if(nrenp.lt.10000)nrenp=10000
        if(nrenp.gt.30000)nrenp=30000
        ipower=0
 105    ipower=ipower+1
        if(2**ipower.lt.mzones)go to 105
        i2p=2**ipower
c
        zrmin2=rmins*rmins
        do 110 i=1,mzones
 110    zb2(i)=xbouni(i+1)**2*zrmin2
        if(i2p.lt.mzones+1)go to 120
        do 115 i=mzones+1,i2p
 115    zb2(i)=zb2(mzones)
c
 120    do 150 jp=1,nrenp
 125    zx=(rex2-rex1)*ranz()+rex1
        zy=(rey2-rey1)*ranz()+rey1
        zr2=zx*zx+zy*zy
        if(zr2.gt.zb2(mzones))go to 125
        iz=0
        isgn=1
        is=1
 135    iz=iz+isgn*2**(ipower-is)
        if(zr2.ge.zb2(iz).and.zr2.le.zb2(iz+1))go to 140
        isgn=1
        if(zr2.le.zb2(iz))isgn=-1
        is=is+1
        go to 135
 140    ip(iz)=ip(iz)+1
 150    continue
c
 160    zp=float(nrenp)
        do 170 i=1,mzones
 170    revols(i)=float(ip(i))/zp
        return
c
 180    do 190 i=1,mzones
 190    revols(i)=dx2i(i+1)
        return
c----------------------------------------------------------------------
c                                                compute weecrh, wiecrh
c
 200    do 205 jz=1,mzones
        weecrh(jz)=0.
 205    wiecrh(jz)=0.
        zrehe=0.
        zrehi=0.
        zt1=tai*10.
        zt2=tbi*10.
        zdt=zt2-zt1
        it=1
 210    zreon=reon(it)*1.e3
        zreoff=reoff(it)*1.e3
        if(zreon.ge.zt2.or.zreoff.le.zt1)go to 240
        if(zreon.ge.zt1)go to 220
        if(zt2.ge.zreoff)go to 215
        zint=1.
        go to 230
 215    zint=(zreoff-zt1)/zdt
        go to 230
 220    if(zreoff.le.zt2)go to 225
        zint=(zt2-zreon)/zdt
        go to 230
 225    zint=(zreoff-zreon)/zdt
        go to 230
 230    zrehe=zrehe+zint*rehe(it)
        zrehi=zrehi+zint*rehi(it)
 240    it=it+1
        if(it.le.20) go to 210
        if(zrehe.eq.0..and.zrehi.eq.0.)return
c
 250    zions=0.
        zelecs=0.
 255    do 260 jz=lcentr,ledge
        zelecs=zelecs+rhoels(2,jz)*revols(jz-1)
 260    zions=zions+rhoins(2,jz)*revols(jz-1)
 265    do 275 jz=lcentr,ledge
        zvols=dx2i(jz) * 2.0 * vols(mzones,1)
        weecrh(jz)=
     &   zrehe*uesp*revols(jz-1)*rhoels(2,jz)/(zelecs*zvols)
 275    wiecrh(jz)=zrehi*uesp*revols(jz-1)*rhoins(2,jz)/(zions*zvols)
        weecrh(1)=weecrh(3)
        wiecrh(1)=wiecrh(3)
        weecrh(mzones)=weecrh(ledge)
        wiecrh(mzones)=wiecrh(ledge)
 290    return
        end
c
c**********************************************************************c
c@rprint
c       aes 15-sep-81 allow both ion and electron heating terms
c       aes 15-sep-81 activate subroutine
c**********************************************************************
        subroutine rprint(ktype,k)
c
      include 'cparm.m'
      include 'cbaldr.m'
c
        dimension ilabel(2,2)
        character *5 ilabel
        data ilabel/'elect','rons ','   io','ns   '/
c----------------------------------------------------------------------
c                                             print r-f heating summary
        if ( rex1 .eq. rex2 ) return
c
        write(nprint,8999)
 8999   format(1h1)
        itype=1
        if(rehe(1).le.epslon) go to 60
   10   continue
        write(nprint,9000)(ilabel(i,itype),i=1,2),rex1,rex2,rey1,rey2
 9000   format(1x,'ecrh heating of ',2a5,'in a window//15x',
     &   'rex1=',f8.2/15x,'rex2=',f8.2/15x,'rey1=',f8.2/15x,'rey2=',
     &   f8.2/1x,'with'/13x,'watts',10x,'time on',10x,'time off')
        it=1
   50   zpower=rehe(it)
        if(itype.eq.2)zpower=rehi(it)
        write(nprint,9001)zpower,reon(it),reoff(it)
 9001   format(10x,1pg10.3,8x,0pf7.2,10x,f7.2)
        it=it+1
        if(it.gt.20) go to 59
        if(reoff(it).gt.epslon)go to 50
   59   write(nprint,9002)
 9002   format(1x//)
c
   60   continue
        if(itype.eq.2) return
        if(rehi(1).le.epslon) return
        itype=2
        go to 10
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@cmpres  .../baldur/code/bald/dauxheat.f
c  rgb 14-aug-96 replaced fuzz with rndeps
c       drm 4-dec-84 add cfutz(41) to control exponent on temperature scaling
c       fgps 29-mar-83 arrays sltis, dnhhs, scroff, and the limit
c                      mxt1 now handled in cbaldr.
c       aes 14-jan-82 allow nfusn=3
c       aes 29-oct-81 array dimensions 52 -->55 in common/tempry/;
c               105 -->111 in zbouni,zzoni
c       aes 16-sep-81 added call to icrf(3)
c       aes 4-nov-80 added call to ecrh(3)
c       fgps 11-oct-79 introduced common/tempry/ --- with mxt1
c                      so as to limit the size of tcompi, redgi,
c                      and rmji when field-ripple effects are active
c       dhfz 5-apr-79 add call to he3
c       amck 18-may-78 fix line overflow in ecompi, wcompi compu.
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 2-jan-78 change zzone, zsmall to .8 and .2
c       amck 14-nov-77 add call alphas(3)
c       amck 30-oct-77 add call heat(3)
c       amck 15-aug-77 no double, add short formula
c       amck 15-aug-77 (temp) double precision
c       amck 13-aug-77 jz+1 -> mzones after #474
c       amck 12-aug-77 set zzoni, zbouni, in 354 loop
c       amck 10-aug-77 change beams(2) -> beams(3)
c       amck 3-aug-77 use adjusted value of rmini if compl adjusted
c       amck 27-jul-77 adjust compl so that outer zone is .gt. zsmall
c       amck 26-jul-77 if dvolx is roundoff, set to 0
c       amck 25-jul-77 add wcompi, fcompi
c       amck 22-jul-77 modify for xbouni = avg. of neighb. xzoni's
c       amck 2-nov-76 remove -1 in ir comp. after #304, add error check
c               in nzcomp, dvolx
c       amck 21-oct-76 add comment about values of comtim vars.
c       amck 3-sep-76 reset cmean at edge if limiter moves in
c       amck 27-aug-76 remove data round
c       amck 26-aug-76 change redgi to min rad. * sqrt(maj. rad.)
c       amck 26-aug-76 fix volume in deletion of plasma (#454)
c       amck 26-aug-76 add if(nzcomp.lt.0) before #468
c       amck 26-aug-76 fix zbint when zxrad0=zxrad1, zradn -> zxradn
c       amck 25-aug-76 pad out cmean, fix adding of xbouni's
c       amck 25-aug-76 fix use of zlimi
c       amck 24-aug-76 fix ze -> zce
c       amck 24-aug-76 remove fuzz, mxt data
c       amck 24-aug-76 add ion energy compresssion
c       amck 24-aug-76 finish generating file
c       amck 19-aug-76 generate file
c**********************************************************************c
c
c
        subroutine cmpres
c
c
cl      2.13    compress the plasma
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c
c       double precision zbint, zx4, zx2, zx22, zxrad0, zxrad1,
c     d z0
        dimension
     r   zbint(3)     , zbouni(111)   , zzoni(111)
c
c
c                               zzone must be .lt.1 and .gt. 0.5
        data    zzone /0.8/,    zsmall /0.20/,  zcross/1.03/
c
c
c------------------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /13/
c
c
        if (.not.nlomt2(isub)) go to 10
      if (nstep .lt. 2)
     &  call mesage(' *** 2.13 subroutine cmpres bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
c
cl      common blocks and variables modified:
c
c       bedgi, compl, compp, dvolx, nlcomp, nzcomp (comtim)
c       bpoli, chi (comsta), rmaji, rmini (commsh), bzi (comtok)
c       acompi, ecompi, ebtoti, totali (comcns),
c       mzones, ledge (comflg)
c
c------------------------------------------------------------------------------
c
c
c       this subroutine performs slow compression by a series
c       of instantaneous compressions.
c
c       the desired major radius is computed from tcompi and rmji.
c       the compression for this timestep (compp) is the present major
c       radius (rmaji) divided by the desired major radius.
c       the zones (and boundaries) are compressed along the minor
c       radius by a factor compp**1/2.
c       bz is increased by a factor compp, b-poloidal by compp**3/2,
c       hence current by a factor compp, volume shrunk by compp**-2,
c       hence density increased by compp**2, energy increased by
c       compp**4/3, and energy density by compp**10/3.
c
c       the limiter radius is independently determined from
c       redgi (and tcompi).
c       if the limiter radius is larger than the compressed
c       rmini, zones are added, if the radius is smaller, they are
c       removed.
c
c
c------------------------------------------------------------------------------
c
c
c       if the plasma is:
c                       compressing     expanding       stationary
c       compp           .gt.1           .lt.1           .eq.1
c
c
c       if the limiter, relative to the plasma, is moving:
c                       in              out             not at all
c       compl           .gt.1           .lt.1           .eq.1
c       nzcomp          .le.0           .ge.0           .eq.0
c       dvolx           .ge.0           .ge.0           .eq.0
c
c       if either compp.ne.1 or compl.ne.1, nlcomp = .true.
c
c
c------------------------------------------------------------------------------
c
c
cl      1)      set zrmaji, zrmini, compp, and compl
c
c
cl      1.1)    set zrmaji, zrmini
c
c
  100   continue
        nlcomp = .false.
c
        it = 1
c
        do 104 jt = 2, mxt1
        if (tcompi(jt).le.0.0) go to 106
        it = jt
        if (tcompi(jt).eq.tbi) go to 106
        if (tcompi(jt).gt.tbi) go to 108
  104   continue
c
  106   continue
        zrmaji = rmji(it)
        zrmini = sqrt(zrmaji) * redgi(it)
        go to 110
c
  108   continue
        z1 = (tcompi(it) - tbi) / (tcompi(it) - tcompi(it-1))
        z2 = (tbi - tcompi(it-1)) / (tcompi(it) - tcompi(it-1))
        zrmaji = rmji(it-1)*z1 + rmji(it)*z2
        zrmini = sqrt(zrmaji)*(z1*redgi(it-1) + z2*redgi(it))
c
  110   continue
c
        compp = rmaji / zrmaji
        compl = rmini / (zrmini * sqrt(compp))
c
        if (abs(compp-1.0) + abs(compl-1.0).gt.rndeps) nlcomp = .true.
                                        call expert(iclass,isub,1)
c
c
c               if no compression, return
c
c
        if (.not.nlcomp) return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
cl      2)      compress
c
c
  200   continue
        zcexp=4./3.
        if(cfutz(41).gt.0.) zcexp=cfutz(41)
        zce = compp**zcexp
        zcn = compp**2
        zcen = zce * zcn
        zcbp = sqrt(compp)**3
        zcbz = compp
        zca = 1.0 / sqrt(compp)
        zceb = compp
c
        do 210 jz = 1, ledge
c
        do 204 jp = lhyd1, limpn
        chi(jp,jz) = chi(jp,jz)*zcn
  204   continue
c
        chi(lelec,jz) = chi(lelec,jz)*zcen
        chi(lion,jz) = chi(lion,jz)*zcen
        bpoli(jz) = bpoli(jz)*zcbp
  210   continue
c
        rmini = rmini*zca
        rmaji = zrmaji
        bpoli(mzones) = bpoli(mzones)*zcbp
        bzi = bzi * zcbz
c
c               conservation variables
c
        zdti = 0.0
        if (tbi-tai.gt.0.0) zdti = 1.0 / (tbi - tai)
        call resetr(fcompi,mxchi,0.0)
        wcompi = 0.0
c
        acompi(lelec) = acompi(lelec) + totali(lelec)*(zce - 1.0)
        fcompi(lelec) = fcompi(lelec) + totali(lelec)*(zce - 1.0)*zdti
        acompi(lion) = acompi(lion) + totali(lion)*(zce - 1.0)
        fcompi(lion) = fcompi(lion) + totali(lion)*(zce - 1.0)*zdti
        totali(lelec) = totali(lelec)*zce
        totali(lion) = totali(lion)*zce
        ecompi = ecompi + ebtoti*(zceb - 1.0)
        wcompi = wcompi + ebtoti*(zceb - 1.0)*zdti
        ebtoti = ebtoti*zceb
                                        call expert(iclass,isub,2)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
cl      3)      add / subtract zones
c
c
c       the profiles are compressed along the minor radius in proportion
c       to the square root of the major radius compression.  if the
c       limiter moves in the same amount, i.e., remains stationary
c       w.r.t. the profiles of the plasma parameters, no zones need be
c       added or removed.
c       if, however, the limiter moves out w.r.t. the plasma profiles,
c       the outer zone must be expanded, new zones added if necessary,
c       and the plasma must be padded out with cool plasma.
c       similarly, if the limiter moves in w.r.t. the plasma profiles,
c       the outer zone must be trimmed, and zones removed if necessary,
c       and the plasma in those zones thrown away.
c       in either case, since the mesh proportions have been changed,
c       the mesh must be rescaled, and the other mesh variables
c       recomputed.
c
c
c
  300   continue
c
c               if limiter radius is proportional to sqrt(rmaji),
c               don't add/subtract plasma
c
        if (abs(compl-1.0).le.rndeps) go to 500
        if (compl.gt.1.0) go to 350
c
c       3.1)    limiter is moving out, relative to compression along
c               the minor radius
c
c
c       since the boundary zone(s) will become larger in comparison
c       with the inner zones (jz .lt. present ledge-1), but xbouni
c       is constant (1.0) at the edge,
c       xbouni for the inner boundaries must be compressed.
c       the compression of the mesh (rescaling of xbouni) times the
c       compression of the limiter will always be equal to the square
c       root of the major radius compression.
c
c       the region between r = zxrad0 * new rmini and
c       r = zxrad1 * new rmini is (in this case) the volume that has
c       added to the (former) outermost real zone.
c       dvolx is analogous to d2xi for that zone.
c       these variables are used to compute how much plasma and b
c       is to be added.
c
c       if the last zone, after either mesh compression or expansion,
c       is too small, compl is modified so that it will not exist.
c
c
c               compress mesh
c
c
        do 304 jz = 1, mzones
        zzoni(jz)  = xzoni(jz)*compl
        zbouni(jz) = xbouni(jz)*compl
  304   continue
c
c
c
c               add additional zones (if possible)
c
c
        zdummy = xbouni(mzones+1)
        i0 = mxzone - 1
c
        do 328 jz = ledge, i0
        ir = jz
        if (radius(2).gt.0.0) go to 308
        go to (310,312),nrfit
c
  308   continue
        i001 = jz+1-lcentr
        i002 = ledge-lcentr
        zzoni(jz) = zzoni(ledge-1)*radius(i001) /
     1                                  radius(i002)
        go to 326
c
  310   continue
        zzoni(jz) = zzoni(ledge-1)*(float(jz-lcentr) + 0.5) /
     1                                  (float(ledge-1-lcentr) + 0.5)
        go to 326
c
  312   continue
        zzoni(jz) = zzoni(ledge-1) *
     1  ((float(jz-lcentr) + 0.75) / (float(ledge-1-lcentr) + 0.75)) *
     2          sqrt(float(ledge-lcentr) / float(jz+1-lcentr))
c
  326   continue
        if (jz.le.ledge) go to 328
        zbouni(jz) = (zzoni(jz) + zzoni(jz-1)) * 0.5
        if (zbouni(jz).ge.(1.0-rndeps)) go to 330
  328   continue
c
        go to 9030
  330   continue
c
c               if last zone is too small, adjust compl so that
c               it will be compressed away, and try again
c
        if ((1.0-zbouni(ir-1)).gt.zsmall*(zbouni(ir-1)-zbouni(ir-2)))
     1          go to 335
        compl = compl / zbouni(ir-1)
        go to 300
c
c               last zone is not too small
c
  335   continue
c
        call copyr(zbouni,1,xbouni,1,ir)
        call copyr(zzoni,1,xzoni,1,ir)
c
        xzoni(ir-1) = min(xzoni(ir-1) ,
     1                  zzone + (1.0-zzone)*xbouni(ir-1))
        xzoni(ir) = 2.0 - xzoni(ir-1)
        xbouni(ir+1) = zdummy
        xbouni(ir) = 1.0
        nzcomp = ir - mzones
        zxrad1 = xbouni(mzones)
        zxrad0 = compl
        go to 370
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
c       3.2)    limiter is moving in, relative to compression along
c               minor radius
c
c
c       since the edge will become smaller in comparison with the inner
c       (zone) boundaries, but xbouni at the edge is constant (1.0),
c       xbouni for the inner boundaries must be expanded.
c       the compression of the limiter divided by the mesh expansion
c       factor (the reciprocal of the mesh compression factor) will
c       always be the square root of the major radius compression
c       factor.
c
c       in this case, the area between r = zxrad0 * new rmini and
c       r = zxrad1 * new rmini is the volume that has been (is to be)
c       trimmed from zone ir - 1 to make it into the new outermost
c       real zone.
c
c
c               expand mesh
c
  350   continue
        ir = lcentr
c
        do 354 jz = 1, mzones
        zzoni(jz) = xzoni(jz) * compl
        zbouni(jz) = xbouni(jz) * compl
        if (zbouni(jz).ge.(1.0-rndeps)) go to 356
        ir = jz + 1
  354   continue
c
        go to 9031
c
  356   continue
c
c               if last zone is too small, adjust compl so that
c               it will be compressed away, and try again
c
        if ((1.0-zbouni(ir-1)).gt.zsmall*(zbouni(ir-1)-zbouni(ir-2)))
     1          go to 360
        compl = compl / zbouni(ir-1)
        go to 300
c
c               last zone is not too small
c
  360   continue
c
        call copyr(zbouni,1,xbouni,1,ir)
        call copyr(zzoni,1,xzoni,1,ir)
c
        xzoni(ir-1) = min(xzoni(ir-1) ,
     1                  zzone + (1.0-zzone)*xbouni(ir-1))
        xzoni(ir) = 2.0 - xzoni(ir-1)
        zdx2i = dx2i(ir)
        zxrad1 = xbouni(ir)
        zxrad0 = 1.0
        xbouni(ir) = 1.0
        xbouni(ir+1) = xbouni(mzones+1)
        nzcomp = ir - mzones
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
c       3.3)    recompute secondary mesh variables
c
c
  370   continue
        if (abs(zxrad0-zxrad1).le.rndeps*zxrad0) zxrad0 = zxrad1
c
        ledge = ledge + nzcomp
        mzones = ledge + 1
        rmini = rmini / compl
c
        dxboui(1) = 0.0
        dxzoni(1) = xbouni(2) - xbouni(1)
c
        do 378 jz = lcentr, mzones
        dxzoni(jz) = xbouni(jz+1) - xbouni(jz)
        dxboui(jz) = xzoni(jz) - xzoni(jz-1)
        dx2i(jz)   = 0.5 * (xbouni(jz+1)**2 - xbouni(jz)**2)
        dx2inv(jz) = 1.0 / dx2i(jz)
c
        if (jz.gt.ledge) go to 378
        if (xbouni(jz+1).le.zcross*xbouni(jz)) go to 374
c
c               zone width > zcross*r, use complicated expression
c
c
        z0 = 0.0
        if (xbouni(jz).gt.epslon) z0 = log(xbouni(jz+1)/xbouni(jz))
        zx2 = xbouni(jz+1)**2 - xbouni(jz)**2
        zx22 = zx2**2
        zx4 = xbouni(jz+1)**4 - xbouni(jz)**4
c
        bint(1,jz) = 2.0*xbouni(jz)**2 / zx22 *
     1  (0.25*zx4 - zx2*xbouni(jz+1)**2 + xbouni(jz+1)**4 * z0)
c
        bint(2,jz) = 2.0 * xbouni(jz)*xbouni(jz+1) / zx22 *
     1  (0.5*zx4 - 2.0*xbouni(jz)**2 * xbouni(jz+1)**2 * z0)
c
        bint(3,jz) = 2.0 * xbouni(jz+1)**2 / zx22 *
     1  (0.25*zx4 - zx2*xbouni(jz)**2 + xbouni(jz)**4 * z0)
c
        go to 378
c
c               zone width .le. zcross*r, use limit expression
c
  374   continue
        zxm = xbouni(jz+1) - xbouni(jz)
        zxp = xbouni(jz+1) + xbouni(jz)
        bint(1,jz) = 2.666666666 * xbouni(jz)**2 * xbouni(jz+1) *
     1          zxm / zxp**2
        bint(2,jz) = 1.333333333 * xbouni(jz) * xbouni(jz+1) *
     1          zxm / zxp
        bint(3,jz) = 2.666666666 * xbouni(jz) * xbouni(jz+1)**2 *
     1          zxm / zxp**2
c
 378    continue
c
        dvolx = 0.5*(zxrad1**2 - zxrad0**2)
c
        z0 = 0.0
        zbint(1) = 0.0
        zbint(2) = 0.0
        zbint(3) = 0.0
c
        if (zxrad1.le.zcross*zxrad0) go to 384
c
c               long formula
c
        if (zxrad0.gt.epslon) z0 = log(zxrad1/zxrad0)
        zx2 = zxrad1**2 - zxrad0**2
        zx22 = zx2**2
        if (zx22.le.epslon) go to 388
c
        zx4 = zxrad1**4 - zxrad0**4
c
        zbint(1) = 2.0*zxrad0**2 / zx22 *
     1  (0.25*zx4 - zx2*zxrad1**2 + zxrad1**4 * z0)
c
        zbint(2) = 2.0 * zxrad0*zxrad1 / zx22 *
     1  (0.5*zx4 - 2.0*zxrad0**2 * zxrad1**2 * z0)
c
        zbint(3) = 2.0 * zxrad1**2 / zx22 *
     1  (0.25*zx4 - zx2*zxrad0**2 + zxrad0**4 * z0)
c
        go to 388
c
c               short formula
c
  384   continue
        zxm = zxrad1 - zxrad0
        zxp = zxrad1 + zxrad0
        zbint(1) = 2.666666666 * zxrad0**2 * zxrad1 *
     1          zxm / zxp**2
        zbint(2) = 1.333333333 * zxrad0 * zxrad1 *
     1          zxm / zxp
        zbint(3) = 2.666666666 * zxrad0 * zxrad1**2 *
     1          zxm / zxp**2
c
  388   continue
c
c
c
                                        call expert(iclass,isub,3)
        if (compl.gt.1.0) go to 450
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
cl      4.1)    add plasma
c
c       pad out former last zone
c
c
c       limiter is moving out w.r.t. the profiles.
c       plasma must be added, the density and energy density of the
c       added plasma being the same as the dummy zone.
c       the current in the added plasma is 0.
c
c
  400   continue
c
c               error checking
c
        if (nzcomp.lt.0.or.dvolx.lt.0.0) go to 9040
c
        iedge = ledge - nzcomp
        z2 = dvolx / dx2i(iedge)
        z1 = 1.0 - z2
        zvoli = vols(mzones,1) * usil**3
        zvols = zvoli * usid
        zebi = zvols * usie * uisb**2 / (8.0*fcpi)
        zdvol = dvolx*2.0 * zvoli
c
        do 404 jp = 1, mchi
        chi(jp,iedge) = chi(jp,iedge)*z1 + chi(jp,iedge+1)*z2
        acompi(jp) = acompi(jp) + zdvol*chi(jp,iedge+1)
        fcompi(jp) = fcompi(jp) + zdvol*chi(jp,iedge+1)*zdti
  404   continue
c
c               pad out cmean
c
        if (mimp.le.0) go to 408
        do 406 ji = 1, mimp
        cmean(ji,2,iedge) = cmean(ji,2,iedge)*z1 + cmean(ji,2,iedge+1)
     1                                                  *z2
  406   continue
c
  408   continue
c
        z0 = compl * bpoli(iedge+1)
        bpoli(iedge+1) = rndup * z0 / xbouni(iedge+1)
        ecompi = ecompi + zebi * 2.0 * z0**2 * log(zxrad1/zxrad0)
        wcompi = wcompi + zebi * 2.0 * z0**2 * log(zxrad1/zxrad0)*zdti
        if (nzcomp.le.0) go to 499
c
c
c               pad out added zones
c
c
        i0 = mzones - nzcomp
c
        do 420 jz = i0, ledge
c
c               pad out cmean
c
        if (mimp.le.0) go to 412
        do 410 ji = 1, mimp
        cmean(ji,2,jz+1) = cmean(ji,2,jz)
  410   continue
c
  412   continue
c
c               pad out chi
c
        do 418 jp = 1, mchi
        chi(jp,jz+1) = chi(jp,jz)
        acompi(jp) = acompi(jp) + chi(jp,jz)*2.0*dx2i(jz)*zvoli
        fcompi(jp) = fcompi(jp) + chi(jp,jz)*2.0*dx2i(jz)*zvoli*zdti
  418   continue
c
        bpoli(jz+1) = bpoli(jz)*xbouni(jz) / xbouni(jz+1)
        ecompi = ecompi + zebi*(bint(1,jz)*bpoli(jz)**2 +
     1  bint(2,jz)*bpoli(jz)*bpoli(jz+1) + bint(3,jz)*bpoli(jz+1)**2)
        wcompi = wcompi + zebi*(bint(1,jz)*bpoli(jz)**2 +
     1  bint(2,jz)*bpoli(jz)*bpoli(jz+1) + bint(3,jz)*bpoli(jz+1)**2)*
     2  zdti
  420   continue
c
        bpoli(mzones+1) = bpoli(mzones)/ xbouni(mzones+1)
        go to 499
c
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c       4.2)    remove plasma
c
c
c       the limiter is moving in w.r.t. the profiles.
c       the plasma in that region is lost, as is the current.
c
c
  450   continue
c
c               error checking
c
        if (nzcomp.gt.0.or.dvolx.lt.0.0) go to 9041
c
        zvoli = vols(mzones,1) * usil**3
        zebi = zvoli*usid * usie * uisb**2 / (8.0*fcpi)
        zdvolx = zvoli * 2.0 * dvolx
        z0 = zvoli * 2.0 * compl**2
        iedge = ledge - nzcomp
c
        if (nzcomp.ge.0) go to 460
c
c
c               delete plasma in zones being removed
c
c
        do 456 jz = mzones, iedge
        if (jz.eq.mzones) z1 = z0 * zdx2i
        if (jz.gt.mzones) z1 = z0 * dx2i(jz)
c
        do 454 jp = 1, mchi
        acompi(jp) = acompi(jp) - z1*chi(jp,jz)
        fcompi(jp) = fcompi(jp) - z1*chi(jp,jz)*zdti
        chi(jp,jz) = 0.0
  454   continue
c
        ecompi = ecompi - zebi*(bint(1,jz)*bpoli(jz)**2 + bint(2,jz)*
     1          bpoli(jz)*bpoli(jz+1) + bint(3,jz)*bpoli(jz+1)**2)
        wcompi = wcompi - zebi*(bint(1,jz)*bpoli(jz)**2 + bint(2,jz)*
     1          bpoli(jz)*bpoli(jz+1) + bint(3,jz)*bpoli(jz+1)**2)*zdti
        bpoli(jz+1) = 0.0
  456   continue
c
  460   continue
c
c
c               delete plasma in last remaining zone (perhaps partial)
c               set dummy zone values
c
c
        do 468 jp = 1, mchi
        acompi(jp) = acompi(jp) - zdvolx*chi(jp,ledge)
        fcompi(jp) = fcompi(jp) - zdvolx*chi(jp,ledge)*zdti
        chi(jp,ledge+1) = chi(jp,iedge+1)
        if (nzcomp.lt.0) chi(jp,iedge+1) = 0.0
  468   continue
c
        if (mimp.le.0) go to 474
c
        do 472 ji = 1, mimp
        cmean(ji,2,ledge+1) = cmean(ji,2,iedge+1)
  472   continue
c
  474   continue
c
        z0 = ((zxrad1**2 - 1.0)*xbouni(mzones-1)*bpoli(mzones-1) +
     1          (1.0 - xbouni(mzones-1)**2)*zxrad1*bpoli(mzones))
     2          / (zxrad1**2 - xbouni(mzones-1)**2)
        ecompi = ecompi - zebi*(zbint(1)*z0**2
     1  + zbint(2)*bpoli(mzones)*z0 + zbint(3)*bpoli(mzones)**2)
        wcompi = wcompi - zebi*(zbint(1)*z0**2
     1  + zbint(2)*bpoli(mzones)*z0 + zbint(3)*bpoli(mzones)**2)*zdti
        bpoli(mzones) = z0
        bpoli(mzones+1) = z0 / xbouni(mzones+1)
c
  499   continue
        bint(1,mzones) = 0.0
        bint(2,mzones) = 0.0
        bint(3,mzones) = 0.0
                                        call expert(iclass,isub,4)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
cl      5)      call other routines which may need to compress
c
c
  500   continue
c
c               neutral beam routines
c
        call beams(3)
c
c               auxiliary heating
c
        call heat(3)
        call ecrh(3)
        call icrf(3)
        if(nfusn.ne.2) call alphas(3)
        if(nfusn.ne.1) call he3(3)
c
c
                                        call expert(iclass,isub,5)
c
        return
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
cl      90)     error handler
c
c
 9030   continue
        call error_olymp(0,iclass,isub,3,
     1          'cannot compress mesh ')
c
 9031   continue
        call error_olymp(0,iclass,isub,3,
     1          'cannot expand mesh ')
 9040   continue
        call error_olymp(0,iclass,isub,4,
     1          'nzcomp or dvolx .lt. 0 ')
 9041   continue
        call error_olymp(0,iclass,isub,4,
     1          'nzcomp .gt. 0  or dvolx .lt. 0 ')
        end

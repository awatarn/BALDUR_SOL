c--------1---------2---------3---------4---------5---------6---------7-c
c@resolv  .../baldur/code/bald/resolv.f
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rgb 28-aug-96 corrected call ivar('timestep,nstep')
c  rgb 24-sep-92 changed calls to harray, each element is a character string
c  rgb 09-jun-89 added label 221 to return to avoid warning
c       dps 07-aug-89 15.12 theta-p center fluxes for IRE code
c       dps 18-jul-89 15.11 remove impurity density check for diagnostic mode
c       dps 07-jul-89 15.11 add neutral impurities to IRE code
c       dps 15-may-89 15.09 insert IRE code into predictor-corrector loop
c       dps 24-aug-88 15.00 skip checks on impurity densities when
c                     using NC code
c       rgb 4-jun-85 do 228 jz = 2, mzones to avoid bpoli(1)
c       fgps 1-apr-83 made sure that all radial zones of index j .gt.
c                     abs(cfutz(izlimt=126)) are always removed from the
c                     timestep control; while zones with  cfutz(iinner)
c                     < j < cfutz(iouter) are only removed from timestep
c                     control when  cfutz(ixdel1) < time < cfutz(ixdel2).
cdoc
c=======================================================================
c
c       ------------
c       sbrtn RESOLV   file DSOLVER
c       ------------
c
c
c      2.11    timestep control
c
c-----------------------------------------------------------------------
c
c
c      common blocks and variables modified:
c
c       bpoli, chi (comsta),
c       bpols, bzs, rmins, rmajs, rhoels, rhoins, rhohs, rhois, tes,
c               tis (comdf2),
c       tai, tbi, dti, nlrpet, nliter (comflg),
c       comext
c
c-----------------------------------------------------------------------
c
c
c       basic notion of extrapolation
c
c       we assume that, if we compute the solution of a partial
c       differential equation by finite timesteps of length dt,
c       the computed value of the solution at a certain time
c       "t" will be a function "f" of dt,
c       and that f(0) will be the ideal solution.
c       if f is taylor expanded, it is:
c
c       f(dt) = f(0) + a1*dt + a2*(dt**2) + ...
c
c       extrapolation is based on the idea of computing f
c       at t with various dt's,  fitting a polynomial to
c       the results (i.e., assuming the taylor expansion is finite),
c       and using the polynomial to guess f(0). (dt = 0).
c
c       this version of this subroutine does a first-order polynomial
c       fit, with  "dti" and 0.5*dti as the two dt's.
c       this means that f(dti) - f(0.5*dti) = 0.5*a1*dti,
c       so
c
c       f(dti) - 2*(f(dti) - f(.5*dti)) = 2*f(.5*dti) - f(dti) = f(0)
c
c
c-----------------------------------------------------------------------
c
c
c       the extrapolation is done in four stages:
c
c       first, (k = 1) the beginning of timestep values of chi and bpoli
c       are saved in xbpol1 and xchi1, as are dti (xdti1),
c       tai (xtai1), tbi (xtbi1).  this is usually done at the
c       end of the previous timestep, but if ilinit = .false.,
c       indicating that they have not been set, or tai is not equal to
c       xtai1, indicating that the values are not up to date,
c       they are set again.
c       second, (k = 2, liter = 1,) after the timestep has been computed
c       with a timestep length of dti, chi and bpoli are reset from
c       xchi1 and xbpol1, and the end-of-timestep values saved
c       in xchi1 and xbpol1.  then nlrpet is set, to cause
c       chi and bpoli to be recomputed without recomputing
c       coefficients, and dti is set to 1/2 its former value.
c       note that this repetition is part of a loop
c       within stepon.
c       third, (k=2,liter=2) after one timestep at 1/2 dti, another timestep
c       is computed, but this time, nliter is set instead of nlrpet,
c       since the coefficients must be recomputed.  the intermediate
c       values of chi and bpoli are not used.
c       fourth, (k=2, liter=3), after we have computed chi(tbi) and bpoli(tbi)
c       in timesteps of 1/2 dti, we use the extrapolation
c       formula above to compute the final values of chi and bpoli.
c       in addition, an error figure is computed by taking the
c       maximum of abs( f(0.5*dti) - f(dti) ) / f(0)
c       (in the terminology of the preceding section)for
c       f representing chi and bpoli for each species and
c       each zone.  if the value is 0/0, it is ignored, otherwise, if
c       f(0) is 0 or negative, it is an error.  the error figure
c       is called xerr, and lperr and lzerr are the parameter and zone
c       indices, respectively, where the error was that value.
c       this error figure, divided by "errmax", the maximum desired error
c       is used to adjust dti for the next timestep.
c       note that dti must lie between dtmini and dtmaxi.
c
c
c       if nlextr is .true., extrapolation is done; otherwise,
c       not.
c
c
c-----------------------------------------------------------------------
c
c
c       outside of the extrapolation loop, a check is made for total
c       change over one timestep, that is (chi(end of t.s.,including
c       extrapolation) - chi(b.of t.s.) ) / chi(end of t.s.,...)
c       this is called "xdel".  xdel/delmax is also used to adjust dti.
c       this check is performed outside of, and independent of,
c       extrapolation (although the change will, of course, depend
c       on whether extrapolation is done).
c
c       if 1) any of the main parameters goes negative or zero,
c       or 2) xdel / delmax is greater than a specified (.gt.1) value
c       (xfutz(idell)), or 3) xerr / errmax is greater than another
c       specified value (xfutz(ierrl)), the timestep is repeated
c       (after statement no. 300), with dt reduced by a specified
c       factor (xfutz(ineg)).
c
c
c
c-----------------------------------------------------------------------
c
c
c       each extrapolation minor-timestep will be "predictor-corrector-ed",
c       if thetap .gt. 0.  predictor-corrector is a single-step
c       iteration, where the timestep is done once with the coefficients
c       based on the beginning of (minor) timestep values, (trial step)
c       then again with the coefficients computed from an average
c       (depending on thetap) of the trial results and the beginning-of-
c       minor-timestep values. (final results)
c       note that each minor timestep for extrapolation is
c       a whole predictor-corrector timestep.
c
c       an option exists to exclude certain prescribed zones from the
c       timestep control (xdel/delmax). there is an on-time for this
c       exclusion and an off-time:  cfutz(ixdel1)="on-time" in secs;
c       cfutz(ixdel2)="off-time" in secs.  when "on", exclusion can be
c       extended over two regions:
c        cfutz(izlimt)=index of radial zone at and beyond which time-
c                      step controls are excluded.  this defines the
c                      first region and default moves it completely
c                      outside the plasma.  moreover, by reprogramming
c                      (1-apr-83), zones in this region are now exclud-
c                      ed from the timestep control for all time.
c        cfutz(iinner)=inner index of second region of exclusion.  if
c                      cfutz(iinner).le.epslon, the second region is
c                      moved outside the plasma; which is the default
c                      condition.
c        cfutz(iouter)=outer radial index of the second region of ex-
c                      clusion from the timestep control.  however, un-
c                      less cfutz(iouter) .gt. cfutz(iinner), the outer
c                      boundary of the second region is located (by de-
c                      fault) outside the plasma.
c
c
cend
c
c-----------------------------------------------------------------------
c
        subroutine resolv(k)
c
        include 'cparm.m'
      include 'cbaldr.m'
        include 'comncr.m'
c
c
        logical ilinit,lcheck
        dimension jlimit(5)
c
c
        data    ilinit /.false./, lcheck /.false./,
     i  idamp/1/, ilim/2/, ierrl/3/, idell/4/, ineg/5/, ihigh/6/
c
c
        data izlimt,ixdel1,ixdel2,iinner,iouter /126,161,162,163,164/
c
c-----------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /11/
c
c
        if (.not.nlomt2(isub)) go to 10
        call mesage('% 2.11 subroutine resolv bypassed')
        return
   10   continue
c
c
c      1.1)    set x----- variables
c
c
        if (k.gt.2.or.k.le.0) go to 9000
        if (k.eq.2) go to 200
c
c  For k = 1: First call to RESOLV; done at top of STEPON
c  before call to COEF.
c
  100   continue
        i1 = mxchi * mzones
c
c
c
        call copyr(aflxii,1,xaii1,1,mchi)
        call copyr(aflxoi,1,xaoi1,1,mchi)
        call copyr(asorci,1,xaci1,1,mchi)
        call copyr(asordi,1,xadi1,1,mchi)
c
        call copyr(bpoli,1,xbpol1,1,mzones)
        call copyr(chi,1,xchi1,1,i1)
c
        call copyr(fflxii,1,xfii1,1,mchi)
        call copyr(fflxoi,1,xfoi1,1,mchi)
        call copyr(fsorci,1,xfci1,1,mchi)
        call copyr(fsordi,1,xfdi1,1,mchi)
c
c..15.09 Save Impurity Rate Equations variables
c
        if ((natomc.eq.3).and.(mimp.ne.0)) then
          do 115 ji=1,mimp
            nk = nkimp(ji)
c
            do 113 j=1,mzones
              xnsold(j,ji) = xns(j,ji)
  113       continue
c
            do 114 ik=0,nk
              do 114 j=1,mzones
                xnold(j,ik,ji) = xn(j,ik,ji)
  114       continue
c
            pflxo(ji) = pflx(ji)
            plosso(ji) = ploss(ji)
            psorco(ji) = psorc(ji)
            flouto(ji) = flout(ji)
            flscro(ji) = flscr(ji)
  115     continue
        end if
c
        do 120 jp=1,mchi
          xccon1(jp) = ccons(jp)
  120   continue
c
        call copyr(totali,1,xtot1,1,mchi)
c
        xebti1 = ebtoti
        xeomi1 = eohmi
        xepoi1 = epoyni
c
        xwbi1 = wbi
        xwomi1 = wohmi
        xwpoi1 = wpoyni
c
        xbcon1 = bcons
c
        xtbi1 = tbi
        xtai1 = tai
        xdti1 = dti
        xdtol1 = dtoldi
c
        liter = 10
        nliter = .false.
        nlrpet = .false.
c
        xerr = 0.0
        xdel = 0.0
c
c
        return
c
c  *** End of k = 1 segment ***
c
c
c      2)      k = 2
c
c  Start here for second (liter = 10) and third (liter = 11) calls
c  to RESOLV from STEPON with predictor-corrector. Both calls are
c  at bottom of STEPON following call to SOLVE.
c
  200   continue
c
c       define limits for dens,te,ti = zdmin,ztemin,ztimin
c       zdmin,ztemin,ztimin are in internal units
c       cfutz(129,130,131) are input in cm**-3, evs, evs
c
c       also set up zones to be excluded from timestep control.
c       default is to exclude outer dummy zone for entire run
c
c
        if(nstep.gt.2) go to 202
        zdmin=0.0
        ztemin=0.0
        ztimin=0.0
        if(cfutz(129).gt.epslon) zdmin =cfutz(129)*usid
        if(cfutz(130).gt.epslon) ztemin=cfutz(130)*evs*usih
        if(cfutz(131).gt.epslon) ztimin=cfutz(131)*evs*usih
        zxdel1=0.
        zxdel2=epsinv
        if(cfutz(ixdel1).gt.epslon) zxdel1=cfutz(ixdel1)*usit
        if(cfutz(ixdel2).gt.cfutz(ixdel1)) zxdel2=cfutz(ixdel2)*usit
        do 201 i=1,5
           jlimit(i)=ninfin
  201   continue
        zcentr=float(lcentr)
        if(nadump(1).gt.lcentr)     jlimit(1)=nadump(1)-1
        if(cfutz(izlimt).gt.zcentr) jlimit(4)=nadump(7)-1
        if(cfutz(izlimt).lt.epslon) jlimit(4)=ledge
        if(cfutz(iinner).gt.zcentr) jlimit(2)=int(cfutz(iinner)+.1)-1
        if(cfutz(iouter).gt.cfutz(iinner))
     .                              jlimit(3)=int(cfutz(iouter)+.1)+1
  202   continue
c
        nliter = .false.
        nlrpet = .false.
        zerr = 0.0
c
c
c           check for negative densities, energy, bpoloidal
c
c               if lcheck=.true., turn off timestep control
c               for  jlimit(2).le.jz.le.jlimit(3)  and
c               for  jz.gt.jlimit(4)
c               if chi<0, set chi(jz)=abs(chi(jz-1))
c
c
        lcheck=.false.
        if(tai.lt.zxdel1)                          go to 2020
        if(tai.gt.zxdel2)                          go to 2020
        lcheck=.true.
 2020   continue
c
c
c
c
c
        do 208 jz = lcentr, mzones
c
        do 204 jp = 1, mchi
        if(jz.gt.jlimit(1)) go to 203
c
c  15.00 Skip checks on impurities when using NC code
c
        if ((natomc.eq.3).and.(jp.ge.limp1)
     1                   .and.(jp.le.limpn)) go to 204
        if (chi(jp,jz).le.epslon) go to 250
  203   continue
c
c
  204   continue
c
c  Calculate electron density - prior to 15.00, done by summing chi
c
        if (versno.gt.14.99) then
          zdsumi = rhoels(2,jz) * usid
        else
          zdsumi = 0.0
          do 205 ip=lhyd1,limpn
            zdsumi = chi(ip,jz)*azzz(ip,jz) + zdsumi
  205     continue
        end if
c
c  Set minimum density
c  15.00 Skip for impurities when using NC code
c
        if (natomc.eq.3) then
          ipmax = lhydn
        else
          ipmax = limpn
        end if
c
        do 206 ip= lhyd1,ipmax
        if((jz.gt.jlimit(1)).and.(chi(ip,jz).lt.zdmin)) chi(ip,jz)=zdmin
  206   continue
c
c
c set minimum values for te,ti=(ztemin,ztimin)*zdsumi
c
        if(jz.le.jlimit(1)) go to 207
        zelmin=ztemin*zdsumi
        ziomin=ztimin*zdsumi
        if(chi(lelec,jz).lt.zelmin) chi(lelec,jz)=zelmin
        if(chi(lion,jz).lt.ziomin) chi(lion,jz)=ziomin
  207   continue
c
c
        if (bpoli(jz).le.epslon.and.jz.gt.lcentr) go to 251
c
c..Check total impurity densities. Exclude point at mzones since it may
c     be < 0 frequently (but small) when using 0 pedestal b. c.
c
        if ((natomc.eq.3).and.(mimp.ne.0)
     1      .and.(jz.ne.mzones).and.(limprd(1).eq.0)) then
          do 210 ji=1,mimp
            if (xns(jz,ji).le.epslon) go to 252
  210     continue
        end if
c
  208   continue
c
c               do extrapolation, predictor-corrector
c
c  If liter =10 then ipred = 0, and if liter = 11 then ipred = 1.
c  So, to finish the first call to RESOLV(2) with predictor-corrector
c  (second call to RESOLV), jump to 2000
c
        ipred = liter - (liter/10)*10
        if (thetap.gt.epslon.and.ipred.le.0) go to 2000
        if (nlextr) go to 1000
c
c
c              compute xdel
c
c  Finish predictor-corrector; i.e., finish second call to RESOLV(2),
c  third call to RESOLV overall.
c
  220   continue
c
        xdel = 0.0
        lpdel = 0
        lzdel = 0
c
        do 228 jz = 2, mzones
c
        do 224 jp = 1, mchi
c
c  15.00 Skip checks on impurities when using NC code
c
        if ((natomc.eq.3).and.(jp.ge.limp1)
     1                   .and.(jp.le.limpn)) go to 224
        zdelc = abs(chi(jp,jz) - xchi1(jp,jz))
        if (zdelc.le.epslon) go to 222
        if((jz.le.jlimit(1)).and.(chi(jp,jz).le.epslon)) go to 250
        zdelc = zdelc / chi(jp,jz)
        if((zdelc.le.xdel).or.(jz.lt.lcentr)) go to 222
        if(jz.gt.jlimit(4))                   go to 222
        if(.not.lcheck)                       go to 2210
        if(jz.le.jlimit(2))                   go to 2210
        if(jz.gt.jlimit(3))                   go to 2210
        go to 222
 2210   continue
        xdel = zdelc
        lpdel = jp
        lzdel = jz
  222   continue
c
  224   continue
c
        zdelb = abs(bpoli(jz) - xbpol1(jz))
        if (zdelb.le.epslon) go to 226
        if (bpoli(jz).le.epslon) go to 251
        zdelb = zdelb / bpoli(jz)
        if (zdelb.le.xdel) go to 226
        xdel = zdelb
        lpdel = mxchi + 1
        lzdel = jz
  226   continue
c
  228   continue
c
c..15.09 Incorporate xnerr into dti calculation; responsible
c        species and zone unknown. Skip for first time-step
c        since extensive ionization from +1 occurs.
c
        if ((natomc.eq.3).and.(xnerr.gt.xdel/delmax)
     1       .and.(tai*uist.gt.tinit).and.(limprd(1).eq.0)) then
          xdel = xnerr * delmax
          lpdel = limp1
          lzdel = 2
c
        end if
c
        zerr = max(zerr , xdel/delmax)
c
c
c              compute new dti
c
c
  230   continue
        if (zerr.gt.1.0) go to 235
c
c
c      2.2.1)  zerr is less than 1.0 -- increase dt
c               here "zerr" is max( xerr/errmax , xdel/delmax , etc.)
c
c
c               new dti is
c               ((1-xfutz(idamp)) * (1/zerr) + xfutz(idamp))
c               * old dti,
c               except new dti .le. xfutz(ilim) * old dti
c
c               note-- 0 .lt. xfutz(idamp) .lt. 1 .lt. xfutz(ilim)
c
        zerr1 = max(zerr ,
     1          (1.0 - xfutz(idamp)) / (xfutz(ilim) - xfutz(idamp)) )
        dti = xdti1 * (xfutz(idamp) + (1.0 - xfutz(idamp)) / zerr1)
        go to 240
c
c
c      2.2.2)  zerr greater than or equal to 1.0
c
c
  235   continue
        dti = xdti1 / zerr
c
c
c      2.3)    test for timestep repetition, restrict range of dti
c
c
  240   continue
c
        if (xfutz(ierrl)*errmax.lt.xerr) go to 260
        if (xfutz(idell)*delmax.lt.xdel) go to 270
c
        if (dti.gt.dtmaxi) dti = dtmaxi
        if (dti.le.dtmini) go to 9024
c
c
        return
c
c  *** End of predictor-corrector loop ***
c
c      2.5)    decided to repeat timestep
c               print message
c
c
c               negative chi
c
  250   continue
        call blines(2)
        call mesage(' *** non-positive density or energy')
        call   ivar('timestep',nstep )
cb        call harray(' species',lhspec(jp),1)
        write (6,*) ' species  = ',lhspec(jp)
        call ivar(' zone   ',jz-1)
        call rvar(' value  ',chi(jp,jz))
        dti = xfutz(ineg) * xdti1
        call mesage(' *** timestep repeated')
        go to 300
c
c               negative bpoli
c
  251   continue
        call blines(2)
        call mesage(' *** non-positive b-poloidal')
        call   ivar('timestep',nstep )
        call ivar(' zone   ',jz-1)
        call rvar(' value  ',bpoli(jz))
        dti = xfutz(ineg) * xdti1
        call mesage(' *** timestep repeated')
        go to 300
c
c               negative impurity density
c
  252   continue
        call blines(2)
        call mesage(' *** non-positive density')
        call   ivar('timestep',nstep )
        jp = ji + limp1 - 1
cb        call harray(' species',lhspec(jp),1)
        write (6,*) ' species  = ',lhspec(jp)
        call ivar(' zone   ',jz-1)
        call rvar(' value  ',xns(ji,jz))
        dti = xfutz(ineg) * xdti1
        call mesage(' *** timestep repeated')
        go to 300
c
c               xerr too high
c
  260   continue
        call blines(2)
        call mesage(' *** extrapolated error too high')
        call   ivar('timestep',nstep )
        if (lperr.le.mxchi) write (6,*) ' species  = ',lhspec(lperr)
cb        call harray(' species',lhspec(lperr),1)
        if (lperr.gt.mxchi)
     1  call mesage(' species = b-poloidal')
        call ivar(' zone   ',lzerr-1)
        if (lperr.le.mxchi) call rvar(' value  ',chi(lperr,lzerr))
        if (lperr.gt.mxchi) call rvar(' value  ',bpoli(lzerr))
        call rvar(' error  ',xerr)
        call rvar(' limit  ',xfutz(ierrl)*errmax)
        call mesage(' *** timestep repeated')
        if (xfutz(ihigh).gt.0.0) dti = xdti1 * xfutz(ihigh)
        go to 300
c
c               xdel too high
c
  270   continue
        call blines(2)
        call mesage(' *** change over one timestep too high')
        call   ivar('timestep',nstep )
        if (lpdel.le.mxchi) write (6,*) ' species  = ',lhspec(lpdel)
cb        call harray(' species',lhspec(lpdel),1)
        if (lpdel.gt.mxchi)
     1  call mesage(' species = b-poloidal')
        call ivar(' zone   ',lzdel-1)
        if (lpdel.le.mxchi) call rvar(' value  ',chi(lpdel,lzdel))
        if (lpdel.gt.mxchi) call rvar(' value  ',bpoli(lzdel))
        call rvar(' change ',xdel)
        call rvar(' limit  ',xfutz(idell)*delmax)
        call mesage(' *** timestep repeated')
        go to 300
c
c
c      3)      repeat the timestep
c
c
  300   continue
c
        i1 = mxchi * mzones
c
c
c
        call copyr(xaii1,1,aflxii,1,mchi)
        call copyr(xaoi1,1,aflxoi,1,mchi)
        call copyr(xaci1,1,asorci,1,mchi)
        call copyr(xadi1,1,asordi,1,mchi)
c
        call copyr(xbpol1,1,bpoli,1,mzones)
        call copyr(xchi1,1,chi,1,i1)
c
        call copyr(xfii1,1,fflxii,1,mchi)
        call copyr(xfoi1,1,fflxoi,1,mchi)
        call copyr(xfci1,1,fsorci,1,mchi)
        call copyr(xfdi1,1,fsordi,1,mchi)
c
c  15.09 Impurity Rate Equations variables
c
        if ((natomc.eq.3).and.(mimp.ne.0)) then
          do 330 ji=1,mimp
            pflx(ji) = pflxo(ji)
            ploss(ji) = plosso(ji)
            psorc(ji) = psorco(ji)
            flout(ji) = flouto(ji)
            flscr(ji) = flscro(ji)
            nk = nkimp(ji)
c
            do 325 j=1,mzones
              xns(j,ji) = xnsold(j,ji)
  325       continue
c
            do 327 ik=0,nk
              do 327 j=1,mzones
                xn(j,ik,ji) = xnold(j,ik,ji)
  327       continue
c
  330     continue
        end if
c
        do 320 jp=1,mchi
          ccons(jp) = xccon1(jp)
  320   continue
c
        call copyr(xtot1,1,totali,1,mchi)
c
        ebtoti = xebti1
        eohmi = xeomi1
        epoyni = xepoi1
c
        wbi = xwbi1
        wohmi = xwomi1
        wpoyni = xwpoi1
c
        bcons = xbcon1
c
        tai = xtai1
        tbi = tai + dti
c
        xtbi1 = tbi
        xdti1 = dti
c
        liter = 10
        nliter = .true.
        nlrpet = .false.
c
        xerr = 0.0
        xdel = 0.0
c
        call getchi(1)
c
c
        return
c
c
c
c
c      10)     extrapolation
c
c
c
 1000   continue
        if (liter.ge.30) go to 1300
        if (liter.ge.20) go to 1200
c
c
c      11)     end of 1st minor timestep: save end-of-timestep values,
c               go back and do 1st half-timestep
c               (from tai to (tai+tbi)/2, 2nd minor timestep)
c
c
 1100   continue
        i1 = mxchi * mzones
c
c               save end-of-timestep values
c
        call copyr(aflxii,1,xaii2,1,mchi)
        call copyr(aflxoi,1,xaoi2,1,mchi)
        call copyr(asorci,1,xaci2,1,mchi)
        call copyr(asordi,1,xadi2,1,mchi)
c
        call copyr(bpoli,1,xbpol2,1,mzones)
        call copyr(chi,1,xchi2,1,i1)
c
        call copyr(fflxii,1,xfii2,1,mchi)
        call copyr(fflxoi,1,xfoi2,1,mchi)
        call copyr(fsorci,1,xfci2,1,mchi)
        call copyr(fsordi,1,xfdi2,1,mchi)
c
        call copyr(ccons,1,xccon2,1,mchi)
        call copyr(totali,1,xtot2,1,mchi)
c
        xebti2 = ebtoti
        xeomi2 = eohmi
        xepoi2 = epoyni
c
        xwbi2 = wbi
        xwomi2 = wohmi
        xwpoi2 = wpoyni
c
        xbcon2 = bcons
c
c               restore old values
c
        call copyr(xaii1,1,aflxii,1,mchi)
        call copyr(xaoi1,1,aflxoi,1,mchi)
        call copyr(xaci1,1,asorci,1,mchi)
        call copyr(xadi1,1,asordi,1,mchi)
c
        call copyr(xbpol1,1,bpoli,1,mzones)
        call copyr(xchi1,1,chi,1,i1)
c
        call copyr(xfii1,1,fflxii,1,mchi)
        call copyr(xfoi1,1,fflxoi,1,mchi)
        call copyr(xfci1,1,fsorci,1,mchi)
        call copyr(xfdi1,1,fsordi,1,mchi)
c
        call copyr(xccon1,1,ccons,1,mchi)
        call copyr(xtot1,1,totali,1,mchi)
c
        ebtoti = xebti1
        eohmi = xeomi1
        epoyni = xepoi1
c
        wbi = xwbi1
        wohmi = xwomi1
        wpoyni = xwpoi1
c
        bcons = xbcon1
c
        liter = 20
        nliter = .false.
        nlrpet = .true.
c
        if (thetap.gt.epslon) nliter = .true.
        if (thetap.gt.epslon) nlrpet = .false.
c
        tai = xtai1
        dti = 0.5 * xdti1
        tbi = xtai1 + dti
c
c
        return
c
c
c      12)     end of 2nd minor timestep: set values to do 2nd
c               half-timestep (from (tai+tbi)/2 to tbi, 3rd minor t.s.)
c
c
 1200   continue
        liter = 30
        nliter = .true.
        nlrpet = .false.
        dti = 0.5 * xdti1
        tai = xtai1 + dti
        tbi = xtbi1
        dtoldi = dti
c
        call getchi(1)
c
c
        if (thetap.gt.epslon) go to 2200
        return
c
c
c      13)     after third time through: compute  f(0), error
c
c
c       13.1)   get extrapolated chi and bpoli
c               get xerr, lperr, lzerr
c
c       error for each param. is (zchi1 - chi) / new chi,
c       and similarly for bpol
c       total error is maximum of individual errors
c
 1300   continue
c
        tai = xtai1
        tbi = xtbi1
        dti = xdti1
        dtoldi = xdtol1
c
        liter = 0
        lperr = 0
        lzerr = 0
        xerr = 0.0
c
        do 1314 jz = 1, mzones
c
        do 1308 jp = 1, mchi
        zdc = chi(jp,jz) - xchi2(jp,jz)
        zc = zdc + chi(jp,jz)
        zdc = abs(zdc)
c
        if (zdc.le.epslon) go to 1304
        if (zc.le.epslon) go to 1304
        ze = 2.0 * zdc / zc
        if (ze.le.xerr.or.jz.lt.lcentr.or.jz.gt.ledge) go to 1304
        xerr = ze
        lperr = jp
        lzerr = jz
 1304   continue
c
        chi(jp,jz) = zc
c
 1308   continue
c
        zdb = bpoli(jz) - xbpol2(jz)
        bpoli(jz) = zdb + bpoli(jz)
        zb = bpoli(jz)
        zdb = abs(zdb)
c
        if (zdb.le.epslon) go to 1310
        if (zb.le.epslon) go to 1310
        ze = 2.0 * zdb / zb
        if (ze.le.xerr.or.jz.lt.lcentr.or.jz.gt.ledge) go to 1310
        xerr = ze
        lperr = mxchi + 1
        lzerr = jz
 1310   continue
c
 1314   continue
c
c
c              set conservation variables
c
c
        do 1318 jp = 1, mchi
c
        aflxii(jp) = 2.0 * aflxii(jp) - xaii2(jp)
        aflxoi(jp) = 2.0 * aflxoi(jp) - xaoi2(jp)
        asorci(jp) = 2.0 * asorci(jp) - xaci2(jp)
        asordi(jp) = 2.0 * asordi(jp) - xadi2(jp)
c
        fflxii(jp) = 2.0 * fflxii(jp) - xfii2(jp)
        fflxoi(jp) = 2.0 * fflxoi(jp) - xfoi2(jp)
        fsorci(jp) = 2.0 * fsorci(jp) - xfci2(jp)
        fsordi(jp) = 2.0 * fsordi(jp) - xfdi2(jp)
c
        ccons(jp) = 2.0 * ccons(jp) - xccon2(jp)
        totali(jp) = 2.0 * totali(jp) - xtot2(jp)
 1318   continue
c
        ebtoti = 2.0 * ebtoti - xebti2
        eohmi  = 2.0 * eohmi  - xeomi2
        epoyni = 2.0 * epoyni - xepoi2
c
        wbi    = 2.0 * wbi    - xwbi2
        wohmi  = 2.0 * wohmi  - xwomi2
        wpoyni = 2.0 * wpoyni - xwpoi2
c
        bcons  = 2.0 * bcons  - xbcon2
c
        zerr = max(zerr , xerr/errmax)
c
c
        go to 220
c
 221    return
c
c
c      20)     predictor-corrector
c
c  Finish first call to RESOLV(2) with predictor-corrector
c  (third call to RESOLV) overall
c
 2000   continue
        if (liter.ge.30) go to 2100
c
c
c              get coefficients at t.s. n+thetap
c               old t.s. values are in x---1
c
c
        zth = 1.0 - thetap
c
        do 2008 jz = 1, mzones
c
        do 2004 jp = 1, mchi
        chi(jp,jz) = xchi1(jp,jz)*zth + chi(jp,jz)*thetap
 2004   continue
c
        bpoli(jz) = xbpol1(jz)*zth + bpoli(jz)*thetap
 2008   continue
c
c  15.09 Do also for Impurity Rate Equations variables
c
        if ((natomc.eq.3).and.(mimp.ne.0)) then
          do 2010 ji=1,mimp
            nk = nkimp(ji)
            do 2010 ik=0,nk
              do 2010 j=1,mzones
                xn(j,ik,ji) = xnold(j,ik,ji)*zth + xn(j,ik,ji)*thetap
 2010     continue
        end if
c
c  Evaluate fundamental physical variables from chi
c
c
        call getchi(1)
c
c
c              repeat timestep (iterate)
c
c
        i1 = mxchi * mzones
c
c
c  Reset chi, etc. to values at end of last time-step so that the
c  chi's appearing explicitly in the solution procedure are NOT the
c  thetap-centered ones (used only to time-center the coefficients),
c  and so that integrated quantities (e.g., aflxii) are not altered
c  by the predictor step.
c
        call copyr(xaii1,1,aflxii,1,mchi)
        call copyr(xaoi1,1,aflxoi,1,mchi)
        call copyr(xaci1,1,asorci,1,mchi)
        call copyr(xadi1,1,asordi,1,mchi)
c
        call copyr(xbpol1,1,bpoli,1,mzones)
        call copyr(xchi1,1,chi,1,i1)
c
        call copyr(xfii1,1,fflxii,1,mchi)
        call copyr(xfoi1,1,fflxoi,1,mchi)
        call copyr(xfci1,1,fsorci,1,mchi)
        call copyr(xfdi1,1,fsordi,1,mchi)
c
c  15.12 Add theta-p centering of flout and flscr
c  15.09 Impurity Rate Equations variables
c
        if ((natomc.eq.3).and.(mimp.ne.0)) then
          do 2040 ji=1,mimp
            pflx(ji) = pflxo(ji)
            ploss(ji) = plosso(ji)
            psorc(ji) = psorco(ji)
            flout(ji) = flouto(ji)*zth + flout(ji)*thetap
            flscr(ji) = flscro(ji)*zth + flscr(ji)*thetap
            nk = nkimp(ji)
c
            do 2035 j=1,mzones
              xns(j,ji) = xnsold(j,ji)
 2035       continue
c
            do 2037 ik=0,nk
              do 2037 j=1,mzones
                xn(j,ik,ji) = xnold(j,ik,ji)
 2037       continue
c
 2040     continue
        end if
c
        do 2025 jp=1,mchi
          ccons(jp) = xccon1(jp)
 2025   continue
c
        call copyr(xtot1,1,totali,1,mchi)
c
        ebtoti = xebti1
        eohmi = xeomi1
        epoyni = xepoi1
c
        wbi = xwbi1
        wohmi = xwomi1
        wpoyni = xwpoi1
c
        bcons = xbcon1
c
        liter = liter + 1
        nliter = .true.
        nlrpet = .false.
c
c
        return
c
c  *** End of first call to RESOLV(2). Return to STEPON, which ***
c  *** next calls COEF.                                        ***
c
c      21)     same as 20), but old vals. are in x--3 vars.
c
c
 2100   continue
c
c
c              get coefficients at t.s. n+thetap
c               old t.s. values are in x---3
c
c
        zth = 1.0 - thetap
c
        do 2108 jz = 1, mzones
c
        do 2104 jp = 1, mchi
        chi(jp,jz) = xchi3(jp,jz)*zth + chi(jp,jz)*thetap
 2104   continue
c
        bpoli(jz) = xbpol3(jz)*zth + bpoli(jz)*thetap
 2108   continue
c
c
c              get vals. for coeffs.
c
c
        call getchi(1)
c
c
c              repeat timestep (iterate)
c
c
        i1 = mxchi * mzones
c
c
c
        call copyr(xaii3,1,aflxii,1,mchi)
        call copyr(xaoi3,1,aflxoi,1,mchi)
        call copyr(xaci3,1,asorci,1,mchi)
        call copyr(xadi3,1,asordi,1,mchi)
c
        call copyr(xbpol3,1,bpoli,1,mzones)
        call copyr(xchi3,1,chi,1,i1)
c
        call copyr(xfii3,1,fflxii,1,mchi)
        call copyr(xfoi3,1,fflxoi,1,mchi)
        call copyr(xfci3,1,fsorci,1,mchi)
        call copyr(xfdi3,1,fsordi,1,mchi)
c
        call copyr(xccon3,1,ccons,1,mchi)
        call copyr(xtot3,1,totali,1,mchi)
c
        ebtoti = xebti3
        eohmi = xeomi3
        epoyni = xepoi3
c
        wbi = xwbi3
        wohmi = xwomi3
        wpoyni = xwpoi3
c
        bcons = xbcon3
c
        liter = liter + 1
        nliter = .true.
        nlrpet = .false.
c
c
        return
c
c
c      22)     save end-of-2nd-minor-timestep values
c               (end of first half-timestep)
c
c
 2200   continue
c
c
        i1 = mxchi * mzones
c
        call copyr(aflxii,1,xaii3,1,mchi)
        call copyr(aflxoi,1,xaoi3,1,mchi)
        call copyr(asorci,1,xaci3,1,mchi)
        call copyr(asordi,1,xadi3,1,mchi)
c
        call copyr(bpoli,1,xbpol3,1,mzones)
        call copyr(chi,1,xchi3,1,i1)
c
        call copyr(fflxii,1,xfii3,1,mchi)
        call copyr(fflxoi,1,xfoi3,1,mchi)
        call copyr(fsorci,1,xfci3,1,mchi)
        call copyr(fsordi,1,xfdi3,1,mchi)
c
        call copyr(ccons,1,xccon3,1,mchi)
        call copyr(totali,1,xtot3,1,mchi)
c
        xebti3 = ebtoti
        xeomi3 = eohmi
        xepoi3 = epoyni
c
        xwbi3 = wbi
        xwomi3 = wohmi
        xwpoi3 = wpoyni
c
        xbcon3 = bcons
c
c
        return
c
c..............................................................................
c
c
c
c      90)     errors
c
c
 9000   continue
        call error_olymp(0,iclass,isub,0,
     1          'resolv called with k .gt.2 or .le.0')
        call error_olymp(3,k,2,1,'  k     ')
        return
c
c      timestep below minimum
c
 9024   continue
        call error_olymp(1,iclass,isub,2,
     1          '? timestep below minimum')
        call error_olymp(2,dti,1,1,'   dt    ')
        return
        end

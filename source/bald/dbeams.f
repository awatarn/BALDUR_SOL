c 14:00 22-Apr-96 .../baldur/code/bald/dbeams.f
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  to obtain this file, type
c cfs get /11040/bald91/wbaldn1
c end
c lib wbaldn1 ^ x dbeams ^ end
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..this file contains the following subroutines:
c
c  beams  - compute neutral beam injection sources and sinks
c              shbems, wibems, webems, shchxs, shblos (comdif)
c              gblosi (comneu)
c  fusigv - fusion cross sections for beam ions with thermal ions
c  hrecrd - read and write restart records
c  hclear - clear variables and arrays
c  dposit - deposit high-energy ions, interface with freya
c  shafs  - calculates analytic shafranov shift for use in beam package
c           NOTE: included only to allow agreement with 1-D code when
c                 leqtyp=0.
c
c  **  monte carlo freya package for neutrals  **
c
c..revised by mikkelsen, april 1987
c
c  gethfr -
c  nhrset -
c  setrmj - height above midplane of beam rays and major radius
c              to the intersections
c  setwgt - beam intensity across its cross section
c  sntgrl - gaussian intgrl across face of ion source to find intensity
c  seghfr - tracks rays through plasma
c               to find the path length in each intersected zone
c  aioniz - ionization cross sections
c  fsig   - cross sections for ionization and charge exchange
c  ttyout - sends character string to output unit 6
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@beams   .../baldur/code/bald/dbeams.f
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 14-aug-96 replaced fuzz and round with rndeps
c  rgb 07-aug-94 move call scaler(gblosi,mxhyd,uist*uisl**2*rmins)
c    after 850 continue end of loop over species
c  rgb 27-jun-94 extrapolate to zone center just outside plasma
c    for rhobis, rhobes, hebems, and ajbs
c       dps 18-dec-87 note bug in tis*uieh; kluge to match 1-D
c       dps 16-dec-87 spliced in following 3 changes from 1-D code;
c                     altered volume elements and radii for eq'm.
c       dps 19-jun-87 fix bugs in d-d beam-beam calculations
c       dps 05-may-87 add beam - beam reactions as diagnostic
c       dps 29-apr-87 allow multiple beam species
c       drm 4-jan-85 correct coding for k .eq. 4 case
c       drm 18-dec-84 add coding for k .eq. 4 case
c       mhr 4-dec-84 use fokbem from mikkelsen
c       drm 10-jan-84 added bpabs and bploss calculation
c       drm 22-feb-83 added zvbeam factor to zddsig(je); added
c       cos(hangle3) to zrtang in beam-driven current
c       mhr 06-aug-82 changed all fcau->fcmp,
c       fxau->fxnucl; add call scaler(halfas, ...)
c       after $200; add idtfus code; hei(je)->(hei(jz)
c       -0.5*hdei(je)) to evaluate hslows, hscats,
c       heion, heelec, zsigcx, hchexs, hsigv, zddsig;
c       moved d-t fusion into do 799 loop and
c       add factor of zdt20; add hfutz(1) and
c       zdt20 to shblos after do 684; use
c       duane's d + d -> he3 + n cross section
c       and remove d+d fusion as a real particle
c       loss; use zdti2 in hddtot; added zloge
c       above do 459; changed approximation for
c       trapped electron correction to beam-driven current;
c       added qm form of zlogi
c       aes 11-feb-82 add hfutz(5);  multiplies slowdown rate hslows
c       aes 26-nov-80 fixed bug near line 622, elliptic integral biz
c       aes 30-sep-80 added elliptic integral business with hrmaj,
c               hrmin as vectors (index = ibeam1)
c       aes ??-aug-80 fixed misplaced parenths. in z0 near line 900
c       fgps 12-aug-80 adjusted hesrci(nhsrc) so as to never exceed
c                      the highest bin energy hei(nhe)
c       drm 30-apr-80 changed zbar to zbarlg, dropped hlogbs,
c               added zlogi and zefflg, corrected hslows,
c               moved 414 and 418 loops inside 459 loop
c       amck 29-mar-78 include uisd in ajbs computation
c       amck 29-mar-78 change zjb -> zbj before 949
c**********************************************************************c
c
c
        subroutine beams(k)
c
c
cl      2.4     compute neutral beam sources and sinks
c
c
c
c
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'cfokkr.m'
        include 'commhd.m'
c
c
        logical ilinit, ilini2
c
c
        dimension
     r  zddsig(20),
     r  ze(10),         zei(20),        zf(10),         zfi(20),
     r  zsigcx(20),     zsrci(10),      zu(20),         zusqrt(20)
c
c
        data    ilinit /.false./,       ilini2 /.false./
c
c
c------------------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /4/
c
c
        if ( nlomt2(isub) ) then
          write (6,*) ' *** 2.4 subroutine beams bypassed'
          return
        endif
c
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       comfok, shbems, wibems, webems, shchxs, shblos (comdif)
c       gblosi (comneu)
c
c------------------------------------------------------------------------------
c
c
cl      local variables:
c
c       ibeam1  - index of first injector with present species.
c                       this is set inside a loop over species.
c                       habeam(ibeam1) is the atomic weight, etc.
c       ie      - backwards e index
c       imu     - backwards mu index
c       je      - forwards e index
c       jmu     - forwards mu index
c       jsrc    - index in hesrci and hfsrci
c       jz      - zone index
c
c       z0
c       z1,
c       z2, z3, - coefficients in equation
c               z1*f(jmu-1) + z2*f(jmu) * z3*f(jmu+1)
c               = f(new timestep) - f(old timestep)
c               where f(jmu..) is thetai* new timestep hfi(jmu..) +
c               (1 - thetai)* old timestep hfi(jmu..)
c       za, zb,
c       zc, zd  - coefficients in equation
c               za*hfi(jmu-1) + zb*hfi(jmu) + zc*hfi(jmu+1) = zd
c               for energy ie, zone jz, and mu jmu, and new timestep hfi.
c               old timesep values of hfi are included in zd.
c       zdti1   - time during this timestep for which
c               beam slowing down and losses have already been
c               computed, and therefor the old values of shbems, etc.
c               should be used.  0 unless a timestep has recently
c               been repeated.
c       zdti2   - time during this timestep for which
c               beam slowing down and losses have not been
c               computed, and hence new plasma source and sink terms
c               must be computed.  this is the timestep
c               for the fast-ion slowing-down calculation this timestep.
c       zdt10   - zdti1 / dti
c       zdt20   - zdti2 / dti
c       ze(jmu),
c       zf(jmu) - coefficients in equation hfi(jmu) + ze(jmu)*hfi(jmu+1)
c                       = zf(jmu), computed from za,... equations.
c       zhfold  - old timestep value of hfi(imu,ie,jz)
c       zhfnew  - new timestep value of hfi(imu,ie,jz)
c       zhf     - value of hfi time-centered according to thetai
c       zeplas, znplas, zeloss, znloss,
c               - used in computation of webems, wibems, shbems, shchxs,
c                       heplas, etc.
c       zint    - interpolation factor used to split up a beam source into
c               beam groups so that the proper energy and number of particles
c               is dumped.   0 < zint <= 1.0
c       zeb     - beam energy in kev, used in d+d reaction comp.
c       zloge   - alog10(zeb)
c       zleb    - log of beam energy je in ev *
c                       (mass hydr.1/mass beam ion)
c       zsigcx(je) - sigma v for charge exchange of beam ions with
c                       neutral gas at energy je
c       zsrci(jmu) - source rate of particles entering group ie.
c               includes slowed-down (from previous ie-step) and
c               injected (from this ie-step).  in dens/time/unit mu
c               in internal units
c       zth     -  = (1.0 - thetai)
c       zvbeam  - v of beam ions at energy je
c       zvoli   - volume in internal units of this (jz) zone
c       zvols   - same as zvoli, but in standard units
c
c
c------------------------------------------------------------------------------
c
c
cl      1)      initialization from input variables
c
c
  100   continue

        if (k.eq.2) go to 200
        if (k.eq.3) go to 1000
        if (k.eq.4) go to 200
        if (k.ne.1) return
c
cl              set mx-- variables
c
        mxhe = 20
        mxhmu = 10
        mxhsrc = 10
        mxhsp = 2
c
        call dposit(-1)
c
c
        if (htbi.ge.epsinv) return
        if (nlres) return
c
c------------------------------------------------------------------------------
c
c
cl              initialization of mesh
c
c
c
cl              clear source terms for transport equations
c
        i1 = mxzone*mxhyd
        i2 = mxzone*mxhsp
        call resetr(shbems,i1,0.0)
        call resetr(shchxs,i1,0.0)
        call resetr(shblos,i1,0.0)
        call resetr(webems,mxzone,0.0)
        call resetr(wibems,mxzone,0.0)
        call resetr(websps,i2,0.0)
        call resetr(wibsps,i2,0.0)
c
c
c
cl      2.1)    set species and e variables
c
        mhsp = 1
        nhbem1(1) = 0
c
c       loop over beams
c
        do 103 jb = 1, mxhbem
        if (nhbeam(jb).eq.0) go to 103
        if (mhsp.ne.1) then
        do 102 jsp = 1, mhsp-1
c
c       identify first occurrence of a particular species
c
        if (nhbeam(jb).eq.nhbeam(nhbem1(jsp))) go to 103
  102   continue
        end if
c
c       allow only mxhsp different species
c
        if (mhsp.gt.mxhsp) go to 103
        nhbem1(mhsp) = jb
        mhsp = mhsp + 1
  103   continue
        mhsp = max( mhsp - 1, 1)
c
c       set up separate energy grids for each species
c
        do 113 jsp = 1, mhsp
        ibeam1 = nhbem1(jsp)
        zton = epsinv
        z0 = 0.0
c
        do 104 jb = 1, mxhbem
        if (nhbeam(jb).eq.nhbeam(ibeam1)) then
c
c       identify max. energy and min. injection time for each species
c
        z0 = max(z0,hebeam(jb))
        zton = min(zton,hton(jb))
        end if
  104   continue
c
        htspon(jsp) = zton
        z0 = z0*uesh*usie / float(nhe)
c
        do 108 je = 1, nhe
        hei(je,jsp) = z0 * float(je)
        hdei(je,jsp) = z0
  108   continue
c
c       set hesrci(nhsrc) .le. hei(nhe)
c
        if(hesrci(nhsrc,jsp).gt.(hei(nhe,jsp)*1.0001)) go to 9060
        hesrci(nhsrc,jsp)=min(hesrci(nhsrc,jsp),hei(nhe,jsp))
c
cl      set lhemin, different for each species
c
        do 112 jz = lcentr, ledge
        ztii = tis(2,jz) * usih
c
        do 110 je = 1, nhe
        if (ztii.lt.hei(je,jsp)) go to 111
  110   continue
c
  111   continue
        lhemin(jz,jsp) = je
  112   continue
c
  113   continue
        lhemax = nhe
c
c
cl              set mu variables
c
c
        z0 = fcpi / float(nhmu)
        hdmub(1) = 0.0
        hmub(1) = -1.0
        hsin(1) = 0.0
        hdmub(nhmu+1) = 0.0
        hsin(nhmu+1) = 0.0
c
        do 124 jmu = 1, nhmu
c
c
cl              if the mu mesh is to be changed, change this statement
c
c
        hmub(jmu+1) = - cos(float(jmu)*z0)
c
c
c
        hmuc(jmu) = (hmub(jmu+1) + hmub(jmu))*0.5
        hdmuc(jmu) = hmub(jmu+1) - hmub(jmu)
        hdmu3(jmu) = 0.333333333 * (hmub(jmu+1)**3 - hmub(jmu)**3)
c
        if(jmu.le.1) go to 124
        hdmub(jmu) = hmuc(jmu) - hmuc(jmu-1)
        hsin(jmu) = (1.0 - hmub(jmu)**2) / hdmub(jmu)
  124   continue
c
        chscat = fcpi * fces**4 / sqrt(2.0*fcau) *
     1          10.0**(fxes*4.0 - 0.5*fxau)
        czgamm = chscat * 2.0
        cznue = 4.0 * sqrt(2.0 * fcpi * fcme) * fces**4 / (3.0 * fcau) *
     1          10.0**(0.5*fxme + 4.0*fxes - fxau)
c
c       calculate cross sections for beam - beam fusion reactions
c
c       set variables for gyro-angle integration
c       zu = cos ( gyro-angle )
c
        imxju = 20
        zdu = 2.0 / imxju
        do 130 ju = 1, imxju
        zuju = -1.0 + (ju-0.5) * zdu
        zu(ju) = zuju
        zusqrt(ju) = 1.0 / sqrt(1.0-zuju**2)
  130   continue
c
c       identify indices for d-t beam-beam reactions, nbbsp(1 & 2)
c       and for d-d, nbbspd.
c
        nbbsp1 = 0
        nbbsp2 = 0
        nbbspd = 0
        do 135 jsp = 1, mhsp
        ibeam1 = nhbem1(jsp)
        if (nbbspd.eq.0.and.nhbeam(ibeam1).eq.-2) nbbspd=jsp
        if (nbbsp1.eq.0) then
        if (nhbeam(ibeam1).eq.-2) then
        nbbsp1 = jsp
        iother = -3
        else if (nhbeam(ibeam1).eq.-3) then
        nbbsp1 = jsp
        iother = -2
        end if
        else if (nhbeam(ibeam1).eq.iother) then
        nbbsp2 = jsp
        end if
  135   continue
c
c       skip all cross section calculations if there is no d,
c       skip just d-t calculation if no t.
c
        if (nbbspd.eq.0) go to 185
        if (nbbsp2.eq.0) go to 150
c
c       begin d-t cross section calculation
c
        call resetr(bbdtsv,mxhmu**2*mxhe**2,0.0)
        bbdt = 0.0
c
        za1 = habeam(nhbem1(nbbsp1))
        za2 = habeam(nhbem1(nbbsp2))
c
c       zamu = mass of 1 amu in grams, z1 & z2 are then 2. / m(1 & 2) in gm**-1
c
        zamu = fcau * 10.0**fxau
        z1 = 2.0 / (zamu*za1)
        z2 = 2.0 / (zamu*za2)
c
c       begin loop over velocity variables
c       zv(1 & 2)sq = velocity **2 in (cm/sec)**2
c
        do 149 je1 = 1, nhe
        zv1sq = (hei(je1,nbbsp1)-0.5*hdei(je1,nbbsp1))*uise*z1
        zv1 = sqrt(zv1sq)
        do 149 je2 = 1, nhe
        zv2sq = (hei(je2,nbbsp2)-0.5*hdei(je2,nbbsp2))*uise*z2
        zv2 = sqrt(zv2sq)
        do 140 jmu1 = 1, nhmu
        zmu1 = hmuc(jmu1)
        zsqmu1 = sqrt(1.0-zmu1**2)
c
c       separate inner mu loop to take advantage of symmetry
c
        do 140 jmu2 = 1, jmu1
        zmu2 = hmuc(jmu2)
        zsqmu2 = sqrt(1.0-zmu2**2)
c
c       average over cos( gyro-angle )
c
        zsigv = 0.0
        do 139 ju= 1, imxju
        zvrel2 = zv1sq + zv2sq - 2.0*zv1*zv2 * ( zmu1*zmu2
     1           + zsqmu1*zsqmu2*zu(ju) )
        zvrel = sqrt(zvrel2)
c
c       zev is equivalent energy of d in ev
c
        zev = 0.5 * zamu * habeam(nhbem1(nbbspd)) * zvrel2 *evsinv
c
c               sigma-v for fusion of energetic deuterium on
c               tritium
c               cross section from b. h. duane in batelle northwest
c               report bnwl-1685
c               sigma-v is considered negligible for e < 7 kev
c
        zsigma = 0.0
        if (zev.ge.7000.0) zsigma =
     1  (5.02e+7 / (1.0 + (1.368e-5*zev - 1.076)**2) + 4.09e+5)
     2          / (zev * (exp(1453.0 / sqrt(zev)) - 1.0))
        zsigv = zsigv + zsigma*zvrel*1.0e-24*zusqrt(ju)
  139   continue
        bbdtsv(jmu2,jmu1,je2,je1) = zsigv * zdu / fcpi
c
c       end average over cos( gyro-angle )
c
  140   continue
        i001 = nhmu-1
        do 145 jmu1 = 1, i001
        i002 = jmu1 + 1
c
c       loop over rest of jmu2, use symmetry of bbdtsv
c
        do 145 jmu2 = i002, nhmu
        bbdtsv(jmu2,jmu1,je2,je1) = bbdtsv(jmu1,jmu2,je2,je1)
  145   continue
  149   continue
c
c       end d-t cross section calculation
c
  150   continue
c
c       begin d-d cross section calculation
c
        call resetr(bbddsv,mxhmu**2*mxhe**2,0.0)
        bbdd = 0.0
c
c       zmd = 2. / mass in gm**-1
c
        zad = habeam(nhbem1(nbbspd))
        zamu = fcau * 10.0**fxau
        zmd = 2.0 / ( zamu * zad)
c
c       begin loop over velocity variables
c       zvd(1 & 2)sq = velocity**2 in (cm/sec)**2
c
        do 170 je1 = 1, nhe
        zvd1sq = (hei(je1,nbbspd)-0.5*hdei(je1,nbbspd))*uise*zmd
        zvd1 = sqrt(zvd1sq)
c
c       now have symmetry in both e and mu, so full calculation
c       needed for only half of matrix.
c
        do 170 je2 = 1, je1
        zvd2sq = (hei(je2,nbbspd)-0.5*hdei(je2,nbbspd))*uise*zmd
        zvd2 = sqrt(zvd2sq)
        do 160 jmu1 = 1, nhmu
        zmu1 = hmuc(jmu1)
        zsqmu1 = sqrt(1.0-zmu1**2)
c
c       separate inner mu loop
c
        do 160 jmu2 = 1, jmu1
        zmu2 = hmuc(jmu2)
        zsqmu2 = sqrt(1.0-zmu2**2)
c
c       do average over cos( gyro-angle )
c
        zsigv = 0.0
        do 159 ju = 1, imxju
        zvrel2 = zvd1sq + zvd2sq - 2.0*zvd1*zvd2* ( zmu1*zmu2
     1           + zsqmu1*zsqmu2*zu(ju) )
        zvrel = sqrt(zvrel2)
c
c       zeb = energy of d in kev.
c
        zeb = zvrel2 * evsinv * 0.001 / zmd
        zsigma = 0.0
c
c       d+d->he3+n cross section from b.h.duane,
c         batelle northwest report bnwl-1685
c
        if (zeb.ge.10.0.and.zeb.le.4000.0)
     1  zsigma = 1.e-24*(482./((1.+(3.08e-4*zeb-1.177)**2)
     1  *zeb*(exp(47.88/sqrt(zeb))-1.)))
        zsigv = zsigv + zsigma * zvrel * zusqrt(ju)
  159   continue
        bbddsv(jmu2,jmu1,je2,je1) = zsigv * zdu / fcpi
c
c       end average over cos( gyro-angle )
c
  160   continue
        i001 = nhmu - 1
        do 165 jmu1 = 1, i001
        i002 = jmu1 + 1
c
c       use symmetry in mu
c
        do 165 jmu2 = i002, nhmu
        bbddsv(jmu2,jmu1,je2,je1) = bbddsv(jmu1,jmu2,je2,je1)
  165   continue
  170   continue
c
c       now use symmetry in e
c
        i001 = nhe - 1
        do 179 je1 = 1, i001
        i002 = je1 + 1
        do 179 je2 = i002, nhe
        do 179 jmu1 = 1, nhmu
        do 179 jmu2 = 1, nhmu
        bbddsv(jmu2,jmu1,je2,je1) = bbddsv(jmu2,jmu1,je1,je2)
  179   continue
  180   continue
c
c       end d-d cross section calculation
c
  185   continue
c
c
        return
c
c------------------------------------------------------------------------------
c
c
cl      2)      compute neutral beam timestep (if any), scale sources
c
c               normally htbi will be  equal to tai at this point (if
c               the beams have turned on.  however, if the last timestep
c               was repeated with a shorter timestep, then the final tbi
c               for that (the previous) timestep, the present tai, will
c               be less than the value of tbi the last time hfi and the
c               sources  were computed.  that means some of the energy
c               and particles that were to have been added to the plasma
c               have not been.   we must advance hfi from time htbi to
c               tbi, and add the new sources to the old sources with
c               the obvious weighting
c
c
  200   continue
c
        if ((tbi.le.htbi).or.(nhbem1(mhsp).eq.0)) return
c
        if (k.eq.4) go to 900
c
        zdti2 = max(0.0 , tbi - htbi)
        zdti1 = max(0.0 ,  dti - zdti2)
        zdt20 = zdti2 / dti
        zdt10 = max(0.0 , 1.0 - zdt20)
c
        zth = 1.0 - thetai
c
        call scaler(halfas,mzones,zdt10)
        call scaler(halsps,mxzone*mxhsp,zdt10)
        call scaler(shbems,mzones*mxhyd,zdt10)
        call scaler(shchxs,mzones*mxhyd,zdt10)
        call scaler(shblos,mzones*mxhyd,zdt10)
        call scaler(webems,mzones,zdt10)
        call scaler(wibems,mzones,zdt10)
        call scaler(websps,mxzone*mxhsp,zdt10)
        call scaler(wibsps,mxzone*mxhsp,zdt10)
        bbdtrs = bbdtrs * zdt10
c
        call resetr(gblosi,  mxhyd ,0.0)
        bpabs = 0.0
        bploss = 0.0
c
c       calculate h(r) for all species
c
        call dposit(0)
c
c       begin main loop over species
c
        do 850 jsp = 1, mhsp
        ibeam1 = nhbem1(jsp)
        if (htspon(jsp)*ueit.gt.tbi) go to 850
c
c       set idtfus for d-t beam-target reactions
c
        idtfus = 0
        if(nhbeam(ibeam1).eq.-2.and.ltrit.gt.0)idtfus=ltrit
        if(nhbeam(ibeam1).eq.-3.and.ldeut.gt.0)idtfus=ldeut
c
c------------------------------------------------------------------------------
c
cl      3)      remap according to variation of ti at each zone,
c               dump ions whose energy <= avg plasma ion thermal energy
c               set lhemin(jz) so that hei(lhemin(jz)) > avg plasma
c               ion thermal energy
c
c
  300   continue
c
        do 349 jz = lcentr, ledge
c
        ztii = tis(2,jz) * usih
c
c
        i001 = lhemin(jz,jsp)
        if (hei(i001,jsp).le.ztii) go to 320
c
c
cl              ti has dropped, add (empty) energy levels
c
c
  304   continue
        if (lhemin(jz,jsp).le.1) go to 349
        if (hei(i001 - 1,jsp).le.ztii) go to 349
        lhemin(jz,jsp) = lhemin(jz,jsp) - 1
        i001 = lhemin(jz,jsp)
        go to 304
c
c
cl              ti has risen, remove energy levels, dump beam ions into plasma
c
c
  320   continue
c
  324   continue
        if (lhemin(jz,jsp).ge.nhe) go to 9030
c
c
        zndump = 0.0
c
        do 326 jmu = 1, nhmu
        zndump = zndump + hfi(jmu,i001,jz,jsp) * hdmuc(jmu)
        hfi(jmu,i001,jz,jsp) = 0.0
  326   continue
c
c
c               set zndump = total no of particles dumped in zone jz
c               znrate = density / sec. (standard units)
c
c
        zndump = zndump * uisd
        znrate = zndump / (dti*uist)
        zndump = zndump * dx2i(jz) *  2.0 * vols(mzones,1)
c
c
c               add dumped energy, particles to transport source terms,
c               beam particle checks
c
c
        hnplas(jz) = hnplas(jz) + zndump
        heplas(jz) = heplas(jz) + zndump * hei(i001,jsp)*uise
        i002 = lhbeam(ibeam1)
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  shbems(i002,jz) = shbems(i002,jz) + znrate
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  addi(i002) = addi(i002) + zndump
        wibems(jz) = wibems(jz) + znrate * hei(i001,jsp)*uise
        wibsps(jz,jsp) = wibsps(jz,jsp)+znrate*hei(i001,jsp)*uise
c
        lhemin(jz,jsp) = lhemin(jz,jsp) + 1
        i001 = lhemin(jz,jsp)
        if (hei(i001,jsp).le.ztii) go to 324
  349   continue
c
c------------------------------------------------------------------------------
c
c
cl      4)      coefficients
c
c      note: all coefficients are evaluated at the energy group
c       midpoints: hei(je)-0.5*hdei(je)
c
  400   continue
c
c
cl              charge exchange cross section
c
c
        hfutz(1) = max(0.0,min(hfutz(1),1.0))
        zm = 2.0 * uise * 10.0**(-fxau) / (fcau * habeam(ibeam1))
c
        do 404 je = 1, lhemax
        zvbeam = sqrt((hei(je,jsp)-0.5*hdei(je,jsp))*zm)
        zleb = log(uise*evsinv * (hei(je,jsp)-0.5*hdei(je,jsp))
     1  / habeam(ibeam1))
        zsigcx(je) = 0.6937 * (1.0 - 0.0673*zleb)**2
     1  / (1.0e+14 + 0.1112*exp( 3.3*zleb )) * zvbeam
  404   continue
c
c
c
        do 479 jz = lcentr, ledge
        zloge=23.9+log(evsinv*tes(2,jz)/sqrt(rhoels(2,jz)))
c
c
c
        i001 = lhemin(jz,jsp)
        do 459 je = i001, lhemax
        zhes = (hei(je,jsp)-0.5*hdei(je,jsp))*uise
        zvbeam = sqrt((hei(je,jsp)-0.5*hdei(je,jsp))*zm)
c
c
cl              sum dens.(js) * z(js) / atomic wt. (js)
c
c
        zefflg = 0.0
        zbarlg = 0.0
        zn0  = 0.0
c
        zlog=sqrt(tes(2,jz)/(fcpi*fces**2*rhoels(2,jz)*
     1  10.**(2.*fxes)))*zhes/(hzbeam(ibeam1)*fces**2*10.**(2.*fxes))
        zlog2=sqrt(tes(2,jz)*fcpi/rhoels(2,jz))
     1  *(2.*zvbeam*habeam(ibeam1)*fcmp/(fch*fces))*
     2  10.**(fxnucl-fxh-fxes)
        do 414 jh = 1, mhyd
        if(zhes*evsinv.lt.habeam(ibeam1)*1.e5)  zlogi=log(zlog
     1  *aspec(jh)/(habeam(ibeam1)+aspec(jh)))
        if(zhes*evsinv.ge.habeam(ibeam1)*1.e5) zlogi=
     1  log(zlog2*aspec(jh)/(aspec(jh)+habeam(ibeam1)))-0.5
        zbarlg = rhohs(jh,2,jz) * zlogi / aspec(jh) + zbarlg
        zefflg = zefflg + zlogi * rhohs(jh,2,jz)/rhoels(2,jz)
        zn0  = zn0 + rhons(jh,jz)
  414   continue
c
        if (mimp.le.0) go to 419
c
        do 418 ji = 1, mimp
        i001 = limp1-1+ji
        if(zhes*evsinv.lt.habeam(ibeam1)*c2mean(ji,2,jz)*1.e5)
     1  zlogi=log(zlog*aspec(i001)/(cmean(ji,2,jz)*
     1  (habeam(ibeam1)+aspec(i001))))
        if(zhes*evsinv.ge.habeam(ibeam1)*c2mean(ji,2,jz)*1.e5)
     1  zlogi=log(zlog2*aspec(i001)/(aspec(i001)+habeam(ibeam1)))
     2  -0.5
        zbarlg = zbarlg + rhois(ji,2,jz) * c2mean(ji,2,jz)
     1  * zlogi / aspec(i001)
        zefflg=zefflg+zlogi*rhois(ji,2,jz)*c2mean(ji,2,jz)/rhoels(2,jz)
  418   continue
c
  419   continue
c
c
cl              pitch-angle scattering rate
c
c
        hscats(je,jz,jsp) = chscat*rhoels(2,jz)*hzbeam(ibeam1)**2*
     1  zefflg / sqrt (zhes**3 * habeam(ibeam1))
c
c
cl              electron collision rate
c               note: this code assumes electron thermal velocity >> beam
c               velocity
c
c
        znue = cznue * rhoels(2,jz) * hzbeam(ibeam1)**2 * zloge /
     1                          (habeam(ibeam1) * sqrt(tes(2,jz))**3)
c
c
cl              ion collision rate * beam energy**3/2
c               note: this code assumes ion thermal velocity << beam velocity
c
c
        zgamma = czgamm * sqrt(habeam(ibeam1)) * hzbeam(ibeam1)**2 *
     1    zbarlg
c
c
cl              slowdown rate
c
c
        if (je.gt.lhemin(jz,jsp)) go to 422
c
c
        hslows(je,jz,jsp) = hfutz(5) * 3.0 * znue /
     1          log((zgamma + znue * sqrt(hei(je,jsp)*uise)**3)
     2          / (zgamma + znue*sqrt(tis(2,jz)/gamin1)**3))
c
        go to 423
c
  422   continue
c
c
        hslows(je,jz,jsp) = hfutz(5) * 3.0 * znue /
     1          log((zgamma + znue*sqrt(hei(je,jsp)*uise)**3) /
     2                  (zgamma + znue*sqrt(hei(je-1,jsp)*uise)**3))
c
c
  423   continue
c
c
cl              slow-down energy distribution (to electrons/to ions)
c
c
        z0 = (zgamma / sqrt(zhes) + znue * zhes)
c
c
cl              to ions:
c
c
        heion(je,jz,jsp) = zgamma / (sqrt(zhes) * z0)
c
c
cl              to electrons:
c
c
        heelec(je,jz,jsp) = znue * zhes / z0
c
c
c               charge exchange loss rate
c
c
        hchexs(je,jz,jsp) = zn0 * zsigcx(je) * hfutz(1)
c
c
  459   continue
c
  479   continue
c
c------------------------------------------------------------------------------
c
c
cl      5)      deposit injected ions
c
c
  500   continue
c
        do 799 jz = lcentr, ledge
c
c
cl              subroutine dposit sets the array hdeps
c               in part./cm**3/sec./unit mu
c
c
        ijz = jz
        call dposit(ijz)
c
c
c               add to heinjs, hninjs
c
c
        zdts2 = zdti2*uist
        zvols = 2.0 * vols(mzones,1) * dx2i(jz)
c
        do 510 jsrc = 1, nhsrc
        do 510 jmu = 1, nhmu
c
        heinjs(jz) = heinjs(jz) + hdeps(jmu,jsrc,jsp)
     1          * hdmuc(jmu)*zvols * zdts2 * uise*hesrci(jsrc,jsp)
        hninjs(jz) = hninjs(jz) + hdeps(jmu,jsrc,jsp)
     1          * hdmuc(jmu)*zvols * zdts2
        bpabs=bpabs+hdeps(jmu,jsrc,jsp)*hdmuc(jmu)*zvols
     1        *uise*hesrci(jsrc,jsp)
  510   continue
c
  512   continue
c
c       compute fusion cross section of beam with thermal ions
c       using function fusigv
c
        do 550 je = 1, nhe
        zebeam = (hei(je,jsp)-0.5*hdei(je,jsp))*uise*evsinv*1.0e-3
        zabeam = habeam(ibeam1)
c
c ...NOTE: bug here; should be tis * useh
c
        zti = tis(2,jz) * uieh * 1.0e-16
        hsigv(je,jsp) = 0.0
        if (zebeam.gt.1.5*zti)
     1     hsigv(je,jsp)=fusigv(3,zebeam,zabeam,zti,ierr)
  550   continue
c
c
c------------------------------------------------------------------------------
c
c
cl      6)      solve equation
c
c
c
  600   continue
c
c
cl              highest energy group--no slowdown into this group
c
c
        call resetr (zsrci,nhmu,0.0)
c
c
        i0 = lhemin(jz,jsp)
        do 749 je = i0, lhemax
        ie = lhemin(jz,jsp) + lhemax - je
c
c
cl              add deposited ions
c
c
        do 638 jsrc = 1, nhsrc
c
        if (ie.lt.lhemax) go to 622
c
c
cl              ie is highest energy group
c               deposited energy cannot exceed this energy
c
c
        if (hesrci(jsrc,jsp).gt.(hei(ie,jsp)*1.0001)) go to 9060
        if (hesrci(jsrc,jsp).gt.hei(ie-1,jsp)) go to 628
        go to 638
c
c
  622   continue
        if (ie.le.lhemin(jz,jsp)) go to 624
c
c
cl              ie is neither highest nor lowest group
c
c
        if (hesrci(jsrc,jsp).gt.hei(ie-1,jsp)
     1  .and.hesrci(jsrc,jsp).le.hei(ie,jsp)) go to 628
        if (hesrci(jsrc,jsp).gt.hei(ie,jsp)
     1  .and.hesrci(jsrc,jsp).lt.hei(ie+1,jsp)) go to 626
        go to 638
c
c
  624   continue
c
c
cl              lowest energy group
c               this energy must not exceed deposited energy
c
c
        if (hesrci(jsrc,jsp).lt.hei(ie,jsp)) go to 638
        if (hesrci(jsrc,jsp).ge.hei(ie+1,jsp)) go to 638
c
c
cl              deposited ions have more energy than group ie
c
c
  626   continue
        zint = (hei(ie+1,jsp) - hesrci(jsrc,jsp)) / hdei(ie+1,jsp)
        go to 630
c
c
cl              deposited ions have less or equal energy than group ie
c
c
  628   continue
        zint = (hesrci(jsrc,jsp) - hei(ie-1,jsp)) / hdei(ie,jsp)
c
c
cl              add deposited ions to source term
c
c
  630   continue
c
        do 634 jmu = 1, nhmu
        zsrci(jmu) = zsrci(jmu) + hdeps(jmu,jsrc,jsp)*usid*
     1                                          uist*zint
  634   continue
c
  638   continue
c
c
c
cl              difference and reduce d/dmu equation
c
c
c
        do 679 jmu = 1, nhmu
c
        z0 = zdti2 * uist * hscats(ie,jz,jsp) / hdmuc(jmu)
        z1 = hsin(jmu) * z0
        z3 = hsin(jmu+1) * z0
        znudt=0.0
        if(idtfus.gt.0) znudt=hsigv(ie,jsp)*rhohs(idtfus,2,jz)
        znutot=hslows(ie,jz,jsp)+hchexs(ie,jz,jsp)+znudt
        z2 = - (z1 + z3 + zdti2*uist*znutot)
c
        zd = 0.0
        za = 0.0
        zc = 0.0
c
c
        if (jmu.le.1) go to 644
c
cl              if not mu = -1.0 boundary
 
c
        za = - z1 * thetai
        zd = zd + zth * z1 * hfi(jmu-1,ie,jz,jsp)
c
  644   continue
c
        zb = 1.0 - z2 * thetai
        zd = zd + (1.0+z2*zth)*hfi(jmu,ie,jz,jsp)+zsrci(jmu)*zdti2
c
        if (jmu.ge.nhmu) go to 648
c
c
cl              if not mu = +1.0 boundary
c
c
        zc = - z3 * thetai
        zd = zd + z3*zth*hfi(jmu+1,ie,jz,jsp)
c
  648   continue
c
c
cl              reduce to  hfi(jmu) + ze(jmu)*hfi(jmu) = zf(jmu)
c
c
        if (jmu.gt.1) go to 652
c
cl              first mu-group - assume za = 0.0
c
        if (abs(zb).le.epslon) go to 9063
        ze(1) = zc / zb
        zf(1) = zd / zb
        go to 679
c
  652   continue
        if (jmu.ge.nhmu) go to 656
c
cl              intermediate mu groups
c
        z0 = zb - za*ze(jmu-1)
        if (abs(z0).le.epslon) go to 9064
        ze(jmu) = zc / z0
        zf(jmu) = (zd - za*zf(jmu-1)) / z0
        go to 679
c
cl              last (nhmu) mu-group
c
  656   continue
c
        z0 = zb - za*ze(jmu-1)
        if (abs(z0).le.epslon) go to 9065
        ze(jmu) = 0.0
        zf(jmu) = (zd - za*zf(jmu-1)) / z0
c
  679   continue
c
c
cl              solve reduced d/dmu equation,
c               add in source terms to transport equation,
c               add to beam conservation check variables,
c               compute slow-down contribution to zsrci
c
c
        do 699 jmu = 1, nhmu
        imu = nhmu + 1 - jmu
c
        zhfold = hfi(imu,ie,jz,jsp)
        zhfnew = zf(imu)
        if (imu.lt.nhmu) zhfnew = zf(imu) - ze(imu)*hfi(imu+1,ie,jz,jsp)
c
cl              use thetai-time-centered hfi
c               to compute sources and sinks
c
        zhf = zth*zhfold + thetai*zhfnew
        zvols = dx2i(jz)* 2.0 * vols(mzones,1)
c
        zeplas = zhf*hslows(ie,jz,jsp)*hdei(ie,jsp)*hdmuc(imu)
     1           *uisd*uise
c
        z1 = zeplas * heelec(ie,jz,jsp) * zdt20
        webems(jz) = webems(jz) + z1
        websps(jz,jsp) = websps(jz,jsp) + z1
        z1 = zeplas * heion(ie,jz,jsp) * zdt20
        wibems(jz) = wibems(jz) + z1
        wibsps(jz,jsp) = wibsps(jz,jsp) + z1
        heplas(jz) = heplas(jz) + zeplas*uist*zdti2*zvols
        zsrci(imu) = zhf * hslows(ie,jz,jsp) * uist
        z1=znudt*zhf*hdmuc(imu)*uisd*zdt20
        if(idtfus.gt.0) shbems(idtfus,jz)=shbems(idtfus,jz)-z1
        halfas(jz)=halfas(jz)+z1
        halsps(jz,jsp)=halsps(jz,jsp)+z1
c
        do 684 jh = 1, mhyd
        shblos(jh,jz) = shblos(jh,jz) +
     1  zdt20*zhf* uisd * hfutz(1)  * zsigcx(ie) * rhons(jh,jz) *
     2  hdmuc(imu)
  684   continue
c
        z0 = zhf * hchexs(ie,jz,jsp) * hdmuc(imu)
     1       * uisd*uist*zdti2*zvols
        z1=zhf*znudt*hdmuc(imu)*uisd*uist*zdti2*zvols
        heloss(jz) = heloss(jz) + (z0+z1) * hei(ie,jsp) * uise
        hnloss(jz) = hnloss(jz) + (z0+z1)
        bploss=bploss+uise*hei(ie,jsp)*(z0+z1)/(uist*zdti2)
c
        hfi(imu,ie,jz,jsp) = zhfnew
c
  699   continue
c
  749   continue
c
c
cl              slowdown out of lowest energy group into plasma
c               if lhemin(jz) = 1, all their energy has been dumped
c               through slowdown, otherwise dump in hei(lhemin(jz)-1)
c               for each particle.
c               also dump particles
c
c
        do 759 jmu = 1, nhmu
c
        znplas = zsrci(jmu) * uisd * usit * hdmuc(jmu)
        i002 = lhbeam(ibeam1)
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  shbems(i002,jz) = shbems(i002,jz) + znplas * zdt20
        hnplas(jz) = hnplas(jz) + znplas * uist * zdti2 * zvols
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  addi(i002) = addi(i002) + znplas * uist * zdti2 * zvols
c
        if (lhemin(jz,jsp).le.1) go to 759
c
        i001 = lhemin(jz,jsp)
        zeplas = znplas * hei(i001-1,jsp)*uise
        wibems(jz) = wibems(jz) + zeplas * zdt20
        wibsps(jz,jsp)=wibsps(jz,jsp)+zeplas*zdt20
        heplas(jz) = heplas(jz) + zeplas * uist * zdti2 * zvols
c
  759   continue
c
c
c               injection groups below hei(lhemin(jz)) get
c               dumped directly into the plasma
c
c
  760   continue
c
        do 769 jsrc = 1, nhsrc
c
        i001 = lhemin(jz,jsp)
        if (hesrci(jsrc,jsp).ge.hei(i001,jsp)) go to 769
        znplas = 0.0
c
        do 763 jmu = 1, nhmu
        znplas = znplas + hdeps(jmu,jsrc,jsp)*hdmuc(jmu)
  763   continue
c
        zeplas = znplas * hesrci(jsrc,jsp)*uise
c
c               add into cons. vars.
c
        hnplas(jz) = hnplas(jz) + znplas * zdti2*uist * zvols
        heplas(jz) = heplas(jz) + zeplas * zdti2*uist * zvols
c
        i002 = lhbeam(ibeam1)
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  shbems(i002,jz) = shbems(i002,jz) + znplas
        if (i002.ge.lhyd1.and.i002.le.lhydn)
     1  addi(i002) = addi(i002) + znplas * uist * zdti2 * zvols
        wibems(jz) = wibems(jz) + zeplas
        wibsps(jz,jsp) = wibsps(jz,jsp) + zeplas
  769   continue
c
c
c               set gblosi
c
c
        if (lhbeam(ibeam1).le.0) go to 775
        i002 = lhbeam(ibeam1)
c
        do 774 jh = 1, mhyd
        gblosi(i002) = gblosi(i002) +
     1          shblos(jh,jz)*dx2i(jz)
  774   continue
c
  775   continue
c
  799   continue
c
        htbi = tbi
c
c
c------------------------------------------------------------------------------
c
c
cl      8)      compute d+d->he3+n fusion rate (for printout only)
c
c
  800   continue
c
        if (nhbeam(ibeam1).ne.-2) go to 850
        if (ldeut.le.0) go to 850
c
        zvols = vols(mzones,1)
c
        ztdd = 0.0
c
        do 849 jz = lcentr, ledge
        zsdd = 0.0
        zeloss = 0.0
c
        do 839 je = 1, nhe
c
c       compute fusion cross section of beam with thermal ions
c       for d - d reaction using function fusigv
c
        zebeam = (hei(je,jsp)-0.5*hdei(je,jsp))*uise*evsinv*1.0e-3
        zabeam = habeam(ibeam1)
c
c ...NOTE: bug here; should be tis * useh
c
        zti = tis(2,jz) * uieh * 1.0e-16
        zddsig(je) = 0.0
        if (zebeam.gt.1.5*zti)
     1     zddsig(je)=fusigv(2,zebeam,zabeam,zti,ierr)
c
        if (ldeut.gt.0) zrfdd = rhohs(ldeut,2,jz) * zddsig(je)
c
        do 829 jmu = 1, nhmu
        zsdd = zsdd + hfi(jmu,je,jz,jsp) * zrfdd * uisd * hdmuc(jmu)
  829   continue
  839   continue
c
        ztdd = ztdd + zsdd * dx2i(jz)
  849   continue
c
        hdds = ztdd * 2.0*zvols
        hddtot = hddtot + hdds*zdti2*uist
c
c
c
  850   continue
c
c
        call scaler(gblosi,mxhyd,uist*uisl**2*rmins)
c
c       calculate beam-beam reaction rates
c       add d-t results to halfas, etc.
c       first, do d-t; skip if don't have d and t
c
        if (nbbsp1.eq.0.or.nbbsp2.eq.0) go to 875
        if (htspon(nbbsp1)*ueit.gt.tbi
     1      .or.htspon(nbbsp2)*ueit.gt.tbi) go to 875
c
        zvols = vols(mzones,1)
c
c       loop over radius for volume integral
c
        do 870 jz = lcentr, ledge
        zbbdt = 0.0
c
c       loop over velocity variables for each species
c
        i001 = lhemin(jz,nbbsp1)
        i002 = lhemin(jz,nbbsp2)
        do 865 je1 = i001, nhe
        do 865 je2 = i002, nhe
        do 865 jmu1 = 1, nhmu
        zf1 = hfi(jmu1,je1,jz,nbbsp1)*hdmuc(jmu1)*uisd
        do 865 jmu2 = 1, nhmu
c
c       compute reaction rate using hfi's and cross section, bbdtsv
c
        zbbdt = zbbdt + zf1*hfi(jmu2,je2,jz,nbbsp2)*hdmuc(jmu2)
     1          * uisd * bbdtsv(jmu2,jmu1,je2,je1) * zdt20
  865   continue
        halfas(jz) = halfas(jz) + zbbdt
        halsps(jz,nbbsp1) = halsps(jz,nbbsp1) + 0.5 * zbbdt
        halsps(jz,nbbsp2) = halsps(jz,nbbsp2) + 0.5 * zbbdt
        bbdtrs = bbdtrs + zbbdt*dx2i(jz)
  870   continue
c
c       bbdtrs = volume-integrated reaction rate in standard units, sec**-1
c
        bbdtrs = bbdtrs*2.0*zvols
c
c       bbdt = total no. of reactions
c
        bbdt = bbdt + bbdtrs*zdti2*uist
c
c       end d-t calculation
c
  875   continue
c
c       begin d-d beam-beam calculation
c       skip if have no d
c
        if (nbbspd.eq.0) go to 899
        if (htspon(nbbspd)*ueit.gt.tbi) go to 899
c
        bbddrs = 0.0
        zvols = vols(mzones,1)
c
c       begin radial loop
c
        do 890 jz = lcentr, ledge
        zbbdd = 0.0
c
c       loop over velocity variables of two reactants
c
        i001 = lhemin(jz,nbbspd)
        do 885 je1 = i001, nhe
        do 885 je2 = i001, nhe
        do 885 jmu1 = 1, nhmu
        zf1 = hfi(jmu1,je1,jz,nbbspd)*hdmuc(jmu1)*uisd
        do 885 jmu2 = 1, nhmu
c
c       compute reaction rate using cross section, bbddsv
c
        zbbdd = zbbdd + zf1*hfi(jmu2,je2,jz,nbbspd)*hdmuc(jmu2)
     1          * uisd * bbddsv(jmu2,jmu1,je2,je1)
  885   continue
c
c       note that a factor of 0.5 is introduced since the
c       deuterium velocity space has been integrated over twice
c
        bbddrs = bbddrs + 0.5 * zbbdd *dx2i(jz)
  890   continue
c
c       bbddrs = volume-integrated reaction rate in standard units, sec**-1
c
        bbddrs = bbddrs*2.0*zvols
c
c       bbdd = total no. of reactions
c
        bbdd = bbdd + bbddrs*zdti2*uist
c
c       end d-d and beam-beam reaction rate calculations
c
c
  899   continue
c
c------------------------------------------------------------------------------
c
c
cl      9)      set rhobis, rhobes, hebems, ajbs
c
c
  900   continue
c
        call resetr (rhobis,2*mxzone,0.0)
        call resetr (rhobes,2*mxzone,0.0)
        call resetr (hebems,mxzone,0.0)
        call resetr (ajbs,mxzone,0.0)
c
c       loop over species
c
        do 950 jsp = 1, mhsp
        ibeam1 = nhbem1(jsp)
        if (htspon(jsp)*ueit.gt.tbi) go to 950
        z0 = fces * sqrt(uise * 2.0 / (fcau * habeam(ibeam1))) *
     1                                  uisd * 10**(fxes - 0.5*fxau)
c
c
        do 949 jz = lcentr, ledge
        zbn = 0.0
        zbe = 0.0
        zbj = 0.0
        i001=lhemin(jz,jsp)
c
        do 939 je = i001, nhe
        zbn0 = 0.0
        zbj0 = 0.0
c
        do 929 jmu = 1, nhmu
        zbn0 = zbn0 + hfi(jmu,je,jz,jsp) * hdmuc(jmu)
        zbj0 = zbj0 + hfi(jmu,je,jz,jsp) * hmuc(jmu) * hdmuc(jmu)
  929   continue
c
        zbn = zbn + zbn0
        zbe = zbe + zbn0 * hei(je,jsp)
        zbj = zbj + zbj0 * sqrt(hei(je,jsp)) * z0
  939   continue
c
        rhobis(2,jz) = rhobis(2,jz) + zbn * uisd
        rhobes(2,jz) = rhobes(2,jz) + zbn * uisd * hzbeam(ibeam1)
        hebems(jz) = hebems(jz) + zbe * uise * uisd
c
c       elliptic integral business
c       switches on when hfutz(4)=1., off when hfutz(4)=0. default:0.
c       assumes hangle(m,n) = default except for hangle(3,ibeam1)
c
c       do the geometry -- zlttlr is 'little r', the radial
c       coordinate, and zrtang is the perpendicular distance
c       from the magnetic axis to the neutral beam path
c       zepsln is the ratio of little r to big r, where big r is
c       the major radius
c
c     rgb 14-nov-85  this is probably where the effect of plasma shift
c  on neutral beam penetration should be included
c
        zlttlr = ahalfs(jz,2)  ! was rmins * xzoni(jz)
        zrtang = sin(fcpi/180.*hangle(3,ibeam1))
     1  *(hrmaj(ibeam1)+hrmin(ibeam1)*cos((fcpi/180.)*hangle(1,ibeam1)))
        zepsln = ahalfs(jz,2) / rmids(jz,2)  ! was zlttlr/rmajs
c
c       zmubsq is mu, the cosine-squared of the angle between the
c       toroidal field and the neutral beam path
        zmub = zrtang / (rmids(jz,2) + ahalfs(jz,2)) ! was (rmajs + zlttlr)
        zmubsq = zmub * zmub
        if (zmubsq .le. epslon) go to 944
c
c       get zm1, the argument for the integral
        zm1 = 1. - (2.*zepsln)/((1.+zepsln)*zmubsq)
        if (zm1 .le. epslon) go to 944
c
c       now get zkm, a polynomial approximation to k(m)
c       ( =k(1.-zm1)), the complete elliptic integral of the second
c       kind
c       multiply its reciprocal by pi over 2. -- so zkm ranges from
c       0. to 1.
        zkm = (((zm1*.07253)+.11197)*zm1+1.38629)
     1          + (((zm1*.02887)+.12135)*zm1+.5)*log(1./zm1)
        zkm = .5 * fcpi / zkm
        go to 945
c
c       zkm is set to zero for beam angles such that zm1 is less than
c       zero
  944   zkm = 0.0
c
c       hfutz(4) switches zkm on/off; when off, zkm = 1.0
  945   continue
        zkm = zkm * hfutz(4)
        if (hfutz(4) .le. epslon) zkm = 1.0
c
c       end of elliptic integral business
c
        zjfact = 1.0 - (1.0 -
     1  (sqrt(ahalfs(jz,2)/rmids(jz,2)) * (1.55 + 0.85/xzeff(2,jz))-
     2  (ahalfs(jz,2)/rmids(jz,2))*(.20+1.55/xzeff(2,jz)))/
     2          (1.0 + hfutz(2)*xnuel(2,jz)))
     3  / xzeff(2,jz)
        zjfact = zjfact * zkm
        zjfact = zjfact * (1. - hfutz(3)*xzeff(2,jz)*sqrt(zepsln))
        ajbs(jz) = ajbs(jz) + zbj * zjfact
  949   continue
c
  950   continue
c
        do 955 jz = lcentr, ledge
        if(rhobis(2,jz).gt.epslon) hebems(jz)=hebems(jz)/rhobis(2,jz)
  955   continue
        rhobis(2,1) = rhobis(2,lcentr)
        rhobes(2,1) = rhobes(2,lcentr)
        hebems(1) = hebems(lcentr)
        ajbs(1) = ajbs(lcentr)
c
c..extrapolate to zone center just outside plasma
c
        z1 = ( xzoni(ledge) - xzoni(ledge+1) )
     &    / ( xzoni(ledge) - xzoni(ledge-1) )
c
        rhobis(2,ledge+1) = max ( epslon, z1 * rhobis(2,ledge-1)
     &                        + ( 1.0 - z1 ) * rhobis(2,ledge) )
        rhobes(2,ledge+1) = max ( epslon, z1 * rhobes(2,ledge-1)
     &                        + ( 1.0 - z1 ) * rhobes(2,ledge) )
        hebems(ledge+1)   = max ( epslon, z1 * hebems(ledge-1)
     &                        + ( 1.0 - z1 ) * hebems(ledge) )
        ajbs(ledge+1)     =  z1 * ajbs(ledge-1)
     &                        + ( 1.0 - z1 ) * ajbs(ledge)

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
        return
c
c
c------------------------------------------------------------------------------
c
c
cl      10)     compress beam particles
c
c
 1000   continue
        if(.not.nlcomp) return
c
        do 1001 jb = 1, mhbeam
        if (hton(jb)*ueit.lt.tbi) go to 1002
 1001   continue
c
        return
c
 1002   continue
        if(abs(compl-1.0).lt.rndeps) go to 1060
        if(compl.gt.1.0) go to 1030
c
c
c               limiter moving out rel. to plasma --
c               pad out old last zone (dilute with vacuum)
c
c               iedge = old ledge .le. (new) ledge
c               dvolx is for zone iedge
c               reset lhemin in new zones, i.e.,
c               zones  jz,  iedge+1 .le. jz .le. ledge
c
c
        iedge = ledge - nzcomp
        zdilut = 1.0 - dvolx * dx2inv(iedge)
c
        do 1008 jsp = 1, mhsp
        do 1008 je = 1, nhe
        do 1008 jmu = 1, nhmu
        hfi(jmu,je,iedge,jsp) = hfi(jmu,je,iedge,jsp) * zdilut
 1008   continue
c
        if (nzcomp.eq.0) go to 1060
c
        i0 = iedge + 1
        do 1010 jsp = 1, mhsp
        do 1010 jz = i0, ledge
        lhemin(jz,jsp) = 1
 1010   continue
c
        go to 1060
c
c
c
c               limiter moving in rel. to plasma --
c               cut away outer zones, note loss of energy,
c               particles in hecmps, hncmps
c
c               iedge = old ledge .ge. (new) ledge
c               dvolx is for zone ledge
c
c
 1030   continue
        iedge = ledge - nzcomp
        i0 = ledge + 1
        zvoli = vols(mzones,1) * usil**3
c
        if (nzcomp.eq.0) go to 1040
c
        do 1038 jz = i0, iedge
c
        do 1034 jsp = 1, mhsp
        do 1034 je = 1, nhe
        do 1034 jmu = 1, nhmu
        hfi(jmu,je,jz,jsp) = 0.0
 1034   continue
c
        hncmps(jz) = hninjs(jz) - hnplas(jz) - hnloss(jz)
        hecmps(jz) = heinjs(jz) - heplas(jz) - heloss(jz)
 1038   continue
c
c
c               partially cut away zone
c
c
 1040   continue
        zv = zvoli * dvolx * 2.0 * compp**2
c
        do 1044 jsp = 1, mhsp
        do 1044 je = 1, nhe
        do 1044 jmu = 1, nhmu
        hncmps(ledge) = hncmps(ledge)
     1                  - hfi(jmu,je,ledge,jsp)*hdmuc(jmu)*zv
        hecmps(ledge) = hecmps(ledge)
     1     - hfi(jmu,je,ledge,jsp)*hdmuc(jmu)*zv*hei(je,jsp) * uise
 1044   continue
c
c
c               compress density, compute theoretical compression
c               energization
c
c
 1060   continue
        zvoli = vols(mzones,1) * usil**3
c
        do 1078 jsp = 1,mhsp
        do 1078 jz = lcentr, ledge
        zepar = 0.0
        zeperp = 0.0
c
        i0 = lhemin(jz,jsp)
        do 1068 je = i0, lhemax
        do 1066 jmu = 1, nhmu
        hfi(jmu,je,jz,jsp) = hfi(jmu,je,jz,jsp) * compp**2
        zepar = zepar + hfi(jmu,je,jz,jsp) * hei(je,jsp) * hdmu3(jmu)
        zeperp = zeperp + hfi(jmu,je,jz,jsp) * hei(je,jsp) *
     1          (hdmuc(jmu) - hdmu3(jmu))
 1066   continue
 1068   continue
c
        hecmps(jz) = hecmps(jz) + uise * zvoli * 2.0 * dx2i(jz) *
     1          (zepar * (compp**2 - 1.0) + zeperp * (compp - 1.0))
 1078   continue
c
c
c               new e-grid for all species
c
c
        do 1149 jsp = 1, mhsp
c
        zei(1) = hei(1,jsp) * compp * (1.0 - rndeps)
        if (nhe.lt.2) go to 9100
c
        z0 = (compp * (1.0 + rndeps))**2
c
        do 1084 je = 2, nhe
        zei(je) = hei(je,jsp) * z0
 1084   continue
c
c
c
c       for some mu, new e = old e *(c**2 * mu**2 + c*(1 - mu**2)).
c       zcomp is the factor of 'old e' averaged over the mu-group.
c       old hfi(e,mu) is assumed to be = new hfi(e*zcomp,new mu).
c       the e-grid is changed to be new e(je) = old e(je) * c**2,
c       except new e(1) = old e(1) * c
c       ('c' is the compression factor)
c       the new mu of a mu,e group is independent of the energy,
c       so the mu-grid is changed so that no remapping is necessary
c       accross mu-groups.
c       we need only map new hfi(e*zcomp,mu) onto the new
c       e-grid. this is done analagous to the mapping of
c       hdeps into hfi.
c       note that repeated remapping does not diffuse particles up
c       the e-grid, only down (except between e-groups 1 and 2)
c       the mu-remapping requires only changing the mu-grid
c       and multiplying hfi * old hdmuc / new hdmuc.
c       this is done by multiplying hfi by old hdmuc when remapping
c       the e-groups (so that hfi(jmu,je,jz) is a density),
c       and dividing by new hdmuc when new hdmuc has been computed.
c
c
c
 1100   continue
c
        do 1148 jz = lcentr, ledge
c
        iemin = lhemax
c
        do 1138 jmu = 1, nhmu
c
c               energization factor
c
        zcomp = (compp**2 - compp) * 0.333333333 *
     1  (hmub(jmu)**2 + hmub(jmu)*hmub(jmu+1) + hmub(jmu+1)**2) + compp
c
        do 1114 je = 1, nhe
        zfi(je) = 0.0
 1114   continue
c
c
c       *       remap hfi(zener,mu)
c
c
        i001 = lhemin(jz,jsp)
        do 1130 je = i001, lhemax
        zener = hei(je,jsp) * zcomp
c
c               find e-groups to map zener group into
c
        do 1124 je2 = 1, nhe
        ie2 = je2
        if(zener.ge.zei(je2).and.zener.lt.zei(je2+1)) go to 1126
 1124   continue
c
        go to 9110
 1126   continue
c
        zfi(ie2) = zfi(ie2) + hfi(jmu,je,jz,jsp) *
     1          (zei(ie2+1) - zener)/(zei(ie2+1) - zei(ie2))
        zfi(ie2+1) = zfi(ie2+1) + hfi(jmu,je,jz,jsp) *
     1          (zener - zei(ie2))/(zei(ie2+1) - zei(ie2))
 1130   continue
c
c               copy new hfi into hfi array, find lowest e-group
c               at this zone
c
        do 1134 je = 1, nhe
        ie = 1 + nhe - je
        hfi(jmu,ie,jz,jsp) = zfi(ie) * hdmuc(jmu)
        if (ie.lt.iemin.and.zfi(ie).gt.epslon) iemin = ie
 1134   continue
c
 1138   continue
c
        lhemin(jz,jsp) = iemin
 1148   continue
c
c               e-variables
c
        do 1149 je = 1, nhe
        hei(je,jsp) = zei(je)
        hdei(je,jsp) = zei(je)
        if(je.gt.1) hdei(je,jsp) = hdei(je,jsp) - hei(je-1,jsp)
 1149   continue
c
c               new e- and mu-mesh variables
c
c
 1150   continue
c
c               new hmub
c
        i = nhmu/2 + 1
        do 1154 jmu = 1, i
        imu = 2 + nhmu - jmu
        hmub(jmu) = - sqrt(compp * hmub(jmu)**2/
     1          ((compp - 1.0) * hmub(jmu)**2 + 1.0))
        hmub(imu) =   sqrt(compp * hmub(imu)**2/
     1          ((compp - 1.0) * hmub(imu)**2 + 1.0))
 1154   continue
c
c               related mu-variables
c
        do 1158 jmu = 1, nhmu
        hmuc(jmu) = (hmub(jmu+1) + hmub(jmu))*0.5
        hdmuc(jmu) = hmub(jmu+1) - hmub(jmu)
        hdmu3(jmu) = (hmub(jmu+1)**3 - hmub(jmu)**3) * 0.333333333
c
        if (jmu.le.1) go to 1158
        hdmub(jmu) = hmuc(jmu) - hmuc(jmu-1)
        hsin(jmu) = (1.0 - hmub(jmu)**2) / hdmub(jmu)
 1158   continue
c
c               convert hfi back to dens. / unit mu
c
c
        do 1178 jsp = 1,mhsp
        do 1178 jmu = 1, nhmu
        z0 = 1.0 / hdmuc(jmu)
c
        do 1176 jz = lcentr, ledge
        do 1174 je = 1, nhe
        hfi(jmu,je,jz,jsp) = hfi(jmu,je,jz,jsp) * z0
 1174   continue
 1176   continue
c
 1178   continue
c
c
c
c               recompute rhobes, rhobis, hebems
c
        go to 900
c
c
c------------------------------------------------------------------------------
c
c
cl      90)     errors
c
c
 9030   continue
        call error_olymp(1,iclass,isub,3,
     1          'mean plasma energy greater than max. beam energy ')
c
        nlend = .false.
        call endrun
        stop
c
 9060   continue
        call error_olymp(1,iclass,isub,6,
     1          'highest beam group energy less than inj. energy ')
        call error_olymp(3,jz,2,0,8h  zone  )

        nlend = .false.
        call endrun
        stop
c
cahk
c 9061   continue
c        call error_olymp(1,iclass,isub,6,
c     1          'lowest beam group energy gr. than inj. energy ')
c        call error_olymp(3,jz,2,0,8h  zone  )
c
c        nlend = .false.
c        call endrun
c        stop
c
 9063   continue
        call error_olymp(1,iclass,isub,6,
     >	     ' cannot reduce d/dmu equation ')
        call error_olymp(3,jz,2,0,'  zone  ')
        call error_olymp(3,ie,2,0,'beam grp')
        call error_olymp(3,1,2,0, 'mu-group')
        call error_olymp(2,zb,1,0,'  zb    ')
        return
c
 9064   continue
        call error_olymp(1,iclass,isub,6,
     1          ' can''t reduce d/dmu equation ')
        call error_olymp(3,jz,2,0,'  zone  ')
        call error_olymp(3,ie,2,0,'beam grp')
        call error_olymp(3,jmu,2,0,'mu-group')
        call error_olymp(2,z0,1,0,'  z0    ')
        return
 9065   continue
        call error_olymp(1,iclass,isub,6,
     1          ' can''t reduce d/dmu equation ')
        call error_olymp(3,jz,2,0,'  zone  ')
        call error_olymp(3,ie,2,0,'beam grp')
        call error_olymp(3,jmu,2,0,'mu-group')
        call error_olymp(2,z0,1,0,'  z0    ')
        return
c
c
c               nhe .lt.2
c
 9100   continue
        call error_olymp(0,iclass,isub,10,
     1          'nhe .lt. 2, therefor can''t compress beam ')
        return
c
c               can't remap zener
c
 9110   continue
        call error_olymp(1,iclass,isub,11,
     1          ' remapped e won''t fit on e grid ')
        call error_olymp(3,jmu,2,0,'mu-group')
        call error_olymp(3,je,2,0,'e-group ')
        call error_olymp(3,zener,1,0,'new e   ')
        call error_olymp(2,zei,5,nhe,'e-grid  ')
        return
c
c
        end
c
c
c******************************************************************************
c******************************************************************************
c/ module fusigv
c
c.              %%%%
c.              %
c.               0.000000usigv          %
c.              %
c.              %%%%
c
       function fusigv(kreac,pebeam,pabeam,ption,kerr)
c
c  kreac denotes the fusion reaction:
c  1       d + d   -> t + p
c  2       d + d   -> he3 + n
c  3       d + t   -> he4 + n
c  4       d + he3 -> he4 + p
c
c  pebeam is the fast-ion energy in kev
c
c  pabeam is the atomic mass of the fast ion
c
c  ption is the temperature (in kev) of the thermal target species
c
c  kerr is an error code set by fusigv if there is an input error
c
c
c  <sigma-v> were calculated from the fits to sigma given by
c  duane, bnwl-1685 (1972); then fit to rational polynomials
c  in tion by using the imsl fitting routine zxssq
c
c
       dimension zredeu(30,4),zrsigv(30,4),zafit(30,4),zbfit(30,4),
     1 zcfit(30,4),zdfit(30,4),zxfit(10),zyfit(10),zdfrnc(10,10),
     2 iremax(4)
c
c  xfit, yfit, and zdfrnc are dimensioned for inewtn up to 10
       data inewtn/4/
c
       data (iremax(l),l=1,4)/25,25,30,25/
c
       data (zredeu(l,1),zrsigv(l,1),zafit(l,1),zbfit(l,1),
     1 zcfit(l,1),zdfit(l,1),l= 1,17)/
     2 10.,  6.855e-22,   45.6243,  491.7580,   -0.0263,    0.1571,
     3 11.,  1.289e-21,   41.6274,  380.2556,    0.2476,    0.0333,
     4 12.,  2.231e-21,   38.2332,  302.5680,    0.4830,   -0.0554,
     5 14.,  5.547e-21,   32.8291,  204.2789,    0.8577,   -0.1579,
     6 17.,  1.576e-20,   27.0591,  126.4253,    1.2344,   -0.2004,
     7 20.,  3.481e-20,   22.9921,   84.7980,    1.4462,   -0.1782,
     8 25.,  9.265e-20,   18.2770,   47.3609,    1.5213,   -0.0953,
     9 32.,  2.397e-19,   13.6178,   17.7425,    0.9749,    0.0216,
     x 40.,  5.087e-19,   39.9706,  290.4696,   30.1380,   -0.6893,
     1 50.,  9.872e-19,   11.6348,   27.2896,    4.0744,    0.1138,
     2 63.,  1.804e-18,    8.6938,   13.3834,    2.9226,    0.1238,
     3 80.,  3.110e-18,    6.6021,    7.0059,    2.2339,    0.0737,
     4 100.,  4.855e-18,    5.0132,    3.5586,    1.6338,   -0.0001,
     5 120.,  6.714e-18,    3.8109,    1.5096,    1.0563,   -0.0752,
     6 140.,  8.614e-18,    2.8719,    0.2407,    0.5392,   -0.1268,
     7 170.,  1.146e-17,    1.7247,   -1.0139,   -0.1914,   -0.1657,
     8 200.,  1.424e-17,    0.8299,   -1.7643,   -0.8211,   -0.1572/
       data (zredeu(l,1),zrsigv(l,1),zafit(l,1),zbfit(l,1),
     1 zcfit(l,1),zdfit(l,1),l=18,25)/
     2 250.,  1.867e-17,   28.9279,   40.4410,   27.5794,    1.7854,
     3 320.,  2.443e-17,    4.0760,    4.0668,    2.9001,    0.2301,
     4 400.,  3.048e-17,    1.9838,    1.7559,    0.8963,    0.3299,
     5 500.,  3.746e-17,    1.5629,    1.6408,    0.5029,    0.5689,
     6 630.,  4.592e-17,    1.7412,    2.3894,    0.6499,    1.0638,
     7 800.,  5.640e-17,    2.2588,    3.9476,    1.0750,    2.1786,
     8 900.,  6.243e-17,    2.5037,    4.7232,    1.2547,    2.9649,
     9 1000.,  6.841e-17,    2.6679,    4.9570,    1.3538,    3.5521/
       data (zredeu(l,2),zrsigv(l,2),zafit(l,2),zbfit(l,2),
     1 zcfit(l,2),zdfit(l,2),l= 1,17)/
     2 10.,  5.271e-22,   48.9829,  606.2766,   -0.2490,    0.2746,
     3 11.,  1.017e-21,   44.7745,  465.8257,    0.0365,    0.1289,
     4 12.,  1.801e-21,   41.1662,  368.7715,    0.2851,    0.0203,
     5 14.,  4.651e-21,   35.3798,  247.1536,    0.6886,   -0.1151,
     6 17.,  1.380e-20,   29.1729,  152.1353,    1.1113,   -0.1914,
     7 20.,  3.152e-20,   24.8010,  102.1260,    1.3721,   -0.1892,
     8 25.,  8.745e-20,   19.7758,   58.1954,    1.5444,   -0.1204,
     9 32.,  2.355e-19,   15.0448,   26.4519,    1.2812,   -0.0046,
     x 40.,  5.162e-19,    9.8904,   -6.8506,   -0.7309,    0.0671,
     1 50.,  1.031e-18,   13.1824,   37.4219,    4.9570,    0.1457,
     2 63.,  1.932e-18,    9.4611,   16.3376,    3.1941,    0.1756,
     3 80.,  3.408e-18,    7.1826,    8.5389,    2.4563,    0.1392,
     4 100.,  5.416e-18,    5.4710,    4.4138,    1.8369,    0.0608,
     5 120.,  7.582e-18,    4.3240,    2.4063,    1.3826,   -0.0004,
     6 140.,  9.816e-18,    3.2418,    0.6647,    0.7747,   -0.0972,
     7 170.,  1.318e-17,    2.0345,   -0.7805,    0.0396,   -0.1815,
     8 200.,  1.647e-17,    2.3218,    0.5971,    0.6277,   -0.0020/
       data (zredeu(l,2),zrsigv(l,2),zafit(l,2),zbfit(l,2),
     1 zcfit(l,2),zdfit(l,2),l=18,25)/
     2 250.,  2.172e-17,    1.4841,   -0.0267,    0.1082,    0.0000,
     3 320.,  2.849e-17,   10.8137,   12.1869,    9.7101,    0.7654,
     4 400.,  3.550e-17,    4.1521,    3.8351,    3.1773,    0.2724,
     5 500.,  4.341e-17,    2.0973,    1.8253,    1.1923,    0.3135,
     6 630.,  5.270e-17,    1.4073,    1.4117,    0.5212,    0.4902,
     7 800.,  6.374e-17,    1.4204,    1.7980,    0.5013,    0.8411,
     8 900.,  6.987e-17,    1.5987,    2.2539,    0.6445,    1.1513,
     9 1000.,  7.582e-17,    1.8181,    2.8108,    0.8225,    1.5564/
       data (zredeu(l,3),zrsigv(l,3),zafit(l,3),zbfit(l,3),
     1 zcfit(l,3),zdfit(l,3),l= 1,17)/
     2 10.,  1.297e-19,   45.8769,  563.8672,   -0.5434,    0.3563,
     3 11.,  2.464e-19,   42.2118,  434.1839,   -0.4015,    0.3364,
     4 12.,  4.312e-19,   39.0039,  345.8783,   -0.2749,    0.3371,
     5 14.,  1.097e-18,   33.7973,  238.2975,   -0.0441,    0.3855,
     6 17.,  3.227e-18,   28.2879,  159.6798,    0.3072,    0.5501,
     7 20.,  7.389e-18,   24.6485,  123.2140,    0.7220,    0.8318,
     8 25.,  2.087e-17,   21.1098,   98.3464,    1.6144,    1.7134,
     9 32.,  5.858e-17,   19.0148,   89.1510,    3.2276,    4.4387,
     x 40.,  1.358e-16,   17.6878,   72.8337,    4.5716,    9.0405,
     1 50.,  2.900e-16,   14.4089,   14.7174,    3.9266,    4.5716,
     2 63.,  5.795e-16,   13.7330,   -0.2729,    7.1026,    0.7890,
     3 80.,  1.034e-15,   18.3527,    3.4768,   16.8053,    7.8477,
     4 100.,  1.474e-15,    8.4112,    2.4827,   11.2036,    8.0189,
     5 120.,  1.663e-15,    6.9368,    1.8853,   11.0703,    8.1596,
     6 140.,  1.642e-15,    3.8003,    0.0828,    7.4776,    2.2102,
     7 170.,  1.442e-15,    0.5874,   -0.9489,    3.0230,   -2.6591,
     8 200.,  1.217e-15,   86.0001,   70.9908,   85.1300,  225.4768/
       data (zredeu(l,3),zrsigv(l,3),zafit(l,3),zbfit(l,3),
     1 zcfit(l,3),zdfit(l,3),l=18,30)/
     2 250.,  9.217e-16,   24.7032,   15.4806,   23.2680,   46.6313,
     3 320.,  6.645e-16,   11.2187,    7.2238,    9.4314,   17.7672,
     4 400.,  4.979e-16,    7.1710,    5.0844,    5.2534,   10.3083,
     5 500.,  3.815e-16,    5.0640,    3.6512,    3.1651,    6.5045,
     6 630.,  3.004e-16,    2.5568,    5.1021,    0.9028,    6.5380,
     7 800.,  2.468e-16,    0.7629,    8.6922,   -0.5642,    8.7898,
     8 1000.,  2.160e-16,    0.4407,    8.6023,   -0.6585,    8.0658,
     9 1200.,  2.006e-16,    0.3937,    5.8982,   -0.5536,    5.5717,
     x 1400.,  1.926e-16,   -0.1359,    5.9812,   -0.8952,    5.5364,
     1 1700.,  1.873e-16,   -0.8553,    9.0282,   -1.3623,    7.9472,
     2 2000.,  1.858e-16,   -0.8982,   13.6016,   -1.2218,   11.6160,
     3 2500.,  1.868e-16,    0.6688,   19.0625,    0.5069,   15.9416,
     4 3000.,  1.894e-16,    1.4983,   16.1376,    1.3781,   13.5456/
       data (zredeu(l,4),zrsigv(l,4),zafit(l,4),zbfit(l,4),
     1 zcfit(l,4),zdfit(l,4),l= 1,17)/
     2 40.,  4.330e-20,   43.6376,  478.1025,   -0.5753,    0.3583,
     3 45.,  9.324e-20,   39.2274,  348.6848,   -0.4425,    0.3588,
     4 50.,  1.784e-19,   35.5424,  267.9410,   -0.3265,    0.3862,
     5 63.,  6.620e-19,   28.6129,  163.6159,   -0.0190,    0.5413,
     6 80.,  2.223e-18,   23.3683,  116.0046,    0.5317,    0.9420,
     7 100.,  6.153e-18,   20.3171,  100.7269,    1.4872,    1.9749,
     8 120.,  1.318e-17,   19.1314,  100.6199,    2.7708,    4.1682,
     9 140.,  2.407e-17,   18.8868,  103.2965,    4.2005,    8.0426,
     x 170.,  4.888e-17,   18.2768,   84.7023,    5.3843,   14.2758,
     1 200.,  8.491e-17,   16.2296,   29.5378,    4.9693,    9.4740,
     2 250.,  1.692e-16,   17.1827,   -1.8880,    9.5490,    0.0000,
     3 320.,  3.145e-16,   35.9119,   22.7345,   34.7617,   37.5282,
     4 400.,  4.415e-16,   10.6303,    7.5539,   14.7738,   19.0596,
     5 500.,  4.723e-16,    6.1354,    2.9240,   10.6891,    9.2953,
     6 560.,  4.468e-16,    2.8817,   -0.2853,    6.4346,   -0.0434,
     7 630.,  4.040e-16,  251.6830,  454.4590,  251.4075, 1113.1820,
     8 700.,  3.609e-16,   61.0598,  100.9300,   60.4976,  231.2654/
       data (zredeu(l,4),zrsigv(l,4),zafit(l,4),zbfit(l,4),
     1 zcfit(l,4),zdfit(l,4),l=18,25)/
     2 800.,  3.079e-16,   29.0890,   45.6210,   27.7902,   93.3121,
     3 1000.,  2.354e-16,   14.2536,   23.2943,   12.3423,   38.1206,
     4 1200.,  1.928e-16,    7.7455,   18.7920,    5.8571,   24.8914,
     5 1400.,  1.669e-16,    4.8014,   19.2218,    3.0570,   21.8785,
     6 1700.,  1.445e-16,    3.4355,   18.8552,    1.9097,   18.9097,
     7 2000.,  1.325e-16,    2.6104,   18.6377,    1.3252,   17.3233,
     8 2500.,  1.233e-16,    1.9658,   36.0058,    1.1457,   30.4580,
     9 3000.,  1.203e-16,    5.9849,   61.3370,    5.4863,   49.3831/
c
c
c  test input
c
       if((kreac.lt.1).or.(kreac.gt.4)) go to 901
c
       zabeam=2.
       zaion=2.
       if(kreac.le.2) go to 50
       zabeam=2.
       zaion=3.
       if(pabeam.gt.2.5) zabeam=3.
       if(pabeam.gt.2.5) zaion=2.
50     continue
c
       zedeu=pebeam*(2./zabeam)
       if(zedeu.lt.zredeu(1,kreac)) go to 903
c
       iemax=iremax(kreac)
       if(zedeu.gt.zredeu(iemax,kreac)) go to 904
c
       if(ption.lt.0.) go to 905
       if(1.5*ption.gt.pebeam) go to 906
c
c  find zredeu entries which bracket zedeu
c
       in=inewtn/2
       ielowr=in
       iend=iremax(kreac)-(inewtn+1)/2
       do 110 je=in,iend
         if(zredeu(je,kreac).gt.zedeu) go to 111
         ielowr=je
110    continue
111    continue
c
c  find sigma-v at nearby values of ebeam by using a
c  rational polynomial fit to the temperature dependence; then
c  fill the first column of the divided difference array
c
       ir=kreac
       zx=(zabeam/zaion)*(ption/pebeam)
       do 220 jn=1,inewtn
         ie=ielowr-inewtn/2+jn
         zxfit(jn)=zredeu(ie,ir)
         zyfit(jn)=zrsigv(ie,ir)*(1.+zafit(ie,ir)*zx+
     1   zbfit(ie,ir)*zx**2)/(1.+zcfit(ie,ir)*zx+zdfit(ie,ir)*zx**2)
         zdfrnc(jn,1)=zyfit(jn)
220    continue
c
c  set up the rest of the divided difference array to apply
c  newton's divided difference formula for the collocation
c  polynomial of degree inewtn-1 through inewtn points
c
       do 240 jdiff=2,inewtn
         do 230 jn=1,inewtn-jdiff+1
           zdfrnc(jn,jdiff)=(zdfrnc(jn+1,jdiff-1)-zdfrnc(jn,jdiff-1))/
     1     (zxfit(jn+jdiff-1)-zxfit(jn))
230      continue
240    continue
c
c  evaluate the collocation polynomial
c
       zf=1.
       zy=0.
       do 280 jn=1,inewtn
         zy=zy+zf*zdfrnc(1,jn)
         zf=zf*(zedeu-zxfit(jn))
280    continue
       fusigv=zy
       return
c
c  error codes
c
c  kreac is out of bounds
901    continue
       kerr=1
       go to 910
c
c  pebeam is too small and out of bounds
903    continue
       kerr=3
       go to 910
c
c  pebeam is too large and out of bounds
904    continue
       kerr=4
       go to 910
c
c  ption is negative
905    continue
       kerr=5
       go to 910
c
c  ption is too large, pebeam is not 'fast'
906    continue
       kerr=6
       go to 910
c
c
910    continue
       fusigv=0.
       return
       end
c@hrecrd  .../baldur/code/bald/dbeams.f
c rgb 05-jun-96 include 'cfokkr.m'
c rgb 22-apr-96 copied contents of ../com/cfokkr.m comfok to write
c rgb 20.26 15-jan-92 removed common/combas/ and common/comddp/
c
         subroutine hrecrd(kledge,kcall,kret)
c
c 0.5  read and write restart records
c
      include 'cfokkr.m'
c
c---------------------------------------------------------------------
         go to (100,200),kcall
c---------------------------------------------------------------------
cl              1.         write record
c
  100    continue
         write(kledge)
     r  bbdd,   bbddrs, bbddsv, bbdt,   bbdtrs, bbdtsv,
     r  chscat, czgamm, cznue , halsps, hchexs, hdei  , hdeps , hdmu3 ,
     r  hdmub , hdmuc , hecmps, heelec, hei   , heinjs, heion , heloss,
     r  heplas, hesrci, hfi   , hlogbs, hmub  , hmuc  , hncmps, hninjs,
     r  hnloss, hnplas, hnsrcs, hscats, hsigv , hsin  , hslows, htbi  ,
     r  htlast, htspon, websps, wibsps,
     i  lhemax, lhemin, libeam, mhsp  , mxhe  , mxhmu , mxhsp , mxhsrc,
     i  nbbsp1, nbbsp2, nbbspd,
     i  nhbem1, nilast
c
         return
c
c---------------------------------------------------------------------
cl              2.         read record
c
  200    continue
         read(kledge,end=202,err=203)
     r  bbdd,   bbddrs, bbddsv, bbdt,   bbdtrs, bbdtsv,
     r  chscat, czgamm, cznue , halsps, hchexs, hdei  , hdeps , hdmu3 ,
     r  hdmub , hdmuc , hecmps, heelec, hei   , heinjs, heion , heloss,
     r  heplas, hesrci, hfi   , hlogbs, hmub  , hmuc  , hncmps, hninjs,
     r  hnloss, hnplas, hnsrcs, hscats, hsigv , hsin  , hslows, htbi  ,
     r  htlast, htspon, websps, wibsps,
     i  lhemax, lhemin, libeam, mhsp  , mxhe  , mxhmu , mxhsp , mxhsrc,
     i  nbbsp1, nbbsp2, nbbspd,
     i  nhbem1, nilast
c
c     success
  201    kret=1
         return
c
c     end of file
  202    kret=2
         return
c
c     error condition
  203    kret=3
         return
c
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@hclear   .../baldur/code/bald/dbeams.f
c pis 11-jun-98 changed 0.0 to 0 in reseti in two places
c rgb 22-apr-96 clear each array in common comfok
c
      subroutine hclear
c
c 1.2  clear variables and arrays
c
      include 'cfokkr.m'
c
c---------------------------------------------------------------------
c
c     block comfok
c
c
      bbdd   = 0.0
      bbddrs = 0.0
      bbddsv = 0.0
      bbdtsv = 0.0
      bbdt   = 0.0
      bbdtrs = 0.0
      chscat = 0.0
      czgamm = 0.0
      cznue  = 0.0
      halsps = 0.0
      hchexs = 0.0
      hdei   = 0.0
      hdeps  = 0.0
      hdmu3  = 0.0
      hdmub  = 0.0
      hdmuc  = 0.0
      hecmps = 0.0
      heelec = 0.0
      hei    = 0.0
      heinjs = 0.0
      heion  = 0.0
      heloss = 0.0
      heplas = 0.0
      hesrci = 0.0
      hfi    = 0.0
      hlogbs = 0.0
      hmub   = 0.0
      hmuc   = 0.0
      hncmps = 0.0
      hninjs = 0.0
      hnloss = 0.0
      hnplas = 0.0
      hnsrcs = 0.0
      hscats = 0.0
      hsigv  = 0.0
      hsin   = 0.0
      hslows = 0.0
      htbi   = 0.0
      htlast = 0.0
      htspon = 0.0
      websps = 0.0
      wibsps = 0.0
      lhemax = 0
      lhemin = 0
      libeam = 0
      mhsp   = 0
      mxhe   = 0
      mxhmu  = 0
      mxhsp  = 0
      mxhsrc = 0
      nbbsp1 = 0
      nbbsp2 = 0
      nbbspd = 0
      nhbem1 = 0
      nilast = 0
c
         return
         end
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@dposit   .../baldur/code/bald/dbeams.f
c  rgb 18-jul-96 set yvols(1) = 0.0 and use vols(iz+1,1)*zrndup
c  rgb 14:28 11-nov-93 changed call abort to call abortb
c  rgb 20.10 25-aug-91 make sure hbti starts no larger than tai
c  rgb 18.81 16-dec-90 vols(iz+1) changed to vols(iz+1,1)
c  rgb 18.41 22-jun-90 lbeams(32) .ge. 0 for uniform FREYA mesh
c     rewrote interpolation from BALDUR to FREYA mesh
c  dps 16-dec-87 splice in subroutine from 1-D code; make changes needed
c                to use eq'm. NOTE: use of shafs is only to allow code
c                  comparisons; definition of zdx2 following 
c                  "314 continue" is still in question.
c  rgb 12.60 25-aug-87 added code to determine HIBEAM from HPOWMW
c  drm 14-may-87 added call to shafs and multiply the shafranov
c               shift by hfutz(6).
c  rgb 05-may-85 changed xzoni(j)*rmins to ahalfs(j,2)
c                changed xbouni(j)*rmins to ahalfs(j,1)
c                implemented arrays yllipt, yvols in common /comloc/
c       drm 21-apr-87 this version calls the h(r) routines from survey
c               instead of the old freya routines.
c       drm 4-dec-84 use mikkelsen version
c       drm 10-jan-84 added bpinjs calculation
c       drm 22-feb-83 changed impurity array transfer to allow 4 imp
c       drm  24-feb-81  set nyprt = 5
c       drm 29-may-79 added yrhoim,yrhob,yqmean,mymp,mlyimp
c       amck 1-feb-78 yc() -> yr aft. 226
c       amck 31-jan-78 fix zdx2 = 0 in 238 loop while compressing
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 12-jan-78 cxxchx -> cxxcex
c       amck 12-jan-78 put in ch.ex. loss deposition, clamda dep. only
c               on ne (not ni)
c       amck 1-dec-77 move (1-yrlose) factor to hrb from hnsrcs
c       amck 1-dec-77 mtbeam -> mybeam in resetl
c       amck 1-dec-77 fix ilprof if libeam = 0
c       amck 30-nov-77 fix setting, testing of libeam() for libeam()=0
c       amck 23-sep-77 use rhoins+rhobis for ion density
c       amck 15-sep-77 set nyprt = 2, set yllipt, yphght, ypwdth, nyaper
c       amck 15-sep-77 set nyprt = 3
c       amck 15-sep-77 fix factor of 2 in hr normalization
c       amck 15-sep-77 cxx-- vars. are 1 / mean free path
c       amck 13-sep-77 check if libeam is all 0
c       amck 9-sep-77 fix subscript order on lybeam
c       amck 7-sep-77 add htlast = , nilast=
c       amck 7-sep-77 message each profile
c       amck 26-aug-77 hndeps/vol. is avg. source rate
c       amck 26-aug-77 iz -> iz1 after "do 318"
c       amck 26-aug-77 set libeam after 300
c       amck 25-aug-77 hrset -> hrset(1)
c       amck 25-aug-77 change nhpart -> nipart, set nyprt
c******************************************************************************
c
c
        subroutine dposit(kz)
c
c
cl      2.17    deposit high-energy ions
c       interface with freya
c
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'commhd.m'
c
c
c
        logical ilprof
c
        dimension
     r  zlksum(10),     zdnsum(10),     zhrsum(10), zshafs(mj)
     & ,ztemp(mj)
c
c------------------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /17/
c
c
        if ( nlomt2(isub) ) then
          write (6,*) ' *** 2.17 subroutine dposit bypassed'
          return
        endif
c                 
c
c------------------------------------------------------------------------------
c
c
c
        if (kz.eq.0) go to 200
        if (kz.gt.0) go to 1000
c
c..initialization
c               flag that no old profile exists
c
  100   continue
        htbi = epsinv
c
        do 104 jb = 1, mhbeam
          if (hton(jb).ge.epsinv) go to 104
          htbi = min(htbi,hton(jb)*ueit)
  104   continue
c
        htbi = max ( htbi, tai )
c
        yrmj = 0.0
        yrp = 0.0
        nilast = 0
        htlast = 0.0
c
        return
c
c
c
cl              compute profile (if necessary)
c
cl      2)      check if profile necessary
c
c
c
  200   continue
c
c  zrndup = 1.0 + a small number to rndeps real numbers up
c
        zrndup = 1.0 + 1.e-7
c
c               check if any beams have gone on or off
c
        ilprof = .false.
        ibeam = 0
c
      do 204 jb = 1, mhbeam
        if ( htbi .gt. htoff(jb)*ueit .or. tbi.le.hton(jb)*ueit) then
          if ( libeam(jb) .ne. 0 ) ilprof = .true.
        else
          ibeam = ibeam + 1
          if (libeam(jb).ne.ibeam) ilprof = .true.
        endif
  204   continue
c
        if (ibeam.le.0) ilprof = .false.
        if (ilprof) go to 254
c
c               if no beams, quit
c
      if ( ibeam .le. 0 ) then
        bpinjs=0.
        call resetr(hrb,mxhsp*6*11*myzone,0.0)
        call reseti(libeam,mxhbem,0)
        call resetr(hnsrcs,mxhsp*mxhsrc,0.0)
        call resetr(hdeps ,mxhsp*mxhsrc*mxhmu,0.0)
        return
      endif
c
c
c               check if we have old profile
c
        if (yrmj.le.0.0) go to 251
c
c               check if profile is too old
c
        if (nstep - nilast.ge.niprof) go to 252
c
c               check if compressed too much
c
        if ( yr(nyzone+1) .le. ahalfs(mzones-1,2)
     1  .or. yr(nyzone+1) .ge. ahalfs(mzones,2)
     2  .or.nzcomp.ne.0) go to 253
c
c
c               check if cross sections have changed much
c
c
      do 245 jspc=1,mhsp
c
        call resetr(zlksum,10,0.0)
        call resetr(zdnsum,10,0.0)
        zdtsum = 0.0
c
        iz = lcentr
        mnengy = min0(mxhsrc,mnengy)
c
        do 238 jc = 1, nyzone
          zten = 0.0
          znen = 0.0
          znin = 0.0
          zdx2 = 0.0
c
          if ( iz .gt. ledge ) go to 238
c
          do 224 jz = iz, ledge
            zten = zten + tes(2,jz)*rhoels(2,jz)*dx2i(jz)
            znen = znen + rhoels(2,jz)*dx2i(jz)
            znin = znin + (rhoins(2,jz)+rhobis(2,jz))*dx2i(jz)
            zdx2 = zdx2 + dx2i(jz)
            iz1 = jz
            if ( yr(jc+1)*zrndup .le. ahalfs(jz+1,2) ) go to 226
  224     continue
c
  226     continue
c
          iz = iz1
          if (yr(jc+1)*zrndup .gt. ahalfs(iz1,2) .and. iz1.lt.ledge)
     &       iz = iz1 + 1
          zdr = yr(jc+1) - yr(jc)                                       
c
        if (znen.gt.epslon) zten = zten / znen
        if (yrhoe(jc).gt.epslon) zdnen = znen / (zdx2*yrhoe(jc)) - 1.0
        if (yrhoi(jc).gt.epslon) zdnin = znin / (zdx2*yrhoi(jc)) - 1.0
        if (yte(jc).gt.epslon)
     1     zdtsum = zdtsum + abs((zten*evsinv / yte(jc) - 1.0)*zdr)
c
          do 230 je = 1, mnengy
            zdnsum(je) = zdnsum(je) + abs(clamda(je,jc,jspc)*zdnen*zdr)
            zlksum(je) = zlksum(je) + clamda(je,jc,jspc)*zdr
  230     continue
c
  238   continue
c
        if (zdtsum.ge.htchek*yr(nyzone+1)) go to 255
c
        do 244 je = 1, mnengy
          if (zlksum(je)*hnchek.lt.zdnsum(je)) go to 256
  244   continue
  245 continue
c
c               if compressing, normalization of hrb changes every t.s.
c
        if (nlcomp) go to 500
        return
c
c
cl              messages to tell why profile computed
c
c
  251   continue
        write (6,*) 'h(r) profile computed for the first time'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
  252   continue
        write (6,*) 'h(r) profile recomputed after niprof timesteps'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
  253   continue
        write (6,*) 'h(r) profile recomputed due to compression'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
  254   continue
        write (6,*) 'h(r) recomputed because beams turned on or off'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
  255   continue
        write (6,*) 'h(r) recomputed because T_e profile changed'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
  256   continue
        write (6,*) 'h(r) recomputed because cross-sections changed'
     &    ,' at timestep ',nstep,' time ',tai
        go to 300
c
c
c
c
cl      3)      get new profile
c
c
c
  300   continue
        ibeam = 0
        call reseti(libeam,mxhbem,0)
c
        do 304 jb = 1, mhbeam
        if (htbi.gt.htoff(jb)*ueit.or.tbi.le.hton(jb)*ueit) go to 304
        ibeam = ibeam + 1
        if (ibeam.gt.mybeam) go to 9030
        libeam(jb) = ibeam
  304   continue
        if (ibeam.le.0) go to 9031
c
c..set up a coarse mesh for the FREYA package
c
        nyzone = min0(min0( myzone-1 , nizone ), ledge-lcentr+1 )
c
c..old method of computing mesh (when lbeams(32) .lt. 0
c
      if ( lbeams(32) .lt. 0 ) then
c
        zfac = (float(ledge-lcentr+1) + 0.125) / float(nyzone)
        if (zfac.le.1.0) go to 9032
c
        iz1 = lcentr - 1
        yr(1) = ahalfs(2,1)
        yrmajr(1) = rmids(2,1)
        if ((hfutz(6).gt.0.).and.(leqtyp.eq.0)) then
          call shafs(zshafs)
          yrmajr(1)=yrmajr(1)+hfutz(6)*zshafs(lcentr)
        end if
        yvols(1) = vols(2,1)
        yllipt(1)=elong(2,1)
c
        do 318 jc = 1, nyzone
        ib = int(float(jc)*zfac) + lcentr
        iz0 = iz1 + 1
        iz1 = ib  - 1
c
        yr(jc+1) = ahalfs(ib,1)
        yrmajr(jc+1) = rmids(ib,1)
        if ((hfutz(6).gt.0.).and.(leqtyp.eq.0))
     1    yrmajr(jc+1)=yrmajr(jc+1)+hfutz(6)*zshafs(ib)
        yllipt(jc+1) = elong(ib,1)
        yvols (jc+1) = vols(ib,1)
c
        zne = 0.0
        zni = 0.0
        zee = 0.0
        znb=0.0
c
c       freya currently is dimensioned for only four impurities
        imimp=min0(mimp,4)
c
        do 308 jimp=1,4
        yrhoim(jimp,jc)=0.
        yqmean(jimp,jc)=0.
308     continue
c
        do 314 jz = iz0, iz1
        zne = zne + rhoels(2,jz) * dx2i(jz)
        zee = zee + tes(2,jz) * rhoels(2,jz) * dx2i(jz)
        znb=znb+rhobis(2,jz)*dx2i(jz)
        if(mhyd.lt.1) go to 311
        do 310 jhyd=1,mhyd
        zni=zni+rhohs(jhyd,2,jz)*dx2i(jz)
310     continue
311     continue
        if(mimp.lt.1) go to 313
        do 312 jimp=1,imimp
        yrhoim(jimp,jc)=yrhoim(jimp,jc)+rhois(jimp,2,jz)*dx2i(jz)
        yqmean(jimp,jc)=yqmean(jimp,jc)+cmean(jimp,2,jz)*dx2i(jz)
312     continue
313     continue
  314   continue
c
        zdx2 = 2.0*vols(mzones,1)/(vols(ib,1)-vols(iz0,1))
c
        yte(jc) = 0.0
        if (zne.gt.epslon) yte(jc) = zee * evsinv / zne
        yrhoe(jc) = zne * zdx2
        yrhoi(jc) = zni * zdx2
        yrhob(jc)=znb*zdx2
c
        if(mimp.lt.1) go to 317
        do 316 jimp=1,imimp
        yrhoim(jimp,jc)=yrhoim(jimp,jc)*zdx2
        yqmean(jimp,jc)=yqmean(jimp,jc)*zdx2
316     continue
317     continue
c
        yiota(jc) = 1.0 / q(iz0)
  318   continue
c
c..new method for computing uniform FREYA mesh (when lbeams .ge. 0)
c
c  make use of OMNILIB routines scopy, saxpy, sdot,...
c  and cubic spline utility routine cubint (see file DINTLIB in WSUTIL)
c
      else
c
c  freya currently is dimensioned for only four impurities
c
        imimp = min0 ( mimp , 4 )
c
        do 330 jc=1,nyzone+1
          yr(jc) = ahalfs(mzones,1) * (jc-1) / real(nyzone)
 330    continue
c
        call scopy (mzones,rmids(1,1),1,ztemp,1)
        if ((hfutz(6).gt.0.).and.(leqtyp.eq.0))
     &    call saxpy (mzones,hfutz(6),zshafs,1,ztemp,1)
c
        call cubint (ahalfs(1,1),ztemp,mzones,0,ceqtmp,kjbal
     &    ,yr,yrmajr,nyzone+1,0,0.0,1
     &    ,'interpolation of major radius in sbrtn dposit')
c
        call cubint (ahalfs(1,1),elong(1,1),mzones,0,ceqtmp,kjbal
     &    ,yr,yllipt,nyzone+1,0,0.0,1
     &    ,'interpolation of elongation in sbrtn dposit')
c
        call cubint (ahalfs(1,1),vols(1,1),mzones,0,ceqtmp,kjbal
     &    ,yr,yvols,nyzone+1,0,0.0,1
     &    ,'interpolation of volume in sbrtn dposit')
c
        yvols(1) = 0.0
c
        do 331 jz=1,mzones
          ztemp(jz) = 1 / q(jz)
 331    continue
c
        call cubint (ahalfs(1,1),ztemp,mzones,0,ceqtmp,kjbal
     &    ,yr,yiota,nyzone+1,0,0.0,1
     &    ,'interpolation of iota in sbrtn dposit')
c
c..volume weighted zone-centered thermodynamic variables
c
        iz = 2
c
        do 346 jc=1,nyzone
c
          zdvols = vols(iz,1) - yvols(jc)
c
          zne = rhoels(2,iz-1) * zdvols
          zee = tes(2,iz-1) * rhoels(2,iz-1) * zdvols
          znb = rhobis(2,iz-1) * zdvols
c
          if ( mhyd .gt. 0 ) then
            do 332 jhyd=1,mhyd
              zni = rhohs(jhyd,2,iz-1) * zdvols
 332        continue
          endif
c
          if ( mimp .gt. 0 ) then
            do 333 jimp=1,imimp
              yrhoim(jimp,jc) = rhois(jimp,2,iz-1) * zdvols
              yqmean(jimp,jc) = cmean(jimp,2,iz-1) * zdvols
 333        continue
          endif
c
 342      continue
          imore = 0
c
          if ( iz+1 .le. mzones
     &         .and. vols(iz+1,1)*zrndup .lt. yvols(jc+1) ) then
c
            zdvols = vols(iz+1,1) - vols(iz,1)
            imore  = 1
c
          else
c
            zdvols = yvols(jc+1) - vols(iz,1)
            imore  = 0
c
          endif
c
          zne = zne + rhoels(2,iz) * zdvols
          zee = zee + tes(2,iz) * rhoels(2,iz) * zdvols
          znb = znb + rhobis(2,iz) * zdvols
c
          if ( mhyd .gt. 0 ) then
            do 334 jhyd=1,mhyd
              zni = zni + rhohs(jhyd,2,iz) * zdvols
 334        continue
          endif
c
          if ( mimp .gt. 0 ) then
            do 335 jimp=1,imimp
              yrhoim(jimp,jc)=yrhoim(jimp,jc)+rhois(jimp,2,iz)*zdvols
              yqmean(jimp,jc)=yqmean(jimp,jc)+cmean(jimp,2,iz)*zdvols
 335        continue
          endif
c
          iz = iz + 1
c
          if ( imore .gt. 0 ) go to 342
c
c..temperature and density at FREYA zone centers
c
          yte(jc) = 0.0
          if ( zne .gt. epslon ) yte(jc) = zee * evsinv / zne
c
          znorm = 1. / ( yvols(jc+1) - yvols(jc) )
c
          yrhoe(jc) = zne * znorm
          yrhoe(jc) = zne * znorm
          yrhoi(jc) = zni * znorm
          yrhob(jc) = znb * znorm
c
          if ( mimp .gt. 0 ) then
            do 336 jimp=1,imimp
              yrhoim(jimp,jc) = yrhoim(jimp,jc) * znorm
              yqmean(jimp,jc) = yqmean(jimp,jc) * znorm
 336        continue
          endif
c
          if ( iz .gt. mzones ) go to 348
c
 346    continue
c
 348  continue
c
      endif
c
c..temporary printout
c
      if ( lbeams(32) .gt. 0 ) then
c
        write (nprint,390)
 390    format (/6x,'radius',6x,'yrmajr',6x,'yllipt',6x,'yvols'
     &    ,7x,'yiota')
c
      do 392 jc=1,nyzone+1
        write (nprint,391) yr(jc),yrmajr(jc),yllipt(jc),yvols(jc)
     &                    ,yiota(jc)
 391    format (1p5e12.4)
 392  continue
c
        write (nprint,394)
 394    format (/6x,'yte   ',6x,'yrhoi ',6x,'yrhob '
     &    ,6x,'yrhoim',6x,'yqmean')
c
      do 396 jc=1,nyzone+1
        write (nprint,395) yte(jc),yrhoi(jc),yrhob(jc)
     &     ,yrhoim(1,jc),yqmean(1,jc)
 395    format (1p5e12.4)
 396  continue
c
      endif
c
      if ( lbeams(32) .gt. 1 )
     &  call abortb(6,'temporary abort in sbrtn dposit')
c
c..end of computing FREYA mesh
c
c               other variables
c
        nyprt = 5
cbate        yllipt = ellipt
cbate        if (yllipt.le.0.0) yllipt = 1.0
        nypart = nipart
        nymu = nhmu
        mymp=mimp
        mlyimp=.false.
        if((extzef.ne.1.).and.(mimp.le.1)) mlyimp=.true.
        yrmj = rmids(mzones,1)
        yrp  = ahalfs(mzones,1)
        call resetl(nlybem,mybeam,.false.)
c
        do 328 jb = 1, mhbeam
        if (libeam(jb).le.0) go to 328
c
        ib = libeam(jb)
c
c  set the fast ion species index for each active beam.
c
        mxysp=mxhsp
        mysp=mhsp
        nyspec(ib)=0
        if ( nhbem1(1) .gt. 0 .and.  nhbem1(1) .le. mhbeam ) then
          if(nhbeam(jb).eq.nhbeam(nhbem1(1))) nyspec(ib)=1
        endif
        if ( nhbem1(2) .gt. 0 .and.  nhbem1(2) .le. mhbeam ) then
          if((mhsp.ge.2).and.(nhbeam(jb).eq.nhbeam(nhbem1(2))))
     &      nyspec(ib)=2
        endif
        if(nyspec(ib).le.0) then
        call ttyout('dposit: problem with the beam species')
        call exit(1)
        endif
c
c..normalize hfract
c
        z0 = 0.0
        do 322 je = 1, mxhfr
        z0 = z0 + hfract(je,jb)
  322   continue
c
        if (z0.gt.epslon) z0 = 1.0 / z0
        z1 = 0.0
        do 324 je = 1, mxhfr
        yfract(ib,je) = hfract(je,jb)*z0
        z1 = z1 + hfract(je,jb) * z0 / real(je)
  324   continue
c
c..determine hibeam [kAmp] from hpowmw [MW] and hebeam [keV]
c
        if ( hpowmw(jb) .gt. epslon .and. hebeam(jb)*z1 .gt. epslon )
     &          hibeam(jb) = hpowmw(jb) / ( hebeam(jb) * z1 )
c
        yabeam(ib) = habeam(jb)
        yzbeam(ib) = hzbeam(jb)
c
        nlybem(ib) = .true.
        yamp(ib) = hibeam(jb)
        yebeam(ib) = hebeam(jb)*uesh*evsinv
c
cbate end of determining hibeam from beam power
c
        yangl1(ib) = hangle(1,jb)
        yangl2(ib) = hangle(2,jb)
        yangl3(ib) = hangle(3,jb)
        yangl4(ib) = hangle(4,jb)
c
        nyshap(ib) = nhshap(jb)
        nyhght(ib) = nhprof(jb)
        nywdth(ib) = nhprfv(jb)
c
        ywdth(ib) = hwidth(jb) * uesl
        yhght(ib) = height(jb) * uesl
c
        yhzfoc(ib) = hfocl(jb) * uesl
        yvtfoc(ib) = hfoclv(jb) * uesl
        ydvghz(ib) = hdiv(jb)
        ydvgvt(ib) = hdivv(jb)
c
        yrmjp(ib)  = hrmaj(jb) * uesl
        yrpiv(ib)  = hrmin(jb) * uesl
        ylengt(ib) = hlenth(jb) * uesl
c
        yphght(ib) = haperv(jb) * uesl
        if (yphght(ib).le.0.0) yphght(ib) = 1.0e+10
        ypwdth(ib) = haper(jb) * uesl
        if (ypwdth(ib).le.0.0) ypwdth(ib) = 1.0e+10
        nyaper(ib) = nhaper(jb)
        if (nyaper(ib).le.0) nyaper(ib) = 1
  328   continue
c
        call copyr(hmub,1,ymub,1,nhmu+1)
        call copyr(hmuc,1,ymu ,1,nhmu  )
        nymu = nhmu
c
c
c
cl              get h(r) profile
c
        call gethfr
        nilast = nstep
        htlast = tai
c
        if (mnengy.gt.mxhsrc) go to 9033
c
c
c
c
cl      5)      normalize hrb on the baldur grid
c
c
c
  500   continue
c
c               compute integrated hrb over mu and volume
c
        do 550 jspc=1,mhsp
        ic0 = 1
        ic1 = 1
        call resetr(zhrsum,10,0.0)
c
        do 538 jz = lcentr, ledge
c
        if (ic0.ge.nyzone) go to 508
        zr = ahalfs(jz,2) * 2.0
        if (zr.le.(yr(ic1) + yr(ic1+1))) go to 508
        ic1 = ic1 + 1
        ic0 = ic1 - 1
        ic1 = min0( ic1,nyzone )
  508   continue
c
        zint1 = (zr - yr(ic0) - yr(ic0+1)) / (yr(ic1+1) - yr(ic0))
        zint0 = 1.0 - zint1
c
c
c
        do 518 je = 1, mnengy
        zhr = 0.0
        do 514 jmu = 1, nhmu
        zhr = zhr + (hrb(je,jmu,ic0,jspc)*zint0 +
     1          hrb(je,jmu,ic1,jspc)*zint1)*hdmuc(jmu)
  514   continue
        zhrsum(je) = zhrsum(je) + zhr*dx2i(jz)
  518   continue
c
c
  538   continue
c
        do 544 je = 1, mnengy
        if (zhrsum(je).gt.epslon)
     1          zhrsum(je) = 0.5*(1 - yrlose(je,jspc)) / zhrsum(je)
  544   continue
c
        do 548 jc = 1, nyzone
        do 548 jmu = 1, nymu
        do 548 je = 1, mnengy
        hrb(je,jmu,jc,jspc) = hrb(je,jmu,jc,jspc) * zhrsum(je)
  548   continue
550     continue
c
c
c
c
cl      6)      set hesrci, hnsrcs
c
c
c
  600   continue
c
        call resetr(hnsrcs,mxhsp*mxhsrc,0.0)
        call resetr(hesrci,mxhsp*mxhsrc,0.0)
c
        do 608 jb = 1, mhbeam
        if (libeam(jb).le.0) go to 608
        ib = libeam(jb)
        ispc=nyspec(ib)
c
        do 604 jfr = 1, mxhfr
c
        i001 = lybeam(ib,jfr)
        hnsrcs(i001,ispc) = hnsrcs(i001,ispc) +
     1          hibeam(jb)*yfract(ib,jfr)
  604   continue
c
  608   continue
c
c
        z0 = uesi * 10.0**(-fxes) / fces
        nhsrc = mnengy
c
        bpinjs=0.
        do 614 jspc=1,mhsp
        do 614 jsrc = 1, nhsrc
        hnsrcs(jsrc,jspc) = hnsrcs(jsrc,jspc) * z0
        hesrci(jsrc,jspc) = yener(jsrc) * evs * usie
        bpinjs=bpinjs+hnsrcs(jsrc,jspc)*yener(jsrc)*evs
  614   continue
c
        return
c
c
c
cl      10)     get hr(kz)
c
c
c
 1000   continue
c
        call resetr(hdeps,mxhsp*mxhsrc*mxhmu,0.0)
        if (yrmj.le.0.0) return
c
        zivols = 1. / vols(mzones,1)
        zr = ahalfs(kz,2) * 2.0
        ic0 = 1
        ic1 = 1
c
        do 1004 jc = 1, nyzone
        ic1 = jc
        if (zr.le.yr(jc) + yr(jc+1)) go to 1006
        ic0 = jc
 1004   continue
c
 1006   continue
c
        zint1 = (zr - yr(ic0) - yr(ic0+1)) / (yr(ic1+1) - yr(ic0))
        zint0 = 1.0 - zint1
        zchex = 0.0
c
        do 1016 jspc=1,mhsp
        do 1015 jsrc = 1, nhsrc
        if(hnsrcs(jsrc,jspc).eq.0.) go to 1015
        zfr0 = 0.0
        if (clamda(jsrc,ic0,jspc).gt.epslon)
     1    zfr0 = cxxcex(jsrc,ic0,jspc) / clamda(jsrc,ic0,jspc) * zint0
        zfr1 = 0.0
        if (clamda(jsrc,ic1,jspc).gt.epslon)
     1    zfr1 = cxxcex(jsrc,ic1,jspc) / clamda(jsrc,ic1,jspc) * zint1
c
        do 1014 jmu = 1, nhmu
        hdeps(jmu,jsrc,jspc) = hnsrcs(jsrc,jspc) * zivols *
     1   ( hrb(jsrc,jmu,ic0,jspc)*zint0 + hrb(jsrc,jmu,ic1,jspc)*zint1)
        zchex = zchex + hnsrcs(jsrc,jspc) * zivols * hdmuc(jmu) *
     1   (hrb(jsrc,jmu,ic0,jspc)*zfr0 + hrb(jsrc,jmu,ic1,jspc)*zfr1)
 1014   continue
 1015   continue
 1016   continue
c
c
c
c
        znh = rhoels(2,kz)
c
        if (znh.gt.epslon) znh = 1.0 / znh
c
        do 1028 jh = 1, mhyd
        shchxs(jh,kz) = rhohs(jh,2,kz) * znh * zchex
 1028   continue
c
c
        return
c
c
cl      90)     errors
c
c
c
 9030   continue
        call error_olymp(0,iclass,isub,3,
     1          'more than mtbeam (mybeam) beams on at once ')
        return
c
 9031   continue
        call error_olymp(0,iclass,isub,3,
     1          'profile being taken with no beams on ')
        return
c
c
c
 9032   continue
        call error_olymp(1,iclass,isub,3,
     1          'no. of injector zones .gt. no. of baldur zones ')
        call error_olymp(2,zfac,1,1,'  zfac  ')
        return
c
c
c
 9033   continue
        call error_olymp(1,iclass,isub,3,
     1          'too many beam energy groups for fokker-plank ')
        call error_olymp(3,mnengy,2,0,'energys ')
        call error_olymp(2,mxhsrc,2,0,'maximum ')
        return
c
        end
c
c
c******************************************************************************
c/ module shafs
c
c.              %%%%
c.              %
c.              hafs           %
c.              %
c.              %%%%
c
       subroutine shafs(pshafs)
c       mhr 15-jun-87 corrected integrations
c       mhr 26-may-87 added alpha pressure to shaf.shift
c       mhr 13-apr-87 shafranov shift for each zone jz
c       calculation of shafranov shift using h. howe's proctor formula
c       corresponding to lao, hirshman and wieland p.f.(1981)
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'cfokkr.m'
        include 'calpht.m'
c
        dimension pshafs(*)
        dimension zs1(mj), zfbprs(2,mj), zspres(2,mj), zintin(2,mj),
     1  zavprs(2,mj), znint(mj), z1(mj), z3(mj), zden(mj), zeperp(mj),
     2   zaps(2,mj)
c
        call resetr(pshafs,mzones,0.0)
c
        call resetr(zs1,mzones,0.0)
        call resetr(zfbprs,2*mzones,0.0)
        call resetr(zspres,2*mzones,0.0)
        call resetr(zintin,2*mzones,0.0)
        call resetr(zavprs,2*mzones,0.0)
        call resetr(zaps,2*mzones,0.0)
        call resetr(znint,mzones,0.0)
        call resetr(z1,mzones,0.0)
        call resetr(z3,mzones,0.0)
        call resetr(zden,mzones,0.0)
        call resetr(zeperp,mzones,0.0)
c
c       start of loop over jz
c
        do 30 jz = lcentr, mzones
c
c       compute fast ion pressure at zone centers
c
        do 10 jsp = 1,mhsp
        do 10 je = 1,nhe
        do 10 jmu=1,nhmu
        z1(jz) = z1(jz) +
     .  hfi(jmu,je,jz,jsp) * hei(je,jsp) *(hdmuc(jmu) - hdmu3(jmu))
        z3(jz) = z3(jz) + hfi(jmu,je,jz,jsp) * hdmuc(jmu)
  10    continue
c       zden is par/vol in each zone
        zden(jz) = z3(jz) * uied
c       zeperp is perp. energy/part in each zone
        zeperp(jz) = z1(jz) /(z3(jz)+epslon) * uise * useh
        zfbprs(2,jz) = zeperp(jz) * zden(jz) * uesh * uesd
c
        zfbprs(2,1) = zfbprs(2,2)
c       interpolate for zone boundaries
        zfbprs(1,jz+1) = .5 * (zfbprs(2,jz) + zfbprs(2,jz+1))
c
c       compute alpha pressure at zone centers
        zaps(2,jz) = alphai(jz) * ealfai(jz) * .6666666666 *uisd*uise
c
c       interpolate alpha pressure for zone boundaries
        zaps(1,jz+1)=.5 * (zaps(2,jz+1) + zaps(2,jz))
        zaps(1,2) = zaps(2,2)
c
c       compute pressure at zone boundaries and centers
        zspres(2,jz) = rhoels(2,jz)*tes(2,jz) + rhoins(2,jz)*tis(2,jz)
     .  +zfbprs(2,jz) + zaps(2,jz)
c
        zspres(1,jz+1) = rhoels(1,jz+1)*tes(1,jz+1) +
     .  rhoins(1,jz+1)*tis(1,jz+1)
     .  +zfbprs(1,jz+1) + zaps(1,jz+1)
c
c       compute internal inductance at zone  boundaries
        zintin(1,jz+1) = zintin(1,jz) +4. * (ellipt**2 +1.) /
     .  ((3. * ellipt**2 +1.))
     .  *dx2i(jz) * bpols(2,jz)**2
c
c       compute average pressure at zone boundaries
        zavprs(1,jz+1) = zavprs(1,jz) +
     .   2.*dx2i(jz)*zspres(2,jz)
c
c
c
 30     continue
        do 35 jz=lcentr,ledge
c
c       internal inductance and average pressure at zone boundaries
c
        zintin(1,jz+1)=zintin(1,jz+1)/(bpols(1,jz+1)*xbouni(jz+1))**2
        zavprs(1,jz+1)=zavprs(1,jz+1)/xbouni(jz+1)**2
c
c       interpolate for zone centers
c
        zintin(2,jz) = .5*(zintin(1,jz) + zintin(1,jz+1))
        zavprs(2,jz) = .5*(zavprs(1,jz) + zavprs(1,jz+1))
        zintin(2,2)=zintin(1,3)
        zavprs(2,2)=zavprs(1,3)
c
c
        zs1(jz+1) = zs1(jz) + rmins**2 *((4. * ellipt**2 /
     .  (3.* ellipt**2 + 1.))
     .  *dx2i(jz) * (zavprs(2,jz) -
     .   zspres(2,jz)) * 8.* fcpi / bpols(2,jz)**2
     .  + dx2i(jz) * .5 * (zintin(2,jz) +
     .   (ellipt**2 -1.) / (3. * ellipt**2 + 1.)))
c
  35    continue
c       compute shafranov shift at center
        znint(lcentr) = zs1(mzones)
        pshafs(lcentr) = rmajs*.5 * (sqrt(1.+
     .  4.* znint(lcentr)/rmajs**2) - 1.)
c
c       compute shafranov shift for other flux surfaces
        do 40 jz = 3,mzones
        znint(jz) = zs1(mzones) -zs1(jz)
        pshafs(jz) =  znint(jz)/(rmajs + pshafs(lcentr))
  40    continue
        return
        end
c
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c
c     monte carlo freya package for neutrals
c
c..the freya package was completely revised by dave mikkelsen april 1987
c
c**********************************************************************c
c
c/ module gethfr
c  rgb 19-oct-90 18.74  changed paramter lh=21 to lh=51
c  rgb 18.40 10:00 22-jun-90 increase ixhofr from 1 to 2
c  rgb 02-jun-89 replaced myimp (not defined) with mymp
c 17-dec-87 dps fix bug is zdvol def'n. (660 loop).
c 16-dec-87 dps splice in 1-D version with multiple beam species.
c 14-may-87 drm changed the call resetr of ziener,zfbeam to use mybeam.
c
c.              %%%%
c.              %
c.               0ethfr          %
c.              %
c.              %%%%
c
       subroutine gethfr
       external case8
c
       include 'cfreya.m'
c
c  parameter index:
c  lh      hofr number of height points in beam cross section.
c  lr      hofr number of rtang points in beam cross section.
c  lz      hofr number of radial zones.
c  lmu     hofr number of cosine of pitch-angle bins.
c  changes here should be duplicated in the parameter statement in seghfr.
       parameter (lh=51, lr=51, lz=26, lmu=11)
       data mhtdim, mrtdim, mzndim, mudim /lh, lr, lz, lmu/
c
       dimension zdvol(lz), zdmuc(lmu)
!cap
       real, allocatable :: ziener(:,:), zfbeam(:,:,:)
       dimension height(lh),artang(lr),
     1 rhrmaj(2,lz,lh),weight(lh,lr),anesig(lz,3),
     2 nhstop(lh,lr),minzon(lh),hrdep(lz,lmu,3),
     3 segmnt(lz,2,lh,lr),nznone(lh,lr),nzntwo(2,lh,lr)
c
c
       data evs,fcpi,fcmp,fxnucl/1.602e-12,3.14159,1.673,-24./
c
!cap
       allocate(ziener(mybeam,mxysp), zfbeam(mybeam,3,mxysp))
c
       if(nymu.gt.mudim) then
         call ttyout('gethfr: too many pitch angle bins for the arrays')
         return
       endif
       icentr=1
       iedge=nyzone
       call nhrset
c  set up zfbeam(jb,jf,jspc): the fraction of the total current for
c  energy index lybeam(jb,jf) and species jspc which is contributed
c  by component jf of beam jb.
!{cap
       ziener = 0.
       zfbeam = 0.
!}
       do 120 jb=1,nybeam
         if(.not.nlybem(jb)) go to 120
        ispc=nyspec(jb)
c  normalize the yfracts.
         znormi=1./(yfract(jb,1)+yfract(jb,2)+yfract(jb,3))
         yfract(jb,1)=znormi*yfract(jb,1)
         yfract(jb,2)=znormi*yfract(jb,2)
         yfract(jb,3)=znormi*yfract(jb,3)
c  add up the contributions.
         do 110 jf=1,3
           ie=lybeam(jb,jf)
           ziener(ie,ispc)=ziener(ie,ispc)+yamp(jb)*yfract(jb,jf)
110      continue
120    continue
c  calculate zfbeam.
       do 140 jb=1,nybeam
         if(.not.nlybem(jb)) go to 140
        ispc=nyspec(jb)
         do 130 jf=1,3
           ie=lybeam(jb,jf)
           if(ziener(ie,ispc).gt.0.) then
             zfbeam(jb,jf,ispc)=yamp(jb)*yfract(jb,jf)/ziener(ie,ispc)
           else
             zfbeam(jb,jf,ispc)=0.
           endif
130      continue
140    continue
c  loop over the beams
       do 800 jb=1,nybeam
         if(.not.nlybem(jb)) go to 800
        ispc=nyspec(jb)
c
c..compute where the chords intersect flux surfaces
c
c  ixhofr = number of chords in each BALDUR zone height
c
c  sbrtn setrmj returns the last four arguments
c
         ixhofr = 2
         call setrmj(ixhofr,mhtdim,mzndim,icentr,iedge,
     1   yr,yrmajr,yllipt,     nhtmax,minzon,height,rhrmaj)
c
c  calculate the distance from the ion source to the minimum
c    minor radius point along the beam path.
c
         zrtang=(yrmjp(jb)+yrpiv(jb))*sin((fcpi/180.)*yangl3(jb))
         if(abs(zrtang).lt.yrmajr(icentr)) then
           zl=sqrt(yrmajr(icentr)**2-zrtang**2)
         else
           zl=0.
         endif
         ztlnth=ylengt(jb)+sqrt((yrmjp(jb)+yrpiv(jb))**2-zrtang**2)-zl
c        convert dvgzs from degrees to distances in the middle of the plasma.
         zvdvgz=ztlnth*ydvgvt(jb)*fcpi/180.
         zhdvgz=ztlnth*ydvghz(jb)*fcpi/180.
         call setwgt(mhtdim,mrtdim,ywdth(jb),yhght(jb),
     1   yphght(jb),ypwdth(jb),zvdvgz,zhdvgz,yvtfoc(jb),yhzfoc(jb),
     2   ztlnth,zrtang,height,     nhtmax,nrtmax,artang,weight)
c
c
c        set sigvbar in all zones for all charge states
c
         izfast=int(yzbeam(jb))
         do 650 jf=1,3
           iq=min0(jf,izfast)
           zef=yebeam(jb)
           if(izfast.eq.1) zef=yebeam(jb)/float(jf)
           zv=sqrt(2.*evs*zef/(fcmp*(10.**fxnucl)*yabeam(jb)))
c
c          start do loop over zones
c
           do 640 jc=1,nyzone
             zei=aioniz(izfast,iq,yte(jc),yabeam(jb),3)
             anesig(jc,jf)=zei*yrhoe(jc)/zv
c
c            loop over impurity species
c
             if(mymp.gt.0) then
               do 630 jimp=1,mymp
                 zzi=yqmean(jimp,jc)
                 zefq=zef/zzi
                 ipel=5
                 if(zefq.gt.yabeam(jb)*2.e5) ipel=1
                 zelimp=zzi*aioniz(izfast,iq,zefq,yabeam(jb),ipel)
                 anesig(jc,jf)=anesig(jc,jf)+zelimp*yrhoim(jimp,jc)
 630           continue
             endif
c            hydrogenic ionization and charge exchange.
c       always use the gryzinski proton ionization formula:
           ipi=1
             zpi=aioniz(izfast,iq,zef,yabeam(jb),ipi)
             zpcx=0.
             if(izfast.eq.1) zpcx=aioniz(izfast,iq,zef,yabeam(jb),4)
             anesig(jc,jf)=anesig(jc,jf)+(zpi+zpcx)*yrhoi(jc)
c            save cross sections calculating the charge exchange fraction.
             ie=lybeam(jb,jf)
             cxxcex(ie,jc,ispc)=zpcx*yrhoi(jc)
             clamda(ie,jc,ispc)=anesig(jc,jf)
640        continue
650      continue
c        prepare to call seghfr.
         do 660 jc=1,nyzone
           zdvol(jc) = (yvols(jc+1) - yvols(jc)) / yvols(nyzone+1)
660      continue
c        calculate the width of the pitch-angle bins.
         do 670 jmu=1,nymu
           zdmuc(jmu)=ymub(jmu+1)-ymub(jmu)
670      continue
         zrcwal=yrmjp(jb)
         zrdwal=yrpiv(jb)
         zffull=yfract(jb,1)
         zfhalf=yfract(jb,2)                                      
         call seghfr(mhtdim,mrtdim,mzndim,mudim,nymu,nhtmax,nrtmax,
     1   icentr,iedge,zrcwal,zrdwal,minzon,height,artang,rhrmaj,
     2   nznone,nzntwo,nhstop,
     3   zffull,zfhalf,anesig,weight,zdvol,ymub,zdmuc,    hrdep)
c        put the h(r) into hrb
         do 730 jf=1,3
           ie=lybeam(jb,jf)
           do 720 jmu=1,nymu
             do 710 jc=1,nyzone
               hrb(ie,jmu,jc,ispc)=hrb(ie,jmu,jc,ispc)+
     1         zfbeam(jb,jf,ispc)*hrdep(jc,jmu,jf)
710          continue
720        continue
730      continue
800    continue
c  set yrlose after integrating hrb.
        do 860 jspc=1,mysp
       do 850 je=1,mnengy
         zsum=0.
         do 840 jmu=1,nymu
           do 830 jc=1,nyzone
             zsum=zsum+hrb(je,jmu,jc,jspc)*zdvol(jc)*zdmuc(jmu)
830        continue
840      continue
         yrlose(je,jspc)=1.-zsum
850    continue
860    continue
!cap
       deallocate(ziener, zfbeam)
c
       return
       end
c******************
      BLOCK DATA case8
       include 'cfreya.m'
       data mybeam, myzone, nyzone/ 12, 26, 20/
       end
c-------------------------------------
c/ module nhrset
c
c.              %%%%
c.              %
c.              hrset          %
c.              %
c.              %%%%
       subroutine nhrset
c
       include 'cfreya.m'
c  index of name changes
c  nzones  nyzone
c  r2      yr
c  r2sq    yrsq
c  rp      yrp
c  rmj     yrmj
c  mtbeam  mybeam
c  rpivot  yrpiv
c  rmjpvt  yrmjp
c  nmu     nymu
c  ntbeam  nybeam
c  nlbeam  nlybem
c  hener   yener
c  ebeam   yebeam
c  lbeam   lybeam
cc-----------------------------------------------------------
c
c
cl      2.1 input parameters
c
       do 210 j=1,nyzone+1
         yrsq(j)=yr(j)*yr(j)
210    continue
c
c  2.7     neutral injector parameters
       do 2720 jb=1,mybeam
         if(.not.nlybem(jb)) go to 2720
         if(yrpiv(jb).le.0.) yrpiv(jb)=yrp
         if(yrmjp(jb).le.0.) yrmjp(jb)=yrmj
 2720  continue
c
c  2.7     neutral injector parameters
c
c-----------------------------------------------------------------------
       nybeam=0
       do 2710 j=1,mybeam
         if(nlybem(j)) nybeam=nybeam+1
2710   continue
c
cl                 3.   numerical scheme
c
c
c  energy boundaries set and ordered
c
c
cl      energy array is first filled
c
       jq=nybeam
       do 3312 j=1,mybeam
         if(.not.nlybem(j)) go to 3312
         jx=3*nybeam-(jq-1)
         yener(jx)=yebeam(j)
         jx=jx-nybeam
         yener(jx)=yebeam(j)/2.
         jx=jx-nybeam
         yener(jx)=yebeam(j)/3.
         jq=jq-1
3312   continue
c
cl      energy array is then arranged in order from low to high energies
c
       ibeam=3*nybeam-1
       if(ibeam.le.0) go to 3350
3320   continue
       ichang=0
       do 3325 j=1,ibeam
         if(yener(j).le.yener(j+1)) go to 3324
         zx=yener(j)
         yener(j)=yener(j+1)
         yener(j+1)=zx
         ichang=1
3324     continue
3325   continue
       if (ichang.ge.1) go to 3320
3350   continue
c
cl      now we throw out duplicate energies
c
       zeps1=.999
       zeps2=1.001
       ix=3*nybeam
       iq=1
3352   continue
       if(.not.(yener(iq).gt.zeps1*yener(iq+1).and.yener(iq).lt.
     1 zeps2*yener(iq+1))) go to 3360
       ix=ix-1
       if(iq.ge.ix) go to 3380
       do 3355 j=iq,ix
         yener(j)=yener(j+1)
3355   continue
       go to 3352
3360   continue
       iq=iq+1
       if(iq.ge.ix) go to 3380
       go to 3352
3380   continue
       mnengy=ix
c
cl      now lybeam array is filled
c  lybeam(ib,ie) is the energy corresponding to beam ib, energy
c  component ie
c
       z3=1./3.
       zeps1=.999
       zeps2=1.001
       do 3391 je=1,mnengy
         do 3390 j=1,mybeam
           if(.not.nlybem(j)) go to 3390
           if((yebeam(j).gt.zeps1*yener(je)).and.
     1     (yebeam(j).lt.zeps2*yener(je))) lybeam(j,1)=je
           zx=yebeam(j)*.5
           if((zx.gt.zeps1*yener(je)).and.(zx.lt.zeps2*yener(je)))
     1     lybeam(j,2)=je
           zx=yebeam(j)*z3
           if((zx.gt.(zeps1*yener(je))).and.(zx.lt.(zeps2*yener(je))))
     1     lybeam(j,3)=je
3390     continue
3391   continue
c
cl      zero out hrb(j1,j2,j3)
c
        do 3400 jspc=1,mysp
       do 3400 j1=1,mnengy
         do 3400 j2=1,nymu
           do 3400 j3=1,nyzone
             hrb(j1,j2,j3,jspc)=0.
3400   continue
       return
       end
c/ module setrmj
c
c.              %%%%
c.              %
c.              etrmj          %
c.              %
c.              %%%%
c
       subroutine setrmj(ksegzn,khtdim,kzndim,kcentr,kedge,
     1 prminr,prmajr,pellip,     khtmax,kminzn,phight,phrmaj)
c
c  setrmj interfaces with the plasma geometry to provide phight,
c  the heights above the midplane for the beam rays, and rhrmaj,
c  the major radius at the intersections of the rays with the
c  plasma zone boundaries. kminzn is the innermost plasma zone
c  index which a ray of height phight(jh) can penetrate.
c
c@setrmj  11040/bald90/wbaldn1  file: DBEAMS
c  rgb 29-jan-89  Increased dimension of pellip and prmajr to kedge+1
c     removed vectorization from do loops 105 and 110
c     otherwise, khtmax is incorrectly computed with value 1 too large
c     changed pellip(jz+1) to pellip(jz) in zdelr after do 125
c  dps 16-dec-87  NOTE: treatment of elongated flux surfaces here is now
c                 fundamentally different than in 1-D code.
c..rgb 05-may-87  For now, the flux surfaces are taken to be
c                 shifted ellipses.
c  pellip has been changed from a scalar to an array.
c
       dimension prminr(kedge+1),prmajr(kedge+1),phight(khtdim),
     & kminzn(khtdim),phrmaj(2,kzndim,khtdim), pellip(kedge+1)
c
c  index of variables used in setrmj  (linear dimensions in cm.)
c
c..input arguments:
c  ksegzn   number of phight elements per plasma zone.
c  khtdim   length of jh dimension in dimension statements.
c  kzndim   length of jz dimension in dimension statements.
c  kcentr   index of central plasma zone in prminr and prmajr.
c  kedge    index of edge plasma zone in prminr and prmajr.
c  prminr(jz)  the height above the midplane of the inner boundary
c  of plasma zone jz (defined for jz=kcentr,kedge+1).
c  prmajr(jz)  the average major radius of the inner boundary
c  of plasma zone jz (defined for jz=kcentr,kedge).
c  pellip(jz)  the ellipticity of magnetic surfaces.
c
c..output arguments:
c  khtmax   number of phight elements used.
c  kminzn(jh)  index of most central plasma zone ray jh penetrates.
c  phight(jh)  height for jh rays (jh=1,khtmax).
c  phrmaj(ib,jz,jh)  the min. (jb=1) and max. (jb=2) major radius of
c  the outer surface of plasma zone jz at phight(jh).
c
       data imess/ 0/
c
c..test for enough room in the plasma zone dimension of rhrmaj.
c
       if(kedge.gt.kzndim) then
         if(imess.lt.100) then
           imess=imess+1
       callttyout(' setseg: too many plasma zones, arrays too small')
         endif
         return
       endif
c
c  set phight array according to the plasma zoning.
c  beam is assumed to be centered around the equatorial plane.
c
c dir$ novector
c
       ih=0
       do 110 jz=kcentr,kedge
c
c        set ksegzn ray heights spaced evenly between plasma zone heights
c        specified by prminr.
c
         zdr=(prminr(jz+1)-prminr(jz))/float(ksegzn)
         do 105 jseg=1,ksegzn
           zh = pellip(jz) * (prminr(jz)+zdr*(float(jseg)-0.5))
           ih=ih+1
           if(ih.le.khtdim) then
             phight(ih)=zh
             kminzn(ih)=jz
             khtmax=ih
           endif
105      continue
110    continue
c
c  set phrmaj array to be the two major radii of the
c  the outer boundary of each zone at a height = phight(jh).
c
       do 130 jh=1,khtmax
         do 115 jz=kcentr,kminzn(jh)-1
           phrmaj(1,jz,jh)=-1.
           phrmaj(2,jz,jh)=-1.
115      continue
         do 125 jz=kminzn(jh),kedge
           zdelr=sqrt(prminr(jz+1)**2-(phight(jh)/pellip(jz))**2)
           phrmaj(1,jz,jh)=prmajr(jz+1)-zdelr
           phrmaj(2,jz,jh)=prmajr(jz+1)+zdelr
125      continue
130    continue
c
c dir$ vector
c
       return
       end
c@setwgt  11040/bald90/wbaldn1  file: DBEAMS
c 18.82 11-jan-91 abs(prtaim) to prtaim to correct beam driven current
c
c.              %%%%
c.              %
c.              etwgt          %
c.              %
c.              %%%%
c
       subroutine setwgt(khtdim,krtdim,pswdth,pshght,
     1 pphght,ppwdth,pvdvgz,phdvgz,pvfocl,phfocl,ptlnth,
     2 prtaim,phight,     khtmax,krtmax,prtang,pweigt)
c
c  uses beam source parameters to set up the arrays which determine
c  the beam intensity across its cross section.
c
       dimension phight(khtdim),prtang(krtdim),pweigt(khtdim,krtdim)
c
c  index of variables used in setwgt  (linear dimensions in cm.)
c  arguments passed on to sntgrl but not used in setwgt:
c  pvfocl   vertical focal length.
c  phfocl   horizontal focal length.
c  ptlnth   distance from ion source to the minimum minor radius
c  along the beam path.
c  input arguments:
c  khtdim   length of jh dimension in dimension statements.
c  krtdim   length of jr dimension in dimension statements.
c  khtmax   number of phight elements used.
c  pswdth   neutral beam source width.
c  pshght   neutral beam source height(assumed to be centered on midplane).
c  pphght   beam port height (full height).
c  ppwdth   beam port width.
c  pvdvgz   vertical divergence, gaussian: exp(-x**2/(2.*pvdvgz**2))
c  phdvgz   horizontal divergence, gaussian: exp(-x**2/(2.*phdvgz**2))
c  prtaim   neutral beam tangency radius (minimum major radius
c  along the centerline of the beam path).
c  phight(jh)  height for jh rays (jh=1,khtmax).
c  output arguments:
c  khtmax   number of phight elements used.
c  krtmax   number of prtang elements used.
c  prtang(jr)  tangency radius for jr rays (jr=1,krtmax).
c  pweigt(jh,jr)  weighting factor of ray jh,jr (sum normalized to 1.)
c
       data imess/ 0/
c
       ihtmax=khtmax
       zhmax=min(0.5*pphght,0.5*pshght+3.*pvdvgz)
       zrmax=min(ppwdth,pswdth+6.*phdvgz)
c
       do 100 jh=1,ihtmax
c        pshght is the full height of the neutral beam source, half
c        is above the midplane and half is below.
         if(phight(jh).le.zhmax) khtmax=jh
100    continue
       if(zhmax.gt.phight(ihtmax)) then
         if(imess.lt.100) then
           imess=imess+1
       callttyout(' setwgt: truncated beam height, khtdim too small')
         endif
       endif
c  set prtang array to have the same average zone spacing
c  as the phight array
       krtmax=max(1,int((zrmax/zhmax)*khtmax + 0.5))
       if(krtmax.gt.krtdim) then
         if(imess.lt.100) then
           imess=imess+1
           callttyout('   setwgt: wider hor. spacing; krtdim too small')
         endif
         krtmax=krtdim
       endif
       zsum=0.
       zdr=zrmax/float(krtmax)
       do 120 jr=1,krtmax
         prtang(jr)=prtaim-zrmax*0.5+zdr*(float(jr)-0.5)
         do 115 jh=1,khtmax
           if(pvdvgz+phdvgz.le.0.) then
c            set pweigt for uniform intensity across the face of the source.
             pweigt(jh,jr)=1./float(khtmax*krtmax)
           else
             zh=phight(jh)
             zr=prtang(jr)-prtaim
c            set pweigt for gaussian intensity distribution due to divergence.
             pweigt(jh,jr)=sntgrl(zr,zh,pswdth,pshght,ptlnth,
     1       phdvgz,pvdvgz,phfocl,pvfocl)
           endif
           zsum=zsum+pweigt(jh,jr)
115      continue
120    continue
c
       if(abs(zsum-1.).gt.0.001) then
c        renormalize pweigt to unit sum over all elements.
         znorm=1./zsum
         do 220 jh=1,khtmax
           do 210 jr=1,krtmax
             pweigt(jh,jr)=znorm*pweigt(jh,jr)
210        continue
220      continue
       endif
c
       return
       end
c/ module sntgrl
c
c.              %%%%
c.              %
c.              ntgrl          %
c.              %
c.              %%%%
c
       function sntgrl(pdr,pdh,pswdth,pshght,ptlnth,
     1 phdvgz,pvdvgz,phfocl,pvfocl)
c
c  evaluates a gaussian integral over the face of the ion source
c  to find the relative intensity of the beam ray specified by
c  pdr and pdh.
c
c  index of variables used in sntgrl  (linear dimensions in cm.)
c  input arguments:
c  pphght   beam port height (full height).
c  ppwdth   beam port width.
c  pvdvgz   vertical divergence, gaussian: exp(-x**2/(2.*pvdvgz**2))
c  phdvgz   horizontal divergence, gaussian: exp(-x**2/(2.*phdvgz**2))
c  pvfocl   vertical focal length.
c  phfocl   horizontal focal length.
c  ptlnth   distance from ion source to the minimum minor radius
c  along the beam path.
c
       data zfcpi/ 3.14159/
c
       ihpnts=max(5,int(2.*pshght/pvdvgz+0.5))
       irpnts=max(5,int(2.*pswdth/phdvgz+0.5))
c  loop over the face of the ion source.
       zsum=0.
       do 180 jr=1,irpnts
         zsr=pswdth*(-0.5 + (float(jr)-0.5)/float(irpnts))
c        advance the point on the source to its focussed position
c        in the plasma at a distance ptlnth from the source.
         zpr=zsr*(1.-ptlnth/phfocl)
         do 170 jh=1,ihpnts
           zsh=pshght*(-0.5 + (float(jh)-0.5)/float(ihpnts))
c          advance the point on the source to its focussed position
c          in the plasma at a distance ptlnth from the source.
           zph=zsh*(1.-ptlnth/pvfocl)
c          calculate the distance of the input point pdh,pdr from the
c          point zph,zpr and divide by the relevant divergence distances.
           zexarg=0.5*(((pdh-zph)/pvdvgz)**2 + ((pdr-zpr)/phdvgz)**2)
           zsum=zsum+exp(-zexarg)
170      continue
180    continue
c  normalize such that in the center of a large source the value
c  of sntgrl will be 1.
       sntgrl=zsum*(pshght*pswdth/(2.*zfcpi*phdvgz*pvdvgz))
     1 /float(ihpnts*irpnts)
       return
       end
c/ module seghfr
c  rgb 19-oct-90 18.74  changed paramter lh=21 to lh=51
c   dps 16-dec-87 replace prtang with zabsrt in if-then following
c                 "140 continue" - as in 1-D code.
c
c.              %%%%
c.              %
c.              eghfr          %
c.              %
c.              %%%%
c
       subroutine seghfr(khtdim,krtdim,kzndim,kmudim,kmu,khtmax,krtmax,
     1 kcentr,kedge,prmjwl,prmnwl,kminzn,phight,prtang,phrmaj,
     2 kznone,kzntwo,khstop,
     3 pffull,pfhalf,pnesig,pweigt,pdvol,pmub,pdmuc,    phrdep)
c
c  seghfr tracks each ray through the plasma to determine which
c  plasma zones are intersected and the path length in each
c  intersected zone. it also determines which rays strike the
c  vacuum vessel wall (or any toroidal surface defined by prmjwl
c  and prmnwl. sets phrdep(jz,jmu,jfrac) using pnesig and zlenth.
c
c  parameter index:
c  lh      hofr number of height points in beam cross section.
c  lr      hofr number of rtang points in beam cross section.
c  lz      hofr number of radial zones.
c  changes here should be made in the parameter statement in gethfr.
       parameter (lh=51, lr=51, lz=26)
       dimension zexpi(lz),zrsrc(lz),zlenth(lz,2),imu(lz,2)
c
       dimension phrmaj(2,kzndim,khtdim),prtang(krtdim),phight(khtdim),
     1 kminzn(khtdim),kznone(khtdim,krtdim),kzntwo(2,khtdim,krtdim),
     2 khstop(khtdim,krtdim)
c
       dimension pweigt(khtdim,krtdim),pnesig(kzndim,3),
     1 phrdep(kzndim,kmudim,3),pdvol(kedge),pmub(kmudim),
     2 pdmuc(kmudim)
c
       data imess/ 0/
c
c  index of variables used in seghfr  (linear dimensions in cm.)
c  input arguments:
c  khtdim   length of jh dimension in dimension statements.
c  krtdim   length of jr dimension in dimension statements.
c  kzndim   length of jz dimension in dimension statements.
c  kmudim   length of jmu dimension in dimension statements.
c  kmu      number of pmub and pdmuc elements used.
c  khtmax   number of phight elements used.
c  krtmax   number of prtang elements used.
c  kcentr   index of central plasma zone.
c  kedge    index of edge plasma zone.
c  prmnwl   minor radius of the toroidal surface (vacuum vessel)
c  which may stop the beam from passing through the
c  plasma twice.
c  prmjwl   major radius of the toroidal surface (vacuum vessel)
c  which may stop the beam from passing through the
c  plasma twice.
c  kminzn(jh)  index of most central plasma zone ray jh penetrates.
c  phight(jh)  height for jh rays (jh=1,khtmax).
c  prtang(jr)  tangency radius for jr rays (jr=1,krtmax).
c  phrmaj(ib,jz,jh)  the min. (jb=1) and max. (jb=2) major radius of
c  the outer surface of plasma zone jz at phight(jh).
c  output arguments:
c  khstop(jh,jr)  1 if ray jh,jr strikes the limiting surface, 0 if not.
c  kznone(jh,jr)   ray jh,jr penetrates to zone kznone on first pass
c  through the plasma to the minimum minor radius.
c  kzntwo(ib,jh,jr)  ray jh,jr runs from zone kzntwo(1,jh,jr) to
c  kzntwo(2,jh,jr) after it passes its minimum minor radius.
c  pdvol(jz)   volume elements for plasma zones (sum over zones = 1).
c  pmub(jmu)   the lower boundary of the cosine of pitch-angle bins.
c  pdmuc(jmu)   the width of the cosine of pitch-angle bins.
c  output arguments:
c  phrdep(jz,jmu,3)  beam deposition function h(r) for each of 3 energies.
c
c  local:
c  zlenth(jz,ip)     length of ray  in zone jz on first
c  pass (ip=1) to minimum minor radius, or the second
c  pass (ip=2) after passing the minimum minor radius.
c
c  zero phrdep before accumulating contributions from beamlets.
       do 20 jfrac=1,3
         do 15 jmu=1,kmu
           do 10 jz=kcentr,kedge
             phrdep(jz,jmu,jfrac)=0.
10         continue
15       continue
20     continue
c
       do 420 jr=1,krtmax
         zzrtng=prtang(jr)
         zabsrt=abs(zzrtng)
         do 410 jh=1,khtmax
c          set default values of kznone, knztwo, and zlenth.
           kznone(jh,jr)=kedge+1
           kzntwo(1,jh,jr)=kedge+1
           kzntwo(2,jh,jr)=kedge
           do 140 jz=kcentr,kedge
             zlenth(jz,1)=0.
             zlenth(jz,2)=0.                      
c            if imu isn't set below an excess of deposition in the counter
c            direction will appear in the h(r,mu).
             imu(jz,1)=1
             imu(jz,2)=1
140        continue
c          check for existence of third and fourth legs (double pass).
           if(zabsrt.le.prmjwl-sqrt(prmnwl**2-phight(jh)**2)) then
             khstop(jh,jr)=1
           else
             khstop(jh,jr)=0
           endif
           do 150 jz=kedge,kcentr,-1
c            find the pitch-angle bin for this zone.
             if(phrmaj(2,jz,jh).ge.zabsrt) then
               zmu=zzrtng/phrmaj(2,jz,jh)
               do 110 jmu=kmu,1,-1
                 iimu=jmu
                 if(pmub(jmu).le.zmu) go to 111
110            continue
111            continue
               imu(jz,1)=iimu
             endif
c            check for prtang greater than the outer edge of the zone.
             if(phrmaj(2,jz,jh).le.zabsrt) then
               kznone(jh,jr)=jz+1
               kzntwo(1,jh,jr)=kedge+1
               go to 151
             endif
c
             if(jz.le.kminzn(jh)) then
c              this is necessarily the innermost zone of the first leg.
               kznone(jh,jr)=jz
c
               if(phrmaj(1,jz,jh).le.zabsrt) then
c                the min. major radius point is reached in this zone.
                 zlenth(jz,1)=sqrt(phrmaj(2,jz,jh)**2-zzrtng**2)
                 kzntwo(1,jh,jr)=kedge+1
                 go to 151
               else
                 zlenth(jz,1)=sqrt(phrmaj(2,jz,jh)**2-zzrtng**2)-
     1           sqrt(phrmaj(1,jz,jh)**2-zzrtng**2)
                 kzntwo(1,jh,jr)=jz+1
                 go to 151
               endif
             elseif(phrmaj(2,jz-1,jh).le.zabsrt) then
c              the min. major radius point is reached in this zone.
               zlenth(jz,1)=sqrt(phrmaj(2,jz,jh)**2-zzrtng**2)
               kznone(jh,jr)=jz
               kzntwo(1,jh,jr)=kedge+1
               go to 151
             else
               zlenth(jz,1)=sqrt(phrmaj(2,jz,jh)**2-zzrtng**2)-
     1         sqrt(phrmaj(2,jz-1,jh)**2-zzrtng**2)
             endif
150        continue
151        continue
c          second pass through the plasma zones to get to the minimum
c          major radius point.
           if(kzntwo(1,jh,jr).gt.kedge) go to 251
           do 250 jz=kzntwo(1,jh,jr),kedge
c            find the pitch-angle bin for this zone.
             if(phrmaj(1,jz-1,jh).ge.zabsrt) then
               zmu=zzrtng/phrmaj(1,jz-1,jh)
               do 210 jmu=kmu,1,-1
                 iimu=jmu
                 if(pmub(jmu).le.zmu) go to 211
210            continue
211            continue
               imu(jz,2)=iimu
             endif
             if(phrmaj(1,jz,jh).le.zabsrt) then
               kzntwo(2,jh,jr)=jz
               zlenth(jz,2)=sqrt(phrmaj(1,jz-1,jh)**2-zzrtng**2)
               go to 251
             else
               zlenth(jz,2)=sqrt(phrmaj(1,jz-1,jh)**2-zzrtng**2)-
     1         sqrt(phrmaj(1,jz,jh)**2-zzrtng**2)
             endif
250        continue
           kzntwo(2,jh,jr)=kedge
251        continue
           if(kznone(jh,jr).lt.kminzn(jh)) then
             if(imess.lt.100) then
               imess=imess+1
               callttyout(' seghfr: kznone(jh,jr).lt.kminzn(jh)')
             endif
           endif
c
           if((kzntwo(1,jh,jr).le.kznone(jh,jr)).and.
     1     (kznone(jh,jr).le.kedge)) then
             if(imess.lt.100) then
               imess=imess+1
               callttyout(' seghfr: kzntwo(1,jh,jr).le.kznone(jh,jr)')
             endif
           endif
c
           if(kzntwo(2,jh,jr).lt.kzntwo(1,jh,jr).and.
     1     (kzntwo(1,jh,jr).le.kedge)) then
             if(imess.lt.100) then
               imess=imess+1
               callttyout(' seghfr: kzntwo(2,jh,jr).lt.kzntwo(1,jh,jr)')
             endif
           endif
c-----------------------------------------------------------------------
c          calculate the neutral beam deposition along the path.
c          loop over the three energy components.
           do 380 jfrac=1,3
             if((jfrac.eq.1).and.(pffull.le.0.)) go to 380
             if((jfrac.eq.2).and.(pfhalf.le.0.)) go to 380
             if((jfrac.eq.3).and.(1.-pffull-pfhalf.le.0.)) go to 380
c            initialize the intensity of the beamlet
             zsorce=pweigt(jh,jr)
c            loop over the four legs of the path through the plasma.
             do 370 jjz=1,4
               if(jjz.eq.1) then
c                does the beamlet intersect the plasma?
                 if(kznone(jh,jr).gt.kedge) go to 371
c..set up for the first leg of the path to the min. minor radius.
                 izstrt=kedge
                 izend=kznone(jh,jr)
                 izdel=-1
                 iseg=1
               endif
               if(jjz.eq.2) then
c                test for the existence of the second leg.
                 if(kzntwo(1,jh,jr).gt.kedge) go to 370
c                set up for the second leg to the min. major radius along
c                the path.                                                    
                 izstrt=kzntwo(1,jh,jr)
                 izend=kzntwo(2,jh,jr)
                 izdel=1
                 iseg=2
               endif
               if(jjz.eq.3) then
c                does the path stop at the inner vacuum vessel wall?
                 if(khstop(jh,jr).gt.0) go to 371
c                if the second leg doesn't exist skip this section.
                 if(kzntwo(1,jh,jr).gt.kedge) go to 370
c                third leg is the reverse of the second.
                 izstrt=kzntwo(2,jh,jr)
                 izend=kzntwo(1,jh,jr)
                 izdel=-1
                 iseg=2
               endif
               if(jjz.eq.4) then
c                does the path stop at the inner vacuum vessel wall?
                 if(khstop(jh,jr).gt.0) go to 371
c                if the first leg doesn't exist, then quit.
                 if(kznone(jh,jr).gt.kedge) go to 371
c                the fourth leg is the reverse of the first leg.
                 izstrt=kznone(jh,jr)
                 izend=kedge
                 izdel=1
                 iseg=1
               endif
               do 340 jz=izstrt,izend,izdel
                 zx=pnesig(jz,jfrac)*zlenth(jz,iseg)
c                approximate exp(x) by 1+x+x**2/2+x**3/6
                 zexpi(jz)=1./(1.+zx*(1.+zx*(0.5+0.166667*zx)))
340            continue
               do 350 jz=izstrt,izend,izdel
c                reduce beamlet intensity.
                 zrsrc(jz)=zsorce
                 zsorce=zsorce*zexpi(jz)
350            continue
               do 360 jz=izstrt,izend,izdel
c                deposit particles in phrdep
                 phrdep(jz,imu(jz,iseg),jfrac)=
     1           phrdep(jz,imu(jz,iseg),jfrac)+
     2           zrsrc(jz)*(1.-zexpi(jz))
360            continue
370          continue
371          continue
380        continue
c
410      continue
420    continue
c  renormalize phrdep to put into standard form.
       do 520 jfrac=1,3
         do 515 jmu=1,kmu
           do 510 jz=kcentr,kedge
             phrdep(jz,jmu,jfrac)=phrdep(jz,jmu,jfrac)/
     1       (pdvol(jz)*pdmuc(jmu))
510        continue
515      continue
520    continue
       return
       end
c/ module aioniz
c
c.              %%%%
c.              %
c.              %       aioniz          %
c.              %
c.              %%%%
c
c  1-apr-82        added xenon; exenon,ixshel,nqxen
c  25-mar-81 changed berylium, added silicon, completed helium
c  2-feb-81 added berylium to ebind in aioniz
c
       function aioniz(kzfast,katom,pze,af,ksig)
       include 'cparm.m'
       dimension ebind(14,5,14),neshel(9,2),exenon(mj),nqxen(4)
c  levels: 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s
       data (neshel(l,1),l=1,9)/2,2,6,2,6,2,10,6,2/,
     1 (neshel(l,2),l=1,5)/2,8,18,18,18/
       data (nqxen(l),l=1,4)/53,45,27,9/
       data ebind(1,1,1)/13.6/
       data (ebind(2,1,jq),jq=1,2)/24.6,54.4/
       data (ebind(3,1,jq),jq=1,3)/68.,75.6,122./,
     1 ebind(3,2,1)/5.4/,
     2 (ebind(4,1,jq),jq=1,4)/129., 140., 154., 218./,
     3 (ebind(4,2,jq),jq=1,2)/9.3, 18.2/,
     4 (ebind(5,1,jq),jq=1,5)/209.,222.,240.,260.,340./,
     5 (ebind(5,2,jq),jq=1,3)/13.4,25.2,37.9/,
     6 ebind(5,3,1)/8.3/,
     7 (ebind(6,1,jq),jq=1,6)/309.,325.,346.,368.,392.,490./,
     8 (ebind(6,2,jq),jq=1,4)/19.6,31.4,47.9,64.5/,
     9 (ebind(6,3,jq),jq=1,2)/11.3,24.4/,
     x (ebind(8,1,jq),jq=1,8)/565.,580.,608.,635.,670.,705.,739.,871./,
     1 (ebind(8,2,jq),jq=1,6)/33.9,48.9,67.5,88.3,114.,138./
       data (ebind(8,3,jq),jq=1,4)/13.6,35.1,54.9,77.4/,
     1 (ebind(10,1,jq),jq=1,10)/890.,915.,940.,975.,1025.,1050.,
     2 1100.,1150.,1195.,1362./,
     3 (ebind(10,2,jq),jq=1,8)/52.5,71.,93.,117.,144.,173.,207.,239./,
     4 (ebind(10,3,jq),jq=1,6)/21.6,41.,63.5,97.,126.,158./,
     5 (ebind(11,1,jq),jq=1,11)/1100.,1110.,1140.,1180.,1220.,1260.,
     6 1310.,1360.,1410.,1465.,1650./,
     7 (ebind(11,2,jq),jq=1,9)/76.,83.,106.,133.,160.,192.,225.,
     8 264.2,300./,
     9 (ebind(11,3,jq),jq=1,7)/41.5,47.3,71.6,98.9,138.4,172.2,
     x 208.5/,ebind(11,4,1)/5.1/
       data (ebind(12,1,jq),jq=1,12)/1320.,1330.,1340.,1380.,1430.,
     1 1480.,1530.,1590.,1640.,1690.,1762.,1963./,
     2 (ebind(12,2,jq),jq=1,10)/103.,111.,122.,150.,180.,212.,
     3 247.,284.,328.,367./,
     4 (ebind(12,3,jq),jq=1,8)/62.,70.,80.,109.,141.,186.,
     5 225.,266./,
     6 (ebind(12,4,jq),jq=1,2)/7.6,15./
       data (ebind(14,1,jq),jq=1,14)/1870.,1880.,1890.,
     1 1910., 1930., 1980., 2030., 2090., 2160., 2230.,
     2 2300., 2370., 2440., 2670./,
     3 (ebind(14,2,jq),jq=1,12)/ 168., 177., 189., 204., 220., 256.,
     4 294., 335., 378., 425., 476., 523./,
     5 (ebind(14,3,jq),jq=1,10)/ 116., 125., 137., 152., 167., 205.,
     6 246., 303., 351., 401./,
     7 (ebind(14,4,jq),jq=1,4)/ 14.7, 22.8, 33.5, 45.1/,
     8 (ebind(14,5,jq),jq=1,2)/ 8.15, 16.3/
c
       data (exenon(jq),jq=1,54)/11.97, 23.54, 35.11, 46.68, 59.69,
     1 71.83, 98.05, 112.3, 170.8, 201.7, 232.6, 263.5, 294.4,
     2 325.3, 358.3, 389.6, 420.9, 452.2, 572.5, 607.7,
     3 642.9, 678.1, 726., 762.4, 852.7, 890.6, 1394., 1491.,1587.,
     4 1684.,1781., 1877., 1987., 2085., 2183., 2291., 2548., 2637.,
     5 2726., 2814., 3001., 3093., 3296., 3386., 7224., 7491., 7758.,
     6 8024., 8617., 8899., 9330., 9569., 39250., 40270./
c
c
c
c  ionization cross sections
c
       ixshel=1
       if(kzfast.eq.54) ixshel=2
       iel= kzfast+1-katom
       ishful= 0
       ii= 0
       do 20 j=1,4
         if((ii+neshel(j,ixshel)) .gt. iel) go to 25
         ii= ii +neshel(j,ixshel)
         ishful= j
  20   continue
  25   continue
       zsig= 0.
       if(ii.eq.iel) go to 29
c
c  do outer shell
c
       ishell= ishful+1
       if(ixshel.eq.1) zeb= ebind(kzfast,ishell,katom)
       if(ixshel.eq.2) zeb=exenon(katom)
       if(zeb .eq. 0.) go to 29
       zsig= float(iel-ii)*fsig(zeb,pze,af,ksig)
  29   continue
       if(ishful .eq. 0) go to 35
c
c  loop over filled shells
c
       do 30 j=1,ishful
         if(ixshel.eq.1) zeb= ebind(kzfast,j,katom)
         if(ixshel.eq.2) zeb=exenon(nqxen(j))
         if(zeb .eq. 0.) go to 30
         zsig= zsig +neshel(j,ixshel)*fsig(zeb,pze,af,ksig)
  30   continue
  35   continue
       aioniz=zsig
       return
       end
c
c/ module fsig
c.              %%%%
c.              %
c.               0.000000sig            %
c.              %
c.              %%%%
c
c  25-mar-81 added ksig=5 option
       function fsig(zeb,ef,af,ksig)
c  27-may-82    corrected gryzinski proton ionization formula
       data za0,za1,za2,za3,za4,za5,za6/-3.17385e1,1.143818e1,
     1 -3.833998,0.7046692,-7.431486e-2,4.153749e-3,-9.486967e-5/
       data zp0,zp1,zp2,zp3,zp4,zp5,zp6/-42.03309,3.557321,
     1 -1.045134,0.3139238,-7.454475e-2,8.459113e-3,-3.495444e-4/
c  corrected fg2: -1/galf before log
       fg2(z)=(z/(z+1.))**1.5*(1./z)*(z/(z+1.)+(2./3.)*(1.-
     1 1./galf)*log(2.7+sqrt(z)))*(1.-1./galf)*(1.-(1./galf)**(1+z))
c
c
       if(ksig .eq. 1) go to 100
       if(ksig.eq.2) go to 200
       if(ksig.eq.3) go to 300
       if(ksig .eq. 4) go to 400
       if(ksig.eq.5) go to 500
c
c  proton ionization from freeman and jones
c
       ze= log(ef/(af*1000.))
       zlgsig= zp0 +ze*(zp1+ze*(zp2+ze*(zp3+ze*(
     1 zp4+ze*(zp5+ze*zp6)))))
       fsig= exp(zlgsig)
       return
  100  continue
c
c  proton ionization from gryzinski
c
       zx=0.0233*sqrt(ef/(af*zeb))
       zz=zx**2
       galf=4.*zz*(1.+1./zx)
       fsig=8.79e-17*746.*fg2(zz)/zeb**2
       if(galf.le.1.) fsig=0.
       return
c
c
200    continue
c
c  proton ionization from garcia,gerjouy, and welker
c
       blimt= (sqrt(2.)-1.)/2.
       ulimt= blimt +1.
       alph= .0233*sqrt(ef/(af*zeb))
       if(alph .lt. blimt) fsig= 0.
       if(alph .lt. blimt) return
       if(alph .gt. ulimt) go to 210
       fsig=(-0.328 +2.*alph**3 +3./(8.*(alph+1.)))/
     1 (3.*(alph**2)*(zeb**2))
       go to 220
  210  continue
       fsig= (5. +3./(8.*(alph+1.)) +alph/(8.*(alph-1.)**2)
     1 -(4.*alph-1.)/(8.*alph**2) -1./(8.*alph*
     2 alph*(alph-1.)**2) -3./(4.*alph*(alph-1.)))/
     3 (3.*(alph**2)*(zeb**2))
  220  continue
       fsig=fsig*740.*8.79e-17
       return
300    continue
c
c  electron ionization (sigma-v-bar) from freeman and jones
c
       ze=log(ef*(13.6/zeb))
       zlgsig=za0+ze*(za1+ze*(za2+ze*(za3+ze*(za4+ze*(za5+ze*za6)))))
       fsig=exp(zlgsig)*(13.6/zeb)**1.5
       return
  400  continue
c
c  hydrogen charge-transfer from freeman and jones
c
       ze= log(ef/af)/2.30258
       fsig= 6.937e-15*(1.-0.155*ze)**2/(1.+1.112e-15*(ef/af)**3.3)
       return
c
c  electron loss by collision with ionized impurity
c
500    continue
c
c  electron loss by impurity ions
c  from olson. et al., phys. rev. letters, 41 (1978) 163
c
       ze=ef*13.6/(1000.*af*zeb*32.)
       fsig=4.6e-16*(13.6/zeb)**2*(1.-exp(-ze))/ze
       return
       end
c/ module ttyout
c
c.              %%%%
c.              %
c.              %       ttyout          %
c.              %
c.              %%%%
c
       subroutine ttyout(cinput)
c  sends a character string to the terminal.
       character cinput*(*),cformt*5
       data iutty/ 6/
c
       ilen=len(cinput)
       write(cformt,'(2h(a,i2,1h))') ilen
       write(iutty,cformt) cinput(1:ilen)
       return
       end

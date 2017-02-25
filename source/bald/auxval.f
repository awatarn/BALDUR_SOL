c--------1---------2---------3---------4---------5---------6---------7-c
c@auxval  .../baldur/code/bald/auxval.f
c  rgb 10-feb-01 Count number of elements that are positive in the
c     rlepwr and rlipwr arrays and use that count to set nrlepwr and
c     nrlipwr in common block comtok in cbaldr.m
c  rgb 12-jun-00 Check rlepwr and rlipwr arrays
c    Negative elements are given value of previous element
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  pis 08-jul-98 moved calculation of nxwexba, ntwexba to beginning  
c  pis 07-may-98 introduced calculation of nxwexba, ntwexba
c  alh 30-oct-97 ensure time dependent tcold & tcoldp value make sense
c  rgb 14-aug-96 replaced zfuzz with rndup = 1.0 + rndeps
c  ahk 24-jun-96 forced izones = nzones after 418
c  rgb 28-may-96 added it = 0 before do 338
c  rgb 04-jan-95 kept ita and itb .le. imax,  xzoni index > 0
c  rgb 07-feb-94 changed edge of parabolic profiles to xzoni(ledge+1)
c  les   nov-90 added d-3he initialization of ions
c            impurity requires nimp(6), wtimp(6) input
c  rgb 18.11 09-mar-90 remove sediti, sploti, tediti, tploti
c  rgb 16.19 27-aug-89 test for hibeam.gt.epslon .or. hpowmw.gt.epslon
c  dps 15.01 19-sep-88 Fix bugs in scrape-off bpols calculation: add 
c            factor of proton mass, replace rmins with rmini, and
c            insert emu0 and twopi factors.
c  rgb 12.54 Corrected TES(2,MZONES) and TIS(2,MZONES) to be consistent
c       with either fitted, tabular, or time-dependent boundary
c       conditions.
c  rgb 12.54 Implemented options to initialize time-dependent edge 
c       values according to bdhyde, bdimpe, bdtie, bdtee
c  dps 23-mar-87 introduce message/check for flimp < 0.
c       rgb 17-may-85 added emu0 in commhd for MKS ISU internal units
c       mhr 15-june-84 added analytic test of bpoloidal diffusion eqn
c       fgps 11-feb-83 modified to handle more than 2 impurity species.
c       aes 16-apr-82 comment out lines with cfutz(139),nadump(20)
c       aes 16-apr-82 eliminate nfusn=3 option
c       aes 13-apr-82 set nadump(11) as predictor-corrector flag
c               for neu-imp influx handled in imprad(2) (same as
c               fgps 18-sep-81)
c       aes 02-mar-82 bypass divertor effects if connection length
c               cfutz(127) .gt. 1.e18 cm (same as fgps 9-sep-81)
c       aes 02-mar-82 move setting of mxt1 to below label 10
c       aes 14-jan-82 allow setting nfusn=3 or "forcing" nfusn=1
c       aes 29-oct-81 array dimension 52-->55 in zjs, 105-->111 in zdens
c       aes 7-aug-81 fixed nfusn bug - always call alphas
c       aes 7-aug-81 set gfract(jt,2) in do 380
c       aes 31-jul-81 added zfract=0. near line 600
c       fgps 24-feb-81 changed "mxt" to "mxt1" in "do 402"
c       aes 1-dec-80 add epslon to test near do 406 so that rminor can
c               be set equal to raidus(jz)
c       aes 19-aug-80 set tcoldp after line 379
c       fgps 11-oct-79 added the capability of prescribing field-
c                      ripple coefficients at specified times
c       fgps 27-jul-79 namelist nadump(*) parameters replaced
c                      by cfutz(*)-factors
c       fgps 12-jul-79 impurity ions included in initial
c                      scrapeoff flow velocity
c       fgps 17-apr-79 modified initial b-poloidal in the
c                      divertor region
c       dhfz 5-apr-79 set lhe3
c       dhfz 5-apr-79 set nfusn
c       dhfz 28-mar-79 fix helium-3 bug at 'do 238'
c       dhfz 30-jan-79 fix bug in bz computation
c       dhfz oct-14-78 put nnhfit,nnifit,nnfit,ntefit,ntifit,
c               nbfit back in
c       dhfz 10-sept-78 modify to avoid neg. bases
c       dhfz 24-july-78 make jz proportional to
c               (1.- (r/a)**eebfit)**ebfit
c       dhfz 20-july-78 change nnhfit to eehfit,nnifit to eeifit,
c               nnfit to eefit, ntefit to eeteft, ntifit to eetift,
c               and nbfit to eebfit
c       dhfz 20-july-78 add exponents ehfit,eifit,efit, etefit, etifit,
c               and ebfit
c       amck 27-jan-78 change subscript exp. which are ill. on icl
c       amck 11-jan-78 don't set nitmax
c       amck 30-oct-77 set lalpha. set ldeut,ltrit,lhbeam only .lt.limp1
c       amck 21-sep-77 jg -> jb in nhbeam
c       amck 25-aug-77 add dynamic default val. for beams
c       amck 16-aug-77 add short formula for bint
c       amck 10-aug-77 initialization call to "beams"
c       amck 8-aug-77 initial setting of beam variables
c       amck 8-aug-77 modify for ngas = z of element
c       amck 2-aug-77 add return after parag. 9
c       amck 2-aug-77 nzones -> izones in paragraph 5
c       amck 29-jul-77 nzones is max no. of zones
c       amck 29-jul-77 move setting of comtim vars. to before comflg
c       amck 15-jul-77 change ihteti (which goes into lhspec) to energy
c       amck 27-may-77 fix tabular te, ti input
c       amck 28-apr-77 fix time point 1 in gflowi
c       amck 21-apr-77 change 975 limit to mimp
c       amck 21-apr-77 check no. of particle species against mxchi
c       amck 20-apr-77 set xbouni acc. to xzoni, use denga0, etc.
c       amck 20-apr-77 make gflowi handle impurities also
c       amck 15-mar-77 call imprad(1) after 3.3)
c       amck 15-mar-77 call imprad(1) after 1004
c       amck 8-mar-77 call errchk
c       amck 18-feb-77 recompute iphyd before setting gflowi
c       amck 17-feb-77 clear tcompi, redgi, rmji, and rcurri before #900
c       amck 15-feb-77 add restart auxiliary values
c       amck 11-feb-77 move mx-- to preset
c       amck 26-jan-77 permute ih index when going from ext. to int.
c       amck 18-jan-77 don't fill in missing points in gflowi
c       amck 18-jan-77 add gtflwi and gflowi setting
c       amck 11-jan-77 add call expert(,,9)
c       amck 6-jan-77 set plotting flags
c       amck 6-oct-76 put "call error"s in
c       amck 6-oct-76 fix indexing bug in packing of fracth/i
c       amck 6-oct-76 dengas=0 or fracth=0 same as ngas=0,
c               ngas not 0-delimited
c       amck 22-sep-76 add tmax(i)
c       amck 2-sep-76 fix do #216, #236 loop limits
c       amck 27-aug-76 fix 'i0.ge.i1'
c       amck 26-aug-76 add setting of tcompi, redgi, rmji, rcurri
c       amck 24-aug-76 set mxt
c       amck 19-aug-76 make bedgi a scalar
c       amck 16-aug-76 fix sign in bint computation
c       amck 10-aug-76 add bint
c       amck 30-jul-76 make eebfit refer to jz fit
c       amck 29-jun-76 mxparm -> mxchi, add mxhyd, mximp, mxzone
c       amck 3-jun-76 add max to zdens computation
c       amck 2-jun-76 change dummy zone centers from on edges to outside
c               of edges.
c       amck 24-may-76 fix refs. to rminai, "nskipi", remove "zfr"
c       amck 20-may-76 change bpolai,bpolbi to bpols, (t)bpoli to (t)bedgi,
c               rmajai, -bi rminai, -bi to rmini, rmaji, rmini and rmaji
c               to redgi and rmji
c       amck 5-may-76 change cbel -> cnuel, cbhyd -> cnuhyd, etc.
c       amck 19-apr-76 #2
cdoc
c******************************************************************************
c
c       -----------
c       sbrtn AUXVAL   file DEFAULT
c       ------------
c
c       1.5     compute auxiliary values from input data
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       comsta, comflg, comdim, comout, comin, comdf2
c
c------------------------------------------------------------------------------
c
c       when field-ripple is active, the array delimiter mxt1 must
c       be reduced from mxt=20 to mxt/2=10, because field-ripple par-
c       ameters are stored in the last 10 elements of tcompi(20),
c       rcurri(20), redgi(20), and rmji(20).
c
c       order of precedence in specification of values:
c
c       curent is ignored if bpoid(1) is not 0 (i.e., is specified)
c
c
c       if dens is specified, i.e., dens(1) is not 0, dens0, dens1, and eefit
c       are ignored.
c
c       If ( dens(1) .eq. 0. .and. bdhyde(1,jh) .gt. epslon )
c               the time dependent edge density is computed at t=tinit
c
c       if radius is specified, i.e., radius(2) is not 0, rminor and nrfit
c       are ignored
c
c       if te is specified, i.e., te(1) is not 0, te0, te1, and eeteft are
c       ignored.
c       similarly, if ti is specified, ti0, ti1, and eetift are ignored.
c       If ( te(1) .eq. 0. .and. bdtee(1) .gt. epslon ), or similarly 
c       if ( ti(1) .eq. 0. .and. bdtie(1) .gt. epslon ), then
c               the time dependent edge temperature is computed
c               at t = tinit
c
c       if wtgas(j) is not specified, i.e. = 0, the weight is taken to be
c       the standard value for that element, as specified by ngas(j).
c       similarly, wtimp(j) is defaulted to the value specified by nimp(j).
c
c       defaulting of boun**(ip) is done by subroutine bounds.
c       see subroutine bounds for details.
c
cend
c
c
        subroutine auxval
c
       include 'cparm.m'
       include 'cbaldr.m'
       include 'cbparm.m'
       include 'commhd.m'
c
        logical ila, ilb
c
        dimension
     i  iphyd(6),       ipimp(idximp),  izdef2(5),
     r  zdens(111),     zadef(2),       zadef2(5),     zjs(mj),
     r  zroot1(5)
c
        character *10 ihspec(5), ihteti(2)
c
        data
     h  ihspec/ 'hydrogen-1',
     h          'deuterium ',
     h          ' tritium  ',
     h          ' helium-4 ',
     h          ' helium-3 '/
        data
     h  ihteti/ 'el. energy',
     h          'ion energy'/
        data    irfit/2/,       izdef2/1,1,1,2,2/
        data    zadef/1.007825,4.00260/
        data    zadef2/1.007825,2.0140,3.01605,4.00260,3.0/
        data    zcross/1.03/
        data    zroot1/3.8317,7.0156,10.1735,13.3237,16.4706/
c
c       cfutz(ikirip) activates field-ripple effects and ikirip
c       must also be initialized in subroutine trcoef
c
        data    ikirip/140/
c
c
c------------------------------------------------------------------------------
c
        data    iclass /1/,     isub /5/
c
c
        if (.not.nlomt1(isub)) go to 10
cap
        call mesage('*** 1.5 subroutine auxval bypassed ***         ')
        return
   10   continue

c
c    Default size lengths (55 and 20) are hardcoded in cbaldr.m
c    1 point minimum is assumed!
c
c... Initialize array size for wexba
c
        if (twexba(1) .eq. twexba(2)) then
          ntwexba  = 0
        else 
          ntwexba = 1
          DO j1 = 2, 20
             if (twexba(j1) .gt. twexba(j1-1)) ntwexba = ntwexba + 1
          END DO
        end if

        if (xwexba(1) .eq. xwexba(2)) then
          nxwexba = 0
        else
          nxwexba = 1
          DO j1 = 2, 55
           if (xwexba(j1) .gt. xwexba(j1-1)) nxwexba = nxwexba + 1
          END DO
        end if
 

c
        mxt1=mxt
        if(cfutz(ikirip).gt.epslon) mxt1=10
        if (nlres) go to 300
c
c
cl      1)      check for errors, set up conversion factors
c
c
c
        call errchk(1)
c
        call units
cdoc
c
c       namelist input parameters cfutz(120) through cfutz(139) are
c       assigned to prescribe properties of any existing scrapeoff
c       layers; in particular:
c        cfutz(120)=innermost radial index of innermost scrapeoff
c                   region or of innermost separatrix
c        cfutz(121)=innermost radial index of second scrapeoff
c                   layer or second separatrix if such exists
c        cfutz(123)=switch for activating charge-exchange friction
c                   in the scrapeoff layers; 0.=off; 1.=on
c        cfutz(126)=outermost radial zone subject to normal timestep
c                   controls
c        cfutz(127)=average path-length (cm) parallel to b in an
c                   innermost scrapeoff layer
c        cfutz(128)=average path-length (cm) parallel to b in a
c                   second scrapeoff layer
c        cfutz(129)=minimum value of electron density (1./cm**3)
c                   in the scrapeoff layers
c        cfutz(130)=minimum value of the electron temperature (ev)
c                   in the scrapeoff layers
c        cfutz(131)=minimum value of the ave ion temperature (ev)
c                   in the scrapeoff layers
c        cfutz(idivrt)=a prescribed fixed value, in terms of *bohm,
c                   for particle diffusion coefficients and, in
c                   terms of *1.5*bohm, for thermal conductivities
c                   in all scrapeoff layers; this overrides all other
c                   determinations of the transport quantities
c
c       the next four parameters involve modeling a straight-line
c       increase in the transport coefficients starting from some
c       prescribed point inside the plasma and rising to fixed
c       values in the scrapeoff layers (hawryluk, jul-79)
c        cfutz(iscrp1)=the fixed value of the hydrogen-ion diffusion
c                   coefficients in the scrapeoff layers
c        cfutz(iscrp2)=the fixed value of the impurity-ion diffusion
c                   coefficients in the scrapeoff layers
c        cfutz(iscrp3)=the fixed value of the electron thermal
c                   conductivity in the scrapeoff layers
c        cfutz(iscrp4)=the starting point for the increase stated
c                   as a fraction multiplying the innermost radius
c                   of the scrapeoff layers
c
c        cfutz(139)=switch for activating a print-out of the
c                   transport model at every radial zone;
c                   de-activated 2-sep-80
c
c
c       as a holdover from earlier versions of the code, some of the
c       input scrapeoff parameters are handled in terms of the integer
c       array nadump(*)
cend
        nadump(1)=0
        if(cfutz(120).gt.epslon) nadump(1)=int(cfutz(120)+0.1)
        nadump(3)=100
        if(cfutz(121).gt.epslon) nadump(3)=int(cfutz(121)+0.1)
        nadump(4)=0
        if(cfutz(123).gt.epslon) nadump(4)=int(cfutz(123)+0.1)
        nadump(5)=1
        nadump(6)=1
        nadump(7)=100
        if(cfutz(126).gt.epslon) nadump(7)=int(cfutz(126)+0.1)
        nadump(8)=1
        nadump(10)=0
cdoc
c
c       nadump(11) is set up as a flag for use in subroutine imprad(2)
c       which handles neu-imp influx:  "0" prevents any neu-imp influx;
c       "1" unlocks predictor mode for neu-imp influx; and, "2" allows
c       use of a corrector mode which presists until nadump(11) is re-
c       set to "1".  to avoid initialization difficulties, nadump(11)=
c       0 here.  in subroutine start, where stepwise advance begins,
c       nadump(11)=1.
cend
        nadump(11)=0
c       nadump(20)=0
c       if(cfutz(139).gt.epslon) nadump(20)=int(cfutz(139)+0.1)
c
c
cl      2)      set ion species variables
c
c
        lhyd1 = 1
        lhydn = 0
        mhyd = 0
        mimp = 0
c
cl      2.1)    set hydrogen variables
c
        do 208 j = 1, mxhyd
        if ((dengas(j,1).le.epslon.and.fracth(j).le.epslon.and.
     1          denga0(j).le.epslon).or.ngas(j).eq.0) go to 208
c
        lhydn = lhydn + 1
        mhyd = mhyd + 1
        if (lhydn+2.gt.mxchi) go to 9020
c
        nspec(lhydn) = ngas(j)
        aspec(lhydn) = wtgas(j)
        iphyd(lhydn) = j
  208   continue
c
c
  210   continue
c
cl      2.2)    set impurity variables
c
        limp1 = lhydn + 1
        limpn = lhydn
c
        do 218 j = 1, mximp
        if ((denimp(j,1).le.epslon.and.fracti(j).le.epslon.and.
     1          denim0(j).le.epslon).or.nimp(j).eq.0) go to 218
c
        limpn = limpn + 1
        mimp = mimp + 1
        if (limpn+2.gt.mxchi) go to 9020
c
        nspec(limpn) = nimp(j)
        aspec(limpn) = wtimp(j)
        ipimp(mimp) = j
  218   continue
c
c
  220   continue
c
c   les  nov-90  d-3he fusion
c      default ngas(1)=d, (2)=t, nimp(1)=p, (2)=3he, (3)=4he,
c       (4)=regular impurity if 6Li, 7Li not used
c      requires nimp(lprotn)=-1 in input
c   treats 3he, 4he specially, not as regular impurity
c      requires input nimp(6) and wtimp(6) for impurity
      if ( cfutz(490) .le. epslon ) go to 221
        lhe3=4
        lalpha=5
      lprotn=3
        ldeut=1
        ltrit=2
        nspec(3)=-1
        nspec(4)=-4
        nspec(5)=-5
        nzspec(1)=1
        nzspec(2)=1
        nzspec(3)=1
        nzspec(4)=2
        nzspec(5)=2
c   otherwise, specify wtgas and wtimp in input
          aspec(1)=2.014
          aspec(2)=3.01605
          aspec(3)=1.007825
          aspec(4)=3.0
          aspec(5)=4.00260
        do 2342 js=1,5
          i002=-nspec(js)
          lhspec(js) = ihspec(i002)
 2342   continue
      nfusn = 4
c   call imprad to initialize regular impurity (mimp=4)
c     tables for impurity 4 have first index = 1, not 4 as usual
c      change call to ofit
      if(mimp.lt.4) go to 2440
      nzspec(limpn)=nspec(limpn)
      limp1z=limp1
      mximpz=mximp
      mimpz=mimp
      limp1=6
      mximp=1
      mimp=1
      call imprad(1)
c
c   this gives the radiation coefficients for impurity 6 with
c      index 1, not 4=mimp
      limp1=limp1z
      mximp=mximpz
      mimp=mimpz
      go to 2440
 221  continue
c
c
c      original baldur
c
        lhe3=0
        ltrit = 0
        ldeut = 0
c
        do 238 js = lhyd1, limpn
        i001 = nspec(js)
        if(i001.eq.-4) go to 224
        if (i001.lt.0) go to 228
        nzspec(js) = i001
        if (i001.gt.2) go to 238
c
        if (aspec(js).le.0.0) aspec(js) = zadef(i001)
        go to (222,224),i001
c
c               ngas = 1 -- decide from wtgas if hyd., deut, or trit.
c
  222   continue
        nspec(js) = -1
        if (aspec(js).gt.1.5) nspec(js) = -2
        if (aspec(js).gt.2.5) nspec(js) = -3
        go to 228
c
c               ngas = 2 -- decide from wtgas if he-3 or he-4
c
  224   continue
        nspec(js) = -5
        if (aspec(js).lt.3.5) nspec(js) = -4
c
  228   continue
c
c
c               hydrogen or helium -- ngas .lt. 0
c
c
        if (nspec(js) .lt.-5) go to 9021
        i002 = - nspec(js)
        if (aspec(js).le.0.0) aspec(js) = zadef2(i002)
        nzspec(js) = izdef2(i002)
c
        if (nspec(js).eq.-5) lalpha = js
        if(nspec(js).eq.-4) lhe3=js
        if (js.gt.lhydn) go to 230
        if (nspec(js).eq.-2) ldeut = js
        if (nspec(js).eq.-3) ltrit = js
  230   continue
c
        lhspec(js) = ihspec(i002)
c
  238   continue
c
        ideut=0
        itrit=0
        ihe3=0
c
        do 240 jh=1,mhyd
        if(ngas(jh).eq.-2) ideut=1
        if(ngas(jh).eq.-3) itrit=1
240     continue
c
        if(mimp.eq.0) go to 242
c
        do 241 ji=1,mimp
        if(nimp(ji).eq.-4) ihe3=1
 241    continue
c
 242    continue
c
        if(nfusn.eq.1) go to 243
        nfusn=1
        if(ideut*ihe3.ne.0) nfusn=2
 243    continue
c
c
c
        call imprad(1)
c
c
 2440  continue
c
        if (mimp.le.0) go to 250
c
        do 244 ji = 1, mimp
        do 244 ji2 = 1, mimp
        ii = ji + limp1 - 1
        ii2 = ji2 + limp1 - 1
        aired(ji,ji2) = aspec(ii) * aspec(ii2) /
     1                  (aspec(ii) + aspec(ii2))
  244   continue
c
  250   continue
c
c
cl      3)      time-dependent conditions
c
c
  300   continue
c
c
c               clear tcompi, redgi, rcurri, rmji
c               (necessary if restart)
c
c
        call resetr(tcompi,mxt,0.0)
        call resetr(rmji  ,mxt,0.0)
        call resetr(redgi ,mxt,0.0)
        call resetr(rcurri,mxt,0.0)
cdoc
c
c               merge tbpoid and tcomp
c
c               jt is the index for tcompi, redgi, rcurri, and rmji
c               ita for tcomp, rmajor, and rminor, and itb for
c               tbpoid, curent, and bpoid.
c               zt, the next time in either tbpoid or tcomp,
c               and hence tcompi(jt),
c               has been already computed at the beginning of
c               the loop.  hence, either tbpoid(itb) or tcomp(ita)
c               or both will be within a factor rndup of zt.
c               if tcomp(ita) is almost zt, redgi is set from rminor
c               and rmji from rmajor (if they are not 0).
c               similarly, if tbpoid(itb) is almost zt, rcurri(jt)
c               is set to bpoid if bpoid is not 0, and to -curent(itb)
c               if bpoid(itb) is 0 and curent(itb) not 0.
c               then, if tcomp(ita) was zt, ita is advanced,
c               and if tcomp(ita) is now less than zt,
c               we set ila=.false. to indicate we have
c               reached the end of the tcomp arrays,
c               if tcomp(ita) almost =zt, check the next rmajor(ita),
c               and if greater, go on to the next jt.
c               similarly for the tbpoid loop (#334).
c
c               then rmji is linearly interpolated (over time),
c               then redgi is set to rminor/sqrt(rmajor) and
c               interpolated, and, finally, rcurri is set
c               to rmajor*curent and interpolated.
c               since before conversion and interpolation,
c               rcurri could be either the current or bpoid,
c               the sign is used to flag which it is.
c               note that the interpolation *must* be done in
c               this order.
c
cend
c
        it = 0
        zt = 0.0
        ita = 1
        itb = 1
        ila = .true.
        ilb = .true.
c
        imax = 20
c
        rndup = 1.0 + rndeps
c
        do 338 jt = 1, mxt1
        tcompi(jt) = zt
        zt = epsinv
c
c               get values near time tcompi(jt)
c
c               get rminor and rmajor
c
  332   continue
        if (.not.ila) go to 334
        if (tcompi(jt).gt.tcomp(ita) * ueit*rndup.or.ita.gt.mxt1)
     1          ila = .false.
        if (.not.ila) go to 334
        if (tcompi(jt)*rndup.lt.tcomp(ita)*ueit) go to 333
c
        if (rmji(jt).eq.0.0.and.rmajor(ita).gt.epslon)
     1          rmji(jt) = rmajor(ita) * ueil
        if (redgi(jt).eq.0.0.and.rminor(ita).gt.epslon)
     1          redgi(jt) = rminor(ita) * ueil
        ita = ita + 1
        if ( ita .gt. imax ) goto 340
        go to 332
c
  333   continue
        zt = min(zt,tcomp(ita)*ueit)
c
c               get curent or bpoid
c
  334   continue
        if (.not.ilb) go to 336
        if (tcompi(jt).gt.tbpoid(itb) * ueit*rndup.or.itb.gt.mxt1)
     1          ilb = .false.
        if (.not.ilb) go to 336
        if (tcompi(jt)*rndup.lt.tbpoid(itb)*ueit) go to 335
c
        if (rcurri(jt).eq.0.0.and.bpoid(itb).gt.epslon)
     1          rcurri(jt) = bpoid(itb) * ueib
        if (rcurri(jt).eq.0.0.and.curent(itb).gt.epslon)
     1          rcurri(jt) = - curent(itb) * ueii
        itb = itb + 1
        if ( itb .gt. imax ) goto 340
        go to 334
c
  335   continue
        zt = min(zt,tbpoid(itb)*ueit)
  336   continue
c
        it = jt
        if (.not.(ila.or.ilb)) go to 340
  338   continue
c
  340   continue
cdoc
c
c       when field-ripple effects are active, prescribed time
c       points are stored in tcompi(11-20), and ripple coefficients
c       at these times are loaded into rcurri(11-20), redgi(11-20),
c       and rmji(11-20)
cend
        if(mxt1.ge.mxt) go to 2340
        jtmin=mxt1+1
        do 1340 jt=jtmin,mxt
        tcompi(jt)=tbpoid(jt)*ueit
        rcurri(jt)=bpoid(jt)
        redgi(jt)=curent(jt)
        rmji(jt)=rmajor(jt)
 1340   continue
 2340   continue
c
c
c               fill in intermediate rmji's
c
c
        z0 = 0.0
        z1 = 0.0
        i0 = 1
        zt0 = 0.0
        zt1 = 0.0
        zdt = 1.0
c
        do 349 jt = 1, it
        if (jt.ne.it.and.rmji(jt).eq.0.0) go to 349
        z0 = z1
        zt0 = zt1
c
        if (rmji(jt).eq.0.0) go to 341
        z1 = rmji(jt)
        zt1 = tcompi(jt)
c
        if (z0.ne.0.0) go to 341
        z0 = z1
        zt0 = zt1
  341   continue
c
        i1 = jt - 1
        if (i0.gt.i1) go to 348
c
        zdt = zt1 - zt0
        if (zdt.le.epslon) go to 343
        zdt = 1.0 / zdt
c
        do 342 jt2 = i0, i1
        rmji(jt2) = (z0*(zt1 - tcompi(jt2)) + z1*(tcompi(jt2) - zt0))
     1                                                  * zdt
  342   continue
c
        go to 348
  343   continue
c
        do 344 jt2 = i0, i1
        rmji(jt2) = z1
  344   continue
c
  348   continue
        if (rmji(jt).eq.0.0) rmji(jt) = z1
        i0 = jt + 1
  349   continue
c
c
c               fill in intermediate redgi
c
c
        z0 = 0.0
        z1 = 0.0
        i0 = 1
        zt0 = 0.0
        zt1 = 0.0
        zdt = 1.0
c
        do 359 jt = 1, it
        if (jt.ne.it.and.redgi(jt).eq.0.0) go to 359
        z0 = z1
        zt0 = zt1
c
        if (redgi(jt).eq.0.0) go to 351
        z1 = redgi(jt) / sqrt(rmji(jt))
        zt1 = tcompi(jt)
c
        if (z0.ne.0.0) go to 351
        z0 = z1
        zt0 = zt1
  351   continue
c
        i1 = jt - 1
        if (i0.gt.i1) go to 358
c
        zdt = zt1 - zt0
        if (zdt.le.epslon) go to 353
        zdt = 1.0 / zdt
c
        do 352 jt2 = i0, i1
        redgi(jt2) = (z0*(zt1 - tcompi(jt2)) + z1*(tcompi(jt2) - zt0))
     1                                                  * zdt
  352   continue
c
        go to 358
  353   continue
c
        do 354 jt2 = i0, i1
        redgi(jt2) = z1
  354   continue
c
  358   continue
        redgi(jt) = z1
        i0 = jt + 1
  359   continue
c
c
c               fill in intermediate rcurri
c
c
        z0 = 0.0
        z1 = 0.0
        i0 = 1
        zt0 = 0.0
        zt1 = 0.0
        zdt = 1.0
c
        do 369 jt = 1, it
        if (jt.ne.it.and.rcurri(jt).eq.0.0) go to 369
        z0 = z1
        zt0 = zt1
c
        if (rcurri(jt).eq.0.0) go to 361
        if (rcurri(jt).gt.0.0)
     1          z1 = rcurri(jt)*rmji(jt)*sqrt(rmji(jt))*redgi(jt)
        if (rcurri(jt).lt.0.0)
     1          z1 = - rcurri(jt)*rmji(jt)
        zt1 = tcompi(jt)
c
        if (z0.ne.0.0) go to 361
        z0 = z1
        zt0 = zt1
  361   continue
c
        i1 = jt - 1
        if (i0.gt.i1) go to 368
c
        zdt = zt1 - zt0
        if (zdt.le.epslon) go to 363
        zdt = 1.0 / zdt
c
        do 362 jt2 = i0, i1
        rcurri(jt2) = (z0*(zt1 - tcompi(jt2)) + z1*(tcompi(jt2) - zt0))
     1                                                  * zdt
  362   continue
c
        go to 368
  363   continue
c
        do 364 jt2 = i0, i1
        rcurri(jt2) = z1
  364   continue
c
  368   continue
        rcurri(jt) = z1
        i0 = jt + 1
  369   continue
c
        tcompi(it+1) = 0.0
c
c
c               fill in gflowi's
c
c
        call resetr(gtflwi,mxt,0.0)
        call resetr(gflowi,mxchi*mxt,0.0)
c
c
cl              set iphyd -- permutation of h index from
c               flgas to gflowi
c
c
        i1 = 1
        do 371 jh = 1, mxhyd
        if ((dengas(jh,1).le.epslon.and.denga0(jh).le.epslon.and.
     1          fracth(jh).le.epslon).or.ngas(jh).eq.0) go to 371
        iphyd(i1) = jh
        i1 = i1 + 1
  371   continue
c
c
        if (mimp.le.0) go to 373
        i1 = 1
        do 372 ji = 1, mximp
        if ((denimp(ji,1).le.epslon.and.denim0(ji).le.epslon.and.
     1          fracti(ji).le.epslon).or.nimp(ji).eq.0) go to 372
        ipimp(i1) = ji
        i1 = i1 + 1
  372   continue
c
  373   continue
c
c
cl              fill in gflowi, gtflwi
c
c
        ith = 1
        iti = 1
        gtflwi(1) = 0.0
        z0 = uiel**2 * uiet
        zt = 0.0
        inegfl=0
c
        do 378 jt = 1, mxt
        if (jt.le.1) go to 374
        zt = epsinv
        if (tgas(ith).gt.epslon) zt = tgas(ith)
        if (timp(iti).gt.epslon.and.mimp.gt.0) zt = min(zt,timp(iti))
        if (zt.ge.epsinv) go to 379
        gtflwi(jt) = zt*ueit
  374   continue
c
        zt0 = 0.0
        if (ith.gt.1) zt0 = tgas(ith-1)
        zint = 0.0
        if (ith.le.1) zint = 1.0
        if (tgas(ith).gt.epslon) zint = (zt - zt0) / (tgas(ith) - zt0)
c
        do 375 jh = 1, mhyd
        zf = 0.0
        i001 = iphyd(jh)
        if (ith.gt.1) zf = flgas(i001,ith-1)
        gflowi(jh,jt) = (zint*flgas(i001,ith) + (1.0-zint)*zf)*z0
  375   continue
c
        if (mimp.le.0) go to 377
        zt0 = 0.0
        if (iti.gt.1) zt0 = timp(iti-1)
        zint = 0.0
        if (iti.le.1) zint = 1.0
        if (timp(iti).gt.epslon) zint = (zt - zt0) / (timp(iti) - zt0)
c
        do 376 ji = 1, mimp
        zf = 0.0
        i001 = ipimp(ji)
        i002 = ji+limp1-1
        if (iti.gt.1) zf = flimp(i001,iti-1)
        gflowi(i002,jt) = (zint*flimp(i001,iti) +
     1                                          (1.0-zint)*zf)*z0
        if (flimp(i001,iti).lt.-epslon) inegfl=1
  376   continue
  377   continue
c
        if (tgas(ith).le.zt*rndup.and.
     1                  (tgas(ith).ge.epslon.or.ith.le.1)) ith = ith + 1
        if (tgas(ith).le.zt*rndup) tgas(ith) = 0.0
c
        if (mimp.le.0) go to 378
        if (timp(iti).le.zt*rndup.and.
     1                  (timp(iti).ge.epslon.or.iti.le.1)) iti = iti + 1
        if (timp(iti).le.zt*rndup) timp(iti) = 0.0
  378   continue
c
c
  379   continue
        if (inegfl.ne.0) then
          if (cfutz(200).eq.0.0) then
            call mesage('  *** negative impurity ion flux specified ')
          else
            go to 9120
          end if
        end if
c
c               tcoldp defaults to tcold
c
        if (tcoldp(1) .le. 0.0) tcoldp(1) = tcold(1)
c
cl              empty elements of the arrays assume last given value
c
        do j = 2, 20
           if (tcold(j) .le. 0.0) tcold(j) = tcold(j-1)
           if (tcoldp(j) .le. 0.0) tcoldp(j) = tcoldp(j-1)
        enddo
c
cl  Check rlepwr and rlipwr arrays
c     nrlepwr = number of positive elements in rlepwr array
c     nrlipwr = number of positive elements in rlipwr array
c
        nrlepwr = 0
        nrlipwr = 0
c
        do j = 1, 20
           if (rlepwr(j) .le. 0.0) exit
           nrlepwr = j
        enddo
c
        do j = 1, 20
           if (rlipwr(j) .le. 0.0) exit
           nrlipwr = j
        enddo
c
        write (6,*) nrlepwr,' = nrlepwr'
        write (6,*) nrlipwr,' = nrlipwr'
c
c
cl              set gfract(jt,2)
c
        do 380 jt=1,mxt
        gfract(jt,2) = (1. - abs(gfract(jt,1))) * sign(1.,gfract(jt,1))
  380   continue
c
c
c
        if (nlres) go to 1000
c
cl      4)      set flags and dimensions
c
cl      4.1)    comdim
c
  400   continue
c
c               compute initial nzones (nzones is max over compression)
c
        zrad = redgi(1)
        it = 1
c
        do 402 jt = 2, mxt1
        if (tcompi(jt).lt.0.0) go to 403
        if (zrad.ge.redgi(jt)) go to 402
        zrad = redgi(jt)
        it = jt
  402   continue
  403   continue
c
        if (zrad.le.epslon) go to 9050
        zrad = zrad / redgi(1)
c
        if (radius(2).le.epslon) go to 410
c
        do 406 jz = 2, nzones
        if (radius(jz)+epslon.ge.rminor(1)) go to 418
        izones = jz
  406   continue
c
        go to 418
c
  410   continue
        go to (412,414), nrfit
c
  412   continue
        izones = int(float(nzones)/zrad)
        go to 418
c
  414   continue
        izones = int(float(nzones)/sqrt(zrad))
c
  418   continue
c
c  ahk 24-jun-96 forced izones = nzones after 418
        izones = nzones
        lcentr = 2
        ledge = izones + 1
        mzones = izones + 2
c
        mchi = mhyd + mimp + 2
        lelec = mchi - 1
        lion = mchi
c
cl      4.2)    comout
c
  420   continue
        nediti = nedit
        nskpi = nskip
        if (nskpi.le.0) nskpi = 1
c
        lhspec(lelec) = ihteti(1)
        lhspec(lion)  = ihteti(2)
c
c
cl      4.3)    comflg
c
c
        dti    = dtinit * ueit
        dtmini = dtmin  * ueit
        dtmaxi = dtmax  * ueit
        tai    = tinit  * ueit
        tbi    = tai
        thetai = theta
        tmaxi  = tmax   * ueit
cdoc
c
cl      5)      set radius
c
c
c       in this code, there are "nzones" "real" zones, the first of
c       which is indexed "lcentr", and the last of which is indexed "ledge".
c       since boundary j is always the inside boundary of zone j, the
c       plasma edge is boundary "ledge+1".
c       in addition, there are two "dummy" zones, zone 1,
c       which is between boundary 1 and boundary "lcentr",
c       and zone "mzones", which lies between boundaries
c       "mzones" and "mzones+1".
c       hence, xzoni(lcentr) = - xzoni(1), and
c       xbouni(ledge+1) = (xzoni(ledge+1) + xzoni(ledge)) / 2
c
c       also note that xzoni and xbouni are scaled so that the plasma
c       boundary is always 1.0, i.e. xbouni(ledge+1) = 1.0.
cend
c
cl      5.1)    tabular input for radii
c
        rmini = rminor(1) * ueil
        if (radius(2).le.0.0) go to 520
        if (rmini.le.epslon) go to 9050
        zradi = 1.0 / rmini
c
        do 500 j = 1, izones
        i001 = j + lcentr - 1
        xzoni(i001) = radius(j) * ueil * zradi
  500   continue
        go to 540
c
  520   continue
c
cl      5.2)    fitted radii
c
        if (nrfit.gt.irfit.or.nrfit.le.0) go to 9051
c
        do 536 j = lcentr, ledge
        go to (522,524), nrfit
  522   continue
        xzoni(j) = (float(j-lcentr) + 0.5) / float(izones)
        go to 536
  524   continue
        xzoni(j) = (float(j-lcentr) + .75) /
     1                          sqrt(float((j-lcentr+1)*izones))
        go to 536
  536   continue
  540   continue
c
cl      5.3)    set xbouni, gx, rmaji, gx, dx2i, dx2inv
c
        xzoni(1) = - xzoni(lcentr)
        xzoni(ledge+1) = 2.0 - xzoni(ledge)
c
        do 542 j = 2, mzones
        xbouni(j) = (xzoni(j) + xzoni(j-1)) * 0.5
  542   continue
c
        xbouni(1) = - xbouni(2*lcentr-1)
        xbouni(mzones+1) = 2.0 - xbouni(mzones-1)
c
c
        do 544 j = 1, mzones
        gx(j) = 1.0
  544   continue
c
        rmaji = rmajor(1) * ueil
        rmji(1) = rmaji
c
        do 546 j = lcentr, mzones
        dx2i(j) = (xbouni(j+1)**2 - xbouni(j)**2) * 0.5
        dx2inv(j) = 1.0 / dx2i(j)
  546   continue
c
        dx2i(1) = dx2i(lcentr)
        dx2inv(1) = dx2inv(lcentr)
c
cl      5.4)    set dxboui and dxzoni
c
        do 562 j = 2, mzones
        dxzoni(j-1) = xbouni(j) - xbouni(j-1)
        dxboui(j) = xzoni(j) - xzoni(j-1)
  562   continue
c
        dxzoni(mzones) = xbouni(mzones+1) - xbouni(mzones)
        dxboui(1) = 0.0
        dxboui(mzones+1) = 0.0
c
c
cl      5.5)    set bint
c
c
c       integral of b**2 in zone j is
c       bint(1,j)*b(j)**2 + bint(2,j)*b(j)*b(j+1) + bint(3,j)*b(j+1)**2
c       times the volume of the torus.
c
c
        do 578 jz = lcentr, ledge
        if (xbouni(jz+1).le.zcross*xbouni(jz)) go to 574
c
c               long formula -- zone width is large w.r.t. radius
c
c
        z0 = 0.0
        if (xbouni(jz).gt.epslon) z0 = log(xbouni(jz+1)/xbouni(jz))
        zx2 = xbouni(jz+1)**2 - xbouni(jz)**2
        zx4 = xbouni(jz+1)**4 - xbouni(jz)**4
c
        bint(1,jz) = 2.0*xbouni(jz)**2 / zx2**2 *
     1  (0.25*zx4 - zx2*xbouni(jz+1)**2 + xbouni(jz+1)**4 * z0)
c
        bint(2,jz) = 2.0 * xbouni(jz)*xbouni(jz+1) / zx2**2 *
     1  (0.5*zx4 - 2.0*xbouni(jz)**2 * xbouni(jz+1)**2 * z0)
c
        bint(3,jz) = 2.0 * xbouni(jz+1)**2 / zx2**2 *
     1  (0.25*zx4 - zx2*xbouni(jz)**2 + xbouni(jz)**4 * z0)
c
        go to 578
c
c               short formula -- zone is narrow w.r.t. radius
c
  574   continue
        zxp = xbouni(jz+1) + xbouni(jz)
        zxm = xbouni(jz+1) - xbouni(jz)
c
        bint(1,jz) = 2.666666666 * xbouni(jz)**2 * xbouni(jz+1) *
     1          zxm / zxp**2
        bint(2,jz) = 1.333333333 * xbouni(jz) * xbouni(jz+1) *
     1          zxm / zxp
        bint(3,jz) = 2.666666666 * xbouni(jz) * xbouni(jz+1)**2 *
     1          zxm / zxp**2
c
  578   continue
c
c
c
cl      6)      set ion densities
c
 600    continue
c
        zfract = 0.0
c
        do 601 jh=1,mhyd
        if(nnhfit(jh).ne.0) eehfit(jh)=nnhfit(jh)
 601    continue
        if(mimp.eq.0) go to 603
        do 602 ii=1,mimp
        if(nnifit(ii).ne.0) eeifit(ii)=nnifit(ii)
 602    continue
 603    continue
        if(nnfit.ne.0) eefit=nnfit
        if(ntefit.ne.0) eeteft=ntefit
        if(ntifit.ne.0) eetift=ntifit
        if(nbfit.ne.0) eebfit=nbfit
c
c
cl      6.1)    individual input -- hydrogen
c
c
        call resetr(zdens,111,0.0)
c
        do 618 jh = 1, mhyd
        i001 = iphyd(jh)
        if (dengas(i001,1).le.0.0) go to 605
c
c               tabular input
c
        do 604 jz = lcentr, mzones
        i002 = jz + 1 - lcentr
        rhohs(jh,2,jz) = dengas(i001,i002) * uesd
  604   continue
c
        if (rhohs(jh,2,ledge+1).le.epslon)
     1          rhohs(jh,2,ledge+1) = rhohs(jh,2,ledge)
        go to 610
c
c               fitted
c
  605   continue
        if (denga0(i001).le.0.0) go to 615
        zd0 = (denga0(i001) - denga1(i001)) * uesd
        zd1 = denga1(i001) * uesd
c
c..option to initialize time-dependent hydrogen edge density
c
        if ( bdhyde(1,i001) .gt. epslon ) then
          call timint (tinit,zd1,bdtime,20,bdhyde(1,i001),1,1)
          zd0 = (denga0(i001) - zd1 ) * uesd
          zd1 = zd1 * uesd
        endif
c
        do 608 jz = lcentr, ledge
          zx = xzoni(jz) / xzoni(ledge+1)
          rhohs(jh,2,jz) = zd1 +
     1           zd0*(1.0 - zx**eehfit(i001))**ehfit(i001)
  608   continue
c
        rhohs(jh,2,ledge+1) = zd1
c
c               note contribution to total density
c
  610   continue
        do 614 jz = lcentr, mzones
          zdens(jz) = zdens(jz) - rhohs(jh,2,jz)
  614   continue
c
        go to 618
c
c               fraction only
c
  615   continue
        zfract = zfract + fracth(i001)
  618   continue
c
  620   continue
c
cl      6.2)    individual input -- impurity
c
        if (mimp.le.0) go to 640
c
        do 638 ji = 1, mimp
        i001 = ipimp(ji)
        if (denimp(i001,1).le.0.0) go to 625
c
c               tabular input
c
        do 624 jz = lcentr, mzones
          i002 = jz + 1 - lcentr
          rhois(ji,2,jz) = denimp(i001,i002) * uesd
  624   continue
c
        if (rhois(ji,2,ledge+1).le.epslon)
     1          rhois(ji,2,ledge+1) = rhois(ji,2,ledge)
        go to 630
c
c               fitted
c
  625   continue
        if (denim0(i001).le.0.0) go to 635
        zd0 = (denim0(i001) - denim1(i001)) * uesd
        zd1 = denim1(i001) * uesd
c
c..option to initialize time-dependent impurity edge density
c
        if ( bdimpe(1,i001) .gt. epslon ) then
          call timint (tinit,zd1,bdtime,20,bdimpe(1,i001),1,1)
          zd0 = (denim0(i001) - zd1 ) * uesd
          zd1 = zd1 * uesd
        endif
c
        do 628 jz = lcentr, ledge
          zx = xzoni(jz) / xzoni(ledge+1)
          rhois(ji,2,jz) = zd1 +
     1           zd0*(1.0 - zx**eeifit(i001))**eifit(i001)
  628   continue
c
        rhois(ji,2,ledge+1) = zd1
c
c               note contribution to total density
c
  630   continue
        do 634 jz = lcentr, mzones
          zdens(jz) = zdens(jz) - rhois(ji,2,jz)
  634   continue
c
        go to 638
c
c               fraction only
c
  635   continue
        zfract = zfract + fracti(i001)
  638   continue
c
  640   continue
c
c
cl      6.3)    total density minus components individually specified
c
c
        if (zfract.le.epslon) go to 660
c
        zfract = 1.0 / zfract
c
        if (dens(1).le.0.0) go to 645
c
c               tabular input
c
        do 644 jz = lcentr, mzones
          i002 = jz+1-lcentr
          zdens(jz) = zdens(jz) + dens(i002)*uesd
          if (zdens(jz).le.epslon.and.jz.lt.mzones) go to 9061
  644   continue
c
        i002 = ledge+1-lcentr
        if (dens(i002+1).le.0.0) zdens(mzones) = zdens(mzones) +
     1          (dens(i002) - dens(i002+1))*uesd
        if (zdens(mzones).le.epslon) go to 9061
        go to 650
c
c               fitted density
c
  645   continue
c
        zd0 = (dens0 - dens1)*uesd
        zd1 = dens1*uesd
c
        do 648 jz = lcentr, ledge
          zx = xzoni(jz) / xzoni(ledge+1)
          zdens(jz) = zdens(jz) + zd1 +
     1           zd0*(1.0 - zx**eefit)**efit
          if (zdens(jz).le.epslon) go to 9061
  648   continue
c
        zdens(mzones) = zdens(mzones) + zd1
        if (zdens(mzones).le.epslon) go to 9061
c
  650   continue
c
c
c               set species densities
c
c
        do 654 jh = 1, mhyd
        i001 = iphyd(jh)
        if (dengas(i001,1).gt.0.0.or.denga0(i001).gt.0.0)
     1                                                  go to 654
c
        do 652 jz = lcentr, mzones
        rhohs(jh,2,jz) = zdens(jz) * fracth(i001) * zfract
  652   continue
  654   continue
c
c
        if (mimp.le.0) go to 660
c
        do 658 ji = 1, mimp
        i001 = ipimp(ji)
        if (denimp(i001,1).gt.0.0.or.denim0(i001).gt.0.0)
     1                                                  go to 658
c
        do 656 jz = lcentr, mzones
        rhois(ji,2,jz) = zdens(jz) * fracti(i001) * zfract
  656   continue
  658   continue
c
  660   continue
c
c               set center values
c
        do 664 jh = 1, mhyd
        rhohs(jh,2,1) = rhohs(jh,2,lcentr)
  664   continue
c
        if (mimp.le.0) go to 670
c
        do 668 ji = 1, mimp
        rhois(ji,2,1) = rhois(ji,2,lcentr)
  668   continue
c
  670   continue
c
c
c              average hydrogen mass
c
c
        z0 = 0.0
        z1 = 0.0
c
        do 674 jz=lcentr,mzones
        do 672 jh=1,mhyd
          z0 = z0 + rhohs(jh,2,jz)
          z1 = z1 + rhohs(jh,2,jz) * aspec(jh)
  672   continue
c
c   les  nov-90  d3he --include protons in "hydrogen"
c
      if ( cfutz(490) .gt. epslon ) then
        z1=z1+rhois(lprotn-lhydn,2,jz)*aspec(lprotn)
        z0=z0+rhois(lprotn-lhydn,2,jz)
      endif
c
        ahmean(2,jz) = z1 / z0
c
  674   continue
c
c
c
cl      7)      set te and ti
c
c
c       electron energy density is computed assuming all impurities
c       have an average z of 1, i.e. ion density=electron density
c
cl      7.1)    tabular te input
c
        if (te(1).le.0.0) go to 720
        do 704 j = lcentr, mzones
        i002 = j+1-lcentr
        tes(2,j) = te(i002) * uesh
  704   continue
c
        if (tes(2,mzones).le.0.0) tes(2,mzones) = tes(2,ledge)
        go to 740
  720   continue
c
cl      7.2)    fitted te
c
        zte1 = te1 * uesh
        zte0 = (te0 - te1) * uesh
c
c..option to initialize time-dependent edge electron temperature
c
        if ( bdtee(1) .gt. epslon ) then
          call timint (tinit,zte1,bdtime,20,bdtee(1),1,1)
          zte0 = (te0 - zte1 ) * uesh
          zte1 = zte1 * uesh
        endif
c
        do 724 j = lcentr, ledge
          zx = xzoni(j) / xzoni(ledge+1)
          tes(2,j) = (zte1 +
     1           zte0 * max((1.0 - zx**eeteft)**etefit,0.0))
  724   continue
        tes(2,mzones) = zte1
  740   continue
        tes(2,1) = tes(2,lcentr)
c
cl      7.3)    tabular ti input
c
        if (ti(1).le.0.0) go to 760
        do 744 j = lcentr, mzones
          i002 = j+1-lcentr
          tis(2,j) = ti(i002) * uesh
  744   continue
c
        if (tis(2,mzones).le.0.0) tis(2,mzones) = tis(2,ledge)
        go to 780
  760   continue
c
cl      7.4)    fitted ti
c
        zti1 = ti1 * uesh
        zti0 = (ti0 - ti1) * uesh
c
c..option to initialize time-dependent edge ion temperature
c
        if ( bdtie(1) .gt. epslon ) then
          call timint (tinit,zti1,bdtime,20,bdtie(1),1,1)
          zti0 = (ti0 - zti1 ) * uesh
          zti1 = zti1 * uesh
        endif
c
        do 764 j = lcentr, ledge
          zx = xzoni(j) / xzoni(ledge+1)
          tis(2,j) = (zti1 +
     1           zti0 * max((1.0 - zx**eetift)**etifit,0.0))
  764   continue
        tis(2,mzones) = zti1
  780   continue
        tis(2,1) = tis(2,lcentr)
c
c
cl      8)      set b-poloidal and bz
c
c
cl      8.1)    set bedgi and bzi
c
        bzi = bz * ueib
c
        if (bpoid(1).le.0.0) go to 810
        bedgi = bpoid(1) * ueib
        go to 816
c
  810   continue
        if((nadump(1).gt.lcentr).and.(cfutz(127).lt.1.e18)) go to 812
        z0 = ueii / rmini
        bedgi = emu0 * curent(1) * z0 / twopi
        go to 816
  812   continue
c
c       restrict divertor current to be less than zjs(j0)
c
c       set zone index for the divertor edge
        idedge=nadump(1)
c 
c  15.01 Bug fix: add proton mass zpmass to factors going into zvs
c
        zpmass = fcau * (10.0**fxau)
c
c       calculate the sheath-limited currents
c
        do 811 j0=idedge,ledge
c
c  inverse average ion mass
c
        z0=0.0
        z1=0.0
        do 8102 jh=1,mhyd
        z0=z0+rhohs(jh,2,j0)
        z1=z1+rhohs(jh,2,j0)*aspec(jh)
 8102   continue
        if(mimp.le.0) go to 8106
        do 8104 ji=1,mimp
        ii=ji+limp1-1
        z0=z0+rhois(ji,2,j0)
        z1=z1+rhois(ji,2,j0)*aspec(ii)
 8104   continue
 8106   continue
        ziionm=z0/(z1*zpmass)
c
c
c       compute parallel flow velocity zvs into divertor limiter
c
        zvs=sqrt((tis(2,j0)+tes(2,j0))*ziionm)
c
c  possible charge-exchange "friction" not included at this point
c  because neutrals have not yet been initialized
c
c  at this point assume unit charge on each impurity ion
c
        zjs(j0)=0.25*z0*zvs*fces*(10.0**fxes)
c
  811   continue
c
c       print out sheath-limited current densities
c
        call rarray('  zjs(j)',zjs,55)
c
c       calculate current in the divertor region
c       (current density set equal to .5*zjs(j0))
c
        zdcurr=0.0
c
c  15.01 Use rmini here since rmins is not defined yet
c
        zxsec=2.*fcpi * (rmini*uisl)**2
        do 814 j0=idedge,ledge
        zdcurr=zdcurr+(.5*zjs(j0))*zxsec*dx2i(j0)
  814   continue   
c
c  15.01 Let emu0 multiply both curent and zdcurr
c
        bedgi=(emu0*(curent(1)*ueii-zdcurr*usii))
     &        / (twopi*rmini*xbouni(idedge))
c
  816   continue
c
cl      8.2)    integrated b-poloidal: eebfit = 0
c
c       solve r * b-poloidal = a * (integral from 0 to r(1 /resistivity) r*dr)
c        = a' * integral from 0 to r of te**(3/2) (r*dr),
c       where a and a' are constants chosen so that b-poloidal at the boundary
c       is the right value.
c
c       z0 is used to normalize the integral to prevent overflows
c
        if (eebfit.ne.0.0) go to 840
        zi = 0.0
        if(tes(2,lcentr).le.epslon) go to 9080
        z0 = 1.0 / tes(2,lcentr)
        bpols(1,lcentr) = 0.0
        jsx=mzones
        if((nadump(1).gt.lcentr).and.(cfutz(127).lt.1.e18)) jsx=idedge-1
        do 824 j = lcentr, jsx
        zi = zi + dx2i(j) * sqrt((z0 * tes(2,j))**3)
        bpols(1,j+1) = zi / xbouni(j+1)
  824   continue
c
        jsy=ledge+1
        if((nadump(1).gt.lcentr).and.(cfutz(127).lt.1.e18)) jsy=idedge
        z0 = bedgi * uisb / bpols(1,jsy)
c
        do 826 j = lcentr, jsx
        bpols(1,j+1) = bpols(1,j+1) * z0
  826   continue
        if((nadump(1).le.lcentr).or.(cfutz(127).gt.1.e18)) go to 850
c
c       calculate b-poloidal in the divertor region
c
        ztcurr=curent(1)*ueii-zdcurr*usii
        do 830 j0=idedge,ledge
        ztcurr=ztcurr+(0.5*zjs(j0))*zxsec*dx2i(j0)*usii
c
c  15.01 Add emu0 and twopi factors
c
        bpols(1,j0+1)=(emu0*ztcurr/(twopi*rmini*xbouni(j0+1)))*uisb
  830   continue
c
        go to 850
c
  840   continue
        if(eebfit.lt.0.) go to 8200
c
cl      8.3)    fitted b-poloidal: eebfit > 0
c
c       jz is proportional to (1 - (r/a)**eebfit)**ebfit.
c       to integrate
c
c          b(r) = (1/r)*integral from 0 to r t*jz(t)dt
c
c       expand the integrand using the binomial expansion
c       and integrate term by term:
c
        icentr=lcentr+1
        do 841  j = icentr, mzones
        zsum = xbouni(j)/2.
        zsgn = 1.
        zx = xbouni(j)**eebfit
        zxx = xbouni(j)
        zcoeff = 1.
        k = 1
c
c       summate:
c
 846    zsgn = -zsgn
        za = 1./(k*eebfit + 2)
        zcoeff = zcoeff*(ebfit-k+1)/k
        zxx = zxx*zx
        zt = zsgn*za*zcoeff*zxx
        zsum = zsum + zsgn*za*zcoeff*zxx
        if(abs(zt).le.0.0001) go to 844
        k = k + 1
        if(k.ge.51) go to 844
        go to 846
 844    continue
        bpols(1,j) = zsum
 841    continue
c
c       normalize so that bpols(1,mzones) = bedgi*uisb
c
        do 847 j = icentr,mzones
        bpols(1,j) = bpols(1,j)*bedgi*uisb/bpols(1,mzones)
 847    continue
        go to 8299
c
c       set up bpols from bessel function
c
 8200   continue
        iroot = max(1,min(5,int(-eebfit)))
        do 8210 jz = lcentr+1,mzones
        zx = zroot1(iroot) * xbouni(jz)
        bpols(1,jz)=bedgi*uisb*(xbouni(jz) + ebfit*fbes(1,zx))
 8210   continue
c
 8299   continue
 850    continue
        bpols(1,1) = bpols(1,lcentr + 1)
        bpols(1,lcentr) = 0.
        bpols(1,mzones+1) = bpols(1,mzones)/xbouni(mzones+1)
c
c
c
c
cl      9)      coefficients for transport coefficients
c
c
        ev50s = evs * 50.0
        ctions = 0.75 * sqrt(fcau / fcpi) / fces**4 *
     1          10.0**(0.5*fxnucl - 4.0*fxes)
        cteles = ctions * sqrt(fcae * 0.5) * 10.0**(0.5 * fxae)
c
        cnuhyd = sqrt(fcau) * 10.0**(0.5 * fxau)
        cnuel = cnuhyd * sqrt(fcae) * 10.0**(0.5 * fxae)
c
        ceta = fcme / fces**2 * 10.0**(fxme - 2.0*fxes)
c
        cdnhs = 8.0 * sqrt(2.0 * fcpi * fcme) * (fces * fcc)**2 / 3.0
     1          * 10.0**(2.0*(fxes + fxc) + 0.5*fxme)
        cdnhis = cdnhs / sqrt(fcae) * 10.0**(-0.5 * fxae)
        cdniis = cdnhis
        cditis = cdnhis * sqrt(0.5)
        cdetes = cdnhs
        ccnu = 3.0 * fcae * 10.0**(fxae)
        cdbohm = 0.0625 * fcc / fces * 10.0**(fxc - fxes)
c
c
c
c
cl      10)     beam variables
c               see also subroutine "beams"
c
c
c
 1000   continue
        if (nlres) go to 1020
c
c               set hzbeam, habeam, nhbeam, lhbeam
c
        do 1018 jb = 1, mxhbem
        if (nhbeam(jb).ne.0) go to 1002
c
c               no species--no injector
c
        hton(jb) = epsinv
        htoff(jb) = 0.0
        go to 1018
c
 1002   continue
        mhbeam = jb
        if (nhbeam(jb).lt.0) go to 1008
        hzbeam(jb) = nhbeam(jb)
        if (nhbeam(jb).gt.2) go to 1012
c
c               h or he -- decide which isotope
c
        i001 = nhbeam(jb)
        if (habeam(jb).le.0.0) habeam(jb) = zadef(i001)
        if (nhbeam(jb).gt.1) go to 1006
c
c               some isotope of hydrogen
c
        nhbeam(jb) = -1
        if (habeam(jb).gt.1.5) nhbeam(jb) = -2
        if (habeam(jb).gt.2.5) nhbeam(jb) = -3
        go to 1008
c
c               some isotope of helium
c
 1006   continue
        nhbeam(jb) = -4
        if (habeam(jb).lt.3.5) habeam(jb) = -5
c
c               explicit specifiation of isotope
c
 1008   continue
        if (nhbeam(jb).lt.-5) go to 9100
        i001 = - nhbeam(jb)
        if (habeam(jb).le.0.0) habeam(jb) = zadef2(i001)
        hzbeam(jb) = float(izdef2(i001))
c
c               find chi index (if any)
c
 1012   continue
        lhbeam(jb) = 0
c
        do 1014 js = lhyd1, limpn
        if (nspec(js).ne.nhbeam(jb)) go to 1014
        if (js.le.lhydn) lhbeam(jb) = js
        if (habeam(jb).le.0.0) habeam(jb) = aspec(jb)
        go to 1016
 1014   continue
cap
        call mesage(' *** beam species does not exist in plasma,    ')
        call mesage('     deposited ions will be discarded          ')
        ijb = jb
        call ivar(8hbeam no.,ijb)
        if (habeam(jb).gt.0.0) go to 1016
c
c               no atomic weight
c
cap
        call mesage(' *** no atomic weight for beam species ')
        call ivar(8hbeam no.,ijb)
        call mesage('     crude estimate used ')
        habeam(jb) = 2.0 * hzbeam(jb)
        call rvar(8ha used  ,habeam(jb))
 1016   continue
c
 1018   continue
c
c
cl              check for "no beam" conditions
c
c
 1020   continue
c
        do 1028 jb = 1, mxhbem
        if (nhbeam(jb).eq.0) go to 1028
        if (hton(jb)*rndup.gt.htoff(jb)) go to 1026
        if (hton(jb).ge.epsinv) go to 1026
        if (hebeam(jb).le.epslon) go to 1026
        if (hibeam(jb).le.epslon .and. hpowmw(jb).le.epslon) go to 1026
c
        zf = 0.0
        do 1024 je = 1, mxhfr
        zf = zf + hfract(je,jb)
 1024   continue
c
        if (zf.le.epslon) go to 1026
        go to 1028
c
c               no beam from this injector
c
 1026   continue
        hton(jb) = epsinv
        htoff(jb) = 0.0
 1028   continue
c
cl              dynamic default values --
c               hrmaj = rmajor(1), hrmin = rminor(1)
c
c
        do 1034 jb = 1, mxhbem
        if (htoff(jb).ge.epsinv) go to 1034
        if (hrmaj(jb).le.0.0) hrmaj(jb) = rmajor(1)
        if (hrmin(jb).le.0.0) hrmin(jb) = rminor(1)
 1034   continue
c
c
c
cl              initialize "beams"
c
        call beams(1)
c
        if (.not.nlres) return
c
c
cl      20)     other restart variables
c
c
 2000   continue
c
c
cl              comout
c
c
        nediti = nedit
        nskpi = nskip
        if (nskpi.le.0) nskpi = 1
c
c
cl              comflg
c
c
        dtmini = dtmin * ueit
        dtmaxi = dtmax * ueit
        nitmax = 25
        thetai = theta
        tmaxi = tmax  * ueit
c
        bzi = bz * ueib
        bzs = bz * uesb
c
c
c
cl              initialize impurity treatment
c
c
        call imprad(1)
c
        return
c
c
cl      90)     errors
c
c       note: these errors should be caught by errchk
c
 9020   continue
        call error_olymp(0,iclass,isub,2,
     1          'number of particle species + 2 .gt. mxchi ')
        return
c
 9021   continue
        call error_olymp(0,iclass,isub,2,
     1          ' *** error *** invalid gas type ')
        return
c
 9023   continue
        call error_olymp(0,iclass,isub,2,
     1          ' *** error *** invalid impurity type ')
        return
c
 9050   continue
        call error_olymp(0,iclass,isub,5,
     1          ' *** error *** minor radius negative or almost 0')
        return
c
 9051   continue
        call error_olymp(0,iclass,isub,5,
     1          ' *** error *** nrfit too large or less than 0 ')
        return
c
 9060   continue
        call error_olymp(1,iclass,isub,6,
     1          ' *** error *** sum of partial densities near 0 ')
        call error_olymp(2,zfract,1,0,'ion dens')
        return
c
 9061   continue
        call error_olymp(1,iclass,isub,6,
     1          'ion dens minus indiv. specified dens. .le.0 ')
        call error_olymp(2,zdens,5,mzones,'density ')
        return
c
 9080   continue
        call error_olymp(0,iclass,isub,8,
     1          ' *** error *** electron energy density near 0 ')
        return
c
 9100   continue
        call error_olymp(0,iclass,isub,10,
     1          ' *** error *** invalid gas type in beam ')
        return
c
 9120   continue
        call error_olymp(0,iclass,isub,3,
     1          ' *** error *** negative neutral impurity influx ')
        return
c
        end












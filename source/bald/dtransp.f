c  19:30 06-May-95 .../baldur/code/bald/dtransp.f
c  BALDUR  file DTRANSP   Bateman, Stotler PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c       Contents of file DTRANSP:
c  TRCOEF - compute transport coefficients, various models
c  XSCALE - gradient scale lengths and smoothing
c  SMOOTH - called by XSCALE to do smoothing
c  FILTER - binomial filter routine used in smoothing
c
c       Subroutines moved to file DNEOCL.TEX:
c  TRNEO1 - simple neoclassical transport coefficients:
c  NCFLUX - Hawryluk-Hirshman neoclassical particle transport
c
c       Subroutines moved to file DRIPPLE:
c  TRIPL1 - compute ripple-induced transport coefficients
c  AVERGE - integrations with respect to poloidal angle for ripple comp
c  INTERP - interpolation on rectangular grid of TF ripple amplitudes
c
c**********************************************************************c
c----------------------------------------------------------------------c
c@trcoef  .../baldur/bald/dtransp.f
c rgb 12-jul-01 removed comtrp
c rgb 27-jun-01 use lneocl(1) = 1 to call nclass_int
c   lneocl(1) replaces cfutz(281) as neoclassical model switch
c rgb 17-jun-01 Added convective velocities
c   veltis, veltes, velnis, velnhs
c   changed from dneo1 to dnneo1(j1,j2,jz)
c rap 29-may-01 Normalization of xineo1 is changed
c rap 02-may-00 Dimension of dneo1 is increased to account different 
c               contribution from different hydrogenic atoms
c rap 23-feb-00 call ifix(...) changed to call int(...)
c rap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c rgb 02-jun-96 removed if(.not.inital) go to 550
c rgb 16-apr-95 change call theory -> call ptheory
c rgb 09-may-94 call prtheory rather than call theory(3)
c rgb 20-mar-94 move common /cmneo1/ to file com/cbaldr.m
c rgb 25-feb-94 call theory(3) each step if lthery(29) < 0
c rgb 13-nov-93 temporarily removed call tripl1 to save space
c rgb 11-oct-93 move the Weiland diffusivity matrix to sbrtn convrt
c   for now until I revise the structure of dnhis, dnihs, ...
c   which have peculiarities of the Rutherford neoclassical transport
c   model built into them.  Sort this out later.
c rgb 11-sep-93 new diffusivity matrix from Weiland model sbrtn etawn6
c rgb 04-feb-93 matrix form of Nordman-Weiland model from sbrtn theory 
c rgb 17-jul-92 added xeneo2(mj) to common
c   electron thermal heat flux due to ion temperature gradient
c  les  2-jan-90  changes for d-3he fusion
c                  add d-3he fast particles to beta
c                  skip banana orbit effects for q<1
c                  ware pinch for ions is veware = electron vware
c                  add subr reorder for ion neoclassical transport to
c                    dd3hfus
c                  no cross terms for hydrogen particle D
c                  want hirshman-hawryluk?
c                  treat protons as a hydrogenic species (check consistency
c                  CMG degraded chi-e and anom pinch in subr cmg, dd3hfus
c                  add anom pinch vxcmg to vxemps (>0 --> inward)
c  les  2-jan-90  change zk(21)->zk(22) and resetr(zk,20,0.0) ->(zk,22,
c                   for new CMG chi-e
c rgb 20.05 28-jul-91 move neoclassical transport to sbrtn TRNEO1
c     move ripple transport calculations to sbrtn TRIPL1
c     call trneo1 and tripl1 just after calls to sbrtn empirc and theory
c rgb 05-may-91 18.88 cthery(25),(26), or (27) .gt. 0.0 turns on ltheor
c     These coefficients are used for the Rebut-Lallia-Watkins model
c rgb 03-dec-90 18.79 changed nxmax=10 to nxmax = 6
c rgb 11-jun-90 18.37 ahalfs(jz,1) to ahalfs(jz,jb) in comp of zdramp
c     make sure zft(jz) = 0.0 at the magnetic axis, zft(jz)=0. after do 1081
c rgb 25-jul-89 remove comment from data statement
c        defining ipdx1,ipdx2,ipdx3,ipdxke
c      removed reference to ikiban, which was never defined
c     ahalfs(jz,j1) changed to ahalfs(jz,1) after zdramp = ...
c  dps 17-oct-88 15.06 make following changes by Singer:
c  elg 11-oct-88 add theory particle and energy fluxes
c  dps 14.03 02-aug-88 add Ware pinch fudge factors mediated by ctrnsp(19)
c       and ctrnsp(20) to do same job on zl13 and zl23 as done on L33 in
c       subroutine getchi; also, fix bug on their zone center values.
c  dps 14.02 20-jun-88 write Ware pinch in terms of loop voltage; remove
c       "go to" that skips Ware pinch if cfutz(19) < epslon.
c  rgb 17-jul-87 inserted  if ( ji2 .eq. ji ) go to 1170 after do 1170
c               to avoid neoclassical imp-imp cross terms on like species
c               Skip loop entirely if mimp .lt. 2.
c       Also:  Skip all cross diffusion terms when cfutz(139) .lt. epslon
c               rather than cfutz(139) .lt. -0.1 as before.
c               Now cfutz(139) will be defaulted to 1.0 rather than 0.0
c               for backward compatability.
c  rgb 4-jun-85 moved sect. 2 "useful quantities" to sbrtn getchi(2)
c  rgb 4-jun-85 if(lmpirc) zk(20)=min(zk(20),cfutz(iblimt)*zdbohm)
c  drm 4-dec-84 use horizontal minor radius to calculate aspect
c       ratio in trapped electron correction to resistivity
c
c***********************************************************************
c
      subroutine trcoef
c
c
cl      2.14    compute transport coefficients
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cbparm.m'
      include 'cd3he.m'
      include 'cmg.m'
c
c..temporary common block for ripple induced transport
c
      common /cmrpl1/ dripl1(mj), dripl2(mj), dripl3(mj), dripl4(mj)
     &  , ykineo(mj), yripl0(mj)
c
c
      logical           inital
c
      dimension
     1   zk(22)       , zahred(2,2)  , yitgrd(mj)   ,
     2   ydbohm(mj)   , ydpsdo(mj)   , ydkp(mj)     ,
     3   ycd(mj)      , ydrift(mj)   , ydtpe(mj)    , ydtpi(mj)    ,
     4   yki6rg(mj)   , ykeneo(mj)   , ykanom(mj)   
!cap     &   aspinv(mj)
c
c
!cap      equivalence   (gx(1),aspinv(1))
c
c-----------------------------------------------------------------------
c
c
      data    iclass,isub,inital /2,14,.true./
c
      data      idhps,idips,ikeps,ikips /1,2,3,4/,
     1          idhbom,idibom,ikebom,ikibom /5,6,7,8/,
     2          idhneo,idineo,ikeneo,ikineo /9,10,11,12/,
     3          iketrp,ikedrf,irdpnd /13,14,15/,
     4          itgrad,iware,isqrt /18,19,20/,
     5          idh,idi,ike,iki /21,22,23,24/,
     6          idhemp,idiemp,izfebp,ikiemp /25,26,27,28/,
     7          ikeem2,ikeem3,ikeemq,ikeem4 /29,30,16,17/,
     8          idhq1,idiq1,ikeq1,ikiq1 /31,32,33,34/,
     9          idhem2,idhem3,idhem4,idiem2,idiem3,idiem4
     1          /35,36,37,38,39,42/
      data    idlimt,ixelmt,ixilmt,iblimt /43,44,52,51/
      data    iplt1,iplt2,iplt3 /45,46,47/,
     1          ipltke,ipltdh,ipltdi /48,49,50/
      data    iply1,iply2,iply3 /345,346,347/,
     1          iplyke,iplydh,iplydi /348,349,350/
      data    ipllmt /351/
      data    incflx,ixdoff /110,139/
      data    irexpn,ikeinr,ikeout,ikeexp /40,57,58,59/
      data    iscrp1,iscrp2,iscrp3,iscrp4 /53,54,55,56/,
     1          idivrt /83/
      data    i6regm,icps6r,ickp6r,iccd6r /91,92,93,94/,
     1          icte6r,icti6r /95,96/,
     2          idhkad,idikad,ikekad /97,98,99/
      data    iki6on /90/
      data  itped1,itped2,itped3,itpex1,itpex2,itpex3,itpea0
     1     /170,171,172,173,174,175,176/
c     where the identifiers "itped1", ---, "itpea0" are available for
c     use in an alternative version of subroutine trcoef.
c
      data    iionxi,ishift /281,282/
c
c             where "iionxi" is used in selecting the
c             formulae for ion-heat conductivity
c             more comments on ion-heat conductivity are given below.
c
      data    impirc,ieware /364,365/
c
c       "cfutz(impirc)" is used in subroutine empirc to adjust the
c       weighting in averages across successive timesteps.
c
c       to restore the treatment of ware pinch on electron-energy flow
c       as used in older versions of the baldur code (prior to 19-may-
c       83), set "cfutz(ieware)=1.0".  default is "cfutz(ieware)=0.0".
c
c       if "cfutz(ispitz)" is .gt. epslon, the spitzer resistivity
c       is multiplied by "cfutz(ispitz)".
c
      data    ispitz /401/
c
c*****note******* in this version dnhhs(ih1,ih2,j)=cfutz(366)* ---, so
c                 that it is necessary to set cfutz(366)=1.0 in order
c                 to restore use of dnhhs(ih1,ih2,j).
c
        data irblmt,irbchi,irbcon,irbalf,irbelp/470,471,472,473,474/
        data igfac,igalfa,igbeta,igamma,igdelt/451, 452, 453, 454, 455/
c
c-----------------------------------------------------------------------
c
c
c
c       *       indices in cfutz of the various futz factors:
c
c               the second and third letters identify
c               the diffusion coefficient to which the futzed number
c               will contribute (d-h, d-i, k-e, k-i),
c               and the rest identify the theory which produced
c               the number to be futzed (pseudo-classical,bohm,neo-
c               classical, ripple-trap, drift-wave, banana-orbit,
c               or trapped-electron.)
c
c      data     idhps,idips,ikeps,ikips /1,2,3,4/,
c     1         idhbom,idibom,ikebom,ikibom /5,6,7,8/,
c     2         idhneo,idineo,ikeneo,ikineo /9,10,11,12/,
c     3         iketrp,ikedrf,irdpnd /13,14,15/,
c     4          itgrad,iware,isqrt /18,19,20/,
c     5         idh,idi,ike,iki /21,22,23,24/,
c     6          idhemp,idiemp,izfebp,ikiemp /25,26,27,28/,
c     7          ikeem2,ikeem3,ikeemq,ikeem4 /29,30,16,17/,
c     8         idhq1,idiq1,ikeq1,ikiq1 /31,32,33,34/,
c     9         idhem2,idhem3,idhem4,idiem2,idiem3,idiem4
c     1         /35,36,37,38,39,42/
 
c               cfutz(izfebp) is the exponent of epsilon**2*beta
c               poloidal factor renormalizing the temperature
c               dependence of q<1 transport
c
c               parameters for limits on diffusion and conduction
c
c      data    idlimt,ixelmt,ixilmt,iblimt /43,44,52,51/
c
c               special empirical parameters for modeling plt
c
c      data    iplt1,iplt2,iplt3 /45,46,47/,
c     1          ipltke,ipltdh,ipltdi /48,49,50/
c
c               parameters for 1978 plt model
c
c      data    iply1,iply2,iply3 /345,346,347/,
c     1          iplyke,iplydh,iplydi /348,349,350/
c
c               soft bohm limit on 1978 and 1980 plt models
c               and on 1983 pdx model
c
c     data    ipllmt /351/
c
c               kaye's empirical pdx model feb-83
c
        data    ipdx1,ipdx2,ipdx3,ipdxke/355,356,357,358/
c
c       note:  cfutz(100) is currently used in subroutine imprad to
c       both activate a calculation of cyclotron losses and set the
c       level of these losses
c
c       incflx is the argument of the cfutz(*)-factor which activates
c       hawryluk-hirshman ion transport.  in subroutine convrt, incflx
c       is also set to 110:
c
c      data    incflx,ixdoff /110,139/
c
c       except for the hawryluk-hirshman neoclassical ion-ion formulae,
c       use of cfutz(ixdoff), where ixdoff=139, allows adjustment of
c       off-diagonal (or cross-ion) diffusion coefficients, namely of
c       the neoclassical effects of one ion species on another:
c        cfutz(ixdoff) <= 0.  this zeroes out all the off-diagonal terms
c        cfutz(ixdoff)  > 0.  this results in multiplying all off-diagon-
c                             al coefficients by cfutz(ixdoff) (default
c                             value is cfutz(ixdoff)=1.).
c
c
c       the neoclassical ion-heat conductivity is calculated by one
c       of three available sets of formulae, the set to be used be-
c       ing selected by means of the switch "cfutz(iionxi)":
c        cfutz(iionxi)=2. causes use of formulae based on the bolton-
c                         ware (1983) analysis.
c        cfutz(iionxi)=1. causes the selection of formulae by chang
c                         and hinton (1982); these include an approx-
c                         imate correction for a shafranov shift.
c        cfutz(iionxi)=0. selects the formula from hinton-hazeltine
c                         (1976). this is the default case.
c
c
c       a special empirical transport model based on plt results
c       (rutherford, feb-80) can be assembled and stored in array
c       element zk(11).  the following cfutz(*)-factors are used:
c        cfutz(idhemp)=the multiplying factor for hydrogenic-ion
c                      diffusivity
c        cfutz(idiemp)=the multiplying factor for impurity-ion
c                      diffusivity
c        cfutz(irdpnd)=the multiplying factor for the radially-de-
c                      pendent term in anomalous particle transport
c        cfutz(irexpn)=the exponent of the radially-dependent term;
c                      i.e., (r/a)**cfutz(irexpn), where "a" denotes
c                      the radius to the plasma edge or separatrix
c        cfutz(ikeinr)=the multiplying factor for electron thermal
c                      conductivity at inner plasma radii
c        cfutz(ikeout)=the multiplying factor for electron thermal
c                      conductivity at outer plasma radii
c        cfutz(ikeexp)=the exponent of the radially-dependent term
c                      in the electron thermal conductivity; i.e.,
c                      (1.0 - (r/a)**2)**cfutz(ikeexp), again "a" de-
c                      notes the radial distance to the plasma edge
c                      or separatrix
c
c note that   idhemp,idiemp,irdpnd /25,26,15/
c      data    irexpn,ikeinr,ikeout,ikeexp /40,57,58,59/
c
c               parameters for prescribing diffusion and conduction
c               in the scrapeoff region (currently idivrt=83)
c
c       namelist input parameters cfutz(120) through cfutz(138) are
c       assigned to prescribe properties of any existing scrapeoff
c       layers; in particular:
c        cfutz(120)=innermost radial index of innermost scrapeoff
c                   region or of innermost separatrix
c        cfutz(121)=innermost radial index of second scrapeoff
c                   layer or second separatrix if such exists;
c                   defaulted to 100
c        cfutz(123)=switch for activating charge-exchange friction
c                   in the scrapeoff layers; 0.=off; 1.=on
c        cfutz(124)=factor adjusting the scrapeoff loss-rate for
c                   hydrogenic ion species; defaulted to 1.0
c        cfutz(125)=factor adjusting the scrapeoff loss-rate for
c                   impurity ion species; defaulted to 1.0
c        cfutz(126)=outermost radial zone subject to normal timestep
c                   controls; defaulted to 100. note that cfutz(ixdel1),
c                   where ixdel1=161, is the on-time (secs) for any
c                   exclusion of certain prescribed radial zones from
c                   the timestep control (xdel/delmax), including that
c                   due to cfutz(126).  see comments in subroutine
c                   resolv.
c        cfutz(127)=average path-length (cm) parallel to b in an
c                   innermost scrapeoff layer (connection length)
c        cfutz(128)=average path-length (cm) parallel to b in a
c                   second scrapeoff layer
c        cfutz(129)=prescribed minimum density (1./cm**3) for any
c                   ionic species in the scrapeoff layers
c        cfutz(130)=minimum value of the electron temperature (ev)
c                   in the scrapeoff layers
c        cfutz(131)=minimum value of the ave ion temperature (ev)
c                   in the scrapeoff layers
c        cfutz(132)=special effective recycling coefficient for all
c                   ions pumped out through the divertor.  defaulted
c                   to grecyc and used in subroutine neugas
c        cfutz(idivrt)=a prescribed fixed value, in terms of *bohm,
c                   for particle diffusion coefficients and, in
c                   terms of *1.5*bohm, for thermal conductivities
c                   in all scrapeoff layers; this overrides all other
c                   determinations of the transport quantities
c
c
c      cfutz indices for sawtooth simulation
       data isawd, isawke, iq1, iq2, iqflag / 402,403,404,405,1/
       data idisdi, idisdo, idiski, idisko / 406,407,408,409/
c
c       isawd   is index of particle diffusion fraction for q=1
c       isawke  is index of electron conductivity fraction for q=1
c       iq1     is index for pressure gradient to trigger q=1
c       iq2     is index for pressure gradient to turn off q=1
c       idisdi  diffusion multiplier inside q=1 during disruption
c       idisdo  outside q=1
c       idiski  chi-e multiplier inside q=1 during disruption
c       idisko  outside q=1
c
c               cfutz(iq1) < (dp/dx)/p
c               cfutz(iq2) > (dp/dx)/p
c               iqflag determines the state of the system
c                1=build-up, 0=disrupted
c
c
c
c
c
c
c       the next four parameters involve modeling a straight-line
c       increase in the transport coefficients starting from some
c       prescribed point inside the plasma and rising to fixed
c       values in the scrapeoff layers (hawryluk, jul-79)
c        cfutz(iscrp1)=the fixed value of the hydrogen-ion diffusion
c                   coefficients (cm**2/sec) in the scrapeoff layers
c        cfutz(iscrp2)=the fixed value of the impurity-ion diffusion
c                   coefficients (cm**2/sec) in the scrapeoff layers
c        cfutz(iscrp3)=the fixed value of the electron thermal
c                   conductivity (cm**2/sec) in the scrapeoff layers
c        cfutz(iscrp4)=the starting point for the increase stated
c                   as a fraction multiplying the innermost radius
c                   of the scrapeoff layers
c
c
c
c       as a holdover from earlier versions of the code, some of the
c       input scrapeoff parameters are handled in terms of the integer
c       array nadump(*)
c
c       as set in subroutine auxval:
c       nadump(1)=0
c       if(cfutz(120).gt.epslon) nadump(1)=int(cfutz(120)+0.1)
c       nadump(3)=100
c       if(cfutz(121).gt.epslon) nadump(3)=int(cfutz(121)+0.1)
c       nadump(4)=0
c       if(cfutz(123).gt.epslon) nadump(4)=int(cfutz(123)+0.1)
c       nadump(5)=1
c       nadump(6)=1
c       nadump(7)=100
c       if(cfutz(126).gt.epslon) nadump(7)=int(cfutz(126)+0.1)
c       nadump(7)=max(nadump(1),nadump(7))
c       nadump(8)=1
c       nadump(10)=0
c       nadump(20)=0
c
c      data    iscrp1,iscrp2,iscrp3,iscrp4 /53,54,55,56/,
c     1          idivrt /83/
c
c               parameters for 6-regime model
c
c      data    i6regm,icps6r,ickp6r,iccd6r /91,92,93,94/,
c     1          icte6r,icti6r /95,96/,
c     2          idhkad,idikad,ikekad /97,98,99/
c
c               index for smoothing parameter set in sub. xscale:
c               isplne /101/
c
c               switch to turn on 6-regime contribution to ki
c
c      data    iki6on /90/
c
c**********************************************************************c
c
      if (.not.nlomt2(isub)) go to 10
      call mesage(' *** 2.14 subroutine trcoef bypassed ')
      return
   10 continue
c
c
      if(cfutz(124).le.epslon) cfutz(124)=1.0  ! scrape-off model
      if(cfutz(125).le.epslon) cfutz(125)=1.0  ! scrape-off model
c
c**********************************************************************c
c
      if(.not.inital) go to 20
c
c     "inital" is set ".false." after subsequent initializations
c
      lmpirc=.false.
      ltheor=.false.
      nxmax = 6
c
c     where "nxmax" is the maximum "n"-dimension of the input arrays
c     "vxsemi(n,ix)", "dxsemi(n,ix)", "xesemi(n)", and "xisemi(n)",
c     which are included in "common/comdf2/".
c
c..the semi-empirical model is turned on (lmpirc = .true.)
c  only if one of the variables
c  xesemi(n), xisemi(n), vxsemi(n,ix), or dxsemi(n,ix) is nonzero
c
          zsemi = 0.
c
      do 14 n=1,nxmax
          zsemi = zsemi + abs(xesemi(n)) + abs(xisemi(n))
        do 14 ix=1,limpn
          zsemi = zsemi + abs(vxsemi(n,ix)) + abs(dxsemi(n,ix))
  14  continue
c
      if (zsemi .gt. epslon) lmpirc = .true.
c
c     where "lmpirc=.true." activates calls to subroutine empirc in
c     which semiempirical contributions to the transport coefficients
c     are assembled.
c
c.. the theoritical model of theory is turned on (ltheor=.true)
c   only if one of the f's is nonzero
c
        zsthe=0.0
      do 15 it=1,3
        zsthe=zsthe+abs(fdr(it))+abs(fdrint)+abs(frm(it))+abs(fig(it))
     &    + abs(fkb(it))+abs(frb(it))+abs(fhf(it))+cthery(24+it)
 15   continue
c
      if(zsthe.gt.epslon) ltheor= .true.
c
c.. this activates subroutine theory
c
   20 continue
c
c---------   end of initialization   -----------------------------------
c
c
c       common blocks and variables modified:
c
c       comdif, comdf2
c
c-----------------------------------------------------------------------
c
c
c       then rhoels is recomputed using the new values of cmean.
c       next, various useful quantities are computed
c       all but q are computed at both boundaries and zones.
c
c       finally, the actual diffusion coefficients themselves are
c       computed.
c
c
c-----------------------------------------------------------------------
c
c
c       the models at present (oct-79) are:
c
c       1 - ntrans=0:  all diffusion and conduction coefficients
c                      whose name begins with "d" are set to zero
c
c       2 - ntrans=2:  includes a family of transport models which
c                      are varied by means of cfutz(*)-factors.  but,
c                      when a variable component of toroidal ripple
c                      is included, the variable-component ripple am-
c                      plitude is prescribed by values assigned to el-
c                      ements ii=11 to ii=20 in the following arrays:
c        tbpoid(ii) -> tcompi(ii)=time points (tbpoid(ii) is in secs.)
c                                 where the following values are set
c        bpoid(ii)  -> rcurri(ii)=del(0)
c        curent(ii) -> redgi(ii) =del(r=a), where "a" corresponds to
c                                 the separatrix
c        rmajor(ii) -> rmji(ii)  =exponential factor (beta) in an em-
c                                 pirical variation of tf ripple with
c                                 poloidal angle "theta", namely,
c                                 exp(-beta*theta*theta)
c
c
c-----------------------------------------------------------------------
c
c       the coefficients generated by trcoef are for the following equations:
c
c       density equations:
c
c       d(rhohs(ih,2,jz))/dt = div ( - flux hyd. ih) + source terms
c               - sink terms - recoms(ih,j)*rhohs(ih,2,jz)
c       d(rhois(ii,2,jz))/dt 
c               = div ( - flux imp. ii) + source terms - sink terms
c
c       energy equations:
c
c       d(elec. energy dens.)/dt = div [ detes(jz)* grad te +
c               denes(jz)* grad elec. dens. + detis(jz)* grad ti
c                       - k * T_e * veltes(jz)
c                       + 3/2 k * T_e * total electron flux ]
c               + ohmic heating term + source terms - sink terms
c               - cnueqs(jz)*k*( T_e - T_i )
c       d(ion energy dens.)/dt = div { 
c                 ditis(jz)* grad ti + dites(jz)* grad te
c               + sum over ih2 [ dinhs(ih2,jz)* grad dens hyd. ih2 ]
c               + sum over ii2 [ dinis(ii2,jz)* grad dens imp. ii2 ]
c               - k * T_i * veltis(jz)
c               + 3/2 k * T_i* total ion flux }
c               + source terms - sink terms
c               + cnueqs(jz) * k * ( T_e - T_i )
c
c       where k is boltzmann's constant (te and ti are in units in which
c               k is 1), ih and ih2 are hydrogen indices, ii and ii2 are
c               impurity species indices.
c       source and sink term variables for hydrogen begin with "sh",
c       for electron energy density, with "we", for ion energy density,
c       with "wi".  at present there are no impurity source or sink terms.
c       jz is the zone or boundary index.
c
c       fluxes are:
c
c     - flux of hydrogen ih =  dnhs(ih,jz) * grad hyd. dens ih
c               + dahs(ih,jz) * grad T_e + dbhs(ih,jz) * grad T_i
c               - velnhs(ih,jz) * hy. dens ih
c               + sum over ii { dnhis(ih,ii,jz) *
c                       [ c2mean(ii,1,jz)*rhois(ii,1,jz)* grad hyd. dens ih
c                       - rhohs(ih,1,jz)* grad (mean z imp ii * dens imp ii) ]
c                       + dbhis(ih,ii,jz)*rhohs(ih,1,jz)*rhois(ii,1,jz)
c                               * grad ti }
c     - flux of impurity ii = dnis(ii,jz) * grad imp. dens ii
c               + dais(ii,jz) * grad T_e + dbis(ii,jz)* grad T_i
c               - velnis(ii,jz) * imp. dens ii
c               - sum over ih { dnihs(ii,ih,jz)*
c                       [ cmean(ii,1,jz)*rhois(ii,1,jz)* grad hyd. dens ih
c                       - rhohs(ih,1,jz)* grad imp. dens ii ]
c                       + dbihs(ii,ih,jz)*rhohs(ih,1,jz)*rhois(ii,1,jz)
c                               * grad ti }
c               + sum over ii2 { dniis(ii,ii2,jz)*
c                       [ c2mean(ii2,1,jz)*rhois(ii2,1,jz)* grad imp. dens. ii
c                       - rhois(ii,1,jz)*cmean(ii,1,jz)*
c                               grad (mean z imp. ii2 * dens imp. ii2) ]
c                       + dbiis(ii,ii2,jz)*rhois(ii,1,jz)*rhois(ii2,1,jz)
c                               * grad ti }
c
c
c
c-----------------------------------------------------------------------
c
c
c
cl      2.)     compute useful quantities   moved to sbrtn getchi(2)
c
c
c..inverse aspect ratio
c  gx(ij) = half-width / major radius to midpoint at zone bndries
c
      do 210 ij=1,mzones+1
        gx(ij) = ahalfs(ij,1) / rmids(ij,1)
  210 continue
c
c       "zsepr0" is a multiplying factor used to obtain "r/a" from
c       r/a = ahalfs(jz,1) / ahalfs(iseprx,1) = zsepr0 * ahalfs(jz,1)
c       where "a" here denotes the halfwidth
c       to the plasma edge or separatrix
c
      iseprx=mzones
      if(nadump(1).gt.lcentr) iseprx=nadump(1)
      zsepr0 = 1.0 / ahalfs(iseprx,1)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
cl      3-8)    transport coefficients
c
c
      i1 = mxhyd * mxzone
      i2 = mximp * mxzone
      i3 = mximp * i1
      i4 = mximp * i2
      i5 = mxhyd * i1
      i6 = i1 + i2
c
      call resetr(dahs,i1,0.0)
      call resetr(dais,i2,0.0)
      call resetr(dbhs,i1,0.0)
      call resetr(dbhis,i3,0.0)
      call resetr(dbihs,i3,0.0)
      call resetr(dbiis,i4,0.0)
      call resetr(dbis,i2,0.0)
c
      call resetr(denes,mxzone,0.0)
      call resetr(detes,mxzone,0.0)
      call resetr(detis,mxzone,0.0)
      call resetr(dites,mxzone,0.0)
      call resetr(dinhs,i1,0.0)
      call resetr(dinis,i2,0.0)
      call resetr(ditis,mxzone,0.0)
c
      call resetr(dnhis,i3,0.0)
      call resetr(dnhs,i1,0.0)
      call resetr(dnihs,i3,0.0)
      call resetr(dnhhs,i5,0.0)
      call resetr(dniis,i4,0.0)
      call resetr(dnis,i2,0.0)
c
      call resetr ( velnhs, i1, 0.0 )
      call resetr ( velnis, i2, 0.0 )
      call resetr ( veltis, mxzone, 0.0 )
      call resetr ( veltes, mxzone, 0.0 )
c
      call resetr(vnwars,mxzone,0.0)
      call resetr(vewars,mxzone,0.0)
      call resetr(vxemps,i6,0.0)
c
c   les  nov-90 d3he
c
      call resetr(vxcmg,mxzone,0.0)
c
c..all diffusion coefs set to zero of ntrans .le. 0
c
      if ( ntrans .le. 0 ) return
c
c..initializations
c
        zrmins=rmins
        zrmajs=rmajs
c
cbate      if(.not.inital) go to 550
c
      inital=.false.
c
c..reduced hydrogen-ion mass in atomic units
c
      do 290 jh=1,mhyd
      do 290 jh2=1,mhyd
        zahred(jh,jh2)=(aspec(jh)*aspec(jh2))/(aspec(jh)+aspec(jh2))
  290 continue
c
c  "zy0" is a possible branch point in an empirical expression
c  for anomalous electron thermal conductivity
c
      zy0=0.0
      if(cfutz(ikeout).le.epslon) go to 300
      if(cfutz(ikeinr).le.epslon) go to 300
      if(cfutz(ikeexp).le.epslon) go to 300
        z0=min(cfutz(ikeinr),cfutz(ikeout))
        z1=1.0/cfutz(ikeexp)
        z2=(z0/cfutz(ikeout))**z1
        zy0=sqrt(1.0-z2)
  300 continue
c
      kscrp1=mzones
      kscrp3=mzones+1
      if(nadump(1).gt.lcentr) kscrp1=nadump(1)
      kscrp2=kscrp1
      if(cfutz(iscrp1).gt.epslon) go to 310
      if(cfutz(iscrp2).gt.epslon) go to 310
      if(cfutz(iscrp3).gt.epslon) go to 310
      cfutz(iscrp4)=0.0
  310 continue
      zdx0=cfutz(iscrp4)*xbouni(kscrp1)
      if((nadump(1).le.lcentr).or.(zdx0.le.epslon)) go to 340
      i1=mzones
  320 continue
      if ( xbouni(i1) .gt. zdx0 ) then
        i1=i1-1
        go to 320
      else
        kscrp2=i1
        kscrp3=kscrp2+1
      endif
  340 continue
      if(kscrp2.lt.kscrp1)
     1  zdx0=1./(xbouni(kscrp1)-xbouni(kscrp2))
c
      zie=1./(fces*(10.0**fxes))
      zrae=sqrt(1.0/fcae)*(10.0**(-0.5*fxae))
      zpmass=fcau*(10.0**fxau)
      zemass=fcme*(10.0**fxme)
      zipmas=1./zpmass
      ziemas=1./zemass
      zibzs=1./bzs
      zketrp=0.0
      zbeta0=0.0
      zms0=0.0
c
  550 continue
c
c---------  end of initialization   -----------------------------------c
c
c**********************************************************************c
c
c       calls to supplementary subroutines xscale and empirc
c
c       subroutine xscale determines the electron-density, electron-
c       temperature, and ion-temperature scale lengths, as well as
c       the thermal-pressure scale length and two forms of magnetic
c       shear.  this subroutine is called whenever any one or more of
c       the following modes is included:
c        (1) a mode involving subroutine empirc in which semiempirical
c            contributions to the transport coefficients (including
c            anomalous pinch) are assembled;
c        (2) the 6-regime model;
c        (3) the ion-temperature-gradient contribution to the ion-
c            thermal conductivity;
c        (4) current-driven drift waves;
c        (5) the tang tpe model (in an alternative version of trcoef).
c
c
c       subroutine empirc (which may use scale-lengths calculated in
c       subroutine xscale) is called if the diagonal transport coeffi-
c       cients [ dnhs(jh,jz), dnis(ji,jz), detes(jz), ditis(jz) ] and/
c       or pinch velocites [ vxemps(ix,jz) ] are enhanced or largely
c       determined by semiempirical contributions.  the contributions
c       from subroutine empirc are added to dnhs(jh,jz), ---, ditis(jz)
c       in the course of executing the main radial do-loop below (i.e.,
c       "do 1549 jz = lcentr, mzones").  in each case use is made of
c       the array element "zk(20)" located under the respective subsec-
c       tion headed by "d-h", "d-i", "x-e", and "x-i".  the generalized
c       pinch velocities "vxemps(ix,jz)" are completely assembled in
c       subroutine empirc except for conventional ware-pinch effects,
c       which are added in below under the subsection entitled "ware
c       pinch".  if semiempirical contributions dominate in the evalu-
c       ation of electron-heat and/or ion-heat conductivity, printed-
c       out columns of "te" and "ti" vs. radius are flagged by "ex".
c
c
      if ( lmpirc .or. ltheor .or. cfutz(i6regm).gt.epslon
     & .or. cfutz(itgrad).gt.epslon .or. cfutz(ikedrf).gt.epslon
     & .or. cfutz(itpex1).gt.epslon .or. lcmg.gt.epslon ) then
        call xscale
        if ( lmpirc ) call empirc(2)
        if ( ltheor ) call ptheory(2)
cbate        if ( ltheor  .and.  lthery(29) .lt. 0 ) call prtheory
cbate        if ( ltheor  .and.  lthery(29) .lt. 0 ) call prethprnt
        if ( ltheor  .and.  lthery(29) .lt. 0 ) call ptheory(3)
c   les  nov-90
        if ( lcmg.gt.epslon ) call cmg(2)
      endif
      if ( .not.lmpirc .and. lcmg.le.epslon )
     &   call resetr(vxemps,i6,0.0)
c
c..simple neoclassical transport models
c  and transport from magnetic ripple
c
      if ( lneocl(1) == 1 ) then 
        call nclass_int
      else 
        call trneo1
      endif	 	
cbate      call tripl1
c
c**********************************************************************c
c
c..get total beta
c
        znelin=0.
        zelec=0.
        zion=0.
        zbeam=0.
        zalpha=0.
      zd3fus=0.
c
        do 576 jzz=lcentr,ledge
          znelin=znelin+rhoels(2,jzz)*dxboui(jzz)
          zelec=zelec+rhoels(2,jzz)*tes(2,jzz)*dx2i(jzz)
          zion=zion+rhoins(2,jzz)*tis(2,jzz)*dx2i(jzz)
          zbeam=zbeam+rhobis(2,jzz)*hebems(jzz)*dx2i(jzz)
          zalpha=zalpha+alphai(jzz)*ealfai(jzz)*dx2i(jzz)
c   les  nov-90
          zd3fus=zd3fus+d3fast(jzz)*dx2i(jzz)
576     continue
        znelin=znelin*zsepr0
c
c..set up for goldston scaling
c
        znt=(zelec+zion)*evsinv*zsepr0**2/1.e16
        zip=(1.+ellipt**2)*5.*rmins**2*bzs/
     1  (2.*ellipt*rmajs*q(ledge+1))
        zgold=cfutz(igfac)*znelin*(rmins/100.)**cfutz(igalfa)*znt**
     1  cfutz(igbeta)/((rmajs/100.)**cfutz(igamma)*(zip/1.e6)**
     1  cfutz(igdelt))
c
c
c.. beta limit model for chi-e  ( tfcd project)
c
        irbuse=0
        if(cfutz(irbchi).lt.0.) go to 579
        if(cfutz(irblmt)*cfutz(irbchi).eq.0.) go to 579
        irbuse=1
        znorm=2.*(8.*fcpi/bzs**2)*zsepr0**2
        zelec=znorm*zelec
        zion=znorm*zion
        zbeam=znorm*zbeam*(2./3.)
        zalpha=zalpha*znorm*(2./3.)*uisd*uise
c   les  nov-90
        zd3fus=zd3fus*znorm*2./3.
        ztbeta=zelec+zion+zbeam+zalpha+zd3fus
c
        zbetal=cfutz(irblmt)*(1.+ellipt**2)*rmins/
     1  (2.*sqrt(ellipt)*rmajs*q(ledge+1))
        zfact=1.
        if(zbetal.ge.0.) zfact=exp((ztbeta/zbetal)**2)
        zrbchi=cfutz(irbchi)*(ahalfs(iseprx,1)/100.)*(1.e6/zip)*
     1  zfact/(ellipt**(cfutz(irbelp)-0.5))
c
c
579     continue
c
c     identify q=1 surface
c
       iqzone=0
       izz=nzones-1
       do 581 jz=lcentr,izz
         if(q(jz).lt.1.)iqzone=jz
  581  continue
c
c     iqzone is index of q=1 surface
c
      if ( iqzone .eq. 0 ) then
        iqflag=0
      else
c
       zpr=2.*(rhoels(1,iqzone-1)*tes(1,iqzone-1)-rhoels(1,iqzone)
     * *tes(1,iqzone))/((rhoels(1,iqzone-1)*tes(1,iqzone-1)
     * +rhoels(1,iqzone)*tes(1,iqzone))*(dxboui(iqzone)+epslon))
c
       if ( iqflag .eq. 0 ) then
c
c  flag is on, transport inhibited for q=1
c
         if ( abs(zpr) .lt. cfutz(iq2) ) iqflag=1
       else
c
c  pressure is tested
c
         if ( abs(zpr) .gt. cfutz(iq1) ) iqflag=0
       endif
      endif
c
c  note: zpr is not used again
c
c
c.. zero the initial separatrix value of electron thermal conductivity
c
      zkes0=0.0
c
c
c
c---------------start of main do-loop over the radial index jz----------
c
c
cpub1
c      write(*,*) 'start of main do-loop over the radial index jz'
c      write(*,*) 'cpub> ,jz,detes(jz),cfutz(ikeneo),xeneo1(jz)'
cpub2 
      do 1540 jz = lcentr, mzones
c
c
cl      3)      useful coefficients
c
c
c..ratio of minor radius to radius
c  of plasma edge or separatrix,
c   namely, "r/a"
c
      zseprx = ahalfs(jz,1) / ahalfs(iseprx,1)
c
c..total hydrogen density
c
c
      znhs = 0.0
      do 590 jh = 1, mhyd
        znhs = znhs + rhohs(jh,1,jz)
  590 continue
c
c
c..inverse electron density
c
      znes=1./rhoels(1,jz)
c
c..ion-gyroradius factor
c
      zgyroi=(fcc*sqrt(2.*ahmean(1,jz)*fcmp)/fces)*
     &  10.**(fxc+.5*fxnucl-fxes)
c
c..powers of the inverse aspect ratio
c
      j0=jz
      if(xbouni(j0).le.epslon) j0=lcentr+1
      zasp12=sqrt(gx(j0))
      zasp22=gx(j0)
      zasp32=zasp12*gx(j0)
      zasp34=sqrt(zasp32)
      zasp42=gx(j0)*gx(j0)
      zasp62=gx(j0)*zasp42
c
c**********************************************************************c
c
c
c       *       particle diffusion coefficients -- some
c               are also used to compute conductivities
c
c
c.. bohm d
c
c
      zdbohm = cdbohm * tes(1,jz) * zibzs
c
c
c.. radially-dependent factor for alternative anomalous
c   particle diffusion (rutherford, feb-80)
c
      if( (cfutz(irdpnd) .gt. epslon) .and. (jz.le.kscrp3) ) then
        z1 = zseprx**cfutz(irexpn)
        zrdpnd = cfutz(irdpnd)*z1+1.0
      else
        zrdpnd=1.0
      endif
c
c
c..pseudo-classical d
c
c
      zdpseu = epsinv
      if (abs(bpols(1,jz)).gt.epslon) zdpseu = gspitz(1,jz) *
     &  cdnhs * rhoels(1,jz) * cloges(1,jz) /(sqrt(tes(1,jz)) *
     &  bpols(1,jz)**2) * (1.0 + tis(1,jz)/tes(1,jz))
c
c..baldur d-ps
c
      if ( 2.0*zdpseu .gt. epsinv ) zdpseu = zdbohm
c
c**********************************************************************c
c
c       *       conductivities not derived from particle diffusion coeffs.
c
c               note -- to get actual conductivities,
c               these coefficients must be multiplied by
c               ne, ni or sum nh (except for "xineo1(jz)").
c
c**********************************************************************c
c
c..k-e due to current-driven drift-waves
c
c
      if ( cfutz(ikedrf) .gt. epslon ) then
c
c  bohm-like term
        zblike=16.0*zdbohm
c  ion sound speed
        zvsdi=sqrt(tes(1,jz)*zipmas/ahmean(1,jz))
c  electron streaming velocity
        zesv=ajzs(1,jz)*znes*zie
c  electron thermal velocity
        zvthe=sqrt(tes(1,jz)*ziemas)
c  velocity ratio
        z0 = zesv/zvthe
        z1 = z0*sqrt(z0)
c  scale-length-shear factor
        if ( jz .gt. lcentr ) then
          z0=abs(slnes(jz))*slbps(jz)
        else
          z0=abs(slnes(jz+1))*slbps(jz+1)
        endif
c
c  c-d drift-wave coefficient
        zkedrf=(z1*zblike*zblike)/(zvsdi*z0)
c
      else
        zkedrf=0.0
      endif
c
c**********************************************************************c
c
c                       k-i due to banana orbits
c
c
      zkiban = 0.0
c
c
c
c..enhancement of transport for q.lt.1 due to mhd
c  (method suggested by p. h. rutherford)
c   as a steady state sawtooth mixing model
c
      if ( q(jz) .lt. 1.0 ) then
        za=q(jz)*q(jz)
        za=(1.-za)/(za+q(jz)+epslon)
        za=za*za
        if(rastar.le.epslon) rastar=(2./3.)*40.0
        epstar=rastar/rmajor(1)
        zeps=gx(jz)/epstar
        zbetap=(rhoels(1,jz)*tes(1,jz)+rhoins(1,jz)*tis(1,jz))/
     &          (epslon+bpols(1,jz)**2/(8.0*fcpi))
 
        zfermi=(1.1*za/(1.+za))*(zbetap*(zeps+epslon)**2)**cfutz(izfebp)
      else
        zfermi=0.
      endif
c
c**********************************************************************c
c
c
cl  kadomtsev model yielding diffusion proportional to
cl    1./[ elec. density ]
c
      if ( cfutz(ikekad) .gt. epslon ) then
c
      if ( nstep .gt. 1 ) go to 900
c
c..this model involves the square of the ratio of the speed
c  of light to the electron plasma frequency, and the square
c  of this ratio is expressed as  zkad0/[ elec. density ]
c  where
c
        zkad0=(fcc*fcc*fcme)*(10.0**(2.0*(fxc-fxes)+fxme))/
     &    (4.0*fcpi*fces*fces)
  900 continue
c
c  electron thermal velocity
c
        zvthe=sqrt(tes(1,jz)*ziemas)
c
c  kadomtsev factor
c
        zkekad=(zkad0*zvthe)/(rhoels(1,jz)*q(jz)*rmajs)
c
      endif
c
c**********************************************************************c
c
cl             six-regime model (rutherford)
c
      if((cfutz(i6regm).gt.epslon).and.(jz.le.kscrp3)) go to 920
        zd6rgm=0.0
        go to 1040
  920 continue
c
c..gyrofrequencies and gyroradii factors
c
      if((nstep.gt.1).or.(jz.gt.lcentr)) go to 930
        zgyrfe=(bzs*fces/(fcc*fcme))*10.**(fxes-fxc-fxme)
        zgyroe=(fcc*sqrt(2.*fcme)/fces)*10.**(fxc+.5*fxme-fxes)
  930 continue
        zgyrfi=((bzs*fces)/(fcc*ahmean(1,jz)*fcmp))*
     &    10.**(fxes-fxc-fxnucl)
c
c..coefficient used in critical densities
c  (includes factors for converting to temperature in "ev")
c
      zcrit0=(gx(jz)*zasp32*evsinv*evsinv)/
     &  (ahalfs(jz,1)*q(jz)*(1.+xzeff(1,jz))*cloges(1,jz)+epslon)
c
c..compute gyroradii
c
c  toroidal and poloidal gyroradii
c
      z0=sqrt(tes(1,jz))
      z1=sqrt(tis(1,jz))
      zrhoet=zgyroe*z0*zibzs
      zrhoep=zgyroe*z0/(bpols(1,jz)+epslon)
      zrhoit=zgyroi*z1*zibzs
      zrhoip=zgyroi*z1/(bpols(1,jz)+epslon)
c
c  collisional-drift regimes
c
c       bohm value = zdbohm
c       pseudoclassical-like value = zdpsdo
c       kadomtsev-pogutse value = zdkp
c
c  electron collision frequency
c
      znueta=gspitz(1,jz)/(telecs(1,jz)+epslon)
c
c  temperature ratio
c
      zratio=1.0+tis(1,jz)/tes(1,jz)
c
c  pseudoclassical-like value
c
      z0=zrhoet/(slbps(jz)+1.e-15)
      zdpsdo=cfutz(icps6r)*znueta*zratio*z0*z0
c
c  kadomtsev-pogutse value
c
      zblike=16.0*zdbohm
      z0=abs(slnes(jz))*abs(slbps(jz))
      z0=z0*z0+epslon
      if ( cfutz(ickp6r) .gt. epslon ) then
        z1=(znueta*zratio*zblike**4)/(z0*zgyrfe*zgyrfi+epslon)
        zdkp=cfutz(ickp6r)*(z1**.33333333333)
      else
        zdkp=0.0
      endif
c
c  combination of all collisional-drift effects
c
      if ( (zdpsdo.gt.epslon) .and. (zdkp.gt.epslon) ) then
        z1 = (zdkp*zdbohm)/(zdkp+zdbohm)
        zdrift = (z1*zdpsdo)/(z1+zdpsdo)
      else
        zdrift=zdpsdo+zdkp
      endif
c
c  additional diffusion due to collisionless drift waves
c
      z1=sqrt(tes(1,jz)/(fcme*z0))*10.**(-.5*fxme)
      zcd=cfutz(iccd6r)*z1*zratio*zrhoet*zrhoet
      zdrift=zdrift+zcd
c
c..trapped-electron effects
c
c  critical density
c
      zncrit=6.0e+13*zcrit0*tes(1,jz)*tes(1,jz)
      z0=rhoels(1,jz)/(zncrit+epslon)
c
      if ( z0 .gt. 13. ) then
c
c..90-degree scattering frequency
c
        znu=(1.0+xzeff(1,jz))/(telecs(1,jz)+epslon)
c
c..trapped-electron diffusion coefficient
c
        z1=slnes(jz)*slnes(jz)+epslon
        z2=(cfutz(icte6r)*zasp32*zblike*zblike)/z1
        z3=1.02*rhoels(2,jz-1)-rhoels(2,jz)
        znu0=gx(jz)*(10.**(-.5*fxnucl))*
     &    sqrt((tes(1,jz)*slbps(jz))/(z1*ahmean(1,jz)*fcmp))
c
        if ( z3 .lt. 0.0 ) then
          zdtpe = 0.1*z2/(znu+5.0*znu0)
        else
          zdtpe = 3.0*z2/(znu+10.0*znu0)
        endif
c
        zdtpe=zdtpe*exp(-z0)
        zdtpe=(zdbohm*zdtpe)/(zdbohm+zdtpe)
c
      else
c
        zdtpe=0.0
c
      endif
c
c     trapped-ion effects
c
c  critical density
c
      zncrit=2.2e+12*zcrit0*tes(1,jz)*tis(1,jz)/
     1  (ahmean(1,jz)**.33333333333)
      z0=rhoels(1,jz)/(zncrit+epslon)
c
      if ( z0 .gt. 13.) then
c
        zdtpi=0.0
c
      else
c
c  trapped-ion diffusion coefficient
c
        z2=(zblike*tis(1,jz))/(tis(1,jz)+tes(1,jz))
        zdtpi=2.0*cfutz(icti6r)*gx(jz)*zasp32*z2*z2/(znu*z1)
        zdtpi=zdtpi*exp(-z0)
        zdtpi=(zdbohm*zdtpi)/(zdbohm+zdtpi)
c
      endif
c
c     ratio of scale lengths
c
      zrnrt=abs(slnes(jz))/(abs(sltes(jz))+1.e-06)
c
c     combined particle-diffusion coefficient
c
      zd6rgm=zrnrt*(1.5*zdtpi+zdtpe)+zdtpi+zdrift
c
c     output arrays
c
      ydbohm(jz)=zdbohm
      ydpsdo(jz)=zdpsdo
      ydkp(jz)=zdkp
      ycd(jz)=zcd
      ydrift(jz)=zdrift
      ydtpe(jz)=zdtpe
      ydtpi(jz)=zdtpi
c
 1040 continue
c
c**********************************************************************c
c
cl      ion-temperature-gradient contribution to the ion thermal
cl      conductivity
c
      if ( cfutz(itgrad) .gt. epslon ) then
c
        if ( cfutz(i6regm) .lt. epslon ) then
c
c  ion gyroradius
c
          zrhoit=zgyroi*zibzs*sqrt(tis(1,jz))
c
        endif
c
c  ion-temperature bohm-like term
c
      zblike=16.0*cdbohm*zibzs*tis(1,jz)
c
c  enhancement factor
c
      z0=1./abs(sltis(jz))
      z1=z0*abs(slnes(jz))
      z2=z1-1.0
      if(z1.gt.2.0) z2=1.0
      if(z1.lt.1.0) z2=0.0
c
c  contribution to the ion thermal conductivity
c
        zitgrd=cfutz(itgrad)*zblike*zrhoit*z0*z1*z2
c
      else
c
        zitgrd=0.0
c
      endif
c
 1070 continue
c
c**********************************************************************c
c
c       kaye's empirical pdx model
c
        zpdx=0.
        if(cfutz(ipdx1).gt.epslon)
     &    zpdx=cfutz(ipdx1) * (ahalfs(jz,2))**cfutz(ipdx2)
     &      *(bpols(2,jz)*1.e-4)**cfutz(ipdx3)
c
c  note:  bpols has been converted to tesla
c
c**********************************************************************c
c
cl              plt empirical modeling as of apr-80
c
      if ( cfutz(iplt1) .gt. epslon ) then
        z0=tes(1,jz)*useh
        zplt=cfutz(iplt1)*znes/(z0**cfutz(iplt3)+epslon)+
     &    cfutz(iplt2)
      else
        zplt=0.0
      endif
c
c               1978 plt modeling
c
        zply=0.0
       if ( cfutz(iply1) .gt. 0. .and. cfutz(iply2) .gt. 0.
     &       .and.  cfutz(iply3) .gt. 0. )
     &  zply = cfutz(iply1)*znes/(1.+(tes(1,jz)/(cfutz(iply2)*uesh))
     &   **cfutz(iply3))**(1./cfutz(iply3))
c
c**********************************************************************c
c
c                       d-h
c
c
      do 1270 jh = 1, mhyd
      if(nadump(1).le.lcentr) go to 1230
      if(cfutz(idivrt).le.epslon) go to 1210
      if(jz.lt.nadump(1)) go to 1230
      dnhs(jh,jz)=cfutz(idivrt)*zdbohm
      go to 1260
 1210 continue
      if(cfutz(iscrp1).le.epslon) go to 1230
      if(jz.lt.kscrp3) go to 1230
      if(jz.le.lcentr) go to 1230
      if(jz.gt.kscrp1) go to 1220
c       zdx0=1./(xbouni(kscrp1)-xbouni(kscrp2)), see 340
      dnhs(jh,jz)=zdx0*(xbouni(jz)-xbouni(kscrp2))*
     1  (cfutz(iscrp1)-dnhs(jh,kscrp2))
      dnhs(jh,jz)=dnhs(jh,jz)+dnhs(jh,kscrp2)      
      go to 1260
 1220 continue
cpub if cfutz(iscrp1) = 1234567, it'll be overided by neocls transport
      if(cfutz(iscrp1).gt.1234566 .and. cfutz(iscrp1).lt.1234568) then
	dnhs(jh,jz)=cfutz(idhneo) * dnneo1(jh, jh, jz)
      else
	dnhs(jh,jz)=cfutz(iscrp1)
      endif
      go to 1260
c
 1230 continue
      call resetr(zk,22,0.0)
      zk(1) = cfutz(idhps)*zdpseu
      zk(2) = cfutz(idhbom) * zdbohm
      zk(3) = cfutz(idhneo) * dnneo1(jh, jh, jz)
      zk(10) = cfutz(idh)
      zk(11)=cfutz(idhemp)*znes*zrdpnd
      zk(8) = cfutz(idhq1) * zdbohm * zfermi
      zk(13)=cfutz(idhem4)*znes
      if((xbouni(jz).gt.epslon).and.(cfutz(idhem2).gt.epslon))
     1  zk(13)=zk(13)+cfutz(idhem2)*(zseprx**cfutz(idhem3))
      zblimt=cfutz(iblimt)*zdbohm
      zblimt=max(zblimt,zdbohm)
      if(zk(11).gt.epslon)
     1  zk(11)=zk(11)*zblimt/(zk(11)+zblimt)
      if(zk(13).gt.epslon)
     1  zk(13)=zk(13)*zblimt/(zk(13)+zblimt)
      if(cfutz(i6regm).gt.epslon) zk(14)=zd6rgm
      if(cfutz(idhkad).gt.epslon) zk(15)=cfutz(idhkad)*zkekad
      if(cfutz(ipltdh).gt.epslon) zk(16)=cfutz(ipltdh)*zplt
      if(cfutz(iplydh).gt.epslon) zk(17)=cfutz(iplydh)*zply
      zpllmt=cfutz(ipllmt)*zdbohm
      zpllmt=max(zpllmt,zdbohm)
      if(zk(16).gt.epslon)
     1  zk(16)=zk(16)*zpllmt/(zk(16)+zpllmt)
      if(zk(17).gt.epslon)
     1  zk(17)=zk(17)*zpllmt/(zk(17)+zpllmt)
c
      if(lmpirc) zk(20) = min( dxemps(jh,jz), cfutz(iblimt) * zdbohm )
c
      if(ltheor)
     &  zk(21) = min( ( dxthes(jh,jz) ) * uisl**2
     &                   ,cfutz(iblimt)*zdbohm)
c
c   les  nov-90
      if(lcmg.gt.epslon) zk(22)=min( dxcmg(jz), cfutz(iblimt)*zdbohm)
c
      zsum = 0.
      do 1240 jk = 1, 22
        zsum = zsum + zk(jk)
 1240 continue
cap
      if ( cfutz(idlimt) .gt. epslon ) then
        if ( lthery(39) .ge. 0 ) then
          zsum = min( zsum, cfutz(idlimt))
        else
          z21 = difthi(2,2,jz) * uisl**2
          zsum = min( zsum+z21, cfutz(idlimt) ) - z21
        endif 
      endif
c
      dnhs(jh,jz) = zsum
 1260 continue
c
c
c     set q=1 transport inhibition
       if(iqflag.eq.1 .and. cfutz(isawd).gt.epslon .and. jz.eq.iqzone)
     & dnhs(jh,jz)=dnhs(jh,jz)*cfutz(isawd)
       if(iqflag.eq.0 .and. cfutz(idisdi).gt.epslon .and. jz.le.iqzone)
     & dnhs(jh,jz)=dnhs(jh,jz)*cfutz(idisdi)
       if(iqflag.eq.0 .and. cfutz(idisdo).gt.epslon .and. jz.gt.iqzone)
     & dnhs(jh,jz)=dnhs(jh,jz)*cfutz(idisdo)
c
c..convective velocity velnhs(jz)
c
ctemp      velnhs(jh,jz) = velnhs(jh,jz)
ctemp     &   - cfutz(idhneo) * vnhneo1(jh,jz)
c
c..off diagonal terms in the hydrogen flux from theory
c
cb      if ( ltheor) then
c
cb        dbhs(jh,jz) = dbhs(jh,jz)
cb     &    + difthi(2,1,jz) * uisl**2 * rhohs(jh,1,jz) / tis(1,jz)
c
cb        dahs(jh,jz) = dahs(jh,jz)
cb     &    + difthi(2,3,jz) * uisl**2 * rhohs(jh,1,jz) / tes(1,jz)
c
cb        do 1262 ji=1,mimp
c
cb          dnhis(jh,ji,jz) = dnhis(jh,ji,jz)
cb     &     + difthi(2,4,jz) * uisl**2 * rhohs(jh,1,jz) / rhois(ji,1,jz)
c
cb 1262   continue
c
cb      endif
c
 1270 continue
c
c**********************************************************************c
c
c                       d-i
c
c
      if (mimp.le.0) go to 1350
c
      do 1340 ji = 1, mimp
       ii=ji+lhydn
      if(nadump(1).le.lcentr) go to 1300
      if(cfutz(idivrt).le.epslon) go to 1280
      if(jz.lt.nadump(1)) go to 1300
      dnis(ji,jz)=cfutz(idivrt)*zdbohm
      go to 1330
 1280 continue
      if(cfutz(iscrp2).le.epslon) go to 1300
      if(jz.lt.kscrp3) go to 1300
      if(jz.le.lcentr) go to 1300
      if(jz.gt.kscrp1) go to 1290
c       zdx0=1./(xbouni(kscrp1)-xbouni(kscrp2)), see 340
      dnis(ji,jz)=zdx0*(xbouni(jz)-xbouni(kscrp2))*
     1  (cfutz(iscrp2)-dnis(ji,kscrp2))
      dnis(ji,jz)=dnis(ji,jz)+dnis(ji,kscrp2)
      go to 1330
 1290 continue
cpub if cfutz(iscrp2) = 1234567, it'll be overided by neocls transport
      if(cfutz(iscrp2).gt.1234566 .and. cfutz(iscrp2).lt.1234568) then
	dnhs(ji,jz)=cfutz(idineo) * dnneo1(ii, ii, jz)
      else 
	dnis(ji,jz)=cfutz(iscrp2)
      endif
      go to 1330
c
 1300 continue
      call resetr(zk,22,0.0)
      zk(1) = cfutz(idips)*zdpseu
      zk(2) = cfutz(idibom) * zdbohm
      zk(3) = cfutz(idineo) * dnneo1(ii, ii, jz)
      zk(10) = cfutz(idi)
      zk(11)=cfutz(idiemp)*znes*zrdpnd
      zk(8) = cfutz(idiq1) * zdbohm * zfermi
      zk(13)=cfutz(idiem4)*znes
      if((xbouni(jz).gt.epslon).and.(cfutz(idiem2).gt.epslon))
     1  zk(13)=zk(13)+cfutz(idiem2)*(zseprx**cfutz(idiem3))
      if(zk(11).gt.epslon)
     1  zk(11)=zk(11)*zblimt/(zk(11)+zblimt)
      if(zk(13).gt.epslon)
     1  zk(13)=zk(13)*zblimt/(zk(13)+zblimt)
      if(cfutz(i6regm).gt.epslon) zk(14)=zd6rgm
      if(cfutz(idikad).gt.epslon) zk(15)=cfutz(idikad)*zkekad
      if(cfutz(ipltdi).gt.epslon) zk(16)=cfutz(ipltdi)*zplt
      if(cfutz(iplydi).gt.epslon) zk(17)=cfutz(iplydi)*zply
      if(zk(16).gt.epslon)
     1  zk(16)=zk(16)*zpllmt/(zk(16)+zpllmt)
      if(zk(17).gt.epslon)
     1  zk(17)=zk(17)*zpllmt/(zk(17)+zpllmt)
c
      if(lmpirc) zk(20) = min(dxemps(ii,jz), cfutz(iblimt) * zdbohm )
c
      if(ltheor)
     &  zk(21) = min( ( dxthes(ii,jz) )*uisl**2
     &                  , cfutz(iblimt)*zdbohm)
c
c   les  nov-90
      if(lcmg.gt.epslon) zk(22) = min(dxcmg(jz),cfutz(iblimt)*zdbohm)
c
      zsum = 0.
      do 1310 jk = 1, 22
        zsum = zsum + zk(jk)
 1310 continue
cap
      if ( cfutz(idlimt) .gt. epslon ) then
        if ( lthery(39) .ge. 0 ) then
          zsum = min( zsum, cfutz(idlimt))
        else
          z21 = difthi(4,4,jz) * uisl**2
          zsum = min( zsum+z21, cfutz(idlimt) ) - z21
        endif 
      endif
c
      dnis(ji,jz) = zsum
 1330 continue
 1340 continue
c
c..convective velocity velnis(jz)
c
ctemp      velnis(ji,jz) = velnis(ji,jz)
ctemp     &   - cfutz(idineo) * vnineo1(ji,jz)
c
c..off diagonal terms in the impurity flux from theory
c
cb      if ( ltheor) then
c
cb        dbis(ji,jz) = dbis(ji,jz)
cb     &    + difthi(4,1,jz) * uisl**2 * rhois(ji,1,jz) / tis(1,jz)
c
cb        dais(ji,jz) = dais(ji,jz)
cb     &    + difthi(4,3,jz) * uisl**2 * rhois(ji,1,jz) / tes(1,jz)
c
cb        do 1342 jh=1,mhyd
c
cb          dnihs(ji,jh,jz) = dnihs(ji,jh,jz)
cb     &     + difthi(4,2,jz) * uisl**2 * rhois(ji,1,jz) / rhohs(jh,1,jz)
c
cb 1342   continue
c
cb      endif
c
 1350 continue
c
c**********************************************************************c
c
c                       k-e
c
c      
c      
      ykeneo(jz)=0.0
      if(nadump(1).le.lcentr) go to 1380
      if(cfutz(idivrt).le.epslon) go to 1360
      if(jz.lt.nadump(1)) go to 1380
      detes(jz)=1.5*cfutz(idivrt)*zdbohm
      go to 1430
 1360 continue
      if(cfutz(iscrp3).le.epslon) go to 1380
      if(jz.lt.kscrp3) go to 1380
      if(jz.le.lcentr) go to 1380
      if(jz.gt.kscrp1) go to 1370
c       zdx0=1./(xbouni(kscrp1)-xbouni(kscrp2)), see 340
      detes(jz)=zdx0*(xbouni(jz)-xbouni(kscrp2))*
     1  (cfutz(iscrp3)-zkes0)
      detes(jz)=detes(jz)+zkes0
      go to 1430
 1370 continue
cpub if cfutz(iscrp3) = 1234567, it'll be overided by neocls transport
      if(cfutz(iscrp3).gt.1234566 .and. cfutz(iscrp3).lt.1234568) then
	detes(jz)=cfutz(ikeneo) * xeneo1(jz)
      else 
	detes(jz)=cfutz(iscrp3)
      endif
      go to 1430
c
 1380 continue
      call resetr(zk,22,0.0)
      zk(1) = cfutz(ikeps) * zdpseu * 1.5
      zk(2) = cfutz(ikebom) * zdbohm * 1.5
      zk(3) = cfutz(ikeneo) * xeneo1(jz)      
c
      ykeneo(jz)=zk(3)
c
      zk(4) = cfutz(iketrp) * zketrp
      zk(5) = cfutz(ikedrf) * zkedrf
      zk(10) = cfutz(ike)
c
      if ( cfutz(ikeout) .gt. epslon ) then
       if ( zseprx .lt. zy0 ) then
        zk(11)=cfutz(ikeinr)*znes/ (-zseprx*zseprx+1.0)**cfutz(ikeexp)
       else
        zk(11)=cfutz(ikeout)*znes
       endif
      endif
c
      zk(8) = cfutz(ikeq1) * zdbohm * 1.5 * zfermi
      zk(13)=cfutz(ikeem4)*znes*q(lcentr)**cfutz(ikeemq)
      if((xbouni(jz).gt.epslon).and.(cfutz(ikeem2).gt.epslon))
     1  zk(13)=zk(13)+cfutz(ikeem2)*(zseprx**cfutz(ikeem3))
c
c       tfcd beta limit model  added to intor scaling
c
        if(zseprx.gt.0.) zshape=(cfutz(irbcon)+zseprx**
     1  cfutz(irbalf))/(cfutz(irbcon)+0.5**cfutz(irbalf))
        if(zseprx.le.0.) zshape=cfutz(irbcon)
     1  /(cfutz(irbcon)+0.5**cfutz(irbalf))
        if(irbuse.eq.1) zk(13)=zk(13)+zrbchi*zshape
        if(cfutz(igfac).gt.0.) zk(13)=zk(13)+zgold*znes
c
      zblimt=1.5*zblimt
      if(zk(11).gt.epslon)
     1  zk(11)=zk(11)*zblimt/(zk(11)+zblimt)
      if(zk(13).gt.epslon)
     1  zk(13)=zk(13)*zblimt/(zk(13)+zblimt)
      if(cfutz(i6regm).gt.epslon) zk(14)=1.5*(5.0*zdtpi+
     1  zdtpe+zdrift)
      if(cfutz(ikekad).gt.epslon) zk(15)=cfutz(ikekad)*zkekad
      if(cfutz(ipltke).gt.epslon) zk(16)=cfutz(ipltke)*zplt
      if(cfutz(ipdxke).gt.epslon) zk(16)=zk(16)+cfutz(ipdxke)*zpdx
      if(cfutz(iplyke).gt.epslon) zk(17)=cfutz(iplyke)*zply
      zpllmt=1.5*zpllmt
      if(zk(16).gt.epslon)
     1  zk(16)=zk(16)*zpllmt/(zk(16)+zpllmt)
      if(zk(17).gt.epslon)
     1  zk(17)=zk(17)*zpllmt/(zk(17)+zpllmt)
c
      if(lmpirc) zk(20) = min(xeemps(jz), cfutz(iblimt) * zdbohm )
c
      if(ltheor) zk(21) = min( ( xethes(jz) )*uisl**2
     &    , cfutz(iblimt) * zdbohm )
c
c   les  nov-90
c
      if(lcmg.gt.epslon) zk(22) = min(xcmg(jz),cfutz(iblimt)*zdbohm)
c
      imax = 0
      zmax = 0.0
      detes(jz) = 0.0
c
      do 1420 jk = 1, 22
        if ( zmax .lt. zk(jk) ) then
          imax = jk
          zmax = zk(jk)
        endif
        detes(jz) = detes(jz) + zk(jk)
 1420 continue
cap
      if ( cfutz(idlimt) .gt. epslon ) then
        if ( lthery(39) .ge. 0 ) then
          detes(jz) = min(detes(jz),cfutz(ixelmt))
          if( abs(detes(jz)-cfutz(ixelmt)) .le. epslon ) imax=10
        else
          z21 = difthi(3,3,jz) * uisl**2
          detes(jz) = min( detes(jz)+z21, cfutz(ixelmt) ) - z21
          if( abs(detes(jz)-cfutz(ixelmt)) .le. epslon ) imax=10
        endif 
      endif
c
 1430 continue
      if ( jz .eq. kscrp2 ) zkes0=detes(jz)
c
c
c
c
c     set q=1 transport inhibition
       if(iqflag.eq.1 .and. cfutz(isawke).gt.epslon .and. jz.eq.iqzone)
     & detes(jz)=detes(jz)*cfutz(isawke)
       if(iqflag.eq.0 .and. cfutz(idiski).gt.epslon .and. jz.le.iqzone)
     & detes(jz)=detes(jz)*cfutz(idiski)
       if(iqflag.eq.0 .and. cfutz(idisko).gt.epslon .and. jz.gt.iqzone)
     & detes(jz)=detes(jz)*cfutz(idisko)
c
c
      ykanom(jz)=detes(jz)-ykeneo(jz)
c
      detepr(jz) = detes(jz)
      detes(jz)  = detes(jz) * rhoels(1,jz)
      lkeflg(jz) = imax
c
c..electron thermal heat flux due to ion temperature gradient
c
      detis(jz) = detis(jz) + xeneo2(jz) * rhoels(1,jz)
c
c..convective velocity veltes(jz)
c
ctemp      veltes(jz) = veltes(jz)
ctemp     &   - cfutz(ikeneo) * veneo1(jz) * rhoels(1, jz)
c
c..cross terms in the electron heat flux
c
cb      if ( ltheor) then
c
cb        detis(jz) = detis(jz)
cb     &    + difthi(3,1,jz) * uisl**2 * rhoels(1,jz)*tes(1,jz)/tis(1,jz)
c
cb        denes(jz) = denes(jz)
cb     &    + difthi(3,2,jz) * uisl**2 * tes(1,jz)
c
c  Note:  This term should be
c    denes(jz) = denes(jz) + uisl**2 * tes(1,jz) *  
c     ( difthi(3,2,jz) * grad n_H / n_H
c     + difthi(3,4,jz) * grad n_Z / n_Z ) * n_e / grad n_e
c
cb      endif
c
c**********************************************************************c
c
c                       k-i
c
c
      ykineo(jz)=0.0
      yitgrd(jz)=0.0
cpub1      
      if(nadump(1).le.lcentr) go to 1440
      if(cfutz(idivrt).le.epslon) go to 1439
      if(jz.lt.nadump(1)) go to 1440
      ditis(jz)=1.5*cfutz(idivrt)*zdbohm*rhoins(1,jz)
      imax=2
      go to 1470
 1439 continue
      if(cfutz(iscrp3).le.epslon) go to 1440
      if(jz.lt.kscrp3) go to 1440
      if(jz.le.lcentr) go to 1440
      if(jz.gt.kscrp1) go to 1438
      ditis(jz)=zdx0*(xbouni(jz)-xbouni(kscrp2))*
     1  (cfutz(iscrp3)-zkes0)
      ditis(jz)=(ditis(jz)+zkes0)*rhoins(1,jz)
      go to 1470
 1438 continue
cpub if cfutz(iscrp3) = 1234567, it'll be overided by neocls transport
      if(cfutz(iscrp3).gt.1234566 .and. cfutz(iscrp3).lt.1234568) then
	ditis(jz)=cfutz(ikineo) * xineo1(jz) * rhoins(1, jz)
      else 
	ditis(jz)=0.5*cfutz(iscrp3)
      endif
      go to 1470
cpub2      

c
 1440 continue
      call resetr(zk,22,0.0)
      zk(1) = cfutz(ikips) * zdpseu * 1.5 * rhoins(1,jz)
      zk(2) = cfutz(ikibom) * zdbohm * 1.5 * rhoins(1,jz)
      zk(3) = cfutz(ikineo) * xineo1(jz) * rhoins(1, jz)
      zk(6) = (dripl2(jz)+dripl3(jz)+dripl4(jz))*rhoins(1,jz)
cbate      zk(7) = cfutz(ikiban) * zkiban * znhs
      zk(10) = cfutz(iki) * rhoins(1,jz)
      zk(11) = cfutz(ikiemp)*znes*rhoins(1,jz)
      yki6rg(jz)=1.5*cfutz(iki6on)*zd6rgm
      if(cfutz(i6regm).gt.epslon) zk(14)=yki6rg(jz)*rhoins(1,jz)
      zk(8) = cfutz(ikiq1) * zdbohm * 1.5 * zfermi * rhoins(1,jz)
      zk(15)=zitgrd*rhoins(1,jz)
c
      if(lmpirc) zk(20) = min( xiemps(jz)*rhoins(1,jz)
     &  , cfutz(iblimt) * zdbohm * rhoins(1,jz) )
c
      if(ltheor) zk(21) = min(
     &  ( xithes(jz) )*uisl**2*rhoins(1,jz)
     &  , cfutz(iblimt) * zdbohm * rhoins(1,jz) )
c
c   les  nov-90
      if(lcmg.gt.epslon) zk(22)=cfutz(496)*xub(jz)
c
      ykineo(jz)=zk(3)
      yitgrd(jz)=zitgrd
c
      imax = 0
      zmax = 0.0
      ditis(jz) = 0.0
c
      do 1460 jk = 1, 22
        if ( zmax .lt. zk(jk) ) then
          imax = jk
          zmax = zk(jk)
        endif
        ditis(jz) = ditis(jz) + zk(jk)
 1460 continue
cap
      z0 = cfutz(ixelmt)
      if ( cfutz(ixilmt) .gt. epslon ) z0=cfutz(ixilmt)
      if ( z0 .gt. epslon ) then
        if ( lthery(39) .ge. 0 ) then
          z0=z0*rhoins(1,jz)
          ditis(jz)=min(ditis(jz),z0)
          if ( ditis(jz) .eq. z0 ) imax=10
        else
          z0=z0*rhoins(1,jz)
          z21 =  difthi(1,1,jz) * uisl**2 * rhoins(1,jz)
          ditis(jz) = min( ditis(jz)+z21, z0 ) - z21
          if ( ditis(jz) .eq. z0 ) imax=10
        endif 
      endif
 1470 continue
c
c
      ditipr(jz) = ditis(jz)/(rhoins(1,jz)+epslon)
      lkiflg(jz) = imax
c
c..convective velocity veltis(jz)
c
      veltis(jz) = veltis(jz)
     &   - cfutz(ikineo) * vineo1(jz)
c
c..cross terms in the ion heat flux
c
cb      if ( ltheor) then
c
cb        dites(jz) = dites(jz)
cb     &    + difthi(1,3,jz) * uisl**2 * rhoins(1,jz)*tis(1,jz)/tes(1,jz)
c
cb        do 1472 jh=1,mhyd
cb          dinhs(jh,jz) = dinhs(jh,jz)
cb     &      + difthi(1,2,jz) * uisl**2 * tis(1,jz)
cb 1472   continue
c
cb        do 1474 jh=1,mimp
cb          dinis(jh,jz) = dinis(jh,jz)
cb     &      + difthi(1,4,jz) * uisl**2 * tis(1,jz)
cb 1474   continue
c
cb      endif
c
c     possible re-adjustments of the composite particle-diffusion
c     coefficients (the diagonal "d's") for each ionic species.
c
c       ion-ion temperature-gradient coefficients
c       (hydrogen-impurities only)
c
      ztis=1.0/(tis(1,jz)+epslon)
c
      do 1480 jh=1,mhyd
      do 1480 ji=1,mimp
        dbhis(jh,ji,jz)=-ztis*dnhis(jh,ji,jz)*
     &    (0.5*c2mean(ji,1,jz)+cmean(ji,1,jz))
        dbihs(ji,jh,jz)=-ztis*dnihs(ji,jh,jz)*
     &    (0.5*cmean(ji,1,jz)+1.)
 1480 continue
c
c       ion-ion temperature-gradient coefficients
c       (impurity-impurity)
c
      do 1510 ji2=1,mimp
      do 1510 ji=1,mimp
        if ( ji .ne. ji2 ) then
          dbiis(ji,ji2,jz)=ztis*dniis(ji,ji2,jz)*
     &      (c2mean(ji2,1,jz)-cmean(ji,1,jz)*cmean(ji2,1,jz))
        else
          dbiis(ji,ji2,jz)=ztis*dniis(ji,ji2,jz)*
     &      (c2mean(ji,1,jz)-cmean(ji,1,jz)*cmean(ji,1,jz))
        endif
 1510 continue
c
c
c
c
c               convection terms
c
c
c    6-regime correction to the electron flux
c
c
      if ( cfutz(i6regm) .gt. epslon ) then
        denes(jz) = denes(jz) + 1.5*tes(1,jz)*((1.0-1.5*zrnrt)*zdtpi)
cbate      else
cbate        denes(jz) =0.0
      endif
c
c
cpub1
c      write(*,*) 'cpub> ',jz,detes(jz),cfutz(ikeneo),xeneo1(jz) 
cpub2
 1540 continue
c
c
c---------------end of main do-loop over the radial index jz------------
c
        return
c**********************************************************************
c
c
      entry trcprt
c
c       edit print-out of principal quantities in six-regime model
c       and ripple model
c
      zzz=uist*1.e+03
      zt=tai*zzz
      zdt=dtoldi*zzz
c
      if((cfutz(i6regm).le.epslon).and.(cfutz(itgrad).le.epslon))
     1  go to 1580
c
      if(cfutz(i6regm).gt.epslon) go to 1560
      ydbohm = 0.0
      ydpsdo = 0.0
      ydkp = 0.0
      ycd = 0.0
      ydrift = 0.0
      ydtpe = 0.0
      ydtpi = 0.0
 1560 continue
c
      lpage = lpage + 1
c
      write (nprint,1610) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
      write(nprint,1620)
      do 1570 jz=lcentr,mzones
      write(nprint,1630) jz,ahalfs(jz,1),slnes(jz),sltes(jz),sltis(jz),
     1  slbps(jz),ydbohm(jz),ydpsdo(jz),ydkp(jz),ycd(jz),ydrift(jz),
     2  ydtpe(jz),ydtpi(jz)
 1570 continue
c
      lpage = lpage + 1
c
      write (nprint,1610) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
      write(nprint,1621)
      do 1571 jz=lcentr,mzones
      write(nprint,1631) jz,ahalfs(jz,1),yki6rg(jz),yitgrd(jz)
 1571 continue
c
c..printout page of ripple model quantities
c
      if(nlpomt(7)) go to 1600
      lpage = lpage + 1
c
      write (nprint,1610) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
      write(nprint,1640)
      do 1590 jz=1,mzones
        ykineo(jz)=ykineo(jz)/(rhoins(1,jz)+epslon)
        write(nprint,1650) ahalfs(jz,1),ykeneo(jz),ykanom(jz),ydtpe(jz),
     1          ykineo(jz),yripl0(jz),dripl1(jz),dripl2(jz),
     2          dripl3(jz),dripl4(jz)
 1590 continue
c
 1600 continue
c
 1580 continue
c
      return
c
 1610   format(1h1,2x,a48,10x,a72/
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.')
c
 1620 format(/29x,' quantities having to do with the six-regime model'
     & //1x,'j',1x,' radius',2x,' slnes(j)',1x,' sltes(j)',1x,
     & ' sltis(j)',1x,' slbps(j)',1x,'  bohm(j)',1x,' pseudo(j)',1x,
     & ' dk-p(j)',1x,'cdless(j)',1x,' drift(j)',1x,'  dtpe(j)',1x,
     & ' dtpi(j)'/)
 1621 format(/29x,' quantities having to do with the six regime model',
     1 ' (continued)'//
     2 1x,'j',1x,' radius',2x,' ki 6rgme',1x,' k itgrad'/)
 1630 format(1x,i2,1x,f6.2,1x,11(1x,1pe9.2))
 1631 format(1x,i2,1x,f6.2,1x,2(1x,1pe9.2))
c
 1640 format(/3x,'radius',5x,'xe neocl',4x,
     1          'xe anoml',4x,'xe trpel',4x,'xi neocl',4x,'fxed rpl',4x,
     2          'vari rpl',4x,'rpl trap',4x,'rpl plat',4x,'ban drft')
 1650 format(3x,f7.2,3x,9(1pe10.3,2x))
c
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@xscale  /11040/baldur/code/bald/dtransp.f
c       dps 27-jan-87 compute derivatives wrt. armins
c       dps 08-jan-87 rewrote xscale
c***********************************************************************
c
      subroutine xscale
c
c       2.17    determine scale lengths and shear
c
c
      include 'cparm.m'
      include 'cbaldr.m'
c
      dimension
     1   zrho(mj),zte(mj),zti(mj),zbponr(mj),zbprsm(mj)
c
c-----------------------------------------------------------------------
c
      data iclass /2/, isub /17/
      if (.not.nlomt2(isub)) go to 10
      call mesage(' *** 2.17 subroutine xscale bypassed ')
      return
   10 continue
c
c-----------------------------------------------------------------------
c
c...Set up local arrays that can be used directly in the
c...argument list of the smoothing routine. Note that Bpol/r <-> zbponr
c...is the variable smoothed for the magnetic shear calculation.
c
      do 20 j=1,mzones
        zrho(j)=rhoels(2,j)
        zte(j)=tes(2,j)
        zti(j)=tis(2,j)
        if (j.ne.lcentr) zbponr(j)=bpols(1,j)/armins(j,1)
   20 continue
c
c...We assume Bp = a*r +b*r**3 near the origin and calculate "a"
c...using Bp at lcentr+1 and lcentr+2. Then, Bp/r = a at lcentr.
c...Since Bp is kept on zone boundaries, an additional step is
c...needed to assign the value at mzones+1.
c
      zx3sq=armins(lcentr+1,1)**2
      zx4sq=armins(lcentr+2,1)**2
      zbponr(lcentr)=(zbponr(lcentr+1)*zx4sq
     1                -zbponr(lcentr+2)*zx3sq)/(zx4sq-zx3sq)
      zbponr(mzones+1)=bpols(1,mzones+1)/armins(mzones+1,1)
c
c...Now call the smoothing routine for each of the variables involved.
c...The first of the explicit arguments indicates grid; only q and Bp
c...are on the boundary grid (kgrid=1). The second indicates if the
c...volume integral of the smoothed quantity is to be conserved
c...(kcons=1); this is not appropriate for q and Bp. The last
c...indicates the symmetry about the origin; all of these
c...variables require even symmetry (ksym=1).
c
      call smooth(zrho,rhosms,mzones,2,1,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(zte,tesms,mzones,2,1,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(zti,tisms,mzones,2,1,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(thrprs,thpsms,mzones,2,1,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(totprs,topsms,mzones,2,1,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(q,qsmth,mzones+1,1,0,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
      call smooth(zbponr,zbprsm,mzones+1,1,0,1,lcentr,dx2i,
     1            smrlow,smlwcy,lsmord)
c
c...Gradient scale lengths are computed for the electron density and
c...the plasma temperatures. Since the gradients go to 0 at lcentr,
c...the scale lengths are infinite. To get the sign of the large
c...number taking the place of infinity correct, determine whether the
c...variable is increasing or decreasing in the immediate vicinity of
c...lcentr. A loop is used in case the zero gradient region extends over
c...more than one zone; it is terminated when all of the signs
c...have been determined.
c
      zsndn=0.
      zsndte=0.
      zsndti=0.
      do 30 j=lcentr+1,mzones
c
        if ((rhosms(j).ne.rhosms(lcentr)).and.(zsndn.eq.0.))
     1     zsndn=(rhosms(j)-rhosms(lcentr))
     2           /abs(rhosms(j)-rhosms(lcentr))
c
        if ((tesms(j).ne.tesms(lcentr)).and.(zsndte.eq.0.))
     1     zsndte=(tesms(j)-tesms(lcentr))
     2            /abs(tesms(j)-tesms(lcentr))
c
        if ((tisms(j).ne.tisms(lcentr)).and.(zsndti.eq.0.))
     1     zsndti=(tisms(j)-tisms(lcentr))
     2            /abs(tisms(j)-tisms(lcentr))
c
        if (zsndn*zsndte*zsndti.ne.0.) go to 40    ! exit loop when done
   30 continue
   40 continue
c
c...Finally, compute radial gradients of the smoothed variables and
c...form the desired derived quantities. All are on zone boundaries.
      do 50 j=lcentr,mzones
c
        zdrs=armins(j,2)-armins(j-1,2)
c
c..."[-d(ln ne)/dr]**-1"
c
        zdnedr=(rhosms(j)-rhosms(j-1))/zdrs
        zrhosm=0.5*(rhosms(j)+rhosms(j-1))
        slnes(j)=-1./(zdnedr/zrhosm+zsndn*epslon)
c
c..."[-d(ln te)/dr]**-1"
c
        zdtedr=(tesms(j)-tesms(j-1))/zdrs
        ztesm=0.5*(tesms(j)+tesms(j-1))
        sltes(j)=-1./(zdtedr/ztesm+zsndte*epslon)
c
c..."[-d(ln ti)/dr]**-1"
c
        zdtidr=(tisms(j)-tisms(j-1))/zdrs
        ztism=0.5*(tisms(j)+tisms(j-1))
        sltis(j)=-1./(zdtidr/ztism+zsndti*epslon)
c
c..."-d(ln Pthermal)/d(ln r)"
c
        zdthp=(thpsms(j)-thpsms(j-1))/zdrs
        zthpsm=0.5*(thpsms(j)+thpsms(j-1))
        slprs(j)=-armins(j,1)*zdthp/zthpsm
c
c..."d(Ptotal)/dr"
c
        slpts(j)=(topsms(j)-topsms(j-1))/zdrs
c
c..."d(ln q)/d(ln r)"
c
        zdr2=armins(j+1,1)-armins(j-1,1)
        zdqdr=(qsmth(j+1)-qsmth(j-1))/zdr2
        shear(j)=armins(j,1)*zdqdr/qsmth(j)
c
c..."abs[ r d/dr(Bp/r) ]/[ Bz d(ln ne)/dr ]" -> magnetic shear
c
        zdbpdr=(zbprsm(j+1)-zbprsm(j-1))/zdr2
        slbps(j)=armins(j,1)*zdbpdr*abs(slnes(j))/bzs
c
   50 continue
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@smooth  /11040/baldur/code/bald/dtransp.f
c  rgb 25-jun-96 return if korder < 1 and abs(prlow*plwcyc) < 1.e-10
c  rgb 28-jan-95 dimension pyin(1) -> dimension pyin(*)
c     dps 09-may-89 15.08 Add kcons=2 option to allow volume integral
c                   to be conserved in a multiplicative sense.
c     dps 6-jan-87
c***********************************************************************
c
      subroutine smooth(pyin,pyout,kmax,kgrid,kcons,ksym,kxsym,pdx,
     1                  prlow,plwcyc,korder)
!cap
      include 'cparm.m'
c
      dimension pyin(*), pyout(*), pdx(*), zpywrk(mj)
c
c...This subroutine uses a binomial filter to smooth the input array, pyin,
c...and return the smoothed values to the array pyout.
c...If kcons is 1, a correction, zcor, is added to the smoothed
c...array so that the volume integral of pyout equals that of pyin
c...This correction is required to have the same sign as the original
c...integral so that the smoothed array is not reduced in amplitude by the
c...correction. If kcons = 2, pyout is scaled multiplicatively to have the
c...same integral as pyin.
c
c...I/O variables:
c      pyin   - array to be smoothed  (input)
c      pyout  - smoothed array (output)
c      kmax   - number of points in pyin and pyout
c      kgrid  - input 1 if pyin is specified on the boundaries of the
c               integration segments (i.e., the elements of pdx),
c               and 2 if pyin is specified at the center of each segment
c      kcons  - input switch to enable conservation of volume integral of py
c               = 1, conserved via addition
c               = 2, conserved via multiplication
c      ksym   - input to set symmetry of pyout about kxsym;
c               if 0, no symmetry; if -1 (+1), use odd (even) symmetry
c      kxsym  - input value of index of grid point about which symmetry
c               is to be imposed
c      pdx    - "length" of the segments used in calculating the integrals
c               of pyin and pyout.
c      prlow  - input argument for subroutine filter
c      plwcyc -   "      "      "      "        "
c      korder -   "      "      "      "        "    ; smoothing can be
c               bypassed by setting korder < 0.
c
c***********************************************************************
c
      do j=kxsym,kmax
        pyout(j)=pyin(j)
      enddo
c
      if ( korder .lt. 1 .and. abs(prlow*plwcyc) .lt. 1.e-10 ) return
c
c...Set lsym (index of grid point used in applying symmetry)
c...for both types of grid.
c
      if (kgrid.eq.1) then
        lsym=kxsym+1
      else if (kgrid.eq.2) then
        lsym=kxsym
c
c...Quit if kgrid not set properly.
      else
        return
      end if
c
c...Calculate number of points to be used in the smoothing and call
c...the binomial filter routine. Note that pyout(kxsym), pyout(kmax-1),
c...and pyout(kmax) are held fixed in this way. The symmetry conditions
c...about kxsym are applied below; the boundary conditions near kmax
c...remain as they were originally by keeping the last two points fixed.
c...The filter routine is not called at all if korder < 0.
c
      if (korder.lt.0) go to 50      ! do no smoothing
      kpts=kmax-kxsym
      call filter(pyout(kxsym),zpywrk,kpts,prlow,plwcyc,korder)
c
c...Calculate integrals for both arrays using the "length" elements pdx.
c...For kgrid= 1 an integration scheme analogous to a trapezoidal rule is
c...appropriate; for kgrid= 2, a rectangular rule is used.
c
      zint=0.
      zintsm=0.
      zx=0.
      if (kgrid.eq.1) then
        do 20 j=kxsym,kmax-2
          zint=zint+0.5*(pyin(j)+pyin(j+1))*pdx(j)
          zintsm=zintsm+0.5*(pyout(j)+pyout(j+1))*pdx(j)
          zx=zx+pdx(j)
   20   continue
      else if (kgrid.eq.2) then
        do 30 j=kxsym,kmax-1
          zint=zint+pyin(j)*pdx(j)
          zintsm=zintsm+pyout(j)*pdx(j)
          zx=zx+pdx(j)
   30   continue
      end if
c
c...Compute and apply the correction.
c
      if (kcons.eq.1) then
c
c...Additive correction; require it to not reduce integral.
c
        zcor = (zint - zintsm) / zx
        if (zcor.gt.0.) then
          do 40 j=kxsym,kmax
            pyout(j) = pyout(j) + zcor
   40     continue
        end if
c
      else if ((kcons.eq.2).and.(zintsm.ne.0.)) then
c
c...Multiplicative correction.
c
        do 45 j=kxsym,kmax
          pyout(j) = pyout(j) * zint / zintsm
   45   continue
      end if
c
   50 continue                       ! come here if no smoothing
c
c...Finally, enforce symmetry conditions; skip if ksym=0 or kxsym < 2.
      if (ksym.eq.0.or.kxsym.lt.2) return
c
c...The values at the points 1 to kxsym-1 are determined according
c...to the values at lsym+kxsym-2 to lsym and ksym. Note that lsym is
c...chosen so that pyout(kxsym-1)=pyout(kxsym+1)... for kgrid=1 and
c...pyout(kxsym-1)=pyout(kxsym)... for kgrid=2.
      do 60 j=1,kxsym-1
        pyout(kxsym-j)=float(ksym)*pyout(lsym+j-1)
   60 continue
c
      return
      end
c***********************************************************************
c@filter  /11040/baldur/code/bald/dtransp.f
c dps 6-jan-87 Use a local variable for KORDER
c
c***********************************************************************
C
         SUBROUTINE FILTER(PY,PYWORK,KPTS,PRLOW,PLWCYC,KORDER)
C
       PARAMETER   (PI=3.1415927)
C
       DIMENSION   PY(*),    PYWORK(*)
C
C---------------------------------------------------------------------
c
CL              1.         NOTES
C
C     DATE: 10/9/86     STEVE SABBAGH (Columbia U.)
C.       THIS SUBROUTINE IS A LOW PASS FILTER CREATED TO ELIMINATE
C     HIGH FREQUENCY NUMERICAL NOISE IN DATA. THE DATA INPUT BY ARRAY
C     PY(KPTS) IS SMOOTHED BY A LOW PASS BINOMIAL FILTER OF ORDER
C     KORDER, IF KORDER.NE.0. IF KORDER.EQ.0, THEN DATA IS SMOOTHED
C     UNTIL AMPLITUDE OF LOWEST FREQUENCY IS REDUCED BY FACTOR PRLOW.
C     THIS ROUTINE WORKS ONLY ON DATA POINTS THAT ARE EQUALLY SPACED.
C     ENDPOINTS ARE FIXED.
C
C     REF: REV. SCI. INSTRUM. 54 (8), AUG 1983.
C
C.    I/O VARIABLES:   PY     - ARRAY OF DATA POINTS
C.                     PYWORK - WORK SPACE
C.                     KPTS   - NUMBER OF DATA POINTS
C.                     PRLOW  - RATIO OF POST-SMOOTHED AMPLITUDE
C.                              OF LOWEST FREQUENCY (BACKGROUND)
C.                              TO PRE-SMOOTHED AMPLITUDE
C.                     PLWCYC - NUMBER OF CYCLES OF LOWEST FRE-
C.                              QUENCY (BACKGROUND) IN DATA
C.                     KORDER - ORDER OF BINOMIAL FILTER
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
CL              2.         DETERMINE ORDER OF FILTER
C
C     NUMBER OF GRID POINTS/CYCLE (LOW FREQUENCY)
         KLORD=KORDER
         IF (prlow*plwcyc.ne.0) THEN
         KPERT = INT(FLOAT(KPTS)/PLWCYC)
         ZHOLD1= LOG(COS(PI*(1.0/FLOAT(KPERT))))
         ZHOLD2= LOG(PRLOW)/2.0
C
         KLORD= INT(ZHOLD2/ZHOLD1)
         ENDIF
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
CL              3.         FILTER DATA USING LOW PASS BINOMIAL FILTER
C
         CALL LOPASS(PY,PYWORK,KPTS,KLORD)
C---------------------------------------------------------------------
C
         RETURN
         END
C
         SUBROUTINE LOPASS(PY,PYWORK,KPTS,KORDER)
C
       DIMENSION   PY(*),    PYWORK(*)
C
C---------------------------------------------------------------------
C
C     LOW PASS BINOMIAL FILTER
         DO 3 K=1,KORDER
         DO 1 J=2,KPTS-1
    1    PYWORK(J)=(PY(J+1)+2.0*PY(J)+PY(J-1))/4.0
C
C     SEND SMOOTHED DATA TO ORIGINAL DATA ARRAY
         DO 2 L=2,KPTS-1
    2    PY(L)=PYWORK(L)
C
    3    CONTINUE
         RETURN
         END
c/ 21:30 22-apr-94 .../baldur/code`/bald/dimprad.f
c  BALDUR  file dimprad.f maintained by Glenn Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type
c cfs get /11040/bald92/wbaldn1
c end
c lib wbaldn1 ^ x dimprad ^ end
c
c**********************************************************************c
c
c     this file contains the following packages:
c  imprad - impurity radiation package
c  neuset - sets default values of neutral impurity influx parameters
c  neudep - calculates neutral impurity influx deposition profile
c  Impurity Rate Equations subroutines:
c  ncdata - sets input data
c  ncsorc - calculates source deposition profile
c  ncrats - calculates reaction and radiation rates
c  ncdifu - calculates diffusion and drift coefficients
c  nccoef - assembles coefficients of rate balance equations
c  ncsolv - solves tridiagonal matrix system for charge state densities
c  nccons - totals conservation checks
c
c  synch  - synchrotron radiation package by S. Tamor
c  syopac - absorbtion coefficients for the synchrotron package
c  syedit - output routine for the synchrotron package
c
c    The following routines are the same synchrotron radiation package
c    modified by Linda Sugiyama at MIT
c
c  syndrv - driver for tamor's synchrotron radiation package by Sugiyama
c  synch2 - corresponds to sbrtn synch above
c  syopac2 - corresponds to sbrtn syopac above
c  syedit2 - corresponds to sbrtn syedit above
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c@imprad  /11040/baldur/code/imprad.f
c  rgb 25-mar-02 changed local dimensions from 2 to numimp
c    where numimp is located in cbparm.m
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 07-feb-00 removed equivalence statement
c  rgb 29-may-96 removed if ( .not. inita1 ) go to 612
c  rgb 02-feb-95 call ncinfl to compute neutral impurity influxing
c  rgb 18-apr-94 implemented ftzeff(it) = Z_eff
c    changed from gtflwi to timp and from gflowi to flimp
c  rgb 17-apr-94 flimp(ii,it) = Z_eff when cfutz(200) = 3.0
c  rgb 16-apr-94 documented neutral impurity influx
c  les 16-jan-91 helium neutral influxing (cfutz(70)=1)
c         change test for 4He index to nimp=-5 to match 1 1/2D Baldur
c         if d-3he (cfutz(490)>0) do not skip helium initializations
c          for subsequent time steps (go to 710)
c  les  dec-90  change 'Englade" synchrotron to corrected version
c             add natomc=2 to d-3he; coronal radiation for impurity 4
c             add cfutz(iheflx=70)=1 to d-3he; helium-3 and 4 recycling
c  les  nov-90  add d3he fusion
c               synchrotron radiation, new bremsstrahlung
c             for d3he require natomc=1; own synchrotron package
c               and inputs in separte namelist - differs from ENGLADE
c  ahk 12-aug-90 combined Englade's radiation calc with standard code
c                changes are set off by ahk comment statements
c                parameter statement and common/insink/ added to connect
c                with synch subroutine
c      LIMPRD(31) = 1 results in using Englade's coding
c      LIMPRD(4)  = 1 results in using Tamor's syncrotron package
c  rgb 02-jun-89 set incflx=110 and ixdoff=139 in data statement
c  dps 06-jun-89 15.09 Fix bug in webrs calculation for natomc=3.
c  dps 15-may-89 15.09 Move NONCOR call to STEPON.
c  dps 09-may-89 15.08 Replace weirs smoothing with calls to SMOOTH (results
c                slightly different, but previous version not working at
c                present); move cmean and c2mean extrapolation to GETCHI.
c  dps 11-jan-89 **.** Change real*8 to real to satisfy CFT77.
c  dps 26-sep-88 15.03 Skip neutral impurity influx and recycling when
c                natomc = 3.
c  dps 24-aug-88 15.00 Make changes for NC code: skip CE calculation,
c                add call to NC code, and fix output bug arising from
c                changes made in previous version.
c  dps 17-aug-88 14.04 Move volume and area evaluations for neutral He
c                recycling out of initialization section.
c  dps 15-aug-88 14.04 Move default settings of neutral impurity influx
c                and recycling switches to separate routine, NEUSET. Move
c                deposition profile calculation to NEUDEP.
c  dps 27-aug-87 to fix impurity influx bug, try instead using GTFLWI and
c                GFLOWI only (see BOUNDS).
c  dps 26-aug-87 switched from GFLOWI to FLIMP for CFUTZ(200)=1. impurity
c                influx specification to fix bug; moved CALL SCALER(DWEIRS,..)
c                near LIMPRD(32) implementation to fix bug.
c  rgb 25-jul-87 impreoved documentation
c  rgb 07-jul-87 LIMPRD(32)=2 computes CMEAN, C2MEAN, WEIRS, and  DWEIRS
c                in the ghost zone just outside the plasma.
c  rgb 02-jul-87 LIMPRD(32)=1 turns on linear extap of cmean and c2mean
c                the default is limprd(32)=0, no extrapolation
c  rgb 12-jun-87 linearly extrapolate cmean and c2mean outside plasma
c  rgb 10-sep-86 the following changes come from BALDP53M:
c       drm 29-april-85 changed znorm just before do 638 to get specified
c       rate of neutral impurity influx to be properly normalized
c       drm 20-feb-85 changed mzones to ledge in do 6050, 6070, 6220
c         and 720, and in setting jedge(ii) in do 610
c       drm 8-feb-85 fixed zt calculation to include factor of uist
c       drm 2-feb-85 change ifrac,imprd1,imprd2 to 327, 328, 329
c            and gtflwi(it) replaced by timp(it)
c     rgb 17-nov-85   1 1/2 D upgrade using common commhd
c        plasma surface area = <|del xi|> V'(xi)
c        changed rmins*xzoni(j) to ahalfs(j,2)
c       drm 20-dec-84 changed defaults for c-205 and c-206 to 0.0
c       drm 7-dec-84 add algebraic prad code
c       mhr 22-june-84 corrected natomc=0 coding
c       mhr 15-jun-84 made natomc=0 turn off bremstralung and impurity
c               line radiation
c       fgps 18-feb-83 made call to subroutine pdx compatible with hav-
c                      ing scroff(mxchi,55) located in common/comdf2/.
c       fgps 1-feb-83 re-programmed for more than one impurity species;
c                     but, only species #1 and #2 are allowed to be
c                     influxed as neutrals.
c       fgps 10-oct-82 allow adjustments at cfutz(270) secs for:
c                      (1) scrapeoff loss rates (cfutz(271) through
c                      cfutz(276)); (2) feedback control of any neu-
c                      tral-impurity influx (cfutz(277) & cfutz(278));
c                      (3) cross-ion transport (cfutz(279) through
c                      cfutz(280)).
c       fgps 10-sep-82 allow for changing max & min edge temperatures
c                      during a run:  cfutz(260), --(261), & --(262).
c       aes 12-aug-82 fix computation of alpha-beta for synch. radiation
c       fgps 29-jul-82 added option to adjust total number of impurity
c                      ions if edge temperature is too low or too high.
c                      involves cfutz(224), --(227), --(228), & --(229).
c       aes 13-apr-82 incorporate fgps' neutral impurity influxes
c       aes 28-jan-82 add note about use of iheflx in subr. sprint
c       aes 07-jan-82 fix spacing in format 800
c       aes 01-dec-81 change .le. to .gt. at call to pdx (fix bug)
c       aes 19-nov-81 move scroff printout to mprint and eliminate
c               arrays zscrof and yscrof; move he printout to
c               entry impprt (called from sprint)
c       aes 29-oct-81 array dimensions 52 --> 55 in common/tempry/,
c               arrays zheli0,yscrof,zscrof; also in equivalence
c       aes 28-oct-81 put all data statements before executable code
c       fgps 29-aug-79 included the option of neutral helium influx
c                      plus its subsequent ionization (update 10-sep-79)
c       fgps 19-apr-79 (like dep 17-apr-79) in trivial model removed
c                      by-pass of bremsstrahlung etc. when mimp=0
c       fgps 28-feb-79 made provision for zero values in the
c                      calculation of wesrs(j)
c       amck 18-may-78 change **-1.388 -> **(-1.388)
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 12-jan-78 convert tes to ev for recombination
c       amck 11-jan-78 add recombination, remove h. rad. below 100 ev
c       amck 29-nov-77 remove calpha
c       dep nov 25-77  skip synchotron losses if bzs=0.
c       dep nov 25-77  insert syncotron losses
c       dep 10-oct-77 fixed bug in not calculating intervals
c               for simple bremstrahllung model
c       dep 10-oct-77 fixed bug in not looping over zones,290
c       dep 11-oct-77   forced c2mean to be 1.0001*cmean**2 for
c                               atomic physics model
c       dep 10-oct-77 added bremstrahlllung to simple model
c       amck 9-aug-77 don't omit if mimp.le.0 and k.gt.1
c       amck 1-aug-77 call otable with z of element,
c                       get atomic weight
c       amck 15-jun-77 scale dweirs by 1 erg in kev, not 1 erg in ev
c       amck 15-jun-77 add dweirs to ofit call
c       amck 15-mar-77 fix indexing of nspec
c       amck 15-mar-77 use ofit to get impurity terms
c       amck 31-oct-76 use nspec to get ispec, not nimp
c       amck 12-oct-76 make cmean(,,mzones)=1.0
c       amck 8-oct-76 remove dcdtes
c       amck 7-oct-76 not-fitted dcdtes, better fit
c       amck 7-oct-76 add weirs = k*te for te < 2 ev, = k*te**.5 for
c               te > 20kev
c       amck 8-sep-76 make sure cmean, c2mean >= 1.0
c       amck 8-sep-76 remove 1.e-7 factor in weirs, webrs
c       amck 3-sep-76 put zatcof data statements into tartar.for
c       amck 1-sept-76 tarter atomic physics
c       amck 20-may-76 remove "(k)"
c       dep  1-sept-76  tarter atomic physics
c       amck 20-may-76 remove "(k)"
c       amck 15-apr-76
c******************************************************************************
c
c
        subroutine imprad(k)
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'cbparm.m'
        include 'commhd.m'
        include 'csynrd.m'
c
cdoc
c
c       ------------
c       sbrtn IMPRAD
c       ------------
c
c   2.6  Sbrtn imprad computes impurity charge states, radiation terms
c               neutral impurity and helium influxes.
c               It interfaces with the coronal radiation package.
c               (call subroutine for scrapeoff losses when active)
c               (prescribe influx of neutral impurities)
c
c       Variables computed:
c
c  CMEAN(II,IZ,JZ)  = mean charge state Z of impurity II
c                       at boundary (IZ=1) or zone (IZ=2) JZ
c
c  C2MEAN(II,IZ,JZ) = mean charge state squared Z**2 of impurity II
c                       at boundary (IZ=1) or zone (IZ=2) JZ
c
c  WEBRS(JZ)        = hydrogen and impurity bremsstrahlung radiation
c                       cooling rate in zone center JZ
c                       ergs / sec / cm**3
c
c  WEIRS(II,JZ)     = radiation cooling rate of electrons in zone JZ
c                       due to impurity species II
c                       ergs / sec / cm**3
c
c  DWEIRS(II,JZ)    = D WEIRS(II,JZ) / D TES(2,JZ)   1. / sec / cm**3
c
c               Statements 400 ...  Synchrotron radiation
c
c  WESRS(JZ)        = synchrotron radiation cooling of electrons
c                       in zone JZ   ergs / sec / cm**3
c
c               Statements 500 ...  Hydrogen recombination
c
c  RECOMS(JH,JZ)    = recombination rate of hydrogen isotope IH
c                       in zone JZ   n_e <\sigma v_e >_recomb   sec-1
c                       Note:  S_a^{REC} = - n_a n_e <\sigma v_e >_recomb
c
c               Statements 600 ... Neutral impurity influx
c                                  only if cfutz(200) > epslon
c
c  SIIONS(II,JZ)    = density source of impurity species II
c                       due to ionization and charge exchange
c                       in zone JZ   particles / cm**3 / sec
c
c                       Changes made to WEIRS(II,JZ)  638 -- 642
c
c               Statements 700 ... Helium recling
c
c                       Changes made to SIIONS and WEIRS for helium.
c
c               Statements 900 ...  Smoothing WEIRS(II,JZ)
c                                   only if CFUTZ(199) .GE. 1.
c
cend
c--------1---------2---------3---------4---------5---------6---------7-c
c
        real
     r   zmxflx       , mhance
        logical
     l   inita0       , inita1       , inita2       ,
     l   lprint       ,
     l   limits       , lchnge
c
        dimension
     r   zweigt(10)    , zatcof(6,4)   , zahe0(9)      , zimpnu(mj)    ,
     r   zpn(36)       , zip(36)       , zinflx(numimp), zeloss(numimp),
     r   ziloss(numimp), zxloss(numimp), zion0(numimp) , ziions(mj)    ,
     r   zsigv(mj)     , znergy(numimp), zstart(numimp), zpatha(numimp),
     r   zslnta(numimp), z0imp(numimp) , z0isum(numimp),
     r   zlevel(numimp), ztme(numimp)  , zmxflx(numimp), zhance(numimp),
     r   mhance(numimp), ztotal(numimp), zcronl(numimp), zsmoo1(mj)    ,
     r   zsmoo2(mj)    , ztotl1(numimp), zfrac(numimp) , zrnorm(numimp),
     x   zrsum(numimp) , ztotl2(numimp), zipsum(numimp,numimp)
c
        character *10 ihdum1(idximp)
c
c   les  nov-90  d3he (bremsstrahlung)
c
      dimension  ztim(mj), zbmi(mj), bremei(mj),bremee(mj),gauntp(7)
c
cgb        equivalence
cgb     x   (nlzzzz(1),zimpnu(1))       , (nlzzzz(56),ziions(1))      ,
cgb     x   (nlzzzz(111),ztotal(1))     , (nlzzzz(113),zcronl(1))     ,
cgb     x   (nlzzzz(115),zsmoo1(1))     , (nlzzzz(170),zsmoo2(1))     ,
cgb     x   (nlzzzz(775),ztotl1(1))     , (nlzzzz(787),ztotl2(1))
c
        data    inita0,inita1,inita2 /.true.,.true.,.true./
        data    lprint /.false./
        data    limits,lchnge /.false.,.true./
c
        data    zmult1 /2.3025850930/
        data    zmult2 /.43429448190/
c
c
c       atomic physics for hydrogen and impurities
c
c       coronal equilibrium
c
c       data taken from b. tarter
c       "radiative losses from impurity ions"
c       ucrl-78119, lll
c
c
c       zatcof is the array of coefficients
c       for the polynomial fit for bremsstrahlung loss rates
c
c       zatcof(i,j)
c
c       i-index of polynomial
c       x=sum over i=1,6 (zatcof(i,j,k,m)*log10(te)**i)
c
c       j-index of temperature interval
c               te in kev
c       j=1, 0.002 kev to 0.02 kev
c       j=2, 0.02 kev to 0.2 kev
c       j=3, 0.2 kev to 2.0 kev
c       j=4, 2.0 kev to 20. kev
c
c
c       coronal-h   loss rate
        data zatcof
     r          /-8.875780e+00, 3.368544e+01, 3.004414e+01,
     r            1.306722e+01, 2.780003e+00, 2.216913e-01,
     r           -2.274558e+01, 2.823806e+00, 4.693808e+00,
     r            4.332644e+00, 1.994706e+00, 3.405031e-01,
     r           -2.322560e+01, 4.390245e-01, 3.994050e-02,
     r           -7.795893e-02, 1.040682e-01, 1.536833e-01,
     r           -2.322962e+01, 4.768998e-01,-9.974274e-02,
     r            2.169734e-01,-1.969521e-01, 6.531130e-02/
c
c       data for impurity ionization rates according to an average ion
c       model (post and hulse 1980):   zpn(i)=number of electrons in
c       valence subshell vs atomic number i; zip(i)=ionization poten-
c       tial of valence subshell vs atomic number i.
        data (zpn(i),i=1,36)
     r         /  1.,2.,1.,2.,1.,2.,3.,4.,5.,6.,1.,2.,1.,2.,3.,4.,
     r            5.,6.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,3.,3.,3.,4.,
     r            5.,6.,7.,6./
        data (zip(i),i=1,36)
     r         / .0136,.0246,.0054,.0094,.0083,.01126,.0145,.0136,
     r           .0174,.0216,.0052,.0076,.0060,.0082,.0105,.0104,
     r           .0130,.0158,.0043,.0061,.0065,.0068,.0067,.0068,
     r           .0074,.0079,.0079,.0076,.0077,.0094,.0060,.0079,
     r           .0098,.0098,.0118,.0140/
c
c       zahe0(1) to zahe0(7) contain coefficients for a polynomial
c       fit to the electron-impact ionization rate (freeman and jones
c       1974) for electron temperatures in the range 2. - 10000. ev.;
c       for temperatures above 10000. ev., zahe0(8) and zahe0(9)
c       define the straight-line extrapolation in log-log space.
c
c       electron-impact he ionization rate
        data zahe0
     r          /-.4450917e+02, .2442988e+02,-.1025714e+02,
     r            .2470931e+01,-.3426362e+00, .2505100e-01,
     r           -.7438675e-03,-.5107924e+00,-.1374394e+02/
c
c   les  nov-90   new bremsstrahlung for d 3he fusion
c
      data gauntp/1.0257e+00,7.2304e-01,-3.7421e-01,3.3907e-02,
     1 1.5562e-02,-3.9483e-03,2.6917e-04/
c
c   bremsstrahlung radiation coeffs
c
        data brmei1/4.7e-24/,brmeef/2.76e-26/
c
c------------------------------------------------------------------------------
c
c
        data    iheflx,icfout,idvout,inpath /70,71,72,73/,
     1          ionlos,islant,itmhe0,istart /74,75,76,77/
        data    incflx /110/, ixdoff /139/, ismrad /199/
        data  impneu,nergy1,nergy2,izlos1,izlos2 /200,201,202,203,204/,
     1        ixlos1,ixlos2,istar1,istar2,inpth1 /205,206,207,208,209/,
     2        inpth2,islnt1,islnt2,iajust,iprint /210,211,212,213,214/,
     3        ifact1,ifact2 /215,216/
        data  ileve1,ileve2,ixflx1,ixflx2,minti1 /220,221,222,223,224/,
     1        ihnce1,ihnce2,maxti1,ilower,iraise /225,226,227,228,229/,
     2        itime2,minti2,maxti2,itime3,ichng1 /260,261,262,270,271/,
     3        ichng2,ichng3,ichng4,ichng5,ichng6 /272,273,274,275,276/,
     4        impnu3,ifboff,ixoff3,ncflx3 /277,278,279,280/
        data ifrac,imprd1,imprd2/ 327, 328, 329/
c
        data    iclass /2/,     isub /6/
c
c--------1---------2---------3---------4---------5---------6---------7-c
cinput
c
c       ------------------------------------------------------
c       Input variables used in the impurity radiation package (sbrtn IMPRAD)
c       ------------------------------------------------------ (file DIMPRAD)
 
c     Tamor's synchrochrotron radiation package was implemented by Ron Englade.
c The physics of this package is described in report LAPS-72 SAI-023-81-189-LJ
c by S. Tamor, June 1981.  To use this new package, set
c
c LIMPRD(4)  = 1   to use Tamor's syncrotron package
c           else   (default) use the synchrotron formula from PPPL-1368
c
c  CFUTZ(100) = coefficient multiplying the local synchrotron radiation loss
c              regardless of the choice of limprd(4)              (default 0.0)
c CFUTZ(103) = fraction of extraordinary reflected as extraordinary wave  (0.0)
c CFUTZ(104) = fraction of extraordinary wave reflected as ordinary wave  (0.0)
c CFUTZ(105) = fraction of ordinary wave reflected as extraordinary wave  (0.0)
c CFUTZ(106) = fraction of ordinary wave reflected as ordinary wave       (0.0)
c
c      A reasonable choice of input variables would be:
c  limprd(4) = 1     cfutz(100) = 1.0
c  cfutz(103) = 0.9  cfutz(104) = 0.05  cfutz(105) = 0.05  cfutz(106) = 0.9
c
c      Results from the actual package implemented in BALDUR were found
c to be about 10% larger than the results plotted in Figures 4 and 5
c of Tamor's report.
c      It should be noted that synchrotron radiation can transport
c energy from one part of the plasma (center) to another (edge).
c This process is represented by a negative loss of radiated power
c near the edge of the plasma.
c
c  LIMPRD(31) = 1 implements Ron Englade's model for Bremsstrahlung
c                 radiation.  These changes are surrounded by 'cahk'.
c            else (default) original Bremsstrahlung model (see below)
c
c       NATOMC  Controls the impurity radiation model used
c               First namelist, default = 2
c               = 3  for Impurity Rate Equations
c               = 2  for coronal equilibrium radiation model
c                       C. B. Tarter, J. Spectros. Radia. Trans
c                                       10 (1977) 531.
c               = 1  for bremsstrahlung only,
c                       impurities are assumed to be fully ionized
c               = 0  for neither bremsstrahlung nor line radiation.
c
c
c            The calculation which treats a partial recycling of helium
c       is activated when nimp(ii)=2 and cfutz(iheflx=70).gt.0. are
c       true simultaneously.  the cfutz(*)-factors listed below are
c       used to prescribe the conditions of such an influx:
c
c        cfutz(iheflx) = flag to activate helium recycling.  default = 0.
c              70       Note, that if the value of iheflx is altered, it
c                       must likewise be altered in subroutine sprint.
c                       Note that even if cfutz(70) > 0., there will be
c                       no helium volume source unless helium-4 is
c                       specified as an impurity species nimp(ii)=2.
c
c        cfutz(icfout) = fraction of helium ions crossing out at r = r_wall
c              71       recycled as a volume source; defaulted to .9.
c
c        cfutz(idvout) = fraction of helium ions lost by parallel flow
c              72       in the scrapeoff (into limiter or divertor)
c                       which is recycled as a volume source; default = 0.9
c
c        cfutz(inpath) = limit on neutral helium penetration expressed
c              73      as a fraction of the overall plasma minor
c                      radius (rmins); defaulted to .9.
c
c        cfutz(ionlos) = electron-energy loss (kev) per initial ion-
c              74      ization of neutral helium; defaulted to 0.1.
c
c        cfutz(islant) = cosine of the angle between the inner normal
c              75      at the outer surface and the direction of helium
c                      influx; defaulted to 1.0.
c
c        cfutz(itmhe0) = incoming kinetic energy (kev) of recycled neu-
c              76      tral helium; defaulted to 0.04 kev.
c
c        cfutz(istart) = identifies the outermost zone allowing neutral
c              77      helium ionization, measured in from the outer
c                      edge of the scrapeoff region; defaulted to 1.
c
c
c
c     note:  cfutz(124) and cfutz(125), which multiply the hydrogenic
c            and impurity-ion scrapeoff loss rate respectively, are
c            defaulted to 1.0 in subroutine trcoef.
c
c
c     note:  in subroutine neugas, cfutz(198) is used to control
c       "tcold", the incoming temperature of recycled hydrogenic neu-
c       trals.  when cfutz(198).gt.epslon, "tcold" is readjusted
c       (if necessary) so as not to exceed the average ion temperature
c       at a radial point equal to  cfutz(198)*plasma minor radius.
c       the default (cfutz(198)=0) is to bypass this control.
c
c     An option is included to smooth sharp peaks in impurity-radi-
c       ation losses which affect the electron energy.
c
c        cfutz(ismrad)=0,1,2.  if non-zero, a binomial-weighted smoothing
c              199     of order cfutz(ismrad)+1 [over 5 (7) adjoining radial
c                      zones for cfutz(ismrad)=1.(2.)] is performed.
c                      Otherwise, this smoothing is bypassed.
c                      Note that this option is independent of
c                      neutral or ionized impurity influx.  default = 0.
c
c     A neutral-impurity influx at the outer boundary is activated
c  by setting cfutz(200) .gt. 0.0, where as always impurity spe-
c  cies are specified by nimp(ii).  However, only impurity species
c  ii=1 and ii=2 are allowed to be influxed as neutrals,
c  and the setting of cfutz(200) affects both ii=1 and ii=2.
c  Otherwise, prescribed impurity influxes are considered
c  to be ionized (see FLIMP(II,JT)).
c
c       When influxing impurities ii=1 and/or ii=2 as neutrals,
c       two alternative outer boundary conditions are available:
c   (1) setting "nbound=0 or 1" leads to fixed pedestal or
c       corrected pedestal boundary conditions; while,
c   (2) setting "nbound=2 or 3"
c       results in zero impurity-ion influx for species ii=1 and ii=2.
c
c       Using the input arrays flimp(ii,it) and timp(it), neutral-im-
c       purity influxing rates are determined in one of two ways:
c   (a) when cfutz(200)=1.0, flimp(ii,it) handles a prescribed
c       neutral influx of impurity ii (per cm*cm per sec) versus time
c       timp(ii) in secs.
c   (b) when cfutz(200).gt.1.0, flimp(ii,it)
c       contains the desired total number of impurity-ii ions at time
c       timp(it) in secs.
c   (c) when cfutz(200) > 2.0, flimp(ii,it) contains Z_eff,
c         which is dimensionless, at time timp(it) is seconds.
c
c       Once a method of influx control has been set,
c       it applies to both ii=1 and ii=2.
c
c  For all influxed impurities to be influxed as ions,
c       both cfutz(200)=0.0 and either nbound=2 or =3.
c
c  However, if either "ii=1" and/or "ii=2" are
c       being influxed as neutrals, regardless of whether nbound=0, 1,
c       2, or 3, "ii=3" and/or "ii=4" can also be influxed across the
c       "wall" boundary, but only as ions; the activator here is that
c       flimp(ii,it) be greator than zero for ii=3 and/or =4 at certain
c       of the prescribed times timp(it).
c
c  Incidentally, zero boundary
c       flux can be imposed on any impurity-ion species ii by making
c       at least one element of flimp(ii,it) slightly greater than
c       epslon (e.g., 10.**-34) at some time timp(it), where timp(it)
c       can be .gt. tmax, where nbound=2 or =3, and where all other
c       elements of flimp(ii,it) are set to zero (see note in subrou-
c       tine bounds).
c
c  Since all boundary conditions are handled in subroutine bounds,
c       an indicator must be passed to this subroutine by inserting
c       "data impneu/200/" there and using cfutz(impneu) as a flag.
c
c  Available with impurity influx, in subroutine convrt
c       there is an adjustable empirical convective drag on
c       impurity ions activated by cfutz(217), modified by cfutz(218)
c       and cfutz(219).
c
c  The following parameters set the conditions for
c       neutral-impurity influx:
c
c
c        cfutz(impneu) = 0.0 bypasses any neutral-impurity influx
c              200     = 1.0 allows only prescribed influxes (/cm**2/sec)
c                        at times specified by timp(it) (sec).
c                      = 2.0 allows only the setting of desired
c                        total numbers of impurity ions at specified
c                        times in secs.
c                        note: flimp(ii,it) = N_imp_tot as a fn of time
c                      = 3.0 use Z_eff to set total number of impurity
c                        note: flimp(ii,it) = Z_eff as a fn of time
c                        default is 0.0.
c
c        cfutz(nergy1) = prescribed incoming kinetic energy (kev) for
c              201     neutral impurity #1; defaulted to .01 kev.
c
c        cfutz(nergy2) = prescribed incoming kinetic energy (kev) for
c              202     neutral impurity #2; defaulted to .01 kev.
c
c        cfutz(izlos1) = electron-energy loss (kev) per initial ion-
c              203     ization of a neutral atom of impurity #1; de-
c                      faulted to .0136 kev.
c
c        cfutz(izlos2) = electron-energy loss (kev) per initial ion-
c              204     ization of a neutral atom of impurity #2; de-
c                      faulted to .0136 kev.
c
c        cfutz(ixlos1) = extra radiation-loss (kev) per impurity atom
c              205     #1 due to nonequilibrium; defaulted to 10. kev.
c
c        cfutz(ixlos2) = extra radiation-loss (kev) per impurity atom
c              206     #2 due to nonequilibrium; defaulted to 10. kev.
c
c        cfutz(istar1) = determines the radial index (jedge(1)) of the
c              207     outermost zone allowing ionization of neutral
c                      impurity #1.  default is 1.
c
c        cfutz(istar2) = determines the radial index (jedge(2)) of the
c              208     outermost zone allowing ionization of neutral
c                      impurity #2.  default is 1.
c
c        cfutz(inpth1) = prescribed penetration limit for neutral im-
c              209     purity #1 expressed as a fraction of the plasma
c                      (or separatrix) radius r(jsepx).  defaulted so
c                      that max penetration equals 0.9*r(jsepx).
c
c        cfutz(inpth2) = determines penetration limit for neutral im-
c              210     purity #2; handled like cfutz(inpth1).
c
c        cfutz(islnt1) = prescribed cosine of the angle between the in-
c              211     ward-directed normal to the outer plasma sur-
c                      face and the influx direction of neutral impur-
c                      ity #1.  used in a rough approx of controlled
c                      gas feed.  default is 1.0.  internally limited
c                      so that 1.0.ge.cfutz(islnt1).ge.0.5.
c
c        cfutz(islnt2) = prescribed cosine parameter like cfutz(islnt1)
c              212     but for neutral impurity #2.
c
c        cfutz(iajust) = 1.0 activates an adjustment of the ionization
c              213     rates determined by an empirical formula due to
c                      d. e. post (1980) so as to approximate ( <15% )
c                      values quoted by freeman and jones (1974).  de-
c                      fault is 0.0 which results in uncorrected use of
c                      the post formula.
c
c        cfutz(ifact1) = a constant factor which for neutral impurity #1
c              215     can be used to increase or decrease the ioniza-
c                      tion rates from the post formula.  defaulted to
c                      1.0.
c
c        cfutz(ifact2) = a constant factor which for neutral impurity #2
c              216     can be used to increase or decrease the ioniza-
c                      tion rates from the post formula.  defaulted to
c                      1.0.
c
c        cfutz(iprint) = time (secs) immediately after which impurity-
c              214     ion source rates vs radial index are printed
c                      out before and after calibration with respect
c                      to the neutral influx.  non-equilibrium radia-
c                      tion enhancement factors are also displayed.
c                      the process is repeated over 3 successive time
c                      steps.  default is 0.0 secs.
c
c        cfutz(ileve1) = the factor available to adjust the influx rate
c              220     which seeks to obtain (flimp(ii,it) vs timp(it))
c                      a prescribed total number of impurity-1 ions.
c                      the default value is 1.0.
c
c        cfutz(ileve2) = an adjustment factor similar to cfutz(ileve1) but
c              221     for impurity ions 2.  defaulted to 1.0.
c
c        cfutz(ixflx1) = maximum recycled neutral influx (per cm*cm per
c              222     sec) of impurity #1.  defaulted to 1.e17 per
c                      cm*cm per sec.
c
c        cfutz(ixflx2) = maximum recycled neutral influx (per cm*cm per
c              223     sec) of impurity #2.  defaulted to 1.e17 per
c                      cm*cm per sec.
c
c        cfutz(ihnce1) = maximum radiation enhancement factor for im-
c              225     purity #1 with respect to coronal equilibrium.
c                      defaulted to 10.
c
c        cfutz(ihnce2) = maximum radiation enhancement factor for im-
c              226     purity #2 with respect to coronal equilibrium.
c                      defaulted to 10.
c
c        cfutz(minti1) = switch for allowing feedback adjustment in the
c              224     desired total number of impurity ions and for
c                      setting a minimum ion temperature (ev) at the
c                      plasma edge (jsepx1).  applies only when cfutz
c                      (impneu) = 2.0; effects all influxed ion species.
c                      defaulted to 0.0 which means no adjustment.
c
c        cfutz(maxti1) = prescribed maximum ion temperature (ev) desired
c              227     at the plasma edge (jsepx1).  defaulted to 1.e34.
c
c        cfutz(ilower) = factor lowering the prescribed neutral-impurity
c              228     influx if ti(jsepx1) is too cold.  default value
c                      is 0.90.
c
c        cfutz(iraise) = factor raising the prescribed neutral-impurity
c              229     influx if ti(jsepx1) is too hot.  default value
c                      is 1.10.
c
c        cfutz(itime2) = prescribed time (secs) at which both the min-
c              260     imum and maximum desired edge ion temperatures
c                      are changed from cfutz(minti1) and cfutz(maxti1)
c                      respectively to cfutz(minti2) and cfutz(maxti2).
c                      defaulted to 1.e34.
c
c        cfutz(minti2) = prescribed minimum edge ion temperature (ev)
c              261     which is to prevail after time cfutz(itime2).
c                      defaulted so that cfutz(minti2) = cfutz(minti1).
c
c        cfutz(maxti2) = prescribed maximum edge ion temperature (ev)
c              262     which is to prevail after time cfutz(itime2).
c                      defaulted so that cfutz(maxti2) = cfutz(maxti1).
c
c********note**********print-out of "zadjst", the influx readjustment
c                      factor, is activated in subroutine sprint when
c                      cfutz(impneu) > 1.1.
c
c       parameters used in subroutine convrt and having to do with
c       enhanced outflux of impurity ions near the boundary:
c             i0drag,i1drag,i2drag /217,218,219/
c
c
c       at time cfutz(itime3=270) in seconds
c       changes can be made in certain  cfutz(*) parameters.
c       In particular, changes can be prescribed
c       for scrapeoff loss rates, for neutral-impurity feedback control,
c       and for cross-ion transport.
c
c        cfutz(ichng1) = the prescribed factor multiplying the value of
c              271     scroff(1,j) obtained from subroutine pdx(*,*,*)
c                      after time cfutz(itime3) in seconds, where "j"
c                      includes all values of the radial index which lie
c                      in the scrapeoff.  apropos default, no change is
c                      made in the output values from pdx(*,*,*) unless
c                      cfutz(itime3) is ".gt." zero and ".lt." infinity;
c                      but, locally cfutz(ichng1) is defaulted to 0.0.
c
c        cfutz(ichng2) thru cfutz(ichng6) are the same as cfutz(ichng1)
c              272     but with respect to scroff(2,j) thru scroff(6,j).
c                      all these are defaulted to 0.0.  note: only 6
c                      readjustment factors are available; hence, when
c                      more than 2 impurity species are prescribed, then
c                      scroff(lion,j) and perhaps scroff(lelec,j) are
c                      not adjustable.
c
c        cfutz(impnu3) = (if greater than 0.0) a prescribed new value of
c              277     cfutz(impneu) which takes affect after cfutz(i
c                      time3) seconds.
c
c        cfutz(ifboff) = if set greater than 0.0, this removes the feed-
c              278     back adjustments (after cfutz(itime3) seconds) on
c                      the prescribed total numbers of impurity ions.
c
c        cfutz(ixoff3) = (if greater than 0.0) a prescribed new value for
c              279     cfutz(ixdoff) which takes affect after cfutz(i-
c                      time3) seconds.  cfutz(ixdoff) adjusts the cross-
c                      ion transport coefficients:  if =0, cross-ion
c                      transport coefficients are left as first calcu-
c                      lated; if >0, cross-ion coefficients are multi-
c                      plied by cfutz(ixdoff) before use; and, if <0,
c                      cross-ion transport is set to zero.
c
c        cfutz(ncflx3) = (if greater than 0.0) a new value of cfutz(incflx)
c              280     which takes affect after cfutz(itime3) seconds.
c                      cfutz(incflx) greater than 0.0 activates hawryluk-
c                      hirshman ion transport.
cend
c--------1---------2---------3---------4---------5---------6---------7-c
c
cbate        if(.not.inita0) go to 5
        inita0=.false.
        i1=mxhyd*mxzone
        i2=mximp*mxzone
        i3=mximp*i1
        i4=mximp*i2
        i5=mxhyd*i1
        max0=i2
        max1=mximp*mzones
        max2=mxchi*mxzone
        mimp0=min0(mimp,2)
        jsmin=lcentr
        if(nadump(1).gt.lcentr) jsmin=nadump(1)
        zhance(1)=1.0
        zhance(2)=1.0
        zchnge=epsinv
        if(cfutz(itime3).gt.epslon) zchnge=cfutz(itime3)*usit
        if(cfutz(ichng1).le.epslon) cfutz(ichng1)=0.0
        if(cfutz(ichng2).le.epslon) cfutz(ichng2)=0.0
        if(cfutz(ichng3).le.epslon) cfutz(ichng3)=0.0
        if(cfutz(ichng4).le.epslon) cfutz(ichng4)=0.0
        if(cfutz(ichng5).le.epslon) cfutz(ichng5)=0.0
        if(cfutz(ichng6).le.epslon) cfutz(ichng6)=0.0
c
    5   continue
c
        if(.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.6 subroutine imprad bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
        if ( tbi .lt. zchnge ) go to 12
        if ( .not. lchnge ) go to 12
        lchnge=.false.
c
c       re-adjust feedback control on neutral-impurity influx
c
        if ( cfutz(impnu3) .gt. epslon ) cfutz(impneu)=cfutz(impnu3)
        if ( cfutz(ifboff) .gt. epslon ) limits=.false.
c
c       modify cross-ion transport coefficients
c
        call resetr(dnhis,i3,0.0)
        call resetr(dnihs,i3,0.0)
        call resetr(dnhhs,i5,0.0)
        call resetr(dniis,i4,0.0)
c
        if ( cfutz(ixoff3) .gt. epslon ) cfutz(ixdoff)=cfutz(ixoff3)
c  where cfutz(ixdoff) is used in subroutine trcoef to adjust
c    cross-ion transport.
        if ( cfutz(ncflx3) .gt. epslon ) cfutz(incflx)=cfutz(ncflx3)
   12   continue
c
        call resetr(scroff,max2,0.0)
c
c  call scrapeoff-loss subroutine if nadump(1).gt.0
c  Note that we can now specify a fixed value of the sound
c  speed / connection length, cimprd(1) (should be > 0).
c
        if ( nadump(1) .gt. lcentr ) call pdx
c
        if ( .not. lchnge ) then
          do js=jsmin,ledge
            scroff(1,js) = cfutz(ichng1)*scroff(1,js)
            scroff(2,js) = cfutz(ichng2)*scroff(2,js)
            scroff(3,js) = cfutz(ichng3)*scroff(3,js)
            scroff(4,js) = cfutz(ichng4)*scroff(4,js)
            scroff(5,js) = cfutz(ichng5)*scroff(5,js)
            scroff(6,js) = cfutz(ichng6)*scroff(6,js)
          enddo
        endif
c
c
cl      1)      initialize coronal physics
c
c
  100   continue
        if (k.gt.1) go to 200
c
        if (mimp.le.0) return
c
      call otable(nzspec(limp1),mximp,mimp,ihdum1,lholab,nounit,zweigt)
c
c   les  nov-90  d3he fusion (?right place?)
c
      if ( cfutz(490) .gt. epslon) then
        aspec(6) = zweigt(1)
        lhspec(6)=ihdum1(1)
        go to 1101
        endif
c
cahk
        IF (LIMPRD(31) .EQ. 0) ihyd = 0
cahk
        do 110 ji = 1, mimp
        ii = ji+lhydn
cahk
        IF (LIMPRD(31) .EQ. 0  .AND.  nimp(ji) .EQ. 3) ihyd = ji
cahk
        if (aspec(ii).le.epslon) aspec(ii) = zweigt(ji)
        if (aspec(ii).le.epslon) aspec(ii) = 1.0
        if (nspec(ii).lt.0) go to 110
c
        lhspec(ii) = ihdum1(ji)
c
  110   continue
cahk
        IF (LIMPRD(31) .EQ. 1) THEN
          if ( ihyd .gt. 0 ) nzspec(ihyd+lhydn)=1
          if ( ihyd .gt. 0 ) aspec(ihyd+lhydn)=1.
          if ( ihyd .gt. 0 ) zweigt(ihyd)=1.
        ENDIF
cahk
c
 1101   continue
c
        return
c
c
cl      2)      trivial model
c
c
  200   continue
c
        if (natomc.ge.3) go to 280
c
c       assume impurities totally ionized.
c
c
c   les  nov-90  d3he
c   d-3he fusion requires natomc=1 (2-jan-90) or natomc=2 (dec-90)
c d 3he fusion --- bremsstrahlung for te.gt.20 keV
c
        if ( cfutz(490) .le. epslon ) go to 280
c
      do 242 jz=1,mzones
        do 240 ji=1,mimp
          ii = ji+lhydn
          cmean(ji,2,jz) = float(nzspec(ii))
 240      c2mean(ji,2,jz) = cmean(ji,2,jz)**2
 242    continue
c
c   bremsstrahlung (erg/cm3/s)
c
      do 258 jz=lcentr,ledge
        ztim(jz)=log10(tes(2,jz)*useh*1000.)
        ztim(jz)=max(1.,ztim(jz))
 258    ztim(jz)=min(5.,ztim(jz))
c
c     horner's rule for gaunt factor = bmi
c
      do 2591 jz=lcentr,ledge
        gaunt=gauntp(7)
        do 259 jk=1,6
 259    gaunt=ztim(jz)*gaunt+gauntp(7-jk)
 2591 zbmi(jz)=gaunt
      zbrei=brmei1*sqrt(useh)
      zbreef=brmeef*useh**1.5
c
c   ztim=rhoels(2,jz)*xzeff(2,jz) - (impurity 4, helium)
c
      do 2592 jz=lcentr,ledge
 2592  ztim(jz)=xzeff(2,jz)*rhoels(2,jz)
c
      if ( natomc .ne. 2  .or.  mimp .lt. 4 ) go to 2596
      do 2593 jz=lcentr,ledge
 2593  ztim(jz)=ztim(jz)-rhois(4,2,jz)*c2mean(4,2,jz)
 2596 do 260 jz=lcentr,ledge
      bremei(jz)=zbrei*ztim(jz)*rhoels(2,jz)*sqrt(tes(2,jz))
     2  *zbmi(jz)*(1.+tes(2,jz)*useh/255.)
 260  bremee(jz)=zbreef*rhoels(2,jz)**2*tes(2,jz)**1.5
c
        do 264 jz=lcentr,ledge
 264    webrs(jz)=bremei(jz)+bremee(jz)
      if ( natomc.ne.2 .or. mimp.lt.4 ) go to 400
        mimpz=1
        mindex=4
        go to 316
 280  continue
c
        if (natomc.ge.2) go to 300
c
c       assume impurities totally ionized.
c
c
c       zmult2=1./log(10.), log10(x)=log(x)*zmult2
c       zmult1=log(10.), 10.**x=exp(x*zmult1)
c
        do 290 jz = 1, mzones
c
        if ( mimp .le. 0 ) go to 221
        do 220 ji = 1, mimp
          ii = ji+lhydn
          cmean(ji,2,jz) = float(nzspec(ii))
          c2mean(ji,2,jz) =  cmean(ji,2,jz)**2
  220   continue
cahk
        IF (LIMPRD(31) .EQ. 1) THEN
          if ( ihyd .gt. 0 ) cmean(ihyd,2,jz)=1.
          if ( ihyd .gt. 0 ) c2mean(ihyd,2,jz)=1.
        ENDIF
cahk
  221   continue
        if ( natomc .ne. 0 ) then
c
c       compute log10(te), te in kev
c
        zte = abs(tes(2,jz)*evsinv*.001)
c
c       zte is the electron temperature in kev
c
c
        ztex = max(min(zte,20.0),0.100)
        ztx=log10(ztex)
c
        iint=1
        if ( ztex .gt. 0.02 ) iint=2
        if ( ztex .gt. 0.2 ) iint=3
        if ( ztex .gt. 2. ) iint=4
c
c       iint is the interval index
cahk
        IF (LIMPRD(31) .EQ. 1) THEN
          zploss=5.34e-24*zte**.5*(1.+zte/255.)*xzeff(2,jz)+
     &      2.76e-26*zte**1.5
          zsum=0.
          webrs(jz)=rhoels(2,jz)*rhoels(2,jz)*zploss
        ELSE
c
          zploss=0.
c
          do 210 j=1,6
            jl=7-j
            zploss=zploss*ztx+zatcof(jl,iint)
  210     continue
c
          zploss=exp(zploss*zmult1)
          if (zte.lt.ztex) zploss = zploss * sqrt(zte/ztex)
          zsum=0.
          webrs(jz)=xzeff(2,jz)*rhoels(2,jz)*rhoels(2,jz)*zploss
        ENDIF
cahk
c-----webrs is hydrogen and impurity bremsstrahlung
c               radiation loss in ergs/sec/cm**3
c
        else
        webrs(jz)=0
        end if
 290    continue
c
c
        go to 400
c
c
c  3)   coronal equilibrium atomic physics
c       ----------------------------------
c
c
  300   continue
c
c       zmult2=1./log(10.), log10(x)=log(x)*zmult2
c       zmult1=log(10.), 10.**x=exp(x*zmult1)
        do 315 jz=1,mzones
c
c       compute log10(te), te in kev
c
        zte = abs(tes(2,jz)*evsinv*.001)
c
c       zte is the electron temperature in kev
c
        ztex = max(min(zte,20.0),0.100)
        ztx=log10(ztex)
c
        iint=1
        if ( ztex .gt. 0.02 ) iint=2
        if ( ztex .gt. 0.2 ) iint=3
        if ( ztex .gt. 2. ) iint=4
c
c       iint is the interval index
c
        zploss=0.
c
c       compute hydrogen radiation first
c
        do 310 j=1,6
          jl=7-j
          zploss=zploss*ztx+zatcof(jl,iint)
  310   continue
c
c       compute hydrogen ion density
c
        zploss=exp(zploss*zmult1)
        if (zte.lt.ztex) zploss = zploss * sqrt(zte/ztex)
        zsum=0.
        do 312 jh=1,mhyd
          zsum=zsum+rhohs(jh,2,jz)
  312   continue
        webrs(jz)=zsum*rhoels(2,jz)*zploss
c-----webrs is hydrogen radiation loss in ergs/sec/cm**3
c
c
cl              impurity radiation, z-bar, and z**2 bar
c
c..15.09 Bug fix: as set up previously, skipped out of webrs
c        calculation for natomc=3 or mimp = 0 (> version 14.99)
c
        if (versno.lt.15.09) then
          if (mimp.le.0) go to 321
          if (natomc.eq.3) go to 320
        else
          if ((mimp.le.0).or.(natomc.eq.3)) go to 315
        end if
c
c   les 1990; separate the do-loop for hydrogen and impurities
c            add entry point 316 for d-3he
c
  315   continue
c
      if ( mimp .le. 0  .or.  natomc .eq. 3 ) go to 319
        mindex=1
        mimpz=mimp
  316   continue
c
      do 318 jz=1,mzones
        zte=abs(tes(2,jz)*evsinv*.001)
c
        call ofit(zte,rhoels(2,jz),rhois(mindex,2,jz),
     1          weirs(mindex,jz),dweirs(mindex,jz),
     2          cmean(mindex,2,jz),c2mean(mindex,2,jz),mimpz)
  318   continue
  319   continue
c
c  ...Rescale dweirs from keV**(-1) to ergs**(-1)
c
        im=mximp*mzones
        call scaler(dweirs,im,evsinv*0.001)
c
 320   continue
 321   continue
c
c  3.3) reset impurity radiation             327
c       with the total power equal to cfutz(ifrac)*(heating power)
c
        if ( nstep .eq. 0 ) go to 370
        if ( cfutz(ifrac) .eq. 0 ) go to 370
c
      if ( cfutz(ifrac) .lt. 0. ) then
        iimp = min0(2,mimp)
c
        if ( cfutz(imprd1) .lt. 0. ) then
          zfrac(1)=-cfutz(imprd1)*rhois(1,2,lcentr)/rhoels(2,lcentr)
        else
          zfrac(1)=rhois(1,2,lcentr)/rhoels(2,lcentr)
        endif
c
        if ( cfutz(imprd1) .gt. 0. ) then
          zrnorm(1)=cfutz(imprd1)
        else
          zrnorm(1)=0.
        endif
c
        if ( cfutz(imprd2) .lt. 0. ) then
          zfrac(2)=-cfutz(imprd2)*rhois(2,2,lcentr)/rhoels(2,lcentr)
        else
          zfrac(2)=rhois(2,2,lcentr)/rhoels(2,lcentr)
        endif
c
        if ( cfutz(imprd2) .gt. 0. ) then
          zrnorm(2)=cfutz(imprd2)
        else
          zrnorm(2)=0.
        endif
      endif
c
c       loop over radius to get total heating power and reset weirs
c   les  nov-90 added  d3he terms
c
        zsum=0.
        znorm=0.
        zrsum(1)=0.
        zrsum(2)=0.
        do 340 jz=lcentr,ledge
          zsum=zsum+dx2i(jz)*(weohms(jz)
     1     + wealfs(jz) + wialfs(jz) + weauxs(jz) + wiauxs(jz)
     2     + webems(jz) + wibems(jz) + weecrh(jz) + wiecrh(jz)
     3     + weicrf(jz) + wiicrf(jz)
     4     + wed3fs(jz) + wid3fs(jz))
          if ( cfutz(ifrac) .gt. 0. ) then
c
c       use  algebraic form for weirs
c
            if ( xzoni(jz) .lt. cfutz(imprd2) ) then
              weirs(1,jz) = cfutz(imprd1) + (1.-cfutz(imprd1))
     x          *xzoni(jz)/cfutz(imprd2)
            else
              weirs(1,jz)=(1.-xzoni(jz)**2)**0.25/
     1          (1.-cfutz(imprd2)**2)**0.25
            endif
            znorm=znorm+weirs(1,jz)*dx2i(jz)
          else
c
c       set weirs as if impurity density were proprtional to
c       electron density
c
        do 335 jimp=1,iimp
          zscale=zfrac(jimp)*rhoels(2,jz)/rhois(jimp,2,jz)
          weirs(jimp,jz)=zscale*weirs(jimp,jz)
          dweirs(jimp,jz)=zscale*dweirs(jimp,jz)
          zrsum(jimp)=zrsum(jimp)+weirs(jimp,jz)*dx2i(jz)
          if(jz.eq.ledge) zrnorm(jimp)=zrnorm(jimp)*zsum/zrsum(jimp)
335     continue
        endif
340     continue
        if ( znorm .ne. 0. ) znorm=zsum*cfutz(ifrac)/znorm
c
c       renormalize weirs to appropriate fractions of heating power
c
        do 350 jz=lcentr,ledge
          if ( cfutz(ifrac) .gt. 0. ) then
            dweirs(1,jz)=0.
            weirs(1,jz)=znorm*weirs(1,jz)
            if(mimp.ge.2) then
              do ji=2,mimp
                dweirs(ji,jz)=0.
                weirs(ji,jz)=0.
              enddo
            endif
          else
            do 347 jimp=1,iimp
              if ( zrnorm(jimp) .eq. 0. ) go to 347
                weirs(jimp,jz)=zrnorm(jimp)*weirs(jimp,jz)
                dweirs(jimp,jz)=zrnorm(jimp)*dweirs(jimp,jz)
347         continue
          endif
350     continue
370     continue
c
c
 400    continue
c
c       4 ) cyclotron losses  wesrs(j)  ( only if CFUTZ(100) > eplson )
c           ----------------
c
c       reflection coefficient is input as cfutz(100)
c
        call resetr (wesrs,mzones,0.0)
c
c   les  nov-90  d3he version of synchrotron rad
c       for 10<Te<120 keV
c
      if ( cfutz(490) .gt. epslon ) then
         nzons1=mzones-lcentr
         call syndrv
         go to 500
      endif
c
        if ( cfutz(100) .le. epslon
     &       .or. bzs .le. epslon
     &       .or. tes(2,lcentr) .le. epslon
     &       .or. rhoels(2,lcentr) .le. epslon ) go to 500
cahk
      IF ( LIMPRD(4) .EQ. 1 ) THEN
c
c use Englade's coding for computing synchrotron radiation
c this uses subroutine synch writen by S. Tamor
c
        ndbug=0
        iedit=0
        nzons=mzones-2
        bmax=0.
        bmin=1.0e30
c
c cfutz(103) and cfutz(104) used to input reflectivity information
c cfutz(105) and cfutz(106) used to input non-refl area eg port area
c
        ree=cfutz(103)
        roo=cfutz(104)
        roe=cfutz(105)
        reo=cfutz(106)
        do 470 jz=1,nzons
          surf(jz)=1.0e04*avi(jz+2,3,1)*avi(jz+2,5,1)
          voq(jz)=vols(jz+2,1)-vols(jz+1,1)
          teqp(jz)=tes(2,jz+1)/1.602e-09
          deny(jz)=rhoels(2,jz+1)
          bavg(jz)=5.*(avi(jz+2,8,1)*avi(jz+2,11,1)+
     &      avi(jz+1,8,1)*avi(jz+1,11,1))
          bmax=max(bavg(jz),bmax)
          bmin=min(bavg(jz),bmin)
  470   continue
        areaw=surf(nzons)
        if ( areaw .gt. epslon ) call synch
        do jz=1,nzons
           wesrs(jz+1) = cfutz(100) * ploss(jz) * 1.602e-09
        enddo
c
c..sychrotron radiation formula taken from PPPL-1368 , #3
c
      ELSE
        z0 = 8. * fcpi / (bzs*bzs)
         do 420 jz = lcentr, ledge
          zbeta = (tes(2,jz) * rhoels(2,jz) + tis(2,jz)* rhoins(2,jz)
     1      + .666666666 * hebems(jz) * rhobis(2,jz)
     2      + .666666666 * alphai(jz) * ealfai(jz) * uisd * uise )
     3      * z0
          wesrs(jz) = ( 2.5e-25 * rhoels(2,jz)**2 ) *
     1      ( tes(2,jz) * evsinv * .001  )**2 * cfutz(100)
     2      *  ((1.-zbeta)/(zbeta + epslon)) *
     3      ( 1. +tis(2,jz) * rhoins(2,jz)/ (tes(2,jz) * rhoels(2,jz) ))
     4      * ( 1. + tes(2,jz)*  evsinv * .001  / 204. )
 420     continue
c
        ENDIF
cahk
c
c
cl      5)      recombination
c               -------------
c
c
  500   continue
c
        if ( .not. nlrcom ) go to 600
c
        do 518 jz = lcentr, ledge
c
        zte = max(tes(2,jz)*evsinv,0.001)
        if ( zte .gt. 400.0 ) then
          zrecom = rhoels(2,jz) * 5.2179e-12 * zte**(-1.388)
        else
          zrecom = rhoels(2,jz) * 1.27e-13 /
     &          ( sqrt(zte*0.073529) * (1.0 + 0.59*(zte*0.073529)) )
        endif
c
        do jh = 1, mhyd
          recoms(jh,jz) = zrecom
        enddo
c
  518   continue
c
c
c  6)   Neutral impurity influx   ( only if CFUTZ(200) > epslon )
c       -----------------------
c
  600   continue
c
cbate        write (nprint,*) 'nstep= ',nstep,'  nadump(11)= ',nadump(11)
cbate     &    ,'  natomc= ',natomc,'  cfutz(impneu)= ',cfutz(impneu)
c
c  initially, nadump(11)=0 in order to prevent premature use of
c  neutral-impurity influx
c
        if ( .not.(nadump(11) .gt. 0) ) go to 1060
c
c  15.03 If using impurity charge state transport, must skip these sections
c  dealing with impurity influx and recycling (15.09: and smoothing)
c
        if (natomc.eq.3) return
c
        if ( cfutz(impneu) .lt. epslon ) go to 700
c
        if ( mimp .le. 0 ) return
c
c..treatment of neutral-impurity influx
c
c
cbate        if ( .not. inita1 ) go to 612
c
c..initializations for neutral-impurity influx
c
        inita1=.false.
c
        call neuset(cfutz,epslon,epsinv,1)
c
        zunit1  = evsinv*1.0e-03
        ztprnt  = cfutz(iprint)*ueit
        nxprnt  = 0
c
        jsepx   = mzones
        if ( nadump(1) .gt. lcentr ) jsepx = nadump(1)
c
        do 610 ii=1,mimp0
          ip=ii+lhydn
          iz=nzspec(ip)
          if (ii.le.1) znergy(ii)=cfutz(nergy1)
          if (ii.gt.1) znergy(ii)=cfutz(nergy2)
          zeloss(ii)=evs*1.0e+03
          if (ii.le.1) zeloss(ii)=zeloss(ii)*cfutz(izlos1)
          if (ii.gt.1) zeloss(ii)=zeloss(ii)*cfutz(izlos2)
          zxloss(ii)=evs*1.0e+03
          if (ii.le.1) zxloss(ii)=zxloss(ii)*cfutz(ixlos1)
          if (ii.gt.1) zxloss(ii)=zxloss(ii)*cfutz(ixlos2)
          if (ii.le.1) zstart(ii)=cfutz(istar1)
          if (ii.gt.1) zstart(ii)=cfutz(istar2)
          if (ii.le.1) zpatha(ii)=cfutz(inpth1)
          if (ii.gt.1) zpatha(ii)=cfutz(inpth2)
          if (ii.le.1) zslnta(ii)=1.0/min(cfutz(islnt1),1.0)
          if (ii.gt.1) zslnta(ii)=1.0/min(cfutz(islnt2),1.0)
          mhance(ii)=1.0
          if (ii.le.1) mhance(ii)=cfutz(ihnce1)
          if (ii.gt.1) mhance(ii)=cfutz(ihnce2)
c
c..gaunt-factor coefficients follow
c
          zvalnc=2.0
          if (iz.lt.3) zvalnc=1.0
          if (iz.gt.10) zvalnc=3.0
          if (iz.gt.28) zvalnc=4.0
        zgaunt=.3035*12.18*(1.+.0335*zvalnc)*exp(-zvalnc/(zvalnc+5.0))
c  later zgaunt will be multiplied by abs(zfy).
          zion0(ii)=zpn(iz)*zgaunt/(zip(iz)*zip(iz))
c  include the normalization factor 3.44e-11
          if (ii.le.1) zion0(ii)=zion0(ii)*cfutz(  ifact1)*3.44e-11
          if (ii.gt.1) zion0(ii)=zion0(ii)*cfutz(ifact2)*3.44e-11
c
  610   continue
c
c  end of initialization
c
  612   continue
c
c..control printout
c
        if(tai.lt.ztprnt) go to 614
        lprint=.true.
        nxprnt=nxprnt+1
        if(nxprnt.lt.10) go to 614
        lprint=.false.
        ztprnt=epsinv
  614   continue
c
c..compute impurity influx
c
c  zsurfs = plasma surface area in standard units
c    = <|del xi|> V'(xi) = avi(mzones,5,1) * avi(mzones,3,1) * uisl**2
c
        zsurfs = avi(mzones,5,1) * avi(mzones,3,1) * uisl**2
        zvols  = vols(mzones,1)
c
        call ncinfl ( zsurfs, zvols, zinflx )
c
c..compute the deposition profiles
c
      do 636 ii = 1,mimp0
        if ( zinflx(ii) .lt. epslon ) then
          do jx=lcentr,mzones
            siions(ii,jx) = 0.0
          enddo
        else
c
          ip = ii+lhydn
          iz = nzspec(ip)
c
c..compute the array of ionization rates zsigv  cm**3 / sec
c
        do 624 j = lcentr,ledge
          zte = tes(2,j)*zunit1
          zxn = abs(zip(iz)/(zte+epslon))
          if ( zxn .ge. 20. ) then
            zsigv(j) = 0.0
            go to 624
          endif
          zexzxn = exp(-zxn)
          z0 = log10(.25*zxn)
          zz0 = z0*z0
          zfy = zz0*(-.01505*zz0-.0459*z0+1.074)+.046*z0+.23
          zsigv(j) = sqrt(zte)*abs(zfy)*zexzxn*(1.-zexzxn)*zion0(ii)
c
c  where zion0(ii) as used here contains the factor 3.44e-11
c
          if ( cfutz(iajust) .gt. epslon )
     &      zsigv(j) = 1.7*(1.3-zexzxn)*zsigv(j)
c
c  cfutz(iajust) allows an adjustment of ionization rates so as
c    to agree (particularly for ne and a) with values quoted by
c    freeman and jones (1974).
c
          zsigv(j) = max(zsigv(j),0.0)
  624   continue
c
c..compute the neutral impurity deposition profile
c
        call neudep ( znergy(ii), ip, zslnta(ii), zpatha(ii), jsepx,
     &              zstart(ii), zsigv, ziions, ztotal(ii) )
c
c..compute impurity source term
c
c  The normalization factor znorm includes factors
c  zsurfs = surface area of zone boundary
c             at which impurities are injected
c  2.*vols(mzones,1) = 2. * total volume of plasma
c  since dx2i(j) used above is d(volume) / (2.*(total volume of plasma))
c
        jmax0 = ledge - int(zstart(ii)+0.1)
        zsurfs = avi(jmax0+1,3,1) * avi(jmax0+1,5,1) * uisl**2
        znorm = zinflx(ii) * zsurfs
     &         / (ztotal(ii) * 2. * vols(jmax0,1) + epslon)
        do j=lcentr,jmax0
          siions(ii,j) = znorm * ziions(j)
        enddo
c
        endif
c
  636 continue
c
c..The impurity-ion radiation is increased over that of coronal
c  equilibrium in two ways:  (1) an empirical non-equilibrium
c  enhancement factor, zhance(ii); and, (2) an empirical
c  initial ionization loss, zeloss(ii)
c
        do 642 ii=1,mimp0
c
          if ( zinflx(ii) .lt. epslon ) then
            do j=lcentr,mzones
              weirs(ii,j) = zhance(ii) * weirs(ii,j)
            enddo
          else
c
            jmax0 = ledge - int(zstart(ii) + 0.1)
            zsurfs = avi(jmax0+1,3,1) * avi(jmax0+1,5,1) * uisl**2
            znorm = zinflx(ii) * zsurfs
     &         / (ztotal(ii) * 2. * vols(jmax0,1) + epslon)
c
            zcronl(ii) = 0.0
            do j=lcentr,jmax0
              zcronl(ii) = weirs(ii,j) * dx2i(j) + zcronl(ii)
            enddo
c
            do 638 j=lcentr,mzones
c
              z0 = 1.0 + abs( zxloss(ii)*znorm*ztotal(ii) /
     &                       (zcronl(ii)+epslon) )
              zhance(ii) = max(z0,zhance(ii))
              z1 = mhance(ii)
              zhance(ii) = min(z1,zhance(ii))
              weirs(ii,j) = zeloss(ii) * siions(ii,j)
     &                    + zhance(ii) * weirs(ii,j)
  638       continue
c
          endif
c
  642   continue
c
c..diagnostic printout
c
c
cbate        if ( lprint ) then
cbate          write(nprint,6016) 'impflx-1'
cbate          write(nprint,6018) ((siions(ii,kj),ii=1,2),kj=1,mzones)
cbate        endif
c
cbate        if ( lprint ) write(nprint,6010) zhance(1),zhance(2)
cbate 6010 format(/2x,
cbate     &' non-equilibrium rad-loss enhancement factor for impurity 1 ='
cbate     & ,1pe10.3/2x,
cbate     &' non-equilibrium rad-loss enhancement factor for impurity 2 ='
cbate     & ,1pe10.3/)
cbate 6016 format(4x,a)
cbate 6018 format((2x,5(2x,1pe11.3,1pe11.3)))
c
  700   continue
c
c  7)   Helium recycling   ( only if CFUTZ(70) > EPSLON )
c       ----------------
c
c       The outer boundary condition normally used with he-4 recycling
c       is that of a pedestal value or pedestal adjusted to prevent
c       a rising density at the "wall" (i.e., nbound=0 or =1).
c
c       Alternatively, a zero he-4 boundary flux can be imposed by
c       following the procedure explained in subroutine bounds, and
c       this involves setting some element of flimp(indxhe,it) slightly
c       greater than epslon at some time timp(it), where "ii=indxhe"
c       is the impurity-ion index for he-4.
c
c
        if(cfutz(iheflx).le.epslon) go to 800
c
        if(.not.inita2) go to 702
        zrmins=rmins
        zrmajs=rmajs
  702   continue
        if(rmins.ne.zrmins) inita2=.true.
        if(rmajs.ne.zrmajs) inita2=.true.
c
        if ( .not. inita2 .and. cfutz(490) .le. epslon ) go to 710
c
c..initializations for helium recycling
c
        inita2=.false.
        zrmins=rmins
        zrmajs=rmajs
        indxhe=0
        do 704 ji=1,mimp
c
c   les jan-91; 4He has nimp=2 or -5 in 1 1/2D BALDUR
c
        if (cfutz(490).le.epslon .and. nimp(ji).eq.2 ) indxhe=ji
        if (cfutz(490).gt.epslon .and. nimp(ji).eq.-5) indxhe=ji
  704   continue
        if(indxhe.eq.0) go to 800
c
c   les  nov-90; for d-3he, repeat with helium-3
c
        istep=0
 705  continue
c
        call neuset(cfutz,epslon,epsinv,2)
c
        zhelos=cfutz(ionlos)*1.0e+03*evs
        iihe0=indxhe+lhydn
        zfactr = 2.0 * vols(mzones,1)
        zfact1 = 2.0 * vols(mzones,1) * usil**3
        jsepx=mzones
        if(nadump(1).gt.lcentr) jsepx=nadump(1)
        zfact2 = avi(jsepx,5,1) * avi(jsepx,3,1) ! surface area
        zslant = min(cfutz(islant),1.0)
        zpath = cfutz(inpath)*rmins/ahalfs(jsepx,1)
c
c 14.04 Standardize v calculation
c
        if (versno.gt.14.03) then
          zhenrg = cfutz(itmhe0)
        else
          zhenrg = cfutz(itmhe0) * fcau / fcmp
        end if
c
  710   continue
c
c  Move volume and area factors out of initialization section
c
        if (versno.gt.14.03) then
          zfactr = 2.0 * vols(mzones,1)
          zfact1 = 2.0 * vols(mzones,1) * usil**3
          zfact2 = avi(jsepx,5,1) * avi(jsepx,3,1) ! surface area
        end if
c
        do 712 jz=lcentr,mzones
          if ( zinflx(indxhe) .lt. epslon ) siions(indxhe,jz)=0.0
  712   continue
c
c  The total neutral helium influx is assembled from the
c  (1) cross-field diffusion at the outer edge of the plasma, and
c  (2) the flow along field lines into the divertor or limiters.
c
c..total cross-field outflux and influx
c
        zheout=max(-fflxoi(iihe0),0.0)*usit
        zhein1=cfutz(icfout)*zheout
c
c..total inflow from divertor or limiters
c
        zhein2=0.0
        do 720 j=jsepx,mzones
        zhein2=zhein2-scroff(iihe0,j)*rhois(indxhe,2,j)*dx2i(j)
  720   continue
        zhein2=cfutz(idvout)*zfactr*zhein2
c
        zinflx(indxhe)=zhein1+zhein2
c
c..outflux of helium ions across the scrapeoff boundary (separatrix)
c
        zflx1=0.0
        do 721 ip=1,mchi
        zflx1=zflx1+aaaa(iihe0,ip,jsepx)*chi(ip,jsepx-1)+
     .  bbbb(iihe0,ip,jsepx)*chi(ip,jsepx)
  721   continue
        zflx1=-zfact2*zflx1*usit
c
c
c  The neutral helium density zimpnu(j) as a function of position
c  is defined at zone boundaries and normalized to 1.0 on the
c  outer edge of the scrapeoff region, r=rmins.
c
c..calculation of the helium ionization rate
c
        do 750 j=lcentr,ledge
        zte=tes(2,j)*evsinv
        if(zte.gt.1.0e+04) go to 740
        if(zte.gt.2.0) go to 730
        zsigv(j)=0.0
        go to 750
  730   continue
        x0=log(zte)
        y0=zahe0(1)+x0*(zahe0(2)+x0*(zahe0(3)+x0*(zahe0(4)+
     .  x0*(zahe0(5)+x0*(zahe0(6)+x0*zahe0(7))))))
        zsigv(j)=exp(y0)
        go to 750
  740   continue
        x0=log(zte)
        y0=x0*zahe0(8)+zahe0(9)
        zsigv(j)=exp(y0)
  750   continue
c
        call neudep(zhenrg,iihe0,zslant,zpath,jsepx,
     1              cfutz(istart)-1.,zsigv,ziions,ztotal(indxhe))
c
        if(nstep.lt.3) call rarray('heions-1',ziions,mzones)
c
c..renormalization so as to conserve total helium influx
c
        znorm=zinflx(indxhe)/(zfactr*ztotal(indxhe)+epslon)
        do 790 j=lcentr,mzones
        siions(indxhe,j)=znorm*ziions(j)+siions(indxhe,j)
        ziions(j)=siions(indxhe,j)
c
c..electron-energy loss rate due to initial helium ionization
c
        weirs(indxhe,j)=weirs(indxhe,j)+zhelos*siions(indxhe,j)
c
  790   continue
c
        if(nstep.lt.3) call rarray('siions-2',ziions,mzones)
c
c   les  nov-90;  for d-3he repeat with helium-3
c
      istep=istep+1
      if ( cfutz(490) .lt. epslon  .or.  istep.ne.1 ) go to 800
        indxhe=lhe3-lhydn
        go to 705
c
  800   continue
c
c  9)   Smooth WEIRS   ( only if CFUTZ(199) > EPSLON )
C       ------------
c
        if(mimp.le.0) go to 1060
c
c..15.08 Replace smoothing with calls to SMOOTH. Note that the order of
c        the smoothing is as before [cfutz(ismrad) = 1 -> smooth over "5
c        zones"], but has been generalized to allow any order of smoothing.
c
      if (cfutz(ismrad).gt.epslon) then
c
        iorder = int(cfutz(ismrad) + 1.01)
        do 950 ii=1,mimp
          do 910 j=lcentr,mzones
            zsmoo1(j) = weirs(ii,j)
  910     continue
c
          call smooth(zsmoo1,zsmoo2,mzones,2,2,0,lcentr,dx2i,0,0,iorder)
c
c..Repeat as in original scheme
c
          call smooth(zsmoo2,zsmoo1,mzones,2,2,0,lcentr,dx2i,0,0,iorder)
c
          do 920 j=lcentr,mzones
            weirs(ii,j) = zsmoo1(j)
  920     continue
  950   continue
c
      end if
c
 1060   continue
        return
c
c**********************************************************************
c
        entry impprt
c
c .. short printout line (called from sprint)
c
        zv0 = sqrt((2.*zhenrg*1.e3*evs)
     1         /(aspec(iihe0)*fcau + epslon))*10.0**(-0.5*fxau)
        write(nprint,1500) zhein1,zhein2,zv0,zflx1
 1500   format(4x,' he-in1=',1pe10.3,2x,' he-in2=',1pe10.3,2x,
     &  ' vel-he=',1pe10.3,2x,' he outflx into scroff=',1pe10.3/)
c
        return
c
        entry fdback
c
c .. short printout of "zadjst" (called from sprint)
c
cbate        write(nprint,1510) zadjst
cbate 1510   format(4x,' neu-imp influx readjustment factor=',1pe10.3/)
c
        return
c
        end

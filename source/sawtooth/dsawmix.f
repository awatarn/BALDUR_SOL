c  12:00 25-aug-04 .../baldur/source/sawtooth/dsawmix.f
c  BALDUR  file DSAWMIX by Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c
c  File DSAWMIX consists of the following subprograms:
c
c sawmix : sawtooth mixing model implemented by Bateman
c sawtst : test for sawtooth crash:  moved to file  DSAWTST.TEX
c recon1 : Kadomtsev magnetic reconnection algorithm
c fpoly1 : polynomial evaluation for zero finder in sbrtn recon1
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@dsawmix  /11040/bald91/wbaldn1 DSAWMIX
c  rgb 02-aug-04 Implemented Porcelli's partial magnetic reconnection
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rap 15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 18:00 11-nov-93 commented out calls to sawmx0 and sawmx1
c  rgb 20.13 14:00 06-oct-91 removed call to sawmx1 (Seidl's sawmix)
c  rgb 20.12 21:00 05-sep-91 changed calls to cubint consistent with
c      imslmath lib
c  dps 07-jul-89 15.11 add neutral impurities to IRE; not flattened by sawteeth
c  dps 15-may-89 15.09 remove ncskip switch from IRE section
c  dps 25-aug-88 15.00 flatten impurity charge state densities for use
c                with NC code.
c  rgb 14.03 29-feb-88 corrected ziota(j) and hence bpoli(j)
c             when partial reconnection model is used
c  rgb 14.01 16-feb-88  partial reconnection model controlled by csawth(1)
c  dps 17-dec-87  allow for multiple beam species
c  rgb 06-jul-87  lsawth(3)=max(1,lsawth(3)) to prevent zero divide
c  rgb 12-jun-87  moved counter increment after test for sawtooth crash
c  rgb 09-jun-87  implemented counter ISAWCT to control frequency of output
c  rgb 26-mar-87  set IMIXED=0 if IHELEV .LT. 4 in sbrtn recon1
c  rgb 25-mar-87  implemented IXINN, IXOUT in sbrtn recon1
c       established starting point in middle of phelev array if necessary
c  rgb 23-mar-87  on input to sbrtn recon1, imixed controls amount of
c       printout from sbrtn recon1:  imixed > 0 for some output,
c       imixed > 4 for full output
c  rgb 6-mar-87 set up cubic spline interpolation in sbrtn recon1
c       for diagnostic printout of coeffs
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@sawmix  .../baldur/code/bald/dsawmix.f
c  rgb 24-nov-99 commented out calls to grafix
c  rgb 16-jun-98 swton --> swton(1)
c  rgb 02-jun-96 initialize and save ztemax, ztemin, znemax, ...
c    data ztlast /0.0/ and save ztlast
c  rgb 18.75 12:00 06-nov-90 changed kjmix from 300 to 900
c  rgb 18:00 17-dec-89 changed do 62 jz=ijz1+1,ibout-2
c                           to do 62 jz=ibinn+1,ibout-2
c  rgb 11:00 18-mar-87 determined ihmnmx from zch(j+1,1)-zch(j,1)
c            rather than from zch(j,2) in sbrtn recon1
c  rgb 10:00 18-mar-87  smoothed chi at edges of mix reg if lsawth(7)>0
c  rgb 18:00 17-mar-87  changed over from quadratic to cubic interpolation
c      in sbrtn recon1 when it became clear that the quadratic
c      interpolation was ringing too much to be useful
c  rgb 11:50 25-feb-87  temp diag output added to sbtrn recon1
c      protection added against itrlev<ihelev etc
c  rgb 25-feb-87 changed zdhdt1 at beginning of sbrtn recon1
c                to reduce ringing of the interpolating polynomial
c  rgb 22-feb-87  differential mixing option for chi(jchi,j)
c  rgb 12-feb-87 implemented completely rewritten version of sbrtns
c       sawmix, sawtst, recon1
c  rgb 01-dec-86 sbrtn grafix(2) called bef and aft crash if lsawth(5)>0
c       Note: NSTEP is increased by 1 by this sbrtn when lsawth(5)>0
c               tai and tbi are increased by 2.*dtmini
c  rgb 30-nov-86 increased long printout, controlled by lsawth(2)
c  rgb 24-nov-86 initialize tbsaws= -1. / epslon
c       new variables lsawth, swton, swtoff, swperd, swqmin, swxmin
c       in common /comsw1/ in cliche clsaw in file dcommon
c  rgb 23-feb-86 implemented nlomt2(24)=.t. to bypass sbrtn sawmix
c                output sawrmx, sawri, sawr1 at end of sbrtn
c  rgb 23-feb-86 implemented nlomt2(24)=.t. to bypass sbrtn sawmix
c                output sawrmx, sawri, sawr1 at end of sbrtn
c  rgb 16-oct-85 initialized tbsaws to - tmaxi*uist
c  rgb 16-oct-85 initialized tbsaws to - tmaxi*uist
c       fgps 24-may-85 modified to be compatible with 1-1/2-d codes.
c
c
c               %%%%%%%%%%%%%
c               %
c               awmix          %
c               %
c               %%%%%%%%%%%%%
c
c
        subroutine sawmix
c
c
cinput
c
c       sawtooth mixing to model internal distruptions, sbrtn SAWMIX
c       ----------------------------------------------
c
c  variables in the third (equilibrium) namelist:
c
c    (note: cfuta(460),... in the first BALDUR namelist
c                          can be used as noted)
c
c  swton  = sawtooth time on [sec]     <-- used to be cfutz(460)
c  swtoff = sawtooth time off [sec]    <-- used to be cfutz(461)
c  swperd = sawtooth period [sec]      <-- used to be cfutz(462)
c
c  swqmin = q-value at magnetic axis needed to trigger sawtooth crash
c                                      <-- used to be cfutz(477)
c
c  swxmin = minimum xbouni(jz) to trigger sawteeth
c           for q < 1., in the case of prescribed period
c           for q < swqmin, in the case of q-value trigger
c
c  lsawth(i), i=1,32, integer-valued switches for the sawtooth model:
c
c  lsawth(1) = 10  for the old (flat-top) Kadomtsev reconnection model
c            = 11  for the differential Kadomtsev model by Seidl (Nov,85)
c                  (removed 6-oct-91)
c            = other, new sawtooth reconnection model
c
c  lsawth(2) controls the amount of long printout each sawtooth crash
c       = 0  for no long printout  (default)
c       = 1  for 1 page of long printout
c       = 2  for some detailed profiles ( q(r) )
c       > 4  for all diagnostic output
c
c  Note: set nlpomt(9)=-1 in the first BALDUR namelist
c       in order to omit long sawtooth printout in sbrtn mprint
c       and just use the long printout contained in sbrtn sawmix.
c
c   Long printout consists of:
c
c  #            = zone boundary number
c  half         = halfwidth [cm] of zone boundary
c  q-pre, q-aft = q-value before and after crash at zone bndries
c  J-pre, J-aft = flux-surface-averaged current density [kA/cm**2]
c  Te-pre, Te-aft       = electron temperature [keV]
c  ne-pre, ne-aft       = electron density [electrons/cm**3]
c  Ti-pre, Ti-aft       = ion temperature [keV]
c  tor-flux    = toroidal magnetic flux [tesla*m**2]
c  helc-pre    = helical magnetic flux before sawtooh crash
c  helcl-aft   = helical magnetic flux after  sawtooh crash
c
c
c  lsawth(3)  controls frequency of long printouts
c       = 0 (default) or 1 to get long printout each time
c                          sbtrn dsawmix is called
c       = n to get long printout every n-th time sbrtn sawmix is called
c
c  lsawth(5) controls the plot output before and after sawtooth crash
c       = 0  for no special calls to sbrtn grafix(2)  (default)
c       = 1  for call to sbrtn grafix(2) before and after crash
c
c  lsawth(6) controls the effect of mixing on chi(jchi,j)
c       = 0  for differential mixing within mixing region
c       = other for the old flat top model for chi(jchi,jz)
c
c  lsawth(7) controls the amount of smoothing
c       = 0  for no smoothing at edges of mixing region
c       = 1  for simple smoothing of bpoli at edges of mixing region
c
c  lsawth(10) used to implement new models of sawtooth period
c             in sbrtn sawtst
c       = 2  for the Park-Monticello sawtooth period PPPL-2601 (1990)
c
c
c  csawth(j), j=1,32, real-valued constants for use in the sawtooth model:
c
c  csawth(1) = fraction of old helical flux added to new helical flux
c              This implements a simple model for partial reconnection
c              After computing the new helical flux with the Kadomtsev
c              model, take a linear combination of this new helical flux
c              with the old helical flux (from just before the crash)
c              new helical flux = csawth(1) * old helical flux
c                               + (1.-csawth(1)) * new helical flux
c              default: csawth(1) = 0.  --> the original Kadomtsev model
c                       csawth(1) = 1.  --> no magnetic reconnection
c
c  csawth(2) = 1.0 - magnetic reconnection fraction using Porcelli's
c              partial magnetic reconnection model
c              from F. Porcelli, D. Boucher, and M.N. Rosenbluth
c              Plasma Phys. Cont. Fusion 38 (1996) 2163-2186
c
c  csawth(5) = fraction of poloidal magnetic energy transferred
c              to the thermal electron energy
c
c
c..Variables in the first BALDUR namelist:
c
c  nlomt2(24) = .t. to bypass sawmix subroutines completely
c
c  nlpomt(9)=-1 (true) to omit long printouts from sbrtn mprint
c           = 0 (false) for long printouts from sbrtn mprint
c                       (it is better to set nlpomt(9)=-1
c                        and lsawth(2) > 0 for long printouts)
c
c  cfutz(468)  determines what quantities are mixed by the sawtooth model
c      = 0 all of the following are mixed:
c      divisible by 2 --> electron energy mixed
c      divisible by 3  --> ion energy mixed
c      divisible by 5  --> hydrogenic densities are mixed
c      divisible by 7  --> impurity densities are mixed
c      divisible by 11 --> fast ion densities are mixed
c      divisible by 13 --> fast alpha particles mixed if NTYPE=1
c
c
cend
c
c  tbsaws = time at last sawtooth crash [sec]
c  sawr1  = minor radius to outermost q=1 surface [cm]
c
c  zrpsi(j) = poloidal - toroidal flux on fine grid
c  zxpsi(j) = x-values on fine grid
c
c**********************************************************************c
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cfokkr.m'
      include 'calpht.m'
      include 'commhd.m'
      include 'clsaw.m'
      include 'comncr.m'
c
        common/comsaw/ tbsaws,sawte,sawne,sawdd,sawate,sawane,sawadd,
     1  sawtau,sawr1,sawri,sawrmx,sawqc,tmguir,njzq1,
     2  sawfus,sawalf,sawuth,estpol,uthavg,tauein,ohptot,
     3  tepre(55),tepost(55),teavg(55),tipre(55),tipost(55),tiavg(55),
     4  ajzpre(55),ajzpst(55),ajzavg(55),ohravg(55)
c
        dimension zbpold(55),zimix(2),zomix(2),
     1  ibjzi(2),ibjzo(2),ibl(2),ibu(2)
c
c..zxpsi, zrpsi, zxl, zxr, zxf, zrfpsi should be dimensioned
c       to at least iint + 5
c
        parameter (kjmix=900, kcfmix=5*kjmix)
c
c..KJMIX is used in dimensions large enough to accommodate the array
c    of toroidal fluxes associated with the helical flux levels
c    within the mixing region
c..KCFMIX is used for dimensions of interpolating coefficients to be
c    used when interpolating arrays dimensioned with KJMIX
c
c        --------------------------------------------------
c
      dimension ztori(kjbal), zpoli(kjbal), zheli(kjbal), zvoli(kjbal)
     & , ziota(kjbal) 
     & , zcf(kcfmix)  , zhelev(kjmix), ztrlev(kjmix), iheltr(kjmix)
     & , zdtnew(kjmix), ztfnew(kjmix), zpfnew(kjmix)
     & , zvoldi(kjmix), zvnewi(kjmix), ztemp(kjmix)
     & , zchtot(kjbal), zchlev(kjmix), zchnew(kjmix)
c
c
c..ZTORI(J) = toroidal flux at BALDUR zone boundaries (MKS units)
c..ZPOLI(J) = poloidal flux at BALDUR zone boundaries (MKS units)
c..ZHELI(J) = helical  flux at BALDUR zone boundaries (MKS units)
c
c               ZHELI(J) = ZQMN * ZPOLI(J) - ZTORI(J)
c
c       Note:  These arrays are updated after each reconnection.
c
c..ZQMN = m / n = q-value at mode rational surface = 1.0 here
c
c..ZCF(J) ... coefficients used in spline fits of any arrays here
c
c  Returned after reconnection from the call to recon1:
c
c..ZHELEV(J) = helical flux levels in the mixing region
c                  sorted in descending order, J=1,IHELEV
c..ZTRLEV(J) = toroidal flux corresponding to each value of ZHELEV
c                  in the mixing region, sorted in ascending order,
c                  J=1,ITRLEV > IHELEV
c..IHELTR(J) = index of each helical flux ZHELEV(IHELTR(J))
c                  correspoding to each toroidal flux ZTRLEV(J), J=1,ITRLEV
c                  1 =< IHELEV(J) =< IHELEV
c..IBINN        = index of BALDUR zone boundary 
c                 at inner edge of mixing region,
c               = 2 if the mixing region includes the magnetic axis
c..IBOUT        = index of BALDUR zone boundary 
c                 beyond outer edge of mixing region
c..IMIXED   = number of separate mixing regions (see sbrtn RECON1)
c
c
c  Computed further in this sbrtn:
c
c..ZDTNEW(J) = toroidal flux difference between ZHELEV(J) and ZHELEV(J+1)
c                  after reconnection, J=1,IHELEV-1
c..ZTFNEW(J) = toroidal flux up to ZHELEV(J) after reconnection, 
c                  J=1,IHELEV
c..ZPFNEW(J) = poloidal flux corresponding to ZTFNEW(J) 
c                  after reconnection, J=1,IHELEV
c..ZVOLDI(J) = volume associated with toroidal flux ZTRLEV(J)
c                  before reconnection, J=1,ITRLEV
c..ZVNEWI(J) = volume associated with toroidal flux ZTFNEW(J)
c                  after reconnection, J=1,IHELEV
c..ZTEMP(J)  = temporary array
c
c..ZCHTOT(J)  = volume integral of chi(i,j) out to BALDUR zone boundaries
c
c..ZCHLEV(J)  = volume integral of chi(i,j) out to ZTRLEV(J) in mix region
c
c..ZCHNEW(J)  = volume integral of chi(i,j) out to ZTFNEW(J) in mix region
c        --------------------------------------------------
c
      dimension  zbtori(kjbal), ztheta(kjbal)
     &  ,zhfpre(kjbal), zqpre(kjbal), znepre(kjbal)
c
c
c  ZBTORI(j) = d (toroidal flux) / d rho  / (2 pi r0ref)
c            = b0ref * rho / r0ref
c              this definition is the toroidal analog of bpoli(j)
c  ZTHETA(j) = r0ref**2 * V'(x) * (d rho / d x)**2 < (del x / R)**2>
c              / ( 4 * emu0 )   at BALDUR zone boundaries
c              to correspond to Seidl's gtheta(j)
c              used to compute poloidal magnetic energy
c  ZHFPRE(J) = helical flux at BALDUR zone boundaries (MKS units)
c                  before reconnection
c        --------------------------------------------------
c
      logical  inital
      data inital /.true./
c
        data iqlimt/477/
        data itimon,itimof,itaus,itim2,itaus2/ 460, 461, 462, 463, 464/
        data itmgr1,itmgr2,itmgr3/ 465, 466, 467/
        data ipmix/468/
        data isawct/0/
        data ztlast /0.0/
c
c  ISAWCT = number of times sbrtn sawmix has been called
c           used to control frequency of long printouts
c
      data  iclass /2/,  isub /24/
c
      save ztemax, ztemin, znemax, znemin, zddmax, zddmin
     &   , zfsmax, zfsmin, zalmax, zalmin, zutmax, zutmin
     &   , ztlast
c
c**********************************************************************c
c
c..choose sawtooth model based on lsawth(1)
c
cbate      if (lsawth(1) .eq. 10) then
cbate        call sawmx0
cbate        return
cbate      elseif (lsawth(1) .eq. 11) then
cbate        call sawmx1
cbate        return
cbate      endif
c
c..bypass sbrtn sawmix if nlomt2(24) = .true.
c
      if (nlomt2(isub)) then
        if (nstep .lt. 2)
     & call mesage(' *** 2.24 subroutine sawmix bypassed ')
        return
      endif
c
c..to ensure backward compatibility
c
        lsawth(3) = max ( 1, lsawth(3) )
        if (cfutz(460) .gt. epslon) swton(1)  = cfutz(460)
        if (cfutz(461) .gt. epslon) swtoff = cfutz(461)
        if (cfutz(462) .gt. epslon) swperd = cfutz(462)
        if (cfutz(477) .gt. epslon) swqmin = cfutz(477)
c
c..initialize
c  ----------
c
      zusifl = usib * usil**2  ! mag flux from standard to internal units
c
        if(inital) then
          write (6,*) ' sbrtn SAWMIX sawtooth mixing by Bateman used'
          njzq1=0
          tmguir=0.
          tbsaws= - 1. / epslon
          sawtau=0.
          ztau=1.
          call sawavg(2,teavg,tiavg,ajzavg,ohravg,ohptot,ztau)
        end if
c
        if(nstep.eq.0) return
c
        call sawavg(1,teavg,tiavg,ajzavg,ohravg,ohptot,ztau)
        call sawamp(zte,zne,zdd,zfus,zalf,zuth)
c
        if(inital) then
          inital = .false.
          ztemax = zte
          ztemin = zte
          znemax = zne
          znemin = zne
          zddmax = zdd
          zddmin = zdd
          zfsmax = zfus
          zfsmin = zfus
          zalmax = zalf
          zalmin = zalf
          zutmax = zuth
          zutmin = zuth
        else
          ztemax = max(ztemax,zte)
          ztemin = min(ztemin,zte)
          znemax = max(znemax,zne)
          znemin = min(znemin,zne)
          zddmax = max(zddmax,zdd)
          zddmin = min(zddmin,zdd)
          zfsmax = max(zfsmax,zfus)
          zfsmin = min(zfsmin,zfus)
          zalmax = max(zalmax,zalf)
          zalmin = min(zalmin,zalf)
          zutmax = max(zutmax,zuth)
          zutmin = min(zutmin,zuth)
        end if
c
        if((nstep.eq.5*(nstep/5)).and.(tbsaws.le.0.)) then
        ztau =uist*tai-ztlast
        call sawavg(2,teavg,tiavg,ajzavg,ohravg,ohptot,ztau)
        ztlast =uist*tai
        uthavg =zuth
        if(ohptot.ne.0.) tauein =uthavg/ohptot
        end if
c
c..Test for sawtooth crash
c  -----------------------
c
        call sawtst (iswtst)
c
        if (iswtst .lt. 1) return
c
c..Increment counter
c
      isawct = isawct + 1
c
c
c..Set up before sawtooth mixing
c  -----------------------------
c
c
c..plot profile before crash
c
c  Note: we have to ensure that the time given to sbrtn grafix
c  does not coincide with the next time sbrtn grafix is called.
c
        if (lsawth(5) .gt. 0) then
          tai = tai + dtmini
          tbi = tbi + dtmini
cbate        call grafix(2)
        endif
c
c
c  set up temporary arrays
c  ZBTORI(j) = d (toroidal flux) / d rho  / (2 pi r0ref)
c            = b0ref * rho / r0ref
c              this definition is the toroidal analog of bpoli(j)
c  ZTHETA(j) = r0ref**2 * V'(x) * (d rho / d x)**2 < (del x / R)**2>
c              / ( 4 * emu0 )   at BALDUR zone boundaries
c              to correspond to Seidl's gtheta(j)
c              used to compute poloidal magnetic energy
c
c  Toroidal and poloidal flux computed from FLTORS(J,1) and FLPOLS(J,1)
c  which were computed in sbrtn MHDBAL on BALDUR zone boundaries.
c  This may need to be changed.
c
      do 10 j=1,mzones
        zbtori(j) = b0ref * avi(j,1,1) / r0ref
        ztheta(j) = r0ref**2 * avi(j,3,1) * avi(j,2,1)**2
     &             * avi(j,7,1) / ( 4. * emu0 )
          ztori(j) = fltors(j,1) * zusifl
          zpoli(j) = flpols(j,1) * zusifl
          zvoli(j) = avi(j,12,1)
        zbpold(j)= bpoli(j)
  10  continue
c
      zbtori(1) = - zbtori(3)
      ztori(1)  = - ztori(3)
      zpoli(1)  = - zpoli(3)
      zvoli(1)  = - zvoli(3)
      zbpold(1) = - zbpold(3)
c
c       get poloidal energy before re-distribution
c
        zepolo=0.
      do 12 jz=lcentr,ledge
        zepolo = zepolo +
     &  (ztheta(jz)+ztheta(jz+1))*bpoli(jz)*bpoli(jz+1)*dxzoni(jz)
  12  continue
c
c
c..Sawtooth Mixing
c  ---------------
c
  20  continue
c
c
c..Compute helical flux ZHELI(J) on BALDUR zone boundaries
c  Here ZQMN = 1.0 = for standard m/n=1 sawtooth reconnection
c
        zqmn = 1.0
c
        do 22 j=1,mzones
          zheli(j)  = zqmn * zpoli(j) - ztori(j)
          zhfpre(j) = zqmn * zpoli(j) - ztori(j)
  22  continue
c
c..Compute Kadomtsev magnetic reconnection
c  Note: On input, ITRLEV = maximum allowed value for itrlev
c                  IMIXED controls amount of printout from sbrtn recon1
c
        imixed = lsawth(2)
        itrlev = kjmix
c
!cap
        ztrlev = 0.
!
        call recon1 (ztori,zheli,mzones,epslon
     & ,zhelev,ihelev,ztrlev,itrlev,iheltr,ibinn,ibout,imixed)
c
c..Check to see if there has been a reconnection
c
        if ( imixed .lt. 1 ) go to 70
c
c..Build up toroidal flux from the inner edge of the mixing region
c
        call resetr (zdtnew,itrlev,0.0)
        call resetr (ztfnew,itrlev,0.0)
c
        do 24 j=1,itrlev-1
c
          jh = min ( iheltr(j), iheltr(j+1) )
          zdtnew (jh) = zdtnew (jh) + ztrlev (j+1) - ztrlev (j)
c
  24  continue
c
        ztfnew (1) = ztrlev (1)
        do 26 j=1,ihelev-1
c
          ztfnew (j+1) = ztfnew (j) + zdtnew (j)
c
  26  continue
c
c..Diagnostic output
c
      if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 4) then
cbate        call prtr1d (6,ztori,mzones,'ztori')
cbate        call prtr1d (6,zheli,mzones,'zheli before reconnection')
cbate        call prtr1d (6,zhelev,ihelev,'zhelev(j),j=1,ihelev')
cbate        call prtr1d (6,ztrlev,itrlev,'ztrlev(j),j=1,itrlev')
cbate        call prtr1d (6,ztfnew,ihelev,'ztfnew(j),j=1,ihelev')
c
c..diagnostic printout
c
         if ( lsawth(2) .gt. 8 ) then
c
           write(nout,*)
           write(nout,*) 'Fundamental arrays before sawtooth crash'
           write (nout,950)
           do j=1,mzones
             write (nout,951) chi(1,j), chi(2,j), chi(3,j), chi(4,j)
     &         , chi(lelec,j), chi(lion,j)
           enddo
c
         endif
      endif
c
c..Compute new helical flux on BALDUR zone bndries within mixing region
c  ZHELI(J) = updated helical flux after reconnection
c
      ib1 = ibinn + 1
      if (ztrlev(1) .lt. 0.5*(ztori(2)+ztori(3)) ) ib1 = 2
      ib2 = ibout - 1
      if (ztrlev(itrlev) .gt. 0.5*(ztori(ibout-1)+ztori(ibout)) )
     &     ib2 = ibout + 1
c
        call cubint (ztfnew,zhelev,ihelev,0,zcf,kjmix
     &  ,ztori(ib1),zheli(ib1),ib2-ib1+1,0,0.,0
     &  ,'new zheli in mixing region on BALDUR grid in sbrtn sawmix')
c
c..Pocelli's model for partial magnetic reconnection
c  which is controlled by csawth(2)
c
      if ( csawth(2) .gt. epslon
     &     .and. csawth(2) .lt. 1.0 + epslon) then
c
c..first, compute the normalized radial extent of the mixing
c  use local linear interpolation of zhfpre(j)
c  assume zhfpre = 0. at the edges of the mixing
c  zxinn = xbouni(j) at inner edge of mixing region
c  zxout = xbouni(j) at outer edge of mixing region
c  iq1   = index just outside the q = 1.0 surface (max zhfpre(j))
c
        zxinn = ( zhfpre(ibinn+1) * xbouni(ibinn)
     &          - zhfpre(ibinn) * xbouni(ibinn+1) )
     &          / ( zhfpre(ibinn+1) - zhfpre(ibinn) )
c
        zxout = ( zhfpre(ibout+1) * xbouni(ibout)
     &          - zhfpre(ibout) * xbouni(ibout+1) )
     &          / ( zhfpre(ibout+1) - zhfpre(ibout) )
c
        do j=ibinn,ibout
          iq1 = j
          if ( zhfpre(j+1) .lt. zhfpre(j) ) exit
        enddo
c
c..compute izbinn and izbout, which are
c  the indicies at the left edge of the zones
c  for which zbout - zbinn = csawth(2) * ( zxout - zxinn )
c  xbouni(izbinn) < zbinn < xbouni(izbinn+1)
c  xbouni(izbout) < zbout < xbouni(izbout+1)
c
        izbinn = ibinn
        izbout = ibout
        do j=ibinn,iq1-1
          izbinn = j
          do j1=ibout,iq1,-1
            izbout = j1
            if ( zhfpre(j1) .gt. zhfpre(j) ) exit
          enddo
          if ( xbouni(izbout) - xbouni(izbinn+1)
     &         .lt. csawth(2) * ( zxout - zxinn ) ) exit
        enddo
c
c..Now use linear interpolation within each of these two zones
c  to solve for the helical flux level zflevel
c  for which zbout - zbinn = csawth(2) * ( zxout - zxinn )
c
       zflevel = (
     & (zhfpre(izbinn+1)*xbouni(izbinn)-zhfpre(izbinn)*xbouni(izbinn+1))
     &  /(zhfpre(izbinn+1)-zhfpre(izbinn)) -
     & (zhfpre(izbout+1)*xbouni(izbout)-zhfpre(izbout)*xbouni(izbout+1))
     &  /(zhfpre(izbout+1)-zhfpre(izbout))
     &  + csawth(2) * ( zxout - zxinn ) ) / ( 
     & (xbouni(izbout+1)-xbouni(izbout))
     &    /(zhfpre(izbout+1)-zhfpre(izbout))
     &  - (xbouni(izbinn+1)-xbouni(izbinn))
     &    /(zhfpre(izbinn+1)-zhfpre(izbinn)) )
c
c..Now compute zbinn and zbout, the range of magnetic reconnection
c
        zbinn = ( zflevel * ( xbouni(izbinn+1) - xbouni(izbinn) )
     & -zhfpre(izbinn)*xbouni(izbinn+1)+zhfpre(izbinn+1)*xbouni(izbinn)
     &  ) / ( zhfpre(izbinn+1)-zhfpre(izbinn) )
c
        zbout = ( zflevel * ( xbouni(izbout+1) - xbouni(izbout) )
     & -zhfpre(izbout)*xbouni(izbout+1)+zhfpre(izbout+1)*xbouni(izbout)
     &  ) / ( zhfpre(izbout+1)-zhfpre(izbout) )
c
c
c..compute range of magnetic reconnection
c  csawth(2) is the fraction of the interval of indexes
c  from ibinn to ibout that is reconnected using the Kadomtsev model
c  Hence, the interval that is reconnected extends from zbinn to zbout
c
c        zbinn = real(ibinn) + 0.5 * csawth(2) * ( ibout - ibinn )
c        izbinn = int ( zbinn + epslon )
c        zbout = real(ibout) - 0.5 * csawth(2) * ( ibout - ibinn )
c        izbout = int ( zbout + epslon )
c
c..interpolate to find zhpinn = zhfpre and zhhinn = zheli
c  Note: zhfpre(j) is the helical flux before sawtooth mixing
c        zheli(j)  is the helical flux after sawtooth mixing
c  at zbinn which is between izbinn and izbinn+1
c  zhpinn = initial helical flux at inner edge of reconnection region
c  zhpout = initial helical flux at outer edge of reconnection region
c
        zhpinn = zflevel
c
        zhhinn = ( ( zbinn - xbouni(izbinn) ) * zheli(izbinn+1)
     &        + ( xbouni(izbinn+1) - zbinn ) * zheli(izbinn) )
     &        / ( xbouni(izbinn+1) - xbouni(izbinn) )
c
        zhpout = zflevel
c
        zhhout = ( ( zbout - xbouni(izbout) ) * zheli(izbout+1)
     &        + ( xbouni(izbout+1) - zbout ) * zheli(izbout) )
     &        / ( xbouni(izbout+1) - xbouni(izbout) )
c
      write(6,*) 'zbinn = ',zbinn,' zbout = ',zbout,
     & ' izbinn = ',izbinn,' izbout = ',izbout,' in file dsawmix.f'
      write(6,*) 'zxinn = ',zxinn,' zxout = ',zxout,' iq1 = ',iq1
      write(6,*) 'zflevel = ',zflevel
      write(6,*) 'zhpinn = ',zhpinn,' zhhinn = ',zhhinn
      write(6,*) 'zhpout = ',zhpout,' zhhout = ',zhhout
c
        do j=ibout,2,-1
c
          if ( j .gt. izbout ) then
c
            zheli(j) = zhfpre(j)
c
          elseif ( j .gt. izbinn .and. j .lt. izbout ) then
c
c            zheli(j) = zheli(j) - zhhout + zhpout
c
            zheli(j) = zflevel
c
          elseif ( j .lt. izbinn ) then
c
            zheli(j) = zhfpre(j) - zhpinn + zhhinn - zhhout + zhpout
c
          endif
c
        enddo
c
        zheli(1) = zheli(3)
c
      else
c
c..Partial reconnection model controlled by csawth(1)
c  Helical flux is taken to be a linear combination of 
c    the helical flux before and after the sawtooth crash
c    as computed by the Kadomtsev model
c
        do j=ibinn,ibout
          zheli(j) = csawth(1) * zhfpre(j) + (1.-csawth(1)) * zheli(j)
        enddo
c
      endif
c
c..Compute new poloidal flux on BALDUR zone bndries within mixing region
c  ZPOLI(J) = updated helical flux after reconnection
c
          zqmnr = 1. / zqmn
          z0 = zheli(2)
        do 32 j=ibinn,ibout
          zpoli(j) = ( ztori(j) + zheli(j) - z0 ) * zqmnr
  32  continue
          zpoli(1) = - zpoli(3)
          zpoli(2) = 0.
c
c..Construct new poloidal flux on fine grid within mixing region
c
        zqmnr = 1. / zqmn
        do 28 j=1,ihelev
          zpfnew(j) = ( ztfnew(j) + zhelev(j) - z0 ) * zqmnr
  28  continue
c
c..Compute new ziota(j) = d poloidal flux / d toroidal flux
c  on BALDUR grid within the mixing region j=ibinn+1,ibout-1
c
        call cubint (ztfnew,zpfnew,ihelev,0,zcf,kjmix
     &  ,ztori(ibinn+1),ziota(ibinn+1),ibout-ibinn-1,1,0.,0
     &  ,'new ziota in mixing region on BALDUR grid in sbrtn sawmix')
c
c..Partial reconnection model controlled by csawth(1)
c
      if ( abs (csawth(1)) .gt. epslon ) then
        zfact = b0ref / r0ref
        do 33 j=ibinn+1,ibout-1
          ziota(j) = csawth(1) * bpoli(j) / ( zfact * avi(j,1,1) )
     &               + ( 1. - csawth(1) ) * ziota(j)
  33    continue
      endif
c
c..Compute new BPOLI = (B0REF/R0REF) * RHO * ZIOTA
c  on BALDUR zone boundaries within the mixing region j=ibinn+1,ibout-1
c
        zfact = b0ref / r0ref
        do 34 j=ibinn+1,ibout-1
          bpoli(j) = zfact * avi(j,1,1) * ziota(j)
  34  continue
c
c..Smooth BPOLI at the edges of the mixing region
c  to avoid excessive current density spikes
c
c  simple smoothing
c
      if ( lsawth(7) .gt. 0 ) then
        bpoli(ibout) = 0.5 * ( bpoli(ibout-1) + bpoli(ibout+1) )
        bpoli(ibinn) = 0.5 * ( bpoli(ibinn-1) + bpoli(ibinn+1) )
      endif
c
        bpoli(2) = 0.
        bpoli(1) = - bpoli(3)
c
c..get poloidal energy after re-distribution
c
          zepoln=0.
      do 36 jz=lcentr,ledge
          zepoln = zepoln +
     &  (ztheta(jz)+ztheta(jz+1))*bpoli(jz)*bpoli(jz+1)*dxzoni(jz)
  36  continue
c
c..find the change in the poloidal magnetic energy
c
          zdbpe = 0.
       do 38 jz=ibinn,ibout
        zdbpe = zdbpe  + (ztheta(jz)+ztheta(jz+1))*
     &  (bpoli(jz)*bpoli(jz+1)-zbpold(jz)*zbpold(jz+1))*dxzoni(jz)
 38    continue
c
      if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 0)
     &   write (6,138) zdbpe, csawth(5)
 138  format (/'  Change in poloidal magnetic energy: ',1pe13.5
     & ,'     and fraction transferred to electrons: ',1pe13.5)
c
c..Diagnostic output
c
      if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 4) then
cbate        call prtr1d (6,zheli,mzones,'zheli(j), j=1,mzones, after recon')
cbate        call prtr1d (6,zpfnew,ihelev,'zpfnew(j),j=1,ihelev')
cbate        call prtr1d (6,ziota,mzones,'ziota(j),j=1,mzones')
c
cbate        call prtr1d (6,zbpold,mzones,'bpoli before')
cbate        call prtr1d (6, bpoli,mzones,'bpoli after')
cbate        call prtr1d (6,ztheta,mzones,'ztheta')
      endif
c
       eohmi = eohmi - zdbpe * csawth(5)
c
c..Flat-top model for profile flattening
c  -------------------------------------
c
c  ZVINNI and ZVOUTI are the volumes at the inner and outer edge
c       of the mixing region
c
        ztemp(1) = ztrlev(1)            ! tor flux at inner edge of mixing
        ztemp(2) = ztrlev(itrlev)       ! tor flux at outer edge of mixing
c
        call cubint (ztori,zvoli,mzones,0,zcf,kjbal
     & ,ztemp,zvnewi,2,0,0.,-1
     & ,'interpolation from tor flux to volume in sbrtn sawmix')
c
        zvinni = zvnewi(1)
        zvouti = zvnewi(2)
c
c  differential volumes just inside the edge of the mixing region
c
        zdvinn = zvoli(ibinn+1) - zvinni
        zdvout = zvouti - zvoli(ibout-1)
c
c  linear interpolation coefficients to be used below
c
        znainn = zdvinn / ( zvoli(ibinn+1) - zvoli(ibinn) )
        zncinn = 1. - znainn
c
        znaout = zdvout / ( zvoli(ibout) - zvoli(ibout-1) )
        zncout = 1. - znaout
c
c  volume within mixing region
c
        zsum0 = zvouti - zvinni
c
c.. loop over the mixing regions :: change to go to 20 below
c
c       get the average particle and energy densities
c
        do 48 jchi=1,mchi
c
c       test for mixing switch
c
        if(jchi.eq.lelec) itest=2
        if(jchi.eq.lion) itest=3
        if((jchi.ge.lhyd1).and.(jchi.le.lhydn)) itest=5
        if((jchi.ge.limp1).and.(jchi.le.limpn)) itest=7
        ipd=int(cfutz(ipmix))
        if(itest*(ipd/itest).ne.ipd) go to 48
c
c..differential model for rearranging chi
c
      if ( lsawth(6) .eq. 0 ) then
c
c  work with the total number of particles or total thermal energy
c  within each BALDUR zone boundary  == ZCHTOT(J)
c
c  For now, do not include the effects of local compression or expansion
c
c..ZCHTOT(J)  = volume integral of chi(i,j) out to BALDUR zone boundaries
c
        zchtot(2) = 0.
      do 92 j=2,mzones-1
        zchtot(j+1) = zchtot(j) + chi(jchi,j) * (zvoli(j+1)-zvoli(j))
  92  continue
        zchtot(1) = - zchtot(3)   ! odd symmetry as a fn of zvoli(j)
c
c  Interpolate from ztori,zchtot to ztrlev,zchlev within mixing region
c..ZCHLEV(J)  = volume integral of chi(i,j) out to ZTRLEV(J) in mix region
c
      call cubint (ztori,zchtot,mzones,0,zcf,kjmix
     & ,ztrlev,zchlev,itrlev,0,0.,-1
     & ,'ztori,zchtot -> ztrlev,zchlev in sbrtn sawmix')
c
c  Rearrange as with toroidal flux within mixing region
c..ZCHNEW(J)  = volume integral of chi(i,j) out to ZTFNEW(J) in mix region
c
      call resetr (ztemp,itrlev,0.0)
      call resetr (zchnew,itrlev,0.0)
c
      do 94 j=1,itrlev-1
        jh = min ( iheltr(j), iheltr(j+1) )
        ztemp (jh) = ztemp (jh) + zchlev (j+1) - zchlev (j)
  94  continue
c
        zchnew (1) = zchlev (1)
      do 96 j=1,ihelev-1
        zchnew (j+1) = zchnew (j) + ztemp (j)
  96  continue
c
c  Interpolate back to BALDUR zone boundaries
c..ZCHTOT(J)  = volume integral of chi(i,j) out to BALDUR zone boundaries
c
      call cubint (ztfnew,zchnew,ihelev,0,zcf,kjmix
     &  ,ztori(ibinn+1),zchtot(ibinn+1),ibout-ibinn-1,0,0.,0
     &  ,'new zchtot in mixing rgn on BALDUR bndries in sbrtn sawmix')
c
c..smooth the edges of the sawtooth mixing region
c
      if (lsawth(7) .gt. 0) then
        zchtot(ibout) = 0.5 * ( zchtot(ibout-1) + zchtot(ibout+1) )
        zchtot(ibinn) = 0.5 * ( zchtot(ibinn-1) + zchtot(ibinn+1) )
      endif
c
      zchtot(2) = 0.
      zchtot(1) = - zchtot(3)
c
c  Compute new chi values
c
      do 98 j=2,mzones-1
        chi(jchi,j) = (zchtot(j+1)-zchtot(j))/(zvoli(j+1)-zvoli(j))
  98  continue
        chi(jchi,1) = chi(jchi,2)
c
c..Flat-top model for chi
c
      else
        zsum = chi(jchi,ibinn) * zdvinn + chi(jchi,ibout-1) * zdvout
        if(ibinn+1.gt.ibout-2) go to 43
       do 42 jz=ibinn+1,ibout-2
          zsum = zsum + chi(jchi,jz) * (zvoli(jz+1) - zvoli(jz))
  42   continue
  43   continue
c
c..increment the energy conservation variables
c
        if(jchi.eq.lelec) then
          zsum = zsum - zdbpe * csawth(5)
          asordi(lelec) = asordi(lelec) - zdbpe * csawth(5)
        endif
c
        zavg = zsum / zsum0
c
c..flatten the densities
c
      if(ibinn+1.gt.ibout-2) go to 45
        do 44 jz=ibinn+1,ibout-2
          chi(jchi,jz) = zavg
  44    continue
  45  continue
c
        chi(jchi,ibinn)   = zavg * znainn + chi(jchi,ibinn) * zncinn
        chi(jchi,ibout-1) = zavg * znaout + chi(jchi,ibout-1) * zncout
        chi(jchi,1) = chi(jchi,2)
c
      endif ! end of rearranging chi
c
  48  continue
c
c  15.00 Flatten impurity charge state densities if present
c  15.11 Add neutrals to IRE code, but do not flatten their densities
c
      if ((natomc.eq.3).and.(mimp.ne.0)) then
c
        itest = 7
        if (itest*(ipd/itest).ne.ipd) go to 159
c
        do 158 ji=1,mimp
c
c  15.00 Prepare also to sum total impurity density
c
          xns(ibinn,ji) = xn(ibinn,0,ji)
          xns(ibout-1,ji) = xn(ibout-1,0,ji)
c
          if (ibinn+1.gt.ibout-2) go to 149
c
          do 148 jz=ibinn+1,ibout-2
            xns(jz,ji) = xn(jz,0,ji)
  148     continue
c
  149     continue
c
          do 157 ik=1,nk
c
            zsum = xn(ibinn,ik,ji)*zdvinn + xn(ibout-1,ik,ji)*zdvout
            if (ibinn+1.gt.ibout-2) go to 152
c
            do 151 jz=ibinn+1,ibout-2
              zsum = zsum + xn(jz,ik,ji)*(zvoli(jz+1) - zvoli(jz))
  151       continue
c
  152       continue
c
            zavg = zsum / zsum0
            if (ibinn+1.gt.ibout-2) go to 155
c
            do 154 jz=ibinn+1,ibout-2
              xn(jz,ik,ji) = zavg
              xns(jz,ji) = xns(jz,ji) + zavg
  154       continue
c
  155       continue
            xn(ibinn,ik,ji) = zavg*znainn
     1                             + xn(ibinn,ik,ji)*zncinn
            xn(ibout-1,ik,ji) = zavg*znaout
     1                             + xn(ibout-1,ik,ji)*zncout
            xn(1,ik,ji) = xn(2,ik,ji)
            xns(ibinn,ji) = xns(ibinn,ji) + xn(ibinn,ik,ji)
            xns(ibout-1,ji) = xns(ibout-1,ji) + xn(ibout-1,ik,ji)
c
  157     continue
c
          xns(1,ji) = xns(2,ji)
  158   continue
  159   continue
c  
      end if
c
c..flatten the fast ion density
c
        itest=11
      if(itest*(ipd/itest).ne.ipd) go to 59

        do 58 jspc=1,mhsp
        do 58 je=1,lhemax
        do 57 jmu=1,nhmu
          zsum = hfi(jmu,je,ibinn,jspc)*zdvinn 
     1           + hfi(jmu,je,ibout-1,jspc)*zdvout
c
          if(ibinn+1.gt.ibout-2) go to 53
            do 52 jz=ibinn+1,ibout-2
              zsum = zsum + hfi(jmu,je,jz,jspc) 
     1                      * (zvoli(jz+1)-zvoli(jz))
  52        continue
  53      continue
c
          zavg=zsum/zsum0
c
        if(ibinn+1.gt.ibout-2) go to 55
          do 54 jz=ibinn+1,ibout-2
            hfi(jmu,je,jz,jspc) = zavg
  54      continue
  55    continue

         hfi(jmu,je,ibinn,jspc)   = zavg*znainn 
     1                          + hfi(jmu,je,ibinn,jspc)*zncinn
         hfi(jmu,je,ibout-1,jspc) = zavg*znaout 
     1                          + hfi(jmu,je,ibout-1,jspc)*zncout
           hfi(jmu,je,1,jspc)       = hfi(jmu,je,2,jspc)
  57    continue
  58   continue
  59  continue
c
c..reset rhobes, rhobes, hebems, ajzbs
c
        call beams(4)
c
c       mix 'slow alphas'
c
       if(ntype.ne.1) go to 69
c
c       if ntype .eq. 1 all the alphas are 'slow alphas'
c
        itest=13
       if(itest*(ipd/itest).ne.ipd) go to 69
c
c       get the average density
c
        zsum = awslow(ibinn)*zdvinn + awslow(ibout-1)*zdvout
c
       if(ibinn+1.gt.ibout-2) go to 63
        do 62 jz=ibinn+1,ibout-2
          zsum = zsum + awslow(jz) * (zvoli(jz+1)-zvoli(jz))
  62    continue
  63   continue
c
          zavg = zsum / zsum0
c
c       now flatten the alpha density
c
       if(ibinn+1.gt.ibout-2) go to 66
        do 65 jz=ibinn+1,ibout-2
          awslow(jz) = zavg
  65    continue
  66   continue
          awslow(ibinn)   = zavg*znainn + awslow(ibinn)*zncinn
          awslow(ibout-1) = zavg*znaout + awslow(ibout-1)*zncout
            awslow(1)       = awslow(2)
  69  continue
c
c..Short printout
c  --------------
c
c..check to see if there are any more reconnection regions to do
c
        if ( imixed .gt. 1 ) go to 20
c
  70  continue
c
c  --------------------------------
c  End of sawtooth mixing algorithm
c  --------------------------------
c
c..save old profiles
c
      do 72 jz=1,mzones
        tepre(jz) = useh * tes(2,jz)
        tipre(jz) = useh * tis(2,jz)
        ajzpre(jz) = usej * ajzs(1,jz)
          zqpre (jz) = q(jz)
          znepre(jz) = rhoels(2,jz)
  72  continue
c
c..recompute ajtori(j,2) = < J-toroidal / R > [Amp/m**3] zone center
c
        do 74 jz=1,mzones-1
          ajtori(jz,2) = r0ref *
     &  ( avi(jz+1,2,1)*avi(jz+1,3,1)*avi(jz+1,7,1)*bpoli(jz+1)
     &  - avi(jz,2,1) * avi(jz,3,1) * avi(jz,7,1) * bpoli(jz) )
     &  / ( emu0 * (avi(jz+1,12,1) - avi(jz,12,1) ) )
  74  continue
          ajtori(1,2) = ajtori(2,2)
c
c
c..update conservation-check bookkeeping 
c  these lines from sbrtn sawmx1 need to be implemented (bateman, 4-dec-86)
c
c       do 76 ip=lhyd1,limpn
c        acompi(ip)=(zn2(ip)-zn1(ip)) + acompi(ip)
c  76   continue
c
c        enadji=(ze2-ze1) + enadji
c        bpadji=(zp2-zp1) + bpadji
c
c       now that the fundamental arrays have been changed, the
c       auxiliary arrays must be updated
c
      if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 1)
     &  call prtr1d (6,q,mzones,'q before reconnection')
c
        call getchi(1)
        call getchi(2)
c
      if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 1)
     &  call prtr1d (6,q,mzones,'q after reconnection')
c
c  these lines from sbrtn sawmx1 need to be implemented (bateman, 4-dec-86)
c
c       do 78 jz=lcentr,mzones
c        jzm1=jz-1
c
c        call fitter(mzones,xbouni,presus,xbouni(jz),prejz,ier)
c        call fitter(mzones,xbouni,presus,xbouni(jzm1),prejzm1,ier)
c
c        presus(jz)=zint0(jz)*prejzm1+zint1(jz)*prejz
c  78   continue
c
c
c       find the inversion radius
c
      do 82 jz=1,mzones
        tepost(jz)=useh*tes(2,jz)
        tipost(jz)=useh*tis(2,jz)
        ajzpst(jz)=usej*ajzs(1,jz)
  82  continue
c
        zxs=0.
        do 84 jz=lcentr,mzones
        if((tepost(jz+1).gt.tepre(jz+1)).and.(tepost(jz).le.tepre(jz)))
     &  iji=jz
        if(zxs.ne.0.) go to 84
        if(jz.ge.mzones) go to 84
        if((tepost(jz+1).gt.tepre(jz+1)).and.(tepost(jz).le.tepre(jz)))
     &  zxs=xbouni(jz)+dxzoni(jz)*(tepre(jz)-tepost(jz))/
     &    (tepost(jz+1)-tepost(jz)-(tepre(jz+1)-tepre(jz)))
  84  continue
        sawri=zxs*rmins
      zxmix = sqrt ( zvouti / zvoli(mzones) )
        sawrmx=zxmix*rmins
c
c
        sawate=0.5*(ztemax+ztemin)
        sawane=0.5*(znemax+znemin)
        sawadd=0.5*(zddmax+zddmin)
        sawne=(znemax-znemin)/sawane
        sawte=(ztemax-ztemin)/sawate
        if(sawadd.ne.0.) sawdd=(zddmax-zddmin)/sawadd
        if(zfsmax.ne.0.) sawfus=(zfsmax-zfsmin)*2./(zfsmax+zfsmin)
        if(zalmax.ne.0.) sawalf=(zalmax-zalmin)*2./(zalmax+zalmin)
        sawuth=(zutmax-zutmin)*2./(zutmax+zutmin)
        uthavg=0.5*(zutmax+zutmin)
c
        call sawamp(ztemin,znemin,zddmin,zfsmin,zalmin,zutmin)
        ztemax=ztemin
        znemax=znemin
        zddmax=zddmin
        zfsmax=zfsmin
        zalmax=zalmin
        zutmax=zutmin
c
        estpol=(zepoln-zepolo)*usee
c
        ztau=tai*uist-tbsaws
        call sawavg(2,teavg,tiavg,ajzavg,ohravg,ohptot,ztau)
        if(ohptot.ne.0.) tauein=uthavg/ohptot
        if((tbsaws.ge.0.).and.(sawtau.eq.0)) sawtau=ztau
        if(tbsaws.ge.0.) sawtau=0.5*(sawtau+ztau)
        tbsaws=tai*uist
c
c
c       diagnostic printout
c       -------------------
c
c..short printout
c
      write (nout  ,9600) tai*uist,sawrmx,sawri,sawr1,rmins
c      write (ntychl,9600) tai*uist,sawrmx,sawri,sawr1,rmins
 9600 format (//' sawtooth crash at time ',f12.6,' sec, r(mix)='
     & ,f6.1,' r(inversion)=',f6.1,' r(q=1)=',f6.1,' cm'
     & ,f6.1,' minor radius'/)
c
c..long printouts
c
       if (mod(isawct,lsawth(3)) .eq. 0 .and. lsawth(2) .gt. 0 ) then
c
c..diagnostic printout
c
         if ( lsawth(2) .gt. 8 ) then
c
           write(nout,*)
           write(nout,*) 'Fundamental arrays after sawtooth crash'
           write (nout,950)
 950       format(t5,'chi(1,j)',t20,'chi(2,j)',t35,'chi(3,j)'
     &       ,t50,'chi(4,j)',t65,'chi(lelec,j)',t80,'chi(lion,j)')
           do j=1,mzones
             write (nout,951) chi(1,j), chi(2,j), chi(3,j), chi(4,j)
     &         , chi(lelec,j), chi(lion,j)
 951         format(1p6e15.7)
           enddo
c
         endif
c
c  ix   = index of zone boundary
c  zbpold = B-poloidal before sawtooth crash 
c         = (d pol flux / d rho)/(2*pi*R0REF)
c  bpoli  = B-poloidal after sawtooth crash
c  tepre  = electron temperature before sawtooth crash
c  tepost = electron temperature after  sawtooth crash
c  ztori  = toroidal magnetic flux
c  zhfpre = helical magnetic flux before sawtooh crash
c  zpoli  = poloidal magnetic flux after sawtooth crash
c  zheli  = helical magnetic flux after  sawtooh crash
c
c
       write (nout,9612)
 9612  format (t4,'#',t7,'half',t13,'q-pre',t20,'q-aft'
     & ,t26,'J-pre',t32,'J-aft',t39,'Te-pre',t47,'Te-aft'
     & ,t55,'ne-pre',t63,'ne-aft',t71,'Ti-pre',t79,'Ti-aft'
     & ,t87,'tor-flux',t97,'helc-pre',t107,'helc-aft')
c
        do 9618 j=1,mzones
       write (nout,9614) j,ahalfs(j,1),zqpre(j),q(j)
     & ,ajzpre(j),ajzs(1,j)*usej,tepre(j),tepost(j)
     & ,znepre(j)*1.e-14,rhoels(2,j)*1.e-14
     & ,tipre(j),tipost(j)
     & ,ztori(j),zhfpre(j),zheli(j)
c
 9614  format (t2,i3,f6.1,2f7.3,2f6.3,6f8.3,3f10.4)
c
 9618   continue
c
       endif
c
c..end of long diagnostic printout
c
c..plot profile before crash
c
c  Note: we have to ensure that the time given to sbrtn grafix
c  does not coincide with the next time sbrtn grafix is called.
c
        if (lsawth(5) .gt. 0) then
          nstep = nstep + 1
          tai = tai + dtmini
          tbi = tbi + dtmini
cbate        call grafix(2)
        endif
c
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@recon1  .../baldur/code/bald/dsawmix.f
c  rgb 25-apr-96 replaced svrgn with sort1
c  rgb 24-apr-96 replace zbren with zero = zeroin(zt1,zt2,fpoly1,zerrel)
c  rgb 17:15 12-nov-93 used x(*) rather then x(1) in dummy dimension
c  rgb 14:28 11-nov-93 changed call abort to call abortb
c  rgb 16.01 20-dec-88 changed do 46 zt=ztrout,ztrinn,-zdtor   to
c      do 46 jt=1,idtor-1 ^ zt=ztorx(ixout)-jt*zdtor
c  rgb 12.84 20-nov-87 include any and all remaining values of zhelx(j)
c         at the end of the phelev(jhelev) array after 47  continue
c  rgb 12.67 04-sep-87 Work from the outside in to find helical flux levels.
c       This change is needed to correctly handle the case with
c       one minimum, one maximum, and the helical flux at the magnetic
c       axis in between these.
c**********************************************************************c
c
c     *********************
c     * SUBROUTINE RECON1 *
c     *********************
c
c       Kadomtsev magnetic reconnection algorithm:
c
c     usage:
c           CALL RECON1 (PTORFL,PHELFL,IBOUND,PEPSLN
c    &   ,PHELEV,IHELEV,PTRLEV,ITRLEV,IHELTR,IBINN,IBOUT,IMIXED)
c
c       Input:
c               PTORFL = toroidal flux on BALDUR zone boundaries
c               PHELFL = helical flux on BALDUR zone boundaries
c                = poloidal flux - (m/n)*toroidal flux
c               IBOUND = number of BALDUR zone boundaries (= MZONES)
c               PEPSLN = machine epsilon
c                        = (smallest number such that 1.+PEPSLN > 1.)
c               ITRLEV = maximum number of elements allowed 
c                               in PTRLRV and PHELEV arrays
c                               (use 300, for example)
c               IMIXED > 4 triggers full diagnostic output
c       Output:
c               PHELEV = helical flux in mixing region sorted 
c                               in descending order, J=1,IHELEV
c               IHELEV = number of helical flux levels
c               PTRLEV = toroidal flux within mixing region, J=1,ITRLEV
c                               in ascending order
c               ITRLEV = number of elements in PTRLEV array =< IMXLEV
c               IHELTR(J) = index of helical flux PHELEV(IHELTR(J)) 
c                   corresponding to toroidal flux PTRLEV(J), J=1,ITRLEV
c               IBINN = index of BALDUR zone boundary at inner edge
c                       of mixing region, corresponding to PTRLEV(1)
c               IBOUT = index of BALDUR zone boundary at outer edge
c                         of mixing region, corresponding to PTRLEV(ITRLEV)
c               IMIXED = number of mixing regions found
c                        = -1 if the mixing region covers the entire plasma
c                        = 0 if no mixing regions were found
c                        = 1 if one and only one mixing region found
c                        = 2 if two or more mixing regions found
c
c**********************************************************************c
c
c  ZCH(KBOUND,3) = cubic spline coefficients
c  ZTMP(KBOUND)  = temporary scratch array
c  ZHELX(j)      = min-max values of helical flux
c  ZTORX(j)      = corresponding values of toroidal flux
c  IBNDX(J)      = BALDUR zone boundary just below (ZTORX(J),ZHELX(J))
c  IHMNMX(J)     = +1 for maximum value
c                  -1 for minimum value
c       IHELX = number of elements on this list
c
c**********************************************************************c
c
        subroutine RECON1 (PTORFL,PHELFL,IBOUND,PEPSLN
     &   ,PHELEV,IHELEV,PTRLEV,ITRLEV,IHELTR,IBINN,IBOUT,IMIXED)
c
        parameter (kbound=100)
c
c  kbound = dimension allowed for BALDUR zone boundaries
c
        dimension  zch(kbound,3), ztmp(kbound)
     & , ztorfl(kbound), zhelfl(kbound)
     & , zhelx(kbound), ztorx(kbound), ibndx(kbound)
     & , ihmnmx(kbound)
     & , ptorfl(*), phelfl(*)
     & , phelev(*), ptrlev(*), iheltr(*)
cap 
	integer, parameter         :: kjmix = 900
        integer, dimension(1)      :: iph
        real, dimension(kjmix)     :: zphelev
c
c..The following common block is to pass cubic interpolation coeffs
c  to sbrtn fpoly1 for the zero finder below
c
      common /cpoly1/ zpoly1(5)
c
      external fpoly1
c
c**********************************************************************c
c
c..No mixing regions found yet
c
        iprint = imixed
        imixed = 0
c
c  Make sure there is odd symmetry across the magnetic axis at j=2
c
c  zsign = + 1. if the outer min-max is a maximum helical flux
c                   (this is the normal case)
c          = - 1. if the outer min-max is a minimum helical flux
c
        zsign  = sign (1. , phelfl(2) - phelfl(ibound) )
        do 6 j=2,ibound
          zhelfl(j) = zsign * ( phelfl(j) - phelfl(2) )
          ztorfl(j) = ptorfl(j) - ptorfl(2)
   6  continue
c
        ztorfl(1) = - ztorfl(3)
        zhelfl(1) = - zhelfl(3)
c
        zdhdt1 = (zhelfl(2)-zhelfl(1)) / (ztorfl(2)-ztorfl(1))
c
c  where zdhdt1 = d helical flux / d toroidal flux at point 1
c    choice made ==> zero second derivative in first zone
c
c  Set up cubic spline interpolation of ZHELFL(J) wrt ZTORFL
c  zsign * helical flux in interval ztorfl(j) =< ztor < ztorfl(j+1)
c       = zhelfl(j) + (ztor-ztorfl(j)) * ( zch(j,1)
c                   + (ztor-ztorfl(j)) * ( zch(j,2)
c                   + (ztor-ztorfl(j)) *   zch(j,3) ) )
c
      call spline (ibound,ztorfl,zhelfl,zch(1,1),zch(1,2),zch(1,3))
c
c      call icsccu (ztorfl,zhelfl,ibound,zch,kbound,ier)
c      if ( ier .gt. 128 ) call abortb (6
c     & ,'problem with call icsccu (ztorfl,zhelfl,... in sbrtn recon1')
c
c..temporary printout
c
      if (iprint .gt. 4) then
       call prtr1d (6,ztorfl,ibound,'ztorfl in recon1')
       call prtr1d (6,zhelfl,ibound,'zhelfl in recon1')
       call prtrv  (6,zch,1,kbound,1,1,1,ibound,3,80
     & ,'cubic spline coefs  zch in sbrtn recon1')
      endif
c
c**********************************************************************c
c
c..Working from the inside out, find all locations where
c       d flhel / d fltor = 0
c  Note: d flhel / d fltor
c        = zch(j,1) + 2. * zch(j,2) * (ztor-ztorfl(j)) 
c                   + 3. * zch(j,3) * (ztor-ztorfl(j))**2 
c       for ztorfl(j) =< ztor < ztorfl(j+1)
c  Make a list of
c       ZHELX(j) = min-max values of helical flux
c       ZTORX(j) = corresponding values of toroidal flux
c       IBNDX(J) = BALDUR zone boundary just below (ZTORX(J),ZHELX(J))
c       IHMNMX(J) = +1 for maximum value
c                   -1 for minimum value
c       IHELX = number of elements on this list
c
c  If any ZHELX(j) falls outside the range
c   (helical flux at magnetic axis, helical flux at edge of plasma)
c   include the helical flux at the magnetic axis on this list
c
c  First find the range of helical flux zhelfl(j), j=2,ibound-1
c
        zhlmin = zhelfl(2)
        zhlmax = zhelfl(2)
      do 8 j=3,ibound-1
        zhlmax = max ( zhlmax, zhelfl(j) )
        zhlmin = min ( zhlmin, zhelfl(j) )
   8  continue
c
        ihelx = 0
c
c  Include magnetic axis?
c
        if (zhlmax .gt. max(zhelfl(2),zhelfl(ibound)) .or.
     &    zhlmin .lt. min(zhelfl(2),zhelfl(ibound)) ) then
          ihelx = ihelx + 1
          ztorx(ihelx) = ztorfl(2)
          zhelx(ihelx) = zhelfl(2)
          ibndx(ihelx) = 2
          ihmnmx(ihelx) = -1
          if ( zhelfl(3) .lt. zhelfl(2) .or. zch(2,1) .lt. 0. )
     &       ihmnmx(ihelx) = 1
        endif
c
        do 10 j=2,ibound-2
c
        if (sign(1.,zch(j,1))*sign(1.,zch(j+1,1)) .lt. 0.) then
          ihelx = ihelx + 1
c
          ztorx(ihelx) = ztorfl(j)
c
c..Find the value of toroidal flux ztorx(ihelx) for which 
c  helical flux goes through a minimax by solving the quadratic equation
c  f(x) = ztmp(2) + (x-ztmp(1)) * ztmp(3) + (x-ztmp(1))**2 * ztmp(4)
c  where ztmp(1)=ztorfl(j), ztmp(2)=zch(j,1), ztmp(3)=2.*zch(j,2),
c  and ztmp(4)=3.*zch(j,3)
c
         ztmp(1) = ztorfl(j)
         ztmp(2) = zch(j,1)
         ztmp(3) = zch(j,2) * 2.
         ztmp(4) = zch(j,3) * 3.
c
         call qdsolv (ztmp,ztorfl(j),ztorfl(j+1),0.,ztorx(ihelx))
c
          ztorx(ihelx) = min ( ztorx(ihelx), ztorfl(j+1) )
          ztorx(ihelx) = max ( ztorx(ihelx), ztorfl(j) )
c
c..Now evaluate the cubic spline at the minmax point to determine ZHELX(IHELX)
c
      zhelx(ihelx) = zhelfl(j) + (ztorx(ihelx)-ztorfl(j)) * ( zch(j,1)
     &               + (ztorx(ihelx)-ztorfl(j)) * ( zch(j,2)
     &               + (ztorx(ihelx)-ztorfl(j)) *   zch(j,3) ) )
c
          ibndx(ihelx) = j
          ihmnmx(ihelx) = - nint ( sign ( 1., zch(j+1,1)-zch(j,1) ) )
c
        endif
  10  continue
c
c..if no mixing regions were found
c
        if (ihelx .lt. 1) then
          imixed = 0
          write (6,*) 
     &  'no values of d flhel / d fltor = 0 found in sbrtn recon1'
          go to 90      ! for long printout
        endif
c
c..if the reconnection region covers the entire plasma
c
        if (ihelx .lt. 2) then
          imixed = -1
          write (6,*)
     & 'reconnection covers the entire plasma in sbrtn recon1'
          go to 90      ! long printout from sbrtn recon1
        endif
c
        if ( ihmnmx(ihelx) .lt. 0 ) call abortb (6
     & ,'the outer min-max is not a maximum in sbrtn recon1')
c
c**********************************************************************c
c
c..determine the extent of the outer mixing region
c
c  find the deepest minimum 
c
        zhmin = zhelx(1)
        ixmin = 1
        do 22 j=2,ihelx
          if ( zhelx(j) .lt. zhmin ) then
            zhmin = zhelx(j)
            ixmin = j
            ztmin = ztorx(j)
          endif
  22  continue
c
c  between the deepest minimum and the outer maximum,
c  find the largest maximum value of zhelx(j)
c
        zhmax = zhelx(ixmin)
        do 24 j=ixmin,ihelx
          if ( zhelx(j) .gt. zhmax ) then
            zhmax = zhelx(j)
            ixmax = j
            ztmax = ztorx(j)
          endif
  24  continue
c
c  This establishes the range of helical flux zhmax - zhmin
c  for this reconnection region.
c..Now establish the range of minimax points IXINN to IXOUT
c  There are  ixout - ixinn + 1  helical flux min-max values
c  within this region.
c
      ixinn = 1
      ixout = ihelx
      if ( ixmin .gt. 1 ) then
        do 26 j=ixmin,1,-1
          if ( zhelx(j) .gt. zhmax ) go to 27
          ixinn = j
  26    continue
  27  continue
      endif
c
c..Check to see if there are more than one mixing regions
c
        imixed = 1
        if ( ixout - ixinn + 1 .lt. ihelx) imixed = 2
        if ( ixout - ixinn + 1 .gt. ihelx) call abortb (6
     & ,'ixout - ixinn + 1 .gt. ihelx in sbrtn recon1')
        if ( ixout - ixinn + 1 .lt. 2    ) call abortb (6
     & ,'ixout - ixinn + 1 .lt. 2     in sbrtn recon1')
c
c**********************************************************************c
c
c..Find the range of toroidal flux for this reconnection region:
c
c  Working from the outside in, find the BALDUR zone boundary IBOUT
c  with helical flux just below zhmin
c
        do 32 j=ibound,2,-1
          if ( zhelfl(j) .gt. zhmin) go to 33
          ibout = j
  32  continue
        call abortb (6
     & ,'could not find outer edge of mixing region in sbrtn recon1')
  33  continue
c
c  Working from the inside out, find the BALDUR zone boundary IBINN
c  with helical flux just above zhmax, or default to magnetic axis (j=2)
c
        do 34 j=2,ibout
          ibinn = j
          if ( zhelfl(j) .lt. zhmax) go to 35
  34  continue
        call abortb (6
     & ,'could not find inner edge of mixing region in sbrtn recon1')
  35  continue
c
c  ZTRNGE = range of toroidal flux within the mixing region
c  ZDTOR  = reference toroidal flux difference used to divide up region
c               between the inner and outer minimax points
c
        ztrnge = abs ( ztorfl(ibout) - ztorfl(ibinn) )
        idtor  = abs ( ( ixout - ixinn + 1 ) * ( ibout - ibinn + 1 ) )
        idtor  = min ( idtor, itrlev - ihelx - 4 )
        zdtor  = abs( ztorx(ixout) - ztorx(ixinn) ) / ( idtor )
c
c**********************************************************************c
c
c..Constrct an array of helical flux levels PHELEV(J), J=1,IHELEV
c
c  Start with the outer maximum and work in to the inner minimum
c  Step through the mixing region 
c    using equal intervals of toroidal flux ZDTOR.
c
c  Let INX be a counter of min-max regions passed.
c  Include points only if INX = 0.
c  Start with INX = 0.
c  Step through the mixing region 
c    using equal intervals of toroidal flux ZDTOR.
c  For each min-max value ZHELX,ZTORX encountered:
c    Include PHELEV(IHELEV) = ZHELX(J)
c    IF (INX .NE. 0) reset INX = 0.
c    IF (INX .EQ. 0)  INX = INX + 1  when passing a maximum
c                           INX = INX - 1  when passing a minimum
c  For each distant minimum passed, INX = INX - 1
c  For each distant maximum passed, INX = INX + 1
c  Do not include the ZHELX values of distant min-max passed.
c  Omit exact duplicates.
c  Include the BLADUR zone boundary at the inner edge of the mixing region.
c  IHELEV = the number of PHELEV(j) values found.
c
c  Let INX be a counter of min-max regions passed.
c  Include points only if INX = 0.
c
c
      ihelev = 1
      phelev(ihelev) = zhelx(ixout)
      inx = 0
      jx = ixout - 1
      jb = ibndx(ixout)
c
        do 46 jt = 1,idtor-1
          zt = ztorx(ixout) - jt*zdtor
c
c  For each min-max value ZHELX,ZTORX encountered:
c    Include PHELEV(IHELEV) = ZHELX(J)
c    IF (INX .NE. 0) reset INX = 0.
c    IF (INX .EQ. 0)  INX = INX + 1  when passing a maximum
c                           INX = INX - 1  when passing a minimum
c
        if (zt .le. ztorx(jx) .and. zt+zdtor .gt. ztorx(jx)) then
          ihelev = ihelev + 1
          phelev(ihelev) = zhelx(jx)
c
          if (inx .eq. 0) then
            inx = ihmnmx(jx)
          else
            inx = 0
          endif
c
          jx = max ( jx-1, ixinn )
c
        else
c
c  No min-max value encountered.
c  Determine the helical flux ZHEL at toroidal flux ZT
c
      zhel = seval (ibound,zt,ztorfl,zhelfl,zch(1,1),zch(1,2),zch(1,3) )
c
c      call icsevu (ztorfl,zhelfl,ibound,zch,kbound,zt,zhel,1,ier)
c      if (ier .gt. 32) call abortb (6
c     & ,'problem with call icsevu (ztorfl,zhelfl,.. in sbrtn recon1')
c
c  Is there a value of zhelx(j) between phelev(ihelev) and zhel?
c  For each distant minimum passed, INX = INX - 1
c  For each distant maximum passed, INX = INX + 1
c  Do not include the ZHELX values of distant min-max passed.
c
          do 42 j=ixinn,ixout
            if (j .ne. jx+1 .and. j .ne. jx .and.
     &          zhelx(j) .ge. min(phelev(ihelev),zhel) .and.
     &          zhelx(j) .lt. max(phelev(ihelev),zhel) )
     &          inx = inx + ihmnmx(j)
  42      continue
c
c  If INX = 0, add ZHEL to the helical flux levels.
c
          if (inx .eq. 0) then
            ihelev = ihelev + 1
            phelev(ihelev) = zhel
          endif
c
        endif
c
  46  continue
c
c  End with the inner minimum helical flux within the mixing region
c
  47  continue
      if ( jx .ge. ixinn ) then
        ihelev = ihelev + 1
        phelev(ihelev) = zhelx(jx)
        jx = jx - 1
        go to 47
      endif
c
      if ( ihelev .lt. 4 ) then
        imixed = 0
        go to 90
      endif
c
c..Sort PHELEV(J), J=1,IHELEV  in DESCENDING order
c
cap      call sort1 ( ihelev, phelev, phelev )
c
cbate        call vsrta (phelev,ihelev)
cbate        call svrgn ( ihelev, phelev, phelev )
c
cap        ihalf = ihelev / 2
cap        do 48 j=1,ihalf
cap          jr = ihelev + 1 - j
cap          zhel = phelev(j)
cap          phelev(j) = phelev(jr)
cap          phelev(jr) = zhel
cap  48  continue
cap {
      iph = minloc(phelev(1:ihelev))
      liph = iph(1)
      do j = 1, ihelev
         iph = maxloc(phelev(1:ihelev))
         jiph = iph(1)
         zphelev(j) = phelev(jiph)
         phelev(jiph) = phelev(liph)
      enddo
      phelev(1:ihelev) = zphelev(1:ihelev)
cap }
c
c**********************************************************************c
c
c..Find PTRLEV(JT), JT=1,ITRLEV =< IMXLEV
c  as the toroidal flux values at the points 
c  where PHELEV(JH), JH=JHMIN,JHMAX, which is sorted in descending order,
c  intersects with the curve
c  ZTORFL(JB), ZHELFL(JB), JB=IBINN,IBOUT
c  JT = index of ptrlev(jt)
c  JB = index of ztorfl(jb) and zhelfl(jb)
c  JX = index of ztorx (jx) and zhelx (jx)
c
        jt = 0
        jb = ibinn
        jx = ixinn
c
c  Find where to start in the PHELEV(JH) array 
c
        if ( zhelfl(jb) .gt. phelev(2) ) then
          jdh = + 1
          jh = 1
        elseif ( zhelfl(jb) .lt. phelev(ihelev-1) ) then
          jdh = - 1
          jh = ihelev
        else
          jdh = - sign (1.,zch(ibinn,1))
          do 52 j=3,ihelev
            jh = j
          if ( zhelfl(jb) .lt. phelev(j-1) .and.
     &         zhelfl(jb) .gt. phelev(j+1) ) go to 53
  52      continue
  53      continue
        endif
c
c
c  Now step through the BALDUR grid within the mixing region
c
cbate        write (6,300)
cbate 300  format (3x,'jb jh jx jt  ztor1',7x,'zhel1',7x
cbate     &  ,'ztor2',7x,'zhel2')
c
        do 58 jb=ibinn,ibout-1
c
        ztor1 = ztorfl(jb)
        zhel1 = zhelfl(jb)
        ztor2 = ztorfl(jb+1)
        zhel2 = zhelfl(jb+1)
c
c  Is there a min-max within the next BALDUR zone?
c
        if ( jx .le. ixout .and. ztorfl(jb) .lt. ztorx(jx)
     &      .and. ztorx(jx) .lt. ztorfl(jb+1) ) then
           ztor2 = ztorx(jx)
           zhel2 = zhelx(jx)
        endif
c
  56  continue
c
cbate        write (6,310) jb,jh,jx,jt,ztor1,zhel1,ztor2,zhel2
cbate 310  format (2x,4i3,1p4e13.5)
c
        if ( jh .gt. ihelev .or. jh .lt. 1 ) go to 60
c
c  Is the next level a min-max?
c
        if ( jx .le. ixout .and. phelev(jh) .eq. zhelx(jx) ) then
          jt = jt + 1
          ptrlev(jt) = ztorx(jx)
          iheltr(jt) = jh
          jdh = ihmnmx(jx)
          jh = jh + jdh
          ztor1 = ztorx(jx)
          zhel1 = zhelx(jx)
          jx = jx + 1
          if (ztor1 .eq. ztorfl(jb+1)) go to 58
          ztor2 = ztorfl(jb+1)
          zhel2 = zhelfl(jb+1)
          go to 56
        endif
c
c
c  is the next level on a BALDUR zone boundary?
c
        if ( phelev(jh) .eq. zhelfl(jb) ) then
          jt = jt + 1
          ptrlev(jt) = ztorfl(jb)
          iheltr(jt) = jh
          jh = jh + jdh
          go to 56
        endif
c
c  Are there any helical flux levels PHELEV(JH) between zhel1 and zhel2?
c
            if ( phelev(jh) .gt. min(zhel1,zhel2) .and.
     &     phelev(jh) .lt. max(zhel1,zhel2) ) then
c
                jt = jt + 1
c
c..Solve the polynomial equation in the interval jb to jb+1
c  isig = 4 = number of significant figures needed for the solution
c  zt1 = ztorfl(jb) - ztor on opposite sides of the solution
c  imaxfn = maximum number of evaluations allowed
c
       isig = 4
       zt1 = ztor1 - ztorfl(jb)
       zt2 = ztor2 - ztorfl(jb)
       zpoly1(1) = zch(jb,1)
       zpoly1(2) = zch(jb,2)
       zpoly1(3) = zch(jb,3)
       zpoly1(4) = zhelfl(jb) - phelev(jh)
       imaxfn = 100
c
       zerrel = 1.e-4
c
       zero = zeroin ( zt1, zt2, fpoly1, zerrel )
c
c       call zbren (fpoly1,0.,zerrel,zt1,zt2,imaxfn)
c
c       call zbrent (fpoly1,0.,isig,zt1,zt2,imaxfn,ier)
c
c       if (ier .gt. 128) call abortb (6
c     & ,'problem with call zbrent (fpoly1,... in sbrtn recon1')
c
       ptrlev(jt) = ztorfl(jb) + zero
c
                iheltr(jt) = jh
                jh = jh + jdh
                go to 56
c
            endif
c
  58  continue
c
  60  continue
c    
c..Restore original sign and zero value to toroidal and helical flux
c
      do 62 j=1,ihelev
        phelev(j) = phelfl(2) + zsign * phelev(j)
  62  continue
c
      do 64 j=1,itrlev
        ptrlev(j) = ptorfl(2) + ptrlev(j)
  64  continue
c
        itrlev = jt     
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c  long printout from sbrtn recon1
c
  90  continue
c
c..Final checks
c
      if (ihelev .lt. 4     ) imixed = 0
      if (itrlev .lt. ihelev) imixed = 0
      if (ibout-2.lt. ibinn ) imixed = 0
c
c..temporary printout
c
      if ( iprint .gt. 0 ) then
        write (6,110) imixed,ibinn,ibout,ihelev,itrlev,ihelx
 110  format (/,'  imixed =',i4,'  ibinn,ibout =',2i4
     & ,'  ihelev,itrlev =',2i4,'  ihelx =',i4,'  in sbrtn RECON1')
c
cbate        call prtis (6,ihelev,'ihelev')
cbate        call prtis (6,itrlev,'itrlev')
cbate        call prtis (6,ibinn, 'ibinn')
cbate        call prtis (6,ibout, 'ibout')
cbate        call prtis (6,imixed,'imixed')
      endif
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@fpoly1  /11040/bald91/wbaldn1 DSAWMIX
c
c     *********************
c     * SUBROUTINE FPOLY1 *
c     *********************
c
c..Evaluate a cubic polynomial for the zero finder used in sbrtn recon1
c  using the coeffs passed through common /cpoly1/
c
       function fpoly1 (px)
c
       common /cpoly1/ zc(5)
c
       fpoly1 = zc(4) + px * ( zc(1) + px * ( zc(2) + px * zc(3) ) )
c
cbate      write (6,100) px,fpoly1,(zc(i),i=1,4)
cbate 100  format (' x=',1pe13.5,' f=',e13.5,' c(i)=',4e13.5)
c
       return
       end

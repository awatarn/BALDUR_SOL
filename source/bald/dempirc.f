c  11:00 09-may-94 .../baldur/code/bald/dempirc.f  Bateman, PPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  To obtain this file, type  (use appropriate date for yymmdd)
c cfs get /11040/bald92/byymmdd/wcode.tar
c tar xf wcode.tar  (this creates directory .../code and subdirectories)
c cd code/bald      (this changes to subdirectory .../code/bald)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c----------------------------------------------------------------------c
      subroutine empirc (knemp)
c
c
cl      2.20   semiempirical contributions to the diffusivities
c                         and to a generalized pinch
c
c
c@empirc /11040/bald92/wbaldn1 DEMPIRC
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 20.26 16-feb-92 zratj is not computed, removed from write statement
c  mhr 19.10 15-mar-91 made zln always positive to keep xi-i(clless) pos
c  rgb 18.81 16-dec-90 statement number 10536 changed to 536
c  rgb 12-aug-90 18.53 lemprc(2) controls averge ion mass za(jz,9)
c  rgb 08-aug-90 18.52 vectorize local computations
c  rgb 14-apr-90 implemented interpolation of density monitor
c     between breakpoint times gftime(jt)
c  rgb 18.17 21-mar-90 moved header from file DIO to after end of init
c  rgb 18.03 06-feb-90 central density control with d0nmon(j) and gainv
c  mhr 17.01 23-oct-89 correct xi-collisionless
c  mhr 17.00 04-oct-89 add collisionless chi-i Tang model and density tuning
c               printout 3 chi-e, 3 chi-i semiempirical terms
c       print out improved values of dxsemi and vxsemi for density tuning
c  dps 14.02 20-jun-88  use cemprc(15) to control critical value of eta-i
c       in Tang model.
c  dps 14.01 09-jun-88  in removing za(jz,25), definition of zaq1 was wiped
c       out; replace at top of subroutine.
c  rgb 14.07 03-apr-88  introduced cemprc(17) as breakpoint value of
c       znustr from collisionless to collisional Tang model
c       default cemprc(17) = 0.15 set in sbrtn PRESET file DEFAULT
c  rgb 14.06 23-mar-88  changed za(jz,31) to 
c       za(jz,31)=(1.+(cemprc(19)/max(xbouni(jz,epslon))**cemprc(20))
c  rgb 14.05 06-mar-88
c       if cemprc(19) .gt. epslon .and. cemprc(20) .gt. espslon
c       za(jz,31) = ( d / dx ) ( x / ( 1. + (cemprc(19)/x)**cemprc(20) )
c       this is used to taper a trapezoidal profile
c       down to a parabolic profile near the magnetic axis
c  rgb 14.02 25-feb-88 Set up LEMPRC(1) to control Tang model
c       define zemaxs = maximum integrated net power to electrons
c              zimaxs = maximum integrated net power to ions
c              zgmaxs = maximum integrated net power to plasma
c       set up za(jz,24) = normalized integrated power to electrons
c              za(jz,25) = normalized integrated power to ions
c           these replace previous definitions
c  dps 16-dec-87  Remove Tang model interpolation.
c  dps 05-oct-87  Switch zqcyl in global factor no. 12 to zqval (match
c                      Bateman); add global factor for central electron
c                      density; add radial geometrical factor for TSC
c                      Tang model.
c  rgb 20:50 18-jun-87 initialized nmax1 and nmax2 to 1 rather than 0
c                      to handle uniform transport contributions
c rgb 20:00 23-mar-87  removed gxp(l,n,11)=1. for Tang model cfutz(485)>0.
c  rgb 22:00 22-mar-87 added two lines of output: various constants
c  rgb 19:30 19-mar-87 impl cemprc(16) to choose zqval btwn q and qcyl
c       rgb 09-mar-87  Removed
c      zalprs=(alphai(jz)*ealfai(jz)
c     &        +alphai(jz-1)*ealfai(jz-1))*0.5*uisd*uise
c      zbeprs=(hebems(jz)*rhobis(2,jz)
c     &        +hebems(jz-1)*rhobis(2,jz-1))*0.5
c      zthprs=(presus(jz)+presus(jz-1))*0.5
c      zpres=zthprs+2.*(cequil(2)*zalprs+cequil(3)*zbeprs)/3.
c                       in preparation for merge with Stotler's version
c       dps 06-mar-87  Made final changes to bring Tang model into line
c                      with Martha Redi's 1D version: change znustr cal'n
c                      to zone boundaries, use abs(zeta), and insert an
c                      amax1 function in volume element for zgtots, etc.
c       dps 27-jan-87  Alter alpha calculations for dP/dr's for ballooning
c                      package to match new XSCALE
c       dps 07-jan-87  Rename presus -> thrprs, pprcri -> aprcri, and
c                      use totprs for total pressure
c       dps 12-dec-86  Place array aprcri in commhd and remove common
c                      ballon. Also, eliminate usage of cfutz(impirc)
c                      and znew, zold factors.
c       dps 09-dec-86  Updated profile consistent Tang model to match
c                      version provided by M. Redi
c       rgb 23-sep-86  Numerically computed ideal ballooning mode
c               stability criterion now implemented if lsweq(4) > 0
c               for empirical factor 30
c       rgb 10-sep-86  Changed sign of zbb in ballooning mode criterion
c               Corrected zone boundary number in write statements
c       rgb 9-sep-86 changed za(jz,27 to 29) to implement TANG model
c               in the same way as the 1-D BALDUR by Martha Redi.
c               Also fixed za(jz,18) to agree with 1-D BALDUR
c               Added do loops 182 and 184 to define zgtots(jz) from Redi.
c               Added zgtots(mj) to dimension statement.
c       dps 5-sep-86 improved analytic ballooning mode implementation
c       lpk 30-jul-86 implemented alalytic ballooning modes dfutz*(n,30)
c       rgb 7-mar-86 intro armins, armajs as minor and major radius
c       rgb 1-mar-86 global premultipliers computed
c       rgb 29-jan-86 major rewrite to simplify and shorten sbrtn
c       mhr 9-dec-85 added central te and total ohmic power variables
c       mhr 21-nov-85 added eta-i threshold function and deleted switch
c       mhr 20-sep-85 changed presur to apresr
c       mhr 6-sep-85 added switch for eta=i mode in xi
c       mhr 10-may-85 added variables for prof. cons. model-Tang
c       mhr 23-apr-85 added high beta factor
c       mhr 15-mar-85 added inverse parabolic dependence
c       mhr 19-feb-85 added nu*e dependence
c       mhr 25-jan-85 added ces empirc with additional printout d and v
c       ces 18-nov-84  add comp. phys. comm. (cpc) eq. numbers
c       ces 12-oct-84  6 terms in xe printout
c       ces 10-oct-84  no sums in xe printout; delete xetang
c       ces 9-oct-84 correct z17 (now contains xnuel)
c       ces 18-sep-84 correct z17
c       ces 11-sep-84 *3.5->**3.5 in xetang
c       ces 10-sep-84 reset z15-z18 and then initialize z17=1.0
c       ces 5-sep-84 added ti/te ratio
c       ces 4-sep-84 added ti scale height, radius and time factors
c       ces 4-sep-84 comment out connor-taylor constraint
c       mhr 6-mar-64 added printout of xe for ohmic,aux and q regimes
c       mhr 18-nov-83 added rastar/rrstar normalization to epsilon
c               term and allowed qmhd1,qmhd2,qstar to be namelist or
c               cfutz input, cfutz input overriding unless zero
c       mhr 24-oct-83 added overflow protection to shear contribution
c       mhr 6-oct-83 added te,scale height of te,and qmhd corrections
c       fgps 8-jul-83 introduced averaging across successive timesteps.
c       fgps 1-jul-83 modified to prevent instances of "0.0**(power)".
c       fgps 6-jun-83 subroutine set up.
c
c***********************************************************************
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'commhd.m'
c
      logical
     1   inita1        , inita2
!cap
      real ::
     &   zv(10)        , zd(10)        , ze(10)        , zi(10)
     & , za(mj,40)     , zzpre(10,10)  , zgpre(10,10)
     & , zzv(10,idxion), zzd(10,idxion), zav(10,idxion), zad(10,idxion)
     & , zgtots(mj) , zetots(mj) , zitots(mj)
 
c
c
c-----------------------------------------------------------------------
c
c
      data               inita1,inita2 /.true.,.true./
      data               impirc /364/
c
c     "cfutz(impirc)" is introduced so as to adjust weights in the
c     averages across timesteps:
c  !!!feature removed 12-dec-86 (dps) - first step towards vectorization
c
      data               iclass /2/    , isub /20/
c
c
      if (.not.nlomt2(isub)) go to  5
      if (nstep .lt. 2)
     & call mesage(' *** 2.20 subroutine empirc bypassed ')
      return
    5 continue
c
c-----------------------------------------------------------------------
c
c
cl      common blocks and variables modified:
c
c       vxemps(ix,jz), dxemps(ix,jz), xeemps(jz), xiemps(jz) in /comdif/
c
c-----------------------------------------------------------------------
c
c
c     the index "n" identifies different sets of exponents (i.e.,
c     "dfutzv(n,k)", "dfutzd(n,k)", etc.) assigned to the various
c     physical parameters (up to "kmax1" in number) which make up a
c     term consisting of the product of these each exponentiated by
c     a prescribed power.  such product terms, which are indexed by
c     "n", are subsequently multiplied by prescribed weighting factors
c     (i.e., "vxsemi(n,ix)", "dxsemi(n,ix)", "xesemi(n)", or "xi-
c     semi(n)") and the results combined additively to yield the
c     semiempirical transport coefficients:  "vxemps(ix,jz)", "dx-
c     emps(ix,jz)", "xeemps(jz)", and "xiemps(jz)".  note that di-
c     verse dependencies with respect to ionic species "ix", as mani-
c     fested by the various physical parameters being raised to dif-
c     ferent exponents, can be represented by making adjustments from
c     species to species in the prescribed weighting factors for par-
c     ticular "n".
c
c..input variables:
c
c  lemprc(n) n=1,20 control switches
c
c    lemprc(1) controls the Tang model (see cfutz(485) below)
c
c    lemprc(2) = 1 for za(jz,9) = ahmean(1,jz) = mean hydrogenic mass
c              else (default) za(jz,9) = aimass(1,jz) = mean ion mass
c
c  cemprc(n) n=1,20 coefficients and constants used in sbrtn empirc
c
c       n=1     used to choose transitions in Tang model
c
c       n=14    coefficient used in parabolic function 23
c                       (1.+cemprc(14)*(armins(j,1)/armins(isep,1))**2)
c
c       n=15    critical eta-i in Tang model; default = 1.5, the value used
c               in calibration by M. Redi.
c
c       n=16    used to choose zqval to represent the q value in the
c               Tang model:
c               < 0.5 -> zqval=q(mzones)
c               > 0.5 -> zqval=zqcyl
c
c       n=17    breakpoint from collisional to collisionless Tang model
c               (default cempirc(17) = 0.15 set in sbrtn PRESET file DEFAULT)
c
c       n=19    normalization radius for modification of ftrap(j1,j)
c       n=20    exponent for modification of ftrap(j1,j) 
c               in file DSOLVER sbrtn GETCHI ftrap(j1,j) is modified by
c               ftrap(j1,j) = 1. + (ftrap(j1,j)-1.)
c                             / ( 1. + ( cemprc(19) / x )**cemprc(20) )
c               where x = xbouni(j) or xzoni(j)
c
c  cgbeta(n), n=1,20, all defaulted to 1.0
c             relative contribution to total beta from
c     n=1          thermal electrons
c     n=2          thermal ions
c     n=3          beam ions
c     n=4          alpha particles
c
c  cgpowr(n), n=1,20, all defaluted to 1.0
c             relative contribution to total input power from
c     n=1          ohmic heating to electrons
c     n=2          ohmic heating to hydrogenic ions (not present)
c     n=3          alpha heating to electrons
c     n=4          alpha heating to hydrogenic ions
c     n=5          auxilliary heating to electrons (not incl sources below)
c     n=6          auxilliary heating to hydrogenic ions
c     n=7          beam heating to electrons
c     n=8          beam heating to hydrogenic ions
c     n=9          ecrh heating to electrons
c     n=10         ecrh heating to hydrogenic ions
c     n=11         icrh heating to electrons
c     n=12         icrh heating to hydrogentic ions
c
c
c     index l in the following arrays refers to type of transport
c       l=1     electron thermal conductivity
c       l=2     ion thermal conductivity
c       l=3     inward particle pinch velocity
c       l=4     hydrogenic particle diffusion coefficient
c       l=5     electron thermal pinch velocity
c       l=6     ion thermal pinch velocity
c
c   index n enumerates the addends contributing to transport coeffs
c
c
c  gcoefs(l,n) = coefficient of global premultiplier [standard units]
c                all defaulted to 1.0
c
c  gxp(l,n,i)  exponentials used in global premultipiers in sbrtn empirc
c       l=1,6  n=1,6  i=1,30  all defaulted to zero in sbrtn preset
c
c       index i enumerates the factors in each addend:
c  i=1  zeff averaged over plasma up to separatrix   dimensionless
c  i=2  mean ion mass averaged over plasma              amu
c  i=3  elongation of flux surface just inside separatrix   dimensionless
c  i=4  toroidal magnetic field (btort(n) at rtort(n)) in tesla
c  i=5  minor radius up to separatrix in cm
c  i=6  major radius of flux surface just inside separatrix
c  i=7  line average density in units of 1.e13 cm**(-3)
c  i=8  total beta as computed above
c  i=9  total plasma current in megampers
c         up to the flux surface just inside the separatrix
c  i=10  total input power as computed above in megawatts
c  i=11  2.*zk**2 / ( 1. + zk**2) where zk = elong(isep,1) dimensionless
c         note this is approximately 1. / < (del rho)**2 >
c  i=12  zqval as determined by cemprc(16).
c           zqcyl = 5.*( ahalfs [cm] )**2*btoroidal [kg]*( 1 + zk**2 )
c                   / ( 2. * rmids [cm] * plasma current [kA] )
c  i=13  atomic mass of the first ion species
c  i=14  atomic mass of the second ion species  amu
c  i=16 to 21: factors for Tang model
c  i=16  power into electrons / electron density
c  i=17  power into ions / ion density
c  i=18  total power / electron density
c  i=19  "trapped electron" mode scaling factor computed with only power
c        into electrons
c  i=20  "trapped electron" mode scaling factor computed with total power
c  i=21  "eta-i" mode scaling factor
c  i=22  central electron density (1.e14 cm**-3)
c
c  addends zzpre(l,n) of the global premultipliers are computed
c       in the following way:
c       zzpre(l,n) = product over i ( factor(i) ** gxp(l,n,i) )
c
c  gcomb(l,m,n)  exponentials used to combine global premultipliers
c        in sbrtn empirc l=1,6  m=1,6  n=1,6
c        defaulted to kronicka delta(m,n) in sbrtn preset
c
c  the premultiplier addends are combined in the following way:
c  zgpre(l,n) = [sum over m (zzpre(l,m)**gcomb(l,m,n)]**(1./gcomb(l,n,n))
c
c..Header for printout
c
      if ( knemp .eq. 3 ) then
c
        write(nprint,10500)
        write(nprint,10502) memprt,memprt
c
10500   format(10x,36hsemiempirical transport coefficients)
10502   format(3x,4hzone,4x,8h radius ,4x,8h xe(1)  ,4x,8h xe(2)  ,
     1  4x,8h xe(3)  ,4x,8h xi(1)  ,4x,8h xi(2)  ,4x,8h xi(3)  ,
     2  4x,4h  v(,i1,3h)  ,4x,4h  d(,i1,3h)  )
c
      endif
c
c..Tang model:
c
c  selection of exponents for Tang model according to value of cfutz(485).
c  lists of exponents; note that input data may be overwritten if cfutz(485)
c  is not 0. interpolation between regimes, etc. is done at end of empirc.
c  note that rastar is set here also if cfutz(485) > 0.
c
c       cfutz(485) = 1 :  tem chi-e and itg chi-i
c                    2 :  tem chi-e and itg chi-i(hardwired on)
c                   3 :  tem chi-e and itg chi-i driven by total heating
c                   4 :  tem chi-e(collisional always) and itg chi-i
c                   5 :  tem chi-e, tem chi-i, itg chi-i
c
c
      if ( cfutz(485) .gt. 0.  .and.  lemprc(1) .eq. 0 )
     &     lemprc(1) = cfutz(485)
c
      if ( lemprc(1) .gt. 0 ) rastar = armins(mzones,1)
c
      if  ((lemprc(1) .gt. 0) .and. (lemprc(1) .ne. 3))  then
c
      dfutze(2,5)=-2.
      dfutze(2,32)=1.
      dfutze(2,35)=-1.   ! collisionless chi-e
      gxp(1,2,16)=.6
      gxp(1,2,21)=1.
c
      dfutze(1,5)=-2.
      dfutze(1,32)=1.
      dfutze(1,35)=-1.   ! collisional chi-e
      gxp(1,1,19)=1.
c
      end if
      if ((lemprc(1) .gt. 0) .and. ( lemprc(1) .ne. 3 )) then
c
      dfutzi(1,5)=-2.
      dfutzi(1,33)=1.
      dfutzi(1,34)=-1.   ! eta-i chi-i
      gxp(2,1,17)=.6
      gxp(2,1,21)=1.
c
      end if
c
c  for collisionless ctem chi-i model inclusion
c
        if ( lemprc(1) .eq. 5 ) then
        dfutzi(3,5)=-2
        dfutzi(3,32)=1
        dfutzi(3,35)=-1
        gxp(2,3,16)=.6
        gxp(2,3,21)=1.
c
        end if
c
      if ( lemprc(1) .eq. 3 ) then
c
      dfutze(3,5)=-2.
      dfutze(3,29)=1.
      dfutze(3,35)=-1.   ! chi-e with total heating
      gxp(1,3,20)=1.
c
      dfutzi(2,5)=-2.
      dfutzi(2,29)=1.
      dfutzi(2,35)=-1.   ! eta-i chi-i with total heating
      gxp(2,2,18)=.6
      gxp(2,2,21)=1.
c
      end if
c
cl             initializations.
c
        itzone=cfutz(120)-2
        if(abs(cfutz(377)).ge.epslon)qstar=cfutz(377)
        if(abs(cfutz(378)).ge.epslon)qmhd1=cfutz(378)
        if(abs(cfutz(379)).ge.epslon)qmhd2=cfutz(379)
c
c
      if(.not.inita1) go to 170
 
      inita1=.false.
c
c     set do-loop limits.
c
c *** note:  !*  identifies lines to be changed
c                when adding new local empirical factors
c
      jzmin=lcentr+1
      nlmax=6      ! number of types of transport coefs
      nxmax=6
      nmax1=1
      nmax2=1
      kmin=3
      kmax = 40    !* second dimension of dfutz's
      kmax1 = 36   !* number of dfutz's actually in use
c
      do 10 n=1,nxmax
      do 10 j=1,mxzone
      za(j,n) = 0.0
  10  continue
c
c
      k1=1
c
c       comment out conor-taylor constraint to avoid unsuspecting use
c       (c. f. cpc paper, ref. no. [36])
c       do 12 n=1,nxmax
c       zpv=-dfutzv(n,7)
c       zpd=-dfutzd(n,7)
c       zpe=-dfutze(n,7)
c       zpi=-dfutzi(n,7)
cc
c       zqv=-dfutzv(n,14)
c       zqd=-dfutzd(n,14)
c       zqe=-dfutze(n,14)
c       zqi=-dfutzi(n,14)
cc
c       zrv=1-dfutzv(n,4)
c       zrd=1-dfutzd(n,4)
c       zre=1-dfutze(n,4)
c       zri=1-dfutzi(n,4)
cc
c12     continue
c
c
      do 110 k=1,kmax
      do 110 n=1,nxmax
      if(abs(dfutzv(n,k)).gt.epslon) nmax1=max0(nmax1,n)
      if(abs(dfutzd(n,k)).gt.epslon) nmax1=max0(nmax1,n)
      if(abs(dfutze(n,k)).gt.epslon) nmax2=max0(nmax2,n)
      if(abs(dfutzi(n,k)).gt.epslon) nmax2=max0(nmax2,n)
c
c       minor radius scaling automatically obeys connor-
c       taylor constraint r=2p+q/2+5r/4 if dfutz?(n,5).le.epslon
c
c       if(abs(dfutzv(n,5)).le.epslon)
c    1          dfutzv(n,5)=2.*zpv+zqv/2.+1.25*zrv
c       if(abs(dfutzd(n,5)).le.epslon)
c    1          dfutzd(n,5)=2.*zpd+zqd/2.+1.25*zrd
c       if(abs(dfutze(n,5)).le.epslon)
c    1          dfutze(n,5)=2.*zpe+zqe/2.+1.25*zre
c       if(abs(dfutzi(n,5)).le.epslon)
c    1          dfutzi(n,5)=2.*zpi+zqi/2.+1.25*zri
  110 continue
c
c     where default values for the elements of arrays "dfutzv(n,k)",
c     "dfutzd(n,k)", "dfutze(n,k)", and "dfutzi(n,k)" are zero as set
c     during the initial zeroing of all common-block variables.
c
c
c     default values for the weighting factors are given below:
c
      do 122 n=1,nxmax
      vx0=vxsemi(n,1)
      dx0=dxsemi(n,1)
      do 122 ix=2,mxions
      absvx=abs(vxsemi(n,ix))
      absdx=abs(dxsemi(n,ix))
      if(absvx.le.epslon) vxsemi(n,ix)=vx0
      if(absdx.le.epslon) dxsemi(n,ix)=dx0
  122 continue
c
c
c     default values for normalization factors.
c
      if(bzstar.le.epslon) bzstar=12.0
c     where "bzstar" is expressed in "kg".
      bzstrs=bzstar*uesb
c
      if(elecd0.le.epslon) elecd0=4.0e13
c     where "elecd0" is given in "cm**-3".
c
      if(apresr.le.epslon) apresr=8.0e13
c     which is expressed in "kev/cm**3".
      apresr=apresr*uesh
c
        if(electe.le.epslon) electe=1.0
c       which is expressed in kev
        elects=uesh*electe
c
      if(rastar.le.epslon) rastar=(2./3.)*40.0
c     which is in "cm".
c
        if(qstar.le.epslon) qstar=2.2
c
        if(rrstar.le.epslon) rrstar=140.
c       which is in cm.
c
  170 continue   ! end of initialization
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..central density control
c
      if ( d0nmon(1) .ge. 1.0 .and. gainv .gt. epslon ) then
c
c  find the appropriate breakpoint time
c  extend the range of reasonable values for d0nmon(jt) if necessary
c
        itime = 1
        do 13 jt=2,20
          if ( d0nmon(jt) .lt. 1.0 ) d0nmon(jt) = d0nmon(jt-1)
          itime = jt
          if ( tai .lt. gftime(jt)*usit ) go to 14
  13    continue
  14    continue
c
c  interpolate
c
        call timint (tai*uist, zd0mon, gftime, itime, d0nmon, 1, 1)
c
c  compute pinch multiplier
c 
        zgain = min ( 1.0, max ( -1.0, dti * uist * gainv *
     &    ( ( zd0mon / rhoels(2,lcentr) ) - 1.0 ) ) )
c
        zmult = ( 2.0 + zgain ) / ( 2.0 - zgain )
c
c  now increase or decrease vxsemi(jn,jx)
c
        do 16 jx=1,mxions
          do 16 jn=1,nxmax
            vxsemi(jn,jx) = vxsemi(jn,jx) * zmult
  16    continue
c
      endif
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..global premultiplier to the semi-empirical transport
c
c..compute beta, total input power, zeff, and mean ion mass
c
      isep = mzones  ! zone bndry at outer edge of plasma
      if (nadump(1) .gt. lcentr) isep = nadump(1)  ! separatrix bndry
      isepm1 = isep - 1
c
      zbeta = max (epslon, cgbeta(1) * betate(isep)
     &          + cgbeta(2) * betati(isep)
     &          + cgbeta(3) * betatb(isep)
     &          + cgbeta(4) * betata(isep) )
c
      zpowr = max (epslon, cgpowr(1)*geohms(isep)
     &   + cgpowr(3)*gealfs(isep) + cgpowr(4)*gialfs(isep)
     &   + cgpowr(5)*geauxs(isep) + cgpowr(6)*giauxs(isep)
     &   + cgpowr(7)*gebems(isep) + cgpowr(8)*gibems(isep)
     &   + cgpowr(9)*geecrs(isep) + cgpowr(10)*giecrs(isep)
     &   + cgpowr(11)*geicrs(isep) + cgpowr(12)*giicrs(isep) )
c
      zpowr = zpowr * usep * 1.e-6     ! convert to MW
c
      znz2 = 0.    ! vol average of zeff(j)*electron density
      zne  = 0.    ! vol average of electron density
      zmass = 0.   ! vol average of ion mass times ion density
      zni  = 0.    ! vol average of ion density
      zee  = 0.
c
      do 26 j=lcentr,isepm1
      znz2 = znz2 + xzeff(2,j)*rhoels(2,j)*2.*dx2i(j)
      zne  = zne  +            rhoels(2,j)*2.*dx2i(j)
      zee  = zee  + chi(lelec,j)*2.*dx2i(j)
c
      do 24 jx=1,mxions
      zmass = zmass + aspec(jx)*chi(jx,j)*2.*dx2i(j)
      zni   = zni   +           chi(jx,j)*2.*dx2i(j)
  24  continue
c
  26  continue
c
      zeff = znz2 / zne      ! mean ion charge of plasma
      zamu = zmass / zni     ! mean ion mass throughout plasma
      ztebar = zee * uisd *uieh / zne
c
c  calculate various global quantities for use with Tang model
c  (including qcyl) in both premultipliers and local factors,
c  from M. Redi (dps 09-dec-86)
c
c  determine zone boundaries where q=1, 1.5 and 2
      jq1=lcentr
      jq2=lcentr
      jqmid=lcentr
      do 27 jqa=lcentr,mzones
      if (q(jqa).lt.1.) jq1=jqa
      if (q(jqa).lt.1.5) jqmid=jqa
      if (q(jqa).lt.2.) jq2=jqa
   27 continue
        if (q(jq1).gt.1.) jqmid=(jq2-jq1)/2
        if (jqmid.eq.0) jqmid=jq2
        if(q(jq1).gt.2.) jqmid=jq2
c
c  calculate eta-i between q=1 and 2
      zeta=((tis(2,jq1)-tis(2,jq2))/tis(2,jqmid))*rhoins(2,jqmid)
     &     /(rhoins(2,jq1)-rhoins(2,jq2)+epslon)
        if(q(2).gt.2) then
        zeta = (tis(2,2)-tis(2,4))/tis(2,3) * rhoins(2,3)
     &  /(rhoins(2,2)-rhoins(2,4) + epslon)
        endif
c
      zeta=abs(zeta)
c
c  calculate nu-star at q=1.5
      znustr=xzeff(1,jqmid)*xnuel(1,jqmid)
c
c  zqval = zqcyl or q(mzones), determined by cemprc(16)
c  zqcyl = 5. * ( ahalfs [cm] )**2 * btoroidal [kg] * ( 1 + zk**2 )
c          / ( 2. * rmids [cm] * plasma current [kA] )
c  zk = elongation, used here and in other premultipliers
c
      zk=elong(isep,1)
      zqcyl = ( 5. * ahalfs(isep,1)**2 * bzs*useb * ( 1. + zk**2 ) )
     &         / ( 2. * rmids(isep,1) * curnts(isep) * usei )
c
c..set zqval according to value of cemprc(16)
c
      zqval = q(mzones)
      if ( cemprc(16) .gt. 0.5 ) zqval = zqcyl
c
c..factor used in Tang model:
c
      zaq1 = zqval + .5
c
c  compute parabolic density exponent using vol avg dens found above
c  note: needs to be generalized for non-circular plasmas
      zalphn=(rhoels(2,lcentr)-rhoels(2,mzones))
     &       /(zne-rhoels(2,mzones)+epslon)-1.
      if (zalphn.le.epslon) zalphn=0.
c
c  compute energy variables (in watts) zgtots, zetots and zitots,
c  cf. M. Redi (dps 09-dec-86)
      call resetr(zgtots,mxzone,0.0)
      call resetr(zetots,mxzone,0.0)
      call resetr(zitots,mxzone,0.0)
c
c       jinvld = zone for warning if model is invalid
c
        jinvld = 0
      do 30 j=2,mzones
      zdvols=2.*vols(mzones,1)*usep*max(dx2i(j-1),0.0)
c
      zgtots(j)=zgtots(j-1)+zdvols*(weohms(j-1)-webrs(j-1)
     &          +webems(j-1)+wibems(j-1)+weauxs(j-1)+wiauxs(j-1)
     &          +wiecrh(j-1)+weecrh(j-1)+wiicrf(j-1)+weicrf(j-1)
     &          +wealfs(j-1)+wialfs(j-1)-wesrs(j-1))
c
      zetots(j)=zetots(j-1)+zdvols*(weohms(j-1)-webrs(j-1)
     &          +webems(j-1)+weauxs(j-1)+weecrh(j-1)+weicrf(j-1)
     &          -cnueqs(j-1)*(tes(2,j-1)-tis(2,j-1))+wealfs(j-1)
     &          -wesrs(j-1))
c
      zitots(j)=zitots(j-1)+zdvols*(wibems(j-1)+wiauxs(j-1)
     &          +wiecrh(j-1)+wiicrf(j-1)
     &          +cnueqs(j-1)*(tes(2,j-1)-tis(2,j-1))+wialfs(j-1))
c
      if (mimp.le.0) go to 30
      do 28 ji=1,mimp
      zgtots(j)=zgtots(j)-zdvols*weirs(ji,j-1)
      zetots(j)=zetots(j)-zdvols*weirs(ji,j-1)
   28 continue
c
        if(zetots(j).le.0) jinvld = j
c
   30 continue
c
c..find maximum integrated power for normalizations below
c
        zemaxs = epslon
        zimaxs = epslon
        zgmaxs = epslon
      do 32 j=2,mzones
        zemaxs = max ( zemaxs, zetots(j) )
        zimaxs = max ( zimaxs, zitots(j) )
        zgmaxs = max ( zgmaxs, zgtots(j) )
  32  continue
c
c..compute global premultiplier
c
      do 34 jn=1,nxmax       ! addend in transport coef
      do 34 jl=1,nlmax       ! type of transporc coef
c
      zpm = gcoefs(jl,jn)     ! global premultiplier
      zk = elong(isep,1)
c
c  zeff averaged over plasma
c
      if(gxp(jl,jn,1).ne.0.) zpm=zpm*(zeff**gxp(jl,jn,1))
c
c  atomic mass averaged over plasma in units of amu
c
      if(gxp(jl,jn,2).ne.0.) zpm=zpm*(zamu**gxp(jl,jn,2))
c
c  elongation of outer flux surface
c
      if(gxp(jl,jn,3).ne.0.) zpm=zpm*(elong(isep,1)**gxp(jl,jn,3))
c
c  vacuum toroidal magnetic field (btort(n) at rtort(n)) in tesla
c
       if (gxp(jl,jn,4) .ne. 0.)
     & zpm = zpm * ( (abs(bzs)*1.e-4)**gxp(jl,jn,4) )
c
c  minor radius of flux surface at separatrix or plasma edge  in cm
c
       if (gxp(jl,jn,5).ne.0.) zpm=zpm*(armins(isep,1)**gxp(jl,jn,5))
c
c  major radius of flux surface at separatrix or plasma edge  in cm
c
       if (gxp(jl,jn,6).ne.0.) zpm=zpm*(armajs(isep,1)**gxp(jl,jn,6))
c
c  line average density in units of 1.e13 cm**(-3)
c
       if (gxp(jl,jn,7) .ne. 0.)
     &   zpm = zpm * ( (enlaes(isep)*1.e-13)**gxp(jl,jn,7) )
c
c  total beta as computed above
c
      if (gxp(jl,jn,8).ne.0.) zpm=zpm*(zbeta**gxp(jl,jn,8))
c
c  total plasma current in megampers
c   up to the flux surface at separatrix or plasma edge
c
      if (gxp(jl,jn,9) .ne. 0.)
     &   zpm = zpm * ( abs(curnts(isep)*usei*1.e-3)**gxp(jl,jn,9) )
c
c  total input power as computed above in megawatts
c
      if (gxp(jl,jn,10) .ne. 0.)  zpm = zpm * (zpowr**gxp(jl,jn,10))
c
c  2. * zk**2 / ( 1. + zk**2) where zk = elong(isep,1) dimensionless
c     note this is approximately 1. / < (del rho)**2 >
c
      if (gxp(jl,jn,11) .ne. 0.) then
      zpm = zpm * ( (2. * zk**2 / (1. + zk**2))**gxp(jl,jn,11) )
      endif
c
c  zqcyl = 5. * ( ahalfs [cm] )**2 * btoroidal [kg] * ( 1 + zk**2 )
c          / ( 2. * rmids [cm] * plasma current [ka] )
c
      if (gxp(jl,jn,12) .ne. 0.) then
      zpm = zpm * abs(zqval)**gxp(jl,jn,12)
      endif
c
c      mass scaling, first species
c
        if( gxp(jl,jn,13) .ne. 0. ) then
        zpm = zpm * (abs(aspec(1)))**gxp(jl,jn,13)
        endif
c
c      second species if any
c
        if( mhyd .eq. 1 .and. gxp(jl,jn,14) .ne. 0. ) then
        zpm = zpm * (abs(aspec(1)))**gxp(jl,jn,14)
        else
        if( mhyd .eq. 2 .and. gxp(jl,jn,14) .ne. 0. ) then
        zpm = zpm * (abs(aspec(2)))**gxp(jl,jn,14)
        endif
        endif
c
      if( gxp(jl,jn,15) .ne. 0. ) then
      zpm = zpm * ( ztebar**gxp(jl,jn,15) )
      endif
c
c  factors i=16 to 21 for Tang model, cf. M. Redi (dps 09-dec-86)
c
c  power into electrons (mw) / electron density (1.e14 cm**-3)
      if (gxp(jl,jn,16).ne.0.) zpm=zpm*
     &  (((max(zetots(mzones),epslon)/1.e6)
     &  /((rhoels(1,lcentr)+epslon)/1.e14))**gxp(jl,jn,16))
c
c  power into ions (mw) / ion density (1.e14 cm**-3)
      if (gxp(jl,jn,17).ne.0.) zpm=zpm*
     &  (((max(zitots(mzones),epslon)/1.e6)
     &  /((rhoins(1,lcentr)+epslon)/1.e14))**gxp(jl,jn,17))
c
c  total power (mw) / electron density (1.e14 cm**-3)
      if (gxp(jl,jn,18).ne.0.) zpm=zpm*
     &  (((max(zgtots(mzones),epslon)/1.e6)
     &  /((rhoels(1,lcentr)+epslon)/1.e14))**gxp(jl,jn,18))
c
c  "trapped electron" mode scaling factor computed with only
c  power into electrons
      if (gxp(jl,jn,19).ne.0.) zpm=zpm*((6.e3*(1.+.25*zalphn)
     &  *(max(zetots(mzones),epslon)/1.e6)**.8/((armajs(isep,1)
     &  /100.)**1.2*zqval**.9*(bzs/1.e4)**.4*xzeff(1,lcentr)**.2
     &  *(armins(isep,1)/100.)**.1*(rhoels(1,lcentr)/1.e14)))
     &  **gxp(jl,jn,19))
c
c  "trapped electron" mode scaling factor computed with total power
      if (gxp(jl,jn,20).ne.0.) zpm=zpm*((6.e3*(1.+.25*zalphn)
     &  *(max(zgtots(mzones),epslon)/1.e6)**.8/((armajs(isep,1)
     &  /100.)**1.2*zqval**.9*(bzs/1.e4)**.4*xzeff(1,lcentr)**.2
     &  *(armins(isep,1)/100.)**.1*(rhoels(1,lcentr)/1.e14)))
     &  **gxp(jl,jn,20))
c
c  "eta-i" mode scaling factor
      if (gxp(jl,jn,21).ne.0.) zpm=zpm*((1.e4/((armajs(isep,1)
     &  /100.)**.8*(bzs/1.e4)**.8*(armins(isep,1)/100.)**.2
     &  *zqval**.8))**gxp(jl,jn,21))
c
c
c  central electron density (1.e14 cm**-3)
c
      if (gxp(jl,jn,22).ne.0.) zpm=zpm*((rhoels(1,lcentr)/1.e14)
     &                                **gxp(jl,jn,22))
      zzpre(jl,jn) = zpm     ! assembled global premultiplier
c
  34  continue
c
c..combine global coefficients using the gcomb array
c  zgpre(l,n) is the final premultiplier
c
      do 38 jl=1,nlmax
      do 38 jn=1,nxmax
c
      if (abs(gcomb(jl,jn,jn)) .lt. epslon) then
      zgpre(jl,jn) = zzpre(jl,jn)
c
      else
        zgpre(jl,jn) = 0.
        do 36 jm=1,nxmax
        if (abs(gcomb(jl,jm,jn)) .gt. epslon)
     &   zgpre(jl,jn) = zgpre(jl,jn) + zzpre(jl,jm)**gcomb(jl,jm,jn)
  36  continue
         zgpre(jl,jn) = zgpre(jl,jn)**(1./gcomb(jl,jn,jn))
c
      endif
c
  38  continue
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      calculations of the semiempirical transport coefficients.
c
c
c     physical parameters independent of radius.
c
c
c     atomic mass (amu).
c
      do 40 ix=1,mxions
      do 40 n=1,nxmax
        zav(n,ix)=1.0
        zad(n,ix)=1.0
        if ( aspec(ix) .gt. epslon ) then
          zav(n,ix)=aspec(ix)**dfutzv(n,2)
          zad(n,ix)=aspec(ix)**dfutzd(n,2)
        endif
  40   continue
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c     definition of some basic physical factors.
c
cpc     @2.1e
c
c
      do 57 jz=jzmin,mzones
c
c   elongation of each flux surface
c
      za(jz,3) = elong(jz,1)
c
c   normalized toroidal magnetic field
c
      za(jz,4) = bzs / bzstrs
c
c   normalized radius.
c
      za(jz,5) = max ( armins(jz,1) / rastar , epslon )
c
c   inverse aspect ratio.
c
      za(jz,6) = max((armins(jz,1)/armajs(jz,1)),epslon)
     &  * rrstar / rastar
c
c   normalized electron density.
c
      za(jz,7)=rhoels(1,jz)/elecd0
c
c   normalized thermal pressure.
c
      za(jz,8)=(thrprs(jz)+thrprs(jz-1))/(2.0*apresr)
c
  57  continue
c
c   mean ion mass aimass(1,jz) or hydrogenic mass ahmean(1,jz) (amu).
c
      if ( lemprc(2) .eq. 1 ) then
c
        do 58 jz=jzmin,mzones
          za(jz,9)=max(ahmean(1,jz),epslon)
  58    continue
c
      else
c
        do 59 jz=jzmin,mzones
          za(jz,9)=max(aimass(1,jz),epslon)
  59    continue
c
      endif
c
      do 75 jz=jzmin,mzones
c
c   z-effective.
c
      za(jz,10)=max(xzeff(1,jz),epslon)
c
c   logarithmic gradient of electron-density.
c
      za(jz,11) =
     &  max( abs((armins(jz,1))/(slnes(jz)+epslon)), epslon)
c
c   logarithmic gradient of thermal pressure.
c
      za(jz,12) = max( abs(slprs(jz)), epslon)
c
c   "[+d(ln q)/d(ln r)]**dfutz---".
cpc     @2.1h
c
        za(jz,13) = max( epslon, abs( cfutz(296) +
     &    shear(jz)+qmhd1/(qmhd2+epslon+abs(q(jz)-q(itzone))) ) )
c
c   normalized electron temperature
c
        za(jz,14) = max( abs(tes(1,jz)/elects), epslon)
c
c   logarithmic gradient of electron temp
c
        za(jz,15) = max( epslon,
     &    cfutz(293) + abs(armins(jz,1)/(sltes(jz)+epslon)) )
c
c   normalized safety factor q
c
        za(jz,16) = max( abs(q(jz)/qstar), epslon)
c
c   logarithmic gradient of ion temperature
c
        za(jz,17) = max(
     &    abs(armins(jz,1)/(sltis(jz)+epslon)), epslon)
c
  75  continue
c
      do 77 jz=jzmin,mzones
c
c   radial factor (primarily for omitting scrapeoff)
c   no time smoothing needed for za(jz,18)
c
        za(jz,18)=1.0
        if(xbouni(jz).gt.cfutz(291)) za(jz,18)=1.0+(cfutz(290)-1.0)*
     &   (xbouni(jz)-cfutz(291))/(abs(cfutz(292)-cfutz(291))+epslon)
        if(xbouni(jz).gt.cfutz(292)) za(jz,18)=cfutz(290)
        za(jz,18)=max(za(jz,18),epslon)
c
  77  continue
c
      do 79 jz=jzmin,mzones
c
c   collisionality-dependent dlnte/dlnr exponent
c
        za(jz,19) = max( epslon,
     &    (max(abs(cfutz(294))+abs((armins(jz,1))/(sltes(jz)+epslon))
     &      ,epslon) )**( xnuel(1,jz) / (cfutz(295) + xnuel(1,jz)) ) )
c
  79  continue
c
c
      do 91 jz=jzmin,mzones
c
c   ion/electron temperature ratio
c
        za(jz,20) = max( cfutz(297)+(tis(1,jz)/tes(1,jz)), epslon)
c
c   function of nu star e
c
        za(jz,21) = max( xnuel(1,jz)/(cfutz(300)+xnuel(1,jz)), epslon)
c
  91  continue
c
      do 93 jz=jzmin,mzones
c
c   high beta cutoff factor
c
          zbetac=max(armins(jz,1)/(armajs(jz,1)*q(jz)*q(jz)),epslon)
          zbeta=8.0*fcpi*thrprs(jz)/(bzs*bzs)
          zratio=max(cfutz(298)*zbeta*za(jz,12)/zbetac,epslon)
          zfact=max(1.0+zratio**cfutz(299),epslon)
        za(jz,22) = max( zfact, epslon)
c
  93  continue
c
      do 99 jz=jzmin,mzones
c
c   inverse parabolic function
c
        za(jz,23) = max( epslon,
     &    1.+cemprc(14)*(armins(jz,1)/armins(isep,1))**2 )
c
c   normalized power to the electrons
c
        za(jz,24) = max ( zetots(jz)/zemaxs, epslon )
c
c   normalized power to the ions
c
        za(jz,25) = max ( zitots(jz)/zimaxs, epslon )
c
c   threshold function
c
        za(jz,26) = max( epslon,
     &    abs(exp(-(za(jz,17)/(za(jz,11)*1.5))**6)) )
c
c..net heating rate
c
        za(jz,27) = max( epslon, zgtots(mzones)*1.e-6)
c
c   central electron temperature
c
        za(jz,28) = max( epslon, tes(1,lcentr)*useh)
c
c..part of Tang model from Martha Redi (rgb 9-sep-86)
c
      za(jz,29) = max( epslon,
     &  (1.-exp(-zaq1))*exp((2./3.)*zaq1*xbouni(jz)**2)
     &      *zgtots(jz)/(zgtots(mzones)+epslon) )
c
  99  continue
c
c
      do 103 jz=jzmin,mzones
c
c..ideal ballooning mode stability criterion
c  factor (1. + ( cfutz(298) * p' / p'crit )**cfutz(299) ) ** dfutz*(n,30)
c
      zbs=shear(jz)
      zslpts=-slpts(jz)
      zbq=q(jz)
      zbbaxi=avi(2,8,1)*avi(2,11,1)    ! b-tor axis, internal units
      zbbaxs=zbbaxi*uisb               ! b-tor axis, standard units
      zbrmins=ators(jz,1)*sqrt(b0ref/zbbaxi)  ! minor radius
      zbrmajs=rbtors(jz,1)/zbbaxs      ! major radius, r*b-tor/b-tor-axis
      zbe=zbrmins/zbrmajs              ! minor radius / major radius
      zbk=elong(jz,1)
      zbt=0.5*triang(jz,1)/zbe         ! q in pcy paper
c
c
      zalps=8.0*fcpi*armins(jz,1)*zbq*zbq
     1      /(zbbaxs*zbbaxs*zbe)
      zalps=zalps*max(0.0,zslpts)   ! alpha in pcy paper
      zalpc=1.0/epslon
c
      if (lsweq(4) .gt. 0) then
c
c..use numerical computation of ideal ballooning mode criterion
c  as computed in sbrtn mhdstb called from sbrtn mhdnew
c  and passed through commhd as aprcri(jz,2)
c  on baldur zone boundaries
c
      if (-aprcri(jz,2).gt.epslon)
     &zalpc = - aprcri(jz,2)
     &        * 2.*emu0*zbq**2*armins(jz,1)*usil/(zbbaxi**2*zbe)
c
      else
c
c.. pogutse, chudin, yurchenko ballooning mode stability criterion
c     from soviet journal of plasma physics, 6 (1980) 341-345.
c
c  note the minor radius is sqrt ( toroidal flux / ( pi * b-tor-axis) )
c       the "triangularity" (q in the pcy paper), remains nonzero
c              in general at the magnetic axis.
c
      zbkk=zbk*zbk
      zbkm1=zbkk-1.
      zbkp1=zbkk+1.
      zb3kp1=3.*zbkk+1.
c
      zbaa=12.0*zbs*zbk**1.5*(zbkk+3.0)/(zbkp1*zbkp1)
      zbab=zbkm1*zbkm1*zbkp1/zbkk/(zbk+1)**3
      zba=-(zbaa+zbab)/zb3kp1/4.0
c
      zbba=1.0-3.0*zbkm1*(zbkk-4.0*zbt)/zbkp1/zb3kp1
      zbb = -zbe*(1.0/(zbq*zbq)-zbba)
c
      zbc=zbs*zbs/2.0
      if( zba .ne. 0.0 ) then
      zbb=zbb/zba
      zbc=zbc/zba
c
      zbd=zbb*zbb-4.0*zbc
      if( zbd .ge. 0.0 ) then
c
c     consider only real and positive roots of alpha
c
      zbd=sqrt(zbd)
      if( zbc .le. 0.0 ) zalpc=0.5*(-zbb+zbd)
      if( zbc .gt. 0.0 .and. zbb .le. 0.0 )
     x                   zalpc=0.5*(-zbb-zbd)
      endif
      else
      if( zbb .ne. 0.0 ) zbd=-zbc/zbb
      if( zbd .ge. 0.0 ) zalpc=zbd
      endif
c
      endif   ! end of critical pressure gradient (zalpc) computation
c
      zratio=max(epslon,cfutz(298)*zalps/max(epslon,zalpc))
      zfact=max( epslon, 1.0+zratio**cfutz(299))
      za(jz,30)=zfact
c
 103  continue
c
c
      do 105 jz=jzmin,mzones
c
c..modification to produce zero slope at the magnetic axis
c  1. + ( cemprc(19) / max(epslon,xbouni(jz)) )**cemprc(20)
c
      if ( cemprc(19) .gt. epslon .and. cemprc(20) .gt. epslon ) then
        za(jz,31) = 1.+ (cemprc(19)/max(epslon,xbouni(jz)))**cemprc(20)
      endif
c
 105  continue
c
c
      do 119 jz=jzmin,mzones
c
c  32 - 35 are four more factors for the Tang model, revised
c  version from M. Redi (dps 09-dec-86)
c
c  electron heating factor
c
      za(jz,32) = max( epslon,
     &  (1.-exp(-zaq1))*exp((2./3.)*zaq1*xbouni(jz)**2)
     &          *zetots(jz)/(zetots(mzones)+epslon) )
c
c  ion heating factor
c
      za(jz,33) = max( epslon,
     &  (1.-exp(-zaq1))*exp((2./3.)*zaq1*xbouni(jz)**2)
     &          *zitots(jz)/(zitots(mzones)+epslon) )
c
c  relative ion density
c
      za(jz,34) = max( epslon, rhoins(1,jz)/rhoins(1,lcentr) )
c
c  relative electron density
c
      za(jz,35) = max( epslon, rhoels(1,jz)/rhoels(1,lcentr) )
c
c  geometrical factor for TSC Tang model; reduces to (r/a)**-2 in
c  cylindrical limit: 4 pi**2 r / ( dv/d(xi) xi (grad(xi))**2 )
c
      za(jz,36) = max( epslon,
     &  4.*fcpi**2*(armajs(isep,1)/100.)/(avi(jz,3,1) *
     &          avxib(jz)*avi(jz,6,1)) )
c
 119  continue
c
c---------start of the main do-loop over the radial index, "jz"---------
c
      do 290 jz=jzmin,mzones
c
c
c     mean charge of ionic species.
c
      do 210 ix=lhyd1,lhydn
      do 210 n=1,nmax1
      zzv(n,ix)=1.0
      zzd(n,ix)=1.0
  210 continue
c     the following do-loop is automatically bypassed if "mimp=0".
      do 220 ii=1,mimp
      ix=ii+lhydn
      z0=max(cmean(ii,1,jz),epslon)
      do 220 n=1,nmax1
      zzv(n,ix)=z0**dfutzv(n,1)
      zzd(n,ix)=z0**dfutzd(n,1)
  220 continue
c
c  construction of terms which are products of physical factors
c  formally independent of ionic species, "ix".
c  note:  first two factors zav,zad,zzv,zzd, depend on ionic species ix
c
      do 240 n=1,nmax1
      zv(n)=1.0
      zd(n)=1.0
c
      do 230 k=3,kmax1
c
      if (dfutzv(n,k) .ne. 0.)  zv(n) = zv(n) * za(jz,k)**dfutzv(n,k)
      if (dfutzd(n,k) .ne. 0.)  zd(n) = zd(n) * za(jz,k)**dfutzd(n,k)
c
  230 continue
c
  240 continue
c
c
      do 260 n=1,nmax2
      ze(n)=1.0
      zi(n)=1.0
c
      do 250 k=3,kmax1
c
      if (dfutze(n,k) .ne. 0.)  ze(n) = ze(n) * za(jz,k)**dfutze(n,k)
      if (dfutzi(n,k) .ne. 0.)  zi(n) = zi(n) * za(jz,k)**dfutzi(n,k)
c
  250 continue
c
c..use premultipliers
c
      ze(n) = zgpre(1,n) * ze(n)
      zi(n) = zgpre(2,n) * zi(n)
c
  260 continue
c
c
c     assembly of the semiempirical contributions to the transport
c     coefficients.
c
      do 270 ix=lhyd1,limpn
      vxemps(ix,jz)=0.0
      dxemps(ix,jz)=0.0
      do 270 n=1,nmax1
cpc     @2.7a
      vxemps(ix,jz)=vxsemi(n,ix)*zgpre(3,n)*zv(n)*zzv(n,ix)*zav(n,ix)
     &             +  vxemps(ix,jz)
cpc     @2.1c
      dxemps(ix,jz)=dxsemi(n,ix)*zgpre(4,n)*zd(n)*zzd(n,ix)*zad(n,ix)
     &               + dxemps(ix,jz)
  270 continue
c
c  Tang model switches to remove various contributions;
c  from M. Redi (dps 09-dec-86)
c
c  remove total heating factors,
c
      if ( (lemprc(1) .gt. 0) .and. (lemprc(1) .ne. 3) ) then
        ze(3)=0.
        zi(2)=0.
c
c  and choose between collisional and collisionless chi-e
c  according to value of nu-star
c  default value cemprc(17) = 0.15 set in sbrtn preset file default
c
        if ( (lemprc(1) .gt. 0) .and. (lemprc(1) .ne. 4)) then
        if (znustr.le.(cemprc(17))) ze(1)=0.
        if (znustr.gt.(cemprc(17))) ze(2)=0.
        if (znustr.gt.(cemprc(17))) zi(3)=0.
c
      end if
        end if
c
c  or remove species heating factors
c
      if ((lemprc(1) .gt. 0) .and. ( lemprc(1) .eq. 3 )) then
        ze(1)=0.
        ze(2)=0.
        zi(1)=0.
        zi(3)=0.
      end if
c
c to ckeep chi-e collisional model always on
c
        if ( lemprc(1) .eq. 4 ) then
          ze(2) = 0
          ze(3) = 0
        end if
c
c  turn off eta-i chi-i if eta-i too small
c
      if ( lemprc(1) .gt. 0  .and.  lemprc(1) .ne. 2 ) then
        if (zeta.le.cemprc(15)) then
          zi(1)=0.
          zi(2)=0.
        end if
      end if
c
c       scale chi-i tem as chi-e tem if plasma is collisonless
c
        if( lemprc(1) .eq. 5 ) then
          zln = min(abs(slnes(jqmid)),rmins)
          xisemi(3)=1.5*(1.+1./(zeta+epslon))*zln*xesemi(2)/rmajs
        end if
c
        if ( znustr .gt. cemprc(17) ) xisemi(3)=0.
c
c
c..prints out three regimes of xe, 3 of xi, 1  of  d and v
c
        if(knemp.ne.3) go to  275
c
        zr = armins(jz,1)   ! minor radius as chosen using nrad
c
 
        if(memprt.le.epslon) then
          write(nprint,10534) jz,zr,(xesemi(np)*ze(np),np=1,3)
     &     ,(xisemi(ns)*zi(ns),ns=1,3)
     &     ,vxemps(1,jz),dxemps(1,jz)
        else
          write(nprint,10534) jz,zr,(xesemi(np)*ze(np),np=1,3)
     &      ,(xisemi(ns)*zi(ns),ns=1,3)
     &      ,vxsemi(memprt,1)*zv(memprt)*zzv(memprt,1)*zav(memprt,1)
     &      ,dxsemi(memprt,1)*zd(memprt)*zzd(memprt,1)*zad(memprt,1)
       end if
10534   format(3x,i4,4x,f8.2,4x,8(1pe10.3,2x))
c
        if(jz.eq.mzones) then
          write(nprint, 10535) znustr, zalphn, zeta, zln
10535    format(10x,6hnustr=,f12.3,3x,8halpha-n=,f12.3,3x,
     &      6heta-i=,f12.3,3x,4hzln=,f12.3)
c
        end if
 
c
        if(lemprc(1).ne.0) then
        if(jz.eq.mzones) then
        if(jinvld.gt.2) then
        write(ntychl,536) jinvld
        write(nprint,536) jinvld
 536   format(3x,"warning: ptote<0 and chi-e formulation invalid",
     &  3x,"j=",i2)
c
        end if
        end if
        end if
c
c..print out improved values for dxsemi and vxsemi for density tuning
c
         if(jz.eq.mzones) then
           if((tbi).gt.tmod) then
             zdxs=0.0
             zvxs=0.
             zvstp=0.
             zdstp=0.
             if(rhoels(2,2).gt.(1.1 * snestr))  then
               zvxs=.5 * vxsemi(1,2) * (snestr/rhoels(2,2))
             else if(rhoels(2,2).lt.(.9 * snestr))  then
               zvxs= 2. * vxsemi(1,2) * (snestr/rhoels(2,2))
             else
               zvxs=vxsemi(1,2)
               zvstp=1.
             endif
c            compute line averaged density
c
             zneb = 0.0
             do 272 jlav=lcentr,ledge
               zneb = zneb + rhoels(2,jlav)*dxzoni(jlav)
 272         continue
             if(snebar.eq.0.) go to 275
c
             if(zneb.gt.(1.1 * snebar)) then
               zdxs= 2. * dxsemi(1,1) * (zneb/snebar)
             else if(zneb.lt.(.9 * snebar))   then
               zdxs = .5 * dxsemi(1,1) * (zneb/ snebar)
             else
               zdxs = dxsemi(1,1)
               zdstp=1.
             endif
c
             if(nstep.gt.1)
     1       write(ntychl,10537) zdxs,zvxs
             write(nprint,10537) zdxs,zvxs
             if((zvstp.eq.1.).and.(zdstp.eq.1.))  then
               write(ntychl,10538)
               write(ntychl,10539) zdstp
             endif
           endif
         endif
10539    format(1x,i1)
10538    format(1x,23hdensity profile matches)
10537    format(1x, 12hdxsemi(1,1)=,1pe8.1,3x,12hvxsemi(1,2)=,
     1   1pe8.1)
c
  275    continue
c
c
      xeemps(jz)=0.0
      xiemps(jz)=0.0
      do 280 n=1,nmax2
cpc     @2.1d
      xeemps(jz)=xesemi(n)*ze(n)+xeemps(jz)
      xiemps(jz)=xisemi(n)*zi(n)+xiemps(jz)
  280 continue
c
  290 continue
c
c---------end of the main do-loop over the radial index, "jz"-----------
c
c
c     values at "r=0".
c
      do 300 ix=lhyd1,limpn
      vxemps(ix,lcentr)=vxemps(ix,jzmin)
      dxemps(ix,lcentr)=dxemps(ix,jzmin)
  300 continue
      xeemps(lcentr)=xeemps(jzmin)
      xiemps(lcentr)=xiemps(jzmin)
c
c..additional printout
c
      if (knemp .eq. 3) then
c
      write (nprint,10540) zalphn,zaq1,znustr,zeta
10540 format (1x,'alphn=',1pe13.5,' zaq1=',1pe13.5
     & ,' nustr=',1pe13.5,' eta_i=',1pe13.5)
c
      write (nprint,10550) zetots(mzones)*1.e-6, zitots(mzones)*1.e-6
     & ,zgtots(mzones)*1.e-6, rastar
10550 format (1x,'zetots=',1pe13.5,' zitots=',1pe13.5,' zgtots='
     & ,1pe13.5,' mw,  rastar=',0pf7.2,' cm')
c
c
       endif
c
      return
      end

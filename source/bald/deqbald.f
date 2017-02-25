c/ 16:00 10 July 2001 .../baldur/code/bald/deqbald.f
c/
c/ subroutines for equilibrium-BALDUR interface
c/
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
cinput
c
c      -------------------------------------
c      Fixed-boundary equilibrium input data
c      -------------------------------------
c
c     One equilibrium variable appears in the first BALDUR namelist:
c
c LEQTYP determines which equilibrium package is to be used
c        in the 1-1/2-D BALDUR code.
c        Variable LEQTYP is stored in common COMMHD in file DCOMMON.
c
c leqtyp = 0   (default) for concentric circular cylinder flux surfaces
c              This is for simulating 1-D BALDUR results.
c              Toroidal effects do not appear in the metric elements.
c
c leqtyp = 1   for concentric circular toroidal flux surfaces
c              toroidal effects are included in the metric elements
c              but there is not shift or elongation, for example.
c
c leqtyp = 8   VMOMS moment solution to the Grad-Shafranov Equation
c              by Lao, Wieland, Houlberg and Hirshman, ORNL/TM-7871 (1982)
c              uses r0bt, ahlfbt, elngbt, trngbt as input
c              for the plasma boundary shape as a function of time
c
c leqtyp = 11  for equilbirum flux surfaces computed with the
c              fixed boundary moments code VMEC from S. P. Hirshman
c
c leqtyp = 12  for equilbirum flux surfaces computed with the
c              fixed boundary moments code VMEC2 from S. P. Hirshman
c
c leqtyp = 20  to read in equilibrium on rectangular grid
c              (must have a file from wbald081 wseqdata
c              which must be renamed ???)
c
c leqtyp = 21  use the free-boundary rectangular-coordinate
c              equilibrium solver from Jardin's TSC code
c              (must have an additional data set inserted
c              just after the usual equilibrium data set)
c
c
c     ---------------------
c     Equilibrium namelist:
c     ---------------------
c
c     All the rest of the equilibrium input variables appear
c            in the third namelist,
c            which was added to the usual BALDUR namelists.
c
c     NOTE:  The equilibrium namelist must be preceded by two lines   !
c            of text, each line up to 80 characters long.             !
c
c
c LSWEQ(n), n=1,32, integer switches for program control:
c
c     lsweq(1) controls verbosity of output in baldur
c              = 0 for brief output; = 9 for verbose output
c                increase these values by 10 for output during
c                the initial equilibrium iteration
c
c     lsweq(2) controls the frequency of printed equilibrium output
c
c     lsweq(3) = 0 for flat extrapolation of flux surface averages
c
c     lsweq(4) > 0 for ballooning mode stability analysis called from mhdNEW
c              = n computed every nth time the equilibrium is computed
c
c     lsweq(5)   controls the frequency with which output is sent to
c                file 'stabil' for Alan Todd's PEST stability analysis
c              =  negative, no ouput is sent to file 'stabil'
c              =  0, output sent each time equilibrium is computed
c              =  n>0, output sent every nth time equilibrium is computed
c
c     lsweq(6) = # of equilibrium iterations during initialization
c
c     lsweq(7) = 0 to initialize equilibrium with a smooth
c                  current profile
c     lsweq(7) .ne. 0 to initialize equil. with uniform vloop
c                  except for region with q = qmin
c
c     lsweq(8) = minimum extent of q=qmin region baldur zone bndry index
c
c     lsweq(10) .gt. 0, reset bpoli(j) to provide a reasonalble
c                       q profile to start the process of
c                       computing an equilibrium
c                otherwise, bpoli(j) remains
c                       as it was computed in sbrtn auxval
c     lsweq(10) is used to control initialization of the
c        poloidal magnetic field bpoli on the first pass.
c     lsweq(10) = 0  (default) starts with bpoli as it was set in
c        sbrtn auxval of the original 1-D BALDUR code
c     lsweq(10) = 1  uses quadratic interpolation of q(xi)
c        between q0 = qmin and qedge = 3. * qmin
c     lsweq(10) = 2  uses cubic interpolation of iota between
c        iota edge as determined by toroidal current after
c             computing V'(xi) and <|del xi|**2/R**2> in sbrtn aveflx
c        iota axis = max ( cequil(10) set by equilibrium namelist,
c             minimum iota consistent with d iota / d xi .le. 0.
c
c     lsweq(11)  is used to experiment with new initialization schemes
c       lsweq(11) = 0 (default) old initialization (before June 1989)
c
c     lsweq(19) .gt. -1 (default) to optimize harmonic representation
c                   of plasma boundary shape in sbrtn eqHR11
c                   (applies only when leqtyp = 11 for VMOMS code)
c
c     lsweq(20) = 0 (default) call eqhrtr with each equilibrium
c                   to compute equilibrium harmonics on BALDUR
c                   zone boundaries
c
c     lsweq(21) = 0 (default) call hrmtht with each equilibrium
c                   to compute (R(theta,xi),Y(theta,xi)) on BALDUR
c                   zone boundaries from the equilibrium harmonics
c
c     lsweq(32) is set equal to the incremental equilibrium counter
c               neqinc in sbrtn mhdnew in order to pass the
c               incremental equilibrium counter through common /commhd/
c
c CEQUIL(J), J=1,32, Equilibrium control variables
c
c       cequil(1) = fraction of fast alpha pressure    (1.0)
c                       to be included in the equilibrium plasma pressure
c
c       cequil(2) = fraction of fast beam particle scalar pressure  (1.0)
c                       to be included in the equilibrium plasma pressure
c
c       cequil(3)   defaulted to 1.0
c
c       WARNING:  Stotler had these definitions shifted up by 1 index.
c                 Hence, I default cequil(3) to 1.0 and avoid using it.
c
c
c MFLXS  = number of equilibrium flux surfaces
c          from magnetic axis (j=1) to edge (j=jr)
c
c MOM    = number of poloidal harmonics computed
c
c MOMBND = number of equilibrium harmonics specified for boundary
c
c NGSKIP = flux surfaces skipped in sbrtn eqgeom
c
c..specify time development of equilibrium:
c
c NTBKMX    = number of breakpoint times to be used (up to 20)
c
c TBK(n)    = breakpoint times in seconds (up to 20)
c
c NEQDT(n)  = number of equilibrium recomputations required
c               between tbk(n) and tbk(n+1) breakpoint time
c
c CURMAT(n) = toroidal plasma current (in megamperes) at time tbk(n)
c
c BTORT(n)  = toroidal vacuum magnetic field (in tesla) at time tbk(n)
c                  measured at rtort(n)
c
c RTORT(n)  = major radius (in meters) at time tbk(n)
c                  at which btort(n) is prescribed
c
c..Description of the plasma boundary shape:
c  -----------------------------------------
c
c     There are three ways to specify the plasma boundary shape.
c These are controlled by what variables are set in the third
c (equilibrium) namelist set.
c
c (LEQBND is a variable set in sbrtn mhdINI, file DEQBALD,
c          according to how the boundary shape is specified on input)
c
c     = 0  (default)  for circular cylinder (1-D BALDUR)
c            the major and minor radius are taken from the 1D BALDUR
c            input data
c
c     = 1  for specification of major radius, half-width,
c               elongation, and triangularity
c
c     = 2  for direct harmonic representation
c               (R0BT, RMBT, YMBT)
c
c    Input used only with leqbnd = 1, or with leqtyp=8:
c
c RMAJBT(N) = major radius (in meters) at time tbk(n)
c
c RMINBT(n) = plasma half-width (in meters) at time tbk(n)
c
c ELNGBT(n) = plasma elongation at time tbk(n)
c           = height to point of peak elongation of plasma boundary
c             divided by the half-width of the plasma boundary
c
c TRNGBT(n) = plasma triangularity at time tbk(n)
c           = ( major radius to geometric center of plasma boundary
c             - major radius to point of peak elongation of boundary )
c             divided by the half-width of the plasma boundary
c
c ELONGA = initial guess for elongation at magnetic axis
c
c VMINIT(j), j=1,6 values used to initialize the VMOMS code
c             see sbrtn vmomeq
c
c    Input used with leqbnd = 2:
c
c R0Bt(n)   = boundary value of zeroth harmonic r0 at time tbk(n)
c
c RMBT(n,m) = bndry value for m th harmonic rmb(m) at time tbk(n)
c
c YMBT(n,m) = bndry value for m th harmonic ymb(m) at time tbk(n)
c
c SHIFT   (not stored)    used in sbrtn EQIN11 (DEQINIT)
c           = Initial guess for outward shift of magnetic axis
c               relative to r0bt(1) in meters
c
c QMIN      = minimum allowed initial q-value used when setting bpoli
c               in sbrtn mhdini
c
c
c     Archaic variable names, no longer used:
c     ---------------------------------------
c
c EQERR(n), n=1,8, array of error tolerances, not yet used
c
c LSW(n), n=1,8, integer switches used for program control
c     lsw(1) = 0 for brief output; = 9 for verbose output
c     lsw(2) = 0 for no plots; = 9 for most plots
c     lsw(3) used in stand alone equilibrium code but not here
c     lsw(4) used in stand alone equilibrium code but not here
c     lsw(5) used to be used to select type of equilibrium computation
c               now LEQTYP in first BALDUR namelist is used for this
c
c JR    = number of equilibrium flux surfaces  = MFLXS
c
c**********************************************************************c
 
c      --------------------------------------------------------
c      Internal variables computed by the equilibrium interface
c      --------------------------------------------------------
 
c ....loop-limits; limits for different grids:
c     mflxs ...# of fluxsurfaces from magnetic axis (j=1) to edge (j=jr)
c              mflxs <-- jr=njeq
c     mhrms ...number of poloidal harmonics   (OLD: mom)
c     mombnd...number of equilibrium harmonics specified for boundary
c     ngskip...flux surfaces skipped in sbrtn AVEglb
c     mjbal ...number of radial gridpoints for BALDUR-code
c              mjbal <-- njav
 
c ....PROFILES:
c
c     pz(j)    = scalar pressure at zone centers
c     eiota(j) = rotational transform at zone boundaries
c     psitor(j)= toroidal flux at zone boundaries
c     torcur(j)= toroidal current (in amps) within flux surface j
c     beta(j)  = plasma beta within flux surface j
c
c..parameters:
c
c pj = maximum number of radial gridpoints from equilibrium code
c pm = maximum number of harmonics from equilibrlium code
c pg = dimension of xg,yg used for graphics
c pnf = number of harmonics used in fast fourier transform routines
c pnfh = dimension of complex arrays used in sbrtns rcfft2 and crfft2
c pt = max number of breakpoints in time for changing bndry shape
 
c ....PARAMETERS:
c
c     Kjflx ... maximum number of radial gridpoints for equilibrium code
c     Kmhrm ... maximum number of harmonics for equilibrlium code
c     Knf   ... number of harmonics used in fast fourier transform routines
c     Knfh  ... dimension of complex arrays
c               used in sbrtns rcfft2 and crfft2
c     Kpbt  ... max number of breakpoints in time for changing bndry shape
 
c======================================================================c
 
c     As far as BALDUR is concerned, the purpose of this equilibrium
c package is to compute needed flux surface averages
c given the appropriate radial profiles and boundary conditions.
c
c  Given pressure and rotational transform as a function of toroidal flux,
c this code computes the shape of the equilibrium flux surfaces in the form:
 
c  R(xi,theta) = r0hr(j) + sum of rmhr(m,j) * cos(m*theta)
c                           for m=1 to mom
c  Y(xi,theta) =           sum of ymhr(m,j) * sin(m*theta)
c                           for m=1 to mom
c
c where r0hr, rmhr, and ymhr are computed at
c
c  xi = sqrt (normalized toroidal flux) = eqxibi(j)
c
c... all spatial dimensions are in meters ...
c
c    Output sent by the equilibrium package to file jobxlpt allows
c    the user to monitor equilibrium convergence:
c FSQR = the square of the covariant component of force imbalance
c      along the major radius
c FSQZ = same along the vertical z-axis
c    After these values are printed out once, their normalized values
c  are printed out after each nskip iterations.  (Multiply the values in
c  the table by the values printed above the table for the true values
c  of FSQR AND FSQZ).
c R00(0) is the major radius of the magnetic axis in meters.
c WMHD = volume integral of (B**2/2. + pressure / (gamma - 1.))
c      gamma = 0.0 is hardwired in at the BALDUR interface
c      (adiabatic compression is handled in the BALDUR transport eqns).
c BETA = volume integral of pressure
c      divided by the volume integral of B**2/2
c      note, the beta values defined elsewhere in the BALDUR code
c      are defined as the volume average of pressure
c      divided by the square of the vacuum toroidal magnetic field
c      at RMID divided by 2*emu0
c <M> is a measure of harmonic optimization, currently using pexp=4.
c
c--------1---------2---------3---------4---------5---------6---------7-c
c     note:  avxiz(j) and avxib(j), j=1,njav, are zone centered and
c zone boundary values of the independent variable xi.
c The equilibrium moments code finds values of  R(xi,theta) and y(xi,theta)
c such that xi = sqrt [ (toroidal flux) / (toroidal flux within boundary)]
c In the future, xi may be any flux surface label.
c
c     The computed flux surface averages are all stored
c in the array  avi(j,n,it), where
c j=1,...,njav = number of radial gridpoints
c n=1,...,nnav = number of flux surface average variables computed
c it = 1,..,5 is a time sequence index explained below
c
c avti(it), it=1,5, is a time sequence variable to be explained below
c
c     In each case
c avi(j,n,1) = the present value of the flux surface average
c     at time avti(1) = 0.5*(tai+tbi) * uist for use in BALDUR
c
c avi(j,n,2) = the present extimate of the time rate of change
c     of the flux surface average
c avti(2) = timestep over which this rate was computed
c
c avi(j,n,3) = extrapolated value of flux surface average
c avti(3) = estimated time at which the equilibrium needs
c     to be recomputed
c
c avi(j,n,4) = most recently computed flux surface average at
c avti(4) = time at which last equilibrium was computed
c
c avi(j,n,5) = previously computed flux surface average at
c avti(5) = time at which equilibrium was previously computed
c
c ^ avi(j,n,-)                              * avi(j,n,3)
c |                                  *  *
c |                avi(j,n,4) *      *
c |                     *         * avi(j,n,1) at time eqti(1)
c |             *               *
c |       *                   * avi(j,n,1) when equilibrium was recomputed
c | * avi(j,n,5)
c  --!------------------------!--------------!---------> time
c   eqti(5)                 eqti(4)-eqti(1)->eqti(3)
c
c     The meaning of each column of the matrix avi(j,n,it) is as follows:
c
c avi(j,1,it) = rho = sqrt [ toroidal flux / (pi * b0ref) ]
c        on zone boundaries
c        = 0. at the magnetic axis j=2
c        = - avi(3,1,it) at the interior ghost point j=1
c        note:  b0ref = a fixed reference magnetic field which should not
c           be changed during a computer run
c
c avi(j,2,it) = d rho / d xi   at zone boundaries
c
c avi(j,3,it)= d V(xi,t) / d xi   on zone boundaries
c           where V(xi,t) is the volume of the flux surface
c           labeled by xi = avxib(j)
c
c avi(j,4,it) = d V(xi,t) / d xi   at zone centers
c
c avi(j,5,it) = <|del xi|>   at zone boundaries
c           here, <...> = the flux surface volume average
c           <...> = d volume integral of ... / d V
c
c           |...| is the absolute value of ...
c
c avi(j,6,it) = < |del xi|**2 >   at zone boundaries
c
c avi(j,7,it) = < |del xi|**2 / R**2 >   at zone boundaries
c           here  R = major radius
c
c avi(j,8,it) = R B-toroidal   at zone boundaries
c
c avi(j,9,it) = R B-toroidal   at zone centers
c
c avi(j,10,it) = < 1. / R**2 >   at zone centers
c
c avi(j,11,it) = < 1. / R >  at zone boundaries
c
c avi(j,12,it) = volume within zone boundary
c
c avi(j,13,it) = area within zone boundary
c
c avi(j,14,it) = rmid = ( rout + rin ) / 2  at zone boundaries
c                major radius midway between outer and inner edge of plasma
c
c avi(j,15,it) = ahalf = (rout - rin) / 2  at zone boundaries
c                half width of flux surface measured to its widest point
c
c avi(j,16,it) = elongation = height / width  at zone boundaries
c
c avi(j,17,it) = trianularity = (rtop - qrmid) / qahalf  at zone boundaries
c
c avi(j,18,it) = indentation  at zone boundaries
c              = (minimum major radius at midplane - rin) / qahalf
c
c avi(j,19,it) = sqrt ( flux surface area / pi ) at zone boundaries
c                a root mean square measure of the minor radius in meters
c
c avi(j,20,it) = volume / ( 2 pi area) [ meters ] at zone boundaries
c                major radius to the geometric center of a flux surface
c
c--------:---------:---------:---------:---------:---------:---------:-c
c
c         The following variables in common block commhd are derived
c         from the avi(j,n,it) array.  In each case
c
c         1)  All are in standard units.
c         2)  Use iz=1 for BALDUR zone boundaries,
c                 iz=2 for BALDUR zone centers
c         3)  All are evaluated at the current BALDUR time step.
c
c  arms(j,iz)   = sqrt (area / pi ) at zone bndries/centers
c               = "root mean square" minor radius
c  rgcs(j,iz)   = volume / ( 2 * pi * area)
c               = major radius to the geometric center of a flux surface.
c
c  areas(j,iz)  = cross-sectional area of flux surface
c  vols(j,iz)   = volume within toroidal foux surface
c
c  ators(j,iz)  = sqrt [ toroidal flux / ( pi * b0ref ) ]
c                 note that b0ref = btort(1) is a fixed reference magnetic
c                 field
c
c  ahalfs(j,iz) = half width of flux surface measured to its widest extent
c  rmids(j,iz)  = major radius midway between widest extent of flux surface
c
c  elong(j,iz)  = elongation ( height / width ) of flux surface
c                 (dimensionless)
c  triang(j,iz) = triangularity of flux surface (dimensionless)
c  dent(j,iz)   = indentation of flux surface (dimensionless)
c
c  fltors(j,iz) = toroidal flux
c  flpols(j,iz) = poloidal flux measured from magnetic axis
c
c  rbtors(j,iz) = major radius times toroidal magnetic field
c
cend
c**** ##################################################################
c**** ##################################################################
c
c  -------------------------
c  contents of file DEQBALD:
c  -------------------------
c
c  (1)  Equilibrium initialization and namelist input for 1-1/2-D BALDUR
c
c                   $$ COMMHD
c                   $$ CL1
c                 < $$ CLintf >
c
c         eqINIT ...initialize equilibrium (RGB)............ (--------)
c
c         eqIN01 ...initialize analytical HR
c         eqIN11 ...initialize RGB/Hirshman- moment's equilibrium code
c         eqIN21 ...initialize wos/SCJ-code
c
c         mhdINI ...initialize MHD.......................... (->BALDUR)
c
c
c**** ##################################################################
c
c  (2) Interface routines between BALDUR & equilibria
c
c                   $$ Cbaldr (BALDUR-common blocks)
c                   $$ COMMHD
c                   $$ CL1
c                   $$ CLintf
c
c
c         mhdINI ...initialize MHD.......................... (->BALDUR)
c                   moved to file deqinit
c
c         mhd    ...control-routine-for-MHD-equilibrium..... (->BALDUR)
c
c         mhdBAL ...transfer equilibrium information to BALDUR->BALDUR)
c         mhdNEW ...new MHD-equilibrium..................... (->BALDUR)
c
c         mhdCYL ...circular cylinder  (1D Baldur).......... (->BALDUR)
c
c         trPEST ...transfer file for PEST.................. (->*PEST*)
c
c         tBALFX ...transfer BALDUR-profiles to FLX-grid.... (->p',ff')
c
c**** ##################################################################
c
c  (3)  Equilibrium control routines;
c
c                   $$ COMMHD
c                   $$ CL1
c                 < $$ CLintf >
c
c
c         eqCTRL ...control routine for equilibria and mapping;averages
c
c         eqHRA1 ...analytical HR
c         eqHR11 ...RGB/Hirshman- moment's equilibrium code
c         eqxy21 ...wos/SCJ-code                ( *** DUMMY *** )
c
c         eq21HR ...mapping from (x,y) to HR    ( *** DUMMY *** )
c         eqHRtr ...transpose from rmhr,ymhr to eqrcjm,eqysjm
c
c         eqout  ...print output of equilibria
c
c**** ##################################################################
c
c  (4)  Flux surface average routines for 1-1/2-D BALDUR transport code
c
c                   $$ COMMHD
c                   $$ CL1
c                 < $$ CLintf >
c
c
c         AVEQHR ...flux surface averages based on HR standard
c         AVEset ...set grid..................................
c         AVEflx ...volume averages on flux surfaces
c         AVEpla ...flusurface-averages of plama quantities
c         AVEglb ...eqgeom-->global geometric paramters.......
c         AVEavi ...store  averages in AVI(j,n,4)
c
c**** ##################################################################
c
c  (5)  MHD stability analysis package
c
c          AVEflx ...add del.xi calculation....................(->mhdSTB)
c
c          mhdNEW ...add call to mhdstb........................(->mhdSTB)
c
c          tfPEST ...retain psi, dp/dpsi, dq/dpsi calculations.(->mhdSTB)
c
c          mhdSTB ...driver for stability analysis.............(->mhdSTB)
c
c          balgam1...Infinite n ideal ballooning mode criterion
c                       routines by Michael W. Phillips 02-oct-87
c
c
c$*** ###################################################################
c
c@eqINIT  .../baldur/code/bald/deqbald.f
c rgb 24-nov-99 no longer open files msg, uout, and stabil
c rgb 24-jul-96 removed rewind 1 and installed read (1,inpt,err=90,end=90)
c rgb 20.31 15-jun-92 leqtyp=11 and leqtyp=12 are now the same option
c     both use the new VMEC2 code
c rgb 20.27 23-feb-92 csawth(9)=1.0 multiplier for Monticello sawtooth period
c rgb 20.11 26-aug-91 update default values of csawth
c rgb 20.04 26-jul-91 leqtyp=12 provided to use VMEC2 code from Hirshman
c rgb 18.80 06-dec-90 removed calls to letout, letini, letter, letend
c     , printx, and prline which were in file dutil0
c rgb 18-dec-89 default cequil(3) to 1.0 since Stotler used cequil(3)
c     to represent the fraction of fast beam particle scalar pressure
c rgb 01-oct-89 removed initialization of qhalf(1)
c rgb 18-jul-89 default lsawth(1)=10 and csawth(5)=1. to old model
c     vminit(1-6) implemented to control initialization of VMOMS package
c rgb 24-may-89 got rid of sbrtn eqin01 through eqin21
c rgb 23-may-89 added rmajbt, rminbt, elngbt, trngbt, leqbnd
c        reworked headers
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqINIT *
c     *********************
c
c     usage:
c           CALL eqINIT
c
c     purpose:
c        Sbrtn eqINIT initializes the MHD equilibrium computation
c        - use only the equilibrium common blocks cl1 and commhd
c     1)  set defaults
c     2)  print headlines
c     3)  read equilibrium namelist
c     4)  establish equilibrium grid eqxibi(j) and eqxizi(j)
c     5)  compute initial equilibrium
c
c***********************************************************************
c
      subroutine eqINIT
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
      include 'clsaw.m'
      include 'clparm.m'
      include 'clislb.m'
c
c------------
c
c..temporary 1-D arrays for Lower Hybrid Current Drive and Power
c  dimensions are 55 x 20 = 1100
c
cbate      real zcdrive(1100), zplhcde(1100)
c
c------------
      namelist /inpt/  lsweq, mflxs, mom, mombnd, ngskip
     &  , ntbkmx, tbk, neqdt, r0bt, rmbt, ymbt, curmat, btort, rtort
     &  , leqbnd, rmajbt, rminbt, elngbt, trngbt, elonga, vminit
     &  , ftol, relerr, deltin, gam, gamma, itmom, nskip
     &  , eqerr, qmin, shift, cequil
     &  , lsawth, csawth, swton, swtoff, swperd, swqmin, swxmin
     & , nmodes, modem,  moden,  qmode,  qb1amp, qb1scl
     & , qdeltx, qdeltq, qxwall, nxiham, nhamad, qpeakc, qshldc, qwexc
     & , cislnd, lislnd
     & , cdcoef, plhcoef, xcd, tcd, cdrive, plhcde
     & , lsw, jr

c     ------------------------------------------------------------------
c
c     ==================================================================
c.....output units used by equilibrium codes deqave and deqmom
c
cbate      open (66,file='msg',status='unknown')
      open (7,file='eqout',status='unknown')
cbate      open (8,file='uout',status='unknown',form='unformatted')
cbate      open (9,file='stabil',status='unknown')
c
c     ==================================================================
c
c
CLLL  0.  set defaults:
c
         nout  = 6                    ! baldur printout
         neqout= 7                    ! printout to file eqout
         ntty  = 5                    ! printout to terminal
         pi = atan2 (0.,-1.)
         twopi = 2. * pi
         emu0 = pi * 4.e-7
	 sqprofile  = ' '
	 scdprofile = ' '
c
      do 3 i= 1,32
         lsweq (i)= 0
         cequil(i)= 0.
         eqerr (i)= 0.
    3 continue
c
c..initialize cequil(1)=1.0 and cequil(2)=1.0
c       cequil(1) = fraction of fast alpha pressure
c                       to be included in the equilibrium plasma pressure
c       cequil(2) = fraction of fast beam particle scalar pressure
c                       to be included in the equilibrium plasma pressure
c
       cequil(1) = 1.0
       cequil(2) = 1.0
       cequil(3) = 1.0
c
      do 4 jt=1,Kpbt
        tbk(jt)    = jt - 1     ! bndry shape time breakpoints
        neqdt(jt)  = 1
        r0bt(jt)   = 0.
        curmat(jt) = 0.
        btort(jt)  = 0.
        rtort(jt)  = 0.
        rmajbt(jt) = 0.
        rminbt(jt) = 0.
        elngbt(jt) = 0.
        trngbt(jt) = 0.
   4  continue
c
      elonga = 0.0
      vminit(1) = 0.1
      vminit(2) = -0.4
      vminit(3) = 1.0
      vminit(4) = 0.0
      vminit(5) = 0.01
      vminit(6) = 2.0
c
      do 5 jt=1,kpbt
        do 5 jm=1,kmhrm
          rmbt(jt,jm) = 0.
          ymbt(jt,jm) = 0.
   5  continue
c
         neqdim = Kjbal     ! dimension of eqtmp, ceqtmp,deqtmp,eeqtmp
         mdim = Kmhrm       ! dimension to be passed as argument
         jdim = Kjflx       ! dimension to be passed as argument
         njav = Kjbal       ! max number of BALDUR radial gridpoints
         nnav = Knave       ! number of types of flux surface averages
c
         mflxs = 21         ! number of radial grid points
         mom = Kmhrm        ! number of harmonics
         ntbkmx = 1         ! maximum number of time breakpoints
c
c.....equ1 (moment's code) only:
c
      leqbnd = 0            ! harmonic representation of moments
      mombnd = 2            ! harmonics specified for bndry shape
      ngskip = 1            ! to skip flux surfaces in sbrtn AVEglb
c
      itmom = 2000          ! max iterations allowed in sbrtn eqmom1
      nskip = 500
      ftol  = 1.e-7         ! abs tolerance used in sbrtn eqmom1
      relerr = 1.e-5
      gam   = 5. / 3.       ! used in sbrtn ave to define entz(j)
      gamma = 0.            ! used in sbrtn eqmom1
      initin = 0            ! to initialize sbrtn eqmom1
      deltin = 2.           ! timesteping scale factor in sbrtn eqmom1
      qmin   = 1.e-4        ! minimum allowed q-value used in mhdini
c
      shift = 0.            !  local default, shift of magnetic axis
c
c..Archaic variables:
c  ------------------
c
      do 2 i=1,8
    2    lsw(i) = 0         ! control switches
c
         jr = -1            ! archaic number of radial grid points
c
c  --------------------------------------------------
c..sawtooth variables used in sbrtn sawmix and sawmx1
c  --------------------------------------------------
!cap
      lsawth = 0
      csawth = 0.
      swton  = 1.e10
!
cbate      swton  = 0.       ! sawtooth time on [sec]        <-- cfutz(460)
      swtoff = 0.       ! sawtooth time off [sec]       <-- cfutz(461)
      swperd = 0.       ! sawtooth period [sec]         <-- cfutz(462)
      swqmin = 0.       ! minimum q-value trigger       <-- cfutz(477)
      swxmin = 0.       ! minimum xbouni(jz) for trigger
c
      lsawth(1) = 10    ! default to old sawtooth model
      csawth(5) = 1.0   ! transfer magnetic energy to electrons by default
      csawth(9)  = 1.0  ! multiplier for Monticello sawtooth period
      csawth(10) = 1.0  ! Rogers-Zakharov sawtooth trigger model
      csawth(11) = 1.0
      csawth(14) = 0.03 ! lower bound on shear for sawtooth crash
c
c  -------------------------------------------------------
c..Saturated tearing mode package initialization and input
c  -------------------------------------------------------
c
      nmodes    = 1             ! number of perturbed helical modes
      modem     = 0
      moden     = 0
      modem(1)  = 2             ! m
      moden(1)  = 1             ! n

      qmode     = 0.
      qpeakc    = 0.
      qshldc    = 0.
      qwexc     = 0.

cbate      qhalf(1)  = 0.1           ! island half-width as fn of q-value
      qb1amp    = 0.
      qb1scl    = 0.
      qdeltx    = 0.01          ! xi starting point for ode's
      qdeltq    = 1.e-3         ! qdeltaq * qedge = smallest island used
      qxwall    = 1.0           ! r_wall / r_plasma
c
      nxiham    = 21            ! number of hamada flux surfaces
      nhamad    = 4             ! number of hamada harmonics
c
      qpeakc(1) = 0.            ! peaking factor within island
c
      do 12 j=1,32
        cislnd(j) = 0.
        lislnd(j) = 0
  12  continue
c
      cislnd(25) = 0.25         ! 1./u**2 part off wings of island
      cislnd(26) = 0.12         ! u**4 part of interior modification
      cislnd(29) = 0.01         ! extent of analytic part of strtn idiftr
      cislnd(31) = 1.e-5        ! tolerance of the shooting integrator
      cislnd(32) = 1.e-2        ! tolerance of the nonlinear eqn solver
c
      lislnd(32) = 11           ! max number of iterations allowed
c
c
c..Defaults for Lower Hybrid Current Drive and Power deposition
c  --------
c
      nxlhcd = 0
      ntlhcd = 0
c
      cdcoef = 1.0
      plhcoef = 1.0
c
      do jx=1,Kjbal
        xcd(jx) = -1.0
        cdprof(jx) = 0.0
        plhprof(jx) = 0.0
      enddo
c
      do jt=1,Knave
        tcd(jt) = -1.0e10
      enddo
c
      do jt=1,Knave
        do jx=1,Kjbal
          cdrive(jx,jt) = 0.0
          plhcde(jx,jt) = 0.0
        enddo
      enddo
c
c
c     ==================================================================
c
c
CLLL  1.  Headlines:
c
CLL   1.1 read headlines:
c
  100 format (a)
  110 format ('  4245 .baldur64 & .bequil64 (+yequilib) RGB,wos 86-04')
  120 format (1x,10i10)
c
crgb      call letout (neqout)
crgb      call letini (1)
crgb      call letter (31,9,3,1,'#',' BALDUR$')
crgb      call letend (0)
crgb      call printx ('                1 1/2 D  version  (may-86)$',4)
crgb      call printx ('                EQUILIBRIUM.......printout$',9)
c
      read (1,100) txtequ(2)
      read (1,100) txtequ(3)
c     ..................................................................
c
CLL   1.2 headline for different equilibrium packages:
c
crgb      call prline       ('leqtyp$$$',5,10)
      write (neqout,120)  leqtyp
c
      if (leqtyp. eq. 21)  then
CL    (21)  (x,y)-equilibrium:
      txtequ(1)= ' Free boundary equilibium - not implemented yet .....'
c
      else if (leqtyp. eq. 20)  then
CL    (20)  (x,y)-equilibrium:
      txtequ(1)= ' Read in equilibrium on rectangular grid ............'
c
      else if (leqtyp. eq. 12)  then
CL    (12)  Hirshman's equilibrium moments code code VMEC2:
      txtequ(1)= ' Hirshmans VMEC2 equilibrium moments code ...........'
c
      else if (leqtyp. eq. 11)  then
CL    (11)  Hirshman's equilibrium moments code code VMEC2:
      txtequ(1)= ' Hirshmans VMEC2 equilibrium moments code ...........'
c
      else if (leqtyp .eq. 8)   then
CL    (8)   VMOMS equilibirum moments code by Lao et al 1982:
      txtequ(1)= ' VMOMS equilibirum moments code by Lao et al 1982....'
c
      else if (leqtyp .eq.  1)  then
CL    ( 1)  analytical MOMENTs specification:
      txtequ(1)= ' analytical specification of equilibrium moments ....'
c
      else
c
        if ( leqtyp .ne. 0 ) then
          leqtyp = 0
      write (nout,*) ' WARNING:  value of LEQTYP changed to 0 (default)'
        endif
c
CL    (00)  concentric circular cylinder:  1D BALDUR
      txtequ(1)= 'leqtyp = 0:1D BALDUR: concentric circular cylinder...'
c
      endif
c     ------------------------------------------------------------------
c
CLL   1.3 print Headlines:
c
   20 continue
         nnnn= nout
      do 21 jjn= 1,2
        if (jjn.gt. 1)  nnnn= neqout
        write (nnnn,110)
        write (nnnn,*) txtequ(1)
        write (nnnn,*) txtequ(2)
        write (nnnn,*) txtequ(3)
        write (nnnn,*)
   21 continue
c     ==================================================================
c
c
CLLL  2.  INPUT:
c
c
c     ------------------------------------------------------------------
c     ..................................................................
c      namelist /inpt/  lsweq, mflxs, mom, mombnd, ngskip
c     &  , ntbkmx, tbk, neqdt, r0bt, rmbt, ymbt, curmat, btort, rtort
c     &  , leqbnd, rmajbt, rminbt, elngbt, trngbt, elonga, vminit
c     &  , ftol, relerr, deltin, gam, gamma, itmax, itmom, nskip
c     &  , eqerr, qmin, shift, cequil, eqerr
c     &  , lsawth, csawth, swton, swtoff, swperd, swqmin, swxmin
c     & , nmodes, modem,  moden,  qmode,  qb1amp, qb1scl
c     & , qdeltx, qdeltq, qxwall, nxiham, nhamad, qpeakc, qshldc, qwexc
c     & , cislnd, lislnd
c     & , lsw, jr
c     ..................................................................
c
cbate960724      rewind 1
c
      read (1,inpt,err=90,end=90)
c
cbate960724 30   continue
cbate960724      rewind 1
c     ------------------------------------------------------------------
      if ( jr .gt. 0 ) mflxs = jr    ! archaic variable jr
         mhrms= mom
c     ------------------------------------------------------------------
c     ==================================================================
c
c
CLLL  3.  transfer information to BALDUR commons {cbaldr,COMMHD}  and
CL             local equilibr. commons {cl1} and interface {clintf}:
c
          mjbal= njav                ! preliminary setting (->commhd)
          mhrms= mom
c
c     ==================================================================
c
c
CLL   (3) establish equilibrium grid; fluxsurface-grid:
CL
CL       eqxibi(j), j=1,mflxs ...normalized equilibrium flux label
CL                                  at equilibrium zone boundaries
CL       eqxizi(j)            ... same zone centered
CL
c
c     These variables are physically equivalent
c  to the variables xbouni and xzoni in the rest of the BALDUR code.
c  They all stand for the square root of the toroidal flux normalized
c  to its value at the edge of the plasma.
c  As such, they are the fundamental independent variable against
c  which all other variables are defined.
c
c..For Hirshman's new equilibrium moments codes MHALF or VMEC2,
c  these arrays are spaced on equal intervals of toroidal flux
c  (that is, eqxibi(j)**2 are equally spaced).
c  This has to be forced when leqtyp = 11 or 12
c
      if ( leqtyp .eq. 12 .or. leqtyp .eq. 11 ) then
        leqxi = 1
        mflxs = 31
        mom   = 5
c
        write (nout,202) leqxi, mflxs, mom
 202    format (// '  WARNING:  The following values have been changed'
     &,/t12,'leqxi changed to ',i3
     &,/t12,'mflxs changed to ',i3
     &,/t12,'mom   changed to ',i3
     &,/t12,'because the Hirshman equilibrium code VMEC2 is being used'
     &,//)
c
      endif
c
cl..eqxibi(j)  equally spaced in toroidal flux
c
      if (leqxi .eq. 1)  then
c
         jaxis = 1
         zfac = 1. / real (mflxs - jaxis)
         eqxibi(1) = 0.0
       do 40 j=2,mflxs
         eqxibi(j) = sqrt ( real (j - jaxis) * zfac )
         eqxizi(j) = sqrt ( real (j - 1.5  ) * zfac )
   40  continue
         eqxibi(1) = 0.0
         eqxizi(1) = - eqxizi(2)
c
cl..eqxibi and eqxizi equally spaced in xi=sqrt(norm toroidal flux)
c
      else   !  (leqxi .ne. 1)
c
         jaxis = 1
         zfac = 1. / real (mflxs - jaxis)
       do 42 j=1,mflxs
         eqxibi(j) = real (j - jaxis) * zfac
         eqxizi(j) = real (j - 0.5 - jaxis) * zfac
   42  continue
c
      endif
c
      do 44 j=1,mflxs
      equixi(j) = eqxibi(j)      ! array in common comisl
  44  continue
c
c
c     ==================================================================
c
c..Set Lower Hybrid Current Drive arrays
c
c..count number of radial grid points
c
      nxlhcd = 0
      do jx=1,Kjbal
        if ( xcd(jx) .gt. -0.9 ) then
          nxlhcd = jx
        else
          go to 51
        endif
      enddo
 51   continue
c
      ntlhcd = 0
      do jt=1,Knave
        if ( tcd(jt) .gt. -99.0 ) then
          ntlhcd = jt
        else
          go to 52
        endif
      enddo
 52   continue
c
c..transfer current drive and power deposition from 1-D to 2-D arrays
c  only if nxlhcd > 0 and ntlhcd > 0
c
cbate      if ( nxlhcd .gt. 0 .and. ntlhcd .gt. 0 ) then
cbate        do jx=1,nxlhcd
cbate          do jt=1,ntlhcd
cbate            cdrive(jx,jt) = zcdrive( jx + (jt-1)*nxlhcd )
cbate            plhcde(jx,jt) = zplhcde( jx + (jt-1)*nxlhcd )
cbate          endif
cbate        endif
cbate      endif
c
c..output to file jobxlpt
c
      if ( nxlhcd .gt. 0 .and. ntlhcd .gt. 0 ) then
c
        write (6,*)
        write (6,*) " Lower Hybrid Current Drive and Heating Power"
        write (6,*) " After namelist inpt in sbrtn eqinit"
        write (6,*)
        write (6,*) " nxlhcd  = ", nxlhcd
        write (6,*) " ntlhcd  = ", ntlhcd
        write (6,*) " cdcoef  = ", cdcoef
        write (6,*) " plhcoef = ", plhcoef
        write (6,*)
        it = ( 1 + ntlhcd ) / 2
        write (6,*) " xcd(jz)    cdrive(jz,",it,")"
        do jz=1,nxlhcd
          write (6,*) xcd(jz), cdrive(jz,it)
        enddo
c
      endif
c
c
c     ==================================================================
c
c..Print out equilibrium namelist
c
      write (nout  ,inpt)
      write (neqout,inpt)
c
c     ==================================================================
c
c
CLLL  4.  compute initial equilibrium:
c
         lprt = lsweq(1)    ! printout controlled by lsweq(1)
         lprt = lprt - 10   ! reduced output frequency during startup
         lsweq(1) = lsweq(1) - 10
c
      call mhdINI
c
         lprt = lprt + 10          ! to control frequency of printout
         lsweq(1) = lsweq(1) + 10
c
      return
c
   90 continue
      call abortb (6,'Error reading equilibrium namelist inpt')
c
      return
      end
c@mhdINI   .../baldur/code/bald/deqbald.f
c  les 01-nov-90 d3he - quartic interpolation of q based on total current
c  rgb 23-aug-90 18.57 new conversion from geometric to harmonic form
c  rgb 08-jul-90 18.46 call eqfbnd before call mhdNEW
c  rgb 01-oct-89 force new initialization scheme if leqtyp = 8
c  rgb 14-jun-89 lsweq(11) .ne. 0 for new initialization scheme
c  rgb 07-jun-89 temporarily remove call mhdBAL, getchi(1) and getchi(2)
c  rgb 06-jun-89 set avi(j,jn,1) and call getchi just after establishing
c     an initial analytic equilibrium shape
c     Removed estimate of toroidal flux and avi(j,1,1) just after that
c     Removed recalculation of avi(j,jn,1) just before equil iteration
c  rgb 05-jun-89 set avi(j,jn,1) before iterating on initial equilibrium
c  rgb 24-may-89 Changed consistency checks to allow specification of
c   geometric representation of boundary with rmajbt, rminbt, elngbt,...
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE mhdINI *
c     *********************
c
c     usage:
c           CALL mhdINI
c
c     purpose:
c           mhdINI initializes the variables subsequently used in
c        sbrtn mhd
c        only common blocks cbaldr and commhd are used here
c     1)  check consistency of input equilibrium data
c     2)  iterate to compute initial equilibrium
c           with desired poloidal magnetic profile
c
c***********************************************************************
c
      subroutine mhdINI
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
      dimension zr2(Kjbal)
      data neqinc /0/
c     ------------------------------------------------------------------
c
      write (nout,101)
  101 format (//' sbrtn mhdINI called to initialize BALDUR '
     & ,'equilibrium variables before calling sbrtn mhd')
c     ------------------------------------------------------------------
c
c..initialize scalars
c
         njav = mzones + 1
         mjbal= njav
         nnav = 20        ! number of types of flux surface averages
      umil = usil * 100.  ! length from meters to baldur internal units
      umie = ueie         ! energy from joules to baldur internal units
      umib = usib * 1.e4  ! magnetic field from tesla to internal units
      umii = usii * 3.e9     ! current from amps to internal units
c
      call blines (1)
      call mesage ('equilibrium initialization')
      call ivar ('mjbal',mjbal)
      call ivar ('nnav',nnav)
      call rvar ('umil',umil)
      call rvar ('umie',umie)
      call rvar ('umib',umib)
c
c..set independent variable avxi* stored in common commhd
c
      do 12 j=1,mjbal
         avxiz(j) = xzoni(j)         ! zone centered values of xi
         avxib(j) = xbouni(j)        ! zone boundary values of xi
  12  continue
c
c..check consistency of equiliblrium data
c
      ntbkmx = max (ntbkmx,2)  ! must have at least two breakpoints
      ntbkmx = min (ntbkmx,20)         ! max 20 elements in array tbk
c
c..the value of leqbnd is set according to what fixed-boundary
c  equilibrium data is given
c
      leqbnd = -1
      if ( rmini .gt. 0. .and. rmaji .gt. 0. ) leqbnd = 0
      if ( rmajbt(1) .gt. 0. .and. rminbt(1) .gt. 0.
     &   .and. elngbt(1) .gt. 0. ) leqbnd = 1
      if ( r0bt(1) .gt. 0. .and. rmbt(1,1) .gt. 0.
     &   .and. ymbt(1,1) .gt. 0. ) leqbnd = 2
c
      if ( leqbnd .lt. 0 ) call abortb (nout
     & ,'major and minor radius not defined in sbtn mhdini')
c
c..if no equilibrium boundary is given in the third namelist
c  use circular plasma boundary defined by the old 1-D BALDUR
c
      if ( leqbnd .eq. 0 ) then
        write ( nout, 102) rmaji, rmini
 102    format (//' Values of rmaji = ',1pe12.5
     &    ,'  and rmini = ',1pe12.5,
     &    /' from the fist BALDUR namelist were'
     &    /' used to define major and minor radius in sbrtn mhdini'/)
c
        do 14 j=1,ntbkmx
          rmajbt(j) = rmaji
          rminbt(j) = rmini
          elngbt(j) = 1.0
          trngbt(j) = 0.0
          r0bt(j)   = rmaji
          rmbt(j,1) = rmini
          ymbt(j,1) = rmini
  14    continue
c
c..If Fourier harmonics of the boundary shape are not given,
c  then use a geometric specification of the boundary shape
c  and compute harmonics at breakpoints as in the VMOMS code by Lao
c
      elseif ( leqbnd .eq. 1 ) then
        write ( nout, 104)
 104    format (//' Values of rmajbt, rminbt, elngbt, and trngbt'
     &    /' used to compute harmonics r0bt, rmbt, and ymbt '
     &    ,' in sbrtn mhdini'/)
c
c..convert from geometric to harmonic representatio
c  R(theta) = R0 + R1 cos(theta) + R2 cos(2 theta)
c  Y(theta) = Y1 sin(theta)
c
        do 16 j=1,ntbkmx
          rmbt(j,1) = rminbt(j)
          rmbt(j,2) = rminbt(j) * 0.5 * trngbt(j)
          ymbt(j,1) = rminbt(j) * elngbt(j)
          r0bt(j)   = rmajbt(j) - rminbt(j) * 0.5 * trngbt(j)
          write(nout,*) j, rminbt(j), rmajbt(j), elngbt(j), trngbt(j)
  16    continue
c
      endif
c
c..plasma current has to be given in the third namelist
c
      if (curmat(1) .eq. 0.)
     &  call abortb (nout,'curmat(1) .eq. 0. in sbrtn mhdini')
c
c..if rtort and btort are not given in the third namelist,
c  then use major radius and toroidal magnetic field
c  from the old 1-D BALDUR code
c
      if (rtort(1) .le. 0.) then
        write (nout,106) rmaji
 106    format (//' Changed rtort(1) to rmaji = ',1pe12.5
     &  ,' in sbrtn mhdini'/)
        do 17 j=1,ntbkmx
          rtort(1) = rmaji
  17    continue
      endif
      if (rtort(1) .le. 0.)  call abortb (nout,'rtort(1) .le. 0.')
c
      if (btort(1) .eq. 0.) then
        write (nout,107) bzi
 107    format (//' Changed btort(1) to bzi = ',1pe12.5
     &  ,' in sbrtn mhdini'/)
        do 18 j=1,ntbkmx
          btort(j) = bzi
  18    continue
      endif
      if (btort(1) .eq. 0.)  call abortb (nout,'btort(1) .le. 0.')
c
      if (tbk(1) .gt. tinit) then
         tbk(1) = tinit
         write (nout,*) ' tbk(1) reset to tinit = ',tbk(1)
      endif
c
c  the array tbk(i) must be in ascending order for i=1,ntbkmx
c  if an element of array tbk(i) is found not in ascending order
c  that value of tbk(i) will be reset to tmax + dtmax
c  and ntbkmx will be reset to that index i
c
      imax = ntbkmx
      do 20 i=2,imax
c
      if (tbk(i) .ge. tmax + dtmax) then
      ntbkmx = i
      go to 21
      endif
c
      if (tbk(i) .le. tbk(i-1)) then   ! test for ascending order
      tbk(i) = tmax + dtmax
      ntbkmx = i
      write (nout,121) i,i,tbk(i),ntbkmx
 121  format (/' tbk(',i3,') found to be out of ascending order'
     &  /' tbk(',i3,') reset to ',1pe12.4,
     &  /' and ntbkmx reset to ',i3)
      go to 21
      endif
c
  20  continue
  21  continue
c
c  make sure last element in tbk(i) is beyond end of run time
c
      if (tbk(ntbkmx) .lt. tmax+dtmax) then
         if (ntbkmx .lt. 20) ntbkmx = ntbkmx + 1
         tbk(ntbkmx) = tmax + dtmax
      write (nout,124) ntbkmx,tbk(ntbkmx)
 124  format (/' found tbk(ntbkmx) .lt. tmax + dtmax'
     &  /' reset ntbkmx =',i3,' and tbk(ntbkmx) = ',1pe12.4)
      endif
c
c  at the plasma boundary, we must have r0bt(i), rmbt(i,m), ymbt(i,m)
c  all greater than zero and curmat .ne. 0.
c  if any of these condtions are violated, reset all of these variables
c  equal to the value they had at the last breakpoint time
c
      do 30 i=2,ntbkmx
      if (r0bt(i) .le. 0. .or. rmbt(i,1) .le. 0. .or. ymbt(i,1) .le. 0.
     &  .or. curmat(i) .eq. 0.) then
      r0bt(i) = r0bt(i-1)
      rmbt(i,1) = rmbt(i-1,1)
      ymbt(i,1) = ymbt(i-1,1)
      curmat(i) = curmat(i-1)
      write (nout,131) i,r0bt(i),rmbt(i,1),ymbt(i,1),curmat(i)
 131  format (/' for breakpoint index i = ',i3,
     & /' r0bt(i) or rmbt(i,1) or ymbt(i,1) was found to be nonpositive'
     & /' or curmat(i) was found to be nonzero'
     & /' r0bt(i) was reset to ',1pe12.4
     & /' rmb(1,i) was reset to ',1pe12.4
     & /' ymb(1,i) was reset to ',1pe12.4
     & /' curmat(i) was reset to ',1pe12.4)
      endif
      if (btort(i) .eq. 0.) btort(i) = btort(i-1)
      if (rtort(i) .le. 0.) rtort(i) = rtort(i-1)
  30  continue
c
c..the vacuum toroidal magnetic field is passed through variable rbtort
c
      do 38 i=1,ntbkmx
      rbtort(i) = rtort(i) * btort(i)  ! R Btoroidal in internal units
  38  continue
c
c=======================================================================
c
c
CLL   establish an initial equilibrium
c      (before any plasma profiles are known, use an analytic estimate)
c
c
c  ntb = present value of breakpoint time interval
c
      ntbk = 1           ! initial time interval index
c
      eqtime = tinit     ! initital time
c
c..vacuum toroidal magnetic field and corresponding major radius
c
      r0ref = rtort(1)
      b0ref = btort(1)
c
      write (nout,140) r0ref, b0ref
 140  format (/' Values set for r0ref and b0ref in sbrtn mhdini'
     & /' r0ref = ',1pe15.5,/' b0ref = ',1pe15.5/)
c
c..establish initial boundary shape and toroidal current (interpolate)
c
      do 39 jm=1,mhrms
        rmb(jm) = 0.
        ymb(jm) = 0.
  39  continue
c
      call eqfbnd
c
c..use a simple analytic harmonic representation to compute an
c  initial equilibrium before the pressure and q(r) profiles are known
c
      if ( lsweq(11) .ne. 0 .or. leqtyp .eq. 8 ) then
c
c..new initialization schemes
c
        if (leqtyp .eq.  0)  then
          call mhdcyl
        else
          call eqhra1
          call aveqhr(6,7)
        endif
c
          do 46 jn=1,nnav
          do 46 j=1,mjbal
            avi(j,jn,1) = avi(j,jn,4)
  46      continue
c
cbate      call mhdBAL   ! to set BALDUR variables in common {cBALDr}
cbate      call getchi(1)
cbate      call getchi(2)
c
      else
c
c..old initialization scheme
c
        if (leqtyp .eq.  0)  then
          call mhdcyl
        else
          call eqhra1
        endif
c
      endif  ! end of new vs old initialization scheme
c
c     ------------------------------------------------------------------
c
c  bpoli(j) remains as it was computed in sbrtn auxval at this point
c
c..estimate toroidal flux and avi(j,1,1) = rho
c
      torflx = pi * rmb(1) * ymb(1) * b0ref   ! cbate
      rhoedg = sqrt (torflx / (pi * b0ref))   ! cbate
c
      do 41 j=1,mjbal  ! cbate
      avi(j,1,1) = xbouni(j) * rhoedg  !  initial estimate for rho  ! cbate
  41  continue  ! cbate
c
c..if lsweq(10) .gt. 0, reset bpoli(j) to provide a reasonalble
c  q profile to start the process of computing an equilibrium
c  otherwise, bpoli(j) remains as it was computed in sbrtn auxval
c
      if (lsweq(10) .eq. 0) then
c
c  Leave bpoli(j) as it was set in sbrtn auxval in the 1-D BALDUR code
c
      elseif (lsweq(10) .eq. 1) then
c
c  Quadratically interpolate q(xi) from q0 = qmin to qedge = 3. * qmin
c
      zbnorm = umib * b0ref / r0ref
      zq0    = qmin
      zqedge = 3. * qmin
      zxnorm = 1. / (mjbal - 2)
      do 42 j=1,mjbal
      zx = xbouni(j)
      bpoli(j) = zbnorm * avi(j,1,1) / ( zq0 + (zqedge-zq0) * zx**2 )
  42  continue
c
      elseif (lsweq(10) .ge. 2) then
c
c..Call aveflx to compute V'(xi) and <|del xi|**2 / R**2>,
c  use these to estimate
c  zf = V'(xi) B0REF <|del xi|**2/R**2> / ( xi 2 pi emu0)
c  Determine the coefficients of
c  iota(xi) = ziota0 + ziota1 * xi**2 + ziota2 * xi**4
c  from the conditions
c  ziotae = iota-edge = plasma current / ( rho**2 zf)
c  d plasma current / d xi = 0 at plasma edge (xi=1)
c  ziota0 = max ( cequil(10) input through the equilibrium namelist,
c     minimum iota-axis consistent with d iota / d xi everywhere)
c  note, zlamf = ( xi / zf )  d zf / d xi
c  ziota2 = 2. * ( (6. + zlamf) * ziotae / 4. - ziota0 ) at plasma edge
c  ziota4 = ziota0 - (4. + zlamf) * ziotae / 2. at plasma edge
c
c  Note:  in order to use this option, the equilibrium harmonics
c     eqrcj0(j), eqrcjm(m,j), and eqysjm(m,j) must be initially set
c     using some reasonable approximation to the flux surface shapes.
c
      call aveqhr(6,7) !  to compute V'(xi) and <|del xi|**2 / R**2>
c
      zf  = avi(mzones,3,4) * b0ref * avi(mzones,7,4)
     &   / ( xbouni(mzones) * twopi * emu0 )
      ziotae = eqcamp / ( rhoedg**2 * zf )
      zf1 = avi(mzones-1,3,4) * b0ref * avi(mzones-1,7,4)
     &             / ( xbouni(mzones-1) * twopi * emu0 )
      zlamf = 2. * (zf - zf1)
     &   / ( (zf + zf1) * (xbouni(mzones) - xbouni(mzones-1)) )
c
      ziota0 = ( 6. + zlamf ) * ziotae / 4.
c
c  this is the min value of iota0 set by d iota / d xi .le. 0. everywhere
c
      ziota0 = max ( ziota0 , cequil(10) )
c
      ziota2 = 2. * ( (6.+zlamf) * ziotae / 4. - ziota0 )
      ziota4 = ziota0 - (4.+zlamf) * ziotae / 2.
c
      zbnorm = umib * b0ref / r0ref
c
      do 44 j=3,mjbal
      zx = xbouni(j)
      bpoli(j) = zbnorm * avi(j,1,1)
     &             * (ziota0 + zx**2 * ( ziota2 + zx**2 * ziota4 ) )
  44  continue
c
      bpoli(2) = 0.
      bpoli(1) = - bpoli(3)
c
c   les - d3he - quartic interpolation of q based on total current
c
      elseif (lsweq(10) .eq. 3) then
c
      zbnorm = umib*b0ref/r0ref
      zq0=qmin
      zqedge=5.*ymb(1)**2/ekappa*b0ref/(r0ref*eqcamp*1.e-06)*
     1  (1.+(1.+2.*edelta)*ekappa**2)/(2.*ekappa)
      zqalph=zqedge/zq0 - 2.
      do 43 j=3,mjbal
        zx=xbouni(j)
        bpoli(j) = zbnorm*avi(j,1,1) / (zq0*(1.+zqalph*zx**2+zx**4))
 43   continue
      bpoli(2)=0.
      bpoli(1)=-bpoli(3)
c
      endif
!
 98   continue
c=======================================================================
c
c
c..now iterate to establish initial equilibrium
c
      itrmax = max ( 1, lsweq(6) )
c
      do 68 iter=1,itrmax
c
      call eqfbnd
c
      call mhdNEW     ! compute equilibrium and flux surface averages
c
      do 49 n=1,nnav  ! copy computed flux surface averages into avi(j,n,1)
      do 49 j=1,mjbal
        avi(j,n,1) = avi(j,n,4)
  49  continue
c
      call mhdBAL   ! to set BALDUR variables in common {cBALDr}
c
      call getchi(1)
      call getchi(2)         ! to compute eta(2,j)
c
c  compute zr2(j) = rho'(xi) * V'(xi) * <|del xi|**2/R**2> * r0ref
c
      do 50 j=1,mzones+1
        zr2(j) = avi(j,2,1)*avi(j,3,1)*avi(j,7,1)*r0ref
  50  continue
c
c..if lsweq(6) < 1, skip rest of bpoli initialization.
c  Then the default values of bpoli(j) set in sbrtn auxval, are left
c
      if (lsweq(6) .lt. 1) go to 60   ! skip rest of bpoli initialization
c
      if (lsweq(7) .eq. 0) then
c
c..initialize bpoli with a smooth current profile
c  matching qmin at the magnetic axis
c  matching current=eqcamp at the plasma edge (j=mzones)
c  requiring d current / d rho = 0 at the plasma edge
c  use a cubic polynomial to join the current in the q=qmin region
c  with the boundary condition at the plasma edge
c  let zx = avi(j,13,1) / avi(mzones,13,1)  where avi(j,13,1) = area
c    be the flux label used in this cubic polynomial fit
c
      bpoli(3) = b0ref * avi(3,1,1) / ( r0ref * qmin )
      bpoli(2) = 0.
      bpoli(1) = - bpoli(3)
c
      zxnorm = 1. / avi(mzones,13,1)
      zcnorm = 1. / ( twopi * emu0 )
c
      zx1 = 0.
      zcur1 = 0.
      zcf1 = eqcamp
c
      zx2 = avi(3,13,1) * zxnorm
      zcur2 = bpoli(3) * zr2(3) * zcnorm
      zcf2 = ( eqcamp - zcur2 ) / (1.-zx2)**2
c
      za3 = ( zcf1 - zcf2 ) / ( zx2 - zx1 )
      za2 = - zcf1 - za3 * ( 2. + zx1 )
      za1 = - 2. * za2 - 3. * za3
      za0 = eqcamp - za1 - za2 - za3
c
      do 52 j=4,mzones-1
c
      zx = avi(j,13,1) * zxnorm
      zcur = za0 + zx * ( za1 + zx * ( za2 + zx * za3 ))
c
c  zbc = bpoli(j) computed using this current profile
c  zbq = bpoli(j) computed using q = qmin
c  choose the smallest of these two and continue
c  note:  bpoli(j) is always assumed to be positive
c  temporarily, lsweq(8) marks the minimum extent of the
c    q = qmin region as a function of baldur zone bndry index
c
      zbc = twopi * emu0 * zcur / zr2(j)
      zbq = b0ref * avi(j,1,1) / ( r0ref * qmin )
c
      if (zbq .lt. zbc .or. j .lt. lsweq(8) ) then
c
      bpoli(j) = zbq
c
c  now we have to compute the coefficients of a new cubic polynomial
c  to extend from the edge of the q = qmin regio to the edge of the plasma
c
      zx1 = zx2
      zcf1 = zcf2
c
      zx2 = avi(j,13,1) * zxnorm
      zcur2 = bpoli(j) * zr2(j) * zcnorm
      zcf2 = ( eqcamp - zcur2 ) / (1.-zx2)**2
c
      za3 = ( zcf1 - zcf2 ) / ( zx2 - zx1 )
      za2 = - zcf1 - za3 * ( 2. + zx1 )
      za1 = - 2. * za2 - 3. * za3
      za0 = eqcamp - za1 - za2 - za3
c
      else ! bpoli computed from fitted current profile
c
      bpoli(j) = zbc
c
      endif
c
  52  continue
c
      do 53 j=mzones,mzones+1
      bpoli(j) = twopi * emu0 * eqcamp / zr2(j)
  53  continue
c
      else
c
c..initialize bpoli(j)
c  so that vloopi(j,2) conforms to a prescribed profile initially
c  beam driven current is not considered for now
c  qmin = the minimum initial q-value allowed
c         with absolute minimum value 1.e-4
c  note:  bpoli(j) must be set positive
c
      zbrmax = b0ref / ( r0ref * max(qmin,1.e-4) )
c
      do 56 j=2,mzones
c
      vloopi(j,2) = 1.   ! change to prescirbed proflile later
c
      bpoli(j+1) = ( zr2(j) * bpoli(j) / avi(j,8,1)
     &        + avi(j,4,1) * dxzoni(j) * vloopi(j,2)
     &        / (twopi * usir * eta(2,j) * avi(j,9,1)) )
     &        * avi(j+1,8,1) / zr2(j+1)
c
  56  continue
      bpoli(2) = 0.
      bpoli(1) = - bpoli(3)
c
c..now compute edge value of bpoli based on total toroidal plasma current
c  and rescale bpoli and vloopi throughout
c
      jb = mzones + 1
      zbedge = twopi * emu0 * eqcamp
     &   / ( avi(mzones,2,1)*avi(mzones,3,1)*r0ref*avi(mzones,7,1) )
c
      zscale = zbedge / bpoli(mzones)
c
      do 58 j=1,jb
      bpoli(j) = min (bpoli(j) * zscale , zbrmax * avi(j,1,1) )
      vloopi(j,2) = vloopi(j,2) * zscale
  58  continue
      bpoli(2) = 0.
      bpoli(1) = - bpoli(3)
c
      endif        ! end of initializing bpoli choices
c
  60  continue
c
c..make sure bpoli(mzones+1) is defined
c  define it so that the current density in zone mzones is zero
c
      bpoli(mzones+1) = bpoli(mzones) * zr2(mzones) / zr2(mzones+1)
c
c..compute initial ajtori(j,2) and vloopi(j,2) arrays
c  note: ajtori(j,2) = <J-tor / R> (Amp/m**3) at zone center j
c        vloopi(j,2) = loop voltage (volts) at zone center j
c
      do 64 j=2,mzones-1
c
      ajtori(j,2) = ( zr2(j+1) * bpoli(j+1) - zr2(j) * bpoli(j) )
     &             / ( emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
c
c
      vloopi(j,2) = twopi * usir * eta(2,j) * ( avi(j,9,1)
     &      * (zr2(j+1)*bpoli(j+1) / avi(j+1,8,1)
     &         - zr2(j) * bpoli(j) / avi(j,8,1) )
     &      / (emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
     &   - cjbeam * ajbs(j) * usij * 0.5*(avi(j+1,11,1)+avi(j,11,1)) )
     &  / avi(j,10,1)
c
  64  continue
c
      vloopi(1,2) = vloopi(2,2)
      ajtori(1,2) = ajtori(2,2)
c
      bpoli(2) = 0.
      bpoli(1) = - bpoli(3)
c
      ajtori(mzones,2) = 0.
      vloopi(mzones,2) = vloopi(mzones-1,2)
c
c..intermediate output
c
      if (lsweq(1) .gt. 5) then
      write (nout,162)
 162  format (/' loop voltage (volts) at BALDUR zone boundaries')
      write (nout,164) (vloopi(j,2),j=1,jb)
 164  format (1p10e12.4)
c
      write (nout,166)
 166  format (/' < toroidal current density (amps/m**3) / R >'
     & /'  at zone centers')
      write (nout,164) (ajtori(j,2),j=1,jb)
c
      endif   !   end of intermediate printout
c
  68  continue     ! end of iteration loop
c
c=======================================================================
c
c..having just computed the initial flux surface averages
c  initialize the flux surface average array avi(j,n,it)
c  with steady state values
c  and set avti(it) to recompute equilibrium after first BALDUR timestep
c
      do 82 n=1,nnav
      do 82 j=1,mjbal
         avi(j,n,1) = avi(j,n,4)
         avi(j,n,3) = avi(j,n,4)
         avi(j,n,5) = avi(j,n,4)
         avi(j,n,2) = 0.           ! time derivative of avi
   82 continue
c
c..set avti(it) to start run
c
         avti(1) = tinit * ueit
         avti(2) = (tbk(2) - tbk(1)) * ueit / max(neqdt(1),1)
         avti(3) = avti(1) + avti(2)
         avti(4) = avti(1)
         avti(5) = avti(1)
c
c..set BALDUR variables found in common cbaldr
c
      call mhdBAL
c     ==================================================================
c
c
c..printout
c
      write (nout,110) neqinc,avti(4)
 110  format (/i5,' new equilibrium computed at time',1pe12.4)
      write (nout,112) (avti(i),i=1,5)
 112  format (/'        avti = ',1p5e12.4)
      write (nout,114) (avi(mjbal,3,i),i=1,5)
 114  format (' V''(edge,it) = ',1p5e12.4//)
c
      return
      end
c@MHD   .../baldur/code/bald/deqbald.f
c  rgb 14-aug-96 avti(1) + zdavti ---> (avti(1) + zdavti)*rndup
c    to control round off error in decision to compute equilibrium
c  rgb 26-may-89 set lprt from lsweq(1)
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE MHD    *
c     *********************
c
c     usage:
c           CALL MHD
c
c     purpose:
c           Sbrtn mhd updates the flux surface average arrays
c           and recomputes a new equilibrium if necessary
c     note that the internal equilibrium common blocks are not present
c
c***********************************************************************
c
      subroutine mhd
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
      data neqinc /0/
c     ------------------------------------------------------------------
          mjbal= njav
          lprt = lsweq(1)
c     ------------------------------------------------------------------
c
CLL   ****** update done every time sbrtn mhd is called
c
c
c..update current values of all flux surface averages
c
         eqtime = 0.5 * (tai + tbi) * uist
c
         zdavti = 0.5*(tai+tbi) - avti(1)  ! time since last update
      do 10 n=1,nnav
      do 10 j=1,mjbal
         avi(j,n,1) = avi(j,n,1) + avi(j,n,2) * zdavti
   10 continue
         avti(1) = 0.5*(tai+tbi)           ! update current time
c
c..find time interval for interpolation of plasma shape and current
c
c  first estimate of eqdts which is the timestep between
c  successive equilibrium computations
c
         ntbk = min (ntbk, ntbkmx-1)
         ntbk = max (ntbk, 1)
         eqdts = (tbk(ntbk+1) - tbk(ntbk)) / max (neqdt(ntbk),1)
c
c..interpolate prescibed boundary shape from values at breakpoint times
c
      call eqfbnd
c
c     ------------------------------------------------------------------
c
c..set BALDUR variables found in common cbaldr
c
      call mhdBAL
c
c=======================================================================
c
CLLL  test to determine if it is time to recompute equilibrium
c
      lequil = 0   ! default, do no recompute equilibrium
c
      if ( (avti(1) + zdavti)*rndup .gt. avti(3)) lequil = 1
c
      if (lequil .eq. 0) return
c
c
CLL   ****** equilibrium recomputed only when necessary
c
      call mhdNEW
c
c     ------------------------------------------------------------------
c
c..determine new equilibrium timestep eqdts and store eqdts in avti(2)
c  eqdts should be no smaller than the BALDUR diffusion timestep tai-tbi
c  the first step into a new time interval should be a short one
c  determined by the minimum of the old and new stepsize.
c
      eqdts=min(eqdts, (tbk(ntbk+1) - tbk(ntbk))/max (neqdt(ntbk),1))
c
c..if needed, eqdts could be adjusted here if the equilibrium
c  is changing too rapidly or too slowly
c     zc = abs( (avi(mjbal,3,4)-avi(mjbal,3,5))
c    &  / ( (abs(avi(mjbal,3,4))+epslon) * errmax ) )
c     eqdts = eqdts / min ( (1.+zc), 2. )
c
      eqdts = abs ( max ( uist * abs(tbi-tai) , eqdts ) )
c
      avti(4) = avti(1)
      avti(2) = eqdts * usit
      avti(3) = avti(4) + avti(2)
c
c..linearly extrapolate from already computed flux surface averages
c  to the time at which we expect to recompute the equilibrium
c
      if (lsweq(3) .eq. 0 .or. avti(4)-avti(5) .eq. 0.) then
c
         zfac = 0.  ! flat extrapolation
c
      else   ! use linear extrapolation
         zfac = (avti(3) - avti(4)) / (avti(4) - avti(5))
      endif
c
      do 42 n=1,nnav
      do 42 j=1,mjbal
         avi(j,n,3) = avi(j,n,4) + (avi(j,n,4) - avi(j,n,5)) * zfac
   42 continue
c
c..decide whether or not these newly extrapolated flux surface average
c  values make sense
c
c
c..estimate time rate of change of flux surface average values
c  based on change from current value avi(j,n,1)
c  to extrapolated values avi(j,n,3)
c
         avti(2) = avti(3) - avti(1)
c
         zfac = 1. / avti(2)
      do 52 n=1,nnav
      do 52 j=1,mjbal
         avi(j,n,2) = (avi(j,n,3) - avi(j,n,1)) * zfac
   52 continue
c     ------------------------------------------------------------------
c
c..printout
c
      if (lprt .lt. 0) return
c
      write (nout,110) neqinc,avti(4)
 110  format (/i5,' new equilibrium computed at time',1pe12.4)
c
      if (lprt .lt. 3) return
c
      write (nout,112) (avti(i),i=1,5)
 112  format ('  avti = ',1p5e12.4)
      write (nout,113) (avi(mjbal,1,i),i=1,5)
 113  format ('rho(edge,it)  = ',1p5e12.4)
      write (nout,114) (avi(mjbal,2,i),i=1,5)
 114  format (' rho''(xi,..)  = ',1p5e12.4)
      write (nout,115) (avi(mjbal,3,i),i=1,5)
 115  format ('  V''(edge,it) = ',1p5e12.4)
      write (nout,116) (avi(mjbal,5,i),i=1,5)
 116  format (' < del xi  >  = ',1p5e12.4)
      write (nout,117) (avi(mjbal,6,i),i=1,5)
 117  format ('<(del xi)**2> = ',1p5e12.4)
      write (nout,118) (avi(mjbal,7,i),i=1,5)
 118  format ('<(dlxi/R)**2> = ',1p5e12.4)
      write (nout,119) (avi(mjbal,10,i),i=1,5)
 119  format ('<(1. / R)**2> = ',1p5e12.4)
c
      return
      end
c@eqfbnd   .../baldur/code/bald/deqbald.f
c  rgb 03-sep-90 18.59  zeroed out rmb and ymb before resetting
c  rgb 07-jun-89 improved logic for time beyond range of breakpoints
c  GB 19-nov-86 created eqfbnd
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqfbnd *
c     *********************
c
c     usage:
c           CALL eqfbnd
c
c     purpose:
c           eqfbmd interpolates fixed boundary shape, rbedge, and eqcamp
c        from values at breakpoint times
c
c     The following variables are set:
c        ntbk, r0b, rmb(jm), ymb(jm), rmajb, rminb, elongb, tringb,
c        rbedge, eqcamp
c
c***********************************************************************
c
      subroutine eqfbnd
c
      include 'cparm.m'
      include 'commhd.m'
      include 'clintf.m'
c
c     ------------------------------------------------------------------
c
      do 5 jm=1,kmhrm
        rmb(jm) = 0.0
        ymb(jm) = 0.0
   5  continue
c
c..special case with only one breakpoint or time before first breakpoint
c
      if ( ntbkmx .le. 2 .or. eqtime .lt. tbk(1) ) then
        ntbk   = 1
        rmajb  = rmajbt(1)
        rminb  = rminbt(1)
        elongb = elngbt(1)
        tringb = trngbt(1)
        rbedge = rbtort(1)
        eqcamp = curmat(1) * 1.e6
        r0b    = r0bt(1)
        do 10 jm=1,mombnd
          rmb(jm) = rmbt(1,jm)
          ymb(jm) = ymbt(1,jm)
  10    continue
        return
c
c..time beyond the last breakpoint
c
      else if ( eqtime .gt. tbk(ntbkmx) ) then
        ntbk   = ntbkmx
        rmajb  = rmajbt(ntbkmx)
        rminb  = rminbt(ntbkmx)
        elongb = elngbt(ntbkmx)
        tringb = trngbt(ntbkmx)
        rbedge = rbtort(ntbkmx)
        eqcamp = curmat(ntbkmx) * 1.e6
        r0b    = r0bt(ntbkmx)
        do 22 jm=1,mombnd
          rmb(jm) = rmbt(ntbkmx,jm)
          ymb(jm) = ymbt(ntbkmx,jm)
  22    continue
        return
c
c..find interpolation interval between breakpoints
c
      else
c
         ntbk = min (ntbk, ntbkmx-1)
         ntbk = max (ntbk, 1)
c
         i1 = ntbk
      do 32 i=i1,ntbkmx-1
         ntbk = i
      if (eqtime .lt. tbk(i+1)) go to 34
   32 continue
   34 continue
c
c..establish linear interpolation parameter ztbint
c
      if (tbk(ntbk+1)-tbk(ntbk) .gt. 0.) then
         ztbint = (eqtime-tbk(ntbk)) / (tbk(ntbk+1)-tbk(ntbk))
         ztbint = max (ztbint,0.)
         ztbint = min (ztbint,1.)
      else
         ztbint = 0.
      endif
c
c..linearly interpolate boundary shape in time
c
       rmajb  = rmajbt(ntbk) + (rmajbt(ntbk+1) - rmajbt(ntbk)) * ztbint
       rminb  = rminbt(ntbk) + (rminbt(ntbk+1) - rminbt(ntbk)) * ztbint
       elongb = elngbt(ntbk) + (elngbt(ntbk+1) - elngbt(ntbk)) * ztbint
       tringb = trngbt(ntbk) + (trngbt(ntbk+1) - trngbt(ntbk)) * ztbint
c
         r0b = r0bt(ntbk) + (r0bt(ntbk+1) - r0bt(ntbk)) * ztbint
      do 16 m=1,mombnd   ! number of boundary moments
         rmb(m) = rmbt(ntbk,m) + (rmbt(ntbk+1,m) - rmbt(ntbk,m))*ztbint
         ymb(m) = ymbt(ntbk,m) + (ymbt(ntbk+1,m) - ymbt(ntbk,m))*ztbint
   16 continue
c
c..vacuum toroidal magnetic field [tesla] time major radius [meter]
c
      rbedge = rbtort(ntbk) + (rbtort(ntbk+1) - rbtort(ntbk)) * ztbint
c
c  eqcamp = total plasma current in amps
c
      eqcamp=(curmat(ntbk)+(curmat(ntbk+1)-curmat(ntbk))*ztbint)*1.e6
c
      endif
c
      return
      end
c@mhdBAL   .../baldur/code/bald/deqbald.f
c  rgb 23-jun-96 compute dvoli(j) from j=1,mjbal-1, not to mjbal
c  les  nov-90  volume element for VOLINT, volume integration subr
c  GB 19-nov-86 after statement label 26,
c               changed symmetry condition in CUBINT to 1 (flpols, even)
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE mhdBAL *
c     *********************
c
c     usage:
c           CALL mhdBAL
c
c     purpose:
c           mhdBAL  resets BALDUR variables found in common block cbaldr
c     which may change each time the flux surface averages are updated
c
c     The following variables are used:
c        mzones, uis*, bpoli(j)
c        mjbal, pi, twopi, r0ref, b0ref, avi(j,i,n),
c                                       avxib(j), avxiz(j), rbedge
c
c***********************************************************************
c
      subroutine mhdBAL
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
      include 'comhd3.m'
c     ------------------------------------------------------------------
c...<<<   mjbal= njav
c     ------------------------------------------------------------------
c
c..set BALDUR variables found in common cbaldr
c
         z2voli = 1. / ( 2. * avi(mzones,12,1) )   ! 1/(2* total volume)
      do 22 j=1,mzones
         dx2i(j)   = ( avi(j+1,12,1) - avi(j,12,1) ) * z2voli
         dx2inv(j) = 1. / dx2i(j)
   22 continue
c
         rmaji = avi(mzones,14,1)
         rmajs = rmaji * uisl
         rmini = avi(mzones,15,1)
         rmins = rmini * uisl
         ellipt = avi(mzones,16,1)
c
c..toroidal magnetic field
c
         bzi = rbedge / rmaji
         bzs = bzi * uisb
c
c..BALDUR flux surface arrays in standard units
c
         zuares = uisl**2
         zuvols = uisl**3
c
      do 20 j=1,mjbal
         arms(j,1) = avi(j,19,1) * uisl
         rgcs(j,1)   = avi(j,20,1) * uisl
c
         ators(j,1)  = avi(j,1,1) * uisl
c
         ahalfs(j,1) = avi(j,15,1) * uisl
         rmids(j,1)  = avi(j,14,1) * uisl
c
         elong(j,1)  = avi(j,16,1)
         triang(j,1) = avi(j,17,1)
         dent(j,1)   = avi(j,18,1)
c
         rbtors(j,1) = avi(j,8,1) * uisl * uisb
         rbtors(j,2) = avi(j,9,1) * uisl * uisb
c
c     ...arrays needed to facilitate using sbrtns from Seidl
         vprimi(1,j) = avi(j,3,1) / rmini ! d vol / d r
         vprimi(2,j) = avi(j,4,1) / rmini ! d vol / d r  center
c
   20 continue
c
      do j=1,mjbal-1
         dvoli(j) = avi(j+1,12,1) - avi(j,12,1)      ! vol within zone j
      enddo
c
c   les  nov-90  volume element for VOLINT, volume integration subr
c
      do 21 j=1,mzones
  21   dvols(j)=dvoli(j)*uisl**3
c
c     ..................................................................
c
c..interpolate to BALDUR zone centers
c
      call cubint (avxib,arms(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,arms(1,2),mjbal,0,0.,-1
     & ,'interpolating arms from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,rgcs(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,rgcs(1,2),mjbal,0,0.,+1
     & ,'interpolating rgcs from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,ators(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,ators(1,2),mjbal,0,0.,-1
     & ,'interpolating ators from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,ahalfs(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,ahalfs(1,2),mjbal,0,0.,-1
     & ,'interpolating ahalfs from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,rmids(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,rmids(1,2),mjbal,0,0.,+1
     & ,'interpolating rmids from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,elong(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,elong(1,2),mjbal,0,0.,+1
     & ,'interpolating elong from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,triang(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,triang(1,2),mjbal,0,0.,-1
     & ,'interpolating triang from BALDUR zone bndries to zone centers')
c
      call cubint (avxib,dent(1,1),mjbal,0,ceqtmp,Kjbal
     & ,avxiz,dent(1,2),mjbal,0,0.,-1
     & ,'interpolating dent from BALDUR zone bndries to zone centers')
c     ..................................................................
c
c..consistently derived variables on both zone boundaries and zone centers
c
         zflt = pi * b0ref * uisb
      do 24 j=1,mjbal
         areas(j,1) = pi * arms(j,1)**2
         areas(j,2) = pi * arms(j,2)**2
c
         vols(j,1)   = twopi * rgcs(j,1) * areas(j,1)
         vols(j,2)   = twopi * rgcs(j,2) * areas(j,2)
c
         fltors(j,1) = zflt * ators(j,1)**2
         fltors(j,2) = zflt * ators(j,2)**2
   24 continue
c     ..................................................................
c
c..compute poloidal flux flpols(j,iz) from bpoli(j) at zone bndries
c  bpoli(j) = d (poloidal flux) / d (rho)  / ( 2 * pi * r0ref )
c             at BALDUR zone boundaries, in internal units
c
      flpols(2,1) = 0.       ! at magnetic axis
c
      call cubint (ators(1,1),bpoli(1),mjbal,0,ceqtmp,Kjbal
     & ,ators(2,1),flpols(2,1),mzones,3,0.,-1
     & ,' integrating poloidal flux flpols(j,1) from bpoli(j)')
c
      flpols(1,1) = flpols(3,1)        ! even symmetry across axis
c
      zf = twopi * r0ref * uisl * uisb
      do 26 j=1,mzones+1
      flpols(j,1) = zf * flpols(j,1)
  26  continue
c
      call cubint (ators(1,1),flpols(1,1),mjbal,0,ceqtmp,Kjbal
     & ,ators(1,2),flpols(1,2),mzones,0,0.,1
     & ,' interpolating flpols from BALDUR zone bndries to centers')
c     ------------------------------------------------------------------
c
c..choice of minor and major radius to be used in sbrtn empirc,...
c
      if (nrad .eq. 0) then
c
      do 30 i=1,2
      do 30 j=1,mjbal
      armins(j,i) = ahalfs(j,i) ! halfwidth of flux surface
      armajs(j,i) = rmids(j,i)  ! major radius to midpoint of flux surface
  30  continue
c
      elseif (nrad .eq. 1) then
c
      do 31 i=1,2
      do 31 j=1,mjbal
      armins(j,i) = arms(j,i)  ! sqrt(area/pi)
      armajs(j,i) = rgcs(j,i)  ! volume / ( 2*pi*area)
  31  continue
c
      elseif (nrad .eq. 2) then
c
      do 32 i=1,2
      do 32 j=1,mjbal
      armins(j,i) = ators(j,i)  ! sqrt [ toroidal flux / (pi*b0ref)]
      armajs(j,i) = rbtors(j,i) / b0ref
  32  continue
c
      else
c
      call abortb (nout,'nrad out of range in sbrtn mhdBAL')
c
      endif  ! end of choice of minor and major radius
c
      return
      end
c@mhdNEW   .../baldur/code/bald/deqbald.f
c  rgb 10-jul-01  delxi2d, delth2d, ejacob2d in common /come2d/
c    interpolated from equilibrium grid to baldur zone boundaries
c  rgb 08:30 14-nov-93 comment out call trpest to save space
c  rgb 14.10 09-apr-88  reinstate original form of isweq6
c        to avoid premature call to sbrtn mhdstb
c  rgb 14.09 07-apr-88  compute zpsi, zpresm, zshear, zq, zqpri, zprpri
c       for use in trpest and mhdstb
c  rgb 12.70 07-oct-87  added section to compute arrays needed for the
c       saturated tearing mode package
c       uses cliches clparm and clislb
c  rgb 12.66 02-sep-87 assign lsweq(32) = equilibrium counter
c    dps 07-jan-87 assign zpresm to totprs calculated in subroutine getchi
c    rgb 28-sep-86 removed sect (3) "establish equilibrium grid..."
c               since the equilibrium grid is already correctly
c               estabilished in sbrtn mhdini and should not be changed.
c    dps 24-sep-86 replace zpresm calculation with thrprs, add shear
c                  to argument list for mhdSTB.
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE mhdNEW *
c     *********************
c
c     usage:
c           CALL mhdNEW
c
c     purpose:
c           recompute equilibrium and flux surface averages
c
c     Sbrtn mhdNEW transfers information from the BALDUR common blocks
c     to the equilibrium common blocks commhd in file dcommon
c     note that the internal equilibrium common blocks are not present
c
c***********************************************************************
c
      subroutine mhdNEW
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
      include 'clparm.m'
      include 'clislb.m'
c
c     ------------------------------------------------------------------
c
      dimension zcav(Kjbal,3),zpresm(Kjbal,2),ziotab(Kjbal,2)
      dimension zpsi(Kjbal,2),zq(Kjbal,2),zqpri(Kjbal,2),zprpri(Kjbal,2)
      dimension zshear(Kjbal)
      dimension  wospp(Kjbal),wospz(Kjbal),wospb(Kjbal)
c
      real :: ztemp2d(kjbal,kfour)     ! temporary 2-D array
c     _______________________________________________________________
c
c     zcav(j,i)= temporary array for spline coefficients on BALDUR grids
c
c       In the following arrays, iz=1 refers to BALDUR zone bndries
c                                iz=2 refers to BALDUR zone centers
c
c     zpresm(j,iz) = total pressure on BALDUR zone bndries/centers
c                        in mks units
c
c     ziotab(j,iz) = rotational transform of BALDUR zone bndaries/cntrs
c
c     zpsi(j,iz)   = poloidal flux / ( 2 pi ) on BALDUR zone bnds/cntrs
c                  = flux stream function used in Grad-Shafranov equ.
c                    in MKS units
c
c     zq(j,iz)     = q-value on BALDUR zone bndries/centers
c
c     zqpri(j,iz)  = d q / d psi
c
c     zprpri(j,iz) = d p / d psi
c
c     zshear(j)    = shear at BALDUR zone bndries
c                  = (r/q) d q / d r   where r = armins(j,1)
c                    note that the array shear(j) was computed in
c                    sbrtn xscale using the smoothed q-values
c                    Here, we use the q(j) array computed in sbrtn getchi
c     ------------------------------------------------------------------
c
      data neqinc /0/
c     ------------------------------------------------------------------
           neqinc = neqinc + 1  ! equilibrium counter
           lsweq(32) = neqinc   ! to pass equilibrium counter
c                               ! through common block commhd
c...<<<   mjbal= njav
c     ------------------------------------------------------------------
c
c     ..................................................................
cbate      call timeused (icp0,ic0,isys0)
c     ..................................................................
c
c
c--   call page
      write (nout,100) neqinc, avti(1)*uiet
      write (ntychl,100) neqinc, avti(1)*uiet
 100  format (/,2x,i3,' update equilibrium in sbrtn mhdNEW at time'
     & ,1pe12.4,' sec')
c
      if (lsweq(1) .gt. 3) then
      write (nout,102) r0b,(rmb(m),m=1,mombnd)
      write (ntychl,102) r0b,(rmb(m),m=1,mombnd)
 102  format (' r0b =',0pf7.3,' rmb =',0p6f7.3)
c
      write (nout,104) (ymb(m),m=1,mombnd)
      write (ntychl,104) (ymb(m),m=1,mombnd)
 104  format (t13,' ymb =',0p6f7.3)
      endif
c     ==================================================================
c
c
CLL   (1) save old flux surface averages
      do 20 n=1,nnav
      do 20 j=1,mjbal
         avi(j,n,5) = avi(j,n,4)
   20 continue
         avti(5) = avti(4)
c     ==================================================================
c
c
CLL   (2) establish grid for flux surface averages;  BALDUR-grid:
c
         zfac = 0.1
      do 32 j=1,mjbal
         avxiz(j) = xzoni(j)         ! zone centered values of xi
         avxib(j) = xbouni(j)        ! zone boundary values of xi
c
CL..compute total pressure on BALDUR grid (zone centers) zpresm(j,2)
c
cl   zpresm(j,2) is set equal to totprs(j) (zone centers)
cl   which is computed in subroutine getchi as the sum of
cl      electron + ion thermal pressure
cl      cequil(2) times the fast alpha particle pressure and
cl      cequil(3) times the fast beam ion pressure
cl              which, for now, is treated as a scalar pressure
cl              (2./3.)*cequil(3)*hebems(j)*rhobis(2,j)*zfac
c
cl    note:  zfac is to convert from pressure in standard (cgs) units
cl                to pressure in mks units  nt / m**2
c
         zpresm(j,2)=totprs(j)*zfac
c
CL    define poloidal flux in MKS-units; at zone-Zenters
c
         pspflx(j)= usib*(usil**2) * flpols(j,2)
   32 continue
c
         avxiz (mjbal)   = 2 * avxiz(mjbal-1) - avxiz(mjbal-2)
c         zpresm(mjbal,2) = zpresm(mjbal-1,2)
c>>?     zpresm(    1,2) = zpresm(      2,2)
         pspflx(    1)  = -pspflx(      2)
c     .................................................................c
c
CL    now compute rotational transform
CL        note: umib converts from tesla to interanal units
c
      do 36 j=3,mjbal
   36    ziotab(j,1) = r0ref * bpoli(j) / ( b0ref * avi(j,1,1))
c
cl    quadratically interpolate to xi=0 with d iota / d xi = 0 there:
c
         ziotab(2,1)= (avxib(4)**2*ziotab(3,1)-avxib(3)**2*ziotab(4,1))
     &               /(avxib(4)**2 - avxib(3)**2)
         ziotab(1,1) = ziotab(3,1)   ! interior guard point
c     ==================================================================
c
c
CLL   (4) interpolate pressure and rotational transorm
CL           from BALDUR grid to equilibrium grid
   45 continue
       call cubint (avxiz,zpresm(1,2),mjbal-1,0,zcav,Kjbal
     & ,eqxibi ,wospb ,mflxs,0,0.,1
     & ,'pressure from baldur grid to equilibrium grid')
c
       call cubint (avxiz,zpresm(1,2),mjbal-1,0,zcav,Kjbal
     & ,eqxizi ,fxp   ,mflxs,0,0.,1
     & ,'pressure from baldur grid to equilibrium grid')
c
       call cubint (avxib,ziotab(1,1),mjbal-1,0,zcav,Kjbal
     & ,eqxizi ,fxiota,mflxs,0,0.,1
     & ,'iota from baldur grid to equilibrium grid')
c     .................................................................c
      do 51 j= 1,mflxs
         fxq   (j)= 1.0/fxiota(j)
         eqprzi(j)= fxp(j)
         eqiotb(j)= fxiota(j)
   51 continue
c     .................................................................c
c
C.WOS calculate  P'  at fluxsurface grid:
c
      call cubint (avxiz,zpresm(1,2),mjbal-1,0,zcav,Kjbal
     & ,eqxizi,fxpp  ,mflxs,1,0.,+1
     & ,'p''  from baldur grid to flux-grid')
c
      call cubint (avxiz,zpresm(1,2),mjbal-1,0,zcav,Kjbal
     & ,eqxibi,wospp ,mflxs,1,0.,+1
     & ,'p''  from baldur grid to flux-grid')
c
c     .................................................................c
      if (leqtyp.lt.20)  go to 50
      call tBALFX (avxib,avxiz, zpresm(1,2),bpoli,2,mjbal-2,
     >                       zcav,Kjbal)
c     .................................................................c
c     ==================================================================
c
c
CLL   (5) compute equilibrium using profiles eqprzi and eqiotb
c
   50 continue
c
Cll   compute equilibrium using profiles eqprzi and eqiotb:
c
      call eqCTRL
c
c     ==================================================================
c
CLL   (6) Prepare arrays for sbrtns mhdstb and trpest
c
      if ( ( lsweq(4) .gt. 0  .or. lsweq(5) .gt. 0 )
     &   .and. flpols(mjbal,1) .ne. 0. ) then
c
      call mhdBAL   ! to make sure flpols(j,iz) is up to date
c
      zpsi0 = usil**2 * usib / twopi
      do 62 j=1,mjbal
      zpsi(j,1) = flpols(j,1) * zpsi0   ! poloidal flux / radian
      zpsi(j,2) = flpols(j,2) * zpsi0
      zq  (j,1) = q(j)                  ! q-value at zone bndries
  62  continue
c
      call cubint (avxiz,zpresm(1,2),mjbal,0,zcav,Kjbal
     &  ,avxib,zpresm(1,1),mjbal,0,0.,1
     &  ,'interpolate pressure from zone centers to boundaries')
c
      do 64 j=3,mjbal-1
        zqpri(j,1)  = (q(j+1) - q(j-1)) / (zpsi(j+1,1) - zpsi(j-1,1))
        zshear(j) = armins(j,1) * ( q(j+1) - q(j-1) )
     &              / ( q(j) * ( armins(j+1,1) - armins(j-1,1) ) )
        zprpri(j,1) = ( zpresm(j,2) - zpresm(j-1,2) )
     &              / ( zpsi(j,2) - zpsi(j-1,2) )
  64  continue
c
c  use symmetry conditions across the magnetic axis
c  and linear extrapolation at the edge of the plasma
c
        zqpri(1,1) = zqpri(3,1)
        zqpri(2,1) = 2 * zqpri(3,1) - zqpri(4,1)
        zqpri(mjbal,1)  = 2 * zqpri(mjbal-1,1)  - zqpri(mjbal-2,1)
        zshear(1) = - zshear(3)
        zshear(2) = 0.
        zshear(mjbal) = 2 * zshear(mjbal-1) - zshear(mjbal-2)
        zprpri(1,1) = zprpri(3,1)
        zprpri(2,1) = 2 * zprpri(3,1) - zprpri(4,1)
        zprpri(mjbal,1) = 2 * zprpri(mjbal-1,1) - zprpri(mjbal-2,1)
c
      endif
c
c     ==================================================================
c
CLL    (6) Ballooning mode stability analysis
c
c  Sbrtn mhdSTB, which computes both numerical and alnalytic ideal
c  ballooning mode stability criteria, is called whenever abs(lsweq(4))>0
c  every lsweq(4) times the equilibrium is computed
c  starting with the second time the equilibrium is computed
c  or the lsweq(6)th time the equilibrium is computed initially
c  (not before flpols(j,1) and other necessary arrays have been set
c  and the initial equilibrium iterations are completed.)
c
      isweq4 = abs(lsweq(4))
      isweq6 = max(lsweq(6),2)
cbate      isweq6 = 0  ! temporarily to agree with Stotler's code
c
      if ( isweq4 .gt. 0 .and. mod(neqinc-isweq6,max(isweq4,1)) .eq. 0
     &  .and. flpols(mjbal,1) .ne. 0. )
     &   call mhdSTB(zpsi,zq,zqpri,zshear,zprpri,zpresm,zcav,kjbal,
     &         bpoli,armins,armajs)
c
c     ==================================================================
c
CLL   (7) Transfer file for Alan Todd PEST stability analysis
c
cbate      if ( lsweq(5) .gt. 0 .and. mod(neqinc,max(lsweq(5),1)) .eq. 0
cbate     &  .and. flpols(mjbal,1) .ne. 0. )
cbate     &   call trPEST (zpsi,zq,zqpri,zprpri,zpresm,ziotab,zcav,Kjbal,
cbate     &                     neqinc,mom)
c
c     ==================================================================
c
CLL    (8) Equilibrium data needed for saturated tearing mode package
c
      do 82 j=1,mflxs
        equixi(j) = eqxibi(j)
        equir0(j) = r0hr(j)
c
        do 82 jm=1,mhrms
          eqrmom(j,jm) = rmhr(jm,j)
          eqymom(j,jm) = ymhr(jm,j)
c
  82  continue
c     _______________________________________________________________
c
c interpolate delxi2d, delth2d, ejacob2d
c   from the equilibrium grid to BALDUR zone boundaries
c
      ztemp2d = delth2d  ! copy 2-D array
c
      do jn = 1 , ntheta
        call cubint(eqxibi,ztemp2d(1,jn),mflxs,0,zcav,kjflx,
     &              avxib ,delth2d(1,jn),njav,0,0.,+1,
     &              'del xi(j,n) ---- BALDUR grids')
      enddo
c
      ztemp2d = delxi2d  ! copy 2-D array
c
      do jn = 1 , ntheta
        call cubint(eqxibi,ztemp2d(1,jn),mflxs,0,zcav,kjflx,
     &              avxib ,delxi2d(1,jn),njav,0,0.,+1,
     &              'del xi(j,n) ---- BALDUR grids')
      enddo
c
      ztemp2d = ejacob2d  ! copy 2-D array
c
      do jn = 1 , ntheta
        call cubint(eqxibi,ztemp2d(1,jn),mflxs,0,zcav,kjflx,
     &              avxib ,ejacob2d(1,jn),njav,0,0.,+1,
     &              'del xi(j,n) ---- BALDUR grids')
      enddo
c     _______________________________________________________________
c
c  compute magnetic field components
c
      call field2d
c
c     _______________________________________________________________
c
cbate      if (lsweq(1) .gt. 3) then
cbate      call timeused (icp1,ic1,isys1)
cbate      cpusec = (icp1 - icp0) * 1.e-6
cbate      write (6,*) cpusec,' = cputime (seconds) used in sbrtn mhdNEW'
cbate      write (7,*) cpusec,' = cputime (seconds) used in sbrtn mhdNEW'
cbate      endif
c
c
      return
      end
c@mhdCYL   .../baldur/code/bald/deqbald.f
c  dps 27-oct-89 15.14 Alter setting of avi(j,4,4) to allow for
c                nonuniform grid.
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE mhdCYL *
c     *********************
c
c     usage:
c           CALL mhdCYL
c
c     purpose:
c        fill in the flux surface average arrays avi(j,n,it)
c        with values appropriate for concentric circular cylinder flux
c        surfaces for an accurate comparison with the old 1-D version
c        of BALDUR
c
c***********************************************************************
c
      subroutine mhdCYL
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clparm.m'
      include 'clislb.m'
c     ------------------------------------------------------------------
c...<<<   mjbal= njav
c     ------------------------------------------------------------------
c
      write (nout  ,120)
      write (ntychl,120)
  120 format (/,' concentric circular cylinder equilibrium flux surf')
c     ------------------------------------------------------------------
c
         zavi3 = rmaji * ( twopi * rmini )**2
         zavi5 = 1. / rmini
         zavi6 = 1. / rmini**2
         zavi7 = 1. / ( rmini * rmaji )**2
         zavi8 = r0ref * b0ref
         zavi10 = 1. / rmaji**2
         zavi11 = 1. / rmaji
         zavi12 = twopi * rmaji * pi * rmini**2
         zavi13 = pi * rmini**2
c
      do 22 j=1,mjbal
         avi(j,1,4) = rmini * avxib(j)
         avi(j,2,4) = rmini
         avi(j,3,4) = zavi3 * avxib(j)
         avi(j,4,4) = zavi3 * 0.5 * (avxib(j) + avxib(j+1))
         avi(j,5,4) = zavi5
         avi(j,6,4) = zavi6
         avi(j,7,4) = zavi7
         avi(j,8,4) = zavi8
         avi(j,9,4) = zavi8
         avi(j,10,4) = zavi10
         avi(j,11,4) = zavi11
         avi(j,12,4) = zavi12 * avxib(j)**2
         avi(j,13,4) = zavi13 * avxib(j)**2
c
         avi(j,14,4) = rmaji
         avi(j,15,4) = rmini * avxib(j)
         avi(j,16,4) = 1.
         avi(j,17,4) = 0.
         avi(j,18,4) = 0.
         avi(j,19,4) = rmini * avxib(j)
         avi(j,20,4) = rmaji
   22 continue
c
c..set harmonics
c
      do 26 j=1,mflxs
        r0hr(j) = rmaji
        do 26 jm=1,mhrms
          rmhr(jm,j) = rmini * equixi(j)
          ymhr(jm,j) = rmini * equixi(j)
  26  continue
c
c..transpose to fill the eqrcjm,eqysjm,... arrays
c
      call eqHRtr
c
      return
      end
c@trPEST   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE trPEST *
c     *********************
c
c     usage:
c           call trPEST (zpsi,zq,zqpri,zprpri,zpresm,ziotab,zcav,Kjb,
c    >                        neqinc,mom)
c
c     purpose:
c        Transfer file for Alan Todd PEST stability analysis
c
c***********************************************************************
c
      subroutine trPEST (zpsi,zq,zqpri,zprpri,zpresm,ziotab,zcav,Kjb,
     >                        neqinc,mom)
c
c    dps 24-sep-86 calculate zqpri and zprpri from shear and slprs,
c                  respectively, instead of using cubint.
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c
      dimension  zcav(Kjb,3),zpresm(Kjb,2),ziotab(Kjb,2)
      dimension  zpsi(Kjb,2),zq(Kjb,2),zqpri(Kjb,2),zprpri(Kjb,2)
      dimension  zdpsidr(Kjbal)
c     ------------------------------------------------------------------
c     zcav(j,i)= temporary array for spline coefficients on BALDUR grids
c
c       In the following arrays, iz=1 refers to BALDUR zone bndries
c                                iz=2 refers to BALDUR zone centers
c
c     zpresm(j,iz) = scalar pressure on BALDUR zone bndries/centers
c                        in mks units
c
c     ziotab(j,iz) = rotational transform of BALDUR zone bndaries/cntrs
c
c     zpsi(j,iz)   = poloidal flux / ( 2 pi ) on BALDUR zone bnds/cntrs
c                  = flux stream function used in Grad-Shafranov equ.
c                    in MKS units
c
c     zq(j,iz)     = q-value on BALDUR zone bndries/centers
c
c     zqpri(j,iz)  = d q / d psi
c
c     zprpri(j,iz) = d p / d psi
c     ------------------------------------------------------------------
c...<<<   mjbal= njav
c     ------------------------------------------------------------------
c
c
c..Transfer file for Alan Todd PEST stability analysis
c
c
      call cubint(armins(1,1),zpsi(1,1),mjbal,1,zcav,Kjbal
     &   ,armins(1,1),zdpsidr,mjbal,1,0.,+1
     &   ,' dpsi/dr ')
c
      if (tbi.le.epslon+tinit) call xscale
      do 54 j=3,mjbal-1
      zden=armins(j,1)*zdpsidr(j)
      zqpri(j,1)=shear(j)*zq(j,1)/zden        !dq/dpsi from shear, dpsi/dr
c
c...Factor of 0.1 needed to convert cgs -> MKS units
c
      zprpri(j,1)=0.1*slpts(j)/zdpsidr(j)     !dp/dpsi from slpts, dpsi/dr
   54 continue
      zdelpsi=zpsi(3,1)-zpsi(2,1)
      zqpri(2,1)=(zq(3,1)-zq(2,1))/zdelpsi    !approx. for value at x=0
      zprpri(2,1)=(zpresm(3,1)-zpresm(2,1))/zdelpsi
      zqpri(1,1)=zqpri(3,1)
      zprpri(1,1)=zprpri(3,1)
      zqpri(mjbal,1)=zqpri(mjbal-1,1)
      zprpri(mjbal,1)=zprpri(mjbal-1,1)
c
c..output to transfer file 'stabil'
c
      if( lsweq(5) .gt. 0 .and. mod(neqinc,max(lsweq(5),1))
     &    .eq. 0 ) then
      write (9,150) neqinc, avti(1), nstep
 150  format (' equilibrium',i5,' at time',1pe15.7,' nstep',i5)
c
      write (9,152) label1(1:48)
      write (9,152) label2(1:48)
      write (9,152) label3(1:48)
      write (9,152) label4(1:48)
      write (9,152) label5(1:48)
 152  format (a)
c
      write (9,154) mjbal, mombnd, mom
 154  format (' mjbal=',i5,' mombnd=',i5,' mom=',i5)
c
      write (9,155) r0b
 155  format (' r0b=',1pe15.7)
c
      write (9,156) (rmb(m),m=1,mombnd)
 156  format(' rmb=',1p5e15.7)
c
      write (9,157) (ymb(m),m=1,mombnd)
 157  format (' ymb=',1p5e15.7)
c
      write (9,158) bzi, rmaji, eqcamp
 158  format (' btor=',1pe15.7,' rtor=',1pe15.7,' eqcamp=',1pe15.7)
c
      write (9,160)
 160  format (t3,'j',t10,'psi',t25,'q-value',t40,'d q / d psi'
     &  ,t55,'pressure',t70,'d p / d psi',t85,'b-pol')
c
      do 56 j=1,mjbal
      write (9,162) j,zpsi(j,1),zq(j,1),zqpri(j,1)
     &  ,zpresm(j,1),zprpri(j,1) ,bpoli(j)
  56  continue
 162  format (i5,1p6e15.7)
c
      endif
c
      return
      end
c$**********************************************************************
c
      BLOCK DATA bd_deqbald
c
      include 'come04.m'
      include 'come14.m'
      include 'come24.m'
      include 'cparm.m'
      include 'commhd.m'
      include 'clintf.m'
c
      data       ndsk/38/
c
      data         ndim/md/,ndima/na/
      data         nin/5/,nout/7/
c      data         cutoff/0.40/,psmxx/1.00/,ncl/3/,elon/1.0/
      data         cutoff/0.40/,psmxx/1.00/,ncl/3/
      data         cutfak/1.00/
c
c
      END
c
c$**********************************************************************
c@tBALFX   .../baldur/code/bald/deqbald.f
c rgb 13-jun-96 removed call link and write(nprt statements
c$**********************************************************************
c
c     *********************
c     * subroutine tBALFX *
c     *********************
c
c     purpose:
c        calculate  p'(psi), ff'(psi), q(psi)
c
c     parameters:
c           ncl    ...number of fluxsurfaces           IN
c
c***********************************************************************
c
      subroutine tBALFX (xsib,xsic,pxsic,bpolib,j1,j2,zcav,ndmj)
c
c     ------------------------------------------------------------------
      dimension  xsib(*),xsic(*),pxsic(*),bpolib(*)
      dimension  zcav(ndmj,3)
      data       it/1/,icount/-1/
      data       nprt/37/
c      data       nprt/37/,ndsk/38/
cbate      data elon/1.0/
c     ------------------------------------------------------------------
      include 'come04.m'
      include 'come14.m'
      include 'come24.m'
c     ------------------------------------------------------------------
      include 'cparm.m'
      include 'commhd.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c
cbate960613      call link ('unit38=(trdqpf,seq), ',
cbate960613     >               'unit37=(wosqpf,create,hardcopy)//')
   37 continue
c     .................................................................c
          icount= icount+1
      if (icount.le. 0)  go to 99
c     ================================================================
c
cll   1.0: calculate p'(psi):  pp(j), defined at zone-centers:
cl         ---------------------------------------------------
c
cl    1.1: calculate ppB(j), def. at zone-bondaries:
      do 11 j= j1+1,j2+1
         zdpsir= r0ref*bpolib(j)
         zdrxsi= avi(j,2,it)
         zdpxsi= (pxsic(j)-pxsic(j-1)) / (xsic(j)-xsic(j-1))
         ppb(j)= zdpxsi / (zdpsir*zdrxsi)
   11 continue
c
cl    1.2: average ppB --> ppC at zone-centers:
      do 12 j= j1,j2
         ppc(j)= 0.5*(ppb(j)+ppb(j+1))
   12 continue
         ppc(j1-1)= -ppc(j1)
         ppb(j1-1)= ppb(j1+1)
c     ------------------------------------------------------------------
c
cll   2.0: calculate ff'(psi):  ffp(j), defined at zone-centers:
cl         -----------------------------------------------------
c
      do 21 j= j1,j2
         zdaxsi= (-avi(j  ,3,it)*r0ref*bpolib(j  )*avi(j  ,7,it)
     >            +avi(j+1,3,it)*r0ref*bpolib(j+1)*avi(j+1,7,it)) /
     >                                  (xsib(j+1)-xsib(j))
         ffpc(j)= (zdaxsi/avi(j,4,it) - emu0*ppc(j)) / avi(j,10,it)
   21 continue
c     ------------------------------------------------------------------
c
cll   3.0: calculate q(psi): qvlc(j), defined at zone-centers:
cl         ---------------------------------------------------
c
cl    3.1: calculate qvlB(j), def. at zone-bondaries:
      do 31 j= j1+1,j2+1
         qvlb(j)= (b0ref*avi(j,1,it)) / (r0ref*bpolib(j))
   31 continue
         qvlb(j1-1)= qvlb(j1+1)
         qvlb(j1  )= 0.5*(qvlb(j1-1)+qvlb(j1+1))
c
cl    3.2: average qvlB --> qvlC at zone-centers:
      do 32 j= j1,j2
         qvlc(j)= 0.5*(qvlb(j)+qvlb(j+1))
   32 continue
         qvlc(j1-1)= qvlc(j1  )
c     ------------------------------------------------------------------
c
cll   4.0: transfer from BALDUR-grid to FLUXsurface-grid:
cl         ----------------------------------------------
c
         psp(1)= 0.
      do 41 j= 2,mflxs
         psp(j)= sqrt(ch(j-1))
   41 continue
c.....note: at this point psp,pspflx and eqxibi should all be the same:
c
      call CUBINT (xsic,pxsic,j2,0,   zcav,ndmj,
     >             eqxibi,fxp ,mflxs,0,0.,0,
     >             'p(tBALFX): p(Bal) --> p(flx)')
      call CUBINT (xsic(j1),ppc(j1),j2-1,0,   zcav,ndmj,
     >             eqxibi,fxpp,mflxs,0,0.,0,
     >             'p''(tBALFX): p''(Bal) --> p''(flx)')
c
      call CUBINT (xsic(j1),ffpc(j1),j2-1,0,   zcav,ndmj,
     >             eqxibi,fxffp,mflxs,0,0.,0,
     >             'ff''(tBALFX): ff''(Bal) --> ff''(flx)')
c
      call CUBINT (xsic(j1),qvlc(j1),j2-1,0,   zcav,ndmj,
     >             eqxibi,fxq  ,mflxs,0,0.,0,
     >             'qvlc(tBALFX): qvlc(Bal) --> fxq (flx)')
c     ------------------------------------------------------------------
c
c
cll   5.0: OUTPUT:
   99 continue
c
cbate960613  211 format (1x,i3,f9.3,1p,6e12.4)
cbate960613  210 format ( /,' tBALFX; icount:',i4)
cbate960613      write (nprt,210)     icount
crgb  call prline ('eqxibi$p(psi)$ p''(psi)$ff''(psi)$q(psi)$$',2,12)
      do 52 j= 1,mflxs
         psp(j)= pspflx(j)
         pp (j)= fxpp (j)
         ffp(j)= fxffp(j)
         qvl(j)= fxq  (j)
cbate960613      write (nprt,211)  j, psp (j),           pp(j),  ffp(j),qvl(j)
cbate960613      write (nout,211)  j, eqxibi(j), fxp(j),fxpp(j),fxffp(j),fxq(j)
   52 continue
c     ------------------------------------------------------------------
c
      return
      end
c@eqCTRL  .../baldur/code/bald/deqbald.f
c rgb 20.26 20-jan-92 leqtyp=11 or 12 calls eqHR12 always
c rgb 20.20 30-nov-91 temporarily removed call eqHR12, reverts to eqHR11
c rgb 20.04 26-jul-91 call eqmom2 if leqtyp=12
c  rgb 26-may-89 revised if statements and added leqtyp=8 option
c     wos 30-apr-86  use LEQTYP:  0,1,  11,  21
c     wos 02-apr-86  use LSWEQ(9):= -1,0,+1,+2;  old: lsw(5): >1,=1,<1
c     wos 01-apr-86  renamed: eqcomp->eqCTRL; different equilbr. control
c     gb  20-apr-85  interpolate flux surface averages out to mzones+1
c     gb  25-mar-85  version of sbrtn eqcomp
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqCTRL *
c     *********************
c
c     usage:
c           CALL eqCTRL
c
c     purpose:
c           subroutine eqCTRL controls 4 different types of equilibria:
c
c      (-1) circular cylinder (==> 1D BALDUR equivalent)
c      ( 0) analytical harmonic representation
c      ( 1) fixed boundary HR equilibrium  (RGB/Hirshman)
c      ( 2) free boundary (x,y) equilibrium  (wos/Jardin-code)
c
c
c           subroutine eqCTRL performs the following steps:
c
c         .1) compute an MHD equilibrium for given profiles p' and ff'
c                                              and boundary conditions
c         .2) mapping (if equilibrium in cartesian coordinates)
c         .3) compute flux surface averages of geometric quantities
c         .4) store these flux surface averages in avi(j,n,4)
c
c     ..........................................................
c     note that only the internal equilibrium common blocks
c     and the transfer common block commhd (dcommon) are present
c     ..........................................................
c
c
c***********************************************************************
c
      subroutine eqCTRL
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c
c..Various options for computing equilibrium flux surface shapes:
c
c
      if (leqtyp .ge. 20)  then
c
cll  free -boundary equilibrium (wos/Jardin-code):
cl   ---------------------------------------------
cl      2.1   calculate equilibrium:
c
        call eqxy21 (leqtyp)
c
cl      2.2   mapping and moment's representation:
c
        call eq21HR
c
c
      elseif (leqtyp .eq. 12)  then
c
cll  fixed-boundary equilibrium (Hirshman VMEC2):
cl   ----------------------------------------------
c
        call eqHR12
c
      elseif (leqtyp .eq. 11)  then
c
cll  fixed-boundary equilibrium (Hirshman VMEC):
cl   ----------------------------------------------
c
        call eqHR12
c
c
      elseif (leqtyp .eq. 8)  then
c
cll  fixed-boundary equilibrium (VMOMS package by Lao):
cl   -------------------------------------------------
c
        call eqHR08
c
c
      elseif (leqtyp .eq.  1)  then
c
cll  analytical moment's representation:
cl   -----------------------------------
c
        call eqHRA1
c
c
      else
c
cll   (-1. )  circular cylinder case; equivalent to 1D BALDUR (default)
cl            -----------------------------------------------
c
cl     -1.1   calculate equilibrium:
cl     -1.2   mapping and moment's representation:
cl     -1.3   flux surface averages:
cl     -1.4   store flux surface averages in  AVI(j,n,4):  then return
c
        call mhdCYL
        return
c
c
      endif
c
c     ==================================================================
c
c     ------------------------------------------------------------------
c
c..compute harmonics and R,Y on BALDUR grid xbouni(j)
c
      if (lsweq(20) .eq. 0) call eqryth
c
c
cll     ...   printed output of equilibria (case: 0,1,2):
c
      call eqout
c     ==================================================================
c
c
cll     x.3   flux surface averages: (for case  0,1,2):
cl      x.4   store flux surface averages in  AVI(j,n,4):
c
      call AVEQHR (nout,neqout)
c     ==================================================================
c
      return
      end
c@eqHRa1   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqHRa1 *
c     *********************
c
c     usage:
c           CALL eqHRa1
c
c     purpose:
c           subroutine eqHRa1 calculates analytical harmonics;
c
c***********************************************************************
c
      subroutine eqHRa1
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c
      write (nout,100)
 100  format (//'  Simple analytic form for shape of flux surfaces'/)
c
c.....analytic harmonics to be used for testing
c
      do 12 j=1,mflxs
      do 12 m=1,mhrms
         RmHR(m,j)= 0.0
         YmHR(m,j)= 0.0
   12 continue
c
      do 14 j=1,mflxs
         R0HR(  j)= r0b + (1.-eqxibi(j)**2) * shift
         zp = 1.
      do 14 m=1,mombnd
         zp = zp * eqxibi(j)
         RmHR(m,j)= rmb(m) * zp
         YmHR(m,j)= ymb(m) * zp
         write (nout,*) 'zp = ',zp
         write (nout,*) 'eqxibi = ',eqxibi(j)
         write (nout,*) 'rmb = ',rmb(m)
         write (nout,*) 'ymb = ',ymb(m)
   14 continue
c
c..transpose from rmhr,ymhr to eqrcjm,eqysjm
c
      call eqHRtr
c
      return
      end
c@eqHR08   .../baldur/code/bald/deqbald.f
c rgb 23-oct-89 implement second chance if vmomeq fails
c rgb 24-jul-89 used mhrms+1 rather than mhrms to set imhdeq for VMOMS
c  rgb 18-jul-89 added imhdeq and vminit to argument list of vmomeq
c  rgb 10-jul-89 added elonga to argument list of sbrtn vmomeq
c  rgb 05-jun-89 output controlled by lsweq(1)
c     print out equilibirum time (eqtime)
c  rgb 01-jun-89 created sbrtn eqHR08 to call the VMOMS package
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqHR08 *
c     *********************
c
c     usage:
c           CALL eqHR08
c
c     purpose:
c           calculate equilibrium in harmonic representation
c           using VMOMS by Lang Lao ORNL/TM-7871 (1982)
c
c***********************************************************************
c
      subroutine eqHR08
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clintf.m'
c
      dimension zhalf1(kjbal), zhalf2(kjbal), zinit0(6)
     &        , zblprs(kjbal), zeqprs(kjbal)
     &        , zblcur(kjbal), zeqcur(kjbal)
     &        , zshift(kjbal), zelong(kjbal), ztring(kjbal)
c
c  zhalf1(j) = half-width of flux surfaces on equilibrium zone bndries
c              eqxibi(j) as input for sbrtn vmomeq
c  zblprs(j) = total plasma pressure from BALDUR
c              on baldur grid avxiz(j) in nt / m**2
c  zeqprs(j) = total plasma pressure interpolated to equilibrium grid
c              eqxibi(j) in nt / m**2
c  qblcur(j) = total toroidal current within flux surface
c              on baldur grid avxib(j) in amperes
c  qeqcur(j) = total toroidal current within flux surface
c              interpolated to equilibrium grid eqxib(j) in amperes
c  zinit0(j) = initial values of vminit to be used in case vmomeq fails
c
c     Output from VMOMS code:
c
c  zhalf2(j) = half-width of flux surfaces on equilibrium zone bndries
c              eqxibi(j) as computed in sbtrn vmomeq
c  zshift(j) = Rmajor(j) - Rmajor(mflxs)
c  zelong(j) = elongation
c  ztring(j) = triangularity
c
c
      data  initin /0/
c     ------------------------------------------------------------------
c
      write (nout,100)
 100  format (//'  Equilibrium computed using VMOMS package by Lao'/)
c
c..save initial values of vminit(j), to be used in case vmomeq fails
c
      if ( initin .eq. 0 ) then
        do 6 j=1,6
          zinit0(j) = vminit(j)
   6    continue
      endif
c
c..compute profiles on equilibrium grid
c
      do 10 j=1,mflxs
        zhalf1(j) = eqxibi(j) * rminb
  10  continue
c
      zprfac = 0.1   ! used to convert pressure to nt / m**2
c
      do 12 j=1,mjbal
        zblprs(j) = totprs(j) * zprfac
        zblcur(j) = bpoli(j) * avi(j,3,1) * avi(j,2,1) * avi(j,7,1)
     &              * r0ref / ( twopi * emu0 )
  12  continue
c
c  interpolate to equilibrium grid
c
      call cubint (avxiz,zblprs,mjbal-1,0,ceqtmp,kjbal
     & ,eqxibi,zeqprs,mflxs,0,0.,1
     & ,'pressure from baldur to equilibrium grid in sbrtn eqhr08')
c
      call cubint (avxib,zblcur,mjbal,0,ceqtmp,kjbal
     & ,eqxibi,zeqcur,mflxs,0,0.,1
     & ,'current from baldur to equilibrium grid in sbrtn eqhr08')
c
c
c..calculate equilibrium:
c
      zbtorb = rbedge / rmajb  ! vacuum toroidal field at R = rmajb
      if ( elonga .lt. 0.8 ) elonga = 1. + (elongb-1.)*0.3
      imhdeq = min ( mhrms+1, 3)
      imhdeq = max ( imhdeq, 2)
      ifail  = 1
c
      call vmomeq ( initin, imhdeq, mflxs
     & ,  rmajb, zbtorb, elongb, elonga, tringb
     & , zhalf1, zeqprs, zeqcur, kjflx, lsweq(1), vminit
     & , zhalf2, zshift, zelong, ztring, eqrcj0, eqrcjm, eqysjm, ifail )
c
c..if vmomeq failed, try again with initial values of vminit
c
      if ( ifail .gt. 0 ) then
c
        do 16 j=1,6
          vminit(j) = zinit0(j)
  16    continue
c
        initin = 0
        if ( elonga .lt. 0.8 ) elonga = 1. + (elongb-1.)*0.3
        imhdeq = min ( mhrms+1, 3)
        imhdeq = max ( imhdeq, 2)
        ifail  = 1
c
        call vmomeq ( initin, imhdeq, mflxs
     & ,  rmajb, zbtorb, elongb, elonga, tringb
     & , zhalf1, zeqprs, zeqcur, kjflx, lsweq(1), vminit
     & , zhalf2, zshift, zelong, ztring, eqrcj0, eqrcjm, eqysjm, ifail )
c
      endif
c
      initin = 1
c
c..transpose harmonics
c
      do 22 j=1,mflxs
        r0hr(j) = eqrcj0(j)
        do 21 jm=1,mhrms
          rmhr(jm,j) = eqrcjm(j,jm)
          ymhr(jm,j) = eqysjm(j,jm)
  21    continue
  22  continue
c
c..output
c
      if ( lsweq(1) .gt. 5 ) then
c
        write (nout,110) eqtime
  110   format (//
     & '     Output from VMOMS equilibrium package at time '
     & ,1pe12.5,'  seconds in sbrtn eqhr08'
     & //' half-width',t18,'shift',t25,'elongation',t41,'triang'
     & ,t52,'current',t63,'pressure')
c
        do 62 j=1,mflxs
          write (nout,112) zhalf2(j), zshift(j), zelong(j), ztring(j)
     &      , zeqcur(j), zeqprs(j)
 112      format ( 0p4f12.4,1p2e15.5)
  62    continue
c
      endif
c
c..call abort if the vmomeq package has failed
c
      if ( ifail .gt. 0 ) then
        write (nout,160) ifail
 160    format (//'  Failure in VMOMS equilibrium code'
     &    ,' called from sbrtn eqhr08 (file DEQBALD)'
     &    //' IFAIL = ',i3,'  see ORNL/TM-7871 (1982) Table III'/)
        call abortb (nout,' Failure in VMOMS equilibrium code')
      endif
c
c
      return
      end
c@eqHR11  .../baldur/code/bald/deqbald.f
c  rgb 26-nov-01 use zy0 = 0.0 as first argument in call realft
c  pis 12-jun-98 changed from do 9 j=3*kmhrm
c  rgb 06-jan-92 20.25 initialize raxin = 0.
c  rgb 08-jul-90 18.45 completely new rescaling as boundary changes
c  rgb 06-jul-90 18.44 changed zftol from 1.e-5 to 1.e-14
c  rgb 28-apr-90 18.27 harmonic optimization (SCRUNCH) from Hirshman
c  rgb 07-aug-88 14.25 r0hr(0) --> r0hr(j) just before 12  continue
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqHR11 *
c     *********************
c
c     usage:
c           CALL eqHR11
c
c     purpose:
c           calculate equilibrium in harmonic representation
c           ======  NEW Bateman/Hirshman version  NEW ======
c
c***********************************************************************
c
      subroutine eqHR11
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c
      parameter  ( keqhrm=6 )
      dimension ztemp1(kmhrm), ztemp2(kmhrm), zr(knf), zy(knf)
      dimension zrbc(3*kmhrm), zybs(3*kmhrm)
      dimension zraxis(0:1), zyaxis(0:1)
c     ------------------------------------------------------------------
c
      write (nout,100)
 100  format (//'  Equilibrium computed using Hirshmans moments code'/)
c
c..harmonic optimization (SCRUNCH) from S. P. Hirshman 28-apr-90
c
      if ( lsweq(19) .gt. -1 ) then
c
c  first represent the boundary as a curve in space
c  zr(j) = R(theta(j)), zy(j) = Y(theta(j))
c
      if ( mhrms .gt. keqhrm-1 ) call abortb (6,
     & ' mhrms .gt. keqhrm-1 in sbrtn eqHR11 before sbrtn scrunch')
c
      neqtht = Knf
      call resetr (ztemp1,mhrms,0.0)
      call resetr (ztemp2,mhrms,0.0)
      zy0 = 0.0
c
      call realft (r0b,rmb,ztemp1,1,mhrms,zr,neqtht,1)
      call realft (zy0,ztemp2,ymb,1,mhrms,zy,neqtht,1)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   before call scrunch'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
 112    format (1p2e17.8)
      endif
c
      zpexp = 4.0
      zqexp = 1.0
      zftol = 1.e-14
      iphi  = 1
      ifp   = 1
      iter  = 1500
      istep = 100
c
c  note: for the moment, the maximum number of harmonics is 6
c  change parameter statement in file DSCRUNCH before changing these
c
      ipol  = mhrms
c
      do 9 j=1, 3*kmhrm ! pis 12-jun-98 changed from do 9 j=3*kmhrm
        zrbc(j) = 0.
        zybs(j) = 0.
   9  continue
c
      call scrunch (zr,zy,zpexp,zqexp,zrbc,zybs,zraxis,zyaxis
     & ,zftol,neqtht,iphi,ipol,ifp,iter,istep)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   after call scrunch'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
      endif
c
c  reset boundary harmonics
c
        r0b = zrbc(mhrms+1)
      do 10 j=1,mhrms
        rmb(j) = zrbc(j+mhrms+1)
        ymb(j) = zybs(j+mhrms+1)
  10  continue
c
c  print out optimized harmonics
c
      write (6,102) r0b
 102  format (/'  After harmonic optimization, r0b = ',0pf9.5)
      write (6,104) (rmb(j),j=1,mhrms)
 104  format (' rmb = ',0p8f9.5)
      write (6,106) (ymb(j),j=1,mhrms)
 106  format (' ymb = ',0p8f9.5)
c
      call resetr (ztemp1,mhrms,0.0)
      call resetr (ztemp2,mhrms,0.0)
      zy0 = 0.0
c
      call realft (r0b,rmb,ztemp1,1,mhrms,zr,neqtht,1)
      call realft (zy0,ztemp2,ymb,1,mhrms,zy,neqtht,1)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   with new harmonics'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
      endif
c
      endif
c
c
cl      1.2   mapping and moment's representation:
c             compute equilibrium harmonics
c
c..rescale r0hr(j), rmhr(m,j), and ymhr(m,j) to accomodate
c  possible changes in the boundary shape r0b, rmb(m), and ymb(m)
c  since the last time eqmom1 was called
c
      zrmin = 1.e-10
      if ( abs(r0hr(mflxs)) .gt. zrmin) then
        zfac0 = r0b / r0hr(mflxs)
        do 11 j=1,mflxs
          r0hr(j) = zfac0 * r0hr(j)
  11    continue
      else
        do 12 j=1,mflxs
          r0hr(j) = r0b + (1.-eqxibi(j)**2) * shift
  12    continue
      endif
c
      do 18 jm=1,mhrms
c
      if ( abs( rmhr(jm,mflxs) ) .gt. zrmin ) then
        zfacr = rmb(jm) / rmhr(jm,mflxs)
        do 13 j=1,mflxs
          rmhr(jm,j) = zfacr * rmhr(jm,j)
  13    continue
      else
        do 14 j=1,mflxs
          rmhr(jm,j) = rmb(jm) * eqxibi(j)**jm
  14    continue
      endif
c
      if ( abs( ymhr(jm,mflxs) ) .gt. zrmin ) then
        zfacy = ymb(jm) / ymhr(jm,mflxs)
        do 15 j=1,mflxs
          ymhr(jm,j) = zfacy * ymhr(jm,j)
  15    continue
      else
        do 16 j=1,mflxs
          ymhr(jm,j) = ymb(jm) * eqxibi(j)**jm
  16    continue
      endif
c
   18 continue
c
c..transpose to fill the eqrcjm,eqysjm,... arrays
c
      call eqHRtr
c
c..recompute flux surface averages
c  in particular, we need the best possible estimate of the toroidal flux
c
         lprt = lprt - 10   ! reduce printout
      call AVEflx
      call AVEpla
         lprt = lprt + 10   ! restore printout
c
         torflx = psitor (mflxs)
      write (6,*)
      write (6,*) torflx,' = toroidal flux within plasma [webers]'
c
          raxin = 0.0
c
c..compute equilibrium harmonics
c
cbate      call eqmom1 (mflxs,mhrms,r0hr,rmhr,ymhr,eqprzi,eqiotb
cbate     & ,torflx,gamma,r0b,rmb,ymb
cbate     &  ,Kmhrm,Kjflx,ftol,relerr,itmom,nskip,initin,raxin,deltin)
c
      call abortb(6
     & ,' sbrtn eqmom1 is not currently available in file deqbald.f')
c
c     ------------------------------------------------------------------
c
      initin = 1   ! so that eqmom1 does not reinitialize from scratch
c     ------------------------------------------------------------------
c
c..transfer harmonics to commhd variables
c
      call eqHRtr
c
      return
      end
c@eqHR12  .../baldur/code/bald/deqbald.f
c  rgb 26-nov-01 use zy0 = 0.0 as first argument in call realft
c rgb 20.26 14:00 21-jan-92 set raxin = r0hr(1) before calling eqmom2
c rgb 20.04 26-jul-91 subroutine eqHR12 copied from eqHR11 and changed
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqHR12 *
c     *********************
c
c     usage:
c           CALL eqHR12
c
c     purpose:
c           calculate equilibrium in harmonic representation
c           using Hirshman VMEC2 code (see file DEQMOM2)
c
c***********************************************************************
c
      subroutine eqHR12
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c
      parameter  ( keqhrm=6 )
      dimension ztemp1(kmhrm), ztemp2(kmhrm), zr(knf), zy(knf)
      dimension zrbc(3*kmhrm), zybs(3*kmhrm)
      dimension zraxis(0:1), zyaxis(0:1)
c     ------------------------------------------------------------------
c
      write (nout,100)
 100  format (//'  Equilibrium computed using Hirshmans VMEC2 code'/)
c
c..harmonic optimization (SCRUNCH) from S. P. Hirshman 28-apr-90
c
      if ( lsweq(19) .gt. -1 ) then
c
c  first represent the boundary as a curve in space
c  zr(j) = R(theta(j)), zy(j) = Y(theta(j))
c
      if ( mhrms .gt. keqhrm-1 ) call abortb (6,
     & ' mhrms .gt. keqhrm-1 in sbrtn eqHR11 before sbrtn scrunch')
c
      neqtht = Knf
      call resetr (ztemp1,mhrms,0.0)
      call resetr (ztemp2,mhrms,0.0)
      zy0 = 0.0
c
      call realft (r0b,rmb,ztemp1,1,mhrms,zr,neqtht,1)
      call realft (zy0,ztemp2,ymb,1,mhrms,zy,neqtht,1)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   before call scrunch'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
 112    format (1p2e17.8)
      endif
c
      zpexp = 4.0
      zqexp = 1.0
      zftol = 1.e-14
      iphi  = 1
      ifp   = 1
      iter  = 1500
      istep = 100
c
c  note: for the moment, the maximum number of harmonics is 6
c  change parameter statement in file DSCRUNCH before changing these
c
      ipol  = mhrms
c
      do 9 j=1,3*kmhrm
        zrbc(j) = 0.
        zybs(j) = 0.
   9  continue
c
      call scrunch (zr,zy,zpexp,zqexp,zrbc,zybs,zraxis,zyaxis
     & ,zftol,neqtht,iphi,ipol,ifp,iter,istep)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   after call scrunch'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
      endif
c
c  reset boundary harmonics
c
        r0b = zrbc(mhrms+1)
      do 10 j=1,mhrms
        rmb(j) = zrbc(j+mhrms+1)
        ymb(j) = zybs(j+mhrms+1)
  10  continue
c
c  print out optimized harmonics
c
      write (6,102) r0b
 102  format (/'  After harmonic optimization, r0b = ',0pf9.5)
      write (6,104) (rmb(j),j=1,mhrms)
 104  format (' rmb = ',0p8f9.5)
      write (6,106) (ymb(j),j=1,mhrms)
 106  format (' ymb = ',0p8f9.5)
c
      call resetr (ztemp1,mhrms,0.0)
      call resetr (ztemp2,mhrms,0.0)
      zy0 = 0.0
c
      call realft (r0b,rmb,ztemp1,1,mhrms,zr,neqtht,1)
      call realft (zy0,ztemp2,ymb,1,mhrms,zy,neqtht,1)
c
      if ( lprt .gt. 5 ) then
        write (6,*) '     zr        zy   with new harmonics'
        write (6,112) (zr(j),zy(j),j=1,neqtht)
      endif
c
      endif
c
c
cl      1.2   mapping and moment's representation:
c             compute equilibrium harmonics
c
c..rescale r0hr(j), rmhr(m,j), and ymhr(m,j) to accomodate
c  possible changes in the boundary shape r0b, rmb(m), and ymb(m)
c  since the last time eqmom1 was called
c
      zrmin = 1.e-10
      if ( abs(r0hr(mflxs)) .gt. zrmin) then
        zfac0 = r0b / r0hr(mflxs)
        do 11 j=1,mflxs
          r0hr(j) = zfac0 * r0hr(j)
  11    continue
      else
        do 12 j=1,mflxs
          r0hr(j) = r0b + (1.-eqxibi(j)**2) * shift
  12    continue
      endif
c
      do 18 jm=1,mhrms
c
      if ( abs( rmhr(jm,mflxs) ) .gt. zrmin ) then
        zfacr = rmb(jm) / rmhr(jm,mflxs)
        do 13 j=1,mflxs
          rmhr(jm,j) = zfacr * rmhr(jm,j)
  13    continue
      else
        do 14 j=1,mflxs
          rmhr(jm,j) = rmb(jm) * eqxibi(j)**jm
  14    continue
      endif
c
      if ( abs( ymhr(jm,mflxs) ) .gt. zrmin ) then
        zfacy = ymb(jm) / ymhr(jm,mflxs)
        do 15 j=1,mflxs
          ymhr(jm,j) = zfacy * ymhr(jm,j)
  15    continue
      else
        do 16 j=1,mflxs
          ymhr(jm,j) = ymb(jm) * eqxibi(j)**jm
  16    continue
      endif
c
   18 continue
c
c..compute equilibrium harmonics
c
      raxin  = r0hr(1)
      iradxi = mflxs
      iradin = 11
c
      call eqmom2 (mflxs,mhrms,r0hr,rmhr,ymhr
     & ,iradxi,eqxizi,eqprzi,eqiotb, raxin, gamma
     & ,r0b,rmb,ymb,Kmhrm,Kjflx,ftol,itmom,nskip
     & ,iradin,initin,ierflag)
c
c     ------------------------------------------------------------------
c
      initin = 1   ! so that eqmom1 does not reinitialize from scratch
c     ------------------------------------------------------------------
c
c..transfer harmonics to commhd variables
c
      call eqHRtr
c
      return
      end
c@eqxy21   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqxy21 *
c     *********************
c
c     usage:
c           CALL eqxy21 (leqtyp)
c
c     purpose:
c           calculate SCJ-equilibrium:  psi(x,y)
c
c           (1) calculate  psi(x,y)
c
c           (2) find magnetic axis and normalize psi:
c
c***********************************************************************
c
      subroutine eqxy21 (leqtyp)
c
c     ------------------------------------------------------------------
c
      if (leqtyp .eq. 20)  then
c
cll   1.0: read psi(x,y) from disk; find geom.ax. and normalize PSI:
cl         ---------------------------------------------------------
  301 format (' ### wos(0): eqxy21(0): call eqjrdd >>>>>')
  333 format (' ### wos.... dBASE7.eqxy21: DUMMY ===> STOP ###')
      write ( 7,301)
      write ( 7,333)
      call EXIT (1)
      return
      endif
c     ------------------------------------------------------------------
c
      if (leqtyp .eq. 21)  then
c
cll   2.0: read psi(x,y) from disk; find geom.ax. and normalize PSI:
cl         ---------------------------------------------------------
  321 format (' ### wos(1): eqxy21(1): call TSCequ >>>>>')
      write ( 7,321)
      write (31,321)
c
CLL   calculate psi(x,y):  ====> TSC-code (transfer psi to /come04/)
c===>>call tscequ
c
CLL   get magnetic axis and renormalize psi: (prepare mapping)
c===>>call eqgexn
      write ( 7,333)
      write (31,333)
      call EXIT (1)
      endif
c     ------------------------------------------------------------------
c
      return
      end
c@eq21HR   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eq21HR *
c     *********************
c
c     usage:
c           CALL eq21HR
c
c     purpose:
c           calculate harmonic representation of fluxsurfaces;
c           MAPPING from psi(x,y) to HR: { rmhr,ymhr }
c
c***********************************************************************
c
      subroutine eq21HR
c
      include 'cparm.m'
      include 'commhd.m'
      include 'come04.m'
      include 'come14.m'
c     ------------------------------------------------------------------
  302 format (' ### wos: eq21HR: call eqxyHR #####')
  309 format (' ### wos: eq21HR: post eqxyHR #####')
c     ------------------------------------------------------------------
cbate      data elon/1.0/
c
c
c
cll   0.   set parameters for fluxsurface representation:
cl         ----------------------------------------------
         NCL = mflxs-1
         IMAX= 25
         MMAX= mhrms
c     ------------------------------------------------------------------
c
cll   1.0: get harmonic representation:
cl         ----------------------------
      write ( 7,302)
      write (31,302)
c===>>call eqxyHR
      call EXIT (1)
      write (31,309)
      write ( 7,309)
c     ------------------------------------------------------------------
c
cll   2.0: redefine coeff. of harmonic representation:
cl         -------------------------------------------
crgb      call prline ('zH$cx$(2)$....$cy$....$$',1,10)
         zsqrt= 0.0
c
c..transpose from rmhr,ymhr to eqrcjm,eqysjm
c
      call eqHRtr
c
      return
      end
c@eqHRtr   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqHRtr *
c     *********************
c
c     Transpose from rmhr,ymhr in common /commhd/
c     to eqrcjm,eqysjm arrays in common /comequ/
c
c***********************************************************************
c
      subroutine eqHRtr
c
      include 'cparm.m'
      include 'commhd.m'
c
c
      do 22 j= 1,mflxs
         eqrcj0(j) = r0hr(j)
      do 21 m= 1,mhrms
         eqrcjm(j,m) = RmHR(m,j)
         eqysjm(j,m) = YmHR(m,j)
   21 continue
   22 continue
c
      return
      end
c@eqout   .../baldur/code/bald/deqbald.f
c rgb 24-nov-99 commented out write (8) ... unformated output to uout
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE  eqout *
c     *********************
c
c     usage:
c           CALL eqout
c
c     purpose:
c        printed output to file 'out'  bateman 25 nov 84
c        added printed output to baldur output file on channel nout
c
c***********************************************************************
c
      subroutine eqout
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c
c
      if (lprt .gt. 7) then
         call prtrv (nout  ,rmhr,1,Kmhrm,1,1,1,mhrms,mflxs,120
     &              ,'rmhr(m,j)')
         call prtrv (nout  ,ymhr,1,Kmhrm,1,1,1,mhrms,mflxs,120
     &              ,'ymhr(m,j)')
         call prtrv (neqout,rmhr,1,Kmhrm,1,1,1,mhrms,mflxs,120
     &              ,'rmhr(m,j)')
         call prtrv (neqout,ymhr,1,Kmhrm,1,1,1,mhrms,mflxs,120
     &              ,'ymhr(m,j)')
      endif
c
      if (lprt .gt. 4) then
         write (nout,130)
         write (neqout,130)
 130  format (t5,'j',t10,'pressure',t23,'iotabar',t35,'r0hr(j)')
c
      do 30 j=1,mflxs
         write (nout,132) j,eqprzi(j),eqiotb(j),r0hr(j)
  30     write (neqout,132) j,eqprzi(j),eqiotb(j),r0hr(j)
 132  format (i5,1p5e13.4)
      endif
c
c
      if (lprt .gt. 8) then
c..print out flux surface averages on BALDUR grid
c
      write (neqout,120)
 120  format (//' Flux surface averages on BALDUR grid')
c
      write (neqout,121)
 121  format (//' j',t6,'avxib',t18,'rho',t30,'d rho/d xi'
     & ,t42,'d V / d xi',t54,'<|del xi|>',t66,'<delxi**2>'
     & ,t78,'<dlxi/R**2>',t90,'R B-tor'/)
c
      do 22 j=1,njav
      write (neqout,122) j, avxib(j), avi(j,1,4), avi(j,2,4), avi(j,3,4)
     & , avi(j,5,4), avi(j,6,4), avi(j,7,4), avi(j,8,4)
  22  continue
c
      write (neqout,123)
 123  format (//' j',t6,'avxic',t18,'d V / d xi',t30,'R B-tor'
     & ,t42,'<1/R**2>'/)
c
      do 24 j=1,njav
      write (neqout,122) j,avxiz(j),avi(j,4,4),avi(j,9,4),avi(j,10,4)
  24  continue
c
 122  format (i3,1p8e12.4)
c
      endif        ! end of flux surface average printout
c     ..................................................................
c
c
c..unformatted output for post processing
c
cbate      if (lprt .gt. 1) then
cbate         write (8) txtequ(1)
cbate         write (8) txtequ(2)
cbate         write (8) txtequ(3)
c
cbat         write (8) mflxs,mhrms
c
cbate         write (8) (r0hr(j),j=1,mflxs)
cbate         write (8) ((rmhr(m,j),m=1,mhrms),j=1,mflxs)
cbate         write (8) ((ymhr(m,j),m=1,mhrms),j=1,mflxs)
c
cbate         write (8) (eqprzi(j),j=1,mflxs)
cbate         write (8) (eqiotb(j),j=1,mflxs)
cbate         write (8) torflx,gamma,b0ref,r0ref,ftol,relerr,itmom,nskip
cbate         write (8) (avti(i),i=1,5)
cbate         write (8) (eqerr(i),i=1,8)
cbate         write (8) (lsweq(i),i=1,8)
cbate      endif
c
      return
      end
c the following was previously in
c filem 11040 .bald86 wbald081(LIB) wseqmap(LIB) deqave       86-09-28
c filem  4245 .baldur65 wbald13(LIB):  dBASE7                 86-05-09
c       11040 .bald86 wbald071(LIB) wseqmap(LIB) dbase7       86-06-15
c     ------------------------------------------------------------------
c@AVEQHR   .../baldur/bald/deqbald.f
c  rgb 18-jul-01 remove call trapf, which is now renamed ftrap_hughes
c  dps 04-mar-87 add call trapf with bootstrap current cf. Mike Hughes.
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE AVEQHR *
c     *********************
c
c     usage:
c           CALL AVEQHR (nout,neqout)
c
c     purpose:
c           calculate flux surface averages (based on harmonic
c       representation of flx srfcs)  and store in  AVI(j,n,4)
c
c***********************************************************************
c
      subroutine AVEQHR (nout,neqout)
c
      include 'cparm.m'
      include 'commhd.m'
c     ------------------------------------------------------------------
c
      write (nout  ,120)
      write (neqout,120)
  120 format (/,' subroutine AVEQHR: flux surface averages from HRs')
c     ------------------------------------------------------------------
c
c
cll   1.   compute flux surface averages (metric quantities):
c
           call AVEflx
c
c
cll   2.   compute plasma averages:
c
           call AVEpla
c
c
cll   3.   compute global geometric quantities (from min/max points):
c
           call AVEglb
c
c
cll   4.   store flux surface averages in  AVI(j,n,4)
c
           call AVEavi
c
C.MHH  -  Calculate trapped particle fraction
c
c         call ftrap_hughes
c
      return
      end
c@eqryth   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE eqryth *
c     *********************
c
c     Compute harmonics and (R,Y) coordinates of flux surfaces
c  on the BALDUR grid xbouni(j).
c     Variables computed:
c  RC0XBI(jm,jx) = RC_0(xbouni(jx))
c  RCMXBI(jm,jx) = RC_jm(xbouni(jx))
c  RSMXBI(jm,jx) = RS_jm(xbouni(jx))
c  YCMXBI(jm,jx) = YC_jm(xbouni(jx))
c  YSMXBI(jm,jx) = YS_jm(xbouni(jx))
c  RTHXBI(jth,jxi) = R ( theta(jth), xbouni(jxi) )
c  YTHXBI(jth,jxi) = Y ( theta(jth), xbouni(jxi) )
c
c***********************************************************************
c
      subroutine eqryth
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c
      dimension   ztemp(Kjbal)
c
c     ------------------------------------------------------------------
c
c..use cubic spline interpolation to find
c
c     cr0(5*j), j=1,mflxs, = interp coef of eqrcj0(j)
c     crj(5*j,m), j=1,mflxs, m=1,mhrms = interp coef of eqrcjm
c     cyj(5*j,m), j=1,mflxs, m=1,mhrms = interp coef of eqysjm
c
c
      call cubint (eqxibi,eqrcj0,mflxs,0,cr0,Kjflx
     & ,avxib,rc0xbi,mjbal,0,0.,1
     & ,'abort after call cubint (eqxibi,r0,... in sbrtn eqryth')
c
      do 22 jx = 1,mjbal
        yc0xbi(jx) = 0.0
  22  continue
c
      do 28 jm=1,mhrms
c
      isym = -1
      if ( mod(jm,2) .eq. 0 ) isym = 1
c
      call cubint (eqxibi,eqrcjm(1,jm),mflxs,0,crj(1,jm),Kjflx
     &  ,avxib,ztemp,mjbal,0,0.,isym
     & ,'abort after call cubint (eqxibi,eqrcjm,... in sbrtn eqryth')
c
      do 24 jx=1,mjbal
        rcmxbi(jm,jx) = ztemp(jx)
        rsmxbi(jm,jx) = 0.0
  24  continue
c
      call cubint (eqxibi,eqysjm(1,jm),mflxs,0,cyj(1,jm),Kjflx
     &  ,avxib,ztemp,mjbal,0,0.,isym
     & ,'abort after call cubint (eqxibi,eqycjm,... in sbrtn eqryth')
c
      do 26 jx=1,mjbal
        ycmxbi(jm,jx) = 0.0
        ysmxbi(jm,jx) = ztemp(jx)
  26  continue
c
  28  continue
c     ------------------------------------------------------------------
c
c..compute R(theta) and Y(theta) on the BALDUR grid avxib(jx)
c
      ntheta = kfour
c
      if ( lsweq(21) .eq. 0 ) then
      do 32 jx=1,mjbal
        call hrmtht (rc0xbi(jx),rcmxbi(1,jx),rsmxbi(1,jx),mhrms
     &      ,rthxbi(1,jx),ntheta)
        call hrmtht (yc0xbi(jx),ycmxbi(1,jx),ysmxbi(1,jx),mhrms
     &      ,ythxbi(1,jx),ntheta)
  32  continue
      endif
c
      return
      end
c@hrmtht   .../baldur/code/bald/deqbald.f
c rgb 22-apr-96 changed Fourier transform routine from crfft2 to realft
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE hrmtht *
c     *********************
c
c     Purpose:  To convert from fourier harmonic to angle representation
c
c  Pfthet(jt) = pc0 + sum_{m=1}^{nharm} [ pcm(m) * cos ( m * theta(jt) )
c                                     + psm(m) * sin ( m * theta(jt) ) ]
c
c  for theta(jt) = 2 * pi * (jt-1) / (nharm-1)  jt = 1, nharm
c
c  Note:  ntheta must be a power of 2
c         ntheta .le. 64
c
c     This routine makes use of the Fast Fourier Transform routine
c  CRFFT2 found in OMNILIB at MFECC.
c
c***********************************************************************
c
      subroutine hrmtht (pc0,pcm,psm,nharm,pfthet,ntheta)
c
      parameter  (knf=64, knfh=1+knf/2, kpw=(3*knf/2)+2 )
c
      dimension  pcm(*), psm(*), pfthet(*)
      complex    zcmplx(knfh), zworkx(kpw)
c
      data  inital/0/
c
c..initialize
c
      if ( inital .ne. ntheta ) then
c
c  check to make sure ntheta is a power of 2 and .le. 64
c
        ipower = 2
        do 10 jp=1,5
          ipower = ipower * 2
          if ( ntheta .eq. ipower ) go to 12
  10    continue
        call abortb (6,
     &          'ntheta not power of 2 .le. 64 in sbrtn hrmtht')
  12    continue
c
c  initialize fast fourier transform
c
c***Comment out the calls to crfft2
c        call crfft2 (1,1,ntheta,zcmplx,zworkx,pfthet)
c
        inital = ntheta
c
      endif
c
c..load complex arrays for fast fourier transform
c
      ihalf = 1 + ntheta/2
c
      do 22 jm=1,ihalf
        zcmplx(jm) = cmplx(0.,0.)
  22  continue
c
      zcmplx(1) = cmplx (pc0,0.)
c
      do 24 jm=1,nharm
        zcmplx(jm+1) = 0.5 * cmplx (pcm(jm), - psm(jm) )
  24  continue
c
c..compute fourier transform
c
      call realft ( pc0, pcm, psm, 1, nharm, pfthet, ntheta, 1 )
c
c***Comment out the calls to crfft2
c      call crfft2 (0,1,ntheta,zcmplx,zworkx,pfthet)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@AVEflx   .../baldur/bald/deqbald.f
c rgb 17-jul-01 delth2d(1,jn) and delxi2d(1,jn) at magnetic axis set to
c   the their averages one grid space away from the axis
c rgb 12-jul-01 removed comtrp
c rgb 11-jul-01 remove raj, rdrdth, and rdydth, remove ../com/comtrp.m
c rgb 10-jul-01 compute delxi2d, delth2d, ejacob2d in common /come2d/
c   on equilibrium flux surfaces
c rgb 01-apr-90 changed computation of r(jt,jz) and y(jt,jz)
c rgb 26-nov-89 removed old (incorrect) Fourier transform calls to crfft2
c rgb 19-jul-89 compute absolute value of jacobian only
c  rgb 18-jul-89 implement calls to sbrtn realft for fourier transforms
c  rgb 29-may-87 added call trapf from Hughes
c  lpk 28-sep-86 added computation of delxie (gradient of xi)
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE AVEflx *
c     *********************
c
c     usage:
c           CALL AVEflx
c
c     purpose:
c             find flux surface averages suitable for 1-1/2 D transport
c
c             only purely geometric surface quantities are computed in
c                  this subroutine
c
c***********************************************************************
c
      subroutine AVEflx
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c
c     ------------------------------------------------------------------
      real :: ztemp1(knf), ztemp2(knf)
c     ------------------------------------------------------------------
cbate      call timeused (icp0,io0,isys0)
c     ------------------------------------------------------------------
c
c..use quasi-Hermite piecewise-cubic interpolation to find
c     cr0(5*j), j=1,mflxs, = interp coef of eqrcj0(j)
c     crj(5*j,m), j=1,mflxs, m=1,mhrms = interp coef of eqrcjm
c     cyj(5*j,m), j=1,mflxs, m=1,mhrms = interp coef of eqysjm
c
c..then use dcsevu to find the radial derivatives
c     dr0dxi(j),j=1,mflxs = d r0 / d xi
c     drmdxi(j,m),j=1,mflxs, m=1,mhrms = d rj(j,m) / d xi at xi=eqxibi(j)
c     dymdxi(j,m),j=1,mflxs, m=1,mhrms = d yj(j,m) / d xi at xi=eqxibi(j)
c
      call cubint (eqxibi,eqrcj0,mflxs,0,cr0,Kjflx
     & ,eqxibi,dr0dxi,mflxs,1,0.,1
     & ,'abort after call cubint (eqxibi,r0,... deriv in sbrtn AVEflx')
c
      do 22 m=1,mhrms
c
      lsym = -1
      if ( mod(m,2) .eq. 0 ) lsym = 1
c
      call cubint (eqxibi,eqrcjm(1,m),mflxs,0,crj(1,m),Kjflx
     &  ,eqxibi,drmdxi(1,m),mflxs,1,0.,lsym
     & ,'abort after call cubint (eqxibi,eqrcjm,... in sbrtn AVEflx')
c
      call cubint (eqxibi,eqysjm(1,m),mflxs,0,cyj(1,m),Kjflx
     &  ,eqxibi,dymdxi(1,m),mflxs,1,0.,lsym
     & ,'abort after call cubint (eqxibi,eqysjm,... in sbrtn AVEflx')
c
  22  continue

      eqrsjm = 0.
      eqycjm = 0.
c     ------------------------------------------------------------------
c
c..initialize work array before using crfft2 sbrtn from omnilib
c
         nf  = Knf
         nfh = Knfh
c
c..start working one flux surface at a time
c
      do 68 j=1,mflxs
c
c.. use sbrtn realft to compute real-valued Fourier transforms
c
c  r(r,j) = major radius as a function of theta(n), eqxibi(j)
c  y(n,j) = vertical distance from midplane (theta(n),eqxibi(j))
c
      do 50 jm=1,mhrms
        ztemp1(jm) = eqrcjm(j,jm)
        ztemp2(jm) = eqrsjm(j,jm)
  50  continue
c
      call realft (eqrcj0(j),ztemp1,ztemp2,1,mhrms ,r(1,j),nf,1)
c
      do 51 jm=1,mhrms
        ztemp1(jm) = eqycjm(j,jm)
        ztemp2(jm) = eqysjm(j,jm)
  51  continue
c
      call realft (eqycj0(j),ztemp1,ztemp2,1,mhrms,y(1,j),nf,1)
c
c  drdth(n) = d r(theta(n),eqxibi(j)) / d theta
c  dydth(n) = d y(theta(n),eqxibi(j)) / d theta
c
        ztemp0 = 0.
      do 52 jm=1,mhrms
        ztemp1(jm) =   jm * eqrsjm(j,jm) ! cos (jm*theta) harmonics
        ztemp2(jm) = - jm * eqrcjm(j,jm) ! sin (jm*theta) harmonics
  52  continue
c
      call realft (ztemp0,ztemp1,ztemp2,1,mhrms,drdth,nf,1)
c
        ztemp0 = 0.
      do 54 jm=1,mhrms
        ztemp1(jm) =   jm * eqysjm(j,jm)   ! cos (jm*theta) harmonics
        ztemp2(jm) = - jm * eqycjm(j,jm)   ! sin (jm*theta) harmonics
  54  continue
c
      call realft (ztemp0,ztemp1,ztemp2,1,mhrms,dydth,nf,1)
c
c  drdxi(n) = d r(theta(n),eqxibi(j)) / d xi
c  dydxi(n) = d y(theta(n),eqxibi(j)) / d xi
c
        ztemp0 = dr0dxi(j)
      do 56 jm=1,mhrms
        ztemp1(jm) = drmdxi(j,jm)
        ztemp2(jm) = 0.
  56  continue
c
      call realft (ztemp0,ztemp1,ztemp2,1,mhrms,drdxi,nf,1)
c
        ztemp0 = 0.
      do 58 jm=1,mhrms
        ztemp1(jm) = 0.
        ztemp2(jm) = dymdxi(j,jm)
  58  continue
c
      call realft (ztemp0,ztemp1,ztemp2,1,mhrms,dydxi,nf,1)
c
c..end of Fourier transform section
c
c..compute the jacobian as a function of angle and flux surface averages
c
c  aj(n) = absolute value of jacobian (theta(n),eqxibi(j))
c             where theta(n) = (n-1)*twopi/(nf)
c  avjac(j) = integral of jacobian over theta on surface eqxibi(j)
c  avir(j)   = flux surface average of (1./major radius)
c  avir2(j)  = flux surface average of (1./major radius)**2
c  avdxi(j)  = flux surface average of | del xi |
c  avdxi2(j) = flux surface average of (del xi)**2
c  dvdxi(j) = d vol / d xi  at zone boundaries
c
c  adhir2(j)   = flux surface average of (del xi)**2 / (major radius)**2
c  here, flux surface average (fsa) is the volume average between
c     differentially close flux surfaces
c  fsa of f(theta,xi) = integral over theta of jacobian * f / integral
c     over theta of jacobian
c
      avjac(j) = 0.
      avir(j)  = 0.
      avir2(j)  = 0.
      avdxi(j)  = 0.
      avdxi2(j) = 0.
      adhir2(j)   = 0.
      dvdxi(j) = 0.
c
      if (j .gt. 1) then
      do 40 n=1,nf
      aj(n) = abs( r(n,j) * (dydth(n)*drdxi(n) - drdth(n)*dydxi(n)) )
      avjac(j) = avjac(j) + aj(n)
      avir(j)  = avir(j)  + aj(n) / r(n,j)
      avir2(j)  = avir2(j)  + aj(n) / (r(n,j)**2)
      avdxi(j) = avdxi(j) + r(n,j) * sqrt (dydth(n)**2 + drdth(n)**2)
      avdxi2(j) = avdxi2(j) + r(n,j)**2*(dydth(n)**2+drdth(n)**2)/aj(n)
      adhir2(j)   = adhir2(j)   + (dydth(n)**2+drdth(n)**2)/aj(n)
c
C.MHH Save metric quantities for later
c
c         raj(j,n)=aj(n)
c         rdrdth(j,n)=drdth(n)
c         rdydth(j,n)=dydth(n)
c
c  store jacobian in 2-D array on the equilibrium grid
c  to be interpolated to BALDUR zone boundaries below
c
         ejacob2d(j,n) = aj(n)
c
  40  continue
c
      avir(j) = avir(j) / avjac(j)
      avir2(j) = avir2(j) / avjac(j)
      avdxi(j)  = avdxi(j) / avjac(j)
      avdxi2(j) = avdxi2(j) / avjac(j)
      adhir2(j) = adhir2(j) / avjac(j)
      avjac(j) = avjac(j) * twopi / nf
      dvdxi(j) = avjac(j) * twopi
      endif
c
c
c     compute ( del xi ) for ballooning mode analysis
c
      do jn = 1 , nf
c
        c1 = dydxi( jn )**2 + drdxi( jn )**2
        c2 = dydth( jn )**2 + drdth( jn )**2
        c3 = dydth( jn ) * drdxi( jn ) - drdth( jn ) * dydxi( jn )
c
        delth2d( j , jn ) = sqrt( c1 ) / max( abs(c3) , 1.0e-30 )
        delxi2d( j , jn ) = sqrt( c2 ) / max( abs(c3) , 1.0e-30 )
c
      enddo
c
  68  continue                    ! end of loop over flux surfaces
c
c..extrapolate to xi=0
c
      c2 = eqxibi(3)**2 / (eqxibi(3)**2 - eqxibi(2)**2)
      c3 = eqxibi(2)**2 / (eqxibi(3)**2 - eqxibi(2)**2)
      avir(1) = c2 * avir(2) - c3 * avir(3)
      avir2(1) = c2 * avir2(2) - c3 * avir2(3)
      avdxi(1) = c2 * avdxi(2) - c3 * avdxi(3)
      avdxi2(1) = c2 * avdxi2(2) - c3 * avdxi2(3)
      adhir2(1)  = c2 * adhir2(2)  - c3 * adhir2(3)
cbate     avjac(1) = c2 * avjac(2) - c3 * avjac(3)
cbate     dvdxi(1) = c2 * dvdxi(2) - c3 * dvdxi(3)
c
c  the value of delth2d(1,jn) and delxi2d(1,jn) at the magnetic axis
c  is the average value of those arrays one grid space away from the axis
c
      zth = 0.0
      zxi = 0.0
      do jn = 1, nf
        zth = zth + delth2d( 2, jn )
        zxi = zxi + delxi2d( 2, jn )
      enddo
      zth = zth / real ( nf )
      zxi = zxi / real ( nf )
c
      do jn = 1, nf
        delth2d( 1, jn ) = zth
        delxi2d( 1, jn ) = zxi
      enddo
c
c
c vol(j)   = volume (m^3) within zone boundary j
c area(j)   = area (m**2) within zone boundary j
c
      vol(1) = 0.
      call cubint (eqxibi,avjac,mflxs,0,ctemp,Kjflx
     & ,eqxibi,vol,mflxs,3,0.,-1
     & ,'abort after call dcsqdu(eqxibi,avjac,.. in sbrtn AVEflx')
c
      do 82 j=1,mflxs
  82  vol(j) = twopi * vol(j)
c
      write (neqout,*)
      write (neqout,*) vol(mflxs),' = plasma volume'
      write (nout,*)
      write (nout,*) vol(mflxs),' = plasma volume'
c
c..compute plasma area
c  area(j) = cross-sectional area within zone boundary j
c
      do 84 j=1,mflxs
  84  temp(j) = avjac(j) * avir(j)
c
      area(1) = 0.
      call cubint (eqxibi,temp,mflxs,0,ctemp,Kjflx
     & ,eqxibi,area,mflxs,3,0.,-1
     & ,'abort after call dcsqdu(eqxibi,temp,... area in sbrtn AVEflx')
c
      rminar(1) = 0.
      rmajgc(1) = eqrcj0(1)
      do 86 j=2,mflxs
      rminar(j) = sqrt ( area(j) / pi )
      rmajgc(j) = vol(j) / ( twopi * area(j) )
  86  continue
      write (neqout,*) area(mflxs),' = plasma area'
      write (nout,*) area(mflxs),' = plasma area'
c     ------------------------------------------------------------------
c
c..interpolate 2-D arrays to BALDUR zone boundaries
c
c     ------------------------------------------------------------------
c
c..output each flux surface
c
      if (lprt .gt. 8) then
c
      write (neqout,160)
      write (nout,160)
 160  format (//t4,'j',t7,'xi',t18,'dvdxi',t29,'avir',t40,'avir2'
     & ,t51,'avdxi2',t62,'adhir2',t73,'volume',t84,'area',t95,'arms[m]'
     & ,t106,'rmajgc[m]')
      do 90 j=1,mflxs
      write (neqout,162) j,eqxibi(j),dvdxi(j),avir(j),avir2(j),avdxi2(j)
     & ,adhir2(j),vol(j),area(j),rminar(j),rmajgc(j)
      write (nout,162) j,eqxibi(j),dvdxi(j),avir(j),avir2(j),avdxi2(j)
     & ,adhir2(j),vol(j),area(j),rminar(j),rmajgc(j)
  90  continue
 162  format (i4,1p10e11.3)
c
      endif ! end of printout in sbrtn AVEflx
c     ------------------------------------------------------------------
c
cbate      call timeused (icp1,io1,isys1)
cbate      cpusec = (icp1-icp0)*1.e-6
cbate      write (7,*) cpusec,' = cputime (seconds) used in sbrtn AVEflx'
c
      return
      end
c@AVEpla   .../baldur/code/bald/deqbald.f
c  rgb 12-jul-01 removed comtrp
c  les   nov-90 store torcur as curtor for C-M-G chi-e
c  rgb 21-dec-89 protect against division by raj(j,n)=0. after do 101
c  rgb 19-dec-89 include 'comtrp.m'
c     to include arrays raj, rdrdth, and rdydth
c  dps 13-oct-87 add calculation of internal inductance as diagnostic.
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE AVEpla *
c     *********************
c
c     usage:
c           CALL AVEpla
c
c     purpose:
c           find flux surface averages suitable for 1-1/2 D transport;
c           physical plasma quantities such as toroidal flux,
c                                              toroidal current,
c           and beta are computed in sbrtn AVEpla
c
c***********************************************************************
c
      subroutine AVEpla
c
      include 'cparm.m'
      include 'commhd.m'
      include 'comhd3.m'
      include 'cl1.m'
      include 'clintf.m'
c
      dimension zbsqav(Kjflx),zbsqia(Kjflx)
c
      do 10 j=1,mflxs
      pz (j) = eqprzi(j)
      eiota(j) = eqiotb(j)  ! changed below
  10  continue
c
c..establish interpolation of pressure and iota profiles
c
      pz(1) = pz(2)
      dxi = eqxibi(mflxs) - eqxibi(mflxs-1)
      do j=mflxs+1,Kjflx
        eqxibi(j)  = eqxibi(j-1)  + dxi
        eqxizi(j) = eqxizi(j-1) + dxi
        pz(j) = 0.
      enddo
c
      c2 = eqxibi(3)**2 / (eqxibi(3)**2 - eqxibi(2)**2)
      c3 = eqxibi(2)**2 / (eqxibi(3)**2 - eqxibi(2)**2)
c
c..interpolate pressure and iota-bar from zone centers to zone boundaries
c
      call cubint (eqxizi,eqiotb,mflxs,0,cpz,Kjflx+4
     &            ,eqxibi,eiota,mflxs,0,0.,1
     & ,'abort after cubint from eiota(j) to ciota(j) in sbrtn AVEpla')
c
      call cubint (eqxizi,pz,mflxs+4,0,cpz,Kjflx+4
     & ,eqxibi,pb,mflxs,0,0.,1
     & ,'abort after cubint(eqxizi,pz,... in sbrtn AVEpla')
      pb(mflxs) = 0.
      pb(1)  = c2 * pb(2) - c3 * pb(3)
c
c..integrate the flux-surface-averaged force balance eqn
c     to obtain an approximate profile for rbtorz(j)
c     given rbtorz(mflxs+1) = rbedge
c rbtorz(j) = major radius * toroidal mag field at zone centers
c rbtorb(j) = major radius * toroidal mag field at zone boundaries
c df2dxi(j) = d f**2 / d xi at zone boundaries
c dpdxi(j) = d pz / d xi at zone boundaries
c..at the same time, compute dsidxi(j) = d psitor / d xi
c
c===< emu0 = twopi * 2.e-7
c
      fac = 1. / twopi**4
      do 70 j=1,mflxs
      temp(j) = eiota(j) * dvdxi(j)**2 * avir2(j)
     &  * adhir2(j) * fac
  70  continue
c
      do 72 j=mflxs+1,mflxs+4
      temp(j) = temp(mflxs)
  72  continue
c
      call cubint (eqxibi,temp,mflxs+4,0,ctemp,Kjflx+4
     & ,eqxibi,dtemp,mflxs,1,0.,1
     & ,'abort after call dcsevu(eqxibi,temp.. in sbrtn AVEpla')
c
      call cubint (eqxizi,pz,mflxs+4,-1,cpz,Kjflx+4
     & ,eqxibi,dpdxi,mflxs,1,0.,1
     & ,'abort after call dcsevu(xiz,pz,... in sbrtn AVEpla')
c
      dpdxi(1) = 0.0
c
      rbtorz(mflxs+1) = rbedge
      do 74 j=mflxs,2,-1
      dxib = eqxizi(j+1) - eqxizi(j) ! dxi at zone bndry
      fac1 = 1. + eiota(j) * (temp(j) + dtemp(j) * dxib)
      fac2 = 1. + eiota(j) * (temp(j) - dtemp(j) * dxib)
c
      rbtorz(j) = sqrt (( fac1 * rbtorz(j+1)**2
     & + 2.*emu0*dxib*dpdxi(j) / avir2(j)) / fac2)
c
      rbtorb(j) = 0.5 * (rbtorz(j+1) + rbtorz(j))
c
      dsidxi(j) = rbtorb(j) * dvdxi(j) * avir2(j) / twopi
c
      df2dxi(j) = -2. * (rbtorb(j)**2 * eiota(j) * dtemp(j)
     & + emu0 * dpdxi(j) / avir2(j)) / (1. + eiota(j) * temp(j))
c
  74  continue
      rbtorz(1) = rbtorz(2)
      rbtorb(1) = c2 * rbtorb(2) - c3 * rbtorb(3)
c
      dsidxi(1) = rbtorb(1) * dvdxi(1) * avir2(1) / twopi
      df2dxi(1) = 0.0
c
c..temporary printout
c
      if (lprt .gt. 9) then
      write (neqout,171) (temp(j),j=1,Kjflx)
 171  format (/'temp = eiota * dvdxi**2 * avir2 * adhir2'
     &  /(1p10e11.3))
c
      write (neqout,172) (dtemp(j),j=1,Kjflx)
 172  format (/'d temp / d xi'
     &  /(1p10e11.3))
c
      write (neqout,173) (dsidxi(j),j=1,Kjflx)
 173  format (/'dsidxi'
     &  /(1p10e11.3))
      endif ! end of temp printout
c
c..integrate to compute psitor(j) as a function of eqxibi(j)
c
c  psitor(j) = troidal flux at zone boundaries
c  rho(j) = sqrt (toroidal flux / (pi * b0ref)) at zone boundaries
c  rhoz(j) = rho at zone centers
c  drhdxi(j) = d rho / d xi  evaluated at zone boundaries
c
      psitor(1) = 0.
      call cubint (eqxibi,dsidxi,mflxs,0,ctemp,Kjflx
     & ,eqxibi,psitor,mflxs,3,0.,-1
     & ,'abort after dcsqdu(eqxibi,dsidxi,... in sbrtn ave')
c
      rho(1) = 0.
      do 76 j=2,mflxs
      rho (j) = sqrt (psitor(j) / (pi * b0ref))
      rhoz(j) = (rho(j) + rho(j-1)) * 0.5
      drhdxi(j) = dsidxi(j) / (2 * pi * b0ref * rho(j))
  76  continue
      rhoz(1) = - rhoz(2)
c
c..extrapolate to ease use of interpolation routines
c
      drho = rho(mflxs) - rho(mflxs-1)
      do 78 j=mflxs+1,Kjflx
      rho(j) = rho(mflxs) + (j-mflxs) * drho
      rhoz(j) = rhoz(mflxs) + (j-mflxs) * drho
  78  continue
c
c..xispi(j) = sqrt of the normalized toroidal flux
c     xierr(j) = difference between eqxibi(j) and xipsi(j)
c     note, under perfect conditions, xipsi(j) should = eqxibi(j)
c
      xipsi(1) = 0.
      xierr(1) = 0.
      do 80 j=2,mflxs
      xipsi(j) = sqrt (psitor(j) / psitor(mflxs))
      xierr(j) = eqxibi(j) - xipsi(j)
  80  continue
c
c..compute volume average beta and plasma current
c dvdrhz(j) = d volume / d rho estimated at zone center
c
      do 88 j=1,mflxs
  88  temp(j) = avjac(j) * pb(j)
c
      beta(1) = 0.     ! at magnetic axis
      call cubint (eqxibi,temp,mflxs,0,ctemp,Kjflx
     & ,eqxibi,beta,mflxs,3,0.,-1
     & ,'abort after call dcsqdu(eqxibi,temp,.. in sbrtn ave')
c
      zvolav = 200. * emu0 * twopi
      zfac = 1. / (twopi * twopi * emu0)
      do 89 j=2,mflxs
      dvdrhz(j) = (vol(j)-vol(j-1)) / (rho(j)-rho(j-1))
      beta(j) = zvolav*beta(j)*rmajgc(j)**2 / (vol(j)*rbtorb(j)**2)
      torcur(j) = dvdxi(j)*dsidxi(j)*eiota(j)*adhir2(j)*zfac
  89  continue
c
c  les  nov-90 - store toroidal current for chi-e cmg
c
      call cubint (eqxibi,torcur,mflxs,0,ctemp,Kjflx,avxib,curtor
     &                      ,mjbal,0,0.,+1
     & ,'transfer tor currnt from equil to BALDUR zone bndries')
c
      dvdrhz(1) = dvdrhz(2)
      beta(1)   = 0.
      torcur(1) = 0.
c     ------------------------------------------------------------------
c
c  ...Compute internal inductance, l_i = volume integral of Bp**2 *
c                                        2 / (( mu0 * Ip )**2 * R )
c  ...Value is alint; used in MPRINT for output.
c
      nf=Knf
      zbsqav = 0.0
      alint  = 0.
c
c  ...Begin radial integration
c
c      do 102 j=2,mflxs
c        zbthp=eiota(j)*dsidxi(j)/twopi
c        zbsqav(j)=0.0
c
c  ...Begin theta integration
c
c        do 101 n=1,nf
c          if ( abs(raj(j,n)) .lt. 1.e-10 ) go to 105
c          zbthi=zbthp*sqrt(rdydth(j,n)**2+rdrdth(j,n)**2)/raj(j,n)
c          zbsqav(j)=zbsqav(j)+zbthi*zbthi*raj(j,n)
c  101   continue
c
c  ...Now have value of Bp**2 integrated over toroidal and poloidal
c  ...angles on this flux surface.
c
c        zbsqav(j)=zbsqav(j)*twopi**2/float(nf)
c  102 continue
c
c  ...Use CUBINT to perform radial integral.
c
      call cubint(eqxibi,zbsqav,mflxs,0,ctemp,Kjflx,eqxibi,zbsqia,mflxs,
     1            3,0.,-1,'abort from cubint call in AVEpla')
c
c  ...Normalize volume integral of Bp**2 by effective Bp at edge for cylinder.
c
      zlimin = (eqrcj0(1)*(emu0*torcur(mflxs))**2) / 2.
c
      if ( abs(zlimin) .gt. 1.e-10 ) then
        alint = zbsqia(mflxs) / zlimin
      else
        alint = 0.0
      endif
c
 105  continue
c
c     ------------------------------------------------------------------
c
c..output each flux surface
c
      if (lprt .gt. 5) then
c
      write (neqout,160)
      write (nout,160)
 160  format (//t4,'j',t7,'xi',t18,'eiota',t29,'psitor',t40,'dpdxi'
     & ,t51,'df2dxi',t62,'xierr',t73,'beta %',t84,'current'
     & ,t95,' ',t106,' ')
      do 90 j=1,mflxs
      write (neqout,162) j,eqxibi(j),eiota(j),psitor(j)
     & ,dpdxi(j),df2dxi(j),xierr(j),beta(j),torcur(j)
      write (nout,162) j,eqxibi(j),eiota(j),psitor(j)
     & ,dpdxi(j),df2dxi(j),xierr(j),beta(j),torcur(j)
  90  continue
 162  format (i4,1p10e11.3)
c
      write (neqout,164)            ! variables on zone centers
      write (nout,164)            ! variables on zone centers
 164  format (/t4,'j',t7,'xiz',t18,'pz',t29,'rbtorz'
     &  ,t40,'dvdrhz')
      do 94 j=1,mflxs+1
      write (neqout,162) j,eqxizi(j),pz(j),rbtorz(j)
     &  ,dvdrhz(j)
      write (nout,162) j,eqxizi(j),pz(j),rbtorz(j)
     &  ,dvdrhz(j)
  94  continue
c
      endif ! end of printout in sbrtn AVEpla
c     ------------------------------------------------------------------
c
cbate      call timeused (icp1,io1,isys1)
cbate      cpusec = (icp1-icp0)*1.e-6
cbate      write (7,*) cpusec,' = cputime (seconds) used in sbrtn AVEpla'
c
      return
      end
c@AVEglb   .../baldur/code/bald/deqbald.f
c rgb 02-apr-90 computed and printed out zshift(jg), jg=1,ngeom
c rgb 22-nov-89 corrected values of th used in interpolations
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE AVEglb *
c     *********************
c
c     usage:
c           CALL AVEglb
c
c     purpose:
c        compute global geometric quantities (elong.,triang., etc.)
c        from points on fluxsurfaces
c
c     rind,yind(jg) = inner point on midplane
c     rin ,yin (jg) = point with smallest major radius on flux surface
c     rtop,ytop(jg) = point with largest y-value on flux surface
c     rbot,ybot(jg) = point with smallest y-value on flux surface
c     rout,yout(jg) = point with largest major radius on flux surface
c
c     this routine uses r(n,j),y(n,j) as functions of theta(n), eqxibi(j)
c     previously computed in sbrtn AVEflx for n=1,nf, j=1,mflxs
c
c     reflection symmetry is assumed across midplane
c     xig(jg) is used as a flux surface label with mag axis at jg=2
c
c***********************************************************************
c
      subroutine AVEglb
c
      include 'cparm.m'
      include 'commhd.m'
      include 'clintf.m'
      include 'cl1.m'
c
      dimension zshift(kjflx)
c     ------------------------------------------------------------------
c
c
         nf  = Knf
         nfh = Knfh
c
      jg = 2
      jgskip = max (ngskip,1)  ! input variable to allow skipping surfaces
      j = 1        ! magnetic axis label used in sbrtn AVEflx
c
  10  continue     ! beginning of loop over j and jg
c
c  let there be increased frequency of flux surfaces considered
c        near the outer boundary of the plasma
c
      if (j + 4*jgskip .gt. mflxs) jgskip = max (jgskip-1,1)
c
      j = j + jgskip
      if (j .gt. mflxs) go to 70
      jg = jg + 1
      xig(jg) = eqxibi(j)        ! flux surface label
      rhoa(jg) = xig(jg)         ! flux surface label rho/a
c
      rind(jg) = eqrcj0(j)       ! compute rind
      ns = 1
      do 14 m=1,mhrms
      ns = - ns
      rind(jg) = rind(jg) + ns * eqrcjm(j,m)
  14  continue
      yind(jg) = 0.          ! assuming midplane symmetry
c
c..find minimum and maximum points on flux surface
c
      rmax = r(1,j)
      nrmax = 1
      ymax = y(1,j)
      nymax = 1
      ymin = y(1,j)
      nymin = 1
      rmin = r(nfh,j)
      nrmin = nfh
c
      do 20 n=2,nf
c
      if (r(n,j) .gt. rmax) then
         rmax = r(n,j)
         nrmax = n
         endif
      if (y(n,j) .gt. ymax) then
         ymax = y(n,j)
         nymax = n
         endif
      if (y(n,j) .lt. ymin) then
         ymin = y(n,j)
         nymin = n
         endif
      if (r(n,j) .lt. rmin) then
         rmin = r(n,j)
         nrmin = n
         endif
  20  continue
c
c..set rout(jg) = rmax
c
      n = nrmax
      if (n .gt. 1) then
         th = 0.5*(r(n-1,j)-r(n+1,j)) / (r(n-1,j)-2.*r(n,j)+r(n+1,j))
         rout(jg) = r(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(r(n-1,j)*(th-1.)+r(n+1,j)*(th+1.))
         yout(jg) = y(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(y(n-1,j)*(th-1.)+y(n+1,j)*(th+1.))
      else
         rout(jg) = r(1,j)
         yout(jg) = y(1,j)
      endif
c
c..set ytop(jg) = ymax
c
      n = nymax
      if (n .gt. 1) then
         th = 0.5*(y(n-1,j)-y(n+1,j)) / (y(n-1,j)-2.*y(n,j)+y(n+1,j))
         rtop(jg) = r(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(r(n-1,j)*(th-1.)+r(n+1,j)*(th+1.))
         ytop(jg) = y(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(y(n-1,j)*(th-1.)+y(n+1,j)*(th+1.))
      else
         call abortb (6,'cannot find ytop in sbrtn AVEglb')
      endif
c
c..set ybot(jg) = ymin
c
      n = nymin
      if (n .gt. 1) then
         th = 0.5*(y(n-1,j)-y(n+1,j)) / (y(n-1,j)-2.*y(n,j)+y(n+1,j))
         rbot(jg) = r(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(r(n-1,j)*(th-1.)+r(n+1,j)*(th+1.))
         ybot(jg) = y(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(y(n-1,j)*(th-1.)+y(n+1,j)*(th+1.))
      else
         call abortb (6,'cannot find ybot in sbrtn AVEglb')
      endif
c
c..set rin(jg) = rmin
c
      n = nrmin
      if (n .eq. nfh) then
         rin(jg) = r(n,j)
         yin(jg) = y(n,j)
      else
      if (n .gt. 1) then
         th = 0.5*(r(n-1,j)-r(n+1,j)) / (r(n-1,j)-2.*r(n,j)+r(n+1,j))
         rin(jg) = r(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(r(n-1,j)*(th-1.)+r(n+1,j)*(th+1.))
         yin(jg) = y(n,j)*(th+1.)*(1.-th)
     &       +0.5*th*(y(n-1,j)*(th-1.)+y(n+1,j)*(th+1.))
      else
         call abortb (6,'cannot find rin in sbrtn AVEglb')
      endif
      endif
c
c..set qrmid, qahalf, qelong, qtrian, and qdent
c
      qrmid(jg)   = 0.5 * (rout(jg) + rin(jg))
      qahalf(jg)  = 0.5 * (rout(jg) - rin(jg))
      qelong(jg)  = abs( (ytop(jg) - ybot(jg)) / (rout(jg) - rin(jg)) )
      qtrian(jg) = (qrmid(jg) - rtop(jg)) / qahalf(jg)
      qdent(jg)   = (rind(jg) - rin(jg)) / qahalf(jg)
c
      go to 10     ! loop back over jg
  70  continue     ! end of loop over flux surfaces
      ngeom = jg   ! total number of flux surfaces analysed in eqgeom
c
c..at the magnetic axis
c
      xig(2)    = 0.
      qrmid(2)   = eqrcj0(1)        ! major radius of magnetic axis
      qahalf(2)  = 0.
      qtrian(2) = 0.
      qdent(2)   = 0.
c
c..extrapolate qelongation to the magnetic axis
c
      c3 = xig(3)**2 / (xig(4)**2 - xig(3)**2)
      c4 = xig(4)**2 / (xig(4)**2 - xig(3)**2)
      qelong(2)  = c4 * qelong(3) - c3 * qelong(4)
c
c..guard point on oposite side of magnetic axis
c
      xig(1)    = - xig(3)
      rhoa(1)    = - rhoa(3)
      qahalf(1)  = - qahalf(3)
      qrmid(1)   =   qrmid(3)
      qelong(1)  =   qelong(3)
      qtrian(1) = - qtrian(3)
      qdent(1)   = - qdent(3)
c
c..compute shift
c
      do 74 jg=1,ngeom
        zshift(jg) = qrmid(jg) - qrmid(ngeom)
  74  continue
c
c..printout
c
      if (lprt .gt. 1) then
c
      write (nout,180)
 180  format (/t5,'geometric variables based on flux surface points'
     & /t4,'j',t9,'xi',t19,'ahalf',t29,'rmid',t39,'elong'
     & ,t49,'triang',t59,'indent',t69,'shift')
c
      do 80 jg=1,ngeom
      write (nout,182) jg,xig(jg),qahalf(jg),qrmid(jg),qelong(jg)
     & ,qtrian(jg),qdent(jg),zshift(jg)
  80  continue
 182  format (1x,i4,0p7f10.5)
c
      endif   !   end of printout
c
      return
      end
c@AVEavi   .../baldur/code/bald/deqbald.f
c$**********************************************************************
c
c     *********************
c     * SUBROUTINE AVEavi *
c     *********************
c
c     usage:
c           CALL AVEavi
c
c     purpose:
c           store flux surface averages in  AVI(j,n,4)
c
c
c        interpolate from fluxsurface-grid  to  BALDUR-grid:
c
c               eqxibi(     )                avxib(     )
c                    (Kjflx) --(CUBINT)--->      (Kjbal)
c               eqxizi(     )                avxiz(     )
c
c***********************************************************************
c
      subroutine AVEavi
c
      include 'cparm.m'
      include 'commhd.m'
      include 'cl1.m'
      include 'clintf.m'
c     ------------------------------------------------------------------
c...<<<  mjbal= njav
c     ------------------------------------------------------------------
c
c
cll   4.   store flux surface averages in  AVI(j,n,4)
c
      call cubint (eqxibi,rho,mflxs,0,ctemp,Kjflx,avxib,avi(1,1,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer rho from equil to BALDUR zone bndries')
c
      call cubint (eqxibi,rho,mflxs,0,ctemp,Kjflx,avxib,avi(1,2,4)
     &                      ,mjbal,1,0.,-1
     & ,'transfer d rho / d xi  from equil to BALDUR zone bndries')
c
      call cubint (eqxibi,dvdxi,mflxs,0,ctemp,Kjflx,avxib,avi(1,3,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer d V / d xi from equil to BALDUR zone bndries')
c
      call cubint (eqxibi,dvdxi,mflxs,0,ctemp,Kjflx,avxiz,avi(1,4,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer d V / d xi from equil to BALDUR zone centers')
c
      call cubint (eqxibi,avdxi,mflxs,0,ctemp,Kjflx,avxib,avi(1,5,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer < del xi > from equil to BALDUR zone bndries')
c
      call cubint (eqxibi,avdxi2,mflxs,0,ctemp,Kjflx,avxib,avi(1,6,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer <(del xi)**2> from equil to BALDUR zone bndries')
c
      call cubint (eqxibi,adhir2,mflxs,0,ctemp,Kjflx,avxib,avi(1,7,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer <(del xi / R)**2> from equil to BALDUR zone bndries')
c
      call cubint (eqxizi,rbtorz,mflxs,0,ctemp,Kjflx,avxib,avi(1,8,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer R * B-tor from equil to BALDUR zone bndries')
c
      call cubint (eqxizi,rbtorz,mflxs,0,ctemp,Kjflx,avxiz,avi(1,9,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer R * B-tor from equil to BALDUR zone centers')
c
      call cubint (eqxibi,avir2,mflxs,0,ctemp,Kjflx,avxiz,avi(1,10,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer < 1 / R**2 > from equil to BALDUR zone centers')
c
      call cubint (eqxibi,avir,mflxs,0,ctemp,Kjflx,avxib,avi(1,11,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer < 1 / R > from equil to BALDUR zone centers')
c
      call cubint (eqxibi,vol,mflxs,0,ctemp,Kjflx,avxib,avi(1,12,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer volume from equil to BALDUR zone boundaries')
c
      call cubint (eqxibi,area,mflxs,0,ctemp,Kjflx,avxib,avi(1,13,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer area from equil to BALDUR zone boundaries')
c
c..geometric variables based on minimax points on flux surfaces
c.....call eqgeom
c
      call cubint (xig,qrmid,ngeom,0,ctemp,Kjflx,avxib,avi(1,14,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer qrmid from equil to BALDUR zone boundaries')
c
      call cubint (xig,qahalf,ngeom,0,ctemp,Kjflx,avxib,avi(1,15,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer qahalf from equil to BALDUR zone boundaries')
c
      call cubint (xig,qelong,ngeom,0,ctemp,Kjflx,avxib,avi(1,16,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer qelong from equil to BALDUR zone boundaries')
c
      call cubint (xig,qtrian,ngeom,0,ctemp,Kjflx,avxib,avi(1,17,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer qtrian from equil to BALDUR zone boundaries')
c
      call cubint (xig,qdent,ngeom,0,ctemp,Kjflx,avxib,avi(1,18,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer qdent from equil to BALDUR zone boundaries')
c
      call cubint (eqxibi,rminar,mflxs,0,ctemp,Kjflx,avxib,avi(1,19,4)
     &                      ,mjbal,0,0.,-1
     & ,'transfer rminar from equil to BALDUR zone boundaries')
c
      call cubint (eqxibi,rmajgc,mflxs,0,ctemp,Kjflx,avxib,avi(1,20,4)
     &                      ,mjbal,0,0.,+1
     & ,'transfer rmajgc from equil to BALDUR zone boundaries')
c
      return
      end
c@mhdSTB   .../baldur/code/bald/deqbald.f
c  rgb 14.12 05-jul-88 Replaced BALCRIT with BALDGAM1 from M. Phillips
c            to compute ideal MHD high-n ballooning modes
c  rgb 14.10 08-apr-88 cequil(4) controls use of PCY pressure gradient
c       as initial guess for numerical ballooning mode computation
c  rgb 12.66 02-sep-87 lsweq(1) .gt. 4 for ballooning mode output
c       frequency of output controlled by lsweq(2)
c       equilibrium counter passed through lsweq(32)
c  dps 27-jan-87  rewrite dP/dr calculation for PCY criterion
c                      to go with new XSCALE
c       dps 07-jan-87  renamed pprcri -> aprcri
c       dps 24-sep-86  replace zshear calculation with shear from xscale
c       rgb 12-sep-86  Changed sign of zbb in ballooning mode criterion
c$**********************************************************************
c
c     compute critical pressure gradient based on either
c     the PCY analytical solution or the infinite n numerical
c     solution from PEST
c
c***********************************************************************
      subroutine mhdSTB(zpsi,zq,zqpri,zshear,zprpri,zpresm,zcav,kjb,
     &                  bpoli,armins,armajs)
c
c     _______________________________________________________________
      include 'cparm.m'
      include 'commhd.m'
      include 'clintf.m'
      include 'cl1.m'
c     _______________________________________________________________
c
c     _______________________________________________________________
      dimension zcav(kjb,3),zpresm(kjb,2),zpsi(kjb,2)
      dimension zq(kjb,2),zqpri(kjb,2),zprpri(kjb,2)
      dimension bpoli(kjb),armajs(kjb,2),armins(kjb,2)
c
      dimension zdpdr(kjbal), zshear(kjbal), zdpsidr(kjbal)
      dimension 
     &          zrjb(kjbal),       zyjb(kjbal),
     &          zrjq(kjflx),       zyjq(kjflx)
      dimension zrth(knf+5), zyth(knf+5), zbp(knf+5)
c     _______________________________________________________________
c
c     aprcri( jz , id ) = critical pressure gradient for ballooning stability
c                  id   = 1        PCY analytical solution
c                  id   = 2        PEST numerical infinite n solution
c
c     iblper = number of integration peroids used to compute ballooning modes
c     _______________________________________________________________
c
      iblper = 20
      epslon=1.0e-30
      uisb = 1.e4
      uisl=100.
      usil=0.01
      jzmin = 3
      jzmax = njav - 1
c
c     _______________________________________________________________
c
      call cubint(armins(1,1),zpsi(1,1),njav,1,zcav,55,
     &            armins(1,1),zdpsidr(1),njav,1,0,+1,
     &            ' dpsi / dr ')
c
      do 40 jz=1,njav
      zdpsidr(jz)= zdpsidr(jz)*uisl
      zdpdr(jz) = zprpri(jz,1)*zdpsidr(jz)
40    continue
c     _______________________________________________________________
c
c      PCY analytical solution
c
       do 60 jz = jzmin , jzmax
c
      zbs=zshear(jz)
      zbq=zq(jz,1)
      zbbaxi=avi(2,8,1)*avi(2,11,1)    ! B-tor axis, internal units
      zbbaxs=zbbaxi*uisb               ! B-tor axis, standard units
      zbrmins=ators(jz,1)*sqrt(b0ref/zbbaxi)  ! minor radius
      zbrmajs=rbtors(jz,1)/zbbaxs      ! major radius, R*B-tor/B-tor-axis
      zbe=zbrmins/zbrmajs              ! minor radius / major radius
      zbk=elong(jz,1)
      zbt=0.5*triang(jz,1)/zbe         ! Q in PCY paper
c
c  note the minor radius is sqrt ( toroidal flux / ( pi * B-tor-axis) )
c       the "triangularity" (Q in the PCY paper), remains nonzero
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
      zbb=-zbe*(1.0/(zbq*zbq)-zbba)
c
      zbc=zbs*zbs/2.0
      zalpc=1.0/epslon
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
      zpric = zalpc * zbbaxi**2 * zbe
     &        /(2.0*emu0*zq(jz,1)**2*armins(jz,1)*usil)
      aprcri(jz,1) = -zpric
c
60    continue
c     _______________________________________________________________
c
c     PEST numerical ballooning stability solver
c
      ibhalf = ntheta / 2 + 1
c
      do 100 jz = jzmin , jzmax
c
c     _______________________________________________________________
c
c      poloidal field at angle theta(jn) 
c
      do 102 jn = 1 , ntheta
        ith = jn + 2
        zbp(ith) = r0ref * bpoli(jz) * avi(jz,2,4) * delxi2d(jz,jn) /
     >          max( rthxbi(jn,jz) , epslon )
        zrth(ith) = rthxbi( jn , jz )
        zyth(ith) = ythxbi( jn , jz )
102   continue
c
      zbp(1)  = zbp(5)
      zbp(2)  = zbp(4)
      zrth(1) = zrth(5)
      zrth(2) = zrth(4)
      zyth(1) = zyth(5)
      zyth(2) = zyth(4)
c
      pguess=zprpri(jz,1)
      if ( cequil(4) .gt. epslon .and. zdpsidr(jz) .gt. epslon )
     &  pguess = cequil(4) * aprcri(jz,1) / zdpsidr(jz)
c
c  zthet0 = angle used for the initiation of integration in balgam1
c
      zthet0 = 0.
c
      aprcri(jz,2)=balgam1(pguess,
     &                     zprpri(jz,1),zq(jz,1),zqpri(jz,1),
     &                     zrth,zyth,zbp,
     &                     ibhalf,iblper,zthet0,1)
c
      aprcri(jz,2)=aprcri(jz,2) * zdpsidr(jz)
c
c
100   continue
c     _______________________________________________________________
c
        if ( lsweq(1) .gt. 4
     &     .and. mod(lsweq(32),max(lsweq(2),1)) .eq. 0 ) then
c
      write(nout,6000)
6000  format(1h1)
      write(nout,6001)
6001  format(5x,'Pressure Gradients for Ballooning Stability Analysis',
     & /)
      write(nout,6002)
6002  format(2x,'j',6x,' r',6x,'psi',8x,'pressure',8x,'q',
     &       10x,'shear',8x,'dp/dr  ',6x,8hp'(PEST),6x,7hp'(PCY)/
     &      ,2x,'-',6x,' -',6x,'---',8x,'--------',8x,'-',
     &       10x,'-----',8x,'-------',6x,8h--------,6x,7h-------
     &       )
c
      j=0
      do 600 jz=jzmin,jzmax
      j=j+1
      zar=armins(jz,1)*usil
      write(nout,6003)j,zar,zpsi(jz,1),zpresm(jz,1),zq(jz,1),
     &                zshear(jz),zdpdr(jz),aprcri(jz,2),aprcri(jz,1)
6003  format(1x,i3,1(1x,f8.4),8(1x,1pe12.5))
600   continue
c
        endif
c     _______________________________________________________________
c
      return
      end
c@balgam1   .../baldur/code/bald/deqbald.f
c  rgb  14.12  05-jul-88  replaced BALCRIT with BALGAM1
c        from Mike Phillips on 05-jul-88
c        which computes the critical pressure gradient for ideal MHD
c        infinite n ballooning modes.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Infinite n ballooning criteria.                                   c
c   Calculates critical pressure gradient for ballooning mode.        c
c   Author:   Michael W. Phillips                                     c
c   Date:     05-jul-88                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Removed cliche blgcom1 to file comblg1.m
c
      function balgam1(guess,ppr,qb,qpb,x,z,bp,nthe,nper
     &,theta0,mode)
c
c-----------------------------------------------------------------------
c
c   Calculates infinite n ballooning stability on a field line.
c   Solves for the critical pressure gradient at which
c   the ballooning mode is marginally stable.
c   Alternatively solves for growth rate.
c
c   Quantities are defined on a flux surface.
c   All units MKS.
c
c   guess    - guess of growth rate or pressure
c   ppr      - dp/dpsi
c   qb       - q on the field line
c   qpb      - dp/dpsi on the field line
c   x        - x coordinate of field line
c   z        - z coordinate of field line. z(3)=0.
c   bp       - poloidal magnetic field along field line
c   nthe     - number of theta points
c   nper     - number of poloidal periods to solve ballooning
c              equation over
c   mode     - 0 solves for growth rate
c              1 solves for critical p'
c
c-----------------------------------------------------------------------
c
      include 'comblg1.m'
      dimension x(*),z(*),bp(*)
c
c   define constants
      pi=4.0*atan(1.0)
      tpi=8.0*atan(1.0)
      u0=pi*4.0e-07
      nthe2=nthe+2
      nful=2*nthe-1
      nful2=nful+2
      dt=pi/(nthe-1)
      rdt=1.0/dt
      the0=theta0
      slen=dt*(nful-1)
c
      if(the0.eq.0.0)then
      npts=nper*(nthe-1)
      i0=1
      ihalf=1
      else
      npts=nper*(nful-1)-1
      i0=(npts+1)/2
      ihalf=0
      endif
c
c   fill in guard points
c   note: x, z, and bp start at i=3
      x(1)=x(5)
      x(2)=x(4)
      x(nthe2+1)=x(nthe2-1)
      x(nthe2+2)=x(nthe2-2)
      z(1)=-z(5)
      z(2)=-z(4)
      z(nthe2+1)=-z(nthe2-1)
      z(nthe2+2)=-z(nthe2-2)
      bp(1)=bp(5)
      bp(2)=bp(4)
      bp(nthe2+1)=bp(nthe2-1)
      bp(nthe2+2)=bp(nthe2-2)
c
c   calculate jacobian, grad(psi)**2 and cost
      do 18 i=3,nthe2
      xte=(8.*(x(i+1)-x(i-1))-x(i+2)+x(i-2))
     &/(12.*dt)
      zte=(8.*(z(i+1)-z(i-1))-z(i+2)+z(i-2))
     &/(12.*dt)
c
c   calc cost
      cost(i)=zte/sqrt(xte*xte+zte*zte)
c
c   calc jacobian
      xjac(i)=sqrt(xte*xte+zte*zte)/bp(i)
18    continue
c
c   calc g
      do 20 i=3,nthe2
      wtemp(i)=xjac(i)/(x(i)*x(i))
20    continue
      call balint(wtemp,ztemp,3,nthe2,dt)
      gb=qb*tpi/(2.0*ztemp(nthe2))
c
c   determine theta location on grid
      ithe=3
      tdiff=0.0
      if(the0.ne.0.0)then
c   map the0 into range 0 to tpi
      if(the0.gt.0.0)then
      the0=the0-int(the0/tpi)*tpi
      else
      the0=the0+(int(abs(the0)/tpi)+1.)*tpi
      endif
      the0=min(the0,tpi)
      the0=max(0.0,the0)
c
c   normalize ztemp to tpi once around
      zcon=pi/ztemp(nthe2)
      do 22 i=3,nthe2
      ztemp(i)=ztemp(i)*zcon
22    continue
      do 24 i=3,nthe2-1
      ztemp(2*nthe2-i)=tpi-ztemp(i)
24    continue
c
c   determine location of theta on grid
      ithe=3
      do 26 i=4,nful2-1
      if(the0.lt.ztemp(i))goto 28
      ithe=i
26    continue
28    continue
      tdiff=(the0-ztemp(ithe))/(ztemp(ithe+1)-ztemp(ithe))
      endif
c
c   calc bp**2 and b**2
      do 30 i=1,nthe2+2
      gpsi(i)=x(i)*bp(i)
      gpsi2(i)=x(i)*x(i)*bp(i)*bp(i)
      bp2(i)=bp(i)*bp(i)
      bsqr(i)=bp(i)*bp(i)+gb*gb/(x(i)*x(i))
30    continue
c
c   calc normal curvature
      do 40 i=3,nthe2
      xte=(8.*(x(i+1)-x(i-1))-x(i+2)+x(i-2))
     &/(12.*dt)
      xtte=(16.*(x(i+1)+x(i-1))-x(i+2)-x(i-2)-30.*x(i))
     &/(12.*dt*dt)
c
      zte=(8.*(z(i+1)-z(i-1))-z(i+2)+z(i-2))
     &/(12.*dt)
      ztte=(16.*(z(i+1)+z(i-1))-z(i+2)-z(i-2)-30.*z(i))
     &/(12.*dt*dt)
      curvn(i)=x(i)*(xtte*zte-ztte*xte)
     &/(bsqr(i)*xjac(i)**3.0)
     &-zte*gb*gb/(bsqr(i)*x(i)*x(i)*xjac(i))
      cappa(i)=(xtte*zte-ztte*xte)/(bp(i)*xjac(i))**3.0
40    continue
c
c  calc geodesic curvature
      do 50 i=3,nthe2
      bte=(8.*(bp2(i+1)-bp2(i-1))-bp2(i+2)+bp2(i-2))
     &/(12.*dt)
      xte=(8.*(x(i+1)-x(i-1))-x(i+2)+x(i-2))
     &/(12.*dt)
c
      curvs(i)=gb*(0.5*bte
     &-gb*gb*xte/x(i)**3.0)/(bsqr(i)*xjac(i))
50    continue
c
c   calculate coefficients needed for O.D.E.
c   that don't depend on galpha
      do 70 i=3,nthe2
      cq2(i)=(qpb*qpb*tpi*tpi/(slen*slen))*gpsi2(i)
     &/(bsqr(i))
      dq1(i)=-2.0*qpb*tpi*xjac(i)*curvs(i)
     &/(slen*bsqr(i))
70    continue
c
      do 72 i=3,nthe2-1
      cq2(2*nthe2-i)=(qpb*qpb*tpi*tpi/(slen*slen))*gpsi2(i)
     &/(bsqr(i))
      dq1(2*nthe2-i)=2.0*qpb*tpi*xjac(i)*curvs(i)
     &/(slen*bsqr(i))
      xjac(2*nthe2-i)=xjac(i)
72    continue
c
      cq2(nful2+1)=cq2(4)
      dq1(nful2+1)=dq1(4)
      xjac(nful2+1)=xjac(4)
c
      cq2(2)=cq2(nful2-1)
      dq1(2)=dq1(nful2-1)
      xjac(2)=xjac(nful2-1)
c
      call bgalpha(ppr,qb,qpb,x,z,bp,nthe)
c
      grate=guess
      if(mode.eq.0)then
      call baldif1(ppr,grate)
      else
      call baldif2(ppr,grate)
      endif
      balgam1=grate
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bgalpha(ppr,qb,qpb,x,z,bp,nthe)
c
c-----------------------------------------------------------------------
c
c   Calculates alpha coefficients.
c
c-----------------------------------------------------------------------
c
      include 'comblg1.m'
      dimension x(*),z(*),bp(*)
c
c   calc grad(alpha).grad(psi)/grad(psi)**2.0 + q'theta
      do 20 i=3,nthe2
      wtemp(i)=xjac(i)*(2.0*cappa(i)+2.0*cost(i)/x(i)
     &-x(i)*u0*ppr/bp(i))/(x(i)**3.0*bp(i))
20    continue
      call balint(wtemp,galpha,3,nthe2,dt)
c
      do 22 i=3,nthe2
      wtemp(i)=xjac(i)/(x(i)**4.0*bp2(i))
22    continue
      call balint(wtemp,ztemp,3,nthe2,dt)
c
      gpb=(tpi*qpb+2.0*gb*galpha(nthe2))
     &/(tpi*qb/gb+gb*gb*2.0*ztemp(nthe2))
c
      do 24 i=3,nthe2
      galpha(i)=gb*galpha(i)-gpb*gb*gb*ztemp(i)
24    continue
c
      do 26 i=3,nthe2
      wtemp(i)=xjac(i)/(x(i)*x(i))
26    continue
      call balint(wtemp,ztemp,3,nthe2,dt)
      do 28 i=3,nthe2
      galpha(i)=galpha(i)-gpb*ztemp(i)
28    continue
      do 30 i=3,nthe2-1
      galpha(2*nthe2-i)=-2.0*pi*qpb-galpha(i)
30    continue
c
      galph0=galpha(ithe)
     &+tdiff*(galpha(ithe+1)-galpha(ithe))
c
      const=pi*qpb/float(nthe-1)
      do 32 i=3,nful2
      galpha(i)=galpha(i)-galph0
     &+const*(i-3.0)
32    continue
c
c
c   calculate coefficients needed for O.D.E.
c   that depend on galpha
      do 50 i=3,nthe2
      cq0(i)=(1.0/gpsi2(i))+gpsi2(i)
     &*galpha(i)*galpha(i)/bsqr(i)
c
      cq1(i)=-2.0*qpb*tpi*gpsi2(i)
     &*galpha(i)/(slen*bsqr(i))
c
      dq0(i)=2.0*xjac(i)*(curvn(i)
     &/gpsi2(i)
     &+galpha(i)*curvs(i)/bsqr(i))
50    continue
c
      do 52 i=3,nthe2-1
      cq0(2*nthe2-i)=(1.0/gpsi2(i))+gpsi2(i)
     &*galpha(2*nthe2-i)*galpha(2*nthe2-i)/bsqr(i)
c
      cq1(2*nthe2-i)=-2.0*qpb*tpi*gpsi2(i)
     &*galpha(2*nthe2-i)/(slen*bsqr(i))
c
      dq0(2*nthe2-i)=2.0*xjac(i)*(curvn(i)
     &/gpsi2(i)
     &-galpha(2*nthe2-i)*curvs(i)/bsqr(i))
52    continue
c
      cq0(nful2+1)=cq0(4)
      cq1(nful2+1)=cq1(4)
      dq0(nful2+1)=dq0(4)
c
      cq0(2)=cq0(nful2-1)
      cq1(2)=cq1(nful2-1)
      dq0(2)=dq0(nful2-1)
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine baldif1(ppres,grate)
c
c-----------------------------------------------------------------------
c
c   Calculate finite difference matrix
c   for solving for growth rate.
c
c-----------------------------------------------------------------------
c
      include 'comblg1.m'
c
      mode=1
      tol=1.0e-7
      pp0=u0*ppres
      nfulm=nful-1
      nful1=nful+1
      num=1+npts/nfulm
      delt=dt
c   lowest range
      tbase=(1-i0)*delt+the0
      tbase=(int(abs(tbase)/tpi)+1)*tpi
c
      do 10 i=1,npts
      a(i)=0.0
      b(i)=0.0
      c(i)=0.0
10    continue
c
c   calc q
      do 20 i=1,npts
      theta=(i-i0)*delt+the0
      xloc=3.0+nfulm*
     &((theta+tbase)/tpi-int((theta+tbase)/tpi))
      ix=min0(int(xloc),nful1)
      delx=xloc-ix
      delx1=1.0-delx
      b(i)=-pp0*(delx1*(dq0(ix)+dq1(ix)*theta)
     &+delx*(dq0(ix+1)+dq1(ix+1)*theta))
20    continue
c
c   calc p
      do 30 i=0,npts+1
      theta=(i-i0)*delt+the0
      xloc=3.0+nfulm*
     &((theta+tbase)/tpi-int((theta+tbase)/tpi))
      ix=min0(int(xloc),nful1)
      delx=xloc-ix
      delx1=1.0-delx
      p(i)=-delx1*((cq0(ix)+(cq1(ix)+cq2(ix)*theta)*theta)
     &/xjac(ix))
     &-delx*((cq0(ix+1)+(cq1(ix+1)+cq2(ix+1)*theta)*theta)
     &/xjac(ix+1))
30    continue
c
c   calc matrix coefficients
      const=1.0/(delt*delt)
      do 36 i=1,npts
      a(i)=0.5*(p(i)+p(i-1))*const
      b(i)=b(i)-0.5*(p(i+1)+2.0*p(i)+p(i-1))*const
      c(i)=0.5*(p(i+1)+p(i))*const
36    continue
c
c   calc r
      do 44 i=1,npts
      theta=(i-i0)*delt+the0
      xloc=3.0+nfulm*
     &((theta+tbase)/tpi-int((theta+tbase)/tpi))
      ix=min0(int(xloc),nful1)
      delx=xloc-ix
      delx1=1.0-delx
      p(i)=delx1*((cq0(ix)+(cq1(ix)+cq2(ix)*theta)*theta)
     &*xjac(ix))
     &+delx*((cq0(ix+1)+(cq1(ix+1)+cq2(ix+1)*theta)*theta)
     &*xjac(ix+1))
44    continue
c
      do 48 i=1,npts
      a(i)=a(i)/p(i)
      b(i)=b(i)/p(i)
      c(i)=c(i)/p(i)
48    continue
c
c   set boundary conditions
      if(ihalf.eq.1)then
      c(1)=c(1)+a(1)
      a(1)=0.0
      c(npts)=0.0
      else
      a(1)=0.0
      c(npts)=0.0
      endif
c
      ierr=0
      call geteig1(a,b,c,npts,grate,mode,tol,ierr)
c
      return
      end
c
c.......................................................................
      subroutine geteig1(a,b,c,n,x,mode,tol,ierr)
c
c   finds closest eigenvalue to x
c
      dimension a(*),b(*),c(*)
c
c   check to see if mode number is good
      if((mode.gt.n).or.(mode.lt.1))then
      ierr=1
      return
      endif
c
c   initialize some stuff
      maxit=100
      nzero=0
      ierr=0
      gval=0.0
      rval=0.0
c
      xnew=x
      do 30 i=1,maxit
      call evpoly1(a,b,c,n,xnew,nzero,pval,gval,rval)
c
      if(nzero.eq.0)then
      valh=gval*gval-rval
      denom=sqrt((n-1.)*(n*valh-gval*gval))
      if(gval.lt.0.0)then
      denom=gval-denom
      else
      denom=gval+denom
      endif
      delta=float(n)/denom
      xnew=xnew-delta
      goto 32
      endif
c
      xold=xnew
      discrim=sqrt(abs(gval*gval-2.0*rval))
      x1=2.0/(-gval+discrim)
      x2=2.0/(-gval-discrim)
      xmin=min(x1,x2)
      xmax=max(x1,x2)
      if(xmin.lt.xnew)then
      xnew=xmin
      else
      xnew=xold-(xmax-xold)
      endif
      if(abs((xnew-xold)/xnew).lt.100.0*tol)go to 32
30    continue
      ierr=2
32    continue
c
c   converge to closest eigenvalue
      do 40 i=1,maxit
      call evpoly1(a,b,c,n,xnew,nzero,pval,gval,rval)
      valh=gval*gval-rval
      denom=sqrt((n-1.)*(n*valh-gval*gval))
      if(gval.lt.0.0)then
      denom=gval-denom
      else
      denom=gval+denom
      endif
      delta=float(n)/denom
      xnew=xnew-delta
      if(abs(delta/xnew).lt.tol)go to 42
40    continue
      ierr=2
42    continue
c
      x=xnew
c
      return
      end
c
c.......................................................................
      subroutine evpoly1(a,b,c,n,x,nzero,pval,gval,rval)
c
c   renormalized so that p0=1.0 to avoid overflow
c   gval and rval are not effected
c
      dimension a(*),b(*),c(*)
c
      jtest=0
c
      nzero=0
      p0=1.0
      p1=b(1)-x
      if((p2*p1.lt.0.0).or.(p1.eq.0.0))nzero=nzero+1
      pp0=0.0
      pp1=-1.0
      ppp0=0.0
      ppp1=0.0
c
      do 10 i=2,n
      p2=(b(i)-x)*p1-a(i)*c(i-1)
      pp2=(b(i)-x)*pp1-a(i)*c(i-1)*pp0-p1
      ppp2=(b(i)-x)*ppp1-a(i)*c(i-1)*ppp0-2.0*pp1
      if((p2*p1.lt.0.0).or.(p1.eq.0.0))nzero=nzero+1
      const=1.0/p1
      p1=p2*const
      pp0=pp1*const
      pp1=pp2*const
      ppp0=ppp1*const
      ppp1=ppp2*const
10    continue
c
      if(jtest.ne.0)then
      write(6,910)nzero,x,p2,pp2,ppp2
910   format(1x,i4,4(1x,1pe12.6))
      endif
c
      pval=p2
      gval=pp2/p2
      rval=ppp2/p2
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine baldif2(ppres,pguess)
c
c-----------------------------------------------------------------------
c
c   Calculate finite difference matrix
c   for solving for critical pressure gradient.
c
c-----------------------------------------------------------------------
c
      include 'comblg1.m'
c
      mode=1
      tol=1.0e-7
      pp0=u0*ppres
      ppg=u0*pguess
      nfulm=nful-1
      nful1=nful+1
      num=1+npts/nfulm
      delt=dt
      tbase=(1-i0)*delt
      tbase=(int(abs(tbase)/tpi)+1)*tpi
c
      do 10 i=1,npts
      a(i)=0.0
      b(i)=0.0
      c(i)=0.0
10    continue
c
c   calc p
      do 30 i=0,npts+1
      theta=(i-i0)*delt+the0
      theta=(i-i0)*delt+the0
      xloc=3.0+nfulm*
     &((theta+tbase)/tpi-int((theta+tbase)/tpi))
      ix=min0(int(xloc),nful1)
      delx=xloc-ix
      delx1=1.0-delx
      p(i)=-delx1*((cq0(ix)+(cq1(ix)+cq2(ix)*theta)*theta)
     &/xjac(ix))
     &-delx*((cq0(ix+1)+(cq1(ix+1)+cq2(ix+1)*theta)*theta)
     &/xjac(ix+1))
30    continue
c
c   calc matrix coefficients
      const=1.0/(delt*delt)
      do 36 i=1,npts
      a(i)=0.5*(p(i)+p(i-1))*const
      b(i)=-0.5*(p(i+1)+2.0*p(i)+p(i-1))*const
      c(i)=0.5*(p(i+1)+p(i))*const
36    continue
c
c   calc q
      do 20 i=1,npts
      theta=(i-i0)*delt+the0
      xloc=3.0+nfulm*
     &((theta+tbase)/tpi-int((theta+tbase)/tpi))
      ix=min0(int(xloc),nful1)
      delx=xloc-ix
      delx1=1.0-delx
      d(i)=delx1*(dq0(ix)+dq1(ix)*theta)
     &+delx*(dq0(ix+1)+dq1(ix+1)*theta)
20    continue
c
c   set boundary conditions
      if(ihalf.eq.1)then
      c(1)=c(1)+a(1)
      a(1)=0.0
      c(npts)=0.0
      else
      a(1)=0.0
      c(npts)=0.0
      endif
c
      ierr=0
      call geteig2(a,b,c,d,npts,ppg,3,tol,ierr)
      pguess=ppg/u0
c
      return
      end
c
c.......................................................................
      subroutine geteig2(a,b,c,d,npts,x,num,tol,ierr)
c
c   finds closest eigenvalues to x
c   returns the minimun eigenvalue
c
      dimension a(*),b(*),c(*),d(*)
      dimension ev(10)
c
c   initialize some stuff
      num=min0(num,10)
      maxit=100
      ierr=0
      gval=0.0
      rval=0.0
c
      do 20 n=1,num
      xnew=x
c   converge to closest eigenvalue
      do 30 i=1,maxit
      call evpoly2(a,b,c,d,npts,xnew,pval,gval,rval)
      valh=gval*gval-rval
c
c   subtract off any modes already found
      do 32 k=1,n-1
      gval=gval-1.0/(xnew-ev(k))
      valh=valh-1.0/((xnew-ev(k))*(xnew-ev(k)))
32    continue
c
      denom=sqrt(max(0.0,(npts-1.)*(npts*valh-gval*gval)))
      denom=gval+sign(denom,gval)
      delta=float(npts)/denom
      xnew=xnew-delta
      if(abs(delta/xnew).lt.tol)go to 38
30    continue
      ierr=2
38    continue
c
      ev(n)=xnew
20    continue
c
c   find maximum negative root or minimum positive root
      x=ev(1)
      do 40 n=2,num
      if(ev(n).lt.x)x=ev(n)
40    continue
      do 42 n=1,num
      if((ev(n).lt.0.0).and.(ev(n).gt.x))x=ev(n)
42    continue
c
      return
      end
c
c.......................................................................
      subroutine evpoly2(a,b,c,d,npts,x,pval,gval,rval)
c
c   renormalized so that p0=1.0 to avoid overflow
c   gval and rval are not effected
c
      dimension a(*),b(*),c(*),d(*)
c
      p0=1.0
      p1=b(1)-x*d(1)
      pp0=0.0
      pp1=-d(1)
      ppp0=0.0
      ppp1=0.0
c
      do 10 i=2,npts
      p2=(b(i)-x*d(i))*p1-a(i)*c(i-1)
      pp2=(b(i)-x*d(i))*pp1-a(i)*c(i-1)*pp0-d(i)*p1
      ppp2=(b(i)-x*d(i))*ppp1-a(i)*c(i-1)*ppp0-2.0*d(i)*pp1
      const=1.0/p1
      p1=p2*const
      pp0=pp1*const
      pp1=pp2*const
      ppp0=ppp1*const
      ppp1=ppp2*const
10    continue
c
      pval=p2
      gval=pp2/p2
      rval=ppp2/p2
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine balint(wtemp,ztemp,nbgn,nend,dt)
c
c----------------------------------------------------------------------
c
c     used by balcrit to calculate integrals
c
c----------------------------------------------------------------------
c
      dimension wtemp(nend),ztemp(nend)
c
      ztemp(nbgn)=0.0
      do 10 i=nbgn+1,nend
      ztemp(i)=0.5*(wtemp(i)+wtemp(i-1))
10    continue
c
      do 12 i=nbgn+2,nend
      ztemp(i)=ztemp(i)+ztemp(i-1)
12    continue
c
      do 20 i=nbgn,nend
      ztemp(i)=ztemp(i)*dt
20    continue
c
      return
      end
c
c----------------------------------------------------------------------
c$**********************************************************************
c@diary   .../baldur/code/bald/deqbald.f
c  rgb 17-jan-88 changed the sign of the imaginary part of zcmplx
c       in sbrtn hrmtht
c       in order to make y(theta,xi) increase with theta from zero
c  rgb 16-aug-88 14.26 replace rthj with rthxbi, ythj with ythxbi,
c        and iblthe with ntheta
c  rgb 16-aug-88 using lsweq(20) = 0 and lsweq(21) = 0 (defaults)
c        to control calls to sbrtns eqryth and hrmtht.
c  rgb 08-aug-88 14.25 changed ntheta from Knfh to Kfour in sbrtn eqhrth
c  rgb 07-aug-88 14.25 call eqHRtr after each time r0hr,rmhr,ymhr is set
c  rgb 07-aug-88 14.25 r0hr(0) --> r0hr(j) in sbrtn eqhr11
c  rgb 25-jul-88 14.22 include clintf in sbrtn aveglb
c      so that eqxibi is available
c  rgb 22-jul-88 14.21 rmhr,ymhr are now the basic harmonic representn
c       sbrtn eHRtr is used to transpose to eqrcjm,eqysjm arrays
c  rgb 21-jul-88 14.20 cleaned up sbrtn headers
c       imlemented leqxi to control spacing of eqxibi(j) and eqxizi(j)
c       improved documentation of equilbirum input variables
c  rgb 20-jul-88 replaced jr and njeq with mflxs
c  rgb 20-jul-88 removed all references to xi(j) and xiz(j)
c       replaced xi(j) with eqxibi(j),  rplaced xiz(j) with eqxizi(j)
c       replaced fxsib with eqxibi,  replaced fxsiz with eqxizi
c  rgb 31-may-87 added shift to sbrtn eqhra1 for analytic flux surfaces
c                called eqIN11 from eqin01 to initialize analytic surfaces
c  rgb 30-apr-87 Set ELLIPT = AVI(MZONES,16,1) in sbrtn MHDBAL
c     GB 20-jun-86  Moved "estabilish equilibrium grid..."
c from sbrtn MHDNEW to sbrtn EQINIT in file deqinit.
c$**********************************************************************
c     GB 19-jun-86  Moved definition of eqxibi, eqxizi
c  to sbrtn EQINIT where it is defined once at the beginning of each run.
c     GB 20-jun-86  Added cequil and eqerr to namelist input.
c  cequil(10) = initial central iota-bar when lsweq(6) .gt. 0
c     GB 21-jun-86
c
c  renamed lsweq( 9) --> lsweq(31)  no longer used, remove later
c
c  lsweq(10) will now be used to control initialization of the
c        poloidal magnetic field bpoli on the first pass.
c     lsweq(10) = 0  (default) starts with bpoli as it was set in
c        sbrtn auxval of the original 1-D BALDUR code
c     lsweq(10) = 1  uses quadratic interpolation of q(xi)
c        between q0 = qmin and qedge = 3. * qmin
c     lsweq(10) = 2  uses cubic interpolation of iota between
c        iota edge as determined by toroidal current after
c             computing V'(xi) and <|del xi|**2/R**2> in sbrtn aveflx
c        iota axis = max ( cequil(10) set by equilibrium namelist,
c             minimum iota consistent with d iota / d xi .le. 0.

cc#! /bin/csh  deqvmom
c#/ 22:00 15-feb-92 /11040/.../baldur/code/equil/deqvmom  Bateman, PPPL
c#/ 23:00 15-aug-91 /11040/bald90/wbaldn1 wsequil DEQVMOM, Bateman, PPPL
c#/ unix version
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c**********************************************************************c
c
c    To obtain this file type:
c cfs get /11040/bald89/wbaldn1 ^ end       (note ^ stands for linefeed)
c lib wbaldn1 ^ x wsequil ^ n wsequil ^ x deqvmom ^ end
c lib wbaldn2 ^ x yequilib yutilib ^ end  ( libraries of compiled sbrtns)
c lib wbaldn3 ^ x xstripx ^ end           ( utility to strip !...)
c
c    To compile this file type:
c csh deqvmom
c
c**********************************************************************c
c@changes  11040/baldur/code/bald/deqvmom
c rgb 23-apr-96 removed function spmpar, not being used
c rgb 20.26 13:00 19-jan-92 changed names: machin to machinb, 
c   d02agE to d02agEb, agzd02 to agzd02b, f03afE to f03afEb, 
c   f04ajE to f04ajEb, p01aaE to p01aaEb, i1mach to i1machb,
c   agyd02 to agyd02b, x03aaE to x03aaEb in order to avoid name
c   conflicts with the correspoding NAG routines they replace.
c   Also changed all hollerith to '...' fields in format statements.
c rgb 20.26 14:00 17-jan-92 moved routines machinb, d02agEb, agzd02b,
c   f03afEb, f04ajEb, p01aaEb, i1machb, agyd02b, and x03aaEb
c   from teqvmom.f  to teqvmain.f
c rgb 18.81 16-dec-90 changed double precision srname to character *8 srname
c     removed multiple declarations of initin and mm1t
c rgb 06-aug-90 failure flag can be set in sbrtn vmaux (ifail=9)
c rgb 29-mar-90 diagnistic printout and loop over input in main
c rgb 20-nov-89 ymb(j,2) = - zeout * rmb(j,2) in sbrtn vmomeq
c  Note:  In the VMOMS code, the angle theta is measured from the inboard
c    edge of the toroidal plasma:  hence R_1 = - rho.
c    In the rest of the BALDUR code, theta is measured from the outboard
c    edge of the plasma:  Hence R_1 = + rho and Y_1 = - elong * R_2.
c rgb 20-nov-89 changed rmb(j,2)=zdout*0.01 to =zdout*0.01*zrho(nrho)
c rgb 18-jul-89 gl(1-4) restored to original in sbrtn vmbcux
c rgb 12-jul-89 correction using qstar to determine zcrtot in main
c rgb 10-jul-89 implemented lprofl=2 profile option
c   Problem with convergence.  Change initialization of param(3) 
c      to a lower value in sbrtn vmhdeq.
c      param(3) is now controlled by input elongation at axis.
c   Implemented nmhdeq = number of equations to be solved
c   Implemented xpinit(j) to input initial values of param from namenlist
c rgb 02-may-89 tightened error bounds:
c      erry from 1.e-2 to 1.e04, errp from 1.e-3 to 1.e-5
c   Note: tightening the error bounds had the effect of bringing the
c         results into much closer agreement with ORNL/TM-7871 (1982)
c      fixed printout of shift, elong, and triang.
c rgb 25-apr-89 test code using parameters from ORNL/TM-7871 (1982)
c rgb 24-apr-89 LPRINT controls amount of printout
c     xe0,... geometric and e0,... moments all computed in sbrtn vmhdeq
c rgb 17-apr-89 corrections from original form of VMOMS program
c rgb 13-apr-89 set nrho,rho,e0,d0 in sbrtn vmhdeq
c     avoided division by pres(1) in sbrtn vmpfun
c     rearranged order of subrtns and added list of sbrtns
c     set xl0 = xl00, which is set by data statement in sbrtn vmhdeq
c     initialize nfl,...,parerr in sbrtn vmhdeq
c     remove variables xa
c rgb 12-apr-89 combine input file into t0pre yielding binary b0
c     remove libraries yequilib and quikdraw
c     Added iptype to argument list of call vmhdeq
c rgb 11-apr-89 implemented implicit none in main and vmomeq
c rgb 31-mar-89 added input/output unit numbers in sbrtn vmhdeq
c**********************************************************************c
c@latex documentation  11040/bald89/weqvmom DEQVMOM
c
c    To extract the LaTex comments between each occurance of 'cbtex'
c    and 'cetex' throughout the code, on the VAX type
c @extex deqvmom. deqvmom.tex
c gtex deqvmom
c
cbtex
c % Latex documentation of file DEQVMOM
c
c \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
c \newcommand{\volave}[1]{\left\langle #1 \right\rangle}
c \newcommand{\leapprox}{\stackrel{<}{\sim}}
c \newcommand{\geapprox}{\stackrel{>}{\sim}}
c
c \newcommand{\beq}{\begin{equation}}
c \newcommand{\beql}[1]{\begin{equation} \label{#1}} 0.000000or labelled equations
c \newcommand{\eeq}{\end{equation}}
c
c \documentstyle[12pt,aip]{article}
c \let\footnotesize=\small
c
c \pagestyle{myheadings}
c \markright{\hfil deqvmom. \hfil Draft by Glenn Bateman \hfil \today}
c
c \begin{document}
c
c \vspace{1 cm} \begin{center} \large \bf{
c DEQVMOM: Driver Program for the \\
c VMOMS Equilibrium Computer Code \\} 0.000000e+00nd of bold face
c \medskip \normalsize
c by L.L. Lao, R.M. Wieland, W.A. Houlberg, and S.P. Hirshman, \\
c ORNL/TM-7871 (Feb, 1982) \\
c  ``VMOMS: A Computer Code for Finding \\
c              Moment Solutions to the Grad-Shafranov Equation'' \\
c \medskip
c DEQVMOM Interface and Documentation \\
c by Glenn Bateman \\
c Princeton Plasma physics Laboratory \\ \medskip \today \end{center}
c
c
c    This file contains a driver program for the equilibrium code VMOMS
c by L.L. Lao, R.M. Wieland, W.A. Houlberg, and S.P. Hirshman,
c ORNL/TM-7871 (Feb, 1982)  ``VMOMS: A Computer Code for Finding
c                  Moment Solutions to the Grad-Shafranov Equation''.
c
c \section{Harmonic Rrepresentation}
c
c     The harmonic representation used in the VMOMS code is written
c \beql{R(rho,theta)} R(\rho,\theta) = R_0(\rho) - \rho \cos(\theta)
c   + R_2(\rho) \cos ( 2 \theta)  \eeq
c \beql{Y(rho,theta)} Y(\rho,\theta) = E(\rho) [ \rho \sin(\theta)
c   + R_2(\rho) \sin ( 2 \theta )  \eeq
c where $\theta$ is measured from the inboard edge of the toroidal plasma.
c  Then the major radius to the midpoint of each flux surface on the 
c midplane is
c \beql{Rmajor} R_{major} = \frac{1}{2} ( R(\rho,0) + R(\rho,\pi) )
c   = R_0(\rho) + R_2(\rho) \eeq
c and the half-width is $\rho$.
c
c  Let $\theta_c$ be the value of $\theta$ at which $Y(\rho,\theta)$
c is maximum --- that is 
c \beq \rho \cos(\theta_c) + 2 R_2(\rho) \cos (2\theta_c) = 0. \eeq
c  Then the geometric elongation is defined as the height to this point
c divided by the half-width
c \beq \kappa = Y(\rho,\theta_c) / \rho, \eeq
c (which is called {\tt elong} in the code.
c  The geometric triangularity is defined as the diference between the
c midpoint major radius and the major radius to this point, all divided
c by the half-width
c \beq \delta = ( R_{major} - R(\rho,\theta_c) ) / \rho \eeq
c (which is called {\tt triang} in the code).
c  The shift is defined as the diference between
c the major radius of each flux surface and the major radius at the edge
c of the plasma (of the computational domain being considered)
c \beq {\tt shift} = R_{major}(\rho) - R_{major}(a)
c   = R_0(\rho)+R_2(\rho)-R_0(\rho)-R_2(\rho) \eeq
c
c  The output harmonic representation is written
c \beq R(\rho,\theta) = R_0(\rho) + R_1(\rho) \cos(\theta)
c    + R_2(\rho) \cos(2 \theta) \eeq
c \beq Y(\rho,\theta) = Y_1(\rho) \sin(\theta)
c    + Y_2(\rho) \sin(2 \theta) \eeq
c
cetex
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c     Subrtns in DEQVMOM:
c     -------------------
c
c MAIN    driver program to test sbrtn vmomeq
c EQVMOM  used to isolate VMOMS package from BALDUR
c VMCOM0  parameter cliche
c VMCOM1  common block cliche
c VMCOM2  common block (counter, io units, machine epsilon) cliche
c VMHDEQ  entry sbrtn for the VMOMS package
c VMGTRN  transforms from geometric to moments representation (GEOTRN)
c VMNROT  
c VMAUX   source functions for ODE's (AUX)
c VMBCUX  boundary conditions for equilibrium moment equations (BCAUX)
c VMRAUX  sets the matching point for the mhd moments (RAAUX)
c VMPRSL  print out iteration parameters for the mhd moments (PRSOL)
c VMAIFN  integrated toroidal current profile (AIFUNR)
c VMPFUN  evaluate plasma pressure and its derivative (PFUN)
c VMTRIG  sets up the sine and cosine arrays
c VMSETX  compute xp array (?)
c VMWGTS  set up integration weighting array (WGTS)
c VMWGTT  set up trapezoidal rule weighting function (WGTT)
c VMWROT  (WROTFN)
c D02AGE  ODE solver from NAG mark 4.5
c AGZD02
c AGYD02  integrates system of ODE's using Merson's method
c F03AFE  for matrix decomposition from NAG mark 2
c F04AJE  linear equation solver from NAG mark 2
c P01AAE  error message from NAG routines from NAG mark 7
c X03AAE  value of a scalar product from NAG mark 4.5
c AAZX03  used in x03aaEb
c SDOT    vector dot product (OMNILIB implementation used instead)
c SPMPAR  machine parameters
c I1MACH  machine parameters
c LQDCMP  Gaussian elimination matrix decomposition (Forsythe,...)
c         (used to be called sbtn DECOMP)
c LQSOLV  linear equation solver from Forsythe, Malcolm and Moler
c         (used to be called sbrtn SOLVE)
c MACHIN  machine epsilon
c
c     Extra sbrtns no longer used:
c
c GPASMA  compute integrated quantities (from TRANSP implementation)
c PLTMEQ  generate output for plots (from TRANSP implementation)
c PRTMEQ  generate printout (from TRANSP implementation)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    This file contains the equilibrium code VMOMS
c by L.L. Lao, R.M. Wieland, W.A. Houlberg, and S.P. Hirshman
c see ORNL/TM-7871 (Feb, 1982)  VMOMS: A Computer Code for Finding
c                  Moment Solutions to the Grad-Shafranov Equation
c
c    This file starts with an interface sbrtn EQVMOM
c which is used to access this equilbirum package 
c entirely through an argument list.
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmomeq  .../baldur/code/bald/deqvmom.f
c  rgb 02-jun-96 initialized zpnorm = zcnorm = 0.0
c
c    Sbrtn vmomeq is used to access the VMOMS package entirely through
c an argument list.
c
      subroutine vmomeq (initin,nmhdeq,nrho,zrmaj,zbtor
     & ,elongb,elonga,tringb,phalfm,ppresm,pcurrm,mjdim,lprint,xpinit
     & ,ahalf,shift,elong,triang,r0b,rmb,ymb,ifail)
c
c    INPUT:
c  initin = 0 on the first call and 1 thereafter
c  nmhdeq = number of equations to be solved (1, 2, or 3, only)
c  nrho   = number of equilibrium flux surfaces
c  zrmaj  = major radius (midpoint at midplane) of plasma boundary (meters)
c  zbtor  = vacuum toroidal field at zrmaj (tesla)
c  elongb = geometric elongation at plasma boundary
c  elonga = elongation at magnetic axis (used for initialization)
c  tringb = geometric triangularity at plasma boundary
c  phalfm = half-width equilibrium flux surface label (meter) array
c  ppresm = plasma pressure at equilbirum flux surfaces (nt/m^2) array
c  pcurrm = plasma current within equilibrium flux surfaces (amp) array
c  mjdim  = first dimension of rmb(j,m) and ymb(j,m)
c  lprint   controls the amount of printout (0 for no printout)
c  xpinit(j) values used to initialize param(j) in sbrtn vmhdeq
c
c    OUTPUT:
c  ahalf(j)  = half-width ( meter )
c  shift(j)  = Rmajor(j) - Rmajor(nrho)  ( meter )
c  elong(j)  = Y(rho,theta_crit) / rho
c  triang(j) = (Rmajor(j) - R(rho,theta_crit) ) / rho
c    where: Rmajor(j) = R0b(j) + Rmb(j,2)
c           theta_crit is the angle at which Y(rho,theta) is maximum
c           that is: 
c       rho(j) cos ( theta_crit) + 2 Rmb(j,2) cos ( 2 theta_crit ) = 0
c  r0b(j)    = zeroth harmonic for equilibrium flux surfaces ( meter)
c  rmb(j,m)  = coeff of cos ( m * theta ) for flux surface j ( meter)
c  ymb(j,m)  = coeff of sin ( m * theta ) for flux surface j ( meter)
c
c**********************************************************************c
c
cbate      implicit none
!cap
      include 'cparm.m'
      parameter (kjdim=mj)
!
c
      integer initin,nmhdeq,nrho,mjdim,ifail,imhdeq,iparam,lprint
     & ,iptype,inrho,imode,j
c
      real zrmaj,zbtor,elongb,tringb,phalfm,ppresm, pcurrm, xpinit
     & ,ahalf,shift,elong,elonga,triang,r0b,rmb,ymb
     & ,zrmaja,zbta,zpres0,zcurnt,zlong0,zrho,zpresn,zcurrn
     & ,zcnorm,zpnorm
     & ,zshift,zelong,ztring,zr7a,zeout,zdout
c
      dimension phalfm(*), ppresm(*), pcurrm(*), xpinit(6)
     & , ahalf(*), shift(*), elong(*), triang(*)
     & , r0b(*), rmb(mjdim,*), ymb(mjdim,*)
c
      dimension zrho(kjdim), zshift(kjdim), zelong(kjdim), ztring(kjdim)
     & ,zpresn(kjdim), zcurrn(kjdim)
c
c..local variables
c
      imhdeq = 3       ! solve R0, E, and R2 equations simultaneously
      if (nmhdeq .lt. 4 .and. nmhdeq .gt. 0 ) imhdeq = nmhdeq
c
      iparam = initin  ! =0 for initial guess, =1 for previous vector p
      iptype = 0       ! geometric form of elongb and tringb given
      inrho  = nrho    ! number of radial equilibrium grid points
                       ! from magnetic axis (j=1) to edge (j=inrho)
      zrmaja = zrmaj * 100.  ! major radius  m to cm
      zbta   = zbtor * 1.e4  ! vacuum tor magnetic field  tesla to Gauss
      zpres0 = ppresm(1) * 6.242e12  ! central pressure  J/m3 to eV/cm3
      zcurnt = pcurrm(inrho)  ! total plasma current  Amps
      zlong0 = max ( elonga, 0.8 )
                       ! elongation at mag axis ( used only if imhdeq = 1 )
      ifail  = 0   ! error flag
c
      zpnorm = 0.0
      zcnorm = 0.0
c
      if ( abs(ppresm(1)) .gt. 1.e-10 ) zpnorm = 1. / ppresm(1)
      if ( abs(pcurrm(inrho)) .gt. 1.e-10 ) zcnorm = 1. / pcurrm(inrho)
c
      do 10 j=1,inrho
        zrho(j) = phalfm(j) * 100.  ! half-width on equilibrium grid
                                    ! zrho(1) = 0   m to cm
                                    ! zrho(inrho) = plasma half-width
        zpresn(j) = ppresm(j) * zpnorm     ! normalized plasma pressure
                                           ! zpresn(1) = 1.0
        zcurrn(j) = pcurrm(j) * zcnorm
                      ! normalized toroidal current within each equilibrium
                      ! flux surface.  zcurrn(1)=0., zcurrn(inrho)=1.
  10  continue
c
c..Call equilibrium solver   !! Beware !!  arguments zpresn and zcurrn
c  have been inserted into the routine provided by Wieland.
c
      call vmhdeq (imhdeq,iparam,iptype,inrho,lprint,xpinit
     & ,zrmaja,zbta,zpres0,zcurnt
     & ,zpresn,zcurrn
     & ,elongb,zlong0,tringb,zrho,zshift,zelong,ztring,ifail)
c
c..convert to harmonic representation and cm to meters
c
      imode = 1  ! to convert from geometric to moments representation
      do 20 j=1,nrho
        zr7a = zrho(j) / zrho(nrho)
        call vmgtrn (imode,zr7a,zelong(j),ztring(j),zeout,zdout)
c
c  convert from cm to meter  and harmonic representation
c
        ahalf(j) = zrho(j) * 0.01
        shift(j) = zshift(j) * 0.01
        elong(j) = zelong(j)
        triang(j) = ztring(j)
        rmb(j,1) = zrho(j) * 0.01
        rmb(j,2) = zdout * zrho(nrho)  * 0.01
        r0b(j)   = shift(j) + zrmaj - rmb(j,2)
        ymb(j,1) = zeout * ahalf(j)
        ymb(j,2) = - zeout * rmb(j,2)
  20  continue
c
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmomsr2  11040/bald89/weqvmom DEQVMOM
c***********************************************************************
c---beware...the version of {d02agEb et al} . has been modified to work
c---with vmoms, thus relying on the structure of the vector defined
c---in that code. the places where it has been patched in such a way
c---are commented as "vmoms patch". other efficiencies are denoted by
c---the comment "efficiency patch"
c***********************************************************************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmhdeq  .../baldur/code/bald/deqvmom.f
c  rgb 02-jun-96 save  neq, pi
c  rgb 10-jul-89 use xpinit(j) to initialize param(j)
c  rgb 10-jul-89 changed to param(3) = e1
c  rgb 07-jul-89 implemented 1D array zc1d for d02agEb
c     replaced m1 with nrho where ever feasible
c  rgb 06-jun-89 call vmgpsm and vmprtq only if lprint .ge. 40
c  rgb 05-jun-89 set m1 back to mm1 again due to error
c  rgb 05-jun-89 set m1=nrho rather than m1=mm1
c
      subroutine vmhdeq(kmhdeq, iparam, iptype, irho, iprint, xpinit,
     &     xr0, xbt0, xp0,  xai0, zpres, zcurt,
     *     zxe0, zxe1, zxd0, xrho, xshift, xelong, xtriang, ifail)
     *     
      external case14
c***********************************************************************
c*****vmhdeq solves for the moments of the grad-shafranov equation.    *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****                                ,phys fluids 24 (1981) aug       *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****               24-apr-89 Glenn Bateman                           *
c*****calculated parameters:                                           *
c*****shift(i)-shift of surface i relative to mag axis-(cm).           *
c*****elong(i)-elongation of surface i-(dimensionless).                *
c*****triang(i)-triangularity of surface i-(dimensionless).            *
c*****input parameters:                                                *
c*****kmhdeq-number of equations to be solved.                         *
c*****      =1 solve shift, approximate elongation and triangularity.  *
c*****      =2 solve shift and elongation, approximate triangularity.  *
c*****      =3 solve shift, elongation, and triangularity.             *
c*****iparam-switch for initial guess at solution.                     *
c*****      =0 construct guess from zxe0, zxe1, and zxd0.              *
c*****      >0 use solution from previous call if available.           *
c*****iptype-switch to specify what "kind" of e&d are coming in        *
c*****      =0 they are "geometrical" and need to be "vmgtrn-ed"       *
c*****      =1 they are "parameterized" and can be used directly       *
c*****nrho-number of radial grid points.                               *
c*****iprint controls amout of printout (argument inserted by Bateman) *
c*****      = 0 for no printout during normal computation              *
c*****      = 10 for printout of full profiles using sbrtn vmprtq      *
cbate xpinit(j) input values used to initialize param(j)               *
c*****r0-major radius to geometric center of outermost surface-(cm).   *
c*****bt0-vacuum toroidal field at r0-(gauss).                         *
c*****p0-plasma pressure at magnetic axis-(ev/cm3).                    *
c*****ai0-total toroidal current-(amps).                               *
c*****   Arguments inserted by Bateman 28-mar-89:                      *
c*****zpres normalized plasma pressure on equilibrium grid points      *
c*****zcurt normalized toroidal current within flux surfaces           *
c*****   Back to original arguments:                                   *
c*****zxe0-elongation of outermost flux surface-(dimensionless).       *
c*****zxe1-elongation on axis for approximate solution-(dimensionless).*
c*****zxd0-triangularity of outermost flux surface-(dimensionless).    *
c*****rho(i)-half-diameter of surface i in midplane-(cm).              *
c*****other comments:                                                  *
c*****let rg(i) be the major radius to the geometric center of flux    *
c*****surface i                                                        *
c*****then shift(i)=rmag-r0.                                           *
c*****let zmax(i) be the maximum height of flux surface i, r(zmax(i))  *
c*****the distance from this point to the major axis then              *
c*****elong(i)=zmax(i)/rho(i).                                         *
c*****triang(i)=(rg(i)-r(zmax(i)))/rho(i).                             *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      dimension xrho(irho), xshift(irho), xelong(irho), xpinit(6),
     *     xtriang(irho), zpres(irho), zcurt(irho), zc1d(mm1neq)
      external vmaux, vmbcux, vmraux, vmprsl
      data erry /1.e-04/, errp /1.e-05/
      data xl00 /0.03/
      data cvel /3.0e+10/, stata /3.0e+9/, everg /1.6022e-12/
c
      save  neq, pi
c
c*****initialization.
c
      if (initm.eq.1) go to 30
c
c..unit numbers added by Bateman, 31-mar-89
c
      nin    = 5
      nout   = 6
      nplot  = 6
      ndebug = 6
      nttyi  = 5
      nttyo  = 6
c
cbate      epsmch = 1.
cbate   2  epsmch = 0.5 * epsmch
cbate      if ( 0.5 * epsmch + 1.0 .gt. 1.0 ) go to 2
c
      call machinb (epslon)
        epsmch = epslon
        epslni = 1. / epslon
c
      ntheta = mtheta
      mtm2 = mtheta - 2
      neq = mneq
cbate      m1     = mm1
      pi = acos(-1.0)
c*****set up sine and cosine arrays.
      call vmtrig
c*****set up the integration weighting array.
      call vmwgts(ntheta,wfcn,cnorm)!   theta weighting fnc
cbate      call vmwgts(m1,wm1,cm1)!  radial weighting fnc
c
      initm  = 1
c
c..end of initializatiion
c
   30 continue
c
      if ( irho .gt. mm1 ) call abortb (6
     &  ,' irho .gt. mm1 in sbrtn vmhdeq in file DEQVMOM')
c
cbate      m1 = mm1
      m1 = irho
      call vmwgts(m1,wm1,cm1)!  radial weighting fnc
c
      lprint = iprint
      nfl    = 0
      nrho   = irho
      xl0    = xl00
      irot   = 0
      hh     = hh0
c
      do 32 j=1,mneq
        error(j)  = erry
        parerr(j) = errp
  32  continue
c
      kmhd = kmhdeq
      r0 = xr0
      bt0 = xbt0
      p0 = xp0
      ai0 = xai0
      initc = 0
      initp = 0
      itpr = 0
c
c..transfer pressure and current...  Bateman 28-mar-89
c
      do 10 j=1,nrho
        rho(j)  = xrho(j)
        pres(j) = zpres(j)
        curt(j) = zcurt(j)
  10  continue
c
c*****set up vmaux/vmgpsm constants
      a0 = xrho(nrho)
      argeom = r0 / a0
      fwall=r0*bt0
      f0 = cvel/2.0*fwall/stata
      ajt0 = (f0*stata)**2/(2.0*pi*a0**2*argeom*(ai0*stata))
      q0=f0/ai0/argeom**2
      betap0 = 2.0*pi*cvel**2*(p0*everg)*a0**2/(ai0*stata)**2
      betat0 = (p0*everg)*4.0*pi/bt0**2
c*****set up dimensionless radial grid.
      call vmsetx (xl0,xrho(2),a0,nrho,xp,dx)
c*****set up rotation form factors
      if (irot.eq.1) call vmnrot
c*****convert to moments representation for elongation and triangularity
c  note:  e0,e1,d0  are moments representation
c         xe0, xe1, and xd0 are geometric representation
c
      if (iptype.eq.0) then
        xe0 = zxe0
        xe1 = zxe1
        xd0 = zxd0
        e1 = zxe1
        call vmgtrn( 1, 1.0, zxe0, zxd0, e0, d0 )
      else
        e0 = zxe0
        d0 = zxd0
        e1 = zxe1
        xe1 = zxe1
        call vmgtrn( 2, 1.0, zxe0, zxd0, xe0, xd0 )
      endif
c
      arc = (r0-d0*a0)/a0
      if (iparam.ne.0 .and. param(3).ne.0.0) go to 50
c
c*****set up the initial guess for the parameters.
c
      param(1) = arc + xpinit(1)  ! was param(1) = arc + 0.1
      param(2) = xpinit(2)        ! was param(2) = -0.4
      param(3) = e1
                    ! was if (kmhd.gt.1) param(3) = 1.0 + (e0-1.0)/2.0
      param(4) = xpinit(4)        ! was param(4) = 0.0
      param(5) = xpinit(5)        ! was param(5) = .01
      param(6) = xpinit(6)*d0     ! was param(6) = 2.*d0
c
c*****call the ode solver.
c
   50 ifail = 1
      call d02agEb(hh, error, parerr, param, zc1d, neq, 2*kmhd, nrho,
     *     vmaux, vmbcux, vmraux, vmprsl, mat, copy, wspace, wspac1,
     *     wspac2, ifail)
      if (iunder.eq.1) ifail = 8
      if (ifail.gt.0) return
c
c..copy from 1D array zc1d to 2D array c(jr,jneq)
c
      do 52 jn=1,neq
      do 52 jr=1,nrho
        i = nrho * (jn-1) + jr
        c(jr,jn) = zc1d (i)
  52  continue
c
c*****interpolate onto the driver grid.
      call splaan(nrho, xp, c(1,1), bsift, csift, dsift)
      call splaan(nrho, xp, c(1,3), belli, celli, delli)
      call splaan(nrho, xp, c(1,5), btria, ctria, dtria)
      do 60 j=2,nrho
           x = xrho(j)/a0
           xshift(j) = seval(nrho,x,xp,c(1,1),bsift,csift,dsift)
           emom = seval(nrho,x,xp,c(1,3),belli,celli,delli)
           trmom = seval(nrho,x,xp,c(1,5),btria,ctria,dtria)
           xshift(j) = (xshift(j)+trmom)*a0 - r0
           call vmgtrn(2, x, emom, trmom, xelong(j), xtriang(j))
   60 continue
      xshift(1) = param(1)*a0 - r0
      xelong(1) = param(3)
      xtriang(1) = 0.0
c
      do 62 j=1,nrho
        shift(j) = xshift(j)
        elong(j) = xelong(j)
        triang(j) = xtriang(j)
  62  continue
c
c..printout
c
      if ( lprint .ge. 40 ) then
        call vmgpsm
        call vmprtq
      endif
c
      return
      end
c******************
      BLOCK DATA case14
c
      include 'tvmcom0.cmm'
      integer initm, initc, initp, itpr, lprint
      common /vcom4/ initm, initc, initp, itpr, lprint
c
c      integer mneq
c      parameter (mneq=6)
      real error(mneq),parerr(mneq),mat(mneq,mneq),wspace(mneq,9),
     1     wspac1(mneq),wspac2(mneq),copy(mneq,mneq),
     2     hh,hh0,epslni

      common /vcom6/ error, parerr, mat, wspace, wspac1,
     *     wspac2, copy, hh, hh0, epslni, m1, ntheta
c
      real xl0,dx,r0,a0,bt0,p0,ai0,e0,e1,d0,arc,slpe,slpr0,betap0
      integer kmhd,iunder
c
      common /vcom8/ xl0, dx, r0, a0, bt0, p0, ai0, e0, e1,
     *     d0, arc, slpe, slpr0, betap0, kmhd, iunder

c      include 'tvmcom1.cmm'
      data slpr0 /1.0e-4/, slpe /1.0e-4/, hh0 /0.01/, initm /0/
      end


c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmgtrn  11040/bald89/weqvmom DEQVMOM
c  (used to be sbrtn geotrn)
      subroutine vmgtrn(mode, xin, ein, din, eout, dout)
c***********************************************************************
c*****vmgtrn transforms elongation and triangularity from a moments to *
c*****a geometrical representation and vice versa.                     *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****calculated parameters:                                           *
c*****eout-output elongation.                                          *
c*****dout-output triangularity.                                       *
c*****input parameters:                                                *
c*****mode-designates direction of transformation.                     *
c*****    =1 geometrical => moments.                                   *
c*****    =2 moments => geometrical.                                   *
c*****xin-input reduced minor radius                                   *
c*****ein-input elongation.                                            *
c*****din-input triangularity.                                         *
c***********************************************************************
      if (xin.le.0.0) go to 30
      if (mode.eq.2) go to 20
c*****geometrical ==> moments representation.
      dout = din/4.0
      ctc = 0.0
      do 10 i=1,10
           ctc = 4.0*dout/(sqrt(xin**2+32.0*dout**2)+xin)
           dout = xin*din/(4.0-6.0*ctc**2)
   10 continue
      eout = xin*ein/(sqrt(1.0-ctc*ctc)*(xin+2.0*dout*ctc))
      return
c*****moments ==> geometrical representation.
   20 ctc = 4.0*din/(sqrt(xin**2+32.0*din**2)+xin)
      stc = sqrt(1.0-ctc*ctc)
      s2tc = 2.0*stc*ctc
      c2tc = 2.0*ctc*ctc - 1.0
      eout = ein*(xin*stc+din*s2tc)/xin
      dout = din*(1.0-3.0*c2tc)/xin
      return
   30 eout = ein
      dout = din
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmnrot  11040/bald89/weqvmom DEQVMOM
c  rgb 05-jun-89 replaced mm1 by m1 in do 30 and call splaan
c
      subroutine vmnrot

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

c***in the form used here, the pressure term looks like
c     p(x,r) = ptilda(x)*exp[u*r**2*w**2/2/t]
c        where u is the amu of the rotating plasma species
c              r = r(x,theta)
c              w = "omega", the angular velocity
c              t is the temperature of the rotating species
c
c***we write here
c     p(x,r) = ptilda(x) * exp[f(x)*(r/r0)**2)]
c       where the form factor f(x) = u*(w(x)*r0)**2/2/t
c       or f(x) = (va0-va1)*(1-x**va2)**va3 + va1
c units-wise: va0/1 = u[amu] * (wr0)[10**7cm/s] / t[kev] * .052
      do 30 i = 1, m1
        formr(i) = (va0-va1)*(1.-xp(i)**va2)**va3
   30 continue
      call splaan(m1,xp,formr,formr1,formr2,formr3)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmaux  11040/bald89/weqvmom DEQVMOM
c  rgb 06-aug-90 added ifail to argument list
c    ifail = 9 if matrix is singular
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  this routine evaluates the derivative of the functions y(x)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vmaux(yp,y,x,paramx,ifail)

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      common/zcoefr0/dr010,dr012,dr021,dr022,dr03,dr040,dr041,dr05
      common/zcoefr2/dr210,dr212,dr221,dr222,dr23,dr240,dr241,dr25
      common/zcoefz1/dz110,dz112,dz121,dz122,dz13,dz140,dz141,dz15
      common/zaaa/al11,al21,al31,al12,al22,al32,al13,al23,al33
      common/zbbb/bl1,bl2,bl3
      dimension abin(3,3),bbin(3),ipvt(3),work(3)
      equivalence (al11,abin), (bl1,bbin)
      dimension yp(mneq),y(mneq),paramx(mneq)
      dimension rwp(mtheta),wp(mtheta)
c
      nfl=nfl+1
c
c
c  get the coeficients in each moment equation
c  then
c  transform the systems of second order odes into first
c  order by solving a system of algebric equations
c
c                     al*yp = bl
c                 *****************
c
      call vmpfun(x,px,ppx)
      call vmaifn(x,aix,aipx)
      call vmwrot(x,w1,w2)
      aipoi=aipx/aix
c
      ai=betap0/2.0/aix**2
c
c
c  vectorize the coeficients in the moment equations
c
      do 50 i=1,ntheta
        r(i)=y(1)-x*c1th(i)+y(5)*c2th(i)
        rx(i)=y(2)-c1th(i)+y(6)*c2th(i)
        rt(i)=x*s1th(i)-2.0*y(5)*s2th(i)
        rxt(i)=s1th(i)-2.0*y(6)*s2th(i)
        rtt(i)=x*c1th(i)-4.0*y(5)*c2th(i)
c
        z(i)=y(3)*(x*s1th(i)+y(5)*s2th(i))
        zx(i)=y(4)*(x*s1th(i)+y(5)*s2th(i))+y(3)*(s1th(i)+y(6)*
     1    s2th(i))
        zt(i)=y(3)*(x*c1th(i)+2.0*y(5)*c2th(i))
        zxt(i)=y(4)*(x*c1th(i)+2.0*y(5)*c2th(i))+y(3)*(c1th(i)+
     1    2.0*y(6)*c2th(i))
        ztt(i)=-y(3)*(x*s1th(i)+4.0*y(5)*s2th(i))
c
        gtt(i)=rt(i)**2+zt(i)**2
        tau(i)=rt(i)*zx(i)-rx(i)*zt(i)
        gsqrt(i)=r(i)*tau(i)
        rtau2(i)=gsqrt(i)*tau(i)
        rtau3(i)=rtau2(i)*tau(i)
c
        gs10(i)=zt(i)*gtt(i)/rtau3(i)
        gs12(i)=gs10(i)*c2th(i)
        gs2(i)=rt(i)*gtt(i)/rtau3(i)
        gs21(i)=gs2(i)*s1th(i)
        gs22(i)=gs2(i)*s2th(i)
        gs4(i)=1.0/r(i)
        gs5(i)=zt(i)/r(i)/gsqrt(i)
        cs13n(i)=gtt(i)*(rx(i)*zxt(i)-zx(i)*rxt(i))
        gttx(i)=rt(i)*rxt(i)+zt(i)*zxt(i)
        gs6(i)=(cs13n(i)+(rx(i)*rt(i)
     1    +zx(i)*zt(i))*(rtt(i)*zx(i)+rt(i)*zxt(i)-rxt(i)*zt(i)
     1    -rx(i)*ztt(i)))/rtau3(i)
        gs7(i)=(gttx(i)-rx(i)*rtt(i)-zx(i)*
     1    ztt(i))/rtau2(i)
        gs8(i)=gtt(i)/rtau2(i)
        ss20(i)=gs5(i)+gs6(i)+gs7(i)
        cs1n(i)=(2.0*gttx(i)-gtt(i)*rx(i)/r(i)+cs13n(i)/tau(i))
     1    /gsqrt(i)
c
        rmom0(i)=zt(i)*wfcn(i)
        tmom(i)=tau(i)*wfcn(i)
        rmom2(i)=(y(3)*rt(i)*s2th(i)-zt(i)*c2th(i))*wfcn(i)
        zmom1(i)=(+x*s1th(i)+y(5)*s2th(i))*rt(i)*wfcn(i)! 01/30/83
c
        rarg = ( r(i) / argeom ) ** 2
        wp(i)=w1**(rarg) * (ppx + px*w2*rarg)
        rwp(i)=wp(i)*r(i)
   50 continue
c
      yp(1)=y(2)
      yp(3)=y(4)
      yp(5)=y(6)
c
c     ****************
c
c     **********   r0 moment package   **********
c
      cr010=sdot(ntheta,rmom0,1,gs10,1)
      cr012=sdot(ntheta,rmom0,1,gs12,1)
      cr021=sdot(ntheta,rmom0,1,gs21,1)
      cr022=sdot(ntheta,rmom0,1,gs22,1)
      cr03=sdot(ntheta,rmom0,1,rwp,1)
      cr04=sdot(ntheta,rmom0,1,gs4,1)
      cr08=sdot(ntheta,rmom0,1,gs8,1)
      cr0ss2=sdot(ntheta,rmom0,1,ss20,1)
c
      cpsi10=sdot(ntheta,tmom,1,gs10,1)
      cpsi12=sdot(ntheta,tmom,1,gs12,1)
      cpsi21=sdot(ntheta,tmom,1,gs21,1)
      cpsi22=sdot(ntheta,tmom,1,gs22,1)
      cpsi3=sdot(ntheta,tmom,1,rwp,1)
      cpsi4=sdot(ntheta,tmom,1,gs4,1)
      cpsi8=sdot(ntheta,tmom,1,gs8,1)
      csiss2=sdot(ntheta,tmom,1,ss20,1)
      ccs1=sdot(ntheta,wfcn,1,cs1n,1)/cpsi8
c
      cr0gs4=cr04/cpsi4
      cr0gs8=cr08/cpsi8
      dr010=cr010-cr0gs8*cpsi10
      dr012=cr012-cr0gs8*cpsi12
      dr021=cr021-cr0gs8*cpsi21
      dr022=cr022-cr0gs8*cpsi22
      dr03=cpsi8**2/cnorm**2*(cr03-cr0gs4*cpsi3)
      dr040=cr0ss2-cr0gs4*csiss2
      dr041=cr08-cr0gs4*cpsi8
      dr05=ccs1*(cr08-cr0gs4*cpsi8)
c
c
c
c
      al11=dr010
      al12=-(x*dr021+y(5)*dr022)
      al13=dr012-y(3)*dr022
      bl1=-ai*dr03-dr040-dr041*aipoi+dr05+2.0*y(4)*(dr021+y(6)*dr022)
c
      if(kmhd.gt.1) go to 100
c
c***  kmhd=1
      bl3=2.*d0
      bl2=2.*(e0-paramx(3))
      bl1=(bl1-al13*bl3-al12*bl2)/al11
      go to 300
c
c     *******************************************
c
c
c
c     **********   z0 moment package   **********
c
100   continue
      cz110=sdot(ntheta,zmom1,1,gs10,1)
      cz112=sdot(ntheta,zmom1,1,gs12,1)
      cz121=sdot(ntheta,zmom1,1,gs21,1)
      cz122=sdot(ntheta,zmom1,1,gs22,1)
      cz13=sdot(ntheta,zmom1,1,rwp,1)
      cz14=sdot(ntheta,zmom1,1,gs4,1)
      cz18=sdot(ntheta,zmom1,1,gs8,1)
      cz1ss2=sdot(ntheta,zmom1,1,ss20,1)
c
      cr210=sdot(ntheta,rmom2,1,gs10,1)
      cr212=sdot(ntheta,rmom2,1,gs12,1)
      cr221=sdot(ntheta,rmom2,1,gs21,1)
      cr222=sdot(ntheta,rmom2,1,gs22,1)
      cr23=sdot(ntheta,rmom2,1,rwp,1)
      cr24=sdot(ntheta,rmom2,1,gs4,1)
      cr28=sdot(ntheta,rmom2,1,gs8,1)
      cr2ss2=sdot(ntheta,rmom2,1,ss20,1)
c
      cz1gs4=cz14/cpsi4
      cz1gs8=cz18/cpsi8
      dz110=cz110-cz1gs8*cpsi10
      dz112=cz112-cz1gs8*cpsi12
      dz121=cz121-cz1gs8*cpsi21
      dz122=cz122-cz1gs8*cpsi22
      dz13=cpsi8**2/cnorm**2*(cz13-cz1gs4*cpsi3)
      dz140=cz1ss2-cz1gs4*csiss2
      dz141=cz18-cz1gs4*cpsi8
      dz15=ccs1*(cz18-cz1gs4*cpsi8)
c
      al21=dz110
      al22=-(x*dz121+y(5)*dz122)
      al23=dz112-y(3)*dz122
      bl2=-ai*dz13-dz140-dz141*aipoi+dz15+2.0*y(4)*(dz121+y(6)*dz122)
c
      if(kmhd.eq.3) go to 150
c
      bl3=2.*d0
      bl1=bl1-al13*bl3
      bl2=bl2-al23*bl3
      go to 200
c
c     *******************************************
c
c
c
c     **********   r2 moment package   **********
c
150    continue
      cr2gs4=cr24/cpsi4
      cr2gs8=cr28/cpsi8
      dr210=cr210-cr2gs8*cpsi10
      dr212=cr212-cr2gs8*cpsi12
      dr221=cr221-cr2gs8*cpsi21
      dr222=cr222-cr2gs8*cpsi22
      dr23=cpsi8**2/cnorm**2*(cr23-cr2gs4*cpsi3)
      dr240=cr2ss2-cr2gs4*csiss2
      dr241=cr28-cr2gs4*cpsi8
      dr25=ccs1*(cr28-cr2gs4*cpsi8)
      al31=dr210
      al32=-(x*dr221+y(5)*dr222)
      al33=dr212-y(3)*dr222
      bl3=-ai*dr23-dr240-dr241*aipoi+dr25+2.0*y(4)*(dr221+y(6)*dr222)
c
c     *******************************************
c
200   call lqdcmp(mnmoms,kmhd,abin,cond,ipvt,work)
      if(cond .gt. epslni) go to 1000
c
      call lqsolv(mnmoms,kmhd,abin,bbin,ipvt)
c
300   yp(2)=bbin(1)
      yp(4)=bbin(2)
      yp(6)=bbin(3)
c
400   continue
c
      return
c
1000  continue
      write(nttyo,1010) x,cond,y
1010  format(' almost-singular system in "vmaux" at x= ',f5.3,/,
     1  ' condition # = ',1pe10.2,/,' y= ',6e10.2)
cbate        call crash
cbate        call bturky
cbate        call abortr
      ifail = 9
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmbcux  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmbcux(gl, gr, paramx)
c***********************************************************************
c*****vmbcux sets up the boundary conditions for the mhd moments.       *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****                                ,phys fluids 24 (1981) aug       *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****calculated parameters:                                           *
c*****gl(i)-value of y(i) at left boundary.                            *
c*****gr(i)-value of y(i) at right boundary.                           *
c*****input parameters:                                                *
c*****param(i)-b.c. parameter for first order ode i.                   *
c*****xl0-lower bound of abscissa.                                     *
c*****slpr0-analytic parameter for left b.c. at x1.                    *
c*****slpe-analytic parameter for left b.c. at x1.                     *
c*****arc-inverse aspect ratio for outermost flux surface.             *
c*****e0-moments representation of elongation at x=1.                  *
c*****d0-moments representation of triangularity at x=1.               *
c*****other comments:                                                  *
c*****vmbcux is required by d02agEb ode solver in nag library.           *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      dimension gl(mneq), gr(mneq), paramx(mneq)

c*****left b.c.'s - center.
c
      gl(1) = paramx(1) - slpr0*xl0**2
      gl(2) = -2.0*slpr0*xl0
c
cbate commented out lines received from wieland in March 1989
cbate      gl(1) = paramx(1) + (arc-paramx(1))*xl0**2
cbate      gl(2) = 2.0*(arc-paramx(1))*xl0
c
      gl(3) = paramx(3) + slpe*xl0**2
      gl(4) = 2.0*xl0*slpe
cbate      gl(3) = paramx(3)
cbate      gl(4) = 0.
      gl(5) = paramx(5)*xl0**2
      gl(6) = 2*paramx(5)*xl0
c*****right b.c.'s - edge.
      gr(1) = arc
      gr(2) = paramx(2)
      gr(3) = e0
      gr(4) = paramx(4)
      gr(5) = d0
      gr(6) = paramx(6)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmraux  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmraux(xl, xr, rmatch, zparam)
c***********************************************************************
c*****vmraux sets the matching point for the mhd moments.               *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****                                ,phys fluids 24 (1981) aug       *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****calculated parameters:                                           *
c*****xl-left (lower) limit on abscissa.                               *
c*****xr-right (upper) limit on abscissa.                              *
c*****rmatch-matching point for shooting method (rmatch=1.0).          *
c*****input parameters:                                                *
c*****xl0-lower bound of abscissa.                                     *
c*****other comments:                                                  *
c*****vmraux is required by d02agEb ode solver in nag library.           *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      dimension zparam(mneq)
      xl = xl0
      xr = 1.0
      rmatch = 1.0
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmprsl  11040/bald89/weqvmom DEQVMOM
c  rgb 06-jun-89 commented out debug output
c
      subroutine vmprsl(zparam, res, n1, err)
c***********************************************************************
c*****vmprsl prints out the iteration parameters for the mhd moments.   *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****                                ,phys fluids 24 (1981) aug       *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****input parameters:                                                *
c*****param(i)-b.c. parameter for first order ode i.                   *
c*****res-sum of squares of err(i).                                    *
c*****n1-number of first order ode's being solved.                     *
c*****err(i)-difference between the two solutions at matching point.   *
c*****other comments:                                                  *
c*****vmprsl is required by d02agEb ode solver in nag library.           *
c*****see d02agEb write-up for optional diagnostics.                    *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      dimension zparam(mneq), err(mneq)
      if (itpr.gt.0) go to 10
cbate      write (ndebug,99999) vtg1,vtg2
cbate      write (ndebug,99998)
   10 itpr = itpr + 1
cbate      write (ndebug,99997) itpr, (zparam(i),i=1,n1)
cbate      write (ndebug,99996) (err(i),i=1,n1), res
      write (nout,99997) itpr, (zparam(i),i=1,n1)
      write (nout,99996) (err(i),i=1,n1), res
      return
cbate99999 format (/, 1x, 20(1h*), 'vmoms output at tg1=',f10.3,
cbate     *     ' tg2= ',f10.3, 3x, 20(1h*) )
cbate99998 format (1x, 20(1h*), ' newton iteration summary ',
cbate     *     20(1h*))
99997 format (1x, ' iteration = ', i3, /, 1x, ' param =',
     *     1h , 6(1pe10.3, 2x))
99996 format (1x, ' err   = ', 6(1pe10.3, 2x),/, '  residue = ',
     *     1pe10.3)
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmaifn  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmaifn(x, ai, aip)
c***********************************************************************
c*****vmaifn provides the integrated toroidal current profile.         *
c*****references:                                                      *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****calculated parameters:                                           *
c*****ai-normalized integrated toroidal current at x-(dimensionless).  *
c*****aip-derivative of ai with respect to x-(dimensionless).          *
c*****input parameters:                                                *
c*****nrho-number of radial grid points.                               *
c*****x-radial position-(dimensionless).                               *
c*****curt(i)-toroidal current integrated to node i.                   *
c*****rho(i)-radial coordinate for node i.                             *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      if (initc.ne.0) go to 10
      call splaan(nrho, rho, curt, bcurt, ccurt, dcurt)
      initc = 1
   10 rloc = x*rho(nrho)
      ai = seval(nrho,rloc,rho,curt,bcurt,ccurt,dcurt)/curt(nrho)
      aip = speval(nrho,rloc,rho,curt,bcurt,ccurt,dcurt)*rho(nrho)/
     *     curt(nrho)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmpfun  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmpfun(x, p, pp)
c***********************************************************************
c*****vmpfun provides the plasma pressure profile.                     *
c*****references:                                                      *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c*****calculated parameters:                                           *
c*****p-normalized pressure profile at x-(dimensionless).              *
c*****pp-derivative of p with respect to x-(dimensionless).            *
c*****input parameters:                                                *
c*****nrho-number of radial grid points.                               *
c*****x-radial position-(dimensionless).                               *
c*****pres(i)-plasma pressure at node i.                               *
c*****rho(i)-radial coordinate for node i.                             *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      if (initp.ne.0) go to 10
      call splaan(nrho, rho, pres, bpres, cpres, dpres)
      initp = 1
   10 rloc = x*rho(nrho)
      zpnorm = 1.
      if ( abs(pres(1)) .gt. 1.e-10 ) zpnorm = 1. / pres(1)
      p = seval(nrho,rloc,rho,pres,bpres,cpres,dpres)*zpnorm
      pp = speval(nrho,rloc,rho,pres,bpres,cpres,dpres)*rho(nrho)*zpnorm
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmtrig  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmtrig
c***********************************************************************
c*****vmtrig sets up the sine and cosine arrays for the mhd moments.   *
c*****references:                                                      *
c*****l.l.lao,s.p.hirshman,r.m.wieland,ornl/tm-7616 (1981).            *
c*****                                ,phys fluids 24 (1981) aug       *
c*****l.l.lao,r.m.wieland,w.a.houlberg,s.p.hirshman,ornl/tm-7871 (1981)*
c*****last revision: 6/81 l.l.lao, r.m.wieland, and w.a.houlberg ornl. *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      pi = acos(-1.0)
      dth = pi/(mtheta-1)
      do 10 i=1,mtheta
           x1 = (i-1)*dth
           xtheta(i) = x1
           x2 = 2.0*x1
           s1th(i) = sin(x1)
           s2th(i) = sin(x2)
           c1th(i) = cos(x1)
           c2th(i) = cos(x2)
   10 continue
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmsetx  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmsetx (xl0,xrho,a0,m1,xp,dx)
      real   xl0,xrho,a0,xp(m1),dx
      integer   m1
      xl0 = min(xl0,xrho/a0)
      dx = (1.0-xl0)/(m1-1)
      do 40 i=1,m1
           xp(i) = xl0 + (i-1)*dx
   40 continue
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmwgts  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmwgts(ntheta,wfcn,cnorm)
      dimension wfcn(ntheta)
c
c set up integration weighting array
c
      wfcn(1)=1.0
      do 160 i=2,ntheta-2,2
        wfcn(i)=4.0
        wfcn(i+1)=2.0
  160 continue
      wfcn(ntheta-1)=4.0
      wfcn(ntheta)=1.0
      cnorm=3.0*(ntheta-1)
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmwgtt  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmwgtt (n,ab,w,c)
c-----trapezoidal rule weighting function
      integer n,i
      real w(n),c,ab
      c = (n-1)/ab
      w(1) = 0.5
      w(n) = 0.5
      if (n.gt.2) then
          do 10 i = 2, n-1
              w(i) = 1.0
10        continue
      end if
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmwrot  11040/bald89/weqvmom DEQVMOM
c  rgb 05-jun-89 replace mm1 by m1 in arg list of speval
c
      subroutine vmwrot(x,w1,w2)

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'
c
      w1=1.
      w2=0.
      if(irot.eq.0) return
c-------------------------
      wx=seval(m1,x,xp,formr,formr1,formr2,formr3)
      wxp=speval(m1,x,xp,formr,formr1,formr2,formr3)
      w1=exp(wx)
      w2=wxp
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@lqdcmp  11040/bald89/weqvmom DEQVMOM
c  rgb 05-jun-89 changed sbrtn name from decomp to lqdcmp to avoid
c     name conflicts with the BALDUR code.
c
      subroutine lqdcmp(ndim,n,a,cond,ipvt,work)
c
      integer ndim,n
      real a(ndim,n),cond,work(n)
      integer ipvt(n)
c
c     decomposes a real matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     use sbrtn lqsolv to compute solutions to linear systems.
c
c     input..
c
c        ndim = declared row dimension of the array containing  a.
c
c        n = order of the matrix.
c
c        a = matrix to be triangularized.
c
c     output..
c
c        a  contains an upper triangular matrix  u  and a permuted
c          version of a lower triangular matrix  i-l  so that
c          (permutation matrix)*a = l*u .
c
c        cond = an estimate of the condition of  a .
c           for the linear system  a*x = b, changes in  a  and  b
c           may cause changes  cond  times as large in  x .
c           if  cond+1.0 .eq. cond , a is singular to working
c           precision.  cond  is set to  1.0e+32  if exact
c           singularity is detected.
c
c        ipvt = the pivot vector.
c           ipvt(k) = the index of the k-th pivot row
c           ipvt(n) = (-1)**(number of interchanges)
c
c     work space..  the vector  work  must be declared and included
c                   in the call.  its input contents are ignored.
c                   its output contents are usually unimportant.
c
c     the determinant of a can be obtained on output by
c        det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n).
c
      real ek, t, anorm, ynorm, znorm
      integer nm1, i, j, k, kp1, kb, km1, m
c
      ipvt(n) = 1
      if (n .eq. 1) go to 80
      nm1 = n - 1
c
c     compute 1-norm of a
c
      anorm = 0.0
      do 10 j = 1, n
         t = 0.0
         do 5 i = 1, n
            t = t + abs(a(i,j))
    5    continue
         if (t .gt. anorm) anorm = t
   10 continue
c
c     gaussian elimination with partial pivoting
c
      do 35 k = 1,nm1
         kp1= k+1
c
c        find pivot
c
         m = k
         do 15 i = kp1,n
            if (abs(a(i,k)) .gt. abs(a(m,k))) m = i
   15    continue
         ipvt(k) = m
         if (m .ne. k) ipvt(n) = -ipvt(n)
         t = a(m,k)
         a(m,k) = a(k,k)
         a(k,k) = t
c
c        skip step if pivot is zero
c
         if (t .eq. 0.0) go to 35
c
c        compute multipliers
c
         do 20 i = kp1,n
             a(i,k) = -a(i,k)/t
   20    continue
c
c        interchange and eliminate by columns
c
         do 30 j = kp1,n
             t = a(m,j)
             a(m,j) = a(k,j)
             a(k,j) = t
             if (t .eq. 0.0) go to 30
             do 25 i = kp1,n
                a(i,j) = a(i,j) + a(i,k)*t
   25        continue
   30    continue
   35 continue
c
c     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
c     estimate obtained by one step of inverse iteration for the
c     small singular vector.  this involves solving two systems
c     of equations, (a-transpose)*y = e  and  a*z = y  where  e
c     is a vector of +1 or -1 chosen to cause growth in y.
c     estimate = (1-norm of z)/(1-norm of y)
c
c     solve (a-transpose)*y = e
c
      do 50 k = 1, n
         t = 0.0
         if (k .eq. 1) go to 45
         km1 = k-1
         do 40 i = 1, km1
            t = t + a(i,k)*work(i)
   40    continue
   45    ek = 1.0
         if (t .lt. 0.0) ek = -1.0
         if (a(k,k) .eq. 0.0) go to 90
         work(k) = -(ek + t)/a(k,k)
   50 continue
      do 60 kb = 1, nm1
         k = n - kb
         t = 0.0
         kp1 = k+1
         do 55 i = kp1, n
            t = t + a(i,k)*work(k)
   55    continue
         work(k) = t
         m = ipvt(k)
         if (m .eq. k) go to 60
         t = work(m)
         work(m) = work(k)
         work(k) = t
   60 continue
c
      ynorm = 0.0
      do 65 i = 1, n
         ynorm = ynorm + abs(work(i))
   65 continue
c
c     solve a*z = y
c
      call lqsolv(ndim, n, a, work, ipvt)
c
      znorm = 0.0
      do 70 i = 1, n
         znorm = znorm + abs(work(i))
   70 continue
c
c     estimate condition
c
      cond = anorm*znorm/ynorm
      if (cond .lt. 1.0) cond = 1.0
      return
c
c     1-by-1
c
   80 cond = 1.0
      if (a(1,1) .ne. 0.0) return
c
c     exact singularity
c
   90 cond = 1.0e+32
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@lqsolv  11040/bald89/weqvmom DEQVMOM
c  rgb 05-jun-89 changed sbrtn name from solve to lqsolv to avoid
c     name conflicts with the BALDUR code.
c
      subroutine lqsolv(ndim, n, a, b, ipvt)
c
      integer ndim, n, ipvt(n)
      real a(ndim,n),b(n)
c
c   solution of linear system, a*x = b .
c   do not use if decomp has detected singularity.
c
c   input..
c
c     ndim = declared row dimension of array containing a .
c
c     n = order of matrix.
c
c     a = triangularized matrix obtained from decomp .
c
c     b = right hand side vector.
c
c     ipvt = pivot vector obtained from decomp .
c
c   output..
c
c     b = solution vector, x .
c
      integer kb, km1, nm1, kp1, i, k, m
      real t
c
c     forward elimination
c
      if (n .eq. 1) go to 50
      nm1 = n-1
      do 20 k = 1, nm1
         kp1 = k+1
         m = ipvt(k)
         t = b(m)
         b(m) = b(k)
         b(k) = t
         do 10 i = kp1, n
             b(i) = b(i) + a(i,k)*t
   10    continue
   20 continue
c
c     back substitution
c
      do 40 kb = 1,nm1
         km1 = n-kb
         k = km1+1
         b(k) = b(k)/a(k,k)
         t = -b(k)
         do 30 i = 1, km1
             b(i) = b(i) + a(i,k)*t
   30    continue
   40 continue
   50 b(1) = b(1)/a(1,1)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@machinb  11040/bald89/weqvmom DEQVMOM
c
      subroutine machinb(u)
c
c   u is the smallest positive number such that (1.0+u) .gt. 1.0 .
c   u is computed approximately as a power of 1./2.
c
c   this code is completely explained and documented in the text,
c   computer solution of ordinary differential equations:  the initial
c   value problem  by l. f. shampine and m. k. gordon.
c
      halfu = 0.5
   50 temp1 = 1.0 + halfu
      if(temp1 .le. 1.0) go to 100
      halfu = 0.5*halfu
      go to 50
  100 u = 2.0*halfu
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmgpsm  11040/bald89/weqvmom DEQVMOM
c
      subroutine vmgpsm

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

c---warning:
c---  the radial integrals calculated here are approximate in the
c---  sense that they are done over the sub-intervel [xl0,1.], rather
c---  than the full interval [0.,1.]. the trapezoidal rule is used.

      dimension wb(mm1),wbd(mm1),gf(mm1),wtrap(mm1)
      dimension pwf(mm1),ppwf(mm1),pmid(mm1t),ppmid(mm1t),
     1   pwmid(mm1t),ppwmid(mm1t)
      dimension agspw(mm1),agsppw(mm1),
     1   gsqpw(mtheta),gsqppw(mtheta),pw(mtheta),ppw(mtheta)
      dimension ppw1(mtheta),ppw2(mtheta),ppw3(mtheta)
      dimension prot(mtheta),pprot(mtheta),rotfac(mtheta)
      dimension rotgsq(mtheta),wfcngs(mtheta),bpsq(mtheta),
     1          btsq(mtheta),bsqinv(mtheta)
      dimension rbfac(mm1),rbfold(mm1),rote(mm1),gr2(mtheta),agr2(mm1)
      dimension ftmp(mm1)
      dimension agr0(mm1),gr0(mtheta)
      data nprint /0/
      data rbfold/mm1*1./
      data amass/1.66e-24/, ev/6.242e+11/, everg/1.6022e-12/
c
      data pi/3.14159/, twopi/6.28318/, cvel/3.0e+10/, stata/3.0e+09/
c
c
c  get the metric elements and various averages
c
c
c
      do 10 j=1,m1
        call vmpfun(xp(j),pxp(j),ppxp(j))
        call vmaifn(xp(j),aispl(j),aippl(j))
   10 continue

      do 200 j=1,m1
c
        call vmwrot(xp(j),w1,w2)
c
        r0m=c(j,1)
        r1m=-xp(j)
        r2m=c(j,5)
        r0mx=c(j,2)
        r1mx=-1.0
        r2mx=c(j,6)
c
        z1m=xp(j)*c(j,3)
        z2=c(j,3)*c(j,5)
        z1x=c(j,3)+xp(j)*c(j,4)
        z2x=c(j,3)*r2mx+c(j,4)*r2m
c
        do 100 i=1,ntheta
          r(i)=r0m+r1m*c1th(i)+r2m*c2th(i)
          rx(i)=r0mx+r1mx*c1th(i)+r2mx*c2th(i)
          rt(i)=-(r1m*s1th(i)+2.0*r2m*s2th(i))
          rxt(i)=-(r1mx*s1th(i)+2.0*r2mx*s2th(i))
          rtt(i)=-(r1m*c1th(i)+4.0*r2m*c2th(i))
c
          z(i)=z1m*s1th(i)+z2*s2th(i)
          zx(i)=z1x*s1th(i)+z2x*s2th(i)
          zt(i)=z1m*c1th(i)+2.0*z2*c2th(i)
          zxt(i)=z1x*c1th(i)+2.0*z2x*c2th(i)
          ztt(i)=-(z1m*s1th(i)+4.0*z2*s2th(i))
c
          tau(i)=rt(i)*zx(i)-rx(i)*zt(i)
          gsqrt(i)=r(i)*tau(i)
          gtt(i)=rt(i)**2+zt(i)**2
          gttgsq(i)=gtt(i)/gsqrt(i)
          tauor(i)=tau(i)/r(i)
          gxt=rx(i)*rt(i)+zx(i)*zt(i)
          gxx=rx(i)**2+zx(i)**2
          grxx(i)=gtt(i)/(gxx*gtt(i)-gxt**2)*wfcn(i)
          r2i(i)=wfcn(i)/r(i)**2
          gr0(i)=gsqrt(i)*r(i)
          gr2(i)=gsqrt(i)*r(i)**2
c
          rarg = (r(i)/argeom)**2
          rotfac(i)=w1**rarg
          rotgsq(i)=rotfac(i)*gsqrt(i)
          prot(i)=pxp(j)*rotfac(i)
          pprot(i)=w1**rarg  * (ppxp(j)+ pxp(j)*w2*rarg)
          gsqpw(i)=gsqrt(i)*prot(i)
          gsqppw(i)=gsqrt(i)*pprot(i)
  100   continue
c
        agspw(j)=sdot(ntheta,wfcn,1,gsqpw,1)/cnorm
        agsppw(j)=sdot(ntheta,wfcn,1,gsqppw,1)/cnorm
        agsqrt(j)=sdot(ntheta,wfcn,1,gsqrt,1)/cnorm
        arif(j)=sdot(ntheta,tau,1,wfcn,1)/cnorm/agsqrt(j)
        agttgs(j)=sdot(ntheta,wfcn,1,gttgsq,1)/cnorm
        atauor(j)=sdot(ntheta,wfcn,1,tauor,1)/cnorm
        ar2if(j)=sdot(ntheta,gsqrt,1,r2i,1)/cnorm/agsqrt(j)
        agrho(j)=sdot(ntheta,gsqrt,1,grxx,1)/cnorm/agsqrt(j)
        agr0(j)=sdot(ntheta,wfcn,1,gr0,1)/cnorm/agsqrt(j)     !<r0>flux avg
        r0avg(j) = agr0(j) * a0
        agr2(j)=sdot(ntheta,wfcn,1,gr2,1)/cnorm
c
        pwf(j)=agspw(j)/agsqrt(j)
        ppwf(j)=agsppw(j)/agsqrt(j)
c
        ffp(j)=-(betap0*agsqrt(j)*ppwf(j)+2.0*aispl(j)*aippl(j)
     1    /agttgs(j))/2.0/q0**2/argeom**4/atauor(j)
        ajtm(j)=-agttgs(j)/aispl(j)*(betat0/argeom*pprot(1)*r(1)
     1    +ffp(j)*argeom/r(1))*ajt0
        ajtm(m1+j)=-agttgs(j)/aispl(j)*(betat0/argeom*pprot(ntheta)
     1    *r(ntheta)+ffp(j)*argeom/r(ntheta))*ajt0
        ajorf(j)=-agttgs(j)/aispl(j)*(betat0/argeom*ppwf(j)+ffp(j)
     1    *argeom*ar2if(j))*ajt0/arif(j)
        rmid(j)=r(1)
        rmid(m1+j)=r(ntheta)
c
        rarg1 = (r(1)/argeom) ** 2
        rargt = (r(ntheta)/argeom) ** 2
c        pmid(m1+1-j)=pdata(j)
c        pmid(m1+j)=pdata(j)
c        ppmid(m1+1-j)=ppdata(j)
c        ppmid(m1+j)=ppdata(j)
c        pwmid(m1+1-j)=pxp(j)*w1**rarg1
c        pwmid(m1+j)=pxp(j)*w1**rargt
c        ppwmid(m1+1-j)=w1**rarg1 *
c     1          (ppxp(j) +pxp(j)*w2*rarg1)
        ppxp(j)=w1**rargt * (ppxp(j) +pxp(j)*w2*rargt)
c
        rotbar(j)=sdot(ntheta,wfcn,1,rotgsq,1)/agsqrt(j)
c
  200 continue
c
c
c--  flow beta
210     continue
c      if(irot.le.1) go to 270
c      do 220 j=1,m1
c        if(irot.eq.2) wrot=wstar(j)/r0
c        if(irot.eq.3.or.irot.eq.4) wrot=wtipb(j)
c        rote(j)=eni(j)*wrot**2*agr2(j)
c220   continue
c      volmet=4.*pi**2*a0**3*sdot(m1,wm1,1,agsqrt,1)/cm1
c      pflow=2.*pi**2*a0**3*awp*amass*sdot(m1,wm1,1,rote,1)/cm1
c      btflow=pflow*8.*pi/(bt*1.e3)**2/volmet/100.
270   continue
c
c  perform integrations
c
      sumf2=0.0
      sumg=agsqrt(1)/2.0*xl0/dx
      sumpsi=aispl(1)/agttgs(1)/2.0*xl0/dx
      sump=pxp(1)*agsqrt(1)/2.0*xl0/dx
cbate      sf=sump/pxp(1)/p0
c      sumpe=sf*presne(1)
c      sumpp=sf*pres(1)
      psi(1)=sumpsi*dx
c
c      call splaan(nd,ro,presne,ps1,ps2,ps3)
c      do 291 i=1,m1
c      pxsne(i)=seval(nd,xp(i),ro,presne,ps1,ps2,ps3)
c291   continue
c
c      call splaan(nd,ro,pres,ps1,ps2,ps3)
c      do 293 i=1,m1
c      pxs(i)=seval(nd,xp(i),ro,pres,ps1,ps2,ps3)
c293   continue
c
      do 300 i=1,m1-1
        sumg=sumg+(agsqrt(i)+agsqrt(i+1))/2.0
        sumf2=sumf2+(ffp(m1-i)+ffp(m1-i+1))/2.0
        sump=sump+(pxp(i)*agsqrt(i)+pxp(i+1)*agsqrt(i+1))/2.0
c        sumpe=sumpe+(pxsne(i)*agsqrt(i)+pxsne(i+1)*agsqrt(i+1))/
c     1 (2.0*p0)
c        sumpp=sumpp+(pxs(i)*agsqrt(i)+pxs(i+1)*agsqrt(i+1))/
c     1 (2.0*p0)
        sumpsi=sumpsi+(aispl(i)/agttgs(i)+aispl(i+1)/agttgs(i+1))/2.0
c
        fpol(m1-i)=sqrt(1.0-sumf2*2.0*dx)
        smq(m1-i)=fpol(m1-i)/aispl(m1-i)*q0*argeom**2*
     1    agttgs(m1-i)*atauor(m1-i)
        psi(i+1)=sumpsi*dx
  300 continue
c
      psic0=4.0*pi*a0*ai0*stata/cvel*sumpsi*dx
      fpol(m1)=1.0
      smq(m1)=q0*argeom**2*agttgs(m1)*atauor(m1)
      betapb=sump/sumg*betap0*agsqrt(m1)*agttgs(m1)
      betatb=sump/sumg*pi*8.0/bt0**2*p0*everg
c      bpbp=sumpp/sumg*betap0*agsqrt(m1)*agttgs(m1)
c      bpbe=sumpe/sumg*betap0*agsqrt(m1)*agttgs(m1)
c...must add contributions to pressure from each beam.
c      do 330 nbeam=mbl,mbu
c      sumpb=sf*presb(1,nbeam)
c      call splaan(nd,ro,presb(1,nbeam),ps1,ps2,ps3)
c      do 310 i=1,m1
c      pxsb(i)=seval(nd,xp(i),ro,presb(1,nbeam),ps1,ps2,ps3)
c310   continue
c      do 320 i=1,m1-1
c      sumpb=sumpb+(pxsb(i)*agsqrt(i)+pxsb(i+1)*agsqrt(i+1))/
c     1 (2.0*p0)
c320   continue
c      avbpb(nbeam)=sumpb/sumg*betap0*agsqrt(m1)*agttgs(m1)/cnorm**2
c330   continue
c
c
c---  toroidal flux
      ftmp(1) = 0.
      ftmp(2) =  fpol(1) * atauor(1)
      call vmwgtt (2,xl0,wtrap,ctrap)
      flxtor(1) = 4.0 * pi * f0 * a0 / cvel * stata *
     1            sdot(2,ftmp,1,wtrap,1) / ctrap
      do 60 i=1,m1
      ftmp(i) = fpol(i) * atauor(i)
   60 continue
      do 62 i = 2,m1
      call vmwgtt (i,xp(i)-xl0,wtrap,ctrap)
      flxtor(i) = 4.0 * pi * f0 * a0 / cvel * stata *
     1            sdot(i,ftmp,1,wtrap,1) / ctrap + flxtor(1)
62    continue
      phic0 = flxtor(m1)
      do 65 i=1,m1
      flxtor(i) = flxtor(i) / phic0
   65 continue

c---  volume terms
      do 400 i=1,m1
        psi(i)=psi(i)/psi(m1)
        dvdx(i)=agsqrt(i)*a0*twopi**2
  400 continue
      dvm(1)=dvdx(1)*a0**2*(xp(2)+xp(1))/2.
      do 405 i=2,m1-1
        dvm(i)=dvdx(i)*a0**2*(xp(i+1)-xp(i-1))/2.
405   continue
      dvm(m1)=dvdx(m1)*a0**2*(xp(m1)-xp(m1-1))/2.
c
c   <b0**2/b**2>
        do 4080 j=1,m1
        r0m=c(j,1)
        r1m=-xp(j)
        r2m=c(j,5)
        r0mx=c(j,2)
        r1mx=-1.0
        r2mx=c(j,6)
c
        z1m=xp(j)*c(j,3)
        z2=c(j,3)*c(j,5)
        z1x=c(j,3)+xp(j)*c(j,4)
        z2x=c(j,3)*r2mx+c(j,4)*r2m
c
        bt00=(2.*fpol(j)*f0*stata/a0/cvel)**2
        bp00=(2.*aispl(j)*ai0*stata/(a0*cvel*agttgs(j)))**2
        do 4050 i=1,ntheta
c
          r(i)=r0m+r1m*c1th(i)+r2m*c2th(i)
          rx(i)=r0mx+r1mx*c1th(i)+r2mx*c2th(i)
          rt(i)=-(r1m*s1th(i)+2.0*r2m*s2th(i))
          rxt(i)=-(r1mx*s1th(i)+2.0*r2mx*s2th(i))
          rtt(i)=-(r1m*c1th(i)+4.0*r2m*c2th(i))
c
          z(i)=z1m*s1th(i)+z2*s2th(i)
          zx(i)=z1x*s1th(i)+z2x*s2th(i)
          zt(i)=z1m*c1th(i)+2.0*z2*c2th(i)
          zxt(i)=z1x*c1th(i)+2.0*z2x*c2th(i)
          ztt(i)=-(z1m*s1th(i)+4.0*z2*s2th(i))
c
          tau(i)=rt(i)*zx(i)-rx(i)*zt(i)
          gsqrt(i)=r(i)*tau(i)
          gtt(i)=rt(i)**2+zt(i)**2
          wfcngs(i)=wfcn(i)*gsqrt(i)
c
          btsq(i)=bt00/r(i)**2
          bpsq(i)=bp00*gtt(i)/gsqrt(i)**2
          bmod(i,j)=sqrt(btsq(i)+bpsq(i))
          bsqinv(i)=bt0**2/(btsq(i)+bpsq(i))
          bpoldl(i,j)=sqrt(bpsq(i))
4050      continue
          bsqinf(j)=sdot(ntheta,bsqinv,1,wfcngs,1)/cnorm/agsqrt(j)
          epsl=xp(j)*a0/r0
          r0p=c(j,2)+c(j,6)
          chb1(j)=(1.+1.5*(epsl**2+epsl*r0p)+0.375*epsl**3*r0p)/
     1              (1.+0.5*epsl*r0p)
          chb2(j)=sqrt(1.-epsl**2)*(1.+.5*epsl*r0p)/
     1          (1.+r0p/epsl*(sqrt(1.-epsl**2)-1.))
          fchhi(j)=0.5*(chb1(j)-chb2(j))/sqrt(epsl)/(epsl*sqrt(epsl))
4080      continue
c
c
c   internal inductance "li"
c
      do 410 i=1,m1
      wb(i)=aispl(i)**2/agttgs(i)
410   continue
      do 420 i=1,m1
      ii=i
      a1=sdot(ii,wm1,1,wb,1)
      a2=sdot(ii,wm1,1,agsqrt,1)
      fli(i)=agsqrt(i)*agttgs(i)*a1/(a2*aispl(i)**2)
420   continue
c
c   diamagnetic integral
c   paramagnetic ==> +
c   diamagnetic  ==> -
      do 460 i=1,m1
      wbd(i)=-(bt0*r0-f0*stata/(cvel/2.)*fpol(i))
460   wb(i)=wbd(i)*atauor(i)
      dmgphi=twopi*a0*sdot(m1,wb,1,wm1,1)/cm1/stata
c
c
c
c   compute aj integral  - needs work
      do 500 i=1,m1
      gf(i)=agsqrt(i)*arif(i)*wm1(i)*cnorm
500   continue
      cumom(1)=0.
      do 510 i=2,m1
      ii=i
      sprod=sdot(ii,gf,1,ajorf,1)
      cumom(i)=twopi*sprod/cm1*a0**2
510   continue
c
      if(nprint.eq.0) go to 730
c      write (ndebug,550) betapb,betatb,bpbp,bpbe
c  550 format (//,'   betapb   = ',f10.3,/,2x,'   betatb   = ',
c     1  f10.3,/,
c     1  '   bpbp= ',f5.3,'  bpbe=',f5.3)
      write (ndebug,600)
  600 format (//,'      x     ','   vp(x)    ','     q      ',
     1  '    jflux   ','    psi     ','  ffp(x)   ',
     1  '   jmid0    ','   jmidpi   ','    fpol    ')
      do 650 i=1,m1
        write (ndebug,660) xp(i),dvdx(i),smq(i),ajorf(i),psi(i),
     1    ffp(i),ajtm(i),ajtm(m1+i),fpol(i)
  650 continue
  660 format (1x,f10.3,1x,8(1x,e10.3,1x))
      write (ndebug,670)
  670 format (//,'      x     ','   agsqrt    ','   agttgs    ',
     1  '    atauor   ','   ar2if    ','  arif   ',
     1  '   agrho    ','   dvm    ','             ','    ')
      do 680 i=1,m1
        write (ndebug,675) xp(i),agsqrt(i),agttgs(i),atauor(i),
     1  ar2if(i),arif(i),agrho(i),dvm(i)
  675 format (1x,f10.3,1x,7(1x,e10.3,1x))
  680 continue
      write (ndebug,690)
  690 format (//,'      x     ','   bsqinf    ',
     1  '    chb1   ','   chb2    ','  fchhi   ')
      do 695 i=1,m1
         write (ndebug,675) xp(i),bsqinf(i),chb1(i),
     1      chb2(i),fchhi(i)
  695 continue
c
c
      write(ndebug,700) (xp(j),j=4,m1,3)
700   format(//,1x,25(1h*),' b poloidal from i(x) ',25(1h*),/,
     1  ' theta',t65,'xp'/
     1  t25,f5.3,t37,f5.3,t49,f5.3,t61,f5.3,t73,f5.3,t85,f5.3,
     2  t97,f5.3,t109,f5.3/)
      do 720 i=1,ntheta
        ang=180.*(i-1)/float(ntheta-1)
        write(ndebug,710) ang,(bpoldl(i,j),j=4,m1,3)
710     format(1x,f6.2,t20,1p8e12.3)
720   continue
c
c
  730 continue
c
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@vmprtq  11040/bald89/weqvmom DEQVMOM
c  rgb 05-jun-89 replace mm1 by m1 in do loops and midpoints
c
      subroutine vmprtq
c***********************************************************************
c*****vmprtq writes an output file summarizing the results.            *
c***********************************************************************

      include 'tvmcom0.cmm'
      include 'tvmcom1.cmm'
      include 'tvmcom2.cmm'

      dimension y(mneq), yp(mneq), cpp(mm1,mnmoms)
      pi = acos(-1.0)
      p0pascal = p0 * 1.6022e-12 / 10.
      phic0t = phic0*1.e-4*1.e-4
      write (ndebug,99999)
      write (ndebug,99998) r0, a0, xe0, xe1, xd0, bt0, p0,
     >  p0pascal,pa, pb
      write (ndebug,99997) ai0, aia, vtg1
      write (ndebug,99996) kmhd, hh0, xl0, parerr, error
      write (ndebug,99995) betap0, f0, q0, psic0, phic0, phic0t,
     >  e0, e1, d0
      write (ndebug,99984) betapb, betatb
      write (ndebug,89984) dmgphi
      write (ndebug,99983) va0,va1,va2,va3,irot
      do 20 i=1,m1
           do 10 j=1,mneq
                y(j) = c(i,j)
   10      continue
           x = xp(i)
           call vmaux(yp, y, x, param,ifail)
      if ( ifail .gt. 1 ) return
           cpp(i,1) = yp(2)
           cpp(i,2) = yp(4)
           cpp(i,3) = yp(6)
   20 continue
      write (ndebug,99994)
      do 30 i=1,m1
           write (ndebug,99993) i, xp(i), pxp(i), ppxp(i), aispl(i),
     *          aippl(i)
   30 continue
      write (ndebug,99992)
      do 40 i=1,m1
           write (ndebug,99991) i, xp(i), c(i,1), c(i,2),
     *          cpp(i,1), c(i,3), c(i,4), cpp(i,2), c(i,5),
     *          c(i,6), cpp(i,3)
   40 continue
      write (ndebug,99990)
      do 50 i=1,m1
           write (ndebug,99989) i, xp(i), smq(i), ajorf(i),
     *          psi(i), fpol(i), ffp(i), ajtm(i), ajtm(m1+i)
   50 continue
      write (ndebug,99988)
      do 60 i=1,m1
           write (ndebug,99989) i, xp(i), agsqrt(i), agttgs(i),
     *          atauor(i), agrho(i), arif(i), ar2if(i), flxtor(i)
   60 continue
      write (ndebug,89988)
      do 65 i=1,m1
           write (ndebug,99989) i, xp(i), fli(i), bsqinf(i)
   65 continue
      write (ndebug,99987)
      write (ndebug,99986)
      do 70 i=1,nrho
           write (ndebug,99985) i, rho(i), pres(i), curt(i),
     *          shift(i), elong(i), triang(i)
   70 continue

c
      write(ndebug,700) (xp(j),j=4,m1,3)
700   format(//,1x,25(1h*),' b poloidal from i(x) ',25(1h*),/,
     1  ' theta',t65,'xp'/
     1  t25,f5.3,t37,f5.3,t49,f5.3,t61,f5.3,t73,f5.3,t85,f5.3,
     2  t97,f5.3,t109,f5.3/)
      do 720 i=1,ntheta
        ang=180.*(i-1)/float(ntheta-1)
        write(ndebug,710) ang,(bpoldl(i,j),j=4,m1,3)
710     format(1x,f6.2,t20,1p8e12.3)
720   continue

      return


99999 format (1h1, /)
99998 format (25x, 20('*'), ' plasma  parameters ', 25('*'),
     *     //, 40x, 'r0 = ', f5.1, t60, 'cm', /, 40x, 'a0 = ',
     *     f5.1, t60, 'cm', /, 40x, 'xe0 = ', f5.2, t60,
     *     'geometric',/, 40x, 'xe1 = ', f5.2, t60,'geometric',
     *     /, 40x, 'xd0 =', f5.2, t60, 'geometric', /,
     *     40x, 'bt0 = ', f7.1, t60, 'gauss', /,
     >     40x, 'p0 = ',1pe9.2, t60, 'ev/cm**3', /,
     >     40x, 'p0 = ',1pe9.2, t60, 'pascals', /,
     >     40x, 'pa = ', 0pf4.2, /,
     *     40x, 'pb = ', f4.2)
99997 format (40x, 'ai0am = ', 1pe8.2, t60, 'amps', /, 40x,
     *     'aia = ', 0pf4.2, /, 40x, 'time = ', 0pf4.2, //)
99996 format (25x, 20('*'), ' d02agEb  parameters ', 20('*'),
     *     //, 40x, 'kmhd   = ', i4, /, 40x, 'hh0 = ', f7.5, /,
     *     40x, 'xl0 = ', f7.5, /, 40x, 'parerr = ', 1p6e10.2,
     *     /, 40x, 'error  = ', 1p6e10.2, ///)
99995 format (25x, 20('*'), ' internal  parameters ', 20('*'),
     *     //, 40x, 'betap0 = ', f6.2, /, 40x, 'f0 = ',
     *     1pe10.2, ' amps', /, 40x, 'q0 = ', 1pe10.2, /, 40x,
     *     'psic0 = ', 1pe10.2, ' gauss-cm**2 (poloidal)', /,
     *     40x,'phic0 = ',1pe10.2, ' gauss-cm**2 (toroidal)', /,
     *     40x,'phic0 = ',1pe10.2, ' tesla       (toroidal)', /,
     *     40x,'e0 = ', 0pf8.6, t60, 'moments', /, 40x, 'e1 = ',
     *     f8.6, t60, 'moments', /, 40x, 'd0 = ', f9.7, t60,
     *     'moments')
99994 format ('1', 25x, 20('*'), ' moments(internal) ',
     *     20('*'), //, '  i', 3x, 'x', t15, 'p-hat', t27,
     *     'pp-hat', t39, 'i-hat', t50, 'ip-hat')
99993 format (i3, f6.2, 1p4e12.3)
99992 format (//, '  i', 3x, 'x', t16, 'r0', t28, 'r0p', t38,
     *     'r0pp', t53, 'e0', t65, 'e0p', t76, 'e0pp', t89,
     *     'd0', t101, 'd0p', t112, 'd0pp')
99991 format (i3, f6.2, 1p9e12.3)
99990 format (//, '  i', 3x, 'x', t16, 'q', t27, 'jflux', t39,
     *     'psi', t51, 'f', t62, 'ffp', t72, 'jmdpl(in)', t84,
     *     'jmdpl(out)', /, t25, '(statamps)', t72,
     *     '(statamps)', t84, '(statamps)')
99989 format (i3, f6.2, 1p8e12.3)
99988 format (//, '  i', 3x, 'x', t14, ' gsqrt', t26, ' gttgs',
     *     t38, ' gsgpp', t51, ' grho', t62, 'arif', t74,
     *     'ar2if', t86, 'toroidal flux')
89988 format (//, '  i', 3x, 'x', t14, ' li', t26, ' <b0**2/b**2>')
99987 format ('1', 25x, 20('*'), ' moments(external) ',
     *     20('*'), ////)
99986 format (//, '  i', 3x, 'r', t17, 'p', t29, 'ai', t39,
     *     'shift', t51, 'elong', t63, 'triang', /, t6, '(cm)',
     *     t13, '(ev/cm**3)', t27, '(amps)', t39, '(cm)', /)
99985 format (i3, f6.2, 1p5e12.3)
99984 format (//,25x, 20('*'), ' volume averaged betas ', 20('*'),
     *     //, 40x, 'beta-p =  ', f6.2, /, 40x, 'beta-t = ',
     *     2pf6.2, ' (%) ', / )
89984 format (//,25x, 20('*'), ' diamagnetic integral ', 20('*'),
     *     /, 40x, 'dmgphi = ', 1pe10.3, / )
99983 format (//,25x, 20('*'), 'rotation parameters', 20('*'),//,
     *     40x, 'va0 = ',1pe15.7, /, 40x, 'va1 = ', 1pe15.7, /,
     *     40x, 'va2 = ',1pe15.7, /, 40x, 'va3 = ', 1pe15.7, /,
     *     40x, 'irot = ',i3)
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@d02agEb  11040/bald89/weqvmom DEQVMOM
c
c***********************************************************************
c---beware...the version of {d02agEb et al} . has been modified to work
c---with vmoms, thus relying on the structure of the vector defined
c---in that code. the places where it has been patched in such a way
c---are commented as "vmoms patch". other efficiencies are denoted by
c---the comment "efficiency patch"
c***********************************************************************

      subroutine d02agEb(h, error, parerr, param, c, n, n1, m1,
     *     vmaux, vmbcux, vmraux, vmprsl, mat, copy, wspace, wspac1,
     *     wspac2, ifail)
c     nag copyright 1975
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 4.5 revised
c     additional comments by r. wieland // ornl // jan.,1981.
c     patched on 8-dec-1987 by r.m. wieland // ppl //
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   the parameters are defined as follows:
c
c     h           input estimate of required integration step length.
c                 changed on return to provide last step length used.
c
c     error(n)    input a real array used as:
c                 a) a local error bound for integration
c                 b) a convergence tolerance on y at the matching point.
c                    on return, the final error.
c
c     parerr(n1)  input a real array used as:
c                 a) a convergence tolerance for the parameters (p)
c                 b) used to approximate delta-p in estimating the
c                    jacobian. on return, the final error.
c
c     param(n1)   input starting values for p;
c                 on return, the corrected values
c
c     c(m1,n)     the solution y(i;j) of the j-th component of y
c                 evaluated at x(i) is returned in c (i,j); the x(i)
c                 spacing is determined by m1, below.
c
c     n           number of ode's
c
c     n1          number of parameters
c
c     m1          the final solution is calculated at m1 equidistant
c                 points.
c
c     vmaux         a user supplied subroutine:
c                 subroutine vmaux (f,y,x,param)
c                 where x,y(n),param(n1) are used to evalute the
c                 n derivatives f(n) at x.
c
c     vmbcux       user supplied subroutine:
c                 subroutine vmbcux(g,g1,param)
c                 where param is used (if necessary) to evaluate y(n)
c                 at the endpoints x and x1 and return them in g(n)
c                 and g1(n).
c
c     vmraux       user supplied subroutine:
c                 subroutine vmraux(x,x1,r,param)
c                 where param is used (if necessary) to evaluate the
c                 endpoints x and x1, and the matching point r.
c
c     vmprsl       a user supplied subroutine:
c                 subroutine vmprsl(param,res,n1,err)
c                 called at each newton iteration; can be used to output
c                 any of the parameters of interest; err(n) are the
c                 errors at r in each component y, and res is the
c                 sum of the squares of these errors.
c
c     mat(n1,n1)   a real work array
c     copy(n1,n1)  a real work array
c     wspace(n,9) a real work array
c     wspac1(n)   a real work array
c     wspac2(n)   a real work array
c
c     ifail       an error flag; on input, enter 0, on output:
c       ifail = 0   normal return
c             = 1   n1 > n
c             = 2   integration failed to converge while calc. jacobian
c             = 3   the condition x <= r <= x1 does not hold.
c             = 4   integration failed to converge
c             = 5   jacobian is singular
c             = 6   after 3/no 6 attempts at halving delta-p, the newton-
c                   raphson loop still does not yield a diminishing
c                   residual.
c             = 7   after 12 n-r iterations, there is still not
c                   sufficient convergence.
c
c    internal flag :
c    vrbose (logical)   .true. to generate debug output on ndebug
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer count, count1, em, one, ct, p01aaEb, n, n1, m1,
     *     ifail, m, i, k, itest, j
      character *8 srname
      real mat(n1,n1), h, error(n), parerr(n1), param(n1),
     *     c(m1,n), copy(n1,n1), wspace(n,9), wspac1(n),
     *     wspac2(n), eps, h0, x, x1, r, dum, resid, d, pert,
     *     oldres, dist, c1, x02aaE
cccccccccccccccccccccccccccccccccccccc
      include 'tvmcom2.cmm'
cccccccccccccccccccccccccccccccccccccc
      external vmaux, vmbcux, vmraux, vmprsl
cbate      external vmaux, %&&&
      logical vrbose
      data vrbose /.true./
      data srname / ' d02agEb'/
c
c     solves a general boundary value
c     problem for n differential equations
c     in n1 parameters using a shooting
c     and matching technique.  eps is the
c     largest real variable such that 1+eps=1
c     all implicitly declared reals may be used double-length
c     the array copy is redundant
c
c     get eps from vmcom2
      eps = epsmch
      m = m1 - 1
      if (n1.le.n) go to 10
c
c     *** ifail = 1 ***
      em = 1
      go to 480
c
c     *** set newton iteration counter ***
   10 count = 0
      h0 = h
      one = 1
      em = -1
c
c     forms the residuals at the
c     matching point
c
c     *** get left-hand pt., right-hand pt., matching pt. ***
   20 call vmraux(x, x1, r, param)
      if ((x-r)*(x1-r).le.0.0) go to 30
c
c     *** ifail = 3 ***
      em = 3
      go to 480
c
   30 if (h0*(x1-x).lt.0.0) h0 = -h0
c     ### g(i)  ==> w1            ###
c     ### g1(i) ==> w2 ==> w(i,8) ###
      call vmbcux(wspac1, wspac2, param)
      h = h0
      do 40 i=1,n
           wspace(i,8) = wspac2(i)
   40 continue
      i = 1
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99999)
      nfl = 0
cccccccccccccccccccccccccccccccccccccc
c
c     *** integrate from x ==> r // on return: y(r-) ==> w1(i) ***
      call agzd02b(x, wspac1, error, one, n, n1, i, r-x, h, vmaux,
     *     wspace, wspac2, param)
      if (i.eq.0) go to 50
c
c     *** ifail = 4 ***
      em = 4
      go to 480
c
c     ###  -y(r-) ==> w(i,8) // g1(i) ==> w1(i) ###
   50 do 60 i=1,n
           dum = wspace(i,8)
           wspace(i,8) = -wspac1(i)
           wspac1(i) = dum
   60 continue
c
c     *** integrate from x1 ==> r // on return: yr(+) ==> w1(i) ***
      h = -h0
      i = 1
      call agzd02b(x1, wspac1, error, one, n, n1, i, r-x1, h,
     *     vmaux, wspace, wspac2, param)
      if (i.eq.0) go to 70
c
c     *** ifail = 4 ***
      em = 4
      go to 480
c
   70 resid = 0.0
      ct = 0
c
c     *** form error residuals at r:                  ***
c     ###   r(i) = y(r+) - y(r-) ==> w(i,8) == s(i;p) ###
c     ###   s(i;p) ==> w1(i)                          ###
      do 80 i=1,n1
           d = wspac1(i)
           dum = wspace(i,8)
           wspace(i,8) = d + dum
           d = 1.0 + abs(dum) + abs(d)
           dum = wspace(i,8)
           if (abs(dum).lt.error(i)*d) ct = ct + 1
           wspac1(i) = dum
           resid = resid + dum*dum
   80 continue
      call vmprsl(param, resid, n1, wspac1)
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99998) ct, nfl
cbate      if (vrbose) write (ndebug,99997) resid
cccccccccccccccccccccccccccccccccccccc
c     *** em = -1 first time thru only ***
      if (em.ne.-1) go to 320
   90 count = count + 1
      if (count.ne.12) go to 100
c
c     *** ifail = 7 ***
      em = 7
      go to 480
c
c     forms the jacobian by numerical
c     differentiation
c     *** del[p(k)] ==> pert ***
c     *** p(k) + del[p(k)] ==> p(k+) ***
  100 continue
cccccccccccccccccccccccccccccccccccccc
      nfl = 0
cccccccccccccccccccccccccccccccccccccc
      do 180 k=1,n1
c*********start of vmoms patch **********************************************
c*********attenberger-houlberg patch :: installed by r. wieland on 8-dec-1987
c*********for the vmoms application, the jacobian for p2,p4 & p6 is already
c*********known. this patch results in a 20 0mprovement in execution time
c*********for a standard case. this is an expensive calculation.
           if (mod(k,2).eq.0) then
             do 105 i=1,n1
               if(i.eq.k) then
                 mat(i,k)= -1.0
               else
                 mat(i,k)= 0.0
               end if
105          continue
             go to 180
           end if
c**********end of vmoms patch ***********************************************
           pert = 10.0*parerr(k)*(1.0+abs(param(k)))
           param(k) = pert + param(k)
           call vmraux(x, x1, r, param)
           if ((x-r)*(x1-r).le.0.0) go to 110
c
c     *** ifail = 3 ***
           em = 3
           go to 480
c
  110      if (h0*(x1-x).lt.0.0) h0 = -h0
           h = h0
           call vmbcux(wspac1, wspac2, param)
c     ### g1(i)[p(k+)] ==> w(i,7) ###
           do 120 i=1,n
                wspace(i,7) = wspac2(i)
  120      continue
           i = 1
c
c     *** integrate from x ==> r using p(k+) ***
           call agzd02b(x, wspac1, error, one, n, n1, i, r-x, h,
     *          vmaux, wspace, wspac2, param)
           if (i.eq.0) go to 130
c
c     *** ifail = 2 ***
           em = 2
           go to 480
c
c     ### m(i,k) = y(r-;p(k+)) ###
  130      do 140 i=1,n1
                mat(i,k) = wspac1(i)
  140      continue
           h = -h0
           i = 1
           do 150 i=1,n
                wspac1(i) = wspace(i,7)
  150      continue
c
c     *** integrate from x1 ==> r using p(k+) ***
           call agzd02b(x1, wspac1, error, one, n, n1, i, r-x1,
     *          h, vmaux, wspace, wspac2, param)
           if (i.eq.0) go to 160
c
c     *** ifail = 2 ***
           em = 2
           go to 480
c
c     ### m(i,k) = (s(p(k+)-s(p(k))/del[p(k)] ###
  160      do 170 i=1,n1
                mat(i,k) = (mat(i,k)-wspac1(i)+wspace(i,8))/pert
                if (abs(mat(i,k)).lt.5.0*eps*abs(wspace(i,8))/
     *               pert) mat(i,k) = 0.0
  170      continue
           param(k) = param(k) - pert
  180 continue
c
c     *** new jacobian ***
      itest = 1
      em = -3
c
c     performs column scaling on the jacobian
c     and forms a triangular decomposition
      do 220 i=1,n1
           d = 0.0
           do 190 j=1,n1
                if (abs(mat(j,i)).gt.d) d = abs(mat(j,i))
  190      continue
           if (d.ne.0.0) go to 200
c
c     *** ifail = 5 ***
           em = 5
           go to 480
c
c     *** normalize m(i,k) ***
  200      do 210 j=1,n1
                mat(j,i) = mat(j,i)/d
  210      continue
           wspace(i,7) = d
  220 continue
      i = 1
c
c     *** lu decomposition of m(i,k) // pivot array: w1(i) ***
      call f03afEb(n1, eps, mat, n1, d, j, wspac1, i)
      if (i.eq.0) go to 230
c
c     *** ifail = 5 ***
      em = 5
      go to 480
c
c     ### pivot array ==> w(i,6) ###
  230 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99996) nfl
cccccccccccccccccccccccccccccccccccccc
      do 240 i=1,n1
           wspace(i,6) = wspac1(i)
  240 continue
c
c     uses a generalised newton raphson
c     technique to solve the nonlinear
c     equations at the matching point
c     *** solve the system :  m(i,k) * del[p(k)] = s(i;p) ***
  250 oldres = resid
c
c     *** set newton - raphson counter ***
      count1 = 0
c
c     ### s(p) ==> w1(i) ###
      do 260 i=1,n1
           wspac1(i) = wspace(i,6)
           wspace(i,1) = wspace(i,8)
  260 continue
c
c     *** the new del[p] = w1 ***
      call f04ajEb(n1, one, mat, n1, wspac1, wspace, n)
c
c     *** normalize del[p] / form new p ***
      do 270 i=1,n1
           wspace(i,1) = wspace(i,1)/wspace(i,7)
           param(i) = param(i) + wspace(i,1)
  270 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99995) param
cccccccccccccccccccccccccccccccccccccc
      if (ct.lt.n1) go to 300
      do 280 i=1,n1
cccccccccccccccccccccccccccccccccccccc
           ipp = i
cccccccccccccccccccccccccccccccccccccc
           if (abs(wspace(i,1)).gt.parerr(i)*(1.0+abs(param(i)))
     *          ) go to 290
  280 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99994)
cccccccccccccccccccccccccccccccccccccc
c
c     *** solution found! ***
      em = -5
      go to 380
  290 continue
cccccccccccccccccccccccccccccccccccccc
      pip = parerr(ipp)*(1.0+abs(param(ipp)))
cbate      if (vrbose) write (ndebug,99993) ipp, wspace(ipp,1), pip
cccccccccccccccccccccccccccccccccccccc
  300 do 310 i=1,n1
           wspace(i,1) = -wspace(i,1)
  310 continue
      go to 20
c
c     ******************************************************
c
  320 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99992) oldres
cccccccccccccccccccccccccccccccccccccc
      if (count1.ne.0) go to 330
      if (resid.ge.oldres/10.0) go to 330
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99991)
cccccccccccccccccccccccccccccccccccccc
c
c     *** go back and form new del[p] using current jacobian, new s ***
      em = -2
      itest = 0
      go to 250
c
c     *** form new jacobian ***
  330 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99990)
cccccccccccccccccccccccccccccccccccccc
      if (resid.lt.oldres) go to 90
c*****vmoms patch :: attenberger & houlberg suggestion
cxxxxxif (count1.ne.3) go to 360
      if (count1.ne.6) go to 360
      if (itest.eq.0) go to 340
c*****end of vmoms patch
c     *** ifail = 6 ***
      em = 6
      go to 480
  340 continue
      do 350 i=1,n1
           param(i) = param(i) + wspace(i,1)
  350 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99989) param
cccccccccccccccccccccccccccccccccccccc
      em = -1
      go to 20
  360 count1 = count1 + 1
      em = -4
c
c     *** scale-down del[p] ***
      do 370 i=1,n1
           wspace(i,1) = wspace(i,1)/2.0
           param(i) = param(i) + wspace(i,1)
  370 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99988) param
cccccccccccccccccccccccccccccccccccccc
      go to 20
c
c     ******************************************************
c
c     calculates the final solution
  380 continue
cccccccccccccccccccccccccccccccccccccc
cbate      if (vrbose) write (ndebug,99987)
cccccccccccccccccccccccccccccccccccccc
      if (m.le.0) go to 480
      call vmraux(x, x1, r, param)
      if ((x-r)*(x1-r).le.0.0) go to 390
      em = 3
      go to 480
  390 if (h0*(x1-x).lt.0.0) h0 = -h0
      h = h0
      call vmbcux(wspac1, wspac2, param)
      do 400 i=1,n
           wspace(i,7) = wspac2(i)
  400 continue
      dist = (x1-x)/float(m)
      j = 1
      c1 = x
      k = 1
c*****vmoms patch :: sea & wh change to preserve edge conditions
c*****--this is only applicable to the vmoms code as currently configured
c*****--08-dec-1987 - old coe has been deleted
  410 i = 1
      call agzd02b(c1, wspac1, error, one, n, n1, i, dist, h,
     *     vmaux, wspace, wspac2, param)
      if (i.eq.0) go to 430
      em = 4
      go to 480
  430 j = j + k
      do 440 i = 1,n
        c(j,i) = wspac1(i)
  440 continue
      if ((r-c1-0.25*dist)*dist.le.0.0) go to 460
      go to 410
  460 call vmbcux(wspac1, wspac2, param)
      do 470 i=1,n
        c(1,i) = wspac1(i) !this is at the "origin" point
  470 continue
c*****end of vmoms patch
  480 if (em.le.0) ifail = 0
      if (em.gt.0) ifail = p01aaEb(ifail,em,srname)
      return
c     end of d02agEb
c99999 format (/' VMOMS equilibrium package'
c     &,' by Lao, Wieland, Houlberg and Hirshman ORNL/TM-7871 1982')
c99998 format (2x, i2, ' y13hs match with ', i6, ' vmaux evals')
c99997 format (' vmprsl: resid = ', 1pe11.2)
c99996 format (' new jacobian with ', i6, ' vmaux evals')
c99995 format (' solve for new p:', /, 5x, 1p6e10.2)
c99994 format (' del p is small enough')
c99993 format ('  del p # ', i1, ' still too large:', 1pe11.3,
c     *     2h >, e11.3)
c99992 format (' compare against old resid- ', 1pe11.3)
c99991 format (' keep the old jac; get new del p using new s')
c99990 format (' new resid is not a factor of 10 better')
c99989 format (' resid is worse & lim=3; restore old p:', /,
c     *     5x, 1p6e10.2)
c99988 format (' resid is worse; back down del p', /, 5x,
c     *     'new p:', /, 5x, 1p6e10.2)
c99987 format (' finished')
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@agzd02b  11040/bald89/weqvmom DEQVMOM
c  rgb 06-aug-90 added ifail to argument list
c
      subroutine agzd02b(x, y, g, t, n, n1, m, h0, h, vmaux,
     *     wspace, wspac1, param)
c
c     uses the logic of nag library routine d02abE
c     eps and dum may be used double length.  eps is the
c     smallest real such that 1+eps>eps.  smax
c     is the largest integer. dum,err and hs maybe declared
c     double precision
c     nag copyright 1975
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 4.5 revised
c     additional comments by r. wieland // ornl // jan.,1981.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     the subroutine parameters are defined as follows:
c
c     x       input inital x value; on return x=x+h0 if successful
c
c     y(n)    input initial values of y(x); on return y(x+h0) if
c             successful
c
c     g(n)    input error bounds on y(n)
c
c     t       input type of error test
c
c     n       number of ode's in the system
c
c     m       input error flag as 0; on return:
c             m = 0   normal return
c               = 1   step length repeatedly halved to < 10**-4 the
c                     initially suggested step length
c               = 2   t .ne. 1 or 2  or 3
c               = 3   the number of steps required exceeds the largest
c                     integer on the machine
c
c     h0      interval of integration
c
c     h       input an estimate of the step length required;
c             on return, the final step length used.
c
c      the remaining parameters are as defined in 'd02agEb'
c
c     on return:
c       wspace(i,2) contains the local error at x0+h for each y(i)
c       wspace(i,5) contains y'(x+h0)
c       wspace(i,3) contains y'(x)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      integer d, s, smax, t, n, n1, m, i, j
      real x, y(n), g(n), h0, h, wspace(n,9), wspac1(n),
     *     param(n1), eps, dum, p, q, hs, err, x02aaE
      include 'tvmcom2.cmm'
      external vmaux
c     *** h0 = interval length // h = step length ***
      eps = epsmch
      smax = i1machb(9)
      m = 0
      dum = max(abs(x),abs(x+h0))
      if (abs(h0).le.eps*dum) return
      dum = h/h0
      if (abs(dum).gt.1.0 .or. h.eq.0.0) go to 10
      dum = abs(h0/h+0.9)
      if (dum.gt.float(smax)) go to 80
      h = h0/aint(dum)
      go to 20
   10 h = h0
c
c     *** t = 1 :: mixed    error test    ***
c     *** t = 2 :: absolute error test    ***
c     *** t = 3 :: relative error test    ***
   20 p = 1.0
      if (t.eq.3) p = 0.0
      q = 1.0
      if (t.eq.2) q = 0.0
      dum = h0/h + 0.1
c     *** s = # of anticipated steps       ***
c     *** hs = minimum step length allowed ***
      s = int(dum)
      hs = 1.0e-4*abs(h)
   30 do 40 i=1,n
           wspace(i,9) = y(i)
   40 continue
c
c     *** integrate from x ==> x+h;                              ***
c     ***   input:  y(x) ==> y // x ==> x;                       ***
c     ***   output: y(x+h) ==> y // x+h ==> x;                   ***
c     ***           local error ==> w(i,2) // y'(x+h) ==> w(i,5) ***
c     ***           y'(x) ==> w(i,3)                             ***
      call agyd02b(y, x, h, n, n1, vmaux, wspace, wspac1, param,ifail)
c
      if ( ifail .gt. 1 ) then
        m = ifail
        return
      endif
c
      d = 0
      do 50 i=1,n
           err = g(i)*(p+q*abs(y(i)))
           if (wspace(i,2).gt.err) go to 60
           if (40.0*wspace(i,2).gt.err) d = 1
   50 continue
      s = s - 1
c
c     *** the only successful way out ***
      if (s.eq.0) return
c
      dum = float(s)/2.0 + 0.1
      j = int(dum)*2
      if (d.ne.0 .or. j.ne.s) go to 30
      h = 2.0*h
      s = s/2
      go to 30
c
c     *** local error too large ***
   60 x = x - h
      do 70 i=1,n
           y(i) = wspace(i,9)
   70 continue
      if (s.gt.smax/2) go to 80
      s = 2*s
      h = 0.5*h
      if (abs(h).gt.hs) go to 30
   80 m = 1
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@agyd02b  11040/bald89/weqvmom DEQVMOM
c  rgb 06-aug-90 added ifail to argument list
c
      subroutine agyd02b(y, x, h, n, n1, vmaux, wspace, wspac1,
     *     param,ifail)
c
c     uses the technique of nag library
c     procedure d02aaE with the
c     specification of vmaux changed
c     all implicitly declared reals may
c     be declared double precision
c     nag copyright 1975
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 4.5 revised
c
c   integrates a system of ode's using merson's meth.
c   cf.  lambert, j.d., computational methods in ordinary differential
c        equations, pp. 130-135, wiley, 1973.
c
      integer n, n1, i
      real y(n), x, h, wspace(n,9), wspac1(n), param(n1), c1,
     *     c2, u, v, w, dum
      external vmaux
c
      call vmaux(wspac1, y, x, param,ifail)
      if ( ifail .gt. 1 ) return
      c1 = 1.0/3.0
      c2 = 1.0/6.0
      do 10 i=1,n
           wspace(i,3) = wspac1(i)
   10 continue
      u = c1*h
      do 20 i=1,n
           wspace(i,4) = y(i)
           y(i) = y(i) + u*wspace(i,3)
   20 continue
      call vmaux(wspac1, y, u+x, param,ifail)
      if ( ifail .gt. 1 ) return
      do 30 i=1,n
           wspace(i,5) = wspac1(i)
   30 continue
      v = h*c2
      do 40 i=1,n
           y(i) = wspace(i,4) + v*(wspace(i,3)+wspace(i,5))
   40 continue
      call vmaux(wspac1, y, u+x, param,ifail)
      if ( ifail .gt. 1 ) return
      do 50 i=1,n
           wspace(i,5) = wspac1(i)
   50 continue
      u = h*0.125
      v = h*0.375
      do 60 i=1,n
           y(i) = wspace(i,4) + wspace(i,3)*u + wspace(i,5)*v
   60 continue
      u = 0.5*h
      v = 1.5*h
      w = 2.0*h
      call vmaux(wspac1, y, u+x, param,ifail)
      if ( ifail .gt. 1 ) return
      do 70 i=1,n
           y(i) = wspace(i,4) + wspace(i,3)*u - wspace(i,5)*v +
     *          wspac1(i)*w
           wspace(i,5) = wspac1(i)
   70 continue
      x = x + h
      call vmaux(wspac1, y, x, param,ifail)
      if ( ifail .gt. 1 ) return
      u = h*c2
      v = 2.0*h*c1
      do 80 i=1,n
           w = wspace(i,4) + u*(wspace(i,3)+wspac1(i)) +
     *          wspace(i,5)*v
           dum = w - y(i)
           wspace(i,2) = 0.2*abs(dum)
           y(i) = w
   80 continue
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@f03afEb  11040/bald89/weqvmom DEQVMOM
c
      subroutine f03afEb(n, eps, a, ia, d1, id, p, ifail)
c     mark 2 release. nag copyright 1972
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 3 revised.
c     mark 4.5 revised
c     sw3 switch revision : r. wieland // ornl //   jan., 1981.
c     .false. gives no extended precision
c     .true. gives extended precision if conditions commented
c     on in aazx03 are fulfilled.
c
c     unsymdet
c     the unsymmetric matrix, a, is stored in the n*n array a(i,j),
c     i=1,n, j=1,n. the decomposition a=lu, where l is a
c     lower triangular matrix and u is a unit upper triangular
c     matrix, is performed and overwritten on a, omitting the unit
c     diagonal of u. a record of any interchanges made to the rows
c     of a is kept in p(i), i=1,n, such that the i-th row and
c     the p(i)-th row were interchanged at the i-th step. the
c     determinant, d1 * 2.0**id, of a is also computed. the
c     subroutine
c     will fail if a, modified by the rounding errors, is singular
c     or almost singular. sets ifail = 0 if successful else ifail =
c     1.
c     1st december 1971
c
      integer isave, ifail, ifail1, i, n, ia, id, k, l, k1, k2,
     *     istart, j, p01aaEb
      character *8 srname
      real y, d2, d1, x, eps, a(ia,n), p(n)
      logical sw3
      data srname / ' f03afEb'/
      data sw3 /.false./
      isave = ifail
      ifail1 = 0
      do 10 i=1,n
           call x03aaEb(a(i,1), ia*n+1-i, a(i,1), ia*n+1-i, n,
     *          ia, ia, 0.0, 0.0, y, d2, sw3, ifail1)
           if (y.le.0.0) go to 100
           p(i) = 1.0/sqrt(y)
   10 continue
      d1 = 1.0
      id = 0
      do 90 k=1,n
           l = k
           x = 0.0
           k1 = k - 1
           k2 = k + 1
           istart = k
           do 20 i=istart,n
                call x03aaEb(a(i,1), n*ia+1-i, a(1,k),
     *               (n-k+1)*ia, k1, ia, 1, -a(i,k), 0.0, y,
     *               d2, sw3, ifail1)
                a(i,k) = -y
                y = abs(y*p(i))
                if (y.le.x) go to 20
                x = y
                l = i
   20      continue
           if (l.eq.k) go to 40
           d1 = -d1
           do 30 j=1,n
                y = a(k,j)
                a(k,j) = a(l,j)
                a(l,j) = y
   30      continue
           p(l) = p(k)
   40      p(k) = l
           d1 = d1*a(k,k)
           if (x.lt.8.0*eps) go to 100
   50      if (abs(d1).lt.1.0) go to 60
           d1 = d1*0.0625
           id = id + 4
           go to 50
   60      if (abs(d1).ge.0.0625) go to 70
           d1 = d1*16.0
           id = id - 4
           go to 60
   70      x = -1.0/a(k,k)
           if (k2.gt.n) go to 90
           do 80 j=k2,n
                call x03aaEb(a(k,1), n*ia+1-k, a(1,j),
     *               (n-j+1)*ia, k1, ia, 1, -a(k,j), 0.0, y,
     *               d2, sw3, ifail1)
                a(k,j) = x*y
   80      continue
   90 continue
      ifail = 0
      return
  100 ifail = p01aaEb(isave,1,srname)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@f04ajEb  11040/bald89/weqvmom DEQVMOM
c
      subroutine f04ajEb(n, ir, a, ia, p, b, ib)
c     mark 2 release. nag copyright 1972
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 4 revised.
c     mark 4.5 revised
c     sw3 switch revision : r. wieland // ornl //   jan., 1981.
c     .false. gives no extended precision
c     .true. gives extended precision if conditions commented
c     on in aazx03 are fulfilled.
c     unsymsol
c     solves ax=b, where a is an unsymmetric matrix and b is an
c     n*ir
c     matrix of ir right-hand sides. the subroutine f04ajEb must by
c     preceded by f03afEb in which l and u are produced in a(i,j),
c     from a, and the record of the interchanges is produced in
c     p(i). ax=b is solved in three steps, interchange the
c     elements of b, ly=b and ux=y. the matrices y and then x are
c     overwritten on b.
c     1st august 1971
c
      integer ifail1, i, n, i1, k, ir, ia, ib, ii
      real x, d1, d2, a(ia,n), p(n), b(ib,ir)
      logical sw3
      data sw3 /.false./
      ifail1 = 0
c     interchanging of elements of b
      do 20 i=1,n
           i1 = p(i) + 0.5
           if (i1.eq.i) go to 20
           do 10 k=1,ir
                x = b(i,k)
                b(i,k) = b(i1,k)
                b(i1,k) = x
   10      continue
   20 continue
      do 50 k=1,ir
c     solution of ly= b
           do 30 i=1,n
                i1 = i - 1
                call x03aaEb(a(i,1), n*ia-i+1, b(1,k),
     *               (ir-k+1)*ib, i1, ia, 1, b(i,k), 0.0, d1,
     *               d2, sw3, ifail1)
                b(i,k) = -d1/a(i,i)
   30      continue
c     solution of ux= y
           b(n,k) = -b(n,k)
           if (n.eq.1) go to 50
           do 40 ii=2,n
                i = n - ii + 1
                i1 = i + 1
                call x03aaEb(a(i,i1), (n-i)*ia-i+1, b(i1,k),
     *               (ir-k+1)*ib-i, n-i, ia, 1, b(i,k), 0.0,
     *               d1, d2, sw3, ifail1)
                b(i,k) = -d1
   40      continue
   50 continue
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@p01aaEb  11040/bald89/weqvmom DEQVMOM
c
      integer function p01aaEb(ifail, error, srname)
c     mark 1 release.  nag copyright 1971
c     mark 3 revised
c     mark 4a revised, ier-45
c     mark 4.5 revised
c     mark 7 revised (dec 1978) .... (apr 1979)
c     nout = ntty :: revised by r. wieland // ornl // jan., 1981.
c     returns the value of error or terminates the program.
c     if a hard failure occurs, this routine calls a fortran auxiliary
c     routine aazp01 which gives a trace, a failure message and halts
c     the program
      integer error, ifail
      character *8 srname
      include 'tvmcom2.cmm'
c     test if no error detected
      if (error.eq.0) go to 20
c     test for soft failure
      if (mod(ifail,10).eq.1) go to 10
c     hard failure
      write (nttyo,99999) srname, error
c     stopping mechanism may also differ
c      call aazp01 (x)
c     stop
c     soft fail
c     test if error messages suppressed
   10 if (mod(ifail/10,10).eq.0) go to 20
      write (nttyo,99999) srname, error
   20 p01aaEb = error
      return
99999 format ('0', 'error detected by nag library routine ',
     *     a8, ' - ifail = ', i5//)
      end
c     auto edit 20 sep 76
c     auto edit 20 sep 76
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@x03aaEb  11040/bald89/weqvmom DEQVMOM
c  pis 12-jun-98 changed declarations double precision --> real
c    and changed casts dble --> real
c
      subroutine x03aaEb(a, isizea, b, isizeb, n, istepa,
     *     istepb, c1, c2, d1, d2, sw, ifail)
c     nag copyright 1975
c     edited by joyce clarke oxford oeg nuclear physics 11th sep 1976
c                   fortran macro version fdia24.tec
c     mark 4.5 release
c
c     calculates the value of a scalar product using basic
c     or additional precision and adds it to a basic or additional
c     precision initial value.
c     x03aaEb calls aazx03 which may be in assembly language.
c
      integer p01aaEb, isave, isizea, isizeb, istepa, istepb,
     *     ifail, is, it, n, i
      real  sum    !cpis double precision -> real
      character *8 srname
      real a(isizea), b(isizeb), c1, c2, d1, d2, x
      logical sw
      data srname / ' x03aaEb'/
      isave = ifail
      ifail = 0
      if (istepa.gt.0 .and. istepb.gt.0) go to 10
      ifail = p01aaEb(isave,1,srname)
      return
   10 is = 1 - istepa
      it = 1 - istepb
      if (sw) go to 40
      x = 0.0
      if (n.lt.1) go to 30
      do 20 i=1,n
           is = is + istepa
           it = it + istepb
           x = x + a(is)*b(it)
   20 continue
   30 d1 = x + (c1+c2)
      d2 = 0.0
      return
   40 sum = 0.0
      if (n.lt.1) go to 60
      do 50 i=1,n
           is = is + istepa
           it = it + istepb
           sum = sum + real(a(is))*b(it) !cpis dble -> real
   50 continue
   60 sum = sum + (real(c1)+c2)  !cpis dble -> real
      call aazx03(sum, d1, d2)
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@aazx03  11040/bald89/weqvmom DEQVMOM
c
      subroutine aazx03(dp, d1, d2)
      real dp     !cpis double precision -> real
      real d1, d2
      d1 = dp
      d2 = dp - d1
      return
c     return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@sdot  11040/bald89/weqvmom DEQVMOM
c
c      function sdot(n, sx, incx, sy, incy)
c      integer   n, incx, incy, ix, iy, i
c      real   sx(1), sy(1), sum
c      sum = 0.0
c      ix = 1
c      iy = 1
c      do 10 i=1,n
c           sum = sum + sx(ix)*sy(iy)
c           ix = ix + incx
c           iy = iy + incy
c   10 continue
c      sdot = sum
c      return
c      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@i1machb  11040/bald89/weqvmom DEQVMOM
c
      integer function i1machb(i)
c
c  i/o unit numbers.
c
c    i1machb( 1) = the standard input unit.
c
c    i1machb( 2) = the standard output unit.
c
c    i1machb( 3) = the standard punch unit.
c
c    i1machb( 4) = the standard error message unit.
c
c  words.
c
c    i1machb( 5) = the number of bits per integer storage unit.
c
c    i1machb( 6) = the number of characters per integer storage unit.
c
c  integers.
c
c    assume integers are represented in the s-digit, base-a form
c
c               sign ( x(s-1)*a**(s-1) + ... + x(1)*a + x(0) )
c
c               where 0 .le. x(i) .lt. a for i=0,...,s-1.
c
c    i1machb( 7) = a, the base.
c
c    i1machb( 8) = s, the number of base-a digits.
c
c    i1machb( 9) = a**s - 1, the largest magnitude.
c
c  floating-point numbers.
c
c    assume floating-point numbers are represented in the t-digit,
c    base-b form
c
c               sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
c
c               where 0 .le. x(i) .lt. b for i=1,...,t,
c               0 .lt. x(1), and emin .le. e .le. emax.
c
c    i1machb(10) = b, the base.
c
c  single-precision
c
c    i1machb(11) = t, the number of base-b digits.
c
c    i1machb(12) = emin, the smallest exponent e.
c
c    i1machb(13) = emax, the largest exponent e.
c
c  double-precision
c
c    i1machb(14) = t, the number of base-b digits.
c
c    i1machb(15) = emin, the smallest exponent e.
c
c    i1machb(16) = emax, the largest exponent e.
c
c  to alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.  also, the values of
c  i1machb(1) - i1machb(4) should be checked for consistency
c  with the local operating system.
c
      integer imach(16), output
c
c      equivalence (imach(4),output)  !use executable stmt---cray trouble
c
c     machine constants for the vax 11/780
c
cv      data imach(1) /    5 /
cv      data imach(2) /    6 /
cv      data imach(3) /    5 /
cv      data imach(4) /    6 /
cv      data imach(5) /   32 /
cv      data imach(6) /    4 /
cv      data imach(7) /    2 /
cv      data imach(8) /   31 /
cv      data imach(9) /2147483647 /
cv      data imach(10)/    2 /
cv      data imach(11)/   24 /
cv      data imach(12)/ -127 /
cv      data imach(13)/  127 /
cv      data imach(14)/   56 /
cv      data imach(15)/ -127 /
cv      data imach(16)/  127 /
c
c     machine constants for the burroughs 1700 system.
c
c     data imach( 1) /    7 /
c     data imach( 2) /    2 /
c     data imach( 3) /    2 /
c     data imach( 4) /    2 /
c     data imach( 5) /   36 /
c     data imach( 6) /    4 /
c     data imach( 7) /    2 /
c     data imach( 8) /   33 /
c     data imach( 9) / z1ffffffff /
c     data imach(10) /    2 /
c     data imach(11) /   24 /
c     data imach(12) / -256 /
c     data imach(13) /  255 /
c     data imach(14) /   60 /
c     data imach(15) / -256 /
c     data imach(16) /  255 /
c
c     machine constants for the burroughs 5700 system.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  48 /
c     data imach( 6) /   6 /
c     data imach( 7) /   2 /
c     data imach( 8) /  39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /   8 /
c     data imach(11) /  13 /
c     data imach(12) / -50 /
c     data imach(13) /  76 /
c     data imach(14) /  26 /
c     data imach(15) / -50 /
c     data imach(16) /  76 /
c
c     machine constants for the burroughs 6700/7700 systems.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  48 /
c     data imach( 6) /   6 /
c     data imach( 7) /   2 /
c     data imach( 8) /  39 /
c     data imach( 9) / o0007777777777777 /
c     data imach(10) /   8 /
c     data imach(11) /  13 /
c     data imach(12) / -50 /
c     data imach(13) /  76 /
c     data imach(14) /  26 /
c     data imach(15) / -32754 /
c     data imach(16) /  32780 /
c
c     machine constants for the cdc 6000/7000 series.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    7 /
c     data imach( 4) /    6 /
c     data imach( 5) /   60 /
c     data imach( 6) /   10 /
c     data imach( 7) /    2 /
c     data imach( 8) /   48 /
c     data imach( 9) / 00007777777777777777b /
c     data imach(10) /    2 /
c     data imach(11) /   48 /
c     data imach(12) / -974 /
c     data imach(13) / 1070 /
c     data imach(14) /   96 /
c     data imach(15) / -927 /
c     data imach(16) / 1070 /
c
c     machine constants for the cray 1
c
c      data imach( 1) /   100 /
c      data imach( 2) /   5 /
c      data imach( 3) /   102 /
c      data imach( 4) /   101 /
c      data imach( 5) /    64 /
c      data imach( 6) /     8 /
c      data imach( 7) /     2 /
c      data imach( 8) /    63 /
c      data imach( 9) /  777777777777777777777b /
c      data imach(10) /     2 /
c      data imach(11) /    48 /
c      data imach(12) / -8192 /
c      data imach(13) /  8191 /
c      data imach(14) /    96 /
c      data imach(15) / -8192 /
c      data imach(16) /  8191 /
c
c     machine constants for the data general eclipse s/200
c
c     data imach( 1) /   11 /
c     data imach( 2) /   12 /
c     data imach( 3) /    8 /
c     data imach( 4) /   10 /
c     data imach( 5) /   16 /
c     data imach( 6) /    2 /
c     data imach( 7) /    2 /
c     data imach( 8) /   15 /
c     data imach( 9) /32767 /
c     data imach(10) /   16 /
c     data imach(11) /    6 /
c     data imach(12) /  -64 /
c     data imach(13) /   63 /
c     data imach(14) /   14 /
c     data imach(15) /  -64 /
c     data imach(16) /   63 /
c
c     machine constants for the harris 220
c
c     data imach( 1) /       5 /
c     data imach( 2) /       6 /
c     data imach( 3) /       0 /
c     data imach( 4) /       6 /
c     data imach( 5) /      24 /
c     data imach( 6) /       3 /
c     data imach( 7) /       2 /
c     data imach( 8) /      23 /
c     data imach( 9) / 8388607 /
c     data imach(10) /       2 /
c     data imach(11) /      23 /
c     data imach(12) /    -127 /
c     data imach(13) /     127 /
c     data imach(14) /      38 /
c     data imach(15) /    -127 /
c     data imach(16) /     127 /
c
c     machine constants for the honeywell 600/6000 series.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /   43 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    6 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -127 /
c     data imach(13) /  127 /
c     data imach(14) /   63 /
c     data imach(15) / -127 /
c     data imach(16) /  127 /
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c     data imach( 1) /   5 /
c     data imach( 2) /   6 /
c     data imach( 3) /   7 /
c     data imach( 4) /   6 /
c     data imach( 5) /  32 /
c     data imach( 6) /   4 /
c     data imach( 7) /   2 /
c     data imach( 8) /  31 /
c     data imach( 9) / z7fffffff /
c     data imach(10) /  16 /
c     data imach(11) /   6 /
c     data imach(12) / -64 /
c     data imach(13) /  63 /
c     data imach(14) /  14 /
c     data imach(15) / -64 /
c     data imach(16) /  63 /
c
c     machine constants for the pdp-10 (ka processor).
c
c     data imach(1) /5/
c     data imach(2) /6/
c     data imach(3) /5/
c     data imach(4) /6/
c     data imach(5) /36/
c     data imach(6) /5/
c     data imach(7) /2/
c     data imach(8) /35/
c     data imach( 9) / "377777777777 /
c     data imach(10) /2/
c     data imach(11) /27/
c     data imach(12) /-128/
c     data imach(13) /127/
c     data imach(14) /54/
c     data imach(15) /-101/
c     data imach(16) /127/
c
c     machine constants for the pdp-10 (ki processor).
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    5 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / "377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   62 /
c     data imach(15) / -128 /
c     data imach(16) /  127 /
c
c     machine constants for pdp-11 fortran's supporting
c     32-bit integer arithmetic.
c
      data imach( 1) /    5 /
      data imach( 2) /    6 /
      data imach( 3) /    5 /
      data imach( 4) /    6 /
      data imach( 5) /   32 /
      data imach( 6) /    4 /
      data imach( 7) /    2 /
      data imach( 8) /   31 /
      data imach( 9) / 2147483647 /
      data imach(10) /    2 /
      data imach(11) /   24 /
      data imach(12) / -127 /
      data imach(13) /  127 /
      data imach(14) /   56 /
      data imach(15) / -127 /
      data imach(16) /  127 /
c
c     machine constants for pdp-11 fortran's supporting
c     16-bit integer arithmetic.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    5 /
c     data imach( 4) /    6 /
c     data imach( 5) /   16 /
c     data imach( 6) /    2 /
c     data imach( 7) /    2 /
c     data imach( 8) /   15 /
c     data imach( 9) / 32767 /
c     data imach(10) /    2 /
c     data imach(11) /   24 /
c     data imach(12) / -127 /
c     data imach(13) /  127 /
c     data imach(14) /   56 /
c     data imach(15) / -127 /
c     data imach(16) /  127 /
c
c     machine constants for the univac 1100 series.
c
c     note that the punch unit, i1machb(3), has been set to 7
c     which is appropriate for the univac-for system.
c     if you have the univac-ftn system, set it to 1.
c
c     data imach( 1) /    5 /
c     data imach( 2) /    6 /
c     data imach( 3) /    7 /
c     data imach( 4) /    6 /
c     data imach( 5) /   36 /
c     data imach( 6) /    6 /
c     data imach( 7) /    2 /
c     data imach( 8) /   35 /
c     data imach( 9) / o377777777777 /
c     data imach(10) /    2 /
c     data imach(11) /   27 /
c     data imach(12) / -128 /
c     data imach(13) /  127 /
c     data imach(14) /   60 /
c     data imach(15) /-1024 /
c     data imach(16) / 1023 /
c
      output=imach(4)
      if (i.lt.1 .or. i.gt.16) go to 10
c
      i1machb = imach(i)
      return
c
   10 write (output,99999)
c
c
      stop
c
99999 format ('1error    1 in i1machb - i out of bounds')
      end

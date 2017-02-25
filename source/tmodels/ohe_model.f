!| %
!| %  This is a LaTeX ASCII file.  To typeset document type: latex theory
!| %  To extract the fortran source code, obtain the utility "xtverb" from
!|  0lenn Bateman (bateman@plasma.physics.lehigh.edu) and type:
!| 0tverb < theory.tex > theory.f
!| %
!| % The following lines control how the LaTeX document is typeset
!| 
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| \begin{document}            0.000000E+00nd of preamble and beginning of text.
!| 
!| \begin{center}
!| \Large {\bf OHE Transport Model} \\
!| \vspace{1pc} \normalsize
!| Matteo Erba, Glenn Bateman, Arnold Kritz \\
!|  Lehigh University Physics Department \\
!| 16 Memorial Drive East, Bethlehem PA 18015 \\
!| erba@plasma.physics.lehigh.edu \\
!| bateman@plasma.physics.lehigh.edu \\
!| kritz@plasma.physics.lehigh.edu
!| \end{center}
!| 
!| This file documents a subroutine which computes
!| plasma transport coefficients using the OHE transport model\cite{ott97}.
!| A complete derivation of the mixed model is given in reference.
!| 
c@ohe_model.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine ohe_model (
     &   rminor,  rmajor,   elong
     & , tekev,    tikev,    q
     & , btor, aimass
     & , grdnz, grdte, grdq, grdnh
     & , wexbs
     & , grdne,   grdni
     & , dense,   densh,    densimp
     & , grdti
     & , thiig,   thdig,    theig,    thzig
     & , thirb,   thdrb,    therb,    thzrb
     & , thikb,   thdkb,    thekb,    thzkb
     & , difthi,   velthi,  vflux
     & , matdim,  npoints,  nprout,   lprint,  nerr
     & , lsuper,  lreset,   lswitch, cswitch, fig,    frb,     fkb)
c
c
c    All the following 1-D arrays are assumed to be defined on flux
c    surfaces called zone boundaries where the transport fluxes are
c    to be computed.  The number of flux surfaces is given by npoints
c    (see below).  For example, if you want to compute the transport
c    on only one flux surface, set npoints = 1.
c
c  Input arrays:
c  -------------
c
c  rminor(jz)   = minor radius (half-width) of zone boundary [m]
c  rmajor(jz)   = major radius to geometric center of zone bndry [m]
c  elong(jz)    = local elongation of zone boundary
c
c  tekev(jz)    = T_e (electron temperature) [keV]
c  tikev(jz)    = T_i (temperature of thermal ions) [keV]
c  q(jz)        = magnetic q-value
c  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
c
c  aimass(jz)   = mean atomic mass of thermal ions [AMU]
c               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
c                 sum_i = sum over all ions, each with mass M_i
c
c  wexbs(jz)    = ExB shearing rate in [rad/s]
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c
c  grdnh(jz) = -R ( d n_h / d r ) / n_h
c  grdnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
c  grdte(jz) = -R ( d T_e / d r ) / T_e
c  grdq (jz) =  R ( d q   / d r ) / q    related to magnetic shear
c
c  where:
c    n_i     = thermal ion density (sum over hydrogenic and impurity)
c    n_h     = thermal hydrogenic density (sum over hydrogenic species)
c    n_Z     = thermal impurity density,  Z = average impurity charge
c                      sumed over all impurities
c
c  Output:
c  -------
c
c    The following effective diffusivities represent contributions
c    to the total diffusivity matrix (difthi and velthi given below)
c    from each of the models that contribute to the Multi-Mode model.
c    Generally, these arrays are used for diagnostic output only.
c
c  thiig(jz) = ion thermal diffusivity from the OHE model
c  thdig(jz) = hydrogenic ion diffusivity from the OHE model
c  theig(jz) = elelctron thermal diffusivity from the OHE model
c  thzig(jz) = impurity ion diffusivity from the OHE model
c
c
c    All of the transport coefficients are given in the following two
c    matricies for diffusion difthi and convection velthi in MKS units.
c    See the LaTeX documentation for difthi and velthi just below.
c
c    NOTE:  difthi and velthi include all of the anomalous transport.
c    There are no additional contributions to the heat fluxes from
c    charged particle convection.
c
c  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c  velthi(j1,jz)    = convective velocities
c  vflux(j1,jz)     = flux matrix
!| 
!| The full matrix form of anomalous transport has the form
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots
!|     \end{array} \right)
!|  = - \nabla \cdot
!| \left( \begin{array}{l} {\rm vFlux}_1 \; n_H T_H \\
!|  {\rm vFlux}_2 \; n_H \\
!|  {\rm vFlux}_3 \; n_e T_e \\
!|  {\rm vFlux}_4 \; n_Z \\
!|  {\rm vFlux}_5 \; n_Z T_Z \\
!|  \vdots \end{array} \right)
!|  + \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  = \nabla \cdot
!| \left( \begin{array}{llll}
!| D_{1,1} n_H & D_{1,2} T_H & D_{1,3} n_H T_H / T_e \\
!| D_{2,1} n_H / T_H & D_{2,2} & D_{2,3} n_H / T_e \\
!| D_{3,1} n_e T_e / T_H & D_{3,2} n_e T_e / n_H & D_{3,3} n_e & \vdots \\
!| D_{4,1} n_Z / T_H & D_{4,2} n_Z / n_H & D_{4,3} n_Z / T_e \\
!| D_{5,1} n_Z T_Z / T_H & D_{5,2} n_Z T_Z / n_H &
!|         D_{5,3} n_Z T_Z / T_e \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!|  \nabla
!|  \left( \begin{array}{c}  T_H \\ n_H \\  T_e \\
!|    n_Z \\  T_Z \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  + \nabla \cdot
!| \left( \begin{array}{l} {\bf v}_1 \; n_H T_H \\ {\bf v}_2 \; n_H \\
!|    {\bf v}_3 \; n_e T_e \\
!|    {\bf v}_4 \; n_Z \\ {\bf v}_5 \; n_Z T_Z \\
!|     \vdots \end{array} \right) +
!|  \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right) $$
!| Note that all the diffusivities are in units of m$^2$/sec while the
!| convective velocities and vfluxes are in units of m/sec.
!| 
!| WARNING:  Do not add separate convective transport terms to this
!| anomalous transport model.  All the anomalous transport
!| predicted by this Multi-Mode model is contained
!| in the diffusion coefficients {\tt difthi} and {\tt velthi} given
!| above.
!| 
c
c  Input integers:
c  ---------------
c
c  matdim  = first and second dimension of transport matricies
c            difthi(j1,j2,jz) and velthi(j1,jz) and the first
c            dimension of gamma and omega.  matdim must be at least 5
c
c  npoints = number of values of jz in all of the above arrays
c
c  nprout  = output unit number for long printout
c
c  nerr    = 0 on input; returning with value .ne. 0 indicates error
c
c  Input switches
c  --------------
c
c  lprint      controls the amount of printout (0 => no printout)
c              higher values yield more diagnostic output
c
c  lsuper  > 0 for supershot simulations
c          = 0 for simulations of all other discharges
c
c  lreset  = 0 to use only internal settings for lswitch, cswitch
c              and for the coefficients fig, frb, and fkb
c
c    Note that when lreset = 0, the values of the switches and
c    coefficients in the argument list are ignored and all the
c    switches and coefficients are set internally.
c
c    WARNING:  use lreset > 0 only if you want to pass all the switches
c    lswitch, cswitch, fig, frb, and fkb through the argument list.
c
c  Internal control variables:
c  ---------------------------
c
c  lswitch(j), j=1,8   integer control variables:
c
c  cswitch(j), j=1,25   general control variables:
c
c  lswitch(5) = 1 to limit magnitude of all normalized gradients
c                    to ( major radius ) / ( ion Larmor radius )
c
c  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
c
c
c***********************************************************************
c
c-----------------------------------------------------------------------
c
c  Compile this routine and routines that it calls with a compiler
c  option, such as -r8, to convert real to double precision when used on
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: ohe_model calls the following routines
c
c  ohe             - Computes diffusivity from the OHE model
c
c-----------------------------------------------------------------------

      implicit none
c
      integer klswitch, kcswitch
c
      parameter ( klswitch = 8, kcswitch = 25 )
c
      integer  matdim,  npoints, nprout,  lprint,   nerr
c
      real
     &   rminor(*),  rmajor(*),   elong(*)
     & , tekev(*),    tikev(*),    q(*),       btor(*)
     & , wexbs(*), aimass(*)
     & , grdti(*)
c
      real
     &   thiig(*),   thdig(*),    theig(*),    thzig(*)
     & , thirb(*),   thdrb(*),    therb(*),    thzrb(*)
     & , thikb(*),   thdkb(*),    thekb(*),    thzkb(*)
     & , densimp(*), grdni(*), grdne(*), densh(*), dense(*)
     & , difthi(matdim,matdim,*), velthi(matdim,*)
     & , vflux(matdim,*)
     & , grdnz(*), grdte(*), grdq(*), grdnh(*)
c
      real     cswitch(*)
c
      integer  lsuper,  lreset,  lswitch(*)
c
      real     fig(*),  fkb(*),  frb(*)
c
c
c-----------------------------------------------------------------------
c
c
c..physical constants
c
      real zckb,  zcme,  zcmp,  zce
c
c  zckb    = energy conversion factor        [Joule/keV]
c  zcme    = electron mass                   [kg]
c  zcmp    = proton mass                     [kg]
c  zce     = electron charge                 [Coulomb]
c
c
c..local variables
c
      integer  jz, j1, j2, jm
c
c..local variables connected to the ohe module
c
      integer  ierr, switch

      real  zq, chibohm,  zeps,    zsound,     zgyrfi,    zrhos,    zrlt
     & , zrltcrit,     zra,    zkappa,   zkapexp,    zbound,   zcoef1
     & , zcoef2, zcoef3, zcoef4, chii, chie, dh, dimp
     & , zrmaj, zshear, zai, zbtor, zrmin
     & , zgth, zgmax, zgte, zgnh, zgnz, zlarpo
     & , zvthi, zlari, zep, zti, zomegde, zkyrho, zelfkb
     & , zgpr, zne, zte, zgne, zni, zgni
     & , zcmu0, zbeta, zepslon, zlgeps, zfbthn, zbprim, zbc1, zdk
     & , zshat, zscyl, zpi, zelonf, zelong, zsmin
c
c.. variables for exb model
c
      real  zwexb
c
c  zwexb    = local copy of ExB shearing rate
c
c
c..initialize arrays
c
      do jz = 1, npoints
        thiig(jz)  = 0.
        thdig(jz)  = 0.
        theig(jz)  = 0.
        thzig(jz)  = 0.
        therb(jz)  = 0.
        thirb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
        thikb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
      enddo
c
      do jz = 1, npoints
        do j1 = 1, matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          do j2 = 1, matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
      enddo
c
c..if lreset < 1, use internal settings for switches and coefficients
c  otherwise, use values passed through the argument list above
c
      if ( lreset .lt. 1 ) then
c
c..initialize switches
c
        do j1=1,kcswitch
          cswitch(j1) = 0.0
        enddo
c
        do j1=1,klswitch
          cswitch(j1) = 0.0
        enddo
c
        cswitch(1)  =  0.5  ! minimum value of shear
        cswitch(2)  =  3.5  ! for fbeta-th in kinetic ballooning
        cswitch(3)  = -4.0  ! elongation scaling for OHE model
        cswitch(6)  =  0.0  ! k_y*rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
        lswitch(5) = 1  ! limit gradients by major radius / ion Larmor radius
c
c  contributions to hydrogenic particle, elec-energy, ion-energy,
c    and impurity ion fluxes

        fig(1) = 1.00
        fig(2) = 1.00
        fig(3) = 1.00
        fig(4) = 1.00
        fkb(1) = 1.00
       	fkb(2) = 0.65
        fkb(3) = 0.65
        fkb(4) = 1.00
c
c
c
      endif

!| 
!| We then enter a loop over the spatial zones, and set the following
!| BALDUR {\tt common} variables to variables local to subroutine {\tt
!| theory}: the mean atomic mass number of the thermal ions, $A_{i}$
!| ({\tt zai}), the electron and ion density and temperature, $n_{e}$
!| ({\tt zne}), $n_{i}=\sum_{a}n_{a}$ ({\tt zni}), $T_{e}$ ({\tt zte}),
!| and $T_{i}$ ({\tt zti}), the safety factor, $q$ ({\tt zq}), the
!| effective charge, $Z_{eff}$ ({\tt zeff}),the midplane halfwidth of a
!| flux surface, $r$ ({\tt zrmin}), the major radius, $R_{o}$ ({\tt
!| zrmaj}), and the toroidal field at $R_0$ major radius, $B_{0}$
!| ({\tt zbtor}).  Therefore, the only variables not defined in
!| subroutine {\tt mmm95} that are needed to complete the rest of the
!| calculation are $\pi$ ({\tt zpi}), the small overflow protection
!| variable {\tt zepslon}.  The points of defining so many local
!| variables are to compact the notation.  The relevant coding for the
!| calculations just described is:
!| 
c-----------------------------------------------------------------------
c
c..physical constants
c
        zpi     = atan2 ( 0.0, -1.0 )
        zcmu0   = 4.0e-7 * zpi
        zckb    = 1.60210e-16
        zcme    = 9.1091e-31
        zcmp    = 1.67252e-27
        zce     = 1.60210e-19
        zkyrho = 0.316
c
c
c.. start the main do-loop over the radial index "jz"..........
c
c
      do 300 jz = 1, npoints
c
c  compute quantities necessary for mixed model
c
      zgte   = grdte(jz)
      zgne   = grdne(jz)
      zgth   = grdti(jz)
      zgni   = grdni(jz)
      zne    = dense(jz)
      zni    = densh(jz) + densimp(jz)
      zte    = tekev(jz)
      zti    = tikev(jz)
      zelong = elong(jz)
      zai    = aimass(jz)
      zrmin  = rminor(jz)
      zrmaj  = rmajor(jz)
      zshear = grdq(jz) * zrmin / zrmaj
      zbtor  = btor(jz)
      zq     = q(jz)
      chibohm = 16.0*62.5*tekev(jz)/btor(jz)
      zeps = rminor(jz) / rmajor(jz)
      zsound = sqrt(zckb * tekev(jz) / (zcmp * aimass(jz)))
      zgyrfi = zce * zbtor / (zcmp * zai)
      zrhos  = zsound / zgyrfi
      zrlt= abs(grdti(jz))
      zrltcrit = 0.0
      zra = rminor(jz) / rminor(npoints)
      zkappa = elong(jz)
      zkapexp = cswitch(3)
      zbound = 0.0
      zcoef1 = 1.00
      zcoef2 = 0.25
      zcoef3 = 1.25
      zcoef4 = 1.25
      switch = 0
c
c
c...Define a local copy of normalized ExB shearing rate : pis
c
        zomegde = 2.0 * zkyrho * zsound / zrmaj
c
        zwexb = cswitch(21) * wexbs(jz) / zomegde
c
c
!| 
!| 
!| In the order in which they are computed, algebraic notation for
!| these quantities is:
!| $$ p=n_e T_e + n_i T_i \ --- \ {\rm (thermal)}\eqno{\tt zprth} $$
!| $$ \omega_{ci}=eB_{o}/(m_{p}A_{i}) \eqno{\tt zgyrfi} $$
!| $$ \beta=(2\mu_{o}k_{b}/B_{o}^{2})(n_{e}T_{e}+n_{i}T_{i})
!|  \eqno{\tt zbeta} $$
!| $$ v_{e}=(2k_{b}T_{e}/m_{e})^{1/2} \eqno{\tt zvthe} $$
!| $$ v_{i}=(2k_{b}T_{i}/m_{p}A_{i})^{1/2} \eqno{\tt zvthi} $$
!| $$ c_{s}=[k_{b}T_{e}/(m_{p}A_{i})]^{1/2} \eqno{\tt zsound} $$
!| $$ \ln (\lambda)=37.8 - \ln (n_{e}^{1/2}T_{e}^{-1}) \eqno{\tt zlog} $$
!| $$ \nu_{ei}=4(2\pi)^{1/2}n_{e}(\ln \lambda)e^{4}Z_{eff}
!|                /[3(4\pi \epsilon_{o})^{2}m_{e}^{1/2}(k_{b}T_{e})^{3/2}]
!|  \eqno{\tt znuei} $$
!| $$ \eta=\nu_{ei}/(2\epsilon_{o}\omega_{pe}^{2}) \eqno{\tt zresis} $$
!| $$ \nu_{eff}=\nu_{ei}/\epsilon \eqno{\tt znueff} $$
!| $$ \nu_{e}^{*}=\nu_{ei}qR_{o}/(\epsilon^{3/2}v_{e}) \eqno{\tt thnust} $$
!| $$ \hat{\nu}=\nu_{eff}/\omega_{De} \eqno{\tt znuhat} $$
!| $$ \rho_{\theta i}=\rho_{i}q/\epsilon \eqno{\tt zlari} $$
!| $$ \rho_{i}=v_{i}/\omega_{ci} \eqno{\tt zlarpo} $$
!| $$ \rho_{s}=c_{s}/\omega_{ci} \eqno{\tt zrhos} $$
!| $$ k_{\perp}=0.3/\rho_{s} \eqno{\tt zwn} $$
!| 
!| 
!| Our formulas for the shear begin with
!| $$ {\hat s}_{cyl}=|(r/q)(\partial q/\partial r)| \eqno{\tt zscyl} $$
!| computed earlier in this subroutine.
!| The minimum prescribed shear is
!| $$ {\hat s}_{min}=max(c_{1},0) \eqno{\tt zsmin} $$
!| where $c_1=$ = {\tt cswitch(1)} so that shear is then given by
!| $$ {\hat s}=max({\hat s}_{min},{\hat s}_{cyl}) \eqno{\tt zshat} $$
c
        thiig(jz) = 0.0
        theig(jz) = 0.0
        thdig(jz) = 0.0
        thzig(jz) = 0.0
c
c
c
c
c
      call ohe(
     &   chibohm,     zq,         zeps,       zrhos,   zrlt
     & , zrltcrit,    zra,        chii,       chie,       dh
     & , dimp,        zkappa,     zkapexp,    zbound,     ierr
     & , zcoef1,      zcoef2,     zcoef3,     zcoef4,     switch)
c
c  If ierr not equal to 0 an error has occured
c
	if (ierr .ne. 0) return
c
c  compute effective diffusivites for diagnostic purposes only
c
c        thdig(jz) = fig(1) * dh
c        theig(jz) = fig(2) * chie
c        thiig(jz) = fig(3) * chii
c        thzig(jz) = fig(4) * dimp
         thdig(jz) = fig(1) * dh
         theig(jz) = fig(2) * chie
         thiig(jz) = fig(3) * chii
         thzig(jz) = fig(4) * dimp
c
c   Define kinetic ballooning model
c
      zgpr = ( zne * zte * ( zgne + zgte )
     &         + zni * zti * ( zgni + zgth ) )
     &         / ( zne * zte + zni * zti )
      zbeta  = (2. * zcmu0 * zckb / zbtor**2) * (zne * zte + zni * zti)
      zbprim = abs( zbeta * zgpr / zrmaj )
      zepslon = 1.0e-34
      zelonf = ( 1. + zelong**2 ) / 2.
      zlgeps  = log ( zepslon )
      zscyl=max(abs(zshear),zepslon)
      zsmin=max(cswitch(1),zepslon)
      zshat=max(zsmin,zscyl)
      zbc1   = abs(zshat)/(1.7*zq**2*zrmaj)
      zfbthn = exp( min(abs(zlgeps),
     &     max(-abs(zlgeps),cswitch(2)*(zbprim/zbc1 - 1.))) )
      zdk = abs( zsound * zrhos**2 * zfbthn * zgpr / zrmaj )
      zelfkb = zelonf**cswitch(3)
      thdkb(jz) = zdk*fkb(1)*zelfkb
      thekb(jz) = zdk*fkb(2)*zelfkb
      thikb(jz) = zdk*fkb(3)*zelfkb
      thzkb(jz) = zdk*fkb(4)*zelfkb
c
c  start computing the fluxes
c
        zep    = zrmin/zrmaj
c        zep    = max( zrmin/zrmaj, zepslon )
        zti    = tikev(jz)
        zvthi  = sqrt(2. * zckb * zti / (zcmp * zai))
        zlari  = zvthi / zgyrfi
c        zlarpo = max(zlari * zq / zep, zepslon)
        zlarpo = zlari * zq / zep
        zgmax = zrmaj / zlarpo
      	zgnh   = grdnh(jz)
      	zgnz   = grdnz(jz)

      	if ( lswitch(5) .eq. 1 ) then
c
        zgnh = sign ( min ( abs ( zgnh ), zgmax ), zgnh )
        zgnz = sign ( min ( abs ( zgnz ), zgmax ), zgnz )
        zgte = sign ( min ( abs ( zgte ), zgmax ), zgte )
        zgth = sign ( min ( abs ( zgth ), zgmax ), zgth )
c
        endif

        vflux(1,jz) = vflux(1,jz) + thiig(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdig(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + theig(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzig(jz) * zgnz / zrmaj
c
        vflux(1,jz) = vflux(1,jz) + thikb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdkb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + thekb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzkb(jz) * zgnz / zrmaj
c
c  compute diffusivity matrix
c
        do j1=1,matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          do j2=1,matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
c
c..set diffthi and velthi
c
c
c  diagonal elements of matrix = effective diffusivities
c
          difthi(1,1,jz) = difthi(1,1,jz) + thiig(jz) + thikb(jz)
          difthi(2,2,jz) = difthi(2,2,jz) + thdig(jz) + thdkb(jz)
          difthi(3,3,jz) = difthi(3,3,jz) + theig(jz) + thekb(jz)
          difthi(4,4,jz) = difthi(4,4,jz) + thzig(jz) + thzkb(jz)
c
c..end of ohe model
c
c
 300  continue
c
c
c   end of the main do-loop over the radial index, "jz"----------
c
      return
      end
!| 
!| %**********************************************************************c
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{ott97}
!| M. Ottaviani, W. Horton, and M. Erba,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 39} (1996) 1461.
!| 
!| \end{thebibliography}
!| 
!| %**********************************************************************c
!| \end{document}              0.000000E+00nd of document.
!| 
!| 
!| 
!| 
!| 
!| 
!| 
!| 
!| 

!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| 
!| \documentstyle[12pt]{article}
!| \headheight 0pt \headsep 0pt
!| \topmargin 0pt  \textheight 9.0 in
!| \oddsidemargin 0pt \textwidth 6.5 in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center}
!| \large {\bf Drift Alfv\'{e}n Transport Model} \\
!| \normalsize  {\tt sda01dif.tex} \\
!| \vspace{1pc}
!| Bruce Scott \\
!| Max Planck Institut f\"{u}r Plasmaphysik \\
!| Euratom Association \\
!| D-85748 Garching, Germany \\
!| \vspace{1pc}
!| Glenn Bateman \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| \today
!| \end{center}
!| 
!| This subroutine computes the transport diffusion matrix driven by
!| drift Alfv\'{e}n turbulence usg a model based on
!| simulations by Bruce Scott.
!| 
c@sda04dif
c rgb 10-oct-98 revised with new model
c rgb 26-feb-98 first draft
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine sda04dif (lswitch, cswitch, nswitch, ndim
     & , lprint, nout
     & , gne, gnh, gnz, gte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , fluxth, fluxnh, fluxte, fluxnz, fluxtz
     & , diffth, diffnh, diffte, diffnz, difftz
     & , condif, convel
     & , nerr )
c
c    input:
c
c  lswitch(j), j=1,nswitch  integer switches
c  cswitch(j), j=1,nswitch  real-valued switches
c  nswitch  = number of switches
c  ndim     = first dimension of condif array (must be 5 or more)
c  lprint    controls amount of printout
c  nout      output unit
c
c  gne      = - R ( d n_e / d r ) / n_e  electron density gradient
c  gnh      = - R ( d n_H / d r ) / n_H  hydrogen density gradient
c  gnz      = - R ( d n_Z / d r ) / n_Z  impurity density gradient
c  gte      = - R ( d T_e / d r ) / T_e  electron temperature gradient
c  gth      = - R ( d T_H / d r ) / T_H  hydrogen temperature gradient
c  gtz      = - R ( d T_Z / d r ) / T_Z  impurity temperature gradient
c  th7te    = T_H / T_e  ratio of hydrogen to electron temperature
c  tz7te    = T_Z / T_e  ratio of impurity to electron temperature
c  fracnz   = n_Z / n_e  ratio of impurity to electron density
c  chargenz = average charge of impurity ions
c  amassh   = average mass of hydrogen ions (in AMU)
c  amassz   = average mass of impurity (in AMU)
c  fracns   = n_s / n_e  ratio of fast ion to electron density
c  chargens = average charge of fast ions
c  betae    = n_e T_e / ( B^2 / 2 mu_o )
c  vef      = collisionality
c  q        = magnetic q value
c  shear    = r ( d q / d r ) / q
c  kappa    = elongation
c
c  zdgth, zdgnz, zdgte, zdgnz, zdgtz
c    are the finite differences in gth, gnh, gte, gnz, and gtz
c
c    output:
c
c  fluxth   = radial velocity of hydrogenic heat / ( n_e T_e c_s )
c  fluxnh   = radial velocity of hydrogenic ions / ( n_e c_s )
c  fluxte   = radial velocity of electron heat / ( n_e T_e c_s )
c  fluxnz   = radial velocity of impurity ions / ( n_e c_s )
c  fluxtz   = radial velocity of impurity heat / ( n_e T_e c_s )
c
c  diffth   = diffusivity of hydrogenic heat / ( n_e T_e c_s gth )
c  diffnh   = diffusivity of hydrogenic ions / ( n_e c_s gnh )
c  diffte   = diffusivity of electron heat / ( n_e T_e c_s gte )
c  diffnz   = diffusivity of impurity ions / ( n_e c_s gnz )
c  difftz   = diffusivity of impurity heat / ( n_e T_e c_s gtz )
c
c  condif   = conductive diffusion matrix with dimension (5,5)
c  convel   = convective velocity vector with dimension (5)
c
c  nerr     = output error condition
c
      implicit none
c
      integer lswitch(*), nswitch, ndim, lprint, nout, nerr
c
      real cswitch(*)
     & , gne, gnh, gnz
     & , gte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , fluxth, fluxnh, fluxte, fluxnz, fluxtz
     & , diffth, diffnh, diffte, diffnz, difftz
c
      real condif(ndim,*), convel(*)
c
      real zgneut, zfracnh
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz
     & , zdgth, zdgnh, zdgte, zdgnz, zdgtz
     & , zfluxth, zfluxnh, zfluxte, zfluxnz, zfluxtz
c
c  zgneut = gne - gnh * n_H/n_e - gnz * Z n_Z / n_e
c  zfracnh = n_H / n_e
c
      real zero, zone, zhalf, zepslon
c
c  zone  = 1.0
c  zhalf = 1.0 / 2
c
      zero  = 0
      zone  = 1
      zhalf = zone / 2
      zepslon = 1.e-10
c
      zfracnh = zone - chargenz * fracnz - chargens * fracns
c
      if ( zfracnh .lt. zero ) then
        zfracnh = zero
        nerr = 2
      endif
c
      zgneut = gne - gnh * zfracnh - gnz * chargenz * fracnz
c
      call sda04flx (lswitch, cswitch, lprint, nout
     & , gne, gnh, gnz, gte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , fluxth, fluxnh, fluxte, fluxnz, fluxtz
     & , nerr )
c
c..compute total effective diffusivities
c
      zgth = gth
      if ( abs(zgth) .lt. zepslon ) zgth = zepslon
      diffth = fluxth / zgth
c
      zgnh = gnh
      if ( abs(zgnh) .lt. zepslon ) zgnh = zepslon
      diffnh = fluxnh / zgnh
c
      zgte = gte
      if ( abs(zgte) .lt. zepslon ) zgte = zepslon
      diffte = fluxte / zgte
c
      zgnz = gnz
      if ( abs(zgnz) .lt. zepslon ) zgnz = zepslon
      diffnz = fluxnz / zgnz
c
      zgtz = gtz
      if ( abs(zgtz) .lt. zepslon ) zgtz = zepslon
      difftz = fluxtz / zgtz
c
c..test to see if ndim is at least 5
c
      if ( ndim .lt. 5 ) then
        nerr = 5
        return
      endif
c
      zdgth = 0.1
      zdgnh = 0.1
      zdgte = 0.1
      zdgnz = 0.1
      zdgtz = 0.1
c
c..derivative wrt gth
c
      zgth = gth + zdgth
c
      call sda04flx (lswitch, cswitch, lprint, nout
     & , gne, gnh, gnz, gte, zgth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , zfluxth, zfluxnh, zfluxte, zfluxnz, zfluxtz
     & , nerr )
c
      condif(1,1) = ( zfluxth - fluxth ) / zdgth
      condif(2,1) = ( zfluxnh - fluxnh ) / zdgth
      condif(3,1) = ( zfluxte - fluxte ) / zdgth
      condif(4,1) = ( zfluxnz - fluxnz ) / zdgth
      condif(5,1) = ( zfluxtz - fluxtz ) / zdgth
c
c..derivative wrt gnh
c
      zgnh = gnh + zdgnh
      zgne = zgnh * zfracnh - gnz * chargenz * fracnz + zgneut
c
      call sda04flx (lswitch, cswitch, lprint, nout
     & , zgne, zgnh, gnz, gte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , zfluxth, zfluxnh, zfluxte, zfluxnz, zfluxtz
     & , nerr )
c
      condif(1,2) = ( zfluxth - fluxth ) / zdgnh
      condif(2,2) = ( zfluxnh - fluxnh ) / zdgnh
      condif(3,2) = ( zfluxte - fluxte ) / zdgnh
      condif(4,2) = ( zfluxnz - fluxnz ) / zdgnh
      condif(5,2) = ( zfluxtz - fluxtz ) / zdgnh
c
c..derivative wrt gte
c
      zgte = gte + zdgte
c
      call sda04flx (lswitch, cswitch, lprint, nout
     & , gne, gnh, gnz, zgte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , zfluxth, zfluxnh, zfluxte, zfluxnz, zfluxtz
     & , nerr )
c
      condif(1,3) = ( zfluxth - fluxth ) / zdgte
      condif(2,3) = ( zfluxnh - fluxnh ) / zdgte
      condif(3,3) = ( zfluxte - fluxte ) / zdgte
      condif(4,3) = ( zfluxnz - fluxnz ) / zdgte
      condif(5,3) = ( zfluxtz - fluxtz ) / zdgte
c
c..derivative wrt gnh
c
      zgnz = gnz + zdgnz
      zgne = zgnh * zfracnh - gnz * chargenz * fracnz + zgneut
c
      call sda04flx (lswitch, cswitch, lprint, nout
     & , zgne, gnh, zgnz, gte, gth, gtz, th7te, tz7te
     & , fracnz, chargenz, amassz, fracns, chargens
     & , betae, vef, q, shear, kappa
     & , zfluxth, zfluxnh, zfluxte, zfluxnz, zfluxtz
     & , nerr )
c
      condif(1,4) = ( zfluxth - fluxth ) / zdgnz
      condif(2,4) = ( zfluxnh - fluxnh ) / zdgnz
      condif(3,4) = ( zfluxte - fluxte ) / zdgnz
      condif(4,4) = ( zfluxnz - fluxnz ) / zdgnz
      condif(5,4) = ( zfluxtz - fluxtz ) / zdgnz
c
c..derivative wrt gtz
c
      condif(1,5) = zero
      condif(2,5) = zero
      condif(3,5) = zero
      condif(4,5) = zero
      condif(5,5) = zero
c
c..convective velocities
c
      convel(1) = fluxth - condif(1,1) * gth - condif(1,2) * gnh
     &  - condif(1,3) * gte - condif(1,4) * gnz - condif(1,5) * gtz
c
      convel(2) = fluxnh - condif(2,1) * gth - condif(2,2) * gnh
     &  - condif(2,3) * gte - condif(2,4) * gnz - condif(2,5) * gtz
c
      convel(3) = fluxte - condif(3,1) * gth - condif(3,2) * gnh
     &  - condif(3,3) * gte - condif(3,4) * gnz - condif(3,5) * gtz
c
      convel(4) = fluxnz - condif(4,1) * gth - condif(4,2) * gnh
     &  - condif(4,3) * gte - condif(4,4) * gnz - condif(4,5) * gtz
c
      convel(5) = fluxtz - condif(5,1) * gth - condif(5,2) * gnh
     &  - condif(5,3) * gte - condif(5,4) * gnz - condif(5,5) * gtz
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      return
      end
!| \end{document}

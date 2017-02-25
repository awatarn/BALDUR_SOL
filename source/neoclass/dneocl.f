c  BALDUR  file DNEOCL   Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  To obtain this file, type           (use appropriate date for yymmdd)
c cfs get /11040/bald94/byymmdd/wcode.tar
c tar xf wcode.tar  (this creates directory .../code and subdirectories)
c cd code/bald      (this moves to subdirectory .../code/bald)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c       Contents of file DNEOCL:
c  TRNEO1 - simple neoclassical transport coefficients
c  NCFLUX - Hawryluk-Hirshman neoclassical particle transport
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@trneo1  .../baldur/bald/dneocl.f
c rgb 18-jul-01 replace function ftrapfn with array trapbr(jz,jb)
c   where trapbr(jz,jb) is computed in subroutine ftrapfn
c   note lneocl(3) controls trapped particle computation option
c   unless cfutz(481) > 0.
c rgb 11-jul-01 replaced function ft with ftrapfn
c rgb 17-jun-01 changed from dneo1 to dnneo1(j1,j2,jz)
c rap 29-may-01 Normalization of xineo1 is changed
c rap 02-may-00 Dimension of dneo1 is increased to account different 
c               contribution from different hydrogenic atoms
c  ap  15-feb-00 changed Hollerith character representation into '...'
c                in call mesage(...) and call error_olymp(...)
c  rgb 07-feb-00 removed equivalence statement
c rgb 20-mar-94 move common /cmneo1/ to file com/cbaldr.m
c rgb 14-jul-92 computed xeneo2(jz) if lneocl(3)=1 after Stringer 1991
c   les  nov-90 d3he; add coppi-sharky anomalous inward pinch
c           vxcmg=0 if lcmg=0
c   les  nov-90 ; d3he ions different
c   impurities have v-ware = veware ???
c   les 90 - note aired (defined in auxval/default) for imp=1-4
c rgb 20.15 07-oct-91  corrected problem with Ware pinch at r=0
c rgb 20.14 06-oct-91  corrections to zalpha in Chang-Hinton \chi_i
c rgb 20.11 27-aug-91  correct Chang-Hinton \chi_i and improve document
c rgb 20.06 05-aug-91  1986 Chang-Hinton implemented for impure plasmas
c rgb 20.05 28-jul-91  sbrtn trneo1 to compute neoclassical transport
c rgb 18.37 11-jun-90 ahalfs(jz,1) to ahalfs(jz,jb) in comp of zdramp
c     make sure zft(jz) = 0.0 at the magnetic axis, zft(jz)=0.
c     after do 1081
c  dps 16-dec-87 reinstate zshift specification by input for
c     Chang-Hinton  chi-i; used only if cfutz(282) > 0.
c  dps 07-may-87 introduce switch on trapped fraction calculation
c                     in ware pinch.
c  dps 04-mar-87 add bootstrap current cf. Mike Hughes; here
c                     install improved trapped fraction calculation.
c  rgb 16-mar-86 ramp trapping factor in ware pinch near axis
c  rgb 28-sep-85 cfutz(282) no longer used to define Shafranov
c                     shift in Chang-Hinton neoclassical ki
c                     instead zshift = d rmids(jz,2) / d ahalfs(jz,2)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c     Input variables:
c
c cfutz(9)   1.0  coeff of simplified neoclassical contribution to
c                   the diagonal hydrogenic paricle diffusivities
c cfutz(10)  1.0  coeff of simplified neoclassical contribution to
c                   the diagonal impurity   paricle diffusivities
c cfutz(11)  1.0  coeff of simplified neoclassical contribution to
c                   the electron thermal diffusivity
c cfutz(12)  1.0  coeff of simplified neoclassical contribution to
c                   the ion      thermal diffusivity
c cfutz(19)  1.0  coeff multiplying Ware pinch of hydrogen isotopes
c                   and electron energy
c cfutz(110) 0.0  activate Hawryluk-Hirshman ion-transport to replace
c                   the simplified off-diagonal neoclassical diffusivities
c cfutz(139) 0.0  multiplier for simplified off-diagonal neoclassical
c                   diffusivities
c cfutz(281) 0.0  switch for neoclassical ion thermal diffusivity model
c                 = 0.0 for Hinton-Hazeltine (1976)
c                 = 1.0 for Chang-Hinton (1982)
c                 = 1.2 for Chang-Hinton (1986) (effect of impurities)
c                 = 2.0 for Bolton-Ware (1983)
c cfutz(282) 0.0  fraction of minor radius that magnetic axis is
c                   shifted outward for computing \chi_i^{CH}
c                   ( = R_0' ).  Maximum value is 0.95.
c cfutz(365) 0.0  to include additional electron energy Ware pinch
c                   and omit convective energy flow controlled by CPVELC
c     (to restore the treatment of ware pinch on electron-energy flow
c      as used in older versions of the baldur code (prior to 19-may-83)
c      set "cfutz(365)=1.0".  default is "cfutz(365)=0.0".)
c cfutz(366) 0.0  multiplier for hydrogen-hydrogen simplified
c                   off-diagonal neoclassical transport diffusivities
c
c lneocl(3)  0  switch for simplified neoclassical contribution to
c                   the electron thermal diffusivity
c          = 0  for the original model in BALDUR
c          = 1  for multiple of neoclassical ion thermal diffusivity
c                 corresponding to theory by Stringer 1991
c
c NCLASS options:
c
c lneocl(1)      neoclassical model 
c         = 1    to call nclass
c         = else to call trneo1
c
c lneocl(2)      version of interface to nclass
c         = 9    simplified interface for testing
c
c lneocl(4)      option for output to nout (-)
c         > 0    errors only
c         = 2    errors and results
c         = else no output
c lneocl(11)      order of v moments to be solved (-)
c         = 2    u and p_q
c         = 3    u, p_q, and u2
c         = else 2
c lneocl(12)      option to include potato orbits (-)
c         = 0    off
c         = else on
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine trneo1
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cbparm.m'
      include 'comtrp.m'
      include 'cmg.m'
c
      dimension zk(21)      , zft(55)           , zaspin(55)
     &  , zifreq(idxion,55) , zgyro2(idxion,55) , znstar(idxion,55)
c
cgb      equivalence   (nlzzzz(1),zgyro2(1,1))
c
c
        ipmax  = mhyd + max0(mimp,0)
c
        z0     = 2.0 * (fxc-fxes) + fxau
        zgyro0 = ( 2.0*fcc*fcc*fcau ) * (10.0**z0) / (fces*fces)
c
      do 10 jz = lcentr, mzones
c
        j0 = jz
        if ( xbouni(j0) .le. epslon ) j0 = lcentr + 1
c
c  zaspin(jz) = half-width / major radius to midpoint at zone bndries
c
          zaspin(jz) = ahalfs(j0,1) / rmids(j0,1)
c
        do 6010 ip=1,ipmax
          ii     = ip - mhyd
          zavez2 = 1.0
          if ( ii .gt. 0 ) zavez2 = c2mean(ii,1,j0)
          zmass0 = sqrt( ahmean(1,j0) / (aspec(ip)+epslon) )
          zifreq(ip,jz) = zmass0 * zavez2 / tions(1,j0)
          zgyro2(ip,jz) = zgyro0 * ( aspec(ip) * tis(1,j0) ) /
     &                  ( zavez2 * bpols(1,j0)**2 )
          znstar(ip,jz) = zavez2 * xnuhyd(1,j0)
 6010   continue
c
  10  continue
c
c
      if ( cfutz(281) .lt. 0.1 ) then
c
      do 12 jz = lcentr, mzones
c
        if ( abs(xbouni(jz)) .gt. epslon ) then
c
c  ion-gyroradius factor
c
          zgyroi = (fcc*sqrt(2.*ahmean(1,jz)*fcmp)/fces)*
     &      10.**(fxc+.5*fxnucl-fxes)
c
          z0 = 0.66 * sqrt(zaspin(jz)) * zgyroi**2 * tis(1,jz)
          z1 = 1.0 + 1.03 * sqrt( xnuhyd(1,jz) + 0.31 * xnuhyd(1,jz) )
          z2 = 1.0 + 0.74 * ( zaspin(jz) * sqrt(zaspin(jz)) )
     &         * xnuhyd(1,jz)
          z3 = 1.76645 * ( zaspin(jz) * sqrt(zaspin(jz)) )
     &           * ( zaspin(jz) * sqrt(zaspin(jz)) ) * xnuhyd(1,jz)
          xineo1(jz) = ( z0 * (z2+z1*z3) ) /
     &      ( z1 * z2 * bpols(1,jz)**2 * tions(1,jz) + epslon )
c
        else
c
          xineo1(jz) = ( 4.25806 * fcc**2 *
     &      sqrt( ahmean(1,jz) * fcmp * tis(1,jz) ) *
     &      tis(1,jz) * ahalfs(jz+1,1) ) /
     &      ( fces**2 * rmajs**2 * bzs * bpols(1,jz+1) ) *
     &      ( 10.0**( 2.0 * (fxc-fxes) + 0.5 * fxnucl ) )
c
        endif
c
  12  continue
c
c
      elseif ( cfutz(281) .lt. 1.9 ) then
c
        zk20 = 0.66
        za20 = 1.03
        zb20 = 0.31
        zc20 = 0.74
c
      do 14 jz = lcentr, mzones
c
        zalpha = 0.0
        zdenh  = 0.0
        zdenzi = 0.0
      if ( cfutz(281) .gt. 1.1  .and.  mimp .gt. 0 ) then
        do 23 jh=1,mhyd
          zdenh = zdenh + rhohs(jh,1,jz)
  23    continue
        do 24 ji=1,mimp
          zdenzi = zdenzi + rhois(ji,1,jz) * c2mean(ji,1,jz)
  24    continue
          zalpha = zdenzi / zdenh
      endif
c
        zshift=0.0
      if (abs(cfutz(282)) .gt. epslon ) then
        zshift = sign ( min(abs(cfutz(282)),.95) , cfutz(282) )
      else
        zshift = (rmids(jz,2)-rmids(jz-1,2))
     &              / (ahalfs(jz,2)-ahalfs(jz-1,2))
      endif
c
        zb1 = ( 1.0 + 1.5 * ( zaspin(jz)**2 + zaspin(jz) * zshift )
     &           + 0.375 * zaspin(jz)**3 * zshift )
     &             / ( 1.0 + 0.5 * zaspin(jz) * zshift )
c
        zb2  = zaspin(jz) * sqrt(1.0-zaspin(jz)**2)
     &    * (1.0 + 0.5*zaspin(jz)*zshift)
     &    / ( zaspin(jz) + zshift*( sqrt(1.0 - zaspin(jz)**2) - 1.0 ) )
c
        zf0 = ( zb1 - zb2 ) / ( 2.0 * sqrt(zaspin(jz)) )
c
        zhp = 1.0 + 1.33*zalpha*(1.0+0.60*zalpha)/(1.0+1.79*zalpha)
c
        zk2hat = zb1 * ( 0.66 * (1.0 + 1.54*zalpha)
     &    + (1.88*sqrt(zaspin(jz)) - 1.54*zaspin(jz))
     &      * (1.0 + 3.75*zalpha) )
c
          xineo1(jz) = 0.0
c
        do 6030 ip=1,ipmax
c
          zmustr = max( abs(znstar(ip,jz) * ( 1.0 + 1.54 * zalpha ))
     &      / ( 1.0 + zalpha) , epslon )
c
          zk2 = zk2hat / ( 1.0 + za20 * sqrt(zmustr) + zb20 * zmustr )
     &          + zk20 * ( zc20**2/zb20 ) * zmustr
     &            * zaspin(jz)*sqrt(zaspin(jz)) * zhp * zf0
     &            / ( 1.0 + zc20 * zmustr*zaspin(jz)*sqrt(zaspin(jz)) )
c
 
          if ( ip .gt. mhyd   ) then
             zdens = rhois(ip - mhyd,1,jz)
          else
            zdens = rhohs(ip,1,jz)
          endif
c
          xineo1(jz) = xineo1(jz)
     &      + zdens * sqrt(zaspin(jz)) * zgyro2(ip,jz)
     &        * zk2 * zifreq(ip,jz) / ( 1.0 + zalpha )
 6030   continue
cap
        xineo1(jz) = xineo1(jz) / rhoins(1, jz) 
c
  14  continue
c
c
      else
c
      do 16 jz = lcentr, mzones
c
        zanc = ( .66 + 2.441 * sqrt(zaspin(jz))
     &    -3.87 * zaspin(jz) + 2.19 * (zaspin(jz)*sqrt(zaspin(jz))) )
        zbnc = ( 0.9362 - 3.109 * zaspin(jz)
     &    + 4.087 * zaspin(jz)**2 )
     &      / sqrt( zaspin(jz) * sqrt(zaspin(jz)) )
        zcnc = ( 0.241 + 3.40 * zaspin(jz)
     &    - 2.54 * zaspin(jz)**2 ) / (zaspin(jz)*sqrt(zaspin(jz)))
        zdnc = ( 0.2664 - 0.352 * zaspin(jz)
     &    + 0.44 * zaspin(jz)**2 ) / (zaspin(jz)*sqrt(zaspin(jz)))
c
        zaps = ( 0.364 - 2.76 * zaspin(jz)
     &    + 2.21 * zaspin(jz)**2 ) * (zaspin(jz)*sqrt(zaspin(jz)))
        zbps = ( 0.553 + 2.41 * zaspin(jz)
     &    - 3.42 * zaspin(jz)**2 ) * (zaspin(jz)*sqrt(zaspin(jz)))
        zcps = ( 1.18 + 0.292 * zaspin(jz) + 1.07 * zaspin(jz)**2 )
        zdps = ( 0.0188 + 0.180 * zaspin(jz) - 0.127 * zaspin(jz)**2 )
c
        xineo1(jz) = 0.0
        za   = 1.3293404 * (zaspin(jz)*sqrt(zaspin(jz)))
c
        do 6020 ip=1,ipmax
          ii  = ip - mhyd
          z0  = za * znstar(ip,jz)
          z1  = sqrt(z0)
          z2  = z0 * z1
          zk1 = zanc / ( 1.0 + zbnc * z1 + zcnc * z0 + zdnc * z0 * z0 )
          zk2 = 1.57 * (zaspin(jz)*sqrt(zaspin(jz)))
     &      + (zaps + zbps*z1) / (1.0 + zcps*z2 + zdps*z0*z2 )
          zka = zk1 + zk2
c
          if ( ip .gt. mhyd   ) then
            zdens = rhois(ii,1,j0)
          else
            zdens = rhohs(ip,1,j0)
          endif
c
          xineo1(jz) = xineo1(jz)
     &      +  sqrt(zaspin(jz)) * zka * zdens
     &           * zgyro2(ip,jz) * zifreq(ip,jz)
 6020   continue
cap
        xineo1(jz) = xineo1(jz) / rhoins(1, jz)
  16  continue
c
      endif
c
c
      do 28 jz = lcentr, mzones
c
        if ( xnuhyd(1,jz) .lt. epsinv ) then
c
c..neither b-poliodal or r are near 0
c
      xeneo1(jz) = cdetes * max ( q(jz)**2, 1.0 )
     &  * cloges(1,jz) * rhoels(1,jz) /
     1  (bzs**2 * sqrt(tes(1,jz))) * (0.73*(1.6 + 2.0*xzeff(1,jz)) /
     2  (1.0 + 0.1*(1.6 + 2.0*xzeff(1,jz))
     3  * xnuel(1,jz) ) * sqrt( rmids(jz,1) / ahalfs(jz,1) )**3 +
     4  1.13 + 0.50*xzeff(1,jz) + 0.55*xzeff(1,jz)/(0.59 + xzeff(1,jz)))
c
        else
c
c..either b-poloidal or r are near 0
c
c  zbr is b/r
c  note--it is assumed that b-poloidal is 0, although
c  r need not be.  i.e., if r=0, b-poloidal must be 0
c
      zbr = bpols(1,jz+1) / (rmins * xbouni(jz+1))
c
c  zetae = r**(3/2) * zbr * xnuel(1,jz)
c
c  i.e., zetah and zetae are xnuhyd and xnuel
c  without the b-poloidal and r**(1/2) terms
c
      zetae = cnuel * rmajs * bzs / telecs(1,jz) *
     &          sqrt(rmajs / tes(1,jz))
c
      xeneo1(jz) = cdetes * q(jz)**2 * cloges(1,jz) * rhoels(1,jz) /
     1          (bzs**2 * sqrt(tes(1,jz))) *
     2  (zbr * 0.73*(1.6 + 2.0*xzeff(1,jz)) / (zetae * 0.1*(1.6 +
     3  2.0*xzeff(1,jz))) * sqrt(rmajs)**3 +
     4  1.13 + 0.50*xzeff(1,jz) + 0.55*xzeff(1,jz)/(0.59 + xzeff(1,jz)))
c
        endif
c
  28  continue
c
c
      call resetr ( xeneo2, mzones, 0.0 )
c
      if ( lneocl(3) .eq. 1 ) then
        do 29 jz = lcentr, mzones
          ztau = tes(1,jz) / tis(1,jz)
cap	  
          xeneo2(jz) = rhoins(1, jz) * xineo1(jz) * ztau**3 / 
     &	   ((1.0 + ztau) * (3.0 + ztau**2 / (1.0 + ztau)))
  29    continue
      endif
c
c
      do jz = lcentr, mzones
c
        do j1=1,mhyd+mimp
          do j2=1,mhyd+mimp
            dnneo1(j1,j2,jz) = cdnhs * rhoels(1,jz) *  q(jz)**2 *
     &          cloges(1,jz) * gspitz(1,jz) * (tes(1,jz) + tis(1,jz)) /
     &          (sqrt(tes(1,jz))**3 * bzs**2)
          enddo
        enddo
c
      enddo
c
c
c
      if ( cfutz(139) .lt. epslon ) go to 1180
      if ( cfutz(110) .gt. epslon ) go to 1180
c
c
        cdnh0=cfutz(139)*cdnhis
        cdni0=cfutz(139)*cdniis
c
      if ( mimp .gt. 0 ) then
c
        do 40 jz = lcentr, mzones
c
          do 1120 jh = 1, mhyd
c
             zdnhis = cdnh0 * q(jz)**2 * clogis(1,jz)
     &                 * sqrt(aspec(jh) / tis(1,jz)) / bzs**2
c
             do 1120 ji = 1, mimp
               dnhis(jh,ji,jz) = zdnhis
               dnihs(ji,jh,jz) = zdnhis
 1120        continue
c
  40    continue
c
      endif
c
c
      do 1160 jh = 1, mhyd
      do 1160 jh2 = 1, mhyd
c
      if ( jh2 .eq. jh ) then
c
        do 42 jz = lcentr, mzones
          dnhhs(jh,jh,jz)=0.0
  42    continue
c
      else
c
        do 44 jz = lcentr, mzones
          dnhhs(jh,jh2,jz) = cfutz(366) * cdni0
     &         * clogis(1,jz) * q(jz)*2
     &         / ( bzs**2 *  sqrt( tis(1,jz) /
     &      (aspec(jh)*aspec(jh2))/(aspec(jh)+aspec(jh2)) )  )
  44    continue
c
      endif
c
 1160 continue
c
c
c  ...For comparison with BALDP86m, needed to use ".le.0" instead of
c  ...".lt.2" in the following line.
c
      if ( mimp .lt. 2 ) go to 1180
c
      do 1170 ji = 1, mimp
      do 1170 ji2 = 1, mimp
c
c  ...Also, for comparison with BALDP86m, remove this if-then-else structure,
c  ...leaving the non-trivial assignment statement (this was a bug).
c
        if ( ji2 .eq. ji ) then
c
          do 46 jz = lcentr, mzones
            dniis(ji,ji2,jz) = 0.0
  46      continue
c
        else
c
          do 48 jz = lcentr, mzones
            dniis(ji,ji2,jz) = cdni0 * clogis(1,jz) * q(jz)**2 /
     &               (sqrt(tis(1,jz)/aired(ji,ji2)) * bzs**2)
  48      continue
c
        endif
c
 1170 continue
c
 1180 continue
c
c
      do 52 jz = lcentr, mzones
c
c       the following relations are taken from an analytic approx-
c       imation for the ware flux worked out by s. p. hirshman and
c       r. j. hawryluk
c
c       the first quantities defined are used in assembling the
c       analytic approximations.  the characters after "z" cor-
c       respond to hirshman and hawryluk's notation
c
C.MHH Dec. 23 1986 - require these coefficients at both cell centres
c     (bootstrap current) and cell boundaries (pinch effect).
c
      do 1081 jb=2,1,-1
c
        zft(jz) = 0.0
c
c..14.03: should only skip calculation of rl13 and rl23 when jb=1
c..and jz=lcentr - allows zone center values at jz=lcentr to be found.
c
      if ( bpols(jb,jz) .le. epslon ) then
        rl13(jz,jb) = 0.0
        rl23(jz,jb) = 0.0
        rly(jz,jb)  = 3.0
      else
c
      zk1    = (0.53+xzeff(jb,jz)) /
     &  (xzeff(jb,jz)*(1.0+1.32*xzeff(jb,jz)))
      zcr1   = zk1 / (1.0+zk1)
c
      zk2    = (0.82*(3.0-xzeff(jb,jz)))/
     &  (xzeff(jb,jz)*(2.57+xzeff(jb,jz)))
      zcr2   = zk2 / (2.5+zk2)
c
      z1beta = 0.4  + 0.11  * xzeff(jb,jz)
      z1zeta = 0.55 + 0.255 * xzeff(jb,jz)
      z2beta = 0.05 + 0.0345* xzeff(jb,jz)
      z2zeta = 0.25 + 0.143 * xzeff(jb,jz)
c
c  zft(jz) is the trapped fraction
c
      zd = abs ( ahalfs(jz,jb) / rmids(jz,jb) )
c
c...Replace this assignment with the following for comparison with
c      BALDP86m.
c      zd = abs ( ahalfs(jz,2) / rmids(jz,2) )
c..........
c
c..smooth cutoff for zft(j), bateman, 16-feb-86
c  use ramp function to smooth ftrap near magnetic axis
c  zdramp = 1.0 for r/a > cemprc(20)
c  zdramp = r / (a * cemprc(20)) for r/a < cemprc(20)
c
      itrap = lneocl(3)
      if ( cfutz(481) - 0.1 .gt. 0.0 )
     &     itrap = int ( cfutz(481) + 0.1 )
c
      if ( itrap .gt. 0 ) then
c
        zft(jz) = trapbr(jz,jb)
c
      else
c
        zdramp = min ( 1.0 ,
     & ahalfs(jz,jb)/(ahalfs(mzones,1)*max(cemprc(20),epslon) ) )
c
        zft(jz) = 1.0 -  (1.0 - zd)**2
     &  / (sqrt(abs(1.0 - zd**2))*(1.0 + 1.46*sqrt( abs(zd*zdramp))))
c
      endif
c
c
      zk11 = zft(jz)
     &  / ( sqrt(z1beta*xnuel(jb,jz))+z1zeta*xnuel(jb,jz)+1.0 )
      zk12 = zft(jz)
     &  / ( sqrt(z2beta*xnuel(jb,jz))+z2zeta*xnuel(jb,jz)+1.0 )
c
c..Introduce factors to "speed up" transition out of banana regime
c..near the axis; i.e., increase rate with which zl13 and zl23 approach
c..0 as r -> 0, analogous to procedure used on neoclassical resistivity
c..in subroutine getchi.
c
            ztrans=1.
        if ((ctrnsp(19).gt.epslon).and.(ctrnsp(20).gt.epslon)) then
          if (jb.eq.1)
     &      ztrans=1./(1.+(ctrnsp(19)/xbouni(jz))**ctrnsp(20))
          if (jb.eq.2)
     &      ztrans=1./(1.+(ctrnsp(19)/xzoni(jz))**ctrnsp(20))
        end if
c
c       principal factors in the analytic expressions
c
C.MHH Dec 23 1986 - save coefficients - use RHH expression for *y*
c
      rl13(jz,jb) = (1.0+zk1) * zk11 * (1.0-zcr1*zk11) * ztrans
      rl23(jz,jb) = (2.5+zk2) * zk12 * (1.0-zcr2*zk12) * ztrans
      rly(jz,jb)  = (1.33+3.0*xnuhyd(jb,jz))/(1.0+xnuhyd(jb,jz))
c
      endif
c
 1081    continue
  52  continue
c
c
      do 61 jz=1,mzones
        vnwars(jz) = 0.0
  61  continue
c
      if ( cfutz(19) .gt. epslon ) then
c
        zcc = fcc * 10.0**fxc
c
        do 62 jz = 3, mzones
c
c..Rewrite Ware pinch in terms of loop voltage
c
        if (versno.gt.14.01) then
          zk15 = uisv * 0.5 * (vloopi(jz,2)+vloopi(jz-1,2))
     &             * zcc / ( 2.*fcpi*r0ref*bpols(1,jz) )
        else
          zk15 = zcc * eta(1,jz) * ajzs(1,jz) / bpols(1,jz)
        end if
c
        vnwars(jz) = cfutz(19) * rl13(jz,1) * zk15
c
        vewars(jz) = cfutz(19) * rl23(jz,1) * zk15 * cfutz(365)
c
  62  continue
c
      endif
c
c     add ware-pinch effects to the generalized pinch velocities
c      (ware-pinch effects on impurities are presently ignored)
c   les  nov-90 d3he; add coppi-sharky anomalous inward pinch
c           vxcmg=0 if lcmg=0
c
      do 64 jh=1,mhyd
        do 63 jz = 1, mzones
          vxemps(jh,jz) = vxemps(jh,jz) + vnwars(jz) - vxcmg(jz)
  63    continue
 64   continue
c
c   les  nov-90 ; d3he ions different
c   impurities have v-ware = veware ???
c
      if ( limpn .gt. mhyd ) then
        do 1103 ji=mhyd+1,limpn
          vxemps(ji,jz) = vxemps(ji,jz) - vxcmg(jz)
 1103   continue
      endif
c
c  les  protons have vnwars
c
      if ( cfutz(490) .gt. epslon )
     &   vxemps(lprotn,jz) = vnwars(jz) + vxemps(lprotn,jz)
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@ncflux  /11040/bald91/wbaldn1 DNEOCL
c  rgb 02-jun-96 remove if(.not.inital) go to 1000
c    save zxtaa, zxtexp, zxlar, zlarex, zxomeg, zxlar2
c       les  nov-90  reorder ion species for d-3he fusion
c       fgps 1-feb-83 revised to handle more than 4 ionic species
c       aes 15-apr-82 fix form feeds in ncfprt formats
c       aes 12-mar-82 make times printed in ncfprt same as in mprint
c       aes 20-jan-82 added label5 to page headers
c       aes 19-nov-81 add 'data incflx' --> control logic for ncfprt
c       aes 19-nov-81 move error lines 9000 to before entry ncfprt
c               -- change label to 8599; clean up labeling order
c       aes 17-nov-81 edit printout --> entry ncfprt
c       aes 29-oct-81 array dimensions 52 -->55 in common/tempry,neoion/
c       aes 28-oct-81 put all data statements before executable code
c       fgps 20-jul-79 adapted subroutine ncflux (hawryluk 1-jun-
c                      79) for inclusion in baldur
c
c**********************************************************************c
c
        subroutine ncflux
c
c
cl      2.22    neoclassical particle transport redefined
c               according to hawryluk and hirshman
c
c
c**********************************************************************c
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cbparm.m'
c
c       do not equivalence arrays; but, when there are 4 impurity
c       species, mxions = 6.
c
c------------------------------------------------------------------------------
c       june 1,1979
c
c       this routine calculates the neoclassical impurity fluxes
c       as well as the hydrogenic fluxes. electron-ion fluxes are
c       not calculated.
c       no mass ratio expansion is used in these calculations.
c       this calculation is based upon pppl-1473.       -appendix a
c
      dimension
     1zmass(idxion)        ,zcharg(idxion)       ,zdens(idxion)        ,
     2zab(idxion,idxion)   ,zc11(idxion,idxion)  ,zc12(idxion,idxion)  ,
     3zxc11(idxion,idxion) ,zxc12(idxion,idxion) ,zxb1b(idxion,idxion) ,
     4zxb2b(idxion,idxion) ,zxi11(idxion,idxion) ,zxi01(idxion,idxion) ,
     5zxi00(idxion,idxion) ,ztaa(idxion)         ,zlar2(idxion)        ,
     6zmu1(idxion)         ,zmu2(idxion)
c
c       june 8,1979
c       the following arrays are required for the pfirsch-schluter fluxes
c       the calculations are based on the work by
c       s.p.hirshman ,phys.fluids, pg 589(1977)
c
      dimension
     1zxn00(idxion,idxion) ,zxn01(idxion,idxion) ,zxn02(idxion,idxion) ,
     2zxn11(idxion,idxion) ,zxn12(idxion,idxion) ,zxn22(idxion,idxion) ,
     3zxm00(idxion,idxion) ,zxm01(idxion,idxion) ,zxm02(idxion,idxion) ,
     4zxm11(idxion,idxion) ,zxm12(idxion,idxion) ,zxm22(idxion,idxion) ,
     5zxhm22(idxion,idxion),zl11(idxion,idxion)  ,zxmmnn(idxion,idxion),
     6zl12(idxion,idxion)  ,zl22(idxion,idxion)  ,zm00(idxion)         ,
     7zm01(idxion)         ,zm11(idxion)         ,zm12(idxion)         ,
     8zalpha(idxion,idxion),zbeta(idxion,idxion) ,zdelta(idxion,idxion),
     9zhalph(idxion,idxion),zhbeta(idxion,idxion),zxa22(idxion,idxion) ,
     &zm20(idxion)         ,zm21(idxion)         ,zvstr(idxion)        ,
     1zw(idxion)           ,zm22(idxion)         ,zm02(idxion)
c
c       the following arrays are related to the strong
c       temperature equilibration approximation
c
      dimension
     1zc2(idxion,idxion)   ,zlh22(idxion)        ,zin2(idxion,idxion)  ,
     2zin22(idxion,idxion) ,zch2(idxion)         ,zs11(idxion,idxion)  ,
     3zs12(idxion,idxion)  ,ziii(idxion,idxion)  ,zzl22(idxion,idxion) ,
     4ipivot(idxion)
c
c
        logical         inital       , ladjst
c
c
c
        data    iclass /2/,     isub /22/
c
cahk
        save zxtaa, zxtexp, zxlar, zlarex, zxomeg, zxlar2
c
        if (.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.22 subroutine ncflux bypassed ')
        return
   10   continue
c
c------------------------------------------------------------------------
c
c       cfutz(incflx) .ne. 0 allows subroutine convrt to call
c       subroutine ncflux.  arguments of cfutz-factors which ad-
c       just the levels of regimes are:  iclflx for classical,
c       ibpflx for banana-plateau, and ipsflx for pfirsch-
c       schluter.
c
        data    incflx,iclflx,ibpflx,ipsflx /110,111,112,113/
c
        data    inital,ladjst /.true.,.false./
c
c
cbate        if(.not.inital) go to 1000
        inital=.false.
        if(cfutz(iclflx).lt.epslon) cfutz(iclflx)=1.0
        if(cfutz(ibpflx).lt.epslon) cfutz(ibpflx)=1.0
        if(cfutz(ipsflx).lt.epslon) cfutz(ipsflx)=1.0
        if(cfutz(iclflx).ne.1.0) ladjst=.true.
        if(cfutz(ibpflx).ne.1.0) ladjst=.true.
        if(cfutz(ipsflx).ne.1.0) ladjst=.true.
c
c       initialize several constants used in the calculation.
c
        zxtaa=3.*sqrt(fcau)/(4.*sqrt(2.*fcpi)*fces**4)
        zxtexp=0.5*fxnucl-4.*fxes
        zxtaa=zxtaa*10.**zxtexp
        zxlar=fcc*sqrt(2.*fcau)/fces
        zlarex=fxc+0.5*fxnucl-fxes
        zxlar=zxlar*10.**zlarex
        zxomeg=sqrt(2./fcau)*10.**(-0.5*fxnucl)
        zxlar2=zxlar*zxlar
c
c       initialize several arrays which depend only upon
c       the atomic mass of the ions.
c
c       limpn=number of ion species
c       mhyd=number of hydrogen species
c       mimp=number of impurity species
c       aspec(ih)=atomic weight of hydrogen, ih=1,mhyd
c       aspec(ii)=atomic weight of impurities, ii=mhyd+1,limpn
c       cmean=mean z of the impurities--the difference between
c               the mean z and the mean z*z has been ignored
c
c       limit set equal to number of active ion species
c
        ihmax=limpn
        if(mimp.le.0) ihmax=mhyd
c
c       inverse aspect ratio at outer! edge of scrapeoff region
c       if such exists, otherwise at outer edge of plasma
c
        zasp0=rmins/rmajs
c
c       store at. wt. of ionic species into zmass(*)
c
c   les  nov-90  d3he fusion -- reorder ions by atomic weigth
c
      ipstrt=1
      if ( cfutz(490).gt.epslon ) then
        zmass(1)=aspec(lprotn)
        zmass(2)=aspec(ldeut)
        zmass(3)=aspec(ltrit)
        ipstrt=4
      endif
c
      do 200 ip=ipstrt,ihmax
        zmass(ip)=aspec(ip)
 200  continue
c
        do 300 ih=1,ihmax
        zxc11(ih,ih)=0.
        zxc12(ih,ih)=0.
c
        do 350 ihh=1,ihmax
c
        zx2=(zmass(ih)/zmass(ihh))
        zxa22(ih,ihh)=1./zx2
        zx4=zx2*zx2
        zx6=zx4*zx2
        zx8=zx6*zx2
        zx=sqrt(zx2)
        zxa1=1./(1.+zx2)
        zx12=sqrt(zxa1)
        zxa12=1./zx12
        zx32=zx12*zxa1
        zx52=zx32*zxa1
        zx72=zx52*zxa1
        zx92=zx72*zxa1
c
c       calculate the constant factor which is part of the classical
c       fluxes  (eq. a5a,b)
c       t(ih)=t(ihh) is assumed throughout.
c
        if(ih .eq. ihh) go to 220
        zxc11(ih,ihh)=zx12
        zxc12(ih,ihh)=-1.5*zx2*zx32
220     continue
c
c       to calculate the anisotropy driven fluxes, several arrays
c       must be calculated-eq a11.
c
c
c       see eq a11a and a11b.
c
        zxb1b(ih,ihh)=zxa12+zx2*log(zx/(1.+zxa12))
        zxb2b(ih,ihh)=zx12
c
c       see eq a12,a13,a14
c
        zxi00(ih,ihh)=(3.+5.*zx2)*zx32
        zxi01(ih,ihh)=(4.5+10.5*zx2)*zx52
        zxi11(ih,ihh)=(35.*zx6+38.5*zx4+46.25*zx2+12.75)*zx72
c
c       the constant coefficients required for the p-s
c       fluxes will be calculated:(sph equation a3)
c
c       the n00 etc. arrays are divided by zx compared with sph
c
        zxm00(ih,ihh)=-zx12
        zxn00(ihh,ih)=zx*zx12
        zxm01(ih,ihh)=-1.5*zx32
        zxn01(ihh,ih)=1.5*zx*zx32
        zxm02(ih,ihh)=-1.875*zx52
        zxn02(ihh,ih)=1.875*zx*zx52
        zxm11(ih,ihh)=-zx52*(7.5*zx4+4.*zx2+3.25)
        zxm12(ih,ihh)=-zx72*(15.75*zx4+6.*zx2+4.3125)
        zxm22(ih,ihh)=-zx92*
     .  (21.875*zx8+28.*zx6+57.375*zx4+17.*zx2+6.765)
        zxn11(ih,ihh)=6.75*zx52*zx2
        zxn12(ih,ihh)=14.0625*zx72*zx4
        zxn22(ih,ihh)=41.016*zx92*zx4
c
c       see equations a19 in ppl 1473
c
        zxhm22(ih,ihh)=-zxm22(ih,ihh)
        zxmmnn(ih,ihh)=3.*zx72*(2.5*zx4+2.*zx2+3.25)*zx2
350     continue
        zxi00(ih,ih)=zxi00(ih,ih)-.707107
        zxi01(ih,ih)=zxi01(ih,ih)-1.06066
        zxi11(ih,ih)=zxi11(ih,ih)-2.65165
c
c       see equations a19 in ppl 1473
c
        zxhm22(ih,ih)=zxhm22(ih,ih)-1.81
        zxmmnn(ih,ih)=zxmmnn(ih,ih)-.33145
c
300     continue
c
c
1000    continue
c
c-----------------start of main do-loop over radial index j---------------
c
        do 10000 j=lcentr,mzones
c
c
c       create a local array composed of the hydrogenic and impurity
c       species and their charge
c
c   les  nov-90  for d3he, reorder ions by atomic weight
c
      if ( cfutz(490).le.epslon ) go to 1009
        zdens(1)=rhois(lprotn-lhydn,1,j)
        zdens(2)=rhohs(ldeut,1,j)
        zdens(3)=rhohs(ltrit,1,j)
        zcharg(1)=1.
        zcharg(2)=1.
        zcharg(3)=1.
        do 1002 ii=4,ihmax
          zdens(ii)=rhois(ii-lhydn,1,j)
 1002  zcharg(ii)=cmean(ii-lhydn,1,j)
      go to 1220
 1009  continue
c
        do 1100 ih=1,mhyd
        zdens(ih)=rhohs(ih,1,j)
        zcharg(ih)=1.
1100    continue
c
c
        if(mimp.le.0) go to 1220
        do 1200 ih=1,mimp
        ii=ih+lhydn
        zdens(ii)=rhois(ih,1,j)
        zcharg(ii)=cmean(ih,1,j)
1200    continue
1220    continue
c
c       calculate the effective charge
c
        do 1300 ih=1,ihmax
        do 1300 ihh=1,ihmax
1300    zab(ih,ihh)=zcharg(ihh)*zcharg(ihh)*zdens(ihh)/
     1  (zcharg(ih)*zcharg(ih)*zdens(ih))
c
c       calculate the mu coefficients eq a8-a9
c
        zmsum=0.
        zeps=xbouni(j)*zasp0
        if(zeps .le. 1.e-10) zeps=1.e-10
        zeps32=zeps**1.5
        z1taa=zxtaa*tis(1,j)**1.5/clogis(1,j)
        z1lar2=zxlar2*tis(1,j)/(bzs*bzs)
        z1omeg=zxomeg*sqrt(tis(1,j))/(q(j)*rmajs+epslon)
c
c
        do 2000 ih=1,ihmax
c
c
        ztaa(ih)=z1taa*sqrt(zmass(ih))/(zdens(ih)*zcharg(ih)**4)
        zlar2(ih)=z1lar2*zmass(ih)/(zcharg(ih)*zcharg(ih))
        zomeg=z1omeg/sqrt(zmass(ih))
        zvstr(ih)=1./(ztaa(ih)*zomeg*zeps32)
c
c       see eq a11-a13
c
        zb1b=0.
        zb2b=0.
        zi11=0.
        zi01=0.
        zi00=0.
        do 2500 ihh=1,ihmax
        zb1b=zb1b + zab(ih,ihh)*zxb1b(ih,ihh)
        zb2b=zb2b + zab(ih,ihh)*zxb2b(ih,ihh)
        zi11=zab(ih,ihh)*zxi11(ih,ihh)+zi11
        zi01=zab(ih,ihh)*zxi01(ih,ihh)+zi01
        zi00=zab(ih,ihh)*zxi00(ih,ihh)+zi00
2500    continue
c
        zb1p=.607/zvstr(ih)
        zb2p=3.*zb1p
        ziot1=2.5*zi11/(zi11*zi00 -zi01*zi01)
        ziot2=2.5*(zi11+3.5*zi01)/(zi11*zi00-zi01*zi01)
        zden=1./(zvstr(ih)*zvstr(ih)*zeps32)
        zb1ps=.514*ziot1*zden
        zb2ps=.514*(ziot2+2.5*ziot1)*zden
        zb1=zb1b/((1.+zb1b/zb1p)*(1.+zb1p/zb1ps))
        zb2=zb2b/((1.+zb2b/zb2p)*(1.+zb2p/zb2ps))
c
c       normalize the original mu matrix coefficients by the proton mass.
c
        zmu1(ih)=zdens(ih)*zmass(ih)*zb1/ztaa(ih)
        zmu2(ih)=zmu1(ih)*(zb2-2.5*zb1)/zb1
        zmsum=zmsum+zmu1(ih)
c
2000    continue
c
c
c       calculate the arrays required for the classical
c       coefficients
c
        do 2600 ih=1,ihmax
        zc11(ih,ih)=0.
        zc12(ih,ih)=0.
        do 2650 ihh=1,ihmax
        if(ih .eq. ihh) go to 2650
        zc11(ih,ihh)=zab(ih,ihh)*zxc11(ih,ihh)
        zc12(ih,ihh)=zab(ih,ihh)*zxc12(ih,ihh)
        zc11(ih,ih)=zc11(ih,ih)-zc11(ih,ihh)
        zc12(ih,ih)=zc12(ih,ih)-zc12(ih,ihh)*zxa22(ih,ihh)
2650    continue
2600    continue
 
c
c       now both the classical (cl11 and cl12)
c       and the anisotropy drive coefficients (bp11,bp12) can be calculated
c
        zb1=.73*q(j)*q(j)/zeps32
        do 3000 ih=1, ihmax
        za1=0.5*zlar2(ih)*zdens(ih)*zcharg(ih)/ztaa(ih)
        do 3500 ihh=1,ihmax
c
c       classical coefficients
c
        za=za1/zcharg(ihh)
        cl11(ih,ihh,j)=za*zc11(ih,ihh)
        cl12(ih,ihh,j)=za*zc12(ih,ihh)
c
c       anisotropy drive fluxes
c
        za=zb1*zcharg(ihh)/zcharg(ih)
     1  *zlar2(ihh)/zmass(ihh)
        z1=zmu1(ih)/zmsum
        if(ih .eq. ihh) z1=z1-1.
c
        bp11(ih,ihh,j)=za*zmu1(ihh)*z1
        bp12(ih,ihh,j)=za*zmu2(ihh)*z1
c
3500    continue
3000    continue
c
c
c       we will now calculate the friction coefficients
c       in the pfirsch-schluter regime.
c
 
c       calculate the m arrays --sph equation 15.
c
        do 4000 ih=1,ihmax
        zm00(ih)=0.
        zm01(ih)=0.
        zm11(ih)=0.
        zm12(ih)=0.
        zm02(ih)=0.
        zm22(ih)=0.
c
        z1=zdens(ih)/ztaa(ih)
        do 4100 ihh=1,ihmax
        zm00(ih)=zm00(ih)+zab(ih,ihh)*zxm00(ih,ihh)
        zm11(ih)=zm11(ih)+zab(ih,ihh)*zxm11(ih,ihh)
        zm12(ih)=zm12(ih)+zab(ih,ihh)*zxm12(ih,ihh)
        zm22(ih)=zm22(ih)+zab(ih,ihh)*zxm22(ih,ihh)
        zm02(ih)=zm02(ih)+zab(ih,ihh)*zxm02(ih,ihh)
        zm01(ih)=zm01(ih)+zab(ih,ihh)*zxm01(ih,ihh)
4100    continue
        zm00(ih)=zm00(ih)*z1
        zm11(ih)=zm11(ih)*z1
        zm12(ih)=zm12(ih)*z1
        zm22(ih)=zm22(ih)*z1
        zm02(ih)=zm02(ih)*z1
        zm01(ih)=zm01(ih)*z1
4000    continue
c
        do 4150 ih=1,ihmax
        zm20(ih)=0.
        zm21(ih)=0.
        do 4175 ihh=1,ihmax
        z1=zdens(ihh)/ztaa(ihh)
        zm20(ih)=zm20(ih)+z1*zab(ihh,ih)*zxm02(ihh,ih)
        zm21(ih)=zm21(ih)+z1*zab(ihh,ih)*zxm12(ihh,ih)
4175    continue
4150    continue
c
c       calculate the alpha,beta,and delta arrays
c       see sph equation 14
c
        do 4200 ih=1,ihmax
        z1=zdens(ih)/ztaa(ih)
        zden=1./(zm22(ih)+z1*zxn22(ih,ih))
c
        do 4250 ihh=1,ihmax
        if(ih .eq. ihh) go to 4250
        z2=zdens(ihh)/ztaa(ihh)*zab(ihh,ih)*zxa22(ih,ihh)*zden
        zalpha(ih,ihh)=zxn02(ihh ,ih)*z2
        zbeta(ih,ihh)=zxn12(ihh,ih)*z2
        zdelta(ih,ihh)=-zxn22(ihh,ih)*z2
4250    continue
        zalpha(ih,ih)=zden*(zm20(ih)+z1*zxn02(ih,ih))
        zbeta(ih,ih)=zden*(zm21(ih)+z1*zxn12(ih,ih))
        zdelta(ih,ih)=0.
4200    continue
c
c       calculate the w arrays
c       see equation a21c in ppl 1473
c
        do 4500 ih=1,ihmax
        zmmnn=0.
        zmm22=0.
        do 4550 ihh=1,ihmax
        zmmnn=zmmnn+zab(ih,ihh)*zxmmnn(ih,ihh)
        zmm22=zmm22+zab(ih,ihh)*zxhm22(ih,ihh)
4550    continue
        zw(ih)=1./(1.+9.57/
     1  (zeps32*zeps32*zvstr(ih)*zvstr(ih)*zmm22*zmmnn))
4500    continue
c
c       calculate the alpha and beta hat arrays
c       see equation a21 in ppl 1473
c
        do 4700 ih=1,ihmax
        do 4750 ihh=1,ihmax
        zsuma=0.
        zsumb=0.
        do 4800 ihk=1,ihmax
        zsuma=zsuma+zw(ihk)*zdelta(ih,ihk)*zalpha(ihk,ihh)
        zsumb=zsumb+zw(ihk)*zdelta(ih,ihk)*zbeta(ihk,ihh)
4800    continue
        zhalph(ih,ihh)=zw(ih)*(zalpha(ih,ihh)+zsuma)
        zhbeta(ih,ihh)=zw(ih)*(zbeta(ih,ihh)+zsumb)
4750    continue
4700    continue
c
c       calculate the friction coefficients
c       see sph eq 20a-c
c
c       the friction coefficients have been normalized by the proton mass
c
        do 5000 ih=1,ihmax
        z1=zdens(ih)/ztaa(ih)
        do 5100 ihh=1,ihmax
        zsuma=0.
        zsumb=0.
        zsumc=0.
c
        do 5150 ihk=1,ihmax
        z2=z1*zab(ih,ihk)
        zsuma=zsuma+z2*zxn02(ih,ihk)*zhalph(ihk,ihh)
        zsumb=zsumb+z2*zxn02(ih,ihk)*zhbeta(ihk,ihh)
        zsumc=zsumc+z2*zxn12(ih,ihk)*zhbeta(ihk,ihh)
5150    continue
        z2=z1*zab(ih,ihh)
        z11=-zm02(ih)*zhalph(ih,ihh)+z2*zxn00(ih,ihh)-zsuma
        z12=zm02(ih)*zhbeta(ih,ihh)-z2*zxn01(ih,ihh)+zsumb
        z22=-zm12(ih)*zhbeta(ih,ihh)+z2*zxn11(ih,ihh)-zsumc
        if(ih .ne. ihh) go to 5160
        z11=z11+zm00(ih)
        z12=z12-zm01(ih)
        z22=z22+zm11(ih)
c
5160    zl11(ih,ihh)=zmass(ih)*z11
        zl12(ih,ihh)=zmass(ih)*z12
        zl22(ih,ihh)=zmass(ih)*z22
5100    continue
5000    continue
c
c       now the coefficients will be evaluated in the strong temperature
c       equilibration approximation
c       see sph section iii.b
c
c       invert the l22 matrix
c
c       create a dummy l22 array and an identity matrix
c
c
        do 6000 ih=1,ihmax
        do 6050 ihh=1,ihmax
        ziii(ih,ihh)=0.
6050    zzl22(ih,ihh)=zl22(ih,ihh)
        ziii(ih,ih)=1.
6000    continue
c
        call matrx1(zzl22,mxions,ihmax,ipivot,ierror)
c
        if(ierror .ne. 0) go to 9000
        call matrx2(zin22,ziii,zzl22,ipivot,mxions,ihmax,1,ihmax,1)
c
c       calculate the c2 and lhat22 matrix
c       see sph equation 35a and b
c
        zlhd=0.
        do 6200 ih=1,ihmax
        zsumb=0.
        do 6250 ihh=1,ihmax
        zsuma=0.
        do 6300 ihk =1,ihmax
        zsuma=zsuma+zl12(ihh,ihk)*zin22(ih,ihk)
6300    continue
        zc2(ih,ihh)=2.5*zdens(ih)/zdens(ihh)*zsuma
        zsumb=zsumb+zin22(ih,ihh)*zdens(ihh)
6250    continue
        zlh22(ih)=2.5*zsumb
        zlhd=zlh22(ih)*zdens(ih)+zlhd
6200    continue
c
c       evaluate zch2--see sph equation 40
c
        do 6400 ih=1,ihmax
        zch2(ih)=0.
        do 6350 ihh=1,ihmax
        zch2(ih)=zch2(ih)+zc2(ihh,ih)
6350    continue
6400    continue
c
c       calculate the effective l matrix (friction coefficients)
c       for strong temperature equilibration
c
        do 7000 ih=1,ihmax
        do 7500 ihh=1,ihmax
        zsum=0.
        do 7600 ihk=1,ihmax
        zsum=zsum+zl12(ih,ihk)*
     1  (zlh22(ihk)*zch2(ihh)/zlhd-zc2(ihk,ihh)/zdens(ihk))
7600    continue
        zs11(ih,ihh)=zl11(ih,ihh)+0.4*zdens(ihh)*zsum
        zs12(ih,ihh)=zdens(ih)*zdens(ihh)*zch2(ih)/zlhd
7500    continue
7000    continue
c
c       calculate the pfirsch-schluter coefficients (ps11 and ps12)
c       in the strong temperature approximation
c
c
c       if the weak temperature approximation is to be used
c       replace zs11 and zs12 respectively with zl11 and zl12 in
c       the following equations.
c
c
        do 8000 ih=1,ihmax
        za1=q(j)*q(j)*zlar2(ih)*zcharg(ih)/zmass(ih)
        do 8100 ihh=1,ihmax
        ps11(ih,ihh,j)=za1/zcharg(ihh)*zs11(ih,ihh)
        ps12(ih,ihh,j)=za1/zcharg(ihh)*zs12(ih,ihh)
8100    continue
8000    continue
c
        if(.not.ladjst) go to 8500
        do 8200 ihh=1,ihmax
        do 8200 ih=1,ihmax
        cl11(ih,ihh,j)=cfutz(iclflx)*cl11(ih,ihh,j)
        cl12(ih,ihh,j)=cfutz(iclflx)*cl12(ih,ihh,j)
        bp11(ih,ihh,j)=cfutz(ibpflx)*bp11(ih,ihh,j)
        bp12(ih,ihh,j)=cfutz(ibpflx)*bp12(ih,ihh,j)
        ps11(ih,ihh,j)=cfutz(ipsflx)*ps11(ih,ihh,j)
        ps12(ih,ihh,j)=cfutz(ipsflx)*ps12(ih,ihh,j)
 8200   continue
c
 8500   continue
c
10000   continue
c
c   les  nov-90  for d3he fusion, reorder ions back to d,t,p,3he,4he
c
      if ( cfutz(490).le.epslon ) then
        call reorder(ihmax,mzones,cl11)
        call reorder(ihmax,mzones,cl12)
        call reorder(ihmax,mzones,bp11)
        call reorder(ihmax,mzones,bp12)
        call reorder(ihmax,mzones,ps11)
        call reorder(ihmax,mzones,ps12)
      endif
c
c-----------------end of main do-loop over the radial index j-------------
c
c
        return
 9000   continue
        call error_olymp(1,iclass,isub,2,
     >         ' *** error *** in solution for l22 ')
        return
c**********************************************************************
c
        entry ncfprt
c
c       edit print-out of neoclassical particle transport coefficients
c
        if(cfutz(incflx).le.epslon) return
c
        lpage = lpage + 1
c
        z0=uist*1.e+03
        zt=tai*z0
        zdt=dtoldi*z0
        write (nprint,9010) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write(nprint,9020)
        do 8600 jz=lcentr,mzones
        write(nprint,9030) cl11(1,1,jz),cl12(1,1,jz),
     1  cl11(1,2,jz),cl12(1,2,jz),bp11(1,1,jz),bp12(1,1,jz),
     2  bp11(1,2,jz),bp12(1,2,jz),ps11(1,1,jz),ps12(1,1,jz),
     3  ps11(1,2,jz),ps12(1,2,jz)
 8600   continue
c
        return
 9010   format(1h1,2x,a48,10x,a72/
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.')
 9020 format(29x,' hawryluk-hirshman neoclassical particle transport'
     & //2x,'cl11(1,1)',1x,'cl12(1,1)',1x,'cl11(1,2)',1x,
     & 'cl12(1,2)',1x,'bp11(1,1)',1x,'bp12(1,1)',1x,'bp11(1,2)',1x,
     & 'bp12(1,2)',1x,'ps11(1,1)',1x,'ps12(1,1)',1x,'ps11(1,2)',1x,
     & 'ps12(1,2)')
 9030 format(1x,12(1x,1pe9.2))
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c

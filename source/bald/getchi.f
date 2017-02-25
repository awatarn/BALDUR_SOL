c--------1---------2---------3---------4---------5---------6---------7-c
c@getchi  .../baldur/bald/getchi.f
c  rgb 18-jul-01 implement call ftrapfn to compute trapbr(jz,js)
c    trapped particle option controlled by either lneocl(3) or
c    by cfutz(481) when cfutz(481) > 0
c  rgb 11-jul-01 function ft replaced with ftrapfn
c  rgb 13:50 12-nov-93 changed ftrap(3,1) to ftrap(1,3)
c  rgb 14:28 11-nov-93 changed call abort to call abortb
c   les  nov-90  trapping correction to resist ft .ge.0 for zeff>3
c         add fast ions from d3he fusion to n_e, zeff,totprs,betat, etc
c         hydrogen mass ameans includes protons (nimp(1))
c  dps 15-may-89 15.09 Add contributions from IRE variables here
c  dps 09-may-89 15.08 Fix old extf bug in rhoels calculation; move
c                cmean and c2mean extrapolation here from IMPRAD.
c  dps 10-nov-88 15.07 Fix bug in enlaes, enlais calculation.
c  dps 24-aug-88 15.00 Make changes required when using NC code:
c                first, skip rhoels, rhois, rhoins calculation in
c                getchi(1); second, skip rhoels and xzeff calculation
c                in getchi(2).
c  rgb 14.04 06-mar-88 Modify ftrap(j1,jz) near the magnetic axis
c       if cnumer(19) .gt. epslon .and. cnumer(20) .gt. epslon)
c       ftrap(j1,j) =
c         1. + (ftrap(j1,j)-1.) / (1. + ( cnumer(19) / x )**cnumer(20) )
c       where x = xbouni(j) for j1=1 or x=xzoni(j) for j1=2
c       this has the effect of flattening the central profile
c       and avoiding a spike of current density at the magnetic axis.
c  rgb 14.02 19-feb-88 Modify ftrap(j1,jz) near the magnetic axis
c       only if  lnumer(20) .gt. 3   .and.   lnumer(19) .lt. lnumer(20)-2
c       ftrap(...) = 1.0  for x .le. xzoni(i1)   i1 = max(lnumer(19),1)
c       ftrap(...) = f1 + a2 * z**2 + a4 * z**4    0. .le. z .le. 1.
c                 for xzoni(i1) .le. x .le. xzoni(i3)   i3 = lnumer(20)
c                 for z = (x-xzoni(i1)) / (xzoni(i3)-xzoni(i1))
c       ftrap(...) matches previous values
c                 at x=xzoni(i3) and x=xzoni(i3-1)
c  rgb 12.62 26-aug-87 enforce symmetry across magnetic axis
c  dps 07-may-87 introduce switch for trapped fraction calculation
c                used in nu-star.
c  dps 04-mar-87 add bootstrap current cf. Mike Hughes.
c  dps 07-jan-87 renamed presus -> thrprs and added totprs for total
c                plasma pressure computed using cequil(2 & 3)
c  rgb 02-sep-86 fixed interpolation of ajtori from centers to bndries
c  rgb 16-mar-86 used ramp fnct to smooth ftrap near mag axis
c  rgb 13-feb-86  bpols(2,1) = - bpols(2,3)
c               expressions for xnuhyd and xnuel improved
c  rgb 13-oct-85 replaced all references to magnetic field with q-value
c     in computation of xnuhyd and xnuel.
c     zd = ahalfs(jz,j1)/rmids(jz,j1)  inverse aspect ratio
c     tested for ftrap(j1,jz) .lt. epslon
c  rgb 25-may-85 changed computation of q(j), ajzs(j,it), bpols(1,mzones+1)
c       jmh 15-jul-82 if cfutz(410).gt.0, the transport coefficient is a
c                     function of the te profile instead of the te value
cdoc
c**********************************************************************c
c
c       ------------
c       sbrtn GETCHI   file DSOLVER
c       ------------
c
c       2.12    set tes, tis, etc.
c
c
c  2)  compute:  rmins, rmajs, bzs, rhohs, rhois, rhoins, rhoels
c                tes, tis, bpols
c
c      given:    rmini, rmaji, bzi, chi, cmean, rhobes, bpoli
c
c
c  3.1)    boundary-center quantities
c
c  linearly interpolate from zone centers to zone boundaries
c  rhohs = hydrogen species density
c  rhois = impurity species density
c  cmean = mean charge state z of impurity
c  c2mean = mean z**2 of impurity
c  tes    = electron temperature in ergs
c  tis    = ion temperature in ergs
c  rhoins = ion density
c  rhobes = beam electron density
c  rhobis = beam ion density
c  bpols  = d poloidal manetic flux / d rho  / 2 pi r0ref
c
c
c  3.2)    compute q, jz, zeff
c
c bateman computing q(j) and ajzs in sbrtn getchi
c
c..compute q(j) = b0ref * rho(j) / ( r0ref * bpoli(j) )  on zone boundaries
c
c
cend
c**********************************************************************c
c
        subroutine getchi(k)
c
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
       include 'comncr.m'
       include 'cd3he.m'
c
c-----------------------------------------------------------------------
c
      dimension  zft(mj)
c
        data    iclass /2/,     isub /12/
c
c       if "cfutz(ispitz)" is .gt. epslon, the spitzer resistivity
c       is multiplied by "cfutz(ispitz)".
c
      data    ispitz /401/
c
c
        if (.not.nlomt2(isub)) go to 10
        call mesage('% *** 2.12 subroutine getchi bypassed')
        return
   10   continue
c
c-----------------------------------------------------------------------
c
c
c      1)      check k
c
c
        if (k.le.0.or.k.gt.2) return
        go to (200,300), k
c
c
c  2)  compute:  rmins, rmajs, bzs, rhohs, rhois, rhoins, rhoels
c                tes, tis, bpols
c
c      given:    rmini, rmaji, bzi, chi, cmean, rhobes, bpoli
c
c
  200   continue
c
        rmins = rmini * uisl
        rmajs = rmaji * uisl
        bzs = bzi * uisb
        znmin = 2.0 * epslon * uisd
c
c
        do 240 j1 = 1, mzones
        znei = 0.0
        znii = 0.0
c
        if (mhyd.le.0) go to 212
        do 210 j2 = lhyd1, lhydn
        rhohs(j2,2,j1) = chi(j2,j1) * uisd
        znei = znei + chi(j2,j1) * extf
        znii = znii + chi(j2,j1)
  210   continue
  212   continue
c
c  15.09 Set chi, cmean, and c2mean for Impurity Rate Equations variables
c  15.11 Do NOT include neutrals here => chi = impurity ion density;
c        code here is UNCHANGED.
c
        if ((natomc.eq.3).and.(mimp.ne.0)) then
          do 215 ji=1,mimp
            jimp = limp1 - 1 + ji
            nk = nkimp(ji)
            znk = 0.0
            znkk = 0.0
            znkk2 = 0.0
c
            do 214 ik=1,nk
              zn = max(znmin,xn(j1,ik,ji))
              znk = znk + zn
              znkk = znkk + float(ik) * zn
              znkk2 = znkk2 + float(ik**2) * zn
  214       continue
c
            chi(jimp,j1) = max(znk*usid,2.0*epslon)
            cmean(ji,2,j1) = znkk / znk
            c2mean(ji,2,j1) = znkk2 / znk
c
c..Diagnostic option
c
            if ((limprd(1).ne.0).and.(tbi.gt.tinit)) then
              chi(jimp,j1) = rhois(ji,2,j1) * usid
              cmean(ji,2,j1) = 1.
              c2mean(ji,2,j1) = 1.
            end if
  215     continue
        end if
c
        if (mimp.le.0) go to 222
        do 220 j2 = limp1, limpn
        i001 = j2+1-limp1
        rhois(i001,2,j1) = chi(j2,j1) * uisd
        znei = znei + chi(j2,j1) * cmean(i001,2,j1)
        znii = znii + chi(j2,j1)
  220   continue
  222   continue
c
        znei = znei + rhobes(2,j1)*usid
        rhoins(2,j1) = znii * uisd
        rhoels(2,j1) = znei * uisd
c
  223   continue
        tes(2,j1) = epsinv
        tis(2,j1) = epsinv
        if (znei.gt.epslon) tes(2,j1) = chi(lelec,j1)
     1                                   / znei * uish
        if (znii.gt.epslon) tis(2,j1) = chi(lion,j1)
     1                                   / znii * uish
        bpols(1,j1) = bpoli(j1)*uisb
  240   continue
c
      bpols(1,mzones+1) = bpoli(mzones+1) * uisb
c
c
        return
c
c
c      3)      call imprad, set boundary points, zeff, q, jz
c
c
  300   continue
c
        call imprad(2)
c
c..15.08 Moved the following from IMPRAD
c..Reset dweirs, cmean, and c2mean at zone center outside plasma
c
c..reset cmean=1. and c2mean=1. as in original 1-D BALDUR
c
      if ( limprd(32) .eq. 0 ) then
c
        call resetr(cmean(1,2,mzones),mimp,1.0)
        call resetr(c2mean(1,2,mzones),mimp,1.0)
c
c..linearly extrapolate cmean and c2mean to zone center outside plasma
c
      elseif ( limprd(32) .eq. 1 ) then
c
        zint1 = (xzoni(mzones  ) - xzoni(mzones-2))
     &        / (xzoni(mzones-1) - xzoni(mzones-2))
        zint2 = (xzoni(mzones-1) - xzoni(mzones  ))
     &        / (xzoni(mzones-1) - xzoni(mzones-2))
c
        do 301 ji=1,mimp
          cmean(ji,2,mzones) = max ( 1. ,
     &  zint1 * cmean(ji,2,mzones-1) + zint2 * cmean(ji,2,mzones-2) )
          c2mean(ji,2,mzones) = max ( 1. ,
     &  zint1 * c2mean(ji,2,mzones-1) + zint2 * c2mean(ji,2,mzones-2) )
 301  continue
c
      endif
c
c      3.1)    boundary-center quantities
c
c  linearly interpolate from zone centers to zone boundaries
c  rhohs = hydrogen species density
c  rhois = impurity species density
c  cmean = mean charge state z of impurity
c  c2mean = mean z**2 of impurity
c  tes    = electron temperature in ergs
c  tis    = ion temperature in ergs
c  rhoins = ion density
c  rhobes = beam electron density
c  rhobis = beam ion density
c  bpols  = d poloidal manetic flux / d rho  / 2 pi r0ref
c
c
        do 308 jb = lcentr, mzones
        zint0 = (xzoni(jb) - xbouni(jb)) / dxboui(jb)
        zint1 = 1.0 - zint0
c
c              hydrogen variables
c
        do 302 jh = 1, mhyd
        rhohs(jh,1,jb) = rhohs(jh,2,jb-1)*zint0 + rhohs(jh,2,jb)*zint1
  302   continue
c
c              impurity variables
c
        if (mimp.le.0) go to 305
c
        do 304 ji = 1, mimp
        rhois(ji,1,jb) = rhois(ji,2,jb-1)*zint0 + rhois(ji,2,jb)*zint1
        cmean(ji,1,jb) = cmean(ji,2,jb-1)*zint0 + cmean(ji,2,jb)*zint1
        c2mean(ji,1,jb) = c2mean(ji,2,jb-1)*zint0 +
     1                                  c2mean(ji,2,jb)*zint1
  304   continue
c
  305   continue
c
c
c
        tes(1,jb) = tes(2,jb-1)*zint0 + tes(2,jb)*zint1
        tis(1,jb) = tis(2,jb-1)*zint0 + tis(2,jb)*zint1
        rhoins(1,jb) = rhoins(2,jb-1)*zint0 + rhoins(2,jb)*zint1
        rhobes(1,jb) = rhobes(2,jb-1)*zint0 + rhobes(2,jb)*zint1
        rhobis(1,jb) = rhobis(2,jb-1)*zint0 + rhobis(2,jb)*zint1
c
c   les  nov-90  d3he fast particle
c
         rh1fst(1,jb)=rh1fst(2,jb-1)*zint0 + rh1fst(2,jb)*zint1
         rh2fst(1,jb)=rh2fst(2,jb-1)*zint0 + rh2fst(2,jb)*zint1
c
c              b-poloidal
c
        bpols(2,jb) = 0.0
        if (xzoni(jb).le.epslon) go to 308
        zint0 = (xbouni(jb+1)**2 - xzoni(jb)**2) /
     1          (xbouni(jb+1)**2 - xbouni(jb)**2)
        zint1 = (1.0 - zint0) * xbouni(jb+1) / xzoni(jb)
        zint0 = zint0 * xbouni(jb) / xzoni(jb)
 
        bpols(2,jb) = bpols(1,jb)*zint0 + bpols(1,jb+1)*zint1
  308   continue
c
c
c
        do 312 jh = 1, mhyd
        rhohs(jh,1,1) = rhohs(jh,1,2*lcentr-1)
  312   continue
c
        if (mimp.le.0) go to 315
        do 314 ji = 1, mimp
        rhois(ji,1,1) = rhois(ji,1,2*lcentr-1)
        cmean(ji,1,1) = cmean(ji,1,2*lcentr-1)
        c2mean(ji,1,1) = c2mean(ji,1,2*lcentr-1)
  314   continue
  315   continue
c
        tes(1,1) = tes(1,2*lcentr-1)
        tis(1,1) = tis(1,2*lcentr-1)
        rhoins(1,1) = rhoins(1,2*lcentr-1)
        rhobes(1,1) = rhobes(1,2*lcentr-1)
        rhobis(1,1) = rhobis(1,2*lcentr-1)
c
c   les  nov-90  fast particles
c
      rh1fst(1,1) = rh1fst(1,2*lcentr-1)
      rh2fst(1,1) = rh2fst(1,2*lcentr-1)
c
        bpols(2,1) = - bpols(2,2*lcentr-1)
c
c
c      3.2)    compute q, jz, zeff
c
c bateman computing q(j) and ajzs in sbrtn getchi
c
c..compute q(j) = r0ref * rho(j) / ( r0ref * bpoli(j) )  on zone boundaries
c
 
      do 320 j=3,mzones+1
      q(j) = b0ref * avi(j,1,1) / ( r0ref * bpoli(j) )
 320  continue
c
      q(2) = ( xbouni(4)**2 * q(3) - xbouni(3)**2 * q(4) )
     &       / ( xbouni(4)**2 - xbouni(3)**2 )
      q(1) = q(3)
c
c..current density on zone centers ajzs(2,j) in standard units
c  defined as the toroidal current through the cross-sectional area
c    between differentially close flux surfaces divided by the
c    cross-sectional area between differentially close flux surfaces
c    = < Jphi / R > / < 1. / R >
c
      ajtori(1,2) = ajtori(2,2)
      do 322 j=1,mzones
      ajzs(2,j) = ajtori(j,2)*uisj / (0.5*(avi(j+1,11,1)+avi(j,11,1)))
 322  continue
c
c..interpolate current density to zone boundaries
c
      call cubint (xzoni,ajtori(1,2),mzones-1,0,ceqtmp, Kjbal
     & ,xbouni,eqtmp,mzones+1,0,0.,1
     & ,'cubic interp of ajtori to zone bndries in sbrtn getchi')
c
      do 324 j=1,mzones+1
      ajzs(1,j) = eqtmp(j) * uisj / avi(j,11,1)
 324  continue
c
cend bateman
c
c  ...To match BALDP86m, need to comput ajzs directly rather than by
c  ...interpolation. For this, use lines below to "326 continue".
c
c      zcc=fcc*10.0**fxc
c      do 326 jz=1,mzones
c        if (jz.le.lcentr) go to 325
c        ajzs(1,jz)=zcc*
c     1  (xzoni(jz)*bpols(2,jz)-xzoni(jz-1)*bpols(2,jz-1))/
c     2  ((xzoni(jz)-xzoni(jz-1))*(xzoni(jz)+xzoni(jz-1))
c     3  *rmins*2.0*fcpi)
c        if (jz.ge.2*lcentr) go to 326
c        i001=2*lcentr-jz
c        ajzs(1,i001)=ajzs(1,jz)
c        go to 326
c  325   continue
c        if (jz.eq.lcentr) ajzs(1,jz)=ajzs(2,jz)
c  326 continue
c  ...............
c
c      rhoels,  zeff, gamma(zeff) (spitzer resistivity term)
c
c  rhoels = electron density computed from charge neutrality
c  xzeff  = z effective at boundary/zone jz
c  gamma  = parameter in definition of canonical angular momentum
c
c  note:  extf  and  extzef  are namelist input variables
c
c
        do 358 jz = 1, mzones
c
        do 357 j1 = 1, 2
c
        zni = 0.0
        do 352 jh = 1, mhyd
        zni = zni + rhohs(jh,j1,jz)
  352   continue
c
        zne = zni*extf
        zne2 = zne
        if (mimp.le.0) go to 355
        do 354 j2 = 1, mimp
        zni = zni + rhois(j2,j1,jz)
        zne = zne + rhois(j2,j1,jz)*cmean(j2,j1,jz)
        zne2 = zne2 + rhois(j2,j1,jz)*c2mean(j2,j1,jz)
  354   continue
  355   continue
c
c   les  nov-90  d3he fast particles
c
         zne = zne + rh1fst(j1,jz) + 2.*rh2fst(j1,jz)
         zne2 = zne2 + rh1fst(j1,jz) + 4.*rh2fst(j1,jz)
c
        rhoels(j1,jz) = zne + rhobes(j1,jz)
        xzeff(j1,jz) = zne2 / zne + (extzef - 1.0)
  356   continue
c
c       21-jan-77 gspitz is due to hirshman, jan. 1977
c
      gspitz(j1,jz) = (xzeff(j1,jz) * (2.67 + xzeff(j1,jz))) /
     &                                  (3.40 * (1.13 + xzeff(j1,jz)))
c
c
  357   continue
c
  358   continue
c
c
c
c      2.)     compute useful quantities
c
c..compute trapped particle fraction
c
c  itrap = lneocl(3) or 
c  itrap = the integer part of cfutz(481) if that is non-zero
c
        itrap = lneocl(3)
        if ( cfutz(481) - 0.1 .gt. 0.0 )
     &       itrap = int ( cfutz(481) + 0.1 )
c
        call ftrapfn( itrap )
c
        if(cfutz(410).le.epslon) cfutz(410)=1.0  ! see cnueqs(jz) below
c
      do 398 jz = 1, mzones
c
c       thermal pressure at zone center
c
      thrprs(jz)=rhoels(2,jz)*tes(2,jz)+rhoins(2,jz)*tis(2,jz)
c
c...totprs(jz) is the sum of
c...electron + ion thermal pressure
c...cequil(2) times the fast alpha particle pressure and
c...cequil(3) times the fast beam ion pressure which, for now,
c...is treated as a scalar pressure.
c   les  nov-90  d3he fusion
c
      totprs(jz)=thrprs(jz)+2.*(cequil(2)*alphai(jz)*ealfai(jz)
     1           *uisd*uise+cequil(3)*hebems(jz)*rhobis(2,jz)
     2           +d3fast(jz))/3.
c
c
c               alpha
c
c
      z0 = 0.0
c
      do 360 j1 = 1, mhyd
      z0 = z0 + rhohs(j1,1,jz)
  360 continue
c
      calph(jz) = rhoels(1,jz) * xzeff(1,jz) / z0   -   1.0
c
c
      do 388 j1 = 1, 2
c
c
c               ion coulomb log
c
c  note:  evsinv and ev50s are set in file default
c
c
      clogis(j1,jz) = log(evsinv * tis(j1,jz) * 1.54e10 /
     1  xzeff(j1,jz) * sqrt((evsinv * tes(j1,jz)) / rhoels(j1,jz)))
      clogis(j1,jz)=max(1.0,clogis(j1,jz))
c
c
c               electron coulomb log
c
c
      if (tes(j1,jz).gt.ev50s) go to 362
c
      cloges(j1,jz) = log(1.54e10 * sqrt((evsinv * tes(j1,jz))**3
     1                                  / rhoels(j1,jz)) / xzeff(j1,jz))
      go to 364
c
  362 continue
c
      cloges(j1,jz) = log(evsinv * tes(j1,jz) * 1.09e11 /
     1                  (sqrt(rhoels(j1,jz)) * xzeff(j1,jz)))
c
  364 continue
      cloges(j1,jz)=max(1.0,cloges(j1,jz))
c
c
c               average hydrogen mass
c
c
      z0 = 0.0
      z1 = 0.0
c
c  aspec = atomic wieght of ion species, set in file default
c
      do 366 j2 = 1, mhyd
        z0 = rhohs(j2,j1,jz) + z0
        z1 = aspec(j2) * rhohs(j2,j1,jz) + z1
  366 continue
        z2=z0
c
c   les  nov-90  for d3he include protons  ???
c
         i2st = 1
      if ( nfusn.eq.4 ) then
         i2st = lprotn +1
        z0 = rhois(lprotn-lhydn,j1,jz) + z0
        z1 = aspec(lprotn) * rhois(lprotn-lhydn,j1,jz) + z1
      endif
c
      ahmean(j1,jz) = z1 / z0
c
c
c               average ion mass
c
c
      if(mimp.gt.0) go to 370
      aimass(j1,jz)=ahmean(j1,jz)
      go to 374
  370 continue
      do 372 i2=1,mimp
      ii=i2+lhydn
      z0=rhois(i2,j1,jz)+z0
      z1=aspec(ii)*rhois(i2,j1,jz)+z1
  372 continue
c
      aimass(j1,jz)=z1/z0
c
  374 continue
c
c
c               ion collision time -- braginski
c       (modified to include the effects of impurities)
c
c  note:  coefficients  ctions  and  cteles  set in file default
c
c
      if (mimp.le.0) go to 378
c
      do 376 ji=1,mimp
      z2=z2+rhois(ji,j1,jz)*c2mean(ji,j1,jz)
  376 continue
c
  378 continue
c
c
      tions(j1,jz) = ctions * tis(j1,jz) * sqrt(tis(j1,jz) *
     1          ahmean(j1,jz)) / (z2 * clogis(j1,jz))
c
c
c               electron collision time -- braginski
c
c
      telecs(j1,jz) = cteles * tes(j1,jz) * sqrt(tes(j1,jz)) /
     1          (rhoels(j1,jz) * cloges(j1,jz))
c
c
c  xnuhyd and xnuel are undefined (infinite) for r = 0 or
c       b-poloidal = 0.
c  ftrap depends on xnuel, but is defined for r or b-poloidal = 0
c  coefficients  cnuhyd  and  cnuel  are set in file default
c
c
      if (abs(ahalfs(jz,j1)).le.epslon) go to 380
c
      if (abs(bpols(j1,jz)).le.epslon) go to 380
c
      zd = abs ( ahalfs(jz,j1) / rmids(jz,j1) )
c
c
c               banana parameter -- hydrogen
c
c
      xnuhyd(j1,jz) = abs ( cnuhyd * rbtors(jz,j1)
     &  * sqrt( abs( rmids(jz,j1)*ahmean(j1,jz)
     &                          / (ahalfs(jz,j1)*tis(j1,jz))) )
     &  / ( bpols(j1,jz) * tions(j1,jz) ) )
c
c
c               banana parameter -- electron
c
c
      xnuel(j1,jz) = abs ( cnuel * rbtors(jz,j1)
     &  * sqrt( abs (rmids(jz,j1) / (ahalfs(jz,j1)*tes(j1,jz))) )
     &  / ( bpols(j1,jz) * telecs(j1,jz) ) )
c
c
c               resistivity enhancement due to trapped electrons
c
c               from:   hirshman, hawryluk, and birge,
c                       "neoclassical conductivity of a tokamak plasma"
c                       princeton plasma physics lab
c                       january, 1977
c
c
c..smooth cutoff for zft(j), bateman, 16-feb-86
c  use ramp function to smooth ftrap near magnetic axis
c  zdramp = 1.0 for r/a > cemprc(20)
c  zdramp = r / (a * cemprc(20)) for r/a < cemprc(20)
c
      zdramp = min ( 1.0 ,
     &  ahalfs(jz,j1)/(ahalfs(mzones,1)*max(cemprc(20),epslon) ) )
c
      zft(jz) = 1.0 -  (1.0 - zd)**2
     &  / (sqrt(abs(1.0 - zd**2))*(1.0 + 1.46*sqrt( abs(zd*zdramp))))
c
      if ( itrap .gt. 0 ) zft(jz) = trapbr(jz,j1)
c
      zxi = 0.58 + 0.20 * xzeff(j1,jz)
c
c   les  nov-90  require ftrap.ge.0
c
      zc = ( 0.56*max(3.0 - xzeff(j1,jz),0.)/
     1                  ( xzeff(j1,jz)*(3.0 + xzeff(j1,jz)) ) )
      z0 = zft(jz) / (1.0 + zxi * xnuel(j1,jz))
c
      ftrap(j1,jz) = (1.0 - z0) * (1.0 - zc*z0)
c
      if (ftrap(j1,jz) .lt. epslon) call abortb (6
     &  ,'ftrap(j1,jz) .lt. epslon in sbrtn getchi')
c
      go to 382
c
  380 continue
c
c
c               values of xnuhyd, xnuel, and ftrap if b-poloidal or r = 0
c
c
      xnuhyd(j1,jz) = epsinv
      xnuel(j1,jz) = epsinv
      ftrap(j1,jz) = 1.0
c
  382 continue
c
c
  388 continue
c
c
c               electron-ion temperature equilibration rate
c
c
      z0 = 0.0
c
c   les nov-90  add protons ???  (not done)
c
      do 390 jh = 1, mhyd
      z0 = z0 + rhohs(jh,2,jz)/aspec(jh)
  390 continue
c
      if (mimp.le.0) go to 394
c
      do 392 ji = 1, mimp
      i001 = ji + limp1 - 1
      z0 = z0 + rhois(ji,2,jz) * c2mean(ji,2,jz) / aspec(i001)
  392 continue
c
  394 continue
c
      cnueqs(jz) = cfutz(410) * ccnu * z0 / telecs(2,jz)
c
c
  398 continue
c
c  end of loop over jz
c
c
c..Modify ftrap(j1,jz) within zone center lnumer(20)
c       only if  lnumer(20) .gt. 3   .and.   lnumer(19) .lt. lnumer(20)-2
c  ftrap as a fn of r should have
c        (1)  ftrap(0.) = f1 = 1.0
c        (2)  d ftrap(0.) / d r = 0.0
c        (3)  ftrap matching the previously computed values
c             at zone centers lnumer(20) and lnumer(20)-1
c             ie  to match in both value and finite difference slope
c       let i1 = lnumer(19)
c           i3 = lnumer(20)   then
c       ftrap(...) = 1.0  for x .le. xzoni(i1)
c       ftrap(...) = f1 + a2 * z**2 + a4 * z**4    0. .le. z .le. 1.
c                 for xzoni(i1) .le. x .le. xzoni(i3)
c                 where   z = (x-xzoni(i1)) / (xzoni(i3)-xzoni(i1))
c       ftrap(...) matches previous values
c                 at x=xzoni(i3) and x=xzoni(i3-1)
c
      if (lnumer(20) .gt. 3) then
        i1  = min ( lnumer(20) - 2 , max ( lnumer(19) , 2 ) )
        i3  = lnumer(20)
        zf1 = 1.
        zf2 = ftrap(2,i3-1)
        zf3 = ftrap(2,i3)
        zinv = 1. / ( xzoni(i3) - xzoni(i1) )
        zrt = ( xzoni(i3-1) - xzoni(i1) ) * zinv
c
        za4 = ( zf3 - zf1 - (zf2-zf1)/zrt**2 ) / ( 1. - zrt**2 )
        za2 = zf3 - zf1 - za4
c
        if ( i1 .gt. 1) then
          do 401 j=1,i1
            ftrap(1,j) = zf1
            ftrap(2,j) = zf1
 401      continue
        endif
c
        do 402 j=i1+1,i3
c
          ftrap(1,j) = zf1 + ( (xbouni(j)-xzoni(i1))*zinv )**2
     &          * ( za2  + za4 * ( (xbouni(j)-xzoni(i1))*zinv )**2 )
          ftrap(2,j) = zf1 + ( (xzoni(j) -xzoni(i1))*zinv )**2
     &          * ( za2  + za4 * ( (xzoni(j) -xzoni(i1))*zinv )**2 )
c
 402    continue
c
      endif
c
c..new modification for ftrap(j1,j):  Bateman 06-mar-88
c
      if ( cnumer(19) .gt. epslon .and. cnumer(20) .gt. epslon ) then
c
        do 403 j=3,mzones
          ftrap(1,j) = 1. + ( ftrap(1,j) - 1. ) /
     &      ( 1. + ( cnumer(19) / xbouni(j) )**cnumer(20) )
          ftrap(2,j) = 1. + ( ftrap(2,j) - 1. ) /
     &      ( 1. + ( cnumer(19) /  xzoni(j) )**cnumer(20) )
 403    continue
c
          ftrap(2,2) = 1. + ( ftrap(2,2) - 1. ) /
     &      ( 1. + ( cnumer(19) /  xzoni(2) )**cnumer(20) )
          ftrap(1,2) = 1.
          ftrap(2,1) = ftrap(2,2)
          ftrap(1,1) = ftrap(1,3)
c
      endif
c
c
c..eta(j1,jz) = electrical resistivity at zone boundary / center jz
c
      if ( cfutz(ispitz) .gt. epslon ) then
c
        do 404 j=1,mzones
c
          eta(1,j) = cfutz(ispitz) * ceta * gspitz(1,j) /
     &               ( rhoels(1,j) * telecs(1,j) )
          eta(2,j) = cfutz(ispitz) * ceta * gspitz(2,j) /
     &               ( rhoels(2,j) * telecs(2,j) )
c
 404    continue
c
      else
c
        do 406 j=1,mzones
c
          eta(1,j) = ceta * gspitz(1,j) /
     &               (rhoels(1,j) * ftrap(1,j) * telecs(1,j))
          eta(2,j) = ceta * gspitz(2,j) /
     &               (rhoels(2,j) * ftrap(2,j) * telecs(2,j))
c
 406    continue
c
      endif
c
c..Enforce symmetry across magnetic axis
c
      if ( versno .gt. 13.99 ) then
c
        thrprs(1)   = thrprs(2)
        totprs(1)   = totprs(2)
c
        calph(1)    = calph(3)
c
        clogis(1,1) = clogis(1,3)
        clogis(2,1) = clogis(2,2)
        cloges(1,1) = cloges(1,3)
        cloges(2,1) = cloges(2,2)
c
        ahmean(1,1) = ahmean(1,3)
        ahmean(2,1) = ahmean(2,2)
        aimass(1,1) = aimass(1,3)
        aimass(2,1) = aimass(2,2)
c
        tions(1,1)  = tions(1,3)
        tions(2,1)  = tions(2,2)
        telecs(1,1) = telecs(1,3)
        telecs(2,1) = telecs(2,2)
c
        xnuhyd(1,1) = xnuhyd(1,3)
        xnuhyd(2,1) = xnuhyd(2,2)
        xnuel (1,1) = xnuel (1,3)
        xnuel (2,1) = xnuel (2,2)
c
        ftrap (1,1) = ftrap (1,3)
        ftrap (2,1) = ftrap (2,2)
c
        eta (1,1)   = eta(1,3)
        eta (2,1)   = eta(2,2)
c
        cnueqs(1)   = cnueqs(2)
c
      endif
c
c
c..integrated variables    Bateman, PPPL, 4-feb-86
c   les  nov-90  add d3he
c
c  set integrated variables to zero at the magnetic axis
c
      curnts(2) = 0.   ! toroidal plasma current
c
      erges(2)  = 0.   ! thermal electron energy
      ergis(2)  = 0.   ! thermal ion energy
      ergbs(2)  = 0.   ! beam particle energy
      ergas(2)  = 0.   ! alpha particle energy
      ergds(2)  = 0.   ! d 3he fusion particle energy
      ergts(2)  = 0.   ! sum of the above
c
      envaes(2) = 0.   ! vol ave electron density
      envais(2) = 0.   ! vol ave ion density
c
      enlaes(2) = 0.   ! line ave electron density
      enlais(2) = 0.   ! line ave ion density
c
      gealfs(2) = 0.   ! alpha particle power to electrons
         ged3fs(2) = 0.    ! d 3he fusion to electrons
         geashs(2) = 0.    ! loss from fast ion ash removal
      geauxs(2) = 0.   ! auxilliary power to electrons
      gebems(2) = 0.   ! beam power to electrons
      gebrs(2)  = 0.   ! bremsstrahlung cooling rate of electrons
      geecrs(2) = 0.   ! ECRH power to electrons
      geicrs(2) = 0.   ! ICRH power to electrons
      geions(2) = 0.   ! hyd. ionization cooling rate of electrons
      geohms(2) = 0.   ! OH heating rate of electrons
      gesrs(2)  = 0.   ! synchrotron radiation cooling rate of electrons
         gesyns(2) = 0.   ! high temp synchrotron rad cooling
      geirs(2) = 0.   ! impurity radiation
c
      gialfs(2) = 0.   ! alpha particle power to ions
         gid3fs(2) = 0.    ! d 3he fusion to ions
      giauxs(2) = 0.   ! auxilliary power to ions
      gibems(2) = 0.   ! beam power to ions
      gichxs(2) = 0.   ! charge exchange cooling rate of ions
      giecrs(2) = 0.   ! ECRH power to ions
      giicrs(2) = 0.   ! ICRH power to ions
      giions(2) = 0.   ! hyd. ionization cooling rate of ions
c
c
c
c..volume integrals out to zone boundaries
c
      do 410 j=2,mzones
c
      zdvols = vols(j+1,1) - vols(j,1)   ! volume in zone j
c
      erges(j+1) = erges(j) + chi(lelec,j)*zdvols*uisd*uise
      ergis(j+1) = ergis(j) + chi(lion ,j)*zdvols*uisd*uise
      ergbs(j+1) = ergbs(j) + hebems(j)*rhobis(2,j)*zdvols
      ergas(j+1) = ergas(j) + alphai(j)*ealfai(j)*zdvols*uisd*uise
      ergds(j+1) = ergds(j) + d3fast(j)*zdvols
      ergts(j+1) = erges(j+1)+ergis(j+1)+ergbs(j+1)+ergas(j+1)
     1  + ergds(j+1)
c
      envaes(j+1) = envaes(j) + rhoels(2,j)*zdvols
      envais(j+1) = envais(j) + rhoins(2,j)*zdvols
c
c  15.07 Bug fix: distance element should be distance between zone boundaries
c  rather than zone centers. At present, this definition is inconsistent
c  with usage elsewhere (e.g., NEUGAS, MPRINT).
c
      if (versno.gt.15.06) then
        zdrs = ahalfs(j+1,1) - ahalfs(j,1)
      else
        zdrs = ahalfs(j+1,2) - ahalfs(j,2)
      end if
c
      enlaes(j+1) = enlaes(j) + rhoels(2,j) * zdrs
      enlais(j+1) = enlais(j) + rhoins(2,j) * zdrs
c
      gealfs(j+1) = gealfs(j) + wealfs(j)*zdvols
       ged3fs(j+1) = ged3fs(j) + wed3fs(j)*zdvols
       geashs(j+1) = geashs(j) + weash(j)*zdvols
      geauxs(j+1) = geauxs(j) + weauxs(j)*zdvols
      gebems(j+1) = gebems(j) + webems(j)*zdvols
      gebrs (j+1) = gebrs (j) + webrs (j)*zdvols
      geecrs(j+1) = geecrs(j) + weecrh(j)*zdvols
      geicrs(j+1) = geicrs(j) + weicrf(j)*zdvols
      geions(j+1) = geions(j) + weions(j)*zdvols
      geohms(j+1) = geohms(j) + weohms(j)*zdvols
      gesrs (j+1) = gesrs (j) + wesrs (j)*zdvols
      gesyns(j+1) = gesyns(j) + wesyn(j)*zdvols
c
      if ( mimp.gt.0 ) then
        zweirs = 0.
        do 409 ji=1,mimp
          zweirs = zweirs + weirs(ji,j)
 409    continue
        geirs(j+1) = geirs(j)  + zweirs*zdvols
      endif
c
      gialfs(j+1) = gialfs(j) + wialfs(j)*zdvols
       gid3fs(j+1) = gid3fs(j) + wid3fs(j)*zdvols
      giauxs(j+1) = giauxs(j) + wiauxs(j)*zdvols
      gibems(j+1) = gibems(j) + wibems(j)*zdvols
      gichxs(j+1) = gichxs(j) + wichxs(j)*zdvols
c
c   les jan-91 change gichxs to follow sprint
c     gichxs(j+1) = gichxs(j) + (wichxs(j)-dwicxs(j)*tis(2,j))*
c    1   zdvols
c
      giecrs(j+1) = giecrs(j) + wiecrh(j)*zdvols
      giicrs(j+1) = giicrs(j) + wiicrf(j)*zdvols
      giions(j+1) = giions(j) + wiions(j)*zdvols
c
 410  continue
c
c..further normalizations
c
      zcur = uisi * r0ref / ( twopi * emu0 )
c
      do 420 j=3,mzones+1
c
      curnts(j) = zcur * bpoli(j) * avi(j,2,1) * avi(j,3,1) * avi(j,7,1)
c
      envaes(j) = envaes(j)  / vols(j,1)
      envais(j) = envais(j) / vols(j,1)
c
      enlaes(j) = enlaes(j) / ahalfs(j,1)
      enlais(j) = enlais(j) / ahalfs(j,1)
c
      zbetas = 16. * fcpi / ( 3. * vols(j,1) * bzs**2 )
c
      betate(j) = erges(j) * zbetas  ! thermal electron toroidal betat
      betati(j) = ergis(j) * zbetas  ! thermal ion      toroidal betat
      betatb(j) = ergbs(j) * zbetas  ! beam ion         toroidal betat
      betata(j) = ergas(j) * zbetas  ! alpha particle   toroidal betat
       betatd(j) = ergds(j) * zbetas   ! d3he fast particle tor betat
      betatt(j) = ergts(j) * zbetas  ! sum of the above
c
 420  continue
c
c
c
c
c
        return
c
        end

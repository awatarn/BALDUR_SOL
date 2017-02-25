c--------1---------2---------3---------4---------5---------6---------7-c
c@solveb  .../baldur/code/bald/solveb.f
c  rgb 07-aug-01 recompute ajboot and ajtpbi only if lneocl(1) < 1
c  rgb 07-feb-00 removed equivalence statement
c     rgb 12-feb-90 18.06 use sfutz(16) to control Chang current ajtpbi
c     ces 28-jul-89 16 fix current drive
c     ces 28-jul-89 16.16 add current drive
c     ces 28-jul-89 16.16 turn "rf" current on and off with heating
c     ces 26-jul-89 16.16 add isep for "rf" current drive profile
c     dps 18-oct-88 15.06 change by Bateman and Chang: add trapped
c                    particle driven current, ajtpbi(j).
c     dps 26-sep-88 15.04 alter scrape-off current calculation to allow
c                   cimprd(1) specification of 1 / tau-parallel
c     dps 11-may-87 fix bug in vloopi(j,2) bootstrap current
c     dps 07-may-87 corrections to bootstrap; add on-off switches
c     dps 04-mar-87 add bootstrap current cf. Mike Hughes
c     rgb 3-sep-86 rearranged section "change resistivity
c                    inside q=1 radius"
c     rgb 18-aug-86 compute vloopi(j,1) and ajtori(j,1) directly from
c                    bpoli(j) at the beginning of the subrtn (do 22 ...)
c     rgb 6-feb-86 OH multiplied by 1. / twopi, error corrected
c     rgb 16-dec-85 bpoli(1) = - bpoli(3)
c       drm 4-dec-84 changed formula for external inductance, add
c       current flattening code, use zbetai
cdoc
c=======================================================================
c       ------------
c       sbrtn SOLVEB   file DSOLVER
c       ------------
c
c      2.8     solve equation for b-poloidal and add ohmic heating term
c
c      common blocks and variables modified:
c
c       bpoli (comsta)
c
c-----------------------------------------------------------------------
c
c
c       in baldur internal units the b-poloidal equation is:
c
c       (d/dt) b  =  (d/dr)( (eta/r)  * (d/dr) (r*b) / emu0 - (eta jb) )
c
c       (where jb is the beam current cjbeam*ajbs )
c       the "b" on the right is taken to be
c       ztheta * b(new) + (1 - ztheta) * b(old)
c
c       we also have boundary conditions
c
c       bpa0 * (d/dr)b + bpb0 * b = bpc0
c
c       at the center, and
c
c       bpa1 * (d/dr)b + bpb1 * b = bpc1
c
c       at the edge
c
c       we difference these equations, boundary by boundary,
c       to form at each zone boundary an equation of the form
c
c       zp*bpoli(j-1) + zq*bpoli(j) + zr*bpoli(j+1) = zs
c
c       where the bpoli are the new timestep values,
c      terms with bpoli at the old timestep being included in zs.
c
c       in effect, what we are doing is solving a matrix of the form
c
c       zq(2)   zr(2)   0       0       .....           bpoli(2)        zs(2)
c       zp(3)   zq(3)   zr(3)   0       0....           bpoli(3)        zs(3)
c       0       zp(4)   zq(4)   zr(4)   0....     *     bpoli(4)   =    zs(4)
c       .. . . . . .... . . . . . . . . . .             ...             ....
c       ......  0       0       zp(m)   zq(m)           bpoli(m)        zs(m)
c
c       this is reduced by row operations to a matrix equation
c
c       1       -ze(2)  0       0       .......         bpoli(2)        zf(2)
c       0       1       -ze(3)  0       0.....          bpoli(3)        zf(3)
c       0       0       1       -ze(4)  0...      *     bpoli(4)   =    zf(4)
c       . .  .  .  .  .  .  .  .  .  .  .  .            . . . .         . . .
c       . . . . 0       0       0       1               bpoli(m)        zf(m)
c
c       this may be solved by back substitution from j = m to 2
c       the row operations may be represented by the recursion formula
c
c       ze(j) = -zr(j) / ( zq(j) + zp(j) * ze(j-1) )
c       zf(j) = ( zs(j) -zp(j) * zf(j-1) ) / ( zq(j) + zp(j) * ze(j-1) )
c
c       and the back substitution by the formula
c
c       bpoli(j) = ze(j) * bpoli(j+1) + zf(j)
c
c
c       in this subroutine, the zp, zq, zr, and zs
c       are computed for each boundary and immediately used to generate
c       ze and zf; the row reduction loop, therefore, takes eta, the old
c        value of bpoli, and various mesh parameters and generates
c       ze and zf.
c       ze and zf are then used to compute bpoli at the new time.
c
c       it should be noted that if we take zp(2) and zr(m) to be 0,
c       we need no special recursion formula at the center (zp=0, therefore
c       the  zp(2)*ze(1) and zp(2)*zf(1) terms are zero),
c       or at the edge (zr(m)=0, therefore ze(m) = 0).
c
c  Boostrap current
c  ----------------
c
c  ajboot(j) = < bootstrap current density / R > [Amps / m**3]
c              at BALDUR zone centers
c  ajtpbi(j) = < current density due to trapped particles in the
c              banana regime / R > [Amps / m**3]
c              at BALDUR zone centers
c              computed with an analytic formula by C. S. Chang, 1988
c
c
c  Notes on Lower Hybrid Current Drive
c  -----------------------------------
c
c      Implemented by Bateman and Voitsekhovitch on 19 November 2000.
c
c      Variables stored in file commhd.m in common block comlhcd:
c
c  xcd(jx) = normalized radial grid points, 0 < jx <= Kjbal == 55
c  tcd(jt) = time slices (sec), 0 < jt <= Knave == 20
c  cdrive(jx,jt) = current drive, A/m^2
c  plhcde(jx,jt) = power density deposited to electrons, W/m^3
c  nxlhcd  = number of normalized radial grid points used
c  ntlhcd  = number of time slices used
c
c
c      Default values:
c
c  nxlhcd  = 0
c  ntlhcd  = 0
c  xcd(jx) = -1.0
c  tcd(jt) = -1.0e10
c  cdrive(jx,jt) = 0.0
c  plhcde(jx,jt) = 0.0
c
c
cend
c**********************************************************************c
c
        subroutine solveb
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c
        dimension
     r   ze(111)      , zf(111)      , zetai(mj), zbetai(mj)
     &  , zr1(mj), zr2(mj), zr3(mj),  zbias(mj)
     &  , zlh(mj), zlhp(mj), zlhup(mj),  zfw(mj), zrfcd(mj)
     &  , ztemp(mj)
c
cgb        equivalence
cgb     x   (nlzzzz(1),ze(1))       , (nlzzzz(112),zf(1))
c
c-----------------------------------------------------------------------
c
cahk
       namelist/lhcurrnt/ zlhp, zlhcu, zfwcu
cahk

        data    iclass /2/,     isub /8/
        data ihton,ihtoff/61,62/
        data ilext,ibboun,itheta,iscrap,igoose / 395,396,397,398,399/
        data iqflat/ 478/
        data zbboun/-1./
c
        if(.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.8 subroutine solveb bypassed')
        return
   10   continue
c
       if ( cfutz(499) .eq. 1 ) then
cahk  Code added to include lower hybrid and fast wave current
cahk  ktest so that name list is read only first time through
        ztau= tai * uist
       if (ztau .lt. 3.0) goto77
       if (ktest .gt. 1) goto 77
c       namelist/lhcurrnt/ zlhp, zlhcu, zfwcu
       open (21, file='lhcurr', status='old')
       read (21, lhcurrnt)
       close (21)
       write(6,*) ztau
       ktest = 2
cahk  The next do loop is used to generate a profile for the driven
cahk  lower hybrid current.  It is only used to produce zlhp(mj) for
cahk  for the namelist lhcurrnt in the file lhcurr
       do 89 j=1, mzones
       zlhup(j)=exp(-((0.5*(avi(j+1,1,1)+avi(j,1,1))-
     +  0.6* avi(mzones,1,1))/(.16*avi(mzones,1,1)))**2)
 89     continue
       open (29, file='zlhprof', status='unknown')
       write(29,*) ztau
       write (29,*) zlhp, zlhcu, zlhup
       if (abs(zlhcu) .lt. 1.e-16 .and. abs(zfwcu) .lt. 1.e-16) goto 77
       zsumlh = 0.0
       zsumfw = 0.0
       do 76 j=1, mzones-1
       zsumlh = zsumlh + zlhp(j)*(avi(j+1,13,1) - avi(j,13,1))
       zsumfw = zsumfw + rleprf(j)*(avi(j+1,13,1) - avi(j,13,1))
 76    continue
       if (zsumlh .lt. 1.0e-16 .and. zsumfw .lt. 1.e-16 ) goto 77
       if (zsumlh .ge. 1.0e-16) zj0lh = zlhcu/zsumlh
       if (zsumfw .ge. 1.0e-16) zj0fw = zfwcu/zsumfw
       do 78 j=1, mzones-1
          zlh(j)=zj0lh  * zlhp(j)
          zfw(j)=zj0fw * rleprf(j)
          zrfcd(j)=zlh(j) + zfw(j)
 78    continue
       do 97 j=1, mzones
          write (29,*) j, avi(j,1,1),zlhup(j),zlh(j),zfw(j),zrfcd(j)
 97    continue
       close(29)
       endif

 77    continue
c-----------------------------------------------------------------------
C.MHH  Dec 23 1986 - fetch bootstrap current
c
c  Note that when lneocl(1) == 1, the bootstrap current is computed
c    in sbrtn nclass_int
c
      if ( lneocl(1) .lt. 1 ) then
        ajboot = 0.0
        ajtpbi = 0.0
        call boots (ajboot,ajtpbi,lcentr,ledge)
      endif
c
c singer add isep and times for current drive profile
c
      isep = mzones
      if (nadump(1) .gt. lcentr) isep = nadump(1)
      zon=cfutz(ihton)*usit
      zoff=cfutz(ihtoff)*usit
c
c zone-independent quantities for flr cutoff of chang current
c
      zb=bzs*usib
      zckb=cfev*10.**(fxk)
      zcmp=fcmp*10.**(fxnucl)*usim
      zce=fce*10.**(fxe)*10.0
c
c bateman  temporary variables needed for magnetic diffusion eqn
c
      do 20 j=1,mzones+1
      zr1(j) = avi(j,4,1) * avi(j,10,1) * r0ref / avi(j,9,1)
      zr2(j) = avi(j,2,1)*avi(j,3,1)*avi(j,7,1)*r0ref / avi(j,8,1)
      zr3(j) = 1.0/(r0ref*avi(j,10,1))
      zbias(j) = 0.5
  20  continue
c
c..compute current drive profile
c
c
c      Variables stored in file commhd.m in common block comlhcd:
c
c  xcd(jx) = normalized radial grid points, 0 < jx <= 55
c  tcd(jt) = time slices (sec), 0 < jt <= 20
c  cdrive(jx,jt) = current drive, A/m^2
c  plhcde(jx,jt) = power density deposited to electrons, W/m^3
c  cdcoef  = coefficient of cdrive (default 1.0)
c  plhcoef = coefficient of plhcde (default 1.0)
c  cdprof(jx) = profile of current drive, A/m^2
c               (computed in sbrtn solveb)
c  plhprof(jx) = profile of power density deposited to electrons, W/m^3
c  nxlhcd  = number of normalized radial grid points used
c  ntlhcd  = number of time slices used
c            (computed in sbrtn eqinit in file deqbald.f)
c
c     Variables xcd, tcd, cdrive, plhcde, cdcoef, and plhcoef
c  are read from the namelist input in sbrtn eqinit in file deqbald.f.
c
c..first, zero out the current drive profile
!cap
      cdprof = 0.0
      ztemp  = 0.0
c
c..interpolate profile only if nxlhcd > 0 and ntlhcd > 0 and
c  abs(cdcoef) > epslon
!cap
      ztime = 0.5 * ( tai + tbi ) * uiet
c
      if ( nxlhcd .gt. 0  .and.  ntlhcd .gt. 0
     &     .and.  abs(cdcoef) .gt. epslon ) then
c
c..interpolate in time
c  The result is ztemp(jz) on the same grid as xcd(jz)
c
        if ( ztime .le. tcd(1) ) then
          do jz=1,nxlhcd
            ztemp(jz) = cdcoef * cdrive(jz,1)
          enddo
        elseif ( ztime .ge. tcd(ntlhcd) ) then
          do jz=1,nxlhcd
            ztemp(jz) = cdcoef * cdrive(jz,ntlhcd)
          enddo
        else
c
          zint = 0.0
          do jt=2,ntlhcd
            it = jt
            if ( tcd(jt) .ge. ztime ) go to 12
          enddo
c
 12       continue
          zint = ( tcd(it) - ztime ) / (  tcd(it) -  tcd(it-1) )
c
          do jz=1,nxlhcd
            ztemp(jz) = cdcoef * ( cdrive(jz,it-1) * zint
     &                + cdrive(jz,it) * ( 1.0 - zint ) )
          enddo
c
        endif
c
c..Now, interpolate in radius to compute cdprof
c
        ix = 1
        call int1d ( 2, 0, 0, nxlhcd, 1, xcd, ztemp
     &    , mzones, xzoni, ix, cdprof )
c
      endif
!cap
 99   continue
c
cdoc
c
c..compute
c  vloopi(j,1) = loop voltage (volts) at zone center j and timestep n
c  vloopi(j,2) = same at timestep n+1
c  ajtori(j,1) = < J-toroidal / R > (Amp/m**3) at zone center j timestep n
c  ajtori(j,2) = same at timestep n+1
c  weohms(j) = Ohmic heating source term [Joules / m**3]
c                  at zone center j and timestep n+1/2
c  for now, the above are in BALDUR internal units and standard units
c
c  zrb = V'(xi)*rho'(xi) * r0ref * bpoli * <|del xi|**2 / R**2 >
c                  at zone boundary j and timestep n+1/2
c  zrbp1 = same at zone boundary j+1
c
cend
      zrbp1 = 0.
c
      do 22 j=lcentr,mzones-1  ! zone center index
c
      zrb = zrbp1
      zrbp1 = avi(j+1,3,1) * avi(j+1,2,1) * r0ref * bpoli(j+1)
     &  * avi(j+1,7,1)
c
       zfrad=max(epslon,armins(j,1)/armins(isep,1))
       zrf=sfutz(18)*(max(epslon,zfrad**sfutz(19)))**sfutz(20)
       if(tai.lt.zon.or.tai.gt.zoff) zrf=0.0
c
cahk Add lhcd current and the fast wave current
c
       if (cfutz(499) .eq. 1 ) then
       zrf=zlh(j)*0.5*(avi(j,11,1)+avi(j-1,11,1)) +
     &     zfw(j)*0.5*(avi(j,11,1)+avi(j-1,11,1))
      else
       zrf=0.0
      endif
c
      vloopi(j,1) = twopi * usir * eta(2,j) * ( avi(j,9,1)
     &  * (zrbp1 / avi(j+1,8,1) - zrb / avi(j,8,1))
     &     / (emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
     &  - cjbeam * ajbs(j) * usij * 0.5*(avi(j+1,11,1)+avi(j,11,1))
     &  - cdprof(j) * 0.5*(avi(j+1,11,1)+avi(j,11,1))
     &  - cfutz(480) * ( ajboot(j) + sfutz(16)*ajtpbi(j) ) - zrf )
     &  / avi(j,10,1)
c
      ajtori(j,1) = (zrbp1 - zrb)
     & / (emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
c
  22  continue
c
c
cend bateman 7-apr-85
c
c
c      1)      b-poloidal boundary conditions
c
c  Use either thetai or cfutz(397) for the implicitness parameter in solveb
c
        ztheta=thetai
        if(cfutz(397).gt.0.) ztheta=min(1.,cfutz(397))
c
c       if the boundary condition switch has been changed to 0 set ibchng
c
        ibchng=0
c
c       1.1)    find interpolated value of rcurri
c
c       find the outermost q=1 surface
c
c
        iq=0
        do 142 jsudo=1,mzones
        jz=mzones+1-jsudo
        if(q(jz).ge.1.0) go to 142
        iq=jz
        go to 144
142     continue
144     continue
        if(cfutz(478).le.0.) iq=0
c
c       change resistivity inside the q=1 radius
c
        do 145 jz=lcentr,mzones
        zbetai(jz)=etai(jz)
 145  continue
c
        if(iq.gt.lcentr) then
        zqx=xbouni(iq)+dxzoni(iq)*(1.-q(iq))/(q(iq+1)-q(iq))
        zqeta=etai(iq)+(etai(iq+1)-etai(iq))*(zqx-xbouni(iq))/
     1  dxzoni(iq)
        zdip=min(cfutz(478),2.)
c
        do 146 jz=lcentr,mzones
        if(jz.le.iq)
     &  zbetai(jz)=zqeta/(1.+zdip*(1.-(xbouni(jz)/zqx)**2))
146     continue
        end if
c
c
c
c       1.2)    set boundary condition coefficients
c
c
cbate        bpa0=0.
cbate        bpb0=1.
cbate        bpc0=0.
c
cbate        bpa1 = 0.0
cbate        bpb1 = 1.0
cbate        bpc1 = zrcurr / (rmini * rmaji)
c
cdoc
c
c  the flux surface averaged poloidal magnetic field at the outer bndry
c  of the plasma is determined here from the total plasma current
c  eqcamp = total toroidal plasma current in amps
c  bpoli = twopi * emu0 * eqcamp
c          / (V'(xi) * rho'(xi) * r0ref * <|del xi / R|**2>
cend
      jb = mzones + 1
      bpoli(jb) = twopi * emu0 * eqcamp
     &  / (avi(jb,3,1) * avi(jb,2,1) * r0ref * avi(jb,7,1))
c
c
c      2)      row reduction loop
c
c
        if(nstep.gt.1) go to 218    ! only on first step
        zpmass=fcmp*(10.0**fxnucl)  ! proton rest mass
        ze0=fces*(10.0**fxes)       ! proton charge e
        zetai = 0.0                 ! zero out zetai
  218   continue
c
        ze(lcentr-1)=0.0  ! set inner ghost points
        zf(lcentr-1)=0.0
c
c  main loop over zone centers
c
        do 280 jz = lcentr, mzones
c
c      2.1)    difference center boundary condition
c
c   difference bpa0(d/dr) bpoli + bpb0*bpoli(lcentr) = bpc0
c   zp * bpoli(1) + zq * bpoli(2) + zr * bpoli(3)  =  zs
c
      if ( jz .eq. lcentr ) then
c
        zp = 0.0   ! bpoli(2) = 0. at magnetic axis
        zr = 0.0
        zq = 1.0
        zs =  0.0
        go to 260
c
      endif
c
c      2.2)    difference edge boundary condition
c
c  bpa1(d/dr) bpoli + bpb1 * bpoli(ledge+1) = bpc1
c  zp * bpoli(mzones-1) + zq * bpoli(mzones) * zr * bpoli(mzones+1) = zs
c
      if ( jz .eq. mzones ) then
c
      zp = 0.
      zq = 1.
      zr = 0.0
c
      zs = twopi * emu0 * eqcamp
     & / (avi(jz,3,1) * avi(jz,2,1) * r0ref * avi(jz,7,1))
cdoc
c
c  eqcamp = total toroidal current within outer zone boundary of plasma
c           computed within equilibrium package for time 0.5 * (tai + tbi)
cend
            go to 260
      endif
c
c      2.3)    difference equation
c
c 1D    (d/dt) b  =  (d/dr)( (eta/r*emu0)  * (d/dr) (r*b) - (eta jb) )
c
c  compute zetai(j), anomalous resistivity in divertor region
c
        zl=0.0
        if(cfutz(398).gt.0.) go to 255
        if(nadump(1).le.lcentr) go to 255
        if(jz.lt.nadump(1)) go to 255
        if(jz.lt.nadump(3)) zl=cfutz(127)
        if(jz.ge.nadump(3)) zl=cfutz(128)
        if(zl.gt.1.e18) go to 255
c
c  ionic charge density
c  15.04 Note that zn = total ion density, zm = total mass density in AU
c
        z0=0.0
        z1=0.0
        zn=0.0
        zm=0.0
        do 242 jh=1,mhyd
        z0=z0+rhohs(jh,2,jz)
        z1=z1+rhohs(jh,2,jz)*scroff(jh,jz)
        zn=zn+rhohs(jh,2,jz)
        zm=zm+rhohs(jh,2,jz)*aspec(jh)
  242   continue
c
c  recover plasma-flow velocity
c
        zvs=-(z1*zl)/(z0+epslon)
c
        if(mimp.le.0) go to 246
        do 244 ji=1,mimp
        z0=z0+rhois(ji,2,jz)*cmean(ji,2,jz)
        ii=ji+lhydn
        zn=zn+rhois(ji,2,jz)
        zm=zm+rhois(ji,2,jz)*aspec(ii)
  244   continue
  246   continue
c
c  15.04 If cimprd(1) is used to specify 1 / tau-parallel = vs / l,
c  calculate l using vs = sound speed
c
        if (cimprd(1).gt.0.0) then
          ziionm=zn/(zm*fcau*(10.0**fxau))
          zvs=sqrt((tis(2,jz)+tes(2,jz))*ziionm)
          zl=zvs/cimprd(1)
        end if
c
c  divertor current zjs
c
        zjs=0.25*ze0*z0*zvs
c
c  anomalous resistivity
c
        zetas=0.0
        if(zl.gt.epslon) zetas=2.*tes(2,jz)/(ze0*zl)
        zcetai=zetas*usir
c
c  check that current is below the sheath limit
c
        if(abs(ajzs(2,jz)).ge.zjs) go to 250
c**there is a question(?) whether to use abs(ajzs(2,jz)) again**********
        if(ajzs(2,jz).le.epslon) ajzs(2,jz)=.5*zjs
        zjsrat=ajzs(2,jz)/zjs
        zarg=(1.+zjsrat)/(1.-zjsrat)
        zetai(jz)=0.5*zcetai*log(zarg)/ajzs(2,jz)
        go to 255
  250   continue
        zetai(jz)=1.e+03*zbetai(jz)
  255   continue
c
        if(jz.lt.nadump(1)) zetai(jz)=0.0
        if((jz-1).lt.nadump(1)) zetai(jz-1)=0.0
c
c       enhance the resistivity used in the current diffusion equation
c       (but not the resistivity used for heating) with cfutz(igoose)
c
        zgoose=1.
        if(cfutz(igoose).gt.0.) zgoose=cfutz(igoose)
        z1eta=zgoose*((1.-cfutz(390))*zetai(jz-1)+zbetai(jz-1))
        z2eta=zgoose*((1.-cfutz(390))*zetai(jz)+zbetai(jz))
c
cdoc
c
c bateman  finite difference magnetic diffusion eqn
c
c      d ( rho'(xi) B-poloidal ) / dt
c      = d ( d(rho)/dt B-poloidal ) / d xi
c     + d ( (eta / zr1) d ( zr2 B-poloidal / emu0 ) / dxi ) / d xi
c     + d ( zr3 J-beam ) / d xi
c
c  zr1(jz) = V' <1/R**2> / (R B-toroidal)  timestep N+1/2  zone center jz
c  zr2(jz) = rho'(xi) V'(xi) < (del xi / R)**2 > / ( (R B-toroidal)
c                  at timestep N+1/2  zone bndry jz
c  zr3(jz) = 1. / (r0ref * <1./r**2> )
c
c  zbias(jz) = 0.5 for now
c                  may be used later for flow biasing
c
c  zdtdrh = dti / ( ( d rho / d xi at step N+1 boundary jz ) )
c                  / (xzoni(jz) - xzoni(jz-1))
c
c  zratio = (d rho / d xi at step N  boundary jz )
c           / ( d rho / d xi at step N+1  boundary jz )
cend
c
      zdhrho = avi(jz,2,2) * 0.5 * (tbi - tai)
      zdtdrh = dti / ( (xzoni(jz)-xzoni(jz-1))
     &            * (avi(jz,2,1) + zdhrho) )
      zratio = ( avi(jz,2,1) - zdhrho ) / ( avi(jz,2,1) + zdhrho )
c
c
      zpp = zdtdrh * ( zr2(jz-1)
     &  *   z1eta / ( emu0 * zr1(jz-1) * dxzoni(jz-1) )
     &                - avi(jz-1,1,2) * (1.-zbias(jz-1)))
c
      zqp = - zdtdrh * (zr2(jz)
     &  * ( z2eta / ( emu0 * zr1(jz) * dxzoni(jz) )
     &    + z1eta / ( emu0 * zr1(jz-1) * dxzoni(jz-1)))
     &          + avi(jz,1,2) * (zbias(jz-1) + zbias(jz) - 1.))
c
      zrp = zdtdrh * ( zr2(jz+1)
     &  *   z2eta / ( emu0 * zr1(jz) * dxzoni(jz) )
     &                 + avi(jz+1,1,2) * zbias(jz) )
c
        zp=-ztheta*zpp
        zq=1.-ztheta*zqp
        zr=-ztheta*zrp
c
C.MHH  Dec 23 1986 - beam current and bootstrap current respectively
c
c..beam current contributions
c
       zjbmp = cjbeam * usij * ajbs(jz)
     &         * 0.5*(avi(jz+1,11,1)+avi(jz,11,1))
       zjbmm = cjbeam * usij * ajbs(jz-1)
     &         * 0.5*(avi(jz,11,1)+avi(jz-1,11,1))
c
       zfrad=max(epslon,armins(jz,1)/armins(isep,1))
       zrf=sfutz(18)*(max(epslon,1.0-zfrad**sfutz(19)))**sfutz(20)
       if(tai.lt.zon.or.tai.gt.zoff) zrf=0.0
c
cahk  Add lower hybrid and fast wave current
c
      jzm2 = max ( jz-2, 1 )
c
      zrfp = 0.0
      zrfm = 0.0
c
      if (cfutz(499) .eq. 1) then
        zrfp=zlh(jz)*0.5*(avi(jz,11,1)+avi(jz-1,11,1)) +
     &     zfw(jz)*0.5*(avi(jz,11,1)+avi(jz-1,11,1))
        zrfm=zlh(jz-1)*0.5*(avi(jz-1,11,1)+avi(jzm2,11,1)) +
     &     zfw(jz-1)*0.5*(avi(jz-1,11,1)+avi(jzm2,11,1))
      endif
c
      zrfp = zrfp + cdprof(jz) * 0.5*(avi(jz,11,1)+avi(jz-1,11,1))
      zrfm = zrfm + cdprof(jz-1)*0.5*(avi(jz-1,11,1)+avi(jzm2,11,1))
c
c..bootstrap current contributions
c
      zjbtp = cfutz(480) * ( ajboot(jz) + sfutz(16)*ajtpbi(jz) )
      zjbtm = cfutz(480) * ( ajboot(jz-1) + sfutz(16)*ajtpbi(jz-1) )
C
      zs = zratio * bpoli(jz) + (1.-ztheta) * 
     &  ( zpp * bpoli(jz-1) + zqp * bpoli(jz) + zrp * bpoli(jz+1) )
     &   - zdtdrh * ( z2eta * zr3(jz) * ( zjbmp + zjbtp + zrfp )
     &      - z1eta * zr3(jz-1) * ( zjbmm + zjbtm + zrfm ) )

c
  260   continue
c
cend bateman 7-apr-85
c
c      2.4)    compute ze and zf
c
        zt=zq+zp*ze(jz-1)
c
        if (abs(zt).le.epslon) go to 9020
        ze(jz)=-zr/zt
        zf(jz)=(zs-zp*zf(jz-1))/zt
  280   continue
c
c
c   3)  back substitution and add ohmic heating term to dddd and weohms
c
c
c       internal curr. dens. to standard
c
c       current density is current / area
c       current = uisi*(r*b at bound. jz+1 - r*b at bound. jz)
c       1 / area = dx2inv / ( 2 pi rmins**2)
c
c
cbate        zvols = 2.0 * (rmins*fcpi)**2 * rmajs
cbate        z0 = uisi*rmini / (rmins**2 * 2.0*fcpi)
cbate        z0jb = cjbeam / z0
cbate        zwohm = zvols * uist * usie
cbate        zebi = zvols * uisb**2 * usie / (8.0*fcpi)
cbate        zwpoyn = fcpi * uisl**3 * uisb**2 * usie * rmaji
c
        zsi = usid * usie / usit
c
        wbi = 0.0
        ebtoti = 0.0
        wohmi = 0.0
        zboldi = bpoli(mzones+1)
c
        do 380 jsudo = lcentr, mzones
        jz = lcentr + mzones - jsudo
c
c              if repeated timestep (dddd not recomputed) subtract out
c               old weohms contribution
c
        if(nlrpet.and..not.nliter)
     1          dddd(lelec,jz) = dddd(lelec,jz) - weohms(jz)*zsi
c
cbate        zdxb1 = dx2inv(jz)*(xbouni(jz+1)*zboldi - xbouni(jz)*bpoli(jz))
cbate        zjolds = zdxb1 * z0
c
c
c              conservation checks
c
c
cbate        if (jz.eq.ledge) zbpoyn = 0.5 * (zboldi + bpoli(jz+1))
cbate        if(jz.eq.ledge) zbpdot=(bpoli(jz+1)-zboldi)/dti
cbate        wbi = wbi - zebi *      ( bint(1,jz) * bpoli(jz)**2
cbate     1                  + bint(2,jz) * bpoli(jz)*zboldi
cbate     2                  + bint(3,jz) * zboldi**2)
c
c
c               solve for b at time nstep+1
c
c
        zboldi = bpoli(jz)
        bpoli(jz) = bpoli(jz+1)*ze(jz) + zf(jz)
c
c
c               curl b
c
c
cbate        zdxb2 = dx2inv(jz) *
cbate     1          (xbouni(jz+1)*bpoli(jz+1) - xbouni(jz)*bpoli(jz))
cbate        zjnews = zdxb2 * z0
c
c               ohmic heating
c
cbate        weohms(jz) = eta(2,jz) * 0.5*(zjnews + zjolds) *
cbate     1  (ztheta*zjnews + (1.0 - ztheta)*zjolds - cjbeam*ajbs(jz))
cbate        dddd(lelec,jz) = dddd(lelec,jz) + weohms(jz)*zsi
c
c
c               conservation checks
c
cbate        if (jz.eq.ledge) wpoyni = zwpoyn * zbpoyn *(zbetai(jz) *
cbate     1  ((1.0 - ztheta)*zdxb1 + ztheta*zdxb2 - z0jb*ajbs(jz))+
cbate     2  zbpdot*rmini**2*(xbouni(mzones)-xzoni(ledge))/zgoose)
cbate        z3 = bint(1,jz)*bpoli(jz)**2 + bint(2,jz)*bpoli(jz)*bpoli(jz+1)
cbate     1          + bint(3,jz)*bpoli(jz+1)**2
cbate        wbi = wbi + zebi * z3
cbate        ebtoti = ebtoti + zebi * z3
cbate        if(jz.le.ledge) wohmi = wohmi + weohms(jz) * dx2i(jz)*2.0*zwohm
c
  380   continue
c
c..compute
c  vloopi(j,1) = loop voltage (volts) at zone center j and timestep n
c  vloopi(j,2) = same at timestep n+1
c  ajtori(j,1) = < J-toroidal / R > (Amp/m**3) at zone center j timestep n
c  ajtori(j,2) = same at timestep n+1
c  weohms(j) = Ohmic heating source term [Joules / m**3]
c                  at zone center j and timestep n+1/2
c  for now, the above are in BALDUR internal units and standard units
c
c  zrb = V'(xi)*rho'(xi) * r0ref * bpoli * <|del xi|**2 / R**2 >
c                  at zone boundary j and timestep n+1/2
c  zrbp1 = same at zone boundary j+1
c
      zfac = 1. / twopi
      zis  = uisd * uise / uist
      zrbp1 = 0.
c
      do 420 j=lcentr,mzones-1  ! zone center index
c
      zrb = zrbp1
      zrbp1 = avi(j+1,3,1) * avi(j+1,2,1) * r0ref * bpoli(j+1)
     &  * avi(j+1,7,1)
c
       zfrad=max(epslon,armins(j,1)/armins(isep,1))
       zrf=sfutz(18)*(max(epslon,zfrad**sfutz(19)))**sfutz(20)
       if(tai.lt.zon.or.tai.gt.zoff) zrf=0.0
c
cahk add lhcd and fast wave current
c
      zrf = 0.0
      if (cfutz(499) .eq. 1) then
        zrf=zlh(j)*0.5*(avi(j,11,1)+avi(j-1,11,1)) +
     &    zfw(j)*0.5*(avi(j,11,1)+avi(j-1,11,1))
      endif
c
      vloopi(j,2) = twopi * usir * eta(2,j) * ( avi(j,9,1)
     &  * (zrbp1 / avi(j+1,8,1) - zrb / avi(j,8,1))
     &     / (emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
     &  - cjbeam * ajbs(j) * usij * 0.5*(avi(j+1,11,1)+avi(j,11,1))
     &  - cdprof(j) * 0.5*(avi(j+1,11,1)+avi(j,11,1))
     &  - cfutz(480) * ( ajboot(j) + sfutz(16)*ajtpbi(j) ) - zrf )
     &  / avi(j,10,1)
c
      ajtori(j,2) = (zrbp1 - zrb)
     & / (emu0 * (avi(j+1,12,1) - avi(j,12,1)) )
c
      zeohmi = zfac * 0.5 * (ajtori(j,1) + ajtori(j,2))
     &  * ((1.-ztheta) * vloopi(j,1) + ztheta * vloopi(j,2))
c
c  add Ohmic heating contribution to electron energy source term
c
      dddd(lelec,j) = dddd(lelec,j) + zeohmi   ! internal units
      weohms(j) = zeohmi * zis                 ! standard units
c
 420  continue
c
c
c               conservation checks
c
c
cbate        wbi = wbi / dti
c
c       adjust eohmi and epoyni so that bcons check works with goosing
c
cbate        wohmi=zgoose*wohmi
cbate        wpoyni=zgoose*wpoyni
cbate        eohmi = eohmi + wohmi * dti
cbate        epoyni = epoyni + wpoyni * dti
cbate        bcons=(ebini+ecompi+epoyni-eohmi)/ebtoti-1.0
c
c
c      4)      set dummy boundary values of b
c
        bpoli(1) = - bpoli(lcentr+1)
        bpoli(mzones+1) = bpoli(mzones)
     &  * avi(mzones,2,1) * avi(mzones,3,1) * avi(mzones,7,1)
     &  / (avi(mzones+1,2,1)*avi(mzones+1,3,1)*avi(mzones+1,7,1))
c              ! zero current density in zone number mzones
c
      vloopi(1,2) = vloopi(2,2)
      ajtori(1,2) = ajtori(2,2)
c
      vloopi(mzones,2) = vloopi(mzones-1,2)
      ajtori(mzones,2) = 0.
      weohms(mzones) = 0.
c
c
c
c
c
        return
c
c
c      90)     errors: singularity in equations
c
c
 9020   continue
c
        call error_olymp(1,iclass,isub,2,
     1          ' *** error *** singularity in b-poloidal')
        call error_olymp(1,iclass,isub,2,
     1          ' *** equations')
        call error_olymp(3,jz,2,1,' zone   ')
        call error_olymp(3,zp,1,1,'   zp   ')
        call error_olymp(3,zq,1,1,'   zq   ')
        call error_olymp(3,zr,1,1,'   zr   ')
        call error_olymp(2,ze(jz-1),1,1,' ze(jz-1) ')
        return
        end

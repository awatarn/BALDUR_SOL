c--------1---------2---------3---------4---------5---------6---------7-c
c@neudep  /11040/baldur/code/bald/dimprad.f
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rgb 31-jan-95 dimension pion(1), piions(1) -> pion(*), piions(*)
c  dps 22-oct-88 15.06 Remove check on zsigv.
c  dps 15-aug-88 Subroutine formed to replace statements in IMPRAD.
c***********************************************************************
c
c
      subroutine neudep(pnergy,kp,pslant,ppath,ksepx,pstart,pion,
     1                  piions,ptotal)
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c     ------------
c     sbrtn NEUDEP
c     ------------
c
c     Input variables:
c
c     PNERGY    = Energy of neutral species in keV.
c
c     KP        = BALDUR index of neutral species.
c
c     PSLANT    = Prescribed cosine of the angle between the inward-
c                 directed normal to the outer plasma surface and the
c                 influx direction of the neutral impurity.
c
c     PPATH     = Prescribed penetration limit for the neutral impurity
c                 expressed as a fraction of the plasam (or separatrix)
c                 radius.
c
c     KSEPX     = Index of plasma edge or separatrix.
c
c     PSTART    = Determines the radial index of the outermost zone
c                 allowing ionization of the neutral impurity.
c
c     PION(JZ)  = Array of ionization rates in cm**3 / sec.
c
c     Output varibles:
c
c     PIIONS(JZ)= Relative impurity-ion source profile.
c
c     PTOTAL    = Total number of particles represented by PIIONS
c                 divided by twice the plasma volume.
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      dimension pion(*), piions(*), zimpnu(mj)
c
c  Velocity in cm/sec
c
      zv0 = sqrt((2.*pnergy*1.e3*evs)/(aspec(kp)*fcau+epslon))
     1            *10.0**(-0.5*fxau)
      ziv0 = 1./(zv0 + epslon)
c
c  Path parameters
c
      zcosin = 1.0/pslant
      zdepth = ppath * ahalfs(ksepx,1)
      zminf = 1.0e-05
      jmax = ledge - int(pstart + 0.1)
      zpath = 0.0
      ptotal = 0.0
      piions(1) = 0.0
c
c  No source for zones outside maximum radius of ionization
c
      do 10 j=mzones,jmax,-1
        piions(j) = 0.0
   10 continue
c
      zimpnu(mzones) = 1.0
      zimpnu(jmax+1) = 1.0
c
      do 50 j=jmax,lcentr,-1
        if (zimpnu(j+1).lt.zminf) go to 30
        if (zpath.gt.zdepth) go to 30
        zdpath = zcosin * (ahalfs(j+1,1) - ahalfs(j,1))
c
c  Calculate path length
c
        zpath = zpath + zdpath
        zsigv = pion(j)
c
c  15.06 The consequences of this check on zsigv are unnecessary and
c  unphysical. The check on zioniz below supersedes it and makes more
c  sense.
c
        if ((versno.le.15.05).and.(zsigv.le.1.e-20)) go to 30
        zioniz = rhoels(2,j) * zsigv * ziv0 * zdpath
        zioniz = max(zioniz,0.0)
        if (zioniz.le.1.e-6) then
          zioniz = 0.0
          zimpnu(j) = zimpnu(j+1)
        else if (zioniz.lt.14.) then
c
c  zimpnu is the ion deposition profile normalized to 1 at the edge
c  (on zone boundaries)
c
          zimpnu(j) = zimpnu(j+1) * exp(-zioniz)
        else
          zimpnu(j) = 0.0
        end if
        go to 40
c
   30   continue
        zsigv = 0.0
        zioniz = 0.0
        zimpnu(j) = 0.0
c
   40   continue
c
c  Calculate relative source rate on zone centers.
c
        piions(j) = ((zimpnu(j+1)*(xzoni(j) - xbouni(j))
     1                + zimpnu(j)*(xbouni(j+1) - xzoni(j)))/dxzoni(j))
     2              * rhoels(2,j) * zsigv
c
c  Integrate to find total number of particles deposited.
c
        ptotal = ptotal + piions(j) * dx2i(j)
c
   50 continue
c
      return
      end

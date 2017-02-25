c--------1---------2---------3---------4---------5---------6---------7-c
c@ncsorc  /11040/baldur/code/bald/ncsorc.f
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rgb 02-feb-95 call ncinfl to compute impurity influx rate
c    to make this routine consistent with sbrtn imprad
c       dps 03-nov-89 15.15 fix "influx instability" by changing flout
c                     and flscr from particle losses to loss rates.
c       dps 16-aug-89 15.13 create subroutine out of lines in NCDATA
c***********************************************************************
c
        subroutine ncsorc(k,ji)
c
c           source calculation for Impurity Rate Equations package
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'commhd.m'
      include 'comncr.m'
c
c-----------------------------------------------------------------------
!cap
      real ::  zn0(mj), znrgfx(kncimp), zslnfx(kncimp),
     1           zpthfx(kncimp), zstrfx(kncimp), zrscfx(kncimp),
     2           zionrc(mj)    , ztotfx(kncimp), zrfxfx(kncimp),
     3           zionfx(mj,kncimp), zinflx(kncimp)
c
      data iheflx, icfout, idvout, inpath /70, 71, 72, 73/,
     1     islant, itmhe0, istart /75, 76, 77/,
     2     impneu, nergy1, nergy2, istar1 /200, 201, 202, 207/,
     3     istar2, inpth1, inpth2, islnt1 /208, 209, 210, 211/,
     4     islnt2 /212/
c
      save jsepx , zunit , zrfxfx, zrscfx, znrgfx, zslnfx, zpthfx,
     1     zstrfx, zrfxrc, zrscrc, znrgrc, zslnrc, zpthrc, zstrrc
c
      if (k.gt.1) go to 100
c
c  Reset cfutz factors for neutral influx calculation where necessary
c
      if (ji.eq.1) then
        call neuset(cfutz,epslon,epsinv,1)
        if (cfutz(iheflx).gt.epslon) call neuset(cfutz,epslon,epsinv,2)
      end if
c
c  Neutral influx parameters
c
      jsepx = mzones
      if (nadump(1).gt.lcentr) jsepx = nadump(1)
      zunit = 1. / (uisl*uisl*uist)
      zrfxfx(ji) = recflx
      zrscfx(ji) = recscr
c
      if (ji.eq.1) then
        znrgfx(ji) = cfutz(nergy1)
        zslnfx(ji) = min(cfutz(islnt1),1.0)
        zpthfx(ji) = cfutz(inpth1)
        zstrfx(ji) = cfutz(istar1)
      else
        znrgfx(ji) = cfutz(nergy2)
        zslnfx(ji) = min(cfutz(islnt2),1.0)
        zpthfx(ji) = cfutz(inpth2)
        zstrfx(ji) = cfutz(istar2)
      end if
c
      if (ji.eq.int(cfutz(iheflx)+0.1)) then
        zrfxfx(ji) = 0.0
        zrscfx(ji) = 0.0
        zrfxrc = cfutz(icfout)
        zrscrc = cfutz(idvout)
        znrgrc = cfutz(itmhe0)
        zslnrc = min(cfutz(islant),1.0)
        zpthrc = cfutz(inpath)
        zstrrc = cfutz(istart)
      end if
c
      return
c
c-----------------------------------------------------------------------
c
  100 continue
c
c  Volume source due to neutral impurity influx. Standard
c  specification via flimp is implemented here as in IMPRAD,
c  but recycling terms are added. A generalization of the cold
c  helium recycling is also included. Instead of just helium, these
c  parameters now refer to impurity number cfutz(iheflx = 70)
c
c  General impurity influx and recycling
c
      ip = ji + lhydn
c
c
c..15.15 Fix "influx instability"; preserve original scheme
c        for backward compatibility
c
      zrflx = zrfxfx(ji)*flout(ji) + zrscfx(ji)*flscr(ji)
      if (versno.lt.15.15) zrflx = zrflx/(delts*ncrept/fcstep)
c
      if (k.eq.2) call neudep(znrgfx(ji),ip,zslnfx(ji),zpthfx(ji),
     1            jsepx,zstrfx(ji),s0(1,ji),zionfx(1,ji), ztotfx(ji))
c
      jmax = ledge - int(zstrfx(ji) + 0.1)
c
      zvols = avi(mzones,12,1) * uisl**3
      zsurfs = avi(jmax+1,3,1) * avi(jmax+1,5,1) * uisl**2
c
      call ncinfl ( zsurfs, zvols, zinflx )
c
c
c  This factor zvols effectively multiplies the normalized volume
c  element dx2i used in computing ztotfx. We shift zvols in time
c  to try to compensate for time variations in the volume without
c  actually altering dx2i.
c  15.02 Multiply specified source (zinflx) by influxing area, zsurfs
c  (determine by source radius), but not recycled flux
c     
      znorm = (zinflx(ji)*zsurfs + zrflx) / (ztotfx(ji) *
     1        2.*zvols + epslon)
c
c  15.07 As done in NEUDEP, hold xn0 constant in zones where
c  s0 = 0. Since there is never ionization at mzones, set xn0 there
c  first and hold constant as we go in until s0 significant.
c
      do 110 j=mzones,lcentr,-1
        dq(j,ji) = zionfx(j,ji) * znorm
c
        if (j.eq.mzones) then
          xn0(j,ji) = znorm
c
        else if ((j.le.jmax).and.(s0(j,ji).gt.epslon)) then
          xn0(j,ji) = dq(j,ji) / (rhoels(2,j) * s0(j,ji))
c
        else
          xn0(j,ji) = xn0(j+1,ji)
        end if
c
  110 continue
c
c  Recycling of impurity number cfutz(iheflx = 70)
c
      if (ji.ne.int(cfutz(iheflx)+0.1)) go to 120
c
      zrflx = zrfxrc*flout(ji) + zrscrc*flscr(ji)
      if (versno.lt.15.15) zrflx = zrflx / (delts*ncrept/fcstep)
c
      if (k.eq.2) call neudep(znrgrc,ip,zslnrc,zpthrc,jsepx,
     1            zstrrc,s0(1,ji),zionrc,ztotrc)
c
      jmax = ledge - int(zstrrc+0.1)
      zvols = avi(mzones,12,1) * uisl**3
      znorm = zrflx / (ztotrc * 2.*zvols + epslon)
c
c  15.07 Calculate contribution to neutral density as above (zdxn0) and
c  add to xn0; contribution from previous zone is zdxnp1.
c
      do 115 j=mzones,lcentr,-1
        dq(j,ji) = dq(j,ji) + zionrc(j)*znorm
c
        if (j.eq.mzones) then
          zdxn0 = znorm
c
        else if ((j.le.jmax).and.(s0(j,ji).gt.epslon)) then
          zdxn0 = zionrc(j) * znorm / (rhoels(2,j) * s0(j,ji))
c
        else
          zdxn0 = zdxnp1
        end if
c
        xn0(j,ji) = xn0(j,ji) + zdxn0
        zdxnp1 = zdxn0
  115 continue
c
  120 continue
c
      return
c
      end

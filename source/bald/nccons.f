c--------1---------2---------3---------4---------5---------6---------7-c
c@nccons  /11040/baldur/code/bald/dimprad.f
c       dps 15-may-89 15.09 remove time-centering of equilibrium quantities.
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 15-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 08/08/84: set sign of "pflx", "ploss" positiv for outflux
c       rhw 03/07/84: reset definition of ccons acc. to baldur
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine nccons(ji)
c
c       2.20.6   conservation checks of non-corona calculation
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'comncr.m'
c
         xntot(ji) = 0.0
         jimp = limp1 - 1 + ji
c
      do 25 j=lcentr,ledge
         zdvs = avi(j,4,1) *dxzoni(j)*uisl**3
         xntot(ji) = xntot(ji) + xns(j,ji)*zdvs
   25 continue
c
         zatot = - pflx(ji) - ploss(ji) + psorc(ji)
      if (xntoti(ji).gt.epslon)
     x           ccons(jimp) = (xntot(ji)-zatot)/xntoti(ji) - 1.0
      if (xntoti(ji).le.epslon.and.xntot(ji).ge.epslon)
     x           ccons(jimp) = (xntot(ji)-zatot)/xntot(ji)
c
      return
      end

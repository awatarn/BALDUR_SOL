c--------1---------2---------3---------4---------5---------6---------7-c
c@noncor  .../baldur/code/bald/noncor.f
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  rgb 19-dec-94 moved initialization to sbrtn ncinit
c       dps 07-jul-89 15.11 add neutral impurities to equations
c       dps 15-may-89 15.09 move IRE code into predictor-corrector loop.
c       dps 10-nov-88 15.07 initialize xntoti just prior to calculation.
c       dps 26-sep-88 15.03 add section that allows NONCOR to work
c                     in a diagnostic mode; activated by limprd(1).
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 12-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 03/08/84: set xn(mzones,1,ji) acc. boundary cond. at begin
c       rhw 12/07/84: insert xntoti, set xns at the beginning
c       rhw 05/06/84: drop out dqrek, dqreko, error condition.
c       rhw 17/05/84: save ilinit, insert xnini and flout, drop ilchek
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine noncor
c
c       2.20    calculate non-corona radiation
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'commhd.m'
      include 'comncr.m'
c
      dimension zweiro(mj), zsmoo1(mj), zsmoo2(mj)
c
      data    ncmini  /2/
      data    ismrad,ifrac /199,327/
c
c
      if ( mimp .le. 0 )    return
c
      if ( nlcomp )   write (nprint,10200)
c
      if ( tes(2,lcentr) .le. 0.0 )   return
      if ( rhoels(2,lcentr) .le. 0.0 )   return
c
         fcstep = 0.5
      if (tai*uist.le.tinit) then
        fcstep = 1.0
c
c  15.07 initialize xntoti prior to calculation in case first
c  time step is repeated.
c
        do 14 ji=1,mimp
          xntoti(ji) = 0.0
          do 14 j=lcentr,ledge
            xntoti(ji) = xntoti(ji) + xn(j,1,ji)
     1                    *avi(j,4,1)*dxzoni(j)*uisl**3
   14   continue
c
      end if
         delts = fcstep * (tbi-tai) * uist / ncrept
c
c  15.03 For diagnostic purposes, can decouple NONCOR
c  from the rest of BALDUR by setting all of the required variables
c  as if the impurities had their initial profile with unit charge.
c  It is activated by setting limprd(1) not = 0.
c
c..Main loop over impurity species
c
      do 30 ji=1,mimp
        nk = nkimp(ji)
c
c..Compute and smooth weirs using old n and T (skip in diagnostic mode)
c
        if (limprd(1).eq.0) then
          call ncrats(3,ji)
c
c..Smooth weirs if requested
c
          if (cfutz(ismrad).gt.epslon) then
            iorder = int(cfutz(ismrad) + 1.01)
            do 22 j=lcentr,mzones
              zsmoo1(j) = weirs(ji,j)
   22       continue
            call smooth(zsmoo1,zsmoo2,mzones,2,2,0,lcentr,
     1                  dx2i,0,0,iorder)
c
c..Repeat, as in original smoothing procedure
c
            call smooth(zsmoo2,zsmoo1,mzones,2,2,0,lcentr,
     1                  dx2i,0,0,iorder)
            do 24 j=lcentr,mzones
              weirs(ji,j) = zsmoo1(j)
   24       continue
          end if
        end if
c
c..Save resulting "old" weirs
c
        do 25 j=1,mzones
          zweiro(j) = weirs(ji,j)
   25   continue
c
c..Advance Impurity Rate Equations
c
        call ncdata (ji)
        call nccons (ji)
c
c..Smooth this new weirs
c
        if (cfutz(ismrad).gt.epslon) then
          do 27 j=lcentr,mzones
            zsmoo1(j) = weirs(ji,j)
   27     continue
          call smooth(zsmoo1,zsmoo2,mzones,2,2,0,lcentr,
     1                dx2i,0,0,iorder)
c
c..Repeat, as in original smoothing procedure
c
          call smooth(zsmoo2,zsmoo1,mzones,2,2,0,lcentr,
     1                dx2i,0,0,iorder)
          do 29 j=lcentr,mzones
            weirs(ji,j) = zsmoo1(j)
   29     continue
        end if
c
   30 continue
c
      if (xnerr.gt.(1.0+epslon))   ncrept = ncrept + 1
      if (xnerr.lt.0.33333)   ncrept = ncrept - 1
         ncrept = max0 (ncrept,ncmini)
c
c..Set dweirs = 0 since it is not used here
c
      call resetr(dweirs,mximp*mzones,0.0)
c
c..Compute theta-centered Prad as done in SOLVEB for weohms
c
      do 40 j=lcentr,mzones
c
        zwrtot = 0.
        do 35 ji=1,mimp
c
          if (limprd(1).eq.0) then
            weirs(ji,j) = (1.-thetai)*zweiro(j) + thetai*weirs(ji,j)
          else
c
c..Diagnostic mode
c
            weirs(ji,j) = 0.
          end if
c
          zwrtot = zwrtot + weirs(ji,j)
   35   continue
c
c..Add weirs to BALDUR source term when specified Prad not operative
c
        if (cfutz(ifrac).eq.0.0) dddd(lelec,j) = dddd(lelec,j)
     1     - zwrtot * uist * usid * usie
   40 continue
c
      return
c
10200 format (//,1x,9('>'),5x,'warning  ----  non-corona-routines ',
     x        'cannot handle compression at the moment',5x,9('<'),//)
      end

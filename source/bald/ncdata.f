c--------1---------2---------3---------4---------5---------6---------7-c
c@ncdata  /11040/baldur/code/bald/dimprad.f
c  rgb 29-dec-94 moved initialization to ncinit
c    removed first argument from sbrtn ncdata(k,ji)
c       dps 16-aug-89 15.13 move source calculation into NCSORC.
c       dps 15-may-89 15.09 eliminate unneeded input variables and
c                     equilibrium time-centering with move of
c                     IRE code into predictor-corrector loop.
c       dps 10-nov-88 15.07 change xn0 calculation to properly treat
c                     regions with s0 = 0 (as in NEUDEP).
c       dps 22-sep-88 15.02 multiply influx rate by area
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 12-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 03/08/84: set sign of "flout", "flscr" positiv for outflux.
c       rhw 12/07/84: preset data and read input data by namelist.
c       rhw 05/06/84: drop out dqrek, dsourc, read output selectors.
c       rhw 17/05/84: change zsion->s0, insert xnini, recycling, beta.
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine ncdata(ji)
c
c       2.20.1   data for non-corona radiation
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'comncr.m'
      include 'comadp.m'
c
c-----------------------------------------------------------------------
c
c    interpolation of rate-coefficients "s0", "sa" and "ra"
c
      call ncrats(2,ji)
c
c    Calculate volume sources due to neutral impurity influxing and recycling
c
      call ncsorc(2,ji)
c
c  Fusion sources
c
      if (nimp(ji).eq.-4) then
c
c  He-3
c
        do 121 j=lcentr,mzones
          dz(j,ji) = aoloss(j)
  121   continue
c
      else if (nimp(ji).eq.2) then
c
c  He-4
c
        do 122 j=lcentr,mzones
          dz(j,ji) = salfs(j)
  122   continue
c
      else
c
        do 123 j=lcentr,mzones
          dz(j,ji) = 0.
  123   continue
c
      end if
c
c    instreaming ionised impurities due to limiter sputtering:
c    ->  volume source for 1. ionisation-stage
c
      if (ji.gt.1)   go to 130
c
c      source for the 1st impurity
c
c           no source at the moment
c
c      resulting  source
c
      do 125 j=lcentr,ledge
         zsrc = 0.0
         dq(j,ji) = dq(j,ji) + zsrc
  125 continue
c
      go to 150
c
  130 continue
c
c      source for the 2nd impurity
c
c           no source at the moment
c
c      resulting  source
c
      do 133 j=lcentr,ledge
         zsrc = 0.0
         dq(j,ji) = dq(j,ji) + zsrc
  133 continue
c
  150 continue
c
c    calculation of diffusion and driftvelocities
c
      call ncdifu(ji)
c
c    calculate final coefficients and solve rate equations
c
      call nccoef(ji)
c
c    calculation of radiation losses
c
      call ncrats(3,ji)
c
      return
      end

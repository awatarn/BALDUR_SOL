c@islbal 13:00 15-jan-92  .../baldur/code/bald/dislbal
c  Glenn Bateman, PPPL, (609) 243 2837
c
c     ------------
c     sbrtn islbal
c     ------------
c
c   Interface between saturated tearing mode package and BALDUR
c
c
c   Input:        3rd (equilibrium) namelist
c   Common block: / comisl /
c   Cliche:       clislb
c   Defaults:     lislnd(*) = 0     cislnd(*) = 0.0
c
c   lislnd(1)     negative    turns off tearing mode package completely
c                 0           performs diagnosic computation only
c                 1           enhances transport coefficients below
c
c   lislnd(2)     1           enhances \chi_e across islands
c
c   lislnd(3)     1           enhances \chi_i across islands
c
c   lislnd(4)     1           enhances D_H    across islands
c
c   lislnd(5)     1           enhances D_{imp} across islands
c
c   lislnd(9)     .lt. 10     incrementally adjusts island width
c                               with two iterations
c                 .eq. 10     uses zero finder to determine
c                               fully converged island width
c                 .gt. 10     incrementally increases xhalf(jm)
c                               until all islands are too wide
c                               then uses zero finder
c
c
c   cislnd(1)     time to start computing tearing modes [sec] def=0.
c
c   cislnd(2)     time between computaiton of saturated tearing modes
c                       must be .ge. epslon for package to be activated
c
c   cislnd(20)    dqfact = cislnd(20) if cislnd(20) .gt. 1.0
c
c   cislnd(21)    ( if > epslon ) = maximum island half-width allowed
c
c           -----------------------------------------------
c
      subroutine islbal
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clparm.m'
      include 'clislb.m'
c
      dimension  zpresi(mj), znhanc(mj)
c
      data inital /0/
c
c  zpresi(j)  = total pressure in MKS (internal) units on zone centers
c  inital     = 0 until sbrtn island is called to set variables
c
c           -----------------------------------------------
c
c..Test whether to call the tearing mode package or not
c
      if ( timisl .gt. tai  .or.  cislnd(1) .gt. tai
     &  .or.  lislnd(1) .lt. 0  .or.  cislnd(2) .lt. epslon)
     &  go to 30
c
c..reset timisl
c
      itime  = ( tai - cislnd(1) ) / cislnd(2)
      timisl = cislnd(1) + cislnd(2) * ( itime + 1 )
c
c..Compute island widths
c
c  only circular cylinder approximation available at present
c
      ieqtyp = 0
c
c  convert from standard to internal units
c
      do 20 j=1,mzones
        zpresi(j) = totprs(j) * usih
  20  continue
c
c..compute island widths
c
      nxiham = mflxs
      nhamad = mhrms
c
      inital = 1
c
      if ( inital .lt. 99 ) then
c
        write (6,*)
     &    'sbrtn island not available until library routines removed'
c
        return
c
      endif
c
c      call island
c     & ( mzones, xzoni,  ajtori(1,2),    zpresi, avi(1,9,1)
c     & , avi(1,10,1),    xbouni, q,      ieqtyp, kjflx
c     & , mflxs,  mhrms,  equixi, equir0, eqrmom, eqymom
c     & , nmodes, modem,  moden,  qmode,  xmode,  xhalf,  xdelta
c     & , dprime, qb1amp, qb1scl
c     & , qdeltx, xhfmin, qxwall, nxiham, nhamad, nithta
c     & , qpeakc, qshldc, qwexc
c     & , cislnd, lislnd)
c
c..printout from BALDUR code
c
      write (6,121) tai*uiet
 121  format (/'   Evaluated at time ',f10.6,' sec')
c
c..modify diffusion coefficients        rgb 16-oct-87
c
  30  continue
c
      if ( inital .gt. 0  .and.  lislnd(1) .gt. 0 ) then
c
        do 38 j=3,mzones
c
c  need q-values at zone centers since the diffusion coefficients
c    are averaged over a control volume from zone center to zone center
c
          zqp = 0.5 * ( xbouni(j+1) + xbouni(j) )
          zqm = 0.5 * ( xbouni(j) + xbouni(j-1) )
c
          do 36 jm=1,nmodes
c
            if ( xhalf(jm) .lt. xhfmin
     &          .or. xmode(jm) + xhalf(jm) .lt. zqm
     &          .or. xmode(jm) - xhalf(jm) .gt. zqp ) go to 36
c
            zuup = ( zqp - xmode(jm) ) / xhalf(jm)
            zuum = ( zqm - xmode(jm) ) / xhalf(jm)
c
            zmup = min ( max ( zuup, -1. ) , 1. )
            zmum = min ( max ( zuum, -1. ) , 1. )
c
c..transport enhancement factor
c
            znhanc(j) = 1. / ( 1. +
     &  ( zmup * ( zmup**2 / 3. - 1. ) - zmum * ( zmum**2 / 3. - 1. ) )
     &  / ( zuup - zuum ) )
c
c  electron thermal diffusivity
c
            if ( lislnd(2) .gt. 0 ) then
              detes(j) = detes(j) * znhanc(j)
              detepr(j) = detes(j) / ( rhoels(1,j) + epslon )
            endif
c
c..ion thermal diffusivity
c
            if ( lislnd(3) .gt. 0 ) then
              ditis(j) = ditis(j) * znhanc(j)
              ditipr(j) = ditis(j) / ( rhoins(1,j) + epslon )
            endif
c
c..hydrogen diffusion coefficients
c
            if ( lislnd(4) .gt. 0 ) then
              do 32 jh=1,mhyd
                dnhs(jh,j) = dnhs(jh,j) * znhanc(j)
  32          continue
            endif
c
c..impurity diffusion coeffients
c
            if ( lislnd(5) .gt. 0 ) then
              do 34 ji=1,mimp
                dnis(ji,j) = dnis(ji,j) * znhanc(j)
  34          continue
            endif
c
  36      continue
c
  38    continue
c
      endif
c
c
c
c
c
c
      return
      end

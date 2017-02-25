c@clcnst.m constants for the island package, Bateman, PPPL
c rgb 12:00 17-jan-92 renamed common /const/ to common /cislcn/
c
      common /cislcn/     pi,     emu0,   epslon
     & , leqtyp, lprofl, maxisl, relerr, abserr, tolzro
     & , cislnd(32),  lislnd(32)
c
      common /cislgr/
     &   grrmin, grrmax, grymin, grymax
     & , text
      character *80 text(3)
c
c  pi     = $ \pi $
c  emu0   = 4 * pi * 1.e-7
c  epslon = machine epsilon ( 1. + epslon = 1. due to round off )
c  leqtyp,... see documentation above
c
c  text(1) = program version
c  text(2) = type of equilibrium

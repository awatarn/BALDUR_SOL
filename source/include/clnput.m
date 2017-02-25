c@clnput for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cinput/
     &   nmodes, modem,  moden,  xhalf,  b1amp,   b1scal, dprime
     & , deltax, xhfmin, qpeakc, qshldc, qwexc,  xwall,  xistop
     & , nxiprf, qaxis,  btaxis, r0ref
     & , nmomeq, nxiequ
     & , rc0bc,  rcmbc, rsmbc, shiftr, yc0bc, ycmbc, ysmbc, shifty
     & , nxiham, nhamad, nithta
     & , nprskp, momskp
     & , beta,   nbeta, betastp, expjt,  nexpjt, expjtstp
     & , praxis, exppr, nexppr, expprstp, nqaxis, qaxisstp, psetbyj
c
      dimension   rcmbc(kpharm),    rsmbc(kpharm)
     & ,          ycmbc(kpharm),    ysmbc(kpharm)
     & ,          qpeakc(kpmode),   qshldc(kpmode)
     & ,          modem(kpmode),    moden(kpmode),    xhalf(kpmode)
     & ,          qmode(kpmode),    xmode(kpmode),    xdelta(kpmode)
     & ,          dprime(kpmode),   b1amp(kpmode),    b1scal(kpmode)
c
c    This common block is intended only for use in the driver program
c  and its sbrtns.  This block contains all the input data except
c  the first line of namelist, which is in common const.
c    For definitions, see input data
c    Dimensions are also given for qmode(j), xmode(j), xdelta(j), and
c  dprime(j).

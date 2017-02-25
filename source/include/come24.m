c@come24.m
c
      parameter  (jbal=55)
c
      common  /come24/
     1        ppb    , ppc    , ffpc   ,          qvlb   , qvlc   ,
     2        psp    , pp     , ffp    ,                   qvl
c     .............................................................
      dimension    ppb(jbal),ppc(jbal),ffpc(jbal),
     >             qvlb(jbal),qvlc(jbal)
      dimension    psp(jm),pp(jm),ffp(jm),qvl(jm)
c
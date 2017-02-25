c@olycom.m
cl                  c1.1.    basic system parameters
c     version 2c           1.8.73      kvr/mhh        culham
       common/combas/
     +   altime,   cptime,   nledge,   nlend,    nlres,    nonlin,
     +   nout,     nprint,   nread,    nrec,     nresum,   nstep,
     +   stime,    label1,   label2,   label3,   label4,   label5,
     +   label6,   label7,   label8,   ndiary,   nin,      npunch,
     +   nrun,     reread,   versno
       logical     nlend,    nlres
       dimension
     h   label1(12),         label2(12),         label3(12),
     h   label4(12),         label5(12),         label6(12),
     h   label7(12),         label8(12)
c/ module comddp
c--------------------------------------------------------------------
cl                  c1.9.    development and diagnostic parameters
c     version 2c           1.8.73      kvr/mhh        culham
       common/comddp/
     i   maxdum,   mxdump,   nadump,   npdump,   npoint,
     i   nsub,     nvdump,
     l   nlched,   nlhead,   nlomt1,   nlomt2,   nlomt3,   nlrept
       logical
     l   nlched,   nlhead,   nlomt1,   nlomt2,   nlomt3,   nlrept
       dimension
     i   nadump(20),         npdump(20),         nvdump(20),
     l   nlhead(9),          nlomt1(50),         nlomt2(50),
     l   nlomt3(50)
c

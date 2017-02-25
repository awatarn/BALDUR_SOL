c@clsaw.m  .../com/clsaw.m
c  rap 28-aug-02 maximum number of sawtooth times is increased to 320 
c  rgb 16-jun-98 added kswmax = 32 and swton --> swton(32)
c  cliche clsaw: for sawtooth input variables
c               and sawtooth variables needed in the rest of BALDUR
c
       integer, parameter :: kswmax = 320
c
       common /comsw1/
     &  lsawth(32), csawth(32), swton(kswmax)
     & , swtoff,swperd, swqmin, swxmin
c
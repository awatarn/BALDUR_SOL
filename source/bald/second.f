c@second   .../baldur/code/bald/second.f
c
c  This routine makes system calls to determine
c  the elapsed time in seconds.
c  Glenn Bateman, Lehigh University, bateman@pppl.gov
c
cap   changed subroitine name (second -> timeInSeconds)
      subroutine timeInSeconds ( ptime )
c
      real ptime, mclock
c
c..IBM workstation system call mclock
c  yeilds elapsed time in milli-seconds
c
      ptime = 0.01 * mclock()
c
      return
      end

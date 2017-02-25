c-----------------------------------------------
 
      subroutine cptimi(ztime)
      implicit none
 
      real ztime
c
c  DMC 5 Sep 1991
c
c    this routine is a more portable FORTRAN version of the old macro
c    CPTIMI.MAR.  It returns the amount of cpu time used so far by the
c    current process, in a floating point word, in units of HUNDREDTHS
c    of seconds *** NOTE FUNNY UNITS which are used for historical
c    reasons
c
c    usage:
c
c       CALL CPTIMI(ZTIME1)
c
c         execute code-to-be-timed
c
c       CALL CPTIMI(ZTIME2)
c
c       ZTIME = 0.01*(ZTIME2-ZTIME1)
c       TYPE 1001,ZTIME
c 1001	FORMAT('  CPU time used = ',1pe11.4,' seconds.')
c
c-------------------------------------------------
c
#ifdef __VMS
      include '($JPIDEF)'
      integer iarg,istat,iret,lib$getjpi
#endif
c
c  executable code:
c
#ifdef __VMS
      iarg = JPI$_CPUTIM
      istat = LIB$GETJPI(iarg,,,iret)
 
      ztime = float(iret)
#endif
 
c  unix version -- return value in funny units for consistency!
 
#ifdef __UNIX
 
#if __HP || __IBM || __RS6000 || __NAGWARE
      integer cpclock
#elif __CRAY
      real second
#elif __ABS
      real etime_
      real tarray(2)
#else
      real etime
      real tarray(2)
#endif
 
C-----------------------
 
#if __HP || __IBM || __RS6000 || __NAGWARE
      ztime = 1.0e-4*cpclock()
#elif __CRAY
      ztime = 100.0*second()
#elif __ABS
      ztime = 100.0*etime_(tarray)
#else
      ztime = 100.0*etime(tarray)
#endif
 
#endif  /* __UNIX */
 
      return
      end

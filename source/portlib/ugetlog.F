      subroutine ugetlog(zbuff)
      character*(*) zbuff               ! for username string, output
c
c  return the current username
c
#ifdef __UNIX
 
      call sget_env('LOGNAME',zbuff)
      if(zbuff.eq.' ') then
         call sget_env('USER',zbuff)
         if(zbuff.eq.' ') then
            zbuff='UNKNOWN'
         endif
      endif
#endif
      return
      end

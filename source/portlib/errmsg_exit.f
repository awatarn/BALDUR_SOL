      subroutine errmsg_exit(msg)
      character*(*) msg
C
C  write message and exit
C
      write(6,*) ' '
      write(6,*) msg
      write(6,*) ' '
C
      call bad_exit
      return
      end

      subroutine bad_exit
#ifdef __NAGWARE
      use f90_unix_proc
#endif
C
      write(6,*) ' %bad_exit:  generic f77 error exit call.'
C
      call err_end
      close(6)
C
#ifdef __VMS
      call lib$stop( %val(2) )
#endif
#ifdef __UNIX
      call c_abort
#endif
C
      STOP                              ! this line not reached.
      END

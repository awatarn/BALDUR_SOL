      program test
c
      character cfile*32, cline*32
c
      cfile = 'prof-exp-TI.dat'
c
      read (*,*) cline
c
      open ( 4, file=cfile, err=90 )
c
      write (*,*) 'opened file ', cfile
      write (*,*) 'cline = ', cline
c
      stop
c
   90 write(*,*) 'unable to open file ', cfile
c
      stop
      end

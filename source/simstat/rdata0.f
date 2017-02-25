c@rdata0.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    Sbrtn to read data from xmgr input file for program simstat.
c
c
      subroutine rdata0 ( radius, data, derror, ndata
     &  , ninput )
c
c
      implicit none
c
      character cfile*32, ctemp*32
      character cline*132
c
c  cfile    = name of the file containing experimental data
c  ctemp    = temporary character string for file name
c  cline    = temporary character string for line of input
c
c.. data points
c
c
      real       radius(*), data(*), derror(*)
c
c  radius(j)   = major radius of data points
c  data(j)     = data points
c  derror(j)   = standard deviations error bars for data
c
      integer ninput, ndata, j
c
c  ndata       = number of data points
c
      real  convradius, convdata
c
c  convradius = conversion factor for radius
c  convdata   = conversion factor for data and derror
c
c
c..determine the name of the file containing experimental data
c
      cfile = '                                '
c
      write (*,*) 'Name of the file containing data:'
      read  (*,*) cfile
c
      if ( cfile(1:1) .eq. ' ' ) then
        write (*,*) 'Abort in sbrtn rdata0'
        write (*,*) ' no name given for cfile = ',cfile
        stop
      endif
c     
c..conversion factors
c     
      convradius = 1.0
      convdata = 1.0
      write (*,*) 'Type conversion factor for radius'
      read (*,*) convradius
      write (*,*) 'Type conversion factor for data'
      read (*,*) convdata
c
      write (*,*) 'Conversion factors are:  convradius = ',convradius
     &  ,'  convdata = ',convdata
c
c
c..Open file and skip headers (lines starting with #)
c
      write (*,*) ' Reading data from file ',cfile
      open ( ninput, file=cfile, status = 'old', err=90)
c
   10 continue
c
      read ( ninput, 100, err=92, end=92 ) cline
c
      if ( cline(1:1) .eq. '#' ) go to 10
c
      backspace ninput
c
c.. read data lines
c
      j = 0
c
   12 continue
c
      j = j + 1
      radius(j) = 0.0
      data(j)   = 0.0
      derror(j) = 0.0
c
      read ( ninput, *, err=14, end=14 ) radius(j), data(j)
      ndata = j
c
c.. use conversion factors
c
      radius(j) = convradius * radius(j)
      data(j)   = convdata * data(j)
c
      go to 12
c
   14 continue
c
c
      close ( ninput )
c
      write (7,*)
      write (7,*) 'Data from file ',cfile
      write (7,*)
c
c..array outputs
c
      write (7,140) ndata
  140 format (t2,'ndata  = ',i4,
     &  /t4,'radius',t17,'data')
c
      do j=1,ndata
        write (7,141) radius(j), data(j)
      enddo
  141 format (1p5e13.4)
c
      return
c
c..error
c
  90  continue
      write (*,*) 'Unable to open file cfile = ', cfile
c
      stop
c
  92  continue
      write (*,*) 'Reached end of experimental data file with no data'
c
      stop
c
 100  format (a132)
 101  format (a32)
c
      end

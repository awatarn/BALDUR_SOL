c@baldata.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    Sbrtn to read experimental data for program simstat.
c
c
      subroutine baldata ( rmsim, tesim, tisim, ensim, nsim
     &  , ninput )
c
c
      implicit none
c
      character csimfile*32, ctemp*32
      character cline*132
c
c  csimfile = name of the file containing BALDUR simulation output
c  ctemp    = temporary character string for file name
c  cline    = temporary character string for line of input
c
c
c..BALDUR simulation profiles
c
      real    rmsim(*), tesim(*), tisim(*), ensim(*)
      real    temp(501)
      integer nsim
c
c  rmsim(j)   = major radii from BALDUR simulation for all profiles
c  tesim(j)   = T_e(R) (keV)
c  tisim(j)   = T_i(R) (keV)
c  ensim(j)   = n_e(R) (m^{-3})
c  nsim       = number of elements in simulation profile arrays
c
c  ninput = unit number for input file
c
      real repsilon
c
c  repsilon = smallest major radius to be accepted
c
c
      integer ninput, input, ndatype, j
c
c
      repsilon = 0.001
c
      input = ninput
      if ( input .lt. 1  .or.  input .gt. 99 ) input = 4
c
c
c..determine the type of input data file (ndatype)
c
      ndatype = 0
      write (*,*) 'Indicate the type of input data file ndatype ='
      write (*,*) '0 for separate xmgr input files'
      write (*,*) '1 jobxlpt file'
c
      read (*,*) ndatype
c
c
c..read separate xmgr input files if ndatype = 0
c
      if ( ndatype .eq. 0 ) then
c
        call rdata0 ( rmsim, tisim, temp, nsim, input )
c
        call rdata0 ( rmsim, tesim, temp, nsim, input )
c
        call rdata0 ( rmsim, ensim, temp, nsim, input )
c
        return
c
      endif
c
c
c..determine the name of the file containing BALDUR simulation output
c
      ctemp = '                                '
      csimfile = 'jobxlpt'
c
      write (*,*)
      write (*,*) 'Name of the file containing BALDUR simulation output'
      read  (*,101) ctemp
c
      if ( ctemp(1:1) .ne. ' ' ) csimfile = ctemp
c
      write (*,*) ' Reading simulation data from file ',csimfile
      open (input,file=csimfile,err=91)
c
c
c..Read simulation profiles
c  continue to the end of the file to read the last profile
c
  30  continue
      read (input,100,err=33,end=33) cline
      if ( cline(4:24) .eq. 'rmajor(m)    ne(m^-3)' ) go to 31
      go to 30
c
  31  continue
      nsim = 0
c
  32  continue
      j = nsim + 1
      read (input,132,err=30,end=33) 
     &  rmsim(j), ensim(j), tesim(j), tisim(j)
 132  format (5(2x,1pe11.4))
      if ( rmsim(j) .lt. repsilon ) go to 32
      nsim = j
      go to 32
c
  33  continue
      if ( nsim .lt. 1 ) go to 93
c
c
c..normal return
c
      write (*,*) ' nsim = ',nsim
c
      close ( input )
c
      return
c
c..error
c
c
  91  continue
      write (*,*) 'Unable to open file csimfile = ', csimfile
c
      stop
c
  93  continue
      write (*,*) 'Reached end of simulation output file with no data'
c
      stop
c
 100  format (a132)
 101  format (a32)
c
      end

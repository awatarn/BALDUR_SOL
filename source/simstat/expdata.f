c@expdata.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    Sbrtn to read experimental data for program simstat.
c
c
      subroutine expdata ( rmteexp, teexp, teerr, nteexp
     &  , rmtiexp, tiexp, tierr, ntiexp
     &  , rmenexp, enexp, enerr, nenexp
     &  , ninput )
c
c
      implicit none
c
      character cexpfile*32, ctemp*32
      character cline*132
c
c  cexpfile = name of the file containing experimental data
c  ctemp    = temporary character string for file name
c  cline    = temporary character string for line of input
c
c..experimental data points
c
c    major radii (m), temperatures (eV) or density (m^{-3}),
c    and error bars
c
      real       rmteexp(*), teexp(*), teerr(*)
     &         , rmtiexp(*), tiexp(*), tierr(*)
     &         , rmenexp(*), enexp(*), enerr(*)
c
c  rmteexp(j) = major radius of experimental T_e data points (m)
c  teexp(j)   = T_e electron temperature data points (keV)
c  teerr(j)   = standard deviations error bars for T_e (keV) ...
c
      integer  nteexp, ntiexp, nenexp, ninput, ndatype
c
c  nteexp = number of T_e experimental data points,
c  ntiexp = number of T_i experimental data points,
c  nenexp = number of n_e experimental data points,
c
c  ninput = unit number for input file
c  ndatype = type of input file
c
      real repsilon, convdens, convtemp
c
c  repsilon = smallest major radius to be accepted
c  convdens = conversion factor for density
c  convtemp = conversion factor for temperatures
c
c
      integer input, iptype, iexprun, j, itemp
c
c  input   input unit number
c  iptype  = 1 for ne, = 2 for Te, = 3 for Ti
c  iexprun = experimental run number
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
      write (*,*) '1 for lists of profiles from Kinsey'
      write (*,*) '2 TRANSP data  prepared using '
     &  ,'[BATEMAN.TRANSP.TOOLS]TEST.EXE'
      write (*,*) '3 from SNAP run using [BATEMAN.BPLOT]COMPARE.COM'
c
      read (*,*) ndatype
c
      if ( ndatype .lt. 0  .or.  ndatype .gt. 3 ) then
        write (*,*) 'natype out of range.  '
     &    ,'Type a number between 1 and 3'
        read (*,*) ndatype
        if ( ndatype .lt. 1  .or.  ndatype .gt. 3 ) then
          write (*,*) 'Abort:  unable to use ndatype = ',ndatype
          stop
        endif
      endif
c
c
c..read separate xmgr input files if ndatype = 0
c
      if ( ndatype .eq. 0 ) then
c
        call rdata0 ( rmtiexp, tiexp, tierr, ntiexp, input )
c
        call rdata0 ( rmteexp, teexp, teerr, nteexp, input )
c
        call rdata0 ( rmenexp, enexp, enerr, nenexp, input )
c
        return
c
      endif
c
c..determine the name of the file containing experimental data
c
      ctemp = '                                '
      cexpfile = 'expdata'
c
      write (*,*) 'Name of the file containing experimental data:'
      read  (*,101) ctemp
c
      if ( ctemp(1:1) .ne. ' ' ) cexpfile = ctemp
c
      write (*,*) ' Reading experimental data from file ',cexpfile
      open (input,file=cexpfile,err=90)
c     
c..conversion factors
c     
      convdens = 1.0
      convtemp = 1.0
      write (*,*) 'Type conversion factor for densities'
      read (*,*) convdens
      write (*,*) 'Type conversion factor for temperatures'
      read (*,*) convtemp
c
      write (*,*) 'Conversion factors are:  convdens = ',convdens
     &  ,'  convtemp = ',convtemp
c
c
      if ( ndatype .eq. 1 ) go to 10
      if ( ndatype .eq. 2 ) go to 20
      if ( ndatype .eq. 3 ) go to 30
c
c
c..list from Kinsey
c
  10  continue
c
        nenexp = 0
        nteexp = 0
        ntiexp = 0
c
  11    continue
        read (input,100,err=92,end=80) cline
c     
cbate        write (*,*) cline
cbate        itemp = index(cline,'ne')
cbate        write (*,*) 'index(cline,ne) = ',itemp
c
        if ( index(cline,'ne') .gt. 0 ) then
          iptype = 1
          go to 11
        elseif ( index(cline,'Te') .gt. 0 ) then
          iptype = 2
          go to 11
        elseif ( index(cline,'Ti') .gt. 0 ) then
          iptype = 3
          go to 11
        elseif ( iptype .eq. 1 ) then
          backspace input
          j = nenexp + 1
          read (input,*,err=92,end=92)
     &      rmenexp(j), enexp(j), enerr(j)
          nenexp = j
          go to 11
        elseif ( iptype .eq. 2 ) then
          backspace input
          j = nteexp + 1
          read (input,*,err=92,end=92)
     &      rmteexp(j), teexp(j), teerr(j)
          nteexp = j
          go to 11
        elseif ( iptype .eq. 3 ) then
          backspace input
          j = ntiexp + 1
          read (input,*,err=92,end=92)
     &      rmtiexp(j), tiexp(j), tierr(j)
          ntiexp = j
          go to 11
        endif
c
        go to 80
c
c
c..for TRANSP data, determine beginning of profiles
c  file prepared using [BATEMAN.TRANSP.TOOLS]TEST.EXE
c
  20  continue
c
      nteexp = 0
      read (input,100,err=92,end=92) cline
      if ( cline(4:26) .eq. 'rmajor       n_e_TRANSP' ) go to 21
      if ( cline(13:14) .eq. 'NE' ) go to 30
      go to 20
c
  21  continue
      j = nteexp + 1
      read (input,*,err=22,end=22) 
     &  rmteexp(j), enexp(j), teexp(j), tiexp(j)
      if ( rmteexp(j) .lt. repsilon ) go to 22
      nteexp = j
      go to 21
c
  22  continue
      if ( nteexp .lt. 1 ) go to 92
      nenexp = nteexp
      ntiexp = nteexp
c
      do j=1,nteexp
        rmenexp(j) = rmteexp(j)
        rmtiexp(j) = rmteexp(j)
      enddo
c
      go to 80
c
c..file prepared from SNAP run using [BATEMAN.BPLOT]COMPARE.COM
c
  30  continue
c
      nenexp = 0
      nteexp = 0
      ntiexp = 0
c
  31  continue
      j = nenexp + 1
      read (input,*,err=32,end=32) 
     &  rmenexp(j), enexp(j), enerr(j)
      if ( rmenexp(j) .lt. repsilon ) go to 32
      nenexp = j
      go to 31
c
  32  continue
      backspace input
      nteexp = 0
  33  continue
      read (input,100,err=92,end=92) cline
      if ( cline(13:14) .eq. 'TE' ) go to 34
      go to 33
c
  34  continue
      j = nteexp + 1
      read (input,*,err=35,end=35) 
     &  rmteexp(j), teexp(j), teerr(j)
      if ( rmteexp(j) .lt. repsilon ) go to 35
      nteexp = j
      go to 34
c
  35  continue
      backspace input
      ntiexp = 0
  36  continue
      read (input,100,err=92,end=92) cline
      if ( cline(13:14) .eq. 'TI' ) go to 37
      go to 36
c
  37  continue
      j = ntiexp + 1
      read (input,*,err=38,end=38) 
     &  rmtiexp(j), tiexp(j), tierr(j)
      if ( rmtiexp(j) .lt. repsilon ) go to 38
      ntiexp = j
      go to 37
c
  38  continue
      go to 80
c
c
c
c
c..normal return
c
  80  continue
c     
c..convert units as needed
c
      if ( convdens .gt. 1.e-34 ) then
        do j=1,nenexp
          enexp(j) = convdens * enexp(j)
          enerr(j) = convdens * enerr(j)
        enddo
      endif
c
      if ( convtemp .gt. 1.e-34 ) then
        do j=1,nteexp
          teexp(j) = convtemp * teexp(j)
          teerr(j) = convtemp * teerr(j)
        enddo
        do j=1,ntiexp
          tiexp(j) = convtemp * tiexp(j)
          tierr(j) = convtemp * tierr(j)
        enddo
      endif
c
c..diagnostic printout
c
      write (*,*) ' nteexp = ',nteexp,'  ntiexp = ',ntiexp
     &  ,'  nenexp = ',nenexp
c
      close ( input )
c
      return
c
c..error
c
  90  continue
      write (*,*) 'Unable to open file cexpfile = ', cexpfile
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

c/ 21:00 21-feb-96 /011040/.../baldur/code/bald/dio.f  Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c  To obtain this file, type           (use appropriate date for yymmdd)
c cfs get /11040/byymmdd/baldyy/wcode.tar
c tar xf wcode.tar  (this creates directory .../code and subdirectories)
c cd code/bald      (this changes to subdirectory .../code/bald)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c  file dio contains input output routines for baldur.  Subrtns:
c
c  stripx  - to strip ! off input data sets
c  data   - input namelist data (not including equilibrium data)
c  redata - reread input data at time REREAD or cfutz(82)
c  output - control output (moved to file DOUTPUT)
c  sprint - generate short printout
c  aprint - alpha particle printout
c  hprint - print out beam data
c  iprint - print out injector data
c  grafix - output to intermediate file for04 for05 for graphics
c  fprint - fusion output
c  gprint - neutral gas printout
c  ncprnt - prints results of impurity charge state transport
c  nclwpr - local radiation printout
c  nciwpr - integrated radiation printout
c  mprint - long printout, many pages
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@stripx   /baldur/code/bald/stripx.f
c  rgb 22-feb-2005 do j=1,kc --> do j=1,kc-1
c  rgb 25-aug-96 clean the data set to conform to namelist requirements
c  rgb 16-aug-96 insert commas as needed in namelist data
c  rgb 06-oct-94 changed line length from 80 to 132 characters
c
      subroutine stripx (nin,ntrans,nout)
c
c  This sbrtn reads input data file logical unit number nin
c  one line at a time up to 80 characters long,
c  prints out the input verbatum to output file logical unit number nout
c  then searches for the first appearance of the character ! on each line
c  and prints out everything to the left of ! to output unit number ntrans
c
c  In this way, everything to the right of ! is treated as a comment
c
      parameter (kc = 132)
c
      character line * 132
c
      logical lcomma, lequil, lspace, lnlist
c
c..lcomma  a comma has been encountered
c  lequil  an "=" has been encountered
c  lspace  the last character was a space
c  lnlist  this line is in a namelist
c
c..iline = line number
c
      iline = 0
      lnlist = .false.
c
  10  read (nin,100,end=20) line
 100  format (a132)
c
      iline = iline + 1
c
c..find number of characters before spaces
c
      ilength = 0
      do j=1,kc
        if ( line(j:j) .ne. ' ' ) ilength = j
      enddo
c
c..echo input data on output file
c
      if ( ilength .gt. 0 ) then
        write (nout,110) line(1:ilength)
      else
        write (nout,*)
      endif
c
c..ixlen = number of nontrivial characters before "!"
c  ieq   = position of "=" in line
c
      ixlen = 0
      ieq   = 1
      do j=1,kc-1
        if ( line(j:j) .eq. '!' ) go to 14
        if ( line(j:j) .ne. ' ' ) ixlen = j
        if ( line(j:j) .eq. '=' ) ieq = j
      enddo
  14  continue
c
c..clean the data set to conform to namelist requirements
c
      i1 = 2
      i2 = min( max(ixlen+1,i1), kc )
c
      if ( index(line(1:ixlen+1),'$end') .gt. 1
     &   .or. index(line(1:ixlen+1),'/') .gt. 1 ) then
        lnlist = .false.
        go to 19
      elseif ( ( index(line(1:ixlen+1),'$') .gt. 1 )
     &    .or. ( index(line(1:ixlen+1),'&') .gt. 1 )
     &    ) then
        lnlist = .true.
        go to 19
      endif
c
c
      if ( lnlist 
     &  .and. ( index(line(1:ixlen+1),'$') .lt. 1 )
     &  .and. ( index(line(1:ixlen+1),'&') .lt. 1 )
     &  ) then
c
      lequil = .false.
      lcomma = .false.
      lspace = .true.
c
      do 18 j=i2,i1,-1
c
        if ( line(j-1:j-1) .eq. ' ' ) then
          lspace = .true.
          go to 18
        endif
c
        if ( line(j-1:j-1) .eq. '=' ) then
          lequil = .true.
          go to 18
        endif
c
        if ( line(j-1:j-1) .eq. ',' ) then
          lcomma = .true.
          go to 18
        endif
c
        if ( lequil .or. lcomma ) then
          lequil = .false.
          lcomma = .false.
          lspace = .false.
          go to 18
        endif
c
        if ( lspace ) then
          line(j:j) = ','
          lspace = .false.
          if ( j .gt. ixlen ) ixlen = j
          go to 18
        endif
c
  18  continue
c
      endif
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
c..echo input data up to "!" on temporary file
c
  19  continue
c
      if ( ixlen .gt. 0 ) then
        write (ntrans,110) line(1:ixlen)
      else
cbate        write (ntrans,*)
      endif
c
 110  format (a)
c
c
      go to 10
c
  20  continue
      rewind ntrans
c
      return
      end

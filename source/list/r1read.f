c@r1read.f   1 December 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine r1read ( crfile, ktunit, knunit, ktmax, crname
     & , ptime, ktime, parray, key, crlabel, crunits, kerr )
c
c   This routine reads 1-D ASCII RPLOT files
c
c  crfile    = base of the RPLOT file names                   (input)
c    Old style:
c      tf<crfile> = file containing labels, units, and array names
c      nf<crfile> = file containing 1-D arrays as a fn of time
c    New style:
c      <crfile>TF.PLN = index file
c      <crfile>YF.PLN = file containing 1-D arrays as a fn of time
c  ktunit    = allowed input unit number for tf<crname>       (input)
c  knunit    = allowed input unit number for nf<crname>       (input)
c
c  ktmax     = maximum number of elements allowed in time     (input)
c  crname    = name of RPLOT variable (up to 10 characters)   (input)
c
c  ptime(jt) = 1-D time array from RPLOT file                 (output)
c  ktime     = number of elements in ptime from RPLOT file    (output)
c  parray(jt) = 1-D array from RPLOT file                     (output)
c
c  key       = element in RPLOT file for crname               (output)
c      key = 0 if crname not found in tf<crname> file
c  crlabel   = labels in tf<crname> file (20 char each)       (output)
c  crunits   = units  in tf<crname> file (10 char each)       (output)
c
c  kerr      = error indicator                          (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the RPLOT variable named <crname> is not there
c            = 2 if the RPLOT file cannot be opened or read
c            > 2 if input is not consistent
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      real parray(*), ptime(*)
c
      integer ktunit, knunit, ktmax, ktime
     &  , key, kerr
c
      real*8 ztemp(1000)
c
      integer   iprint, ierr, iarray(32)
     &  , ichmax, ic1, ic2, islash
     &  , ir1var, ikeys
     &  , j, jr, jt
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c  ic1     = first non-blank character
c  ic2     = last contiguous non-blank character
c  islash  = end of subdirectory name (last '/' character in crfile)
c
c  ir1var  = number of 1-D variables in 1-D RPLOT files
c  ikeys   = number of variable names that match variabel name crname
c  ierr    = input value of kerr
c
      character  cline*80, crfile*(*)
     &  , crname*(*), crlabel*(*), crunits*(*)
     &  , ctfile*64, cnfile*64, ctemp*20
c
c
c..test input
c
      if ( ktunit .lt. 1  .or.  ktunit .gt. 99 ) then
        write (*,*)
        write (*,*) ' ktunit = ',ktunit,'  out of range in sbrtn r1read'
        kerr = 9
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( knunit .lt. 1  .or.  knunit .gt. 99 ) then
        write (*,*)
        write (*,*) ' knunit = ',knunit,'  out of range in sbrtn r1read'
        kerr = 9
        if ( ierr .lt. 0 ) stop
        return
      endif
c
      if ( ktunit .eq. knunit ) then
        write (*,*)
        write (*,*) ' ktunit = knunit = ',knunit,'  in sbrtn r1read'
        kerr = 9
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..set defaults
c
      ierr    = kerr
      ktime   = 0
      ir1var  = 0
c
      kerr    = 0
c
      iprint  = 0
c
c..clear arrays
c
      do j=1, ktmax
       ptime(j)  = 0.0
       parray(j) = 0.0
      enddo 
c
c
c..determine filename crfile(ic1:ic2)
c
      ic1 = 0
      ic2 = 0
      islash = 0
c
      ichmax = len ( crfile )
      if ( ichmax .lt. 1 ) go to 92
c
      do j=1,ichmax
        if ( crfile(j:j) .ne. ' ' ) then
          if ( crfile(j:j) .eq. '/' ) islash = j
          if ( ic1 .lt. 1 ) then
            ic1 = j
          else
            ic2 = j
          endif
        else
          if ( ic1 .gt. 0 ) go to 12
        endif
      enddo
c
  12  continue
c
c
c..open the input ASCII RPLOT files
c
c  New style of file names
c
      ctfile = ' '
      cnfile = ' '
      if ( ic1 .eq. 0 ) then
        ctfile = 'TF.PLN'
        cnfile = 'YF.PLN'
      else
        ctfile = crfile(ic1:ic2) // 'TF.PLN'
        cnfile = crfile(ic1:ic2) // 'YF.PLN'
      endif
c
      open (unit=ktunit,file=ctfile,status='old',err=14)
      open (unit=knunit,file=cnfile,status='old',err=14)
      go to 18
c
c  Old style of file names
c
  14  continue
c
      ctfile = ' '
      cnfile = ' '
      if ( ic1 .eq. 0 ) then
        ctfile = 'tf'
        cnfile = 'nf'
      else if ( islash .gt. 0 ) then
        if ( islash .lt. ic2 ) then
          ctfile = crfile(ic1:islash) // 'tf'//crfile(islash+1:ic2)
          cnfile = crfile(ic1:islash) // 'nf'//crfile(islash+1:ic2)
        else
          ctfile = crfile(ic1:islash) // 'tf'
          cnfile = crfile(ic1:islash) // 'nf'
        endif
      else
        ctfile = 'tf' // crfile(ic1:ic2)
        cnfile = 'nf' // crfile(ic1:ic2)
      endif
c
      open (unit=ktunit,file=ctfile,status='old',err=92)
      open (unit=knunit,file=cnfile,status='old',err=92)
c
  18  continue
c
c..Read tf<crfile> file and set up keys
c
      ctemp = ' '
      read (ktunit,102,err=92,end=92) (iarray(j),j=1,6), ctemp
 102  format (i4,i6,4i5,a)
c
      ir1var = iarray(4)
      if ( ir1var .lt. 0 ) then
        write (*,*)
        write (*,*) ' Error: number of 1-D variables = ',ir1var,
     &    '  in sbrtn r1read'
        kerr = 4
c
        close ( ktunit )
        close ( knunit )
c
        return
      endif
c
      read (ktunit,*)
      read (ktunit,*)
c
      ikeys = 0
      do j=1,ir1var
        cline = ' '
        read (ktunit,100) cline
 100    format (a)
c
          if ( crname(1:5) .eq. cline(31:35) ) then
            key = j
            ikeys = ikeys + 1
            crlabel = cline(1:20)
            crunits = cline(21:29)
            !exit
          endif
c
      enddo
c
      if ( ikeys .lt. 1 ) then
        write (*,*)
        write (*,*) ' ikeys = 0 in sbrtn r1read'
        write (*,*) ' There is no match to  RPLOT variable name'
c
        write (*,*)
        write (*,*) crname(1:5),'= crname(1:5)'
        write (*,*) ctfile,' = ctfile'
        write (*,*) cnfile,' = cnfile'
        write (*,*) cline(31:35),'= cline(31:35)'
        write (*,*) key,' = key'
        write (*,*) ir1var,' = ir1var'
c
        kerr = 1
c
        close ( ktunit )
        close ( knunit )
c
        return
      endif
c
c
c..read nf<crfile> file and fill time-dependent arrays
c
      do jt=1,ktmax 
c
        read (knunit,110,end=30,err=92) (ztemp(jr),jr=1,ir1var)
!         read (knunit,110,end=30) (ztemp(jr),jr=1,ir1var)
  110   format (1x,(1p6e12.5))
c
        ktime = jt
        ptime(jt) = ztemp(1)
c
        parray( jt ) = ztemp( 1 + key )
c
      enddo
c
c
c..normal end of routine
c
c
  30  continue
c
      close ( ktunit )
      close ( knunit )
c
      return
c
c
c..error conditions
c
c
  92  continue
c
      write (*,*)
     & 'error reading line in RPLOT file in sbrtn r1read'
      write (*,*) 'ctfile = ',ctfile
      write (*,*) 'at timestep',ktime
c
      if ( ktime .gt. 0 ) then
        kerr = 0
        write (*,*) 'continuing with kerr = 0'
      else
        kerr = 2
      endif
c
      close ( ktunit )
      close ( knunit )
c
      if ( ierr .lt. 0 ) stop
c
      return
c
  94  continue
c
      write (*,*)
     & 'error reading line in RPLOT file in sbrtn r1read'
      write (*,*) 'cnfile = ',cnfile
      write (*,*) 'at timestep ',ktime
c
      kerr = 2
      if ( ierr .lt. 0 ) stop
c
      close ( ktunit )
      close ( knunit )
c
      return
c
c
      end

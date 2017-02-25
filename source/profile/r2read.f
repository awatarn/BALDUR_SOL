c@r2read.f   1 December 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine r2read ( crfile, kmunit, krdim, ktdim, cname
     & , ptime, ktime, kradius, parray, kerr )
c
c   This routine reads 2-D ASCII RPLOT file
c
c  crfile    = base of the RPLOT file names                   (input)
c    mf<crname> = file containing 2-D array
c  kmunit    = allowed input unit number for tf<crname>       (input)
c
c  krdim     = maximum number of elements allowed in radius   (input)
c  ktdim     = maximum number of elements allowed in time     (input)
c  cname     = name of array  (up to 10 characters)           (input)
c
c  ptime(jt) = 1-D time array from RPLOT file                 (output)
c  ktime    = number of elements in ptime from RPLOT file     (output)
c  kradius   = number of radial elements from RPLOT file      (output)
c  parray(jrt) = 2-D array from RPLOT file                    (output)
c    stored as parray( 1...ktime ) for jr=1
c              parray( ktmax + 1...ktime ) for jr=2 ...
c              parray( (jr-1)*ktmax + 1...ktime for jr
c
c  kerr      = error indicator                          (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the RPLOT file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if kru > krdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      real parray(*), ptime(*), zu(6)
c
      integer  kmunit, krdim, kradius, ktdim, ktime, kufile, kerr
c
      real ztemp(1000)
c
      integer   iprint, ierr, iline, itemp, index
     &  , ichmax, ic1, ic2, islash
     &  , j
c
c  iprint  = integer to control amount of diagnostic output
c  ichmax  = maximum allowed length of file and directory names
c
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
      character  crfile*64, cmfile*64, cline*72, cname*(*), ctemp*10
c
c
c..test input
c
      if ( kmunit .lt. 1  .or.  kmunit .gt. 99 ) then
        write (*,*)
        write (*,*) ' kmunit = ',kmunit,'  out of range in sbrtn r1read'
        kerr = 9
        if ( ierr .lt. 0 ) stop
        go to 98
      endif
c
c
c..set defaults
c
      ierr    = kerr
      iline   = 0
      ktime   = 0
      kradius = 0
      index   = 0
c
      kerr    = 0
c
      iprint  = 0
c
c..clear arrays
c
      do j=1,ktdim
        ptime(j)    = 0.0
      enddo
c
      do j=1,ktdim*krdim
        parray(j) = 0.0
      enddo
c
      write (*,*) 'cname = ',cname
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
c..open the input ASCII RPLOT file
c
c  New style of file names
c
      cmfile = ' '
      if ( ic1 .eq. 0 ) then
        cmfile = 'XF.PLN'
      else
        cmfile = crfile(ic1:ic2) // 'XF.PLN'
      endif
c
      open (unit=kmunit,file=cmfile,status='old',err=14)
      go to 18
c
c  Old style of file names
c
  14  continue
c
c
      cmfile = ' '
      if ( ic1 .eq. 0 ) then
        cmfile = 'mf'
      else if ( islash .gt. 0 ) then
        if ( islash .lt. ic2 ) then
          cmfile = crfile(ic1:islash) // 'mf'//crfile(islash+1:ic2)
        else
          cmfile = crfile(ic1:islash) // 'mf'
        endif
      else
        cmfile = 'mf' // crfile(ic1:ic2)
      endif
c
      open (unit=kmunit,file=cmfile,status='old',err=92)
c
  18  continue
c
c..read next set of records
c
  10  continue
c
c..read header line
c
      itemp = 0
      ctemp = ' '
      read (kmunit,102,err=92,end=94) itemp, ctemp
 102  format (i4,1x,a)
      iline = iline + 1
c
c
      if ( itemp .lt. 1 ) then
        write (*,*) 'error -- number of elements < 1 for',ctemp
        write (*,*) 'on line ',iline,'  of RPLOT file'
        kerr = 1
        go to 98
      endif
c
c..is the next element part of the time array ?
c
      if ( ctemp(1:4).eq.'time' .or. ctemp(1:4).eq.'TIME' ) then
        ktime = ktime + 1
        read (kmunit,*,err=92,end=94) ptime(ktime)
        iline = iline + 1
        go to 10
c
c..is the next element part of the array to be read ?
c  if so, there are six real numbers per line
c
      elseif ( ctemp .eq. cname ) then
        read (kmunit,110,err=92,end=94) (parray(index+j),j=1,itemp)
  110   format(1x,(1p6e12.5))
       kradius = itemp
        index = index + itemp
        iline = iline + ( itemp + 5 ) / 6
        go to 10
c
c..if not, then skip lines
c
      else
        itemp = ( itemp + 5 ) / 6
        do j=1,itemp
          read (kmunit,*,err=92,end=94)
          iline = iline + 1
        enddo
        go to 10
c
      endif
c
c
c
c..end of routine
c
c
      go to 98
c
c
c..error conditions
c
c
  92  continue
c
      write (*,*)
     & 'error reading line in RPLOT file in sbrtn r2read'
      write (*,*)
      write (*,*) 'on line ',iline
c
c..temporary printout
c
      write (*,*) 'cname = ',cname
      write (*,*) 'ctemp = ',ctemp
      write (*,*) 'itemp = ',itemp
      write (*,*) 'iline = ',iline
      kerr = 1
      if ( ierr .lt. 0 ) stop
      go to 98
c
  94  continue
c
      write (*,*)
     & 'reached the end of RPLOT file in sbrtn r2read'
      write (*,*) 'after line ',iline
c
      kerr = 0
      if ( ktime .lt. 1  .or.  kradius .lt. 1 ) then
         write (*,*) 'set kerr = 2 because ktime or kradius .lt. 1'
         write (*,*) 'ktime = ',ktime
         write (*,*) 'kradius = ',kradius
         kerr = 2
      else
         go to 98
      endif
c
      if ( ierr .lt. 0 ) stop
      go to 98
c
c..close file and return
c
  98  continue
      close ( kmunit )
      return
c
      end

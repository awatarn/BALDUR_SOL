c@br2read.f   04 Dec 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine br2read ( cvarname, kunit, kxdim, kydim
     & , pxu, kxu, cxname, pyu, kyu, cyname, parray, kerr )
c
c   This routine reads concatinated 2-D ASCII U-file (Boucher files)
c  or 2-D ASCII RPLOT file.
c
c  cvarname  = name of the independent variable parray        (input)
ccc  crname    = name of radius array                           (input)
c  kunit     = input unit number to read file (default 70)    (input)
c  kxdim     = maximum number of elements allowed in x array  (input)
c              (first dimension of parray)
c  kydim     = maximum number of elements allowed in y array  (input)
c              (second dimension of parray)
c
c  pxu(jx)   = 1-D independent variable x array from file     (output)
c  kxu       = number of elements in x array from file        (output)
c  cxname    = name of independent variable x                 (output)
c  pyu(jx)   = 1-D independent variable y array from file     (output)
c  kyu       = number of elements in y array from file        (output)
c  cyname    = name of independent variable x                 (output)
c  parray(jx,jy) = 2-D array from file                        (output)
c              (treated as a 1-D array in this routine)
c
c  kerr      = error indicator                          (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the file cannot be opened or read
c            = 2 reached the end of file before finding cvarname
c            = 3 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer kxdim, kxu, kydim, kyu, kunit, kdchar, kfchar, kerr
c
      real parray(*), pxu(*), pyu(*), zu(6)
c
      integer  iftype, ifdim
     &  , ixname, iyname, imaxname, iprint, ierr
     &  , iunit, iendfile, iscalar, iu, iumax, iskip
     &  , itemp, iline
     &  , ilow1, iupp1, iuppt
     &  , j, jx, jy, ju
c
c  iftype  = 1 for RPLOT file,  = 2 for U-file,  = 3 for concatenated
c
c  iprint  = integer to control amount of diagnostic output
c
c  ifdim   = 1 if a 1-D U-file is found, = 2 otherwise
c
c  iunit  = unit number for reading input file
c  iendfile = counter number of times end of input file has been reached
c
c  ixname  = number of characters in cxname
c  iyname  = number of characters in cyname
c  imaxname = maximum number of characters allowed in cvarname
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
c  ilow#, iupp# = lower and upper index of names in character strings
c
      character cline*72, cline1*72, cline2*72, ctemp*72
      character cvarname*(*), cxname*(*), cyname*(*)
c
c..set defaults
c
      ierr     = kerr
      kerr     = 0
      imaxname = 20
      iendfile = 0
c
      iprint   = 0
c
      iftype   = 0
      ifdim    = 2
      kxu      = 1
      kyu      = 1
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.0
      enddo
c
      do jx=1,kxdim*kydim
        parray(jx) = 0.0
      enddo
c
      cline  = ' '
      cline1 = ' '
      cline2 = ' '
      ctemp  = ' '
c
c..find first nonblank characters in cvarname
c
      iumax = max ( 1, min ( imaxname, len ( cvarname ) ) )
      do j=1,iumax
        ilow1 = j
        if ( cvarname(j:j) .ne. ' ' ) go to 5
      enddo
c
      write (*,*)
      write (*,*) 'Abort: no independent variable name given '
     &    ,cvarname
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
 5    iupp1 = ilow1
      do j=ilow1,iumax
        if ( cvarname(j:j) .eq. ' ' ) go to 6
        iupp1 = j
      enddo
 6    continue
c
      ctemp = ' ' // cvarname(ilow1:iupp1) // ' '
      iuppt = iupp1 - ilow1 + 3
c
c..test input unit number
c
      iunit = kunit
      if ( kunit .lt. 1  .or.  kunit .gt. 99 ) goto 90
c
c
c..test to see if cfile is an ASCII RPLOT file
c
      cline = ' '
      read (iunit,100,err=92,end=94) cline
 100  format (a)
c
      if ( cline(6:9).eq.'time' .or. cline(6:9).eq.'TIME' ) then
c
        iftype = 1
        rewind iunit
c
       kxu = 0
       kyu = 0
c
c..read next set of records in 2-D ASCII RPLOT file
c
  10    continue
c
c..read header line
c
        cline = ' '
        read (iunit,102,err=16,end=16) itemp, cline
 102    format (i4,1x,a)
        iline = iline + 1
c
c
        if ( itemp .lt. 1 ) then
          write (*,*) 'error -- number of elements < 1 for',cline
          write (*,*) 'on line ',iline,'  of input file'
          kerr = 1
          if ( ierr .lt. 0 ) stop
          rewind iunit
          return
        endif
c
c..is the next element part of the time array ?
c
        if ( cline(1:4).eq.'time' .or. cline(1:4).eq.'TIME' ) then
          kyu = kyu + 1
          read (iunit,*,err=14,end=14) pyu(kyu)
          iline = iline + 1
          go to 10
c
c..is the next element part of the array to be read ?
c  if so, there are six real numbers per line
c
        elseif ( index( cline, cvarname(ilow1:iupp1) ) .gt. 0 ) then
          kxu = itemp
          read (iunit,110,err=14,end=14) 
     &      (parray(j+(kyu-1)*kxdim),j=1,kxu)
  110     format(1x,(1p6e12.5))
          iline = iline + ( kxu + 5 ) / 6
          go to 10
c
c..is the next element part of the radius array to be read ?
c  if so, there are six real numbers per line
c
c        elseif ( cline .eq. crname ) then
c          read (iunit,110,err=14,end=14) (pradius(j),j=1,kxu)
c          iline = iline + ( itemp + 5 ) / 6
c          go to 10
c
c..if not, then skip lines
c
        else
          kxu = itemp
          itemp = ( kxu + 5 ) / 6
          do j=1,itemp
            read (iunit,*,err=14,end=14)
            iline = iline + 1
          enddo
          go to 10
c
        endif
c
c..bug out due to read error, try to recover
c
 14     continue
c
        kyu = kyu - 1
c
c..reached the end of the RPLOT file
c
 16     continue
c
        rewind iunit
c
        if ( kxu .lt. 1  .or.  kyu .lt. 1 ) then
          write (*,*) 'error in br2read'
          write (*,*) kxu,' = kxu'
          write (*,*) kyu,' = kyu'
          kerr = 3
          if ( ierr .lt. 0 ) stop
        endif
c
c..fill the pxu array
c  (For a better treatment, the index table would need to be read
c   to find out which independent variable should actually be used.)
c
        do j=1,kxu
          pxu(j) = j
        enddo
c
        return
c
c..read input file as an ASCII 2-D U-file
c
      else
c
        iftype = 2
        rewind iunit
c
c..read concatenated U-file till we get to the requested cvarname
c
  20  continue
c
      cline1 = cline2
      cline2 = cline
c
      read (iunit,100,err=92,end=28) cline
c
c..find dependent variable line
c
      if ( index( cline, '-DEPENDENT VARIABLE LABEL-' )
     &     .gt. imaxname ) then
c
c..check to see if this is a 1-D or 2-D U-file
c
        if ( index( cline1, '-INDEPENDENT VARIABLE LABEL' )
     &       .gt. imaxname ) then
          ifdim = 2
        else
          ifdim = 1
        endif
c
c..check to see if correct dependent variable name has been found
c
        if ( index( cline, ctemp(1:iuppt) ) .gt. 0
     &    .and. index( cline, ctemp(1:iuppt) )
     &    .lt. imaxname ) go to 21
c
c..skip lines as needed
c
      read (iunit,100,err=92,end=93) cline
      read (iunit,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iunit,104,err=94,end=93) kyu
      endif
c
      iskip = kxu/6 + kyu/6 + (kxu*kyu)/6
c
      do j=1,iskip
        read (iunit,*,err=92,end=93)
      enddo
c
      endif
c
c..go back and read more lines
c
      go to 20
c
c..cvarname has been found, now set cxname and cyname
c
 21   continue
c
      write (*,*) cline
c
      iumax = len ( cline1 )
      cxname = ' '
      iu = 0
      if ( ifdim .ne. 1  .and.   iumax .gt. 0 ) then
        do j=1,iumax
          if ( cline1(j:j) .ne. ' ' ) then
            iu = iu + 1
            cxname(iu:iu) = cline1(j:j)
          else
            if ( iu .gt. 0 ) go to 23
          endif
        enddo
      endif
 23   continue
c
      iumax = len ( cline2 )
      cyname = ' '
      iu = 0
      if ( iumax .gt. 0 ) then
        do j=1,iumax
          if ( cline2(j:j) .ne. ' ' ) then
            iu = iu + 1
            cyname(iu:iu) = cline2(j:j)
          else
            if ( iu .gt. 0 ) go to 24
          endif
        enddo
      endif
 24   continue
c
c..skip line for PROC CODE
c
      read (iunit,100,err=92,end=93) cline
c
c..read number of x and y points
c
      read (iunit,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iunit,104,err=94,end=93) kyu
      endif
c
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (*,*) kxu,' = kxu'
        write (*,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (*,*) 'Too many elements in 2-D file arrays'
        write (*,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (*,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn br2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iunit, *, err=96,end=93 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
      enddo
c
      if ( ifdim .ne. 1 ) then
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iunit, *, err=96,end=93 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
      enddo
      endif
c
c
c..read 2-D profile
c
      iu = 7
      do jy=1,max(kyu,1)
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
            read ( iunit, 122, err=96,end=93) ( zu(ju),ju=1,iumax )
  122       format(1x,6e13.6)
cbate            read ( iunit, *, err=96,end=93 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx+(jy-1)*kxdim) = zu(iu)
        enddo
      enddo
c
      return
c
 28   continue
c
c..reached the end of the input file
c  rewind input file and try again, or bug out
c
      if ( iendfile .lt. 1 ) then
        iendfile = iendfile + 1
        rewind iunit
        go to 20
      else
        go to 93
      endif
c
      endif
c
c
c..end of routine
c
      return
c
c
c..error conditions
c
  90  continue
c
      write (*,*)
     & 'error reading unit number ', kunit
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cline
      write (*,*) 'cvarname = ',cvarname
      write (*,*)
     & 'error reading line in UFILE in sbrtn br2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  93  continue
c
c..reached the end of the input file
c
      write (*,*)
     & 'reached the end of U-file in sbrtn br2read'
      write (*,*) ' before finding cvarname = ',cvarname
      kerr = 2
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (*,*)
      write (*,*) kxu,' = kxu'
      write (*,*) kyu,' = kyu'
      write (*,*)
     & 'error reading kxu or kyu in UFILE in sbrtn br2read'
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn br2read'
      write (*,*) 'cxname = ', cxname
      write (*,*) 'cyname = ', cyname
      write (*,*) 'cvarname = ',cvarname
      write (*,*) 'zu = ',(zu(j),j=1,iumax)
      write (*,*) 'jx = ',jx,'  jy = ',jy
      write (*,*) 'parray = ',(parray(j+(jy-1)*kxdim),j=1,jx)
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
      end

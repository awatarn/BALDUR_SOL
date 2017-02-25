c@b2read.f   02 Dec 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine b2read ( cvarname, kufile, kxdim, kydim
     & , pxu, kxu, cxname, pyu, kyu, cyname, parray, kerr )
c
c   This routine reads concatinated 2-D ASCII U-files (Boucher files).
c
c  cvarname  = name of the independent variable parray        (input)
c  kufile    = input unit number to read file (default 70)    (input)
c  kxdim     = maximum number of elements allowed in x array    (input)
c              (first dimension of parray)
c  kydim     = maximum number of elements allowed in y array    (input)
c              (second dimension of parray)
c
c  pxu(jx)   = 1-D independent variable x array from file     (output)
c  kxu       = number of elements in x array from file        (output)
c  cxname    = name of independent variable x                 (output)
c  pyu(jx)   = 1-D independent variable y array from file     (output)
c  kyu       = number of elements in y array from file        (output)
c  cyname    = name of independent variable x                 (output)
c  parray(jx,jy) = 2-D array from file                        (output)
c
c  kerr      = error indicator                            (input/output)
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
      integer kxdim, kxu, kydim, kyu, kufile, kdchar, kfchar, kerr
c
      real parray(kxdim,*), pxu(*), pyu(*), zu(6)
c
      integer  ixname, iyname, imaxname, iprint, ierr, ifdim
     &  , iufile, iendfile, iscalar, iu, iumax, iskip
     &  , ilow1, iupp1
     &  , j, jx, jy, ju
c
c
c  iprint  = integer to control amount of diagnostic output
c
c  ifdim   = 1 if a 1-D U-file is found, = 2 otherwise
c
c  iufile  = unit number for reading input file
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
      character cline*72, cline1*72, cline2*72
      character cvarname*(*), cxname*(*), cyname*(*)
c
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
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.0
        enddo
      enddo
c
      cline  = ' '
      cline1 = ' '
      cline2 = ' '
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
        if ( cvarname .eq. ' ' ) go to 6
        iupp1 = j
      enddo
 6    continue
c
c..open concatinated U-file
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) goto 90
c
c
c..read till we get to the requested cvarname
c
  10  continue
c
      cline1 = cline2
      cline2 = cline
c
      read (iufile,100,err=92,end=93) cline
 100  format (a)
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
        if ( index( cline, cvarname(ilow1:iupp1) ) .gt. 0
     &    .and. index( cline, cvarname(ilow1:iupp1) )
     &    .lt. imaxname ) go to 12
c
c..skip lines as needed
c
      read (iufile,100,err=92,end=93) cline
      read (iufile,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iufile,104,err=94,end=93) kyu
      endif
c
      iskip = kxu/6 + kyu/6 + (kxu*kyu)/6
c
      do j=1,iskip
        read (iufile,*,err=92,end=93)
      enddo
c
      endif
c
c..go back and read more lines
c
      go to 10
c
c..cvarname has been found, now set cxname and cyname
c
 12   continue
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
            if ( iu .gt. 0 ) go to 21
          endif
        enddo
      endif
 21   continue
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
            if ( iu .gt. 0 ) go to 22
          endif
        enddo
      endif
 22   continue
c
c..skip line for PROC CODE
c
      read (iufile,100,err=92,end=93) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iufile,104,err=94,end=93) kyu
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
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn b2read'
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
          read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
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
          read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
      enddo
      endif
c
c
c..read 1-D or 2-D profile
c
      if ( ifdim .eq. 1 ) then
c
        iu = 7
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, kxu + 1 - jx )
            read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
          parray(jx,1) = zu(iu)
        enddo
c
      else
c
        iu = 7
        do jy=1,max(kyu,1)
          do jx=1,kxu
            iu = iu + 1
            if ( iu .gt. 6 ) then
              iumax = min ( 6, kyu*kxu + 1 - (jy-1)*kxu - jx )
              read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
              iu = 1
            endif
          parray(jx,jy) = zu(iu)
          enddo
        enddo
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
     & 'error reading unit number ', kufile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (*,*)
      write (*,*) cline
      write (*,*)
     & 'error reading line in UFILE in sbrtn b2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  93  continue
c
c..reached the end of the input file
c  rewind input file and try again, or bug out
c
      if ( iendfile .lt. 1 ) then
        iendfile = iendfile + 1
        rewind iufile
        go to 10
      endif
c
      write (*,*)
     & 'reached the end of U-file in sbrtn b2read'
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
     & 'error reading kxu or kyu in UFILE in sbrtn b2read'
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (*,*)
      write (*,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn b2read'
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
      end

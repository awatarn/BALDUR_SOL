c*/ 17:00 30-sep-89 /11040/bald89/wbaldn1 wsutil DPRTUTIL, Bateman
c BALDUR  file DPRTUTIL library of printout routines by Bateman, PPPL
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c**********************************************************************c
c
c contents:
c abortb:  terminate run with a message and leave dropfile
c          (moved to separate file abortb.f)
c header:  creates column headers out of the character string head
c prtlns:  this subroutine skips nlines blank lines on output nout
c prtrs:   this subroutine prints a real number on output number nout
c prtis:   this subrtn prints integer number n on output number nout
c prtiv:   print out 2D array of integers out of a 3D array
c prtrv:   print out 2D array of real values out of a 3D array
c prtr1d:  print out text followed by 1D real array
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@header
c
       subroutine header (nout,head)
c
c   This sbrtn creates column headers out of the character string head
c  which must be of the form text1,#,text2,#,...
c  where text* is the text to appear at the top of each column and
c  # is an integer representing the width of each column.
c  The total width must be no more than 132 columns wide.
c
      character head *(*), txtout*132, txtpos*1
c
c  head   = input character string
c  txtout = ouput character string printed out
c  txtpos = 'l' for left justified text
c         = 'c' for centerd text (default)
c         = 'r' for right justified text
c
      length = len (head)
      if ( length .lt. 1 ) return
      length = min ( length, 132)
c
cd      txtout = ' txtout example' // ' ' // head
cd      write (nout,*) txtout,' = txtout'
c
      itext  = 1             ! position of next character in 'head'
      txtout = ' '           ! first output character is blank
      iadjst = 0             ! adjust width of first column
      istout = 1             ! position of next character in txtout
      istmax = 132           ! maximum value allowed for istout
      txtpos = 'c'           ! centered text is initial default
c
cd      write (nout,*) ' After initialization'
cd      write (nout,*) head,' = head'
cd      write (nout,*) length,' = length'
c
  10  continue               ! search for next character substring
c                            ! find next occurance of ','
      inc = index ( head( min(itext,length) :length), ',' )
c
cd      write (nout,*) head( min(itext,length) :length)
cd     &  ,' = head( min(itext,length) :length)'
cd      write (nout,*) itext,' = itext  ',inc,' = inc'
cd     &  ,istout,' = istout'
c
      if ( inc .gt. 1 ) then ! handle text or column width
c
        itext1 = itext       ! position of first character in substring
c                          ! position of last character in substring
        itext2 = inc + itext1 - 2  
        iwdtxt = itext2 - itext1 + 1  ! width of text in field
c
cd        write (6,*) itext1,' = itext1,  ',itext2,' = itext2  '
cd     &   ,inc,' = inc'
c
        itext = itext + inc
c
c                            ! find next occurance of ','
        inc = index ( head( min(itext,length) :length), ',' )
c
        iwidth = 0
        inum1  = itext  ! first position of number field
        inum2  = inum1+inc-2  ! last position of number field
        if ( inc .lt. 1 ) inum2 = length  ! last number in string
c
        if ( inum2 .ge. inum1 ) then
          do 20 j=inum1,inum2
c
            if ( head(j:j) .eq. 'l' ) txtpos = 'l'
            if ( head(j:j) .eq. 'c' ) txtpos = 'c'
            if ( head(j:j) .eq. 'r' ) txtpos = 'r'
c
            if ( ichar( head(j:j) ) .ge. ichar ('0') .and.
     &           ichar( head(j:j) ) .le. ichar ('9') )
     &        iwidth = 10 * iwidth
     &                 + ichar(head(j:j)) - ichar('0')
  20      continue
        endif
c
        iwidth = max ( 1, iwidth)    ! field must be at least 1 col wide
        iwidth = min ( istmax - istout, iwidth)  ! max num col remaining
        iwidth = iwidth + iadjst   ! adjust column width if needed
c
        if ( iwidth .gt. 0 .and. iwidth .lt. istmax ) then
c
          iwdtxt = min ( iwidth-1, iwdtxt )  ! width of text
c
          if ( txtpos .eq. 'l' ) ifill1 = 0  ! left justified
          if ( txtpos .eq. 'c' ) ifill1 = ( iwidth - iwdtxt ) / 2
          if ( txtpos .eq. 'r' ) ifill1 = iwidth - iwdtxt ! right justified
c
          txtout(istout+ifill1:istout+ifill1+iwdtxt)
     &        = head(itext1:itext1+iwdtxt-1)
c
          istout = istout + iwidth     ! move to next part of txtout
c
cd          write (nout,*) inum1,' inum1  ',inum2,' = inum2  '
cd     &      ,inc,' = inc  '
cd          write (nout,*) iwidth,' = iwidth  ',iwdtxt,' iwdtxt  '
cd     &      ,ifill1,' = ifill1  '
cd          write (nout,*) head(itext1:itext1+iwdtxt-1)
cd     &      ,' = head(itext1:itext1+iwdtxt-1)'
cd          write (nout,*) txtout,' = txtout'
c
        endif                ! next text or column width
c
        itext = itext + inc  ! next character string 'head' after ','
        iadjst = 0
c
      else                   ! no more occurances of ',' in head
c
        write (nout,*) txtout
        return
c
      endif
c
      go to 10               ! go on to next field
c
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtlns
c
      subroutine prtlns (nout,nlines)
c
cx prtlns:  this subroutine skips nlines blank lines on output nout
c**
      do 10 l=1,nlines
      write (nout,100)
 100  format (1h )
  10  continue
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtrs
c
        subroutine prtrs (nout,x,textin)
c
cx prtrs:  this subroutine prints a real number on output number nout
c    together with a message
c
      character *(*) textin
c
      write (nout,100) x,textin
 100  format (1x,1pe12.5,' = ',a)
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtis
c
        subroutine prtis (nout,n,textin)
c
cx prtis:  this subrtn prints integer number n on output number nout
c    together with a message
c
      character *(*) textin
c
      write (nout,100) n,textin
 100  format (1x,i12,' = ',a)
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtiv
c
      subroutine prtiv (nout,n,ldim,idim,l1,i1,j1,i2,j2,ncol,textin)
cx prtiv
c  *****
c  nout is the unit number for the output
c  n is the three-dimensional integer array to be printed out
c   the last two dimensions are printed as a 2-d array
c  ldim,idim  are the dimensions of the first two indicies
c  i1,i2,j1,j2  are the lower and upper bounds of the portion of the
c   array to be printed out
c  ncol is the number of array columns printed across the page
c   the widest printout takes 10 + 10*ncol printed columns
c  textin  called with "text" in the argument list, used as header
c  ----
      dimension n(ldim,idim,1)
c
      character *(*) textin
c
      write (nout,100) n,textin
 100  format (1x,i12,' = ',a)
c
c
      iend = i1 - 1
      do 40 k=1,20
      iin = iend + 1
      iend = min0 (iend+ncol,i2)
      icol = iend - iin + 1
      if (icol .lt. 1) go to 990
c
      jform = mod (k-1,5) + 1
      go to (2,4,6,8,10), jform
c
 990  write (nout,*) 'error in prtiv'
      return
c
   2  continue
      do 12 j=j1,j2
  12  write (nout,902) (n(l1,i,j),i=iin,iend)
      go to 30
   4  continue
      do 14 j=j1,j2
  14  write (nout,904) (n(l1,i,j),i=iin,iend)
      go to 30
   6  continue
      do 16 j=j1,j2
  16  write (nout,906) (n(l1,i,j),i=iin,iend)
      go to 30
   8  continue
      do 18 j=j1,j2
  18  write (nout,908) (n(l1,i,j),i=iin,iend)
      go to 30
  10  continue
      do 20 j=j1,j2
  20  write (nout,910) (n(l1,i,j),i=iin,iend)
c
  30  continue
      if (j2-j1+1 .ne. 1)  call prtlns (nout,2)
      if (iend .ge. i2) return
  40  continue
c
 902  format (2x,10(i11))
 904  format (4x,10(i11))
 906  format (6x,10(i11))
 908  format (8x,10(i11))
 910  format (10x,10(i11))
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtrv
c
c     *********************
c     * SUBROUTINE prtrv  *
c     *********************
c
c     usage:
c           CALL prtrv (nout,a,ldim,idim,l1,i1,j1,i2,j2,maxcol,title)
c
c     purpose:
c           prtrv prints out any 2-D part of a multidimensional array
c
c     parameters:
c           nout  ...device # for printout
c           A     ...3D array  A(ldim,idim,1)
c           title...any character string
c
c     nout is the unit number for the output
c     a is the three-dimensional real array to be printed out
c          the last two dimensions are printed as a 2-d array
c     ldim,idim  are the dimensions of the first two indicies
c     i1,i2,j1,j2  are the lower and upper bounds of the portion of the
c           array to be printed out
c     maxcol = max number of printed columns across the page
c           (typically 80 or 120)
c     title is a character string which is printed out as a header
c
c***********************************************************************
c
      subroutine prtrv (nout,a,ldim,idim,l1,i1,j1,i2,j2,maxcol,title)
c
      dimension a(ldim,idim,1)
      character *(*) title
c     ------------------------------------------------------------------
c
      write (nout,*)
      write (nout,*) title
c
      ncol = max((maxcol - 8)/11,1)
      ncol = min(ncol,10)
      iend = i1 - 1
      do 40 k=1,20
      iin = iend + 1
      iend = min0 (iend+ncol,i2)
      icol = iend - iin + 1
      if (icol .lt. 1) go to 990
c
      jform = mod (k-1,5) + 1
      go to (2,4,6,8,10), jform
c
 990  write (nout,*) ' error in prtrv'
      return
c
   2  continue
      write (nout,912) (i,i=iin,iend)
      do 12 j=j1,j2
  12  write (nout,902) j,(a(l1,i,j),i=iin,iend)
      go to 30
   4  continue
      write (nout,914) (i,i=iin,iend)
      do 14 j=j1,j2
  14  write (nout,904) j,(a(l1,i,j),i=iin,iend)
      go to 30
   6  continue
      write (nout,916) (i,i=iin,iend)
      do 16 j=j1,j2
  16  write (nout,906) j,(a(l1,i,j),i=iin,iend)
      go to 30
   8  continue
      write (nout,918) (i,i=iin,iend)
      do 18 j=j1,j2
  18  write (nout,908) j,(a(l1,i,j),i=iin,iend)
      go to 30
  10  continue
      write (nout,920) (i,i=iin,iend)
      do 20 j=j1,j2
  20  write (nout,910) j,(a(l1,i,j),i=iin,iend)
c
  30  continue
      if (j2-j1+1 .ne. 1)  write (nout,*)
      if (iend .ge. i2) return
  40  continue
c
 902  format (i4,10(1pe11.3))
 904  format (1x,i4,10(1pe11.3))
 906  format (2x,i4,10(1pe11.3))
 908  format (3x,i4,10(1pe11.3))
 910  format (4x,i4,10(1pe11.3))
c
 912  format (4x,10(i4,7x))
 914  format (5x,10(i4,7x))
 916  format (6x,10(i4,7x))
 918  format (7x,10(i4,7x))
 920  format (8x,10(i4,7x))
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@prtr1d
c
c     *********************
c     * SUBROUTINE PRTR1D *
c     *********************
c
c     usage:
c           CALL PRTR1D (nout,pa,nmax,title)
c
c     purpose:
c           prtr1d prints out 1-D array pa(j), j=1,nmax
c
c     parameters:
c           nout  ...device # for printout
c           pa    ...real-valued 1-D array
c           nmax  ...number of elements to be printed out
c           title...any character string
c
c***********************************************************************
c
      subroutine prtr1d (nout,pa,nmax,title)
c
      dimension pa(1)
      character *(*) title
c     ------------------------------------------------------------------
c
      write (nout,*)
      write (nout,*) title
c
      write (nout,110) (pa(j),j=1,nmax)
 110  format (2x,1p5e13.5)
c
      return
      end



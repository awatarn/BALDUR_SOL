c@strim.f   6 Jan 1995   Glenn Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine strim ( string, strout, nout)
c
c     This routine finds the first non-blank string of characters
c
c  string  = input character string
c  strout  = output character string
c  nout    = number of non-blank characters in strout
c
      implicit none
c
      character string*(*), strout*(*)
c
      integer nout
c
      integer ichar1, ichar2, imax, j
c
c  imax  = number of characters in the input string
c
      imax = len ( string )
c
      strout = ''
      nout = 0
c
      if ( imax .lt. 1 ) return
c
c..find the number of non-blank characters in string
c
      do j=1,imax
        if ( string(j:j) .ne. ' ' ) then
          nout = nout + 1
          strout(nout:nout) = string(j:j)
        else
          if ( nout .gt. 0 ) return
        endif
      enddo
c
c
      return
      end

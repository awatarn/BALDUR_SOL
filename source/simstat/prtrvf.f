c@prtrvf.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c   Print a real variable in  0pf13.4 format
c     followed by a line of text to output unit and to the screen.
c
      subroutine prtrvf ( nout, px, text )
c
      real px
c
      character *(*) text
c
      write ( nout, 110 )  px, text
      write ( *,    110 )  px, text
 110  format (0pf13.4,2x,(a))
c
      return
      end


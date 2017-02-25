c@prtext.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c   Print a line of text to output unit and to the screen.
c
      subroutine prtext ( nout, text )
c
      character *(*) text
c
      write ( nout, * ) text
      write ( *, * ) text
c
      return
      end


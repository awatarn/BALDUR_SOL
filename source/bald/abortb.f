c@abortb  Bateman 24 July 1992
c
      subroutine abortb (nout,text)
c
c   Subrtn abortb prints out text and stops run
c
      character *(*) text
c
      write (nout,*) text
c
      stop
      end

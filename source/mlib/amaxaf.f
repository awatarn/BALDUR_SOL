c@amaxaf.f  Glenn Bateman, PPPL
c
c  Find the largest value in an array.
c
c     usage:
c           xmax = amaxaf (array,ifirst,ilast)
c
c     parameters:
c           array:  real array
c           ifirst: first element of the array to be searched
c           ilast:  last element of the array to be searched
c           xmax:   maximum value of the array
c
c  Note:  This routine replaces the routine by the same name in LIBMATH.
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      function amaxaf (array,ifirst,ilast)
c
      real amaxaf, array(*)
c
      integer ifirst, ilast, istep, j
c
      amaxaf = array(ifirst)
c
      istep = isign ( 1, ilast - ifirst )
c
      do j=ifirst,ilast,istep
c
        amaxaf = max ( amaxaf, array(j) )
c
      enddo
c
      return
      end

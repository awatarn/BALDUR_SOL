c@lug.f  .../baldur/code/bald/lug.f  Glenn Bateman, Lehigh
c rgb 04-sep-96 removed return after go to 10
c
c  Find the location in an table with table(index-1) .le. x .lt. table(index).
c  index = 1 if x < table(1).
c  index = n+1 if x .ge. table(n)
c
c     usage:
c           index = lug (x,table,n,iguess)
c
c     parameters:
c           table:  real array with strictly increasing values
c                   table(1) < table(2) < ... < table(n)
c           n:      length of the array table to be searched
c           iguess: guess for the index
c           index:  maximum value of the table
c
c  Note:  This routine replaces the routine by the same name in LIBMATH.
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      function lug (x,table,n,iguess)
c
      real x, table(1)
c
      integer lug, iguess, imn, imx
c
      lug = 1
c
      if ( n .lt. 1  .or.  x .lt. table(1) ) then
        lug = 1
        iguess = lug
        return
      endif
c
      if ( x .ge. table(n) ) then
        lug = n+1
        iguess = lug
        return
      endif
c
      imn = 2
      imx = n
c
      iguess = max ( 2, min ( n, iguess ) )
c
      if ( table(iguess-1) .le. x ) imn = iguess
      if ( x .lt. table(iguess) )   imx = iguess
c
  10  continue
c
      if ( imn .ge. imx ) then
        lug = imx
        iguess = lug
        return
      endif
c
      iguess = ( imn + imx ) / 2
c
      if ( x .lt. table(iguess) ) then
        imx = iguess
      else
        imn = iguess + 1
      endif
c
      go to 10
c
      end

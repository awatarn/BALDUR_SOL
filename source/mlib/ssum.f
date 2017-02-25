c--------1---------2---------3---------4---------5---------6---------7-c
      function ssum(n,sx,incx)
c
c     ******************************************************************
c     *  ssum adds up the elements of the array sx.                    *
c     *  last revision: 3/81 w.a.houlberg and s.e.attenberger ornl.    *
c     *  other comments:                                               *
c     *  use  omnilib version for cray optimization.                   *
c     ******************************************************************
c
      dimension sx(*)
c
      sum=0.0
      do i=1,n,max(incx,1)
        sum = sum + sx(i)
      enddo
      ssum=sum
c
      return
      end

c@ismax  .../baldur/code/util/ismax.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      function ismax ( n, sx, incx )
c
c    Find index of maximum value
c  n     = number of elements to be searched
c  sx(j) = vector
c  incx  = skip distance between elements in sx(j)
c        = 1 for contiguous elements
c
      integer   n, incx, ix, imax, j
      real      sx(*), zmax
c
      ix = 1
      zmax = sx(1)
      imax = 1
      do j=1,n
        if ( sx(ix) .gt. zmax ) then
           zmax = sx(ix)
           imax = j
        endif
        ix = ix + max(incx,1)
      enddo
      ismax = imax
c
      return
      end

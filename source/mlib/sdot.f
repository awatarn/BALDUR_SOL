c@sdot  .../baldur/code/util/sdot.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      function sdot (n, sx, incx, sy, incy)
c
c    Inner product of two vectors
c  n     = number of elements in the product
c  sx(j) = vector
c  incx  = skip distance between elements in sx(j)
c        = 1 for contiguous elements
c  sy(j) = vector
c  incy  = skip distance between elements in sy(j)
c        = 1 for contiguous elements
c
      integer   n, incx, incy, ix, iy, j
      real   sx(*), sy(*), sum
c
      sum = 0.0
      ix = 1
      iy = 1
      do j=1,n
           sum = sum + sx(ix)*sy(iy)
           ix = ix + max(incx,1)
           iy = iy + max(incy,1)
      enddo
      sdot = sum
c
      return
      end

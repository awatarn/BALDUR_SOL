c@scopy  .../baldur/code/util/scopy.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine scopy (n, sx, incx, sy, incy)
c
c    Copy a real vector
c  n     = number of elements to be copied
c  sx(j) = vector to be copied
c  incx  = skip distance between elements in sx(j)
c        = 1 for contiguous elements
c  sy(j) = result vector
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
           sy(iy) = sx(ix)
           ix = ix + max(incx,1)
           iy = iy + max(incy,1)
      enddo
c
      return
      end

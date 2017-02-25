c@saxpy  .../baldur/code/util/saxpy.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine saxpy (n, a, sx, incx, sy, incy)
c
c    Compute sy(j) = sy(j) + a * sx(j), where
c  n     = number of elements in the operation
c  a     = scalar constant
c  sx(j) = vector
c  incx  = skip distance between elements in sx(j)
c        = 1 for contiguous elements
c  sy(j) = vector
c  incy  = skip distance between elements in sy(j)
c        = 1 for contiguous elements
c
      integer   n, incx, incy, ix, iy, j
      real      a, sx(*), sy(*)
c
      ix = 1
      iy = 1
      do j=1,n
           sy(iy) = sy(iy) + a * sx(ix)
           ix = ix + max(incx,1)
           iy = iy + max(incy,1)
      enddo
c
      return
      end

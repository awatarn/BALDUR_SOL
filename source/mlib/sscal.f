c@sscal  .../baldur/code/util/sscal.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine sscal (n, a, sx, incx )
c
c    Compute sx(j) = a * sx(j), where
c  n     = number of elements in the operation
c  a     = scalar constant
c  sx(j) = vector
c  incx  = skip distance between elements in sx(j)
c        = 1 for contiguous elements
c
      integer   n, incx, ix, j
      real      a, sx(*)
c
      ix = 1
      do j=1,n
           sx(ix) = a * sx(ix)
           ix = ix + max(incx,1)
      enddo
c
      return
      end

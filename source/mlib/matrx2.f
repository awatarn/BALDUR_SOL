c@matrx2  .../baldur/code/bald/matrx2.f
c  pis 12-jun-98 removed kdim = kdim
c******************************************************************************
c
        subroutine matrx2(px,py,plu,kpivot,kdim,kn,kx1,kx2,ky)
c
c
cl      z.12    solve l-u equations
c
c
c------------------------------------------------------------------------------
c
c
c       a . x(kx1+j) = l . u . x(kx1+j)  =  y(ky+j)
c
c       where  kx1 <= kx1+j <= kx2, and x(i) and y(i) refer to the
c       i-th column vectors in  x and y, respectively.
c
c
c       we solve  l.u.x = y by first solving  l.x' = y and then solving
c       u.x = x' .
c       because of the indexing, we may solve a.x = y where
c       x and y are n by n matrices by setting
c       kx1 = ky1 = 1, and kx2 = kn.
c
c                       *********
c                       * note  *
c                       *********
c
c       a  *must*  have been decomposed into a l-u matrix by matrx1
c       *first*.
c
c
c-----------------------------------------------------------------------
c
c
        dimension
     i  kpivot(kdim),
     r  plu(kdim,kdim), px(kdim,kdim),  py(kdim,kdim),
     r  zx(20)
c
c
cl      1)      initialize internal variables
c
c
c        kdim = kdim   ! this (strange) line may cause segmentation faults!!!!
        in = kn
        ix1 = kx1
        ix2 = kx2
        iy = ky - kx1
        i0 = in - 1
c
c
cl      2) solve for px, using a major do-loop
c
c
        do 238 j0 = ix1, ix2
c
c
cl      2.1)    solve  l.zx = y(j0+iy)
c
c       note that iy = ky - ix1, so when j0 = kx1, j0 + iy = ky
c
c
        zx(1) = py(kpivot(1),j0+iy)
c
        do 210 j1 = 2, in
                ipiv = kpivot(j1)
                i3 = j1 - 1
                zsum = 0.0
                do 208 j2 = 1, i3
                        zsum = zsum + plu(ipiv,j2)*zx(j2)
  208           continue
                zx(j1) = py(ipiv,j0+iy) - zsum
  210   continue
c
c
cl      2.2)    solve  u.px(j0) = zx
c
c       same substitution, except working from row n back to row 1
c
c
        px(in,j0) = zx(in) / plu(kpivot(in),in)
c
        do 230 j1 = 1, i0
                i2 = in - j1
                ipiv = kpivot(i2)
                i3 = i2 + 1
                zsum = 0.0
                do 228 j2 = i3, in
                        zsum = zsum + plu(ipiv,j2)*px(j2,j0)
  228           continue
                px(i2,j0) = (zx(i2) - zsum) / plu(ipiv,i2)
  230   continue
  238   continue
        return
        end

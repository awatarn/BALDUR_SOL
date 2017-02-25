c@matrx1  .../baldur/code/bald/dolymp.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c
        subroutine matrx1(pmatrx,kdim,kn,kpivot,kerror)
c
c
cl      z.11)   decomposes pmatrx into l-u matrix
c
c
c
c
c----------------------------------------------------------------------------
c
c
c       matrx1 decomposes pmatrx in place into a l-u decomposed matrix.
c       the method is gaussian elimination with partial pivoting,
c       using kpivot to hold the permutation vector indicating the pivots.
c       hence, pmatrx(kpivot(i),j) will be a l-u matrix,
c       and kpivot indicates the row switches that were done in pivoting.
c
c       this algorithm was modified from the algorithm
c       in %%%?????
c
c
c------------------------------------------------------------------------------
c
c
c
c       this subroutine does the first part of an algorithm to
c       solve equations of the type:
c
c               a . x = y
c
c       where  a  is a given  n by n matrix,  y  a given n by m
c       matrix,  x  an unknown n by m matrix, and
c       "." denotes matrix multiplication.
c
c
c
c       the scheme used is l-u decomposition, done by gaussian elimination.
c       the essence of the scheme is to decompose matrix  a  into  l.u,
c       where "." denotes matrix multiplication, and
c       l is of the form:                       and u:
c
c               1 0 0 . . . 0 0 0               * * * . . . . * *
c               * 1 0 0 . . . 0 0               0 * * . . . . * *
c               * * 1 0 . . . 0 0               0 0 * * . . . * *
c               * * * . . . . . 0               0 0 0 * . . . * *
c               . . . . . . . . .               . . . . . . . . .
c
c               * * * . . * 1 0 0               0 0 . . . 0 * * *
c               * * * . . * * 1 0               0 0 . . . . 0 * *
c               * * * . . . . . 1               0 0 . . . . 0 0 *
c
c
c       i.e.,  l is a lower triangular matrix with 1's along the diagonal,
c       and u is an upper triangular matrix.
c       once these matrices are found,  l.x' = y and  u.x = x' are easily
c       solved ( since x'(1) = y(1), x'(2) = y(2) - x'(1)*l(2,1),
c       x(n) = x'(n) / u(n,n), x(n-1) = (x'(n-1) - x(n)*u(n-1,n))/ u(n-1,n-1),
c       etc.).  note that this requires that none of the diagonal entries
c       in u be near or equal to zero.
c
c       in this system, u and l may be squeezed into the same matrix,
c       since l(j,j) is always 1, hence need not be stored.
c       l(i,j) = 0 for i < j, and u(i,j) = 0 for i > j.
c       hence lu(i,j) = l(i,j) for i > j, and = u(i,j) otherwise.
c
c       to get this l-u decomposition matrix, gaussian elimination is used.
c       the idea is to do column operations of  a  to generate l, and keep
c       record of what you have done in u,  so that l.u will reconstruct  a.
c       when we are working on column j, we have l(i,k) = 0 for
c       i < j and i < k,  so we are trying to set row j to 0 in all
c       columns after j.
c       first we divide column j by l(j,j) so l(j,j) becomes 1.
c       meanwhile, we fill in  u(j,j) to be equal to the old value of l(j,j).
c       next, we subtract from column i column j times l(j,i).
c       i.e.,  l(k,i) = l(k,i) - l(j,i)*l(k,j) for i > j.  note that the
c       subtraction is unnecessary for k < j.
c       hence, l(j,i) will become 0.  we set u(j,i) equal to
c       the old value of l(j,i).
c       later, when we finish, we find that a(i,j) = sum for i>=k and
c       j>=k of  l(i,k)*u(k,j).
c       this turns out to be just the inversion of the column operations
c       that were performed on a to produce l.
c
c       pivoting:
c
c       the process is actually more complex, since "pivoting" is used.
c       remember that  u.x = x' cannot be solved if any of the diagonal
c       elements of u is 0 or near 0.
c       pivoting is simply performing row interchanges so that this danger
c       is avoided.  specifically,  just before dividing the column by
c       l(j,j), one searches the remaining rows (since the previous rows
c       have been set to 0) for the row  k with the maximum
c       abs(l(k,j)).  then rows k and j are interchanged, and the column
c       is divided by the new l(j,j).  then, when solving l.x' = y,
c       y must be appropriately permuted, and similarly x' when solving
c       u.x = x'.
c       this is done in practice by indexing  lu by a vector kpivot,
c       so one decomposes the matrix  lu(kpivot(i),j).
c       when a row interchange is desired, one simply interchanges elements
c       of kpivot.  this vector is later used to untangle the pivoting
c       when solving for x' and x.
c
c
c
c==============================================================================
c
c
c
c
       include 'cparm.m'
       include 'cbaldr.m'
c
        dimension
     i  kpivot(kdim),
     r  pmatrx(kdim,kdim),
     r  zscale(20)
c
c
cl      1)      initialize kpivot and zscale
c
c       1.1) kpivot:
cc
        if(kdim.gt.20) go to 9010
        do 114 j = 1, kn
        kpivot(j) = j
c
c       1.2) zscale:
c
        z0 = 0.0
        do 110 j1 = 1, kn
        if(z0.le.abs(pmatrx(j,j1))) z0 = abs(pmatrx(j,j1))
  110   continue
        if (z0.le.epslon) go to 9011
        zscale(j) = 1.0 / z0
  114   continue
c
c
cl      2)       pivoting
c
c
        i0 = kn - 1
        do 218 j = 1, i0
        z0 = 0.0
c
c       2.1) set i1 to the index of the row with the greatest
c        scaled magnitude element in lower column j
c
                do 204 j1 = j, kn
                        z1 = abs(pmatrx(kpivot(j1),j)) *
     1                          zscale(kpivot(j1))
                        if (z1.le.z0) go to 204
                                z0 = z1
                                i1 = j1
  204           continue
c
c       2.2) error exit if maximum is almost 0
c
                if (z0.le.epslon) go to 9020
c
c       2.3) interchange rows j and i1
c
                if (i1.eq.j) go to 208
                i2 = kpivot(j)
                kpivot(j) = kpivot(i1)
                kpivot(i1) = i2
  208           continue
c
c       gaussian elimination: subtract column j from remaining columns,
c        with appropriate factor zpivot thrown in
c
c       3.1) set zpivot
c
        i1 = kpivot(j)
        zpivot = 1.0 / pmatrx(i1,j)
        i2 = j + 1
c
c       3.2) scaling and subtracting:
c
        do 214 j1 = i2, kn
                i3 = kpivot(j1)
                z0 = pmatrx(i3,j)*zpivot
                pmatrx(i3,j) = z0
                do 214 j2 = i2, kn
                        pmatrx(i3,j2) = pmatrx(i3,j2) - z0*pmatrx(i1,j2)
  214   continue
  218   continue
c
c       4) set kerror
c
        if (abs(pmatrx(kpivot(kn),kn)).le.epslon) go to 9030
        kerror = 0
        return
c
c
 9010   continue
        call error_olymp(1,99,11,1,
     1          ' *** error *** call to matrx1 to invert matrix ')
        call error_olymp(1,99,11,1,
     1          ' *** larger than 20 by 20 ')
        call error_olymp(2,kdim,2,1,'  kdim  ')
        return
c
 9011   continue
        kerror = 1
        return
c
 9020   continue
        kerror = 2
        return
c
 9030   continue
        kerror = 3
        return
        end

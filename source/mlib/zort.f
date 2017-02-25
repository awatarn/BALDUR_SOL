      SUBROUTINE ZORT(x,ieq,Indx)
C
C Sort the elements of x and return the pointer in Indx
C
      REAL    x(*), xmin
      INTEGER temp, Indx(10), j, jj
C== Sort the values by x
C==
      do j = 1,ieq
         Indx(j)=j
      enddo
      do jj = 1,ieq-1
         xmin=x(Indx(jj))
         do j=jj,ieq
            IF(x(Indx(j)).GT.xmin)THEN
                temp    = Indx(j)
                Indx(j) = Indx(jj)
                Indx(jj)= temp
                xmin    = x(Indx(jj))
            ENDIF
         enddo
      enddo
      RETURN
      END

c--------1---------2---------3---------4---------5---------6---------7-c
c@ncsolv  /11040/baldur/code/bald/dimprad.f
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 15-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine ncsolv (ik,ji)
c
c       2.20.5   solver-routine for non-corona radiation
c
      include 'comncr.m'
c
         n1 = n2 - 1
         eee(1) = aaa(1)/bbb(1)
         fff(1) = ddd(1)/bbb(1)
c
      do 10 j=2,n2
         i = j - 1
         zqq = 1.0/(bbb(j)-ccc(j)*eee(i))
         eee(j) = zqq * aaa(j)
         fff(j) = zqq * (ddd(j)+ccc(j)*fff(i))
   10 continue
c
         xn(n2,ik,ji) = fff(n2)
c
      do 20 j=1,n1
         jj = n2 - j
         xn(jj,ik,ji) = eee(jj)*xn(jj+1,ik,ji) + fff(jj)
   20 continue
c
      return
      end

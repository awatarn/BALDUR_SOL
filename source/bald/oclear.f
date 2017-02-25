c@oclear  .../baldur/code/bald/oclear.f 
c rgb 22-nov-99 oclear generated using cleargen 
c fgps 31-mar-83 added "mximp" to argument list for "otable(---)", 
c                    where inside subroutine otable  kximp=mximp. 
c fgps 11-mar-83 file "corona.for" modified to handle 4 impurities. 
c 
c***note*** 
c       to change the impurity-species capacity do the following, where 
c       mximp=(the maximum number of impurity species allowed): 
c      1.re-dimension:  ocoef(6,5*mximp,3), orange(2,5*mximp,3) in 
c                       /comcor/. 
c      2.re-dimension locally in subroutine otable(---): 
c                       zcperm(10,5*mximp), zrperm(2,5*mximp), 
c                       irange(mximp,3); and, correspondingly 
c                       icp /50*mximp/, irp /10*mximp/, ira /3*mximp/. 
c      3.re-dimension locally in subroutine oclear: 
c                       real81(120*mximp), and 
c                  call resetr(real81,120*mximp,0.0). 
c      3.re-define in subroutine otable(---): 
c                       mxorng=5*mximp 
c                       if(mimp.gt.0) mxorng=(5*mximp)/mimp. 
c      4.re-adjust in subroutine olist(---): 
c                  call raray3(8hocoef   ,6,5*mximp,mocoef,5*mximp,3) 
c                  call raray3(8horange  ,2,5*mximp,2,5*mximp,3). 
c 
c       incidentally, the use of "klabel(15)" rather than "klabel(12)" 
c       prevents the truncation of part of a title line.  to avoid 
c       possible array overlap, "lholab(15)" in "/comflg/" must also 
c       be dimensioned out to "15".  moreover, to be consistent "10141 
c       format(---)" in subroutine mprint must contain "15a4". 
c 
c--------1---------2---------3---------4---------5---------6---------7-c 
c 
         subroutine oclear 
c 
c 1.2  clear variables and arrays 
c 
c
c  cleargen 8:17 22 Nov 99
c
c sbrtn clear was generated using the cleargen program
c
      ocoef = 0.0
      orange = 0.0
      mocoef = 0
      mxorng = 0

c 
      return 
      end 

c@blgcom.m  common blocks for ballooning mode package balgam1
c  by Michael W. Phillips  05-jul-88
c
      parameter (LTHE=200,LTHE2=2*LTHE,LTOT=2500)
c
      common /bal1/ pnew,pold,nful,nful2,nthe2,nodes,xnode,dt,the0
     &,u0,pi,tpi,rdt,gb,gpb,slen,npts,i0,ihalf,ithe,tdiff
c
      common /bal2/ cost(LTHE),gpsi2(LTHE)
     &,bp2(LTHE),cappa(LTHE)
     &,bsqr(LTHE),curvn(LTHE),curvs(LTHE),galpha(LTHE2)
     &,gpsi(LTHE),xjac(LTHE2),cq0(LTHE2),cq1(LTHE2),cq2(LTHE2)
     &,dq0(LTHE2),dq1(LTHE2)
c
      common /bal3/ wtemp(LTHE),ztemp(LTHE)
c
      common /bal4/ a(LTOT),b(LTOT),c(LTOT),p(0:LTOT)
     &,d(LTOT)
c

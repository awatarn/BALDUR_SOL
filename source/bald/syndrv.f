cahk
c--------1---------2---------3---------4---------5---------6---------7-c
c@syndrv  /11040/baldur/code/bald/dimprad.f
c rgb 03-jun-92 replace aminaf and amaxaf with do loop to find bmin and bmax
c
      subroutine syndrv
c
c   called from dimprad
c   driver for tamor's synchrotron radiation package
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'csynrd.m'
c
c      common/insink/ nzons,ree,roo,reo,roe,bmin,bmax,pwall,ptot,areaw,
c     2  taucrt,nw1,ndbug,iedit,
c     2  stemp(jx),sdens(jx),bavg(jx),svol(jx),surf(jx),wsyn(jx)
c   input in comin2
      iedit=iedit1
      ndbug=ndbug1
c   nzons = 25-->mzones=52
      nzons=nzons1
      taucrt=tacrt
      areaw=areaw1
      ree=ree1
      roo=roo1
      reo=reo1
      roe=roe1
c   kludge to avoid blowing up at t=0 when avi=0 - skip calc
      if(tbi.lt.tinit*ueit+epslon) return
c
c   synch radial nzons = number of baldur zones =mzones-lcentr
c      r=0 corresponds to synch index k = 0
c      synch zone bdy/center k --> baldur zone bdy/center k-lcentr
      do 11 k=1,nzons
      jz=k+lcentr
      jz1=k+lcentr-1
      surf(k)=avi(jz,5,1)*avi(jz,3,1)*uisl**2
      svol(k)=(avi(jz,12,1)-avi(jz1,12,1))*uisl**3
      stemp(k)=tes(2,jz1)*useh
      sdens(k)=rhoels(2,jz1)*used
  11  bavg(k)=sqrt(avi(jz1,9,1)**2*avi(jz1,10,1)+0.5*(avi(jz,7,1)+
     1 avi(jz1,7,1))*(r0ref*0.5*(bpoli(jz)*avi(jz,2,1)+
     2 bpoli(jz1)*avi(jz1,2,1)))**2)*uieb
c       use nzons=25 to ledge (r=wall) zones for baldur nzones=50
c     assume lcentr=2 (r=0);  nzons=nzones/2 + 1
c         synch zone centers index k --> BALDUR zone bdy j=2k+1
c        synch zone bdy index k --> BALDUR zone bdy j=2k+2
c     jz=2*k+1
c     jz2=2*k+2
c     surf(k)=avi(jz2,5,1)*avi(jz2,3,1)*uisl**2
c     svol(k)=(avi(jz2,12,1)-avi(jz2-2,12,1))*uisl**3
c     stemp(k)=tes(1,jz)*useh
c     sdens(k)=rhoels(1,jz)*used
c   ave b is toroidal only ??
c 11  bavg(k)=sqrt(avi(jz,8,1)**2*(0.5*(avi(jz,10,1)+
c    1 avi(jz-1,10,1)))+
c    2 avi(jz,7,1)*(r0ref*bpoli(jz)*avi(jz,2,1))**2)*uieb
c   Max and min B are from average
c
      bmin = bavg(1)
      bmax = bavg(1)
      do 12 j=2,nzons
        bmin = min ( bmin, bavg(j) )
        bmax = max ( bmax, bavg(j) )
 12   continue
c      bmin=aminaf(bavg,1,nzons)
c      bmax=amaxaf(bavg,1,nzons)
c
      call synch2 (nstep,nzesv,nzosv,nwsv,wstops)
c
c         wsyn = rad in keV/cm3;  wsyn(k) on j=2*jz+1
      do 30 k=1,nzons
  30  wsyn(k)=wsyn(k)*1.6022e-9
c   power density for BALDUR on zone centers
c      linear interpolation
      do 40 k=1,nzons
  40  wesyn(k+lcentr-1)=wsyn(k)
c     do 40 k=1,nzons-1
c     jz=2*k+1
c     jz2=2*k+2
c     zw=(wsyn(k+1)-wsyn(k))/(xbouni(jz+2)-xbouni(jz))
c     wesyn(jz)=wsyn(k)+zw*(xzoni(jz)-xbouni(jz))
c 40  wesyn(jz2)=wsyn(k)+zw*(xzoni(jz2)-xbouni(jz))
c     wesyn(mzones-1)=wsyn(nzons)
      wesyn(lcentr)=wsyn(1)
      wesyn(1)=wsyn(1)
      wesyn(mzones)=wsyn(nzons)
c   find radii where wsyn changes sign
c  47  do 48 j=1,n
c  48  rsyn(j)=0.
c      rsyn(1)=ra
c      j1=n
c      qs=qsyn(2)
c      k=1
c      do 50 j=3,n
c      if(qs*qsyn(j).gt.0.) go to 50
c      rsyn(k)=(r1(j-1)+qsyn(j-1)*dr1(j-1)/(qsyn(j-1)-qsyn(j)))
c      if(k.eq.1) j1=j-1
c      qs=qsyn(j)
c      k=k+1
c  50  continue
c  52  frj=(rsyn(1)-r1(j1))/dr1(j1)
c      rsyn1=rsyn(1)
c      tsync=vint(j1,rsyn1,frj,r1,dr1,qsyn)
      return
      end

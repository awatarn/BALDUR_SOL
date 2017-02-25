c     function zncnt(jy)
c     parameter(jx=55)
c     zncnt=55
c     return
c     end
c@syedit2  /11040/baldur/code/bald/dimprad.f
c
      subroutine syedit2
c
      include 'csynrd.m'
c
c     common/insink/ nzons,ree,roo,reo,roe,bmin,bmax,pwall,ptot,areaw,
c    2  taucrt,nw1,ndbug,iedit,
c    2  stemp(mj),sdens(mj),bavg(mj),svol(mj),surf(mj),wsyn(mj)
      write(6,100)
  100 format(1h1,1x,5h zone,13h       volume,13h  temperature,
     *13h      density,9x,4hbavg,13h    kev/cm**3,
     *13h       edot/e  )
      do 200 j=1,nzons
      ts1=wsyn(j)
      ts2=ts1/(1.5*sdens(j)*stemp(j))
      write(6,220) j,svol(j),stemp(j),sdens(j),bavg(j),ts1,ts2
  200 continue
  220 format(i7,6e13.2)
      write(6,300) ptot,pwall,nw1
  300 format(///19h  total volume loss ,e13.5,/
     * 19h total surface loss ,e13.5,//
     * 22h number of frequencies, i5 )
      return
      end

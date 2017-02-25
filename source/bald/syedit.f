c@syedit  /11040/baldur/code/bald/dimprad.f
      subroutine syedit
!cap
      include 'cparm.m'
      parameter ( kq=mj )
        include 'csynrd.m'
c
      write(6,100)
  100 format(1h1,1x,' zone       volume  temperature',
     & '      density',9x,'bavg','    kev/cm**3',
     & '       edot/e  ')
      do 200 j=1,nzons
      ts1=ploss(j)
      ts2=ts1/(1.5*deny(j)*teqp(j))
      write(6,220) j,voq(j),teqp(j),deny(j),bavg(j),ts1,ts2
  200 continue
  220 format(i7,6e13.2)
      write(6,300) ptot,pwall,nw1
  300 format(///'  total volume loss ',e13.5,/
     & ' total surface loss ',e13.5,//
     & ' number of frequencies', i5 )
      return
      end

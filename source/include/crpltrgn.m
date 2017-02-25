c@pltrgn
c
c   specifications for r plot generation
c
c  NUMBER OF x-AXES DEFINED IN xplot5.for
c
      integer t1imx, t1isw
c
      parameter (t1imx=50)
c
c  NUMBER OF LOGICAL OUTPUT GROUPS
c
      parameter (t1isw=50)
c
c
      common /rplgen/ lunpl1, lunpl2,
     >                nwr(2,2,t1imx), nset(t1isw)
      character*32 filpl1,  filpl2
c
      common /rp1gen/  filpl1, filpl2
c
c--------1---------2---------3---------4---------5---------6---------7-c

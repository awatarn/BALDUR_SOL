c--------1---------2---------3---------4---------5---------6---------7-c
c@ftrapfn  .../baldur/bald/ftrapfn.f
c  rgb 11-jul-01 changed ft.f to ftrapfn.f
c  rgb 11-jul-01 removed comtrp.m and inserted cbaldr.m
c  rgb 26-jul-88  removed cliches cl1 and clintf
c   dps 04-mar-87 inserted with bootstrap current cf. Mike Hughes
c
         function ftrapfn( kpt, ktyp, kopt)
c
c. return trapped particle fraction to baldur
c
C. KPT     refers to BALDUR mesh index
C. KTYP    refers to zone centre(2) or zone bound. (1)
C. KOPT    pointer to new (1) or old(0) calculation
C
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
C---------------------------------------------------------------------
CL              1.         New Calculation
C
         if ( kopt .eq. 0 ) go to 200
c
         ft = trapbr(kpt,ktyp)
         return
c
c-----------------------------------------------------------------------
CL              2.         Old Calculation
c
  200    continue
         zd = abs( ahalfs(kpt,ktyp) / rmids(kpt,ktyp) )
         ft=1.0-((1.0-zd)**2/(sqrt(1.0-zd*zd)*(1.0+1.46*sqrt(zd))))
c
         return
         end

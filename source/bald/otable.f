c--------1---------2---------3---------4---------5---------6---------7-c
c@otable  .../baldur/code/bald/default.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c rgb 20.26 09-jan-92 klabels type character format (a)
c       amck 23-jan-78 give error mess. if too many te ranges.
c       amck 4-aug-77 use z as key for element, read atomic weight
c       amck 15-mar-77 check if klen.le.0
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine otable(kelem,kximp,klen,khname,klabel,kunit,pweigt)
c
c
cl              set up tables of impurity fits
c
c
c---------------------------------------------------------------------
cl                  c8.1     comcor--coronal equilibrium fit coef
c
       common/comcor/
     r  ocoef , orange,
     i  mocoef, mxorng
       dimension
     r   ocoef(6,20,3),      orange(2,20,3)
c
c
        dimension       kelem(klen), irange(10,3),
     r  pweigt(klen),   zcperm(10,20),  zrperm(2,20),   zcoef(10)
c
cahk changed 1 to * in khname
        character *(*) khname(*), klabel
        character *10  ihname
        character *2   ihtype
c
        data    irp /40/,       icp /200/,      ira /30/
c
c------------------------------------------------------------------------------
c
c
c       calling format:
c
c       dimension lzelem(mimp),weight(mimp)
c       character lableo*60, lhspec(mchi)*10
c       data    lzelem /6,8,5,...etc./
c
c       ncoron = unit to read fits
c    call otable(lzelem,mximp,mimp,lhspec(limp1),labelo,ncoron,weight)
c       write (nprint,100) labelo, lhspec(ji+limp1-1)
c 100   format(4x,a60,(/10x,8helement ,a10))
c
c------------------------------------------------------------------------------
c
c
c       1)      initialize
c
c
  100   continue
        if (klen.le.0) return
c
c       call oclear
c
        mocoef = 6
        mxorng = 5*kximp
        if (klen.gt.0) mxorng = mxorng / klen
c       where  mxorng=5*(maximum allowed number of impurities)
c
        mxorbr = 5
c
        call reseti(irange,ira,0)
        irbrem = 0
c
        read (kunit,10100) klabel
c
c
        do 108 ji = 1, klen
          pweigt(ji) = 0.0
          khname(ji) = ' unknown  '
  108   continue
c
c
c
cl      2)      read fits
c
c
  200   continue
c
c               read element symbol, quantity, and name (if any)
c
        read (kunit,10200) ielem, ihtype, ihname, zweigt
c
c               check eof (blank card), name card (cols. 13-14 =na)
c
        if ( ielem .eq. 0 ) go to 300
        if ( ihtype .eq. 'na' ) go to 210
c
c               read one fit for one te range
c
        read (kunit,10201) zte0, zte1, (zcoef(j),j=1,mocoef)
c
c               check zte0 .lt. zte1
c
        if (zte1.ge.zte0) go to 206
        call mesage(' *** temperatures for temp. range for fit in ')
        call mesage('     wrong order ')
        call ivar(8helement ,ielem)
        call hvar(8hfit for ,ihtype)
        z = zte1
        zte1 = zte0
        zte0 = z
  206   continue
c
c
        itype = 0
        if ( ihtype .eq. 'lo' ) itype = 1
        if ( ihtype .eq. 'zb' ) itype = 2
        if ( ihtype .eq. 'zs' ) itype = 3
        if ( itype  .eq.  0 ) go to 200
c
c               stuff coefficients into appropriate part of ocoef
c
  210   continue
c
        if (klen.le.0) go to 200
        do 238 ji = 1, klen
        if ( ielem .ne. kelem(ji) ) go to 238
        if ( ihtype .eq. 'na' ) go to 230
c
  220   continue
        if (irange(ji,itype).lt.mxorng) go to 222
        call mesage(' *** too many te ranges. extra ranges ignored ')
        i = kelem(ji)
        call ivar(8helement ,i)
        call ivar(8hmxorng  ,mxorng)
        go to 238
c
  222   continue
        irange(ji,itype) = irange(ji,itype) + 1
        ir = irange(ji,itype) + (ji - 1)*mxorng
c
        orange(1,ir,itype) = zte0
        orange(2,ir,itype) = zte1
c
        do 226 jc = 1, mocoef
        ocoef(jc,ir,itype) = zcoef(jc)
  226   continue
c
        go to 238
c
c               set khname
c
  230   continue
c
        khname(ji) = ihname
c
        pweigt(ji) = zweigt
  238   continue
        go to 200
c
c
cl      3)      arrange fits in ascending order of temperature
c               in particular, lower temps. of fits for each
c               element and quantity will be in strict increasing
c               order.  if two fits have the same lower limit,
c               the second (in the deck) will be thrown away
c
c
  300   continue
        zinfin = 1.0e+34
        if (klen.le.0) go to 340
c
        do 338 jt = 1, 3
c
        do 318 ji = 1, klen
c
c               if only one range, don't permute
c
        if (irange(ji,jt).le.1) go to 318
c
        call resetr(zrperm,irp,0.0)
        call resetr(zcperm,icp,0.0)
c
c               each time through the do 308 loop,
c               the next lower limit above the last one
c               is found. the last one is in zte0.
c               the index of the next one is in ir2
c               after the do 304.   the correct ordered
c               set of fits is in zcperm and zrperm.
c
c
        zr1 = 0.0
        ir0 = (ji-1)*mxorng + 1
        ir1 = ir0 - 1 + irange(ji,jt)
c
        do 308 jr1 = ir0, ir1
c
c               find next lower limit (orange(1,x,jt)
c
        zr2 = zinfin
        ir2 = 0
c
        do 304 jr2 = ir0, ir1
        if(orange(1,jr2,jt).le.zr1.or.orange(1,jr2,jt).ge.zr2) go to 304
        ir2 = jr2
        zr2 = orange(1,jr2,jt)
  304   continue
c
        if (ir2.eq.0) go to 310
c
        zr1 = zr2
        zrperm(1,jr1) = orange(1,ir2,jt)
        zrperm(2,jr1) = orange(2,ir2,jt)
c
        do 306 jc = 1, mocoef
        zcperm(jc,jr1) = ocoef(jc,ir2,jt)
  306   continue
c
  308   continue
c
  310   continue
c
cl              copy correctly ordered fits back into ocoef, orange
c
        ir11 = ir0 + mxorng - 1
        do 314 jr = ir0, ir11
        orange(1,jr,jt) = zrperm(1,jr)
        orange(2,jr,jt) = zrperm(2,jr)
c
        do 312 jc = 1, mocoef
        ocoef(jc,jr,jt) = zcperm(jc,jr)
  312   continue
c
  314   continue
c
  318   continue
  338   continue
c
  340   continue
c
c
c
cl      4)      check for problems in fits
cl              check for: (1) no fits at all, (2) overlapping
cl              temperature ranges, (3) gaps between ranges
c
c
  400   continue
        zfuzz = 1.001
        if (klen.le.0) go to 460
c
cl      4.1)    impurity fits
c
        do 418 jt = 1, 3
c
        do 414 ji = 1, klen
c
        ir = (ji-1)*mxorng + 1
c
c               check for no fits at all
c
        if (orange(1,ir,jt).gt.0.0) go to 402
c
        if (jt.eq.1) call mesage(
     1          ' *** no loss term fit coefficients in deck for ')
        if (jt.eq.2) call mesage(
     1          ' *** no   z-bar   fit coefficients in deck for ')
        if (jt.eq.3) call mesage(
     1          ' *** no z-sq. bar fit coefficients in deck for ')
        call ivar(8helement ,kelem(ji))
        if (jt.eq.1) call mesage(
     1          '     0 used ')
        if (jt.gt.1) call mesage('     1 used ')
        go to 414
c
  402   continue
c
        ir2 = ir + mxorng - 2
        if (mxorng.le.1) go to 414
        do 410 jr = ir, ir2
c
        if (orange(1,jr+1,jt).le.0.0) go to 414
        if (orange(2,jr,jt).le.orange(1,jr+1,jt)*zfuzz) go to 406
c
        call mesage(' *** temp. range overlap in deck for fit for  ')
        if (jt.eq.1) call mesage('     radiation loss for ')
        if (jt.eq.2) call mesage('     z-bar  for ')
        if (jt.eq.3) call mesage('     z-sq.bar for ')
        call ivar(8helement ,kelem(ji))
        call rarray('range 1 ',orange(1,jr,jt),2)
        call rarray('range 2 ',orange(1,jr+1,jt),2)
        call mesage('     range 1 used ')
        go to 410
c
c               check missing range
c
  406   continue
c
        if (orange(2,jr,jt)*zfuzz.ge.orange(1,jr+1,jt)) go to 408
        call mesage(' *** missing temp. range in deck for fit for ')
        if (jt.eq.1) call mesage('     radiation loss for ')
        if (jt.eq.2) call mesage('     z-bar  for ')
        if (jt.eq.3) call mesage('     z-sq.bar for ')
        call ivar('element ',kelem(ji))
        call rarray('range   ',orange(2,jr,jt),2)
        call mesage('     range will be interpolated ')
c
  408   continue
c
  410   continue
c
  414   continue
  418   continue
c
  440   continue
c
c
cl              check atomic weight
c
c
        do 444 ji = 1, klen
        if (pweigt(ji).gt.0.0) go to 444
c
        call mesage(' *** no atomic weight in deck for ')
        call ivar('atom.no.',kelem(ji))
        call mesage('     1.0 used ')
        pweigt(ji) = 1.0
  444   continue
c
  460   continue
c
        return
c
c
c
cl      100)    format statements
c
c
10100   format(a)
10200   format(i3,9x,a2,10x,a10,10x,f10.5)
10201   format(2e10.3,3(2x,e13.6)/3(2x,e13.6))
        end

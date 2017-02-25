c--------1---------2---------3---------4---------5---------6---------7-c
c@preset  .../baldur/code/bald/preset.f
c  ton 13-aug-01 add lbound and cbound
c  rgb 08-jan-01 zeroed out all lthery, cthery, and theory coeffs
c  rgb 12-jun-00 initialized arrays rlepwr(jt) and rlipwr(jt)
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  pis 02-jul-98 added initialization of vrota array
c  pis 11-jun-98 initialized xwexba, wexba, twexba arrays
c    and nxwexba, ntwexba scalars
c  alh 30-oct-97 make tcold & tcoldp time dependent arrays
c  rgb 14-aug-96 round ---> rndeps, fuzz removed
c  jek 08-dec-95 defaulted cthery(123) = 1 for coeff of k_parallel
c  rgb 04-jan-95 fixed d0nmon(jt) = 0.0 index
c  rgb 16-dec-94 initialize cimprd(j) = 1.0, j=10,19 as coefficients for the
c    sources and sinks in the non-equilibrium impurity radiation model
c  rgb 21-apr-94 initialize cimprd(2) = 1.0 for zeff control
c  rgb 17-apr-94 initialize ftzeff(jt) = 1.0, jt=1,20
c  rgb 26-dec-93 added tcoef(jt) and coeft(jt,jc) general purpose coefs
c    removed lstart(j) and cstart(j)
c  rgb 14-jul-92 change lpelet to lneocl, change cpelet to cneocl
c  rgb 14:00 14-feb-92 set fmh(j)=0.0 and increase j=1,3 to j=1,5
c  rgb 20.11 26-aug-91 update cthery defaults
c  rgb 06-may-90 18.34 default cpvelc and cpvion to 1.5
c  rgb 12-feb-90 initialize sfutz(js) to 0.0, js=1,20
c  rgb 06-feb-90 d0nmon(j)=0.0, j=1,20 and gainv=0.0
c  rgb 03-feb-90 cthery(30)=1.0 and cthery(31)=1.9
c  rgb 29-jan-90 new defaults for cthery(20-22)
c  rgb 14-jan-90 replaced sfutz(j) with cthery(j), defaulted lthery(j)
c  rgb 15-dec-89  default tmod=0., snestr=0., snebar=0., used in dempirc
c  rgb 27-oct-89  default sfutz(19) to -1.5
c  rgb 09-oct-89  default sfutz(17) and sfutz(18) to 1.0
c  ces 27-jul-89 sfutz(12--20)=0.0
c  ces 24-jul-89 set sfutz(11)=1.0 for eta-i-crit
c  dps 17-oct-88 15.06 Make the following change from Singer:
c  elg 11-oct-88 set the default values for subroutine theory's variables
c  dps 14.02 20-jun-88 set default of cemprc(15)=1.5.
c  rgb 14.07 03-apr-88 initialized cemprc(17) to 0.15
c      cemprc(17) is used as breakpoint in Tang model in sbrtn EMPIRC
c  rgb 12.60 25-aug-87 initialized HPOWMW to 0.0
c  rgb 12.55 12-aug-87 initialized HTON(j)=epsinv, HTOFF(j)=-epsinv
c  rgb 12.54 01-aug-87 initialized bdhyde, bdimpe, bdtee, bdtie, bdtime
c  rgb 17-jul-87  change default of cfutz(139) from 0. to 1.
c       This turns on off-diag neocl diff for backward compatability
c  rgb 01-jul-87 default for control switches
c  wos 08-may-86 default for equilibrium:  1D BALDUR
c  rgb 5-mar-86 preset cemprc, cgbeta, cgpowr, gxp, gcomb for empirc
c       fgps 9-feb-83 modified to handle more than 2 impurity species.
c       aes 19-apr-82 change default of cfutz(198) to 0.0
c       aes 03-mar-82 default nlgpin=.t. --> proton ioniz. in monte *on*
c       aes 02-mar-82 set cfutz(198)=.8 --> default tcold readjustment
c       aes 16-feb-82 set afslow=1.0
c       aes 11-feb-82 set hfutz(5)=1.0
c       aes 07-dec-81 cfutz(351)=1.e10 --> default plt bohm limit off
c       aes 04-dec-81 default npelou = 1
c       aes 03-dec-81 nlpomt(8) = .true.
c       aes 26-oct-81 cfutz(51)=1.e10 --> default soft bohm limit off
c       aes 12-aug-81 nlpomt(7) = .true.
c       aes 7-aug-81 gwmin = 1.e-6, nhskip = 2
c       aes 2-aug-81 nlpomt(20) = .true.
c       aes 2-aug-81 nlpomt(5),nlpomt(6) = .true.
c       aes 17-feb-81 gxdome --> gxdome(2)
c       aes 16-feb-81 ngxene --> 31; gxmaxe --> 100 kev
c       aes 6-nov-80 changed ngxene to 25; elminated ngxmod
c       aes 4-nov-80 updated comments near line 180
c       aes 4-nov-80 added rlfreq,rlipwr,rlepwr,nrlmod,nrldmp
c       aes 28-oct-80 set ngxmod = 2
c       aes 27-oct-80 set gxthet, gxphi
c       aes 22-oct-80 added neutral outflux spectrum variables
c       aes 8-oct-80 set ntgraf=4,nrgraf=3
c       dhfz 26-apr-79 set cfutz's to 1 x neoclassical + 1 x ware pinch
c       drm  19-april-79 added ntype, rcwals, and rdwals
c       dhfz 5-apr-79 npel=1,npel2=1
c       dhfz 5-apr-79 set ntychl=5
c       dhfz 5-apr-79 add gftime and gflmax
c       dhfz 2-sept-78 set fracti(1) = 0.0
c       dhfz 2-sept-78 set nimp(1) = 0
c       dhfz 20-july-78 change nnhfit to eehfit,nnifit to eeifit,
c               nnfit to eefit, ntefit to eeteft, ntifit to eetift,
c               and nbfit to eebfit
c       dhfz 20-july-78 add ehfit,eifit,efit,etefit,etifit,
c               and ebfit
c       amck 18-may-78 nglmon -> nlgmon
c       amck 28-mar-78 turn off beam current page (nlpomt(16))
c       amck 22-mar-78 set hfutz's
c       amck 31-jan-78 turn off 1st hprint page (nlpomt(13))
c       amck 23-jan-78 set tspare
c       amck 19-jan-78 change extz -> extf
c       amck 11-jan-78 set nitmax, nlrcom, extz
c       amck 28-sep-77 set xfutz(6)
c       amck 16-sep-77 set haper, haperv, nhaper
c       amck 7-sep-77 nsedit = 1
c       amck 25-aug-77 set niprof, nizone, nipart
c       amck 25-aug-77 set ellipt, grecyc
c       amck 22-aug-77 set hfocl(v), hfract
c       amck 11-aug-77 oxygen is now no. 8
c       amck 9-aug-77 set neutral beam stuff, mxhbem, mxhfr
c       amck 29-jul-77 set nhfit, nifit
c       amck 15-jul-77 add xfutz(1-5)
c       amck 25-may-77 add delmax
c       amck 21-apr-77 change dengas -> fracth, etc.
c       amck 14-apr-77 all bohm cfutz's 0 except ke
c       amck 23-mar-77 make nrgraf not conflict with nounit
c       amck 15-mar-77 add nounit
c       amck 16-feb-77 55 zone, 2hyd, 2imp, 20t
c       amck 12-feb-77 105 zone, 2 hyd, 2 imp
c       amck 11-feb-77 45 zone, 1 hyd, 1 imp, set mx-- vars.
c       amck 29-jan-77 more aurora defaults
c       amck 27-jan-77 add defaults for aurora (monte-carlo neutrals)
c       amck 18-jan-77 modify defaults to usual choices
c       amck 10-jan-77 default plotting flags
c       amck 2-nov-76 change cfutz for bohm to 0.03
c       amck 26-oct-76 set some beam vars.
c       amck 12-oct-76 change nlsorb, nlsorc -> nlsorc, nlsord
c       amck 5-oct-76 set cfutz
c       amck 27-sep-76 change dtinit, dtmin
c       amck 22-sep-76 add tmax
c       amck 17-sep-76 add nbound=0
c       amck 7-sep-76 change nedit to "10"
c       amck 3-sep-76 add natomc, eioniz
c       amck 26-aug-76 make rmajor, rminor arrays
c       amck 17-aug-76 add presets of new namelist vars.
cdoc
c******************************************************************************
c
c       ------------
c       sbrtn PRESET   file DEFAULT
c       ------------
c
c       1.3     set default input values
c
cend
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine preset

c
        include 'cparm.m'
      include 'cbaldr.m'
        include 'cbparm.m'
c
        data    iclass /1/,     isub /3/
c
c
c------------------------------------------------------------------------------
c
c
        cfev = 1.60210
        cxev = -12.0
        cfevdk = 1.16049
        cxevdk = 4.0
        fca0 = 5.29167
        fxa0 = -9.0
        fcae = 5.48597
        fxae = -4.0
        fcalfa = 7.29720
        fxalfa = -3.0
        fcan = 1.0086654
        fcap = 1.00727663
        fcc = 2.997925
        fxc = 10.0
        fcc1 = 3.7405
        fxc1 = -5.0
        fcc2 = 1.43879
        fxc2 = 0.0
        fce = 1.60210
        fxe = -20.0
        fces = 4.80298
        fxes = -10.0
        fcf = 9.64870
        fxf = 3.0
        fcg = 6.670
        fxg = -8.0
        fch = 6.6256
        fxh = -27.0
        fck = 1.38054
        fxk = -16.0
        fclnae = 2.42621
        fxlnae = -10.0
        fclnap = 1.32140
        fxlnap = -13.0
        fcme = 9.1091
        fxme = -28.0
        fcmn = 1.67482
        fxmn = -24.0
        fcmp = 1.67252
        fxnucl = -24.0
        fcna = 6.02252
        fxna = 23.0
        fcpi = 3.1415926536
        fcr = 8.3143
        fxr = 7.0
        fcre = 2.81777
        fxre = -13.0
        fcrinf = 1.0973731
        fxrinf = 5.0
        fcsb = 5.6697
        fxsb = -5.0
        fcsigt = 6.6516
        fxsigt = -25.0
        fcv0 = 2.24136
        fxv0 = 4.0
        fcwien = 2.8978
        fxwien = -1.0
        epsinv = 1.0e34
        epslon = 1.0e-34
        ninfin = 1000000000
        rndeps = 1.5e-08
c
c------------------------------------------------------------------------------
c
c
c
        if (.not.nlomt1(isub)) go to 10
        call mesage(' *** 1.3 subroutine preset bypassed *** ')
        return
   10   continue
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       5.1--comin
c
c------------------------------------------------------------------------------
c
c
c       since all parameters to the code are entered through
c       comin, combas, and comddp, it is only necessary to preset these.
c       combas is already preset in modify.
c
c
cl      0)      output
c
        ntychl=5
c
cl      1)      maximum dimensions
c
c
        mxzone = mj
        mxhyd = 2
        mximp = idximp
c       where "mximp" denotes the maximum number of impurity species.
        mxions=mxhyd+mximp
        mxt = 20
        mxchi = 2 + mxhyd + mximp
        mxhbem = 12
        mxhfr = 3
c
cl      2)      magnetic fields,        in kilo-gauss
c
        bz = 35.0
c
cl      3)      current, in kilo-amperes
c
        curent(1) = 600.0
c
cl      4)      densities, in part/cu cm
c
        dens0 = 5.0e+13
        dens1 = 1.0e+13
c
cl      5)      partial densities, on any scale
c
c
        fracth(1) = 1.0
        fracti(1) = 0.0
c
cl      6)      ion species
c
        ngas(1) = 1
        nimp(1) = 0
c
cl      7)      delta t, in seconds
c
        dtinit = 1.0e-06
        dtmax = 2.0e-03
        dtmin = 1.0e-09
c
cl      8)      fitting exponents or flags
c               (see auxval for use)
c
c       0 <= eebfit
c
        eebfit = 0.
c
c       1 <= nrfit <= 2
c
        nrfit = 1
c
c       1 <= eefit, eeteft, eetift
c
        eefit = 3.
        call resetr(eehfit,mxhyd,eefit)
        call resetr(eeifit,mximp,eefit)
        eeteft = 2.
        eetift = 2.
c
c       1 <= ehfit,eifit,efit,etefit,etifit,ebfit
c
        efit = 1.
        call resetr(ehfit,mxhyd,efit)
        call resetr(eifit,mximp,efit)
        etefit = 1.
        etifit = 1.
        ebfit = 1.
c
cl      9)      output flags
c
        nedit = 10
        nskip = 1
        sedit = 0.00
        nsedit = 1
        nlpomt(5) = .true.
        nlpomt(6) = .true.
        nlpomt(7) = .true.
        nlpomt(8) = .true.
        nlpomt(13) = .true.
        nlpomt(16) = .true.
        nlpomt(20) = .false.
c
        nplot = 10
        splot = 0.00
        nrgraf = 3
        ntgraf = 4
c
cl      10)     code flags
c
        nitmax = 25
        nzones = 20
        theta = 1.0
        tmax  = 1.0
        nrun = 0
        tspare = 60.0
c
        nlextr = .true.
        errmax = 0.01
        delmax = 0.1
        xfutz(1) = 0.5
        xfutz(2) = 1.5
        xfutz(3) = 2.0
        xfutz(4) = 2.0
        xfutz(5) = 0.333
        xfutz(6) = 0.333
c
        nldiff = .true.
        nlsorc = .true.
        nlsord = .true.
c
        nlrcom = .true.
c
cl      11)     te and ti, in kev
c
        te0 = 0.2
        ti0 = 0.2
        te1 = 0.01
        ti1 = 0.01
c
cl      12)     major and minor radii, in cm.
c
        rmajor(1) = 130.0
        rminor(1) = 40.0
        ellipt = 1.0
c
cl      13)     diffusion coefficients
c
        cfutz(1) = .0
        cfutz(2) = .0
        cfutz(3) = .0
        cfutz(4) = .0
        cfutz(5)=.0
        cfutz(6)=.0
        cfutz(7)=.0
        cfutz(8)=.0
        cfutz(9)=1.
        cfutz(10)=1.
        cfutz(11)=1.
        cfutz(12)=1.
c
        cfutz(19)=1.
c
        cfutz(139) = 1.0
c
        cfutz(51)=1.e10
        cfutz(351)=1.e10
c
c
        cpvelc = 1.5
        cpvion = 1.5
c
        extf   = 1.0
        extzef = 1.0
c
        ntrans = 1
c
cl      14)     impurities:
c
        nounit = 22
        natomc = 2
c
cl      15)     neutrals:
c
        eioniz = 0.040
        grecyc = 1.0
        gwmin  = 1.e-6
        ngpart = 500
        ngprof = 10
        ngsplt = 5
        ngzone = 20
        nlgmon = .true.
        nlgpin = .true.
        ngxene = 31
        gxmaxe = 1.e5
        gxmine = 1.
        gxdome(1)=4.
        gxdome(2)=10.
        gxphi(2) =15.
        gxphi(3) =30.
        gxphi(4) =60.
        gxthet(2)=15.
        gxthet(3)=30.
        gxthet(4)=60.
c
        do jj = 1, 20
           tcold(jj) = -1.0
           tcoldp(jj) = -1.0
        enddo
        tcold(1)=0.003
        tcoldp(1)=0.003
c
cl      16)     boundary conditions:
c
        nbound = 1
c
cl      17)     neutral beam variables
c
        nhsrc = 1
        nhmu = 10
        nhe = 10
        nipart = 10000
        niprof = 10
        nizone = 20
        hnchek = 0.05
        htchek = 0.20
        nhskip = 2
c
        call reseti(nhshap,mxhbem,1)
        call reseti(nhprof,mxhbem,1)
        call reseti(nhprfv,mxhbem,1)
        call reseti(nhaper,mxhbem,1)
c
        call resetr(hton ,mxhbem, epsinv)
        call resetr(htoff,mxhbem,-epsinv)
        call resetr(hr,mxzone,1.0)
        call resetr(hfocl ,mxhbem,1.0e+10)
        call resetr(hfoclv,mxhbem,1.0e+10)
        call resetr(haper ,mxhbem,1.0e+07)
        call resetr(haperv,mxhbem,1.0e+07)
c
        call resetr(hpowmw,mxhbem,0.0)
c
        do 170 jb = 1, mxhbem
        hfract(1,jb) = 1.0
  170   continue
c
        hfutz(1) = 1.0
        hfutz(2) = 1.0
        hfutz(5) = 1.0
c
cl      18)     miscellaneous
c
c
c               gas puffing
                do jt=1,20
                  gftime(jt)=1.e34
                  gflmax(jt)=1.e16
                  d0nmon(jt) = 0.0
                enddo
c
                   gainv = 0.0
c
c               lower hybrid heating
c
c  Note:  the default values of rlepwr(1) and rlipwr(1) are the same
c    as the scalars rlepwr and rlipwr used to be set before 12-jun-00
c    The rest of the arrays rlepwr(jt) and rlipwr(jt) are set to a
c    negative number to indicate that they need to be reset later
c    in sbrtn auxval.
c
                nrldmp = 3
                nrlmod = 3
                rlepwr(1) = 1.0e+06
                rlipwr(1) = 1.0e+06
                do jt=2,20
                  rlepwr(jt) = -99.0
                  rlipwr(jt) = -99.0
                enddo
                rlfreq = 5.0e+08
c
c               pellets
                npel=1
                npel2=1
                npelou=1
c
c..alpha particle treatment
c
                ntype=2
                afslow=1.
                rcwals=0.
                rdwals=0.
c
c..smoothing parameters for sbrtn xscale
c
                smrlow=0.
                smlwcy=0.
                lsmord=3
c
cl      19)     for sbrtn empirc, global transport, bateman, 6-mar-86
c
      do 190 j=1,mxt
      cgbeta(j) = 1.0
      cgpowr(j) = 1.0
 190  continue
c
      do 192 jn=1,6
      do 192 jl=1,6
      gcoefs(jl,jn) = 1.0
 192  continue
c
      do 194 ji=1,30
      do 194 jn=1,6
      do 194 jl=1,6
      gxp(jl,jn,ji) = 0.
 194  continue
c
      do 196 jn=1,6
      do 196 jm=1,6
      do 196 jl=1,6
      gcomb(jl,jm,jn) = 0.0
      gcomb(jl,jn,jn) = 1.0
 196  continue
c
      tmod   = 0.
      snestr = 0.
      snebar = 0.
c
      nrad = 0  ! used to control armins and armajs in sbrtn balnew
c
c..default 'equilibrium': concentric circ.cylinder --> 1D BALDUR
c
         leqtyp= 0
c
c..control switches in common / cswtch /
c
      do 200 j=1,32
        lnumer(j) = 0
        cnumer(j) = 0.0
c
        linout(j) = 0
        cinout(j) = 0.0
c
        ltrnsp(j) = 0
        ctrnsp(j) = 0.0
c
        lemprc(j) = 0
        cemprc(j) = 0.0
c
        lbeams(j) = 0
        cbeams(j) = 0.0
c
        ldivrt(j) = 0
        cdivrt(j) = 0.0
c
        lauxil(j) = 0
        cauxil(j) = 0.0
c
        limprd(j) = 0
        cimprd(j) = 0.0
c
        lnugas(j) = 0
        cnugas(j) = 0.0
c
        lneocl(j) = 0
        cneocl(j) = 0.0
c
        lfsion(j) = 0
        cfsion(j) = 0.0
c
        lstabl(j) = 0
        cstabl(j) = 0.0
c
        lbound(j) = 0
        cbound(j) = 0.0
c
 200  continue
c
      cimprd(2) = 1.0
      cimprd(4) = 1.0
c
      do jc=10,19
        cimprd(jc) = 1.0
      enddo
c
      do 203 js=1,20
        sfutz(js) = 0.0
 203  continue
c
c  15.06 for subroutine theory
c
      do 201 jc=1,150
        cthery(jc) = 0.0
 201  continue
c
      do 202 jl=1,50
        lthery(jl) = 0
 202  continue
c
cbate      cthery(1)=0.0
cbate      cthery(2)=1.0
cbate      cthery(3)=0.5
cbate      cthery(4)=0.0
cbate      cthery(5)=6.0
cbate      cthery(6)=1.0
cbate      cthery(7)=1.0
cbate      cthery(8)=3.0
cbate      cthery(9)=6.0
cbate      cthery(10)=0.0
cbate      cthery(11)=1.0
cbate      cthery(12)=-4.0
cbate      cthery(13)=-4.0
cbate      cthery(14)=-4.0
cbate      cthery(15)=-4.0
cbate      cthery(16)=0.0
cbate      cthery(17)=0.0
cbate      cthery(18)=0.0
cbate      cthery(19) = -1.5
cbate      cthery(20) = 1.0
cbate      cthery(21) = 1.0
cbate      cthery(22) = 0.25
cbate      cthery(23) = 1.0
c
cbate      cthery(30) = 1.0
cbate      cthery(31) = 5.0
cbate      cthery(34) = 0.25
cbate      cthery(35) = 5.0
cbate      cthery(36) = 4.0
c
cbate      cthery(41) = 3.0
cbate      cthery(42) = 1.0
cbate      cthery(43) = 0.16667
cbate      cthery(44) = 1.0
cbate      cthery(45) = 1.0
cbate      cthery(47) = 8.0
cbate      cthery(48) = 0.0
cbate      cthery(49) = 1.0
cbate      cthery(50) = 10.0
cbate      cthery(51) = 10.0
cbate      cthery(52) = 10.0
cbate      cthery(53) = 10.0
cbate      cthery(54) = 10.0
cbate      cthery(60) = 10.0
cbate      cthery(70) = 1.0
cbate      cthery(71) = -4.0
cbate      cthery(72) = 0.0
cbate      cthery(73) = 1.0
cbate      cthery(74) = 0.0
cbate      cthery(75) = 1.0
cbate      cthery(76) = 0.16667
cbate      cthery(77) = 10.0
cbate      cthery(80) = 1.0
cbate      cthery(81) = -4.0
cbate      cthery(82) = 0.0
cbate      cthery(83) = -4.0
cbate      cthery(85) = 2.0
cbate      cthery(86) = 0.15
cbate      cthery(111) = 0.0
cbate      cthery(112) = 0.0
cbate      cthery(113) = 0.0
cbate      cthery(114) = 0.0
cbate      cthery(119) = 1.0
cbate      cthery(121) = 1.0
cbate      cthery(123) = 1.0
cbate      cthery(124) = 1.0
cbate      cthery(135) = 2.0
cbate      cthery(136) = 1.0
c	   cthery(137) = 1.0
c
c   set the contributions to the fluxes and interchange
c
        fdrint=0.0
      do 220 j=1,5
        fdr(j)=0.0
        fig(j)=0.0
        fec(j)=0.0
	feg(j)=0.0
        fti(j)=0.0
        frm(j)=0.0
        fkb(j)=0.0
        frb(j)=0.0
        fhf(j)=0.0
        fmh(j)=0.0
 220  continue
cap { 
cbate      feg(3) = 1.0
cbate      fig = 0.8
cbate      fkb = 0.01
cbate      frb = 1.0
cap }       
c..14.02: Default value for critical eta-i in Tang model; same as used
c..in calibration by M. Redi.
c
        cemprc(15)=1.5
      cemprc(17) = 0.15
c
c..Time-dependent boundary conditions
c  bdhyde(20,2), bdimpe(20,4), bdtee(20), bdtie(20), bdtime(20)
c
        call resetr (bdhyde,40,0.0)
        call resetr (bdimpe,80,0.0)
        call resetr (bdtee ,20,0.0)
        call resetr (bdtie ,20,0.0)
        call resetr (bdtime,20,0.0)
        call resetr (ftzeff,20,1.0)
c
c  les  nov-90  d3he new namelist
c
      ekappa = 1.
      edelta = 0.
      ash4 = 0.
      ashp = 0.
      fploss=1.
      f4loss=1.
      ftloss=1.
      f3loss=1.
      cfutz(490) = 0.
      cfutz(499) = 0.
c   cmg thermal and particle transport
      alcmg = 0.25
      gcmg=0.
      epscmg = 0.2
      alpcmg = 0.7
      tohcmg=1000.
      lcmg = 0
c   synchrotron radiation (tamor)
      roo1=0.90
      ree1=0.90
      roe1=0.05
      reo1=0.05
c         tamor uses tacrt=0.7
      tacrt=0.6
      areaw1=0.
c         grid spacing nzons hardwired
      nzons1=nzones     !for nzones=50 mzones to wall = 52
      iedit1=0
      ndbug1=0
c   analytic auxiliary heatin and particle source
      lpaux=0
      lspaux=0
      lrepld=0
      lreplt=0
      lrepl3=0
      do 102 j=1,20
      atim(j)=0.0
      pauxe(j)=0.0
      pauxi(j)=0.0
      stim(j)=0.0
      spd(j)=0.0
      spt(j)=0.0
      spp(j)=0.0
      sp3(j)=0.0
      sp4(j)=0.0
 102  spimp4(j)=0.0
      apaux=2.
      spaux=3.
      apauxe=2.
      apauxi=2.
      do 104 j=1,6
 104  spauxi(j)=3.0
c
c..general purpose time-dependent coefficients
c
      do jt=1,20
        tcoef(jt) = 0.0
        do jc=1,20
          coeft(jt,jc) = 1.0
        enddo
      enddo
c pis: 29-may-1998
c preset values associated with wexb shearing rate in etaw17split
!cap
      xwexba = 0.0
      wexba  = 0.0
      vrota  = 0.0
      nxwexba = 0
      ntwexba = 0
      twexba = 0.0
c
      return
      end

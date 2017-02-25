c/ 21:00 21-feb-96 /011040/.../baldur/code/bald/dio.f  Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c  To obtain this file, type           (use appropriate date for yymmdd)
c cfs get /11040/byymmdd/baldyy/wcode.tar
c tar xf wcode.tar  (this creates directory .../code and subdirectories)
c cd code/bald      (this changes to subdirectory .../code/bald)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c  file dio contains input output routines for baldur.  Subrtns:
c
c  stripx  - to strip ! off input data sets
c  data   - input namelist data (not including equilibrium data)
c  redata - reread input data at time REREAD or cfutz(82)
c  output - control output (moved to file DOUTPUT)
c  sprint - generate short printout
c  aprint - alpha particle printout
c  hprint - print out beam data
c  iprint - print out injector data
c  grafix - output to intermediate file for04 for05 for graphics
c  fprint - fusion output
c  gprint - neutral gas printout
c  ncprnt - prints results of impurity charge state transport
c  nclwpr - local radiation printout
c  nciwpr - integrated radiation printout
c  mprint - long printout, many pageszte
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@stripx   /baldur/code/bald  file DIO
c  rgb 06-oct-94 changed line length from 80 to 132 characters
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@data   /baldur/code/bald  file DIO
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  pis 02-jul-98 added vrota to nurun2 namelist
c  pis  8-jun-98 compressed nurun2 namelist def to 15 continuation lines
c  pis 15-may-98 added wexba, xwexba, twexba to nurun2 namelist
c  alh 16-jan-98 interpolate value for tcold
c  rgb 14-aug-96 added rndeps to first namelist
c  rgb 17-apr-94 added ftzeff
c  rgb 26-dec-93 added tcoef(jt) and coeft(jt,jc) general purpose coefs
c  rgb 16-oct-93 added fec
c  rgb 14-jul-92 change lpelet to lneocl, change cpelet to cneocl
c  rgb 03-jun-92 use linetest to bypass nurun3 as needed
c  les    nov-90 added namelist nurun3 for d-3he fusion
c  rgb 17-mar-92 inserted fmh and rearranged namelist nurun1
c  rgb 11-feb-92 changed nurun3 append D-3He data to 2nd namelist nurun2
c  les    nov-90 added namelist nurun3 for d-3he fusion
c  rgb 06-feb-90 added d0nmon and gainv to /nurun2/
c  rgb 01-feb-90 added tmod, snestr, snebar to /nurun1/
c  rgb 14-jan-90 18.00 added lthery and cthery
c  dps 17-oct-88 15.06 Add following change from Singer:
c  elg 11-oct-88 add subroutine theory's parameters to namelist nrun1
c  rgb 02-jul-87  Added control switches LNUMER,..., CSTABL to 1st namelist
c       ces 27-jun-86  nlpomt(20) true skips nurun1, nurun2 printout
c
         subroutine data
c
c 1.4  define data specific to run
c
c**********************************************************************c
c**********************************************************************c
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'clsaw.m'
c---------------------------------------------------------------------
       namelist/nurun1/ tmod, snestr, snebar,lstabl,cstabl, rndeps,      !cpis
     .  nledge, nlres , nonlin, nout  , nprint, nread , nrec  , nresum,
     .  ndiary, nin   , npunch, nrun  , nadump, npdump, npoint,
     .  nsub  , nvdump, nlched, nlhead, nlomt1, nlomt2, nlomt3, nlrept,
     .  ellipt, rehe  , rehi  , reoff , reon  , rex1  , rex2  , rey1  ,
     .  rey2  , rle0  , rle1  , rleprf, rlepwr, rlfreq, rli0  , rli1  ,
     .  rliprf, rlipwr, rlpara, rlpowr, nrldmp, nrleex, nrliex, nrlmod,
     .  rlhepwr, rlheprf,
     .  nrlpar, cfutz , dfutzd, dfutze, dfutzi, dfutzv, dxsemi, extf  ,
     .  extzef, vxsemi, xesemi, xisemi, gnchek, grecyc, gtchek, gwmin ,
     .  gxdome, gxmaxe, gxmine, gxphi , gxthet, ngpart, ngprof, ngsplt,
     .  ngxene, ngzone, nlglim, nlgmon, nlgpin, nlgref, nlgspt, delmax,
     .  errmax, theta , thetap, tspare, natomc, nbound, nfusn , nitmax,
     .  nounit, ntrans, nldiff, nlextr, nlrcom, nlsorc, nlsord, xfutz ,
     .  nbdump, nplot , nrgraf, nsedit, ntgraf, nlpomt, reread, versno
     & ,cemprc, cgbeta, cgpowr, gcoefs, gxp   , gcomb , nrad  , leqtyp
     & ,smrlow, smlwcy, lsmord, sfutz
     &  , lnumer, cnumer, linout, cinout, tcoef, coeft, ltrnsp, ctrnsp
     &  , lemprc, lbeams, cbeams, ldivrt, cdivrt, lauxil, cauxil
     &  , limprd, cimprd, lnugas, cnugas, lneocl, cneocl, lfsion, cfsion
     &  ,fdr,fig,fec,feg,fti,fdrint,frm,fkb,frb,fmh,fhf,lthery,cthery
     &  , lbound,cbound
c---------------------------------------------------------------------
       namelist/nurun2/ npelgc,frout,frcor,rpelc,ypa,hpowmw,qmhd1,qmhd2
     & , qstar,rrstar,cjbeam,habeam,hangle,haper,haperv,hdiv,hdivv
     & , hebeam,height,hfocl,hfoclv,hfract,hfutz,hibeam,hlenth,hnchek
     & , hpbeam,hprof,hprofv,hr,hrmaj,hrmin,htchek,htoff,hton,hwidth
     & , nhaper,nhbeam,nhe,nhmu,nhprfv,nhprof,nhshap,nhskip,nhsrc 
     & , nipart,niprof,nizone,afslow,rcwals,rdwals,ntype,bpoid,bz,bzstar
     & , cpvelc,cpvion,curent,denga0,denga1,dengas,denim0,denim1,denimp
     & , denmon,dens,dens0,dens1,dtinit,dtmax,dtmin,ebfit,eebfit,eefit
     & , eehfit,eeifit,eeteft,eetift,efit,ehfit,eifit,eioniz,elecd0
     & , electe,etefit,etifit,flgas,flimp,fracth,fracti,gflmax,gfract
     & , gftime,apresr,radius,rastar,rmajor,rminor,rpa,rpela,rqs,sedit
     & , splot,tbpoid,tcold,tcoldp,tcomp,te,te0,te1,tedit,tgas,ti,ti0
     & , ti1,timp,tinit,tmax,tpela,tplot,vpela,wtgas,wtimp,nbfit,nedit
     & , ngas,nimp,nnfit,nnhfit,nnifit,npel,npel2,npelga,npelou,npuff
     & , nrfit,nskip,ntefit,ntifit,ntty,nzones,bdhyde,bdimpe,bdtee
     & , bdtie ,bdtime,d0nmon,gainv,ftzeff,wexba,vrota,xwexba,twexba
     & , cfz230,itvpel
c---------------------------------------------------------------------
       namelist/nurun3/lcmg,alcmg,epscmg,alpcmg,gcmg,tohcmg,
     .   ash4,ashp,
     .   roo1,ree1,roe1,reo1,tacrt,areaw1,ndbug1,iedit1,
     .   fploss,f4loss,f3loss,ftloss,
     .   ekappa,edelta,
     .   lpaux,pauxe,pauxi,atim,apaux,apauxe,apauxi,
     .   lspaux,lrepld,lreplt,lrepl3,spd,spt,spp,sp3,sp4,spimp4,
     .   stim,spaux,spa
c---------------------------------------------------------------------
       namelist/reset/  npelgc, frout , frcor , rpelc ,
     .  nledge, nlres , nonlin, nout  , nprint, nread , nrec  , nresum,
     .  ndiary, nin   , npunch, nrun  , nadump, npdump, npoint,
     .  nsub  , nvdump, nlched, nlhead, nlomt1, nlomt2, nlomt3, nlrept,
     .  cfutz , dfutzd, dfutze, dfutzi, dfutzv, dxsemi, extf  , extzef,
     .  vxsemi, xesemi, xisemi, gnchek, grecyc, gtchek, gwmin , gxdome,
     .  gxmaxe, gxmine, gxphi , gxthet, ngpart, ngprof, ngsplt, ngxene,
     .  ngzone, nlglim, nlgmon, nlgpin, nlgref, nlgspt, delmax, errmax,
     .  theta , thetap, tspare, natomc, nbound, nitmax, nounit, ntrans,
     .  nldiff, nlextr, nlrcom, nlsorc, nlsord, xfutz , nbdump, nplot ,
     .  nrgraf, nsedit, ntgraf, nlpomt, cjbeam, hangle, haper , haperv,
     .  hdiv  , hdivv , hebeam, height, hfocl , hfoclv, hfract, hfutz ,
     .  hibeam, hlenth, hnchek, hpbeam, hprof , hprofv, hr    , hrmaj ,
     .  hrmin , htchek, htoff , hton  , hwidth, nhaper, nhprfv, nhprof,
     .  nhshap, nhskip, bpoid , bz    , bzstar, cpvelc, cpvion, curent,
     .  denmon, dtmax , dtmin , eioniz, elecd0, electe, flgas , flimp ,
     .  gflmax, gfract, gftime, apresr, rastar, rmajor, rminor, rpa   ,
     .  rpela , rqs   , sedit , splot , tbpoid, tcold , tcoldp, tcomp ,
     .  tedit , tgas  , timp  , tmax  , tpela , tplot , vpela , nedit ,
     .  npel  , npel2 , npelga, npelou, npuff , nskip , ntty  , cfz230,
     .  itvpel
c---------------------------------------------------------------------
cl              1.         new run
c
         if(nlres) go to 200
         read(nread,nurun1)
         read(nread,nurun2)
c
          if(.not.nlpomt(20)) write(nout,nurun1)
        call rarray('cfutz   ',cfutz,500)
          if(.not.nlpomt(20)) write(nout,nurun2)
c
c..test to see if namelist nurun3 is there and, if so, read it
c
      if ( linetest(nread,'nurun3') .eq. 1 ) then
          read(nread,nurun3)
          if(.not.nlpomt(20)) write(nout,nurun3)
      endif
c     
         return
c
c---------------------------------------------------------------------
cl              2.         run restarted
c
  200    continue
         read(nread,reset)
         write(nout,reset)
c
         return
         end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@linetest
c
      function linetest (lunit, string)
c
c  function linetest tests to see if the character string 'string'
c  appears anywhere on the next line and then backs up so that the
c  line can be read over again.
c
c    On output:
c  linetest = 1 if the character string does appear on the line
c  linetest = 0 if not
c
c  lunit  = unit number to be read (input integer)
c  string = character string with no more than 8 characters
c
      character *(*) string
      character *80 line
c
      linetest = 0
c
      linelen = 80
c
c..find the length of the character string
c
      istring = len ( string )
      if ( istring .lt. 1 ) return
      istring = min ( istring, 8 )
c
c..read one line and then backspace to the beginning of that line
c
      read (lunit, 100, end=90 ) line
 100  format (a)
      backspace lunit
c
c..test to find string on line.  If so, set linetest = 1
c
      do 10 j=1, linelen-istring
        if ( string(1:istring) .eq. line(j:j+istring-1) ) linetest = 1
 10   continue
c
 90   continue
c
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@redata   /baldur/code/bald  file DIO
c rgb 14-oct-90 18.73 moved reread from sbrtn stepon file dsolver to here
c
      subroutine redata
c
c     Reread input data.  Called from sbrtn stepon in file dsolver.
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clsaw.m'
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
        namelist /nmread/ reread,leqtyp,nrun,nlomt1,nlomt2,nlomt3,cfutz
     &,rehe,rehi,reoff ,reon,rex1,rex2,rey1,rey2,rle0,rle1,rleprf,rlepwr
     &,rlfreq,rli0,rli1,rliprf,rlipwr,rlpara,rlpowr,nrldmp,nrleex,nrliex
     &,rlhepwr,rlheprf
     &,nrlmod,nrlpar,extf,extzef,afslow,rcwals,rdwals,ntype,nfusn
     &,dfutzd,dfutze,dfutzi,dfutzv,dxsemi,vxsemi,xesemi,xisemi
     &,cemprc,cgbeta,cgpowr,gcoefs,gxp ,gcomb ,nrad
     &,grecyc,gwmin,gxmaxe,gxmine,gxphi,gxthet,ngpart,ngprof
     &,nbound,nlglim,nlgmon,nlgpin,nlgref,nlgspt,natomc
     &,nitmax,theta,thetap,delmax,errmax,smrlow,smlwcy,lsmord,rqs
     &,ntrans,nldiff,nlextr,nlrcom,nlsorc,nlsord,xfutz
     &,nedit,nplot,splot,nsedit,sedit,nskip,ntty,nlpomt
     &,lnumer,cnumer,linout,cinout,ltrnsp,ctrnsp,lemprc
     &,lbeams,cbeams,ldivrt,cdivrt,lauxil,cauxil,limprd,cimprd
     &,lnugas,cnugas,lneocl,cneocl,lfsion,cfsion,lstabl,cstabl,versno
     &,lsawth,csawth,swton,swtoff,swperd,swqmin,swxmin,lsweq,cequil
     &,hton,htoff,hangle,habeam,haper,haperv,hdiv,hdivv,hebeam,height
     &,hfocl,hfoclv,hfract,hfutz,hibeam,hlenth,hnchek,hpbeam,hprof,hr
     &,hprofv,hrmaj,hrmin,htchek,hwidth,nhaper,nhbeam,niprof,nizone
     &,nhe,nhmu,nhprfv,nhprof,nhshap,nhskip,nhsrc,nipart,cjbeam,hpowmw
     &,npuff,eioniz,tcold,tcoldp,cpvelc,cpvion,dtmax,dtmin,tmax
     &,fdr,fig,fec,feg,fti,fdrint,frm,fkb,frb,fmh,fhf,lthery,cthery
cjk add fdr, etc to reread namelist
c
c..items left out of the reread namelist
c
c    &,tpela,vpela,rpa,rpela,npelga,npelou,npelgc,frout,frcor,rpelc,ypa
c    &,bdhyde,bdimpe,bdtee ,bdtie ,bdtime
c    &,flgas,flimp,gflmax,gfract,gftime
c
c
      cfutz(82) = 0.
      reread  = 1. / epslon
c
c
      write (nout,120)  nstep, tai*uist
 120  format (///5x,70('*')/
     & 5x,'REREAD NAMELIST AT STEP ',i4,'   and time ',f8.3,' seconds'
     & /,5x,70('*')//)
c
      read (nread,nmread)
c
      if ( .not. nlpomt(20) ) write (nout,nmread)
c
      dtmaxi=ueit*dtmax
      thetai = theta
c
c
c  for backward compatability,  CFUTZ(82) may be used instead
c
      if ( cfutz(82) .gt. epslon ) reread = cfutz(82)
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@sprint   /baldur/code/bald  file DIO
c  rgb 25-mar-02 commented out call impprt
c  les  nov-90  add synchrotron radiation to electon losses
c  les  nov-90  d-3he mass includes helium, protons
c  les  nov-90  d-3he printout from sbrtn d3prin
c  rgb 16-feb-92 for the moment, comment out call to d3prin
c  les 11-feb-92 if cfutz(490) .gt. epslon, include D-3He mass and print
c  rgb 20.21 21-dec-91 do 21 jz=2,isep changed to do 21 jz=2,isep-1
c  rgb 20.08 08-aug-91 print total energy and power balance lines
c  rgb 21-aug-90 18.46 ITER89-P scaling printed out
c  rgb 18-aug-89 16.18 changed the mass dependence in Goldston scaling
c    to agree with recent standard by Kaye sqrt(mass (AMU))
c    Removed Kaye-Goldston scaling altogether from short output
c    Print out inverse square combination of Goldston and neoalcator
c    Slightly shortened format of confinement time line
c    Added magnetic field in tesla to short output
c  dps 10-nov-88 15.07 generalize expressions for volume and line average
c                densities for 1-1/2-D and scrape-off; disable calls to
c                IMPPRT and FDBACK when using IRE code.
c  rgb 19:00 21-mar-87  corrected Goldston scaling and changed 5x -> 1x
c  rgb 22:00 18-mar-87 changed tau-e comparisons, added Goldston scaling
c******************************************************************************
c
c
        subroutine sprint(kunit,k)
c
c
cl      3.3     generate short printout
c
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
!cap
        character*(2) :: ihcomp(2) 
	data ihcomp/'  ','**'/
c
c
c------------------------------------------------------------------------------
c
c
        data    iclass /3/,     isub /3/
        data    iheflx /70/,    impneu /200/
c
c
        if (.not.nlomt3(isub)) go to 10
        call mesage(' *** 3.3 subroutine sprint bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
c
        if (k.lt.1.or.k.gt.3) return
        if (kunit.lt.0) return
        iunit = kunit
        if (kunit.eq.0) iunit = nprint
        go to (100,200,300), k
c
c
cl      1)      first timestep edit
c
c
  100   continue
c
c
cl      2)      intermediate printouts
c
c
  200   continue
c
c
        if(nadump(1).le.lcentr) then
          isep   = mzones
          isepm1 = ledge
          zscale = 2.0 * used
        else
          isep   = nadump(1)
          isepm1 = isep - 1
          zscale = zscale * avi(mzones,12,1) / avi(isep,12,1)
        endif
c
cl              average quantities
c
c
          zee = 0.0
          zei = 0.0
          znz2 = 0.0
          zne = 0.0
          zni = 0.0
c
        do 204 jz = lcentr, isepm1
          zee = zee + chi(lelec,jz) * dx2i(jz)
          zei = zei + chi(lion,jz) * dx2i(jz)
          znz2 = znz2 + xzeff(2,jz) * rhoels(2,jz) * dx2i(jz)
          zne = zne + rhoels(2,jz) * dx2i(jz)
          zni = zni + rhoins(2,jz) * dx2i(jz)
  204   continue
c
        ztebar = zee * uisd*uieh / zne
        ztibar = zei * uisd*uieh / zni
        zzbar =  znz2 / zne
        znebar = zne * zscale
        znibar = zni * zscale
c
c
cl              central values
c
c
        zte0 = tes(2,lcentr) * useh
        zti0 = tis(2,lcentr) * useh
        zne0 = rhoels(2,lcentr) * used
        zni0 = rhoins(2,lcentr) * used
        zzeff0 = xzeff(2,lcentr)
        zq0 = q(lcentr)
c
c..geometry
c
        za = rmins * usel * 0.01
        zr = rmajs * usel * 0.01
        zelong = elong(isep,1)
        ztring = triang(isep,1)
c
        zt = tai * uist*1000.0
        zdt = dtoldi * uist*1000.0
c
c  zi     = total toroidal plasma current in MA
c  zbtor  = toroidal field in tesla
c  zvl    = loop voltage in volts
c
      zi = avi(mzones,2,1) * avi(mzones,3,1) * r0ref * bpoli(mzones)
     &  * avi(mzones,7,1) * 1.e-6 / (twopi * emu0)
c
      zbtor = bzs * usib
c
cl              loop voltage
c
      zvl = vloopi(ledge,2)  ! loop voltage at edge zone center
c
cbate   zvl = 2.0*fcpi * rmajs * fcc*100.0 *
cbate1   (ajzs(2,ledge)-cjbeam*ajbs(ledge))*eta(2,ledge)
c
c
        write (iunit,10100) nstep, zte0, zti0, zne0, zni0, zzeff0,
     &          zt, ztebar, ztibar, znebar, znibar, zzbar
c
10100   format(/1x,'time-step=',i5,'  te-axis=',1pe10.3,
     1          '  ti-axis=',1pe10.3,'  ne-axis=',1pe10.3,
     2          '  ni-axis=',1pe10.3,'   z-axis=',1pe10.3/1x,
     3          't=',0pf10.3,' ms','  te-avg.=',1pe10.3,
     4          '  ti-avg.=',1pe10.3,'  ne-avg.=',1pe10.3,
     5          '  ni-avg.=',1pe10.3,'   z-avg.=',1pe10.3)
c
      write (iunit,112)
     &          zdt, za, zr, zelong, ztring, zbtor, zi, zvl, zq0
c
 112  format (' dt=',0pf9.3,' ms  a=',0pf6.3,' R=',0pf6.3
     & ,' elong=',0pf6.3,' delta=',0pf6.3,' B=',0pf7.3,' I=',0pf7.3
     & ,' vloop=',1pe10.3,' q-axis=',1pe10.3)
c
c..total plasma energy and power balance lines
c  use zergmj to convert energies to mega joules
c
        zergmj = usee * 1.e-6
      write (iunit,121) tai*uist, ergts(isep)*zergmj
     & ,erges(isep)*zergmj, ergis(isep)*zergmj
     & ,ergbs(isep)*zergmj, ergas(isep)*zergmj
c
 121  format (' t=',0pf10.6,' s Wtot =',0pf12.6
     & ,' Weth =',0pf12.6,' With   =',0pf12.6
     & ,' Wbeam=',0pf12.6,' Walpha=',0pf12.6
     & ,' MJ')
c
c  zpeirs = impurity radiation power
c
            zpeirs = 0.0
      if ( mimp .gt. 0 ) then
        do 23 ji = 1, mimp
          do 21 jz = 2, isep-1
            zpeirs = zpeirs + weirs(ji,jz)*(vols(jz+1,1)-vols(jz,1))
  21      continue
  23    continue
      endif
c
c  zeiequ = equipartition power from electrons to ions in watts
c
        zeiequ = 0.0
      do 25 jz = 2, isep
        zeiequ = zeiequ + ( vols(jz+1,1)-vols(jz,1) ) * ( weiths(jz)
     &   + cnueqs(jz) * ( tes(2,jz) - tis(2,jz) ) ) * usip
  25  continue
c
c  use zpowmw to convert from watts to mega watts
c      zpsmw  to convert from ergs/sec to MW
c
        zpowmw = 1.e-6
        zpsmw  = usip * 1.e-6
c
c  zeheat = total heating power to electrons in watts
c  ziheat = total heating power to ions
c  zeloss = total power lost by electrons
c  ziloss = total power lost by ions
c  zeptot = total power to and from electrons
c  ziptot = total power to and from ions
c  zpwtot = total power to and from plasma
c
        zeheat = usip * ( geohms(isep) + gealfs(isep) + geauxs(isep)
     &                  + gebems(isep) + geecrs(isep) + geicrs(isep) )
c
        ziheat = usip * ( gialfs(isep) + giauxs(isep) + gibems(isep)
     &                  + giecrs(isep) + giicrs(isep) )
c
        zeloss = usip * ( gebrs (isep) + geions(isep) + gesrs (isep)
     &                  + zpeirs )
c
        ziloss = usip * ( - gichxs(isep) - giions(isep) )
c
        zeptot = zeheat - zeloss - zeiequ
c
        ziptot = ziheat - ziloss + zeiequ
c
        zpwtot = zeptot + ziptot
        zpheat = zeheat + ziheat
        zploss = zeloss + ziloss
c
cl              d(W_total) / dt
c
        zedot = fflxii(lelec) + fflxoi(lelec) + fsorci(lelec) +
     1          fsordi(lelec) + fcompi(lelec) + fflxii(lion) +
     2          fcompi(lion)  + fflxoi(lion)  + fsorci(lion) +
     3          fsordi(lion)
c
      write (iunit,123) tai*uist, zeptot*zpowmw
     & , (geauxs(isep)+gebems(isep)+geecrs(isep)+geicrs(isep))*zpsmw
     & , gealfs(isep)*zpsmw, geohms(isep)*zpsmw
     & , (gebrs(isep)+geions(isep)+gesrs(isep)+zpeirs)*zpsmw
c
 123  format (' t=',0pf10.6,' s Petot=',0pf12.6
     & ,' Peaux=',0pf12.6,' Pealpha=',0pf12.6,' Peohm=',0pf12.6
     & ,' Peloss=',0pf12.6,' MW')
c
      write (iunit,125) tai*uist, ziptot*zpowmw
     & , (giauxs(isep)+gibems(isep)+giecrs(isep)+giicrs(isep))*zpsmw
     & , gialfs(isep)*zpsmw, zeiequ*zpowmw
     & , ziloss*zpowmw
c
 125  format (' t=',0pf10.6,' s Pitot=',0pf12.6
     & ,' Piaux=',0pf12.6,' Pialpha=',0pf12.6
     & ,' Piequ=',0pf12.6,' Piloss=',0pf12.6,' MW')
c
      write (iunit,127) tai*uist, zpwtot*zpowmw
     & , zpheat*zpowmw, zploss*zpowmw, zedot*zpowmw
c
 127  format (' t=',0pf10.6,' s Ptot =',0pf12.6
     & ,' Pheat=',0pf12.6,' Ploss  =',0pf12.6,' dW/dt=',0pf12.6
     & ,' MW')
c
c       if helium recycling is active, print relevant quantities
c  15.07 this is not applicable to the IRE code.
c
cbate        if (cfutz(iheflx).gt.epslon.and.(natomc.ne.3)) call impprt
c
c       if neutral-impurity influx is subject to feedback readjustment,
c       print readjustment factor
c  15.07 this is not applicable to the IRE code.
c
        if (cfutz(impneu).gt.1.1.and.(natomc.ne.3)) call fdback
c
c
cl..confinement times
c
        zvols = 2.0 * avi(isep,12,1) * uisl**3   ! plasma volume
        zsurfi = avi(isep,3,1) * avi(isep,5,1)   ! plasma surface area
c
      zlosei = 0.0  ! eloctron power loss due to radiation and ionization
      zlosii = 0.0  ! ion power loss due to ionization and charge exchange
c
        do 718 jz = lcentr, isepm1
c
        zeirs = 0.0   ! electron power lost due to impurity radiation
        if (mimp.gt.0) then
        do 704 ji = 1, mimp
          zeirs = zeirs + weirs(ji,jz)
  704   continue
      endif
c
c   les  nov-90  add synchr rad
c
        zlosei = zlosei + usip*zvols*dx2i(jz)*
     1    (weions(jz) + webrs(jz) + zeirs + wesyn(jz) + wesrs(jz) )
        zlosii = zlosii - usip*zvols*dx2i(jz)*
     1                          (wiions(jz) + wichxs(jz))
c
  718   continue
c
      zflosi = 0.0   ! loss due to heat flux
c
      do 716 jp = lelec,lion
c
          zfluxi = 0.0
        do 714 jp2 = 1, mchi
          zfluxi= zfluxi- aaaa(jp,jp2,isep)*chi(jp2,isepm1)
     &                  - bbbb(jp,jp2,isep)*chi(jp2,isep)
  714   continue
c
          zflosi = zflosi + zfluxi * zsurfi * ueit
  716 continue
c
      zloss = (zlosei + zlosii + zflosi) * uisp
c
      if (zloss .gt. epslon)
     &  ztaue = (erges(isep) + ergis(isep)) / zloss
c
c..q-cylindrical
c
      zk = elong(isep,1)
      zqcyl = ( 5. * ahalfs(isep,1)**2 * bzs*useb * ( 1. + zk**2 ) )
     &         / ( 2. * rmids(isep,1) * curnts(isep) * usei )
c
c..confinement scalings
c
      zbeta = max (epslon, cgbeta(1) * betate(isep)
     &          + cgbeta(2) * betati(isep)
     &          + cgbeta(3) * betatb(isep)
     &          + cgbeta(4) * betata(isep) ) * 100.   ! vol ave beta %
c
      zpowr = max (epslon, cgpowr(1)*geohms(isep)
     &   + cgpowr(3)*gealfs(isep) + cgpowr(4)*gialfs(isep)
     &   + cgpowr(5)*geauxs(isep) + cgpowr(6)*giauxs(isep)
     &   + cgpowr(7)*gebems(isep) + cgpowr(8)*gibems(isep)
     &   + cgpowr(9)*geecrs(isep) + cgpowr(10)*giecrs(isep)
     &   + cgpowr(11)*geicrs(isep) + cgpowr(12)*giicrs(isep) )
c
      zpowr = zpowr * usep * 1.e-6     ! convert to MW
c
c  effective mass in AMU
c   les  nov-90  d-3he mass includes helium, protons
c
      zmeff = aspec(1)
      if (mhyd .eq. 2) zmeff = ( aspec(1) + aspec(2) ) / 2.
      if ( cfutz(490) .gt. epslon )
     & zmeff=(aspec(1) + aspec(2)+aspec(3)+aspec(4)+aspec(5))/5.
c
c  Neoalcator scaling [in seconds]
c
      ztoh = 7.e-22 * enlaes(isep) * ahalfs(isep,1)
     &  * (rmids(isep,1))**2 * zqcyl
c
      ztoh = min (ztoh,1./epslon)
c
      zftoh = ztaue / ztoh  ! ratio of confinement times
c
c..Goldston scaling [Plasma Phys and Cont. Fus. 26 (1984) 87]
c  with mass correction by Kaye, Aug 1988
c  as a function of heating power
c
       ztauxg = 0.0302 * (abs(curnts(isep))*usei*1.e-3)
     &  * (ahalfs(isep,1)*0.01)**(-0.37) * zk**0.5
     &  * (rmids(isep,1)*0.01)**1.75 * zpowr**(-0.5)
     &  * (abs(zmeff))**0.5
c
      ztauxg = max ( epslon, min (ztauxg,1./epslon) )
c
      zftaug = ztaue / ztauxg  ! ratio of taue / tau_Goldston
c
c..ITER89-P scaling
c
      ziter  = 0.048 * (abs(zmeff))**0.5
     &  * ( abs(curnts(isep))*usei*1.e-3)**0.85
     &  * (rmids(isep,1)*0.01)**1.2 * (ahalfs(isep,1)*0.01)**0.3
     &  * zk**0.5 * (enlaes(isep)*1.e-14)**0.1 * (bzs*1.e-4)**0.2
     &  * zpowr**(-0.5)
c
      ziter  = max ( epslon, min ( ziter, 1./epslon ) )
c
c  ratio of tau to ITER89-P scaling
c
      zfiter = ztaue / ziter
c
c  inverse squared combination
c
cbate      ztcomb = sqrt (1. / (1./ztoh**2 + 1./ztauxg**2) )
c
cbate      zfcomb = ztaue / ztcomb  ! ratio of taue / tau-Combination
c
      write (iunit,131) tai*uist,ztaue,zftaug,ztauxg,zfiter,ziter
     &  ,zftoh,ztoh
 131  format (' t=',0pf10.6,' s taue =',0pf9.5
     &  ,' = (',0pf8.3,')*',0pf9.5,' GL, = (',0pf8.3,')*',0pf9.5
     &  ,' ITER89-P,  (',0pf8.3,')*',0pf9.5,' NA')
c
      write (iunit,132) tai*uist,enlaes(isep),zbeta,zqcyl
 132  format (' t=',0pf10.6
     &  ,' s ne-bar=',1pe11.3,' cm-3, beta = ',0pf8.3,'%'
     &  ,' qcyl=',0pf6.3)
c
c   les  nov-90  d-3he printout
c
      if ( cfutz(490) .gt. epslon ) call d3prin(iunit,isep,zbeta)
c
        return
c
c
cl      3)      final print-out
c
c
  300   continue
        go to 200
c
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@aprint   /baldur/code/bald  file DIO
c       rgb 4-sep-85 1 1/2 D modification of volume, surface area, radius
c                volume within zone bndry j = avi(j,12,1) in internal units
c                surface area of zone bndry j = avi(j,3,1) * avi(j,5,1)
c                minor radius (half width) to zone bndry j = avi(j,15,1)
c******************************************************************************
c
c
        subroutine aprint(ktype,k)
c
c
c..sbrtn aprint prints out information about alpha particle effects
c
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'calpht.m'
c
c
        logical
     l  ilend
c
c------------------------------------------------------------------------------
c
        data    iclass/3/,      isub/9/
c
        if (.not.nlomt3(isub))  go to 10
        call mesage(' *** 3.9 subroutine aprint bypassed ')
        return
   10   continue
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       lpage (comout)
c
c------------------------------------------------------------------------------
c
        if (ktype.gt.0) return
        if (k.ne.1)     go to 200
c
        return
c
c
c
cl      2.)     print first page of large time-step edit
c
c
  200   continue
        if (nlpomt(19)) return
c
c               if no alphas, don't print out
c
        if (acons(2).lt.1.0) return
c
        zt = atbi * uist * 1000.0
        zdt = dtoldi * uist * 1000.0
        lpage = lpage + 1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write (nprint,10201) lhlen, lhdens, lheden, lhtime, lhdenr,
     1          lhdenr, lhdens
c
        zivole = used / (avi(mzones,12,1) * uisl**3)
c
c
c       loop begins here
c
c
        i1 = nskip + 1
        do 235 jz = i1, ledge, nskip
c
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uiel ! zone center radius
c
        zna = alphai(jz) * uied
        zea = ealfai(jz) * alphai(jz) * uied*uiee
        ztslow = 0.0
        if (aslows(jz).gt.epslon) ztslow = uset / aslows(jz)
        ztherm = afuses(jz) * used*uest
        zfusn  = ztherm + halfas(jz) * used*uest
c
        z0 = (weohms(jz) + webems(jz) + wibems(jz))
        zqe = 0.0
        if (z0.gt.epslon) zqe = zfusn * 17.5e+06*evs / z0 * uset*uesd
        ztfuse = (abfuse(jz) + atfuse(jz)) * zivole*dx2inv(jz)*0.5
        if(wealfs(jz).ne.0.) zfalfi=wialfs(jz)/(wealfs(jz)+wialfs(jz))
c
c
c
        ijz = jz - 1
        write (nprint,10202) ijz, zrad, zna, zea, ztslow, ztherm,
     1                          zfusn, aoloss(jz), zqe, ztfuse,zfalfi
c
c
  235   continue
c
c
c       loop has ended
c
c
c
c       compute and print totals
c
c
        zzna = 0.0
        zzea = 0.0
        zztsl = 0.0
        zztthr = 0.0
        zztfus = 0.0
        zzloss = 0.0
        zzqeb = 0.0
        zzqebs = 0.0
        zzqtot = 0.0
        zzftot = 0.0
c
        do 258 jz = lcentr, ledge
        zzna = zzna + dx2i(jz)*alphai(jz)
        zzea = zzea + dx2i(jz)*ealfai(jz)*alphai(jz)
        zztsl = zztsl + dx2i(jz)*aslows(jz)*ealfai(jz)*alphai(jz)
        zztthr = zztthr + dx2i(jz)*afuses(jz)
        zztfus = zztfus + dx2i(jz)*(afuses(jz) + halfas(jz))
        zzloss = zzloss + dx2i(jz)*aoloss(jz)*(afuses(jz) + halfas(jz))
        zzftot = zzftot + (atfuse(jz) + abfuse(jz))
        zzqtot = zzqtot +
     1                  dx2i(jz)*(webems(jz) + wibems(jz) + weohms(jz))
        zzqebs = zzqebs + dx2i(jz)*(webems(jz) + wibems(jz))
  258   continue
c
c
c               average values
c
c
        zna = zzna * 2.0*uied
        zea = zzea * 2.0 * uied*uiee
        ztslow = 0.0
        if (zztsl.gt.epslon) ztslow = zzea / zztsl * uset
        ztherm = zztthr * 2.0 * used * uest
        zfusn = zztfus * 2.0 * used * uest
        zoloss = 0.0
        if (zztfus.gt.epslon) zoloss = zzloss / zztfus
        zqe = 0.0
        if (zzqtot.gt.epslon) zqe = (zztfus * 17.5e+06*evs) / zzqtot
        ztfuse = zzftot * zivole
c
        write (nprint,10203) zna, zea, ztslow, ztherm, zfusn, zoloss,
     1          zqe, ztfuse
c
c
c               volume totals
c
c
        zvole = avi(mzones,12,1) * ueid
c
        zna = zna * zvole
        zea = zea * zvole
        ztherm = ztherm * zvole
        zfusn = zfusn * zvole
        zoloss = zoloss * zfusn
        ztfuse = ztfuse * zvole
c
        write (nprint,10204) lhener, lhtime, zna, zea, ztherm, zfusn,
     1          zoloss, ztfuse
c
c
c               write out status information
c
        write(nprint,10210) npactv,npdead,nldead,ntsmpl,ntperm,nulong,
     1  nsqeze,nbounc,nstart,narrow,ndzerr,nslost,mskip
        zz=actvfe*usee*evs
        zz1=actvse*usee*evs
        write(nprint,10211) actvfw,zz,actvsw,zz1
        return
c
c
cl      format statements
c
c
10200   format(1h1,2x,a48,10x,a72//
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     2  'time =',f12.3,'  millisecs.',12x,'dt =',f12.6,'  millisecs.')
10201   format(/4x,'zone',3x,'radius',2x,'hi-e alphas',3x,'energy',4x,
     1          'slowdown t.',1x,'therm-fusion rate-total',1x,
     2          '1st orb loss',4x,'qe',3x,'total fusion',3x,'ion'
     2          /8x,a10,1x,a10,2(2x,a10)
     3          ,2(3x,a10),2x,'per-cent',12x,a10,3x,'fract')
10202   format(4x,i4,2x,0pf7.2,5(2x,1pe10.3),2x,2pf8.3,3x,0pf8.3,
     1  3x,1pe10.3,2x,0pf5.3)
10203   format(20x,4('-------------------------')/
     1          '   ** averages ** ',1x,5(1pe10.3,2x),2pf8.3,3x,0pf8.3,
     2          3x,1pe10.3)
10204   format(20x,'particles',2x,a10,14x,'---- particles per ',a10,
     1          ' ----',13x,'particles'/'   **  totals  **',2(2x,
     2          1pe10.3),12x,2(2x,1pe10.3),1x,1pe10.3,13x,1pe10.3)
10210   format(/,/,5x,'npactv,npdead,nldead=',i4,',',i4,',',i4,
     1  4x,'ntsmpl,ntperm=',i4,',',i4,4x,'nulong,nsqeze=',i5,',',i3,
     2  /,5x,'nbounc,nstart,narrow,ndzerr=',i6,3(',',i4),4x,
     3  'nslost,mskip=',i5,',',i4)
10211 format(5x,'# fast particles=',1pe10.3,5x,'fast energy=',1pe10.3,
     1  5x,'# slow particles=',1pe10.3,5x,'slow energy=',1pe10.3)
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
c@hprint   /baldur/code/bald  file DIO
c  rgb 27-jun-94 print rhobis, rhobes, and hebems with ajbs
c  rgb 07-jan-92 changed zdarea from surface area to cross-sectional area
c  rgb 05-jun-89 commented out unused format statements
c       dps 16-dec-87 match to 1-D code and allow for more than one
c                     beam species.
c       rgb 4-sep-85 1 1/2 D modification of volume, surface area, radius
c                volume within zone bndry j = avi(j,12,1) in internal units
c                surface area of zone bndry j = avi(j,3,1) * avi(j,5,1)
c                minor radius (half width) to zone bndry j = avi(j,15,1)
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c******************************************************************************
c
c
        subroutine hprint(ktype,k)
c
c
cl      3.5     print out beam data
c
c
c
c
        include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
        include 'cfokkr.m'
c
c
        dimension
     r  zd(20),znb(100),zmui(20), zhei(20),zshei(20),zsmu(20),
     r  zznb(2),zzez(2),zzep(2),zzepp(2),zzezp(2),zzmub(2),
     r  zzet(2),zzetp(2)
c
c
c
c------------------------------------------------------------------------------
c
c
        data    iclass  /3/,    isub    /5/
c
c
        if (.not.nlomt3(isub)) go to 10
        call mesage(' *** 3.5 subroutine hprint bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
c
c               local variables:
c
c       ibeam1  - index of first injector w/ species specified
c       zchksm  -- test for hangle warning message
c       zd      -- scratch array for printout
c       zden    -- part. / vol. for this zone
c       zdt     - timestep length in etu
c       zeden   -- total energy / vol. at this zone
c       zepar   -- parallel energy / part. at this zone
c       zeperp  -- perpendicular energy / part. at this zone
c       zetot   -- total energy / part. at this zone
c       zfuse   -- no. of fusions / time / vol.
c       zfutot  -- no. of fusions / time
c       zib     -- total beam-driven current
c       zit     -- average beam-driven current
c       zmubar  -- avg. mu / part.
c       zmudev  -- stand. dev. of mu among part.
c       zmui(jmu)-- sin of pitch ang. * dmu
c       zrad    -- radius in elu
c       zt      -- end of timestep time in etu
c       zvoli   -- volume of torus = zvoli * sum dx2i
c       zzchex  -- part. / sec. lost from beam d.t. charge exchange
c       zzechk  -- energy check for beam
c       zzep    -- total perp. energy in beam
c       zzepp   -- avg. perp. energy / part. in beam
c       zzet    -- total energy in beam
c       zzetp   -- avg. energy / part. in beam
c       zzez    -- total paral. energy in beam
c       zzezp   -- avg. par. energy / part. in beam
c       zzheat  -- energy / time beam heating of plasma
c       zzloss  -- part. / time lost from beam d.t. loss cones (not used)
c       zzmub   -- avg. mu / part. in beam
c       zzmudv  -- std. dev. of mu / part. in beam
c       zznb    -- number of part. in beam
c       zznchk  -- part. check for beam
c
c
        if (ktype.gt.0) return
        if (k.ne.1) go to 200
c
c
c       here insert title page generator
c
c
  100   continue
c
        return
c
c
cl      first page
c
c
c
  200   continue
        if (nlpomt(8).and.nlpomt(13).and.nlpomt(14).and.nlpomt(15)
     1          .and.nlpomt(16)) return
        if (nstep.eq.0) return
c
        do 202 jb = 1, mhbeam
        if (hton(jb)*ueit.lt.tai) go to 204
  202   continue
        return
c
  204   continue
c
        if (nlpomt(8)) go to 210
c
c       19-oct-79 code removed for convenience due to increase
c       in no. of beam species allowed.
c
  210   continue
        zt = tai * uist*1000.0
        zdt = dtoldi * uist*1000.0
        zvoli = 2.0 * avi(mzones,12,1)
c
        if (nlpomt(13)) go to 237
c
c       this section also removed due to increase in no. of
c       beam species allowed.
c
  237   continue
c
c
c       totals
c
c
        zznbt = 0.0
        zzett = 0.0
        call resetr (zd,mxhe,0.0)
c
        do 260 jsp = 1, mhsp
        zznb(jsp) = 0.0
        zzep(jsp) = 0.0
        zzez(jsp) = 0.0
        z4 = 0.0
c
        do 250 je = 1, nhmu
        z0 = 0.0
        z1 = 0.0
        z2 = 0.0
c
        do 245 jz = lcentr, ledge
        do 245 jmu = 1, nhmu
        z0 = z0 + hfi(jmu,je,jz,jsp) * dx2i(jz) * hdmuc(jmu)
        z1 = z1 + hfi(jmu,je,jz,jsp) * dx2i(jz) * hdmu3(jmu)
        z2 = z2 + hfi(jmu,je,jz,jsp) * dx2i(jz) * (hdmuc(jmu)
     1       - hdmu3(jmu))
        z4 = z4 + hfi(jmu,je,jz,jsp) * dx2i(jz) * hmuc(jmu)
     1       * hdmuc(jmu)
  245   continue
c
        zd(je) = zd(je) + z0 * zvoli
        zznb(jsp) = zznb(jsp) + z0*zvoli
        zzez(jsp) = zzez(jsp) + hei(je,jsp)*uiee * z1 * zvoli
        zzep(jsp) = zzep(jsp) + hei(je,jsp)*uiee * z2 * zvoli
  250   continue
c
        zzepp(jsp) = 0.0
        zzezp(jsp) = 0.0
        zzmub(jsp) = 0.0
        if (zznb(jsp).le.epslon) go to 255
        z0 = uese*useh / zznb(jsp)
        zzepp(jsp) = zzep(jsp) * z0
        zzezp(jsp) = zzez(jsp) * z0
        zzmub(jsp) = z4*zvoli / zznb(jsp)
  255   continue
c
        zzetp(jsp) = zzepp(jsp) + zzezp(jsp)
        zzet(jsp) = zzep(jsp) + zzez(jsp)
        zznbt = zznbt + zznb(jsp)
        zzett = zzett + zzet(jsp)
  260   continue
c
c
        if (.not.nlpomt(13)) write (nprint,10205) (zd(je),je=1,nhe)
c
c
c
cl              second page
c
c
c
  300   continue
        if (nlpomt(14).and.nlpomt(15)) go to 395
c
        do 390 jsp = 1, mhsp
        if (htspon(jsp)*ueit.gt.tai) go to 390
        ibeam1 = nhbem1(jsp)
c
        lpage = lpage + 1
c
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
        if (nlpomt(14)) go to 336
        write (nprint,10302) jsp
        write (nprint,10301) lhlen, lhtemp, lhtemp, lhtime, lhdens,
     1                                  lheden, lhdenr, lhdenr
c
c
        i1 = nhskip+1
        do 335 jz = i1, ledge, nhskip
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uiel ! zone center radius
        z0 = 0.0
        z1 = 0.0
        z2 = 0.0
        z3 = 0.0
        zchex = 0.0
        zfuse = 0.0
        ztaus= 0.0
        z11 = 0.0
c
        do 310 je = 1, nhe
        z12 = hchexs(je,jz,jsp) * hei(je,jsp)
        if(je.ge.lhemin(jz,jsp)) ztaus = ztaus + 1./hslows(je,jz,jsp)
c
        do 310 jmu = 1, nhmu
        z0 = z0 + hfi(jmu,je,jz,jsp) * hei(je,jsp) * hdmu3(jmu)
        z1 = z1 + hfi(jmu,je,jz,jsp) * hei(je,jsp)
     1       * (hdmuc(jmu) - hdmu3(jmu))
        z2 = z2 + hfi(jmu,je,jz,jsp) * hmuc(jmu) * hdmuc(jmu)
        z3 = z3 + hfi(jmu,je,jz,jsp) * hdmuc(jmu)
        zchex = zchex + hfi(jmu,je,jz,jsp) * hdmuc(jmu) * z12
  310   continue
c
        zepar = 0.0
        zeperp = 0.0
        zetot = 0.0
        ztaus = ztaus * uset
        zeden = 0.0
        zmubar = 0.0
        zden = z3 * uied
        zfuse = halsps(jz,jsp) * used * uest
        zchex = zchex * uied * uise * usep
        zsourc = 0.0
        i001 = lhbeam(ibeam1)
        if(i001.gt.0)
     1          zsourc = shbems(i001,jz) * used*uset
        zifrac = 0.0
        z4 = websps(jz,jsp) + wibsps(jz,jsp)
        if (z4.gt.epslon) zifrac = wibsps(jz,jsp) / z4
        if (z4.gt.epslon) zchex = zchex / (z4 * usep + zchex)
c
        if (zden.le.epslon) go to 312
c
        zepar = z0 / z3 * uise * useh
        zeperp = z1 / z3 * uise*useh
        zetot = zeperp + zepar
        zeden = (z0 + z1) * uiee * uied
        zmubar = z2 / z3
  312   continue
c
        z4 = 0.0
c
        do 315 je = 1, nhe
        do 315 jmu = 1, nhmu
c
        z4 = z4 + (zmubar - hmuc(jmu))**2 * hfi(jmu,je,jz,jsp)
     1       * hdmuc(jmu)
  315   continue
c
        zmudev = 0.0
        if (z3.gt.epslon) zmudev = sqrt(z4/z3)
c
        ijz = jz - 1
        write (nprint,10304) ijz, zrad, zepar, zetot, ztaus,
     1  zden, zeden, zmubar, zmudev, zchex, zsourc, zfuse, zifrac
c
  335   continue
c
c
c
cl      totals
c
c
c
  336   continue
        z0 = 0.0
        z1 = 0.0
        z2 = 0.0
        z3 = 0.0
        z4 = 0.0
        z5 = 0.0
        z6 = 0.0
        z7 = 0.0
        z8 = 0.0
        z9 = 0.0
        z11 = 0.0
        zzalfa = 0.0
        zzsorc = 0.0
        zfutot = 0.0
c
        do 345 jz = lcentr, ledge
        z0 = z0 + hninjs(jz) - hnloss(jz) - hnplas(jz) + hncmps(jz)
        z1 = z1 + heinjs(jz) - heloss(jz) - heplas(jz) + hecmps(jz)
        zzalfa = zzalfa + halfas(jz) * dx2i(jz)
        zfutot = zfutot + halsps(jz,jsp) * dx2i(jz)
        z5 = z5 + webems(jz) * dx2i(jz)
        z6 = z6 + wibems(jz) * dx2i(jz)
        i001 = lhbeam(ibeam1)
        if (i001.gt.0)
     1          zzsorc = zzsorc + shbems(i001,jz)*dx2i(jz)
c
        do 342 je = 1, nhe
c
        z10 = hdei(je,jsp) * hslows(je,jz,jsp) * dx2i(jz)
        if(je.eq.lhemin(jz,jsp)) z10 = hei(je,jsp)
     1      * hslows(je,jz,jsp) * dx2i(jz)
        if (je.lt.lhemin(jz,jsp)) z10 = 0.0
c
        do 340 jmu = 1, nhmu
c
        z2 = z2 + hfi(jmu,je,jz,jsp) * hdmuc(jmu) * dx2i(jz) *
     1                             (zzmub(jsp) - hmuc(jmu))**2
        z3 = z3 + hfi(jmu,je,jz,jsp) * hdmuc(jmu) * z10
        z4 = z4 + hfi(jmu,je,jz,jsp) * hdmuc(jmu) * dx2i(jz) *
     1                             hchexs(je,jz,jsp) * hei(je,jsp)
  340   continue
c
  342   continue
c
  345   continue
c
c
        zzalfa = zzalfa * zvoli * usid
        zfutot = zfutot * zvoli * usid
c
        z10 = 10.0**(-fxes) / fces
        zznchk = 1.0 - z0 / zznbt
        zzechk = 1.0 - z1*usee / zzett
c
        zzchex = z4 * zvoli * uest * uiee
        zzloss = 0.0
        zzheat = z3 * zvoli * uest * uiee
        zzcxfr= 0.0
        if(zzchex+zzheat.gt.epslon) zzcxfr=zzchex/(zzchex+zzheat)
        zzmudv = 0.0
        if(zznb(jsp).gt.epslon) zzmudv = sqrt(z2*zvoli / zznb(jsp))
        zzfrac = 0.0
        if (z5 + z6.gt.epslon) zzfrac = z6 / (z5 + z6)
        zzsorc = zzsorc * zvoli*usid * uest
c
        write (nprint,10305) lhtime, zzezp(jsp), zzetp(jsp),
     1  zznb(jsp), zzet(jsp), zzmub(jsp), zzmudv, zzcxfr, zzsorc,
     2  zfutot, zzfrac
        if (nlpomt(15)) go to 390
c
c
        write (nprint,10306) jsp, lhener, jsp, lhpowr
        write (nprint,10307) zzep(jsp), zzchex, zzechk, zzez(jsp),
     1                   zzheat, zznchk, zzet(jsp), bbdtrs, zzalfa
c
390     continue
  395   continue
c
c
cl      4)      third page -- beam currents
c
c
  400   continue
        if (nlpomt(16)) go to 495
c
        lpage = lpage + 1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
c
        write (nprint,10401) lhlen, lhcden, cjbeam
c
        i1 = nskip + 1
        do 435 jz = i1, ledge+1, nskip
        zr = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uiel ! zone center radius
        zj = ajbs(jz) * usej
        zfr = ajbs(jz) / (ajzs(2,jz)+epslon)
        ijz = jz - 1
        write (nprint,10402) ijz, zr, zj, zfr, rhobis(2,jz)
     &    , rhobes(2,jz), hebems(jz)
  435   continue
c
c               totals
c
        zit = 0.0
        zib = 0.0
c
        do 445 jz = lcentr, ledge
          zdarea = ( avi(jz,13,1) - avi(jz-1,13,1) ) * uisl**2
          zib = zib + ajbs(jz) * zdarea
          zit = zit + ajzs(2,jz) * zdarea
  445   continue
c
        zafr = zib / zit
        zib = zib * usei
        write (nprint,10403) lhcurr, zib, zafr
c
  495   continue
c
        if (hfutz(3) .eq. 0.) go to 505
c
           zchksm = 0.0
           do 500 i2 = 1,4
           if (i2 .eq. 3) go to 500
                do 499 jz = 1,10
  499           zchksm = zchksm + hangle(i2,jz)
  500      continue
           if (zchksm .ne. 0.0) write (nprint,10405)
c
  505   continue
c
        return
c
c
cl      format statements
c
c
10200   format(1h1,2x,a48,10x,a72//
     1          2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     2          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,
     3          2x,'millisecs.')
c10201   format(/14x,a10,1x,'energy levels in ',a10,' per particle')
c10202   format(12x,8(2x,1pe9.2),2(1x,1pe9.2))
c10203   format(1x,'zone radius '/3x,a10,21x,
c     1          'densities at each energy level in ',a10)
c10204   format(1x,i3,1x,f7.2,8(2x,1pe9.2),2(1x,1pe9.2))
10205   format(12x,3('------------------------------------')/
     1          1x,'total part.',8(2x,1pe9.2),2(1x,1pe9.2))
c
cl      second page format statements
c
10301   format(/14x,'avg. energy / ptcle',5x,'taus',3x,
     x  '  beam dens',
     1  '  ener dens',' cos pitch ang','  cx power     part   ',
     2  '   fusions   ion'/' zone radius','  parallel     total  ',
     3  9x,25x,'mean std.dev.    over      source',15x,'heat'
     3  /1x,6(1x,a10),14x,'heat+cx powr',2(1x,a10),' fract')
10302   format(/44x,'energetic species no.',i1)
10304   format(1x,i3,1x,f7.2,2(1x,1pe10.3),0pf10.5,1x,2(1x,1pe10.3),
     x  2(1x,0pf6.4),0pf10.5,1x,2(1x,1pe10.3),
     1  1x,0pf5.3)
10305   format(20x,4('-------------------------')/47x,'particles',36x,
     1  '    parts.'/,a10/5x,'totals',1x,2(1x,1pe10.3),11x,
     2  2(1x,1pe10.3),2(1x,0pf6.4),3(1x,1pe10.3),1x,0pf5.3)
10306   format(1x/6x,'species no. ',i1,17x,a10,4x,'species no. ',
     1  i1,13x,a10,4x,'all species')
10307   format(10x,'total perp. energy  . . . ',1pe10.3,4x,
     1  'charge exchange loss rate ',1pe10.3,4x,'energy check  . ',
     2  1pe10.3/10x,'total parallel energy . . ',1pe10.3,4x,
     3  'plasma heating rate . . . ',1pe10.3,4x,'particle check  ',
     4  1pe10.3/10x,'total beam energy . . . . ',1pe10.3,4x,
     5  'd-t beam-beam rate  . . . ',1pe10.3,4x,'beam targt rate ',
     6  1pe10.3)
c
10401   format(/' zone radius     beam j   beam j/'
     &    
     1    ,t100,'fraction used in b'
     &    /3x,a10,2x,a10,t36,'rhobis',t49,'rhobes',t62,'hebems'
     2    ,t100,' total j             equation:',2pf6.2,' %')
10402   format(1x,i3,1x,0pf7.2,2x,1pe10.3,1x,0pf7.4
     &    ,1p3e13.4)
10403   format(/14x,a10,5x,'average'/4x,'total',5x,1pe10.3,5x,0pf7.4)
10405   format(//14x,'***warning:  beam-driven current assumes'
     1  /1x,'hangle(n,ibeam1) = 0.0 unless n=3.  values of'
     2  /1x,'hangle are given at the beginning of this printout'
     3  /1x,'(before timestep 0).'/)
        end
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@iprint   /baldur/code/bald  file DIO
c  rgb 28-jan-95 nhprof(i002) -> nhprof(i001)
c       dps 16-dec-87 splice in changes for multiple beam species
c       mhr 17-apr-86 changed column heading for active injector no
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c******************************************************************************
c
c
        subroutine iprint(kunit,k)
c
c
cl      3.8     print out injector data
c
c
c
        include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
        include 'cfokkr.m'
        include 'cfreya.m'
c
c
        dimension       zpowr(4,12),    zener(12),      zhr(12),
     r  zebeam(12),     zibeam(12),     zfract(3,12),   zton(12),
     r  ztoff(12),      zrmaj(12),      zrmin(12),      zwidth(12),
     r  zheigt(12),     zlenth(12),     zangle(4,12),   zfocl(12),
     r  zfoclv(12),     zdiv(12),       zdivv(12),      zaper(12),
     r  zaperv(12)
c
        character ihtitl(12)*10, ihname(12)*10, ihblnk(1)*4
        character ihshap(3)*10, ihprof(3)*10
c
        dimension ibperm(12)
c
        data    ihblnk /' '/
        data
     &  ihshap  /' circular ','rectangle ','**********'/,
     &  ihprof  /'   flat   ','**********','**********'/
c
c------------------------------------------------------------------------------
c
c
        data    iclass /3/,     isub /8/
c
c
        if (.not.nlomt3(isub)) go to 10
        call mesage(' *** 3.8 subroutine iprint bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
c
c
        if (kunit.gt.0) go to 1000
        if (k.gt.1) go to 200
c
c
c       here insert title page generator (if any)
c
c
  100   continue
        if (mhbeam.le.0) return
c
        write (nprint,10100) nipart, nizone, niprof, htchek, hnchek
c
c               beam parameters
c
        ib = 0
        do 104 jb = 1, mhbeam
          call reseth(ihname(jb),10,ihblnk(1))
          if (nhbeam(jb).eq.0) go to 104
          ib = ib + 1
          i001 = lhbeam(jb)
          if (i001.gt.0)  ihname(jb) = lhspec(i001)
cbat  call copyi(lhspec(1,i001),1,ihname(1,jb),1,10)
          ibperm(ib) = jb
          zebeam(ib) = hebeam(jb)
          zibeam(ib) = hibeam(jb)
          zfract(1,ib) = hfract(1,jb)
          zfract(2,ib) = hfract(2,jb)
          zfract(3,ib) = hfract(3,jb)
          zton(ib) = hton(jb)*uest
          ztoff(ib) = htoff(jb)*uest
          zrmaj(ib) = hrmaj(jb)
          zrmin(ib) = hrmin(jb)
          zwidth(ib) = hwidth(jb)
          zheigt(ib) = height(jb)
          zlenth(ib) = hlenth(jb)
          zangle(1,ib) = hangle(1,jb)
          zangle(2,ib) = hangle(2,jb)
          zangle(3,ib) = hangle(3,jb)
          zangle(4,ib) = hangle(4,jb)
          zfocl(ib)  = hfocl(jb)
          zfoclv(ib) = hfoclv(jb)
          zdiv(ib)   = hdiv(jb)
          zdivv(ib)  = hdivv(jb)
          zaper(ib)  = haper(jb)
          zaperv(ib) = haperv(jb)
          i001 = nhshap(jb)
c          call copyi(ihshap(1,i001),1,ihtitl(1,ib),1,5)
          ihtitl(ib) = ihshap(i001)
  104   continue
c
c
        write (nprint,10101) (ibperm(jb),jb=1,ib)
        write (nprint,10102) (ihblnk,jb=1,ib)
c
        write (nprint,10103) (ihname(jb),jb=1,ib)
        write (nprint,10104) (habeam(jb),jb=1,ib)
c
        write (nprint,10105) lhtemp, (zebeam(jb),jb=1,ib)
        write (nprint,10106) lhcurr, (zibeam(jb),jb=1,ib)
        write (nprint,10107) (zfract(1,jb),jb=1,ib)
        write (nprint,10108) (zfract(2,jb),jb=1,ib)
        write (nprint,10109) (zfract(3,jb),jb=1,ib)
c
        write (nprint,10110) (zton(jb),jb=1,ib)
        write (nprint,10111) (ztoff(jb),jb=1,ib)
c
        write (nprint,10112) (zrmaj(jb),jb=1,ib)
        write (nprint,10113) lhlen, (zrmin(jb),jb=1,ib)
c
        write (nprint,10114) (ihtitl(jb),jb=1,ib)
        write (nprint,10115) lhlen, (zwidth(jb),jb=1,ib)
        write (nprint,10116) lhlen, (zheigt(jb),jb=1,ib)
        write (nprint,10117) lhlen, (zlenth(jb),jb=1,ib)
c
        write (nprint,10118) (zangle(1,jb),jb=1,ib)
        write (nprint,10119) (zangle(2,jb),jb=1,ib)
        write (nprint,10120) (zangle(3,jb),jb=1,ib)
        write (nprint,10121) (zangle(4,jb),jb=1,ib)
c
        write (nprint,10122) (zfocl(jb),jb=1,ib)
        write (nprint,10123) lhlen, (zfoclv(jb),jb=1,ib)
c
        write (nprint,10124) (zdiv(jb),jb=1,ib)
        write (nprint,10125) (zdivv(jb),jb=1,ib)
c
        do 114 jb = 1, ib
          i001 = ibperm(jb)
          i002 = nhprof(i001)
c         call copyi(ihprof(1,i002),1,ihtitl(1,jb),1,5)
          ihtitl(jb) = ihprof(i002)
  114   continue
c
        write (nprint,10126) (ihtitl(jb),jb=1,ib)
c
        do 118 jb = 1, ib
          i001 = ibperm(jb)
          i002 = nhprfv(i002)
c         call copyi(ihprof(1,i002),1,ihtitl(1,jb),1,5)
          ihtitl(jb) = ihprof(i002)
  118   continue
c
        write (nprint,10127) (ihtitl(jb),jb=1,ib)
c
c
        do 124 jb = 1, ib
          i001 = ibperm(jb)
          i002 = nhaper(i002)
c         call copyi(ihshap(1,i002),1,ihtitl(1,jb),1,5)
          ihtitl(jb) = ihshap(i002)
  124   continue
c
        write (nprint,10128) (ihtitl(jb),jb=1,ib)
        write (nprint,10129) lhlen, (zaper(jb),jb=1,ib)
        write (nprint,10130) lhlen, (zaperv(jb),jb=1,ib)
c
        return
c
c
cl      first page of timestep printout
c
c
c
  200   continue
        if (nlpomt(17).and.nlpomt(18)) return
        do 500 jspc=1,mhsp
c
c               check if beams have started
c
        do 202 jb = 1, mhbeam
        if(nhbeam(jb).ne.nhbeam(nhbem1(jspc))) go to 202
        if (hton(jb).lt.tai*uiet) go to 203
  202   continue
c
        go to 500
  203   continue
c
c
c
        zt = tai * uist*1000.0
        zdt = dtoldi * uist*1000.0
        ztold = htlast * uist*1000.0
c
        lpage = lpage + 1
        write (nprint,10200) label1(1:48),label5(1:72),
     1           lpage,nstep,zt,zdt,  nilast,ztold,nipart
c
c
c               check if any beams on
c
c
        do 204 jb = 1, mhbeam
        if (libeam(jb).gt.0) go to 205
  204   continue
c
        write (nprint,10202)
        return
c
  205   continue
c
c
c               describe the beams that are on
c
c
c
        z0 = 10.0**(-fxes) / fces * uesi * uesh * usep
c
        ib=0
        do 214 jb = 1, mhbeam
        if (libeam(jb).le.0) go to 214
        if(nhbeam(jb).ne.nhbeam(nhbem1(jspc))) go to 214
        ib = ib+1
        ibperm(ib) = jb
        zener(ib) = hebeam(jb)
c
        do 212 jfr = 1, mxhfr
        zpowr(jfr,ib) =
     1  yfract(ib,jfr) * z0 * hibeam(jb)*hebeam(jb) / float(jfr)
  212   continue
  214   continue
c
        it = 1
        if (ib.gt.1) it = 2
c
        write (nprint,10203) (ihblnk,lhpowr,j=1,it)
        write (nprint,10204) (ihblnk,lhtemp,j=1,it)
        write (nprint,10205)
     1  (ibperm(jb),zener(jb),(zpowr(jf,jb),jf=1,mxhfr),jb=1,ib)
c
c
cl      3)      h(r) for each e
c
c
  300   continue
        if (nlpomt(17)) go to 395
c
        call copyr(yener,1,zener,1,mnengy)
        call scaler(zener,mnengy,evs*useh)
c
        write (nprint,10300) lhtemp, lhlen, (zener(je),je=1,mnengy)
c
        do 318 jc = 1, nyzone
        call resetr(zhr,mnengy,0.0)
c
        do 308 je = 1, mnengy
        do 304 jmu = 1, nymu
        zhr(je) = zhr(je) + hrb(je,jmu,jc,jspc)*hdmuc(jmu)
  304   continue
  308   continue
c
        zr1 = yrmajr(jc+1) * usel
        zr2 = yr(jc+1) * usel
c
        write (nprint,10301) jc, zr1, zr2, (zhr(je),je=1,mnengy)
  318   continue
c
        write (nprint,10302) (yrlose(je,jspc),je=1,mnengy)
c
  395   continue
c
c
cl      4)      total h(r), avg. mu and standard dev. of mu,
c               particle deposition rate, and charge-exchange deposition
c               rate
c
c
c
  400   continue
        if (nlpomt(18)) go to 495
c
        write (nprint,10400) lhlen, lhdenr, lhdenr
c
        zipowr = 0.0
        do 404 je = 1, mnengy
        zipowr = zipowr + yener(je)*hnsrcs(je,jspc)
  404   continue
c
        if (zipowr.gt.epslon) zipowr = 1.0 / zipowr
c
        zvols = avi(mzones,12,1)*uisl**3
        zivols = 1.0 / zvols
        zir2 = 1.0 / yr(nyzone+1)**2
c
        zz1 = 0.0
        zz2 = 0.0
        zz3 = 0.0
        zzdn = 0.0
        zzcex = 0.0
        zloss=0.0
c
c
c
        do 418 jc = 1, nyzone
c
        z1 = 0.0
        z2 = 0.0
        z3 = 0.0
        z4 = 0.0
        z5 = 0.0
c
        do 414 je = 1, mnengy
        z7 = 0.0
        z8 = 0.0
        z9 = 0.0
c
        do 408 jmu = 1, nymu
        z7 = z7 + hrb(je,jmu,jc,jspc)*hdmuc(jmu)
        z8 = z8 + hrb(je,jmu,jc,jspc)*hdmuc(jmu)*hmuc(jmu)
        z9 = z9 + hrb(je,jmu,jc,jspc)*hdmuc(jmu)*hmuc(jmu)**2
  408   continue
c
        z1 = z1 + yener(je)*hnsrcs(je,jspc)*z7
        z2 = z2 + yener(je)*hnsrcs(je,jspc)*z8
        z3 = z3 + yener(je)*hnsrcs(je,jspc)*z9
        z4 = z4 + hnsrcs(je,jspc)*z7
        z5 = z5 + hnsrcs(je,jspc)*z7* (cxxcex(je,jc,jspc) /
     1            clamda(je,jc,jspc))
  414   continue
c
        zhrt = z1 * zipowr
        zloss = zloss + zhrt*(yr(jc+1)**2 - yr(jc)**2)
        zmu = 0.0
        zdev = 0.0
        if (z1.gt.epslon) zmu = z2 / z1
        if (z1.gt.epslon) zdev = sqrt(max(0.0 , z3/z1 - zmu**2))
        zdn = z4 * zivols * used*uest
        zcex = z5 * zivols * used*uest
c
        zr = (yr(jc) + yr(jc+1))*0.5*usel
        if(jc.gt.mhbeam) go to 415
        if(libeam(jc).le.0) go to 415
        write (nprint,10403) jc, zr, zhrt, zmu, zdev, zdn, zcex,
     x  jc, yaumin(jc), yaumaj(jc)
 
        go to 416
415     continue
c
        write (nprint,10401) jc, zr, zhrt, zmu, zdev, zdn, zcex
416     continue
c
        zdx2 = zir2 * (yr(jc+1)**2 - yr(jc)**2)
c
        zz1 = zz1 + z1*zdx2
        zz2 = zz2 + z2*zdx2
        zz3 = zz3 + z3*zdx2
        zzdn = zzdn + z4*zdx2
        zzcex = zzcex + z5*zdx2
c
  418   continue
c
        zloss = 100.*(1. - zloss*zir2)
        if (zz1.gt.epslon) zmu = zz2 / zz1
        if (zz1.gt.epslon) zdev = sqrt(max(0.0 , zz3/zz1 - zmu**2))
        zzdn = zzdn * uest
        zzcex = zzcex * uest
        zdn = zzdn * zivols*used
        zcex = zzcex * zivols*used
c
        write (nprint,10402) zloss, zmu, zdev, zdn, zcex, jspc,
     1   lhtime, zzdn, zzcex
c
  495   continue
500     continue
c
c
c
        return
c
c    .    .    .    .    .    .    .    .    .    .    .    .
c
c
c
cl      10)     short printout
c
c
c
 1000   continue
c
        return
c
c
c
cl      100)    format statements
c
c
c
10100   format(/5x,'the deposition profiles for high-energy hydrogen',
     1  ' beams are computed by a monte-carlo method, using',i8,
     2  ' test'/5x,'particles and',i5,' zones.  the profiles are ',
     3  '(re)computed every time an injector is turned on or off',/5x,
     4  'every',i5,' timesteps, and every time the te profile',
     5  ' changes by more than',2pf6.2,' %, or the mean-free-path'/5x,
     6  'profile changes by more than',2pf6.2,'%.')
10101   format(/5x,'injector number:',2x,8(i11))
10102   format(27x,8(a1,'----------'))
10103   format(27x,8(1x,a10))
10104   format(5x,'atomic weight:',6x,8(0pf11.1))
10105   format(/5x,'energy (',a10,'):',2x,8(1pe11.3))
10106   format(5x,'current (',a10,'):',1x,8(1pe11.3))
10107   format(5x,'fraction: 1/1 energy',2x,8(0pf11.3))
10108   format(15x,'1/2 energy',2x,8(0pf11.3))
10109   format(15x,'1/3 energy',2x,8(0pf11.3))
10110   format(/5x,'on at  (msec.):',7x,8(3pf11.3))
10111   format(5x,'off at (msec.):',7x,8(3pf11.3))
10112   format(/5x,'pivot point: major r.',8(0pf11.2))
10113   format(5x,'(',a10,') minor r.',8(0pf11.2))
10114   format(/5x,'beam shape:',11x,8(1x,5a2))
10115   format(5x,'width  (',a10,'):',1x,8(0pf11.2))
10116   format(5x,'height (',a10,'):',1x,8(0pf11.2))
10117   format(5x,'length (',a10,'):',1x,8(0pf11.2))
10118   format(/5x,'angle      1:',8x,8(0pf11.2))
10119   format(5x,'(degrees)  2:',8x,8(0pf11.2))
10120   format(16x,'3:',8x,8(0pf11.2))
10121   format(16x,'4:',8x,8(0pf11.2))
10122   format(/5x,'focal length: horiz.',1x,8(0pf11.2))
10123   format(5x,'(',a10,')   vert.',1x,8(0pf11.2))
10124   format(/5x,'divergence: horiz.',4x,8(0pf11.3))
10125   format(5x,'(degrees)    vert.',4x,8(0pf11.3))
10126   format(/5x,'beam profile: horiz.',2x,8(1x,5a2))
10127   format(20x,'vert.',2x,8(1x,5a2))
10128   format(/5x,'aperture     shape:',3x,8(1x,5a2))
10129   format(5x,'(',a10,') width:',2x,8(0pf11.2))
10130   format(5x,'(',a10,') height:',1x,8(0pf11.2))
c
10200   format(1h1,2x,a48,10x,a72,
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.'/
     3  10x,'h(r) profile at timestep ',i5,',     time  ',0pf12.3,
     4  '  millisecs., ',i5,' test particles per injector')
10202   format(//20x,'***************',5x,'active injector no. ',5x,
     1  '***************'//1x)
10203   format(/a1,'active  ', 2x,
     1  ' no.   energy  power in energy compon. (',a10,')',
     2  a1,6x,' no.   energy  power in energy compon. (',a10,')')
10204   format(a1,'injectors:',5x,
     1  a10,'  1/1 energy  1/2 energy  1/3 energy',
     2  a1,11x,a10,'  1/1 energy  1/2 energy  1/3 energy')
10205   format((11x,i3,4(2x,1pe10.3),7x,i3,4(2x,1pe10.3)))
c
10300   format(/'  zone    rmajor  rminor',8x,
     1  'h(r) at energy levels (',a10,1h)/11x,a10,10x,1p8e11.2)
10301   format(2x,i4,1x,0p2f8.2,4x,0p20(f11.2:))
10302   format(27x,8(6x,'-----')/'   ** 0ss **',12x,2p8f11.2)
c
10400   format(/2x,'zone',2x,'avg.rad.',2x,'total h(r)',3x,
     1  'cos of pitch ang',1x,'deposition  charge ex.',
     x  14x,'beam no. ',2x,'optical depth to minimum:',
     2  /7x,a10,14x,'mean    std.dev.',1x,a10,2x,a10,
     x  25x,'minor rad.',2x,'major rad.')
10401   format(2x,i4,2x,0pf7.2,0pf11.2,0p2f10.4,1p2e12.3)
10402   format(30x,'------',4x,'------',2(2x,'----------')
     1  /'  ** 0ss / avg. ',0pf6.1,' %',0p2f10.4,1p2e12.3,
     1  15x,'species',i3,
     2  /'  ** total part./',a10,'** ',16x,1p2e12.3)
10403   format(2x,i4,2x,0pf7.2,0pf11.2,0p2f10.4,1p2e12.3,
     x  18x,i2,4x,2(0pf8.4,4x))
c
        end
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@grafix   /baldur/code/bald  file DIO
c  rgb 30-may-96 save inold and set ix0 just before using
c       les   -nov-90  add d-3he fusion  to 'alpha' heating
c              add d-3he fast particles to zalfap and ztbalf
c              add synchrotron rad to total losses; note that
c                   loss from burning of thermal ions is subtracted from
c                   total fusion power, not from losses
c              add d-3he fusion radial plots
c              add synchrotron radiation radial plots
c              add total powers, convection and conduction (from dd3hfus)
c              add total fusion heating, neutron power from d-3he
c              note .6666666 is cutoff in ztbalf in orig BALDUR!!
c             note no synchrotron rad in zloss in original!!
c       dps 10-nov-88 15.07 alter line-averaged electron density to be
c                     consistent with other calculations.
c       dps 12-may-87 use vloopi(jz,2) for loop voltage and Etor,
c                     which is now actually parallel electric field
c       rgb 4-sep-85 1 1/2 D modification of volume, surface area, radius
c                volume within zone bndry j = avi(j,12,1) in internal units
c                surface area of zone bndry j = avi(j,3,1) * avi(j,5,1)
c                minor radius (half width) to zone bndry j = avi(j,15,1)
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c
         subroutine grafix(k)
c
c 3.7 write out variables for radial profiles
c
c  ***note***   in order to obtain correct plotting, a version of
c               "cdrplt" must be used in which all read statements
c               involve only variables that match precisely the write
c               statements of subroutine grafix, both in calling se-
c               quence and in array dimensioning.
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'cbparm.m'
       include 'cd3he.m'
       include 'cmg.m'
       include 'cpower.m'
c
c---------------------------------------------------------------------
c
         dimension
     1   zzni(idxchi)  , zrconf(5)     , ztconf(idchp2,5)        ,
     2   zloss(idxchi) , zlne(6)       , ztbalf(mj),
     3   zetors(2,mj)
c
         data ipts /4/
c
      save inold
c
c---------------------------------------------------------------------
cl               1. initialize
c
  100    continue
         if (nrgraf .le. 0) go to 500
         if (k .ne. 1) go to 200
c
         rewind nrgraf
         inold=-1
c
         write (nrgraf) label1,label2,label3,label4,label5,label6,
     *   lhyd1,lhydn,limp1,limpn,mhyd,mimp,mchi,lhspec
c
c---------------------------------------------------------------------
cl               2. intermediate output
c
  200    continue
         if (inold .eq. nstep) go to 300
         ztime=tai*uist
         zdelt=dti*uist
c
c calculate confinement times
c initial values
         zzei = 0.
         zzee = 0.
         zzeb = 0.
         do 205 jz=lcentr,ledge
         zzee = zzee + chi(lelec,jz) * dx2i(jz)
         zzei = zzei + chi(lion,jz) * dx2i(jz)
         zzeb = zzeb + hebems(jz) * rhobis(2,jz) * dx2i(jz)
  205    continue
         z0 = 2. * avi(mzones,12,1)
         zzee = zzee * z0 * uiee
         zzei = zzei * z0 * uiee
         zzeb = (zzeb * usee) * usid * z0
c
        zvols = avi(mzones,12,1) * uisl**3
         zsurfi = avi(mzones,3,1) * avi(mzones,5,1)  ! plasma surface area
c
         zohme = 0.
         zalfe = 0.
         zexte = 0.
         do 210 jz = lcentr,ledge
         zohme = zohme + dx2i(jz) * weohms(jz)
         zalfe = zalfe + dx2i(jz) * (wealfs(jz) + wialfs(jz) +
     &     wed3fs(jz) + wid3fs(jz))
         zexte = zexte + dx2i(jz) * (webems(jz) + wibems(jz))
  210    continue
         zohme = zohme * zvols * usep
         zalfe = zalfe * zvols * usep
         zexte = zexte * zvols * usep
c
c confinement times
c
        zpts = 1.0 / float(ipts)
        i = 1
        call resetr(zzni  ,mxchi,0.0)
         ix0=5*(mxchi+2)
        call resetr(ztconf,ix0,0.0)
        call resetr(zloss ,mxchi,0.0)
        call resetr(zlne  ,6    ,0.0)
c
        zzlne = 0.0
        do 518 jz = lcentr, ledge
        zeirs = 0.0
        if (mimp.le.0) go to 505
        do 504 ji = 1, mimp
        zeirs = zeirs + weirs(ji,jz)
  504   continue
  505   continue
c
        zloss(lelec) = zloss(lelec) + usip*zvols*dx2i(jz)*
     1        (weions(jz) + webrs(jz) + zeirs + wesrs(jz) + wesyn(jz))
        zloss(lion) = zloss(lion) - usip*zvols*dx2i(jz)*
     1                          (wiions(jz) + wichxs(jz))
c
c
c
        z0 = zvols * uisd * dx2i(jz)
        do 512 jp = 1, mchi
        zzni(jp) = zzni(jp) + z0*chi(jp,jz)
  512   continue
c
        zzlne = zzlne + rhoels(2,jz) * dx2i(jz) * 2.
c
        if (jz.lt.ledge.and.(xbouni(jz+1)+rndeps).lt.(float(i)*zpts))
     1          go to 518
c
        zrconf(i) = avi(jz+1,15,1) * uiel
        zeloss = 0.0
        zlne(i) = zzlne
c
        do 516 jp = 1, mchi
c
        zflux = 0.0
        do 514 jp2 = 1, mchi
        zflux = zflux - aaaa(jp,jp2,jz+1)*chi(jp2,jz)
     1                - bbbb(jp,jp2,jz+1)*chi(jp2,jz+1)
  514   continue
c
        z0 = (zloss(jp) + zflux*zsurfi*xbouni(jz+1)) * ueit
        if (z0.gt.epslon) ztconf(jp,i) = zzni(jp) / z0
        if (jp.eq.lelec.or.jp.eq.lion) zeloss = zeloss + z0
  516   continue
c
        if (zeloss.gt.epslon) ztconf(mchi+1,i) =
     1  (zzni(lelec) + zzni(lion)) / zeloss
        zlne(i) = zlne(i) * ztconf(mchi+1,i) / xbouni(jz+1)**2
        i = i + 1
  518   continue
c
c
c       experimental times
c
        ztenr1 = 0.0
        if (zohme.gt.epslon)
     1          ztenr1 = (zzei + zzee) / zohme * usep*uset*uese
        ztenr2 = 0.0
        z0 = zohme + zexte + zalfe
        if (z0.gt.epslon)
     1          ztenr2 = (zzei + zzee + zzeb) / z0 * usep*uset*uese
c
c compute alpha beta, line average density
c  15.07 calculate line average density as done elsewhere; note that
c  we explicitly prevent a contribution from the outer ghost zone.
c
         zlined = 0.0
         zalfap = 0.0
      do 250 jz = lcentr,mzones
         zdrs = ahalfs(jz+1,1) - ahalfs(jz,1)
         if (jz.eq.mzones) zdrs = 0.
         zlined = zlined + rhoels(2,jz) * zdrs
         zalfap = zalfap + (alphai(jz) * ealfai(jz) +d3fast(jz)) *
     &     dx2i(jz)
         ztbalf(jz) = (d3fast(jz) +
     &     alphai(jz) * ealfai(jz)) * uisd * uise * .666666666
     &                *8. * fcpi / bzs**2
  250 continue
c
c  15.07 in particular, we must allow for use of the scrape-off; note
c  that the volume average density is computed externally and does not
c  account for use of the scrape-off.
c
      if (nadump(1).gt.lcentr) then
        isep = nadump(1)
      else
        isep = mzones
      end if
      zlined = zlined / ahalfs(isep,1)
c
         zalfap = zalfap * uisd * uise * .666666666 * 2.
     *            * 8. * fcpi / bpols(1,mzones)**2
         zalfat = zalfap * (bpols(1,mzones)/bzs)**2
c
c       loop voltage and toroidal electric field
c         note that zetors is now actually parallel electric field
c
        do 255 jz=lcentr,mzones
         zetors(1,jz)=0.25*(vloopi(jz,2)+vloopi(jz-1,2))*(avi(jz,10,1)
     1               +avi(jz-1,10,1))/(2.0*fcpi*avi(jz,11,1))*uisv
         zetors(2,jz)=vloopi(jz,2)*avi(jz,10,1)/(fcpi*(avi(jz,11,1)
     1               +avi(jz+1,11,1)))*uisv
255     continue
         zv = vloopi(ledge,2)
c
c
         write (nrgraf) nstep,zdelt
c
c   les  nov-90 include d3he fusion
c
         write (nrgraf) ztime,lcentr,ledge,mzones,xzoni,rmajs,rmins,
     *   ajzs,rhoins,rhoels,tes,tis,xzeff,q,zetors,wealfs,webrs,webems,
     *   weohms,weirs,wialfs,wibems,cmean,rhohs,rhois,tns,rhons,
     *   wed3fs,weash,wesyn,wid3fs,wid3fl,d3fast,
     *   detepr,ditipr,dnhs,dnis,vxcmg,vnwars,
     *   bzs,bpols(1,mzones),hebems,rhobis,dx2i,ztconf,ztenr2,zlne,
     *   zalfap,zalfat,zlined,zv,ztbalf
c
c---------------------------------------------------------------------
cl               3. close files
c
  300    continue
         inold=nstep
         if (k .ne. 3) go to 500
c
         zdelt=-5.
         write (nrgraf) nstep,zdelt
c
         endfile nrgraf
c
c
  500    return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@tgraf   /baldur/code/bald  file DIO 
c  rgb 30-may-96 save inold
c
         subroutine tgraf(k)
c
c 3.8 write out variables for time profiles
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'cpower.m'
c
      save inold
c
c---------------------------------------------------------------------
cl               1. initialize
c
  100    continue
         if (ntgraf .le. 0) go to 500
         if (k .ne. 1) go to 200
c
         rewind ntgraf
         inold=-1
c
         if (nrgraf .gt. 0) go to 200
         write (ntgraf) label1,label2,label3,label4,label5,label6,
     *   lhyd1,lhydn,limp1,limpn,mhyd,mimp
c
c---------------------------------------------------------------------
cl               2. intermediate output
c
  200    continue
         if (inold .eq. nstep) go to 300
         ztime=tai*uist
         zdt=dti*uist
c
         write (ntgraf) ztime,tes(2,lcentr),tis(2,lcentr),
     *   rhoels(2,lcentr),rhoins(2,lcentr),xzeff(2,lcentr),zdt,
     *   poht,pbrem,psyn,pirad,pfust,pfusd3,pfusdt,plosst,pauxt,
     *   pcond,pconde,pcondi,pconve,pconvi,pntn
c
c---------------------------------------------------------------------
cl               3. close file
c
  300    continue
         inold=nstep
         if (k .ne. 3) go to 500
c
         ztime=-5.
         write (ntgraf) ztime,ztime,ztime,ztime,ztime,ztime,ztime,
     *   ztime,ztime,ztime,ztime,ztime,ztime,ztime,ztime,ztime,
     *   ztime,ztime,ztime,ztime,ztime,ztime
c
         endfile ntgraf
c
c
  500    return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@fprint   /baldur/code/bald  file DIO
c  rgb 05-jun-89 commented out unused format statements
c       rgb 4-sep-85 1 1/2 D modification of volume, surface area, radius
c                volume within zone bndry j = avi(j,12,1) in internal units
c                surface area of zone bndry j = avi(j,3,1) * avi(j,5,1)
c                minor radius (half width) to zone bndry j = avi(j,15,1)
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c*************************************************************
        subroutine fprint(ktype,k)
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        common/comtst/zll(8,mj),zione(8,mj),zelece(8,mj)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
        dimension zwealf(mj),zwialf(mj),zsh1(mj),zsh2(mj),
     1          zsum(12),ztot(12),zvols(mj)
c
        logical ilend
c
        data  iclass/3/,    isub/9/
c
        if(.not.nlomt3(isub)) go to 10
        return
c
 10     continue
c
        if(ktype.gt.0) return
        if(k.ne.1) go to 200
c
 100    continue
c
        return
c
c       2.)  print fusion page of large time step edit
c
 200    continue
c
        if(nlpomt(19)) return
c
        ihyd = 0
        ideut = 0
        ihe3 = 0
        ihe4 = 0
c
        do 205 ih = 1,mhyd
        if(ngas(ih).eq.1) ihyd = ih
        if(ngas(ih).eq.-2) ideut = ih
 205    continue
c
        if(mimp.eq.0) go to 215
c
        do 210 ii = 1,mimp
        if(nimp(ii).eq.-4) ihe3 = ii
        if(nimp(ii).eq.2) ihe4 = ii
 210    continue
c
 215    continue
c
c
c               if no fusions don't print out
c
        if(ideut*ihe3.eq.0) return
c
        zshfus = 0.
c
c
c
        do 220 jz = lcentr,ledge
      zvols(jz) = ( avi(jz+1,12,1) + avi(jz,12,1) ) * uisl**3
        zshfus = zshfus + (shfus(ideut,jz)+shfus(ihe3,jz))*zvols(jz)
 220    continue
c
        zdts=dti*uist
c
        zshfus = zshfus*zdts
c
        if(zshfus.le.1.) return
c
c--------------------------------------------------------------
c
c               page 1
c
        lpage=lpage+1
        zt=tai*uist*1000.
        zdt=dtoldi*uist*1000.
c
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write(nprint,10201)
c
c               loop begins here
c
        i1 = nskip + 1
c
        do 222 jz = i1,ledge,nskip
c
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
        ijz = jz - 1
        zwealf(jz) = wealfs(jz)*usee
        zwialf(jz) = wialfs(jz)*usee
        zsh1(jz) = -shfus(ideut,jz)
        zsh2(jz) = aoloss(jz)
c
        write(nprint,10202)ijz,zrad,afuses(jz),aslows(jz),
     1          atfuse(jz),zwealf(jz),zwialf(jz),zsh1(jz),zsh2(jz)
c
 222    continue
c
c               loop has ended
c
c               compute and print totals
c
        write(nprint,10203)
        call resetr(zsum,7,0.)
        call resetr(ztot,7,0.)
c
        do 224 jz=lcentr,mzones
        zsum(1)=zsum(1)+afuses(jz)
        zsum(2)=zsum(2)+aslows(jz)
        zsum(3)=zsum(3)+atfuse(jz)
        zsum(4)=zsum(4)+zwealf(jz)
        zsum(5)=zsum(5)+zwialf(jz)
        zsum(6)=zsum(6)+zsh1(jz)
        zsum(7)=zsum(7)+zsh2(jz)
c
 224    continue
c
        zdenom=float(mzones-lcentr)+1.
c
        do 226 ind=1,7
        zsum(ind)=zsum(ind)/zdenom
 226    continue
c
        write(nprint,10204)(zsum(j),j=1,7)
c
        do 228 jz=lcentr,ledge
        z1=rhohs(ideut,2,jz)*rhois(ihe3,2,jz)*afuses(jz)
        z2=.5*rhohs(ideut,2,jz)*rhohs(ideut,2,jz)*aslows(jz)
        z3=.5*rhohs(ideut,2,jz)*rhohs(ideut,2,jz)*atfuse(jz)
        ztot(1)=ztot(1)+zvols(jz)*z1
        ztot(2)=ztot(2)+zvols(jz)*z2
        ztot(3)=ztot(3)+zvols(jz)*z3
        ztot(4)=ztot(4)+zvols(jz)*zwealf(jz)
        ztot(5)=ztot(5)+zvols(jz)*zwialf(jz)
        ztot(6)=ztot(6)+zvols(jz)*zsh1(jz)
        ztot(7)=ztot(7)+zvols(jz)*zsh2(jz)
 228    continue
c
        do 230 ind=1,7
        ztot(ind)=ztot(ind)*dti*uist
 230    continue
c
        write(nprint,10206)(ztot(j),j=1,7)
c
        abouni(1,1)=abouni(1,1)+ztot(1)
        abouni(2,1)=abouni(2,1)+ztot(2)
        abouni(3,1)=abouni(3,1)+ztot(3)
c
        write(nprint,10205)(abouni(i,1),i=1,3)
c
c------------------------------------------------------
c
c               page 2
c
        lpage=lpage+1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
        write(nprint,10300)
c
        i1 = nskip+1
c
c               loop begins here
c
        do 240 jz = i1,ledge,nskip
c
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
        ijz = jz-1
c
        write(nprint,10301) ijz,zrad,(zll(j,jz),j=1,6)
c
 240    continue
c
c               loop ends
c
c---------------------------------------------------------
c
c               page 3
c
        lpage=lpage+1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write(nprint,10400)
c
c               loop begins
c
        do 250 jz = i1,ledge,nskip
c
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
        ijz = jz-1
        zzz = zione(7,jz)+zione(8,jz)
c
        write(nprint,10401)ijz,zrad,(zione(j,jz),j=1,6),zzz
c
 250    continue
c
c               loop ends
c
        write(nprint,10203)
        call resetr(zsum,7,0.)
        call resetr(ztot,7,0.)
c
c               compute and print totals
c
        do 254 jz=lcentr,ledge
c
        do 252 ind=1,6
        zsum(ind)=zsum(ind)+zione(ind,jz)
        ztot(ind)=ztot(ind)+zvols(jz)*zione(ind,jz)*zdts
 252    continue
        zsum(7)=zione(7,jz)+zione(8,jz)+zsum(7)
        ztot(7)=ztot(7)+zvols(jz)*(zione(7,jz)+zione(8,jz))*zdts
 254    continue
c
        do 256 ind=1,7
        zsum(ind)=zsum(ind)/zdenom
 256    continue
c
        write(nprint,10302)(zsum(j),j=1,7)
c
        write(nprint,10303)(ztot(j),j=1,7)
c
c----------------------------------------------------
c
c               page 4
c
        lpage=lpage+1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
c
        write(nprint,10500)
c
c               loop begins here
c
        do 260 jz = i1,ledge,nskip
        zrad = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
        ijz = jz - 1
        zzw = zelece(7,jz)+zelece(8,jz)
c
        write(nprint,10401)ijz,zrad,(zelece(j,jz),j=1,6),zzw
 260    continue
c
c               loop ends
c
c
        write(nprint,10203)
        call resetr(zsum,7,0.)
        call resetr(ztot,7,0.)
c
c               compute and print totals
c
        do 264 jz=lcentr,ledge
c
        do 262 ind=1,6
        zsum(ind)=zsum(ind)+zelece(ind,jz)
        ztot(ind)=ztot(ind)+zvols(jz)*zelece(ind,jz)*zdts
 262    continue
        zsum(7)=zelece(7,jz)+zelece(8,jz)+zsum(7)
        ztot(7)=ztot(7)+zvols(jz)*(zelece(7,jz)+zelece(8,jz))*zdts
 264    continue
c
        do 266 ind=1,7
        zsum(ind)=zsum(ind)/zdenom
 266    continue
c
        write(nprint,10302)(zsum(j),j=1,7)
c
        write(nprint,10303)(ztot(j),j=1,7)
c
        return
c
c-----------------------------------------------------
c
c               format statements
c
c10100   format(1h1,4x,'zone',4x,'radius',7x,'te',11x,'nhe3',10x,
c     1          'nd'/14x,'cm',8x,'kev',10x,'part/cc',4x,'part/cc')
c10101   format(6x,i2,3x,f7.2,3(3x,e10.3))
10200   format(1h1,2x,a48,10x,a72//
     1          2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     2          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     3          'millisecs.',/)
10201   format(/4x,'zone',4x,'radius',6x,
     1  'thermal-fusion rates(part/cc*s)',6x,
     2  'energy rate(joules/cc*s)',2x,
     3  'source/sink terms(part/cc*s)'/
     4  13x,'cm',9x,'ii',11x,'iii',10x,'iv',9x,
     5  'electrons',7x,'ions',8x,'deuterium',6x,'he3')
10202   format(6x,i2,3x,f7.2,5(3x,e10.3),2x,2(3x,e10.3))
10203   format(12x,4('-------------------------'))
10204   format(1x,'** averages **',6x,5(e10.3,3x),2x,2(e10.3,3x))
10205   format(1h1,1x,'** total no. of reactions since t=0:'/
     1  /5x,'ii: d + he3 -> he4 + p',15x,e10.3/
     2  5x,'iii: 3d -> p + he4 + n + 1.01mev',5x,e10.3/
     3  5x,'iv: d + d -> he3 + n',17x,e10.3)
10206   format(36x,'-- reactions --',18x,
     1          '-- joules --',12x,'-- particles --'/
     2       1x,'** totals **',8x,5(e10.3,3x),2x,2(e10.3,3x))
10300   format(44x,
     1  'fraction of energy given to ions'//5x,'zone',4x,'radius',
     2  6x,'he4(3.5)mev)',5x,'p(3.03mev)',6x,'he3(.82mev)',4x,
     3  'he4(3.67mev)',4x,'p(14.67mev)',4x,'t(1.01mev)')
10301   format(6x,i2,3x,f7.2,6(9x,f7.3))
10302   format(1x,'** averages **',6x,7(e11.3,3x))
10303   format(55x,'-- joules --'/
     2       1x,'** totals **',8x,7(e11.3,3x))
10400   format(39x,'energy rate to ions(joules/sec)'//
     1  4x,'zone',4x,'radius',4x,'he4(3.52mev)',
     2  3x,'p(3.03mev)',3x,'he3(.82mev)',3x,
     3  'he4(3.67mev)',3x,'p(14.67mev)',4x,'t(1.01mev)',7x,'other')
10401   format(6x,i2,3x,f7.2,7(4x,e10.3))
10500   format(40x,'energy rate to electrons(joules/sec)'//
     1  4x,'zone',4x,'radius',4x,'he4(3.52mev)',
     2  3x,'p(3.03mev)',3x,'he3(.82mev)',3x,
     3  'he4(3.67mev)',2x,'p(14.67mev)',3x,'t(1.01mev)',7x,'other')
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@gprint   /baldur/code/bald  file DIO
c  rgb 30-may-96 set isurf = nsurf - 1 before do 278
c    if ( nsurf .lt. 1 ) return
c  rgb 28-jan-95 used max ( 1, min ( 3, ... ) ) for index of ihtype
c  rgb 20.26 18-jan-92 in write (nprint,10101) change lhgas to lhgas(1)
c  rgb 05-jun-89 put label 490 after go to 499
c     and commented out unused format statement 10140
c       dps 23-oct-86 add ypa array to pellet summary printout
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c******************************************************************************
c
c
        subroutine gprint(kunit,kcall)
c
c
cl      3.4     neutral gas printouts
c
c
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'cmonte.m'
c
c
        dimension
     i  j0(9),          izonet(4),
     i  iperm(6),       isrc(4),        intrvl(20),
     r  zifl(4,4),      zperm(6),       ztifl(4,4),     zofl(4),
     r  znhtot(4),      zn0(4),         zn0tot(4),      zs0(4),
     r  zs0tot(4),      zt0(4),         zt0tot(4),
     r  zsvcfe(9),      zsvcfp(9),      ztboun(5),      zpboun(5),
     r  zspect(100,5,4),zenrgy(100),    zyvalu(100,4),  ztimax(4,4),
     r  znorml(4,4),    zentot(5,4),    zcorel(4,4),    znhall(mj),
     r  zn0all(mj),     zntrvl(20),     zperc(20),      zpero(20)
c
        logical
     l  ilfeed
c
        character *1  ihblnk, ihs, ih0
        character *10 iheu, ihpart, ihtit1(2)
        character *6  iharnt(2)
        character *5  ihdont(2), ihisnt(2)
        character *14 ihmono(2)
        character *10 ihspec(4), ihtype(3), ihpel, ihpelc
c
        data ihblnk /' '/, ihs /'s'/,
     &       ihpart /'part/cm**3'/, iheu /'erg/cc-sec'/,
     &       ihtit1 /'ch.ex.loss','variation '/,
     &       iharnt /'aren''t',' are  '/,
     &       ihdont /'don''t',' do  '/,
     &       ihisnt /'isn''t',' is  '/,
     &       ihmono /'  maxwellian  ','mono-energetic'/
     &       ihtype/'hydrogen  ','deuterium ','tritium   '/
c
        data ipellt,(intrvl(i),i=1,20),ihance,noprnt
     i      /230,231,232,233,234,235,236,237,238,239,240,
     i       241,242,243,244,245,246,247,248,249,250,251,252/
c       for definitions of ipellt, intrvl, ihance, and noprnt
c       see subroutine pdrive
c
c                       data for neutral spectrum business
        data j0/0,1,2,3,4,5,6,7,8/
        data zsvcfe/-0.3173850e02, 0.1143818e02 ,-0.3833998e01 ,
     1               0.7046692e00,-0.7431486e-01, 0.4153749e-02,
     2              -0.9427185e-04,0.           , 0.            /
        data zsvcfp/-0.1490861e03, 0.7592575e02 ,-0.2209281e02 ,
     1               0.3909709e01,-0.4402168e00 , 0.3209047e-01,
     2              -0.1493409e-02,0.4094151e-04,-0.5069777e-06 /
c
c------------------------------------------------------------------------------
c
c
        data    iclass /3/,     isub /4/
c
c
        if (.not.nlomt3(isub)) go to 10
        call mesage(' *** 3.4 subroutine gprint bypassed ')
        return
   10   continue
c
c------------------------------------------------------------------------------
c
c
cl      common blocks and variables modified:
c
c       lpage (comout)
c
c
c------------------------------------------------------------------------------
c
c
        if ( nsurf .lt. 1 ) return
c
        if (kunit.gt.0) go to 1000
        if (kcall.gt.1) go to 200
c
c
cl      1)      initial large printout
c
c
  100   continue
c
        write (nprint,10100) ngprof, ngpart, ngzone, ngsplt
c
        if (nlglim) then
          write (nprint,10102) grecyc
        else
          zt = tai * uist
          call timint (zt,zztcold,bdtime,20,tcold,1,1)
          write (nprint,10101) grecyc, zztcold, lhgas(1)
        endif
c
        iref = 1
        if (nlgref) iref = 2
        ispt = 1
        if (nlgspt) ispt = 2
        iinc = 1
        do 108 ji = limp1, limpn
          if (nzspec(ji).eq.26) iinc = 2
  108   continue
c
        write (nprint,10103) iharnt(iref),
     1                  ihdont(ispt), ihisnt(iinc)
c
c               gas feed
c
        ilfeed = .false.
        do 118 jt = 1, mxt
        do 114 jh = 1, mhyd
        if (gflowi(jh,jt).gt.epslon) ilfeed = .true.
  114   continue
c
        if ((jt.gt.1.and.gtflwi(jt).le.epslon).or.ilfeed) go to 120
  118   continue
c
  120   continue
c
        if (.not.ilfeed) go to 139
c
        imono = 1
        if (nlgmon) imono = 2
c
        zt = tai*uist
        call timint (zt,ztcoldp,bdtime,20,tcoldp,1,1)
        write (nprint,10120) ihmono(imono), ztcoldp, lhtemp
c
c
cl              influx printout
c
c
        i = 0
        do 128 js = lhyd1, limpn
c
        do 124 jt = 1, mxt
          if (jt.gt.1.and.gtflwi(jt).le.0.0) go to 126
          if (gflowi(js,jt).le.epslon) go to 124
          i = i + 1
          iperm(i) = js
cbate        call copyi(lhspec(1,js),1,ihspec(1,i),1,10)
          ihspec(i) = lhspec(js)
          go to 126
  124   continue
  126   continue
  128   continue
c
c
        write (nprint,10121) lhnflx,(ihspec(js),js=1,i)
c
        do 138 jt = 1, mxt
        if (jt.gt.1.and.gtflwi(jt).le.0.0) go to 139
        zt = gtflwi(jt)*uist*1000.0
c
        do 134 js = 1, i
        i001 = iperm(js)
        zperm(js) = gflowi(i001,jt) * ueil**2 * ueit
  134   continue
c
        write (nprint,10122) zt, (zperm(j),j=1,i)
  138   continue
c
  139   continue
c
c               density monitoring
 140    continue
c
        if(gftime(1).ge.1.e10.or.npuff.eq.0) go to 180
c
        do 142 jh = 1,mhyd
          igpuf = max ( 1, min ( 3, iabs(ngas(jh)) ) )
          lhgas(jh) = ihtype(igpuf)
 142    continue
c
        go to (145,146),npuff
c
 145    write(nprint,10130)
        if(denmon(1).le.epslon) go to 150
        go to 160
c
 146    continue
c
        write(nprint,10131)
        if(denmon(1).lt.epslon) go to 150
        go to 160
c
 150    continue
c
        write(nprint,10132) lhgas(1),lhgas(2)
c
        itimes=1
c
 151    continue
c
        if(itimes.eq.11) go to 170
        if(gftime(2*itimes-1).ge.1.e10) go to 170
        zftime=gftime(2*itimes)
        if(gftime(2*itimes).ge.1.e10) zftime=tmax
        write(nprint,10133)gftime(2*itimes-1),
     1          zftime,gfract(itimes,1),gfract(itimes,2),gflmax(itimes)
        itimes=itimes+1
        go to 151
c
 160    continue
c
        write(nprint,10135) lhgas(1),lhgas(2)
        itimes=1
        if(gftime(itimes).gt.epslon) go to 165
c
        zzla = 0.
        isep = mzones
        if ( nadump(1) .gt. lcentr ) isep=nadump(1)
        isepm1 = isep - 1
c
        if ( npuff .eq. 2 ) then
c
          do 163 jz = lcentr,isepm1
            zzla = zzla + rhoels(2,jz)*dx2i(jz)
 163      continue
c
          zzla = 2. * zzla / xbouni(isep)**2
c
        else
c
          do 161 jz = lcentr,ledge
            zzla = zzla + rhoels(2,jz)*dxzoni(jz)
 161      continue
c
          zzla = zzla / xbouni(isep)
c
        endif
c
        write(nprint,10136)gftime(itimes),zzla,
     1          gfract(itimes,1),gfract(itimes,2),gflmax(itimes)
c
        itimes = itimes + 1
c
 165    continue
c
        if(itimes.eq.21) go to 170
        if(gftime(itimes).ge.1.e10) go to 166
        write(nprint,10136)gftime(itimes),denmon(itimes),
     1          gfract(itimes,1),gfract(itimes,2),gflmax(itimes)
c
        itimes=itimes + 1
        go to 165
c
 166    continue
c
        itz = itimes-1
        if(gftime(itz).ge.(tmax-epslon)) go to 170
        write(nprint,10136)tmax,denmon(itz),gfract(itz,1),
     1          gfract(itz,2),gflmax(itz)
c
c
 170    continue
c
        imono=1
        if(nlgmon)imono=2
        zt = tai*uist
        call timint (zt,ztcoldp,bdtime,20,tcoldp,1,1)
        write(nprint,10120) ihmono(imono),ztcoldp,lhtemp
c
c       pellet-injection title page
c
 180    if(npelga(1).eq.0) go to 1900
        if(cfutz(ipellt).le.epslon) cfutz(ipellt)=0.0
        do 1810 it=1,20
        iii=intrvl(it)
        if(cfutz(iii)   .le.epslon) cfutz(iii)=0.0
 1810   continue
        if(cfutz(ihance).le.epslon) cfutz(ihance)=50.0
        if(cfutz(noprnt).le.epslon) cfutz(noprnt)=0.0
        do 1830 it=1,20
        zntrvl(it)=10000.
        if(cfutz(ipellt).le.epslon) go to 1820
        iii=intrvl(1)
        if(cfutz(ipellt).gt.1.1) iii=intrvl(it)
        zntrvl(it)=cfutz(iii)*1000.
 1820   continue
 1830   continue
        write(nprint,9000)
 9000 format(1h1//5x,'pellet injections:'/
     * /5x,'pellet injections are set up to occur in clusters involving'
     1    ,' one or more consecutive pellet firings.  pellet properties'
     2 /5x,'are uniform in each cluster, but can be changed from cluste'
     3    ,'r to cluster.  the number of firings per cluster equals the'
     4 /5x,'difference in "on-time" for this cluster and the next divid'
     5    ,'ed by a time interval, which is the maximum of a prescribed'
     6 /5x,'firing interval (secs) and the latest timestep (dti).  when'
     7    ,' the prescribed firing interval is less than "dti", a depo-'
     8 /5x,'sition-enhancement factor is automatically applied; namely,'
     9 /36x,'(enhancement factor) = (dti)/(prescribed firing interval)'
     * /5x,'firing intervals are prescribed by "cfutz" factors; see com'
     1    ,'ment statements in subroutine pdrive.  "on-times" displayed'
     2 /5x,'below correspond to the release time for the first pellet '
     3    ,'firing in a cluster.'//6x,'on',8x,'hydrogen',14x,'core'
     4 ,18x,'pellet',4x,'core',6x,'pellet',9x,'pellet path',5x,'firing'
     5 /5x,'time',7x,'species(%)',12x,'species(%)',12x,'radius'
     6 ,4x,'radius',4x,'velocity',5x,'radius',2x,'height',4x,'interval'
     7 /5x,'(sec)',50x,'(cm)',6x,'(cm)',6x,'(cm/sec)'
     8 ,6x,'(cm)',4x,'(cm)',4x,'(millisecs)'/)
c
        it=1
 1840   ipelg  = max ( 1, min ( 3, iabs(npelga(it)) ) )
        ipelgc = max ( 1, min ( 3, iabs(npelgc(it)) ) )
c
        ihpel = ihtype(ipelg)
        ihpelc = ihtype(ipelgc)
c
        zpero(it)=100.0*(1.0 - frcor(it))
        zperc(it)=100.0*(1.0 - frout(it))
        write(nprint,9001)tpela(it),ihpel,zpero(it),
     &   ihpelc,zperc(it),
     &   rpela(it),rpelc(it),vpela(it),rpa(it),ypa(it),zntrvl(it)
 9001   format(4x,f7.4,5x,a10,'(',f5.1,')',5x,a10,'(',f5.1,')',5x,
     1  f5.3,6x,f5.3,4x,1pe9.2,4x,0pf5.1,4x,0pf5.1,4x,0pf8.3)
        it=it+1
        if(npelga(it).ne.0.and.it.le.20) go to 1840
 1900   continue
        write (nprint,10123)
c
c
        return
c
c
c
c
cl      2)      large printout (intermediate)
c
c
c
c
  200   continue
        if (rplsma.le.0.0) return
        if (nlpomt(11).and.nlpomt(12)) return
c
c
cl      2.1)    first lines
c
c
        lpage = lpage + 1
        zt = tai * uist * 1000.0
        zdt = dtoldi * uist * 1000.0
        ztold = gtprfi * uist * 1000.0
        zzt = tai * uist
        call timint (zzt,zztcold,bdtime,20,tcold,1,1)
        ztcold = zztcold * uesh * evsinv
c
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt,ngsprf, ztold, ngpart, ztcold
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
cl      2.2)    rsurf, te, ti, ne, nh, wchex, variation by zone, and avg
c
c
  220   continue
        if (nlpomt(11)) go to 250
        write (nprint,10201)
     &         (lhspec(js),js=lhyd1,lhydn),ihtit1(1),ihtit1(2)
        write (nprint,10202) (ihpart,js=lhyd1,lhydn), iheu
c
        ztetot = 0.0
        ztitot = 0.0
        znetot = 0.0
        call resetr(znhtot,ngases,0.0)
        zcxtot = 0.0
        zvtot  = 0.0
c
c
        isurf = nsurf - 1
        do 238 js = 1, isurf
c
c               charge exchange loss (w. ti and nh from old values)
c
        znh = 0.0
        do 222 jg = 1, ngases
        znh = znh + denion(jg,js)
  222   continue
c
        zchex = 0.0
        zchexi = 0.0
        do 224 jsrc = 1, nsrces
        zchex = zchex - fluxn(jsrc)*spi(js,jsrc)
        zchexi = zchexi - fluxn(jsrc)*spii(js,jsrc)
  224   continue
c
        zchex = (zchex + zchexi*tiin(js)) * znh * evs * 1.5
c
c               variation (crude estimate)
c
        zvar = 0.0
        znn = 0.0
        do 228 jg = 1, ngases
        zv = 0.0
        zng = 0.0
        inpt = 0
c
        do 226 jsrc = 1, nsrces
        zv = zv + (fluxn(jsrc) * den0(jg,js,jsrc))**2 *
     1  (sigma(jg,js,jsrc)**2 + 1.0/float(npts(jsrc)))
        zng = zng + fluxn(jsrc)*den0(jg,js,jsrc)
        inpt = inpt + npts(jsrc)
  226   continue
c
        if (zng.gt.epslon) zv = sqrt(zv/zng**2 - 1.0/float(inpt))
        znn = znn + denion(jg,js)
        zvar = zvar + zv * denion(jg,js)
  228   continue
c
        if (znn.gt.epslon) zvar = zvar / znn
        if (znn.le.epslon) zvar = 0.0
c
c               is this a splitting surface?
c
        ih0 = ihblnk
        do 230 jsplt = 1, maxrad
        if (rsplit(jsplt).le.epslon) go to 231
        if ((rsplit(jsplt).gt.(rsurf(js+2)+rsurf(js+1))*0.5).or.
     1  (rsplit(jsplt).le.(rsurf(js)+rsurf(js+1))*0.5)) go to 230
        ih0 = ihs
        go to 231
  230   continue
c
  231   continue
c
c
c
        write (nprint,10220) js, rsurf(js), rsurf(js+1), ih0, tein(js),
     1          tiin(js), denein(js), (denion(jg,js),jg=1,ngases),
     2          zchex, zvar
c
c               totals
c
        zdv = rsurf(js+1)**2 - rsurf(js)**2
        ztetot = ztetot + denein(js)*tein(js)*zdv
        ztitot = ztitot + znn * tiin(js)*zdv
        znetot = znetot + denein(js)*zdv
c
        do 232 jg = 1, ngases
        znhtot(jg) = znhtot(jg) + denion(jg,js)*zdv
  232   continue
c
        zcxtot = zcxtot + zchex*zdv
        zvtot = zvtot + zvar*zdv
c
  238   continue
c
c
c               print averages
c
c
        znn = 0.0
        do 242 jg = 1, ngases
        znn = znn + znhtot(jg)
  242   continue
c
        if (znetot.gt.epslon) ztetot = ztetot / znetot
        if (znn.gt.epslon) ztitot = ztitot / znn
c
        zv = 1.0 / rsurf(nsurf)**2
        znetot = znetot*zv
c
        do 244 jg = 1, ngases
        znhtot(jg) = znhtot(jg)*zv
  244   continue
c
        zcxtot = zcxtot * zv
        zvtot = zvtot * zv
c
c
c
        write (nprint,10221) (ihblnk,jg=1,ngases)
c
        write (nprint,10222) ztetot, ztitot, znetot,
     1                  (znhtot(jg),jg=1,ngases), zcxtot
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
cl      2.3)    n0, t0, and ionization for each h, by zone and avg.
c
c
c
  250   continue
        if (nlpomt(12)) go to 300
c
c               title line
c
        write (nprint,10250) (ihblnk,lhspec(jg),jg=lhyd1,lhydn)
        write (nprint,10251) (ihblnk,jg=lhyd1,lhydn)
c
c
        call resetr(zn0tot,ngases,0.0)
        call resetr(zt0tot,ngases,0.0)
        call resetr(zs0tot,ngases,0.0)
c
        isurf = nsurf - 1
        do 278 js = 1, isurf
c
        call resetr(zn0,ngases,0.0)
        call resetr(zt0,ngases,0.0)
        call resetr(zs0,ngases,0.0)
        zdv = rsurf(js+1)**2 - rsurf(js)**2
c
        znh = 0.0
        do 252 jg = 1, ngases
        znh = znh + denion(jg,js)
  252   continue
c
        do 258 jg = 1, ngases
c
        do 256 jsrc = 1, nsrces
        zn0(jg) = zn0(jg) + fluxn(jsrc)*den0(jg,js,jsrc)
        zt0(jg) = zt0(jg) +
     1          fluxn(jsrc)*den0(jg,js,jsrc)*eneut(jg,js,jsrc)
        zs0(jg) = zs0(jg) + fluxn(jsrc)*
     1          (sne(jg,js,jsrc)*denein(js) + sni(jg,js,jsrc)*znh)
  256   continue
c
c
        zn0tot(jg) = zn0tot(jg) + zn0(jg)*zdv
        zt0tot(jg) = zt0tot(jg) + zt0(jg)*zdv
        zs0tot(jg) = zs0tot(jg) + zs0(jg)*zdv
c
        if (zn0(jg).gt.epslon) zt0(jg) = zt0(jg) / zn0(jg)
        if (zn0(jg).le.epslon) zt0(jg) = 0.0
  258   continue
c
c
        zrad = (rsurf(js+1) + rsurf(js)) * 0.5
c
        write (nprint,10255) js, zrad,
     1                  (zn0(jg),zt0(jg),zs0(jg),jg=1,ngases)
c
  278   continue
c
c               average n0, t0, s0
c
        zv = 1.0 / rsurf(nsurf)**2
c
        do 284 jg = 1, ngases
        zzt0 = 0.0
        if (zn0tot(jg).gt.epslon) zzt0 = zt0tot(jg) / zn0tot(jg)
        zt0tot(jg) = zzt0
        zn0tot(jg) = zn0tot(jg) * zv
        zs0tot(jg) = zs0tot(jg) * zv
  284   continue
c
c
        write (nprint,10281) (ihblnk,jg=1,ngases)
        write (nprint,10282)
     1                  (zn0tot(jg),zt0tot(jg),zs0tot(jg),jg=1,ngases)
c
c    .    .    .    .    .    .    .    .    .    .    .    .    .    .
c
c
cl      3)      fluxes
c
c
c
  300   continue
c
        call reseti(isrc,4,0)
        call resetr(zofl,4,0.0)
        zzifl = 0.0
        zzofl = 0.0
        zzsfl = 0.0
c
        do 306 jsrc = 1, nsrces
        igas = ngasi(jsrc)
        isrc(igas) = isrc(igas) + 1
        if (isrc(igas).gt.4) go to 9030
        i001 = isrc(igas)
        zifl(i001,igas) = fluxin(jsrc)*fluxn(jsrc)
        ztifl(i001,igas) = e0in(jsrc)
        zzifl = zzifl + zifl(i001,igas)
c
        do 304 jg = 1, ngases
        zofl(jg) = zofl(jg) + fluxn(jsrc)*outflx(jg,jsrc)
        zzofl = zzofl + fluxn(jsrc)*outflx(jg,jsrc)
  304   continue
c
        zzsfl = zzsfl + sflux(jsrc)*fluxn(jsrc)
  306   continue
c
c
c
        write (nprint,10309)
        if(npuff.eq.0) go to 310
c
        write(nprint,10310)
     1  (gflmon(jh),lhgas(jh),jh=1,mhyd)
c
 310    continue
c
        do 314 jg = 1, ngases
        is = isrc(jg)
c
        if(is.eq.1 .or. (is.eq.2 .and. ztifl(2,jg).gt.epslon))
     1          write (nprint,10311) lhspec(jg),
     2          (zifl(js,jg),ztifl(js,jg),js=1,is)
        if(is.eq.2 .and. ztifl(2,jg).le.epslon)
     1          write (nprint,10312) lhspec(jg),
     2          zifl(1,jg),ztifl(1,jg),zifl(2,jg)
        if(is.eq.3)
     1          write (nprint,10311) lhspec(jg),
     2          (zifl(js,jg),ztifl(js,jg),js=1,2),zifl(3,jg)
c
  314   continue
c
c
c
        write (nprint,10314) (lhspec(jg),zofl(jg),jg=1,ngases)
        write (nprint,10315) zzsfl, zzofl, zzifl
c
c
cl      4)      neutral outflux energy spectrum
c
c       ******* warning:  this coding has not been fully tested. *******
c
  400   continue
c      the following used to say "if(nlpomt(20)) go to 449"
c      if neutral efflux is ever debugged, similar control will be needed
        go to 499
 490    lpage = lpage + 1
        write (nprint,10400) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
c
c***********************************************
c       charge exchange neutrals outflux spectrum -- chord line avg.
c       restrictions on current version: tokamak has circular cross-
c       section and is assumed cylindrical;
c       energy spectrum is logarithmic; max. energy .le. 1.e5 (ev)
c
c               based on m.h. hughes and d.e. post, a monte carlo
c               algorithm for calculating neutral gas transport in
c               cylindrical plasmas, journal of computational physics
c               28, pp. 43-55, july 1978.  the integral performed
c               here is found on p. 53 of this article, and the
c               ion temperature computation on pp. 53-54.
c
c               units:  energies and temperatures
c                       are in ev, angles in degrees, densities
c                       in (cm)**-3
c                       velocities in cm/sec, masses in ev
c
c
c               variables:
c                       imaxen = maximum energy-bin used in temperature
c                               computation
c                       iminen = minimum energy-bin used in temperature
c                               computation
c                       iphi   = number of bins in phi-integration
c                       itheta = number of bins in theta-integration
c                       izones = number of radial zones
c                       izonet(4) = saves izone0 values for later use
c                               in temperature computation printout
c                       izone0 = innermost zone intersected by chord,
c                               for given angle
c                       izone1 = izone0 + 1
c                       jg     = index of neutral species
c                       j,j0(9),j1,j2,j3,j5,j6,j7,j9=index variables
c                       j8     = index to mark program flow
c                       zaabb = a trigonometric factor
c                       zc    = speed of light in cm/sec
c                       zcorel(4,4) = correlation coefficient in
c                               least-squares fitting
c                       zcp    = cosine of zphi
c                       zct    = cosine of ztheta
c                       zdenrg = 'differential' energy increment
c                       zdgrad= converts from degrees to radians
c                       zdx   = 'differential' length element of chord
c                       zenrgy(j) = energy of incident neutral
c                       zentot(5,4) =
c                       zeta = log of the attenuation factor
c                       zfixit = numerical factor for integral
c                               normalization
c                       zgmau = converts grams to atomic mass units
c                       zintgl = (partial sum for) the line integral
c                       zlambd = reciprocal of the mean free path in
c                               current zone
c                       zlgrle = log of electron temperature (in ev)
c                       zlgrlp = log of relative energy of charge-
c                               exchange particles
c                       zn0all(j) = total neutral density, zone j,
c                               summed over all species
c                       znhall(j) = total ion density, zone j, summed
c                               over all species
c                       znnn   = number of terms in least-squares sum
c                       znorml(4,4) = normalization factor for line integral
c                       zohboy = a (rather messy) trigonometric factor
c                       zpboun
c                       zphi   = gxphi, converted to radians
c                       zrelnp = relative energy of charge-exchange
c                               particles wrt ions
c                       zr00  = shortest distance from center line
c                               of minor cross-section to the chord
c                       zspect(100,4,4) = final result: normalized
c                               energy spectrum ;
c                               also used in least-squares fitting
c                       zst    = sine of ztheta
c                       zst2  = sine-squared of ztheta
c                       zsumx,zsumxx,zsumxy,zsumy,zsumyy = accumulator
c                               variables for least-squares fitting
c                       zsvbre = electron ionization sigma-v-bar
c                       zsvbrp = proton ionization sigma-v-bar
c                       zsvbrx = charge-exch. sigma (not -v-bar)
c                       zsvcfe(9) = polynomial  coefficients for zsvbre
c                       zsvcfp(9) = polynomial  coefficients for zsvbrp
c                       ztboun(5)
c                       ztemp = temporary variable
c                       ztheta = gxthet, converted to radians
c                       ztimax(4,4) = calculated ion temperature value,
c                               using least-squares fitting routine
c                       ztotxx,ztotxy,ztotyy = subtotals for
c                               least-squares fitting
c                       zveloc = lab-frame velocity of outgoing neutral
c                       zyvalu(100,4) = y-values for least-squares fitting
c                       z4bypi = 4./ pi
c
c
c       *** initialize *********************************************
c
c       number of zones may not exceed 55, else dimensions
c        of znhall, zn0all will be exceeded
                izones = min0(nsurf-1,55)
c
c       number of energy bins may not exceed 100,
c        else dimensions of zspect, zenrgy will be exceeded
                ngxene = min0(ngxene,100)
c
c       numerical factors
                zc =  fcc * 10.**fxc
                zdgrad = fcpi / 180.
                zfixit = 4./sqrt(fcpi)
                zgmau  = fcau * 10.**fxau
                z4bypi = 4./fcpi
c
c       energy bins
                zdenrg = (gxmaxe/gxmine)**(1./float(ngxene-1))
                zenrgy(1) = gxmine
                do 401 j = 2,ngxene
                   zenrgy(j) = zenrgy(j-1) * zdenrg
  401           continue
c
c       total ion density
                znhall = 0.0
                do 402 j  = 1,izones
                do 402 jg = 1,ngases
                   znhall(j) = znhall(j) + denion(jg,j)
  402           continue
c
c       total neutral density
                zn0all = 0.0
                do 403 jsrc = 1,nsrces
                do 403 j    = 1,izones
                do 403 jg   = 1,ngases
                   zn0all(j) = zn0all(j) + fluxn(jsrc)*den0(jg,j,jsrc)
  403           continue
c
c       determine number of arguments in gxthet and gxphi
                itheta = 1
                iphi   = 1
                do 404 j = 2,4
                   if (gxthet(j).gt.epslon) itheta = itheta + 1
                   if (gxphi(j).gt.epslon)  iphi   = iphi   + 1
  404           continue
c
c       set up table of trig integrals for normalization
                ztboun(1) = 0.
                zpboun(1) = 0.
                ztboun(itheta+1) = fcpi/2.
                zpboun(iphi  +1) = fcpi/2.
                do 406 j2 = 2,itheta
                   ztboun(j2) = zdgrad * (gxthet(j2)-gxthet(j2-1))/2.
  406           continue
                do 407 j3 = 2,iphi
                   zpboun(j3) = zdgrad * (gxphi(j3)-gxphi(j3-1))/2.
  407           continue
                do 408 j3 = 1,iphi
                do 408 j2 = 1,itheta
                   znorml(j2,j3) = zfixit * (ztboun(j2+1)-ztboun(j2)
     1                  +.5*(sin(2.*ztboun(j2+1))-sin(2.*ztboun(j2))))
     2                  *(sin(zpboun(j3+1))-sin(zpboun(j3)))
  408           continue
c
c       reset working variables
                j2 = 1
                j3 = 1
                j9 = 1
                ztheta = 0.0
                zphi   = 0.0
                call resetr(zspect,500*ngases,0.0)
c
c
c
c       *** loops begin here ***************************************
c
c       *** loop over theta
  410   continue
c
c       convert theta to radians
                ztheta = gxthet(j2) * zdgrad
c
c       get useful trig coefficients
                zst   = sin(ztheta)
                zst2  = zst**2
                zct   = cos(ztheta)
c
c       *** loop over phi
  411   continue
c
c       convert phi to radians
                zphi  = gxphi(j3) *  zdgrad
c
c       get useful trig coefficients
                zcp   =  cos(zphi)
                zaabb =  zst2 + (zct*zcp)**2
                zohboy= zst2 * (rmins/zaabb)**2
c
c       get izone0,  innermost zone for given angle
                izone1 = 1
                zr00  = rsurf(nsurf) * zst / sqrt (zaabb)
  412           izone1 = izone1 + 1
                if(zr00.gt.rsurf(izone1)) go to 412
                izone0 = izone1 - 1
                if (j3.eq.1) izonet(j2) = izone0
c
c       *** loop over energies
  414   continue
c
c       *** loop over species
        do 450 jg = 1,ngases
c
c
c       *** computation loop
c
c       reinitialize
                zeta = 0.0
                zintgl =0.0
                zveloc = sqrt(2.* zenrgy(j9)*evs/rmass(jg))
c
  415   continue
c
c       set indices -- see computed goto after line 440
c              *first pass
                j5 = izones
                j6 = izone1
                j7 = -1
                j8 = 1
                if (izone1 .le. izones) go to 420
c              *second pass -- innermost zone
  418           j5 = izone0
                j6 = izone0
                j7 = 1
                j8 = 2
                go to 420
c              *third pass
  419           j5 = izone1
                j6 = izones
                j7 = 1
                j8 = 3
                if (izone1 .gt. izones) go to 445
c
  420   continue
c
        do 440 j = j5,j6,j7
c
c       geometry - find next zdx
                if (j8 .eq. 2) go to 421
c                      *case where j8 .ne. 2
                        zdx = sqrt (rsurf(j+1)**2/zaabb - zohboy)
     1                          - sqrt (rsurf(j)**2/zaabb - zohboy)
                        go to 425
c                      *case where j8 .eq. 2 (innermost zone)
  421                   zdx = 2.*sqrt (rsurf(izone1)**2/zaabb - zohboy)
                        if (zdx .le. epslon) go to 440
  425           continue
c
c       find sigma-v-bars in this zone
c       formulae from freeman and jones, and hughes and post
c       first get relative energy -- this is a fudge
c       strictly, should select from an energy distribution
                zrelnp = (zenrgy(j9) + z4bypi*tiin(j))
     1                  * zgmau / rmass(jg)
                zlgrle = log(tein(j))
                zlgrlp = log(zrelnp)
c
c       next get ch. x. sigma
                zsvbrx = .6937e-14 * (1. - 0.155*log10(zrelnp))**2
     1                  /(1. + .112e-14 * zrelnp**3.3)
     2                  * sqrt (2.*zrelnp*evs/zgmau )
c
c       next get electron and proton ionization sigma-v-bars
                zsvbre = 0.0
                zsvbrp = 0.0
                do 430 j1 = 1,9
                   zsvbre = zsvbre + zsvcfe(j1) * zlgrle**j0(j1)
                   zsvbrp = zsvbrp + zsvcfp(j1) * zlgrlp**j0(j1)
  430           continue
                zsvbre = exp(zsvbre)
                zsvbrp = exp(zsvbrp)
c
c       find the attenuation factor
                zlambd = (denein(j)*zsvbre + znhall(j)*zsvbrp
     1                  + znhall(j)*zsvbrx)     /zveloc
                zeta = zeta - zlambd*zdx
c
c       integrate
c       ztemp is for overflow/underflow check -- 70 is a somewhat
c       arbitrary number
c       .5*zlambd*zdx --> get attenuation factor from center of zone
                ztemp = zeta + .5*zlambd*zdx  - zenrgy(j9)/tiin(j)
                if (abs(ztemp) .gt. 70.) go to 439
                zintgl = zintgl +
     1                  (denion(jg,j) * zn0all(j) *zsvbrx*
     2                  sqrt(zenrgy(j9))* exp(ztemp)
     3                  / tiin(j)**1.5)                       * zdx
  439           continue
c
c
  440   continue
c
c       *go from outermost zone to innermost, and then
c       back out
                go to (418,419,445), j8
c
c       *** end of computation loop
  445   continue
c
c
c
c       store results
                if (j3.eq.1) zspect(j9,j2,jg) = zintgl
                zspect(j9,itheta+1,jg) = zspect(j9,itheta+1,jg) +
     1                                   znorml(j2,j3)*zintgl
c
c       *** end of loop over species
  450   continue
c
c       *** end of loop over energies
                j9 = j9 + 1
                if (j9 .le. ngxene) go to 414
c
c       *** end of loop over phi
                j9 = 1
                j3 = j3 + 1
                if (j3 .le. iphi) go to 411
c
c       *** end of loop over theta
                j3=1
                j2 = j2 + 1
                if (j2.le.itheta) go to 410
c
c       *** end of loops *******************************************
c
c
c       *** totals -- integrate over energy
c          energy intervals are delimited by geometric mean of
c          zenergy(j9), zenergy(j9+1)
                ztemp = sqrt(zdenrg) - sqrt(1./zdenrg)
                do 455 jg = 1,ngases
                do 455 j2 = 1,itheta+1
                do 455 j9 = 1,ngxene
                   zentot(j2,jg) = zentot(j2,jg) +
     1                  zspect(j9,j2,jg) * zenrgy(j9) * ztemp
  455           continue
c
c
c       *** temperature computations
c
            do 469 jg=1,ngases
            do 469 j2=1,itheta
c
c               compute least-squares fitting domain
c               domain contains at least four points
c
c                    exclude zeroes --> ceiling on imaxen
                        do 461 j = ngxene,1,-1
                           if (zspect(j,j2,jg).gt.epslon) go to 462
  461                   continue
  462                   imaxen = j
c
c                    minimum energy
                        izone0= izonet(j2)
                        ztemp = gxdome(1) * tiin(izone0)
                        do 463 j = imaxen-4,1,-1
                           if (zenrgy(j).le.ztemp) go to 464
  463                   continue
  464                   iminen = j+1
c
c                    maximum energy
                        ztemp = gxdome(2) * tiin(izone0)
                        do 465 j = iminen+4,imaxen
                           if (zenrgy(j).ge.ztemp) go to 466
  465                   continue
  466                   imaxen = j-1
c
c
c               compute y-values for least-squares fitting
                   call resetr(zyvalu,400,0.0)
                   do 467 j = iminen, imaxen
                      zyvalu(j,jg)=log(zspect(j,j2,jg)/sqrt(zenrgy(j)))
  467              continue
c
c               least-squares fitting -- fit straight line to zyvalu
c                as a function of zenrgy
                   zsumx = 0.
                   zsumy = 0.
                   zsumxx= 0.
                   zsumxy= 0.
                   zsumyy= 0.
                   do 468 j = iminen, imaxen
                      zsumx = zsumx + zenrgy(j)
                      zsumy = zsumy + zyvalu(j,jg)
                      zsumxx= zsumxx+ zenrgy(j)**2
                      zsumxy= zsumxy+ zenrgy(j)*zyvalu(j,jg)
                      zsumyy= zsumyy+ zyvalu(j,jg)**2
  468              continue
                   znnn = float(imaxen-iminen+1)
                   ztotxx = zsumxx - zsumx*zsumx /znnn
                   ztotxy = zsumxy - zsumx*zsumy /znnn
                   ztotyy = zsumyy - zsumy*zsumy /znnn
c
c               correlation coeficient -- shows goodness of fit
                   zcorel(jg,j2) = ztotxy / sqrt(ztotxx*ztotyy)
c
c
c
c               temperature = negative reciprocal of line slope
                   ztimax(jg,j2) = - ztotxx/ztotxy
c
  469       continue
c
c
c
c
c       *** output
c
c          print spectra
c
                write(nprint,10480)
                write(nprint,10481)(gxphi(j3),j3=1,iphi)
                if (ngases.eq.1)
     1             write(nprint,10482) (gxthet(j2),j2=1,itheta),
     2             (lhspec(1),j2=1,itheta+1)
                if (ngases.gt.1)
     1             write(nprint,10483) (gxthet(j2),j2=1,itheta),
     2             ((lhspec(jg),jg=1,ngases),j2=1,itheta+1)
                do 484 j=1,ngxene
                write(nprint,10484) zenrgy(j),
     1             ((zspect(j,j2,jg),jg=1,ngases),j2=1,itheta+1)
  484           continue
c
c          print totals
                write(nprint,10485)
                write(nprint,10486)
     1             ((zentot(j2,jg),jg=1,ngases),j2=1,itheta+1)
c          print temperature, correlation coefficient
                if (ngases.eq.1) write(nprint,10490)
                if (ngases.gt.1) write(nprint,10491)
                write(nprint,10492) (lhspec(jg),jg=1,ngases)
                do 495 j2 = 1,itheta
                   j = izonet(j2)
                   write(nprint,10495) gxthet(j2),rsurf(j),
     1              tiin(j),(ztimax(jg,j2),zcorel(jg,j2),jg=1,ngases)
  495           continue
c
c************************************************
c
c
  499   continue
c
        return
c
c
c.......................................................................
c
c
c
cl      10)     short printouts
c
c
c
 1000   continue
c
c
        zn0(1) = 0.0
        ze0 = 0.0
        zna = 0.0
        zea = 0.0
        zin = 0.0
        zout = 0.0
        zspt = 0.0
c
        icells = nsurf - 1
c
        do 1008 jsrc = 1, nsrces
c
        do 1004 jg = 1, ngases
        zn0(1) = zn0(1) + den0(jg,1,jsrc)*fluxn(jsrc)
        ze0 = ze0 + den0(jg,1,jsrc)*eneut(jg,1,jsrc)*fluxn(jsrc)
        zna = zna + den0(jg,1,jsrc)*fluxn(jsrc)
        zea = zea +
     1          den0(jg,icells,jsrc)*eneut(jg,icells,jsrc)*fluxn(jsrc)
        zin = zin + fluxr(jg,jsrc)*fluxn(jsrc)
        zout = zout + outflx(jg,jsrc)*fluxn(jsrc)
 1004   continue
c
        zin = zin + fluxin(jsrc)*fluxn(jsrc)
        zspt = zspt + sflux(jsrc)*fluxn(jsrc)
 1008   continue
c
c               charge exchange loss (ergs/sec)
c
c
        zchex = 0.0
        if (rplsma.le.epslon) go to 1040
c
        iz0 = 1
        iz1 = 1
        icells = nsurf - 1
c
        do 1038 jz = lcentr, ledge
c
        zni = 0.0
        do 1024 jh = 1, mhyd
        zni = zni + rhohs(jh,2,jz)
 1024   continue
c
        if (iz0.ge.icells) go to 1028
        zr = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
        if (zr.le.(rsurf(iz1) + rsurf(iz1+1))) go to 1028
        iz1 = iz1 + 1
        iz0 = iz1 - 1
        iz1 = min0(iz1,icells)
 1028   continue
c
        zint1 = (zr - rsurf(iz0) - rsurf(iz0+1)) /
     1                                  (rsurf(iz1+1) - rsurf(iz0))
        zint0 = 1.0 - zint1
c
c
 1030   continue
c
        zwx = 0.0
        zwxi = 0.0
c
        do 1034 jsrc = 1, nsrces
        zwx = zwx +
     1  (spi(iz0,jsrc)*zint0 + spi(iz1,jsrc)*zint1)*fluxn(jsrc)
        zwxi = zwxi +
     1  (spii(iz0,jsrc)*zint0 + spii(iz1,jsrc)*zint1)*fluxn(jsrc)
 1034   continue
c
        zchex = zchex - (zwx*evs + zwxi*tis(2,jz)) * zni*1.5*dx2i(jz)
c
c
 1038   continue
c
c
        zchex = zchex * (usep * 2.0 * avi(mzones,12,1) * uisl**3)
 1040   continue
c
        zt0(1) = 0.0
        if (zn0(1).gt.epslon) zt0(1) = ze0*evs / zn0(1) * useh
        if (zna.gt.epslon) zea = zea*evs / zna * useh * 1.5
c
        z0 = uesl**2 * uest
        zspt = zspt * z0
        zout = zout * z0
        zin  = zin  * z0
c
        zn0(1) = zn0(1) * used
c
        write (kunit,11000) zn0(1), zt0(1), zchex, zin, zout, zea, zspt
c
        return
c
c
cl      90)     errors
c
c
 9030   continue
        call error_olymp(1,iclass,isub,3,
     1          'more than 4 sources per gas species ')
        call error_olymp(2,igas,2,1,'gas num. ')
        return
c
c
c
cl      f)      format statements
c
c
c
10100   format(5x,'neutral gas:'//5x,'neutral hydrogen profiles are ',
     1  'computed using a multi-species, multi-influx, monte-carlo ',
     2  'method.'/5x,'profiles are computed every   ',i5,
     3  ' (ngprof) timesteps, with   ',i5,
     4  ' (ngpart) test particles per influx,'/5x,
     5  'on  a  radial  grid,   with   ',i5,
     6  ' (ngzone) cells, and with   ',i5,
     7  ' (ngsplt) splitting surfaces.')
10101   format(5x,2pf6.2,' 0f ','escaped plasma ions are returned',
     1  ' with a  ',1pe10.3,1x,a10,
     2  ' (tcold), mono-energetic distribution.'/10x,'(nlglim=f)')
10102   format(5x,2pf6.2,' 0f escaped plasma ions have a maxwellian',
     1  ' distribution at the edge ti, and are reflected off the',
     2  ' limiter.'/10x,'(nlglim=t)')
10103   format(5x,'escaping neutrals ',a6,
     1  ' reflected (nlgref), and ',a5,
     2  ' cause fe sputtering (nlgspt) which ',a5,' included.')
10120   format(/5x,'gas feed has a ',a14,' distribution at ',1pe10.3,
     1  1x,a10,' (tcoldp)')
10121   format(/15x,'time',5x,'neutral influx (',a10,')'/14x,'(msec)',
     1  4(5x,a10))
10122   format(8x,0pf12.3,4(5x,1pe10.3))
10123   format(1h1)
10130   format(1x/
     1  5x,'line average density monitoring parameters are:')
10131   format(1x/
     1  5x,'vol. ave. density monitoring parameters are:')
10132   format(1x/15x,'time on',5x,'time off',15x,
     1  7hpercent,17x,'percent',5x,
     2  'max. neutral influx (part/cm2*s)'/
     3  16x,'(sec)',8x,'(sec)',2(15x,a10)/)
10133   format(13x,f9.5,4x,f9.5,14x,f6.3,14x,f6.3,12x,e10.3)
10135   format(1x/
     1  17x,'time',6x,'density',16x,
     2  'percent',14x,'percent',9x,
     3  'max. neutral influx (part/cm2*s)'/16x,'(sec)',5x,
     4  '(part/cm3)',13x,a10,13x,a10/)
10136   format(13x,f9.5,3x,e10.3,14x,f6.3,14x,f6.3,12x,e10.3)
c10140   format(5x,'pellet(',a10,')injection at t=',x,f7.4,
c     1  x,'(sec):'/5x,'radius=',e9.3,' (cm)','velocity=',
c     2  x,e9.3,' (cm/sec)')
c
10200   format(1h1,2x,a48,10x,a72/
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.'/
     3  10x,'neutral profile at timestep ',i5,',  time  ',0pf12.3,
     4  '  millisecs., ',i5,' particles  tcold is ',1pe8.1,' ev')
10201   format(1x,'  zone   in-surface-out',8x,'te',10x,'ti',7x,
     1          'electrons',5(2x,a10))
10202   format(10x,'centimeters',10x,'ev',10x,'ev',6x,'part/cu cm',
     1          5(2x,a10))
10220   format(2x,i4,2x,0pf7.2,1x,0pf7.2,1x,a1,7(2x,1pe10.3),2x,1pe9.2)
10221   format(25x,4(2x,'----------'),3(a1,1x,'----------'),a1,1x,
     1          '---------')
10222   format('  **** average value ****',7(2x,1pe10.3),2x,1pe9.2)
c
10250   format(/'  zone  avg.rad',3(a1,4x,'neutral ',a10,
     1          ' -- source  '))
10251   format(15x,3(a1,'  part/cu cm',5x,'ev',4x,'part/cc-sec'))
10255   format(2x,i4,2x,0pf7.2,3(2x,3(1x,1pe10.3)))
10281   format(15x,3(a1,2x,'--------------------------------'))
10282   format('  ** average **',3(2x,3(1x,1pe10.3)))
c
10309   format(1x)
10310   format(1x,'dens. monitoring influx - ',
     1  2(e10.3,'  (',a10,')'))
10311   format(1x,a10,' influx--',2(1pe9.2,' at',1pe9.2,' ev,'),
     1          1pe9.2,' at local ion temperature (volume source)')
10312   format(1x,a10,' influx--',1pe9.2,' at',1pe9.2,' ev,',
     1          1pe9.2,' at local ion temperature (volume source)')
10314   format(' outflux of neutral ',a10,'=',1pe9.2,
     1          3(',    ',a10,'=',1pe9.2))
10315   format(' influx of sputtered   iron   =',1pe9.2,
     1          ', total outflux=',1pe9.2,', total influx =',1pe9.2)
10400   format(1h1,2x,a48,10x,a72//
     1  '  -',i2,'-  *** time step ',i5,' ***',14x,'time =',
     2  0pf12.3,'  millisecs.',12x,'dt =',0pf12.6,'  millisecs.')
10480   format(//1x,
     &  'line-integrated neutral spectrum (hughes-post)')
10481   format(1x,'phi =',4(f5.2,2x),'degrees'//)
10482   format(1x,'energy(ev)',2x,'theta=',3(1f5.2,5x),1f5.2,2x,
     &                  'integrated'/16x,5(a10))
10483   format(2x,'energy (ev)',5x,4('theta =',1f5.2,8x),2x,
     &                  'integrated'/16x,10(a10))
10484   format(1x,1pg12.3,2x,10(1pg10.3))
10485   format(1x,'------------------------------------------------',
     &  '------------------------------------------------------------')
10486   format(4x,'totals',5x,10(1pg10.3)//)
10490   format(//3x,'theta',3x,'radius',2x,'known ti',2x,
     &          'computed ti',2x,'correlation',4x)
10491   format(//3x,'theta',3x,'radius',2x,'known ti',2x,
     &          2('computed ti',2x,'correlation',4x))
10492   format(2x,'degrees',4x,'cm',7x,'ev',12x,
     &          a10,18x,a10/)
10495   format(2x,1f6.2,3x,1f6.2,1pg11.3,2(2(1pg12.3),4x))
11000   format(' n0-axis=',1pe8.2,' t0-axis=',1pe8.2,' cx.loss=',1pe10.2
     1 ,' influx=',1pe8.2,'  outflux=',1pe8.2,' e-edge=',1pe8.2
     2 ,'  fe.flux=',1pe8.2)
        end
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@ncprnt
c       dps 07-jul-89 15.11 add neutral impurities
c       dps 15-may-89 15.09 eliminate unneeded variables with move of
c                     IRE code into predictor-corrector loop
c       dps 24-aug-88 15.00 incorporate routines into DIO
c       rhw 08/08/84: change print format 11000.
c       rhw 05/06/84: change print of radiation.
c       rhw 17/05/84: ins. print xnini, recycling, beta, local radiation
c       rhw 20/03/84: insert print of ccons on short printout
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine ncprnt(kout,kcall)
c
c       2.20.7   print of non-corona radiation results
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'comncr.m'
      include 'comadp.m'
c
      dimension id(28)
c
      if ((mimp.le.0).or.(natomc.ne.3))   return
c
      if (kout.gt.0)   go to 100
      if (kcall.gt.1)   go to 20
c
c    initial print-out of non-corona data read in by sr. ncdata
c
      write (nprint,10100) (lhspec(js),js=limp1,limpn)
c
      write (nprint,10111)  recflx,recscr,beta
      if (ifneo.le.0)   write (nprint,10112)
      if (ifneo.gt.0)   write (nprint,10113)
c
c  Write out ADPAK input data
c
      write (nprint,10120) (naden(ji),ji=1,mimp)
      write (nprint,10121) (nadtip(ji),ji=1,mimp)
      write (nprint,10122) (neci(ji),ji=1,mimp)
      write (nprint,10123) (ncxx(jh),jh=1,mhyd)
      write (nprint,10124) (cermlt(ji),ji=1,mimp)
      write (nprint,10125) (cizmlt(ji),ji=1,mimp)
      write (nprint,10126) (cxrmlt(ji),ji=1,mimp)
c
      return
c
c-----------------------------------------------------------------------
c
   20 continue
c
c    print-out of non-corona impurity densities at various timesteps
c
      if (nlpomt(21))   go to 99
c
      do 90 ji=1,mimp
         nk = nkimp(ji)
         jimp = limp1 - 1 + ji
         lpage = lpage + 1
         nk2 = min0(nk,8)
         zx0tot = 0.0
         zt = tai * uist * 1000.0
         zdt = delts * 1000.0
c
      do 21 i=1,nk
         id(i) = i
   21 continue
c
      write (nprint,10200)  label1(1:48),label5(1:72),
     1                      lpage,nstep,zt,zdt
      write (nprint,10201)  lhspec(jimp),(id(i),i=1,nk2)
c
      do 25 jz=lcentr,mzones
         ijzp = jz + 1
         zrb1 = avi(jz,15,1) * uisl
         zrb2 = avi(ijzp,15,1) * uisl
         zr = 0.5 * (zrb1 + zrb2)
c
      if (jz.le.ledge)   zx0tot = zx0tot + xn0(jz,ji) * avi(jz,4,1)
     1                   *dxzoni(jz)*uisl**3
c
c..15.11 Include transport neutrals in column output, but total at
c        bottom is still the total number of source neutrals.
c
      write (nprint,10202)  jz, zr, dq(jz,ji),
     x                  xns(jz,ji), (xn(jz,ik,ji),ik=0,nk2)
   25 continue
c
      write (nprint,10203) xntot(ji), zx0tot
c
c    print-out final stages, if nk2.lt.nk
c
      if (nk2.ge.nk)   go to 30
c
   26 continue
         nk1 = nk2 + 1
         nk2 = nk1 + 10
         nk2 = min0(nk2,nk)
c
      write (nprint,10204)  lhspec(jimp),(id(i),i=nk1,nk2)
c
      do 28 jz=lcentr,mzones
        ijzp = jz + 1
        zrb1 = avi(jz,15,1) * uisl
        zrb2 = avi(ijzp,15,1) * uisl
        zr = 0.5 * (zrb1 + zrb2)
        write (nprint,10205)  jz, zr, (xn(jz,ik,ji),ik=nk1,nk2)
   28 continue
c
      if (nk2.lt.nk)   go to 26
c
c    printout of impurity radiation
c
   30 continue
c
c    1. cooling rate
c
         iwd = nwcool/1000
      if (iwd.le.0)   go to 32
c
      call ncwlpr (weir,nk,ji,-1,nstep,zt,zdt,
     1             '  cooling rate  ')
c
   32 continue
         iwcool = nwcool - 1000*iwd
         iwd = iwcool/100
      if (iwd.le.0)   go to 40
c
      call ncwipr (weir,nk,ji,-1,nstep,zt,zdt,
     1             '  cooling rate  ')
c
   40 continue
c
c    2. line radiation
c
         iwd = nwline/1000
      if (iwd.le.0)   go to 42
c
      call ncwlpr (wline,nk,ji,-1,nstep,zt,zdt,
     1             ' line radiation ')
c
   42 continue
         iwline = nwline - 1000*iwd
         iwd = iwline/100
      if (iwd.le.0)   go to 50
c
      call ncwipr (wline,nk,ji,-1,nstep,zt,zdt,
     1             ' line radiation ')
c
   50 continue
c
c    3. ionisation loss
c
         iwd = nwioni/1000
      if (iwd.le.0)   go to 52
c
      call ncwlpr (wioni,nk,ji,-1,nstep,zt,zdt,
     1             'ionization loss ')
c
   52 continue
         iwioni = nwioni - 1000*iwd
         iwd = iwioni/100
      if (iwd.le.0)   go to 60
c
      call ncwipr (wioni,nk,ji,-1,nstep,zt,zdt,
     1             'ionization loss ')
c
   60 continue
c
c    4. recombination loss
c
         iwd = nwreko/1000
      if (iwd.le.0)   go to 62
c
      call ncwlpr (wreko,nk,ji,0,nstep,zt,zdt,
     1             'recombinat. radn')
c
   62 continue
         iwreko = nwreko - 1000*iwd
         iwd = iwreko/100
      if (iwd.le.0)   go to 70
c
      call ncwipr (wreko,nk,ji,0,nstep,zt,zdt,
     1             'recombinat. radn')
c
   70 continue
c
c    5. charge exchange loss
c
         iwd = nwchex/1000
      if (iwd.le.0)   go to 72
c
      call ncwlpr (wchex,nk,ji,0,nstep,zt,zdt,
     1             'charge exch radn')
c
   72 continue
         iwchex = nwchex - 1000*iwd
         iwd = iwchex/100
      if (iwd.le.0)   go to 80
c
      call ncwipr (wchex,nk,ji,0,nstep,zt,zdt,
     1             'charge exch radn')
c
   80 continue
c
c    6. bremsstrahlung
c
         iwd = nwbrem/1000
      if (iwd.le.0)   go to 82
c
      call ncwlpr (wbrem,nk,ji,0,nstep,zt,zdt,
     1             ' bremsstrahlung ')
c
   82 continue
         iwbrem = nwbrem - 1000*iwd
         iwd = iwbrem/100
      if (iwd.le.0)   go to 90
c
      call ncwipr (wbrem,nk,ji,0,nstep,zt,zdt,
     1             ' bremsstrahlung ')
c
   90 continue
c
   99 return
c
c-----------------------------------------------------------------------
c
  100 continue
c
c    short printout at every "nsedit" timesteps
c
         zint = (xzoni(mzones)-xbouni(mzones)) / dxboui(mzones)
c
      do 102 ji=1,mimp
         jimp = limp1 - 1 + ji
         zna = xns(ledge,ji)*zint + xns(mzones,ji)*(1.0-zint)
c
      write (nprint,11000)  nstep,ncrept,ji,xns(lcentr,ji),zna,
     x                      ccons(jimp),pflx(ji),ploss(ji),psorc(ji)
  102 continue
c
      return
c
10100 format (///,5x,'Data used by the impurity charge state ',
     x        'transport model',/,5x,54('*'),///,64x,a10,3x,a10)
10111 format (/,' Recycling factor for diffusive and drift outflux =',
     x        f8.3,//,' Recycling factor for scrape-off losses =',f8.3,
     x        //,' Acceleration factor beta for the rate coeffs =',
     x        f10.3)
10112 format (/,' Without neoclassical terms in diffusion and drift',/)
10113 format (/,' With neoclassical terms in diffusion and drift',/)
c
10120 format (/,' Atomic energy level calculation (0=Mayer, 1=More):',
     x        15x,i1,12x,i1)
10121 format (/,' Ionization potentials are (=1), are not (=0) ',
     x        'tabulated:',10x,i1,12x,i1)
10122 format (/,' Electron coll''al ionization rate (1=XSNQ, ',
     x        '2=Belfast, 3=Ygr):',5x,i1,12x,i1)
10123 format (/,' Cross sections for CX, Hyd. indices (1=OSAS, ',
     x        '2=GJ, 3=OSCT):',6x,i1,12x,i1)
10124 format (/,' Multiplier for electron collisional recombination ',
     x        'rates:',4x,1pe11.3,2x,1pe11.3)
10125 format (/,' Multiplier for electron collisional ionization ',
     x        'rate:',8x,1pe11.3,2x,1pe11.3)
10126 format (/,' Multiplier for charge exchange recombination ',
     x        'rate:',10x,1pe11.3,2x,1pe11.3)
c
10200   format (1h1,2x,a48,10x,a72/2x,1h-,i2,1h-,2x,
     1          '*** time step',i5,' ***',14x,'time =',f12.3,
     2          2x,'millisecs.',12x,'dt =',f12.6,2x,'millisecs.')
10201 format (/,'  zone   radius    source',40x,
     x        a10,'  densities (cm**-3)',/,11x,'cm    cm-3 sec-1',
     x        2x,'n-imp-total',3x,'n-0',3x,8(4x,'n-',i1,3x),//)
10202 format (1x,i4,f10.2,2x,1pe10.3,1pe11.3,1x,1p9e10.3)
10203 format (/,' particles:',16x,1p2e11.3)
10204 format ('1',//,' zone   radius',43x,a10,'  densities (cm**-3)',
     x        /,10x,'cm',5x,11(4x,'n-',i2,2x),//)
10205 format (1x,i4,f10.2,2x,1p11e10.3)
c
11000 format (' step=',i5,' substeps=',i2,' imp.#',i1,' n0=',
     x        1pe9.2,' na=',1pe9.2,' cons=',1pe9.2,' outflx=',1pe9.2,
     x        ' scroff=',1pe9.2,' vol.srce=',1pe9.2)
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@ncwlpr
c       dps 15-may-89 15.09 remove time-centering of equilibrium quantities
c       dps 24-aug-88 15.00 incorporate routines into DIO
c       dps 24-aug-88 adapt to 1-1/2-D BALDUR
c       rhw 04/06/84: first write-up
c***********************************************************************
c
        subroutine ncwlpr (we,nk,ji,in,iter,zt,zdt,text)
c
c       2.20.71  print of local non-corona radiations
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
      dimension  we(mj,28,2), id(28), zrad(28)
      character*1 text(16)
c
         lpage = lpage + 1
         nk2 = min0(nk,10)
         jimp = limp1 - 1 + ji
c
      do 10 i=1,nk
         id(i) = i + in
   10 continue
c
      write (nprint,100)  label1(1:48),label5(1:72),
     1                    lpage,iter,zt,zdt
      write (nprint,110)  lhspec(jimp),text,(id(i),i=1,nk2)
c
      do 15 jz=lcentr,mzones
        ijzp = jz + 1
        zrb1 = avi(jz,15,1) * uisl
        zrb2 = avi(ijzp,15,1) * uisl
        zr = 0.5 * (zrb1 + zrb2)
        zrads = 0.
c
      do 12 ik=1,nk
         zrad(ik) = we(jz,ik,ji) * usep
         zrads = zrads + zrad(ik)
   12 continue
c
      write (nprint,112)  jz, zr, zrads, (zrad(i),i=1,nk2)
   15 continue
c
c    print-out final stages, if nk2.lt.nk
c
      if (nk2.ge.nk)   go to 30
c
   20 continue
         nk1 = nk2 + 1
         nk2 = nk1 + 10
         nk2 = min0(nk2,nk)
c
      write (nprint,111)  lhspec(jimp),text,(id(i),i=nk1,nk2)
c
      do 25 jz=lcentr,mzones
        ijzp = jz + 1
        zrb1 = avi(jz,15,1) * uisl
        zrb2 = avi(ijzp,15,1) * uisl
        zr = 0.5 * (zrb1 + zrb2)
c
      do 22 ik=nk1,nk2
         zrad(ik) = we(jz,ik,ji) * usep
   22 continue
c
      write (nprint,113)  jz, zr, (zrad(ik),ik=nk1,nk2)
   25 continue
c
      if (nk2.lt.nk)   go to 20
c
   30 continue
c
      return
c
  100 format (1h1,2x,a48,10x,a72/2x,1h-,i2,1h-,2x,
     x        '*** time step',i5,' ***',14x,'time =',f12.3,
     x        2x,'millisecs.',12x,'dt =',f12.6,2x,'millisecs.')
  110 format (/,' zone     radius',20x,a10,2x,16a1,
     x        ' local (watt*cm**-3)',/,12x,'cm',4x,'sum over nk',
     x        10(4x,'w-',i2,2x),//)
  111 format ('1',//,' zone    radius',20x,a10,2x,16a1,
     x        ' local (watt*cm**-3)',/,12x,'cm',4x,
     x        11(4x,'w-',i2,2x),//)
  112 format (1x,i4,f10.2,2x,1pe11.3,1x,1p10e10.3)
  113 format (1x,i4,f10.2,2x,1p11e10.3)
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@ncwipr
c       dps 15-may-89 15.09 remove time-centering of equilibrium quantities
c       dps 24-aug-88 15.00 incorporate routines into DIO
c       dps 24-aug-88 adapt to 1-1/2-D BALDUR
c       rhw 04/06/84: first write-up
c***********************************************************************
c
        subroutine ncwipr (we,nk,ji,in,iter,zt,zdt,text)
c
c       2.20.72  print of volume-integrated non-corona radiations
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
      dimension  we(mj,28,2), id(28), zrad(28)
      character*1 text(16)
c
         lpage = lpage + 1
         nk2 = min0(nk,10)
         jimp = limp1 - 1 + ji
c
      do 10 i=1,nk
         zrad(i) = 0.0
         id(i) = i + in
   10 continue
c
      write (nprint,100)  label1(1:48),label5(1:72),
     1                    lpage,iter,zt,zdt
      write (nprint,110)  lhspec(jimp),text,(id(i),i=1,nk2)
c
      do 15 jz=lcentr,ledge
         ijzp = jz + 1
         zr = avi(ijzp,15,1) * uisl
         zrads = 0.0
c
      do 12 ik=1,nk
         zrad(ik) = zrad(ik) + we(jz,ik,ji)*usep*avi(jz,4,1)
     1              *dxzoni(jz)*uisl**3
         zrads = zrads + zrad(ik)
   12 continue
c
      write (nprint,112)  ijzp, zr, zrads, (zrad(i),i=1,nk2)
   15 continue
c
c    print-out final stages, if nk2.lt.nk
c
      if (nk2.ge.nk)   go to 30
c
   20 continue
         nk1 = nk2 + 1
         nk2 = nk1 + 10
         nk2 = min0(nk2,nk)
c
      do 21 ik=nk1,nk2
         zrad(ik) = 0.0
   21 continue
c
      write (nprint,111)  lhspec(jimp),text,(id(i),i=nk1,nk2)
c
      do 25 jz=lcentr,ledge
         ijzp = jz + 1
         zr = avi(ijzp,15,1) * uisl
c
      do 22 ik=1,nk
         zrad(ik) = zrad(ik) + we(jz,ik,ji)*usep*avi(jz,4,1)
     1              *dxzoni(jz)*uisl**3
   22 continue
c
      ijz = jz - 1
      write (nprint,113)  ijz, zr, (zrad(ik),ik=nk1,nk2)
   25 continue
c
      if (nk2.lt.nk)   go to 20
c
   30 continue
c
      return
c
  100 format (1h1,2x,a48,10x,a72/2x,1h-,i2,1h-,2x,
     x        '*** time step',i5,' ***',14x,'time =',f12.3,
     x        2x,'millisecs.',12x,'dt =',f12.6,2x,'millisecs.')
  110 format (/,' zone     radius',20x,a10,2x,16a1,
     x        ' volume integrated (watts)',/,12x,'cm',4x,'sum over nk',
     x        10(4x,'w-',i2,2x),//)
  111 format ('1',//,' zone    radius',20x,a10,2x,16a1,
     x        ' volume integrated (watts)',/,12x,'cm',4x,
     x        11(4x,'w-',i2,2x),//)
  112 format (1x,i4,f10.2,2x,1pe11.3,1x,1p10e10.3)
  113 format (1x,i4,f10.2,2x,1p11e10.3)
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@mprint   /baldur/code/bald  file diof.
c  rgb 07-aug-01 printed vnwars(jz)+vnneo1(1,jz), vewars(jz)+vnneo1(2,jz)
c  rgb 07-feb-00 removed equivalence (nlzzzz(1),zsig(1)), ...
c  rgb 14-aug-96 replaced fuzz with rndeps
c  rgb 30-may-96 set ix0 just before using
c  rgb 31-may-94 allw ihreg to have 24 elements, fill blanks with 'zz'
c  rgb 09-may-94 call prtheory instead of theory(3)
c  rgb 21-feb-94 more complete long printout when nstep = 0
c  rgb 14-dec-93 improved formatting of k-e totl, ...
c  rgb 13-dec-93 removed zone-centering parameter (bound(i,j), i=1,6)
c      and added ion source rates vs half-width
c  rgb 18-sep-92 corrected format 12304
c  rgb 21-jul-92 changed from zjb7b to zjtor for toroidal current dens
c  rgb 15-jul-92 added total currents to magnetics page
c  rgb 11-feb-92 limit range of index for ihreg in 1st page of long printout
c  les 02-jan-91 turned off vectorization of do 310 for cft77
c  les   nov-90  added d-3he fusion fast particles to zalfap,ztbalf
c  rgb 23-jun-90 18.41 temporary printout of variables in FREYA
c  rgb 21-jun-90 18.39 temporary printout of rhobis(jb,jz)
c      only when lbeams(32) .gt. 0
c  rgb 08-may-90 18.35 new conduction and convection linout(32) .ge. 0
c  rgb 20-mar-90 moved header for output to sbrtn THEORY
c  rgb 18-dec-89 replaced logical variable lthery with variable ltheor
c  mhr 05-oct-89 17.00 printout 3 terms in chi-e and in chi-i from empirc
c  rgb 05-jun-89 commented out unused format statements
c       dps 10-nov-88 15.07 alter average density calculations to be
c                     consistent with others and with 1-1/2-D and
c                     scrape-off; fix bug in printout of magnetic diffusion.
c       dps 18-oct-88 15.06 change by Bateman: remove <1/R**2> and put in
c                     ajtpbi(j) on magnetic diffusion page
c       dps 17-oct-88 15.06 make following change from Singer:
c       elg 11-oct-88 add subroutine theory's output
c       dps 08-feb-87 finish alterations to beta-poloidal calculation with
c                     changes to alpha component; include old calculation
c                     of internal inductance for leqtyp=0.
c       dps 16-dec-87 add to dd reaction printout on summary page
c       dps 14-oct-87 use new internal inductance calculation; improve
c                     poloidal beta computation.
c       dps 09-oct-87 fix bug in grad T computation of conducted power
c  rgb 12.66 02-sep-87 print out eta-neoclassical / eta-Spitzer
c               ( = 1/ftrap) rather than ftrap
c  rgb 12.65 29-aug-87 implemented the following changes from
c               ux1:[baldur.balsav]mprint.for  from Martha Redi
c       mhr 8-oct-86 included zeff in nu-el printout
c  rgb 08-jun-89  For previous changes, see file CHNGSOLD.
c******************************************************************************
c
c     nlpomt(j) (j=1,25) controls the following pages of long output
c     nlpomt(j) = .F prints output,  nlpomt(j) = .T omits output
c
c nlpomt(1)   T_e, T_i, n_e, n_i, v_loop, J_tor, q, beta  vs half-width
c nlpomt(2)   nu_el, nu_h, hydrogen and impurity density, Z_eff
c               as a function of half-width
c nlpomt(3)   energy losses from plasma within radius r(j) in watts
c nlpomt(4)   confinement times, profile averages
c nlpomt(5)   total diffusion and pinch coefficients vs half-width
c               also, diffusivities from sbrtn empirc
c               and n_e, T_e, and T_i profiles vs major radius
c nlpomt(6)   ion source rates vs half-width
c nlpomt(7)   ripple diffusivities vs half-width from sbrtn trcoef
c nlpomt(8)   neutral beam distribution in energy and pitch angle
c               in each zone and fits to the distribution
c nlpomt(9)   sawtooth profiles
c nlpomt(10)  magnetic diffusion equation variables vs half-width
c               (resistivity, driven currents, metric elements)
c nlpomt(11)  first half of neutrals printout
c               (temperatures, densities, charge-exchange losses)
c nlpomt(12)  second half of neutrals printout
c               (neutral particle source rates)
c nlpomt(13)  neutral beam energy levels and densities at each
c               energy level
c nlpomt(14)  neutral beam average energy per particle,
c               particle source rate, etc.
c nlpomt(15)  summary at bottom of main neutral beam printout
c nlpomt(16)  neutral beam driven current
c nlpomt(17)  neutral beam H(r) for 1/1, 1/2, 1/3 energy components
c nlpomt(18)  neutral beam total H(r), charge exchange rates, etc 
c nlpomt(19)  alpha particle or D-D fusion from sbrtn aprint
c nlpomt(20)  complete namelist output from sbrtn data
c               and neutral outflux spectrum
c nlpomt(21)  non-equilibrium coronal impurity radiation
c               from sbrtn ncprnt
c nlpomt(22)  diffusivities from sbrtn theory
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
        subroutine mprint(k)
c
c       3.2     generate large edits on unit nprint
c
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'cfokkr.m'
       include 'cfreya.m'
      include 'commhd.m'
       include 'cbparm.m'
       include 'cd3he.m'
c
c
        logical
     l  ilend,  ila,    ilb,
     l  inital
c
c
c
        common/comsaw/ tbsaws,sawte,sawne,sawdd,sawate,sawane,sawadd,
     1  sawtau,sawr1,sawri,sawrmx,sawqc,tmguir,njzq1,
     2  sawfus,sawalf,sawuth,estpol,uthavg,tauein,ohptot,
     3  tepre(mj),tepost(mj),teavg(mj),tipre(mj),tipost(mj),tiavg(mj),
     4  ajzpre(mj),ajzpst(mj),ajzavg(mj),ohravg(mj)
c
        dimension
     x   zd(idxchi)    , zloss(idxchi) , zlne(6)       ,
     r   zzni(idxchi)  , zrconf(5)     , ztconf(idchp2,5)              ,
     r   zscflx(idxion), zdvflx(idxion), zwaflx(idxion), zsig(mj)      ,
     r   zmfp(mj)      , ztcx(mj)      , ztis(mj)      , zrats(mj)     ,
     r   zpdiv(mj)     , zndiv(mj)     , zlden(mj)     , zei(mj)       ,
     r   yscrof(mj)    , zscrof(mj)    , zarray(9)     , zwimp(6)
c
c
cgb        equivalence
cgb     x   (nlzzzz(1),zsig(1))       , (nlzzzz(56),zmfp(1))      ,
cgb     x   (nlzzzz(111),ztcx(1))     , (nlzzzz(166),ztis(1))     ,
cgb     x   (nlzzzz(221),zrats(1))    , (nlzzzz(276),zpdiv(1))    ,
cgb     x   (nlzzzz(331),zndiv(1))    , (nlzzzz(386),zlden(1))    ,
cgb     x   (nlzzzz(441),zei(1))      , (nlzzzz(496),yscrof(1))   ,
cgb     x   (nlzzzz(551),zscrof(1))   , (nlzzzz(606),ztconf(1,1)) ,
cgb     x   (cfutz(380),zarray(1))
c
        data    inital /.true./
c
c       ion-ion neoclassical indicators (see subroutine trcoef).
c
        data    ikineo,incflx,ixdoff,iionxi,ishift /12,110,139,281,282/
c
c
        character *1 ihsign(3), ihs(idxchi)
        character *2 ihte, ihti, ihni, ihreg(24)
        character *10 ihbpol, ihpart, ihzeff
        character ihh(2)*6, ihi(3)*7, ihdont(2)*5
c
        data    ihsign /' ','o','-'/, ihte/'te'/, ihti/'ti'/, ihni/'ni'/
        data    ihbpol /'b-poloidal'/,
     h          ihpart /'particles '/, ihzeff /'   zeff   '/
        data    ihreg
     h          /'pc','bm','nc','te','cd','rp','bd','q ','sk','lm',
     h           'em','bo','mp','zz','dr','te','ti','tp','zz','ex',
     &           'th','cm','zz','zz'/
        data    ihh /' even ','except'/,
     h          ihi /'density','influx ','source '/
     h          ihdont /'don''t',' do  '/
        data    ipts /4/
c
c------------------------------------------------------------------------------
c
        data    iclass/3/,      isub/2/
c
        if (.not.nlomt3(isub))  go to 10
        call mesage(' *** 3.2 subroutine mprint bypassed ')
        return
   10   continue
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       nout (combas)
c
c------------------------------------------------------------------------------
c
        if (k.ne.1)     go to 200
c
c
c
c
cl      1.1)
c
c
  100   continue
        lpage = 1
c
c
cl      1.2)    title page
c
c
  120   continue
        nout = nprint
        call page
        call daytim ( nout )
        call blines(2)
c
        write (nprint,10100)
        write (nprint,10101)
c
c               major / minor radius
c
        zrmn = rmini * uiel
        zrmj = rmaji * uiel
c
        if (radius(2).le.epslon) go to 122
        write (nprint,10102) zrmj, lhlen, zrmn, lhlen
        go to 125
c
  122   continue
        if (nrfit.eq.1) write (nprint,10103) zrmj, lhlen, zrmn, lhlen
        if (nrfit.eq.2) write (nprint,10104) zrmj, lhlen, zrmn, lhlen
        if (nrfit.gt.2) write (nprint,10105) zrmj, lhlen, zrmn, lhlen
  125   continue
c
c               b, current
c
      zcurr = avi(mzones,2,1) * avi(mzones,3,1) * r0ref * bpoli(mzones)
     &  * avi(mzones,7,1) * uiei / (twopi * emu0)
        zbz = bzi * uieb
c
        if (eebfit.eq.0) write (nprint,10110) zbz, lhb, zcurr, lhcurr
        if(eebfit.gt.0) write(nprint,10111)
     1           zbz,lhb,zcurr,lhcurr,eebfit,ebfit
c
c               te
c
        zte0 = tes(2,lcentr) * useh
        zte1 = tes(2,mzones) * useh
c
        if (te(1).le.epslon) write(nprint,10115)
     1  ihte, zte0, lhtemp, ihte, zte1, lhtemp, ihte, ihte,
     2          eeteft,etefit
        if (te(1).gt.epslon) write(nprint,10116)
     1  ihte, zte0, lhtemp, ihte, zte1, lhtemp, ihte
c
c               ti
c
        zti0 = tis(2,lcentr) * useh
        zti1 = tis(2,mzones) * useh
c
        if (ti(1).le.epslon) write(nprint,10115)
     1  ihti, zti0, lhtemp, ihti, zti1, lhtemp, ihti, ihti,
     2          eetift,etifit
        if (ti(1).gt.epslon) write(nprint,10116)
     1  ihti, zti0, lhtemp, ihti, zti1, lhtemp, ihti
c
c               ne
c
        zne0 = rhoels(2,lcentr) * used
        zne1 = rhoels(2,mzones) * used
c
        write(nprint,10117) zne0, lhdens, zne1, lhdens
c
c               ni
c
        ila = .false.
c
        do 127 jh = 1, mxhyd
        if (ngas(jh).ne.0.and.fracth(jh).gt.epslon.and.
     1          .not.(denga0(jh).gt.epslon.or.dengas(jh,1).gt.epslon))
     2          ila = .true.
  127   continue
c
        do 128 ji = 1, mximp
        if (nimp(ji).ne.0.and.fracti(ji).gt.epslon.and.
     1          .not.(denim0(ji).gt.epslon.or.denimp(ji,1).gt.epslon))
     2          ila = .true.
  128   continue
c
        zni0 = rhoins(2,lcentr) * used
        zni1 = rhoins(2,mzones) * used
c
        if (dens(1).le.epslon.and.ila) write(nprint,10115)
     1  ihni, zni0, lhdens, ihni, zni1, lhdens, ihni, ihni,
     2  eefit,efit
        if (dens(1).gt.epslon.and.ila) write(nprint,10116)
     1  ihni, zni0, lhdens, ihni, zni1, lhdens, ihni
        if (.not.ila) write(nprint,10118)
     1  ihni, zni0, lhdens, ihni, zni1, lhdens
c
c
c               individual species
c
c
        write (nprint,10120)
c
        z0 = 0.0
        do 132 js = lhyd1, limpn
        z0 = z0 + totali(js)
  132   continue
c
        do 134 js = lhyd1, limpn
        zn0 = chi(js,lcentr) * uied
        zn1 = chi(js,mzones) * uied
        zfract = totali(js) / z0
        write (nprint,10122) lhspec(js), zn0, lhdens, zn1,
     1          lhdens, zfract, aspec(js), nzspec(js)
  134   continue
c
c
c               t's, dt's, and print flags
c
c
        z0 = uist * 1000.0
        ztin = tai * z0
        zdtin = dti * z0
        zdtmn = dtmini * z0
        zdtmx = dtmaxi * z0
c
        z0 = uest * 1000.0
        zspr = sedit * z0
        zspl = splot * z0
        ztpr = tedit(1) * z0
        ztpl = tplot(1) * z0
c
        write (nprint,10125) ztin, nedit, nsedit, nplot, zdtin, zspr,
     1          zspl, zdtmn, nskip, zdtmx, ztpr, ztpl
c
        ila = .false.
        ilb = .false.
c
        do 138 jt = 2, mxt
        if (tedit(jt).le.0.0) ila = .true.
        if (tplot(jt).le.0.0) ilb = .true.
        if (ila.and.ilb) go to 139
        ztpr = tedit(jt) * z0
        ztpl = tplot(jt) * z0
c
        if (.not.(ila.or.ilb)) write(nprint,10126) ztpr, ztpl
        if (ilb) write(nprint,10127) ztpr
        if (ila) write(nprint,10128) ztpl
c
  138   continue
c
  139   continue
c
c
c
c               transport coefficient fudge factors
c
        call page
c
        if(cfutz(incflx).gt.epslon) go to 146
        write(nprint,10135) ntrans
        go to 148
  146   continue
        write(nprint,11135) ntrans
  148   continue
c
        if (ntrans.eq.2) write (nprint,10136)
     1          cfutz(3),cfutz(7),cfutz(11),cfutz(13),cfutz(14),
     2          cfutz(4),cfutz(8),cfutz(12),cfutz(140),cfutz(16),
     3          cfutz(1),cfutz(5),cfutz(9),cfutz(18),cfutz(17),
     4          cfutz(2),cfutz(6),cfutz(10),
     5          cfutz(19)
c
        if(cfutz(ixdoff).lt.-0.1) write(nprint,10137)
c
        zshift=0.0
        if(cfutz(ikineo).le.epslon) go to 154
        if(cfutz(iionxi).le.epslon) go to 156
        if(cfutz(iionxi).gt.1.1   ) go to 152
        z0=abs(cfutz(ishift))
        if(z0.le.epslon) go to 150
        z1=min(z0,0.95)
        zshift=sign(z1,cfutz(ishift))
        write(nprint,12138) zshift
        go to 158
  150   continue
        write(nprint,11138)
        go to 158
  152   continue
        write(nprint,13138)
        go to 158
  154   continue
        write(nprint,14138)
        go to 158
  156   continue
        write(nprint,10138)
  158   continue
c
c
c               radiation/impurity radiation
c
c
        if (natomc.eq.1) write (nprint,10140)
        if (natomc.eq.2) write (nprint,10141) lholab
c
c               boundary conditions
c
        ii = nbound / 2 + 1
        ih = nbound - ii*2 + 3
        write (nprint,10142) ihi(ii), ihh(ih)
c
c               extrapolation
c
        if (.not.nlextr) write (nprint,10143) delmax
        if (nlextr) write (nprint,10144) errmax, delmax
c
c               equation flags
c
        iab = 1
        ic = 1
        id = 1
        if (nldiff) iab = 2
        if (nlsorc) ic = 2
        if (nlsord) id = 2
c
        write (nprint,10145) ihdont(iab),ihdont(ic),ihdont(id)
c
c               theta
c
        write(nprint,10146) theta
c
        return
c
c
c
cl      2.)     print first page of large time-step edit
c
c       in internal units:
c
c       zeff(r) = sum of (dens. of spec. j)*(mean z**2 of spec.j)
c               / sum of (dens. of spec. j)*(mean z of spec. j)
c
c..dx2i(j) = (zone bndry volume(j+1) - zone bndry volume(j))
c            / ( 2. * zone bndry volume(mzones) )
c  dx2inv(j) = 1. / dx2i(j)
c
  200   continue
        zt = tai*uist * 1000.0
        zdt = dtoldi * uist * 1000.0
        if (nlpomt(1)) go to 235
        lpage = lpage + 1
        write (nprint,10200) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write (nprint,10201)
        write (nprint,10202)
c
c
c       loop begins here
c
c
        j = 1
        ilend = .false.
  205   continue
        zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uiel ! zone center radius
c
c
        zte = tes(2,j) * useh
        zti = tis(2,j) * useh
        zne = rhoels(2,j) * used
        zni = rhoins(2,j) * used
c
        zj = 0.0
        if (j.le.ledge) zj = ajzs(2,j)*usej
        ij = min0(j,ledge)
cbate   zez = (ajzs(2,ij) - cjbeam*ajbs(ij))*eta(2,ij)*usev
c
        zbeta = (rhoels(2,j)*tes(2,j) + rhoins(2,j)*tis(2,j)) *
     1          8.0*fcpi / bzs**2
        ike = min (24, max (1, lkeflg(j) ) )
        iki = min (24, max (1, lkiflg(j) ) )
c
cbate        ij = j - 1
c
        write (nprint,10204) j, zrad, zte, ihreg(ike), zti,
     1          ihreg(iki), zne, zni, vloopi(j,2), zj, q(j), zbeta
c
c
        if (ilend) go to 235
        j = j + nskpi
        if (j.lt.mzones) go to 205
        ilend = .true.
        j = mzones
        go to 205
  235   continue
c
c       loop has ended
c
c
c
c       compute and print totals
c
c bateman  differential volume computed using avi(j,4,1) from sbrtn mhd
c
c
        idim = mhyd + mimp
        zzei = 0.0
        zzee = 0.0
        zzeb = 0.0
        zzne = 0.0
        zznion = 0.0
        call resetr(zzni,mxchi,0.0)
      zvols = 0.
      z0 = uisl**3  ! to convert volume from internal to standard units
c
        do 249 jz = lcentr, ledge
      zdvols = avi(jz,4,1)*dxzoni(jz)*z0  ! dV = V'(xi) * d xi at zone cent
      zvols  = zvols + zdvols
        zzee = zzee + chi(lelec,jz) * zdvols
        zzei = zzei + chi(lion,jz) * zdvols
        zzeb = zzeb + hebems(jz)*rhobis(2,jz) * zdvols
        zzne = zzne + rhoels(2,jz) * zdvols
        zznion = zznion + rhoins(2,jz) * zdvols
c
        do 244 js = lhyd1, limpn
        zzni(js) = chi(js,jz) * uisd * zdvols  + zzni(js)
  244   continue
c
  249   continue
c
        zzbeta = (zzee + zzei) * 5.333333*fcpi * uise*uisd
     &  / ( zvols * bzs**2 )
c
        zzee = zzee * uiee * uisd  ! total electron energy
        zzei = zzei * uiee * uisd  ! total ion energy
        zzeb = (zzeb * usee)         ! total beam energy
        zzne = zzne         ! total number of electrons
        zznion = zznion         ! total number of ions
cbate call scaler (zzni,idim,z0)     ! total number of each impurity ion
c
c bateman total plasma current for 1 1/2 D BALDUR
c
      jb = ledge + 1
      zi = uiei * avi(jb,2,1) * avi(jb,3,1) * r0ref * bpoli(jb)
     &  * avi(jb,7,1) / (twopi * emu0)
c
cend bateman
c
        if (nlpomt(1)) go to 295
        write (nprint,10205)    lhener, lhener, lhcurr
        write (nprint,10206)    zzee, zzei, zzne, zznion, zi, zzbeta
  295   continue
c
c
cl      3)      print second page of mprint edit
c
        if (nlpomt(2)) go to 395
c
        lpage = lpage + 1
c
        write (nprint,11300) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        if(mimp.le.0) go to 303
        if(mimp.lt.3) go to 302
        if(mimp.lt.5) go to 301
c       case of more than 4 impurities not yet implemented
        go to 304
  301   continue
        write(nprint,12301) (lhspec(js),js=lhyd1,limpn),
     1  ihzeff
        write(nprint,12302) lhlen,(lhdens,js=lhyd1,limpn)
        go to 304
  302   continue
        write(nprint,11301) (lhspec(js),js=lhyd1,limpn),
     1  ihzeff
        write(nprint,10302) lhlen,(lhdens,js=lhyd1,limpn)
        go to 304
  303   continue
        write(nprint,11301) (lhspec(js),js=lhyd1,lhydn)
        write(nprint,10302) lhlen,(lhdens,js=lhyd1,lhydn)
  304   continue
c
c
c       loop begins here
c
c
        j = 2
c       "j=2" in order to hold displays of ion densities vs. r to
c       one page.
c
        ilend = .false.
  305   continue
c
        zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uiel ! zone center radius
        i1 = 0
c
        if (mhyd.le.0) go to 312
        do 310 j1 = 1, mhyd
        i1 = i1 + 1
        zd(i1) = rhohs(j1,2,j) * used
  310   continue
  312   continue
c
        if (mimp.le.0) go to 322
        do 320 j1 = 1, mimp
        i1 = i1 + 1
        zd(i1) = rhois(j1,2,j) * used
  320   continue
  322   continue
c
c
c               compute sign of fluxes
c
c
c
        do 328 jp = lhyd1, limpn
        z0 = 0.0
        z1 = 0.0
        ihs(jp) = ihsign(1)
        if (j.le.lcentr) go to 328
c
        do 326 jp2 = 1, mchi
        z0 = z0 - aaaa(jp,jp2,j)*chi(jp2,j-1) -
     1                          bbbb(jp,jp2,j)*chi(jp2,j)
        z1 = z1 + abs(  aaaa(jp,jp2,j)*chi(jp2,j-1)  ) +
     1                          abs(  bbbb(jp,jp2,j)*chi(jp2,j)  )
  326   continue
c
        if (z0.le.rndeps*z1) ihs(jp) = ihsign(2)
        if (z0.lt.-rndeps*z1) ihs(jp) = ihsign(3)
  328   continue
c
c
c
c
        ij = j - 1
c
        if ( mimp .le. 0 ) then
          write(nprint,10304) ij,zrad,xnuel(1,j)*xzeff(1,j),xnuhyd(1,j),
     1                      (zd(j2),ihs(j2),j2=1,i1)
        elseif ( mimp .lt. 3 ) then
          write(nprint,10304) ij,zrad,xnuel(1,j)*xzeff(1,j),xnuhyd(1,j),
     1                      (zd(j2),ihs(j2),j2=1,i1),xzeff(2,j)
        elseif ( mimp .lt. 5 ) then
          write(nprint,12304) ij,zrad,xnuel(1,j)*xzeff(1,j),xnuhyd(1,j),
     1                      (zd(j2),ihs(j2),j2=1,i1),xzeff(2,j)
c       case of more than 4 impurities not yet implemented
        endif
c
c
        if (ilend) go to 335
        j = j + nskpi
        if (j.lt.mzones) go to 305
        ilend = .true.
        j = mzones
        go to 305
  335   continue
c
c       loop has ended
c
c
c       print totals
c
        if(mimp.lt.3) go to 340
        write(nprint,12305)     (ihpart,js=lhyd1,limpn)
        write(nprint,12306)     (zzni(js),js=lhyd1,limpn)
        go to 395
  340   continue
        write (nprint,10305)    (ihpart,js=lhyd1,limpn)
        write (nprint,10306)    (zzni(js),js=lhyd1,limpn)
c
  395   continue
c
c
cl              third page--energy fluxes
c
c
  400   continue
cbate        if (tai.le.0.0) return
        if (nlpomt(3)) go to 401
c
        lpage = lpage + 1
c
        write (nprint,10400) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
        write (nprint,10401) lhpowr, lhlen
  401   continue
        if(nadump(1).le.lcentr) then
        isep = mzones
        isepm1 = ledge
        else
        isep = nadump(1)
        isepm1 = isep - 1
        endif
c
c
        zohme = 0.0
        zchxe = 0.0
        zalfe = 0.0
        zrade = 0.0
        zexte = 0.0
        zeion = 0.0
c
        zvols = 2.0 * avi(mzones,12,1) * uisl**3
         zsurfi = avi(mzones,3,1) * avi(mzones,5,1)  ! plasma surface area
c
c
c               loop begins here
c
c
        jb = 2 + nskpi
        ibold = lcentr
        ilend = .false.
  405   continue
c
        zrad = avi(jb,15,1) * uiel
c
c
cl              compute fluxes
c
c  note <grad f> = (<|grad xi|**2>/<|grad xi|>) * d f / d xi,
c  as in sbrtn CONVRT (dps 09oct87)
c
        z0 = avi(jb,6,1) / (avi(jb,5,1)*dxboui(jb) * uisl)
c
        zgrdte = (tes(2,jb) - tes(2,jb-1))*z0
        zgrdti = (tis(2,jb) - tis(2,jb-1))*z0
c
c
c               energy fluxes (total)
c
      if ( linout(32) .lt. 0 ) then
c
c     This way of computing energy fluxes was in the original BALDUR code
c  It suffers from using the updated chi values which are out of synch
c  with the values of aaaa and bbbb.
c
        zeeflx = 0.0
        zeiflx = 0.0
c
        do 422 jp = 1, mchi
        zeeflx = zeeflx - aaaa(lelec,jp,jb) * chi(jp,jb-1) -
     1                          bbbb(lelec,jp,jb) * chi(jp,jb)
        zeiflx = zeiflx - aaaa(lion,jp,jb) * chi(jp,jb-1) -
     1                          bbbb(lion,jp,jb) * chi(jp,jb)
  422   continue
c
        zeeflx = zeeflx * avi(jb,3,1)*avi(jb,5,1) * uiep
        zeiflx = zeiflx * avi(jb,3,1)*avi(jb,5,1) * uiep
c
c
c               conductivities
c
c
        zicond=-ditis(jb) *zgrdti*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
        zecond=-detes(jb) *zgrdte*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
        zineut=-ditins(jb)*zgrdti*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
c
        zeconv = zeeflx - zecond
        ziconv = zeiflx - zicond - zineut
c
      else
c
c     What follows will be the new default mode.
c  The arrays condei(jb), condii(jb), cvctei(jb), cvctii(jb), ctotei(jb)
c  and ctotii(jb) are computed in sbrtn convrt and cnvcof in file DSOLVER.
c
        zeeflx = ctotei(jb) * uiep
        zeiflx = ctotii(jb) * uiep
        zecond = condei(jb) * uiep
        zicond = condii(jb) * uiep
        zeconv = cvctei(jb) * uiep
        ziconv = cvctii(jb) * uiep
c
        zineut=-ditins(jb)*zgrdti*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
c
      endif
c
  430   continue
c
c
c               cumulative sources/sinks
c               (sources/sinks integrated from 0 to r(jb))
c
c
        z1ohm = 0.0
        z1chx = 0.0
        z1alf = 0.0
        z1rad = 0.0
        z1ext = 0.0
        z1eion = 0.0
c
        iz = jb - 1
        if (ibold.gt.iz) go to 440
c
        do 438 jz = ibold, iz
c
        z1ohm = z1ohm + dx2i(jz) * weohms(jz)
        z1chx = z1chx - dx2i(jz)*(wiions(jz) + wichxs(jz))
        z1alf = z1alf + dx2i(jz) * (wealfs(jz) + wialfs(jz) +
     1     wed3fs(jz) + wid3fs(jz))
c
c   les  jan-91 added d-3he fusion
c
        z1ext = z1ext + dx2i(jz) * (webems(jz) + wibems(jz) +
     1      weauxs(jz)+wiauxs(jz)+weecrh(jz)+wiecrh(jz)+
     2      weicrf(jz)+wiicrf(jz))
        z1eion=z1eion+dx2i(jz)*cnueqs(jz)*(tes(2,jz)-tis(2,jz))
c
c   les  nov-90 added synchr rad
c
        z1rad = z1rad+dx2i(jz) * (webrs(jz)+weions(jz)+wesrs(jz)
     1    + wesyn(jz))
c
        if (mimp.le.0) go to 438
c
        do 436 ji = 1, mimp
        z1rad = z1rad + dx2i(jz) * weirs(ji,jz)
  436   continue
c
  438   continue
c
        zohme = zohme + z1ohm*zvols*usep
        zchxe = zchxe + z1chx*zvols*usep
        zalfe = zalfe + z1alf*zvols*usep
        zexte = zexte + z1ext*zvols*usep
        zrade = zrade + z1rad*zvols*usep
        zeion = zeion + z1eion*zvols*usep
c
        ibold = jb
c
  440   continue
c
        ztote = zohme + zalfe + zexte - zecond - zicond - zeconv -
     1          ziconv - zineut - zchxe - zrade
        if(iz.le.isepm1) zthdot=ztote
c
        zineut = zineut + zchxe
c
        ijb = jb - 2
        if (.not.nlpomt(3))
     1  write (nprint,10403) ijb, zrad, zecond, zeconv, zicond, ziconv,
     1                  zineut, zrade, zohme, zalfe, zexte, ztote,zeion
c
c
        if (ilend) go to 435
        jb = jb + nskpi
        if (jb.lt.mzones) go to 405
        ilend = .true.
        jb = mzones
        go to 405
  435   continue
c
cl      loop has ended
c
c
c
cl              fourth page-- usually omitted
c
c
  500   continue
        if(nlpomt(5)) go to 595
        lpage = lpage + 1
c
        write (nprint,10500) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
c
c       print-out of total impurity radiation and
c       current recycling coefficient
c
        z0 = avi(mzones,12,1) * uisl**3 * usep
        ztot=0.0
        do 515 jz=lcentr,ledge
        do 510 ji=1,mimp
        if(jz.eq.lcentr) zwimp(ji)=0.
        zwimp(ji)=zwimp(ji)+weirs(ji,jz)*dx2i(jz)
        if(jz.eq.ledge) zwimp(ji)=z0*zwimp(ji)
        if(jz.eq.ledge) ztot=ztot+zwimp(ji)
  510   continue
  515   continue
c
        ie=max0(1,mimp)
        zt = tai * uist
        call timint (zt,zztcold,bdtime,20,tcold,1,1)
        write(nprint,10515) zztcold,grecyc,(zwimp(l),l=1,ie),ztot
c
c
c      print out power flux from electrons to ions due to collisional
c       equilibration. define zei = cnueqs*(te-ti)
c
        zei(1)=0.0
        do 520 i=lcentr,mzones
        zei(i)=cnueqs(i)*(tes(2,i)-tis(2,i))
  520   continue
c
c       print out diffusion coefficients
c
cbate        write(nprint,10531)
c
        write(nprint,10532)
        do jz=lcentr,mzones
          write(nprint,10534) jz,ahalfs(jz,1),detepr(jz),ditipr(jz),
     &      vnwars(jz) + vnneo1(1,jz),vewars(jz) + vnneo1(2,jz),
     &          (dnhs(jh,jz),jh=1,mxhyd),(dnis(ji,jz),ji=1,mimp)
        enddo
c
c
  595   continue
c
c
cl      page of semi-empirical coefficients
c
      if ( lmpirc .and. ( .not. nlpomt(5) )
     &    .and. ( nstep .gt. 0 ) ) then
        lpage = lpage+1
c
        write(nprint,10700) label1(1:48),label5(1:72),
     &    lpage,nstep,zt*1000.0,zdt
c
        call empirc(3)
      endif
c
cl... 15.06 page of theoretical coefficients of theory
c
      if ( ltheor .and. ( .not. nlpomt(22) ) ) then
        lpage=lpage+1
c       call prtheory
        call ptheory(3)
      endif
c
c
cbate      call prtcfn
cbate      call prtsfrm
c
c
cl  fifth page -- particle source rates
c
c
  600   continue
c
      if ( .not. nlpomt(6) ) then
c
        lpage = lpage + 1
c
        write (nprint,10600) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
        write (nprint,'(1A)') 'Particle source rates 1'
        write (nprint,10606)
        do j=lcentr,mzones
          zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uiel
          zsrch1 = shions(1,j) + shbems(1,j) - shchxs(1,j)
     &      - shfus(1,j) + shblos(1,j)
          write (nprint,10607) j, zrad
     &      , zsrch1
     &      , siions(1,j)
     &      , shions(1,j), shbems(1,j), shchxs(1,j)
     &      , shfus(1,j),  shblos(1,j)
        enddo
!ap
        write (nprint,'(1A)') 'Particle source rates 2'
        write (nprint,10606)
        do j=lcentr,mzones
          zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uiel
          zsrch2 = shions(2,j) + shbems(2,j) - shchxs(2,j)
     &      - shfus(2,j) + shblos(2,j)
          write (nprint,10607) j, zrad
     &      , zsrch2
     &      , shions(2,j)
     &      , shions(2,j), shbems(2,j), shchxs(2,j)
     &      , shfus(2,j),  shblos(2,j)
        enddo
c
cbate        write(nprint,10610)
cbate        do 620 j=lcentr,mzones
cbate        write(nprint,10620) j, (bound(i,j),i=1,6)
cbate  620   continue
c
      endif
c
c..Long printout for magnetic diffusion equation
c
        if ( nlpomt(10) ) go to 695
c
        lpage = lpage + 1
c
        write (nprint,10600) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
        write (nprint,10650)
        write (nprint,10652)
c
c  ajboot(j)    = < bootstrap current density / R >     [A/m**3]
c  ajtpbi(j)    = < current density due to trapped
c                 particles in the banana regime
c                 / R > at zone center j                [A/m**3]
c  zrad         = plasma half-width                     [m]
c  eta(2,j)     = resistivity at zone centers           [sec]
c  zjtor        = <J_tor>                               [MA/m**2]
c  zcdrive      = beam and RF driven current density    [MA/m**2]
c  zboot        = bootstrap current density             [MA/m**2]
c  zjtpb        = <current density due to trapped
c                 particles
c  zconv        = convective contribution to
c                 current density due to motion
c                 of toroidal flux                      [Ma/m**2]
c               = \Partial{\psi_{pol}}{\rho}
c                 \Partial{\rho}{t}
c                 \langle 1/R^2 \rangle
c                 / 2 \pi \eta \langle 1/R \rangle
c  avi(j,9,1)   = R*B_tor                               [m*tesla]
c  z7r          = <1/R>                                 [1./m]
c  avi(j,10,1)  = <1/R**2>                              [1./m**2]
c  avi(j,7,1)   = < |del xi|**2 / R**2 >                [1./m**4]
c
c  zitorr = total toroidal current                      [MA]
c  zibeam = total beam-driven current                   [MA]
c  ziboot = total bootstrap current                     [MA]
c  zijtpb = total current due to trapped particles      [MA]
c
c  Note:  if ( lneocl(1) == 1 ) then
c    ajboot and ajtpbi are computed in sbrtn nclass_int
!cap
        if ( lneocl(1) .lt. 1 ) then
          ajboot = 0.0
          ajtpbi = 0.0
          call boots  ( ajboot,ajtpbi,lcentr,ledge)
        endif
c
        zibeam = 0.0
        ziboot = 0.0
        zijtpb = 0.0
c
c  zdarea = differential cross-sectional area [m**2]
c
        do 680 j=1,mzones-1
c
          zdarea = avi(j+1,13,1) - avi(j,13,1)
c
        zrad = 0.5 * ( avi(j+1,15,1) + avi(j,15,1) )
c
c        zjb7b = 1.e-6 *  avi(j,9,1) *
c     & ( avi(j+1,2,1) * avi(j+1,3,1) * r0ref * bpoli(j+1)
c     &          * avi(j+1,7,1) / avi(j+1,8,1)
c     & - avi(j,2,1) * avi(j,3,1) * r0ref * bpoli(j)
c     &          * avi(j,7,1)/avi(j,8,1) )
c     & / ( emu0 * ( avi(j+1,12,1) - avi(j,12,1) ) )
c     & / ( 0.5  * ( avi(j+1,11,1) + avi(j,11,1) ) )
c
c        zib7bt = zib7bt + zjb7b * zdarea
c
      zjtor = 1.e-6 * r0ref *
     & ( avi(j+1,2,1) * avi(j+1,3,1) * bpoli(j+1) * avi(j+1,7,1)
     &   - avi(j,2,1) * avi(j,3,1) * bpoli(j) * avi(j,7,1) )
     &   / ( twopi * emu0 * zdarea )
c
        zcdrive = 1.e-6 * ( cjbeam * ajbs(j) * usij
     &            + cdprof(j) )
        zibeam = zibeam + zcdrive * zdarea
c
        zboot = 1.e-6 * cfutz(480) * ajboot(j)
     &          / ( 0.5 * ( avi(j+1,11,1) + avi(j,11,1) ) )
        ziboot = ziboot + zboot * zdarea
c
        zjtpb = 1.e-6 * cfutz(480) * sfutz(16) * ajtpbi(j)
     &          / ( 0.5 * ( avi(j+1,11,1) + avi(j,11,1) ) )
        zijtpb = zijtpb + zjtpb * zdarea
c
        zconv = 1.e-6 * avi(j,10,1) * r0ref * (1. / eta(2,j))
     & * 0.5 * ( bpoli(j+1) * avi(j+1,1,2) / avi(j+1,11,1)
     &         + bpoli(j)   * avi(j,1,2)   / avi(j,11,1) )
c
        z7r = 0.5 * ( avi(j+1,11,1) + avi(j,11,1) )
c
        write (nprint,10658) j, zrad, eta(2,j)*usir, 1./ftrap(2,j)
     & , zjtor, zcdrive, zboot, zjtpb, zconv, avi(j,9,1),  z7r
     & , avi(j,7,1)
c
 680  continue
c
c  zitorr = total toroidal plasma current in MA
c
      zitorr = avi(mzones,2,1) * avi(mzones,3,1) * r0ref * bpoli(mzones)
     &  * avi(mzones,7,1) * 1.e-6 / (twopi * emu0)
c
        write (nprint,10659) zitorr, zibeam, ziboot, zijtpb
c
        lpage = lpage + 1
c
        write (nprint,10600) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
        iin = ( ntheta / 2 ) + 1   ! index of theta at inboard edge
c
        write (nprint,10660)
        do jz=3,mzones
          write (nprint,10662) avi(jz,15,1), trapbr(jz,1)
     &      ,  bpol2di(jz,1), bpol2di(jz,iin)
     &      ,  btor2di(jz,1), btor2di(jz,iin)
     &      ,  rthxbi(1,jz),  rthxbi(iin,jz)
        enddo
c
 695  continue
c
c
c
cl              sixth page--misc. quantities
c                         -- usually appears as fourth page
c
c
  700   continue
        if (nlpomt(4)) go to 799
        lpage = lpage + 1
c
        write (nprint,10700) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
c
        isep = mzones
        isepm1 = ledge
c
c       if scrapeoff is on, redefine some quantities
c
        if(nadump(1).le.lcentr) go to 702
        isep = nadump(1)
        isepm1 = isep - 1
        zzee = 0.0
        zzei = 0.0
        zzeb = 0.0
        zzne = 0.0
        zznion = 0.0
        zohme = 0.0
        zalfe = 0.0
        zexte = 0.0
c
        do 701 jz = lcentr,isepm1
        zzee = zzee + chi(lelec,jz) * dx2i(jz)
        zzei = zzei + chi(lion,jz) * dx2i(jz)
        zzeb = zzeb + hebems(jz)*rhobis(2,jz) * dx2i(jz)
        zzne = zzne + rhoels(2,jz) * dx2i(jz)
        zznion = zznion + rhoins(2,jz) * dx2i(jz)
        zohme = zohme + dx2i(jz) * weohms(jz)
c
c   les  nov-90 add d3he fusion power
c
        zalfe = zalfe + dx2i(jz) * (wealfs(jz) + wialfs(jz)
     1    + wed3fs(jz) + wid3fs(jz))
        zexte = zexte + dx2i(jz) * (webems(jz) + wibems(jz) +
     1      weauxs(jz)+wiauxs(jz)+weecrh(jz)+wiecrh(jz)+
     2      weicrf(jz)+wiicrf(jz))
  701   continue
c
        zzee = zzee * uisd * uiee
        zzei = zzei * uisd * uiee
        zzeb = zzeb * usee
        zohme = zohme * usep
        zalfe = zalfe * usep
        zexte = zexte * usep
c
c
cl              confinement times
c
c
  702   continue
        zpts = 1.0 / float(ipts)
        i = 1
        call resetr(zzni,mxchi,0.0)
        ix0=5*(mxchi+2)
        call resetr(ztconf,ix0,0.0)
        call resetr(zloss ,mxchi,0.0)
        call resetr(zlne  , 6,0.0)
c
        zzlne = 0.0
        do 718 jz = lcentr, isepm1
        zeirs = 0.0
        if (mimp.le.0) go to 705
        do 704 ji = 1, mimp
        zeirs = zeirs + weirs(ji,jz)
  704   continue
  705   continue
c
c   les  nov-90  add synchrotron loss
c
        zloss(lelec) = zloss(lelec) + usip*zvols*dx2i(jz)*
     1    (wesyn(jz) +wesrs(jz) +weions(jz) + webrs(jz) + zeirs)
        zloss(lion) = zloss(lion) - usip*zvols*dx2i(jz)*
     1                          (wiions(jz) + wichxs(jz))
c
c
c
        z0 = zvols * uisd * dx2i(jz)
        do 712 jp = 1, mchi
        zzni(jp) = zzni(jp) + z0*chi(jp,jz)
  712   continue
c
        zzlne = zzlne + rhoels(2,jz) * dx2i(jz) * 2.
c
        if (jz.lt.isepm1.and.
     1  (xbouni(jz+1)/xbouni(isep)+rndeps).lt.(float(i)*zpts)) go to 718
c
        zrconf(i) = avi(jz+1,15,1) * uiel  !  halfwidth in external units
        zeloss = 0.0
        zlne(i) = zzlne
c
        do 716 jp = 1, mchi
c
        zflux = 0.0
        do 714 jp2 = 1, mchi
        zflux = zflux - aaaa(jp,jp2,jz+1)*chi(jp2,jz)
     1                - bbbb(jp,jp2,jz+1)*chi(jp2,jz+1)
  714   continue
c
c                     surface area V'(xi)*<|del xi|> internal units
        z0 = (zloss(jp) + zflux*avi(jz+1,3,1)*avi(jz+1,5,1)) * ueit
        if (z0.gt.epslon) ztconf(jp,i) = zzni(jp) / z0
        if (jp.eq.lelec.or.jp.eq.lion) zeloss = zeloss + z0
  716   continue
c
        if (zeloss.gt.epslon) ztconf(mchi+1,i) =
     1  (zzni(lelec) + zzni(lion)) / zeloss
        zlne(i) = zlne(i) * ztconf(mchi+1,i) / xbouni(jz+1)**2
        i = i + 1
  718   continue
c
c
c       experimental times
c
        ztenr1 = 0.0
        if (zohme.gt.epslon)
     1          ztenr1 = (zzei + zzee) / zohme * usep*uset*uese
        ztenr2 = 0.0
        z0 = zohme + zexte + zalfe
        if (z0.gt.epslon)
     1          ztenr2 = (zzei + zzee + zzeb) / z0 * usep*uset*uese
c
c       average ni, ne, ti, te
c  15.07 First, clarify the scaling factor z01: without scrape-off, it
c  is just 1/volume; it multiplies the number of particles computed for
c  p.1. With scrape-off, just have to scale the density by the volume
c  within the separatrix.
c
        z01 = 1. / (avi(mzones,12,1) * uiel**3)
        if(nadump(1).gt.lcentr)
     1      z01 = 2.0 * used * avi(mzones,12,1) / avi(isep,12,1)
        znibar = zznion * z01
        znebar = zzne * z01
        ztibar = zzei / zznion * ueie*uieh
        ztebar = zzee / zzne * ueie*uieh
c
c       loop voltage
c
        zv = vloopi(ledge,2)
c
c       beta -- e, i, beam, total
c
        zbps=emu0*zi*ueii/(twopi*avi(isep,15,1))*uisb
        zpb = zbps**2 / (8.0*fcpi) * used*usee
        z1 = z01 / zpb
        zbetae = 0.66666666 * zzee * z1
        zbetai = 0.66666666 * zzei * z1
        zbetab = 0.66666666 * zzeb * z1
c
        zbeta = zbetae + zbetai + zbetab
c
c               total beta (include bz)
c
        z1 = (zbps / bzs)**2
        ztbete = zbetae * z1
        ztbeti = zbetai * z1
        ztbetb = zbetab * z1
        ztbeta  = zbeta  * z1
c
c       compute alpha beta
c
        zalfap = 0.0
        do 730 jz = lcentr,isepm1
        zalfap = zalfap + (alphai(jz) * ealfai(jz) + d3fast(jz)) *
     &    dx2i(jz)
  730   continue
c
        zalfap = zalfap * uisd * uise * .6666666666 * 2.
     1  * 8. * fcpi /zbps**2
        zalfat = zalfap * z1
        zbeta = zbeta + zalfap
        ztbeta = ztbeta + zalfat
c
c       internal inductance and lambda
c  ...Note: lambda here is just beta-pol + l-i/2 !
c  ...l-i = alint is computed in subroutine AVEPLA.
c
        zsum=0.
        zflux=0.
        do 732 jz=lcentr,isepm1
        zsum=zsum+bpols(2,jz)**2*dx2i(jz)    ! volume integral of bpols**2
        zflux = zflux + bpols(2,jz)*(avi(jz + 1,1,1) - avi(jz,1,1))
732     continue
        zlint=2.*zsum/bpols(1,isep)**2
        if (leqtyp.eq.0) alint=zlint
        zlambd=zbeta+alint/2.
c
c..1.e-8 converts poloidal flux from cm**2*gauss to m**2*tesla
c
        zflux = 1.e-8 * 2. * fcpi * r0ref * zflux * uisl**2
c
c       compute line averaged density
c  15.07 Make the definition consistent with uses elsewhere.
c
        zlined=0.0
        do 735 jz= lcentr , ledge
        zlined = zlined + rhoels(2,jz)
     1                     * (ahalfs(jz+1,1) - ahalfs(jz,1))
 735    continue
c
        zlined = zlined / ahalfs(isep,1)
c
c               print out results
c
        write (nprint,10701)
        write (nprint,10702) lhlen, (zrconf(j),j=1,ipts)
c
        do 737 jp = 1, mchi
        write (nprint,10703) lhspec(jp), lhtime,
     1                                          (ztconf(jp,j),j=1,ipts)
  737   continue
c
        write (nprint,10704) lhtime, (ztconf(mchi+1,j),j=1,ipts)
        write (nprint,10720) (zlne(i),i=1,ipts)
        write (nprint,10705) ztenr1, lhtime, ztenr2, lhtime
        write (nprint,10706)  znebar, lhdens,
     1          znibar, lhdens, zbetae, ztbete,
     2          ztebar, lhtemp, zbetai, ztbeti,
     3          ztibar, lhtemp, zbetab, ztbetb,
     4          zlined, lhdens, zalfap, zalfat,
     5           zbeta,  ztbeta
        write (nprint,10707) zv,zlambd,alint,zflux
c
c               d+d reaction neutrons
c
        zdds = adds + hdds + bbddrs
        zddtot = addtot + hddtot + bbdd
c
        write (nprint,10708) adds, hdds, bbddrs, zdds, addtot,
     1        hddtot, bbdd, zddtot
c
c
c               qtech (for beam heating only)
c
c       compute the total injected beam power
        zpowr = 0.0
        do 752 jb = 1, mhbeam
        if (libeam(jb).le.0) go to 752
        ib = libeam(jb)
        do 751 jfr = 1, mxhfr
        zpowr = zpowr +
     1          yfract(ib,jfr) * hibeam(jb)*hebeam(jb) / float(jfr)
  751   continue
  752   continue
        zpowr = zpowr * 10.0**(-fxes)/fces * uesi * uesh * usep
c
c       integrate out to separatrix to get total powers
c
        zphbem=0.
        zphrf=0.
        zphalf=0.
        zpath=0.
        zpathl=0.
        zpabt=0.
        zpabtl=0.
        do 754 jz=lcentr,isepm1
        zphbem=zphbem+dx2i(jz)*(webems(jz)+wibems(jz))
        zphrf=zphrf+dx2i(jz)*(weauxs(jz)+wiauxs(jz)+
     &    weicrf(jz)+wiicrf(jz)+weecrh(jz)+wiecrh(jz))
c
c   les  nov-90  add d3he fusion power
c
        zphalf=zphalf+dx2i(jz)*(wealfs(jz)+wialfs(jz) +
     &    wed3fs(jz) + wid3fs(jz))
        zpath=zpath+dx2i(jz)*afuses(jz)
        zpathl=zpathl+aoloss(jz)*dx2i(jz)*afuses(jz)
        zpabt=zpabt+dx2i(jz)*halfas(jz)
        zpabtl=zpabtl+aoloss(jz)*dx2i(jz)*halfas(jz)
754     continue
        zpbinj=bpinjs*usep
        zpbabs=bpabs*usep
        zpblos=bploss*usep
        zvols = 2.0 * avi(mzones,12,1) * uisl**3
        zphbem=zphbem*zvols*usep
        zphrf=zphrf*zvols*usep
        zphalf=zphalf*zvols*usep
        zpath=zpath*zvols*usep*efusei*uise
        zpathl=zpathl*zvols*usep*efusei*uise
        zpabt=zpabt*zvols*usep*efusei*uise
        zpabtl=zpabtl*zvols*usep*efusei*uise
c
c       calculate loss fractions
c
        if(zpbinj.ne.0.) then
        zlbem=(zpbinj-zpbabs+zpblos)/zpbinj
        else
        zlbem=0.
        if(zpblos+zphbem.ne.0.) zlbem=zpblos/(zpblos+zphbem)
        endif
        zlrf=0.
        if(cfutz(60).lt.1.) zlrf=max(0.,cfutz(60))
        zprfin=zphrf/(1.-zlrf)
        if(zpath.ne.0.) then
        zlath=zpathl/zpath
        else
        zlath=0.
        endif
        if(zpabt.ne.0.) then
        zlabt=zpabtl/zpabt
        else
        zlabt=0.
        endif
c
c       calculate stored energy time derivatives
c
        if(zpbinj.ne.0.) then
        zebdot=zpbinj*(1.-zlbem)-zphbem
        else
        zebdot=-(zpblos+zphbem)
        endif
        zeadot=zpath*(1.-zlath)+zpabt*(1.-zlabt)-zphalf
c
c       calculate the steady state beam power
c
        if(zpbinj.ne.0.) then
        zpbnow=zpbinj
        else
        zpbnow=zpblos+zphbem
        endif
        if(zphbem.gt.0.) then
        zpbstr=(zphbem+zphrf-zthdot-zeadot+zpabt*(1.-zlabt)
     1  -zprfin*(1.-zlrf))
     1  /((1.-zlbem)*(1.+zpabt*(1.-zlabt)/zphbem))
        else
        zpbstr=-1.
        endif
c
c       calculate the steady state rf power
c
        if(zpbstr.gt.0.) then
        zprstr=zprfin
        else
        zpbstr=0.
        zprstr=(zphbem+zphrf-zthdot-zeadot+zpabt*(1.-zlabt))
     1  /(1.-zlrf)
        if(zprfin.le.0.) zprstr=0.
        endif
c
c       calculate the q values
c
        zpss=zpbstr+zprstr
        zpcomp=zphbem+zpblos+zprfin-zthdot-zeadot
        if(zpbnow+zprfin.ne.0.) then
        zqinst=5.*(zpath+zpabt)/(zpbnow+zprfin)
        else
        zqinst=-1.
        endif
        zbfact=0.
        if(zphbem.ne.0.) zbfact=zpbstr*(1.-zlbem)/zphbem
        if(zpss.ne.0.) then
        zqss=5.*(zpath+zpabt*zbfact)/zpss
        else
        zqss=-1.
        endif
        if(zpcomp.ne.0.) then
        zqcomp=5.*(zpath+zpabt)/zpcomp
        else
        zqcomp=-1.
        endif
        zheat=zphbem+zphrf
        if(zheat.ne.0.) then
        zqdenm=(zpowr+zprfin)*(1.-zthdot/zheat)
        else
        zqdenm=0.
        endif
        if(zqdenm.ne.0.) then
        zqtech=5.*(zpath+zpabt)/zqdenm
        else
        zqtech=-1.
        endif
c
c       write out all the powers and losses
c
        write(nprint,756) zpbstr,zpbnow,zpbabs,zphbem,
     x  zebdot,zlbem
756     format(/,5x,'pbstr=',1pe9.2,3x,'pbnow=',1pe9.2,
     x  3x,'pbabs=',1pe9.2,3x,'phbem=',1pe9.2,3x,
     x  'ebdot=',1pe9.2,5x,'lbem=',0pf6.3)
        write(nprint,757) zprstr,zprfin,zphrf,zlrf,
     x  zthdot
757     format(5x,'prfstr=',1pe9.2,3x,'prfinj=',1pe9.2,
     x  3x,'phrf=',1pe9.2,5x,'lrf=',0pf6.3,15x,'ethdot=',
     x  1pe9.2)
        write(nprint,758) zpath,zpabt,zphalf,zeadot,
     x  zlath,zlabt
758     format(5x,'path=',1pe9.2,3x,'pabt=',1pe9.2,
     x  3x,'phalf=',1pe9.2,3x,'eadot=',1pe9.2,5x,
     x  'lath=',0pf6.3,3x,'labt=',0pf6.3)
c
c       write out the q values
c
        write(nprint,759) zqinst,zqss,zqcomp,zqtech
759     format(/,10x,'qinst=',0pf7.3,5x,'qss=',0pf7.3,5x,'qcomp=',
     x  0pf7.3,20x,'qtech=',0pf7.3)
c
c
c
c               extrapolation results
c
  760   continue
        if (lzerr.le.0.or.lperr.le.0.or.k.ne.2) go to 790
        izerr = lzerr - 1
        if (lperr.le.mxchi)
     1  write (nprint,10711) xerr, lhspec(lperr), izerr
        if (lperr.eq.mxchi+1)
     1  write (nprint,10711) xerr, ihbpol, izerr
  790   continue
c
c               maximum change over one timestep
c
        if (lzdel.le.0.or.lpdel.le.0.or.k.ne.2) go to 795
        izdel = lzdel - 1
        if (lpdel.le.mxchi)
     1  write (nprint,10712) xdel, lhspec(lpdel), izdel
        if (lpdel.eq.mxchi+1)
     1  write (nprint,10712) xdel, ihbpol, izdel
  795   continue
c
c
cl              conservation checks
c
c
        write (nprint,10701)
        write (nprint,10715) (lhspec(js),ccons(js),js=1,mchi)
        write (nprint,10716) bcons
c
  799   continue
c
c
cl       seventh page -- print out divertor region quantities
c                     -- usually appears as fifth page (if at all)
c
c
c
  800   continue
c
c       divertor region data
c
c       if divertor is turned off skip divertor print-out
        if(nadump(1).le.lcentr) go to 899
c
        lpage = lpage + 1
        write (nprint,10800) label1(1:48),label5(1:72),
     1          lpage,nstep,zt,zdt
c
c
        if(.not.(nlcomp.or.inital)) go to 810
        inital=.false.
        isep=nadump(1)
        zpmass=fcau*(10.0**fxau)
        zf0=avi(mzones,3,1) * avi(mzones,5,1)   ! plasma surface area
        zf1=2. * avi(mzones,12,1)               ! 2 * plasma volume
        zf2=zf1 * uisl**3                       ! 2*volume standard units
        zf3=avi(isep,3,1) * avi(isep,5,1) *usit
        zf4=zf0*xbouni(mzones)*usit
c
  810   continue
        call resetr(zscflx,mxions,0.0)
        call resetr(zdvflx,mxions,0.0)
        call resetr(zwaflx,mxions,0.0)
        zndiv = 0.0
c
c  flux of ions into the scrapeoff region and flux onto wall
c
        do 820 ip=1,mchi
        do 815 jh=1,mhyd
        zscflx(jh)=-(aaaa(jh,ip,isep)*chi(ip,isep-1)
     .              +bbbb(jh,ip,isep)*chi(ip,isep))*zf3+zscflx(jh)
        zwaflx(jh)=-(aaaa(jh,ip,mzones)*chi(ip,mzones-1)
     .              +bbbb(jh,ip,mzones)*chi(ip,mzones))*zf4+zwaflx(jh)
  815   continue
        if(mimp.le.0) go to 818
        do 817 ji=1,mimp
        ii=ji+lhydn
        zscflx(ii)=-(aaaa(ii,ip,isep)*chi(ip,isep-1)
     .              +bbbb(ii,ip,isep)*chi(ip,isep))*zf3+zscflx(ii)
        zwaflx(ii)=-(aaaa(ii,ip,mzones)*chi(ip,mzones-1)
     .              +bbbb(ii,ip,mzones)*chi(ip,mzones))*zf4+zwaflx(ii)
  817   continue
  818   continue
  820   continue
c
c  start of loop over scrapeoff region
c
        do 850 j=isep,ledge
        if(j.lt.nadump(3)) zl=cfutz(127)
        if(j.ge.nadump(3)) zl=cfutz(128)
        zhydm=ahmean(2,j)*zpmass
        zihydm=1./zhydm
        zvoli = avi(j+1,12,1) - avi(j,12,1)     ! vol(j+1) - vol(j)
c
c  inverse average ion mass
c
        z0=0.0
        z1=0.0
        do 821 jh=1,mhyd
        z0=z0+rhohs(jh,2,j)
        z1=z1+rhohs(jh,2,j)*aspec(jh)
  821   continue
        if(mimp.le.0) go to 825
        do 823 ji=1,mimp
        ii=ji+lhydn
        z0=z0+rhois(ji,2,j)
        z1=z1+rhois(ji,2,j)*aspec(ii)
  823   continue
  825   continue
        ziionm=z0/(z1*zpmass)
c
c
c  if nadump(4)=0, use sound speed model for plasma
c  flow into limiter/divertor
c  if nadump(4)=1, use model with neutral friction included
c
c  compute parallel flow velocity zvs of plasma into limiter/div.
c
c  sound speed model
c
        zvs=sqrt((tis(2,j)+tes(2,j))*ziionm)
        if(nadump(4).lt.epslon) go to 835
c  include charge exchange "friction" on plasma.
c
c  total density, average mass and temperature of the neutrals
c
        zsneu=0.0
        zmneu=0.0
        ztneu=0.0
        do 830 jh=1,mhyd
        zsneu=zsneu+rhons(jh,j)
        zmneu=zmneu+rhons(jh,j)*aspec(jh)
        ztneu=ztneu+rhons(jh,j)*tns(jh,j)
  830   continue
        z0=1./(zsneu+epslon)
        zmneu=z0*zmneu*zpmass
        ztneu=z0*ztneu
c
c  find average neutral velocity zv0
c
        zv0=0.0
        if(zsneu.gt.epslon) zv0=sqrt(3.*ztneu/zmneu)
c
c  compute average hydrogen ion velocity zvi
c
        zvi=sqrt(2.*tis(2,j)*zihydm)
c
c  compute average collision energy ze in terms of h-1
c
        zrelvv=zv0*zv0+zvi*zvi
        ze=0.5*zpmass*zrelvv
c
c  convert ze in units of ergs to units of evs
c
        ze=ze*evsinv
c
c  compute c-x cross-section zsig in cm**2
c
        zsig(j)=0.6937e-14*(1.-0.155*log10(ze))**2/
     x  (1.+0.1112e-14*ze**3.3)
c
c  compute charge exchange frequency zfreq
c
        zfreq=zsneu*zsig(j)*sqrt(zrelvv)
c
c  compute flow velocity of plasma along field zvs
c
        z00=sqrt(2.)*zvs
        zvs=z00*zvs/(z00+zl*zfreq)
c
c  compute mean free path for c-x
c
        zmfp(j)=zvi/(zfreq+epslon)
c
c  compute the confinement time due to c-x ztcx=zl**2*zfreq/(2zvi)
c
        ztcx(j)=zfreq*zl**2/(2.*zvi**2)
c
c  compute the confinement time for the ion sound speed model
c
        ztis(j)=sqrt(2.)*zl/zvi
c
c  compare the two times zrats(j)=ztcx(j)/ztis(j).  if zrats > or
c  = 1, neutral friction is important
c
        zrats(j)=ztcx(j)/ztis(j)
  835   continue
c       power flux to divertor from zone j
        zpdiv(j)=-zvoli*uiee*
     .  (chi(lelec,j)*scroff(lelec,j)+chi(lion,j)*scroff(lion,j))
c       change from watts to kwatts
        zpdiv(j)=zpdiv(j)*1.e-3
c
c       particle fluxes to the divertor from zone j
c
        z0=-zf2*dx2i(j)
        do 841 jh=1,mhyd
        z1=z0*rhohs(jh,2,j)*scroff(jh,j)
        zdvflx(jh)=z1+zdvflx(jh)
        zndiv(j)=z1+zndiv(j)
  841   continue
        if(mimp.le.0) go to 844
        do 843 ji=1,mimp
        ii=ji+lhydn
        z1=z0*rhois(ji,2,j)*scroff(ii,j)
        zdvflx(ii)=z1+zdvflx(ii)
        zndiv(j)=z1+zndiv(j)
  843   continue
  844   continue
c
c       line density  *** needs to be fixed  ***
c
        zlden(j)=rhoels(2,j) * (avi(j+1,15,1) - avi(j,15,1))
c
  850   continue
c
        if(mhyd.gt.1) go to 848
        do 846 ix=limpn,limp1,-1
        ixp1=ix+1
        ixp1=min0(ixp1,mxions)
        zdvflx(ixp1)=zdvflx(ix)
        zscflx(ixp1)=zscflx(ix)
        zwaflx(ixp1)=zwaflx(ix)
  846   continue
        zdvflx(2)=0.0
        zscflx(2)=0.0
        zwaflx(2)=0.0
  848   continue
c
        zpower=0.
        zrsep = avi(isep,15,1) * uisl
        if(nadump(3).lt.100) zrlim2 = avi(nadump(3),15,1) * uisl
c       zpdiv(isep)=0.
        zdflux=0.
        zlined=0.
c
        do 855 j=isep,ledge
 
        zpower=zpower+zpdiv(j)
        zdflux=zdflux+zndiv(j)
        zlined=zlined+zlden(j)
 
  855   continue
        write(nprint,10850)
        write(nprint,10851)zpower
        if(mimp.gt.4) go to 857
        if(mimp.gt.2) go to 856
        write(nprint,10852) (zdvflx(i),i=1,4)
        write(nprint,10853) zdflux
        write(nprint,10854) (zscflx(i),i=1,4)
        write(nprint,10870) (zwaflx(i),i=1,4)
        go to 858
  856   continue
        write(nprint,11852) (zdvflx(i),i=1,6)
        write(nprint,10853) zdflux
        write(nprint,11854) (zscflx(i),i=1,6)
        write(nprint,11870) (zwaflx(i),i=1,6)
        go to 858
  857   continue
c       although not needed in present versions, the case of more
c       than 4 impurity species has been implemented here.
        write(nprint,12852) (zdvflx(i),i=1,8)
        write(nprint,10853) zdflux
        write(nprint,12854) (zscflx(i),i=1,8)
        write(nprint,12870) (zwaflx(i),i=1,8)
  858   continue
        if(nadump(3).eq.100) write(nprint,10855) zlined,zrsep
        if(nadump(3).lt.100) write(nprint,10856) zlined,zrsep,zrlim2
        write(nprint,10857)
        zarray(7)=cfutz(386)
        zarray(8)=cfutz(387)
        zarray(9)=cfutz(388)
        zgas=float(ngas(1))
        do 860 j=isep,ledge
        if(cfutz(383).eq.0.0) go to 859
        call diver(j,tes(2,j),zt2,rhoels(2,j),zd2,zu1,zamach1,ztau2,
     x  rmajs,xzeff(1,j),q(j),zgas,zarray)
 859    continue
        idj=j-1
        zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uisl ! zone center radius
        zted=tes(2,j)*evsinv
        ztid=tis(2,j)*evsinv
        zt2=zt2*evsinv
        write(nprint,10858)idj,zrad,zpdiv(j),zndiv(j),zted,ztid,
     x  rhoels(2,j),rhoins(2,j),zu1,zamach1,ztau2,zt2,zd2
 
  860   continue
c
c  print out neutral friction model quantities. if nadump(4)=0,
c  skip this section
c
        if(nadump(4).eq.0) go to 899
c
        write(nprint,10860)
        do 870 j=isep,ledge
        idj=j-1
        zrad = 0.5*(avi(j+1,15,1)+avi(j,15,1))*uisl ! zone center radius
        write(nprint,10865)idj,zrad,zmfp(j),ztcx(j),ztis(j),
     x  zrats(j),zsig(j)
  870     continue
c
c       print scrapeoff quantities (first timesteps only)
c
        if(nstep.gt.2) go to 885
        do 880 j=1,mzones
        yscrof(j)=scroff(lhyd1,j)
        if(mimp.gt.0) zscrof(j)=scroff(limp1,j)
  880   continue
        call rarray('scroff-h',yscrof,mzones)
        if(mimp.gt.0) call rarray('scroff-i',zscrof,mzones)
  885   continue
c
c
  899   continue
c
c
cl              ninth page--sawtooth profiles and summary
c
c
  900   continue
        if(nlpomt(9)) go to 995
        if(ohravg(lcentr).eq.0.) go to 995
        lpage = lpage + 1
c
        write (nprint,10600) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
c       prepare to normalize oh(r)
c
        zsum=0.
        do 910 jz=lcentr,ledge
        zsum=zsum+dx2i(jz)*ohravg(jz)
910     continue
        zsum=2.*zsum
        if(zsum.eq.0.) go to 995
c
c       write column headers
c
        write(nprint,911)
911     format(/,1x,'jz  radius te-pre tepost te-avg ti-pre tipost',
     1  ' ti-avg jz-pre jzpost jz-avg oh(r)',/,7x,
     2  '(cm)  (kev)  (kev)  (kev)  (kev)  (kev)  (kev)',
     3  ' ka/cm2 ka/cm2 ka/cm2')
        ztdif=0.
        ijmix=0
c
        do 920 jz=lcentr,mzones
          iz=jz-1
          zrad = avi(jz,15,1) * uisl
          zraz = 0.5*(avi(jz+1,15,1)+avi(jz,15,1))*uisl ! zone center radius
          if ( zraz .le. 0.5*sawr1 ) then
            ztdif=ztdif+dx2i(jz)*(teavg(jz)-tiavg(jz))
            ijmix=jz
          end if
           zoh=ohravg(jz)/zsum
c
c       write out columns
c
        write(nprint,915) iz,zrad,tepre(jz),tepost(jz),teavg(jz),
     1  tipre(jz),tipost(jz),tiavg(jz),ajzpre(jz),ajzpst(jz),
     2  ajzavg(jz),zoh
915     format(1x,i2,1x,0pf7.2,9(0pf7.3),0pf6.2)
920     continue
        if ( ijmix .gt. 1 ) ztdif=ztdif*2./xbouni(ijmix+1)**2
c
c       write summary
c
        ztcdif=teavg(lcentr)-tiavg(lcentr)
        write(nprint,925)
925     format(/,2x,'tlast   q(0) period  tau-m   rq1  rinv  rmix',
     1  ' tauein poh  dub   te0   ti0  te-ti :avg  del:te   ne neut',
     2  ' pfus  alf  uth',/,3x,
     3  'sec',10x,'msec   msec    cm    cm    cm    sec',
     4  '   kw   kj   kev   kev   ev    ev,6x,17h %',
     5  '    %')
        write(nprint,926) tbsaws,sawqc,sawtau,tmguir,sawr1,sawri,sawrmx,
     1  tauein,ohptot,estpol,teavg(lcentr),tiavg(lcentr),ztcdif,ztdif,
     2  sawte,sawne,sawdd,sawfus,sawalf,sawuth
926     format(1x,0pf7.3,f6.3,3pf7.2,3pf7.1,3(0pf6.1),0pf6.3,
     1  -3pf6.0,-3pf5.0,2(0pf6.2),2(3pf6.0),2x,6(2pf5.1))
995     continue
c
c..temporary printout of neutral beam injection fast ion density (Bateman)
c
      if ( lbeams(32) .gt. 0 ) then
c
        lpage = lpage + 1
c
        write (nprint,10600) label1(1:48),label5(1:72),
     1          lpage,nstep,zt*1000.0,zdt
c
      write (nprint,19980)
19980 format (/2x,'radius',2x,'rhobis'/)
c
      do 990 jz=1,mzones
        do 990 jb=1,2
          write (nprint,19982) ahalfs(jz,jb),rhobis(jb,jz)
19982     format (0pf8.3,1pe14.4)
 990  continue
c
c..temporary printout of variables in FREYA computed in sbrtn DPOSIT
c
        write (nprint,10390)
10390   format (/6x,'radius',6x,'yrmajr',6x,'yllipt',6x,'yvols'
     &    ,7x,'yiota')
c
      do 392 jc=1,nyzone+1
        write (nprint,10395) yr(jc),yrmajr(jc),yllipt(jc),yvols(jc)
     &                    ,yiota(jc)
 392  continue
c
        write (nprint,10394)
10394   format (/6x,'yte   ',6x,'yrhoi ',6x,'yrhob '
     &    ,6x,'yrhoim',6x,'yqmean')
c
      do 396 jc=1,nyzone+1
        write (nprint,10395) yte(jc),yrhoi(jc),yrhob(jc)
     &     ,yrhoim(1,jc),yqmean(1,jc)
10395   format (1p5e12.4)
 396  continue
c
      endif
c
        return
c
c..@@@return
c
c
cl      format statements
c
10100   format(40x,'****************************************'//54x,
     1  'B A L D U R'/
     &  /39x,'1-1/2-Dimensional Tokamak Transport Code'/
     &  /49x,'G.Bateman, D.P.Stotler,'
     2  /19x,'D.Post, P.H.Rutherford, A.McKenney, G.Lister, M.Hughes,',
     3  ' L.Foote, R.Jensen, F.Seidl'
     4  /38x,'D.Heifetz, R.Hawryluk, W.Houlberg, D.Mikkelsen'
     5  /38x,'C.Singer, A.Silverman, J.Ogden, M.Redi'
     6  /30x,'Princeton Plasma physics',
     7  ' Laboratory -- Princeton, New Jersey'/
     8  40x,'****************************************')
10101   format(2x,'initial conditions:')
10102   format(//5x,'major radius - ',0pf9.2,1x,a10,5x,
     1          'minor radius - ',0pf9.2,1x,a10,5x,
     2          'mesh input')
10103   format(//5x,'major radius - ',0pf9.2,1x,a10,5x,
     1          'minor radius - ',0pf9.2,1x,a10,5x,
     2          'mesh equal spaced (nrfit=1)')
10104   format(//5x,'major radius - ',0pf9.2,1x,a10,5x,
     1          'minor radius - ',0pf9.2,1x,a10,5x,
     2          'mesh volume weighted (nrfit=2)')
10105   format(//5x,'major radius - ',0pf9.2,1x,a10,5x,
     1          'minor radius - ',0pf9.2,1x,a10,5x,
     2          'special mesh')
c
10110   format(/5x,'toroidal b - -',1pe10.3,1x,a10,5x,'toroidal i',
     1          ' - -',1pe10.3,1x,a10,5x,'j proportional to te**1.5')
10111   format(/5x,'toroidal b - -',1pe10.3,1x,a10,5x,'toroidal i',
     1  ' - -',1pe10.3,1x,a10,5x,'j proportional to (1-(r/a)**',
     2  0pf6.2,')**',0pf6.2)
c
10115   format(/5x,a2,'(0)  - - - -',1pe10.3,1x,a10,5x,a2,
     1          '(a)  - - - -',1pe10.3,1x,a10,5x,a2,'(r)-',a2,
     2          '(a) prop. to (1-(r/a)**',0pf6.2,')**',0pf6.2)
10116   format(/5x,a2,'(0)  - - - -',1pe10.3,1x,a10,5x,a2,
     1          '(a)  - - - -',1pe10.3,1x,a10,5x,a2,'(r) input')
10117   format(/5x,'ne(0)  - - - -',1pe10.3,1x,a10,5x,
     1          'ne(a)  - - - -',1pe10.3,1x,a10,5x,
     2          'computed from ion densities and z''s')
10118   format(/5x,a2,'(0)  - - - -',1pe10.3,1x,a10,5x,a2,
     1          '(a)  - - - -',1pe10.3,1x,a10,5x,
     2          'computed from species densities')
c
10120   format(1x)
c10121   format(5x,a10,' n(0)=',1pe10.3,1x,a10,', n(a)=',1pe10.3,1x,
c     1          a10,',',2pf6.2,' 0f plasma, atomic wt=',0pf6.1)
10122   format(5x,a10,' n(0)=',1pe10.3,1x,a10,', n(a)=',1pe10.3,1x,
     1          a10,',',2pf6.2,' 0f plasma, atomic wt=',0pf6.1,
     2          ', charge=',i3)
c
10125   format(/13x,'millisecs.',7x,'profile printouts every:',6x,
     1          'short printouts every:',4x,'profile plots every:'/
     2          2x,'initial t =',0pf9.3,19x,i5,' timesteps',10x,i5,
     3          ' timesteps',15x,i5,' timesteps'/2x,'initial dt=',
     4          0pf12.6,12x,0pf9.3,' millisecs',36x,f9.3,' millisecs'/
     5          2x,'minimum dt=',0pf12.6,17x,i4,' zones',/2x,
     6          'maximum dt=',0pf12.6,5x,'and at:',0pf9.3,' millisecs'
     7          ,29x,'and at:',0pf9.3,' millisecs')
10126   format(37x,0pf9.3,' millisecs',36x,0pf9.3,' millisecs')
10127   format(37x,0pf9.3,' millisecs')
10128   format(92x,0pf9.3,' millisecs')
c
cb 10130   format(/13x,'time',7x,'major radius',3x,'minor radius',5x,
cb      1          'current'/10x,'millisecs.',3(5x,a10)/1x)
cb 10131   format(10x,0pf10.3,' --- ',0pf9.2,'  ---  ',0pf8.2,'  --- ',
cb      1          1pe10.3)
c
10135   format(/2x,'transport model',i3,
     1          ' includes neoclassical ion-ion terms, plus:')
11135   format(/2x,'transport model',i3,
     1          ' includes neoclassical ion-ion terms according to '
     2 ,'hawryluk-hirshman (1978), plus:')
10136   format(5x,'xe = ',f6.2,' pseudocl(3) +',f6.2,' bohm(7) +',
     3          f6.2,' neoclass(11) +',f6.2,' trapped elec(13) +',
     4          f6.2,' drift wave(14)'/
     5          5x,'xi = ',f6.2,' pseudocl(4) +',f6.2,' bohm(8) +',
     6          f6.2,' neoclass(12) +',f6.2,' toroid rpl (140) +',
     7          f6.2,' banana  drift(16)'/
     8          5x,'dh = ',f6.2,' pseudocl(1) +',f6.2,' bohm(5) +',
     9          f6.2,' neoclass(9)  +',f6.2,' trapped elec(18) +',
     1          f6.2,' drift wave(17)'/
     2          5x,'di = ',f6.2,' pseudocl(2) +',f6.2,' bohm(6) +',
     3          f6.2,' neoclass(10)'/5x,f6.2,' ware pinch(19)')
10137   format(/5x
     &,'off-diagonal particle-diffusion coefficients zeroed out.'/)
10138   format(5x,'xi (ion-thermal conductivity) based on formulae by:'
     1,2x,'bolton-ware (1983).')
11138   format(5x,'xi (ion-thermal conductivity) based on formulae by:'
     1,2x,'chang-hinton (1982) with zero shafranov shift.')
12138   format(5x,'xi (ion-thermal conductivity) based on formulae by:'
     1,2x,'chang-hinton (1982) with mean shafranov shift of',2x,
     2f6.3)
13138   format(5x,'xi (ion-thermal conductivity) based on formulae by:'
     1,2x,'hinton-hazeltine (1976).')
14138   format(5x,'xi (ion-thermal conductivity) based on formulae by:'
     1,2x,'some prescribed empirical fit.')
c
10140   format(/5x,'radiation and ionization states are not computed')
10141   format(/5x,'radiation and ionization states are coronal equili'
     1          ,'brium values --',a60)
10142   format(/5x,'boundary impurity ',a7,' is fixed, other ',
     1          'densities and temperatures fixed ',a6,
     2          ' when gradients are negative')
10143   format(/5x,'dt is limited so change over one timestep is les',
     1          's than ',2pf6.2,' %')
10144   format(/5x,'dt is limited so extrapolated error is less than',
     1  2pf6.2,' %',', and change over one timestep is less than ',
     2  2pf6.2,' %')
10145   format(/5x,'equations ',a5,' include flux',a5,
     1          ' include density/energy sources/sinks',a5,
     2          ' include fixed sources/sinks')
10146   format(/5x,'equations solved are ',2pf4.0,'  0mplicit'/1h1)
c
10200   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10201   format (4x,'zone',4x,'radius',6x,'te',4x,'reg',6x,'ti',4x,
     1          'reg',6x,'ne',10x,'ni',9x,'vloop',8x,'jz',10x,'q',
     2          6x,'thermal beta')
cbate     2          7x,'total beta')
10202    format (t16,'cm',t24,'kev',t40,'kev',t52,'part/cu cm'
     & ,t64,'part/cu cm',t79,'volts',t88,'kamp/sq cm')
c10203   format(1x)
10204   format(4x,i4,3x,f7.2,2(2x,1pe10.3,1x,a2),6(2x,1pe10.3))
10205   format(18x,2(2x,a10,3x),2(3x,'particles'),14x,a10)
10206   format(6x,'* totals ** ',2(2x,1pe10.3,3x),2(2x,1pe10.3),
     &         2(14x,1pe10.3))
c
c10300   format(1h1,2x,a48,10x,a72//
c     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
c     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
c     2          'millisecs.')
11300   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
c10301   format(/16x,'ion densities'/2x,'zone',2x,'radius',2x,
c     &          ' nu-el',10x,'nu-h',1x,6(5x,a10))
11301   format(16x,'ion densities'/2x,'zone',2x,'radius',2x,
     &          ' nu-el',10x,'nu-h',1x,6(5x,a10))
12301   format(16x,'ion densities'/1x,'zone',2x,'radius',2x,
     &          ' nu-el',10x,'nu-h',1x,6(5x,a10))
cbate     1          '    nu-el     nu-h   ',2x,7(3x,a10))
10302   format(6x,a10,22x,5(4x,a10,1x))
12302   format(5x,a10,23x,6(2x,a10,1x))
c10303   format(1x)
10304   format(1x,i4,1x,f7.2,1x,1p2e12.4,6(3x,1pe10.3,1x,a1))
12304   format(1x,i3,1x,f7.2,1x,1p2e12.4,7(2x,1pe10.3,a1))
10305   format(37x,5(5x,a10))
12305   format(38x,6(3x,a10))
10306   format(4x,'** totals **',20x,5(5x,1pe10.3))
12306   format(4x,'** totals **',21x,6(3x,1pe10.3))
c10308   format (2x,'conservation checks ',5(5x,1pe10.3))
c
10400   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10401   format(/30x,'energy losses from plasma within radius r(j) in'
     1  ,a10/5x,'rad.',2x,'electron',2x,'electron',5x,'ion',7x,
     2  'ion',6x,'neutral',2x,'radiative',3x,'ohmic',5x,'alpha',5x,
     3  'other',5x,'total',5x,'e-i',/1x,'j',a10,'conduct',2x,
     4  'convect.',2x,'conduct.',2x,'convect.',3x,'losses',4x,'losses',
     5  4x,'heating',3x,'heating',3x,'heating',4x,'gain',5x,
     6  'coupling')
10403   format(1x,i2,1x,f5.1,1x,11(1x,1pe9.2))
c
10500   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10515 format(2x,
     1  ' tcold(kev)=',1pe9.2,3x,' recycling=',0pf6.3,5x,
     2  'total imp. rad.=',5(1pe10.2)/)
c10531   format(2x,' power flux zei from electrons to ions,',
c     1  1x,'and diffusion coefficients')
10532   format(/2x,'zone',4x,'radius',5x,
     1  'k-e totl',4x,'k-i totl'
     2  ,6x,'vnware',6x,'veware',4x,'d-h totl',16x,'d-i totl')
10534   format (2x,i4,2x,0pf7.2,1p9e12.3)
c
cbate 10600   format(1h1,2x,a48,10x,a72/
10600   format(/1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10606 format (/3x,'particle source rates  particles / (cm^3 sec)'
     &  /' zone',t8,'radius',6x,'S_Hyd'
     &  ,6x,'S_imp'
     &  ,6x,'shions',6x,'shbems',6x,'shchxs',6x,' shfus',6x,'shblos')
10607 format (1x,i4,1x,f7.3,1p9(e12.3:))
10610 format(/3x,' zone-centering parameter bound(ip,j), where ip=1,-'
     &,'---,6 refers to the ion and energy densities'/)
10620 format(3x,i4,10x,6(2x,1pe10.3))
c
10650 format (/4x,'zone',3x,'ahalf',8x,'eta',1x,'neocl/Spitz'
     & ,5x,'<J_tor>',5x,'<J_drive>',4x,'<J_boot>',5x,'<J_tpb>'
     & ,5x,'convect'
     & ,2x,'RB_tor',2x,'<1/R>',2x,'<dxi**2/R**2>')
10652 format (13x,'m',8x,'ohm*m',13x,5(5x,'MA/m**2'),3x,'m*tesla'
     & ,2x,'1./m',3x,'1/m**4')
10658 format (3x,i4,2x,0pf7.3,1pe12.3,1pe14.5,1p5e12.3
     & ,0pf8.3,0pf7.3,1pe11.3)
10659   format (/4x,'** totals **',26x,1p5e12.3)
10660   format (/t2,'radius',t12,'trapbr',t25,'bpol_out'
     &  ,t38,'bpol_in',t51,'btor_out',t64,'btor_in'
     &  ,t77,'R_out',t90,'R_in')
10662   format (t1,0pf7.3,t10,1p8e13.4)
c
10700   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10701   format(/1x)
10702   format(32x,'radius (',a10,'):',5(0pf8.2,5x))
10703   format(10x,a10,'  confinement time (',a10,'):',5(3x,1pe10.3))
10704   format(10x,'total energy confinement time(',a10,'):',
     1  5(3x,1pe10.3))
10705   format(/10x,'experimental energy confinement times:',
     1  /15x,29hthermal/ohmic  . . . . . . . ,1pe10.3,1x,a10,
     2  /15x,'total  . . . . . . . . . . . ',1pe10.3,1x,a10)
10706   format(/10x,'mean electron density . . . ',1pe10.3,2x,a10,22x,
     1  'beta-poloidal',2x,'beta-toroidal'/10x,
     2  'mean ion density  . . . . . ',1pe10.3,2x,a10,15x,
     3  'electron: ',1pe10.3,5x,1pe10.3/10x,
     4  'mean electron temperature',':  ',1pe10.3,2x,a10,15x,
     5  'ion:  . . ',1pe10.3,5x,1pe10.3/10x,
     6  'mean ion temperature  . . . ',1pe10.3,2x,a10,15x,
     7  'beam ion: ',1pe10.3,5x,1pe10.3/10x,
     8  'line avg. electron density.   ',1pe10.3,2x,a10,15x,
     9  'alpha:  . ',1pe10.3,5x,1pe10.3/75x,
     x  'total:  .   ',1pe10.3,5x,1pe10.3)
10707   format(/10x,'loop voltage: ',1pe10.3,1x,'volts',
     1          5x,'lambda=',0pf6.3,5x,'l-i=',0pf6.3,5x,
     2          'internal flux=',0pf6.3,2x,'volt-sec')
10708   format(/10x,'d+d reaction neutrons   ','thermonuclear',3x,
     1          'beam-plasma',5x,'beam-beam',8x,'total'/15x,
     2          'neutrons / sec.',4(5x,1pe10.3)/15x,
     3          'neutrons (total)',4x,1pe10.3,3(5x,1pe10.3))
10711   format(/29x,'maximum corrected error was', 2pf7.3,
     1          ' %, for ',a10,', at zone',i4)
10712   format(/29x,'maximum change over one timestep was', 2pf7.3,
     1          ' %, for ',a10,', at zone',i4)
10715   format(17x,a10,' conservation: ',1pe10.3)
10716   format(10x,'b-poloidal energy conservation: ',1pe10.3)
10720   format(37x,'n-tau(energy) :',5(3x,1pe10.3))
c
10800   format(1h1,2x,a48,10x,a72/
     1  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     1          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     2          'millisecs.')
10850   format(/6x,'scrapeoff region')
10851   format(/5x,'power(kw)=',1pe10.2)
10852 format(2x,' ion fluxes into divertor  (per sec):',1x,' h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3)
11852 format(2x,' ion fluxes into divertor  (per sec):',1x,' h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx =',
     3  1pe10.3)
12852 format(2x,' ion fluxes into divertor  (per sec):',1x,9h h1 flx =,
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx  ',
     3  1pe10.3,2x,' i5 flx =',1pe10.3,2x,' i6 flx =',1pe10.3)
10853   format(2x,' total ion flux into the divertor (per sec):',1x,
     1  ' tot flx =',1pe10.3)
10854 format(2x,' ion fluxes into scrapeoff (per sec):',1x,' h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3)
11854 format(2x,' ion fluxes into scrapeoff (per sec):',1x,' h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx =',
     3  1pe10.3)
12854 format(2x,' ion fluxes into scrapeoff (per sec):',1x,' h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx =',
     3  1pe10.3,2x,' i5 flx =',1pe10.3,2x,' i6 flx =',1pe10.3)
10855   format(/5x,'line den(ions/cm**2)=',1pe10.2,45x,
     1 '     separatrix=',1pe10.3,' cm ')
10856   format(/5x,'line den(ions/cm**2)=',1pe10.2,10x,'lim. radius=',
     1  e12.4,10x,' cm ',5x,'virt. lim. radius=',e12.4,' cm ')
10857   format(/5x,' zone','    r(cm) ',' power(kw)','  ions/sec',
     1  '   te(ev) ','   ti(ev) ','      ne  ','      ni  ',
     1  'u1(cm/sec)','  amach1  ','   tau2   ','   t2(ev) ',
     1  '     d2   ')
10858   format(5x,i5,f10.2,1p11e10.2)
10860   format(/5x,' zone','  r(cm)   ','c-x mfp cm',
     x  ' tau c-x s','  tau in s',' tcx/tis  ','sig cx cm**2')
10865   format(5x,i5,f10.2,1p5e10.2)
10870 format(2x,' ion fluxes onto the wall  (per sec):,1x,9h h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3)
11870 format(2x,' ion fluxes onto the wall  (per sec):,1x,9h h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx =',
     3  1pe10.3)
12870 format(2x,' ion fluxes onto the wall  (per sec):,1x,9h h1 flx =',
     1  1pe10.3,2x,' h2 flx =',1pe10.3,2x,' i1 flx =',1pe10.3,2x,
     2  ' i2 flx =',1pe10.3/40x,' i3 flx =',1pe10.3,2x,' i4 flx =',
     3  1pe10.3,2x,' i5 flx =',1pe10.3,2x,' i6 flx =',1pe10.3)
        end
c
c--------1---------2---------3---------4---------5---------6---------7-c

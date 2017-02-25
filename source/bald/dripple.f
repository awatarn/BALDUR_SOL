c  20:00 16-feb-92 .../baldur/code/bald/dripple.f
c/ 20.05 13:00 29-jul-91 /11040/bald91/wbaldn1 DRIPPLE, Bateman, PPPL
c  BALDUR  file DRIPPLE   Bateman, Stotler PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type           (use appropriate date for yymmdd)
c cfs get /11040/bald92/byymmdd/wcode.tar
c tar xf wcode.tar  (this creates directory .../code and subdirectories)
c cd code/bald      (this changes to subdirectory .../code/bald)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c       Contents of file DRIPPLE:
c  TRIPL1 - compute ripple-induced transport coefficients
c  AVERGE - integrations with respect to poloidal angle for ripple comp
c  INTERP - interpolation on rectangular grid of TF ripple amplitudes
c
c----------------------------------------------------------------------c
c@tripl1  /11040/bald92/wbaldn1 DRIPPLE
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c rgb 16-feb-92 defined zasp32, zpmass, zseprx, zsepr0, and iseprx
c   and zaspinv, ported over from dtransp.f
c rgb 20.05 28-jul-91 write sbrtn tripl1 to compute ripple transport
c
c**********************************************************************c
c
      subroutine tripl1
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
c..temporary common block for ripple induced transport
c
      common /cmrpl1/ dripl1(mj), dripl2(mj), dripl3(mj), dripl4(mj)
     &  , ykineo(mj), yripl0(mj)
c
      common/apolar/ yripla,yms00,yriple,xpoint,ypoint,galpha,
     1               x00,y00,imax0,jmax0
c
      dimension
     &   yripla(mj)   , yms00(mj)    , yriple(21,21),
     &   xpoint(21)   , ypoint(21)   , galpha(mj) 
c
      logical           inital       , lriple       , lgalfa
c
      data    inital,lriple /.true.,.false./
c
      data    ikirip,ncoil1,idelt0,idelta /140,141,142,143/,
     1          n1rple,ibeta1,ncoil2,n2rple /144,145,146,147/,
     2          infile,igtrpl,iaveth /148,149,150/
c
c               parameters for activating and adjusting
c               the level of field-ripple effects are
c               stored in cfutz(140) through cfutz(150)
c
c       two ripple components are handled:  (1) a steady component due
c       to the discreteness of the tf coils, and (2) a prescribed var-
c       iable component which can be introduced to control a burn.  all
c       ripple amplitudes are peak to average
c
c        cfutz(ikirip)=0. ensures that all the toroidal-ripple effects
c                      are set to zero
c        cfutz(ikirip)=1. activates the three principal effects aris-
c                      ing from the discrete number of main tf coils
c        cfutz(ikirip)=2. activates only the ripple-trapping effect
c        cfutz(ikirip)=3. activates only the ripple-plateau transport
c        cfutz(ikirip)=4. activates only the banana-drift transport
c        cfutz(ikirip)=5. activates only ripple-trapping and ripple-
c                      plateau effects, not banana-drift transport
c        cfutz(ncoil1)=the number of main toroidal-field coils; de-
c                      faulted to 12.0
c        cfutz(idelt0)=the amplitude at the plasma center (r=0) of the
c                      steady ripple component; no default but ignored
c                      when input ripple-data file is used
c        cfutz(idelta)=the amplitude at the plasma edge or separatrix
c                      (r=a) of the steady component; no default but
c                      ignored when input ripple-data file is used
c        cfutz(n1rple)=exponent of "r/a" in an empirical fit of the
c                      change in amplitude with minor radius for the
c                      steady-ripple component; defaulted to 2.0
c        cfutz(ibeta1)=the factor "beta" in an empirical approximation
c                      for the variation in the steady component with
c                      poloidal angle "theta", namely, exp(-beta*theta*
c                      theta); defaulted to 0.2
c        cfutz(ncoil2)=the effective number of tf coils generating the
c                      prescribed variable ripple, which can be intro-
c                      duced to control a burn; if set to 0.0, the var-
c                      iable ripple component is neglected
c        cfutz(n2rple)=exponent of "r/a" in an empirical fit of the
c                      change in amplitude with minor radius for the
c                      variable-ripple component; defaulted to 2.0
c        cfutz(infile)=1. activates the acquisition from an input data
c                      file of the steady-component ripple amplitude as
c                      a function of position; defaulted to 0.0 which
c                      means reverting to the empirical formula.  when
c                      a data file is used, cfutz elements with argu-
c                      ments idelt0, idelta, n1rple, and ibeta1 are ig-
c                      nored; however, the variable ripple component is
c                      always accessible in the same way
c        cfutz(igtrpl)=1. activates the goldston-towner method of hand-
c                      ling the spatial variations of the tf ripple
c                      amplitude
c        cfutz(iaveth)=time (secs) between successive re-evaluations
c                      of weakly time-dependent integrals over the
c                      poloidal angle theta which determine the ripple-
c                      trapping coefficient, g(alpha,--); defaulted
c                      to 5.*dtmax
c
c      data    ikirip,ncoil1,idelt0,idelta /140,141,142,143/,
c     1          n1rple,ibeta1,ncoil2,n2rple /144,145,146,147/,
c     2          infile,igtrpl,iaveth /148,149,150/
c
c       note:   ikirip must also be initialized in subroutine auxval
c
c       note:   whenever the variable tf ripple component is active, the
c       array delimiter mxt1 must be reset from mxt=20 to mxt/2=10, be-
c       cause variable ripple parameters are stored in the last 10 ele-
c       ments of tcompi(20), rcurri(20), redgi(20), and rmji(20).  these
c       parameters define the time variation of the variable component
c       and are set by prescribing values in the main namelist for tb-
c       poid(11-20), bpoid(11-20), curent(11-20), and rmajor(11-20)
c
c**********************************************************************c
c
cap   data inital /.true./
c
      if ( inital ) then
c
        inital = .false.
c
        z0=2.0*(fxc-fxes)+fxau
        zgyro0=(2.0*fcc*fcc*fcau)*(10.0**z0)/(fces*fces)
        zneu0=1.3293404
        mxt1=mxt
        if(cfutz(ikirip).gt.epslon) mxt1=10
        if(cfutz(ncoil1).le.epslon) cfutz(ncoil1)=12.0
        if(cfutz(n1rple).le.epslon) cfutz(n1rple)=2.0
        if(cfutz(ibeta1).le.epslon) cfutz(ibeta1)=0.2
        if(cfutz(n2rple).le.epslon) cfutz(n2rple)=2.0
        if((cfutz(iaveth).le.epslon).and.
     &    (cfutz(ikirip).gt.epslon)) cfutz(iaveth)=5.*dtmax
        if(cfutz(infile).gt.epslon .and. nstep.eq.1)
     &    call link("unit20=(ripple,open)//")
c
c
      if(cfutz(ikirip).gt.epslon) lriple=.true.
      coils=cfutz(ncoil1)
      exbeta=cfutz(ibeta1)
      nriple=0
      if(cfutz(ikirip).gt.epslon) nriple=int(cfutz(ikirip)+0.1)
      kriple=0
      if(cfutz(igtrpl).gt.epslon) kriple=int(cfutz(igtrpl)+0.1)
      zdti=cfutz(iaveth)*usit
      zti=tbi
      x00=0.0
      y00=0.0
      imax0=20
      jmax0=20
c       zshrnk and zexpnd are used to map data for a plasma of
c       elliptic cross-section onto an equivalent circle
      zmap=1.0
      if(ellipt.gt.epslon) zmap=sqrt(ellipt)
      zshrnk=1.0/zmap
      zexpnd=zmap
c
      i1 = mxhyd * mxzone
      i2 = mximp * mxzone
      i3 = mximp * i1
      i4 = mximp * i2
      i5 = mxhyd * i1
      i6 = i1 + i2
c
      call resetr(yripla,mzones,0.0)
      call resetr(yms00,mzones,0.0)
      call resetr(galpha,mzones,0.0)
c
      do 270 jj=1,21
        xpoint(jj)=0.0
        ypoint(jj)=0.0
      do 270 ii=1,21
        yriple(ii,jj)=0.0
  270 continue
c
c
c       initialization of time-independent quantities
c       involved in toroidal-field ripple effects
c
      if(.not.lriple) go to 540
      lriple=.false.
      lgalfa=.true.
      crple1=(fcc/fces)*(10.0**(fxc-fxes))
      crple2=32.880
      crple3=3.0*sqrt(fcpi/2.0)
      crple4=40./fcpi
      zmxrho = ahalfs(iseprx,1)
      jtmin=mxt1+2
      eps=1.e-04
      mmax=1000
c
c       branch depending on whether the tf ripple amplitude is acquir-
c       ed from an input data file or set up from an empirical formula
c
      if(cfutz(infile).le.epslon) go to 490
c
c       the input data file "ripple" is assumed to contain:
c       (1) yriple(i,j), the peak-to-average tf ripple amplitude
c       in percent at x-y points indexed by i and j;
c       (1) xpoint(i), ypoint(j), the x-y coordinates in meters;
c       (2) rcentr, the distance (meters) of the x-y origin from
c       the center-line of the toroidal magnetic field;
c       (3) y00, the vertical displacement (meters) of the x-y
c       origin above a horizontal plane, usually the plane of mag-
c       netic symmetry; and
c       (4) imax0, jmax0, the maximum number of x and y points re-
c       spectively
c
      rewind 20
      read(20,410) rcentr,y00,imax0,jmax0
      if((imax0.lt.22).and.(jmax0.lt.22)) go to 350
      write(nprint,420) imax0,jmax0
      nlend=.true.
      call endrun
  350 continue
c
c       load arrays yriple(i,j), xpoint(i), and ypoint(j) from input
c       file ripple.  note:  the format statement 430 must correspond
c       with the format of the input file
c
      read(20,430) (xpoint(ii),ii=1,imax0)
      read(20,430) (ypoint(ii),ii=1,jmax0)
      do 360 jj=1,jmax0
      read(20,430) (yriple(ii,jj),ii=1,imax0)
  360 continue
c
c       re-express distances in centimeters and fractionize yriple(i,j)
c
      rcentr=100.*rcentr
      y00=100.*y00
      do 370 i=1,imax0
      xpoint(i)=100.*xpoint(i)
  370 continue
      do 380 j=1,jmax0
      ypoint(j)=100.*ypoint(j)
  380 continue
      do 390 j=1,jmax0
      do 390 i=1,imax0
      yriple(i,j)=1.e-02*yriple(i,j)
  390 continue
c
c       printout of ripple file information
c
      write (nprint,401) label1(1:48)
      call rarray(' xpoint ',xpoint,21)
      call rarray(' ypoint ',ypoint,21)
      write(nprint,440) rcentr,y00
      do 400 jj=1,21
      write(nprint,450) (yriple(ii,jj),ii=1,21)
  400 continue
c
  401 format(1h1,2x,a48//
     1 29x,'quantities having to do with ripple data file')
  410 format(1x,2(e13.6,2x),2(i3,2x))
  420 format(/2x,' ripple input file overflows allotted storage',
     & ' in yriple(21,21)',3x,' imax0=',i3,1x,' jmax0=',i3)
  430 format(1x,8(e13.6,2x))
  440 format(/3x,' yriple(i,j), i increases across,',
     & ' j increases down:',3x,' rcentr=',1pe11.3,3x,' y00=',
     & 1pe11.3)
  450 format(/3(1x,1p10e12.3/))
c
c
c       interpolation on yriple(i,j) to load yripla(i).  the array
c       yripla(i) is loaded with peak-to-average values of the
c       steady ripple component evaluated along a line pointing
c       in the theta-zero direction
c
      rdx=float(imax0-1)/(xpoint(imax0)-xpoint(1))
      rdy=float(jmax0-1)/(ypoint(jmax0)-ypoint(1))
      x00=rcentr-rmajs
      j0=int((y00-ypoint(1))*rdy)+1
      if((y00.ge.ypoint(j0)).and.(y00.lt.ypoint(j0+1))) go to 460
      if(y00.lt.ypoint(j0)) j0=max0(j0-1,1)
      if(y00.ge.ypoint(j0+1)) j0=min0(j0+1,jmax0)
  460 continue
      j0p1=min0(j0+1,jmax0)
      zb=(y00-ypoint(j0))*rdy
      do 480 iz=lcentr,mzones
      x0 = zshrnk * ahalfs(iz,1) - x00
      i0=int((x0-xpoint(1))*rdx)+1
      if((x0.ge.xpoint(i0)).and.(x0.lt.xpoint(i0+1))) go to 470
      if(x0.lt.xpoint(i0)) i0=max0(i0-1,1)
      if(x0.ge.xpoint(i0+1)) i0=min0(i0+1,imax0)
  470 continue
      i0p1=min0(i0+1,imax0)
      za=(x0-xpoint(i0))*rdx
      zc=za*zb
      yripla(iz)=(1.0+zc-(za+zb))*yriple(i0,j0)
     1                   +(za-zc)*yriple(i0p1,j0)
     2                   +(zb-zc)*yriple(i0,j0p1)
     3                        +zc*yriple(i0p1,j0p1)
  480 continue
c
c       the array yms00(j) is used to store, as a function of radial
c       index j, the values of mean-square tf ripple amplitude.  these
c       values are calculated in subroutine averge(eps,mmax), where eps
c       specifies the relative error being sought and mmax denotes the
c       maximum number of iterations allowed.  the mean is taken with
c       respect to the poloidal angle theta
c
      call averge(eps,mmax)
c
      go to 540
c
  490 continue
c
c       midplane amplitude of the steady-ripple component using an
c       empirical formula
c
      do 500 i=lcentr,mzones
      zradal=1.0
      zexpon=cfutz(n1rple)
      if(abs(zexpon).gt.epslon) zradal=(zsepr0*ahalfs(i,1))**(zexpon)
      yripla(i)=cfutz(idelt0)+(cfutz(idelta)-cfutz(idelt0))*zradal
  500 continue
c
c       factor for the poloidal averaging of the squared steady-
c       component amplitude
c
      z0=fcpi*sqrt(2.0*cfutz(ibeta1))
      if(z0.gt.2.0) go to 510
      zms00=1.0-.278443269*z0
      go to 520
  510 continue
      zms00=.8862269255/z0
  520 continue
      do 530 jy=1,mzones
      yms00(jy)=zms00*yripla(jy)*yripla(jy)
  530 continue
  540 continue
c
      do 542 j = 1, mzones
        yripl0(j) = yripla(j)
 542  continue
c
      endif
c
c ---- end of initialization ---
c
c
c       when ripple-trapping effects are included (ripple-plateau and
c       banana-drift are not involved here), the ripple-trapping coef-
c       ficient g(alpha,--) for each radial zone is determined at this
c       point ahead of the main radial do-loop
c
      go to (552,552,554,554,552), nriple
      go to 554
c
  552 continue
c
c       g(alpha,--)=galpha(j) is determined by one of two methods:
c       (1) the uckan-tsang-callen method (ornl/tm-5438, june 1976)
c       which is activated when kriple=0; or, (2) the goldston-
c       towner method (pppl-1638r, feb. 1980) activated when kriple
c       =1.  note that initially lgalfa is set to be .true.
c
      if ( tbi .gt. zti ) lgalfa=.true.
      if ( lgalfa ) then
        lgalfa = .false.
        zti    = tbi+zdti
        call intgrl(kriple,coils,exbeta,0.0,0.0,eps,infile,mmax)
        if(nstep.lt.3) call rarray('galpha-j',galpha,mzones)
      endif
c
 554  continue
c
c
c---------------start of main do-loop over the radial index jz----------
c
c
      do 1540 jz = lcentr, mzones
c
c
c..k-i due to toroidal-field ripple
c
c
      if((cfutz(ikirip).gt.epslon)) go to 700
c
      yripla(jz)=0.0
      dripl1(jz)=0.0
      dripl2(jz)=0.0
      dripl3(jz)=0.0
      dripl4(jz)=0.0
      go to 870
  700 continue
c
c
      dripl1(jz)=0.0
      if(cfutz(ncoil2).le.epslon) go to 790
c
c       midplane amplitude of the variable-ripple component
c
      i=mxt1+1
      do 710 jt=jtmin,mxt
      if(tcompi(jt).le.0.0) go to 720
      i=jt
      if(tcompi(jt).eq.tbi) go to 720
      if(tcompi(jt).gt.tbi) go to 730
  710 continue
c
  720 continue
      zdel0=rcurri(i)
      zdela=redgi(i)
      zbeta1=rmji(i)
      zexpon=cfutz(n2rple)
      go to 750
c
  730 continue
      if(rcurri(i).ne.rcurri(i-1)) go to 740
      if(redgi(i).ne.redgi(i-1)) go to 740
      go to 720
c
  740 continue
      z0=1.0/(tcompi(i)-tcompi(i-1))
      z1=tcompi(i)-tbi
      z2=tbi-tcompi(i-1)
      zdel0=z0*(z1*rcurri(i-1)+z2*rcurri(i))
      zdela=z0*(z1*redgi(i-1)+z2*redgi(i))
      zbeta1=z0*(z1*rmji(i-1)+z2*rmji(i))
      zexpon=cfutz(n2rple)
c
  750 continue
c
c       "zsepr0" is a multiplying factor used to obtain "r/a" from
c       r/a = ahalfs(jz,1) / ahalfs(iseprx,1) = zsepr0 * ahalfs(jz,1)
c       where "a" here denotes the halfwidth
c       to the plasma edge or separatrix
c
      iseprx=mzones
      if(nadump(1).gt.lcentr) iseprx=nadump(1)
      zsepr0 = 1.0 / ahalfs(iseprx,1)
      zseprx = ahalfs(jz,1) / ahalfs(iseprx,1)
c
      zradal=1.0
      if(abs(zexpon).gt.epslon) zradal=zseprx**(zexpon)
      dripl1(jz)=zdel0+(zdela-zdel0)*zradal
c
c       poloidal average of the square of the variable-ripple
c       amplitude, <dripl1*dripl1>, assuming an approximate pol-
c       oidal dependence exp(-2.*zbeta1*theta*theta), where "theta"
c       is the poloidal angle.  note that initially zms0=0.0 and
c       zbeta0=0.0
c
      if(zms0.le.epslon) go to 760
      if(zbeta1.eq.zbeta0) go to 780
  760 continue
      zbeta0=zbeta1
      z0=fcpi*sqrt(2.0*zbeta0)
      if(z0.gt.2.0) go to 770
      zms0=1.0-.278443269*z0
      go to 780
  770 continue
      zms0=.8862269255/z0
  780 continue
c
  790 continue
c
      zmsqr0=yms00(jz)
      zmsqr1=0.0
      if(cfutz(ncoil2).gt.epslon) zmsqr1=zms0*dripl1(jz)*dripl1(jz)
c
      dripl2(jz)=0.0
      dripl3(jz)=0.0
      dripl4(jz)=0.0
      go to (810,810,830,830,810), nriple
c
c       ripple-trapping transport
c
  810 continue
c
      z1=crple1*tis(1,jz)/(rmajs*bzs+epslon)
      z1=z1*z1
      z1=crple2*z1*tions(1,jz)
c
c       the ripple-trapping coefficient g(alpha,--) is determined
c       either by the goldston-towner method or the uckan-tsang-
c       callen method.  the former is described in pppl-1638r, feb.,
c       1980; the latter is described in ornl/tm-5438, june 1976.
c       the call to execute the determination is made ahead of the
c       main do-loop, namely, ahead of "do 1540 --"
c
      dripl2(jz)=galpha(jz)*z1*yripla(jz)*sqrt(yripla(jz))
      if(cfutz(ncoil2).le.epslon) go to 820
      zalpha=bpols(1,jz)/(cfutz(ncoil2)*bzs*dripl1(jz)+epslon)
c
c       g(alpha,--) in the variable-ripple case is determined by the
c       uckan-tsang-callen method only
c
      xx=zalpha
      yy=zbeta1
      yy1=sqrt(yy)
      if(xx.gt.2.0) go to 816
      xx2=xx*xx
      zg=-1.641*xx+0.234*xx2-
     1  yy1*(1.000-0.243*xx-0.13*xx2)
      zg=10.0**zg
      go to 818
  816 continue
      xx2=xx-2.0
      z0=2.346+0.462*yy1
      zg=-z0*((1.0+0.64370*xx2)**.3127)
      zg=10.0**zg
  818 continue
      dripl2(jz)=zg*z1*dripl1(jz)*sqrt(dripl1(jz))+dripl2(jz)
  820 continue
      if(nriple.eq.2) go to 870
c
c       factors shared by ripple-plateau and banana-drift transport
c
  830 continue
c       ave ion mass
c
      zpmass=fcau*(10.0**fxau)
      z0=zpmass*aimass(1,jz)
c
c       ave ion thermal speed
c
      z1=sqrt(tis(1,jz)/z0)
c
c       effective poloidal ion gyroradius
c
      z2=crple1*z0*z1/(bpols(1,jz)+epslon)
      z2=min(z2,zmxrho)
      if(nriple.eq.4) go to 850
c
c       ripple-plateau transport (boozer 1979)
c       (assuming ambipolarity)
c
      z0=crple3*z1*z2*z2/rmajs
      dripl3(jz)=z0*cfutz(ncoil1)*zmsqr0
      if(cfutz(ncoil2).le.epslon) go to 840
      dripl3(jz)=z0*cfutz(ncoil2)*zmsqr1+dripl3(jz)
  840 continue
      if(nriple.ne.1) go to 870
c
c       ripple banana-drift transport (boozer 1980)
c       (assuming ambipolarity)
c
  850 continue
      z3=z1*z2
      z4=q(jz)*q(jz)*q(jz)*rmajs*rmajs
c
c..inverse aspect ratio
c  zaspinv = half-width / major radius to midpoint at zone bndries
c
      j0=jz
      if(xbouni(j0).le.epslon) j0=lcentr+1
        zaspinv = abs( ahalfs(j0,1) / rmids(j0,1) )
      zasp32 = sqrt(zaspinv) * zaspinv
c
      z0=crple4*zasp32*z3*z3*tions(1,jz)
      dripl4(jz)=z0*zmsqr0/(z4*cfutz(ncoil1)+epslon)
      if(cfutz(ncoil2).le.epslon) go to 860
      dripl4(jz)=z0*zmsqr1/(z4*cfutz(ncoil2)+epslon)+dripl4(jz)
  860 continue
  870 continue
c
 1540 continue
c
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@averge  /11040/bald92/wbaldn1 DRIPPLE
c       aes 29-apr-82 corrected use of zshrnk, zexpnd
c       aes 29-oct-81 array dimensions 52 --> 55 in common/apolar/
c       fgps 3-sep-80 subroutine set up (update 29-oct-80)
c
c**********************************************************************c
c                     
        subroutine averge(eps,mmax)
c
c
cl      2.21 integrations with respect to the poloidal angle
c
c       this subroutine carries out two kinds of integration with re-
c       spect to the poloidal angle, theta.  (1) the first has to do
c       with calculating the poloidal average of the squared tf ripple
c       amplitude and is performed only during the initialization of
c       a run or a re-initialization.  each active radial zone, j, is
c       handled and the results are loaded into the array yms00(j).
c       (2) the second integration, executed again for all active zones
c       but at periodic intervals, is used for approximating goldston-
c       towner's generalization of connor and hastie's ripple-trapping
c       coefficient, g(alpha,--).  this result is stored in galpha(j).
c       both integration procedures involve linear interpolations on a
c       given ripple-amplitude grid yriple(i,j) defined on a rectangu-
c       lar coordinate system, xpoint(i), ypoint(j).  two interpolations
c       are executed in succession by calling interp(deltf1,deltf2,x1,
c       y1,x2,y2), where deltf1 and deltf2 contain the results
c
c******************************************************************************
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
        logical         inital
c
c
        dimension
     r   yripla(mj)   , yms00(mj)    , yriple(21,21), xpoint(21)   ,
     r   ypoint(21)   , galpha(mj)   , zx(3)        , zw(3)
c
        common/apolar/ yripla,yms00,yriple,xpoint,ypoint,galpha,
     .                 x00,y00,imax0,jmax0
c
c******************************************************************************
c
c       the quadratures follow a program written by h. h. towner of
c       ppl, and use a 6-point gauss-legendre procedure, where zx(j)
c       are the roots of the legendre polynomial and zw(j) contains
c       weight factors.
c
        data
     1   zx /.238619186083197,.661209386466265,.932469514203152/,
     2   zw /.467913934572691,.360761573048139,.171324492379170/,
     3   z0 /6.283185307/,
     4   inital / .true. /
c
c
c       zshrnk and zexpnd are used to map data for a plasma of ellip-
c       tical cross-section onto an equivalent circle
c
        zmap=1.0
        if(ellipt.gt.epslon) zmap=sqrt(ellipt)
        zshrnk=1.0/zmap
        zexpnd=zmap
c
c       determination of the poloidal average of the squared ripple
c       amplitude.  the integration range is from zero to 2.*pi and
c       i equals the number of subintervals into which the range is
c       divided
c
      do 240 jr=lcentr,mzones
        i=1
        test2=0.0
  200   continue
        xlow=0.0
        test1=test2
        sum=0.0
      do 220 k=1,i
        xhigh=xlow+z0/float(i)
        c1=.5*(xhigh-xlow)
        c0=c1+xlow
      do 210 j=1,3
        z=c1*zx(j)
        r = ahalfs(jz,1)
        x1=zshrnk*r*cos(+z+c0)-x00
        y1=zexpnd*r*sin(+z+c0)-y00
        x2=zshrnk*r*cos(-z+c0)-x00
        y2=zexpnd*r*sin(-z+c0)-y00
        call interp(deltf1,deltf2,x1,y1,x2,y2)
        sum=c1*zw(j)*(deltf1*deltf1+deltf2*deltf2)+sum
210   continue
        xlow=xhigh
220   continue
        test2=sum
        re=abs((test2-test1)/(test2+epslon))
        if(re.le.eps) go to 230
        i=i*2
        if(i.le.mmax) go to 200
        i=i/2
        write(nprint,250) jr,i,re
  230   continue
c
        yms00(jr)=test2/z0
c
240   continue
c
        call rarray('msq rple',yms00,mzones)
c
        return
c
  250 format(1x,' yms00(j) at radial index',i3,2x,' after',
     & i5,' attempts the integral has not converged',3x,
     & ' relative error=',1pe10.3/)
c
c******************************************************************************
c
        entry intgrl(k0,coils,exbeta,a,b,eps,infile,mmax)
c
c       branch depending on how g(alpha,--)=galpha(j) is determined:
c       (1) by the uckan-tsang-callen method (ornl/tm-5438, june 1976);
c       or (2) by the goldston-towner method (pppl-1638r, feb. 1980)
c
c       alfrpl is the ratio of toroidal drift to ripple drift in
c       mid-plane, and in its calculation (inverse aspect ratio/q)
c       is replaced by the equivalent (b-poloidal)/(b-toroidal)
c
        if(k0.gt.0) go to 300
c
        yy1=sqrt(exbeta)
      do 290 jz=lcentr,mzones
        alfrpl=bpols(1,jz)/(coils*bzs*yripla(jz)+epslon)
        xx=alfrpl
        if(xx.gt.2.0) go to 270
        xx2=xx*xx
        zg=-1.641*xx+0.234*xx2-
     .  yy1*(1.000-0.243*xx-0.013*xx2)
        go to 280
  270   continue
        xx2=xx-2.0
        z0=2.346+0.462*yy1
        zg=-z0*((1.0+0.64370*xx2)**.3127)
  280   continue
        galpha(jz)=10.0**zg
290   continue
c
      return
c
300   continue
c
c       the following algorithm is set up to calculate goldston and
c       towner's generalization of the connor-hastie ripple-trapping
c       coefficient g(alpha,--)=galpha(j).  the algorithm involves an
c       integration from theta=a to theta=b, default limits being a=0.0
c       and b=z0=2.*pi.  intermediate values of theta, as well as sub-
c       intervals, are determined in the algorithm itself.  the index
c       i counts the number of subintervals
c
c       initializations
c
        if(.not.inital) go to 310
        inital=.false.
        if(a.le.epslon) a=0.0
        if(b.le.epslon) b=z0
        c=b-a
        z1=4.220232732e-02
        z2=4.80168702
        z3=z1*z2
c
c       precautionary re-initialization of zshrnk and zexpnd
c
        zmap=1.0
        if(ellipt.gt.epslon) zmap=sqrt(ellipt)
        zshrnk=1.0/zmap
        zexpnd=zmap
  310   continue
c
      do 400 jz=lcentr,mzones
        alfrpl=bpols(1,jz)/(coils*bzs*yripla(jz)+epslon)
        i=1
        test2=0.0
  320   continue
        xlow=a
        test1=test2
        sum=0.0
      do 350 k=1,i
        xhigh=xlow+c/float(i)
        c1=.5*(xhigh-xlow)
        c0=c1+xlow
      do 340 j=1,3
        z=c1*zx(j)
        theta1=+z+c0
        theta2=-z+c0
        sine1=sin(theta1)
        cosin1=cos(theta1)
        sine2=sin(theta2)
        cosin2=cos(theta2)
        r = ahalfs(jz,1)
        x1=zshrnk*r*cosin1-x00
        y1=zexpnd*r*sine1-y00
        x2=zshrnk*r*cosin2-x00
        y2=zexpnd*r*sine2-y00
        if(cfutz(infile).gt..9) go to 3202
        if(theta1.gt.fcpi) theta1=z0-theta1
        if(theta2.gt.fcpi) theta2=z0-theta2
        ratio1=exp(-exbeta*theta1*theta1)
        ratio2=exp(-exbeta*theta2*theta2)
        go to 3204
 3202   continue
        call interp(deltf1,deltf2,x1,y1,x2,y2)
        za=1.0/(yripla(jz)+epslon)
        ratio1=za*deltf1
        ratio2=za*deltf2
 3204   continue
        astar1=(alfrpl*abs(sine1))/(ratio1+epslon)
        astar2=(alfrpl*abs(sine2))/(ratio2+epslon)
        astar1=abs(astar1)
        astar2=abs(astar2)
        arcos1=0.0
        if(astar1.lt.1.0) arcos1=acos(astar1)
        arcos2=0.0
        if(astar2.lt.1.0) arcos2=acos(astar2)
c
c       note:  the inverse aspect ratio has been set up and stored
c       in array gx(mj), see trcoef and common/commsh/
c
        if(arcos1.gt.0.0) go to 324
  322   continue
        zf1=0.0
        go to 326
  324   continue
        za=sqrt(1.0-astar1*astar1)
        zb=za-astar1*arcos1
        if(zb.le.0.0) go to 322
        zb=zb*ratio1
        zc=zb**1.5
        zf1=zc*sine1*sine1*arcos1*(1.0+gx(jz)*cosin1)
  326   continue
        if(arcos2.gt.0.0) go to 330
  328   continue
        zf2=0.0
        go to 332
  330   continue
        za=sqrt(1.0-astar2*astar2)
        zb=za-astar2*arcos2
        if(zb.le.0.0) go to 328
        zb=zb*ratio2
        zc=zb**1.5
        zf2=zc*sine2*sine2*arcos2*(1.0+gx(jz)*cosin2)
  332   continue
        sum=c1*zw(j)*(zf1+zf2)+sum
340   continue
        xlow=xhigh
350   continue
        test2=sum
        re=abs((test2-test1)/(test2+epslon))
        if(re.le.eps) go to 380
        i=i*2
        if(i.le.mmax) go to 320
        i=i/2
        write(nprint,360) jz,i,re
  380   continue
c
        galpha(jz)=z3*test2
c
400   continue
c
        return
c
  360 format(1x,' galpha at radial index=',i3,2x,' after',
     & i5,' attempts the integral has not converged',3x,
     & ' relative error=',1pe10.3/)
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
        subroutine interp(del1,del2,x1,y1,x2,y2)
c
c
cl      2.22 interpolation on a rectangular grid of tf ripple amplitudes
c
c
c       aes 29-oct-81 array dimensions 52 -->55 in common/apolar/
c       fgps 5-sep-80 subroutine set up
c
c**********************************************************************c
!cap
	include 'cparm.m'
        logical         inital
c
c
        real ::
     r   yripla(mj)   , yms00(mj)    , yriple(21,21), xpoint(21)   ,
     r   ypoint(21)   , galpha(mj)   , zx(3)        , zw(3)
c
        common/apolar/ yripla,yms00,yriple,xpoint,ypoint,galpha,
     .                 x00,y00,imax0,jmax0
c
c******************************************************************************
c
c
        data            inital / .true. /
c
c       initializations
c
        if(.not.inital) go to 190
        inital=.false.
        rdx=float(imax0-1)/(xpoint(imax0)-xpoint(1))
        rdy=float(jmax0-1)/(ypoint(jmax0)-ypoint(1))
  190   continue
c
c       double linear interpolation
c
        i1=int((x1-xpoint(1))*rdx)+1
        j1=int((y1-ypoint(1))*rdy)+1
        i2=int((x2-xpoint(1))*rdx)+1
        j2=int((y2-ypoint(1))*rdy)+1
        if((x1.ge.xpoint(i1)).and.(x1.lt.xpoint(i1+1))) go to 200
        if(x1.lt.xpoint(i1)) i1=max0(i1-1,1)
        if(x1.ge.xpoint(i1+1)) i1=min0(i1+1,imax0)
  200   continue
        if((y1.ge.ypoint(j1)).and.(y1.lt.ypoint(j1+1))) go to 210
        if(y1.lt.ypoint(j1)) j1=max0(j1-1,1)
        if(y1.ge.ypoint(j1+1)) j1=min0(j1+1,jmax0)
  210   continue
        if((x2.ge.xpoint(i2)).and.(x2.lt.xpoint(i2+1))) go to 220
        if(x2.lt.xpoint(i2)) i2=max0(i2-1,1)
        if(x2.ge.xpoint(i2+1)) i2=min0(i2+1,imax0)
  220   continue
        if((y2.ge.ypoint(j2)).and.(y2.lt.ypoint(j2+1))) go to 230
        if(y2.lt.ypoint(j2)) j2=max0(j2-1,1)
        if(y2.ge.ypoint(j2+1)) j2=min0(j2+1,jmax0)
  230   continue
        i1p1=min0(i1+1,imax0)
        j1p1=min0(j1+1,jmax0)
        i2p1=min0(i2+1,imax0)
        j2p1=min0(j2+1,jmax0)
        za=(x1-xpoint(i1))*rdx
        zb=(y1-ypoint(j1))*rdy
        zc=za*zb
        del1=(1.0+zc-(za+zb))*yriple(i1,j1)
     .               +(za-zc)*yriple(i1p1,j1)
     .               +(zb-zc)*yriple(i1,j1p1)
     .                    +zc*yriple(i1p1,j1p1)
        za=(x2-xpoint(i2))*rdx
        zb=(y2-ypoint(j2))*rdy
        zc=za*zb
        del2=(1.0+zc-(za+zb))*yriple(i2,j2)
     .               +(za-zc)*yriple(i2p1,j2)
     .               +(zb-zc)*yriple(i2,j2p1)
     .                    +zc*yriple(i2p1,j2p1)
c
        return
        end

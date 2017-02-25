c@synch2  /11040/baldur/code/bald/dimprad.f
c
      subroutine synch2 (nstep,nzesv,nzosv,nwsv,wstops)
c
c       routine for rapid estimation of energy transport by
c       synchrotron radiation in a magnetized plasma
c       for details of the physical model and the equations
c       see report SAI-023-81-189-lj.
c
c       written by:  s. tamor, sai   phone (619)456-6466
c
c        date of this revision:  9/17/82
c
c
c    This routine provides an approximate solution of the steady
c     state transport equation for synchrotron radiation in a
c     plasma assuming plasma properties to be constant on a set of
c     nested surfaces.
c     Input variables are the volume, surface area,
c     temperature, density, and average magnetic field in each
c     zone, and the surface area and reflection coefficients
c     of the surrounding wall.
c     The routine returns the energy loss per unit
c     volume for each zone. Data are supplied and results
c     returned via the common block "COMSYN".
c
c
c
c        Input variables:
c          nzons  ..............number of zones in the problem (.le. jx)
c          ree,roo,roe,reo .........matrix of reflection coefficients
c                                    giving the reflected ordinary
c                                    and extraordinary mode intensities
c                                    in terms of the incident intensities
c          bmin  ................... smallest value of magnetic field
c          bmax .....................largest value of magnetic field
c          areaw ....................area of reflecting wall
c                                    (default=surf(nzons))
c         stemp(i) ..................zone temperatures
c         sdens(i) ..................zone densities
c          bavg(i) ..................average magnetic field in each zone
c          svol(i) ...................zone volumes
c          surf(i) ..................surface area of outer boundary
c                                    of each zone
c          ndbug .........=0  no debug edits
c                         =1  edits soutine synch for each frequency
c                         =2  also edits opacty for each frequency
c          iedit..........=0  no final edit
c                         =1  causes summary edit of output data
c
c        Output variables
c          wsyn(i) .................zone energy loss rate per unit volume
c          ptot .....................volume integrated energy loss
c          pwall ....................energy loss to wall (energy check)
c
c
c         Units:
c           length.................cm
c           energy & temperature...kev
c           magnetic field.........kilogauss
c
c
c     common/insink/ nzons,ree,roo,reo,roe,bmin,bmax,pwall,ptot,areaw,
c    2  taucrt,nw1,ndbug,iedit,
c    2  stemp(jx),sdens(jx),bavg(jx),svol(jx),surf(jx),wsyn(jx)
!cap
      include 'cparm.m'
      include 'csynrd.m'
c
      common/syntst/ wsyn1(mj),wsyne(mj),wsyno(mj)
      common/alsync/  wfreq,alphae(mj),alphao(mj)
      dimension delta(mj)
      data pi/ 3.14159265/
      istop=0
      testi=0.0
      testi0=0.0
      pwall=0.0
      dwght=1.0
      fwght=3.0
      ts1=surf(1)
      wsyn(1)=0.0
      wsyn1(1)=0.0
      wsyne(1)=0.0
      wsyno(1)=0.0
      delta(1)=2.*svol(1)/surf(1)
c
c        we here calculate the average thickness of each zone.
c
c
      do 100 nz=2,nzons
      ts2=surf(nz)
      delta(nz)=2.*svol(nz)/(ts1+ts2)
      ts1=ts2
      wsyn(nz)=0.0
      wsyn1(nz)=0.0
      wsyne(nz)=0.0
      wsyno(nz)=0.0
  100 continue
      if(areaw.lt.surf(nzons))  areaw=surf(nzons)
c
c        we now choose the lowest frequency and the interval for the
c        frequency integration.  units = radians/sec.
c         wmin=2.8 min omega_ce
c            wmax=14 min omega_ce
c
      wmin=5.0e10*bmin
      wmax=2.4e11*bmax
      delw=0.08*(wmax-wmin)
c
c        we now start main frequency loop.
c      save diagnostic params
      nzesv=0
      nzosv=0
      nwsv=0
      wstops=0.
c
      do 500 nw=1,100
      nw1=nw
c     xe=0.0
c     xo=0.0
c     srce=0.0
c     srco=0.0
c     emass=0.0
c     omass=0.0
c     areae=0.0
c     areao=0.0
      wfreq=wmin+delw*(nw-1)
      dwght=-dwght
      wght=fwght+dwght
      if(nw.eq.1) wght=1.
      wght=wght*0.333333*delw
c
c        les 1990 - test to remove jumps in wsyn
      nze0=0
      nzo0=0
      istep=0
c
c        blackb is the rayleigh-jeans intensity divided by the
c        stemperature. =omega**2/((2.*pi)**3*c**2)
c
      blackb=4.48e-24*wfreq*wfreq
c
c        we will now call for a table of values of the absorption
c        coefficients in all zones. the quantities returned need
c        to be normalized by wmult=(omega^p)**2/(c*omega^c).
c
      call syopac2 (nzons,stemp,bavg,ndbug)
      taue=0.0
      tauo=0.0
      nze=nzons
      nzo=nzons
c
c        we now march in from the outer wall calculating the accumulated
c        optical thickness  for each plasma mode. the
c        edge of the thick region is defined by the zone at which
c        the thickness exceeds taucrt.
c
      do 150 ii=1,nzons
      nz=nzons-ii+1
      wmult=6.03e-12*sdens(nz)/bavg(nz)
      alphae(nz)=alphae(nz)*wmult
      alphao(nz)=alphao(nz)*wmult
      taue=taue+alphae(nz)*delta(nz)
      if(taue.lt.taucrt) nze=nz-1
      tauo=tauo+alphao(nz)*delta(nz)
      if(tauo.lt. taucrt) nzo=nz-1
  150 continue
c
c        we will now calculate the coefficients and source terms
c        in the equations for the self consistent e wave and o wave
c        intensities.
c
c  les 1990
      testi0=testi
 151  continue
c
      xe=0.0
      emass=0.0
      srce=0.0
      areae=0.0
      if(nze.eq.0) go to 155
      areae=surf(nze)
      srce=areae*stemp(nze)
      do 153 nz=1,nze
      emass=emass+sdens(nz)*svol(nz)
  153 continue
  155 continue
      xo=0.0
      omass=0.0
      srco=0.0
      areao=0.0
      if(nzo.eq. 0) go to 160
      areao=surf(nzo)
      srco=areao*stemp(nzo)
      do 158 nz=1,nzo
      omass=omass+sdens(nz)*svol(nz)
  158 continue
  160 continue
      nze1=nze+1
      nzo1=nzo+1
      if(nze1.gt.nzons) go to 175
      do 170 nz=nze1,nzons
      ts=4.*svol(nz)*alphae(nz)
      xe=xe+ts
      srce=srce+ts*stemp(nz)
  170 continue
  175 continue
      if(nzo1.gt.nzons) go to 185
      do 180 nz=nzo1,nzons
      ts=4.*svol(nz)*alphao(nz)
      xo=xo+ts
      srco=srco+ts*stemp(nz)
  180 continue
  185 continue
c         save params  - added 10/88
      nzesv=max0(nzesv,nze)
      nzosv=max0(nzosv,nzo)
c
c        we now solve the pair of algebraic equations for the self
c        consistent intensities eray and oray. the quantities
c        qee, qoo, qeo, qoe are the elements of the transposed
c        cofactor matrix.
c
      qee=areao+areaw*(1.-roo)+xo
      qoo=areae+areaw*(1.-ree)+xe
      qoe=roe*areaw
      qeo=reo*areaw
c   les 1990; fnorm sensitive to ree=roo=0.95, reo=roe=0.05, low T
c     fnorm=blackb/(qee*qoo-qeo*qoe)
      fnorm=blackb/(qoo*((qee-qeo)+(qeo/qoo)*(qoo-qoe)))
      eray=fnorm*(qee*srce+qeo*srco)
      oray=fnorm*(qoe*srce+qoo*srco)
c
c        testi keeps track of the maximum value of eray. the
c        frequency integration is terminated when eray .lt. 0.03*testi.
c
      testi=max(testi0,eray)
      if(eray.lt.0.03*testi) istop=1
c
c   les 1990; check optically thick region to avoid jumps in wsyn
c
      if(nze.eq.0) go to 199
      if(istep.ge.1) go to 189
      nze0=nze
      nzo0=nzo
 189  istep=istep+1
      nzeol=nze
      nzool=nzo
c     test eray
      if(nze.ge.nzons) go to 190
      psyntn=4.*alphae(nze)
      psyntk=areae*sdens(nze)/emass
      if(psyntn.gt.psyntk*1.0) nze=nze+1
      if((psyntk-psyntn).gt.psyntk*0.3) nze=nze+1
c     test oray
 190  if(nzo.eq.0) go to 192
      if(nzo.ge.nzons) go to 192
      psyntn=4.*alphao(nzo)
      psyntk=areao*sdens(nzo)/omass
      if(psyntn.gt.psyntk*1.0) nzo=nzo+1
      if((psyntk-psyntn).gt.psyntk*0.3) nzo=nzo+1
c     did nze or nzo change?
 192  if(nze.eq.nzeol.and.nzo.eq.nzool) go to 199
c        allow change of one grid spacing only
      if(istep.lt.21) go to 151
      nze=nzeol
      nzo=nzool
c      write() 'nz should be larger'
      go to 199
c
 199  continue
c
c        accumulate contribution from current frequency to power
c        loss from each zone in optically thin region
c
      if(nze1.gt.nzons) go to 205
      do 200 nz=nze1,nzons
      ts=4.*pi*alphae(nz)*(blackb*stemp(nz)-eray)*wght
      wsyn(nz)=wsyn(nz)+ts
      wsyn1(nz)=wsyn1(nz)+ts
  200 continue
  205 continue
      if(nzo1.gt.nzons) go to 225
      do 220 nz=nzo1,nzons
      ts=4.*pi*alphao(nz)*(blackb*stemp(nz)-oray)*wght
      wsyn(nz)=wsyn(nz)+ts
      wsyn1(nz)=wsyn1(nz)+ts
  220 continue
  225 continue
c
c        add contribution from thick region
c
      tse=pi*areae*(blackb*stemp(nze)-eray)
      tso=pi*areao*(blackb*stemp(nzo)-oray)
      ptot=0.0
      do 350 nz=1,nzons
c     if(nz.lt.nze1) wsyn(nz)=wsyn(nz)+tse*wght*sdens(nz)/emass
c     if(nz.lt.nzo1) wsyn(nz)=wsyn(nz)+tso*wght*sdens(nz)/omass
      if(nz.lt.nze1) then
      wsyn(nz)=wsyn(nz)+tse*wght*sdens(nz)/emass
      wsyne(nz)=wsyne(nz)+tse*wght*sdens(nz)/emass
      endif
      if(nz.lt.nzo1) then
      wsyn(nz)=wsyn(nz)+tso*wght*sdens(nz)/omass
      wsyno(nz)=wsyno(nz)+tso*wght*sdens(nz)/omass
      endif
      ptot=ptot+wsyn(nz)*svol(nz)
  350 continue
      pwall=pwall+pi*areaw*((1.-ree-roe)*eray+(1.-roo-reo)*oray)*wght
      if(ndbug.eq.0) go to 9999
      if(mod(nstep,100).ne.0) go to 9999
         write(6,9010) nw1,nze,nzo,wfreq,wght,blackb,areae,areao,
     $   emass,omass,qee,qoo,qoe,qeo,xe,xo,srce,srco,eray,oray,
     $   pwall,ptot
 9010    format(1x,19h  synch edit....nw= ,i5//
     $   ,3x,3hnze,2x,3hnzo,5x,5hwfreq,6x,4hwght,4x,6hblackb,5x,5hareae
     $   ,5x,5hareao,
     $   /,1x,2i5,5e10.2,//6x,5hemass,5x,5homass,7x,3hqee,7x,3hqoo,
     $   7x,3hqoe,7x,3hqeo,/1x,6e10.2//9x,2hxe,8x,2hxo,6x,4hsrce,6x,
     $   4hsrco,4x,6h  eray,4x,6h  oray,/1x,6e10.2//,
     $   10x,5hpwall,11x,4hptot,/2e15.5)
         write(6,9015)
 9015    format(1x,81(1h-),/
     $   4x,'nz',6x,'sdens',7x,3hvol,'    alphae','    alphao',5x,
     $   ' wsyn',/1x,81(1h-))
         do 800 nz=1,nzons
         write(6,9020)nz,sdens(nz),svol(nz),alphae(nz),alphao(nz),
     $   wsyn(nz)
  800    continue
 9020    format(1x,i5,5e10.2)
 9999 continue
      if(istop.eq.1) go to 600
  500 continue
  600 continue
         nwsv=nw
         wstops=wfreq
c
      if(iedit.gt.0) call syedit2
      return
      end

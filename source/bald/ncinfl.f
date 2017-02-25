c--------1---------2---------3---------4---------5---------6---------7-c
c@ncinfl  .../baldur/code/bald/ncinfl.f
c  rgb 25-mar-97 zgain2/(dtmax*uest) --> zgain2
c  rgb 24-mar-97 use zflimp(ii) = interpolated value of flimp(ii,it)
c    as maximum influx rate for zeff monitor
c  rgb 02-feb-95 allow routine to run when natomc = 3
c  rgb 22-apr-94 replace zzdtis with dtmax in expression for z0
c    implement zgain4
c  rgb 21-apr-94 zgain = cimprd(2) = gain factor for feedback loop
c    make sure this routine is done only once each step
c    zgain3 = cimprd(3) = gain factor for d N / d t
c  rgb 20-apr-94 created sbrtn ncinfl from part of sbrtn imprad
c    allowed z0 to be negative
c
      subroutine ncinfl (psurfs, pvols, pinflx )
c
c   Compute the neutral impurity influx
c
c  psurfs = plasma surface area (cm**2)
c  pvols  = plasma volume (cm**3)
c  pinflx = neutral impurity influx ( particles / ( cm**2 sec ) )
c
      real psurfs, pvols, pinflx(*)
c
c     A neutral-impurity influx at the outer boundary is activated
c  by setting cfutz(200) .gt. 0.0, where as always impurity spe-
c  cies are specified by nimp(ii).  However, only impurity species
c  ii=1 and ii=2 are allowed to be influxed as neutrals,
c  and the setting of cfutz(200) affects both ii=1 and ii=2.
c
c     When cfutz(200) = 3.0, Z_eff controlled by ftzeff(it)
c  pinflx(ii) = cimprd(4) * pinflx(ii) + 
c    cimprd(2) * [ z0imp(ii) - zipsim(ii,2)
c      - cimprd(3) * ( zipsim(ii,2) - zipsum(ii,1) ) ]
c
c   ***  See documentation in subroutine imprad  ***
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'commhd.m'
c
      real ztme(2), z0imp(4), z0isum(4), zinflx(4), ziloss(4)
     &  , zlevel(2), zmxflx(4), zipsum(4,2), z0(4), zflimp(4)
cap {
      integer impnum
cap }
c
      data zinflx / 0.0, 0.0, 0.0, 0.0 /
c
      save zinflx, zipsum
c
c  ztargt    = target value for Z_eff
c  zeff      = computed value for Z_eff
c  zne       = total number of electrons
c  z0imp(ii) = the target number of impurity ions (N_target)
c  zipsum(ii,2) = the actual number of impurity ions N_new
c  zipsum(ii,1) = number of impurity ions last step  N_old
c  z0isum(ii)   = initial values of zimpsum(ii,1)
c  zinflx(ii)   = flux / cm^2 crossing plasma edge
c  ziloss(ii) is a loss rate (inverse confinement time)
c  z0 = ( N_target - N_new ) / dt  if  N_target > N_new
c  zadjst
c
c    ziloss is computed in the following way:
c  area * d N / d t = flux of ions out - flux of neutral in
c  ziloss = ( area *  (flux of neutrals in ) - d N / d t ) / N
c         = area * ( flux of ions out ) / N
c
      logical inita1, lcorec, limits, lprint, inflx2(4)
c
c  inflx2(ji) = true if there is any influx for impurity species ji
c
      data istep / 0 /
      data inita1, lcorec, limits /.true., .false., .false./
      data lprint /.false./
c
c       data  impneu,nergy1,nergy2,izlos1,izlos2 /200,201,202,203,204/,
c    1        ixlos1,ixlos2,istar1,istar2,inpth1 /205,206,207,208,209/,
c    2        inpth2,islnt1,islnt2,iajust,iprint /210,211,212,213,214/,
c    3        ifact1,ifact2 /215,216/
c
c       data  ileve1,ileve2,ixflx1,ixflx2,minti1 /220,221,222,223,224/,
c    1        ihnce1,ihnce2,maxti1,ilower,iraise /225,226,227,228,229/,
c    2        itime2,minti2,maxti2,itime3,ichng1 /260,261,262,270,271/
c
      save inita1, lcorec, limits, lprint, inflx2
cap {
      save impnum
cap }
c
      save istep
c
      save zoldts, ztme, zadjst, zunit2, ztprnt, zchnge
     &  , zshift, zmint1, zmint2, zmaxt1, zmaxt2, zf1, zf2, zf3, zf4
c
c
c  6)   Neutral impurity influx   ( only if CFUTZ(200) > epslon )
c       -----------------------
c
cbate        write (nprint,*) 'nstep= ',nstep,'  nadump(11)= ',nadump(11)
cbate     &    ,'  natomc= ',natomc,'  cfutz(200)= ',cfutz(200)
c
c..set influx to 0.0
c
      do ji=1,mimp
        pinflx(ji) = 0.0
      enddo
c
c  initially, nadump(11)=0 in order to prevent premature use of
c  neutral-impurity influx
c
        if ( nadump(11) .lt. 1 ) return
c
        if ( cfutz(200) .lt. epslon ) return
c
        if ( mimp .le. 0 ) return
c
c..treatment of neutral-impurity influx
c
        zoldts  = dtmax
        ztme(1) = tai*uist
        if ( nstep .gt. 2 ) ztme(1)=ztme(2)
c
        if ( .not. inita1 ) go to 612
c
c..initializations for neutral-impurity influx
c
        inita1=.false.
        lcorec=.false.
c
        zchnge=epsinv
        if(cfutz(270).gt.epslon) zchnge=cfutz(270)*usit
c
        write (nprint,6020)
c
        impnum = min ( mimp, 2 )
c
        call neuset(cfutz,epslon,epsinv,1)
c
        zadjst  = 1.0
cbate        zunit2  = 1./(uiel*uiel*uiet)
        zunit2  = 1.0
        ztprnt  = cfutz(214)*ueit
        nxprnt  = 0
c
        jsepx1  = mzones - 1
        if ( nadump(1) .gt. lcentr ) jsepx1 = nadump(1) - 1
c
        zshift  = epsinv
        zmint1  = 0.0
        zmint2  = 0.0
        zmaxt1  = epsinv
        zmaxt2  = epsinv
c
        zf1     = 0.0
        zf2     = 0.0
        zf3     = 0.0
        zf4     = 0.0
        if ( cfutz(224) .gt. epslon) limits=.true.
        if ( ( tbi .ge. zchnge )  .and.
     &       ( cfutz(279) .gt. epslon ) ) limits=.false.
        if ( limits ) then
          zshift = cfutz(260)*usit
          zmint1 = cfutz(224)*evs
          zmint2 = cfutz(261)*evs
          zmaxt1 = cfutz(227)*evs
          zmaxt2 = cfutz(262)*evs
          if(cfutz(224).gt.epslon) zf1=10./cfutz(224)
          if(cfutz(261).gt.epslon) zf3=10./cfutz(261)
          if(cfutz(227).gt.cfutz(224)) zf2=10./cfutz(227)
          if(cfutz(262).gt.cfutz(261)) zf4=10./cfutz(262)
        endif
c
        do ii=1,impnum
          z0imp(ii)  = 0.0
          z0isum(ii) = 0.0
          zinflx(ii) = 0.0
          ziloss(ii) = 0.0
          zlevel(ii) = 0.0
          zmxflx(ii) = 0.0
        enddo
c
        do ii=1,impnum
          inflx2(ii) = .false.
cbate          if ( cfutz(200) .gt. 2.9 ) then
cbate            do it=1,mxt
cbate              if ( ftzeff(it) .gt. 0.99 ) inflx2(ii) = .true.
cbate            enddo
cbate          else
            do it=1,mxt
              if ( flimp(ii,it) .gt. epslon ) inflx2(ii)=.true.
            enddo
cbate          endif
        enddo
c
c  If there is no impurity influx, set cfutz(200)=0.0
c    and skip the rest of this section
c
          if ( inflx2(1) ) go to 608
          if ( inflx2(2) ) go to 608
            cfutz(200) = 0.0
            return
  608   continue
c
        if ( inflx2(1) ) then
          zlevel(1)   = cfutz(220)
          zmxflx(1)   = cfutz(222)
          zipsum(1,1) = 0.0
          do jx=lcentr,ledge
            zipsum(1,1) = rhois(1,2,jx)*dx2i(jx)+zipsum(1,1)
          enddo
          zipsum(1,1) = 2.0 * pvols * zipsum(1,1)
          z0isum(1)   = zipsum(1,1)
        endif
c
        if ( inflx2(2) ) then
          zlevel(2) = cfutz(221)
          zmxflx(2) = cfutz(223)
          zipsum(2,1)=0.0
          do jx=lcentr,ledge
            zipsum(2,1) = rhois(2,2,jx)*dx2i(jx)+zipsum(2,1)
          enddo
          zipsum(2,1) = 2.0 * pvols *zipsum(2,1)
          z0isum(2)   = zipsum(2,1)
        endif
c
c  end of initialization
c
  612   continue
c
c..skip routine if nstep = istep
c
      if ( nstep .eq. istep ) then
        do ji=1,mimp
          pinflx(ji) = zinflx(ji)
        enddo
        return
      endif
c
      istep = nstep
c
c..control printout
c
        if ( tai .ge. ztprnt ) then
          lprint = .true.
          nxprnt = nxprnt+1
          if ( nxprnt .ge. 10 ) then
            lprint = .false.
            ztprnt = epsinv
          endif
        endif
c
c..time-interpolation factors
c
        it=1
        zint=0.0
        zt=.5 * uist * (tai+tbi)
        do jt=2,mxt
          if(timp(jt).le.0.0) go to 620
          it=jt
          if(timp(jt).eq.zt) go to 620
          if(timp(jt).gt.zt) go to 618
        enddo
  618   continue
c
        zint=(timp(it)-zt)/(timp(it)-timp(it-1))
c
  620   continue
c
c.."zadjst" is an adjustment factor on the neutral-impurity influx
c       control so as to approach certain prescribed limits on
c       the plasma-edge temperature.
c
        zt = zt*uist
        ztme(2) = tbi*uist
        if ( (.not.inflx2(1)) .and. (.not.inflx2(2)) ) go to 6090
        zdtime = ztme(2)-ztme(1)
        zdtime = max(zdtime,dtinit)
        zzdtis = max ( dti*uist, epslon )
        zinvdt = 1.0 / zzdtis
        zadjst = 1.0
        if ( .not. limits ) go to 6090
        if ( zt .gt. zshift ) then
          if ( tis(2,jsepx1) .gt. zmint2 ) then
            if ( tis(2,jsepx1) .lt. zmaxt2 ) go to 6090
            z2 = zf4 * max( tis(2,jsepx1)*evsinv - cfutz(262), 0.0)
            z2 = min(z2,7.0)
            zadjst = cfutz(229)
            if ( z2 .gt. 1.0 ) zadjst = cfutz(229)**(z2)
          else
            z2 = zf3 * max( cfutz(261) - tis(2,jsepx1)*evsinv, 0.0)
            z2 = min(z2,7.0)
            zadjst = cfutz(228)
            if ( z2 .gt. 1.0 ) zadjst = cfutz(228)**(z2)
          endif
        elseif ( tis(2,jsepx1) .gt. zmint1 ) then
          if ( tis(2,jsepx1) .lt. zmaxt1 ) go to 6090
          z2 = zf2 * max( tis(2,jsepx1)*evsinv - cfutz(227), 0.0)
          z2 = min(z2,7.0)
          zadjst = cfutz(229)
          if ( z2 .gt. 1.0 ) zadjst = cfutz(229)**(z2)
        else
          z2 = zf1 * max( cfutz(224) - tis(2,jsepx1)*evsinv, 0.0)
          z2 = min(z2,7.0)
          zadjst = cfutz(228)
          if ( z2 .gt. 1.0 ) zadjst = cfutz(228)**(z2)
        endif
 6090   continue
c
c..compute impurity influx and deposition profile
c
      do 636 ii = 1,impnum
c
        pinflx(ii) = 0.0
        zinflx(ii) = 0.0
c
        if ( .not. inflx2(ii) ) go to 636
c
c
c..time-interpolated neutral influx for impurity species iz
c
          zflimp(ii) = 
     &      (flimp(ii,it-1)*zint+flimp(ii,it)*(1.-zint))*zunit2
c
        if ( cfutz(200) .lt. 1.9 ) then
c
c..time-interpolated neutral influx for impurity species iz
c
          zinflx(ii) = zflimp(ii)
c
        else
c
c..here influx is determined so as to approach a prescribed
c       total number of impurity ions in the plasma
c
c..compute the actual number of impurity ions in the plasma
c
        zipsum(ii,2) = 0.0
        do jx = lcentr,ledge
          zipsum(ii,2) = rhois(ii,2,jx)*dx2i(jx) + zipsum(ii,2)
        enddo
        zipsum(ii,2) = 2.0 * pvols *zipsum(ii,2)
c
c..compute the target number of impurity ions
c
        if ( flimp(ii,it) .le. epslon ) then
c
c  revert to the original number of impurities
c
          z0imp(ii) = z0isum(ii)
c
        elseif ( cfutz(200) .gt. 2.9 ) then
c
c  compute average Z_eff and interpolate flimp to set target Z_eff
c
          zne  = 0.0
          znz2 = 0.0
          isepm1 = ledge - 1
          if ( nadump(1) .gt. lcentr ) isepm1 = nadump(1) - 1
c
          do jz=lcentr,isepm1
            zne  = zne + rhoels(2,jz) * dx2i(jz)
            znz2 = znz2 + xzeff(2,jz) * rhoels(2,jz) * dx2i(jz)
          enddo
c
          zeff = max ( znz2 / zne, 1.0 )
c
          ztime = 0.5 * ( tai + tbi )
          call timint ( ztime, ztargt, timp, mxt, ftzeff, 1, 1 )
          if ( ztargt .lt. 1.0 ) ztargt = zeff
c
          z0imp(ii) = zipsum(ii,2) * ztargt / zeff
c
        else
c
c  interpolate flimp to set target number of impurities
c
          z0imp(ii) = zadjst * zflimp(ii)
c
        endif
c
c  compute loss rate
c
      if ( zlevel(ii) .gt. epslon ) then
c
cbate        if ( lcorec ) ziloss(ii) = 
cbate     &    (psurfs*zinflx(ii)*zoldts-(zipsum(ii,2)-zipsum(ii,1)))/
cbate     &    (zoldts*zipsum(ii,2)+epslon)
c
c  compute influx
c
        zgain2 = cimprd(2)
        zgain3 = cimprd(3)
        zgain4 = cimprd(4)
c
        ziloss(ii) = - zgain3 * ( zipsum(ii,2) - zipsum(ii,1) )
c
        z0(ii) = ( z0imp(ii)-zipsum(ii,2) )
c
        zinflx(ii) = zgain4 * zinflx(ii) +
     &      zgain2 * ( z0(ii) + ziloss(ii) ) / psurfs
c
        zinflx(ii) = max ( 0.0,
     &    zlevel(ii) * min ( zflimp(ii), zinflx(ii) ) )
c
cbate        if ( .not. lcorec ) ziloss(ii) = 
cbate     &    (psurfs*zinflx(ii)*zzdtis-(zipsum(ii,2)-zipsum(ii,1)))/
cbate     &    (zzdtis*zipsum(ii,2)+epslon)
c
        endif
c
c..diagnostic printout
c
        write (nprint,6021) nstep, ii, zt, ztargt, zeff, zne
     &    , zipsum(ii,1), zipsum(ii,2), zinflx(ii), z0(ii)
     &    , ziloss(ii)
c
      endif
c
 6020 format ('#i nstep   ii',t17,'time',t29,'ztargt',t41,'zeff'
     &  ,t53,'zne',t65,'zipsum(1)',t77,'zipsum(2)',t89,'zinflx'
     &  ,t101,'z0',t113,'ziloss')
 6021 format ('#i ',2i5,1p10e12.3)
c
c  end of compuation of zinflx
c
        pinflx(ii) = zinflx(ii)
c
  636   continue
c
c
        zoldts = zzdtis
        lcorec = .true.
          do ii = 1, impnum
            zipsum(ii,1) = zipsum(ii,2)
          enddo
          ztme(1) = tbi*uist
        if ( nadump(11) .le. 1 ) then
          lcorec = .false.
        endif
c
c..diagnostic printout
c
c
        if ( lprint .and. inflx2(1) )
     &    write(nprint,6012) ziloss(1),nstep,nadump(11)
        if ( lprint .and. inflx2(2) )
     &    write(nprint,6014) ziloss(2),nstep,nadump(11)
 6012 format(/2x,
     &' loss-rate per sec as used for recycling impurity 1 ='
     & ,1pe10.3,2x,' nstep=',i4,2x,' nadump(11)=',i3/)
 6014 format(/2x,
     &' loss-rate per sec as used for recycling impurity 2 ='
     & ,1pe10.3,2x,' nstep=',i4,2x,' nadump(11)=',i3/)
 6016 format(4x,a)
 6018 format((2x,5(2x,1pe11.3,1pe11.3)))
c
      return
      end










c--------1---------2---------3---------4---------5---------6---------7-c
c@bounds.../baldur/code/bald/dsolver.f
c  ton 01-aug-01  rewrote ion and electron temperature boundary conditions
c  ajr 15-jun-95  Fixed some confusion on the temperature boundaries
c  rgb 10-aug-94  Installed a patch to skip pedestal boundary conditions
c    when natomc > 2
c  rgb 10-aug-94 save local variables that are initializ
c  rgb 17-jul-94 rewrote impurity boundary conditions
c  rgb 10:15 07-jun-94  rewrote temperature boundary conditions
c    temporarily remove uine and replace usnd -> usid
c  rgb 18.54 13-aug-90  clean up and allow use of bd... when nbound=0
c  rgb 12.52 30-jul-87  purely diffusive form of matched boundary cond.
c  rgb 21-jul-87  extensive rewrite of impurity boundary conditions.
c       Implemented section to match transport coefficients to the
c       prescribed impurity influx if LNUMER(32) > 0 and CNUMER(32) > 0.
c  dps 23-mar-87 allow for "flimp(ii,it) < 0" to pump out impurities;
c                    enabled only if cfutz(200)=0.
c      dps 08-dec-86 additional alterations as made on 14-nov but
c                    for cfutz(200) > 0 case.
c      dps 14-nov-86 improved flux at n+theta for nbound>1
c       fgps 20-sep-83 made consistent with older version in which
c                      zero boundary flux prevails when "nbound=2 or
c                      3" and "flimp(ii,it)=0.0".
cdoc
c=======================================================================
c
c       ------------
c       sbrtn BOUNDS   file DSOLVER
c       ------------
c
c       2.7     compute boundary condition coefficients
c
c  alpha1.(chi at r-(dr/2)) + (beta1 + gamma1).(chi at r+(dr/2)) = delta1
c
c  where alpha1, beta1, and gamma1 are matrices, and delta1 is a vector
c
c  sum over js  {    alpha1(jr,js) * chi(js,mzones-1)
c                 + [ beta1(jr,js) + gamma1(jr,js) ] * chi(js,mzones)  }
c                 = delta1(jr)
c
c       NOTES:
c
c  All units are internal units in this sbrtn.
c
c  CPEDST(JR)  is initialized in sbrtn START
c  NZCOMP = 0  unless sbrtn CMPRES is called
c               ( this option should not be used in the 1-1/2-D code).
c
cend
cinput
c=======================================================================
c       --------------------------------
c       Input variables for sbrtn BOUNDS
c       -------------------------------- (file DSOLVER)
c
c  NBOUND    determines the type of boundary condition used:
c  ------
c
c       0 -- edge densities and temperatures fixed
c
c       1 -- edge densities and temps. set to the minimum of the
c               fixed value of nbound=0 and less than the value of
c               the next zone in
c
c       2 -- same as 0, but impurity-ion influxes can be prescribed
c                       and/or specified by gsputs
c
c       3 -- same as 1, but impurity-ion influxes can be prescribed
c                       and/or specified by gsputs
c
c       4 -- same as 2, but influx is done as a volume source
c
c       5 -- same as 3, but influx is done as a volume source
c
c Corresponding values of internal variables:
c
c       nbound  iiboun  ihboun
c       ------  ------  ------
c         0       1       0
c         1       1       1
c         2       2       0
c         3       2       1
c         4       3       0
c         5       3       1
c
c       note:   volume sources are only implemented for impurity
c               species #1 and #2.  when species #1 and/or #2 are sub-
c               ject to neutral influxing, the setting "nbound=0 or 1"
c               assigns pedestal or corrected pedestal boundary condi-
c               tions to these species; while, "nbound=2 or 3" assigns
c               zero-influx boundary conditions.
c
c       note:   "lflxng(ii)" is a tag on impurity species ii.  if
c               species ii is subject to any kind of influxing, then
c               lflxng(ii)=.true.; otherwise, lflxng(ii)=.false..
c
c  Time-dependent boundary conditions:
c  -----------------------------------
c
c       All defaults are 0.0 which turns this option off.
c
c       These options change the value of the pedestal boundary
c               conditions at the edge of the plasma.
c       They have no effect when the pedestal boundary conditions
c               are not in use (see NBOUND).
c
c       Pedestal densities and temperatures are interpolated linearly
c               as a function of time between the time breakpoints.
c
c  BDTIME(JT)    = breakpoint times [sec]               JT=1,...,20
c
c  BDHYDE(JT,JH) = density of hydrogen isotope JH       JH=1,2
c                       [cm-3] in the ghost zone at the
c                       edge of the plasma.
c                       Implemented only if bdhyde(1,jh) > 0.
c                       Must be positive when used.
c
c  BDIMPE(JT,JI) = density of impurity species JI       JI=1,...,4
c                       [cm-3] in the ghost zone at the
c                       edge of the plasma.
c                       Implemented only when bdimpe(1,ji) > 0.
c                       Must be positive when used.
c
c  BDTEE(JT)     = electron temperature [keV] in the ghost zone
c                       at the edge of the plasma.
c                       Implemented only when bdtee(1) > 0.
c                       Must be positve when used.
c
c  BDTIE(JT)     = ion temperature [keV] in the ghost zone
c                       at the edge of the plasma.
c                       Implemented only when bdtie(1) > 0.
c                       Must be positve when used.
c
cend
c**********************************************************************c
c
        subroutine bounds
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'commhd.m'
        include 'cd3he.m'
        include 'clintf.m'
        include 'clsaw.m'
        include 'clparm.m'
        include 'clislb.m'
c
        logical   inital,lozint
c
        data      inital /.true./
	data      mode /0/
c
        save inital, ix1, zfuzz, iiboun, ihboun
c
        save mode
c
        integer      isep,     jz,     jk,      jj,
     &       jjj,    jjjj,     kk,     kkk,     pk,
     &       pkk,    pkkk,     type,   charge
c
        real zeheat,  ziheat,     zpheat,    zrmajor,  zrminor,
     &       zbtor,   zhydmass,   zne,       znebar,   pth,
     &       nped,    diff,       ratio,     zcurrent, zkappa,
     &       zdelta,  zzeff,      tped_e,    tped_i,   ratio_e,
     &       ratio_i, tped_e_old, tped_i_old,neped,    neped_mod,
     &       zpedst_old(8), nped_old, zcharge(8),      zshear,
     &       zq,      zbpol 
c
c-----------------------------------------------------------------------
cinput
c
c     new (24-feb-81):
c  if cfutz(200) .eq. 0.0, any prescribed impurity influx
c       is considered ionized; this is the default case.
c  if cfutz(200) .gt. 0.0, any prescribed impurity influx
c       of species #1 and #2 is treated as neutral impurities,
c       which ionize and radiate on penetration.  the numerical
c       simulation is worked out in subroutine imprad so as to
c       yield volume-distributed sources for impurity species #1
c       and #2.
c
cend
        data    iclass /2/,   isub /7/
c
        if ( nlomt2(isub) ) then
        call mesage(' *** 2.7 subroutine bounds bypassed')
          return
        endif
c
c
c-----------------------------------------------------------------------
c
c..temporary normalized units conversion factors
c
      uine = 1.0
      usnd = usid
c
c               initializations
c
        if ( inital ) then
          inital=.false.
          ix1=mxchi*mxchi
          zfuzz=1.0+0.001*dxboui(mzones)
          iiboun=(10*nbound+21)/20
          ihboun=nbound+2-iiboun*2
c
          if ( cfutz(200) .gt. epslon ) iiboun=3
c
c       tag the impurity species subject to any kind of influxing
c       the "abs" in the do 8 loop allows specified impurity outflow
c
          call resetl(lflxng,mximp,.false.)
          do 8 it=1,mxt
            do 8 ii=1,mximp
              if(abs(flimp(ii,it)).gt.epslon) lflxng(ii)=.true.
    8     continue
          do 12 ii=1,mximp
            if(lflxng(ii)) go to 20
   12     continue
            if(nbound.lt.2) iiboun=1
   20     continue
      endif
c
c
c      1)      pedestal boundary conditions -- no changes
c
c
cbate        if (nbound.le.0) return
cbate        if((mimp.le.0).and.(ihboun.le.0)) return
c
        call resetr(alpha1,ix1,0.0)
        call resetr(beta1 ,ix1,0.0)
        call resetr(gamma1,ix1,0.0)
        call resetr(delta1,mxchi,0.0)
c
c..calculate charge state for each species
c
        do 100 pk = 1, mhyd        
           zcharge(pk) = 1.0
 100     continue
c
        do 110 pkk = 1, mimp
           pkkk = mhyd+pkk
           zcharge(pkkk) = cmean(pkk,1,mzones)
 110     continue     
c
c..safe pedestal density for smoothing in the next time step
c
        nped_old = 0.
        jk = mhyd + mimp
c
        do 150 jj = 1,jk
            nped_old = nped_old + cpedst(jj)*zcharge(jj)
 150     continue
c
c      2)      hydrogen
c
c
  200   continue
c
      do 218 jh = 1, mhyd
c
c..option to reset cpedst(jh) as a function of time
c
          if ((bdhyde(1,jh) .gt. epslon)) then
            ztime = 0.5 * ( tai + tbi ) * uiet
            call timint (ztime,zdene,bdtime,20,bdhyde(1,jh),1,1)
            cpedst(jh) = zdene * ueid
          endif
c
          if ( ihboun .gt. 0  .and.
     &       rhohs(jh,2,mzones-1)*usid .lt. zfuzz*cpedst(jh)) then
            alpha1(jh,jh) = - 1.0
            beta1(jh,jh) = zfuzz
c
          else
c
            gamma1(jh,jh) = 1.0
            delta1(jh) = cpedst(jh)
c
          endif
c
  218  continue
c
c
c      3)      impurity
c
c
  300   continue
c
      if (mimp.le.0) go to 395
c
c
      lozint=.true.
      do 390 ji=1,mimp
        ii=ji+lhydn
        jix=ji
c
c
c       note:
c       "iiboun=1" indicates no prescribed impurity influxing of any
c                  kind as well as pedestal boundary conditions (i.e.,
c                  "nbound=0 or 1").
c       "iiboun=2" indicates that one or more impurity species are sub-
c                  ject to prescribed influxing in the ionized state.
c       "iiboun=3" indicates that impurity species #1 and/or #2 are in-
c                  fluxed as neutral atoms, and that species #3 and/or
c                  #4 may or may not be subject to prescribed influx-
c                  ing as ions.
c
c  rgb 10-aug-94  Installed a patch to skip pedestal boundary conditions
c    when natomc > 2
c
c..pedestal boundary conditions
c
        if ( ( iiboun .eq. 1  .or.  iiboun .eq. 3 )
     &    .and. natomc .lt. 3 ) then
c
c..option to reset cpedst(ii) as a function of time
c
          if (bdimpe(1,ji) .gt. epslon) then
            ztime = 0.5 * ( tai + tbi ) * uiet
            call timint (ztime,zdene,bdtime,20,bdimpe(1,ji),1,1)
            cpedst(ii) = zdene * ueid
          endif
c
c..check status of plasmas
c
           isep   = mzones
c
           zeheat = usip * ( geohms(isep) + gealfs(isep) + geauxs(isep)
     &                     + gebems(isep) + geecrs(isep) + geicrs(isep))
c
           ziheat = usip * ( gialfs(isep) + giauxs(isep) + gibems(isep)
     &                     + giecrs(isep) + giicrs(isep) )
c
c..calculate the total heating power
c
           zpheat = 1.e-6 * (zeheat + ziheat)
c
c..setup the input for L-H model 
c
           zrmajor   = max (rmajb, 0.01)
           zrminor   = max (rminb, 0.01) 
           zbtor     = max (bzs * 1.0E-4, 0.01)  
           zhydmass  = max (ahmean(1,mzones), 1.00)
c
c..calculate average electron density
c
           if (lbound(6) .eq. 2) then
              znebar = denmont * 1.e6
           else
              znebar = enlaes(mzones) * 1.e6
           endif
c
c..calculate pedestal density using models
c..determine L-H transition by using power threadhold (P > P_LH)
c
c  mode = 0  for L-mode
c       = 1  for H-mode 
c
           ztime = 0.5 * ( tai + tbi ) * uiet	  
c
           if (((cbound(2)>epslon) .and. (nstep>50)) 
     &         .and. ((ztime .ge. cbound(2)))) then
c
             if (lthery(47) == 1) then
                call bdlhmode (
     &           lbound,     cbound,     zpheat,     zrmajor,
     &           zrminor,    zbtor,      zhydmass,   znebar,
     &           pth,        mode)
             endif
c
           else   
             mode = 0
           endif
c
c..calculate pedestal density
c
        nped = 0.
        jk = mhyd + mimp
c
           if ( mode .eq. 0 ) then
             do j=1,jk
               zpedst_old(j) = cpedst(j)
             enddo
           endif
c
        do 397 jj = 1,jk
            nped = nped + cpedst(jj)*zcharge(jj) 
 397    continue
c
        neped = nped
c
c..calculate boundary density by model
c
        if ((lbound(2) .ge. 2) .and. (mode .eq. 1)) then
c 
           zcurrent = max(eqcamp*1E-6,1.0E-6)
           zbtor    = max(bzs*1.0E-4,1.0E-4)
c
           call bdhden (lbound,   cbound,    znebar,   zcurrent,
     &                  zbtor,    neped_mod)
c
        write(nout,*)	
        write(nout,*)' At time t = ', ztime, ' sec'
c
        if (lbound(3) .eq. 1) then
           write(nout,*) ' Pedestal Density Model 1'
        elseif (lbound(3) .eq. 2) then
           write(nout,*) ' Pedestal Density Model 2'
        endif
c 
        write(nout,*)' The predicted pedestal density is ', 
     &               neped_mod, ' particles/m^3.'
c
c..calculate ratio density between this time step and last time step
c
           diff = (neped_mod - nped_old)/nped_old
c
c..apply smoothing for L-H transion for hydrogenic densities
c
           zdiffmax = max ( 0.05, cbound(3) )
c
           if ( abs(diff) .gt. zdiffmax ) then 
              diff = sign(min(abs(diff), zdiffmax ), diff)
           endif   
c
c..pedestal density target from model
c
           neped_mod = nped_old + diff*nped_old
c
c..calculate ratio for each species
c
           ratio = neped_mod/neped
c
           do  kk = 1,jk
              cpedst(kk) = cpedst(kk) * ratio
              alpha1(kk,kk) = 0.0
              beta1(kk,kk)  = 0.0
              gamma1(kk,kk) = 1.0
              delta1(kk)     = cpedst(kk)
           enddo
c
           nped = 0.0
c
          do kkk = 1,jk
            nped = nped + cpedst(kkk)*zcharge(kkk)
          enddo
c
        endif   
c
          if ( ihboun .gt. 0.  and.
     &         rhois(ji,2,mzones-1)*usid .lt. zfuzz*cpedst(ii)) then
c
            alpha1(ii,ii) = - 1.0
            beta1(ii,ii) = zfuzz
c
          else
c
            gamma1(ii,ii) = 1.0
            delta1(ii) = cpedst(ii)
c
          endif
c
        else
c
c..ionized-impurity influx
c
cbate        if ( ( iiboun .eq. 2 .or. iiboun .eq. 3 )
cbate     &       .and. ( lflxng(ji) .or. nbound.gt.1 ) ) then
c
c..zdelta = influx of impurity ii at time tai [ particles / ( m**2 sec ) ]
c
c
            zdelta = gsputs(ji)*uisl*uisl*uist
c
          if ( iiboun .eq. 2 .or.
     &         cfutz(200) .le. epslon .or. jix .gt. 2 ) then
c
c..interpolate prescribed impurity-ion influx as a function of time
c
        if ( lozint ) then
          lozint=.false.
          it = 1
          zint = 0.0
          zt = 0.5*(tai + tbi)
c
          do 320 jt = 2, mxt
            if (gtflwi(jt).le.0.0) go to 324
            it = jt
            if (gtflwi(jt).eq.zt) go to 324
            if (gtflwi(jt).gt.zt) go to 322
  320     continue
  322     continue
c
          zint = (gtflwi(it) - zt) / (gtflwi(it) - gtflwi(it-1))
c
        endif
c
  324   continue
c
            zdelta = gsputs(ji) * uisl*uisl * uist +
     &          gflowi(ii,it-1)*zint + gflowi(ii,it)*(1.0 - zint)
c
          endif
c
cinput
c
c..Match transport coefficients to the prescribed impurity influx zdelta.
c  Add equal amounts to aaaa(ii,ii,jz) and bbbb(ii,ii,jz)
c  to match (or partially match) the inpurity influx at the boundary.
c  Smoothly match this flux to the transport driven flux in the
c  interior of the plasma.
c
c  LNUMER(32) = number of zone boundaries over which aaaa and bbbb
c               are adjusted
c               must be positive to implement this procedure
c             = 0 for no adjustment (default)
c             = 1 to adjust aaaa and bbbb at jz=mzones only
c             > 0 to adjust aaaa and bbbb from jz=mzones-lnumer(32)
c                                         out to jz=mzones
c             = 3 recommended value
c
c  CNUMER(32) = fraction by which flux is adjusted to match influx
c             = 0. no adjustment (default)
c             = 1. full adjustment (recommended value)
c
c  CNUMER(31) = measure of how flat chi(ii,jz) profile is prescribed
c               across plasma boundary before new values of chi(ii,jz)
c               are computed.
c             =  sgn(1.,zfluxi) * ( chi(ii,mzones) - chi(ii,mzones-1) )
c                               / ( chi(ii,mzones) + chi(ii,mzones-1) )
c               inversely proportional to imposed diffusion coef
c                               at plasma boundary
c             = 0.05 (recommended value)
c             = 0.0  (default)
c
c               IF ( ABS(CNUMER(31)) .LT. EPSLON ) CNUMER(31) = 0.05
cend
c
c
        if ( versno .gt. 12.51  .and.  lnumer(32) .gt. 0  .and.
     &       abs( cnumer(32) ) .gt. epslon ) then
c
c  add up the prescribed influx minus the sum of the off-diagonal fluxes
c
        zfluxi = zdelta
c
        do 350 js=1,mchi
          if ( js .ne. ii ) zfluxi = zfluxi
     &          - aaaa(ii,js,mzones) * chi(js,mzones-1)
     &          - bbbb(ii,js,mzones) * chi(js,mzones)
 350    continue
c
c  compute zaddi and zbddi = contributions to be added to aaaa and bbbb
c
        if ( abs(cnumer(31)) .lt. epslon ) then
          z31 = cnumer(31)
          cnumer(31) = 0.05
          write (nout,*) '   cnumer(31) changed from'
     &          ,z31,' to ',cnumer(31)
        endif
c
        z31max = 0.80
        if ( abs(cnumer(31)) .gt. z31max ) then
          z31 = cnumer(31)
          cnumer(31) = sign ( z31max, cnumer(31) )
          write (nout,*) '   cnumer(31) changed from'
     &          ,z31,' to ',cnumer(31)
        endif
c
        zd = sign (1., zfluxi) * cnumer(31)
        chi(ii,mzones) = chi(ii,mzones-1) * (1.+zd) / (1.-zd)
c
        zaddi = - aaaa(ii,ii,mzones)   - abs(zfluxi)
     &        / ( cnumer(31) * ( chi(ii,mzones) + chi(ii,mzones-1) ) )
c
        zbddi = ( zfluxi
     &   - ( aaaa(ii,ii,mzones) + zaddi ) * chi(ii,mzones-1)
     &          - bbbb(ii,ii,mzones) * chi(ii,mzones) )
     &          / ( chi(ii,mzones) )
c
c  match this contribution to the transport in the interior of the plasma
c  zfract = cnumer(32)  at jz = mzones
c         = 0.0         at jz = mzones - lnumer(32) - 1
c
        do 354 jz=mzones-lnumer(32), mzones
c
          zfract = cnumer(32) * real( jz - (mzones-lnumer(32)-1) )
     &                        / real( lnumer(32) + 1 )
c
          aaaa(ii,ii,jz) = aaaa(ii,ii,jz) + zaddi * zfract
          bbbb(ii,ii,jz) = bbbb(ii,ii,jz) + zbddi * zfract
c
 354    continue
c
      endif
c
c
c..set up delta1 so that flux at t.s. n+theta is specified, put
c       chi at t.s. n terms into delta1
c
        z1=0.
        do 386 js = 1, mchi
          alpha1(ii,js) = aaaa(ii,js,mzones)
          beta1(ii,js)  = bbbb(ii,js,mzones)
          z1=z1 + alpha1(ii,js)*chi(js,ledge)
     &          + beta1(ii,js) *chi(js,mzones)
  386   continue
          delta1(ii) = ( zdelta - (1.0-thetai)*z1) / thetai
c
c
        endif
c
c..end of loop over impurity species ji
c
 390  continue
c
c..check status of plasma
c
 395  continue
c
c      4)      electron temperature
c
c
 400  continue
c
c..save tped for smoothing in the next time step
c
        tped_e_old = cpedst(lelec) / ueih
c
c..option to reset cpedst(lelec) as a function of time
c
        if (bdtee(1) .gt. epslon) then
c
             ztime = 0.5 * ( tai + tbi ) * uiet
             call timint (ztime,ztemp,bdtime,20,bdtee,1,1)
             cpedst(lelec) = ztemp * ueih
c
        endif  
c
        if ((lbound(4) .ge. 2) .and. (mode .eq. 1)) then
c
c..input for t_ped calculation
c
          zrminor  = max (rminb, 0.01)   
          zrmajor  = max (rmajb, 0.01)
          zkappa   = max (elongb, 1.00)     
          zdelta   = tringb
          zcurrent = max (eqcamp * 1.E-6, 1.E-6)
          zbtor    = max (bzs * 1.E-4, 1.E-4)
          zhydmass = max (ahmean(1,mzones), 1.00)
          zzeff    = max (xzeff(1,mzones), 1.00)
          zshear   = max (shear(mzones), 0.01)
          zq       = max (q(mzones), 0.01)
          zbpol    = max (bpoli(mzones) * 1.E-4, 0.01)
c
          zdiffmin = max ( 0.001, cbound(4) )
c
          type = 1
	  charge = 1 	
c
          call bdhtemp (
     &       lbound,       cbound,     zrminor,      zrmajor,
     &       zkappa,       zdelta,     zcurrent,     zbtor,
     &       nped,         zhydmass,   zzeff,        zshear,
     &       zq,           zbpol,      znebar,       zpheat,
     &       tped_e)
c
          write(nout,*)	
          write(nout,*)' At time t = ', ztime, ' sec'
c
        if (lbound(5) .eq. 1) then
           write(nout,*) ' Pedestal Temperature Model based on model 1'
        elseif (lbound(5) .eq. 2) then
           write(nout,*) ' Pedestal Temperature Model based on model 2'
        elseif (lbound(5) .eq. 3) then
           write(nout,*) ' Pedestal Temperature Model based on model 3'
        elseif (lbound(5) .eq. 4) then
           write(nout,*) ' Pedestal Temperature Model based on model 4'
        elseif (lbound(5) .eq. 5) then
           write(nout,*) ' Pedestal Temperature Model based on model 5'
        elseif (lbound(5) .eq. 6) then
           write(nout,*) ' Pedestal Temperature Model based on model 6'
        elseif (lbound(5) .eq. 11) then
           write(nout,*) ' Pedestal Temperature Model based on model 11'
        elseif (lbound(5) .eq. 12) then
           write(nout,*) ' Pedestal Temperature Model based on model 12'
        endif
c
          write(nout,*)' The electron pedestal temperature is ', 
     &               tped_e, ' keV.'
c
c..apply smoothing for electron temperature
c
          diff = tped_e - tped_e_old
c
          if ( abs(diff) .gt. zdiffmin ) then 
             diff = sign(min(abs(diff), zdiffmin ), diff)
          endif   
c
          tped_e = tped_e_old + diff
c
          cpedst(lelec) = tped_e * ueih
c
        endif   
c
      if ( ihboun .gt. 0  .and.
     &    tes(2,mzones-1)*usih .le. cpedst(lelec)*zfuzz ) then
c
        alpha1(lelec,lelec) = - 1.0 / (rhoels(2,mzones-1) * usnd)
        beta1(lelec,lelec) = 1.0 / (rhoels(2,mzones) * usnd) * zfuzz
c
      else
c
        do 404 js = lhyd1, limpn
          gamma1(lelec,js) = cpedst(lelec)*azzz(js,mzones) * uine
  404   continue
c
        gamma1(lelec,lelec) = - 1.0
c
        endif
c
c
c      5)      ion temperature
c
c
  500   continue
c
c..save tped for smoothing in the next time step
c
        tped_i_old = cpedst(lion) / ueih
c
c..option to reset cpedst(lion) as a function of time
c
        if (bdtie(1) .gt. epslon) then
c
             ztime = 0.5 * ( tai + tbi ) * uiet
             call timint (ztime,ztemp,bdtime,20,bdtie,1,1)
             cpedst(lion) = ztemp * ueih
c
        endif
c
        if ((lbound(4) .ge. 2) .and. (mode .eq. 1)) then
c
          type = 1
          charge = 2	
c 
          call bdhtemp (
     &       lbound,       cbound,     zrminor,      zrmajor,
     &       zkappa,       zdelta,     zcurrent,     zbtor,
     &       nped,         zhydmass,   zzeff,        zshear,
     &       zq,           zbpol,      znebar,       zpheat,
     &       tped_i)
c
          write(nout,*)	
          write(nout,*)' At time t = ', ztime, ' sec'
          write(nout,*)' The ion pedestal temperature is ', 
     &               tped_i, ' keV.'
c
c..apply smoothing for L-H transion for impurity densities
c
          diff = tped_i - tped_i_old
c
          if (abs(diff) .gt. zdiffmin ) then 
            diff = sign(min(abs(diff), zdiffmin ), diff)
          endif   
c
          tped_i = tped_i_old + diff
c
          cpedst(lion) = tped_i * ueih 	
c
        endif   
c
      if (ihboun.gt.0.and.
     &      tis(2,mzones-1)*usih .le. cpedst(lion)*zfuzz ) then
c
          alpha1(lion,lion) = - 1.0 / (rhoins(2,mzones-1) * usnd)
          beta1(lion,lion) = 1.0 / (rhoins(2,mzones) * usnd) * zfuzz
c
      else
c
        do 504 js = lhyd1, limpn
          gamma1(lion,js) = cpedst(lion) * uine
  504   continue
c
        gamma1(lion,lion) = - 1.0
c
      endif
c
      return
      end

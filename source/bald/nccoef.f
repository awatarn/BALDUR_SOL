c--------1---------2---------3---------4---------5---------6---------7-c
c@nccoef  /11040/baldur/code/bald/dimprad.f
c  rgb 31-jan-95 equilibrium ionization states used for boundary cond
c    removed renormalization of 19-dec-94
c  rgb 19-dec-94 boundary condition renormalized after each step
c  rgb 14:40 11-nov-93 replaced call pdx(scroff,6,52) with call pdx
c  rgb 20.18 31-oct-91 allow time-varying impurity boundary condition
c    cimpde(1,ji) is reset by interpolating bdimpe wrt time
c       dps 03-nov-89 15.15 fix "influx instability" by converting flout
c                     and flscr from particle losses to loss rates.
c       dps 26-oct-89 15.14 fix bug in zdelp1 arising with nonuniform grid.
c       dps 16-aug-89 15.13 iterate on main loop during corrector step
c                     to improve accuracy of recycling.
c       dps 07-jul-89 15.11 add neutral impurities to equations
c       dps 15-may-89 15.09 remove equilibrium time-centering factors.
c       dps 10-nov-88 15.07 insert beta factors in bbb (bug fix).
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 15-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 09/08/84: sign of flout, flscr, pflx, ploss to "+": outflux
c       rhw 03/08/84: new power-law-scheme, for boundary condition too.
c       rhw 25/07/84: change theta to thetr at boundary condition.
c       rhw 05/06/84: insert beta, calculate flout, flscr.
c       rhw 19/03/84: integration sceme for large drift velocities
c       rhw 14/03/84: change time step control
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine nccoef(ji)
c
c       2.20.4   coefficients for non-corona radiation solution
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'commhd.m'
      include 'comncr.m'
c
      real zk(kzmax), zn(kzmax)
!cap
      dimension zscrot (mj),
     1          zxna(mj,0:kzmax,kncimp), zdeln(0:kzmax)
c
      data  thetr  /1.0/
      data  iheflx, icfout, idvout, impneu /70, 71, 72, 200/
c
c    calculate scrape-off loss rates
c
      call resetr (scroff,312,0.0)
      if (nadump(1).gt.lcentr)   call pdx
c (scroff,6,52)
c
c  15.13 Reset counter for iterations on corrector step
c
         ipc = 0
c
         jimp = ji + lhydn
         zdelti = 1.0/delts
         nk1 = nk - 1
c
c  Zero flux central boundary conditions
c
         aaa(1) = 1.0
         bbb(1) = 1.0
         ccc(1) = 0.0
         ddd(1) = 0.0
c
c  aaa = 0. for all allowed outer boundary conditions
c
         aaa(mzones) = 0.0
c
      if (ji.le.1)   xnerr = 0.0
c
      do 2 j=lcentr,ledge
         zscrot(j) = scroff(jimp,j)
    2 continue
c
c  Interpolate ionized impurity influx
c  Note that zdelta = outwardly directed flux at boundary
c
      zfuzz = 1. + 0.001 * dxboui(mzones)
      if (nbound.gt.1) then
        zdelta = - gsputs(ji)
        if (cfutz(impneu).le.epslon) then
          it = 1
          zint=0.0
          zt = 0.5 * (tai + tbi)
c
          do 3 jt=2,mxt
            if (gtflwi(jt).le.0.0) go to 5
            it = jt
            if (gtflwi(jt).eq.zt) go to 5
            if (gtflwi(jt).gt.zt) go to 4
    3     continue
c
    4     continue
          zint = (gtflwi(it) - zt) / (gtflwi(it) - gtflwi(it-1))
c
    5     continue
          zdelta = - gsputs(ji) - (gflowi(jimp,it-1) * zint
     1             + gflowi(jimp,it) * (1.0 - zint))/(uisl**2*uist)
c
        end if
      end if
c
    7 continue
c
c
c..option to reset cimped(1,ji) as a function of time
c     Note: cimped is in the same units as bdimpe
c
      if ( bdimpe(1,ji) .gt. epslon ) then
        ztime = 0.5 * ( tai + tbi ) * uiet
        call timint (ztime,zdene,bdtime,20,bdimpe(1,ji),1,1)
c
c..Solve the equilibrium ionization rate equations
c  assuming steady state with no transport, sources, or sinks
c  -S(k) n(k) + S(k-1) n(k-1) - R(k) n(k) + R(k+1) n(k+1) = 0
c  where n(k) = density of ionization state k
c  S(k) = ionization rate
c  R(k) = recombination rate
c
c  Let n(k) = zk(k) n(k-1).
c  Compute zk(k) = S(k-1) / ( S(k) + R(k) - zk(k+1) R(k+1) )
c
        zepsln = 1.e-6
c
        nk = nkimp(ji)
c
        jz = mzones
c
        zk(nk) = sa(jz,nk-1,ji) / ( sa(jz,nk,ji) + ra(jz,nk,ji) )
        do jk=nk-1,1,-1
          zk(jk) = sa(jz,jk-1,ji) /
     &      max ( zepsln * sa(jz,jk-1,ji), 
     &      (sa(jz,jk,ji) + ra(jz,jk,ji) - zk(jk+1) * ra(jz,jk+1,ji)) )
        enddo
c
c..find the index where zk(jk) goes from > 1 to < 1
c
c
        iz = 1
        do jk=1,nk
          if ( zk(jk) .gt. 1.0 ) iz = jk
        enddo
c
c..compute density factors
c
        zn(iz) = 1.0
        if ( iz .lt. nk ) then
          do jk=iz+1,nk
            zn(jk) = zn(jk-1) * zk(jk)
            if ( zn(jk) .lt. zepsln ) zn(jk) = 0.0
          enddo
        endif
c
        if ( iz .gt. 1 ) then
          do jk=iz-1,1,-1
            zn(jk) = zn(jk+1) / zk(jk+1)
            if ( zn(jk) .lt. zepsln ) zn(jk) = 0.0
          enddo
        endif
c
        zsum = 0.0
        do jk=1,nk
          zsum = zsum + zn(jk)
        enddo
c
c..normalize the density of each ionization state so that the sum
c  of these densities equals the total impurity density
c
        do jk=1,nk
          cimped(jk,ji) = zdene * zn(jk) / zsum
        enddo
c
c
      endif
c
c..15.13 Store fluxes to test at end of each recycling-corrector iteration;
c        reset fluxes and local (this ji) error criterion
c
      zflout = flout(ji)
      zflscr = flscr(ji)
      flout(ji) = 0.0
      flscr(ji) = 0.0
      zxnerr = 0.0
c
      do 99 nc=1,ncrept
c
      do 98 k=1,2
c
c    first calculate "downwards in ik" = recombination implicit (k=1)
c    second calculate "upwards in ik" = ionisation implicit (k=2)
c    the first step is only an "ionisation"-step
c..15.11 Note: neutrals do not have scrape-off losses; source of +1 due
c        to impurity influx now appears only in central do-loop.
c
      if (tai*uist.le.tinit.and.k.eq.1)   go to 98
c
c  Define time differentials for use in extrapolating equilibrium
c  quantities. They are given at 0.5 * (tai + tbi).
c  The V' factors on the left side of the
c  equation are evaluated at the same point in time as the densities.
c  On the right side, everyting is at the centered time, 0.5 * (tai + tbi).
c
      zdtot = (delts * ncrept / fcstep) * usit
      zdtn = - 0.5 * zdtot + ((float(nc) - 1.) / fcstep
     1       + (float(k) - 1.) - (2. - 1./fcstep)) * delts * usit
      zdtnp1 = zdtn + delts * usit
c
c     save old densities
c
      do 8 ik=0,nk
      do 8 j=1,mzones
         zxna(j,ik,ji) = xn(j,ik,ji)
    8 continue
c
c     calculate final coefficients and solve rate equations
c
c      stages  "nk" (k=1) and "0" (k=2) resp.
c
         ik = (2-k)*nk
         ikm1 = nk1
         ikp1 = 1
c
      do 10 j=lcentr,ledge
         jj = j + 1
         jm = j - 1
c
         zdxi = avi(j,5,1) / uisl
         zdxip = avi(jj,5,1) / uisl
         zdr = avi(j,6,1) / (uisl**2 * zdxi * dxboui(j))
         zdrp1 = avi(jj,6,1) / (uisl**2 * zdxip * dxboui(jj))
         zvmag = avi(j,1,2) / (uist * avi(j,2,1) * zdxi)
         zvmagp = avi(jj,1,2) / (uist * avi(jj,2,1) * zdxip)
         zvpnew = avi(j,4,1) + zdtnp1 * avi(j,4,2)
         zdel = avi(j,3,1) * zdxi / (dxzoni(j) * zvpnew)
         zdelp1 = avi(jj,3,1) * zdxip / (dxzoni(j) * zvpnew)
         zvdvp1 = avi(j,4,1) / zvpnew
         zvtrat = (avi(j,4,1) + zdtn * avi(j,4,2)) / (zvpnew * delts)
c
         zxd1 = da(j,ik,ji) * zdr
         zxd2 = da(jj,ik,ji) * zdrp1
         zp1 = (va(j,ik,ji) - zvmag) / zxd1
         zp2 = (va(jj,ik,ji) - zvmagp) / zxd2
         zq1 = (1.0 - 0.1*abs(zp1))**5
         zq2 = (1.0 - 0.1*abs(zp2))**5
         zf1 = max (0.0,-zp1) + max (0.0,zq1)
         zf2 = max (0.0,-zp2) + max (0.0,zq2)
         zf11 = max (0.0,zp1) + max (0.0,zq1)
         zf22 = max (0.0,zp2) + max (0.0,zq2)
         za = zdelp1 * zxd2 *zf2
         zc = zdel * zxd1 * zf11
         zb = zdelp1 * zxd2 * zf22 + zdel * zxd1 * zf1
     1         - (2-k) * zvdvp1 * zscrot(j)
         aaa(j) = theta * za
         ccc(j) = theta * zc
         bbb(j) = zdelti + theta * zb
     1          + zvdvp1 * ((k-1)*sa(j,ik,ji) + (2-k)*ra(j,ik,ji))
     2                      * beta * rhoels(2,j)
         ddd(j) =  (1.-theta) * za * xn(jj,ik,ji)
     1         + (1.-theta) * zc * xn(jm,ik,ji)
     2         - ( zvdvp1*((k-1)*ra(j,ik,ji) + (2-k)*sa(j,ik,ji))
     3                      * beta * rhoels(2,j)
     4                  + (1.-theta)*zb - zvtrat ) * xn(j,ik,ji)
     5         + zvdvp1*( (k-1)*ra(j,ikp1,ji)
     6                      * beta * rhoels(2,j) * xn(j,ikp1,ji)
     7                  + (2-k)*sa(j,ikm1,ji)
     8                      * beta * rhoels(2,j) * xn(j,ikm1,ji)
     9                  + (2-k)*dz(j,ji) )
   10 continue
c
c  Set outer boundary condition. Use variables defined in last pass
c  through above loop to obtain quantities evaluated at mzones.
c
      if (nbound.gt.1) then
        bbb(mzones) = theta * zxd2 * zf2
        ccc(mzones) = theta * zxd2 * zf22
        zdd1 = (1.-theta) * zxd2 * zf2
        zdd2 = (1.-theta) * zxd2 * zf22
        ddd(mzones) = - zdd1 * xn(mzones,ik,ji)
     1                + zdd2 * xn(ledge,ik,ji)
      else if ((nbound.eq.1) .and.
     1           (zxna(ledge,ik,ji).lt.zfuzz*cimped(ik,ji))) then
        bbb(mzones) = zfuzz
        ccc(mzones) = 1.0
        ddd(mzones) = 0.0
      else
        bbb(mzones) = 1.0
        ccc(mzones) = 0.0
        ddd(mzones) = cimped(ik,ji)
      end if
c
      call ncsolv (ik,ji)
c
c      stage  "nk-1" to "1" (k=1)  and  "1" to "nk-1" (k=2)  resp.
c
      do 20 iik=1,nk1
         ik = (k-1)*iik + (2-k)*(nk-iik)
         ikm1 = ik - 1
         ikp1 = ik + 1
c
      do 22 j=lcentr,ledge
         jj = j + 1
         jm = j - 1
c
         zdxi = avi(j,5,1) / uisl
         zdxip = avi(jj,5,1) / uisl
         zdr = avi(j,6,1) / (uisl**2 * zdxi * dxboui(j))
         zdrp1 = avi(jj,6,1) / (uisl**2 * zdxip * dxboui(jj))
         zvmag = avi(j,1,2) / (uist * avi(j,2,1) * zdxi)
         zvmagp = avi(jj,1,2) / (uist * avi(jj,2,1) * zdxip)
         zvpnew = avi(j,4,1) + zdtnp1 * avi(j,4,2)
         zdel = avi(j,3,1) * zdxi / (dxzoni(j) * zvpnew)
         zdelp1 = avi(jj,3,1) * zdxip / (dxzoni(j) * zvpnew)
         zvdvp1 = avi(j,4,1) / zvpnew
         zvtrat = (avi(j,4,1) + zdtn * avi(j,4,2)) / (zvpnew * delts)
c
         zxd1 = da(j,ik,ji) * zdr
         zxd2 = da(jj,ik,ji) * zdrp1
         zp1 = (va(j,ik,ji) - zvmag) / zxd1
         zp2 = (va(jj,ik,ji) - zvmagp) / zxd2
         zq1 = (1.0 - 0.1*abs(zp1))**5
         zq2 = (1.0 - 0.1*abs(zp2))**5
         zf1 = max (0.0,-zp1) + max (0.0,zq1)
         zf2 = max (0.0,-zp2) + max (0.0,zq2)
         zf11 = max (0.0,zp1) + max (0.0,zq1)
         zf22 = max (0.0,zp2) + max (0.0,zq2)
         za = zdelp1 * zxd2 *zf2
         zc = zdel * zxd1 * zf11
         zb = zdelp1 * zxd2 * zf22 + zdel * zxd1 * zf1
     1         - zvdvp1 * zscrot(j)
         zdd0 = 0.0
         if (ik.eq.1) zdd0 = dq(j,ji)
         aaa(j) = theta * za
         ccc(j) = theta * zc
         bbb(j) = zdelti + theta * zb
     1          + zvdvp1 * ((k-1)*sa(j,ik,ji) + (2-k)*ra(j,ik,ji))
     2                      * beta * rhoels(2,j)
         ddd(j) =  (1.-theta) * za * xn(jj,ik,ji)
     1         + (1.-theta) * zc * xn(jm,ik,ji)
     2         - ( zvdvp1*((k-1)*ra(j,ik,ji) + (2-k)*sa(j,ik,ji))
     3                      * beta * rhoels(2,j)
     4                  + (1.-theta)*zb - zvtrat ) * xn(j,ik,ji)
     5         + zvdvp1*( ra(j,ikp1,ji)
     6                      * beta * rhoels(2,j)*xn(j,ikp1,ji)
     7                  + sa(j,ikm1,ji)
     8                      * beta * rhoels(2,j)*xn(j,ikm1,ji) + zdd0)
   22 continue
c
      if (nbound.gt.1) then
        bbb(mzones) = theta * zxd2 * zf2
        ccc(mzones) = theta * zxd2 * zf22
        zdd0 = 0.0
        if (ik.eq.1) zdd0 = zdelta
        zdd1 = (1.-theta) * zxd2 * zf2
        zdd2 = (1.-theta) * zxd2 * zf22
        ddd(mzones) = - zdd1 * xn(mzones,ik,ji)
     1                + zdd2 * xn(ledge,ik,ji) - zdd0
      else if ((nbound.eq.1) .and.
     1           (zxna(ledge,ik,ji).lt.zfuzz*cimped(ik,ji))) then
        bbb(mzones) = zfuzz
        ccc(mzones) = 1.0
        ddd(mzones) = 0.0
      else
        bbb(mzones) = 1.0
        ccc(mzones) = 0.0
        ddd(mzones) = cimped(ik,ji)
      end if
c
      call ncsolv (ik,ji)
   20 continue
c
c      stages  "0" (k=1) and "nk" (k=2) resp.
c
         ik = (k-1)*nk
         ikm1 = nk1
         ikp1 = 1
c
      do 30 j=lcentr,ledge
         jj = j + 1
         jm = j - 1
c
         zdxi = avi(j,5,1) / uisl
         zdxip = avi(jj,5,1) / uisl
         zdr = avi(j,6,1) / (uisl**2 * zdxi * dxboui(j))
         zdrp1 = avi(jj,6,1) / (uisl**2 * zdxip * dxboui(jj))
         zvmag = avi(j,1,2) / (uist * avi(j,2,1) * zdxi)
         zvmagp = avi(jj,1,2) / (uist * avi(jj,2,1) * zdxip)
         zvpnew = avi(j,4,1) + zdtnp1 * avi(j,4,2)
         zdel = avi(j,3,1) * zdxi / (dxzoni(j) * zvpnew)
         zdelp1 = avi(jj,3,1) * zdxip / (dxzoni(j) * zvpnew)
         zvdvp1 = avi(j,4,1) / zvpnew
         zvtrat = (avi(j,4,1) + zdtn * avi(j,4,2)) / (zvpnew * delts)
c
         zxd1 = da(j,ik,ji) * zdr
         zxd2 = da(jj,ik,ji) * zdrp1
         zp1 = (va(j,ik,ji) - zvmag) / zxd1
         zp2 = (va(jj,ik,ji) - zvmagp) / zxd2
         zq1 = (1.0 - 0.1*abs(zp1))**5
         zq2 = (1.0 - 0.1*abs(zp2))**5
         zf1 = max (0.0,-zp1) + max (0.0,zq1)
         zf2 = max (0.0,-zp2) + max (0.0,zq2)
         zf11 = max (0.0,zp1) + max (0.0,zq1)
         zf22 = max (0.0,zp2) + max (0.0,zq2)
         za = zdelp1 * zxd2 *zf2
         zc = zdel * zxd1 * zf11
         zb = zdelp1 * zxd2 * zf22 + zdel * zxd1 * zf1
     1         - (k-1) * zvdvp1 * zscrot(j)
         aaa(j) = theta * za
         ccc(j) = theta * zc
         bbb(j) = zdelti + theta * zb
     1          + zvdvp1 * ((k-1)*sa(j,ik,ji) + (2-k)*ra(j,ik,ji))
     2                      * beta * rhoels(2,j)
         ddd(j) =  (1.-theta) * za * xn(jj,ik,ji)
     1         + (1.-theta) * zc * xn(jm,ik,ji)
     2         - ( zvdvp1*((k-1)*ra(j,ik,ji) + (2-k)*sa(j,ik,ji))
     3                      * beta * rhoels(2,j)
     4                  + (1.-theta)*zb - zvtrat ) * xn(j,ik,ji)
     5         + zvdvp1*( (2-k)*ra(j,ikp1,ji)
     6                      * beta * rhoels(2,j) * xn(j,ikp1,ji)
     7                  + (k-1)*sa(j,ikm1,ji)
     8                      * beta * rhoels(2,j) * xn(j,ikm1,ji)
     9                  + (k-1)*dz(j,ji) )
   30 continue
c
      if (nbound.gt.1) then
        bbb(mzones) = theta * zxd2 * zf2
        ccc(mzones) = theta * zxd2 * zf22
        zdd1 = (1.-theta) * zxd2 * zf2
        zdd2 = (1.-theta) * zxd2 * zf22
        ddd(mzones) = - zdd1 * xn(mzones,ik,ji)
     1                + zdd2 * xn(ledge,ik,ji)
      else if ((nbound.eq.1) .and.
     1           (zxna(ledge,ik,ji).lt.zfuzz*cimped(ik,ji))) then
        bbb(mzones) = zfuzz
        ccc(mzones) = 1.0
        ddd(mzones) = 0.0
      else
        bbb(mzones) = 1.0
        ccc(mzones) = 0.0
        ddd(mzones) = cimped(ik,ji)
      end if
c
      call ncsolv (ik,ji)
c
c..renormalize boundary densities
c
c  The sum of the ionization states should equal the prescribed impurity
c  density at the boundary.
c
cbate      if ( nbound .lt. 2 ) then
c
cbate        zsum = 0.0
cbate        do jk=1,nk
cbate          zsum = zsum + xn(mzones,jk,ji)
cbate        enddo
c
cbate        if ( zsum .gt. epslon  .and.  cimped(1,ji) .gt. epslon ) then
cbate          do jk=1,nk
cbate            xn(mzones,jk,ji) = xn(mzones,jk,ji) * cimped(1,ji) / zsum
cbate          enddo
cbate        endif
c
cbate      endif
c
c     save old values of xns
c
      do 41 j=1,mzones
         xnas(j) = xns(j,ji)
         xns(j,ji) = 0.0
   41 continue
c
c..15.11 Make xns = total impurity density, but use rhois, chi,
c        etc. for impurity ion density.
c
      do 44 ik=0,nk
      do 44 j=lcentr,mzones
         xns(j,ji) = xns(j,ji) + xn(j,ik,ji)
   44 continue
c
c    calculate accuracy criterium
c
      do 48 j=lcentr,ledge
         imax = ismax(nk+1,xn(j,0,ji),mxzone)-1
         zxsmax = xn(j,imax,ji)
      if (zxsmax.le.epslon)   go to 48
c
      do 47 ik=0,nk
         zdeln(ik) = abs(xn(j,ik,ji)-zxna(j,ik,ji)) / (zxsmax*delmax)
   47 continue
c
         imax = ismax(nk+1,zdeln,1)-1
         zxnerr = max(zxnerr,zdeln(imax))
   48 continue
c
c    calculate conservation check variables
c
c       surface outflux: zflx
c
c  Again, we make use of some of the variables calculated in
c  the last pass through the above loop.
c
      zflx = 0.0
c
      do 50 ik=0,nk
        zxd2 = da(mzones,ik,ji) * zdrp1
        zp2 = (va(mzones,ik,ji) - zvmagp) / zxd2
        zq2 = (1.0 - 0.1 * abs(zp2))**5
        zf2 = max(0.0,-zp2) + max(0.0,zq2)
        zf22 = max(0.0,zp2) + max(0.0,zq2)
        znbi = theta*xn(ledge,ik,ji) + (1.-theta)*zxna(ledge,ik,ji)
        znba = theta*xn(mzones,ik,ji) + (1.-theta)*zxna(mzones,ik,ji)
        zflx = zflx + znbi*zxd2*zf22 - znba*zxd2*zf2
   50 continue
c
         zsurfs = avi(mzones,3,1) * zdxip * uisl**3
c
c  Use area and time-integrated diffusive losses for computing recycling
c..15.15 To avoid influx instability, store flout and flscr as rates.
c
         if (versno.gt.15.14) then
           flout(ji) = flout(ji) + zflx*zsurfs/(ncrept/fcstep)
         else
           flout(ji) = flout(ji) + zflx*delts*zsurfs
         end if
         zflx = zsurfs * zflx
c
c       volume-losses by scrape-off: zflx1
c..15.11 Neutrals not affected by scrape-off => subtract out
c       volume-sources: zflx2
c
         zflx1 = 0.0
         zflx2 = 0.0
c
      do 55 j=lcentr,ledge
         zdvs = avi(j,4,1) *dxzoni(j)*uisl**3
         znbi = theta * (xns(j,ji) - xn(j,0,ji))
     1          + (1.0-theta) * (xnas(j) - zxna(j,0,ji))
         zflx1 = zflx1 - scroff(jimp,j)*znbi*zdvs
         zflx2 = zflx2 + (dq(j,ji) + dz(j,ji))*zdvs
   55 continue
c
c  Use volume and time-integrated scrape-off losses for computing recycling
c
         if (versno.gt.15.14) then
           flscr(ji) = flscr(ji) + zflx1/(ncrept/fcstep)
         else
           flscr(ji) = flscr(ji) + zflx1*delts
         end if
         zloss = zflx1
         zsorc = zflx2
c
c       time integrated fluxes
c
         pflx(ji) = pflx(ji) + zflx * delts
         ploss(ji) = ploss(ji) + zloss * delts
         psorc(ji) = psorc(ji) + zsorc * delts
c
   98 continue
   99 continue
c
c..15.13 Install mechanism for iterating on loop 99 to improve
c        accuracy of recycling
c
c  Verify that this is a corrector step and that recycling is
c        active (at least one type)
c
      if ((liter-(liter/10)*10.ne.0) .and.
     1    (     ((ji.ne.cfutz(iheflx))
     2          .and.(recflx*zflout + recscr*zflscr.gt.0.))
     3     .or.
     4          ((ji.eq.cfutz(iheflx))
     5          .and.(cfutz(icfout)*zflout
     6                 + cfutz(idvout)*zflscr.gt.0.)))
     7         ) then
c
        zth = 1. - thetap
        zflots = flouto(ji)*zth + flout(ji)*thetap
        zflsts = flscro(ji)*zth + flscr(ji)*thetap
c
c  Compare with previous iteration (also theta-centered)
c
        zdelfx = abs(zflout - zflots) / max(abs(zflout),epslon)
        zdelsc = abs(zflscr - zflsts) / max(abs(zflscr),epslon)
c
c  Test to see if we need another iteration; check counter also
c  Note that the tolerance (0.01) and max. iterations (10) are hardwired
c
        if ((max(zdelfx,zdelsc).gt.0.01).and.(ipc.lt.10)) then
c
          ipc = ipc + 1
c
c  Reset impurity variables as if we were repeating a time-step
c
          do 110 ik=0,nk
            do 110 j=1,mzones
              xn(j,ik,ji) = xnold(j,ik,ji)
  110     continue
c
          do 120 j=1,mzones
            xns(j,ji) = xnsold(j,ji)
  120     continue
c
          pflx(ji) = pflxo(ji)
          ploss(ji) = plosso(ji)
          psorc(ji) = psorco(ji)
c
c  Put new fluxes into common variables and recalculate source
c
          flout(ji) = zflots
          flscr(ji) = zflsts
          call ncsorc(3,ji)
c
          go to 7
        end if
      end if
c
         xns(1,ji) = xns(lcentr,ji)
c
c..15.13 Check for overall maximum (all ji) error
c
      xnerr = max(xnerr,zxnerr)
c
      return
      end

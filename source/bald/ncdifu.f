c--------1---------2---------3---------4---------5---------6---------7-c
c@ncdifu  /11040/baldur/code/bald/ncdifu.f
c  rgb 02-feb-95  corrected units conversion for velthi
c  rgb 12-jul-94  added velthi(4,j) to vxemps  convective velocity
c       dps 07-jul-89 15.11 add D & v settings for neutrals.
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 12-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine ncdifu(ji)
c
c       2.20.3   diffusion- and drift-coefficients
c
c                note: units in this routine are mks !
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'comncr.m'
c
      real :: zgfps  (mj), zgfcl  (mj), zcolfa (mj), zxsih (mj)
      real :: zll   (mj), zgaa   (mj), zmass (mj), znhs (2,mj)
      real :: zxnzs (mj), zxnz1s(mj), zxnsw (mj), zxnzsw(mj)
      real :: zxnziw(mj), zxnzqs(mj), zxnziq(mj)
      real :: zie1  (mj), zie2  (mj), zia   (mj), zia0  (mj)
      real :: zia1  (mj), zia2  (mj)
      real :: ziahxx(mj), ziaaxx(mj), ziahyx(mj), ziaayx(mj)
      real :: ziahxy(mj), ziaaxy(mj), ziahyy(mj), ziaayy(mj)
c
      data  zdielk /8.8543e-12/, zevws /1.6021e-19/, zumass /0.001/
      data  zuchar /10.0/, zumagn /1.0e-04/, zudens /1.0e+06/
      data  zulen /0.01/
      data  idlimt /43/
c
c       functions for neoclassical coefficients
c
         c1st(z) = 1.0 - 3.515625*z / (3.977476+6.765625*z)
         c2st(z) = 1.5 - 1.875*(1.06066+4.3125*z)/(3.977476+6.765625*z)
         c3st(z) = 1.414214 + 3.25*z - (1.06066+4.3125*z)**2
     x                                           /(3.977476+6.765625*z)
         fd3(z) = 2.607456 + 3.25*z - (2.303621+4.3125*z)**2
     x                                           /(5.790127+6.765625*z)
         fd4(z) = (1.0+2.612*z+1.343*z**2+0.1719*z**3)
     x                         /(1.683+2.852*z+0.5828*z**2)
         fd5(z) = -(0.3965+1.121*z+0.3545*z**2)
     x                         /(1.683+2.852*z+0.5828*z**2)
         falf(z) = (1.0+1.198153*z+0.222222*z**2)
     x                         / (1.0+2.96592*z+0.753472*z**2)
c
      if (nstep.le.0)   return
c
c
c    calculation of diffusion-coefficients and drift-velocities
c
c  note: da and va are defined as boundary-values; e.g.
c        da(lcentr) at xbouni(lcentr), da(mzones) at xbouni(mzones).
c
      do ik=1,nk
c
        zfactd = 1.0
        if ( abs(daqexp(ji)) .gt. epslon )
     &     zfactd = (float(ik))**daqexp(ji)
        zfactv = 1.0
        if ( abs(vaqexp(ji)) .gt. epslon )
     &     zfactv = (float(ik))**vaqexp(ji)
c
        do j=lcentr,mzones
          da(j,ik,ji) = dnis(ji,j) * zfactd
          va(j,ik,ji) = - ( vxemps(ji+lhydn,j) + velthi(4,j)*uisl )
     1                    * zfactv
        enddo
c
      enddo
c
      if (ifneo.gt.0)   go to 10
      go to 100
c
c    neoclassical part of diffusion-coefficients and drift-velocities
c
   10 continue
         jimp = limp1 - 1 + ji
         zaimp = aspec(jimp)
         zhmass = fcmp * 10.0**fxnucl * zumass
         zimass = zaimp * zhmass
         zemass = fcme * 10.0**fxme * zumass
         zcharg = fce * 10.0**fxe * zuchar
         zbtor = bzs * zumagn
         lcent1 = lcentr + 1
         zxnmin = 1.0
c
      do 11 j=lcentr,mzones
         zxnzs(j) = 0.0
         zxnz1s(j) = 0.0
         zxnsw(j) = 0.0
         zxnzsw(j) = 0.0
         zxnziw(j) = 0.0
         zxnzqs(j) = 0.0
         zxnziq(j) = 0.0
         zmass(j) = 0.0
   11 continue
c
      do 12 ik=1,nk
      do 12 j=lcentr,mzones
         zxnzs(j) = zxnzs(j) + ik*xn(j,ik,ji)*zudens
         zxnz1s(j) = zxnz1s(j) + xn(j,ik,ji)*zudens/float(ik)
   12 continue
c
      do 14 ik=1,nk
         zik = ik
c
      do 14 j=lcent1,mzones
         zint0 = (xzoni(j)-xbouni(j)) / dxboui(j)
         zint1 = 1.0 - zint0
         zxn = (xn(j-1,ik,ji)*zint0 + xn(j,ik,ji)*zint1) * zudens
         zxnsw(j) = zxnsw(j) + zxn
         zxnzsw(j) = zxnzsw(j) + zxn*zik
         zxnziw(j) = zxnziw(j) + zxn/zik
         zxnziq(j) = zxnziq(j) + zxn/zik**2
         zxnzqs(j) = zxnzqs(j) + zxn*zik**2
   14 continue
c
      do 15 j=lcent1,mzones
         zxnmin = min(zxnmin,zxnsw(j))
   15 continue
c
      do 16 i=1,2
      do 16 j=lcentr,mzones
         znhs(i,j) = 0.0
   16 continue
c
      do 17 jh=1,mhyd
      do 17 j=lcentr,mzones
         znhs(2,j) = znhs(2,j) + chi(jh,j)*uisd*zudens
         znh = (chi(jh,j-1)*(xzoni(j)-xbouni(j))
     x                  + chi(jh,j)*(xbouni(j)-xzoni(j-1))) / dxboui(j)
         znhs(1,j) = znhs(1,j) + znh*uisd*zudens
   17 continue
c
      do 18 jh=1,mhyd
      do 18 j=lcent1,mzones
         znh = (chi(jh,j-1)*(xzoni(j)-xbouni(j))
     x                  + chi(jh,j)*(xbouni(j)-xzoni(j-1))) / dxboui(j)
         zmass(j) = zmass(j)
     x                + zhmass*aspec(jh) * znh*uisd*zudens / znhs(1,j)
   18 continue
c
c    if no density at one zone-boundary, return
c
      if (zxnmin.le.0.0)    return
c
      do 20 j=lcent1,mzones
         zrne = rhoels(1,j) * zudens
         zte = tes(1,j) * evsinv
         zti = tis(1,j) * evsinv
         zeps = xbouni(j)*rmins/rmajs
         zgfps(j) = (q(j)/zeps)**2
     x                * (1.0+1.5*zeps**2-sqrt(1.0-zeps**2)) / zbtor**2
         zgfcl(j) = (1.0+1.5*zeps**2) / zbtor**2
c
         zarg = 1.0e-03 * sqrt(zrne/zte**3)
         zlnla = 23.0 - log(zarg)
         zcolfa(j) = 2.116455e-02 * zlnla * zcharg**4
     x                           / (zdielk**2 * zevws * sqrt(zevws))
         zxsifa = 2.116455e-02 * zlnla / (zdielk**2*zevws*sqrt(zevws))
         zxsih(j) = zxsifa * sqrt(zmass(j)/zti**3)
         zxsie = zxsifa * sqrt(zemass/zte**3)
c
         zze = (znhs(1,j)+zxnzqs(j)) / zrne
         zzh = zxnzqs(j) / znhs(1,j)
         zza = 0.0
         zfak = zcolfa(j) * zrne**2 * sqrt(zemass/zte**3) / zcharg**2
         zie1(j) = -zfak * (c1st(zze)*zgfps(j)+zgfcl(j))
         zie2(j) = zfak * (c2st(zze)*zgfps(j)+1.5*zgfcl(j))
c
         zcolhh = zcolfa(j) * sqrt(zmass(j)/zti**3) * znhs(1,j)**2
         zcolaa = zcolfa(j) * sqrt(zimass/zti**3) * zxnzqs(j)**2
         zcolha = zcolfa(j) * sqrt(zmass(j)/zti**3)
     1              * znhs(1,j)*zxnzqs(j)
         zcolea = zcolfa(j) * sqrt(zemass/zte**3) * zrne*zxnzqs(j)
         zfak = -zrne / (zcharg**2 * zxnzqs(j))
         zsum1 = zcolaa * zgfcl(j) * (zza + 0.707107)
         zsum2 = (zgfps(j)+zgfcl(j)) * (zcolha+zcolea)
         zia(j) = zfak * (zcolaa*zgfps(j)*fd4(zza) + zsum1 + zsum2)
c
         zphia = zxnzqs(j)*zxnziq(j)/zxnsw(j)**2  -  1.0
         zll(j) = -zxnsw(j)**2/zcolaa *(1.0/c3st(zza)+zphia/fd3(zza))
     x                              - znhs(1,j)**2/(zcolhh*c3st(zzh))
         zsum1 = zza * zxnsw(j) * c2st(zza) / (zxnzqs(j) * c3st(zza))
         zsum2 = zxnsw(j) * fd5(zza) / zxnzqs(j)
         zgaa(j) = (-zsum1 + zsum2 + c2st(zzh)/c3st(zzh)) / zcharg**2
         zgah = (-zzh*c2st(zzh)/c3st(zzh) + fd5(zzh)) / zcharg**2
         zfak = zgfps(j) * (zrne*zcharg)**2
         zfak2 = zgfcl(j) * (zrne*zcharg)**2
         zsum1 = zxsie * (1.0-c1st(zze)) / zze
         zsum2 = zgaa(j) * zgah / zll(j)
         ziahxx(j) = zfak * (falf(zzh)*zxsih(j) + zsum1 + zsum2)
     x                                          + zfak2 * zxsih(j)
         zsum3 = -zxsih(j) * sqrt(zaimp) * (zza*falf(zza) - fd4(zza))
         ziaaxx(j) = zfak * (zsum3 + zxsih(j)*(1.0-falf(zzh))/zzh
     x    +zsum1+zgaa(j)*2/zll(j))+zfak2*0.707107*zxsih(j)*sqrt(zaimp)
         ziahyx(j) = -zgfps(j)*zrne*zxnziw(j)*fd5(zza)*zgah/zll(j)
   20 continue
c
      do 25 j=lcent1,mzones
         zrne = rhoels(1,j) * zudens
         zti = tis(1,j) * evsinv
         zza = 0.0
         zzh = zxnzqs(j) / znhs(1,j)
         zcolaa = zcolfa(j) * sqrt(zimass/zti**3) * zxnzqs(j)**2
         ziaayx(j) = -zgfps(j)*zrne*zxnziw(j)*fd5(zza)*zgaa(j)/zll(j)
         ziahxy(j) = -zgfps(j)*zrne*zxnziw(j)*fd5(zzh)*zgaa(j)/zll(j)
         ziaaxy(j) = ziaayx(j)
         ziahyy(j) = zgfps(j) * zxnziw(j)**2 * fd5(zza) * fd5(zzh)
     x                                          / (zcharg**2 * zll(j))
         ziaayy(j) = zgfps(j)*zxnziw(j)**2
     1                 *fd5(zza)**2/(zcharg**2*zll(j))
c
         zfak = zcolaa * (0.353553+zza) / (zcharg**2 * zxnzqs(j))
         zia0(j) = zgfcl(j) * 1.5 * zrne * zfak
         zsum1 = zgfps(j) * zxnziw(j) * zgaa(j) / (zcharg * zll(j))
         zsum2 = 0.53033*zxsih(j)*sqrt(zaimp)*zxnzsw(j)
     x                                    + 1.5*zxsih(j)*znhs(1,j)
         zia1(j) = zcharg * zrne * (zsum1 - zgfcl(j)*zcharg*zsum2)
         zia2(j) = -zgfps(j)*zxnziw(j)**2
     1                *fd5(zza)/(zcharg**2 * zll(j))
c
         do 24 ik=1,nk
           da(j,ik,ji) = da(j,ik,ji)
     1                - zudens*zulen*zia(j)*zevws*zti/zrne
   24    continue
   25 continue
c
      do 27 j=lcent1,mzones
         zdr = evsinv*zevws / (dxboui(j)*rmins*zulen**2)
         zrhoe = rhoels(1,j)*zudens
         zgse = znhs(1,j) + zxnziw(j)
         zdtis = tis(2,j) - tis(2,j-1)
         zvsu1 = zie1(j)
     1            * (rhoels(2,j)*tes(2,j)-rhoels(2,j-1)*tes(2,j-1))
     x            / (zrhoe*rhoels(1,j))
     3            + zie2(j)*(tes(2,j)-tes(2,j-1))/zrhoe
         zvsu2 = zia(j) * zdtis / zrhoe
c
         zfak = tis(1,j) * (znhs(2,j)-znhs(2,j-1)) + znhs(1,j) * zdtis
         zfak1 = tis(1,j)*(zxnzs(j)-zxnzs(j-1)) + zxnzsw(j)*zdtis
         zfak2 = tis(1,j)*(zxnz1s(j)-zxnz1s(j-1)) + zxnziw(j)*zdtis
c
         do 26 ik=1,nk
         zvsu3 = ik*ziahxx(j)/zrhoe**2 + ziahyx(j)/(ik*zgse*zrhoe)
         zvsu4 = ik*ziaaxx(j)/zrhoe**2 + ziaayx(j)/(ik*zgse*zrhoe)
c?     x                             + ziaayx(j)/(ik*zgse*zrhoe)
         zvsu5 = (ik*ziahxy(j)/zrhoe + ziahyy(j)/(ik*zgse)) / zgse
         zvsu6 = (ik*ziaaxy(j)/zrhoe + ziaayy(j)/(ik*zgse)) / zgse
         zvsu7 = (zia0(j)/zrhoe+ik*zia1(j)/zrhoe
     x                                 +zia2(j)/(ik*zgse))  *  zdtis
         va(j,ik,ji) = va(j,ik,ji)
     x            + zdr * (ik*zvsu1 + zvsu2 + zfak*zvsu3 + zfak1*zvsu4
     x                            + zfak*zvsu5 + zfak2*zvsu6 + zvsu7)
   26 continue
   27 continue
c
  100 continue
c
c  Enforce limit on D
c
      if (cfutz(idlimt).gt.epslon) then
        do 110 ik=1,nk
          do 110 j=lcentr,mzones
            zd = da(j,ik,ji)
            da(j,ik,ji) = min(zd,cfutz(idlimt))
  110   continue
      end if
c
c..15.11 Set values for neutrals
c
      do 120 j=lcentr,mzones
        da(j,0,ji) = 1.e6
        va(j,0,ji) = 0.
  120 continue
c
      return
      end

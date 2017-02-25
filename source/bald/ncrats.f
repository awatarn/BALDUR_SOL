c--------1---------2---------3---------4---------5---------6---------7-c
c@ncrats  /11040/baldur/code/bald/dimprad.f
c  rgb 16-dec-94 use cimprd(10) as coefficient for sources and sinks
c  rgb 14:28 11-nov-93 changed call abort to call abortb
c       dps 07-jul-89 15.11 use ik=0 for neutral impurities.
c       dps 15-may-89 15.09 replace usage of nimp with nkimp
c       dps 29-sep-88 15.05 add conversion from W to erg to radbcx.
c       dps 24-aug-88 15.00 incorporate routines into DIMPRAD
c       dps 12-aug-88 Begin adapting routines to 1-1/2-D BALDUR.
c       rhw 13.06.88: update calc. of ionis. + recomb. rates acc. to
c                     behringer: jet report jet-r(87)08
c       rhw 10/08/84: change format "9001".
c       rhw 12/07/84: preset alm in ncdata, read alm by namelist ncodat
c       rhw 03/07/84: change "5" to "nin" in read statements.
c       rhw 02/07/84: select species, read aspec, separate radiations.
c       rhw 17/05/84: save local variables, calculate s0
c       rhw 12/03/84: interpolate rate-coeffs. exponentially
c       rhw 02/03/84: first write-up
c***********************************************************************
c
        subroutine ncrats(k,ji)
c
c       2.20.2   ionisation- and recombination-rates for non-corona
c                radiation
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'comncr.m'
      include 'comadp.m'
c
      common  /comrat/
     r    ztel(knte)      , znel(knne)  , zsal(kzmax,knte,kncimp)  ,
     r    zral(kzmax,knne,knte,kncimp)  , zbradl(kzmax,knte,kncimp),
     r    zsradl(kzmax,knte,kncimp)     , zlradl(kzmax,knte,kncimp),
     r    zrradl(kzmax,knne,knte,kncimp), zen(100,10,kncimp)       ,
     r    zdtlg  , ztealg , zdnlg  , znealg ,
     i    nanspc(2)       , navaln(100,kncimp),
     i    nk1    , nte    , nne
!cap
      real ::  zwbrem(mj), znv(mj), zrcx(mj,28,2),
     1           zradcx(mj,kzmax,kncimp)
c
      equivalence (wchex(1,1,1),zrcx(1,1,1)) , (zwbrem(1),znv(1))
c
      save  /comrat/
c
      go to (100,200,300), k
c
  100 continue
c
c    data for ionisation- and recombination-rates
c
c  Initialize ADPAK routines for this impurity
c
      laden = naden(ji)
      ladtip = nadtip(ji)
      ldrmlt = ndrmlt(ji)
      leci = neci(ji)
      yhnna = yhnn(1,ji)
      yhnnb = yhnn(2,ji)
      yhnnc = yhnn(3,ji)
      yhnma = yhnm(1,ji)
      yhnmb = yhnm(2,ji)
      yhnmc = yhnm(3,ji)
c
c..15.11 Offset indices here to match MIST convention
c
      do 5 ik=2,kzmax+1
        cdnn(ik) = dnnmlt(ik-1,ji)
        cdnm(ik) = dnmmlt(ik-1,ji)
    5 continue
c
      call adset(nkimp(ji))
c
c  Save information needed to compute charge exchange
c  recombination so that we don't have to reinitialize
c  entire package later on
c
      nanspc(ji) = nspc
      do 10 j=1,100
        navaln(j,ji) = nvalnc(j)
        do 10 j1=1,10
          zen(j,j1,ji)=en(j,j1)
   10 continue
c
      nk1 = nk - 1
c
c    data to set up tables of ionisation- and recombination-rates
c
         ztea = 1.0
         ztee = 1.0e+05
         nte = knte
         ztealg = log10(ztea)
         zteelg = log10(ztee)
         zdtlg = (zteelg - ztealg) / (float(nte) -1.)
c
         znea = 1.0e+10
         znee = 1.0e+15
         nne = knne
         znealg = log10(znea)
         zneelg = log10(znee)
         zdnlg = (zneelg - znealg) / (float(nne) - 1.)
c
      do 12 i=1,nte
         ztel(i) = ztealg + (i-1)*zdtlg
   12 continue
c
      do 13 j=1,nne
         znel(j) = znealg + (j-1)*zdnlg
   13 continue
c
c  Call ADPAK routines to set up tables of rates. Note that
c  these routines assume T is in keV and n is in cm**(-3).
c
      do 30 i=1,nte
        zte = 10.0**ztel(i) / 1.e3
c
c  Set n = 1. so that rates are not functions of the
c  electron density. The required overall multiplier
c  is added during the actual solution procedure.
c
        zne = 1.0
        call adeci(zte,zne)
        call adrrec(zte,zne)
        call adecex(zte,zne)
        call adbrem(zte,zne)
c
c  Note that reaction rates have units of cm**3/sec, and radiation
c  rates have units of W*cm**3. The latter are converted to
c  erg*cm**3/sec by adding 7.0 to the logarithm.
c
        do 15 ik=1,nk
          zsal(ik,i,ji)   = log10( cizmlt(ji)*rclion(ik) + epslon )
          zsradl(ik,i,ji) =
     1          log10( cizmlt(ji)*cizlos(ik) + epslon ) + 7.0
          zlradl(ik,i,ji) = log10( radclx(ik) + epslon ) + 7.0
          zbradl(ik,i,ji) = log10( radbrm(ik+1) + epslon ) + 7.0
   15   continue
c
c  Now calculate recombination rates. The density
c  dependence is intrinsic now in the dielectronic contribution,
c  but we divide out the overall factor for consistency.
c
        do 25 j=1,nne
c
          zne = 10.0**znel(j)
          call addrec(zte,zne)
c
          do 20 ik=1,nk
            zral(ik,j,i,ji)   = log10( cermlt(ji)*(rrarec(ik+1)
     1                           + rdirec(ik+1)/zne) + epslon )
            zrradl(ik,j,i,ji) = log10( cermlt(ji)*(radrrc(ik+1)
     1                           + raddrc(ik+1)/zne) + epslon ) + 7.0
   20     continue
   25   continue
   30 continue
c
      return
c
c-----------------------------------------------------------------------
c
  200 continue
c
      nk1 = nk - 1
c
c    interpolation of rate-coefficients s0, sa, and ra
c
      do 202 ik=0,nk
         sa(1,ik,ji) = 0.0
         ra(1,ik,ji) = 0.0
  202 continue
c
      do 204 j=lcentr,mzones
         ztev = tes(2,j)*evsinv
         ztelg = log10(ztev)
         znelg = log10(rhoels(2,j))
         ite = (ztelg-ztealg)/zdtlg + 1.99999
         ine = (znelg-znealg)/zdnlg + 1.99999
         ite = max0 (ite,2)
         ine = max0 (ine,2)
         ite = min0 (ite,nte)
         ine = min0 (ine,nne)
         zint = (ztel(ite)-ztelg) / (ztel(ite)-ztel(ite-1))
         zint1 = (znel(ine)-znelg) / (znel(ine)-znel(ine-1))
c
c  Neutral ionization rates
c
         zslg = zsal(1,ite-1,ji)*zint + zsal(1,ite,ji)*(1.0-zint)
         s0(j,ji) = 10.0**zslg
c
c  Ionization rates of higher stages
c
         do 203 ik=0,nk1
           zslg = zsal(ik+1,ite-1,ji)*zint
     1             + zsal(ik+1,ite,ji)*(1.0-zint)
           sa(j,ik,ji) = 10.0**zslg
  203    continue
c
c  No ionization from stage nk to nk+1
c
         sa(j,nk,ji) = 0.0
c
c  Recombination rates
c
c..15.11 No recombination of neutrals
c
         ra(j,0,ji) = 0.0
c
         do 204 ik=1,nk
           zra1 = zral(ik,ine-1,ite-1,ji)*zint
     x                           + zral(ik,ine-1,ite,ji)*(1.0-zint)
           zra2 = zral(ik,ine,ite-1,ji)*zint
     x                           + zral(ik,ine,ite,ji)*(1.0-zint)
           zrlg = zra1*zint1 + zra2*(1.0-zint1)
           ra(j,ik,ji) = 10.0**zrlg
  204 continue
c
c  Call ADPAK routines to calculate charge
c  exchange recombination rates
c
      nspc = nanspc(ji)
c
      do 210 j=1,100
        nvalnc(j) = navaln(j,ji)
        do 210 j1=1,10
          en(j,j1) = zen(j,j1,ji)
  210 continue
c
      ivunit = 1
      incxb = 1
c
c  Must add contributions from various hydrogenic
c  species. So, will initialize these arrays and
c  then sum the recombination rates in them.
c
      do 212 ik=1,nk
        do 212 j=lcentr,mzones
          zrcx(j,ik,ji) = 0.0
          zradcx(j,ik,ji) = 0.0
  212 continue
c
      do 220 jh=1,mhyd
        do 220 j=lcentr,mzones
          zhn=rhons(jh,j)
c
c  Neutral thermal velocity is in cm/sec
c
          zhv = sqrt( 2.0*tns(jh,j)/(aspec(jh)*fcau*10.**fxnucl) )
          icxx = ncxx(jh)
          if ( zhn.lt.1.0 .or. zhv.eq.0.0 ) go to 220
c
          call adbcxr(zhn,zhv,incxb,icxx,icxerr,ivunit)
c
          if (icxerr.ne.0) call abortb(nprint,'ADBCXR error')
c
          do 215 ik=1,nk
            zdcx = cxrmlt(ji) * rcxrec(ik+1) / rhoels(2,j)
            zrcx(j,ik,ji) = zrcx(j,ik,ji) + zdcx
            ra(j,ik,ji) = ra(j,ik,ji) + zdcx
c
c  15.05 The factor of 1.e7 multiplying radbcx is required to
c  convert from W/cm**3 to ergs / s cm**3.
c
            zradcx(j,ik,ji) = zradcx(j,ik,ji)
     1                           + cxrmlt(ji)*radbcx(ik+1)*1.e7
  215     continue
c
  220 continue
c
      return
c
c-----------------------------------------------------------------------
c
  300 continue
c
c    interpolation of radiation and cooling rates
c
      do 310 j=lcentr,mzones
         ztev = tes(2,j)*evsinv
         ztelg = log10(ztev)
         znelg = log10(rhoels(2,j))
         ite = (ztelg-ztealg)/zdtlg + 1.99999
         ine = (znelg-znealg)/zdnlg + 1.99999
         ite = max0 (ite,2)
         ine = max0 (ine,2)
         ite = min0 (ite,nte)
         ine = min0 (ine,nne)
         zint = (ztel(ite)-ztelg) / (ztel(ite)-ztel(ite-1))
         zint1 = (znel(ine)-znelg) / (znel(ine)-znel(ine-1))
c
c  Neutral impurity
c
c  Line radiation due to electron collisional excitation of neutrals
c..15.11 Include now both types of neutrals
c
        zradlg = zlradl(1,ite-1,ji)*zint + zlradl(1,ite,ji)*(1.0-zint)
        wline(j,1,ji) = cimprd(10) * (10.0**zradlg) * rhoels(2,j)
     1                  * (xn0(j,ji) + xn(j,0,ji))
c
c  Radiation due to recombination of stage 1 ions to neutrals
c  (excluding charge exchange)
c
        zradl1 = zrradl(1,ine-1,ite-1,ji)*zint
     1             + zrradl(1,ine-1,ite,ji)*(1.0-zint)
        zradl2 = zrradl(1,ine,ite-1,ji)*zint
     1             + zrradl(1,ine,ite,ji)*(1.0-zint)
        zradlg = zradl1*zint1 + zradl2*(1.0-zint1)
        wreko(j,1,ji) = cimprd(10) * 
     &    (10.0**zradlg) * rhoels(2,j) * xn(j,1,ji)
c
c  Radiation due to charge exchange of stage 1 ions to neutrals
c
        wchex(j,1,ji) = cimprd(10) * zradcx(j,1,ji) * xn(j,1,ji)
c
c  Electron cooling by ionization of source neutrals plus
c  ionization of neutrals formed by recombination of stage 1 ions
c..15.11 The latter replaces the previous "auto-ionization" assumption
c
        zradlg = zsradl(1,ite-1,ji)*zint + zsradl(1,ite,ji)*(1.0-zint)
        wioni(j,1,ji) = cimprd(10) * (10.0**zradlg)
     1                  * (dq(j,ji) / max(s0(j,ji),epslon)
     2                     + rhoels(2,j) * xn(j,0,ji))
c
c  Electron cooling bremsstrahlung radiation off of stage 1 impurity
c
        zradlg = zbradl(1,ite-1,ji)*zint + zbradl(1,ite,ji)*(1.0-zint)
        wbrem(j,1,ji) = cimprd(10) * 
     &    (10.0**zradlg) * rhoels(2,j) *xn(j,1,ji)
c
c  Total electron cooling power due to presence of neutral impurity
c
        weir(j,1,ji) = wline(j,1,ji) + wioni(j,1,ji)
        weirs(ji,j)  = weir(j,1,ji)
c
c  Higher stages
c
        do 305 ik=2,nk
c
c  Line radiation due to electron collisional excitation of
c  stage ik-1 ions
c
          zradlg = zlradl(ik,ite-1,ji)*zint
     1               + zlradl(ik,ite,ji)*(1.0-zint)
          wline(j,ik,ji) = cimprd(10) * (10.0**zradlg)
     1               * rhoels(2,j) * xn(j,ik-1,ji)
c
c  Radiation due to recombination of stage ik ions
c  (excluding charge exchange)
c
          zradl1 = zrradl(ik,ine-1,ite-1,ji)*zint
     1               + zrradl(ik,ine-1,ite,ji)*(1.0-zint)
          zradl2 = zrradl(ik,ine,ite-1,ji)*zint
     1               + zrradl(ik,ine,ite,ji)*(1.0-zint)
          zradlg = zradl1*zint1 + zradl2*(1.0-zint1)
          wreko(j,ik,ji) = cimprd(10) * 
     &      (10.0**zradlg) * rhoels(2,j) * xn(j,ik,ji)
c
c  Radiation due to charge exchange of stage ik ions
c
          wchex(j,ik,ji) = cimprd(10) * zradcx(j,ik,ji) * xn(j,ik,ji)
c
c  Electron cooling by ionization of stage ik-1 ions
c
          zradlg = zsradl(ik,ite-1,ji)*zint
     1               + zsradl(ik,ite,ji)*(1.0-zint)
          wioni(j,ik,ji) = cimprd(10) * (10.0**zradlg)
     1               * rhoels(2,j) * xn(j,ik-1,ji)
c
c  Electron cooling bremsstrahlung radiation off of stage ik ions
c
          zradlg = zbradl(ik,ite-1,ji)*zint
     1               + zbradl(ik,ite,ji)*(1.0-zint)
          wbrem(j,ik,ji) = cimprd(10) * 
     &      (10.0**zradlg) * rhoels(2,j) * xn(j,ik,ji)
c
c  Electron cooling power due to presence of stage ik-1 impurity
c
          weir(j,ik,ji) = 
     &      wline(j,ik,ji) + wioni(j,ik,ji) + wbrem(j,ik-1,ji)
          weirs(ji,j) = weirs(ji,j) + weir(j,ik,ji)
c
  305   continue
c
c  For stage nk ions (fully stripped), only bremsstrahlung
c  gives rise to electron energy losses. We could include the
c  following statement:
c       weir(j,nk+1,ji) = wbrem(j,nk,ji)
c  but weir is not dimensioned properly. For now, just add
c  to total impurity radiation.
c
        weirs(ji,j) = weirs(ji,j) + wbrem(j,nk,ji)
c
  310 continue
c
      return
c
      end

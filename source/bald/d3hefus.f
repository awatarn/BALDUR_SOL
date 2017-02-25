c 20:00 27-sep-92 .../baldur/code/bald/d3hefus.f
c/ mar-1989 /wbaldn1/DD3HFUS, Sugiyama, MIT
c
c  rgb 31-may-92 removed sbrtns syndrv, synch, opacty, and editor
c    to file baldur/code/new/dsynch_sugiyama.f since the Tamor
c    synchrotron package is already installed in BALDUR
c  rgb 31-may-92 removed function seval
c     which already exists in the utility library yutil.a
c  ahk 3-jun-96 replace zf(1),zpf(1) with zf(*),zpf(*)
c  pis 18-jun-98 trucation real -> integer resolved for idum
c************************************************************************
c
c..This file contains the following subroutines:
c
c   fusdrv    ! driver for d-3he fusion; called in coef
c   fusion    ! d-3he cycle fusion, including simple slowing down
c             !fusion(0,..) in STEPON to update fast particles
c             ! called in FUSDRV to calculate fusion reactions
c  block data dt3fus   ! fusion reaction coefficients; spline fit
c              ! to MILEY, et al.
c   igntst    ! test for ignition at each time step; called in stepon
c             ! also tests beta_T(r=0) = 1 and
c             ! savestotal powers for printout to plots
c   reorder   ! reorder ions for neoclassical transport
c   cmg       ! CMG transport chi-e and D-e, V-e
c   auxanl   ! simple analytic auxiliarty heating and particle sources
c             ! called in coef
c   bussac     ! calc bussac's ideal MHD m=1 mode critical beta-p
c   d3prin    ! short printout for d-3He, called in OUTPUT
c   dprint    ! long radial printout for d-3He
c   volint    ! utility volume integration
c   syndrv    ! driver for Tamor's synchrotron radiation
c             ! called in IMPRAD
c   synch      ! Tamor's synchrotron calc., modified
c   opacty     ! Tamor
c   editor     ! Tamor
c
c************************************************************************
c@fusdrv   ---   les version for d 3he   --- .../baldur/code/d3hefus.f
c  rgb  04-mar-96  introduce integer argument for sbrnt fusdrv
c    fusdrv(1) -> initialize arrays rh1fst and rh2fst
c  les  16-jan-91  remove first orbit electron energy losses from
c        wed3fs
c     input march 1989
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine fusdrv(k)
c     driver for fusion d-3he
c
c
c************************************
c@cd3he
c
c     d-3he fusion
c
c     cliche cd3he
c     common/testf/tppd(n9),tpp3(n9),tpt(n9),tp3(n9),tp4t(n9),
c    1 tp43(n9),tp42t(n9),tepd(n9),tep3(n9),tet(n9),te3(n9),te4t(n9),
c    2 te43(n9),te42t(n9),fpde(n9),fp3e(n9),fte(n9),f3e(n9),f4te(n9),
c    3 f43e(n9),f42te(n9),fddt(n9),fdd3(n9),fdt(n9),fdtb(n9),fd3(n9),
c    4 fd3b(n9),ftt4(n9),fttb(n9),
c    r   epde(n9),ep3e(n9),etre(n9),e3e(n9),e4te(n9),e43e(n9),
c    r   e42te(n9),epdi(n9),ep3i(n9),etri(n9),e3i(n9),e4ti(n9),
c    r   e43i(n9),e42ti(n9),
c    r   drnp(n9),drnd(n9),drnt(n9),drn3(n9),drn4(n9)
c       common/d3fus/q1fus(n9),q2fus(n9),q1lfus(n9),q2lfus(n9),
c    r   qefus(n9),wntn(n9),d3fast(n9)
c      common/d3fstp/fnpd0(n9),fnp30(n9),fnt0(n9),fn30(n9),fn4t0(n9),
c    1 fn430(n9),fn42t0(n9),fnpd(n9),fnp3(n9),fnt(n9),fn3(n9),
c    2 fn4t(n9),fn43(n9),fn42t(n9),rh1fst(2,n9),rh2fst(2,n9),
c    3 epd0(n9),ep30(n9),et0(n9),e30(n9),e4t0(n9),e430(n9),e42t0(n9),
c    4 epd(n9),ep3(n9),et(n9),e3(n9),e4t(n9),e43(n9),e42t(n9)
c     endcliche
c**********************************************************************c
c   common block cmdif2 should appear in cliche cbaldr
c     common/cmdif2/ shfuel, sifuel,
c    r wed3fs, wid3fs, wid3fl, weash, shd3fs, sid3fs, wesyn,
c    r wefus, wifus
c       dimension
c    r shfuel(2,55), sifuel(4,55), shd3fs(2,55), sid3fs(4,55),
c    r wed3fs(mj), wid3fs(mj), wid3fl(mj), weash(mj), wesyn(mj),
c    r wefus(mj), wifus(mj)
c
c************************************
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cd3he.m'
c   transfer values to fusion subr
      common/d3hef/rnp(n9),rnd(n9),rnt(n9),rn3(n9),rn4(n9),rn(n9),
     a  tefus(n9),tifus(n9),pfols(n9),wefol(n9),
     a   dtf
      data id3he/490/
c
c..initialize
c
      i1 = 2*n9
      call resetr (rh1fst,i1,0.0)
      call resetr (rh2fst,i1,0.0)
c
      if ( k .lt. 2 ) return
c
c      time step (sec)
      dtf=dti*uist
c   densities in cm-3, temps in keV
c     note that protons are a BALDUR impurity species
c     deuterium and tritium are the 'hydrogenic' species
c     lhydn=2; ldeut=1, ltrit=2, lprotn=3, lhe3=4,lalpha=5, for example
c
      do 4 jz=lcentr,ledge
        rn(jz)    = rhoels(2,jz)
        rnp(jz)   = rhois(lprotn-lhydn,2,jz)
        rnd(jz)   = rhohs(ldeut,2,jz)
        rnt(jz)   = rhohs(ltrit,2,jz)
        rn3(jz)   = rhois(lhe3-lhydn,2,jz)
        rn4(jz)   = rhois(lalpha-lhydn,2,jz)
        tifus(jz) = tis(2,jz)*useh
   4    tefus(jz) = tes(2,jz)*useh
c
c   fusion slowing down model
c      cfutz(id3he)=1  direct depostion
c      cfutz(id2he)=2  central time (t+dt/2) weighted slowing down
c
      kfusn=cfutz(id3he)+1.0e-06
      call fusion(kfusn,lcentr,ledge,fploss,f4loss)
c
c   source and sink terms in cgs units
c     to transfer these variables to the rest of BALDUR,
c      use common block ctmp03 inside cbaldr and cliche cd3he
c
      do 20 jz=lcentr,ledge
c  power deposition from slowing down of fusion fast particles
c   wed3fs - to electrons (erg/cm3/s)
c   wid3fs - to ions
c
        wed3fs(jz)=qefus(jz)
c
c   wefol(jz)=pfols(jz)*1.5*tefus(jz)*uesh
c
        wid3fs(jz)=(q1fus(jz)+q2fus(jz))
c
c  miscellaneous - power lost through particle losses
c   wid3fl - from ions due to fusion burnup of thermal ions
c   weash - from electrons due to imposed ash removal of ions
c       (maintain charge neutrality)
c
        wid3fl(jz)=(q1lfus(jz)+q2lfus(jz))
        weash(jz)=1.5*tefus(jz)*(ashp*drnp(jz)+ash4*drn4(jz))*
     &    uesh
c
c   net thermal particle gain from fusion for each ion species
c      (number/cm3/sec)
c
        shd3fs(ldeut,jz)        = drnd(jz)
        shd3fs(ltrit,jz)        = drnt(jz)
        sid3fs(lprotn-lhydn,jz) = (1.-ashp)*drnp(jz)
        sid3fs(lhe3-lhydn,jz)   = drn3(jz)
  20    sid3fs(lalpha-lhydn,jz) = (1.-ash4)*drn4(jz)
      return
      end
c@fusion
c  rgb 27-sep-92 removed (j,...) from seval argument list
c  rgb 27-sep-92 bounded zti by sx(1) and sx(nsv)
c  rgb 31-may-92 removed parameter(nsv=12) which is already in cd3he.m
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine fusion(ifusn,lc,le,fploss,f4loss)
c
c       D,T,3He fusion reactions; subr works for any inital thermal ion
c         densities
c       Values from current time t have suffix 0, from t+dt no suffix
c   fusion power and particle source rate densities
c   must be updated by call fusion(0) after sucessful time step
c
c      parameter(nsv=12)
      include 'cd3he.m'
      common/d3hef/rnp(n9),rnd(n9),rnt(n9),rn3(n9),rn4(n9),rn(n9),
     a  tefus(n9),tifus(n9),pfols(n9),wefol(n9),
     a   dtf
      dimension sx(nsv),syddt(nsv),sbddt(nsv),
     a  scddt(nsv),sdddt(nsv),
     b  sydd3(nsv),sbdd3(nsv),scdd3(nsv),sddd3(nsv),
     c  sydt(nsv),sbdt(nsv),scdt(nsv),sddt(nsv),
     d  syd3(nsv),sbd3(nsv),scd3(nsv),sdd3(nsv),
     e  sytt4(nsv),sbtt4(nsv),sctt4(nsv),sdtt4(nsv),
     f  sxb(7),sydtb(6),sbdtb(6),scdtb(6),sddtb(6),syttb(7),sbttb(7),
     g  scttb(7),sdttb(7),syd3b(6),sbd3b(6),scd3b(6),sdd3b(6)
c
      dimension svd3(n9),svd3b(n9),svddt(n9),svdd3(n9),svdt(n9),
     1 svtt4(n9),svttb(n9),svdtb(n9),zyc(n9),ztc(n9),
     2 ztemp(n9)
      dimension zstarb(n9),zxlme(n9)
c
c  atomic mass am_i = m_i/m_p
c
      data amd/2./,amt/3./,am3/3./,am4/4./,zii1/1./,zii2/2./
c
c     block data dt3fus
c   sx(b) = ln(T_i(keV))
c   sy    = ln(<sigma-v> cm3/s)
c
      data sx/0.0000e+0,6.9310e-1,1.3863e+0,1.9459e+0,2.3026e+0,
     a  2.9957e+0,3.6889e+0,4.2485e+0,4.6052e+0,5.2983e+0,
     b  5.9915e+0,6.5511e+0/
      data sxb/ 4.6050e+0, 5.2980e+0, 5.9910e+0, 6.5510e+0, 6.9080e+0,
     a 7.6010e+0, 8.2940e+0/
c   Maxwellian-Maxwellian rates
      data syddt/-5.0843e+1,-4.7318e+1,-4.4611e+1,-4.2896e+1,-4.1989e+1,
     a -4.0557e+1,-3.9472e+1,-3.8777e+1,-3.8393e+1,-3.7724e+1,
     b -3.7170e+1,-3.6920e+1/
      data sbddt/5.7677e+0,4.4487e+0,3.4104e+0,2.7330e+0,2.3648e+0,
     a  1.7906e+0,1.3669e+0,1.1305e+0,1.0311e+0,9.0514e-1,
     b  6.4162e-1,2.2587e-1/
      data scddt/-1.0482e+0,-8.5484e-1,-6.4293e-1,-5.6768e-1,-4.6448e-1,
     a -3.6404e-1,-2.4718e-1,-1.7520e-1,-1.0340e-1,-7.8359e-2,
     b -3.0178e-1,-4.4116e-1/
      data sdddt/9.2980e-2,1.0190e-1,4.4821e-2,9.6440e-2,4.8303e-2,
     a  5.6195e-2,4.2875e-2,6.7100e-2,1.2042e-2,-1.0744e-1,
     b -8.3023e-2,-8.3023e-2/
      data sydd3/-5.1025e+1,-4.7399e+1,-4.4610e+1,-4.2849e+1,-4.1915e+1,
     a -4.0444e+1,-3.9335e+1,-3.8635e+1,-3.8255e+1,-3.7613e+1,
     b -3.7050e+1,-3.6713e+1/
      data sbdd3/5.9212e+0,4.5843e+0,3.5063e+0,2.8101e+0,2.4353e+0,
     a  1.8359e+0,1.3878e+0,1.1290e+0,1.0062e+0,8.7269e-1,
     b  7.1843e-1,4.7227e-1/
      data scdd3/-1.0563e+0,-8.7258e-1,-6.8264e-1,-5.6133e-1,-4.8963e-1,
     a -3.7515e-1,-2.7131e-1,-1.9114e-1,-1.5319e-1,-3.9371e-2,
     b -1.8316e-1,-2.5674e-1/
      data sddd3/8.8351e-2,9.1331e-2,7.2261e-2,6.7000e-2,5.5061e-2,
     a  4.9932e-2,4.7755e-2,3.5458e-2,5.4741e-2,-6.9142e-2,
     b -4.3830e-2,-4.3830e-2/
      data sydt/-4.6653e+1,-4.2783e+1,-3.9727e+1,-3.7763e+1,-3.6756e+1,
     a -3.5396e+1,-3.4759e+1,-3.4644e+1,-3.4703e+1,-3.5004e+1,
     b -3.5414e+1,-3.5714e+1/
      data sbdt/6.2405e+0,4.9588e+0,3.9008e+0,3.1030e+0,2.5325e+0,
     a  1.4069e+0,4.8331e-1,-4.0715e-2,-2.7796e-1,-5.5175e-1,
     b -5.9224e-1,-4.5459e-1/
      data scdt/-9.9417e-1,-8.5512e-1,-6.7115e-1,-7.5456e-1,-8.4459e-1,
     a -7.7947e-1,-5.5289e-1,-3.8353e-1,-2.8159e-1,-1.1343e-1,
     b  5.5020e-2,1.9096e-1/
      data sddt/6.6872e-2,8.8462e-2,-4.9681e-2,-8.4132e-2,3.1317e-2,
     a  1.0895e-1,1.0088e-1,9.5260e-2,8.0875e-2,8.1001e-2,
     b  8.0975e-2,8.0975e-2/
      data syd3/-5.8761e+1,-5.2609e+1,-4.7781e+1,-4.4643e+1,-4.2928e+1,
     a -4.0114e+1,-3.7973e+1,-3.6824e+1,-3.6362e+1,-3.5950e+1,
     b -3.5969e+1,-3.6134e+1/
      data sbd3/9.9788e+0,7.8462e+0,6.1590e+0,5.0941e+0,4.5387e+0,
     a  3.5932e+0,2.5343e+0,1.5788e+0,1.0264e+0,2.2413e-1,
     b -2.2184e-1,-3.3003e-1/
      data scd3/-1.6962e+0,-1.3806e+0,-1.0533e+0,-8.4959e-1,-7.0767e-1,
     a -6.5640e-1,-8.7124e-1,-8.3625e-1,-7.1225e-1,-4.4527e-1,
     b -1.9808e-1,4.7406e-3/
      data sdd3/1.5176e-1,1.5740e-1,1.2136e-1,1.3263e-1,2.4657e-2,
     a -1.0331e-1,2.0847e-2,1.1587e-1,1.2840e-1,1.1886e-1,
     b 1.2081e-1,1.2081e-1/
      data sytt4/-4.9468e+1,-4.6395e+1,-4.4042e+1,-4.2559e+1,-4.1772e+1,
     a -4.0526e+1,-3.9559e+1,-3.8901e+1,-3.8498e+1,-3.7698e+1,
     b -3.7108e+1,-3.6986e+1/
      data sbtt4/5.0308e+0,3.8750e+0,2.9536e+0,2.3667e+0,2.0534e+0,
     a  1.5691e+0,1.2484e+0,1.1267e+0,1.1523e+0,1.0767e+0,
     b  5.5697e-1,-1.6708e-1/
      data sctt4/-9.1696e-1,-7.5071e-1,-5.7839e-1,-4.7046e-1,-4.0790e-1,
     a -2.9087e-1,-1.7169e-1,-4.5795e-2,1.1746e-1,-2.2651e-1,
     b -5.2326e-1,-7.7061e-1/
      data sdtt4/7.9957e-2,8.2858e-2,6.4291e-2,5.8462e-2,5.6282e-2,
     a  5.7310e-2,7.4991e-2,1.5256e-1,-1.6543e-1,-1.4270e-1,
     b -1.4734e-1,-1.4734e-1/
c   beam-Maxwellian rates
      data sydtb/-3.4579e+1,-3.4550e+1,-3.4802e+1,-3.5312e+1,-3.5694e+1,
     a-3.6148e+1/
      data sbdtb/ 1.5964e-1,-1.1150e-1,-6.7900e-1,-1.0682e+0,-1.0205e+0,
     a-1.2003e-1/
      data scdtb/-1.1868e-1,-2.7259e-1,-5.4631e-1,-1.4875e-1, 2.8240e-1,
     a 1.0170e+0/
      data sddtb/-7.4031e-2,-1.3166e-1, 2.3664e-1, 4.0258e-1, 3.5334e-1,
     a 3.5334e-1/
      data syttb/-3.9625e+1,-3.8854e+1,-3.8173e+1,-3.7599e+1,-3.7155e+1,
     a-3.6257e+1,-3.7657e+1/
      data sbttb/ 1.2395e+0, 1.0171e+0, 9.7794e-1, 1.0990e+0, 1.4687e+0,
     a 3.5713e-1,-5.0704e+0/
      data scttb/-2.2858e-1,-9.2391e-2, 3.5933e-2, 1.8025e-1, 8.5545e-1,
     a-2.4595e+0,-5.3724e+0/
      data sdttb/ 6.5507e-2, 6.1724e-2, 8.5901e-2, 6.3044e-1,-1.5945e+0,
     a-1.4011e+0,-1.4011e+0/
      data syd3b/-3.7247e+1,-3.6536e+1,-3.5853e+1,-3.5568e+1,-3.5633e+1,
     a-3.6236e+1/
      data sbd3b/ 8.7888e-1, 1.0864e+0, 8.1017e-1, 1.0844e-1,-4.4450e-1,
     a-1.2665e+0/
      data scd3b/ 3.3735e-1,-3.7896e-2,-3.6070e-1,-8.9239e-1,-6.5647e-1,
     a-5.2961e-1/
      data sdd3b/-1.8049e-1,-1.5527e-1,-3.1648e-1, 2.2028e-1, 6.1021e-2,
     a 6.1021e-2/
c     end of block data
c
      gfun(y)=1.-(log((1.-sqrt(y)+y)/(1.+sqrt(y))**2)/3.+1.15470*
     1 (atan((2.*sqrt(y)-1.)/1.732051)+0.5235987))/y
      tsfun(y)=log(1.+y**1.5)/3.
c
      if ( ifusn ) 200,300,1
c
c  1  dtf=dt*1.0e-6
c     dtfi=1./dtf
c
  1   continue
c
      do 4 jz=lc,le
c zxlme(j)=cloge(j) ?
        zte=tefus(jz)*1000.
        if ( zte .lt. 10.) then
          zxlme(jz)= 23.-log(sqrt(rn(jz)/zte)/zte)
        else
          zxlme(jz)=24.-log(sqrt(rn(jz))/zte)
        endif
   4  continue
c
c  fusion cross section sigma v: S V ab (c)(b), (b=suprathermal 'beam')
c     input ln T_i (keV)
c
      do 9 j=lc,le
        zti = log(tifus(j))
cbate        zti = max ( sx(1), min ( sx(nsv), log(tifus(j)) ) )
        svdt(j) = exp(seval(nsv,zti,sx,sydt,sbdt,scdt,sddt))
        svddt(j)= exp(seval(nsv,zti,sx,syddt,sbddt,scddt,sdddt))
        svdd3(j)= exp(seval(nsv,zti,sx,sydd3,sbdd3,scdd3,sddd3))
        svtt4(j)= exp(seval(nsv,zti,sx,sytt4,sbtt4,sctt4,sdtt4))
        svd3(j)= exp(seval(nsv,zti,sx,syd3,sbd3,scd3,sdd3))
   9  continue
c
c   fast 3He + D
c
      do 7 j=lc,le
        if ( fn3(j) .le. 0. ) then
          svd3b(j)=0.
        else
          zti=6.6666667e-4*e3(j)/fn3(j)
          if ( zti .lt. 100. ) then
            svd3b(j)=0.
          else if ( zti .lt. 500. ) then
            svd3b(j) = exp(-35.694+(500.-zti)*(35.694-41.551
     &        * exp(-0.015464 * sqrt(tifus(j))))/400.)
          else
            zti=log(zti)
            svd3b(j)=exp(seval(6,zti,sxb,syd3b,sbd3b,scd3b,sdd3b))
          endif
        endif
   7  continue
c
      do 8 j=lc,le
        if ( fnt(j) .le. 0. .or. et(j) .le. 0. ) then
          ztemp(j)=0.
          svdtb(j)=0.
          svttb(j)=0.
        else
          ztemp(j)=log(6.666667e-4*et(j)/fnt(j))
        endif
   8  continue
c
c   fast T + D
c
      do 2 j=lc,le
        if ( ztemp(j) .gt. 0. )
     &   svdtb(j)=exp(seval(6,ztemp(j),sxb,sydtb,sbdtb,scdtb,sddtb))
   2  continue
c
c   fast T + T
c
      do 3 j=lc,le
        if ( ztemp(j) .gt. 0. )
     &   svttb(j)=exp(seval(7,ztemp(j),sxb,syttb,sbttb,scttb,sdttb))
   3  continue
c
c  zyc = E_c factors indep of the fusion reaction -- E_c(eV)
c  ztc = tau_s factors indept of the fusion reaction
c  zstarb = [Z]; should have rn,fn *Coul logs if T_e<10* Zf**2(eV)
c
       do 12 j=lc,le
         zstarb(j)=((rnp(j)+fnpd(j)+fnp3(j)+rnd(j)/amd+
     &   (rnt(j)+fnt(j))/amt)+4.*((rn3(j)+fn3(j))/am3+(rn4(j)+
     &   fn43(j)+fn4t(j))/am4))/rn(j)
  12  continue
c
      do 14 j=lc,le
        zyc(j)=14.804e03*tefus(j)*zstarb(j)**0.6666667
        ztc(j)=1.99815e13*tefus(j)**1.5/(rn(j)*zxlme(j))
  14  continue
c
c      slowing down times for particle a
c           T P a
c             |
c             [P=particle , E=energy]
c
c   D + D --> T + p for p
c
       do 20 j=lc,le
  20  ztemp(j)=3.20e6/zyc(j)
      do 21 j=lc,le
  21  fpde(j)=gfun(ztemp(j))
      do 22 j=lc,le
      tepd(j)=0.5*ztc(j)*fpde(j)
  22  tppd(j)=ztc(j)*tsfun(ztemp(j))
c
c        D + He3 --> He4 + p for p
c
      do 30 j=lc,le
  30  ztemp(j)=1.47e7/zyc(j)
      do 31 j=lc,le
  31  fp3e(j)=gfun(ztemp(j))
      do 32 j=lc,le
      tep3(j)=0.5*ztc(j)*fp3e(j)
  32  tpp3(j)=ztc(j)*tsfun(ztemp(j))
c
c        D + D --> T + p for T
c
      do 40 j=lc,le
  40  ztemp(j)=1.01e6/(amt*zyc(j))
      do 41 j=lc,le
  41  fte(j)=gfun(ztemp(j))
      do 42 j=lc,le
      tet(j)=0.5*ztc(j)*amt*fte(j)
  42  tpt(j)=ztc(j)*amt*tsfun(ztemp(j))
c
c        D + D --> He3
c
      do 50 j=lc,le
  50  ztemp(j)=8.2e5/(am3*zyc(j))
      do 51 j=lc,le
  51  f3e(j)=gfun(ztemp(j))
      do 52 j=lc,le
      te3(j)=0.5*ztc(j)*am3/4.*f3e(j)
  52  tp3(j)=ztc(j)*am3/4.*tsfun(ztemp(j))
c
c        D + T --> He4
c
      do 60 j=lc,le
  60  ztemp(j)=3.5e6/(am4*zyc(j))
      do 61 j=lc,le
  61  f4te(j)=gfun(ztemp(j))
      do 62 j=lc,le
      te4t(j)=0.5*ztc(j)*am4/4.*f4te(j)
  62  tp4t(j)=ztc(j)*am4/4.*tsfun(ztemp(j))
c
c       D + He3 --> He4 + p for He4
c
      do 70 j=lc,le
  70  ztemp(j)=3.6e6/(am4*zyc(j))
      do 71 j=lc,le
  71  f43e(j)=gfun(ztemp(j))
      do 72 j=lc,le
      te43(j)=0.5*ztc(j)*am4/4.*f43e(j)
  72  tp43(j)=ztc(j)*am4/4.*tsfun(ztemp(j))
c
c       T + T --> He4 + 2n
c
      do 75 j=lc,le
  75  ztemp(j)=1.85e6/(am4*zyc(j))
      do 76 j=lc,le
  76  f42te(j)=gfun(ztemp(j))
      do 77 j=lc,le
      te42t(j)=0.5*ztc(j)*am4/4.*f42te(j)
  77  tp42t(j)=ztc(j)*am4/4.*tsfun(ztemp(j))
c
c      fusion rate (cm-3/s); (b) = beam-Maxwellian rates for T, 3He
c   FDDT      D + D --> T + p
c   FDD3       D + D --> He3
c   FDT(b)      D + T --> He4
c   FD3(b)     D + He3 --> He4 + p
c   FTT4(b)    T + T --> He4 + 2n
c
      do 80 j=lc,le
       fddt(j) = 0.5*svddt(j)*rnd(j)**2
       fdd3(j) = 0.5*svdd3(j)*rnd(j)**2
       fdt(j)  = svdt(j)*rnd(j)*rnt(j)
       fd3(j)  = svd3(j)*rnd(j)*rn3(j)
  80   ftt4(j) = 0.5*svtt4(j)*rnt(j)**2
c
      if(ifusn.eq.1) go to 84
      do 82 j=lc,le
       fdtb(j) = svdtb(j)*rnd(j)
       fd3b(j) = svd3b(j)*rnd(j)
  82   fttb(j) = svttb(j)*rnt(j)
  84  continue
c
c   fn-(-) is the particle density (cm-3) for suprathermal 'fast' ions
c   fnpd         p from D + D --> T + p
c   fnp3         p from D + He3 --> He4 + p
c   fnt          T from  D + D --> T + p
c   fn3          He3 from D + D --> He3
c   fn4t         He4 from D + T --> He4
c   fn43         He4 from D + He3 --> He4 + p
c   fn42t        He4 from T + T --> He4 + 2n  (average neutron power used)
c
c    fploss = 1 - first orbit loss of 14.7 Mev protons (fraction)
c
      if(ifusn.ge.2) go to 100
c   trial - instantaneous slowing down
      do 86 j=lc,le
        epde(j)  = 3.02e06*fddt(j)*fpde(j)
        ep3e(j)  = 1.47e07*fd3(j)*fp3e(j)
        etre(j)  = 1.01e06*fddt(j)*fte(j)
        e3e(j)   = 8.2e05*fdd3(j)*f3e(j)
        e4te(j)  = 3.5e06*fdt(j)*f4te(j)
        e43e(j)  = 3.6e6*fd3(j)*f43e(j)
        e42te(j) = 1.85e06*ftt4(j)*f42te(j)
        epdi(j)  = 3.02e06*fddt(j)*(1.-fpde(j))
        ep3i(j)  = 1.47e07*fd3(j)*(1.-fp3e(j))
        etri(j)  = 1.01e06*fddt(j)*(1.-fte(j))
        e3i(j)   = 8.2e05*fdd3(j)*(1.-f3e(j))
        e4ti(j)  = 3.5e06*fdt(j)*(1.-f4te(j))
        e43i(j)  = 3.6e6*fd3(j)*(1.-f43e(j))
  86    e42ti(j) = 1.85e06*ftt4(j)*(1.-f42te(j))
c
c   thermal particle sources (#/cm3/s)
c
      do 144 j=lc,le
        drnp(j) = fddt(j)+fd3(j)*fploss
        drnd(j) = -2.*(fddt(j)+fdd3(j))-fdt(j)-fd3(j)
        drnt(j) = fddt(j)-fdt(j)-2.*ftt4(j)
        drn3(j) = fdd3(j)-fd3(j)
 144    drn4(j) = fdt(j)*f4loss+fd3(j)+ftt4(j)
      go to 112
c
c   simple slowing down model
c
c
c   ifusn=cfutz(id3he)=2;  slowing down weigted to t+dt/2
c
 100  dtf2=0.5*dtf
      do 101 j=lc,le
      fnpd(j)=exp(-dtf/tppd(j))*(fnpd0(j)+dtf*fddt(j)*exp(dtf2/
     1  tppd(j)))
      fnp3(j)=exp(-dtf/tpp3(j))*(fnp30(j)+dtf*(fd3(j)+fd3b(j)*fn30(j))
     1 *exp(dtf2/tpp3(j))*fploss)
      fnt(j)=exp(-dtf*(1./tpt(j)+fdtb(j)+fttb(j)))*(fnt0(j)+dtf*
     1 fddt(j)*exp(dtf2*(1./tpt(j)+fdtb(j)+fttb(j))))
      fn3(j)=exp(-dtf*(1./tp3(j)+fd3b(j)))*(fn30(j)+dtf*fdd3(j)*
     1  exp(dtf2*(1./tp3(j)+fd3b(j))))
      fn4t(j)=exp(-dtf/tp4t(j))*(fn4t0(j)+dtf*(fdt(j)+fdtb(j)*fnt0(j))*
     1  f4loss*exp(dtf2/tp4t(j)))
      fn43(j)=exp(-dtf/tp43(j))*(fn430(j)+dtf*(fd3(j)+fd3b(j)*fn30(j))*
     1  exp(dtf2/tp43(j)))
 101  fn42t(j)=exp(-dtf/tp42t(j))*(fn42t0(j)+dtf*(ftt4(j)+fttb(j)*
     1 fnt0(j))*exp(dtf2/tp42t(j)))
c
c   Ea_ = (3/2) n_a T_a = energy density of fast particle a (eV/cm3)
c
      do 104 j=lc,le
      epd(j)=exp(-dtf/tepd(j))*(epd0(j)+dtf*3.02e6*fddt(j)*
     1  exp(dtf2/tepd(j)))
      ep3(j)=exp(-dtf/tep3(j))*(ep30(j)+dtf*1.47e7*(fd3(j)+fd3b(j)*
     1 fn30(j))*exp(dtf2/tep3(j))*fploss)
      et(j)=exp(-dtf*(1./tet(j)+fdtb(j)+fttb(j)))*(et0(j)+dtf*1.01e6*
     1 fddt(j)*exp(dtf2*(1./tet(j)+fdtb(j)+fttb(j))))
      e3(j)=exp(-dtf*(1./te3(j)+fd3b(j)))*(e30(j)+dtf*8.2e5*fdd3(j)*
     1  exp(dtf2*(1./te3(j)+fd3b(j))))
      e4t(j)=exp(-dtf/te4t(j))*(e4t0(j)+dtf*3.5e6*(fdt(j)+fdtb(j)*
     1 fnt0(j))*exp(dtf2/te4t(j))*f4loss)
      e43(j)=exp(-dtf/te43(j))*(e430(j)+dtf*3.6e6*(fd3(j)+fd3b(j)*
     1 fn30(j))*exp(dtf2/te43(j)))
 104  e42t(j)=exp(-dtf/te42t(j))*(e42t0(j)+dtf*1.85e6*(ftt4(j)+
     1 fttb(j)*fnt0(j))*exp(dtf2/te42t(j)))
c
c     et(j)=(et0(j)*(dtfi-0.5*(1./tet(j)+fdtb(j)+fttb(j)))+
c    1 1.01e6*fddt(j))/(dtfi+0.5*(1./tet(j)+fdtb(j)+fttb(j)))
c
c
c   kludge for slowing down time te__ lt dt, when solution decays expon.
c     et(j)=max(et(j),0.)
c   particles not strictly conserved (eg, edge)
      do 108 j=lc,le
        fnpd(j)  = max(1.e-2,fnpd(j))
        fnp3(j)  = max(1.e-2,fnp3(j))
        fnt(j)   = max(1.e-2,fnt(j))
        fn3(j)   = max(1.e-2,fn3(j))
        fn4t(j)  = max(1.e-2,fn4t(j))
        fn43(j)  = max(1.e-2,fn43(j))
 108    fn42t(j) = max(1.e-2,fn42t(j))
      do 109 j=lc,le
        epd(j)  = max(1.e-2,epd(j))
        ep3(j)  = max(1.e-2,ep3(j))
        et(j)   = max(1.e-2,et(j))
        e3(j)   = max(1.e-2,e3(j))
        e4t(j)  = max(1.e-2,e4t(j))
        e43(j)  = max(1.e-2,e43(j))
 109    e42t(j) = max(1.e-2,e42t(j))
c
c   fusion power density sources q1fus, q2fus (W/cm3/microsec)/punit
c   ztemp = total power density to ions from collis. slowing down
c      to species j,  * n_j*Z_j**2*(m_p/m_j)/n_e/[Z]
c   fusion heating to e (i) from each fast particle - epde (epdi)
c
      do 110 j=lc,le
        epde(j)  = 0.5*(epd(j)+epd0(j))*fpde(j)/tepd(j)
        ep3e(j)  = 0.5*(ep3(j)+ep30(j))*fp3e(j)/tep3(j)
        etre(j)  = 0.5*(et(j)+et0(j))*fte(j)/tet(j)
        e3e(j)   = 0.5*(e3(j)+e30(j))*f3e(j)/te3(j)
        e4te(j)  = 0.5*(e4t(j)+e4t0(j))*f4te(j)/te4t(j)
        e43e(j)  = 0.5*(e43(j)+e430(j))*f43e(j)/te43(j)
        e42te(j) = 0.5*(e42t(j)+e42t0(j))*f42te(j)/te42t(j)
        epdi(j)  = 0.5*(epd(j)+epd0(j))*(1.-fpde(j))/tepd(j)
        ep3i(j)  = 0.5*(ep3(j)+ep30(j))*(1.-fp3e(j))/tep3(j)
        etri(j)  = 0.5*(et(j)+et0(j))*(1.-fte(j))/tet(j)
        e3i(j)   = 0.5*(e3(j)+e30(j))*(1.-f3e(j))/te3(j)
        e4ti(j)  = 0.5*(e4t(j)+e4t0(j))*(1.-f4te(j))/te4t(j)
        e43i(j)  = 0.5*(e43(j)+e430(j))*(1.-f43e(j))/te43(j)
 110    e42ti(j) = 0.5*(e42t(j)+e42t0(j))*(1.-f42te(j))/te42t(j)
c
      do 150 j=lc,le
      drnp(j)=0.5*((fnpd(j)+fnpd0(j))/tppd(j)+(fnp3(j)+fnp30(j))/
     1 tpp3(j))
      drnd(j)=-2.*(fddt(j)+fdd3(j))-fdt(j)-fd3(j)-0.5*(fdtb(j)*
     1 (fnt0(j)+fnt(j))+fd3b(j)*(fn30(j)+fn3(j)))
      drnt(j)=0.5*(fnt(j)+fnt0(j))*(1./tpt(j)-fttb(j))-fdt(j)-2.*
     1 ftt4(j)
      drn3(j)=0.5*(fn3(j)+fn30(j))/tp3(j)-fd3(j)
 150  drn4(j)=0.5*((fn4t(j)+fn4t0(j))/tp4t(j)+(fn43(j)+
     1 fn430(j))/tp43(j)+(fn42t(j)+fn42t0(j))/tp42t(j))
c
c   total fast particle densities and thermal energy (standard units)
c
      do 159 j=lc,le
        d3fast(j)=(epd(j)+ep3(j)+et(j)+e3(j)+e4t(j)+e43(j)+e42t(j))*
     &    1.6022e-12
        rh1fst(2,j)=fnpd(j)+fnp3(j)+fnt(j)
 159    rh2fst(2,j)=fn3(j)+fn4t(j)+fn43(j)+fn42t(j)
c
c
 112  continue
c
c
c   volume integrated fusion powers from slowing down
c
      do 116 j=lc,le
 116  ztemp(j)=epdi(j)+ep3i(j)+etri(j)+e3i(j)+e4ti(j)+e43i(j)+e42ti(j)
c     pepd=volint(lc,le,epde,zdum)
c
      do 120 j=lc,le
 120  qefus(j)=epde(j)+ep3e(j)+etre(j)+e3e(j)+e4te(j)+
     1 e43e(j)+e42te(j)
c
      do 130 j=lc,le
      q1fus(j)=(rnp(j)+rnd(j)/amd+rnt(j)/amt)*zii1**2/(rn(j)*
     1 zstarb(j))*ztemp(j)
 130  q2fus(j)=(rn3(j)/am3+rn4(j)/am4)*zii2**2/
     1 (rn(j)*zstarb(j))*ztemp(j)
c
c  power density loss from burnup of thermal, fusing ions
c
      do 140 j=lc,le
      q1lfus(j)=1.5*tifus(j)*(2.*(fddt(j)+fdd3(j)+fdt(j)+ftt4(j))+
     1 fd3(j)+0.5*(fdtb(j)+fttb(j))*(fnt(j)+fnt0(j))+fd3b(j)*0.5*
     2 (fn3(j)+fn30(j)))
 140  q2lfus(j)=1.5*tifus(j)*fd3(j)
c
c   standard units (ergs/cm3/s)
c
      do 142 j=lc,le
        qefus(j)  = qefus(j)*1.6022e-12
        q1fus(j)  = q1fus(j)*1.6022e-12
        q2fus(j)  = q2fus(j)*1.6022e-12
        q1lfus(j) = q1lfus(j)*1.6022e-12
 142    q2lfus(j) = q2lfus(j)*1.6022e-12
c
c
c  energetic neutron production (ergs/cm3/s)
c
      wntn(1)=0.
      do 170 j=lc,le
 170  wntn(j)=((fdt(j)+fdtb(j)*0.5*(fnt(j)+fnt0(j)))*1.4e7+fdd3(j)
     1 *2.45e6+(ftt4(j)+fttb(j)*0.5*(fnt(j)+fnt0(j)))*9.45e6)*
     2  1.6022e-12
c
c   electrons lost to maintain charge neutrality for first orbit loss
c
      do 180 j=lc,le
 180  pfols(j)=(fd3(j)+fd3b(j)*0.5*(fn30(j)+fn3(j)))*(1.-fploss)
     1 +2.*(fdt(j)+fdtb(j)*0.5*(fnt0(j)+fnt(j)))*(1.-f4loss)
c
      return
c
c   initialize fast particles = 0 at t=0
c      ifusn=-1
c
 200  do 210 j=1,385
        fnpd(j)=0.
 210    epd(j)=0.
c
      do 220 j=1,n9
        fdtb(j)=0.
        fd3b(j)=0.
 220    fttb(j)=0.
c
c   update time-dependent arrays after successful time step
c      ifusn=0
c
 300  continue
      do 310 j=lc,le+1
        fnpd0(j)  = fnpd(j)
        fnp30(j)  = fnp3(j)
        fnt0(j)   = fnt(j)
        fn30(j)   = fn3(j)
        fn4t0(j)  = fn4t(j)
        fn430(j)  = fn43(j)
        fn42t0(j) = fn42t(j)
        epd0(j)   = epd(j)
        ep30(j)   = ep3(j)
        et0(j)    = et(j)
        e30(j)    = e3(j)
        e4t0(j)   = e4t(j)
        e430(j)   = e43(j)
 310    e42t0(j)  = e42t(j)
c
      return
      end
c@voling
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      function volint(jstart,jend,zf,zpf)
c
c   volume integral for 2D geometry; dvols =vol element in mhdbal
c     =vols(j+1,1)-vols(j,1);  integrates zf(j), zone-centers;
c      zpf(j) = volume-
c      integrated zf with zpf(1)=0 on axis and zpf(jend+1)=total
c      on zone boundaries
c
      include 'cparm.m'
      include 'commhd.m'
      include 'comhd3.m'
cahk960603
      dimension zf(*),zpf(*)
c
      zpf(1)=0.
      do 10 jz=jstart,jend
  10  zpf(jz+1)=zpf(jz)+zf(jz)*dvols(jz)
      volint=zpf(jend+1)
      return
      end
c@igntst
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine igntst
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cd3he.m'
      include 'cpower.m'
      dimension zdum(mj)
c
      if(nstep.le.1) return
      poht0=poht
      pradt0=pradt
      pfust0=pfust
      pauxt0=pauxt
      plost0=plosst
      pcond0=pcond
      jza=mzones
      poht=geohms(jza)
      pradt=gesyns(jza)+gebrs(jza)+geirs(jza)+gesrs(jza)
      pfust=ged3fs(jza)+gid3fs(jza)+gialfs(jza)+gealfs(jza)
      pauxt=geauxs(jza)+giauxs(jza)+gebems(jza)+gibems(jza)+
     1 geecrs(jza)+giecrs(jza)+geicrs(jza)+giicrs(jza)
      plosst=gebrs(jza)+gesyns(jza)+geirs(jza)+geions(jza)-
     1 giions(jza)-gichxs(jza)
c   convert to watts
      poht=poht*usep
      pfust=pfust*usep
      pradt=pradt*usep
      pauxt=pauxt*usep
      plosst=plosst*usep
c
      pe4t=volint(lcentr,ledge,e4te,zdum)*1.6022e-19
      pi4t=volint(lcentr,ledge,e4ti,zdum)*1.6022e-19
      pfusdt=pe4t+pi4t
      pfusd3=(ged3fs(jza)+gid3fs(jza))*usep-pfusdt
      psyn=gesyns(jza)*usep
      pbrem=gebrs(jza)*usep
      pirad=geirs(jza)*usep
      pntn=volint(lcentr,ledge,wntn,zdum)*usep
c
c   conduction and convection in Watts
      z0=avi(jza,6,1)/(avi(jza,5,1)*dxboui(jza)*uisl)
      zgrdti=(tis(2,jza)-tis(2,jza-1))*z0
      pconde=condei(jza)*uiep
      pcondi=condii(jza)*uiep
      pconve=cvctei(jza)*uiep
      pconvi=cvctii(jza)*uiep
        zineut=-ditins(jza)*zgrdti*avi(jza,3,1)*avi(jza,5,1)*
     1   uisl**2*usep
c
      pcond=pconde+pcondi+pconve+pconvi
      if(nstep.le.5) return
c
      if(pfust.lt.poht.or.pfust0.ge.poht0) go to 150
      write(nprint,101)
 101     format(' ign-preburn')
      call sprint(nprint,2)
      call d3prin(nprint)
 150  if(pfust.lt.poht+pauxt.or.pfust0.ge.poht0+pauxt0) go to 200
      write(nprint,151)
 151  format(' ign-preburn with aux heating')
      call sprint(nprint,2)
      call d3prin(nprint)
 200  if(pfust.lt.poht+pradt.or.pfust0.ge.poht0+pradt0) go to 250
      write(nprint,201)
 201    format(' ign-burn - fusion > OH + rad')
      call sprint(nprint,2)
      call d3prin(nprint)
 250  if(pfust.lt.poht+pradt+pauxt.or.pfust0.ge.poht0+pauxt0+
     1 pradt0) go to 300
      write(nprint,251)
 251    format(' ign-burn with aux heating')
      call sprint(nprint,2)
      call d3prin(nprint)
 300    if(pfust.lt.pcond+plosst.or.pfust0.ge.pcond0+plost0)
     1 go to 400
      write(nprint,301)
 301    format(' ignition - fusion > conduction + convection + rad')
      call sprint(nprint,2)
      call d3prin(nprint)
 400    return
      end
c@reorder
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c  reorder species p,d,t for output of ncflux if d-3he fusion
        subroutine reorder(ih,mz,w)
        dimension w(6,6,55), wd(6,6,55)
        do 100 jz=1,mz
        wd(1,1,jz)=w(2,2,jz)
        wd(1,2,jz)=w(2,3,jz)
        wd(1,3,jz)=w(2,1,jz)
        wd(2,1,jz)=w(3,2,jz)
        wd(2,2,jz)=w(3,3,jz)
        wd(2,3,jz)=w(3,1,jz)
        wd(3,1,jz)=w(1,2,jz)
        wd(3,2,jz)=w(1,3,jz)
        wd(3,3,jz)=w(1,1,jz)
        do 90 ji=4,6
        wd(1,ji,jz)=w(2,ji,jz)
        wd(2,ji,jz)=w(3,ji,jz)
  90    wd(3,ji,jz)=w(1,ji,jz)
 100   continue
        do 200 jz=1,mz
        do 200 j=1,ih
        do 200 k=1,ih
  200   w(j,k,jz)=wd(j,k,jz)
        return
        end
c***********************************************************************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bussac
c   compute Bussac's beta-poloidal for stability; also <B_p**2>
c     beta_pB=(8*pi/<B_p**2>_o)*(<<p>>_o-p_o), <<.>> = volume average,
c      subscript o means q=1 surface
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cd3he.m'
      include 'comhd3.m'
c
      dimension zdum(mj),zprsth(mj),zebms(mj)
c   <b_p^2> = bpolsq(j) (T**2)
      do 10 jz=lcentr,mzones
      bpolsq(jz)=avi(jz,7,1)*(r0ref*bpoli(jz)*avi(jz,2,1))**2
      zprsth(jz)=rhoels(2,jz)*tes(2,jz)+rhoins(2,jz)*tis(2,jz)
  10  zebms(jz)=hebems(jz)*rhobis(2,jz)
c  beta-poloidal-Bussac thermal, hot-isotropic, hot-anisotropic (beams)
      buscth=0.
      buschi=0.0
      buscha=0.0
c   find largest q=1 surface
      if(q(mzones).lt.1.0) go to 22
      do 20 jz=mzones,lcentr+1,-1
      if(q(jz).le.1.0) go to 30
  20  continue
  22  buscth=-1.0
      buschi=-1.0
      buscha=-1.0
      go to 100
  30  jzq1=jz
      jzq2=jzq1+1
      zfracq=(1.-q(jzq1))/(q(jzq2)-q(jzq1))
c   rho(j) = avi(j,1,it)
c   interpolation, arbitrarily, is based on toroidal flux radius
      rhoqz1=avi(jzq1,1,1)
      rhoqz2=avi(jzq2,1,1)
      rhoq1=avi(jzq1,1,1)+(avi(jzq2,1,1)-avi(jzq1,1,1))*zfracq
      rhzq1=ahalfs(jzq1,1)+(ahalfs(jzq2,1)-ahalfs(jzq1,1))*zfracq
      prs0th=zprsth(jzq1)+(zprsth(jzq2)-zprsth(jzq1))*zfracq
      d3f0=d3fast(jzq1)+(d3fast(jzq2)-d3fast(jzq1))*zfracq
      ebms0=zebms(jzq1)+(zebms(jzq2)-zebms(jzq1))*zfracq
      zvolq1=(avi(jzq1,12,1)+(avi(jzq2,12,1)-avi(jzq1,12,1))*
     1 zfracq)*uisl**3
      vprsth=volint(lcentr,jzq1,zprsth,zdum)
      vprsth=(vprsth+0.5*(zprsth(jzq1)+prs0th)*(vols(jzq2,1)-
     1 vols(jzq1,1))*
     2  (rhoq1**2-rhoqz1**2)/(rhoqz2**2-rhoqz1**2))/zvolq1
      vprshi=volint(lcentr,jzq1,d3fast,zdum)
      vprshi=(vprshi+0.5*(d3fast(jzq1)+d3f0)*(vols(jzq2,1)-
     1  vols(jzq1,1))*
     2  (rhoq1**2-rhoqz1**2)/(rhoqz2**2-rhoqz1**2))/zvolq1
      vprsha=volint(lcentr,jzq1,zebms,zdum)
      vprsha=(vprsha+0.5*(zebms(jzq1)+ebms0)*(vols(jzq2,1)-
     1 vols(jzq1,1))*
     2  (rhoq1**2-rhoqz1**2)/(rhoqz2**2-rhoqz1**2))/zvolq1
      bpsq0=(bpolsq(jzq1)+(bpolsq(jzq2)-bpolsq(jzq1))*zfracq)*uisb**2
c   bussac beta-p
      buscth=8.0*fcpi*(vprsth-prs0th)/bpsq0
      buschi=8.0*fcpi*(vprshi-d3f0)/bpsq0
      buscha=8.0*fcpi*(vprsha-ebms0)/bpsq0
  100  return
      end
c***********************************************************************
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c@cmg
c   coppi-mazzucato-gruber chi-e
c      generalized to 2D geometry
      subroutine cmg(icmg)
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'comhd3.m'
      include 'cd3he.m'
      include 'cmg.m'
      dimension znue(n9),zbetap(n9),zdpdi(n9)
c   computed at zone boundaries
c   units cm2/s; alcmg = 0.25, epscmg = 0.2 standard
c     alpcmg = 0.7; gcmg=4.5 H-mode, 9.0 L-mode
      jza=mzones
c
c   toroidal current curtor in (amps) from torcur
c   tes(.,.,) *useh in keV, lengths in (m) from deqbald
      zeps = alcmg*2.68e09*avi(jza,12,1)**2/(avi(jza,13,1)**1.5
     1 *useh)
      do 10 jz=lcentr,mzones
      zstara(jz)=0.
      zbetap(jz)=0.
      do 4 ih=1,mhyd
    4    zstara(jz)=zstara(jz)+rhohs(ih,1,jz)/aspec(ih)
      if(mimp.eq.0) go to 8
      do 6 ii=1,mimp
       ji=ii+lhydn
    6   zstara(jz)=zstara(jz)+c2mean(ii,1,jz)*rhois(ii,1,jz)/aspec(ji)
    8   zstara(jz)=zstara(jz)/rhoels(1,jz)
        zbetap(jz)=(zbetap(jz)+rhoels(2,jz)*tes(2,jz))+
     1  d3fast(jz)
   10  continue
      if(geohms(jza).eq.0.) then
      zf=1.
      else
      zf=ged3fs(jza)+geauxs(jza)+gealfs(jza)+geecrs(jza)+geicrs(jza)
     1  + gid3fs(jza)+giauxs(jza)+gialfs(jza)+giecrs(jza)+giicrs(jza)
      zf=zf/(geohms(jza)+zf)
      endif
c
      if(cfutz(492).lt.epslon) go to 38
      zwtot=volint(lcentr,ledge,zbetap,znue)
      zwtot=zwtot*4.*twopi**5/vols(mzones,1)
c   <1/R**2> at boundary
        zint0 = (xzoni(mzones) - xbouni(mzones)) / dxboui(mzones)
        zint1 = 1.0 - zint0
        zr2=avi(mzones-1,10,1)*zint0 + avi(mzones,10,1)*zint1
c   factor 1. + gcmg*betap_e*((Pf+Pa)/(Pf+Pa+Poh))**2
      zbetpe=zwtot*(q(jza)/(avi(jza,3,1)*avi(jza,8,1)*
     1  zr2*uisb))**2/avi(jza,7,1)
      zgcmg=1.+gcmg*zf**2*zbetpe
c      special xcmg(lcentr) - contains curtor=0/(dV/d xi=0)
      do 20 jz=lcentr+1,mzones
  20  xcmg0(jz)=zeps*curtor(jz)*(cloges(1,jz)*zstara(jz))**0.4*
     1 zgcmg/(rhoels(1,jz)**0.8*tes(1,jz)*avi(jz,3,1)**2*avi(jz,6,1))
c   tohcmg (keV) is cutoff temp for purely ohmic chi-e (5 keV) for
c      Ignitor-like)
      if(tes(1,lcentr)*useh.le.tohcmg) go to 22
      do 21 jz=lcentr+1,mzones
  21  xcmg0(jz)=xcmg0(jz)*useh*tes(1,lcentr)/tohcmg
      xcmg0(lcentr)=xcmg0(lcentr+1)
  22  continue
c
c   ubiquitous mode (B. Coppi, Comm. Pl. Phys. Cont. Fusion 12 (1989)
c      312).
c
  38  if(cfutz(493).le.epslon) go to 58
c   coh=cfutz(492)=1, cub=cfutz(493)=1.4, cfutz(494)=1/4-1/2 (3/8),
c     alphaub=cfutz(495)=1-3, cfutz(496)=0-? where chi_i=f_i*xub+chi_neo
c   tohcmg=1000. (large), gcmg=0. with ubiquitous mode
      do 50 j=lcentr+1,mzones-1
      zdpdi(j)=(rhoels(1,j)*(tes(2,j)-tes(2,j-1))+rhoins(1,j)*
     1 (tis(2,j)-
     1 tis(2,j-1))+tes(1,j)*(rhoels(2,j)-rhoels(2,j-1))+tis(1,j)*(
     2 rhoins(2,j)-rhoins(2,j-1)))/(curtor(j+1)-curtor(j-1))
      zdpdi(j)=zdpdi(j)*(fltors(j+1,1)-fltors(j-1,1))/(fltors(j,2)-
     1 fltors(j-1,2))
      zai=0.
      do 52 ih=1,mhyd
  52  zai=zai+rhohs(ih,1,j)*aspec(ih)
      do 53 ii=1,mimp
  53  zai=zai+rhois(ii,1,j)*aspec(ii+lhydn)
      zai=zai/rhoins(1,j)
  50  xub0(j)=5.265e25*(ahalfs(jza,1)*zdpdi(j)/
     1 rhoels(1,j))**2*sqrt(ahalfs(jza,1)/rgcs(j,1))*
     2 (avi(j,12,1)/avi(jza,12,1))**cfutz(494)*
     1 xzeff(1,j)**0.5/(q(j)*tes(1,j)**0.5*rgcs(j,1)*zai)
c   xub(mzones) - use current derivative curtor(j)-curtor(j-1)
      j=mzones
      zdpdi(j)=(rhoels(1,j)*(tes(2,j)-tes(2,j-1))+rhoins(1,j)*
     1 (tis(2,j)-
     1 tis(2,j-1))+tes(1,j)*(rhoels(2,j)-rhoels(2,j-1))+tis(1,j)*(
     2 rhoins(2,j)-rhoins(2,j-1)))/(curtor(j)-curtor(j-1))
      zdpdi(j)=zdpdi(j)*(fltors(j,1)-fltors(j-1,1))/(fltors(j,2)-
     1 fltors(j-1,2))
      zai=0.
      do 54 ih=1,mhyd
  54  zai=zai+rhohs(ih,1,j)*aspec(ih)
      do 55 ii=1,mimp
  55  zai=zai+rhois(ii,1,j)*aspec(ii+lhydn)
      zai=zai/rhoins(1,j)
      xub0(j)=5.265e25*(ahalfs(jza,1)*zdpdi(j)/
     1 rhoels(1,j))**2*sqrt(ahalfs(jza,1)/rgcs(j,1))*
     2 (avi(j,12,1)/avi(jza,12,1))**cfutz(494)*
     1 xzeff(1,j)**0.5/(q(j)*tes(1,j)**0.5*rgcs(j,1)*zai)
c
      xub0(lcentr)=xub0(lcentr+1)
c
  58  do 60 j=lcentr,mzones
  60  xub(j)=cfutz(493)*xub0(j)*zf**cfutz(495)
      do 62 j=lcentr,mzones
  62  xcmg(j)=cfutz(492)*xcmg0(j)+xub(j)
c
c   particle transport
c
      do 70 jz=lcentr,mzones
  70  dxcmg(jz)=epscmg*xcmg(jz)
c   inward pinch vxcmg; note vxcmg<0 inward
      do 72 jz=lcentr,mzones
  72  znue(jz)=1.2589e-37*rgcs(jz,1)*q(jz)*xzeff(1,jz)*
     1 cloges(1,jz)*rhoels(1,jz)/tes(1,jz)**2
c     zalp=0.1*alpcmg/vols(mzones,1)
      zalp=alpcmg/vols(mzones,1)*uisl**2
      do 80 jz=lcentr,mzones
  80    vxcmg(jz)=-zalp*dxcmg(jz)*(1.+4.1*sqrt(znue(jz)))*avi(jz,3,1)/
     1  avi(jz,2,1)
         vxcmg(mzones)=0.
      return
      end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c@auxanl(iaux)
c   analytic forms for heatin and particle source
c   les  mar-89
c
      subroutine auxanl(iaux)
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      dimension zdum(mj),zdum1(mj)
c
c   fuel only inside scrapeoff
c     isep = cfutz(120)-1.+epslon
      isep=ledge
      if(lpaux.lt.epslon) go to 500
c   auxiliary heating with given profile (ergs/cm3/s)
c      pauxe, etc given as total MW integrated over plasma
c   note vols(2,1)=0, vols(2,2).gt.0
      zpaux=-apaux/vols(isep+1,1)
c   compute volume integral each time step
      do 5 jz=lcentr,isep
   5  zdum(jz)=exp(zpaux*vols(jz,2))
      zvaux=volint(lcentr,isep,zdum,zdum1)
      do 6 j=1,19
      if(tai.ge.atim(j).and.tai.lt.atim(j+1)) go to 8
   6  continue
      go to 500
  8   zpxe1=pauxe(j)*1.0e13/zvaux
      zpxi1=pauxi(j)*1.0e13/zvaux
        do 20 jz=lcentr,isep
        weauxs(jz)= zpxe1*exp(zpaux*vols(jz,2))
  20    wiauxs(jz)= zpxi1*exp(zpaux*vols(jz,2))
c
c   replacement fuelling for d-3he fusion
 500    znd=0.
      znt=0.
      zn3=0.
      zsaux=-spaux/vols(isep+1,1)
      do 505 jz=lcentr,isep
 505  zdum(jz)=exp(zsaux*vols(jz,2))
      zvaux=volint(lcentr,isep,zdum,zdum1)
c
c
      if(lrepld.eq.0) go to 520
c   replace fusion-reduced deuterium (no/cm3/s)
        do 510 jz=lcentr,ledge
 510  zdum(jz)=shd3fs(ldeut,jz)
      znd=volint(lcentr,ledge,zdum,zdum1)
      znd=-min(0.,znd)
c    replace fusion-reduced tritium
 520  if(lreplt.eq.0) go to 550
         do 522 jz=lcentr,ledge
 522  zdum(jz)=shd3fs(ltrit,jz)
      znt=volint(lcentr,ledge,zdum,zdum1)
      znt=-min(0.,znt)
c   replace fusion-reduced 3He
 550    if(lrepl3.eq.0) go to 600
      do 560 jz=lcentr,ledge
 560  zdum(jz)=sid3fs(lhe3-lhydn,jz)
      zn3=volint(lcentr,ledge,zdum,zdum1)
      zn3=-min(0.,zn3)
 600  continue
c   sp.. is additional fuelling (volume averaged)
      znewd = znd/zvaux
      znewt = znt/zvaux
      znew3 = zn3/zvaux
      znewp=0.
      znew4=0.
      znewi4=0.
      if(lspaux.eq.0) go to 700
      do 602 j=1,19
      if(tai.ge.stim(j).and.tai.lt.stim(j+1)) go to 605
 602  continue
      go to 700
 605  znewd=znewd+spd(j)*vols(isep+1,1)/zvaux
      znewt=znewt+spt(j)*vols(isep+1,1)/zvaux
      znewp=znewp+spp(j)*vols(isep+1,1)/zvaux
      znew3=znew3+sp3(j)*vols(isep+1,1)/zvaux
      znew4=znew4+sp4(j)*vols(isep+1,1)/zvaux
      znewi4=znewi4+spimp4(j)*vols(isep+1,1)/zvaux
 700    do 720 jz=lcentr,isep
      shfuel(ldeut,jz)=znewd*exp(zsaux*vols(jz,2))
      shfuel(ltrit,jz)=znewt*exp(zsaux*vols(jz,2))
      sifuel(lprotn-lhydn,jz)=znewp*exp(zsaux*vols(jz,2))
      sifuel(lhe3-lhydn,jz)=znew3*exp(zsaux*vols(jz,2))
      sifuel(lalpha-lhydn,jz)=znew4*exp(zsaux*vols(jz,2))
 720  sifuel(mimp,jz)=znewi4*exp(zsaux*vols(jz,2))
      return
      end
c@d3prin
c  rgb 27-sep-92 removed r29 from write statement (not defined)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine d3prin(iunit)
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'comhd3.m'
      include 'cd3he.m'
      include 'cpower.m'
      dimension zdu(mj),zduma(mj),idum(mj)
c
c   isep and zbeta as from mprint
      isep=mzones
      if(nadump(1).gt.lcentr) isep=nadump(1)
c   short printout has total fast particles, total fusion heating due
c     to each species, neutron production
c
c   convert erg -> MJ
      zusep=usep*1.0e-06
      zp1l=volint(lcentr,ledge,q1lfus,zdu)*zusep
      zp2l=volint(lcentr,ledge,q2lfus,zdu)*zusep
      zp1=volint(lcentr,ledge,q1fus,zdu)*zusep
      zp2=volint(lcentr,ledge,q2fus,zdu)*zusep
      zpe=volint(lcentr,ledge,qefus,zdu)*zusep
      zeftot=volint(lcentr,ledge,d3fast,zdu)*zusep
      zwntn=volint(lcentr,ledge,wntn,zdu)*zusep
      write(iunit,2003) zpe,zp1,zp2,zp1l,zp2l
 2003   format('fus pwr MW: e=',1pe12.4,' hy=',1pe12.4,' he=',
     1   1pe12.4,'hyloss=',1pe12.4,' he l=',1pe12.4)
        write(iunit,2004) zwntn,zeftot
 2004   format('neutron pwr=',1pe12.4,'fast part ke=',1pe12.4)
      zdum1=geohms(isep)*zusep
      zdum2=ged3fs(isep)*zusep
      zdum3=geauxs(isep)*zusep
      zdum4=gebems(isep)*zusep
      zdum5=geecrs(isep)*zusep
      zdum6=geicrs(isep)*zusep
        zpwr=zdum1+zdum2+zdum3+zdum4+zdum5+zdum6
        write(iunit,400) zdum1,zdum2,zdum3,zdum4,zdum5,zdum6
        zdum1=0.
        zdum2=gid3fs(isep)*zusep
        zdum3=giauxs(isep)*zusep
        zdum4=gibems(isep)*zusep
        zdum5=giecrs(isep)*zusep
        zdum6=giicrs(isep)*zusep
        zpwr=zpwr+zdum2+zdum3+zdum4+zdum5+zdum6
        write(iunit,401) zdum1,zdum2,zdum3,zdum4,zdum5,zdum6
 400   format('poh (MW)    d3fus       aux         beam        ech',
     1 '         ich',/,6(1pe12.4))
 401   format(6(1pe12.4))
        zdum1=gebrs(isep)*zusep
        zdum2=gesyns(isep)*zusep
        zdum3=geions(isep)*zusep
        zdum4=geirs(isep)*zusep
        zdum5=-giions(isep)*zusep
        zdum6=-gichxs(isep)*zusep
        zls=zdum1+zdum2+zdum3+zdum4+zdum5+zdum6
        write(iunit,402) zdum1,zdum2,zdum3,zdum4,zdum5,zdum6
 402  format('      brem        syn        ioniz       ir',
     1 '          ioniz-i     cx-i',/,'brem=',6(1pe12.4))
c   print W-dot losses due to electrons corresponding to ash removal
c      and burnup of thermal ions in fusion reactions
      zdum1=zpwr-zls
      zdum5=volint(lcentr,ledge,wid3fl,zdu)*zusep
      zdum6=geashs(isep)*zusep
 4021 format('dw/dt loss - tot',1pe12.4,'e-ash removal',1pe12.4,
     1 ' fusion burnup',1pe12.4)
c   print conduction losses, from igntst
c     convert W-> MW
      zcond=pcond*1.e-06
      zconde=pconde*1.e-06
      zcondi=pcondi*1.e-06
      zconve=pconve*1.e-06
      zconvi=pconvi*1.e-06
      write(iunit,4020) zcond,zconde,zcondi,zconve,zconvi
 4020  format('cond+conv',1pe12.4,' cond e,i',1p2e12.4,' conv e,i',
     1 1p2e12.4)
c   convert eV -> MJ
      pepd=volint(lcentr,ledge,epde,zdu)*1.6022e-25
      pep3=volint(lcentr,ledge,ep3e,zdu)*1.6022e-25
      pet=volint(lcentr,ledge,etre,zdu)*1.6022e-25
      pe3=volint(lcentr,ledge,e3e,zdu)*1.6022e-25
      pe4t=volint(lcentr,ledge,e4te,zdu)*1.6022e-25
      pe43=volint(lcentr,ledge,e43e,zdu)*1.6022e-25
      pe42t=volint(lcentr,ledge,e42te,zdu)*1.6022e-25
      pipd=volint(lcentr,ledge,epdi,zdu)*1.6022e-25
      pip3=volint(lcentr,ledge,ep3i,zdu)*1.6022e-25
      pit=volint(lcentr,ledge,etri,zdu)*1.6022e-25
      pi3=volint(lcentr,ledge,e3i,zdu)*1.6022e-25
      pi4t=volint(lcentr,ledge,e4ti,zdu)*1.6022e-25
      pi43=volint(lcentr,ledge,e43i,zdu)*1.6022e-25
      pi42t=volint(lcentr,ledge,e42ti,zdu)*1.6022e-25
      write(iunit,2000)
      write(iunit,2001) pepd,pep3,pet,pe3,pe4t,pe43,pe42t
      write(iunit,2001) pipd,pip3,pit,pi3,pi4t,pi43,pi42t
 2000   format('fusion heating (MW): pd,p3,t,3,4dt,4d3,4tt, e/i')
 2001   format(1p7e12.4)
      zdum1=(ged3fs(isep)+gid3fs(isep))*zusep
      if(abs(zdum1).gt.epslon) zdum1=(pe4t+pi4t)/zdum1
        write(iunit,403) zdum1,zp1l,zp2l
 403    format('dt/d3 power=',1pe12.4,'hy fus loss (MW)',1pe12.4,
     1  'he fus ls',1pe12.4)
        zdum1=rhois(lprotn-lhydn,1,lcentr)*used
        zdum2=rhohs(ldeut,1,lcentr)*used
        zdum3=rhohs(ltrit,1,lcentr)*used
        zdum4=rhois(lhe3-lhydn,1,lcentr)*used
        zdum5=rhois(lalpha-lhydn,1,lcentr)*used
        zdum6=rhois(4,1,lcentr)*used
        write(iunit,405) zdum1,zdum2,zdum3,zdum4,zdum5,zdum6
 405   format('ni0',6(1pe12.4))
c   beta-poloidal Bussac
      call bussac
        zdum1=avi(lcentr,8,1)*uieb/avi(lcentr,20,1)
        zdum2=avi(lcentr,20,1)*uiel
         zbp=bpoli(isep)
      zbpsq=sqrt(bpolsq(isep))
      zdum4=5.*uieb*uiel*(avi(isep,8,1)-avi(lcentr,8,1))
        write(iunit,406) zdum1,zdum2,zbp,zbpsq,zdum4
 406  format('bt',1pe11.4,'r',1pe11.4,' bpbldr',1pe11.4,' bpsq',
     1 1pe11.4,' Ipol(kA)',1pe11.4)
      zdum1=volint(lcentr,ledge,ajboot,zdu)*1.0e-03*uiej
      zdum2=volint(lcentr,ledge,ajtpbi,zdu)*1.0e-03*uiej
      do 219 j=lcentr,ledge
 219  zduma(j)=ajbs(j)*usij*0.5*(avi(j+1,11,1)+
     1 avi(j,11,1))/uiel
      zdum3=volint(lcentr,ledge,zduma,zdu)*1.0e-03*uiej
      do 220 j=lcentr,ledge
 220    zduma(j)=ajtori(j,2)
      zdum4=volint(lcentr,ledge,zduma,zdu)*1.0e-09/twopi
      write(iunit,407) zdum1,zdum2,zdum3,zdum4
 407  format('tor curr (kA) - bs',1pe12.4,'trap p',1pe12.4,'beam',
     1 1pe12.4,'tot',1pe12.4)
c
c   central betas - thermal and fast particle
      zprs0=rhoels(2,lcentr)*tes(2,lcentr)+rhoins(2,lcentr)*
     1  tis(2,lcentr)
      zprs0h=d3fast(lcentr)+hebems(lcentr)*rhobis(2,lcentr)
      if(zbpsq.lt.epslon) then
        zdum1=0.
        zdum2=0.
      else
      zdum1=4.*twopi*(zprs0+zprs0h)/(bpolsq(isep)*uisb**2)
      zdum2=4.*twopi*zprs0h/(bpolsq(isep)*uisb**2)
      endif
      zbeta0=4.*twopi*(zprs0+zprs0h)
     1 /(avi(lcentr,8,1)*uisl*uisb/rgcs(1,lcentr))**2
       zbetf0=zprs0h/(avi(lcentr,8,1)*uisl*uisb/rgcs(1,lcentr))**2
      zbetf0=zbetf0*4.*twopi
      write(iunit,411) vloopi(isep-1,1),vloopi(lcentr,1),zbeta0,
     1  zbetf0,zdum1,zdum2
 411  format('va,v0',1p2e11.4,' bta0',1pe11.4,' fst',1pe11.4,
     1  ' btpa0,f',1p2e11.4)
c   total volume integrated betas
      zbeta = max (epslon, cgbeta(1) * betate(isep)
     &          + cgbeta(2) * betati(isep)
     &          + cgbeta(3) * betatb(isep)
     &          +             betatd(isep)
     &          + cgbeta(4) * betata(isep) ) * 100.   ! vol ave beta %
      zdum1=zbeta
      zdum2=betatb(isep)+betatd(isep)+betata(isep)
      if(zbpsq.lt.epslon) then
        zdum3=0.
        zdum4=0.
      else
        zdum3=0.01*zbeta*bzi**2/bpolsq(isep)
        zdum4=0.01*(betatb(isep)+betatd(isep)+betata(isep))*
     1 bzi**2/bpolsq(isep)
      endif
      write(iunit,4112) zdum1,zdum2,zdum3,zdum4
 4112  format('total betator(%)',1pe11.4,' fst t',1pe11.4,' pol',
     1  1pe11.4,' fst p',1pe11.4)
      write(iunit,4111) buscth,buschi,buscha,rhzq1,rhoq1
 4111  format('btap-bussac - therml',1pe11.4,' fusion',1pe11.4,
     1 ' beam',1pe11.4,' rhz,rho q=1', 1p2e11.4)
      write(iunit,4113) prs0th,vprsth
 4113  format('q=1 p0',1pe12.4, ' vol-ave p',1pe12.4,' erg/cm3')
        zdum2=avi(isep,13,1)*uiel**2
        zdum3=avi(isep,12,1)*uiel**3
      write(iunit,412) zdum2,zdum3
 412  format('area (cm2)', 1pe12.4,'vol',1pe12.4)
        if(tis(1,lcentr)*useh.le.29.) go to 500
        do 240 j=lcentr+1,ledge
        if(tis(2,j)*useh.lt.29.) go to 241
 240   continue
      go to 500
 241   j29=j
        do 250 j=j29-1,ledge
 250    zdu(j)=1.5*(tes(2,j)*rhoels(2,j)+tis(2,j)*rhoins(2,j))
        zke=volint(j29,ledge,zdu,zduma)
        tau29=zke/(zpwr-zls)
        tau29s=zke/zpwr
        write(iunit,408) tau29,tau29s,j29
 408   format('tau29',1pe12.4,'tau-s',1pe12.4,'r29',2(1pe12.4))
c   ballooning stability, following dempirc
 500  do 259 j=lcentr,ledge
 259   idum(j)=0
      do 260 j=lcentr,ledge
      zalpc=-2.*emu0*aprcri(j,2)*usil
      zalps=8.0*fcpi*max(0.0,-slpts(j))/uisb**2
      if(zalps.ge.zalpc.and.-aprcri(j,2).gt.epslon) idum(j)=1
 260   continue
      zdum1=tai*uiet
      zdum2=tes(1,lcentr)*useh
c     write(iunit,409) (idum(j),j=lcentr,ledge)
c409  format('balloon',(50i1))
      write(iunit,409) (idum(j),j=lcentr,ledge),zdum1,zdum2
 409  format('balloon',(50i1),' t=',f5.2,' T=',f5.1)
        return
        end
c@dprint
c     les 21-mar-89 new
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine dprint(ktype,k)
c
c
cl    print d-3he fusion and cmg transport data
c
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'cd3he.m'
      include 'cmg.m'
c
      dimension zdu(mj),zdu2(mj),zdu1(mj),zdu3(mj),zdu4(mj),
     1 zdu5(mj),zdu6(mj),zrad(mj)
c
c-----------------------------------------------------------------------
c   no isub
c-----------------------------------------------------------------------
      if(ktype.gt.0) return
      if(k.ne.1) go to 200
c   no title page
  100   continue
c                                      call expert ?
c
cl    first page of large time-step edit
cl    radial profiles
c   print fast particle density, energy, slowing down times
c     fusion reaction rates, fusion heatingdue to each species,
c     total powers, neutron production
c
c   zone center radius (cm)
 200    do 102 j=1,mzones
 102    zrad(j)= 0.5*(avi(j+1,15,1)+avi(j,15,1))*uiel ! zone center
      i1 = nskip+1
c
      write(nprint,1001) ash4,ashp
 1001  format('ash removal fractions - 4he',f8.2, ' protons',f8.2)
        do 120 j=i1,ledge,nskip
      zdu1(j)=wed3fs(j)*1.0e-13
      zdu2(j)=weash(j)*1.0e-13
      zdu3(j)=wid3fs(j)*1.0e-13
      zdu4(j)=wid3fl(j)*1.0e-13
 120    zdu5(j)=wntn(j)*1.0e-13
c     write(nprint,1002) (zrad(j),wed3fs(j),weash(j),wid3fs(j),
c    1 wid3fl(j),wntn(j),j=i1,ledge,nskip)
      write(nprint,1002) (zrad(j),zdu1(j),zdu2(j),zdu3(j),
     1 zdu4(j),zdu5(j),j=i1,ledge,nskip)
      do 130 j=i1,ledge,nskip
      zdu1(j)=epde(j)*1.6022e-25
      zdu2(j)=ep3e(j)*1.6022e-25
      zdu3(j)=etre(j)*1.6022e-25
      zdu4(j)=e3e(j)*1.6022e-25
      zdu5(j)=e4te(j)*1.6022e-25
      zdu6(j)=e43e(j)*1.6022e-25
 130    zdu(j)=e42te(j)*1.6022e-25
      write(nprint,1003) (zrad(j),zdu1(j),zdu2(j),zdu3(j),
     1 zdu4(j),zdu5(j),zdu6(j),zdu(j),j=i1,ledge,nskip)
c     write(nprint,1003) (zrad(j),epde(j),ep3e(j),etre(j),e3e(j),
c    1  e4te(j),e43e(j),e42te(j),j=i1,ledge,nskip)
      write(nprint,1004) (zrad(j),fnpd(j),fnp3(j),fnt(j),fn3(j),
     1   fn4t(j),fn43(j),fn42t(j),j=i1,ledge,nskip)
 1004   format(3x,'r fast n: p (d+d)     p (d+3he)  trit        3he',
     1   '         4he (d+T)   4he (d+3he) 4he (t+t)',/,
     2   (1p8e12.4))
        write(nprint,1005) (zrad(j),tepd(j),tep3(j),tet(j),te3(j),
     1 te4t(j),te43(j),te42t(j),j=i1,ledge,nskip)
 1005   format(3x,'r tau-e-sd:p (d+d)    p (d+3he)   trit       3he',
     1    '         4h3 (d+t)   4he (d+3he) 4he (t+t)',/,
     2    (1p8e12.4))
      write(nprint,1006) (zrad(j),fddt(j),fdd3(j),fdt(j),fd3(j),
     1   ftt4(j),fdtb(j),fd3b(j),j=i1,ledge,nskip)
 1006   format(3x,'r fus rate: d+d        d+d->3he    d+t',
     1   '         d+3he       t+t         d+t (fast)  d+3 (fast)',/,
     2   (1p8e12.4))
        write(nprint,1007) (zrad(j),drnp(j),drnd(j),drnt(j),drn3(j),
     1   drn4(j),j=i1,ledge,nskip)
 1007   format(3x,'r dn/dt,th: p        d           t           3he',
     1  '         4he',/,(1p6e12.4))
c
 1002   format(3x,'r  powers: e heat',5x,'e ashrem',4x,'i heat',6x,
     1 'i loss-burn',1x,'neutron',/,(1p6e12.4))
 1003   format(3x,'r  powers->e: pd        p3          trit',
     1  '        3he         4 (d+t)     4 (d+3he)  4 (t+t)',/,
     2  (1p8e12.4))
c
c   transport coefficients
      if(lcmg.le.epslon) return
      write(nprint,1008) (zrad(j),xcmg0(j),xub(j),xcmg(j),vxcmg(j),j=
     1 1,ledge)
 1008  format(//,3x,'r chi(cm2/s)    cmg        ubiq         tot     ',
     1 ' v-an (cm/s)',/,(1p5e12.4))
      return
      end

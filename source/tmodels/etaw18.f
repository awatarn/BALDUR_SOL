!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| 
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt
!| \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center}
!| {\bf {\tt etaw18.tex} \\
!| Toroidal Ion Temperature Gradient Mode \\
!| \vspace{1pc}
!| Glenn Bateman \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| Jan Weiland and Hans Nordman \\
!| Chalmers University of Technology \\
!| G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Jon Kinsey \\
!| General Atomics \\
!| P.O. Box 85608, San Diego, CA 92186} \\
!| \vspace{1pc}
!| \today
!| \end{center}
!| 
!| 
c
c@etaw18.f
c
c  ***************************************************
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. THE MODIFICATIONS CONSIST OF THE INCLUSION
c  OF IMPURITIES IN THE DILUTION APPROX. IN THE SYSTEM WITH
c  4 EQUATIONS AND THE INCLUSION OF COLLISIONS ON TRAPPED
c  ELECTRONS IN THE FULL SYSTEM. THIS SYSTEM THEN CONSISTS
c  OF 7 EQUATIONS. WHEN PARALLEL ION MOTION IS INCLUDED THERE
c  ARE 8 EQUATIONS WITH COLLISIONS. IN ETAW18 ALSO ELECTROMAGNETIC
c  EFFECTS ARE INCLUDED.
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c  ********************************************************
cc
c   The variables are ordered as: e phi/Te, Th, Nh, Te, Nz, Tz, F, Vp,
c   Av, K where F is due to collisions, Vp is due to parallel ion motion,
c   Av is the vector potential and K is its time derivative.
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c rap 31-jul-02 dimag function is replaced by aimag
c
      subroutine etaw18 (letain, cetain, lprintin, neq
     & , epsnein, epsnhin, epsnzin
     & , epstein, epsthin, epstzin, tauzin, tauhin
     & , fnzin, czin, azin, zfs, betae, ftrapein
     & , vef, q, S, aspinv, kappa, hydmass, ekyrhoin, ekparlin
     & , ndim, omega, gamma, zamrt, zbmrt, zamit, zbmit, zvr, zvi
     & , zalfr, zalfi, zbeta, nmodes, nerr, ITS,ITERA, em)
c
      implicit none
c
      INTEGER idp
      PARAMETER ( idp = 15 )
c
      LOGICAL inital, lmatv
      data inital /.true./, lmatv /.true./
c
      REAL  cetain(32), letain(32)
     &  , omega(idp), gamma(idp)
      COMPLEX ZZ(10)
c
      REAL epsnhin, epsnzin, epstein, epsthin, epstzin, tauhin, tauzin
     & , fnzin, czin, azin, ftrapein, epsnein
     & , ekyrhoin, g, etai, etae
     & , zb, si, eq, kiq, kxq, bt, bt1, vef, tvr, ftr, ftrt
c
      REAL KAPPA, RAV, GAV,ALA,SH,SH2,WZI,HCN,He,WMHD,EPS,AVC
      REAL H1,XH,EN,EI,EE,TAU,TAUI,FL,FT,GM,BTA,XT,R,HQR,HQI,WM
      REAL RFL,fh,fs,ftt,ftf,fzf,alff,WZIMAX,WIMAX(100)
      COMPLEX ALPC,ALPK,ALPHA,HQ,WZ,WZ1,WZ2,IU,H2,E,E1,WZJ(100),WZH
      COMPLEX WZP
c
      SAVE WZJ
      SAVE WIMAX
c
      INTEGER lprintin, lprint, neq, idim, ndim
     & , ieq, j1, j2, j, iret, IK, IST, IM
c
c ndim  = first dimension of the 2-D array difthi
c           and the maximum number of unstable modes allowed
c ieq   = number of equations
c
      REAL zamr(idp,idp),zami(idp,idp),zbmr(idp,idp),zbmi(idp,idp)
     &  ,zamrt(idp,idp), zbmrt(idp,idp), zamit(idp,idp), zbmit(idp,idp)
     &  ,zalfr(idp),zalfi(idp),zbeta(idp),zvr(idp,idp),zvi(idp,idp),ztol
c
      INTEGER iter(idp), ifail
c
c zamr(i,j) = matrix A
c zbmr(i,j) = matrix B
c   Note that the eigenvalues are
c omega(j) = ( zalfr(j) + i zalfi(j) ) / zbeta(j)
c where beta(j) will be 0.0 in the case of an infinite eigenvalue
c zvr(j) = eigenvector real part
c zvi(j) = eigenvector imag part
c
      REAL wr,wi,H,fft,fzft
c
      REAL zepsilon, zepsmach, zepsqrt
     & , zetah, zetaz, zetae
     & , zepsne, zepsnh, zepste, zepsth
     & , ztauh, ztauz, zep2nh, zep2nz, zep2ne, zft
     & , zimp, zfnz, zmass, zflh, zflz,zfs,A,em1
      REAL q,S,Cs,alp,k1,k2,kps,betae,eni,alf,kpc
      INTEGER ITC,ITL,ITS,ITERA, nmodes, nerr, em
      REAL TOL, aspinv, hydmass, ekparlin
c
      REAL zzz,zzzz,zzzzz, c1, c2
c
c
      bt           =   1.5
      IK           =   1
      IST          =   1
      ITC          =   1
      ITL          =   100
      TOL          =   1.e-6
c
      HCN          = ( 1 / vef ) * ( 1836 * hydmass ) * (1/aspinv)
     & * ( q * ekyrhoin ) ** ( -2 )
c
      WMHD         = 1 / ( epsnein * ekyrhoin )
c
c
c ...  STRUCTURE OF MATRIX EQUATION ...
c
c ...  omega*zbmr(i,j) = zamr(i,j)
c
c    variables i=1,6: efi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
c    variables j=1,6 same as for i
c
c  ...................................................
c
c
c..initialize variables
c
      if ( inital ) then
c
        idim = idp
c
        tvr = 2./3.
        ftr = 5./3.
        ftrt = ftr/tauhin
        em1=em
c        em1=1.D0
c
        A = azin
        zepsilon = 1.e-4
c
        zepsmach = 0.5
  2     if ( 0.5 * zepsmach + 1.0 .gt. 1.0 ) then
          zepsmach = 0.5 * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
        inital = .false.
c
      endif
c
      lprint = lprintin
c
      ftrt=ftr/tauhin
      ieq = max ( 2, neq )
c
      IU=(0.,1.)
c
      ftrt = ftr/tauhin
      bt1 = bt - 2.5
c
      iret = 0
c..print header
c
      if ( lprint .gt. 8 ) then
c
        write (6,*)
        write (6,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine etaw18'
        write (6,*) '(all frequencies normalized by omega_{De})'
        write (6,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
      endif
c
c
c..check validity of input data
c
      if ( neq .lt. 2 ) call abortb (6
     & ,' neq .lt. 2 in sbrtn etaw18')
c
      if ( abs(epstein) .lt. zepsqrt ) call abortb (6
     & ,' abs(epstein) .lt. zepsqrt in sbrtn etaw18')
c
      if ( abs(epsthin) .lt. zepsqrt ) call abortb (6
     & ,' abs(epsthin) .lt. zepsqrt in sbrtn etaw18')
c
      if ( abs(epstzin) .lt. zepsqrt ) call abortb (6
     &   ,' abs(epstzin) .lt. zepsqrt in sbrtn etaw18')
c
      if ( abs(epsnzin) .lt. zepsqrt ) call abortb (6
     &   ,' abs(epsnzin) .lt. zepsqrt in sbrtn etaw18')
c
      if ( czin .lt. 1.0 ) call abortb (6
     &   ,' czin .lt. 1.0 in sbrtn etaw18')
c
      zepsth = epsthin
      zepsnh = epsnhin
      zepste = epstein
c
c..compute the rest of the dimensionless variables needed
c
      zetah  = zepsnh / zepsth
      zetae  = zepsne / zepste
c
c  *******  NOTE THE INVERSE DEFINITION OF ZTAUH ! ******
      ztauh  =1./tauhin
c
      zep2nh=epsnhin
      zep2ne=epsnein
      zetae = zep2ne/zepste
      zft    = ftrapein
      zflh   = ekyrhoin**2
      eni = 1./zep2ne
c
      zimp   = czin
      zmass  = azin
c
        zfnz   = fnzin * zimp
c        zetaz  = zepsnz / zepstz
        zetaz=epsnzin/epstzin
c
c  ******  NOTE THE INVERSE DEFINITION OF ZTAUZ ! ******
        ztauz = 1./(czin*tauzin)
        zep2nz=epsnzin
c        zep2nz = 2.0 * zepsnz
        zflz   = zmass * zflh / zimp**2
c
c..diagnostic output
c
      if ( lprint .gt. 8 ) then
        write (6,*)
        write (6,*) '--------------------------------------'
        write (6,*)
        write (6,*)
        write (6,*) zetah,' = zetah'
        write (6,*) zetaz,' = zetaz'
        write (6,*) zetae,' = zetae'
        write (6,*) ztauh,' = ztauh'
        write (6,*) ztauz,' = ztauz'
        write (6,*) zep2nh,' = zep2nh'
        write (6,*) zep2nz,' = zep2nz'
        write (6,*) zep2ne,' = zep2ne'
        write (6,*) zft,' = zft'
        write (6,*) zimp,' = zimp'
        write (6,*) zmass,' = zmass'
        write (6,*) zfnz,' = zfnz'
        write (6,*) zflh,' = zflh'
        write (6,*) zflz,' = zflz'
        write (6,*) zepsqrt,' = zepsqrt'
        write (6,*) zepsmach,' = zepsmach'
        endif
c
c..set matricies for eigenvalue equation
c
      do j1=1,idim
        zalfr(j1) = 0.0
        zalfi(j1) = 0.0
        zbeta(j1) = 0.0
        do j2=1,idim
          zamr(j1,j2) = 0.0
          zami(j1,j2) = 0.0
          zbmr(j1,j2) = 0.0
          zbmi(j1,j2) = 0.0
          zvr(j1,j2) = 0.0
          zvi(j1,j2) = 0.0
        enddo
      enddo
c
      if ( ieq .eq. 2 ) then
c
c..two equations when trapped particles and FLR effects omitted
c
        zamr(1,1) = ( 1.0 / zep2nh ) - ztauh - 1.0
        zamr(1,2) = - ztauh
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(1,1) = 1.0
        zbmr(1,2) = 0.0
        zbmr(2,1) = - tvr
        zbmr(2,2) = 1.0
c
      elseif ( ieq .eq. 4 ) then
c
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
      if  ( lprint .gt. 5 ) then
        write (6,*)
        write (6,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
      endif
c
c       ion continuity
c
        zamr(1,1) = 1.0 - zep2nh - zflh * ztauh * ( 1.0 + zetah )
        zamr(1,2) = - zep2nh * ztauh
        zamr(1,3) = - zep2nh * ztauh
c
        zbmr(1,1) = zflh * zep2nh
        zbmr(1,3) = zep2nh
c
c  ion energy
c
        zamr(2,1) = zetah - tvr
        zamr(2,2) = - zep2nh * ztauh * ftr
c
        zbmr(2,2) =   zep2nh
        zbmr(2,3) = - zep2nh * tvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation. The dilution factor 1-zfnz has now been
c   added.
c
        zamr(3,1) = zft - zep2ne
        zamr(3,3) = zep2ne * (1. - zfnz - zfs)
        zamr(3,4) = zft * zep2ne
c
        zbmr(3,1) = ( zft - 1.0 ) * zep2ne
        zbmr(3,3) = zep2ne * (1. - zfnz -zfs)
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr )
        zamr(4,4) = zft * zep2ne * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * zep2ne * tvr
        zbmr(4,3) = - zep2ne * tvr
        zbmr(4,4) = zft * zep2ne
c
      else if ( ieq .eq. 6 ) then
c
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c  hydrogen density
c
        zamr(1,1) = - 1.0
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  trapped electron density
c
        zamr(3,1) = - 1.0 + zft / zep2ne
        zamr(3,3) = 1.0 - zfnz -zfs
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zetae - tvr ) / zep2ne
        zamr(4,4) = zft * ftr
c
        zbmr(4,1) = ( 1.0 - zft ) * tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs) * tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0
     &    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz * ftr
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
      else if ( ieq .eq. 7 ) then
c
c..Seven equations with impurities, trapped electrons, parallel ion motion
c and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here F is defined as F = GM*e phi/T_e where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (6,*)
      write (6,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and Vp'
      endif
c
c********
      H=0.5*ABS(S)/q
c********
c
c  hydrogen density
c
      zamr(1,1) = -1.
     & + (1. - zflh * ztauh * (1. + zetah ))/zep2nh
      zami(1,1) = -H
      zamr(1,2) = - ztauh
      zami(1,2) = -ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = -ztauh*H
C
      zbmr(1,1) = zflh
      zbmr(1,3) = 1.
c
c hydrogen energy
c
      zamr(2,1) = (zetah - tvr )/zep2nh
      zamr(2,2) = - ztauh * ftr
c
      zbmr(2,2) = 1.
      zbmr(2,3) = -tvr
c
c trapped electron density
c
      zamr(3,1) = -1. + zft/zep2ne
      zamr(3,3) = 1. - zfnz -zfs
      zamr(3,4) = zft
      zamr(3,5) = zfnz
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz -zfs
      zbmr(3,5) = zfnz
c
c trapped electron energy
c
      zamr(4,1) = zft * (zetae - tvr )/zep2ne
      zamr(4,4) = zft * ftr
c
      zbmr(4,1) = (1. - zft) * tvr
      zbmr(4,3) = - (1. - zfnz -zfs) * tvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * tvr
c
c impurity density
c
      zamr(5,1) = -1.
     & + ( 1. - zflz * ztauz * (1. + zetaz ))/zep2nz
      zami(5,1)=-zimp*H/A
      zamr(5,5) = - ztauz
      zami(5,5) = -zimp*ztauz*H/A
      zamr(5,6) = - ztauz
      zami(5,6) = -zimp*ztauz*H/A
c
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
c
c impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = - ztauz * ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
c
c
      else if ( ieq .eq. 8 ) then
c
c..Eight equations with impurities, trapped electrons, parallel ion
c  motion collisions and FLR
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Vp
c  Here F = omega*gamma (collisions).
c
      k1 = 0.25*q*q*zflh*SQRT(ABS((1.+zetah)*ztauh/((1.-zft)*zep2nh)))
      k2 = q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
      alp = 0.5*(k1+SQRT(k1+S*S*k2))
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
c
c
c hydrogen density
c
        zamr(1,1) = - 1.0
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,8) = kps
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.0
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr ) / zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.0
        zbmr(2,3) = - tvr
c
c  total electron density expressed in ion and imp densities
c
        zamr(3,1) = - 1.0 + zft / zep2ne
        zami(3,1) = vef*(1.-zft)
        zamr(3,3) = 1.0 - zfnz -zfs
        zami(3,3) = -vef*(1. - zfnz -zfs)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - 1.0
        zbmr(3,3) = 1.0 - zfnz -zfs
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*(zetae - tvr)/zep2ne
        zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
        zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
        zamr(4,4) = zft * ftr
        zami(4,5) = -vef*tvr*bt1*zfnz
        zami(4,7) = -ftr*vef*zft
c
        zbmr(4,1) = ( 1.0 - zft ) *tvr
        zbmr(4,3) = - ( 1.0 - zfnz -zfs ) *tvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * tvr
c
c  impurity density
c
        zamr(5,1) = - 1.0
     &    + ( 1.0 - zflz * ztauz * ( 1.0 + zetaz ) ) / zep2nz
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = 1.0
c
c  impurity energy
c
        zamr(6,1) = ( zetaz - tvr ) / zep2nz
        zamr(6,6) = - ztauz * ftr
c
        zbmr(6,5) = - tvr
        zbmr(6,6) = 1.0
c
c  variable F
c
         zamr(7,1) = zetae/zep2ne - 1.
         zami(7,1) = vef
         zamr(7,7) = 1.
         zami(7,7) = -vef
c
         zbmr(7,1) = -1.
         zbmr(7,7) = 1.
c
c     Parallel ion motion Vpi/Cs
c
         zamr(8,1) = kps
         zamr(8,2) = kps*ztauh
         zamr(8,3) = kps*ztauh
c
         zbmr(8,8) = 1.
c
c
      else if ( ieq .eq. 9 ) then
c
c..Nine  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta and parallel motion of impurities
c
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
c   magnetic vector potential.
c

      EN=zep2nh
      ENI=1./EN
      EI=zetah
      EE=zetae
      TAU=1./ztauh
      TAUI=ztauh
      FL=zflh
      FT=zft
c
      GM=1./(1.-FT)
      BTA=FTRT*(GM+TAUI)
      H1=4.*q*q*FL
      XT=1./(1.+TAUI)
c
      ALA=em*q*q*betae*(1.+EE+TAUI*(1.+EI))/EN
      AVC=EPS*(1.-1./q**2)
      SH2=2.*S-1.+(KAPPA*(S-1.))**2
      SH=SQRT(SH2)
      H=0.5*ABS(SH)/q
      H2=IU*H
      GAV=1.
c
      IF(IST.NE.1) GO TO 800
      E1=FTRT*(1.+FL)-ENI+FL*TAUI*ENI*(1.+EI)+(GM+FTRT)*GAV
     &+H2*(1.+FTRT)
      E1=0.5*E1/(1.+FL)
      E=(TAUI*ENI*GM*(EI-TVR)+BTA)*(GAV+H2)-FTRT*ENI*(1.-FL*TAUI*(1.
     &+EI))
      E=E/(1.+FL)
      WZ1=-E1+SQRT(E1*E1-E)
      WZ2=-E1-SQRT(E1*E1-E)
      WZ=WZ1
c
      IF(aimag(WZ2).GT.aimag(WZ1)) WZ=WZ2
c      WZ=WZ-IU*ROT*WEXB
      WZI = aimag(WZ)
      IF(WZI.LT.0.001) WZ=WZ+IU*(0.001-WZI)
      HQ=-IU*H
      ALPHA=-IU*ABS(SH)*q*FL
      R=2.*ABS(REAL(WZ*ALPHA))
c
      WZJ(IK)=WZ

      WZH=EN*WZ
c      WRITE(*,10001) WZH
10001 FORMAT(2X,'WZ=',2G11.3)
C
  800 CONTINUE
      WZ=WZJ(IK)
      WZP=WZ
      ALPK=0.5*SH*SQRT(H1*XT*FL*(1.+TAUI*(1.+EI)/(EN*WZ)))
      IF(REAL(ALPK).GE.0.) GOTO 801
      ALPK=-ALPK
  801 CONTINUE
      ALPC=-IU*ALPK
      ALPHA=-IU*ABS(SH)*q*FL
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(REAL(WZ*ALPC))
      HQ=2.*ALPC/H1
  802 GAV=(1.+0.5*S/R)*EXP(-0.25/R)
      GAV=GAV-0.5*ALA*(1.-EXP(-1./R))
      IF(GAV.LT.0.001) GAV=0.001
      GAV=GAV-AVC
c
c      k1 = 0.25*q*q*zflh*SQRT((1.+zetah)*ztauh/((1.-zft)*zep2nh))
c      k2 = q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
c      alp = 0.5*(k1+SQRT(k1+S*S*k2))
      alp=0.5*R
      IF(alp.LT.0.1) alp=0.1
      alf = alp/(2.D0*zflh*q*q*betae*(1.D0 - zft))
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*SH2/alp
c      WRITE(*,10002) alp,SH2,RAV
10002 FORMAT(2X,'alp=',G11.3,' SH2=',G11.3,' RAV=',G11.3)
c      WRITE(*,10003) XH,GAV,alf
10003 FORMAT(2X,'XH=',G11.3,' GAV=',G11.3,' alf=',G11.3)
c
c
c  *********
      HQR = REAL(HQ)
      HQI = aimag(HQ)
c  *********
c hydrogen density
c
        zamr(1,1) = - GAV+HQR
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zami(1,1) = HQI
        zamr(1,2) = (HQR-GAV)*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-GAV)*ztauh
        zami(1,3) = ztauh*HQI
        zamr(1,8) = -em*ztauh*HQR*(1.+zetah)/(kpc*zep2nh)
        zami(1,8) = -em*ztauh*HQI*(1.+zetah)/(kpc*zep2nh)
        zamr(1,9) = -em*HQR/kpc
        zami(1,9) = -em*HQI/kpc
c
        zbmr(1,1) = zflh
        zbmr(1,3) = 1.
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr )/zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr
c
c  total electron density expressed in ion density and imp density
c
       zamr(3,1) = -1. + zft/zep2ne
       zami(3,1) = vef*(1.-zft)
       zamr(3,3) = 1. - zfnz -zfs
       zami(3,3) = -vef*(1. - zfnz - zfs)
       zamr(3,4) = zft
       zamr(3,5) = zfnz
       zami(3,5) = -vef*zfnz
       zami(3,7) = vef*zft
       zamr(3,8) = -em*(1. - zft)/(kpc*zep2ne)
       zami(3,8) = em*(1.-zft)*vef/(kpc*zep2ne)
       zamr(3,9) = em*(1. - zft)*(1.+eni)/kpc
       zami(3,9) = -em*(1.-zft)*vef/kpc
c
      zbmr(3,1) = zft - 1.
      zbmr(3,3) = 1. - zfnz - zfs
      zbmr(3,5) = zfnz
      zbmr(3,9) = em*(1. - zft)/kpc
c
c  trapped electron energy
c
      zamr(4,1) = zft*(zetae - tvr)/zep2ne
      zami(4,1) = vef*tvr*(bt-2.5*(1.-zft))
      zami(4,3) = -vef*tvr*bt1*(1.-zfnz -zfs)
      zamr(4,4) = zft*ftr
      zami(4,5) = -vef*tvr*bt1*zfnz
      zami(4,7) = -ftr*vef*zft
c
      zbmr(4,1) = (1. - zft)*tvr
      zbmr(4,3) = -(1. - zfnz -zfs)*tvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*tvr
c
c  impurity density
c
      zamr(5,1) = - GAV +zimp*HQR/A
     & +(1. -zflz*ztauz*(1.+zetaz))/zep2nz
      zami(5,1) = zimp*HQI/A
      zamr(5,5) = (HQR*zimp/A-GAV)*ztauz
      zami(5,5) = zimp*ztauz*HQI/A
      zamr(5,6) = (HQR*zimp/A-GAV)*ztauz
      zami(5,6) = zimp*ztauz*HQI/A
      zamr(5,8) = -em*HQR*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zami(5,8) = -em*HQI*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zamr(5,9) = -em*HQR*zimp/(kpc*A)
      zami(5,9) = -em*HQI*zimp/(kpc*A)
c
      zbmr(5,1) = zflz
      zbmr(5,5) = 1.
c
c  impurity energy
c
      zamr(6,1) = (zetaz - tvr)/zep2nz
      zamr(6,6) = -ztauz*ftr
c
      zbmr(6,5) = -tvr
      zbmr(6,6) = 1.
c
c  variable F
c
      zamr(7,1) = zetae/zep2ne - 1.
      zami(7,1) = vef
      zamr(7,7) = 1.
      zami(7,7) = -vef
c
      zbmr(7,1) = -1.
      zbmr(7,7) = 1.
c
c
c  electromagnetic parallel vectorpotential Av = e A_par/Te
c
      fft=(1.-zfnz)/(1.-zft)
      fzft=zfnz/(1.-zft)
      zamr(8,1) = em1*kpc*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,1) = em1*HQI*(fft+zimp*fzft/A)*kpc
      zamr(8,2) = em1*HQR*ztauh*fft*kpc
      zami(8,2) = em1*HQI*ztauh*fft*kpc
      zamr(8,3) = em1*HQR*ztauh*fft*kpc
      zami(8,3) = em1*HQI*ztauh*fft*kpc
      zamr(8,5) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,5) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,6) = em1*HQR*zimp*ztauz*fzft*kpc/A
      zami(8,6) = em1*HQI*zimp*ztauz*fzft*kpc/A
      zamr(8,8) = em1*((1.+zetae)*eni - alf*zflh*RAV)
     &-em1*HQR*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zami(8,8) = -em1*HQI*(fft*ztauh*(1.+zetah)/zep2nh
     &+zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*A))
      zamr(8,9)= -em1*(eni+HQR*(fft+zimp*fzft/A))
      zami(8,9) = -em1*HQI*(fft+zimp*fzft/A)
c
      zbmr(8,1) = em1*kpc
      zbmr(8,8) = em1
      zbmr(8,9)= -em1
c
c     K = omega*Av
c
      zamr(9,9) = em1
c
      zbmr(9,8) = em1
c
c
c
c ----------------------------------------------------------------
      else if ( ieq .eq. 10 ) then
c
c  Ten  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta,  parallel motion of impurities
c  and edge physics
c  Equations for e phi/Te, Ti, net, nef, Tet, Tef, n_z, T_z, F, Av
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
c   magnetic vector potential.
c

      EN=zep2nh
      ENI=1./EN
      EI=zetah
      EE=zetae
      TAU=1./ztauh
      TAUI=ztauh
      FL=zflh
      RFL=SQRT(FL)
      FT=zft
c
      GM=1./(1.-FT)
      BTA=FTRT*(GM+TAUI)
      H1=4.*q*q*FL
      XT=1./(1.+TAUI)
c
      ALA=em*2.*q*q*betae*(1.+EE+TAUI*(1.+EI))/EN
      AVC=EPS*(1.-1./q**2)
      SH2=2.*S-1.+(KAPPA*(S-1.))**2
      SH=SQRT(SH2)
      H=0.5*ABS(SH)/q
      H2=IU*H
      GAV=1.
      ITS=0
      ITERA=1
c
      IF(IST.NE.1) GO TO 900
      E1=FTRT*(1.+FL)-ENI+FL*TAUI*ENI*(1.+EI)+(GM+FTRT)*GAV
     &+H2*(1.+FTRT)
      E1=0.5*E1/(1.+FL)
      E=(TAUI*ENI*GM*(EI-TVR)+BTA)*(GAV+H2)-FTRT*ENI*(1.-FL*TAUI*(1.
     &+EI))
      E=E/(1.+FL)
      WZ1=-E1+SQRT(E1*E1-E)
      WZ2=-E1-SQRT(E1*E1-E)
      WZ=WZ1
      IF(aimag(WZ2).GT.aimag(WZ1)) WZ=WZ2
c      WZ=WZ-IU*ROT*WEXB
      WZI = aimag(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
c
      WZJ(IK)=WZ
      WZH=EN*WZ
c      WRITE(*,20001) WZH
20001 FORMAT(2X,'WZH=',2G11.3)
      zzz    = 0.001
      zzzz   = 0.2*WMHD
      WIMAX(IK)=MAX(zzz,zzzz)
c      IF (0.2*WMHD .ge. 0.001 ) THEN
c      WIMAX(IK) = 0.2*WMHD
c      ELSE
c      WIMAX(IK) = 0.001
c      ENDIF
C
  900 CONTINUE
c      WRITE(*,911) ITERA
  911 FORMAT(' ITER=',I5)
      WZ=WZJ(IK)
      WZP=WZ
      WZH=ABS(EN)*WZ
      WZIMAX=WIMAX(IK)
      ALPK=0.5*SH*SQRT(H1*XT*FL*(1.+TAUI*(1.+EI)/(EN*WZ)))
      IF(REAL(ALPK).GE.0.) GOTO 901
      ALPK=-ALPK
  901 CONTINUE
      ALPC=-IU*ALPK
      ALPHA=-IU*ABS(SH)*q*FL
      XH=ABS(ALPHA/ALPC)
      ALPC=XH*ALPC
      R=2.*ABS(REAL(WZ*ALPC))
      HQ=2.*ALPC/H1
  902 GAV=(1.+0.5*S/R)*EXP(-0.25/R)
      GAV=GAV-0.5*ALA*(1.-EXP(-1./R))
      IF(GAV.LT.0.001) GAV=0.001
      GAV=GAV-AVC
c
      fh=1.-zfnz-zfs
C ***   Note Change !!!!
c      k1 = 0.25*q*q*zflh*SQRT((1.+zetah)*ztauh/((1.-zft)*zep2nh))
c      k2 = q*q*zflh*zflh*(1.+zetah)*ztauh/((1.-zft)*zep2nh)
c      alp = 0.5*(k1+SQRT(k1+S*S*k2))
c      R=2.*ALP
c      GAV=(1.+0.5*S/R)*EXP(-0.25/R)
c      GAV=GAV-0.5*ALA*(1.-EXP(-1./R))
C *****************************
      alp=0.5*R
      IF(alp.LT.0.1) alp=0.1
      alff = alp/(2.D0*zflh*q*q*betae*fh)
      kps = 0.5*SQRT(alp/zflh)/q
      kpc = 1.D0
      RAV=1.+0.25*SH2/alp
c      WRITE(*,20002) alp,SH2,RAV
20002 FORMAT(2X,'alp=',G11.3,' SH2=',G11.3,' RAV=',G11.3)
c      WRITE(*,20003) XH,GAV,alff
20003 FORMAT(2X,'XH=',G11.3,' GAV=',G11.3,' alff=',G11.3)
20004 FORMAT(2X,' DISP10: He=',G11.3,' WZIMAX=',G11.3)
c
c
c  *********
      HQR = REAL(HQ)
      HQI = aimag(HQ)
C *******
      He=ABS(SH)*RFL*SQRT(HCN*WZIMAX)
c      WRITE(*,20004) He,WZIMAX
c      WRITE(*,20005) HQ
20005 FORMAT(' HQ=',2G11.3)
c  *********
      ftt=zft/fh
      ftf=(1.-zft)/fh
      fzf=zfnz/fh
      alf=alff/ftf
c  **********
c
c hydrogen density
c
        zamr(1,1) = - GAV+HQR
     &   + ( 1.0 - zflh * ztauh * ( 1.0 + zetah ) ) / zep2nh
        zami(1,1) = HQI
        zamr(1,2) = (HQR-GAV)*ztauh
        zami(1,2) = ztauh*HQI
        zamr(1,3) = (HQR-GAV)*ztauh*ftt
        zami(1,3) = ztauh*HQI*ftt
        zamr(1,4)=(HQR-GAV)*ztauh*ftf
        zami(1,4)=ztauh*HQI*ftf
        zamr(1,7)=-(HQR-GAV)*ztauh*fzf
        zami(1,7)=-ztauh*HQI*fzf
        zamr(1,10) = -em*ztauh*HQR*(1.+zetah)/(kpc*zep2nh)
        zami(1,10) = -em*ztauh*HQI*(1.+zetah)/(kpc*zep2nh)
c
        zbmr(1,1) = zflh
        zbmr(1,3) = ftt
        zbmr(1,4)=ftf
        zbmr(1,7)=-fzf
        zbmr(1,10)=em*HQR/kpc
        zbmi(1,10)=em*HQI/kpc
c
c  hydrogen energy
c
        zamr(2,1) = ( zetah - tvr )/zep2nh
        zamr(2,2) = - ztauh * ftr
c
        zbmr(2,2) = 1.
        zbmr(2,3) = - tvr*ftt
        zbmr(2,4)= -tvr*ftf
        zbmr(2,7)= tvr*fzf
c
c trapped electron cont.
c
       zamr(3,1) = -1. + 1./zep2ne
       zamr(3,3) = 1.
       zami(3,3) = -vef
       zamr(3,5) = 1.
       zami(3,9) = vef
c
      zbmr(3,3) = 1.
c
c  free electon cont.
c
      zamr(4,1) =1./zep2ne - 1.
      zami(4,1) =He
      zamr(4,4) =1.
      zami(4,4) =-He
      zamr(4,6) =1.
      zami(4,6) =-He
      zami(4,10) = em*He*(1.+zetae)/(zep2ne*kpc)
c
      zbmr(4,4) =1.
      zbmi(4,10) = em*He/kpc
c

c  trapped electron energy
c
      zamr(5,1) = (zetae - tvr)/zep2ne
      zami(5,1) = vef*tvr*bt
      zami(5,3) = vef*tvr
      zamr(5,5) = ftr
      zami(5,9) = -ftr*vef
c
      zbmr(5,3) = -tvr
      zbmr(5,5) = 1.
c
c     free electron energy
c
      zamr(6,1) = (zetae-tvr)/zep2ne
      zamr(6,6) = ftr
      zami(6,6) = -1.06*He
      zami(6,10) = em*1.06*He*zetae/(zep2ne*kpc)
c
      zbmr(6,4) =-tvr
      zbmr(6,6) =1.
c
c  impurity density
c
      zamr(7,1) = - GAV +zimp*HQR/A
     & +(1. -zflz*ztauz*(1.+zetaz))/zep2nz
      zami(7,1) = zimp*HQI/A
      zamr(7,7) = (HQR*zimp/A-GAV)*ztauz
      zami(7,7) = zimp*ztauz*HQI/A
      zamr(7,8) = (HQR*zimp/A-GAV)*ztauz
      zami(7,8) = zimp*ztauz*HQI/A
      zamr(7,10) = -em*HQR*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
      zami(7,10) = -em*HQI*zimp*ztauz*(1.+zetaz)/(kpc*zep2nz*A)
c
      zbmr(7,1) = zflz
      zbmr(7,7) = 1.
      zbmr(7,10) = em*HQR*zimp/(kpc*A)
      zbmi(7,10) = em*HQI*zimp/(kpc*A)
c
c  impurity energy
c
      zamr(8,1) = (zetaz - tvr)/zep2nz
      zamr(8,8) = -ztauz*ftr
c
      zbmr(8,7) = -tvr
      zbmr(8,8) = 1.
c
c  variable F
c
      zamr(9,1) = zetae/zep2ne - 1.
      zami(9,1) = vef
      zamr(9,9) = 1.
      zami(9,9) = -vef
c
      zbmr(9,1) = -1.
      zbmr(9,9) = 1.
c
c
c  Amperes law
c
      fft=(1.-zfnz)/(1.-zft)
      fzft=zfnz/(1.-zft)
c
      zamr(10,1) = em1*HQR*(1.+zimp*fzf/A)
      zami(10,1) = em1*(HQI*(1.+zimp*fzf/A)-ftf*He)
      zamr(10,2) = em1*HQR*ztauh
      zami(10,2) = em1*HQI*ztauh
      zamr(10,3) = em1*HQR*ztauh*ftt
      zami(10,3) = em1*HQI*ztauh*ftt
      zamr(10,4) = em1*HQR*ztauh*ftf
      zami(10,4) = em1*He*ftf+em1*HQI*ztauh*ftf
      zami(10,6) = em1*He*ftf
      zamr(10,7) = -em1*HQR*fzf*(ztauh-zimp*ztauz/A)
      zami(10,7) = -em1*HQI*fzf*(ztauh-zimp*ztauz/A)
      zamr(10,8) = em1*HQR*fzf*zimp*ztauz/A
      zami(10,8)= em1*HQI*fzf*zimp*ztauz/A
      zamr(10,10) = -em1*(zflh*RAV*alff+HQR*((1.+zetah)*ztauh/zep2nh
     &+fzf*zimp*ztauz*(1.+zetaz)/(A*zep2nz)))/kpc
      zami(10,10)=-em1*(HQI*((1.+zetah)*ztauh/zep2nh+
     &fzf*zimp*ztauz*(1.+zetaz)/(A*zep2nz))+ftf*He*(1.+zetae)/zep2ne)
     &/kpc
c
      zbmr(10,10) = em1*HQR*(1.+zimp/A)/kpc
      zbmi(10,10) = em1*(HQI*(1.+zimp/A)-He*ftf)/kpc
c
      else
c
      write(6,*)
      write(6,*) ieq,' = ieq in sbrtn etaw18'
      call abortb(6
     &,'the value of ieq is wrong in sbrtn etaw18')
c
       endif
c
c
c..save copy of matrix which is over-written by sbrtn f02bjf
c
      do j2=1,ieq
        do j1=1,ieq
           zamrt(j1,j2) = zamr(j1,j2)
           zbmrt(j1,j2) = zbmr(j1,j2)
           zamit(j1,j2) = zami(j1,j2)
           zbmit(j1,j2) = zbmi(j1,j2)
        enddo
      enddo
c
c..diagnostic output
c
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (6,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (6,*)
c
        write (6,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write(6,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
        write (6,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (6,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
       write(6,*) ' zbmi(j1,j2)  j2->'
       do j1=1,ieq
         write(6,192) (zbmi(j1,j2),j2=1,ieq)
       enddo
c
 192  format (1p5e12.4)
      endif
c
c      write(*,193) zep2ne,zep2nh,zep2nz
  193 format(2X,'zep2n=',G12.5,' zep2nh=',G12.5,' zep2nz=',G12.5)
c      zgne=2./zep2ne
c      zgnh=2./zep2nh
c      zgnz=2./zep2nz
c
c      write(*,194) zgne,zgnh,zgnz
  194 format(2X,'zgne=',G12.5,' zgnh=',G12.5,' zgnz=',G12.5)
c

c..find the eigenvalues using NAG14 routine f02bjf or f02gjf
c
       zzzzz = 0.0
       ztol=max(zzzzz,cetain(32))
c      if (cetain(32) .ge. 0.0) then
c      ztol = cetain(32)
c      else
c      ztol = 0.0
c      endif
c
      ifail = 1
c
c      if( vef .gt. 0.0001 ) go to 201
c
c      call f02bjf ( ieq, zamr, idim, zbmr, idim, ztol
c     & , zalfr, zalfi, zbeta, lmatv, zvr, idim, iter, ifail )
c
c      go to 202
c
c  201 continue
c
c      call f02gjf ( ieq,zamr,idim,zami,idim,zbmr,idim,zbmi,idim,ztol
c     & ,zalfr,zalfi,zbeta,lmatv,zvr,idim,zvi,idim,iter,ifail )
c
c...onjun add tomsqz eigen solver
c
c
      call tomsqz (idim, ieq, zamr, zami, zbmr, zbmi, zalfr, zalfi
     & , zbeta, zvr, zvi, ifail)
c
  202 continue
c
        nerr = ifail
c
       if( ifail .le. 0 ) goto 210
       iret = 1
c
      go to 215
c
  210 continue
c..compute the complex eigenvalues
c
      do j=1,ieq
c
        zb = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsilon ) zb = zepsilon
        omega(j) = zalfr(j) / zb
        gamma(j) = zalfi(j) / zb
        ZZ(j)=omega(j)+(0.D0,1.D0)*gamma(j)
c
      enddo
c
      if(.not.lmatv) go to 215
      if(lprint.NE.2) go to 215
      do 214 j=1,ieq
      if(gamma(j).LT.0.001) go to 214
      write(*,220) omega(j), gamma(j)
      do 212 j1=1,ieq
  212 write(*,211) zvr(j1,j),zvi(j1,j),j1
  214 continue
  211 format(" Vr=",G11.3," Vi=",G11.3," j=",I5)
c
  220 format(2X,'wr=',G11.4,' wi=',G11.4)
c
  215 continue
c
        if ( lprint .gt. 6 ) then
        do j=1,ieq
        write(*,220) omega(j),gamma(j)
      enddo
      endif
c
      WM=0.
      DO 00062 j=1,ieq
      IF(gamma(j).LT.WM) GOTO 00062
      WM=gamma(j)
00062 CONTINUE
c      WM=WM-IU*ROT*WEXB
      IF(WM.LT.0.01) WM=0.01
      WIMAX(IK)=WM
c
      WM=0.
      IM=0.
      DO 00082 j=1,ieq
      IF(omega(j).GE.0.) GO TO 00082
      IF(gamma(j).LT.WM) GO TO 00082
      WM=gamma(j)
      IM=j
00082 CONTINUE
c
      IF(IM.EQ.0) GO TO 00083
c      WZ=ZZ(IM)-IU*ROT*WEXB
      WZ=ZZ(IM)
      IF(ABS((WZ-WZP)/WZ).LT.TOL) ITS=1
      WZ=(WZ+WZP)/2.
      WZI = aimag(WZ)
      IF(WZI.LT.0.01) WZ=WZ+IU*(0.01-WZI)
      WZJ(IK)=WZ
00083 CONTINUE
c
c..onjun adding
c
      if (IM .eq. 0 .and. ITERA .gt. 1) then
          ITS = 1
          return
      endif
c
c..end of Onjun adding
c
      IF(ITC.NE.1) GO TO 00085
      IF(ITS.EQ.1) GO TO 00085
      IF(ITERA.GE.ITL) GO TO 00085
c  *** A new iteration will be made  ***
      DO j1=1,idim
      zalfr(j1)=0.
      zalfi(j1)=0.
      zbeta(j1)=0.
c
      DO j2=1,idim
      zamr(j1,j2)=0.
      zami(j1,j2)=0.
      zbmr(j1,j2)=0.
      zbmi(j1,j2)=0.
      zvr(j1,j2)=0.
      zvi(j1,j2)=0.
       enddo
      enddo
c
      ITERA=ITERA+1
      GO TO 900
00085 CONTINUE
c
      return
      end

!| \end{document}
!| 
!| 
!| 
!| 
!| 
!| 
!| 
!| 
!| 

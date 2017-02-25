!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| 
!| \documentclass{article}
!| \usepackage{longtable}
!| 
!| \headheight 0pt \headsep 0pt \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center}
!| {\bf {\tt weiland18flux} \\
!| Weiland Model for Drift Mode Transport} \\
!| \vspace{1pc}
!| Jan Weiland \\
!| Chalmers University of Technology \\
!| G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Thawatchai Onjun, Glenn Bateman and Jon Kinsey \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| P\"{a}r Strand \\
!| Oak Ridge National Laboratory \\
!| Oak Ridge, TN 37830
!| \vspace{1pc} \today
!| \end{center}
!|
      subroutine weiland18flux (
     &   letain,    cetain,    lprint,    neq,     nout,    gnein
     & , gnhin,     gnzin,     gtein,     gthin,   gtzin,   tauhin
     & , tauzin,    fnzin,     czin,      azin,    fnsin,   betaein
     & , ftrapein,  vef,       q,         shear,   aspinv,  kappa
     & , hydmass,   ekyrhoin,  ekparlin,  wexb
     & , hcn,       rot,       tol,       itl
     & , ndim,      omega,     gamma
     & , chieff,    flux
     & , nmodes,    nerr,      ITS,       ITERA,    em)
c
      implicit none
c
      integer idp
c
      parameter ( idp = 12 )
c
      integer neq, ieq, imatrx, j, j1, j2, jd, nerr, nout, ndim
     & , lprint, zzz, nmodes, idim, iter(idp), ifail, ITERA, ITS
     & , itl, ist, ik, itc
c
      integer letain(*)
c
      real cetain(*)
c
      real omega(*), gamma(*), chieff(*), flux(*), em

      real  zalfr(idp), zalfi(idp)
     & , zvr(idp,idp), zvi(idp,idp)
     & , zeigenv(idp,idp), ztempa(idp), ztempb(idp), zbeta(idp)
     & , zamrt(idp,idp), zbmrt(idp,idp), zamit(idp,idp)
     & , zbmit(idp,idp), zamr(idp,idp), zbmr(idp,idp), zami(idp,idp)
     & , zbmi(idp,idp), zomega(idp), zgamma(idp)
c
      real ztemp1, zepsmach, zepsqrt, zdg, zgamax, gtein, gthin
     & , gtzin, gnein, gnhin, gnzin, gnsin, zerreal, zerimag, zerrmax
     & , zfnz, zfns, ztol, zone, zflxph, zflxnh, zflxpe, zflxnz
     & , zflxpz, zphsph, zphsnh, zphspe, zphsnz, zphspz, zft
     & , zreal, zimag, ztemp2, ztauh, ztauz, ekyrhoin
     & , ekparlin, wexb, zflz, czin, azin
     & , fnsin, q, shear, hydmass, vef, zmass, zimp
     & , kappa, aspinv, fnzin, betaein, zeni, zbetae
     & , zepstz, zepsts, zepsth, zepste, kiq, kxq, zgts, eq
     & , zepsnz, zepsne, zepsnh, zetaz, zetah, zetae, zepsmin
     & , tauhin, tauzin, zgth,zgte, zgnz, zgns, zgtz, zgnh, zgne
     & , ztwo, zhalf, ftrapein
     & , thte, tzte
c
      real  hcn, rot, tol
c
      real zfnh, zphsnet, zflxnet, zphsnef, zflxnef
     &  , zphstet, zflxtet, zphstef, zflxtef
     &  , zphsth, zflxth, zphstz, zflxtz, zflxpi
c
      complex zalpha(idp), zevec(idp,idp)
c
c
c
c...Calculation starts
c
        zepsmin = 1.e-10
c
        idim = idp
c
        zone   = 1.0
        ztwo   = 2.0
c
        zhalf = zone / ztwo
c
        zepsmach = zhalf
c
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
      ieq = max ( 2, neq )
c
      imatrx = min ( ieq - 1, 4 )
      zzz = letain(2)-1
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, zzz )
c
c
c..Print Header
c
c        write (nout,*)
c        write (nout,*)
c     &   'Weiland-Nordman eigenvalue equations, subroutine weiland18'
c        write (nout,*)
c
c        if ( letain(6) .gt. 1 ) then
c          write (nout,*)
c     &     ' Eigenvalues computed using IMSL routine gvcrg'
c        else
c          write (nout,*)
c     &     ' Eigenvalues computed using TOMSQZ routine'
c        endif
c
c        if ( ieq .eq. 11) then
c          write (nout,*)
c          write (nout,*)
c     &     'Eleven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
c     &     F, Vp_H, Av, K, and Vp_Z'
c        elseif ( ieq .eq. 10) then
c          write (nout,*)
c          write (nout,*)
c     &     'Ten eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c     &      , F, Vp, Av, and K'
c        elseif ( ieq .eq. 9) then
c          write (nout,*)
c          write (nout,*)
c     &     'Nine eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
c     &     F, Av, K, and Vp in strong ballooning limit'
c        elseif ( ieq .eq. 8) then
c          write (nout,*)
c          write (nout,*)
c     &     'Eight eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, and Vp'
c        elseif ( ieq .eq. 7) then
c          write (nout,*)
c          write (nout,*)
c     &     'Seven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, and F'
c        elseif ( ieq .eq. 6 ) then
c          write (nout,*)
c          write (nout,*)
c     &     'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
c           elseif ( ieq .eq. 5 ) then
c          write (nout,*)
c          write (nout,*)
c     &     ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
c        elseif ( ieq .eq. 4 ) then
c          write (nout,*)
c          write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
c        elseif ( ieq .eq. 2 ) then
c           write (nout,*)
c           write (nout,*) ' Two eqns for e phi/T_e and T_H'
c        else
c          write (nout,*)
c          write (nout,*) ieq,' = ieq in sbrtn weiland18'
c          call abortb(nout
c     &     ,'the value of ieq is wrong in sbrtn weiland18')
c        endif
c
c
      nerr = 0
c
c..print header
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine weiland18flux'
        write (nout,*) '(all frequencies normalized by omega_{De})'
        write (nout,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
        write (6,*) '----------------------------------------'
c
        write (nout,108) 'gnein', gnein
        write (nout,108) 'gnhin', gnhin
        write (nout,108) 'gnzin', gnzin
        write (nout,108) 'gtein', gtein
        write (nout,108) 'gthin', gthin
        write (nout,108) 'gtzin', gtzin
        write (nout,108) 'tauhin', tauhin
        write (nout,108) 'tauzin', tauzin
        write (nout,108) 'ftrapein', ftrapein
        write (nout,108) 'fnzin', fnzin
        write (nout,108) 'czin', czin
        write (nout,108) 'azin', azin
        write (nout,108) 'fnsin', fnsin
        write (nout,108) 'betaein', betaein
        write (nout,108) 'vefin', vef
        write (nout,108) 'qin', q
        write (nout,108) 'shearin', shear
        write (nout,108) 'aspinv', aspinv
        write (nout,108) 'kappain', kappa
        write (nout,108) 'hydmass', hydmass
        write (nout,108) 'ekyrhoin', ekyrhoin
        write (nout,108) 'ekparlin', ekparlin
 108    format (1x,a8,' = ',1pe14.6,',')
c
      endif
c
c
c..check validity of input data
c
      if ( neq .lt. 2 ) call abortb (nout
     & ,' neq .lt. 2 in sbrtn Weiland18flux')
c
      if ( ndim .gt. idim ) call abortb (nout
     & ,' ndim .gt. idim in sbrtn Weiland18flux')
c
c..initialize arrays
c
      do j1=1,ndim
        omega(j1)   = 0.0
        gamma(j1)   = 0.0
        chieff(j1)  = 0.0
      enddo
c
c..set up initial gradients
c
      zgne   = gnein
      zgnh   = gnhin
      zgnz   = gnzin
      zgte   = gtein
      zgth   = gthin
      zgtz   = gtzin
c
      zimp   = max ( czin, zone )
      zmass  = max ( azin, zone )
      zfnz   = fnzin * zimp
      zfns   = fnsin
c
c..Save zgns = - R ( d n_s / d r ) / n_e
c    where n_s is the fast ion density
c
      zgns = zgne - zgnh * ( zone - zfnz - zfns ) - zgnz * zfnz
c
      if ( zfns .lt. zepsmach) zgns = 0
c
c      write (nout,109) 'gnsin', zgns
c 109  format (1x,a8,' = ',1pe14.6,', computed from quasi-neutrality')
c
c..loop over gradient matrix
c
      jd = 1
c
c..set variables
c
      do j1=1,ndim
        zomega(j1)  = 0.0
        zgamma(j1)  = 0.0
        zalpha(j1)  = ( 0.0, 0.0 )
        zalfr(j1)   = 0.0
        zalfi(j1)   = 0.0
        zbeta(j1)   = 0.0
        do j2=1,ndim
          zevec(j1,j2)   = ( 0.0, 0.0 )
          zeigenv(j1,j2) = 0.0
        enddo
      enddo
c
      do j1=1,5
        flux(j1)    = 0.0
      enddo
c
c..Since the density gradients may have changed,
c  reconstruct the normalized electron density gradient.
c
      zgne = zgnh * ( zone - zfnz - zfns ) + zgnz * zfnz + zgns
c
c  check for incompatibilities:  Weiland's definitions
c
      zft    = ftrapein
c
      zbetae = betaein
c
c      zbetae = max ( cetain(20) * betaein, zepsqrt )
c
c..end of initailization
c
        if (abs(zgnh) .gt. zepsmin) zetah      = zgth / zgnh
        if (abs(zgnz) .gt. zepsmin) zetaz      = zgtz / zgnz
        if (abs(zgne) .gt. zepsmin) zetae      = zgte / zgne
c
        if (abs(tauhin) .gt. zepsmin)  ztauh   = 1 / tauhin
        if (abs(tauzin) .gt. zepsmin)  ztauz   = 1 / tauzin
c
        if (abs(zgne) .gt. zepsmin) zepsne     = 2.0 / zgne
        if (abs(zgnh) .gt. zepsmin) zepsnh     = 2.0 / zgnh
        if (abs(zgnz) .gt. zepsmin) zepsnz     = 2.0 / zgnz
        if (abs(zgth) .gt. zepsmin) zepsth     = 2.0 / zgth
        if (abs(zgte) .gt. zepsmin) zepste     = 2.0 / zgte
        if (abs(zgtz) .gt. zepsmin) zepstz     = 2.0 / zgtz
        if (abs(zgts) .gt. zepsmin) zepsts     = 2.0 / zgts
c
        if (abs(zgnz) .gt. zepsmin) eq         = zgtz / zgnz
        if (abs(zgnz) .gt. zepsmin) kiq        = zgnh / zgnz
        if (abs(zgne) .gt. zepsmin) kxq        = zgnz / zgne
c
c..diagnostic output
c
      if ( lprint .gt. 2 ) then
        write (nout,*)
        write (nout,*) '--------------------------------------'
        write (nout,*)
        write (nout,*) ' jd = ',jd
        write (nout,*)
c
       if ( lprint .gt. 2 .or. ( lprint .gt. 0 .and. jd .eq. 1 ) ) then
c
        write (nout,*) zepsnh,' = epsnh'
        write (nout,*) zepsne,' = epsne'
        write (nout,*) zepsnz,' = epsnz'
        write (nout,*) zepsth,' = epsth'
        write (nout,*) zepste,' = epste'
        write (nout,*) zepstz,' = epstz'
        write (nout,*) zetae,' = etae'
        write (nout,*) zetah,' = etah'
        write (nout,*) ztauh,' = tauh'
        write (nout,*) ztauz,' = tauz'
        write (nout,*) zft,' = zft'
        write (nout,*) zimp,' = zimp'
        write (nout,*) zmass,' = zmass'
        write (nout,*) zfnz,' = zfnz'
        write (nout,*) zfns,' = zfns'
        write (nout,*) zflz,' = zflz'
        write (nout,*) zepsqrt,' = zepsqrt'
        write (nout,*) zepsmach,' = zepsmach'
       endif
      endif
c
c..always start with analytic eigenvalue calculation
c
      ik  = 1   ! grid index
      ist = 1   ! 1st iteration
      itc = 1   ! muller's method
c
c..compute eigenvalues and eigenvectors
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7--
C  The following call is for Vers6
c
      call weiland18disp (lprint,neq,zgne,zgnh,zgnz,zgte,zgth,zgtz,
     &      tauhin, tauzin,
     &      fnzin, zimp, zmass, zfns, zbetae, ftrapein,
     &      vef, q, shear, aspinv, kappa, ekyrhoin, hcn, wexb, rot, tol,
     &      ITL, ndim, 
     &      zomega, zgamma, zamrt, zbmrt, zamit, zbmit, zvr, 
     &      zvi, zalfr, zalfi, zbeta, nerr, IK, ITS, ITERA, em)
c
c  The following call is for original
c       call weiland18disp (letain, cetain, lprint, neq, nout
c    & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
c    & , fnzin, zimp, zmass, zfns, zbetae, ftrapein
c    & , vef, q, shear, aspinv, kappa, hydmass, ekyrhoin, ekparlin
c    & , hcn, wexb, rot, tol, itl
c    & , ndim, omega, gamma, zamrt, zbmrt, zamit, zbmit, zvr, zvi
c    & , zalfr, zalfi, zbeta, nmodes, nerr
c    & , ik, ist, itc, ITS, ITERA, em)

c
c..diagnostic output
c
      if ( lprint .gt. 0 ) then
        write (6,*)
        write (6,*) ' after calling weiland18disp, itera = ', itera
        write (6,*) ' its = ', its
        write (6,*) ' nerr = ', nerr
        write (6,*) ' zomega   zgamma'
        do j1=1,ieq
          write (6,*) zomega(j1), zgamma(j1)
        enddo
      endif
c
      if ( lprint .gt. 6 ) then
        write (6,*)
        write (6,*) ' zamrt(j1,j2)  j2 ->'
        do j1=1,ieq
          write (6,192) (zamrt(j1,j2),j2=1,ieq)
        enddo
c
        write (6,*)
c
        write (6,*) ' zamit(j1,j2)  j2 ->'
        do j1=1,ieq
          write(6,192) (zamit(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
        write (6,*) ' zbmrt(j1,j2)  j2->'
        do j1=1,ieq
          write (6,192) (zbmrt(j1,j2),j2=1,ieq)
        enddo
c
       write (6,*)
       write(6,*) ' zbmit(j1,j2)  j2->'
       do j1=1,ieq
         write(6,192) (zbmit(j1,j2),j2=1,ieq)
       enddo
c
 192  format (1p5e12.4)
      endif
c
c..compute the complex eigenvalues from results of NAG routines
c
        do j=1,ieq
          zalpha(j) = cmplx( zalfr(j), zalfi(j) )
        enddo
c
        do j2=1,ieq
          do j1=1,ieq
              zevec(j1,j2) = cmplx ( zvr(j1,j2), zvi(j1,j2) )
          enddo
        enddo
c
c..save growth rates and frequencies
c
        zgamax = 0.0
        do j=1,ieq
          omega(j) = zomega(j)
          gamma(j) = zgamma(j)
          zgamax = max ( zgamax, zgamma(j) )
        enddo
        if ( zgamax .lt. zepsqrt ) go to 90
c
c..check the eigenfunctions
c
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
c
c  Real and imaginary parts
c
      zerrmax = 0.0
c
        do j=1,ieq
c
          do j1=1,ieq
c
            ztempa(j1) = 0.0
            ztempa(j1) = 0.0
c
            do j2=1,ieq
c
              zerreal =
     &            zamrt(j1,j2) * real(zevec(j2,j))
     &          - zamit(j1,j2) * aimag(zevec(j2,j))
     &          - zbmrt(j1,j2) * real(zevec(j2,j))  * zomega(j)
     &          + zbmrt(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * real(zevec(j2,j))  * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
              zerimag =
     &            zamrt(j1,j2) * aimag(zevec(j2,j))
     &          + zamit(j1,j2) * real(zevec(j2,j))
     &          - zbmrt(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
     &          - zbmrt(j1,j2) * real(zevec(j2,j))  * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          - zbmit(j1,j2) * real(zevec(j2,j))  * zomega(j)
c
              ztempa(j1) = ztempa(j1) + zerreal
              ztempb(j1) = ztempb(j1) + zerimag
c
            enddo
c
            zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
c
          enddo
c
          if ( lprint .gt. 12 ) then
            write (nout,*)
            write (nout,*) ' LHS - RHS for j =  ',j
            do j1=1,ieq
              write (nout,142) ztempa(j1), ztempb(j1)
            enddo
          endif
c
        enddo
c
        ztemp1 = max ( zepsqrt, 10.0*ztol )
        if ( zerrmax .gt. ztemp1 ) nerr = max ( nerr, 1 )
        if ( lprint .gt. 0 ) then
          write (nout,*)
          write (nout,*) zerrmax,' = zerrmax'
          write (nout,*) nerr
     &      ,' = nerr, error in eigenvalue in weiland18flux'
        endif
c
 142  format (1p10e12.4)
c
c..end of checking eigenfunctions
c
      endif
c
!| 
!| The heat and particle fluxes are computed directly from the
!| eigenvectors in the following way:
!| Consider, for example, the flux of hydrogen particles produced by
!| the perturbed $ {\mathbf E} \times {\mathbf B} $ motion of the plasma
!| \[ \Gamma_H = 2 {\rm Re} ( \tilde{n}_H \tilde{v}_E^* )
!| = 2 ( {\rm Re} \tilde{n}_H {\rm Im} \tilde{\phi}
    !| - {\rm Im} \tilde{n}_H {\rm Re} \tilde{\phi} ) k_y / B. \]
!| Using the quasi-linear assumption, the saturation level is
!| approximated by
!| \[ \hat{\phi} \equiv \frac{e \tilde{\phi}}{T_e}
!| \approx \frac{1}{k_x \rho_{sH}} \frac{ \gamma }{ k_y c_{sH} }
!| = \frac{2}{R k_x} \frac{\gamma}{\omega_{De}}. \]
!| Hence, the hydrogen particle flux is given by
!| \[ \frac{ \Gamma_H }{ n_H } \frac{ R k_x^2 }{ \omega_{De} }
 !| = 2 \frac{ \gamma^2 }{ \omega_{De}^2}
 !| \frac{ {\rm Re} \hat{n}_H {\rm Im} \hat{\phi}
      !| - {\rm Im} \hat{n}_H {\rm Re} \hat{\phi} }{ | \hat{\phi} |^2 }. \]
!| Correspondingly, the flux of heat flow through the hydrogen ions with
!| pressure $ n_H T_H $ is given by
!| \[ \frac{ \Gamma_{p_H} }{ n_H T_H } \frac{ R k_x^2 }{ \omega_{De} }
 !| = 2 \frac{ \gamma^2 }{ \omega_{De}^2}
 !| \frac{ {\rm Re} ( \hat{n}_H + \hat{T}_H ) {\rm Im} \hat{\phi}
 !| - {\rm Im} ( \hat{n}_H + \hat{T}_H ) {\rm Re} \hat{\phi} }{
   !| | \hat{\phi} |^2 }. \]
!| Four channels of transport are used in our simulations --- electron and
!| ion heat flux, as well as hydrogenic and impurity charged particle flux.
!| 
c
c..compute effective diffusivities directly from eigenvalues
c
        nmodes = 0
c
        do j=1,ieq
c
          ztemp1 =  real(zevec(1,j)) *  real(zevec(1,j))
     &           + aimag(zevec(1,j)) * aimag(zevec(1,j))
c
c..Use only those eigenvalues with positive growth rate
c  and non-zero perturbed potential (ztemp1 > 0)
c  nmodes counts the number of such modes
c
          if ( zgamma(j) .gt. zepsqrt
     &      .and. ztemp1 .gt. zepsqrt
     &       ) then
c
            nmodes = nmodes + 1
c
!| 
!| 
c
            if ( ieq .eq. 10 ) then
c
c..ieq .eq. 10 in the routine disp10 from Weiland
c  in this case Weiland changed the order of the perturbed variables to
c  e\phi/T_e, \hat{T_i}, \hat{n}_{et}, \hat{n}_{ef}, \hat{T}_{et}
c    \hat{T}_{ef}, \hat{n}_{eZ}, \hat{T}_{Z}, F, e A_\parallel / T_e
c
c..fraction of hydrogenic ions
c
              zfnh = zone - zfnz - zfns
c
              if ( zfnh .lt. zepsmin ) then
                write ( nout, * ) 'zfnh = ', zfnh
                call abortb ( nout, 'abort from sbrtn weiland18flux')
              endif
c
c..trapped electron particle flux
c
              zphsnet = - ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxnet = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..passing electron particle flux
c
              zphsnef = - ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxnef = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..impurity particle flux
c
              zphsnz = - ( aimag(zevec(7,j)) * real(zevec(1,j))
     &          - real(zevec(7,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxnz = - 2.0 * (abs( zgamma(j) ))**2
     &        * ( aimag(zevec(7,j)) * real(zevec(1,j))
     &          - real(zevec(7,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              flux(4) = flux(4) + zflxnz
c
c..hydrogenic particle flux
c
              zflxnh = ( zft * zflxnet + ( zone - zft ) * zflxnef
     &                   - zfnz * zflxnz ) / zfnh
c
              flux(2) = flux(2) + zflxnh
c
c..trapped electron temperature flux
c
              zphstet = - ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxtet = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..passing electron temperature flux
c
              zphstef = - ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxtef = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(4,j)) * real(zevec(1,j))
     &          - real(zevec(4,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..electron heat flux
c
              zflxpe = zflxnet + zft * zflxtet
     &          + ( zone - zft ) * zflxtef
c
              flux(3) = flux(3) + zflxpe
c
c..hydrogenic temperature flux
c
              zphsth = - ( aimag(zevec(2,j)) * real(zevec(1,j))
     &          - real(zevec(2,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxth = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(2,j)) * real(zevec(1,j))
     &          - real(zevec(2,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..impurity temperature flux
c
              zphstz = - ( aimag(zevec(8,j)) * real(zevec(1,j))
     &          - real(zevec(8,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxtz = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(8,j)) * real(zevec(1,j))
     &          - real(zevec(8,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c..ion heat flux, including hydrogenic and impurity ions
c
              zflxpi = ( zfnh * ( zflxnh + zflxth )
     &                 + zfnz * ( zflxnz + zflxtz ) )
     &                 / ( zfnh + zfnz )
c
              flux(1) = flux(1) + zflxpi
c
c
            else
c
c..ieq .ne. 10
c  compute fluxes as before August 2000
c  assume eigenvectors are arranged in the order of
c  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
            zreal =  real(zevec(2,j)) +  real(zevec(3,j))
            zimag = aimag(zevec(2,j)) + aimag(zevec(3,j))
c
            zphsph = - ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
            zflxph = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
            flux(1) = flux(1) + zflxph
c
c            if ( lprint .gt. 7 ) then
c              write (nout,*)
c              write (nout,*) zflxph,' = zflxph for j = ',j
c            endif
c
            zphsnh = - ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c

            zflxnh = - 2.0 * ( abs ( zgamma(j) ) )**2
     &        * ( aimag(zevec(3,j)) * real(zevec(1,j))
     &          - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
            flux(2) = flux(2) + zflxnh
c
c            if ( lprint .gt. 7 ) then
c              write (nout,*)
c              write (nout,*) zflxnh,' = zflxnh for j = ',j
c            endif
c
            if ( ieq .gt. 3 ) then
c
              zreal =  real(zevec(4,j))
     &          + ( zone - zfnz - zfns ) * real(zevec(2,j))
     &          + zfnz * real(zevec(5,j))
              zimag = aimag(zevec(4,j))
     &          + ( zone - zfnz - zfns ) * aimag(zevec(2,j))
     &          + zfnz * aimag(zevec(5,j))
c
c  Note, the electron heat flux is reduced by the fraction of
c    trapped electrons
c
              zphspe = - zft * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxpe =
     &          - 2.0 * zft * (abs( zgamma(j) ))**2
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              flux(3) = flux(3) + zflxpe
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxpe,' = zflxpe for j = ',j
c              endif
            endif
c
            if ( ieq .gt. 4 ) then
c
              zphsnz = - ( aimag(zevec(5,j)) * real(zevec(1,j))
     &          - real(zevec(5,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              zflxnz = - 2.0 * (abs( zgamma(j) ))**2
     &        * ( aimag(zevec(5,j)) * real(zevec(1,j))
     &          - real(zevec(5,j)) * aimag(zevec(1,j)) ) / ztemp1
c
              flux(4) = flux(4) + zflxnz
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxnz,' = zflxnz for j = ',j
c              endif
            endif
c
            if ( ieq .gt. 5 ) then
c
              zreal =  real(zevec(5,j)) +  real(zevec(6,j))
              zimag = aimag(zevec(5,j)) + aimag(zevec(6,j))
c
              zphspz = - ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              zflxpz = - 2.0 * (abs( zgamma(j) ))**2
     &        * ( zimag * real(zevec(1,j))
     &          - zreal * aimag(zevec(1,j)) ) / ztemp1
c
              flux(5) = flux(5) + zflxpz
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxpz,' = zflxpz for j = ',j
c              endif
            endif
c
c..end of choice of number of equations
c
            endif
c
c..header for diagnostic printout of frequencies, fluxes, and phases
c
        if ( letain(29) .gt. 0 .and.  jd .eq. 1 ) then
          write (nout,134)
c
 134  format (//' Diagnostic printout of frequencies, fluxes and phases'
     &  /' Note: fluxph = flux of hydrogenic thermal energy'
     &  ,' = 2.0 * gamma^2 * phaseph'
     &  /9x,' (phaseph is related to the phases of the perturbations)'
     &  /7x,'fluxnh = flux of hydrogenic ions = 2.0 * gamma^2 * phasenh'
     &  /7x,'fluxpe = flux of electron thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepe'
     &  /7x,'fluxnz = flux of impurity ions = 2.0 * gamma^2 * phasenz'
     &  /7x,'fluxpz = flux of impurity thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepz'
     &  //1x,'radius',t10,'omega',t20,'gamma'
     &  ,t30,'fluxph',t40,'phaseph',t50,'fluxnh',t60,'phasenh'
     &  ,t70,'fluxpe',t80,'phasepe',t90,'fluxnz',t100,'phasenz'
     &  ,t110,'fluxpz',t120,'phasepz  #m')
c
c..diagnostic printout of frequencies and fluxes mode by mode
c
        write (nout,135) cetain(29), zomega(j), zgamma(j)
     &    , zflxph, zphsph, zflxnh, zphsnh, zflxpe, zphspe
     &    , zflxnz, zphsnz, zflxpz, zphspz
        endif
c
 135  format (0pf7.3,1p12e10.2,' #m')
c
          endif
c
        enddo
c
c..compute effective total diffusivities
c
        chieff(1) = flux(1)
     &   / sign ( max ( abs ( zgth ), zepsqrt ),  zgth )
        chieff(2) = flux(2)
     &   / sign ( max ( abs ( zgnh ), zepsqrt ),  zgnh )
        chieff(3) = flux(3)
     &   / sign ( max ( abs ( zgte ), zepsqrt ),  zgte )
        chieff(4) = flux(4)
     &   / sign ( max ( abs ( zgnz ), zepsqrt ),  zgnz )
        chieff(5) = flux(5)
     &   / sign ( max ( abs ( zgtz ), zepsqrt ),  zgtz )
c
      if ( lprint .gt. 2 ) then
c
c..print eigenvalues and eigenfunctions
c
        write (nout,121)
        do j=1,ieq
          write (nout,122) zomega(j), zgamma(j)
     &     , zalfr(j), zalfi(j), zbeta(j)
        enddo
 121    format (/' Solution of the eigenvalue equations'
     &   /t4,'zomega',t18,'zgamma',t32,'zalfr',t46,'zalfi',t60,'zbeta')
 122    format (1p5e14.5,i5)
c
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,ieq-1)
c
      endif
c
      if ( lprint .gt. 99 ) then
c
        write (nout,*)
        write (nout,*) ' Eigenvectors zeigenv(j1,j2) j2->'
        do j1=1,ieq
          write (nout,124) (zeigenv(j1,j2),j2=1,ieq)
 124      format (1p12e11.3)
        enddo
c
        write (nout,*)
        write (nout,*) ztol,' = tolerance used in sbrtn f02gjf'
c
      endif
c
c..end of loop over diffusivity matrix
c
  50  continue
c
  90  continue
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn weiland18flux'
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,imatrx)
c
      endif
c
c
      return
c
 130    format (/t3,'chi_i',t16,'dif_h',t29,'chi_e'
     &    ,t42,'dif_Z',t55,'chi_Z')
c
 132    format (1p11e12.4)
c
      end

!| 
!| \end{document}

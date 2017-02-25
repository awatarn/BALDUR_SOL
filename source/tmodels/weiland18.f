      subroutine weiland18 (
     &   letain,    cetain,    lprint,    neq,     nout,    gnein
     & , gnhin,     gnzin,     gtein,     gthin,   gtzin,   tauhin
     & , tauzin,    fnzin,     czin,      azin,    fnsin,   betaein
     & , ftrapein,  vef,       q,         shear,   aspinv,  kappa
     & , hydmass,   ekyrhoin,  ekparlin,  ndim,    omega,   gamma
     & , difthi,    velthi,    chieff,    nmodes,  nerr,    ITS
     & , ITERA,     em)
c
      implicit none
c
      integer idp
c
      parameter ( idp = 15 )
c
      integer neq, ieq, imatrx, j, j1, j2, jd, nerr, nout, ndim
     & , lprint, zzz, nmodes, idim, iter(idp), ifail, ITERA, ITS
     & , em
c
      real cetain(32), letain(32), perform(idp)
c
      real omega(idp), gamma(idp), difthi(ndim,idp), velthi(idp)
     & , chieff(idp), zalfr(idp), zalfi(idp), zflxm(idp,idp)
     & , zchim(idp,idp), zgm(idp,idp), zvr(idp,idp), zvi(idp,idp)
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
     & , ekparlin, zflz, czin, azin
     & , fnsin, q, shear, hydmass, vef, zmass, zimp
     & , kappa, aspinv, fnzin, betaein, zeni, zbetae
     & , zepstz, zepsts, zepsth, zepste, kiq, kxq, zgts, eq
     & , zepsnz, zepsne, zepsnh, zetaz, zetah, zetae, zepsmin
     & , tauhin, tauzin, zgth,zgte, zgnz, zgns, zgtz, zgnh, zgne
     & , ztwo, zhalf, ftrapein
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
c     &   'Weiland-Nordman eigenvalue equations, subroutine etaw18'
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
c          write (nout,*) ieq,' = ieq in sbrtn etaw18'
c          call abortb(nout
c     &     ,'the value of ieq is wrong in sbrtn etaw18')
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
     & 'Weiland-Nordman eigenvalue equations, subroutine etaw18'
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
c..check validity of input data
c
      if ( neq .lt. 2 ) call abortb (nout
     & ,' neq .lt. 2 in sbrtn Weiland18')
c
      if ( ndim .gt. idim ) call abortb (nout
     & ,' ndim .gt. idim in sbrtn Weiland18')
c
c..initialize arrays
c
      do j1=1,ndim
        omega(j1)   = 0.0
        gamma(j1)   = 0.0
        chieff(j1)  = 0.0
        velthi(j1)  = 0.0
        perform(j1) = 0.0
        do j2=1,ndim
          difthi(j1,j2) = 0.0
          zchim(j1,j2)  = 0.0
          zflxm(j1,j2)  = 0.0
        enddo
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
c..set up gradient matrix zgm(j1,jd)
c
c  zdg = the finite difference used to construct the zgm matrix
c
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 0.01
c
c  input values for zgm matrix
c
      do jd=1,imatrx+1
        zgm(1,jd) = zgth
        zgm(2,jd) = zgnh
        zgm(3,jd) = zgte
        if ( ieq .gt. 4 ) then
          zgm(4,jd) = zgnz
          zgm(5,jd) = zgtz
        endif
      enddo
c
c  incremental values
c  Note:  if zgm(jd,jd) is negative, increment in negative direction
c
      do jd=1,imatrx+1
        zgm(jd,jd+1) = zgm(jd,jd+1) + sign ( zdg, zgm(jd,jd+1) )
      enddo
c
      if ( lprint .gt. 8 ) then
c
        write (nout,*)
        write (nout,*) '  zgm(j1,jd) matrix,  -> jd'
        do j1=1,ieq-1
          write (nout,132) (zgm(j1,jd),jd=1,ieq)
        enddo
c
      endif
c
c..loop over gradient matrix
c
      do 50 jd=1,imatrx+1
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
c..compute gradient scale lengths from the zgm(j1,jd) matrix
c
      zgth = zgm(1,jd)
      zgnh = zgm(2,jd)
      zgte = zgm(3,jd)
      zgnz = zgm(4,jd)
      zgtz = zgm(5,jd)
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
      zbetae = max ( cetain(20) * betaein, zepsqrt )
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
      call etaw18 (letain, cetain, lprint, neq
     & , zepsne, zepsnh, zepsnz
     & , zepste, zepsth, zepstz, ztauz, ztauh
     & , fnzin, zimp, zmass, zfns, zbetae, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass, ekyrhoin, ekparlin
     & , ndim, zomega, zgamma, zamrt, zbmrt, zamit, zbmit, zvr, zvi
     & , zalfr, zalfi, zbeta, nmodes, nerr, ITS, ITERA, em)
c
c..diagnostic output
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
c..save growth rates and frequencies during first element of jd loop
c
      if ( jd .eq. 1 ) then
        zgamax = 0.0
        do j=1,ieq
          omega(j) = zomega(j)
          gamma(j) = zgamma(j)
          zgamax = max ( zgamax, zgamma(j) )
        enddo
        if ( zgamax .lt. zepsqrt ) go to 90
      endif
c
c..check the eigenfunctions
c
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
      endif
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
          write (nout,*) nerr,' = nerr, error in eigenvalue in etaw18'
        endif
c
 142  format (1p10e12.4)
c
c..compute effective diffusivities directly from eigenvalues
c  assume eigenvectors are arranged in the order of
c  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
        nmodes = 0
c
        do j=1,ieq
c
          ztemp1 =  real(zevec(1,j)) *  real(zevec(1,j))
     &           + aimag(zevec(1,j)) * aimag(zevec(1,j))
c
          if ( zgamma(j) .gt. zepsqrt
     &      .and. ztemp1 .gt. zepsqrt
     &       ) then
c
            nmodes = nmodes + 1
c
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
            zflxm(1,jd) = zflxm(1,jd) + zflxph
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
            zflxm(2,jd) = zflxm(2,jd) + zflxnh
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
              zflxm(3,jd) = zflxm(3,jd) + zflxpe
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
              zflxm(4,jd) = zflxm(4,jd) + zflxnz
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
              zflxm(5,jd) = zflxm(5,jd) + zflxpz
c
c              if ( lprint .gt. 7 ) then
c                write (nout,*)
c                write (nout,*) zflxpz,' = zflxpz for j = ',j
c              endif
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
        zchim(1,jd) = zflxm(1,jd)
     &   / sign ( max ( abs ( zgth ), zepsqrt ),  zgth )
        zchim(2,jd) = zflxm(2,jd)
     &   / sign ( max ( abs ( zgnh ), zepsqrt ),  zgnh )
        zchim(3,jd) = zflxm(3,jd)
     &   / sign ( max ( abs ( zgte ), zepsqrt ),  zgte )
        zchim(4,jd) = zflxm(4,jd)
     &   / sign ( max ( abs ( zgnz ), zepsqrt ),  zgnz )
        zchim(5,jd) = zflxm(5,jd)
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
        write (nout,132) (zchim(j1,jd),j1=1,ieq-1)
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
c..save effective diffusivities
c
      do j1=1,min(ieq-1,5)
        chieff(j1) = zflxm(j1,1) /
     &    sign( max( abs(zgm(j1,1)), zepsqrt), zgm(j1,1) )
      enddo
c
      if ( imatrx .lt. 1 ) go to 90
c
c..computation of difthi(j1,j2) from zflxm(j1,jd) and zgm(j1,jd)
c
      do j1=1,imatrx
c
        ztemp2 = 0.0
c
        do jd=1,imatrx
c
          ztemp1 = zgm(jd,jd) - zgm(jd,1)
          if ( abs(ztemp1) .lt. zepsqrt )
     &       ztemp1 = zdg * sign(zone,zgm(jd,1))
c
          difthi(j1,jd) = ( zflxm(j1,jd+1) - zflxm(j1,1) ) / ztemp1
          ztemp2 = ztemp2 + difthi(j1,jd) * zgm(jd,1)
c
        enddo
c
cap        letain(7) = 0
c
        if ( letain(7) .eq. 0 ) then
c
c..compute convective velocities
c
          velthi(j1) = zflxm(j1,1) - ztemp2
c
        else
c
c..alternatively, set the convective velocities to zero
c    and rescale the diffusivity matrix
c
          velthi(j1) = 0.0
          if ( abs(ztemp2) .lt. zepsqrt )
     &      ztemp2 = sign ( zepsqrt, ztemp2 )
          do jd=1,imatrx
            difthi(j1,jd) = difthi(j1,jd) * zflxm(j1,1) / ztemp2
          enddo
c
        endif
c
      enddo
c
c..diagnostic printout
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ' matrix zflxm(j1,jd)  -> jd'
        do j1=1,ieq-1
          write (nout,132) (zflxm(j1,jd),jd=1,ieq-1)
        enddo
c
      endif
c
  90  continue
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn etaw18'
c
        write (nout,*)
        write (nout,*) ' diffusion matrix from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
        write (nout,*)
        write (nout,*) ' difthi(j1,j2)  -> j2'
        do j1=1,imatrx
          write (nout,132) (difthi(j1,j2),j2=1,imatrx)
        enddo
c
        write (nout,*)
        write (nout,*) ' convective velocities'
     &    ,' normalized by omega_{De} / k_y'
        do j1=1,imatrx
          write (nout,132) velthi(j1)
        enddo
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


c@clprof for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cprofl/
     &      nxiprf,           xizoni(kpgrid),   xibndi(kpgrid)
     & ,    uj0zen(kpgrid),   cuj0zn(kpgrid*5)
     & ,    dpnorm(kpgrid),   cdpnrm(kpgrid*5)
     & ,    b0uzet(kpgrid),   cb0uzt(kpgrid*5)
     & ,    rinv2(kpgrid),    crinv2(kpgrid*5)
     & ,    qprofl(kpgrid),   cqprfl(kpgrid*5)
c
c
c  nxiprf         number of grid points on which profiles are prescribed
c                       including ghost points
c                       = mzones in BALDUR
c
c  xizoni(j)      independent variable xi on BALDUR zone centers
c                       j=1,nxiprf
c                       the magnetic axis (xi=0.) lies between j=1 and 2
c
c  xibndi(j)      independent variable xi on BALDUR zone boundaries
c                       j=1,nxiprf
c                       xi(1) is a ghost point.
c                       The magnetic axis is at xi(2) = 0.0
c
c  uj0zen(j)      \mu_0 J^{0 \zeta} / B^{0 \zeta}
c                       on xizoni(j) grid
c                       not explicitly including the effects of magnetic
c                       islands
c
c  cuj0zn(j,js)   cubic spline interpolation coefficients
c                 for uj0zen(j) with respect to xizoni(j)
c                 computed in sbrtn island.
c
c                 At any flux surface labeled xi,
c                       with xizoni(j) .le. xi .lt. xizoni(j+1)
c                 \mu_0 J^{0 \zeta} / B^{0 \zeta} =
c                       uj0zen(j) + (xi-xizoni(j)) * ( cuj0zn(j,1)
c                                 + (xi-xizoni(j)) * ( cuj0zn(j,2)
c                                 + (xi-xizoni(j)) * ( cuj0zn(j,3) ) ) )
c
c  dpnorm(j)      2 \mu_0 dp^0(\xi)/d\xi  / ( (B^{0\zeta))^2 )
c                       on xibndi(j) grid
c
c  cdpnrm(j,js)   cubic spline interpolation coefficients
c                 for dpnorm(j) with respect to xizoni(j).
c
c  b0uzet(j)      B^{0 \zeta} = <1/R^2> R B_toroidal
c                       on xizoni(j) grid
c
c  cb0uzt(j,js)   cubic spline interpolation coefficients
c                 for b0uzet(j) with respect to xizoni(j).
c
c  rinv2(j)       < 1 / R^2 > on xizoni(j) grid
c
c  crinv2(j,js)   cubic spline interpolation coefficients
c                 for rinv2(j) with respect to xizoni(j)
c
c  qprofl(j)      q-value = B^{0 \zeta} / B^{0 \theta}
c                       on xibndi(j) grid
c
c  cqprfl(j,js)   cubic spline interpolation coefficients
c                 for qprofl(j) with respect to xibndi(j).

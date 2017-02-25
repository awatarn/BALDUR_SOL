c@clmetr for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /hametr/
     & nxiham,     nhamad,   nithta
     &,xihamd(kpgrid),       hamjac(kpgrid)
     &,harcj0(kpgrid),harcjm(kpgrid,kpharm),harsjm(kpgrid,kpharm)
     &,haycj0(kpgrid),haycjm(kpgrid,kpharm),haysjm(kpgrid,kpharm)
     &,hamadr(kpgrid,kpthta),hamady(kpgrid,kpthta),hamadz(kpgrid,kpthta)
     &,hdrdth(kpgrid,kpthta),hdydth(kpgrid,kpthta),hdzdth(kpgrid,kpthta)
     &,hdrdxi(kpgrid,kpthta),hdydxi(kpgrid,kpthta),hdzdxi(kpgrid,kpthta)
     &,hgxixi(kpgrid,kpthta),hgthth(kpgrid,kpthta),hgzeze(kpgrid,kpthta)
     &,hgxith(kpgrid,kpthta),hgxize(kpgrid,kpthta),hgthze(kpgrid,kpthta)
     &,hcxixi(kpgrid, 0:kpharm) ,hcthth(kpgrid, 0:kpharm)
     &,hczeze(kpgrid, 0:kpharm) ,hsxith(kpgrid, 0:kpharm)
     &,hsxize(kpgrid, 0:kpharm) ,hcthze(kpgrid, 0:kpharm)
c
c  nxiham         number of flux surfaces on which
c                       Hamada coordinate metric elements are computed
c                       = nxiequ at present
c
c  nhamad         number of hamada harmonics desired
c
c  nithta         number of angles or Fourier haromonics
c                       around each flux surface
c
c  xihamd(j)      xi values for those flux surfaces on which
c                       Hamada coordinate metric elements are computed
c
c  hamjac(j)      Hamada Jacobian = xi * hamjac(j)
c
c  harcj0(jx)     Harmonic representation of Hamada coordinates
c  harcjm(jx,jm)
c  haysjm(jx,jm)
c
c  hamadr(jx,jt)  major radius R as a function of xihamd(jx), theta(jt)
c
c  hamady(jx,jt)  Y              as a function of xihamd(jx), theta(jt)
c
c  hamadz(jx,jt)  Z=\zeta - \phi as a function of xihamd(jx), theta(jt)
c                 where \zeta is the Hamada variable corresponding to
c                       \phi  which is the toroidal angle
c
c  hdrdth(jx,jt)  d R / d theta  as a function of xihamd(jx), theta(jt)
c
c  hdydth(jx,jt)  d Y / d theta  as a function of xihamd(jx), theta(jt)
c
c  hdzdth(jx,jt)  d Z / d theta  as a function of xihamd(jx), theta(jt)
c
c  hdrdxi(jx,jt)  d R / d xi     as a function of xihamd(jx), theta(jt)
c
c  hdydxi(jx,jt)  d Y / d xi     as a function of xihamd(jx), theta(jt)
c
c  hdzdxi(jx,jt)  d Z / d xi     as a function of xihamd(jx), theta(jt)
c
c
c       The following metric elements in Hamada coordinates are stored
c  as a function of xihamd(jx), theta(jt):
c
c  hgxixi(jx,jt) =   ( d R / d xi )^2 + ( d Y / d xi )^2 
c
c  hgthth(jx,jt) = [ ( d R / d theta )^2 + ( d Y / d theta )^2 ] / xi^2
c
c  hgzeze(jx,jt) = 1 / R^2
c
c  hgxith(jx,jt) = [  (d R / d xi)(d R / d theta)
c                   + (d Y / d xi)(d Y / d theta) ] / xi
c
c  hgxize(jx,jt) = d Z / d xi
c
c  hgthze(jx,jt) = ( d Z / d theta ) / xi
c
c       These variables are used to transform from
c  B^{1\xi}, B^1_\theta, B^1_{\zeta}
c  to  B^1_{\xi}, B^{1 \theta}, B^{1 \zeta}
c  as in Eqs. (21 - 23) of Phys. Fluids 29 (1986) 753--761.
c
c       The transformation equations are:
c
c  [ ( d R / d theta )^2 + ( d Y / d theta )^2 ]  B^{\theta}
c            = - [  (d R / d xi)(d R / d theta)
c                 + (d Y / d xi)(d Y / d theta) ] B^{\xi}
c              + B_{\theta} + ( d Z / d \theta )  B_{\zeta}
c
c  B^{\zeta} = ( d Z / d \xi ) B^{\xi} + ( d Z / d \theta ) B^{\theta}
c              + ( 1 / R^2 ) B_{\zeta}
c
c  B^{\xi}   = [ ( d R / d xi )^2 + ( d Y / d xi )^2 ] B^{\xi}
c              + [  (d R / d xi)(d R / d theta)
c                 + (d Y / d xi)(d Y / d theta) ] B^{\theta}
c              - ( d Z / d \xi ) B_{\zeta}
c
c  hcxixi(jx,jm) are the cos harmonics of hgxixi(jx,jt)
c  hcthth(jx,jm) are the cos harmonics of hgthth(jx,jt)
c  hczeze(jx,jm) are the cos harmonics of hgzeze(jx,jt)
c  hsxith(jx,jm) are the sin harmonics of hgxith(jx,jt)
c  hsxize(jx,jm) are the sin harmonics of hgxize(jx,jt)
c  hcthze(jx,jm) are the cos harmonics of hgthze(jx,jt)
c
c       These variables are harmonics of the matrix elements needed
c  to transform from  B^{1v}, B^1_\theta, B^1_{\zeta}
c  to  B^1_v, B^{1 \theta}, B^{1 \zeta}
c  as in Eqs. (21 - 23) of Phys. Fluids 29 (1986) 753--761.
c  They are computed in sbrtn metevl.
c
c  xma(j,m)             : 1.0 / grad v dot grad v
c
c  xmb(j,m)             : grad v dot grad htheta / grad v dot grad v
c
c  xmc(j,m)             : grad v dot grad hzeta  / grad v dot grad v
c
c  xmd(j,m) * xi**2     : [ grad htheta dot grad htheta
c               - ( grad htheta dot grad v )**2 / grad v dot grad v ]
c
c  xme(j,m)             : [ grad htheta dot grad hzeta
c               - ( grad htheta dot grad v ) * ( grad v dot grad hzeta)
c               / grad v dot grad v ]
c
c  xmf(j,m)             : [ grad hzeta dot grad hzeta
c               - ( grad hzeta dot grad v )**2 / grad v dot grad v

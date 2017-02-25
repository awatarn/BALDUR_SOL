c@clpert for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cbpert/ b1uxi(kpmode),  b1utr(kpmode),  b1uzr(kpmode)
     & ,              b1lxi(kpmode),  b1ltr(kpmode),  b1lzr(kpmode)
     & ,              sgnbux(kpmode), xib1a(kpgrid),  ndiftr(kpgrid)
     &,b1uxia(kpgrid,kpmode),b1utra(kpgrid,kpmode),b1uzra(kpgrid,kpmode)
     &,b1lxia(kpgrid,kpmode),b1ltra(kpgrid,kpmode),b1lzra(kpgrid,kpmode)
     & ,bmxprt(2*kpmode), nxib1a, ndifco
     &,trm201(kpgrid,kpmode),trm202(kpgrid,kpmode),trm203(kpgrid,kpmode)
     &,trm204(kpgrid,kpmode), tp1mn(kpgrid,kpmode)
     &,tdlbzi(kpgrid),tdj0dx(kpgrid),tdp0dx(kpgrid),thamji(kpgrid)
     &,tuj0zi(kpgrid)
     &,tbpert(kpgrid,2*kpmode),tdbpdx(kpgrid,2*kpmode)
c
c       Contravariant and covariant components of the perturbed
c  magnetic field for each harmonic m.
c  b1lzr, b1lxi, b1utr, and b1uzr are computed in sbrtn difco.
c  b1uxi, b1ltr,          are computed in sbrtn diftr.
c
c  b1uxi(m) = real part of [ - i B^{1 \xi}_{m n} ]
c  b1utr(m) = real part of [ B^{1 \theta }_{m n} ]
c  b1uzr(m) = real part of [ B^{1 \zeta }_{m n} ]
c
c  b1lxi(m) = real part of [ - i B^{1}_{\xi m n} ]
c  b1ltr(m) = real part of [ B^{1}_{ \theta m n} ]
c  b1lzr(m) = real part of [ B^{1}_{ \zeta m n} ]
c
c  sgnbux(m) = initial sign of b1uxi(m) computed in sbrtn iwidth
c
c  xib1a(jx) = \xi array at which perturbed magnetic field stored
c
c  ndiftr(jx) = number of times sbrtn idiftr has been called by ode
c                 solver from \xi = deltaxi up to xib1a(jx)
c
c  b1uxia(jx,jm) = b1uxi(jm) at \xi = xib1a(jx) in sbrtn idiftr
c  b1utra(jx,jm) = b1utr(jm) ...
c    These arrays are stored for printout and diagnostic
c
c  bmxprt(jm) = max abs perturbation amplitude within q .lt. qmode(jm)
c     here: bmxprt(jm) refers to the b^1_{\theta m n} part and
c           bmxprt(jm+nmodes) refers to the B^{1x}_{mn} part
c             of the perturbation as defined in sbrtn idiftr.
c
c  nxib1a   = index of the perturbed field array, used as a counter
c               set in sbrtn iwidth
c
c  ndifco   = number of times sbrtn idiftr has been called from ode
c
c  trm201,... are terms in Eq. (20) of PF 29 (1986) 754
c
c  tbpert(jx,jm) = pb(jm) perturbation from sbrtn idiftr
c  tdbpdx(jx,jm) = pdb(jm) derivative of perturbation from sbrtn idiftr
c                    here jm=1,2*nmodes

c
      parameter(mu=6,nu=128,nv=1,nv2=1+nv/2,mnd=mu*(nv+1),n2d=4*mu+nu)
        common/scmgue/r0n(nv),z0n(nv),raxis(nv),zaxis(nv)
        common/scmind/mpol,mpol2,mpol3,mpol4,nfp,ntheta,nphi,nphi2,n2
        common/scmmna/mm(mu),nn(nv),m1(mnd),n1(mnd),
     &           dm1(mu),faccon(mu),xmpq(mu,4),mpnt
        common/scmspc/gnorm,specw,delt,deltf
        common/scmtrg/angle(nu,nv),dnorm,elongate,twopi
        common/scmxst/xvec(n2d),gvec(n2d),xdot(n2d),xstore(n2d),
     &                result(nv2,mu,4)
c

cat  > ceq2n2.m << 'EOF'
        real iota,iotas,mass,mscale,nscale,jtheta,jzeta
        common/e2current/bpco(nsd),bzco(nsd),jtheta(nsd),jzeta(nsd)
        common/e2extfld/potvac(nznt),xmpot(nznt),xnpot(nznt),           vac
     >         brv(nznt),bphiv(nznt),bzv(nznt),mpmax,bscale,ivac,ivac2  vac
        common/e2fsqu/fnorm,fsqr,fsqz,fsqr1,fsqz1,fsql,fsq,
     >         fedge,wb,wp,fsqt(100),wdot(100),equif(nsd)
        common/e2inputdat/am(0:5),ai(0:5),raxis(0:nmax),zaxis(0:nmax),
     >         rb(0:nmax,0:mpol1,2),zb(0:nmax,0:mpol1,2),ftol,gam,
     >         ncurr,nlam,nfp,niter,nstep,nvacskip
       common/e2mnarray/xrz3(0:nmax,0:mpol1),xrz4(0:nmax,0:mpol1),
     >         xlam3(0:nmax,0:mpol1),xlam4(0:nmax,0:mpol1),
     >         xmpq(0:mpol1,3),mscale(0:mpol2),nscale(0:nmax3)
        common/e2profs/phip(nrztd),iota(nrztd),torcur(nrztd),
     >         vp(nsd),mass(nsd),pres(nsd),iotas(nsd),phips(nsd)
        common/e2scalars/dnorm,hs,ohs,twopi,voli,ijacob,itfsq,iequi,
     >         irst,iter1,iter2,isigng,meven,modd,mnsold,
     >         ndamp,ns,ns4,neqs,nrzt,mns
        common/e2scalefac/faclam(2*mnd*nsd),shalf(nrztd),sqrts(nrztd),
     >         scalxc(4*mnd*nsd),wint(nrztd)
        common/e2spectra/faccon(0:nmax,0:mpol1),specw(nsd),tcon(nsd)
        common/e2time/delt,otav,otau(2*nsd),timer(0:10)
        common/e2trign/cosmu (ntheta2,0:mpol1),sinmu (ntheta2,0:mpol2),
     >                 cosmum(ntheta2,0:mpol1),sinmum(ntheta2,0:mpol1),
     >               cosmui (ntheta2,0:mpol2),cosmumi(ntheta2,0:mpol1),
     >                 cosnv (nzeta,0:nmax2),sinnv (nzeta,0:nmax2),
     >                 cosnvn(nzeta,0:nmax),sinnvn(nzeta,0:nmax)
        common/e2magfield/bsqsav(nznt,3),dbsq(nznt),bsqvac(nznt),       vac
     >         rbsq(nznt),curpol,curtor,cpol,ctor,delbsq,mvac           vac
        common/e2xstuff/xcdot(neq),xstore(neq),xc(neq),gc(neq)
        common /e2eqpr/ nradeq, eqradx(nsd), eqpres(nsd), eqiota(nsd)
     &    , eqcoef(nsd+5,3)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c

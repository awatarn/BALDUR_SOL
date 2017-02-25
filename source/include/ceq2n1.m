cat  > ceq2n1.m << 'EOF'
        parameter(nsd=31, mpol=6, nmax=0, ntheta=18, nzeta=01, nvac=0,
     >  nsd1=nsd+1, ntheta1=2*(ntheta/2), ntheta2=1+ntheta1/2, nznt=
     >  nzeta*ntheta2,nrztd = nznt*nsd+1,mpol1 = mpol-1,mpol2 = 2*mpol1,
     >  nmax1=1+nmax, nmax2=2*nmax, nmax3=nmax2+1, mnd=mpol*nmax1, mnd1=
     >  mnd-1, mnd2=2*mnd, mnmax=nmax1+mpol1*(1+2*nmax), neq=6*nsd*mnd)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c

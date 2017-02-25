c 22:00 25-apr-96 .../baldur/code/bald/deqmom2.f
c */ 23:00 27-jul-91 /11040/bald91/wbaldn1(LIB) wsequil DEQMOM2, Bateman
c */ from cfs get /20330/vmec/0190/vmec0190.2d
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c**********************************************************************c
c
c     To obtain this file, type:       (Note: ^ stands for linefeed)
c cfs get /11040/bald91/wbaldn1
c lib wbaldn1 ^ x wsequil ^ n wsequil ^ x deqmom2 ^ end
c
c     To compile these souce routines:
c First obtain the common block file DCOMMON from wbaldn1
c and obtain library YBALDLIB from wbaldn2
c and obtain controlle XTRIPX from wbaldn3, then type:
c cosmos deqmom2
c
c     The following files will be produced:
c t0pre  - all but the cosmos statements in this file
c t1pre  - same, with inline comments (!...) stripped out by XTRIPX
c t0     - same, after being precompiled
c b0     - compiled binary file
c ybaldlib - library of compiled routines
c log    - cosmos log output
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@changes
c  rgb 21-jan-92 protected agianst (...)**gam wherever gam = 0.0
c    if max pressure = 0.0, set eqpres(j) = 1.e-10 * (1 - eqradx(j)**2)
c    change test program to run eqmom2 twice, changing pressure
c    implement changes from vmec2.for from Wieland
c  rgb 11-jun-91 changing names of subroutines to all start with e2
c  rgb 26-mar-90 changed from rcft to cft77 and ldr
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c From 000251 hirshman steve@mfe (orn) on 02/02/90 at 05:20:34  e
c Subject :
c I fixed up the bcovar subroutine, as well as a few other little
c odds-and-ends.  The attached threed1 file shows that it works
c as expected.  I will finish commenting out the rest of the 3d stuff,
c and will send the new VMEC.2D file to you in the next few minutes.
c -Steve
c From 000251 hirshman steve@mfe (orn) on 02/02/90 at 06:05:27  e
c Subject : vmec.2d should be attached
c Dick:
c I commented out a bunch of 3d lines in TOTZSP, FORCES, AND TOMNSP.
c It is important to comment them ALL or NOT AT ALL, but do not
c eliminate some and not others. (Some vectors are used as scratch
c arrays and can be non-zero even when you think they are zero...)
c The file worked fine for the 2d test case (I sent earlier).
c
c Let me know if there are any problems. I also corrected the BCOVAR
c stuff - yesterday, the overall PHIP factor was missing, and also
c because of the way the angle goes only from 0 to pi, angle integrations
c MUST be done with the WINT weight function and dnorm, rather than
c the simple SSUM we used yesterday to compute psiav.
c
c -STEVE
c       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
c       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
c       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
c       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
c       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION LAMBDA,
c       WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED VARIATIONALLY
c       ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS DETERMINED BY
c       MINIMIZING <M> = m**2 S(m) , WHERE S(m) =  Rm**2 + Zm**2 .
c       AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE NO. OF R AND Z IS
c       USED TO IMPROVE RADIAL RESOLUTION. A FREE-BOUNDARY OPTION IS
c       AVAILABLE (FOR nvac > 0), WITH A USER-SUPPLIED SUBROUTINE
c       "VACUUM" NEEDED TO COMPUTE THE PLASMA BOUNDARY VALUE OF B**2.
c
c       Added features since last edition
c       1.  Implemented preconditioning algorithm for R,Z
c       2.  The physical (unpreconditioned) residuals are used
c           to determine the level of convergence
c       3.  The lambda force is not solved at js = 2 (for nlam=1), but
c           rather a two point extrapolation law is used for improved accuracy
c
c       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
c       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
c       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
c       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c     This file contains a collection of subroutines used to
c     set up and call the equilibrium moments code
c     to compute flux surface averages and
c     to interface with the 1 1/2 D BALDUR upgrade
c
c sbrtn   originally  comment
c -----   -------     -------
c eqmom2              interface routine
c eqfixa  fixaray
c e2solv  eqsolve
c e2prof  profil3d
c e2volv  evolve
c e2bcov  bcovar
c e2prec  precondn
c e2lamc  lamcal
c e2boun  boundary
c e2curr  current
c e2qfor  eqfor
c e2scal  scalfor
c e2trid  trid
c e2getf  getfsq
c e2forc  forces
c e2fnct  funct3d
c e2xtrp  extrap
c e2alia  alias
c e2ntrp  interp
c e2bss   bss
c e2jacb  jacobian
c e2half  mhalf
c e2polr  polar
c e2pres  pressure
c e2prnt  printout
c e2rsid  residue
c e2strt  restart
c e2spct  spectrum
c e2tomn  tomnsp
c e2totz  totzsp
c e2trig  triggf
c e2rout  wrout
c e2conv  convert
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@eqmom2 /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
c
      subroutine eqmom2 ( nrad, nmom, r0, rm, ym
     &  , nradxi, eqxi, epress, eiota, raxin, gamma
     &  , r0b, rmb, ymb, kmdim, krdim, ftolin, nitmax, nskip
     &  , nradin, ninit, ierflag)
c
      include 'ceq2n1.m'
      include 'ceq2n2.m'
      include 'ceq2n3.m'
      include 'ceq2n4.m'
c
        real gtt(nrztd),gtz(nrztd),gzz(nrztd),rmnc(mnmax),
     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),rax(nznt),
     >  zax(nznt),work1(12*nrztd),work2(12*nrztd)
c
      dimension  r0(krdim), rm(kmdim,krdim), ym(kmdim,krdim)
     &  , eqxi(*)         , epress(*)      , eiota(*)
     &  , rmb(kmdim)      , ymb(kmdim)
     &  , ztempr(mpol)    , ztempy(mpol)   , ztempl(mpol)
c
        character*80 werror(0:6)
cbate        character*20 date,mach,timeloc
        data intflag/0/
c
        external case9
c
c
c  Argument list:
c  --------------
c
c  nrad   = number of radiad grid points from the magnetic axis (j=1)
c             to the edge (j=nrad)
c  nmom   = number of poloidal harmonics
c  r0(jr), rm(jm,jr), ym(jm,jr) = output values of harmonics
c  nradxi = number of radial points in equilbrium profiles (input)
c  eqxi(jr)   = equilibrium profile radial grid (input)
c  epress(jr) = equilibrium scalar pressure (input)
c  eiota(jr)  = equilibrium rotational transform (input)
c  torflx = total toroidal flux surface within last flux surface
c             in webers (input)
c  gamma  = ratio of specific heats (input)
c  r0b, rmb(jm), ymb(jm) = boundary values of r0(nrad), rm(jm,nrad)
c             and ym(jm,nrad) (input)
c  kmdim  = first dimension of rm, ym, rmb, and ymb (input)
c  krdim  = second dimension of rm, ym,
c             and dimension of epress, eiota (input)
c  ftolin = value of fsq = fr**2 + fq**2 at which iteration ends (input)
c  nitmax = maximum number of iterations allowed (input)
c  nskip  = number of iterations between calls to convergence check (input)
c  nradin = number of radial grid points for the initial iteration
c  ninit    must be 0 on first call to initialize, reset to 1 thereafter
c  ierflag = error flag (output)
c
c  -----------------------------
c
c  Given pressure and rotational transform as a function of toroidal flux,
c this code computes the shape of the equilibrium flux surfaces in the form:
c
c  R(xi,theta) = r0(j) + sum of rm(m,j)*cos(m*theta) for m=1 to mom
c  Y(xi,theta) =         sum of ym(m,j)*sin(m*theta) for m=1 to mom
c
c where r0, rm, and ym are found on equally spaced intervals of
c
c  xi = sqrt (normalized toroidal flux)
c
c... all spatial dimensions are in meters ...
c%
c
c..assign input data to variables in common/e2inputdat/
c
      ncurr  = 0
      nlam   = 1
      nfp    = 1
      niter  = nitmax
      nstep  = nskip
      nvacskip = 3
      intflag = 0
! cap      
      ierflag = 0
c
      ftol   = ftolin
      gam    = gamma
c
        nradeq = nradxi
        zmaxpr = 0.0
c
      do 22 j=1,nradeq
        eqradx(j) = eqxi(j)
        eqpres(j) = epress(j)
        zmaxpr = max ( zmaxpr, epress(j) )
        eqiota(j) = eiota(j)
  22  continue
c
      if ( zmaxpr .eq. 0.0 ) then
        do 21 j=1,nradeq
          eqpres(j) = 1.e-10 * ( 1. - eqradx(j)**2 )
  21    continue
      endif
c
      raxis(0) = raxin
      zaxis(0) = 0.0
      rb(0,0,1) = r0b
      rb(0,0,2) = 0.0
c
      do 24 jm=1,mpol1
        rb(0,jm,1) = rmb(jm)
        rb(0,jm,2) = 0.0
        zb(0,jm,1) = 0.0
        zb(0,jm,2) = ymb(jm)
  24  continue
c
c               1.3  COMPUTE INVARIANT ARRAYS (IF VMEC IS USED
c               AS A SUBROUTINE, CALL FIX_ARAY ONLY ONCE AT SETUP)
c
        if ( ninit .lt. 1 )  call e2fixa(ierflag)
        if ( ierflag .ne. 0 ) return
c
c               1.4  MAKE CALL TO EQUILIBRIUM SOLVER
c
        if ( nradin .lt. 2 .or. nradin .ge. nsd ) nradin = nsd
        call e2solv(nradin,intflag,ierflag)
c
c               1.4a PERFORM COARSE-TO-FINE MESH INTERPOLATION
c
        if ( nradin .lt. nsd .and. ierflag .eq. 0 ) then
          intflag = nradin
          call e2solv(nsd,intflag,ierflag)
        end if
c
c..convert to harmonic form
c
      do 32 js=1,ns
        call e2conv ( ztempr, ztempy, ztempl, xm, xn, js
     &    , xc, xc(1+mns)
     &    , xc(1+2*mns), xc(1+3*mns), xc(1+4*mns), xc(1+5*mns) )
        r0(js) = ztempr(1)
        do 31 jm=1,nmom
          rm(jm,js) = ztempr(jm+1)
          ym(jm,js) = ztempy(jm+1)
  31    continue
  32  continue
c
c..reset initialization flag
c
      ninit = 1
c
c..Diagnostic printout
c
      write (6,108) mns,neqs,nrad,nmom,kmdim,krdim
 108  format (/'  mns =',i5/,'  neqs = ',i5/,'  nrad = ',i5/
     & ,'  nmom = ',i5/,'  kmdim = ',i5/,'  krdim = ',i5/)
c
      write (6,120) (r0(j),j=1,nrad)
 120  format (' r0 = ',/(1pe14.5) )
      write (6,121)
 121  format (/' rm(jm,js) = ')
      do 84 js=1,nrad
        write (6,122) (rm(jm,js),jm=1,nmom)
  84  continue
 122  format(1p5e14.5)
      write (6,123)
 123  format (/' ym(jm,js) = ')
      do 85 js=1,nrad
        write (6,122) (ym(jm,js),jm=1,nmom)
  85  continue
c
c..error messages
c
      if ( ierflag .eq. 0 ) then
        return
      elseif ( ierflag .eq. 1 ) then
        write(6,*)' INITIAL JACOBIAN CHANGED SIGN (NEED A BETTER GUESS)'
        stop
      elseif ( ierflag .eq. 2 ) then
        write(6,*)' MORE THAN 100 JACOBIAN ITERATIONS (DECREASE DELT)'
        stop
      elseif ( ierflag .eq. 3 ) then
        write(6,*)' m IN BOUNDARY ARRAY EXCEEDS mpol-1 (INCREASE MPOL)'
        stop
      elseif ( ierflag .eq. 4 ) then
        write(6,*)' n IN BOUNDARY ARRAYS OUTSIDE nmax,(-nmax) RANGE'
        stop
      elseif ( ierflag .eq. 5 ) then
        write(6,*)' ASSUMED SIGN OF JACOBIAN (ISIGNG) IS WRONG'
        stop
      else
        write(6,*)' IERFLAG = ',ierflag,' IN SBRTN EQMOM2'
        stop
      endif
c
        return
        end

c******************
      BLOCK DATA case9
      include 'ceq2n1.m'
      include 'ceq2n2.m'
      include 'ceq2n3.m'
      include 'ceq2n4.m'
c
c        real gtt(nrztd),gtz(nrztd),gzz(nrztd),rmnc(mnmax),
c     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),rax(nznt),
c     >  zax(nznt),work1(12*nrztd),work2(12*nrztd)
        data meven/0/, modd/1/, mvac/0/
c
        data rb/mnd2*0./, zb/mnd2*0./, nlam/1/, nvacskip/3/
        data am/6*0./, ai/6*0./, raxis/nmax1*0./, zaxis/nmax1*0./
c
        data mscale/.707,mpol2*1.0/,nscale/.707,nmax3*1.0/
c
        data bscale/1.0/, delbsq/1.0/, mns/0/, ndamp/10/, ns4/25/,
     >  isigng/-1/, itfsq/0/, ivac/0/, timer/11*0./
c
        data curtor/0./, curpol/0./,
     >  iotas(1)/0./,phips(1)/0./
c
c        real gtt(nrztd),gtz(nrztd),gzz(nrztd),rmnc(mnmax),
c     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),rax(nznt),
c     >  zax(nznt),work1(12*nrztd),work2(12*nrztd)
c
        end

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2fixa /11040/bald91/wbaldn1  wsequil  DEQMOM2
c rgb 30-may-96 do j=1,mpol2+1 --> do j=1,mpol2, etc
c
        subroutine e2fixa(ierflag)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        data rtest /0./, ztest/0./
c
c               3.0   INDEX OF LOCAL VARIABLES
c
c         mscale  array for norming trig functions (internal use only)
c         nscale  array for norming trig functions (internal use only)
c
c               3.1  COMPUTE TRIGONOMETRIC FUNCTION ARRAYS
c
        twopi=8.*atan(1.)
        dnorm=1./float(nzeta*(ntheta2-1))
c
      isigng = -1
      mscale(0) = 0.707
      do j=1,mpol2
        mscale(j) = 1.0
      end do
c
      nscale(0) = 0.707
      do j=1,nmax3
        nscale(j) = 1.0
      end do
c
      dmu0 = 1.256637e-6
c
      do 10 i = 1,ntheta2
        argi = twopi*float(i-1)/float(ntheta1)
        do 10 m = 0,mpol2
          arg = argi*float(m)
          cosmui(i,m) = cos(arg)*mscale(m)
          if(i.eq.1.or.i.eq.ntheta2)cosmui(i,m) = .5*cosmui(i,m)
          sinmu (i,m) = sin(arg)*mscale(m)
          if(m.gt.mpol1) go to 10
          cosmu (i,m) = cos(arg)*mscale(m)
          sinmum(i,m) =-sinmu(i,m)*float(m)
          cosmum(i,m) = cosmu(i,m)*float(m)
          cosmumi(i,m)= cosmui(i,m)*float(m)
 10   continue
c
      do 20 j = 1,nzeta
        argj = twopi*float(j-1)/float(nzeta)
        if (nmax2 .eq. 0) then
          n = 0
          arg = argj*float(n)
          cosnv (j,n) = cos(arg)*nscale(n)
          sinnv (j,n) = sin(arg)*nscale(n)
          if ( n .gt. nmax ) go to 20
          cosnvn(j,n) = cosnv(j,n)*float(n*nfp)
          sinnvn(j,n) =-sinnv(j,n)*float(n*nfp)
          else
        do n = 0,nmax2
          arg = argj*float(n)
          cosnv (j,n) = cos(arg)*nscale(n)
          sinnv (j,n) = sin(arg)*nscale(n)
          if ( n .gt. nmax ) go to 20
          cosnvn(j,n) = cosnv(j,n)*float(n*nfp)
          sinnvn(j,n) =-sinnv(j,n)*float(n*nfp)
          enddo
          endif
 20       continue
c
      do 30 m = 0,mpol1
        power = -.5*float(m)
        xmpq(m,1) = max0(0,m**2-1)
        xmpq(m,2) = m**4
        xmpq(m,3) = m**5
        do 30 n = 0,nmax
          t1 = -.20*isigng*dnorm/((1.+ m)**3*max0(m,1))
          if ( m .eq. 0  .or.  n .eq. 0 )
     &      faccon(n,m)=.5*t1/(mscale(m)*nscale(n))**2
          if ( m .ne. 0  .and.  n .ne. 0 )
     &      faccon(n,m) = t1/(mscale(m)*nscale(n))**2
          xlam3(n,m) = 2.*3.**power
          xlam4(n,m) = -5.**power
          xrz3(n,m) = 2.*2.**power
          xrz4(n,m) = -  3.**power
  30  continue
c               3.3  CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS ISIGNG)
      rtest = 0.0
      ztest = 0.0
      do 40 n = 0,nmax
        rtest = rtest + rb(n,1,1)
 40     ztest = ztest + zb(n,1,2)
c
        if( (rtest*ztest*float(isigng)) .ge. 0.) ierflag = 5
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2solv /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2solv(nsval,intflag,ierflag)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real :: mclock
c
c       ____________________________________________________________
c
c               2.0   INDEX OF LOCAL VARIABLES
c
c         hs      radial mesh size increment
c         iequi   counter used to call EQFOR at end of run
c         ijacob  counter for number of times jacobian changes sign
c         irst    counter monitoring sign of jacobian; resets R, Z, and
c                 Lambda when jacobian changes sign and decreases time step
c         isigng  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)
c         iterj   stores position in main iteration loop (j=1,2)
c         itfsq   counter for storing FSQ into FSQT for plotting
c         ivac    counts number of free-boundary iterations
c         ndamp   number of iterations over which damping is averaged
c         meven   parity selection label for even poloidal modes of R and Z
c         modd    parity selection label for odd poloidal modes of R and Z
c         nvac    fixed boundary (=0) or free boundary (>0)
c         rb      boundary coefficient array for R (rcc,rss)
c         zb      boundary coefficient array for Z (zcs,zsc)
c         gc      stacked array of R, Z, Lambda Spectral force coefficients
c         xc      stacked array of scaled R, Z, Lambda Fourier coefficients
c         STACKING ORDER:
c           1:mns   => rmncc Fourier coefficients
c           1+mns,2*mns => rmnss Fourier coefficients
c           1+2*mns,3*mns => zmncs Fourier coefficients
c           1+3*mns,4*mns => zmnsc Fourier coefficients
c           1+4*mns,5*mns => lmncs Fourier coefficients
c           1+5*mns,neqs => lmnsc Fourier coefficients
c       ________________________________________________________________
c % (check to see if these should be reinitialized each call)
c
      bscale = 1.0
      delbsq = 1.0
      mns    = 0
      ndamp  = 10
      ns4    = 25
      isigng = -1
      itfsq  = 0
      ivac   = 0
! cap
      w1     = 0.
      r01    = 0.
      res1   = 0.
!      
      do j=0,10
        timer(j)  = 0.0
      end do
c
        fsql = 0.
        ns = nsval
        hs = 1./float(ns-1)
        ohs = 1./hs
        mnsold = mns
        mns = ns*mnd
        neqs = 6*mns
        nrzt = nznt*ns
        iter2=1
        irst=1
        iequi = 0
        write(6,*)'NS = ',ns,' NO. FOURIER MODES = ',mnmax
c
c               2.1  COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
c               AND STORE XC, XCDOT FOR POSSIBLE RESTART
c
        call e2prof(xc,xc(1+2*mns),intflag)
        delt = 1.0
        ijacob=-1
 10     call e2strt
        ijacob=ijacob+1
        iter1=iter2
c
c               2.2  FORCE ITERATION LOOP
c
        do 30 miter=iter1,niter
c
c               2.3  ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
c
        iter2=miter
c3d     if(miter.le.2)iter2=miter-1
        call e2volv(ierflag)
        if(ierflag.ne.0)return
c
c               2.4  COMPUTE ABSOLUTE STOPPING CRITERION
c
        r00 = xc(1)*mscale(0)*nscale(0)
        w0=wb+wp/(gam-1.)
        wdota=abs(w0-w1)/w0
        r0dot=abs(r00-r01)/r00
        r01=r00
        w1=w0
        if( (fsqr.le.ftol.and.fsqz.le.ftol.and.fsql.le.ftol
     >      .and.miter-iter1.gt.10)
     >      .or.(ns.lt.nsd.and.(fsqr+fsqz+fsql).le.2.e-8))goto 40
        if(ijacob.ge.100.or.miter.eq.niter)go to 40
        if(ivac.eq.1)then                                              ! vac
           print 110, miter,1./bscale                                  ! vac
           write(6,110)miter,1./bscale                                 ! vac
           ivac = ivac+1                                               ! vac
        endif                                                          ! vac
        if(mod(miter,(niter/100+1)).ne.0.or.(ns.lt.nsd))go to 15
        itfsq=itfsq+1
        fsqt(itfsq) = fsqr + fsqz
        wdot(itfsq)=max(wdota,1.e-13)
c
c               2.5  TIME STEP CONTROL
c
 15     if(miter.eq.iter1)res0=fsq
        res0=min(res0,fsq)
        if( (fsq.le.res0) .and. ((iter2-iter1).gt.nsd/2) )then
          call e2strt
        else if(fsq.gt.100.*res1.and.miter.gt.iter1)then
          irst = 2
        else if(fsq.gt.(1.5*res1).and.mod(miter-iter1,5).eq.0.and.
     >  (miter-iter1).gt.2*ns4)then
          irst = 3
        endif
        if(irst.ne.1)goto 10
        res1=res0
c
c               2.6  PRINTOUT EVERY NSTEP ITERATIONS
c
 20     if(mod(miter,nstep).ne.0.and.miter.gt.1)go to 30
        call e2prnt(iter2,w1,r00)
 30     continue
 40     call e2prnt(iter2,w0,r00)
cbate        call second(timer(0))
cap
        timer(0) = 0.01 * mclock()
        if(ijacob.ge.100)ierflag = 2
        if(ns.lt.nsd)write(6,100)timer(0),ijacob
        write(6,60)wdota,r0dot
 60     format(/,' d(ln W)/dt = ',1pe10.3,' d(ln R0)/dt = ',1pe10.3/)
 100    format(/,'  COMPUTATIONAL TIME FOR INITIAL INTERPOLATION = ',
     >  1pe10.2,' SECONDS ',' JACOBIAN RESETS = ',i4)
 110    format(/'  VACUUM PRESSURE TURNED ON AT ',i4, ' ITERATIONS'/,   ! vac
     >  '  MAGNETIC FIELD IN PLASMA IS SCALED AT END BY ',1pe10.2,      ! vac
     >  ' (VAC. UNITS)'/)                                               ! vac

        return
        end

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2prof /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2prof(rmn,zmn,intflag)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real rmn(ns,0:nmax,0:mpol1,2),zmn(ns,0:nmax,0:mpol1,2)
c
        data phiedge/9.370000/, dmu0/1.256637e-6/
c
        dmass(x) = am(0)+x*(am(1)+x*(am(2)+x*(am(3)+x*(am(4)+x*am(5)))))
        diota(x) = ai(0)+x*(ai(1)+x*(ai(2)+x*(ai(3)+x*(ai(4)+x*ai(5)))))
c       ________________________________________________________________
c
c             4.0   INDEX OF LOCAL VARIABLES
c
c         ai        vector of coefficients in phi-series for iota
c         am        vector of coefficients in phi-series for mass
c         iota      rotational transform, two-dimensional array on half-grid
c         mass      mass profile on half-grid
c         sqrts     sqrt(s), two-dimensional array on full-grid
c         phip      radial derivative of phi/(2*pi) on half-grid
c         phiedge  value of real toroidal flux at plasma edge (s=1)
c         phips(iotas)  same as phip (iota), one-dimensional array
c         pressure  pressure profile on half-grid, mass/phip**gam
c         shalf     sqrt(s) ,two-dimensional array on half-grid
c         jzeta   toroidal current inside flux surface (vanishes like s)
c       ________________________________________________________________
c
c             4.0   VERIFY SIGN OF CURPOL CONSISTENT WITH TOROIDAL FLUX
c
        if( phiedge*curpol .lt. 0 )then
        print *,' CHANGING SIGN OF PHIEDGE '
        write(6,*)' CHANGING SIGN OF PHIEDGE '
        phiedge = -phiedge
        endif
c
c             4.1   COMPUTE PHIP, MASS, CHIP PROFILES ON HALF-GRID
c
      do 10 js=2,ns
        phij = hs*(js-1.5)
        shalf(js) = sqrt(hs*abs(js-1.5)) * sign(1.0,js-1.5)
        phips(js) = isigng * phiedge / twopi
        jzeta(js)=curtor*(1. - (1.-phij)**2)
cbate        mass(js) = dmu0*dmass(phij)*(abs(phips(js))*rb(0,0,1))**gam
cbate        iotas(js) = diota(phij)
  10  continue
c
c..interpolate profiles
c
        call cubint (eqradx,eqpres,nradeq,0,eqcoef,nradeq+5
     &    ,shalf,mass,ns,0,0.0,1
     &    ,'cubint (eqradx,eqpres,... in sbrtn e2prof')
c
        call cubint (eqradx,eqiota,nradeq,0,eqcoef,nradeq+5
     &    ,shalf,iotas,ns,0,0.0,1
     &    ,'cubint (eqradx,eqiota,... in sbrtn e2prof')
c
c..convert from pressure to mass array
c
      if ( gam .eq. 0.0 ) then
        do 11 js=2,ns
          mass(js) = dmu0 * mass(js)
  11    continue
      else
        do 12 js=2,ns
          mass(js) = dmu0 * mass(js)
     &      * (abs(phips(js))*rb(0,0,1))**gam
  12    continue
      endif
        mass(1) = mass(2)
c
        do 15 js=1,ns
          if(ncurr.eq.1)iotas(js)=0.
          do 15 lk=1,nznt
            loff = js + ns*(lk-1)
            shalf(loff) = sqrt(hs*abs(js-1.5)) * sign(1.0,js-1.5)
            sqrts(loff) = sqrt(hs*(js-1))
            torcur(loff)=isigng*jzeta(js)/twopi
            phip(loff)=phips(js)
            iota(loff)=iotas(js)
  15    continue
c
        lk = 0
        do 20 lt=1,ntheta2
          do 20 lz=1,nzeta
            lk = lk+1
            do 20 js=2,ns
              wint(js+ns*(lk-1)) = cosmui(lt,0)/mscale(0)
  20    continue
c
        do 25 lk=1,nznt
          wint(1+ns*(lk-1)) =  0.
  25    continue
c
c..for testing purposes, print out profile arrays
c
      write (6,110)
 110  format (/t3,'js',t9,'shalf(js)',t23,'mass(js)'
     &  ,t35,'iotas(js) in sbrtn e2prof from interpolation')
c
      do 14 js=1,ns
        write (6,112) js, shalf(js), mass(js), iotas(js)
 112    format (2x,i3,1p3e13.4)
  14  continue
c
c             4.2   COMPUTE INITIAL R & Z PROFILES FROM SCALED BOUNDARY DATA
c
      do 30 l=1,neqs
        gc(l) = xc(l)
        xc(l)=0.
        xcdot(l)=0.
  30  continue
c
      do 35 l = 1,4*mns
        gc(l) = gc(l) * scalxc(l)
  35  continue
c
c             4.3 COMPUTE INITIAL VALUES FOR R AND Z FOURIER COEFFICIENTS
c
      do 40 ntype = 1,2
        do 40 m=0,mpol1
          do 40 n=0,nmax
            t1 = 1./(mscale(m)*nscale(n))
            do 40 js=1,ns
              sm0 = 1. - (hs*(js-1))
              l = js + ns*(n + nmax1*m)
              if(mod(m,2).eq.1)scalxc(l)=1./max(sqrts(js),sqrts(2))
              if(mod(m,2).eq.0)scalxc(l)=1.
              if(m.eq.0.and.ntype.eq.1)then
                rmn(js,n,m,1) = (rb(n,m,1)+(raxis(n)-rb(n,m,1))*sm0)*t1
                zmn(js,n,m,1) = (zb(n,m,1)+(zaxis(n)-zb(n,m,1))*sm0)*t1
              else if(m.eq.0.and.ntype.eq.2)then
                rmn(js,n,m,2) = 0.
                zmn(js,n,m,2) = 0.
              else if(m.ne.0)then
                facj = t1*sqrts(js)**m
                rmn(js,n,m,ntype) = rb(n,m,ntype)*facj
                zmn(js,n,m,ntype) = zb(n,m,ntype)*facj
              endif
 40   continue
c
      do 50 l = 1+mns,4*mns
 50   scalxc(l) = scalxc(l-mns)
c
c             4.4   INTERPOLATE FROM COARSE TO FINE RADIAL GRID
c
        if ( intflag .ne. 0 ) call e2ntrp(xc,gc,scalxc,intflag)
c
        return
        end

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2volv /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2volv(ierflag)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
c
c             6.1   COMPUTE MHD FORCES
c
        call e2fnct
        if(iter2.ne.0.or.irst.ne.2)go to 10
        ierflag = 1
        return
c
c             6.2   COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
c                   R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
c
 10     if(iter2.ne.iter1)go to 15
        do 35 i=1,ndamp
 35     otau(i)=.15/delt
 15     fsq1 = fsqr1 + fsqz1
        if(iter2.gt.iter1)dtau=min(abs(log(fsq/fsq1)),0.15)
        fsq = fsq1
        if(iter2.le.1)return
        do 25 i=1,ndamp-1
 25     otau(i)=otau(i+1)
        if(iter2.gt.iter1)otau(ndamp)=dtau/delt
        otav=ssum(ndamp,otau,1)/ndamp
        dtau=delt*otav
        b1=1.-.5*dtau
        fac=1./(1.+.5*dtau)
        do 20 l=1,neqs
        xcdot(l)=fac*(xcdot(l)*b1+delt*gc(l))
 20     xc(l)=xc(l)+xcdot(l)*delt
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2strt /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2strt
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
c
      go to (10,20,20) irst
c
  10  call scopy(neqs,xc,1,xstore,1)
        return
c
  20  continue
        do 30 l=1,neqs
          xcdot(l)=0.
          xc(l)=xstore(l)
  30    continue
c
        delt = delt*(0.96*(irst-2)+0.90*(3-irst))
        irst = 1
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2fnct /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2fnct
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        include 'ceq2n3.m'
        include 'ceq2n4.m'
        real gtt(nrztd),gtz(nrztd),gzz(nrztd),rmnc(mnmax),
     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),rax(nznt),
     >  zax(nznt),work1(12*nrztd),work2(12*nrztd)
        real :: mclock
        equivalence(work1,azodd),(work2,r)
c
c             8.0   EXTRAPOLATE M>2 MODES AT JS = 2
c
        call e2xtrp(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns),
     >  xc(1+4*mns),xc(1+5*mns),xrz3,xrz4,xlam3,xlam4,ns,nlam)
c
      do 5 l = 1,4*mns
 5    gc(l) = xc(l) * scalxc(l)
c
c             8.1   INVERSE FOURIER TRANSFORM TO S,THETA,ZETA SPACE
c                   R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE
c
        call e2totz(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),xc(1+4*mns),
     >  xc(1+5*mns),r,rt,rz,z,zt,zz,lt,lz,rcon,zcon,work1,meven)
c
        call e2totz(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),xc(1+4*mns),
     >  xc(1+5*mns),rodd,rtodd,rzodd,zodd,ztodd,zzodd,lt,lz,
     >  gtt,gzz,work1,modd)
c
c             8.2   COMPUTE CONSTRAINT FORCE (GCON) AND SCALE LT,LZ
c
      do 10 l=1,nrzt
        lt(l)   = phip(l)*lt(l)
        lz(l)   = phip(l)*lz(l)
        rt0(l)  = rt(l) + rtodd(l)*sqrts(l)
        zt0(l)  = zt(l) + ztodd(l)*sqrts(l)
        rcon(l) = rcon(l) + gtt(l)*sqrts(l)
        zcon(l) = zcon(l) + gzz(l)*sqrts(l)
        gcon(l) = 0.
 10   armn(l) = (rcon(l)-rcon0(l))*rt0(l)+(zcon(l)-zcon0(l))*zt0(l)
c
        if(iter2.eq.1.and.nvac.eq.0) call scopy(nrzt,rcon,1,rcon0,1)
        if(iter2.eq.1.and.nvac.eq.0) call scopy(nrzt,zcon,1,zcon0,1)
        if(iter2.gt.1) call e2alia(gcon,armn,work1,gc,gc(1+mns))
c
c             8.3a  COMPUTE S AND THETA DERIVATIVE OF R AND Z
c                   AND JACOBIAN ON HALF-GRID
c
        call e2jacb(r,z,rt,zt,rodd,zodd,rtodd,ztodd,armn,azmn,
     >  brmn,bzmn,crmn,arodd,brodd,wint,shalf,ohs,nrzt,nznt,irst)
c
        if ( irst .eq. 2 ) return
c
c             8.4   COMPUTE PRESSURE AND VOLUME ON HALF-GRID
c
        call e2pres(crmn,czmn,wint,wp,dnorm,mass,vp,pres,
     >  gam,nznt,ns)
c
        if ( iter2 .eq. 1 ) voli = (twopi**2)*hs*ssum(ns-1,vp(2),1)
c
c             8.5   COMPUTE COVARIANT COMPONENTS OF B, MAGNETIC
c                   PRESSURE, AND METRIC ELEMENTS ON HALF-GRID
c
        call e2bcov(clmn,blmn,crmn,czmn,arodd,bzmn,brmn,
     >  azmn,armn,brodd,gtt,gtz,gzz,azodd,bzodd,czodd,work1)
c
        bz0 = sdot(nznt,blmn(ns),ns,wint(ns),ns)
c
c             8.6   COMPUTE VACUUM MAGNETIC PRESSURE AT PLASMA EDGE
c
c         NOTE: FOR FREE BOUNDARY RUNS, THE PLASMA VOLUME CAN BE INCREASED
c               BY DECREASING CURPOL AT FIXED PHIPS.  THE VALUE FOR CURPOL
c               GIVEN BELOW SHOULD BE USED AS AN INITIAL GUESS, WHICH CAN BE
c               CHANGED BY READING IN CURPOL FROM DATA STATEMENT

      if ( nvac .eq. 0  .or.  iter2 .le. 1 ) go to 55
cbate        call second(timeon)
cap
        timeon = 0.01 * mclock()
        ivac2=mod(iter2-iter1,nvacskip)                                 ! vac
        if( (fsqr+fsqz).le.5.e-1.and.ivac2.eq.0 )ivac=ivac+1            ! vac
        if( (fsqr+fsqz).lt.10.*ftol .and. ivac.gt.1)mvac=1              ! vac
        cpol = curpol                                                   ! vac
        if(ivac.eq.1.and.cpol.eq.0.)cpol = twopi*dnorm*bz0              ! vac
          if(ivac2.eq.0.and.ivac.ge.1.and.mvac.eq.0)then                ! vac
        if(mod(ivac-1,10).eq.0.and.(fsqr+fsqz).gt.2.e-7)                ! vac
     >  ctor=twopi*dnorm*( 1.5*sdot(nznt,clmn(ns),ns,wint(ns),ns)       ! vac
     >                  -  0.5*sdot(nznt,clmn(ns-1),ns,wint(ns-1),ns) ) ! vac
        call e2conv(rmnc,zmns,lmns,xm,xn,ns,xc,xc(1+mns),               ! vac
     >  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))                ! vac
        do 20 lk = 1,nznt                                               ! vac
        rax(lk) = r(1+ns*(lk-1))                                        ! vac
 20     zax(lk) = z(1+ns*(lk-1))                                        ! vac
cbate        call vacuum(rmnc,zmns,xm,xn,mnmax,ctor,cpol,bsqvac,rax,zax)!     vac
          endif                                                         ! vac
        do 40 lk=1,nznt                                                 ! vac
        bsqsav(lk,3) = 1.5*(czmn(ns*lk)   - pres(ns))                   ! vac
     >               - 0.5*(czmn(ns*lk-1) - pres(ns-1))                 ! vac
        rbsq(lk) = bsqvac(lk)*ohs*(r(ns*lk) + rodd(ns*lk))              ! vac
 40     dbsq(lk)=abs(bsqvac(lk)-bsqsav(lk,3))                           ! vac
        if(ivac.eq.1)call scopy(nznt,czmn(ns),ns,bsqsav(1,1),1)         ! vac
        if(ivac.eq.1)call scopy(nznt,bsqvac,1,bsqsav(1,2),1)            ! vac
cbate        call second(timeoff)
cap
        timeoff = 0.01 * mclock()
        timer(1) = timer(1) + (timeoff-timeon)
c
c             8.7   COMPUTE AVERAGE EQUILIBRIUM PROPERTIES AT END OF RUN
c
  55  if ( iequi .eq. 1 ) then
c
        call e2qfor(clmn,blmn,bpco,bzco,jtheta,jzeta,equif,
     >  pres,phips,vp,iotas,czmn,specw,ohs,voli,twopi,bscale,
     >  xc,xc(1+2*mns),wint,dnorm,mscale,nscale,ns,isigng)
c 
        call e2rout(czmn,crmn,clmn,blmn,brodd,azodd,bzodd,czodd)
c
        return
      endif
c
c             8.8  COMPUTE MHD FORCES ON INTEGER-MESH
c
        call e2forc(rbsq,gtt,gtz,gzz,sqrts,shalf,ohs,nrzt,ns,ivac)
c
c             8.9  FOURIER-TRANSFORM MHD FORCES TO (M,N)-SPACE
c
      do 75 l = 1,neqs
 75   gc(l) = 0.
c
        call e2tomn(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),
     >  gc(1+5*mns),armn,brmn,crmn,azmn,bzmn,czmn,blmn,clmn,
     >  rt0,zt0,rcon,zcon,work2,meven)
c
      do 80 l=2,nrzt
        rt0(l) =  rt0(l)*sqrts(l)
        zt0(l) =  zt0(l)*sqrts(l)
        rcon(l) = rcon(l)*sqrts(l)
 80     zcon(l) = zcon(l)*sqrts(l)
c
        call e2tomn(gc,gc(1+mns),gc(1+2*mns),gc(1+3*mns),gc(1+4*mns),
     >  gc(1+5*mns),arodd,brodd,crodd,azodd,bzodd,czodd,blmn,clmn,
     >  rt0,zt0,rcon,zcon,work2,modd)
c
      do 90 l = 1,4*mns
 90   gc(l) = gc(l) * scalxc(l)
c
c             8.10  COMPUTE FORCE RESIDUALS
c
        call e2rsid(gc,gc(1+2*mns),gc(1+4*mns),bz0,work2)
c
        return
        end
c******************

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2xtrp /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2xtrp(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,
     >  x3,x4,xl3,xl4,ns,nlam)
c
        include 'ceq2n1.m'
        real rmncc(ns,0:mnd1),rmnss(ns,0:mnd1),x3(0:mnd1),
     >       zmncs(ns,0:mnd1),zmnsc(ns,0:mnd1),x4(0:mnd1),
     >       lmncs(ns,0:mnd1),lmnsc(ns,0:mnd1),xl3(0:mnd1),xl4(0:mnd1)
c
      do 10 mn = 2*nmax1,mnd1
        rmncc(2,mn) = x3(mn)*rmncc(3,mn) + x4(mn)*rmncc(4,mn)
        rmnss(2,mn) = x3(mn)*rmnss(3,mn) + x4(mn)*rmnss(4,mn)
        zmncs(2,mn) = x3(mn)*zmncs(3,mn) + x4(mn)*zmncs(4,mn)
        zmnsc(2,mn) = x3(mn)*zmnsc(3,mn) + x4(mn)*zmnsc(4,mn)
  10  continue
c
      if( nlam.eq.0 )return
      do 20 mn = 0,mnd1
        lmncs(2,mn) = xl3(mn)*lmncs(3,mn) + xl4(mn)*lmncs(4,mn)
        lmnsc(2,mn) = xl3(mn)*lmnsc(3,mn) + xl4(mn)*lmnsc(4,mn)
  20  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2alia /11040/bald91/wbaldn1  wsequil  DEQMOM2
c  pis 12-jun-98 do 15 n= ntheta2 --> do 15 n = 1, ntheta2
c    do 14 m = nzeta --> do 14 m = 1, nzeta
c
        subroutine e2alia(gcon,zcon,work,gcs,gsc)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real gcon(ns,nzeta,ntheta2),zcon(ns,nzeta,ntheta2),
     >  gcs(ns,0:nmax,0:mpol1),gsc(ns,0:nmax,0:mpol1),
     >  work(ns,nzeta,0:mpol1,4)

c        do 10 l = 1,ns*mnd
c        gcs(l,0,0) = 0.
c 10     gsc(l,0,0) = 0.
c        do 15 l = 1,nrzt
c 15     gcon(l,1,1) = 0.
c        do 20 l = 1,4*ns*nzeta*mpol
c 20     work(l,1,0,1) = 0.
        do 12 n = 0, mpol1      
        do 11 m = 0, nmax
        do 10 l = 1,ns
        gcs(l,m,n) = 0.
        gsc(l,m,n) = 0.
 10     continue
 11     continue
 12     continue
        do 15 n = 1, ntheta2   !cpis do15n= ntheta2  !!!
        do 14 m = 1, nzeta     !cpis do14m= nzeta
        do 13 l = 1,ns
        gcon(l,m,n) = 0.
 13   continue
 14     continue
 15     continue
        do 20 n = 1,4
        do 19 m = 0, mpol1
        do 18 l = 1, nzeta
        do 17 k = 1, ns
        work(k,l,m,n) = 0.
 17     continue
 18     continue
 19     continue
 20     continue

c
c       BEGIN DE-ALIASING (TRUNCATION OF GCON IN FOURIER-SPACE)
c
        do 60 m = 1,mpol1-1
        do 25 i = 1,ntheta2
        do 25 jk = 1,ns*nzeta
        work(jk,1,m,01) = work(jk,1,m,01) + zcon(jk,1,i)*cosmui(i,m)
 25     work(jk,1,m,02) = work(jk,1,m,02) + zcon(jk,1,i)*sinmu (i,m)
        do 30 n = 0,nmax
        fm = faccon(n,m)
        do 30 k = 1,nzeta
        do 30 js= 2,ns
        gcs(js,n,m) =gcs(js,n,m) +fm*tcon(js)*work(js,k,m,01)*sinnv(k,n)
 30     gsc(js,n,m) =gsc(js,n,m) +fm*tcon(js)*work(js,k,m,02)*cosnv(k,n)
c
c       RECONSTRUCT DE-ALIASED GCON
c
        do 40 n = 0,nmax
        do 40 k = 1,nzeta
        do 40 js= 2,ns
        work(js,k,m,03) = work(js,k,m,03) + gcs(js,n,m)*sinnv(k,n)
 40     work(js,k,m,04) = work(js,k,m,04) + gsc(js,n,m)*cosnv(k,n)
        do 50 i = 1,ntheta2
        do 50 jk= 1,ns*nzeta
 50     gcon(jk,1,i) = gcon(jk,1,i) + work(jk,1,m,03)*cosmu(i,m)
     >                              + work(jk,1,m,04)*sinmu(i,m)
 60     continue
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2totz /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2totz(rmncc,rmnss,zmncs,zmnsc,lmncs,lmnsc,
     >  r,rt,rz,z,zt,zz,lt,lz,rcon,zcon,work,mparity)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
c
        real rmncc(ns,0:nmax,0:mpol1),rmnss(ns,0:nmax,0:mpol1),
     >       zmncs(ns,0:nmax,0:mpol1),zmnsc(ns,0:nmax,0:mpol1),
     >       lmncs(ns,0:nmax,0:mpol1),lmnsc(ns,0:nmax,0:mpol1),
     >  r(ns,nzeta,ntheta2),rt(ns,nzeta,ntheta2),
     >  rz(ns,nzeta,ntheta2),z(ns,nzeta,ntheta2),
     >  zt(ns,nzeta,ntheta2),zz(ns,nzeta,ntheta2),
     >  lt(ns,nzeta,ntheta2),lz(ns,nzeta,ntheta2),
     >  rcon(ns,nzeta,ntheta2),zcon(ns,nzeta,ntheta2)
        dimension work(ns,nzeta,0:mpol1,12),jmin(0:mpol)
cbate        data jmin/1,1,mpol1*2/
c
      jmin(0) = 1
      jmin(1) = 1
      do j=2,mpol1+1
        jmin(j) = 2
      end do
c
c       THIS PROGRAM ASSUMES THE FOLLOWING STACKING OF R, Z, L ARRAYS:
c       rmncc(ns,0:nmax,0:mpol1),rmnss,zmncs,zmncc,lmncs,lmnsc
c
cahk      do 10 l=1,nrzt
cahk        r(l,1,1)  = 0.
cahk        rt(l,1,1) = 0.
c3d     rz(l,1,1) = 0.
cahk        z(l,1,1)  = 0.
cahk        zt(l,1,1) = 0.
c3d     zz(l,1,1) = 0.
cahk        rcon(l,1,1)=0.
cahk        zcon(l,1,1)=0.
cahk  10  continue
c
      do 12 l=1,ns
        do 11 m=1,nzeta
          do 10 n=1,ntheta2 
            r(l,m,n)  = 0.
            rt(l,m,n) = 0.
c3d         rz(l,m,n) = 0.
            z(l,m,n)  = 0.
            zt(l,m,n) = 0.
c3d         zz(l,m,n) = 0.
            rcon(l,m,n)=0.
            zcon(l,m,n)=0.
  10  continue
  11  continue
  12  continue

      if ( mparity .eq. meven ) then
c        do 15 l = 1,12*ns*nzeta*mpol
c 15     work(l,1,0,1) = 0.
c        do 20 l=1,nrzt
c        lt(l,1,1) = 1.
c 20     lz(l,1,1) = iota(l)
        do 16 k = 1, ns
        do 15 l = 1, nzeta
        do 14 m = 0, mpol1
        do 13 n = 1, 12
        work(k,l,m,n) = 0.
 13     continue
 14     continue
 15     continue
 16     continue
c       
        do 21 n=1, ntheta2
        do 20 m=1, nzeta
        do 19 l=1, ns
        lt(l,m,n) = 1.
        kz = l+ns*(m-1)+ns*nzeta*(n-1)
        lz(l,m,n) = iota(kz)
 19     continue
 20     continue
 21     continue
        do 30 n = 0,nmax
        rmncc(1,n,1) = rmncc(2,n,1)
c3d     rmnss(1,n,1) = rmnss(2,n,1)
c3d     zmncs(1,n,1) = zmncs(2,n,1)
 30     zmnsc(1,n,1) = zmnsc(2,n,1)
      endif
c
c            10.2   COMPUTE R, Z, AND LAMBDA IN REAL SPACE
c               a.  INVERSE TRANSFORM IN N-ZETA
c
      do 70 m = mparity,mpol1,2
        do 50 n = 0,nmax
          do 50 k = 1,nzeta
CDIR$ IVDEP
            do 50 js= jmin(m),ns
        work(js,k,m,1) = work(js,k,m,1) + rmncc(js,n,m)*cosnv (k,n)
c3d     work(js,k,m,2) = work(js,k,m,2) + rmnss(js,n,m)*sinnv (k,n)
c3d     work(js,k,m,3) = work(js,k,m,3) + rmncc(js,n,m)*sinnvn(k,n)
c3d     work(js,k,m,4) = work(js,k,m,4) + rmnss(js,n,m)*cosnvn(k,n)
c3d     work(js,k,m,5) = work(js,k,m,5) + zmncs(js,n,m)*sinnv (k,n)
        work(js,k,m,6) = work(js,k,m,6) + zmnsc(js,n,m)*cosnv (k,n)
c3d     work(js,k,m,7) = work(js,k,m,7) + zmncs(js,n,m)*cosnvn(k,n)
c3d     work(js,k,m,8) = work(js,k,m,8) + zmnsc(js,n,m)*sinnvn(k,n)
c3d     work(js,k,m,9) = work(js,k,m,9) + lmncs(js,n,m)*sinnv (k,n)
c3d     work(js,k,m,10)= work(js,k,m,10)+ lmnsc(js,n,m)*cosnv (k,n)
c3d     work(js,k,m,11)= work(js,k,m,11)+ lmncs(js,n,m)*cosnvn(k,n)
c3d     work(js,k,m,12)= work(js,k,m,12)+ lmnsc(js,n,m)*sinnvn(k,n)
 50       continue
c
c               b.  INVERSE TRANSFORM IN M-THETA
c
        do 60 i = 1,ntheta2
          cosmux = xmpq(m,1)*cosmu(i,m)
          sinmux = xmpq(m,1)*sinmu(i,m)
            do 60 jk = 1, nzeta*ns
              r( jk,1,i) = r( jk,1,i) +
     >          work(jk,1,m,01)*cosmu (i,m)
c3d  >      +   work(jk,1,m,02)*sinmu (i,m)
              rt(jk,1,i) = rt(jk,1,i) +
c3d  >          work(jk,1,m,02)*cosmum(i,m) +
     >          work(jk,1,m,01)*sinmum(i,m)
c3d           rz(jk,1,i) = rz(jk,1,i) +
c3d  >      work(jk,1,m,03)*cosmu (i,m) + work(jk,1,m,04)*sinmu (i,m)
              z( jk,1,i) = z( jk,1,i) +
c3d  >          work(jk,1,m,05)*cosmu (i,m) +
     >          work(jk,1,m,06)*sinmu (i,m)
              zt(jk,1,i) = zt(jk,1,i) +
     >          work(jk,1,m,06)*cosmum(i,m)
c3d  >       +  work(jk,1,m,05)*sinmum(i,m)
c3d           zz(jk,1,i) = zz(jk,1,i) +
c3d  >      work(jk,1,m,07)*cosmu (i,m) + work(jk,1,m,08)*sinmu (i,m)
c3d           lt(jk,1,i) = lt(jk,1,i) +
c3d         work(jk,1,m,10)*cosmum(i,m) + work(jk,1,m,09)*sinmum(i,m)
c3d           lz(jk,1,i) = lz(jk,1,i) -
c3d        (work(jk,1,m,11)*cosmu (i,m) + work(jk,1,m,12)*sinmu (i,m))
              rcon(jk,1,i) = rcon(jk,1,i) +
     >          work(jk,1,m,01)*cosmux
c3d  >       +  work(jk,1,m,02)*sinmux
              zcon(jk,1,i) = zcon(jk,1,i) +
c3d  >          work(jk,1,m,05)*cosmux   +
     >          work(jk,1,m,06)*sinmux
  60    continue
  70  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2jacb /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2jacb(r,z,rt,zt,rodd,zodd,rtodd,ztodd,zt12,
     >  rt12,zs,rs,gsqrt,r12,tau,wint,shalf,ohs,nrzt,nznt,irst)
        real r(*),zs(*),rs(*),gsqrt(*),z(*),rodd(*),zodd(*),r12(*),
     >  rtodd(*),ztodd(*),rt12(*),zt12(*),zt(*),rt(*),shalf(*),
     >  tau(*),wint(*)
c****** (RS, ZS)=(R, Z) SUB S, (RT12, ZT12)=(R, Z) SUB THETA,
c****** AND GSQRT=SQRT(G) ARE DIFFERENCED ON HALF MESH
c
      do 10 l=2,nrzt
        lm = l-1
        rt12(l) = .50*(rt(l)+rt(lm) + shalf(l)*(rtodd(l)+rtodd(lm)))
        zt12(l) = .50*(zt(l)+zt(lm) + shalf(l)*(ztodd(l)+ztodd(lm)))
        rs(l)   = ohs*( r(l)- r(lm) + shalf(l)*( rodd(l)- rodd(lm)))
        zs(l)   = ohs*( z(l)- z(lm) + shalf(l)*( zodd(l)- zodd(lm)))
        r12(l)  = .50*( r(l)+ r(lm) + shalf(l)*( rodd(l)+ rodd(lm)))
        gsqrt(l)= r12(l) * ( rt12(l)*zs(l) - rs(l)*zt12(l) + 0.25*
     >  ( rtodd(l)*zodd(l) + rtodd(lm)*zodd(lm)
     >  - ztodd(l)*rodd(l) - ztodd(lm)*rodd(lm) + (rt(l)*zodd(l)
     >  + rt(lm)*zodd(lm) - zt(l)*rodd(l) - zt(lm)*rodd(lm))/shalf(l)) )
        tau(l) = wint(l)*gsqrt(l)
  10  continue
c
c             TEST FOR SIGN CHANGE IN JACOBIAN
c
        taumax = 0.
        taumin = 0.
      do 20 l=2,nrzt
        taumax = max(tau(l),taumax)
        taumin = min(tau(l),taumin)
  20  continue
c
        if ( taumax*taumin .lt. 0. ) irst = 2
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2pres /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2pres(gsqrt,bsq,wint,wp,dnorm,mass,
     >  vp,pres,gam,nznt,ns)
        real mass(*),gsqrt(ns,*),bsq(ns,*),vp(*),pres(*),wint(ns,*)
c
      do 10 js=2,ns
        vp(js)=sdot(nznt,gsqrt(js,1),ns,wint(js,1),ns)
  10  continue
c
      do 20 js=2,ns
        vp(js)=dnorm*abs(vp(js))
        pres(js) = mass(js)
  20  continue
c
      if ( gam .gt. 0.0 ) then
        do 22 js=2,ns
          pres(js)=mass(js)/vp(js)**gam
  22    continue
      endif
c
        wp=sdot(ns,vp,1,pres,1)/float(ns-1)
c
      do 30 js=2,ns
        do 30 lk=1,nznt
          bsq(js,lk)=pres(js)
  30  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2bcov /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2bcov(btheta,bzeta,gsqrt,bsq,r12,
     >  rs,zs,rt12,zt12,bsubs,gtt,gtz,gzz,br,bphi,bz,worka)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        include 'ceq2n3.m'
        common/e2precond/
     &                 ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),
     >        cr(nsd1),azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)
        real lmatrix((2*mnd+1)*(2*mnd+1)),worka(*),ar(nsd),gsqrt(*),
     >  az(nsd),bsq(*),btheta(*),bzeta(*),r12(*),rt12(*),zt12(*),
     >  rs(*),zs(*),bsubs(*),gtt(*),gtz(*),gzz(*),br(*),bphi(*),bz(*)
c
c       ON ENTRY, BSQ IS THE KINETIC PRESSURE
c
        do 5 l = 1,nrzt
        gtt(l) = 0.
        gtz(l) = 0.
        gzz(l) = 0.
 5      continue
        do 15 is = -1,0
          do 10 l=2,nrzt
          lm = l + is
          br(l)  = .5*sqrts(lm)*sqrts(lm)
          gtt(l) = gtt(l) + .5*(rt(lm)*rt(lm) + zt(lm)*zt(lm)) +
     >     br(l)*(rtodd(lm)*rtodd(lm) + ztodd(lm)*ztodd(lm))   +
     >     shalf(l)*(rt(lm)*rtodd(lm) + zt(lm)*ztodd(lm))
c3d       gtz(l) = gtz(l) + .5*(rt(lm)*rz(lm) + zt(lm)*zz(lm)) +
c3d  >     br(l)*(rtodd(lm)*rzodd(lm) + ztodd(lm)*zzodd(lm))   +
c3d  >  .5*shalf(l)*(rt(lm)*rzodd(lm) + rz(lm)*rtodd(lm)
c3d  >             + zt(lm)*zzodd(lm) + zz(lm)*ztodd(lm))
 10       gzz(l) = gzz(l) +
c3d  >    .5*(rz(lm)*rz(lm) + zz(lm)*zz(lm)) +
c3d  >     br(l)*(rzodd(lm)*rzodd(lm) + zzodd(lm)*zzodd(lm))   +
c3d  >     shalf(l)*(rz(lm)*rzodd(lm) + zz(lm)*zzodd(lm))      +
     >    .5*r(lm)**2 + br(l)*rodd(lm)**2 + shalf(l)*r(lm)*rodd(lm)
 15     continue
c
c       COMPUTE IOTA PROFILE IF ZERO NET TOROIDAL CURRENT IS PRESCRIBED
c
        if(ncurr.eq.1)then
        do 25 l=2,nrzt
        btheta(l)=(lz(l)*gtt(l)+lt(l)*gtz(l))/gsqrt(l) - torcur(l)
 25     bzeta(l)=-phip(l)*gtt(l)/gsqrt(l)
        call e2curr(iotas,phips,lz,bzeta,btheta,wint,ns,nznt)
         endif
c
c       COMPUTE LAMBDA FOR TWO DIMENSIONAL CASE
c
        do 30 l =2,nrzt                                                 ! c2d
 30     bzeta(l) = gsqrt(l)/gzz(l)                                      ! c2d
        do 35 js = 2,ns                                                 ! c2d
        psiav = sdot(nznt,bzeta(js),ns,wint(2),ns) * dnorm              ! c2d
        do 35 l = js,nrzt,ns                                            ! c2d
 35     lt(l) = phip(l)*bzeta(l)/psiav                                  ! c2d

 50     if(iequi.eq.1)  call e2bss(r12
     >  ,rs,zs,rt12,zt12,shalf,bsubs,gsqrt,br,bphi,bz,nrzt)

        wb=-wp
        do 55 l=2,nrzt
        lt(l)=lt(l)/gsqrt(l)
        lz(l)=lz(l)/gsqrt(l)
        bzeta(l)=lt(l)*gzz(l)+lz(l)*gtz(l)
        btheta(l)=lt(l)*gtz(l)+lz(l)*gtt(l)
        gzz(l)=lt(l)*lt(l)*gsqrt(l)
        bsq(l)=.5*(lz(l)*btheta(l)+lt(l)*bzeta(l))+bsq(l)
 55     wb = wb + hs*dnorm*wint(l)*abs(gsqrt(l))*bsq(l)
c
C       COMPUTE PRE-CONDITIONING MATRIX ELEMENTS AND FORCE NORMS
c
        if(mod(iter2-iter1,ns4).eq.0)then
        call e2prec(gzz,bsq,gsqrt,r12,zs,zt12,zt,ztodd,zodd,shalf,
     >  wint,pres,arm,ard,brm,brd,cr,ohs,ns,iter2)
        call e2prec(gzz,bsq,gsqrt,r12,rs,rt12,rt,rtodd,rodd,shalf,
     >  wint,pres,azm,azd,bzm,bzd,cr,ohs,ns,iter2)
        do 60 l=2,nrzt
 60     gtt(l) = gtt(l)*r12(l)**2
        volume = hs*ssum(ns-1,vp(2),1)
        fnorm = dnorm/(sdot(nrzt,gtt,1,wint,1)*(wb/volume)**2)
c
c       COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
c
        do 65 js=2,ns-1
        ar(js) = 0.
 65     az(js) = 0.
        do 70 js=2,ns-1
        do 70 lk = 1,nrzt,ns
        ar(js) = ar(js) + wint(js+lk-1)*rt0(js+lk-1)**2
 70     az(js) = az(js) + wint(js+lk-1)*zt0(js+lk-1)**2
        do 75 js=2,ns-1
 75     tcon(js) = min(abs(ard(js,1)/ar(js)),abs(azd(js,1)/az(js)))
        tcon(ns) = tcon(ns-1)
        endif
c
        do 80 l=2,nrzt
        gtt(l)=lz(l)*lz(l)*gsqrt(l)
        gtz(l)=lt(l)*lz(l)*gsqrt(l)
        lz(l)=bsq(l)*gsqrt(l)/r12(l)
 80     lt(l)=bsq(l)*r12(l)
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2prec /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2prec(gzz,bsq,gsqrt,r12,xs,xt12,xte,xto,xodd,
     >  shalf,wint,pres,axm,axd,bxm,bxd,cx,ohs,ns,iter2)
c
        include 'ceq2n1.m'
        real bsq(ns,*),gsqrt(ns,*),r12(ns,*),gzz(ns,*),xt12(ns,*),
     >  xte(ns,*),xto(ns,*),xodd(ns,*),xs(ns,*),wint(ns,*),shalf(*),
     >  pres(*),ax(nsd1,4),bx(nsd1,4),sm(nsd1),sp(0:nsd1),ptau(nznt),
     >  axm(nsd1,2),axd(nsd1,2),bxm(nsd1,2),bxd(nsd1,2),cx(*)
c
      if(iter2.eq.1)then
!ap {
        sm = 0.0
        sp = 0.0  ! }
        do 5 js=2,ns
          sm(js) = sqrt( (js-1.5)/(js-1.) )
          sp(js) = sqrt( (js-0.5)/(js-1.) )
  5     continue
        sp(1)  = sm(2)
      endif
c
c       COMPUTE PRECONDITIONING MATRIX ELEMENTS (ALL ARE MULTIPLIED BY 0.5)
c
      do 10 i = 1,4
        do 10 js = 1,ns+1
          ax(js,i) = 0.
          bx(js,i) = 0.
  10  continue
      do 15 js = 1,ns+1
 15     cx(js) = 0.

      do 20 js = 2,ns
c
c       COMPUTE DOMINANT (1/DELTA-S)**2 PRECONDITIONING MATRIX ELEMENTS
c
        do 30 lk = 1,nznt
          ptau(lk) = r12(js,lk)**2*(bsq(js,lk)-pres(js))
     >         * wint(js,lk)/gsqrt(js,lk)
          t1 = xt12(js,lk)*ohs
          t2 = .25*(xte(js  ,lk)/shalf(js) + xto(js  ,lk))/shalf(js)
          t3 = .25*(xte(js-1,lk)/shalf(js) + xto(js-1,lk))/shalf(js)
          ax(js,1) = ax(js,1) + ptau(lk)*t1*t1
          ax(js,2) = ax(js,2) + ptau(lk)*(-t1+t3)*(t1+t2)
          ax(js,3) = ax(js,3) + ptau(lk)*(t1+t2)**2
          ax(js,4) = ax(js,4) + ptau(lk)*(-t1+t3)**2
  30    continue
c
c       COMPUTE ORDER M**2 PRECONDITIONING MATRIX ELEMENTS
c
        do 35 lk=1,nznt
          t1 = .5*(xs(js,lk) + .5*xodd(js  ,lk)/shalf(js))
          t2 = .5*(xs(js,lk) + .5*xodd(js-1,lk)/shalf(js))
          bx(js,1) = bx(js,1) + ptau(lk)*t1*t2
          bx(js,2) = bx(js,2) + ptau(lk)*t1**2
          bx(js,3) = bx(js,3) + ptau(lk)*t2**2
          cx(js) = cx(js) + .25*gzz(js,lk)*wint(js,lk)
  35    continue
  20  continue
c
      do 40 js = 1,ns
        axm(js,1) =-ax(js,1)
        axd(js,1) = ax(js,1) + ax(js+1,1)
        axm(js,2) = ax(js,2) * sm(js) * sp(js-1)
        axd(js,2) = ax(js,3)*sm(js)**2 + ax(js+1,4)*sp(js)**2
        bxm(js,1) = bx(js,1)
        bxm(js,2) = bx(js,1) * sm(js) * sp(js-1)
        bxd(js,1) = bx(js,2) + bx(js+1,3)
        bxd(js,2) = bx(js,2)*sm(js)**2 + bx(js+1,3)*sp(js)**2
        cx(js) = cx(js) + cx(js+1)
  40  continue
        axd(ns,1) = 1.25*axd(ns,1)
        axd(ns,2) = 1.25*axd(ns,2)
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2lamc   .../baldur/code/bald/deqmom2.f
c rgb 25-apr-96 replace call sgefa with call matrx1 in BALDUR
c   replace sgesl with matrx2
c
        subroutine e2lamc(phip1,gtt,gtz,gzz,gsqrt,lt,lz,
     >  lmatrix,matrd,worka,workb,lmncs,lmnsc)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real lt(ns,nzeta,*),lz(ns,nzeta,*),gtt(ns,nzeta,*),
     >  gtz(ns,nzeta,*),gzz(ns,nzeta,*),ltemp(2*mnd+1),
     >  lmatrix(matrd,matrd),worka(nzeta,0:mpol2,6),
     >  gttc(0:nmax2,0:mpol2),gtts(0:nmax2,0:mpol2),
     >  gtzc(0:nmax2,0:mpol2),gtzs(0:nmax2,0:mpol2),
     >  gzzc(0:nmax2,0:mpol2),gzzs(0:nmax2,0:mpol2),
     >  lmncs(ns,0:nmax,0:mpol1),lmnsc(ns,0:nmax,0:mpol1)
        real gsqrt(ns,nzeta,*),tsign(0:nmax,0:nmax),ym(mnd),yn(mnd)
        real phip1(ns,nzeta,*),i1(0:nmax,0:mpol1),i2(0:nmax,0:mpol1),
     >  i3(0:nmax,0:mpol1),i4(0:nmax,0:mpol1),workb(ns,nzeta,0:mpol1,4)
        integer ipvt(2*mnd+1),iyn(mnd),iym(mnd)
        data isetup/1/, mn /0/

        if(isetup.eq.1)then
        isetup=2
        do 1 n = 0,nmax
        do 1 np= 0,nmax
        if(n.ge.np)tsign(n,np) = 1.
 1      if(n.lt.np)tsign(n,np) = -1.
        do 2 m = 0,mpol1
        do 2 n = 0,nmax
        mn = mn + 1
        iym(mn) = m
        iyn(mn) = n
         ym(mn) = mscale(m)*nscale(n)*float(m)
 2       yn(mn) = mscale(m)*nscale(n)*float(n*nfp)
        endif
c
c       COMPUTE MATRIX ELEMENTS FOR INDIVIDUAL RADIAL SURFACES
c       (FACTOR OF .25 IS SUPPRESSED - SOURCE TERMS ARE MULTIPLIED BY 4)
c
        do 130 js = 2,ns

        do 10 l = 1,6*nzeta*(mpol2+1)
 10     worka(l,0,1) = 0.
        do 20 l = 0,(nmax2+1)*(mpol2+1)-1
        gttc(l,0) = 0.
        gtts(l,0) = 0.
        gtzc(l,0) = 0.
        gtzs(l,0) = 0.
        gzzc(l,0) = 0.
 20     gzzs(l,0) = 0.
c
c       PERFORM M-THETA TRANSFORMS FIRST ON 0 <= m <= MPOL2
c
        do 30 k = 1,nzeta
        do 30 m = 0,mpol2
CDIR$ IVDEP
        do 30 i = 1,ntheta2
        og = 1./(gsqrt(js,k,i)*mscale(m))
        worka(k,m,1) = worka(k,m,1) + gtt(js,k,i)*cosmui(i,m)*og
        worka(k,m,2) = worka(k,m,2) + gtt(js,k,i)*sinmu (i,m)*og
        worka(k,m,3) = worka(k,m,3) + gtz(js,k,i)*cosmui(i,m)*og
        worka(k,m,4) = worka(k,m,4) + gtz(js,k,i)*sinmu (i,m)*og
        worka(k,m,5) = worka(k,m,5) + gzz(js,k,i)*cosmui(i,m)*og
 30     worka(k,m,6) = worka(k,m,6) + gzz(js,k,i)*sinmu (i,m)*og
c
c       DO N-ZETA TRANSFORM
c
        do 40 m = 0,mpol2
        do 40 n = 0,nmax2
        do 40 k = 1,nzeta
        gttc(n,m) = gttc(n,m) + worka(k,m,1)*cosnv(k,n)/nscale(n)
        gtts(n,m) = gtts(n,m) + worka(k,m,2)*sinnv(k,n)/nscale(n)
        gtzc(n,m) = gtzc(n,m) + worka(k,m,3)*cosnv(k,n)/nscale(n)
        gtzs(n,m) = gtzs(n,m) + worka(k,m,4)*sinnv(k,n)/nscale(n)
        gzzc(n,m) = gzzc(n,m) + worka(k,m,5)*cosnv(k,n)/nscale(n)
 40     gzzs(n,m) = gzzs(n,m) + worka(k,m,6)*sinnv(k,n)/nscale(n)
c
c       COMPUTE COEFFICIENTS OF PHIP AND CHIP FOR SOURCE TERMS (RHS)
c       (THE FACTOR OF 4 ARISES FROM WRITING COS COS ... AS SUMS ON LHS)
c
        mn = 0
        do 50 m = 0,mpol1
        do 50 n = 0,nmax
        mn = mn + 1
        i1(n,m) =-4.*(yn(mn)*gttc(n,m) + ym(mn)*gtzs(n,m))
        i2(n,m) = 4.*(yn(mn)*gtzc(n,m) + ym(mn)*gzzs(n,m))
        i3(n,m) = 4.*(yn(mn)*gtts(n,m) + ym(mn)*gtzc(n,m))
 50     i4(n,m) =-4.*(yn(mn)*gtzs(n,m) + ym(mn)*gzzc(n,m))

        do 60 j = 1,mnd
CDIR$ IVDEP
        do 62 i = j,mnd
        mp = iym(i) + iym(j)
        np = iyn(i) + iyn(j)
        mm = iym(i) - iym(j)
        nm = iabs(iyn(i) - iyn(j))
        sgn = tsign(iyn(i),iyn(j))
        t1 = gttc(np,mp) + gttc(np,mm) + gttc(nm,mp) + gttc(nm,mm)
        t2 = gtzs(np,mp) + gtzs(np,mm) +(gtzs(nm,mp) + gtzs(nm,mm))*sgn
        t3 = gtzs(np,mp) - gtzs(np,mm) +(gtzs(nm,mm) - gtzs(nm,mp))*sgn
        t4 = gzzc(np,mp) + gzzc(nm,mm) -(gzzc(nm,mp) + gzzc(np,mm))
        lmatrix(i,j) = yn(i)*(yn(j)*t1 + ym(j)*t3)
     >               + ym(i)*(yn(j)*t2 + ym(j)*t4)
        t1 = gttc(np,mp) + gttc(nm,mm) -(gttc(nm,mp) + gttc(np,mm))
        t4 = gzzc(np,mp) + gzzc(np,mm) + gzzc(nm,mp) + gzzc(nm,mm)
        lmatrix(i+mnd,j+mnd) = yn(i)*(yn(j)*t1 + ym(j)*t2)
     >                       + ym(i)*(yn(j)*t3 + ym(j)*t4)
        t1 = gtts(np,mp) - gtts(np,mm) +(gtts(nm,mm) - gtts(nm,mp))*sgn
        t2 = gtzc(np,mp) + gtzc(nm,mm) -(gtzc(nm,mp) + gtzc(np,mm))
        t3 = gtzc(np,mp) + gtzc(np,mm) + gtzc(nm,mp) + gtzc(nm,mm)
        t4 = gzzs(np,mp) + gzzs(np,mm) +(gzzs(nm,mm) + gzzs(nm,mp))*sgn
 62     lmatrix(i,j+mnd) =-(yn(i)*(yn(j)*t1 + ym(j)*t3)
     >                   +  ym(i)*(yn(j)*t2 + ym(j)*t4))
        do 66 i = 1,j-1
        mp = iym(i) + iym(j)
        np = iyn(i) + iyn(j)
        mm =-iym(i) + iym(j)
        nm = iabs(iyn(i) - iyn(j))
        sgn = tsign(iyn(i),iyn(j))
        t1 = gtts(np,mp) + gtts(np,mm) -(gtts(nm,mm) + gtts(nm,mp))*sgn
        t2 = gtzc(np,mp) + gtzc(nm,mm) -(gtzc(nm,mp) + gtzc(np,mm))
        t3 = gtzc(np,mp) + gtzc(np,mm) + gtzc(nm,mp) + gtzc(nm,mm)
        t4 = gzzs(np,mp) - gzzs(np,mm) +(gzzs(nm,mp) - gzzs(nm,mm))*sgn
 66     lmatrix(i,j+mnd) =-(yn(i)*(yn(j)*t1 + ym(j)*t3)
     >                   +  ym(i)*(yn(j)*t2 + ym(j)*t4))
 60     continue
c
c       COMPUTE SYMMETRIC BLOCKS
c
        do 70 i = 1,mnd
CDIR$ IVDEP
        do 70 j = 1,mnd
 70     lmatrix(j+mnd,i) = lmatrix(i,j+mnd)
        do 74 i = 2,mnd
CDIR$ IVDEP
        do 74 j = 1,i-1
        lmatrix(j,i) = lmatrix(i,j)
 74     lmatrix(j+mnd,i+mnd) = lmatrix(i+mnd,j+mnd)
c
c       COMPUTE EIGENVALUE MATRIX (TLAM)
c       FOR SCALING LAMBDA FORCES IN (M,N)-SPACE
c
        do 210 ntype=1,2
        mnoff = mnd*(ntype-1)
        do 210 mn = 1,mnd
        power = .6*float(iym(mn))/(1.+float(iym(mn)))
        tlam = 0.
        do 200 mnp = 1,mnd
 200    tlam = tlam + .25*
     >  (abs(lmatrix(mnp,mn+mnoff)) + abs(lmatrix(mnp+mnd,mn+mnoff)))
 210    if( tlam*phips(js).ne.0. )faclam(js +ns*(mn-1) +(ntype-1)*mns) =
     >  -.666*(isigng/tlam/phips(js))*shalf(js)**power
        if(iter2.ne.iter1)goto 130
c
c       COMPUTE SOURCE (RHS)
c
        diota =-iotas(js)*float(1-ncurr)
CDIR$ IVDEP
        do 80 i = 1,mnd
        ltemp(i    ) = diota*i1(i-1,0) + i2(i-1,0)
 80     ltemp(i+mnd) = diota*i3(i-1,0) + i4(i-1,0)
c
c       COMPUTE MATRIX ELEMENTS, SOURCE IF NCURR = 1 (CURRENT PRESCRIBED)
c
        if(ncurr.eq.1)then
CDIR$ IVDEP
        do 90 i = 1,mnd
        lmatrix(i,matrd) = i1(i-1,0)
 90     lmatrix(i+mnd,matrd) = i3(i-1,0)
        lmatrix(matrd,matrd) = 4.*gttc(0,0)
        fac =4.*isigng/(twopi*dnorm*phips(js))
        ltemp(matrd) = -4.*gtzc(0,0) + fac*jzeta(js)
CDIR$ IVDEP
        do 100 i = 1,matrd
 100    lmatrix(matrd,i) = lmatrix(i,matrd)
        endif
c
c       REMOVE SINGULAR MATRIX ELEMENTS (SIN[M=0,N=0] TERMS)
c
        const = lmatrix(1+nmax1+mnd,1+nmax1+mnd)
c       (N = 0 TERMS)
CDIR$ IVDEP
        do 105 l = 1,1+mpol1*nmax1,nmax1
 105    lmatrix(l,l) = const
c       (M = 0 TERMS)
CDIR$ IVDEP
        do 107 l = 1+mnd,nmax1+mnd
 107    lmatrix(l,l) = const

c       FACTOR MATRIX AND SOLVE FOR LAMBDA

        call matrx1 ( lmatrix, matrd, matrd, ipvt, info )
c
cbate        call sgefa(lmatrix,matrd,matrd,ipvt,info)
c
        if(info.eq.0)go to 120
        print 110, info
        write(6,110)info
 110    format(' info is nonzero: info = ',i4)
        stop
 120    continue
c
        call matrx2 ( ltemp, ltemp, lmatrix, ipvt, matrd, matrd
     &    , 1, 1, 1 )
c
cbate       call sgesl(lmatrix,matrd,matrd,ipvt,ltemp,0)
c
c       STORE SOLUTION FOR LAMBDA IN LMNCS,LMNSC ARRAYS

        call scopy(mnd,ltemp,1,lmncs(js,0,0),ns)
        call scopy(mnd,ltemp(1+mnd),1,lmnsc(js,0,0),ns)
        if(ncurr.eq.1)then
        iotas(js) = ltemp(matrd)
        endif
 130    continue
c
c       RECONSTRUCT LAMBDA IN REAL SPACE
c
        if(iter2.ne.iter1)return
        do 140 lk = 1,nznt
        do 140 js = 2,ns
        lt(js,lk,1) = phips(js)
 140    lz(js,lk,1) = iotas(js)*phips(js)
        do 150 l = 1,4*ns*nzeta*mpol
 150    workb(l,1,0,1) = 0.
        do 160 m = 0,mpol1
        do 170 n = 0,nmax
        do 170 k = 1,nzeta
        do 170 js= 2,ns
        workb(js,k,m,1) = workb(js,k,m,1) + lmncs(js,n,m)*sinnv (k,n)
        workb(js,k,m,2) = workb(js,k,m,2) + lmnsc(js,n,m)*cosnv (k,n)
        workb(js,k,m,3) = workb(js,k,m,3) + lmncs(js,n,m)*cosnvn(k,n)
 170    workb(js,k,m,4) = workb(js,k,m,4) + lmnsc(js,n,m)*sinnvn(k,n)
        do 180 i = 1,ntheta2
        do 180 jk = 1,nzeta*ns
        lt(jk,1,i) = lt(jk,1,i) + phip1(jk,1,i)*
     >     (workb(jk,1,m,2)*cosmum(i,m) + workb(jk,1,m,1)*sinmum(i,m))
 180    lz(jk,1,i) = lz(jk,1,i) - phip1(jk,1,i)*
     >     (workb(jk,1,m,3)*cosmu (i,m) + workb(jk,1,m,4)*sinmu (i,m))
 160    continue
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2curr /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2curr(iota,phip,lz,bze,bth,wt,ns,nznt)
c
        real lz(ns,*),iota(*),phip(*),bze(ns,*),bth(ns,*),wt(*)
c
      do 10 js=2,ns
        iota(js) = sdot(nznt,bth(js,1),ns,wt(2),ns)
     >          /sdot(nznt,bze(js,1),ns,wt(2),ns)
        do 10 lk=1,nznt
          lz(js,lk)=lz(js,lk) + iota(js)*phip(js)
  10  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2forc /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2forc(rbsq,gtt,gtz,gzz,sqrts,shalf,ohs,nrzt,ns,ivac)
        include 'ceq2n1.m'
        include 'ceq2n3.m'
        include 'ceq2n4.m'
        real rbsq(*),gtt(*),gtz(*),gzz(*),gtts(nrztd),
     >  gtzs(nrztd),gzzs(nrztd),shalf(*),sqrts(*)  
        real aredge(nznt),azedge(nznt),bredge(nznt),bzedge(nznt),
     >  credge(nznt),czedge(nznt),aroedge(nznt),azoedge(nznt)
        equivalence (z,gtzs),(crmn,gzzs),(czodd,gtts)
c
c       ON ENTRY, ARMN=ZT,BRMN=ZS,AZMN=RT,BZMN=RS,LT=R*BSQ
c       IT IS ESSENTIAL THAT LT,LZ AT j=1 ARE ZERO INITIALLY
c
      do 5 l = 1,nrzt,ns
        lt (l) = 0.
 5      lz (l) = 0.
c
      do 10 l=1,nrzt
        gtts(l)  = gtt(l)*shalf(l)
        gtzs(l)  = gtz(l)*shalf(l)
        gzzs(l)  = gzz(l)*shalf(l)
        armn(l)  = ohs*armn(l)*lt(l)
        azmn(l)  =-ohs*azmn(l)*lt(l)
        arodd(l) = armn(l)*shalf(l)
        azodd(l) = azmn(l)*shalf(l)
        brmn(l)  = brmn(l)*lt(l)
        bzmn(l)  =-bzmn(l)*lt(l)
        brodd(l) = brmn(l)*shalf(l)
        bzodd(l) = bzmn(l)*shalf(l)
 10     czmn(l)  = lt(l)/shalf(l)
CDIR$ IVDEP
!ap {
      czmn(nrzt+1) = 0.
      lt(nrzt+1)   = 0.
      gtt(nrzt+1)  = 0.
      gtts(nrzt+1) = 0.
      gtz(nrzt+1)  = 0.
      gtzs(nrzt+1) = 0.
      gzz(nrzt+1)  = 0.
      gzzs(nrzt+1) = 0. !}
      do 15 l = 1,nrzt
        czmn(l) = .25*(czmn(l) + czmn(l+1))
        lt(l)   = .25*(lt(l)   +   lt(l+1))
        gtt(l)  = .50*(gtt(l)  +  gtt(l+1))
        gtts(l) = .50*(gtts(l) + gtts(l+1))
        gtz(l)  = .50*(gtz(l)  +  gtz(l+1))
        gtzs(l) = .50*(gtzs(l) + gtzs(l+1))
        gzz(l)  = .50*(gzz(l)  +  gzz(l+1))
 15     gzzs(l) = .50*(gzzs(l) + gzzs(l+1))
c
c       CONSTRUCT AND SAVE EDGE FORCES
c
        if( ivac.ge.1 )then
        do 20 lk=1,nznt                                                 ! vac
        l1     = ns*lk                                                  ! vac
        ared = zt0(l1)*rbsq(lk)                                         ! vac
        azed =-rt0(l1)*rbsq(lk)                                         ! vac
        rzedge = rz(l1) + rzodd(l1)                                     ! vac
        zzedge = zz(l1) + zzodd(l1)                                     ! vac
        aredge(lk) = ared - armn(l1) + .5*lz(l1)                        ! vac
     >  - (gzz(l1)*r(l1) + gzzs(l1)*rodd(l1))                           ! vac
        azedge(lk) = azed - azmn(l1)                                    ! vac
        aroedge(lk)= ared - arodd(l1) + .5*lz(l1)*shalf(l1)             ! vac
     >  - (gzzs(l1)*r(l1) + gzz(l1)*rodd(l1) + czmn(l1)*zt0(l1))        ! vac
        azoedge(lk)= azed - azodd(l1) + czmn(l1)*rt0(l1)                ! vac
        bredge(lk) = .5*brmn(l1) + zodd(l1)*czmn(l1)                    ! vac
     >  - (gtt(l1)*rt0(l1) + gtz(l1)*rzedge)                            ! vac
        bzedge(lk) = .5*bzmn(l1) - rodd(l1)*czmn(l1)                    ! vac
     >  - (gtt(l1)*zt0(l1) + gtz(l1)*zzedge)                            ! vac
        credge(lk) = gtz(l1)*rt0(l1) + gzz(l1)*rzedge                   ! vac
 20     czedge(lk) = gtz(l1)*zt0(l1) + gzz(l1)*zzedge                   ! vac
        endif                                                           ! vac
c
c       CONSTRUCT CYLINDRICAL FORCE KERNELS
c
CDIR$ IVDEP
      do 30 l=1,nrzt-1
        s2 = sqrts(l)**2
        gtts2   = gtt(l)*s2
        gtzs2   = gtz(l)*s2
        gzzs2   = gzz(l)*s2
        armn(l) = armn(l+1) - armn(l) + .50*(lz(l) + lz(l+1))
     >  -(gzz(l)*r(l) + gzzs(l)*rodd(l))
        azmn(l) = azmn(l+1) - azmn(l)
        arodd(l)= arodd(l+1) - arodd(l)-(ztodd(l)*lt(l) + zt(l)*czmn(l))
     >          + .50*(lz(l)*shalf(l) + lz(l+1)*shalf(l+1))
     >          -     (gzzs(l)*r(l) + gzzs2*rodd(l))
        azodd(l)= azodd(l+1) - azodd(l) + rtodd(l)*lt(l) + rt(l)*czmn(l)
        brmn(l) = .50*(brmn(l) + brmn(l+1)) + zodd(l)*czmn(l)
     >  - (gtt(l)*rt(l)+gtts(l)*rtodd(l))
c3d  >  - (gtz(l)*rz(l)+gtzs(l)*rzodd(l))
        bzmn(l) = .50*(bzmn(l) + bzmn(l+1)) - rodd(l)*czmn(l)
     >  - (gtt(l)*zt(l)+gtts(l)*ztodd(l))
c3d  >  - (gtz(l)*zz(l)+gtzs(l)*zzodd(l))
        brodd(l)= .50*(brodd(l) + brodd(l+1)) + zodd(l)*lt(l)
     >  - (gtts(l)*rt(l)+gtts2*rtodd(l))
c3d  >  - (gtzs(l)*rz(l)+gtzs2*rzodd(l))
        bzodd(l)= .50*(bzodd(l) + bzodd(l+1)) - rodd(l)*lt(l)
     >  - (gtts(l)*zt(l)+gtts2*ztodd(l))
c3d  >  - (gtzs(l)*zz(l)+gtzs2*zzodd(l))
c3d        crodd(l)= gtzs(l) * rt(l) + gtzs2 * rtodd(l)
c3d     >          + gzzs(l) * rz(l) + gzzs2 * rzodd(l)
c3d        czodd(l)= gtzs(l) * zt(l) + gtzs2 * ztodd(l)
c3d     >          + gzzs(l) * zz(l) + gzzs2 * zzodd(l)
c3d        czmn(l) = gtz(l)  * zt(l) + gtzs(l) * ztodd(l)
c3d     >          + gzz(l)  * zz(l) + gzzs(l) * zzodd(l)
c3d        crmn(l) = gtz(l)  * rt(l) + gtzs(l) * rtodd(l)
c3d     >          + gzz(l)  * rz(l) + gzzs(l) * rzodd(l)
 30     continue
c
c       COMPUTE CONSTRAINT FORCE KERNELS
c
      do 40 l = 1,nrzt
        rcon(l) = (rcon(l) - rcon0(l)) * gcon(l)
        zcon(l) = (zcon(l) - zcon0(l)) * gcon(l)
        rt0(l) =  rt0(l) * gcon(l)
 40     zt0(l) =  zt0(l) * gcon(l)
c
c       ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
c
        if( ivac.lt.1 )return                                           ! vac
        do 50 lk = 1,nznt                                               ! vac
        l1 = ns*lk                                                      ! vac
        armn(l1) = aredge(lk)                                           ! vac
        arodd(l1)= aroedge(lk)                                          ! vac
        azmn(l1) = azedge(lk)                                           ! vac
        azodd(l1)= azoedge(lk)                                          ! vac
        brmn(l1) = bredge(lk)                                           ! vac
        brodd(l1)= bredge(lk)                                           ! vac
        bzmn(l1) = bzedge(lk)                                           ! vac
        bzodd(l1)= bzedge(lk)                                           ! vac
        crmn(l1) = credge(lk)                                           ! vac
        crodd(l1)= credge(lk)                                           ! vac
        czmn(l1) = czedge(lk)                                           ! vac
 50     czodd(l1)= czedge(lk)                                           ! vac
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2tomn /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2tomn(frcc,frss,fzcs,fzsc,flcs,flsc,armn,brmn,crmn,
     >  azmn,bzmn,czmn,blmn,clmn,arcon,azcon,brcon,bzcon,work,mparity)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        dimension frcc(ns,0:nmax,0:mpol1), frss(ns,0:nmax,0:mpol1),
     >  fzcs(ns,0:nmax,0:mpol1),fzsc(ns,0:nmax,0:mpol1),jmin(0:mpol ),
     >  flcs(ns,0:nmax,0:mpol1),flsc(ns,0:nmax,0:mpol1),
     >  armn(ns,nzeta,*),brmn(ns,nzeta,*),crmn(ns,nzeta,*),
     >  azmn(ns,nzeta,*),bzmn(ns,nzeta,*),czmn(ns,nzeta,*),
     >  blmn(ns,nzeta,*),clmn(ns,nzeta,*),arcon(ns,nzeta,*),
     >  azcon(ns,nzeta,*),brcon(*),bzcon(*),work(ns,nzeta,0:mpol1,12)
cbate        data jmin/1,2,mpol1*3/
c
      jmin(0) = 1
      jmin(1) = 2
      do j=2,mpol1+1
        jmin(j) = 3
      end do
c
      if ( mparity .eq. meven ) then
c        do 10 l = 1,12*ns*nzeta*mpol
c          work(l,1,0,1) = 0.
c  10    continue
        do 13 k = 1, ns
        do 12 l = 1, nzeta
        do 11 m = 0, mpol1
        do 10 n = 1, 12
        work(k,l,m,n) = 0.
 10     continue
 11     continue
 12     continue
 13     continue


      endif
c
c       DO M-THETA TRANSFORM FIRST
c
c      do 15 l = 1,nrzt
c        brmn(l,1,1) = brmn(l,1,1) + brcon(l)
c        bzmn(l,1,1) = bzmn(l,1,1) + bzcon(l)
c  15  continue
        do 16 n=1, ntheta2
        do 15 m=1, nzeta
        do 14 l=1, ns
        kz = l+ns*(m-1)+ns*nzeta*(n-1)
        brmn(l,m,n) = brmn(l,m,n) + brcon(kz)
        bzmn(l,m,n) = bzmn(l,m,n) + bzcon(kz)
 14     continue
 15     continue
 16     continue
c
      do 40 m = mparity,mpol1,2
        jmax = ns
        if(ivac.lt.1.or.mvac.eq.1)jmax = ns-1
        do 20 i = 1,ntheta2
          do 25 jk= 1,ns*nzeta
            temp1 = armn(jk,1,i) + xmpq(m,1) * arcon(jk,1,i)
            temp3 = azmn(jk,1,i) + xmpq(m,1) * azcon(jk,1,i)
        work(jk,1,m,01) = work(jk,1,m,01) + temp1        * cosmui(i,m)
     >                                    + brmn(jk,1,i) * sinmum(i,m)
c3d     work(jk,1,m,02) = work(jk,1,m,02) - crmn(jk,1,i) * cosmui(i,m)
        work(jk,1,m,03) = work(jk,1,m,03) + temp1        * sinmu (i,m)
     >                                    + brmn(jk,1,i) * cosmumi(i,m)
c3d     work(jk,1,m,04) = work(jk,1,m,04) - crmn(jk,1,i) * sinmu (i,m)
        work(jk,1,m,05) = work(jk,1,m,05) + temp3        * cosmui(i,m)
     >                                    + bzmn(jk,1,i) * sinmum(i,m)
c3d     work(jk,1,m,06) = work(jk,1,m,06) - czmn(jk,1,i) * cosmui(i,m)
        work(jk,1,m,07) = work(jk,1,m,07) + temp3        * sinmu (i,m)
     >                                    + bzmn(jk,1,i) * cosmumi(i,m)
c3d     work(jk,1,m,08) = work(jk,1,m,08) - czmn(jk,1,i) * sinmu (i,m)
 25       continue
c
c       NOTE FOR 2D: ELIMINATE 30 LOOP and ALL LINE WITH WORK9-WORK12
c        ONCE LAMBDA FORCE IS DEBUGGED
c
      do 30 jk = 1,ns*nzeta
        work(jk,1,m,09) = work(jk,1,m,09) + blmn(jk,1,i) * sinmum(i,m)
        work(jk,1,m,10) = work(jk,1,m,10) - clmn(jk,1,i) * cosmui(i,m)
        work(jk,1,m,11) = work(jk,1,m,11) + blmn(jk,1,i) * cosmumi(i,m)
        work(jk,1,m,12) = work(jk,1,m,12) - clmn(jk,1,i) * sinmu (i,m)
  30  continue
  20    continue
c
c       DO N-ZETA TRANSFORM
c
        do 45 n = 0,nmax
          do 45 k = 1,nzeta
            do 50 js= jmin(m),jmax
        frcc(js,n,m) = frcc(js,n,m) + work(js,k,m,01)*cosnv (k,n)
c3d  >                              + work(js,k,m,02)*sinnvn(k,n)
c3d     frss(js,n,m) = frss(js,n,m) + work(js,k,m,03)*sinnv (k,n)
c3d  >                              + work(js,k,m,04)*cosnvn(k,n)
c3d     fzcs(js,n,m) = fzcs(js,n,m) + work(js,k,m,05)*sinnv (k,n)
c3d  >                              + work(js,k,m,06)*cosnvn(k,n)
        fzsc(js,n,m) = fzsc(js,n,m) + work(js,k,m,07)*cosnv (k,n)
c3d  >                              + work(js,k,m,08)*sinnvn(k,n)
  50        continue
c
c       NOTE: CAN ELIMINATE 55 LOOP WHEN LAMBDA FORCE IS DEBUGGED
c
            do 55 js= 2+nlam,ns
        flcs(js,n,m) = flcs(js,n,m) + work(js,k,m,09)*sinnv (k,n)
     >                              + work(js,k,m,10)*cosnvn(k,n)
        flsc(js,n,m) = flsc(js,n,m) + work(js,k,m,11)*cosnv (k,n)
     >                              + work(js,k,m,12)*sinnvn(k,n)
  55        continue
  45    continue
  40  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2rsid /11040/bald91/wbaldn1  wsequil  DEQMOM2
c  rgb 11-aug-96 if ( nmax .lt. 1 ) then  n = nmax ...
        subroutine e2rsid(gcr,gcz,gcl,bz0,work)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        common/e2precond/
     &                 ard(nsd1,2),arm(nsd1,2),brd(nsd1,2),brm(nsd1,2),
     >        cr(nsd1),azd(nsd1,2),azm(nsd1,2),bzd(nsd1,2),bzm(nsd1,2)
        real gcr(ns,0:nmax,0:mpol1,2),gcz(ns,0:nmax,0:mpol1,2),
     >  gcl(ns,0:nmax,0:mpol1,2),work(mns,6)
c
c       IMPOSE M=1 MODE CONSTRAINT (ON Z CS, NTYPE = 1)
c
      if ( fsqz .lt. 1.e-8  .and.  (iter1 .ne. 1) ) then
        if ( nmax .lt. 1 ) then
          n = nmax
          do js= 2,ns
            gcz(js,n,1,1) = 0.
          enddo
        else
          do 10 n = 1,nmax
          do 10 js= 2,ns
 10       gcz(js,n,1,1) = 0.
        endif
      endif
c
c       CONSTRUCT INVARIANT RESIDUALS
c
        bz0  = 2.*hs/bz0**2
        fsql = bz0*sdot(2*mns,gcl,1,gcl,1)
        call e2getf(gcr,gcz,fsqr,fsqz,fnorm,meven)
        fedge = fnorm*(sdot(2*mnd,gcr(ns,0,0,1),ns,gcr(ns,0,0,1),ns)
     >        +        sdot(2*mnd,gcz(ns,0,0,1),ns,gcz(ns,0,0,1),ns))
c
c       PERFORM PRECONDITIONING AND COMPUTE RESIDUES
c
        call e2scal(gcr,arm,brm,ard,brd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
c
        call e2scal(gcz,azm,bzm,azd,bzd,cr,work,work(1,2),
     >  work(1,3),work(1,4),work(1,5))
c
        call e2getf(gcr,gcz,fsqr1,fsqz1,hs,modd)
c
c       CONSTRUCT PRECONDITIONED (SCALED) LAMBDA FORCES
c
      do 30 l=1,2*mns
 30   gc(l+4*mns)=faclam(l)*gc(l+4*mns)
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2scal /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2scal(gcx,axm,bxm,axd,bxd,cx,bx,dx,ax,gm,alf)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        dimension gcx(ns,0:nmax,0:mpol1,2),ax(ns,0:nmax,0:mpol1),gm(*),
     >  bx(ns,0:nmax,0:mpol1),dx(ns,0:nmax,0:mpol1),in(0:mpol),alf(*),
     >  axm(nsd1,2),axd(nsd1,2),cx(*),bxm(nsd1,2),bxd(nsd1,2)
cbate        data in/1,2,mpol1*3/
c
      in(0) = 1
      in(1) = 2
      do j=2,mpol1+1
        in(j) = 3
      end do
c
      jmax = ns
      if ( ivac .lt. 1  .or.  mvac .eq. 1 ) jmax = ns-1
      do 10 m = 0,mpol1
        mp = mod(m,2)+1
        do 10 n = 0,nmax
          do 20 js=in(m),jmax
            ax(js,n,m) = axm(js+1,mp) + bxm(js+1,mp)*m**2
            bx(js,n,m) = axm(js  ,mp) + bxm(js  ,mp)*m**2
            dx(js,n,m) = axd(js  ,mp) + bxd(js  ,mp)*m**2
     &                     + cx(js)*(n*nfp)**2
  20      continue
          if ( m .eq. 1 ) dx(2,n,1)  = dx(2,n,1) + bx(2,n,1)
          if ( n .ne. 0 ) dx(ns,n,m) = 8.*dx(ns,n,m)
          if ( n .ne. 0 ) bx(ns,n,m) = 8.*bx(ns,n,m)
  10  continue
c
      do 30 ntype=1,2
        call e2trid(ax,dx,bx,gcx(1,0,0,ntype),gm,alf,jmax,ns)
  30  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2trid /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2trid(a,d,b,c,gam,alf,nn,ns)
c
        include 'ceq2n1.m'
cnomo        parameter(mup2=nmax+nmax1,mlo3=mup2+1)
c
c       SOLVES B(I)*X(I-1)+D(I)*X(I)+A(I)*X(I+1)=C(I), I=IN,NN
c       AND RETURNS ANSWER IN C(I)
c
        dimension a(ns,0:mnd1),b(ns,0:mnd1),c(ns,0:mnd1),d(ns,0:mnd1),
     >  alf(0:mnd1,ns),gam(0:mnd1,ns),mupper(3),mlower(3),js(3)
cbate        data mupper/nmax,mup2,mnd1/, mlower/0,nmax1,mlo3/
cbate        data js /1,2,3/
c
      mup2 = nmax + nmax1
      mlo3 = mup2 + 1
      mupper(1) = nmax
      mupper(2) = mup2
      mupper(3) = mnd1
      mlower(1) = 0
      mlower(2) = nmax1
      mlower(3) = mlo3
      js(1) = 1
      js(2) = 2
      js(3) = 3
c
c       SEPARATE M=0 (IMODES=1), M=1 (IMODES=2), M>1 (IMODES=3)
c
      do 100 imodes = 1,3
        in = js(imodes)
        in1 = in + 1
        do 10 mn = mlower(imodes),mupper(imodes)
          gam(mn,in) = d(in,mn)
  10    continue
c
        do 20 i=in1,nn
          do 20 mn = mlower(imodes),mupper(imodes)
            alf(mn,i-1) = a(i-1,mn)/gam(mn,i-1)
            gam(mn,i)   = d(i,mn) - b(i,mn)*alf(mn,i-1)
  20    continue
c
        do 30 mn = mlower(imodes),mupper(imodes)
          c(in,mn) = c(in,mn)/gam(mn,in)
  30    continue
c
        do 40 i=in1,nn
          do 40 mn = mlower(imodes),mupper(imodes)
            c(i,mn) = (c(i,mn) - b(i,mn)*c(i-1,mn))/gam(mn,i)
  40    continue
c
        n2 = nn + in
        do 50 i=in1,nn
          i1 = n2 -i
          do 50 mn = mlower(imodes),mupper(imodes)
            c(i1,mn) = c(i1,mn) - alf(mn,i1)*c(i1+1,mn)
  50    continue
 100  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2getf /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2getf(gcr,gcz,gnormr,gnormz,gnorm,mprecon)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real gcr(ns,0:mnd1,2),gcz(ns,0:mnd1,2)
c
        gnormr = 0.
        gnormz = 0.
        jsmax = (ns-1) + mprecon
      do 10 mn = 0,mnd1
        do 10 js = 1,jsmax
          gnormr = gnormr + gnorm*(gcr(js,mn,1)**2 + gcr(js,mn,2)**2)
          gnormz = gnormz + gnorm*(gcz(js,mn,1)**2 + gcz(js,mn,2)**2)
  10  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2ntrp /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2ntrp(xnew,xold,scale,nsin)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real xnew(ns,0:nmax,0:mpol1,4),xold(nsin,0:nmax,0:mpol1,4)
        real scale(ns,0:nmax,0:mpol1)
c
        hsold=1./float(nsin-1)
c********   INTERPOLATE R,Z ON FULL GRID   ****************
      do 30 ntype = 1,4
        do 20 n = 0,nmax
 20     xold(1,n,1,ntype) = 2.*xold(2,n,1,ntype) - xold(3,n,1,ntype)
        do 10 js=1,ns
          sj=(js-1)*hs
          js1 = 1+int(sj/hsold+1.e-10)
          js2 = min0(js1+1,nsin)
          s1=(js1-1)*hsold
cornl          do 10 mn = 0,mnd1
          do 10 mn = 0,mpol1
            xint=(sj-s1)/hsold
            xnew(js,0,mn,ntype) = ((1.-xint)*xold(js1,0,mn,ntype)
     >                      + xint*xold(js2,0,mn,ntype))/scale(js,0,mn)
  10    continue
c
        do 40 n = 0,nmax
          xnew(1,n,1,ntype) = 0.
 40     continue
 30   continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2bss /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2bss(r12,rs,zs,rt12,zt12,shalf,bsubs,
     >  gsqrt,br,bphi,bz,nrzt)
        include 'ceq2n1.m'
        include 'ceq2n3.m'
        real r12(*),rs(*),zs(*),rt12(*),zt12(*),shalf(*),
     >  bsubs(*),br(*),bphi(*),bz(*),gsqrt(*)
        do 10  l=2,nrzt
        rz12 =      .5*(   rz(l)+   rz(l-1)
     >       +shalf(l)*(rzodd(l)+rzodd(l-1)))
        zz12 =      .5*(   zz(l)+   zz(l-1)
     >       +shalf(l)*(zzodd(l)+zzodd(l-1)))
        gst= rs(l)*rt12(l) + zs(l)*zt12(l)
     >        +.25*(rodd(l)*rtodd(l)+rodd(l-1)*rtodd(l-1)
     >             +zodd(l)*ztodd(l)+zodd(l-1)*ztodd(l-1)
     >        +(    rodd(l)*rt(l)   +rodd(l-1)*rt(l-1)
     >        +     zodd(l)*zt(l)   +zodd(l-1)*zt(l-1))/shalf(l))
        gsz= rs(l)*rz12    + zs(l)*zz12
     >        +.25*(rodd(l)*rzodd(l)+rodd(l-1)*rzodd(l-1)
     >             +zodd(l)*zzodd(l)+zodd(l-1)*zzodd(l-1)
     >        +(    rodd(l)*rz(l)   +rodd(l-1)*rz(l-1)
     >        +     zodd(l)*zz(l)   +zodd(l-1)*zz(l-1))/shalf(l))
        br(l)   = (lz(l)*rt12(l) + lt(l)*rz12 )/gsqrt(l)
        bphi(l) =                  lt(l)*r12(l)/gsqrt(l)
        bz(l)   = (lz(l)*zt12(l) + lt(l)*zz12 )/gsqrt(l)
 10     bsubs(l)= (lz(l)*gst     + lt(l)*gsz )/gsqrt(l)
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2rout /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2rout(bsq,gsqrt,bsubt,bsubz,bsubs,br,bphi,bz)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        include 'ceq2n3.m'
        real bmod(nznt),bsq(ns,*),gsqrt(ns,*),bsubt(ns,*),bsubz(ns,*),
     >  bsubs(ns,*),br(*),bphi(*),bz(*),phi(nsd),rmnc(mnmax),
     >  zmns(mnmax),lmns(mnmax),xm(mnmax),xn(mnmax),bmodmn(mnmax)
c       THIS SUBROUTINE CREATES THE FILE WOUT WHICH CONTAINS
c       THE CYLINDRICAL COORDINATE SPECTRA RMNC,ZMNS,LMNS
cbate        write(8)voli,gam,1./nfp,mnmax,nmax1,1,itfsq,niter/100+1
        do 20 js=1,ns
        call e2conv(rmnc,zmns,lmns,xm,xn,js,xc,xc(1+mns),
     >  xc(1+2*mns),xc(1+3*mns),xc(1+4*mns),xc(1+5*mns))
        mn = 0
        do 25 lk=1,nznt
 25     bmod(lk)=sqrt(2.*abs(bsq(js,lk)/bscale**2-pres(js)))
        do 20 m = 0,mpol1
        nmin0 = -nmax
        if(m.eq.0)nmin0 = 0
        do 20 n = nmin0,nmax
        mn = mn+1
        dmult=2.*dnorm/(mscale(m)*nscale(n))
        if(m.eq.0.and.n.eq.0)dmult=.5*dmult
        if(js.eq.1)go to 20
        gmn = 0.
        bmn = 0.
        bsubtmn = 0.
        bsubzmn = 0.
        bsubsmn = 0.
        do 200 j = 1,ntheta2
        do 200 k = 1,nzeta
        lk = k + nzeta*(j-1)
        if(n.ge.0)then
        tcosi = dmult*(cosmui(j,m)*cosnv(k,n)+sinmu(j,m) *sinnv(k,n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,n)-cosmui(j,m)*sinnv(k,n))
        else
        tcosi = dmult*(cosmui(j,m)*cosnv(k,-n)-sinmu(j,m) *sinnv(k,-n))
        tsini = dmult*(sinmu (j,m)*cosnv(k,-n)+cosmui(j,m)*sinnv(k,-n))
        endif
        bmn = bmn + tcosi*bmod(lk)
        gmn = gmn + tcosi*gsqrt(js,lk)
        bsubtmn = bsubtmn + tcosi*bsubt(js,lk)
        bsubzmn = bsubzmn + tcosi*bsubz(js,lk)
        bsubsmn = bsubsmn + tsini*bsubs(js,lk)
 200    continue
        if(js.eq.ns/2)bmodmn(mn) = bmn
 20     continue
cbate        write(8)xm(mn),xn(mn),rmnc(mn),zmns(mn),lmns(mn)/bscale,
cbate     >  bmn,gmn,bsubtmn/bscale,bsubzmn/bscale,bsubsmn/bscale
        phi(1)=0.
c       NOTE:  CHIPS, PHIPS, PRES, & CURRENTS WERE SCALED BY 1/BSCALE IN EQFOR
        do 30 js=2,ns
 30     phi(js)=twopi*hs*ssum(js-1,phips(2),1)
        fac=abs(bscale)**(gam-2.)
cbate        write(8)(iotas(js)*phips(js),mass(js)*fac,pres(js),phips(js),
cbate     >  bpco(js),bzco(js),phi(js),vp(js),jtheta(js),jzeta(js),
cbate     >  specw(js),js=2,ns)
cbate        write(8)(fsqt(i),wdot(i),i=1,100)
        if(nvac.eq.0)return                                             ! vac
cbate        open(unit=9,file='bextn',status='unknown')                      ! vac
        write(6,60)cpol/bscale,ctor/bscale,bscale                       ! vac
 60     format(/' IPOL = ',1pe16.8,' ITOR = ',1pe16.8,' BSCALE = ',     ! vac
     >  1pe16.8)                                                        ! vac
        do 80 iprint=1,2                                                ! vac
        if(iprint.eq.1)write(6,70)                                      ! vac
        if(iprint.eq.2)write(6,75)                                      ! vac
        ntskip=1+ntheta1/12                                             ! vac
        nzskip=1+nzeta/6                                                ! vac
        do 80 l=1,nzeta,nzskip                                          ! vac
        zeta =360*(l-1)/float(nzeta)                                    ! vac
        do 80 k=1,ntheta2,ntskip                                        ! vac
        lk=l+nzeta*(k-1)                                                ! vac
        if(iprint.eq.1)write(6,90)zeta,r(ns*lk)+rodd(ns*lk),            ! vac
     >  z(ns*lk)+zodd(ns*lk),(bsqsav(lk,n)/bscale**2,n=1,3),            ! vac
     >  bsqvac(lk)/bscale**2                                            ! vac
        if(iprint.eq.2)write(6,95)zeta,r(ns*lk)+rodd(ns*lk),            ! vac
     >  z(ns*lk)+zodd(ns*lk),br(ns*lk)/bscale,                          ! vac
     >  bphi(ns*lk)/bscale,bz(ns*lk)/bscale,brv(lk)/bscale,             ! vac
     >  bphiv(lk)/bscale,bzv(lk)/bscale                                 ! vac
 80     continue                                                        ! vac
 70     format(/4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,                       ! vac
     >  'BSQMHDI',5x,'BSQVACI',5x,'BSQMHDF',5x,'BSQVACF'/)              ! vac
 75     format(/4x,'ZETA',8x,' Rb ',8x,' Zb ',6x,                       ! vac
     >  'BR',8x,'BPHI',6x,'BZ',8x,'BRv',7x,'BPHIv',5x,'BZv',/)          ! vac
 90     format(1pe10.2,1p6e12.4)                                        ! vac
 95     format(1pe10.2,1p2e12.4,1p6e10.2)                               ! vac
        write(6,100)                                                    ! vac
 100    format(/,'    mb   nb     rbc         zbs     |B|(s=.5)',/)     ! vac
        do 110 mn=1,mnmax                                               ! vac
 110    write(6,115)nint(xm(mn)),nint(xn(mn)/nfp),rmnc(mn),zmns(mn),    ! vac
     >  bmodmn(mn)                                                      ! vac
 115    format(i5,i4,1p3e12.4)                                          ! vac
cbate        write(9,125)nfp,mpmax,mnmax                                     ! vac
cbate        write(9,120)(rmnc(mn),zmns(mn),xm(mn),xn(mn),mn=1,mnmax)        ! vac
cbate        write(9,120)(potvac(mn)/bscale,xmpot(mn),xnpot(mn),mn=1,mpmax)  ! vac
 120    format(1p6e12.4)                                                ! vac
 125    format(3i6)                                                     ! vac
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2conv /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2conv(rmnc,zmns,lmns,xm,xn,js,rmncc,rmnss,
     >  zmncs,zmnsc,lmncs,lmnsc)
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real rmnc(*),zmns(*),lmns(*),xm(*),xn(*),
     >  rmncc(ns,0:nmax,0:mpol1),rmnss(ns,0:nmax,0:mpol1),
     >  zmncs(ns,0:nmax,0:mpol1),zmnsc(ns,0:nmax,0:mpol1),
     >  lmncs(ns,0:nmax,0:mpol1),lmnsc(ns,0:nmax,0:mpol1)
c
c       CONVERTS INTERNAL MODE REPRESENTATION TO STANDARD
c       FORM FOR OUTPUT (COEFFICIENTS OF cos(mu-nv), sin(mu-nv))
c
        mn = 0
        do 10 m = 0,mpol1
        nmin0 = -nmax
        if(m.eq.0)nmin0 = 0
        do 10 n = nmin0,nmax
        n1 = iabs(n)
        t1 = mscale(m)*nscale(n1)
        mn = mn + 1
        xm(mn) = float(m)
        xn(mn) = float(n*nfp)
        if( n.eq.0 )then
                rmnc(mn) = t1*rmncc(js,n,m)
                zmns(mn) = t1*zmnsc(js,n,m)
                lmns(mn) = t1*lmnsc(js,n,m)
        else if( m.eq.0 )then
                rmnc(mn) = t1*rmncc(js,n,m)
                zmns(mn) =-t1*zmncs(js,n,m)
                lmns(mn) =-t1*lmncs(js,n,m)
        else if( js.gt.1 )then
                sign = float(n/n1)
                rmnc(mn) = .5*t1*(rmncc(js,n1,m) + sign*rmnss(js,n1,m))
                zmns(mn) = .5*t1*(zmnsc(js,n1,m) - sign*zmncs(js,n1,m))
                lmns(mn) = .5*t1*(lmnsc(js,n1,m) - sign*lmncs(js,n1,m))
        endif
 10     continue
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2qfor /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2qfor(bp,bz,bpco,bzco,jp,jz,equif,
     >  pres,phips,vp,iotas,bsq,specw,ohs,voli,twopi,bscale,
     >  raxis,zaxis,wint,dnorm,mscal,nscal,ns,isigng)
        include 'ceq2n1.m'
        real bp(ns,*),bz(ns,*),bsq(ns,*),bpco(*),bzco(*),wint(*),
     >  equif(*),pres(*),phips(*),vp(*),iotas(*),specw(*),jp(*),jz(*),
     >  raxis(ns,0:nmax),zaxis(ns,0:nmax),mscal(0:mpol2),nscal(0:nmax2)
        data phi1/0./,chi1/0./
        write(6,5)
 5      format(/,'     S       EQUIF      PHI       CHI      JTHETA ',
     >  '    IOTA      JZETA     P''/V''    V''/PHIP   <M>',/)
        betaxis =
     >   1.5*pres(2)/(dnorm*sdot(nznt,bsq(2,1),ns,wint(2),ns)-pres(2))
     >  -0.5*pres(3)/(dnorm*sdot(nznt,bsq(3,1),ns,wint(2),ns)-pres(3))
        do 15 i=1,ns
        bpco(i)=dnorm*sdot(nznt,bp(i,1),ns,wint(2),ns)/bscale
        bzco(i)=dnorm*sdot(nznt,bz(i,1),ns,wint(2),ns)/bscale
        pres(i)=pres(i)/bscale**2
 15     phips(i)=phips(i)/bscale
        do 20 js=2,ns-1
        t0=.5*(vp(js+1)+vp(js))
        jz(js)=  isigng*ohs*(bpco(js+1)-bpco(js))/t0
        jp(js)= -isigng*ohs*(bzco(js+1)-bzco(js))/t0
        aiotaf=.5*(iotas(js+1)+iotas(js))
        t1=jz(js)*aiotaf
        t2=ohs*(pres(js+1)-pres(js))/t0
        t3=.5*(vp(js+1)/phips(js+1) + vp(js)/phips(js))
        phi1=phi1+phips(js)/ohs
        chi1=chi1+iotas(js)*phips(js)/ohs
        equif(js)=(jp(js)-t1-t2*t3)/(abs(t1)+abs(jp(js))+abs(t2*t3))
        es=(js-1)/float(ns-1)
 20     write(6,30)es,equif(js),isigng*twopi*phi1,isigng*twopi*chi1,
     >  jp(js),aiotaf,jz(js),t2,t3,specw(js)
 30     format(1p9e10.2,0pf7.3)
        volf=(twopi**2)*ssum(ns,vp,1)/float(ns-1)
        write(6,40)voli,volf,betaxis
 40     format(/,'  INITIAL VOLUME = ',1pe20.10,'  FINAL VOLUME = ',
     >  1pe20.10/'  BETA ON AXIS   = ',2x,1pe12.4/)
        write(6,50)
 50     format(2x,'MAGNETIC AXIS COEFFICIENTS'/,
     >  '    n     raxis       zaxis'/)
        do 60 n=0,nmax
        t1 = mscal(0)*nscal(n)
 60     write(6,70)n,t1*raxis(1,n),-t1*zaxis(1,n)
 70     format(i5,1p2e12.4)
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2spct /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2spct(rcc,rss,zcs,zsc)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
        real rcc(ns,0:nmax,0:mpol1),rss(ns,0:nmax,0:mpol1),
     >  zsc(ns,0:nmax,0:mpol1),zcs(ns,0:nmax,0:mpol1)
        real t1(nsd),dnumer(nsd),denom(nsd)
c
      do 5 js=2,ns
        dnumer(js) = 0.
        denom (js) = 0.
   5  continue
c
      do 10 n = 0,nmax
        do 10 m = 1,mpol1
          scale = (mscale(m)*nscale(n))**2
          do 15 js=2,ns
            t1(js) = (rcc(js,n,m)**2 + rss(js,n,m)**2
     >              + zcs(js,n,m)**2 + zsc(js,n,m)**2)*scale
  15      continue
          do 25 js=2,ns
            dnumer(js) = dnumer(js) + t1(js)*xmpq(m,3)
            denom (js) = denom (js) + t1(js)*xmpq(m,2)
  25      continue
  10  continue
c
      do 30 js=2,ns
        specw(js) = dnumer(js)/denom(js)
  30  continue
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@e2prnt /11040/bald91/wbaldn1  wsequil  DEQMOM2
c
        subroutine e2prnt(i,w0,r00)
c
        include 'ceq2n1.m'
        include 'ceq2n2.m'
c
        betav=wp/wb
        w=w0*twopi*twopi/nfp
        avm=0.
        den=0.
        specw(1)=1.
        call e2spct(xc,xc(1+mns),xc(1+2*mns),xc(1+3*mns))
      do 10 j=2,ns
        den=den+phips(j)
        avm=avm+phips(j)*(specw(j)+specw(j-1))
  10  continue
        avm = 0.5*avm/den
c
        if( ivac.ge.1 )delbsq=sdot(nznt,dbsq,1,wint(1+nznt),1)/         ! vac
     >                 sdot(nznt,bsqsav(1,3),1,wint(1+nznt),1)          ! vac
c
cbate        if(i.eq.1.and.nvac.ne.0)print 20
        if(i.eq.1.and.nvac.ne.0)write(6,15)
cbate        if(i.eq.1.and.nvac.eq.0)print 30
        if(i.eq.1.and.nvac.eq.0)write(6,25)
 15     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     ',
     >  '<M>   DEL-BSQ   FEDGE'/)
 20     format(/,' ITER    FSQR      FSQZ      FSQL      ',
     >  'R00(0)     WMHD      DEL-BSQ'/)
 25     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz      DELT     R00(0)        WMHD        BETA     <M>'/)
 30     format(/,' ITER    FSQR      FSQZ      FSQL      fsqr      ',
     >  'fsqz     R00(0)     WMHD'/)
c
        if(nvac.ne.0)go to 35
cbate        print 45, i,fsqr,fsqz,fsql,fsqr1,fsqz1,r00,w
        write(6,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,betav,avm
        return
c
 35     continue
cbate        print 50, i,fsqr,fsqz,fsql,r00,w,delbsq
        write(6,40)i,fsqr,fsqz,fsql,fsqr1,fsqz1,delt,r00,w,
     >  betav,avm,delbsq,fedge
 40     format(i5,1p6e10.2,1pe11.4,1pe15.8,1pe9.2,0pf7.3,1p2e9.2)
 45     format(i5,1p5e10.2,1pe10.3,1p2e11.4)
 50     format(i5,1p3e10.2,1pe10.3,1p2e11.4)
c
        return
        end

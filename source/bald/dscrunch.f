c 12:00 25-aug-96  .../baldur/code/bald/dscrunch.f   Bateman, Lehigh
c SCRUNCH routines from S.P. Hirshman, April 1990
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c**********************************************************************c
c
c    To obtain this file type:
c cfs get /11040/bald90/wbaldn1 ^ end     (note ^ stands for linefeed)
c lib wbaldn1 ^ x wsequil ^ n wsequil ^ x dscrunch ^ end
c lib wbaldn2 ^ x yequilib yutilib ^ end  ( libraries of compiled sbrtns)
c lib wbaldn3 ^ x xstripx ^ end           ( utility to strip !...)
c
c    To compile this file type:
c cosmos dscrunch
c
c**********************************************************************c
c@changes  baldur/code/bald/dscrunch.f
c rgb 25-aug-96 force more than one iteration with if (iter .gt. 1 )...
c rgb 02-jun-92 nresets initialized in sbrtn scrunch and passed to
c   sbrtn scrstrt where it is incremented and tested
c rgb 06-jul-90 changed ftol from 1.e-5 to 1.e-14 in data statement
c rgb 05-jul-90 changed parameter nu from 32 to 128
c rgb 30-apr-90 changed parameters nu=ntheta from 30 to 32
c     initialize r0n(i)=0.0 and z0n(i)=0.0 in sbrtn scgues
c    changed write(3,.. to write(6,..., where unit 6 is file scrout
c**********************************************************************c
c
c       THIS IS PROGRAM DESCUR - WHICH USES A STEEPEST DESCENT
c       ALGORITHM TO FIND A LEAST SQUARES APPROXIMATION TO AN 
c       ARBITRARY 3-D SPACE CURVE. ANGLE CONSTRAINTS, BASED ON A
c       MINIMUM SPECTRAL WIDTH CRITERION, ARE APPLIED.
c       THE CONSTRAINT IS SATISFIED BY TANGENTIAL VARIATIONS
c       ALONG THE CURVE.
c
c       THE MAIN PROGRAM SETS UP THE INITIAL POINTS TO BE FITS AND
c       THEN CALLS THE SUBROUTINE SCRUNCH, WHICH DOES THE ACTUAL
c       CURVE-FITTING.
c***********************************************************************
c       REVISION 1:  January 26, 1989
c                    Compute Fit to individual toroidal planes separately
c                    and then Fourier decompose in phi, using optimized 
c                    angle representation
c From HIRSHMAN@ORFE.ESNET on April 16, 1990 at 7:32 EDT
c Subject & Phone: Let me know if this works all right...
c***********************************************************************
c
c     Subroutines in driver program:
c  printit
c  plotter
c  bubble
c  setpag0
c  flsur
c  totz
c
c     Subroutines in SCRUNCH package:
c
c  scrunch
c  scfixa
c  scgues
c  scorder
c  scgeta
c  scampl
c  scrstrt
c  scevol
c  scfnct
c  scfftr
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scrunch
        subroutine scrunch(rin,zin,pexp,qexp,rbc,zbs,rmnaxis,zmnaxis,
     &                     ftol,nthetain,nphin,mpolin,nfpin,niter,nstep)
        external case17
        include 'name1.m'
        real rin(nthetain,nphin),zin(nthetain,nphin),rbc(*),zbs(*),
     &  rmnaxis(*),zmnaxis(*)
        real :: mclock

c**********************************************************************
c       DATA USED INTERNALLY (USUALLY NOT TO BE CHANGED BY USER)
c**********************************************************************
        data irst /1/, tottime/0./
        
c**********************************************************************
c       CHECK THAT INPUT DIMENSIONS ARE CONSISTENT WITH NAME1 PARAMETERS
c**********************************************************************
        if( (mod(nphin,2).ne.0) .and. nphin.ne.1 .or. nphin.eq.0 )then
        write (6,*) ' NPHI & 0 must be EVEN for non-symmetric systems'
        stop
        endif
        if( nthetain.gt.nu )then
        write (6,*) ' NTHETA input to SCRUNCH is too large'
        stop
        endif
        if( nphin.gt.nv )then
        write (6,*) ' NPHI input to SCRUNCH is too large'
        stop
        endif
        if( mpolin.gt.mu )then
        write (6,*) ' MPOL input to SCRUNCH is too large'
        stop
        endif
        if( nfpin.le.0 )then
        write (6,*) ' NFP & 0 must be positive and exceed zero'
        stop
        endif

c**********************************************************************
c       INITIALIZE FIXED m,n ARRAYS AND FILL INPUT COMMON BLOCK
c**********************************************************************
cbate        open(unit=6,file='scrout',status='unknown')
        twopi=8.*atan(1.)
        call scfixa(pexp,qexp,nthetain,mpolin,nphin,nfpin)

c**********************************************************************
c       COMPUTE INITIAL GUESSES (MUST ORDER ANGLES MONOTONICALLY FOR
c       NON-STARLIKE DOMAINS)
c**********************************************************************
        call scgues(rin,zin)

c**********************************************************************
c       BEGIN MAIN INTERATION LOOP
c**********************************************************************
cap
csgi    call second(timein)
        timein = 0.01 * mclock()
c
        write (6,10)
 10     format(/' ITERATIONS    RMS ERROR    FORCE GRADIENT    <M>',
     &  '    m #''s <= ON')
c
      nresets = 0
c 
        do 1000 nplane = 1,nphi2
c
        write (6,200)nplane
 200    format(/'                  Fitting toroidal plane # ',i3)
c**********************************************************************
c       INITIALIZE M=0 and M=1 MODE AMPLITUDES
c
c       STACKING OF XVEC ARRAY (FOR FIXED TOROIDAL PLANE)
c       XVEC(1,mpol):Rmcos            XVEC(1+mpol,2*mpol):Rmsin  
c       XVEC(1+2*mpol,3*mpol):Zmcos   XVEC(1+3*mpol,4*mpol):Zmsin
c       XVEC(4*mpol,n2): Theta angle
c**********************************************************************
        do 20 n = 1,n2
        xstore(n)= 0.
        xdot(n)  = 0.
 20     xvec(n)  = 0.

        call scampl(r0n(nplane),z0n(nplane),angle(1,nplane),
     &  xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3),
     &  xvec(1+mpol4),rin(1,nplane),zin(1,nplane))
        r10 = .5*(abs(xvec(2      )) + abs(xvec(2+ mpol))
     &      +     abs(xvec(2+mpol2)) + abs(xvec(2+mpol3)))
        imodes = mpol
        delt=deltf

        do 30 iter=1,niter
        call scfnct(xvec,xvec(1+mpol),xvec(1+mpol2),xvec(1+mpol3),
     &  xvec(1+mpol4),gvec,gvec(1+mpol),gvec(1+mpol2),
     &  gvec(1+mpol3),gvec(1+mpol4),fsq,r10,rin(1,nplane),
     &  zin(1,nplane),imodes)

        go to (50,40)iter
        gmin = min(gmin,gnorm)
        call scevol(g11)
c**********************************************************************
c       RUDIMENTARY TIME STEP CONTROL
c**********************************************************************
        if(gnorm/gmin.gt.1.e6)irst = 2
        if ( (irst.eq.2) .or. (gmin.eq.gnorm) )
     &    call scrstrt(irst,nresets)
 
        if(gnorm.lt.1.e-4.and.imodes.lt.mpol)then
        imodes = imodes + 1
        call scrstrt(irst,nresets)
        irst = 2
        delt = delt/.95
        endif
        if(mod(iter,nstep).eq.0.or.(gnorm.lt.ftol**2))go to 60
        go to 30
 40     g11=gnorm
        go to 30
 50     gmin = gnorm
        imodes = 3
 60     gout = sqrt(gnorm)
        modeno = imodes - 1
        if(iter.eq.1)modeno = mpol-1
cap
csgi    call second(timeou)
        timeou = 0.01 * mclock()
        tottime = tottime + (timeou - timein)
c
        write (6,110) iter,fsq,gout,specw,modeno
cap
csgi    call second(timein)
        timein = 0.01 * mclock()
        if ( (iter .gt. 1 ) .and.
     &   ( (gnorm.lt.ftol**2.and.imodes.eq.mpol)
     &     .or.fsq.lt.ftol) ) go to 80
 30     continue
 110    format(i8,1p2e16.3,0pf10.2,i8)
c**********************************************************************
c       STORE FINAL ANSWER FOR FINAL PHI-TRANSFORM
c**********************************************************************
 80     do 120 ntype = 1,4
        ntoff = mpol*(ntype-1)
        do 120 m = 1,mpol
 120    result(nplane,m,ntype) = xvec(m+ntoff)
 1000   continue
c**********************************************************************
c       OUTPUT LOOP
c**********************************************************************

        write (6,90) tottime
 90     format(/,' COMPUTATIONAL TIME = ',1pe12.3,' SECONDS'/)
 320    continue
        write (6,330) pexp,qexp
 330    format(' ANGLE CONSTRAINTS WERE APPLIED ',/,
     &  ' BASED ON RM**2 + ZM**2 SPECTRUM WITH P = ',
     &  f8.2,' AND Q = ',f8.2,/)
c**********************************************************************
c       PERFORM PHI-FOURIER TRANSFORM
c**********************************************************************
        call scfftr(result(1,1,1),result(1,1,2),result(1,1,3),
     &  result(1,1,4),rbc,zbs,rmnaxis,zmnaxis)
        return
        end
c******************
      BLOCK DATA case17
        common/scmspc/gnorm,specw,delt,deltf
        data deltf/1.0/
      end

c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scfixa
        subroutine scfixa(pexp,qexp,nthetain,mpolin,nphin,nfpin)
        include 'name1.m'
c**********************************************************************
c       This routine stores toroidal and poloidal mode number arrays
c
c       MPOL = NUMBER OF POLOIDAL MODES USED IN CURVE-FIT
c       NPHI = NUMBER OF EQUALLY-SPACED PHI PLANES PER FIELD
c              PERIOD IN WHICH THE CURVE-FITTING IS PERFORMED
c       IMPORTANT: NPHI MUST BE EVEN FOR NON-SYMMETRIC SYSTEMS

c       MPNT = NUMBER OF R AND Z MODES APPEARING IN FIT
c       NTHETA=NUMBER OF THETA POINTS IN EACH PHI PLANE
c       N2   = TOTAL NUMBER OF RESIDUALS IN FSQ (PER PHI PLANE)
c       NFP  = NUMBER OF TOROIDAL FIELD PERIODS (IRRELEVANT)
c**********************************************************************
        ntheta = nthetain
        mpol = mpolin
        mpol2= 2*mpol
        mpol3= mpol + mpol2
        mpol4= mpol + mpol3
        nphi = nphin
        nphi2= 1 + nphi/2
        nfp  = nfpin
        n2 = 4*mpol + ntheta

        l = 0
        ntor = max0(1,nphi-1)
        nn0 = 1 - (ntor+1)/2
        do 10 n=1,ntor
 10     nn(n) = (nn0 + (n-1))*nfp
        do 40 m=1,mpol
        mm(m)=m-1
        do 40 n=1,ntor
        if( mm(m).eq.0 .and. nn(n).lt.0 )go to 40
        l=l+1
        go to 30
 20     l=l+1
 30     continue
        m1(l)=mm(m)
        n1(l)=nn(n)
 40     continue
        mpnt=l
        dnorm=2./ntheta

        do 50 m = 1,mpol
        dm1(m) = float(m-1)
        faccon(m) = .125*dnorm/(1.+dm1(m))**pexp
        if( m.eq.1 )go to 50
        xmpq(m,1) = dm1(m)**pexp
        xmpq(m,2) = dm1(m)**qexp
        xmpq(m,4) = xmpq(m,1)
        if( m.le.2 )xmpq(m,4) = 0.
        xmpq(m,3) = sqrt(xmpq(m,4))
 50     continue
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scgues
        subroutine scgues(rin,zin)
        include 'name1.m'
        real rin(ntheta,nphi),zin(ntheta,nphi)
c**********************************************************************
c       This subroutine obtains initial guesses for the centroid at each
c       toroidal plane. By default, the polar axis is set equal to this
c       centroid.  It is imperative that the polar axis be enclosed by
c       the surface (otherwise the computed theta angles will not span
c       [0,2pi]). For certain complex cross-sections, this simple estimate
c       of the polar axis may fail, and the user must supply values for
c       raxis(i),zaxis(i).  In addition, for non-starlike domains, the
c       points along the curve must be monitonically increasing as a
c       function of the arclength from any fixed point on the curve. This
c       ordering is attempted by the subroutine ORDER, but may fail in
c       general if too few theta points are used.
c**********************************************************************
c       COMPUTE CENTROID
c**********************************************************************
        do 10 i = 1,nphi
          r0n(i) = 0.0
          z0n(i) = 0.0
          do 20 j = 1,ntheta
            r0n(i) = r0n(i)+rin(j,i)/float(ntheta)
            z0n(i) = z0n(i)+zin(j,i)/float(ntheta)
  20      continue
        raxis(i) = r0n(i)
 10     zaxis(i) = z0n(i)
c**********************************************************************
c       ORDER POINTS ON FLUX SURFACE AT EACH TOROIDAL ANGLES
c**********************************************************************
        write (6,*) 'ORDERING SURFACE POINTS'
        do 30 i = 1,nphi
 30     call scorder(rin(1,i),zin(1,i),raxis(i),zaxis(i))
c**********************************************************************
c       COMPUTE OPTIMUM ANGLE OFFSETS FOR M=1 MODES
c**********************************************************************
        call scgeta(rin,zin,angle,r0n,z0n)
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scorder
        subroutine scorder(rval,zval,xaxis,yaxis)
        include 'name1.m'
        real rval(*),zval(*),newdist,tempr(nu),tempz(nu)
c**********************************************************************
c       Program ORDER : orders points on a magnetic surface at a
c       fixed toroidal plane and assigns right-handed circulation
c       around flux surface. XAXIS, YAXIS:  Polar-type axis (must lie
c       inside curve to check sign of rotation)
c**********************************************************************
        olddist = 1.e20
        do 10 i = 1,ntheta-1
        ip1 = i + 1
        i1 = i
        shortest = 1.e20
 15              do 20 j = ip1,ntheta
                 if(i1.gt.1)olddist = (rval(i1-1)-rval(j))**2
     &                              + (zval(i1-1)-zval(j))**2
                 newdist =(rval(i1)-rval(j))**2 + (zval(i1)-zval(j))**2
                 if((newdist.le.olddist).and.
     &           (newdist.lt.shortest))then
                         next = j
                         shortest = newdist
                 endif
 20     continue
c**********************************************************************
c       Swap nearest point (next) with current point (ip1)
c**********************************************************************
        if(shortest.ge.1.e10)then
        saver = rval(i-1)
        rval(i-1) = rval(i)
        rval(i)= saver
        savez = zval(i-1)
        zval(i-1) = zval(i)
        zval(i)= savez
        i1 = i1 - 1
        ip1 = ip1 - 1
        go to 15
        endif
        saver = rval(ip1)
        rval(ip1) = rval(next)
        rval(next)= saver
        savez = zval(ip1)
        zval(ip1) = zval(next)
        zval(next)= savez
 10     continue
c**********************************************************************
c       Check that xaxis,yaxis is inside surface and
c       ascertain that the angle rotates counterclockwise
c       using Cauchy's theorem in "complex"-plane
c**********************************************************************
        residue = 0.
        do 30 i = 1,ntheta-1
        x = .5*(rval(i)+rval(i+1)) - xaxis
        y = .5*(zval(i)+zval(i+1)) - yaxis
        dx= rval(i+1) - rval(i)
        dy= zval(i+1) - zval(i)
        residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.e-10)
 30     continue
        x = .5*(rval(1)+rval(ntheta)) - xaxis
        y = .5*(zval(1)+zval(ntheta)) - yaxis
        dx= rval(1) - rval(ntheta)
        dy= zval(1) - zval(ntheta)
        residue = residue + (x*dy - y*dx)/(x**2 + y**2 + 1.e-10)

        if( residue .lt.(-.90*twopi))then
        do 40 i = 2,ntheta
        j = ntheta - i + 2
        tempr(i) = rval(j)
        tempz(i) = zval(j)
 40     continue
        do 50 i = 2,ntheta
        rval(i) = tempr(i)
        zval(i) = tempz(i)
 50     continue
        else if( abs(residue) .lt. (.90*twopi) )then
c
        write (6,*)' The magnetic axis is not enclosed by boundary '
        stop
        endif

        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scgeta
        subroutine scgeta(rval,zval,xpts,rcenter,zcenter)
        include 'name1.m'
        real rval(ntheta,*),zval(ntheta,*),xpts(ntheta,*),
     &  rcenter(*),zcenter(*),rcos(nv),rsin(nv),
     &  zcos(nv),zsin(nv),phiangle(nv)
c**********************************************************************
c       Compute angle offset consistent with constraint Z1n = Z1,-n
c       Note: This is done iteratively, since elongation is unknown
c**********************************************************************
        do 10 i = 1,nphi
        do 10 j = 1,ntheta
 10     xpts(j,i) = twopi*(j-1)/float(ntheta)
        do 100 iterate = 1,5
        do 20 i = 1,nphi
        rcos(i) = 0.
        rsin(i) = 0.
        zcos(i) = 0.
        zsin(i) = 0.
        do 20 j = 1,ntheta
        xc = rval(j,i) - rcenter(i)
        yc = zval(j,i) - zcenter(i)
        rcos(i) = rcos(i) + cos(xpts(j,i))*xc
        rsin(i) = rsin(i) + sin(xpts(j,i))*xc
        zcos(i) = zcos(i) + cos(xpts(j,i))*yc
        zsin(i) = zsin(i) + sin(xpts(j,i))*yc
 20     continue
c**********************************************************************
c       Compute new angles starting from offset phiangle(i)
c**********************************************************************
        elongate = ssum(nphi,zsin,1)/ssum(nphi,rcos,1)
        delangle = 0.
        do 30 i = 1,nphi
        phiangle(i) = atan2(elongate*zcos(i)-rsin(i),
     &                     elongate*zsin(i)+rcos(i))
        delangle = max(delangle,abs(phiangle(i)))
        do 30 j = 1,ntheta
 30     xpts(j,i) = xpts(j,i) + phiangle(i)
        if(delangle.lt.0.02)go to 40
 100    continue
 40     write (6,*) ' Average elongation = ',elongate
c
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scampl
        subroutine scampl(rcenter,zcenter,angin,rmc,rms,
     &  zmc,zms,xpts,xin,yin)
        include 'name1.m'
        real rmc(*),rms(*),zmc(*),zms(*),angin(*),xpts(*),xin(*),yin(*)
c*****************************************************************
c       This subroutine assigns initial guesses for angles and
c       Fourier mode amplitudes to the appropriate components of
c       the xvec array
c*****************************************************************
        rmc(1) = rcenter
        zmc(1) = zcenter
        call scopy(ntheta,angin,1,xpts,1)
        xmult = 2./float(ntheta)
        do 10 j = 1,ntheta
        arg = angin(j)
        xi = xmult*(xin(j) - rcenter)
        yi = xmult*(yin(j) - zcenter)
        do 10 m = 2,mpol
        rmc(m) = rmc(m) + cos((m-1)*arg)*xi
        rms(m) = rms(m) + sin((m-1)*arg)*xi
        zmc(m) = zmc(m) + cos((m-1)*arg)*yi
 10     zms(m) = zms(m) + sin((m-1)*arg)*yi
        if( abs(elongate*zmc(2)-rms(2)).gt.
     &  1.e-5*(abs(elongate*zms(2))+abs(rms(2))) )then
        write (6,*) ' Incorrect initial angle assignment'
        stop
        endif
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scrstrt
        subroutine scrstrt(irst,nresets)
        include 'name1.m'
c
c**********************************************************************
c       This routine either stores an accepted value of the local solution
c       (irst=1) or reset the present solution to a previous value (irst=2)
c**********************************************************************
        go to (15,25)irst
 
 15     call scopy(n2,xvec,1,xstore,1)
        return
 
 25     do 20 l = 1,n2
        xdot(l) = 0.
 20     xvec(l) = xstore(l)
        delt = .95 * delt
        irst = 1
        nresets = nresets + 1
        if ( nresets .ge. 100 )then
        write (6,*) ' Time step reduced 100 times without convergence'
        stop
        endif
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scevol
        subroutine scevol(g11)
        include 'name1.m'
        data bmax /0.15/
 
        ftest=gnorm/g11
        dtau=abs(1.-ftest)
        g11=gnorm
        otav=dtau/delt
        dtau=delt*otav+.001
        if(dtau.gt.bmax)dtau=bmax
        b1=1.-.5*dtau
        fac=1./(1.+.5*dtau)
        do 10 i=1,n2
        xdot(i)=fac*(xdot(i)*b1-delt*gvec(i))
 10     xvec(i)=xvec(i)+xdot(i)*delt
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scfnct
        subroutine scfnct(rmc,rms,zmc,zms,xpts,grc,grs,
     &  gzc,gzs,gpts,fsq,r10,xin,yin,mmax)
        include 'name1.m'
        real yc(mu),ys(mu),gcon(nu),gtt(nu),r1(nu),
     &  z1(nu),rt1(nu),zt1(nu),arg(nu),rcon(nu),
     &  zcon(nu),cosa(nu,mu),sina(nu,mu),rmc(*),rms(*),
     &  zmc(*),zms(*),xpts(*),grc(*),grs(*),gzc(*),gzs(*),gpts(*),
     &  xin(*),yin(*)
        fsq=0.
        denom=0.
        gnorm=0.
        specw = 0.
        do 10 j=1,ntheta
        r1(j)=-xin(j)
        z1(j)=-yin(j)
        rcon(j) = 0.
        zcon(j) = 0.
        rt1(j)  = 0.
 10     zt1(j)  = 0.
        do 20 l=1,n2
 20     grc(l)=0.
c**********************************************************************
c       COMPUTE SPECTRAL WIDTH OF CURVE
c**********************************************************************
        do 30 m=2,mmax
        t2=(rmc(m)*rmc(m)+zmc(m)*zmc(m)+
     &      rms(m)*rms(m)+zms(m)*zms(m))*xmpq(m,1)
        denom=denom+t2
 30     specw=specw+xmpq(m,2)*t2
        specw=specw/denom
c**********************************************************************
c       COMPUTE CURVE AND CONSTRAINT FORCES
c**********************************************************************
        do 40 m=1,mmax
        do 40 l=1,ntheta
        arg(l)    = dm1(m) * xpts(l)
        cosa(l,m) = cos(arg(l))
        sina(l,m) = sin(arg(l))
        gtt(l)  = rmc(m)*cosa(l,m) + rms(m)*sina(l,m)
        gcon(l) = zmc(m)*cosa(l,m) + zms(m)*sina(l,m)
        r1(l)   = r1(l)  + gtt(l)
        z1(l)   = z1(l)  + gcon(l)
        rt1(l)  = rt1(l) + dm1(m)*(rms(m)*cosa(l,m) - rmc(m)*sina(l,m))
        zt1(l)  = zt1(l) + dm1(m)*(zms(m)*cosa(l,m) - zmc(m)*sina(l,m))
        rcon(l) = rcon(l) + xmpq(m,4)*gtt(l)
        zcon(l) = zcon(l) + xmpq(m,4)*gcon(l)
 40     continue
        do 50 l=1,ntheta
        gtt(l)  = rt1(l)*rt1(l) + zt1(l)*zt1(l)
        gpts(l) = .25*(r1(l)*rt1(l)+z1(l)*zt1(l))/gtt(l)
        gcon(l) = rcon(l)*rt1(l) + zcon(l)*zt1(l)
 50     continue
c**********************************************************************
c       COMPUTE MEAN-SQUARE DEVIATION BETWEEN POINTS AND FIT
c**********************************************************************
        do 60 l=1,ntheta
 60     fsq = fsq + .5*(r1(l)**2 + z1(l)**2)
        fsq = sqrt(dnorm*fsq)/r10
c**********************************************************************
c       FILTER CONSTRAINT FORCE TO REMOVE ALIASING
c**********************************************************************
        tcon = 1./amaxaf(gtt,1,ntheta)/sqrt(1.+xmpq(mmax,3))
        do 70 m=2,mmax-1
        yc(m) = sdot(ntheta,cosa(1,m),1,gcon,1)*faccon(m)*tcon
 70     ys(m) = sdot(ntheta,sina(1,m),1,gcon,1)*faccon(m)*tcon
        do 80 l = 1,ntheta
 80     gcon(l) = 0.
        do 90 m=2,mmax-1
        do 90 l = 1,ntheta
 90     gcon(l) = gcon(l) + yc(m)*cosa(l,m) + ys(m)*sina(l,m)
c**********************************************************************
c       ADD CURVE AND CONSTRAINT FORCES
c**********************************************************************
        do 100 m=1,mmax
        do 100 l=1,ntheta
        rcon(l) = r1(l)  + gcon(l)*rt1(l)*xmpq(m,3)
        zcon(l) = z1(l)  + gcon(l)*zt1(l)*xmpq(m,3)
        grc(m)  = grc(m) + cosa(l,m)*rcon(l)
        gzc(m)  = gzc(m) + cosa(l,m)*zcon(l)
        grs(m)  = grs(m) + sina(l,m)*rcon(l)
 100    gzs(m)  = gzs(m) + sina(l,m)*zcon(l)
c**********************************************************************
c       COMPUTE m=1 CONSTRAINT
c**********************************************************************
        gzc(2) = (elongate*grs(2) + gzc(2))/(1.+elongate**2)
        grs(2) = gzc(2)*elongate

        call sscal(mpol4,dnorm,grc,1)
        grc(1) = 0.5*grc(1)
        do 110 j=1,mmax
 110    gnorm = gnorm + grc(j)*grc(j) + gzc(j)*gzc(j)
     &                + grs(j)*grs(j) + gzs(j)*gzs(j)
        gnorm = gnorm/r10**2
        do 120 j=1,ntheta
 120    gnorm = gnorm + dnorm*gpts(j)*gpts(j)
        do 130 m = 2,mmax
        grc(m) = grc(m)/(1.+tcon*xmpq(m,3))
        gzc(m) = gzc(m)/(1.+tcon*xmpq(m,3))
        grs(m) = grs(m)/(1.+tcon*xmpq(m,3))
 130    gzs(m) = gzs(m)/(1.+tcon*xmpq(m,3))
        return
        end
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c@scfftr
        subroutine scfftr(rmc,rms,zmc,zms,rbc,zbs,rmnaxis,zmnaxis)
        include 'name1.m'
        real rmc(nphi2,*),rms(nphi2,*),zmc(nphi2,*),zms(nphi2,*),
     &  intgrate(nv2),argi(nv2),rbc(0:mpol-1,-nphi2:nphi2),
     &  zbs(0:mpol-1,-nphi2:nphi2),rmnaxis(0:nphi2),zmnaxis(0:nphi2)
c**********************************************************************
c       PERFORM FOURIER TRANSFORM IN phi
c**********************************************************************
        delphi = 2./float(nphi)
        do 20 i = 1,nphi2
        intgrate(i) = delphi
        argi(i) = twopi*(i-1)/float(nphi*nfp)
 20     if(i.eq.1.or.i.eq.nphi2)intgrate(i) = .5*delphi
        do 70 mn=1,mpnt
        mreal = m1(mn)
        nreal = n1(mn)/nfp
        m = mreal + 1
        dn= float(n1(mn))
        rbc(mreal,nreal) = 0.
        zbs(mreal,nreal) = 0.
        do 40 i = 1,nphi2
        arg = dn*argi(i)
        tcosn = cos(arg)
        tsinn = sin(arg)
        rbc(mreal,nreal) = rbc(mreal,nreal)
     &                   + intgrate(i)*(tcosn*rmc(i,m) + tsinn*rms(i,m))
 40     zbs(mreal,nreal) = zbs(mreal,nreal)
     &                   + intgrate(i)*(tcosn*zms(i,m) - tsinn*zmc(i,m))
        if(rbc(mreal,nreal).eq.(0.).and.zbs(mreal,nreal).eq.0.)go to 70
        if(mreal.eq.0.and.nreal.eq.0)then
        rmnaxis(0) = sdot(nphi2,intgrate,1,raxis,1)
        zmnaxis(0) = 0.
        else if(mreal.eq.0.and.nreal.gt.0)then
        rmnaxis(nreal) = 0.
        zmnaxis(nreal) = 0.
        do 50 i = 1,nphi2
        rmnaxis(nreal) = rmnaxis(nreal) 
     &                 + 2.*intgrate(i)*raxis(i)*cos(dn*argi(i))
 50     zmnaxis(nreal) = zmnaxis(nreal)
     &                 - 2.*intgrate(i)*zaxis(i)*sin(dn*argi(i))
        endif
 70     continue
 80     format(i5,i4,1p4e12.4)
        return
        end

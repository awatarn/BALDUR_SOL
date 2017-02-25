c@realft  Glenn Bateman, 20 April 1996
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c     *********************
c     * SUBROUTINE realft *
c     *********************
c
c     Usage: 
c  call realft ( pc0, pcm, psm, kcsdim, nharm, pfthet, ntheta, ixfour)
c
c     Purpose:  Real-valued Fourier analysis.
c
c  ixfour .ge. 0:  Convert from fourier harmonic to angle representation
c                    pc0, pcm(m), psm(m) --> pfthet(jt)
c
c  ixfour .lt. 0:  Convert from angle to Fourier harmonic representation
c                    pfthet(jt) --> pc0, pcm(m), psm(m)
c
c     where:
c
c  pfthet(jt)
c       = pc0 + sum_{m=1}^{nharm} [   pcm(1,m) * cos ( m * theta(jt) )
c                                   + psm(1,m) * sin ( m * theta(jt) ) ]
c
c  for theta(jt) = 2 * pi * (jt-1) / (ntheta-1)  jt = 1, ntheta
c
c  kcsdim = first dimension of the 2-D arrays pcm(j,m) and psm(j,m)
c
c  Note:  ntheta must be a power of 2
c         must have  8 .le. ntheta .le. 64
c
c         nharm must be .le. ( ntheta / 2 ) - 1
c
c     This routine makes use of the Fast Fourier Transform routine
c  fft1 found in Numerical Recipies as sbrtn four1.
c
c
c     Glenn Bateman  16-sep-88  Checked using stand-alone program.
c
c**********************************************************************c
c
      subroutine realft (pc0,pcm,psm,kcsdim,nharm,pfthet,ntheta,ixfour)
c
      implicit none
c
      integer knf, kcsdim, nharm, ntheta, ixfour
c
      parameter  ( knf=64 )
c
      real  pcm(kcsdim,*), psm(kcsdim,*), pfthet(*), ztemp(2*knf)
c
      real  pc0, zepslon, zpfmax, zpfeps, znorm, one
c
      integer inital, ipower, inf, jm, jt, ir
c
      data  inital /0/
c
      save inital, zepslon, one
c
c..initialize
c
      if ( inital .ne. ntheta ) then
c
c  check to make sure ntheta is a power of 2 and .le. knf
c
        ipower = 4
        inf = knf / 4
        do while ( inf .gt. 1 )
          ipower = ipower * 2
          inf    = inf / 2
          if ( ntheta .eq. ipower ) go to 12
        enddo
        call abortb (6
     &    ,'ntheta not power of 2 .le. knf in sbrtn realft')
  12    continue
c
c..compute machine epslon
c
        one = 1.0
        zepslon = one
        do while ( one + zepslon .gt. one )
          zepslon = zepslon * 0.1
        enddo
        zepslon = zepslon * 10.0
c
        inital = ntheta
c
      endif
c
c..check to make sure nharm is .le. ( ntheta / 2 ) - 1
c
      if ( nharm .gt. ( ntheta / 2 ) - 1 )
     & call abortb (6,'nharm .gt. ntheta/2 in sbrtn realft')
c
c..initialize zxmplx(jm) and zymplx(jm)
c
      do jm=1,2*ntheta
        ztemp(jm) = 0.0
      enddo
c
c..IXFOUR .ge. 0 ==> Compute pfthet(jt) from pc0, pcm(m), psm(m)
c
      if ( ixfour .ge. 0 ) then
c
          ztemp(1) = pc0
c
        ir = 2 * ntheta + 1
        do jm=1,nharm
          ztemp(2*jm+1)    =   0.5 * pcm(1,jm)
          ztemp(2*jm+2)    = - 0.5 * psm(1,jm)
          ztemp(ir-2*jm)   =   0.5 * pcm(1,jm)
          ztemp(ir-2*jm+1) =   0.5 * psm(1,jm)
        enddo
c
        call fft1 (ztemp,ntheta,-1)
c
          zpfmax = abs ( ztemp(1) )
        do jt=1,ntheta
          pfthet(jt) = ztemp(2*jt-1)
          zpfmax = max ( zpfmax, abs ( pfthet(jt) ) )
        enddo
c
          zpfeps = zpfmax * zepslon
        do jt=1,ntheta
          if ( abs ( pfthet(jt) ) .lt. zpfeps ) pfthet(jt) = 0.0
        enddo
c
c..IXFOUR .lt. 0 ==> Compute pc0, pcm(1,m), psm(1,m) from pfthet(jt)
c
      else
c
          znorm = 1. / ntheta
        do jt=1,ntheta
          ztemp(2*jt-1) = znorm * pfthet(jt)
        enddo
c
        call fft1 (ztemp,ntheta,1)
c
          pc0 = ztemp(1)
c
          zpfmax = abs ( pc0 )
        do jm=1,nharm
          pcm(1,jm) =  2. * ztemp(2*jm+1)
          psm(1,jm) = -2. * ztemp(2*jm+2)
          zpfmax = max ( zpfmax, abs(pcm(1,jm)), abs(psm(1,jm)) )
        enddo
c
          zpfeps = zpfmax * zepslon
          if ( abs(pc0) .lt. zpfeps ) pc0 = 0.0
        do jm=1,nharm
          if ( abs(pcm(1,jm)) .lt. zpfeps ) pcm(1,jm) = 0.0
          if ( abs(psm(1,jm)) .lt. zpfeps ) psm(1,jm) = 0.0
        enddo
c
      endif
c
      return
      end

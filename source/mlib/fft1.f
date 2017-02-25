c@fft1  Fast Fourier Transform,  Glenn Bateman
c  pis 18-jun-98 Actually implemented the changes described below
c  rgb 03-jun-96 dsin --> sin, sngl(wr) --> wr, sngl(wi) --> wi
c    real*8 --> real  (use -r8 option on 32 bit workstations)
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine fft1 ( data, nn, isign )
c
c  data(j) = real array with 2*nn elements on input
c            is replaced by its discrete Fourier transform on output
c
c  nn      = number of harmonics.  MUST be an integer power of two.
c  isign   = +1 to replace data by its discrete Fourier transform
c          = -1 to replace data by nn times its inverse discrete
c               Fourier transform
c
c  based on sbrtn four1 from Numerical Recipies, 1986, p394
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      integer nn, isign, n, j, i, m, mmax, istep
c
      real data(2*nn), tempr, tempi
c
c  use double precision for the trigonometric recurrences
c
      real  wr, wi, wpr, wpi, wtemp, theta
c
      n = 2*nn
      j = 1
c
c..bit reversal
c
      do i=1,n,2
        if ( j .gt. i ) then
          tempr = data(j)
          tempi = data(j+1)
          data(j) = data(i)
          data(j+1) = data(i+1)
          data(i) = tempr
          data(i+1) = tempi
        endif
c
        m = n / 2
   1    continue
        if ( ( m .ge. 2 )  .and.  ( j .gt. m ) ) then
          j = j - m
          m = m / 2
          go to 1
        endif
        j = j + m
      enddo
c
c..Danielson-Lanczos algorithm
c
      mmax = 2
   2  continue
c
      if ( n .gt. mmax ) then
        istep = 2 * mmax
        theta = 6.28318530717959 / ( isign * mmax )
        wpr   = - 2.0 * sin ( 0.5 * theta )**2
        wpi   = sin ( theta )
        wr    = 1.0
        wi    = 0.0
        do m=1,mmax,2
          do i=m,n,istep
            j = i + mmax
            tempr = wr * data(j) - wi * data(j+1)
            tempi = wr * data(j+1) + wi * data(j)
            data(j) = data(i) - tempr
            data(j+1) = data(i+1) - tempi
            data(i) = data(i) + tempr
            data(i+1) = data(i+1) + tempi
          enddo
          wtemp = wr
          wr = wr * wpr - wi * wpi + wr
          wi = wi * wpr + wtemp * wpi + wi
        enddo
        mmax = istep
        go to 2
      endif
c
      return
      end
      

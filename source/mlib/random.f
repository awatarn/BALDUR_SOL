c   File random.f  Random number generators from Numerical Recipies
c
c@ran1   .../baldur/code/util/random.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      real function ran1( idum )
c
c   Returns a uniform random deviate between 0.0 and 1.0.
c  Set idum < 0 to initialize or reinitialize the sequence.
c
      implicit none
c
      integer idum, m1, m2, m3, ia1, ia2, ia3, ic1, ic2, ic3, iff
     &  , ix1, ix2, ix3, j
c
      real r(97), rm1, rm2
c
      parameter ( m1=259200, ia1=7141, ic1=54773, rm1=1.0/m1 )
      parameter ( m2=134456, ia2=8121, ic2=28411, rm2=1.0/m2 )
      parameter ( m3=243000, ia3=4561, ic3=51349 )
c
      data iff / 0 /
c
      save iff, ix1, ix2, ix3
      save r
c
c..initialize on first call or if idum < 0
c
      if ( idum .lt. 0  .or.  iff .eq. 0 ) then
        iff = 1
c
c..seed the first routine
c
        ix1 = mod ( ic1 - idum, m1 )
        ix1 = mod ( ia1 * ix1 + ic1, m1 )
c
c..use ix1 to seed the second and third routines
c
        ix2 = mod ( ix1, m2 )
        ix1 = mod ( ia1 * ix1 + ic1, m1 )
        ix3 = mod ( ix1, m3 )
c
c..fill table with sequential uniform deviates generated
c  by the first two routines
c
        do j=1,97
          ix1 = mod ( ia1 * ix1 + ic1, m1 )
          ix2 = mod ( ia2 * ix2 + ic2, m2 )
c
c..low and high order pieces combined here
c
          r(j) = ( float(ix1) + float(ix2)*rm2 ) * rm1
        enddo
        idum = 1
      endif
c
c..start here when not initializing.  Generate next number in sequence
c
      ix1 = mod ( ia1 * ix1 + ic1, m1 )
      ix2 = mod ( ia2 * ix2 + ic2, m2 )
      ix3 = mod ( ia3 * ix3 + ic3, m3 )
c
c  use the third sequence to get an integer between 1 and 97
c
      j = 1 + ( 97 * ix3 ) / m3
c
      if ( j .gt. 97  .or.  j .lt. 1 ) then
        write (*,*) 'Abort:  error in fnc ran1, j = ',j
        stop
      endif
c
c  return that table entry
c
      ran1 = r(j)
c
      if ( ran1 .lt. 0.0  .or.  ran1 .gt. 1.0 ) then
        write (*,*) 'Abort: ran1 < 0.0 or ran1 > 1.0'
        stop
      endif
c
c  then refill table entry
c
      r(j) = ( float(ix1) + float(ix2)*rm2 ) * rm1
c
      return
      end
c@ran2   .../baldur/code/util/random.f
c rgb 17-jun-96 removed diagnostic printout
c rgb 17-jun-96 modified range of results to  rm .le. ran2 .le. 1.0
c rgb 14-jun-96 added diagnostic printout
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      real function ran2( idum )
c
c  Faster routine
c  Returns a uniform random deviate between 0.0 and 1.0.
c  Set idum < 0 to initialize or reinitialize the sequence.
c
      implicit none
c
      integer idum, m, ia, ic, iy, iff, icount, j, ir(97)
c
      real rm
c
      parameter ( m=714025, ia=1366, ic=150889, rm=1.0/(m+1) )
c
      data iff / 0 /, icount /0/
      save iff, iy, ir
c
c..initialize on first call or if idum < 0
c
      if ( idum .lt. 0  .or.  iff .eq. 0 ) then
        iff = 1
c
c..seed the first routine
c
        idum = mod ( ic - idum, m )
c
c..Initialize the shuffle table
c
        do j=1,97
          idum = mod ( ia * idum + ic, m )
          ir(j) = idum
        enddo
c
        idum = mod ( ia * idum + ic, m )
        iy   = idum
      endif
c
c..start here when not initializing.  Generate next number in sequence
c
      j = 1 + ( 97 * iy ) / m
c
      if ( j .gt. 97  .or.  j .lt. 1 ) then
        write (*,*) 'Abort:  error in fnc ran1, j = ',j
        stop
      endif
c
c  return that table entry
c
      iy = ir(j)
      ran2 = ( iy + 1 ) * rm
c
      icount = icount + 1
cbate      write (6,110) icount, idum, ran2
cbate 110  format ('@r',2x,i12,2x,i6,2x,0pf20.16)
c
      if ( ran2 .lt. 0.0  .or.  ran2 .gt. 1.0 ) then
        write (*,*) 'Abort: ran2 < 0.0 or ran2 > 1.0'
        write (*,*) ' ran2 = ',ran2
        write (*,*) ' j    = ',j
        write (*,*) ' iy   = ',iy
        write (*,*) ' idum = ',idum
        stop
      endif
c
c  then refill table entry
c
      idum = mod ( ia * idum + ic, m )
      ir(j) = idum
c
      return
      end
c@ran3   .../baldur/code/util/random.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      real function ran3( idum )
c
c  Knuth's portable random number generator.
c  Returns a uniform random deviate between 0.0 and 1.0.
c  Set idum < 0 to initialize or reinitialize the sequence.
c
      implicit none
c
      integer idum, i, ii, k, iff
      integer inext, inextp
c
      real fac
c
      integer  mbig, mseed, mz, mj, mk, ma(55)
c
      parameter ( mbig = 1000000000, mseed = 161803398
     &  , mz = 0, fac = 1.0/mbig )
c
c..substitute the following two lines for completely floating point
c
c      real  mbig, mseed, mz, mj, mk, ma(55)
c      paramter ( mbig=4000000., mseed=1618033., mz=0., fac=1.0/mbig )
c
c..any large mbig and somewhat smaller mseed may be substituted
c
      data iff / 0 /
c
      save iff, inext, inextp
      save ma
c
c..initialize on first call or if idum < 0
c
      if ( idum .lt. 0  .or.  iff .eq. 0 ) then
        iff = 1
c
c..initialize ma(55) using the seed idum and large number mseed
c
        mj = mseed - iabs(idum)
        mj = mod ( mj, mbig )
        ma(55) = mj
c
c..initialize the rest of the table
c
        mk = 1
        do i=1,54
          ii = mod ( 21*i, 55 )
          ma(ii) = mk
          mk = mj - mk
          if ( mk .lt. mz ) mk = mk + mbig
          mj = ma(ii)
        enddo
c
c..randomize them by warming up the generator
c
        do k=1,4
          do i=1,55
            ma(i) = ma(i) - ma( 1 + mod ( i + 30, 55 ) )
            if ( ma(i) .lt. mz ) ma(i) = ma(i) + mbig
          enddo
        enddo
c
c..prepare indicies for first generated number (note, 31 is special)
c
        inext = 0
        inextp = 31
        idum = 1
      endif
c
c..start here when not initializing.  Generate next number in sequence
c  Increment inext wrapping around 56 to 1
c
      inext = inext + 1
      if ( inext .eq. 56 ) inext = 1
      inextp = inextp + 1
      if ( inextp .eq. 56 ) inextp = 1
c
c..generate a new random number subtractively
c
      mj = ma(inext) - ma(inextp)
c
c..be sure it is in range, then store it and output uniform deviate
c
      if ( mj .lt. mz ) mj = mj + mbig
      ma(inext) = mj
      ran3 = mj * fac
c
      if ( ran3 .lt. 0.0  .or.  ran3 .gt. 1.0 ) then
        write (*,*) 'Abort: ran3 < 0.0 or ran3 > 1.0'
        stop
      endif
c
      return
      end

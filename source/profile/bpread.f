c@bpread.f   02 Dec 1996   Glenn Bateman, Lehigh, bateman@pppl.gov
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine bpread ( cvarname, kunit, kxdim, ptime
     & , prad, krad, cradname, profile, kerr )
c
c   Extract a radial profile time-slice from a 2-D ASCII input file.
c
c  cvarname  = name of the independent 2-D variable parray      (input)
c  kunit     = input unit number to read input file (default 70)(input)
c  kxdim     = maximum number of elements allowed in x array    (input)
c  ptime     = time at which profile is to be extracted         (input)
c
c  prad(jx)  = 1-D independent variable x array from input file (output)
c  krad      = number of elements in x array from input file    (output)
c  cradname  = name of radial 1-D array                         (output)
c  profile(jx) = 1-D profile array from input file              (output)
c
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the input file cannot be opened or read
c            = 2 if there are too many elements in the arrays
c                (that is, if krad > kxdim)
c            = 3 if there are no elements in the arrays
c
c  Routines used:
c
c    br2read
c    int1d
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c     implicit none
c
      integer kr, kt, krt
c
      parameter (kr=200, kt=4000, krt=kr*kt)
c
c  kr  = max number of radial points
c  kt  = max number of time points
c
      real, INTENT(out)         :: profile(*), prad(*)
      integer, INTENT(out)      :: kerr, krad
      real, INTENT(in)          :: ptime
      integer, INTENT(in)       :: kxdim, kunit
      character*(*), INTENT(in) :: cvarname
      character*(*), INTENT(out):: cradname

c
      real z2d(kr,kt), zrad(kr), ztime(kt), zt(1), zp(1)
c
      integer kydim, kyu, kdchar, kfchar
c
      integer ierr, ixdim, iydim, irad, itime, ixlast, jx
c
c  ierr    = input value of kerr
c  ixdim   = first dimension of z2d and zrad arrays
c  iydim   = second dimension of z2d and ztime arrays
c  irad    = actual number of radial points in 2-D input file
c  itime   = actual number of time values in 2-D input file
c  ixlast  = last index used in interpolation routine
c
c
      character ctimname*32
c
c
cbate      write ( 6, *)
cbate      write ( 6, *) 'kunit = ',kunit
cbate      write ( 6, *) 'kxdim  = ',kxdim
c
c..set defaults
c
      kerr   = 0
      ierr   = kerr
c
c..the following array dimensions are fixed above
c
      ixdim  = kr
      iydim  = kt
c
c..clear arrays
c
      do jx = 1, kxdim
        prad(jx) = 0.
        profile(jx) = 0.
      enddo
c
c..read the 2-D input file
c
      call br2read ( cvarname, kunit, ixdim, iydim
     &  , zrad, irad, cradname, ztime, itime, ctimname, z2d, kerr )
c
c      write (*,*)
c      write (*,*) 'itime = ',itime
c      write (*,*) 'ztime(itime) = ',ztime(itime)
c      write (*,*) 'j ',cradname, cvarname
c      do jx=1,irad
c         write (*,102) jx, zrad(jx), z2d(jx,itime)
c      enddo
c 102  format (i4,1p2e15.5)
c
c..error checks
c
      if ( irad .lt. 1  .or.  itime .lt. 1 ) kerr = 3
c
      if ( kerr .gt. 0 ) return
c
c..interpolate profile at ptime
c
        krad   = irad
        zt(1)  = ptime
        ixlast = 1
c
      if ( itime .lt. 2 ) then
        do jx=1,krad
          prad(jx) = zrad(jx)
          profile(jx) = z2d(jx,1)
        enddo
      else
c
        do jx=1,krad
c
          call int1d ( 1, 0, 0, itime, ixdim, ztime, z2d(jx,1)
     &      , 1, zt, ixlast, zp )
c
          prad(jx)    = zrad(jx)
          profile(jx) = zp(1)
c
        enddo
c
      endif
c
c
      if ( irad .lt. 9999) return
c
      return
      end


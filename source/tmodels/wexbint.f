      subroutine wexbint (kxdim, nxwexba, xwexba, ntwexba, twexba 
     & , wexba, nxb, xb, time, wexbxb)
c-----------------------------------------------------------------------
c Subroutine for interpolating a 2-D array defined on a (xwexba, twexba)
c grid to a single vector (xb) at a specified time. The routine uses an 
c implementation of Akima's shape preserving 1-D algorithm for 
c interpolation by Glenn Bateman with some extrapolation features added.
c-----------------------------------------------------------------------
c Wexbint: Written by P. Strand 29-may-1998
c-----------------------------------------------------------------------
c Dependencies: int1d.f by G. Bateman
c-----------------------------------------------------------------------

      Implicit None

c-----------------------------------------------------------------------
c Declaration of input variables
c-----------------------------------------------------------------------

      integer kxdim         !Leading order of input array wexba(kxdim,*)
      integer nxwexba       !Length of input x vector xwexba(nxwexba)
      integer ntwexba       !Length of input x vector xwexba(nxwexba)
      integer nxb           !Length of xb -vector (interpolates to this)

      real xwexba(nxwexba)  !Real values on radial grid
      real twexba(ntwexba)  !Real values on time grid
      real wexba (kxdim,*)  !Dependent variable wexb
      real xb    (nxb)      !The routine interpolates to this grid
      real time, tm(1)      !Interpolates to this time

c-----------------------------------------------------------------------
c Declaration of output variables
c-----------------------------------------------------------------------    

      real wexbxb (nxb)     !The interpolated result

c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      integer j, j1, j2     !Loop variables

      integer ikintrp, ikextrp, ikderiv, ikskip, ikxlast, iku
c
      parameter (ikintrp = 2) !Use Akimas method
      parameter (ikextrp = 0) !Extrapolate using flat profiles
      parameter (ikderiv = 0) !Calculate function values rather than der.
      parameter (ikskip  = 1) !send in vector rather than array to int1d
      parameter (iku     = 1) ! Only one time point is requested
      real zloc(1)
      real zwtmed(20)       !Intermediate vector for interpolation in t
      real zwxmed (55)      !Intermediate vector for interpolation in x
c
c-----------------------------------------------------------------------
c Start of actual code 
c-----------------------------------------------------------------------
c
c... set output initially to 0
c
 
         DO j1 = 1, nxb
           wexbxb(j1) = 0.0
         End do
c
c... check to se if we should exit with wexba(j) = 0.0  j = 1,...,nxb
c
      if (ntwexba .eq. 0 .or. nxwexba .eq.0 ) Return

c-----------------------------------------------------------------------
c Interpolate in time!
c-----------------------------------------------------------------------

      ikxlast = 0
      Do j1 = 1, nxwexba            
        DO j2 = 1, ntwexba
           zwtmed(j2) = wexba(j1,j2)
        END DO
	! ap{
	tm(1) = time
	!}
        call int1d (ikintrp, ikextrp, ikderiv, ntwexba, ikskip,
     &              twexba, zwtmed, iku, tm, ikxlast, zloc)
        ! ap {
        zwxmed(j1) = zloc(1)
	!}
      END DO

c-----------------------------------------------------------------------
c Interpolate to requested radial grid
c-----------------------------------------------------------------------

      ikxlast = 0
      call int1d (ikintrp, ikextrp, ikderiv, nxwexba, ikskip,
     &            xwexba, zwxmed, nxb, xb, ikxlast, wexbxb)

c-----------------------------------------------------------------------
c Return
c-----------------------------------------------------------------------
      END
















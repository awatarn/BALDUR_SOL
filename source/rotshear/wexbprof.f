      subroutine wexbprof(nx,rho,ni,ti,bpol,btor,vtor,rminor,rmajor
     &    ,wexba)
c-----------------------------------------------------------------------
      integer nx
      INTEGER, PARAMETER :: idp = KIND(1.0D0)
c-----------------------------------------------------------------------
      real(kind=idp)    rho   (nx)               
      real(kind=idp)    ni    (nx)               
      real(kind=idp)    Ti    (nx)               
      real(kind=idp)    Bpol  (nx)               
      real(kind=idp)    Btor  (nx)               
      real(kind=idp)    vtor  (nx)               
      real(kind=idp)    rminor(nx)               
      real(kind=idp)    rmajor(nx)               
      real(kind=idp)    wexba (nx)        
c-----------------------------------------------------------------------     
      real(kind=idp), parameter  :: e = 1.60219e-19   !unit charge

      integer :: intrp  = 2   !Use akima splines
      integer :: kextrp = 1   !Flat extrapolaation
      integer :: kderiv = 1   !Caclulate derivatives
      integer :: kskip  = 1   !Input data is 1-d
      integer :: kxlast = 0   ! Extrapolation point unknown
 
c-----------------------------------------------------------------------

      real(kind=idp), dimension (:), allocatable ::x      
      real(kind=idp), dimension (:), allocatable ::dTi      
      real(kind=idp), dimension (:), allocatable ::dpi      
      real(kind=idp), dimension (:), allocatable ::vpol
      real(kind=idp), dimension (:), allocatable ::Erad
      real(kind=idp), dimension (:), allocatable ::Epsr
      real(kind=idp), dimension (:), allocatable ::zx

c-----------------------------------------------------------------------
c
c... Allocate some local variables
c
 
      allocate (x(nx))
      allocate (dti(nx))
      allocate (dpi(nx))
      allocate (vpol(nx))
      allocate (epsr(nx))
      allocate (Erad(nx))
      allocate (zx(nx))
 
      DO I = 1,nx
           x(i) = rho(i)
         dpi(i) = ni(i) * Ti(i) *e
         dti(i) = Ti(i)
         epsr(i) = rminor(i)/(rmajor(i)+rminor(i))
         write(6,*) ni(i),ti(i),rho(i), epsr(i)
      end do

      call int1d(intrp,kextrp,kderiv,nx,kskip,x,dpi,nx,x,kxlast,dpi)
      call int1d(intrp,kextrp,kderiv,nx,kskip,x,dti,nx,x,kxlast,dTi)

      DO i = 2, nx
         fc = 1.0-1.46*sqrt(epsr(i))+0.46*epsr(i)*sqrt(epsr(i))
         k1 = 0.8839*fc/(0.3477+0.4058*fc)
         Vpol(i) = K1*dti(i)*Btor(i)/(BTOR(i)**2+ Bpol(i)**2)
         Erad(i) = dpi(i)/e/ni(i)-vpol(i)*Btor(i)+vtor(i)*bpol(i)
         zx  (i) = Erad(i)/(rmajor (i) + rminor(i)) /Bpol(i)
      end do
      zx(1) = zx(2) 
      call int1d(intrp,kextrp,kderiv,nx,kskip,x,zx,nx,x,kxlast,zx)

      DO i = 1,nx
          wexba(i) = ((Rmajor(i)+rminor(i))*bpol(i))**2/
     &      sqrt(Btor(i)**2 + Bpol(i)**2) *zx(i)
      END DO
      END
  
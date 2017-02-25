      subroutine wexbprof(nx,rho,ni,ti,bpol,btor,vtor,rminor,rmajor
     &    ,wexba)
c-----------------------------------------------------------------------
      integer nx
c-----------------------------------------------------------------------
      real    rho   (nx)               
      real    ni    (nx)               
      real    Ti    (nx)               
      real    Bpol  (nx)               
      real    Btor  (nx)               
      real    vtor  (nx)               
      real    rminor(nx)               
      real    rmajor(nx)               
      real    wexba (nx)        
c-----------------------------------------------------------------------     
      real e 
      parameter( e = 1.60219e-19)   !unit charge
      
      integer  intrp  
      integer  kextrp 
      integer  kderiv 
      integer  kskip  
      integer  kxlast 
 
c-----------------------------------------------------------------------
      include 'dim_readvec'
      real  x(kxdim)      
      real dTi(kxdim)      
      real dpi(kxdim)      
      real vpol(kxdim)
      real Erad(kxdim)
      real Epsr(kxdim)
      real zx(kxdim)

c-----------------------------------------------------------------------
c
c... Allocate some local variables
c
      
        intrp  = 2   !Use akima splines
        kextrp = 1   !Flat extrapolaation
        kderiv = 1   !Caclulate derivatives
        kskip  = 1   !Input data is 1-d
        kxlast = 0   ! Extrapolation point unknown      allocate (x(nx))
 
      DO I = 1,nx
           x(i) = rho(i)
         dpi(i) = ni(i) * Ti(i) *e
         dti(i) = Ti(i)
         epsr(i) = rminor(i)/(rmajor(i)+rminor(i))
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
  
c
c This subroutine is used to calculate the temperature at the top of the 
c pedestal used the model based on Thermal Conduction concept. For more 
c details, see section 3.4 in pedestal_models.ps.
c 
c@bdheight11.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bdheight11 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      zploss,     znebar,
     & neped,      tped)
c
c  Input:
c  ------
c lbound(j)   integer control array
c cbound(j)   general control array
c zrminor     minor radius [m]
c zrmajor     major radius [m]
c zkappa      elongation
c zdelta      triangularity
c zcurrent    plasma current [MA]
c zbtor       toroidal field at rmajor [Tesla]
c zhydmass    average hydrogenic mass [AMU]
c zzeff       effective charge at the edge of plasma
c zploss      power across separatrix [MW]
c znebar      line averaged electron density [particles/m^3]
c neped       density at the top of the pedestal [particles/m^3]
c     
c  Output:
c  -------
c tped      temperature at the top of the pedestal [keV]
c                                                               
      implicit none                                             
c 
      integer    
     &  lbound(32)
c
      real       
     &  cbound(32),  zrminor,  zrmajor,  zkappa,
     &  zdelta,      zcurrent, zbtor,    zhydmass, 
     &  zzeff,       zploss,   znebar,   neped,
     &  pi,          kb,       esp,      kappa95,
     &  delta95,     q95,      qcyl,     Fq,
     &  wped,        vol,      pped,     tped     
c
c..constant
c
      pi = atan2 ( 0.0, -1.0 )
      kb = 1.6022E-16
c
c..calculate inverse aspect ratio
c
      esp  = zrminor/zrmajor
c
c..calculate safety factor at 95% flux surface
c
         kappa95 = 0.914*zkappa
         delta95 = 0.85*zdelta
         q95 = ((5.*(zrminor**2)*zbtor)/(2.*zcurrent*zrmajor)           
     &        *(1+(kappa95**2)*(1+2.*(delta95**2)-1.2*(delta95**3))) 
     &        *(1.17-0.65*esp)/((1-(esp**2))**2))
c
c..calculate cylindrical safety factor
c
      qcyl = 5.*zkappa*(zrminor**2)*zbtor/(zrmajor*zcurrent)
c
c..calculate shaping factor
c
      Fq = q95/qcyl
c
c..calculate the stored energy [J] in the pedestal regime
c
      wped = (1.E+6)*0.00064*(zcurrent**1.58)*(zrmajor**1.08)
     &      *(zploss**0.42)*((znebar*1.E-19)**(-0.08))*(zbtor**0.06)
     &      *(zkappa**1.81)*(esp**(-2.13))*(zhydmass**0.20)
     &      *(Fq**2.09)
c
c..calculate plasma volume [m^3]
c
      vol = 2.*(pi**2)*(zrminor**2)*zrmajor*zkappa
c
c..calculate total pedestal pressure (pe_ped + pi_ped = wped/0.92 vol)
c  Note that constant of 0.92 is the fraction of the total volume 
c  occupied by tge pedestal 
c
      pped = wped/(0.92*vol)
c
c..calcalate pedestal tempeature [keV] (pped = 3nktped)
c
      tped = pped/(3.*neped*kb)
c
      return
      end















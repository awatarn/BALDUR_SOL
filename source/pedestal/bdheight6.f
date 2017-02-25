c
c This subroutine is used to calculate the temperature at the top of the 
c pedestal used the model that pressure gradient is limited by the first 
c ballooning limit and pedestal width is based on normalized pressure. 
c For more details, see section 3.3 in pedestal_models.ps.
c 
c@bdheight6.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bdheight6 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      zshear,
     & zq,         zbpol,      tped,       dpdr_c, 
     & width,      iter_t)
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
c neped       density at the top of the pedestal [particles/m^3]
c     
c  Output:
c  -------
c tped      temperature at the top of the pedestal [keV]
c dpdr_c    critical pressure gradient [Pa/m]
c width     width of the pedestal [m]
c iter_t    number of iteration used in the temperature calculation                                                              
c                                                               
      implicit none                                             
c 
      integer    
     &  lbound(32),  num,      iter_t,   iter_t_max,
     &  ierr(2)
c
      real       
     &  cbound(32),  zrminor,  zrmajor,  zkappa,
     &  zdelta,      zcurrent, zbtor,    zhydmass, 
     &  zzeff,       neped,    pi,       mu,
     &  kb,          esp,      kappa95,  delta95,
     &  t_guess,     zq,       zshear,   width,
     &  dpdr_c,      pr_c,     t_pre,    tped,
     &  alpha_c,     q95,      zbpol,    diff_t_max,
     &  Cw,          beta_p      
c
c..define constant
c
      pi = atan2 ( 0.0, -1.0 )
      mu = 4*pi*1.E-7
      kb = 1.6022E-16
c
c..set control for iteration
c
      if (lbound(7) .gt. 0) then
         iter_t_max = lbound(7)
      else
         iter_t_max = 1000
      endif
c
      if (cbound(7) .gt. 0.) then
         diff_t_max = cbound(7)
      else
         diff_t_max = 1.E-2
      endif   
c
c..calculate inverse aspect ratio
c
      esp  = zrminor/zrmajor
c
c..approximate elongation, triangularity and safety factor at 95% flux surface
c
c      kappa95 = 0.914*zkappa
c      delta95 = 0.85*zdelta
      kappa95 = zkappa
      delta95 = zdelta
c
c..Start iteration for calculating the height of pedestal
c                
      t_guess = 1.0
      iter_t  = 0
      Cw      = cbound(5)       
c
      do 
c
c..calculate poloidal beta
c
      beta_p = (8.0384E-22)*neped*t_guess/(zbpol**2)
c
c..calculate the width of the pedestal and position of the top of the pedestal
c
      width = Cw*zrmajor*sqrt(beta_p)
c
c..calculate the maximum normalized pressure gradient from
c  ballooning instability
c
      alpha_c = 0.4*zshear*(1+(kappa95**2)*(1+5.0*(delta95**2)))
c
c..calculate critical pressure gradient of the first ballooning mode
c
      dpdr_c = (zbtor**2)*alpha_c/(2*mu*zrmajor*(zq**2))
c
c..calculate critical pedestal pressure by assuming pressure gradient 
c  within the pedestal is constant
c
      pr_c = dpdr_c*width
c
c..calcalate the temperature at the top of pedestal using pedestal width 
c  based on normalized poliodal pressure
c
      t_pre = pr_c/(2.*kb*neped) 
c
c..iterate to solve non linear equation of pedestal height
c
      if(abs(t_pre - t_guess) .le. diff_t_max) then
         tped = t_pre
         return
      endif
c
      iter_t = iter_t+1
c
c..calculate T_guess for next step
c
      if (iter_t .le. iter_t_max) then
         t_guess = (t_pre+t_guess)/2
c
      else
         write(6,*) 'Iteration on temperature exceeds', iter_t_max,
     &              ' steps.'
         ierr(2) = 1
         tped    = t_pre   
         return
      endif
c
      enddo
c
      return
      end















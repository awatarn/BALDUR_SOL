c
c This subroutine is used to calculate the magnetich shear at the top 
c of the pedestal. The magnetic shear is a function of the pedestal 
c width, which is based on narmalized poloidal pressure.
c
c@bdshear3.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine bdshear3 (
     & lbound,    cbound,     zrminor,    zrmajor,
     & zbtor,     zhydmass,   zzeff,      neped,
     & kappa95,   delta95,    q95,        t_guess,
     & q,         shear,      width)
c
c  Input:
c  ------
c lbound(j)   integer control array
c cbound(j)   general control array
c zrminor     minor radius [m]
c zrmajor     major radius [m]
c zbtor       toroidal field at rmajor [Tesla]
c zhydmass    average hydrogenic mass [AMU]
c zzeff       effective charge at the edge of plasma
c neped       density at the top of the pedestal [particles/m^3]
c kappa95     elongation at the 95% flux surface
c delta95     triangularity at the 95% flux surface
c q95         safety factor at the 95% flux surface
c t_guess     guessing pedestal temperature [keV]
c     
c  Output:
c  -------
c q           safety factor at the top of the pedestal
c shear       magnetic shear at the top of the pedestal
c width       pedestal width [m]
c
      IMPLICIT NONE
c
      integer    
     &  lbound(32)
c
      real       
     &  cbound(32),  zrminor,      zrmajor,  zbtor,
     &  zhydmass,    zzeff,        neped,    kappa95,
     &  delta95,     q95,          t_guess,  esp,
     &  rho,         sh,           width,    x,
     &  q,           nuhat_e,      nuhat_i,  c_b,
     &  shear,       b_pol,        beta_p,   Cw
c
c..calculate aspect ratio
c
      esp = zrminor/zrmajor
c
c..calculate shaping factor
c
      sh  = (1+(kappa95**2)*(1+2*(delta95**2)-1.2*(delta95**3))) 
     &      *(1.17-0.65*esp)/(2*(1-(esp**2))**2)
c
c..calculate poloidal magnetic field
c
      b_pol = ((5*sh)/(3.14*(1+kappa95)))*zrminor*zbtor/
     &        (zrmajor*q95)
c
c..calculate poloidal beta
c
      beta_p = (8.0384E-22)*neped*t_guess/(b_pol**2)
c
c..calculate the width of the pedestal and position of the top of the pedestal
c
      Cw    = cbound(5)
      width = Cw*zrmajor*sqrt(beta_p)
      x     = 1 - (width/zrminor)
c
c..calculate safety factor at the top of the pedestal
c
      q=(q95/2.93)*(((1+((x/1.4)**2))**2) + (0.267*abs(LOG(1-x))))
c
c..calculate collisionality at the top of the pedestal
c	
      call bdcoll(
     &   zrminor,   zrmajor,   zhydmass,   zzeff,
     &   neped,     q,         t_guess,    nuhat_e,
     &   nuhat_i)
c
c..calculate collisional effect in boostrap current
c
      call bdboots (
     &   zrminor,   zrmajor,   nuhat_i, nuhat_e,
     &   c_b)
c
c..calculate magnetic shear
c
      shear = x*((4*x/(1.4*1.4))*(1+((x/1.4)**2)) 
     &        + 0.267*abs(1/(1-x)))
     &        /(((1+((x/1.4)**2))**2) + (0.267*abs(LOG(1-x))))
c
c..calculate magnetic shear modified by bootstrap current
c
      shear = shear/(1+((0.1*c_b*shear
     &        *(1+(kappa95**2)*(1+5.0*(delta95**2))))
     &        /sqrt(esp)))
c
      return
      end
c

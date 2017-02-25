c
c This subroutine is used to calculate the magnetich shear at the top 
c of the pedestal. The magnetic shear is a function of the pedestal 
c width, which is based on magnetic and flow shear stabilization
c
c@bdshear1.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine bdshear1 (
     & lbound,    cbound,     zrminor,    zrmajor,
     & zbtor,     zhydmass,   zzeff,      neped,
     & kappa95,   delta95,    q95,        t_guess,
     & q,         shear,      width,      iter_s,
     & ierr)
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
c iter_s      number of iterations used in the magnetic shear calculation                                                              
c
      IMPLICIT NONE
c
      integer    
     &  lbound(32),  iter_s,       iter_s_max, ierr(2)
c
      real       
     &  cbound(32),  zrminor,      zrmajor,    zbtor,
     &  zhydmass,    zzeff,        neped,      kappa95,
     &  delta95,     q95,          t_guess,    esp,
     &  rho,         shear_guess,  width,      x,
     &  q,           nuhat_e,      nuhat_i,    c_b,
     &  shear,       diff_s_max,   Cw
c
c..set the maximum number of iteration in shear calculation
c
      if (lbound(8) .gt. 0) then
         iter_s_max = lbound(8)
      else
         iter_s_max = 1000
      endif
c
c..set the accuracy in shear calculation
c
      if (cbound(8) .gt. 0.) then
         diff_s_max = cbound(8)
      else
         diff_s_max = 1.E-2
      endif   
c
c..calculate aspect ratio
c
      esp = zrminor/zrmajor
c
c..calculate ion gyro radius
c
      rho = (4.57E-3)*sqrt(zhydmass*t_guess)/(zbtor)
c
c..iterate shear
c
      shear_guess = 2.0
      iter_s      = 0
      Cw          = cbound(5)
c
      do
c
c..calculate the width of the pedestal and position of the top of the pedestal
c
         width = Cw*rho*(shear_guess**2)
         x     = 1 - (width/zrminor)
c
c..calculate safety factor at the top of the pedestal
c
         q=(q95/(((1+((0.95/1.4)**2))**2)+(0.267*abs(LOG(1-0.95)))))
     &      *(((1+((x/1.4)**2))**2) + (0.267*abs(LOG(1-x))))
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
c..calculate magnetic shear from s=(r/q)(dq/dr)
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
c..iterate to solve non linear equation of magnetic shear
c
         if (abs(shear - shear_guess) .le. diff_s_max) then
            return
         endif
c
c..calculate shear_guess for next step
c
         if (iter_s .le. iter_s_max) then
            shear_guess = 0.5*shear_guess + 0.5*shear
            iter_s = iter_s+1
c
         else 
            write(6,*) 'Iteration on magnetic shear exceeds', 
     &                 iter_s_max,' steps.'
            ierr(2) = 12
            shear_guess = 0.5*shear_guess + 0.5*shear
            return
         endif   
c            
      enddo
c
      return
      end
c

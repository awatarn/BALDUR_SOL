c
c This subroutine is used to calculate the collisionality at the top 
c of the pedestal.
c
c@bdcoll.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
        subroutine bdcoll (
     &   zrminor,   zrmajor,   zhydmass,   zzeff,
     &   neped,     q,         t_guess,    nuhat_e,
     &   nuhat_i)
c
c  Input:
c  ------
c zrminor     minor radius [m]
c zrmajor     major radius [m]
c zhydmass    average hydrogenic mass [AMU]
c zzeff       effective charge at the edge of plasma
c neped       density at the top of the pedestal [particles/m^3]
c q           safety factor at the top of the pedestal
c t_guess     guessing pedestal temperature [keV]
c     
c  Output:
c  -------
c nuhat_e     electron normalized collisionality 
c nuhat_i     ion normalized collisionality
c  
      implicit none
c
      real 
     &   zrminor,   zrmajor,   zhydmass,  zzeff,
     &   neped,     q,         t_guess,   esp,
     &   ln_e,      ln_i,      nu_e,      nu_i,
     &   w_be,      w_bi,      nuhat_e,   nuhat_i    
c
c..calculate inverse aspect ratio
c
        esp  = zrminor/zrmajor
c
c..calculate value of ln(lambda) for collision
c
        ln_e = 15.2 - 0.5*LOG(neped*1E-20) + LOG(t_guess) 
        ln_i = 17.3 - 0.5*LOG(neped*1E-20) + 1.5*LOG(t_guess) 
c 
c..calculate value of collision frequency
c
        nu_e = (1.09E+16)*(t_guess**(3./2.))
        nu_e = nu_e/(neped*(zzeff**2)*ln_e) 
        nu_e = 1/nu_e
c
        nu_i = (6.6E+17)*sqrt(zhydmass)*(t_guess**(3./2.))
        nu_i = nu_i/(neped*(zzeff**4)*ln_i) 
        nu_i = 1/nu_i
c 
c..calculate value of trapped particle bounce frequency
c
        w_be = sqrt(esp)*((1.33E+7)*sqrt(t_guess))/(zrmajor*q)
        w_bi = sqrt(esp)*((3.09E+5)*sqrt(t_guess/zhydmass))/(zrmajor*q)
c 
c..calculate value of normalized collisionality
c
        nuhat_e = nu_e/(esp*w_be)
        nuhat_i = nu_i/(esp*w_bi)
c
      return
      end






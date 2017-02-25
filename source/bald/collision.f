c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
        subroutine collision (
     &   rminor,   rmajor,   dense,    hydmass,  
     &   zeff,     q,        t_ped,    znuhat_e,
     &   znuhat_i)
c
      implicit none
c
      real 
     &   rminor,   rmajor,   dense,    hydmass,  
     &   q,        zeff,     t_ped    
c
      real
     &   zep,      ln_e,     ln_i,     nu_e,     
     &   nu_i,     w_be,     w_bi,     znuhat_e, 
     &   znuhat_i    
c
c..calculate inverse aspect ratio
c
        zep  = rminor/rmajor
c
c..calculate value of ln(lambda) for collision
c
        ln_e = 15.2 - 0.5*LOG(dense*1E-20) + LOG(t_ped) 
        ln_i = 17.3 - 0.5*LOG(dense*1E-20) + 1.5*LOG(t_ped) 
c 
c..calculate value of collision frequency
c
        nu_e = (1.09E+16)*(t_ped**(3./2.))
        nu_e = nu_e/(dense*(zeff**2)*ln_e) 
        nu_e = 1/nu_e
c
        nu_i = (6.6E+17)*sqrt(hydmass)*(t_ped**(3./2.))
        nu_i = nu_i/(dense*(zeff**4)*ln_i) 
        nu_i = 1/nu_i
c 
c..calculate value of trapped particle bounce frequency
c
        w_be = sqrt(zep)*((1.33E+7)*sqrt(t_ped))/(rmajor*q)
        w_bi = sqrt(zep)*((3.09E+5)*sqrt(t_ped/hydmass))/(rmajor*q)
c 
c..calculate value of normalized collisionality
c
        znuhat_e = nu_e/(zep*w_be)
        znuhat_i = nu_i/(zep*w_bi)
c
      return
      end






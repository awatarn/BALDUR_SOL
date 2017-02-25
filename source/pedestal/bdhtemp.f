c
c This subroutine is used to call the routines that compute the 
c temperature at the top of the pedestal
c
c@bdhtemp.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bdhtemp (                                         
     &       lbound,       cbound,     zrminor,      zrmajor,
     &       zkappa,       zdelta,     zcurrent,     zbtor,
     &       neped,        zhydmass,   zzeff,        zshear,
     &       zq,           zbpol,      znebar,       zploss,
     &       tped)
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
c zzeff       effective charge at the top of the pedestal
c zshear      magnetic shear at the top of the pedestal
c zq          safety factor at the top of the pedestal
c neped       density at the top of the pedestal [particles/m^3]
c     
c  Output:
c  -------
c tped      temperature at the top of the pedestal [keV]
c dpdr_c    critical pressure gradient [Pa/m]
c width     width of the pedestal [m]
c iter_t    number of iterations used in the temperature calculation                                                              
c iter_s    number of iterations used in the magnetic shear calculation                       
c iter_q    number of iterations used in the safety factor calculation
c                                                              
      implicit none                                             
c 
      integer    
     &  lbound(32),  iter_t,   iter_s,   iter_q,
     &  ierr(2)
c
      real       
     &  cbound(32),  zrminor,  zrmajor,  zkappa,
     &  zdelta,      zcurrent, zbtor,    zhydmass, 
     &  zzeff,       zshear,   zq,       zbpol,
     &  neped,       tped,     dpdr_c,   width,
     &  znebar,      zploss
c
         if (lbound(5) .le. 1) then
c
c..call the pedestal temperature model using the pedestal width based on 
c  magnetic and flow shear stabilization. For more details, see section 3.1 
c  in pedestal_models.ps.
c
         if ((cbound(5) .lt. 0.01) .or. (cbound(5) .gt. 5.00)) then
            write(6,*)'Cw is out of range. Use the default value.'
            cbound(5) = 2.42
         endif
c
         call bdheight1 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      tped,
     & dpdr_c,     width,      iter_t,     iter_s,
     & ierr)
c
         elseif (lbound(5) .eq. 2) then
c
c..call the pedestal temperature model using the pedestal widhth based on 
c  flow shear stabilization. For more details, see section 3.2 in 
c  pedestal_models.ps.
c
         if ((cbound(5) .lt. 0.01) .or. (cbound(5) .gt. 1.00)) then
            write(6,*)'Cw is out of range. Use the default value.'
            cbound(5) = 0.22
         endif
c
         call bdheight2 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      tped,
     & dpdr_c,     width,      iter_t,     iter_q,
     & ierr)
c
         elseif (lbound(5) .eq. 3) then
c
c..call the pedestal temperature model using the pedestal width based on 
c  normalized pressure. For more details, see section 3.3 in 
c  pedestal_models.ps.
c
         if ((cbound(5) .lt. 0.001) .or. (cbound(5) .gt. 0.1)) then
            write(6,*)'Cw is out of range. Use the default value.'
            cbound(5) = 0.021
         endif
c
         call bdheight3 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      tped,
     & dpdr_c,     width,      iter_t,     ierr)
c
         elseif (lbound(5) .le. 4) then
c
c..call the pedestal temperature model using the pedestal width based on
c  magnetic and flow shear stabilization.
c
c         if ((cbound(5) .lt. 0.01) .or. (cbound(5) .gt. 5.00)) then
c            write(6,*)'Cw is out of range. Use the default value.'
c            cbound(5) = 2.42
c         endif
c
         call bdheight4 (
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      zshear,
     & zq,         tped,       dpdr_c,     width,
     & iter_t)
c
         elseif (lbound(5) .eq. 5) then
c
c..call the pedestal temperature model using the pedestal widhth based on
c  flow shear stabilization. For more details, see section 3.2 in
c  pedestal_models.ps.
c
c         if ((cbound(5) .lt. 0.01) .or. (cbound(5) .gt. 1.00)) then
c            write(6,*)'Cw is out of range. Use the default value.'
c            cbound(5) = 0.22
c         endif
c
         call bdheight5 (
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      zshear,
     & zq,         tped,       dpdr_c,     width,
     & iter_t)
c
         elseif (lbound(5) .eq. 6) then
c
c..call the pedestal temperature model using the pedestal width based on
c  normalized pressure. For more details, see section 3.3 in
c  pedestal_models.ps.
c
c         if ((cbound(5) .lt. 0.001) .or. (cbound(5) .gt. 0.1)) then
c            write(6,*)'Cw is out of range. Use the default value.'
c            cbound(5) = 0.021
c         endif
c
         call bdheight6 (
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      neped,      zshear,
     & zq,         zbpol,      tped,       dpdr_c,     
     & width,      iter_t)
c
         elseif (lbound(5) .eq. 11) then
c
c..call the pedestal temperature model based on thermal 
c  conduction model I. For more details, see section 3.4 in 
c  pedestal_models.ps.
c
         call bdheight11 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      zploss,     znebar,
     & neped,      tped)
c
         elseif (lbound(5) .eq. 12) then
c
c..call the pedestal temperature model based on thermal 
c  conduction model II. For more details, see section 3.5 
c  in pedestal_models.ps.
c
         call bdheight12 (                                         
     & lbound,     cbound,     zrminor,    zrmajor,
     & zkappa,     zdelta,     zcurrent,   zbtor,
     & zhydmass,   zzeff,      zploss,     znebar,
     & neped,      tped)
c
         endif
c
      return
      end















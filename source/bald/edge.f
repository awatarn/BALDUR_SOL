c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine edge (                                         
     &   rminor,   rmajor,   kappa,    delta,                   
     &   current,  btor,     dense,    hydmass,                 
     &   zeff,     t_ped,    ierr,     iwarn )                  
c                                                               
      implicit none                                             
c                                                               
c..declare and set constant variables                           
c                                                               
      real
     &   delta_max, ratio_min, iter_max, t_guess_int,
     &   diff_min
c
      parameter ( delta_max   = 2.23293,      
     &            ratio_min   = 0.2,        
     &            iter_max    = 100,               
     &            t_guess_int = 1.0,
     &            diff_min    = 1.0E-5)                      
c                                                               
c..declare variable for input data                              
c                                                               
      real                                                      
     &   rminor,   rmajor,   kappa,    delta,                      
     &   current,  btor,     dense,    hydmass,                 
     &   zeff                                                   
c                                                               
c..declare variable                                             
c                                                               
      integer                                                
     &   nout,     ierr,     iwarn,    j,                       
     &   iter                                                   
c                                                               
      real                                                 
     &   zepsilon, zpi,      zrminor,  zrmajor,                 
     &   zkappa,   zdelta,   zcurrent, zbtor,                   
     &   zdense,   zhydmass, zzeff,    zep,                     
     &   kappa95,  delta95,  q,        n_gr,                    
     &   ratio,    t_guess,  t_pre,    t_ped
c
c..Inputs:
c
c    rminor    :   Minor radius  [m]                                
c    rmajor    :   Major radius  [m]                                 
c    kappa     :   Elongation at the separatrix                      
c    delta     :   Triangularity at the separatrix                  
c    current   :   Plasma current [MA]                               
c    btor      :   Magnetic field [T]
c    dense     :   Pedestal density [particle/m^3]
c    hydmass   :   Average hydrogenic mass  [amu]
c    zeff      :   Effective charge
c
      nout = 9
c
      zepsilon = 1.e-10	
c
c..initial output variables
c
      t_ped = 0.0
      ierr  = 0
      iwarn = 0
c
c..physical constants
c
      zpi = atan2 ( 0.0, -1.0 )
c
c..change name of input variable
c
      zrminor  = rminor  
      zrmajor  = rmajor
      zkappa   = kappa     
      zdelta   = delta
      zcurrent = current
      zbtor    = btor
      zdense   = dense
      zhydmass = hydmass
      zzeff    = zeff        
c
c..check input for validity
c
      if (zrminor .lt. zepsilon) then
         ierr = 1
         return
      elseif (zrmajor .lt. zepsilon) then
         ierr = 2
         return
      elseif (zkappa .lt. zepsilon) then
         ierr = 3
         return
      elseif (zcurrent .lt. zepsilon) then
         ierr = 4
         return
      elseif (zbtor .lt. zepsilon) then
         ierr = 5
         return
      elseif (zdense .lt. zepsilon) then
         ierr = 6
         return
      elseif ((zhydmass .lt. zepsilon) .or. (zhydmass .gt. 3.0)) then
         ierr = 7
         return
      elseif (zzeff .lt. zepsilon) then
         ierr = 8
         return
      endif
c
c..If delta is higher than 2.23293, the prediction will be negative. 
c  The prediction will be set to be at delta = 2.2329 and give warning
c
      if (zdelta .ge. delta_max) then   
         zdelta = 2.2329
         iwarn = 1
      endif
c
c..calculate Greenwald density
c
      n_gr = (current/(zpi*rminor**2))*1.0E+20
c
c..check the validation of model prediction capability. 
c  If ratio of pedestal density to the Greenwald density is lower than
c  0.20, the ratio will set to be 0.20 and give the warning.
c
      ratio = zdense/n_gr
      if (ratio .lt. ratio_min) then
         zdense = ratio_min*n_gr
         iwarn = 2
      endif
c
c..calculate inverse aspect ratio
c
      zep  = rminor/rmajor
c
c..calculate elongation and triangularity at 95 0.000000lux surface
c
      kappa95 = 0.914*zkappa
      delta95 = 0.85*zdelta
c
c..calculate the safety factor
c
      q = (5*(zrminor**2)*zbtor)/(2*zcurrent*zrmajor)           
     &      *(1+(kappa95**2)*(1+2*(delta95**2)-1.2*(delta95**3))) 
     &      *(1.17-0.65*zep)/((1-(zep**2))**2)
c
c..calculate the height of pedestal
c                
      t_pre = t_guess_int 
      t_guess = t_guess_int
      iter = 1      
c
      do 
c
      call height1 (                                             
     &   zrminor,  zrmajor,  kappa95,  delta95,                
     &   zcurrent, zbtor,    zdense,   zhydmass,             
     &   zzeff,    q,        t_guess,  t_pre)
c
c..iterate to solve non linear equation of pedestal height
c
      if(abs(t_pre - t_guess) .le. diff_min) then
         t_ped = t_pre
         return
      endif
c
      t_guess = t_pre
c
      if(iter .ge. iter_max) then
         iwarn = 3
         return
      endif
c
      enddo
c
      return
      end















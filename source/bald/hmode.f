c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine hmode ( ztime, T_ped )
c
          include 'cparm.m'
          include 'cbaldr.m'
          include 'commhd.m'
          include 'cd3he.m'
          include 'clintf.m'
          include 'clsaw.m'
          include 'clparm.m'
          include 'clislb.m'
c 
          integer    ierr,        iwarn
c
          real       nepedestal,  n_gr,       pi, 
     &   current,    hydmass,    nhpedestal,
     &   r_maj,      r_min,       kappa,      delta, 
     &   b_tor,      nzpedestal, zeff,
     &   tekev,      tikev,       T_ped,      zeheat,
     &   ziheat,     P_heat,      zne,        znebar,
     &   zs,         P_th,        ratio,      ztime
c
         r_min   = max (rminb, 0.01) 
         r_maj   = max (rmajb, 0.01)
         kappa   = elongb   
         delta   = tringb
         current = max (eqcamp * 1.E-6, 1.0E-6)
         b_tor   = max (bzs * 1.0E-4, 1.0E-4)  
         nepedestal = rhoels(1,mzones) * 1.E+6
         hydmass    = ahmean(1,mzones)
         zeff     = xzeff(1,mzones)
c
        call edge (
     &   r_min,    r_maj,    kappa,         delta,
     &   current,  b_tor,    nepedestal,    hydmass,  
     &   zeff,     T_ped,    ierr,          iwarn)
c
	      write(nout,*)	
	      write(nout,*)' At time t = ', ztime, ' sec'
      	      write(nout,*)' The Height of Pedestal is ',T_ped, ' keV.'
c
      return
      end

c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine lhmode ( ztime, mode )
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
          integer    :: ierr,        iwarn,      mode
	  integer    :: ij(1)
c
          real       nepedestal,  n_gr,       pi,
     &   current,    hydmass,     nhpedestal,
     &   r_maj,      r_min,       kappa,      delta,
     &   b_tor,      nzpedestal, zeff,
     &   tekev,      tikev,       T_ped,      zeheat,
     &   ziheat,     P_heat,      zne,        znebar,
     &   zs,         P_th,        ratio,      ztime
c
         pi      = atan2 ( 0.0, -1.0 )
c
c..initail mode to be L-mode
c
c  mode = 0  for L-mode
c       = 1  for H-mode
c
         mode = 0
c
         iaxis = lcentr
         iedge = mzones
         isep = mzones
         if (nadump(1) .gt. lcentr) isep = nadump(1)
c
         current = max (eqcamp * 1.E-6, 0.01)
         r_min   = max (rminb, 0.01) 
         b_tor   = max (bzs * 1.0E-4, 0.01)  
c
c..calculate the Greenwald limit
c
         n_gr = (current / (pi * r_min * r_min))*(1.E+20)
c
         nepedestal  = rhoels(1,mzones) * 1.E+6
c 
c..calculate the total power apllied
c
        zeheat = usip * ( geauxs(isep) + gebems(isep) +
     &                   geecrs(isep) + geicrs(isep) )
c
        ziheat = usip * ( giauxs(isep) + gibems(isep) +
     &                   giecrs(isep) + giicrs(isep) )
c
        P_heat = (zeheat + ziheat) * 1.E-6 
c
c..calculate average electron density
c
          zne = 0.0
c
        do 204 jz = lcentr, mzones
          zne = zne + rhoels(1,jz) * 1.E+6
  204   continue
c
        znebar = zne / (mzones - lcentr + 1)
        if (lthery(47)==1) then
	  mode = lhte()
          !
          ij  = maxloc(slnes(mj-15:mj))
          !
          te_crit = 0.45*(b_tor*b_tor*xzeff(2,ij(1))/
     &              sqrt(rmajb*ngas(1)))**(1./3.)*sqrt(slnes(ij(1)))
	  if (te_crit>te(ij(1))) mode = 1
	  return
          !
	endif
c
c..calculate threshold power
c
	zs       = avi(iedge,4,1) * avi(iedge,5,1)            
	P_th     = 0.04 * (znebar * 1.E-20) * B_tor * zs
c
        ratio    = nepedestal / n_gr
c
	if ((P_heat .gt. P_th)) then
           mode = 1
        endif
c
      return
      end

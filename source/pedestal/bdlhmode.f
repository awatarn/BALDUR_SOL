c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine bdlhmode (
     &       lbound,     cbound,     P_heat,     rmajor,     
     &       rminor,     b_tor,      hydmass,    znebar,
     &       P_th,       mode)
c
      IMPLICIT NONE
c 
      integer    lbound(32), mode    
c
      real       cbound(32), P_heat,   rmajor,     rminor,
     &           b_tor,      hydmass,  znebar,     P_th 
c
c Model for L-H transition is obtained from Y. Shimomura, et. al.
c This empirical mode was presented at the IAEA 2000
c
c  mode = 0  for L-mode
c       = 1  for H-mode 
c
c..use model of power threshold from an empirical formula by Y. Shimomura
c  IAEA 2000
c
      if (lbound(1) .le. 1) then
c
c        write(6,*) hydmass,' ',b_tor,' ',znebar,' ',rmajor,' ',rminor, 
c     &  ' ',mode
c
	P_th     = 2.84 * (1/hydmass) * (b_tor**0.82) 
     &             * ((znebar * 1.E-20)**(0.58))
     &             * (rmajor) * (rminor**0.81)
c
      endif  
c
c..use cbound(1) is the fraction of power from H-L
c
        if (mode .eq. 1) then 
           P_th = cbound(1) * P_th
        endif   
c
c..check status of plasma
c  if P_heat > P_th, plasma is in H-mode (mode = 1)
c
	if ((P_heat .ge. P_th)) then
           mode = 1
        else
           mode = 0
        endif        
!
        write(6,'(1A, i2, 1A, e13.5, 1A, e13.5, 1A)')
     &      'LH-Mode = ', mode, '  (P_cr-P_heat=', P_th, ' -', P_heat,
     &       ')'
c
      return
      end

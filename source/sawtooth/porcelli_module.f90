!| \documentstyle{article}
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \title{Porcelli model}
!| \author{Canh N. Nguyen and Glenn Bateman}
!| 
!| \begin{document}           % End of preamble and beginning of text.
!| 
!| \maketitle                 % Produces the title.
!| 
!| The model for testing whether a sawtooth occurs presented below is
!| by Porcelli \cite{porcelli96a}.  The conditionss for a sawtooth crash
!| are given by equation
!| (13)-(15) in the reference.  Namely, sawteeth crashes are triggered when
!| one of the following conditions is met:\\
!| \indent $-\delta \hat{W}_{\rm core} > c_{\rm h}\omega_{Dh}\tau_{\rm A}$\\
!| \indent $-\delta \hat{W}  > 0.5\omega_{\ast \rm i}\tau_{\rm A}$\\
!| \indent $-c_{\rho}\hat{\rho} < -\delta \hat{W}< 0.5\omega_{\ast \rm i}\tau_{\rm A}$
!| and $\omega_{\ast \rm i} < c_{\ast} \gamma_{\rho}$\\
!| where $\delta\hat{W}$ is the nomalized potential energy functional and
!| $\tau_{\rm A}$ is the Alfven time.
!| The numerical factor $c_{\rm h}$, $c_{\rm\rho}$, and $c_{\ast}$ are of
!| order unity (their precise evaluation requires knowledge of the plasma
!| equilibrium).
!| From the Appendix of the \cite{porcelli96a}:\\
!| \indent $\delta\hat{W}=\delta\hat{W}_{\rm Busac}
!|   +\delta\hat{W}_{\rm el}+\delta\hat{W}_{\rm KO}
!| +\delta\hat{W}_{\rm fast}$\\
!| With this definition of $\delta\hat{W}$, the internal kink growth rate (disregarding
!| neoclassical effects) is $\gamma=\delta\hat{W}\tau_{\rm A}^{-1}$ when layer physics
!| effects are not important.
!| 
!| \begin{enumerate}
!| 
!| \item The Alfven time is: $\tau_{\rm A}=\sqrt{3}R/v_{\rm A}$\\
!| $\tau_{\rm A}(\rm sec)=0.8\frac{R(\rm m)}{B(\rm T)}\sqrt{\mu n_{0}(10^{20} m^{-3})}$\\
!| where $\mu$ is the ratio of fast ion mass and mass of proton.
!| 
!| \item Fast particle precession frequency:
!| $\omega_{\rm Dh}=cE_{\rm f}/4eBR\bar{r}_{1}$\\
!| $\omega_{\rm Dh}(\rm sec^{-1})=250.0*E_{\rm f}(\rm keV)/z_{\rm i}B(\rm T)R(\rm m)
!| \bar{r}_{1}(\rm m)$\\
!| where $z_{\rm i}$ is the charge of the fast ion.
!| 
!| \item Ion diamagnetic frequency: $\omega_{\ast i}=cT/eB\bar{r}_{1}r_{\rm p}$\\
!| $\omega_{\ast \rm i}(\rm sec^{-1})
!|   = 1000T_{\rm p0}(keV)/B(\rm T)\bar{r}_{1}(m)r_{\rm p}(m)$\\
!| where $r_{\rm p}=|dp/dr|^{-1}p_{\rm i}$ is the pressure scale length.
!| 
!| \item Normalized ion larmor radius $\hat{\rho}=\rho_{\rm i}/\bar{r}_{1}$\\
!| $\hat{\rho}= 0.00323*\sqrt{T_{\rm p0}(keV)}/B(T)\bar{r}_{1}(m)$.
!| 
!| \item Semi-collisional m = 1 growth rate:
!|   $\gamma_{\rho}=1.1\hat{\rho}^{4/7}S^{-1/7}s_{1}^{6/7}\tau_{\rm A}^{-1}$\\
!| where the Lundquist number ($\tau_{\rm R}/\tau_{\rm A}$) is:\\
!| $S=43.7T_{\rm p0}^{3/2}(keV)\bar{r}_{1}^{2}/\tau_{\rm A}$.
!| 
!| \end{enumerate}
!| 
!
      module porcelli_constants

      REAL, PARAMETER :: &
       ctrho  = 1.0 & ! multiplier for rho_star in Eq. (15)
     , ctstar = 3.0 & ! multiplier for gamma_rho in Eq. (15)
     , cwrat  = 1.0 & ! multiplier 
     , chfast = 0.4 & ! multiplier c_h in Eq. (13)
     , zpi    = 3.141592653589793 &     ! pi
     , zmu0   = 1.256637061435917E-06 & ! magnetic permeability
     , zepslon = 1.0e-12
!
      INTEGER, PARAMETER :: &
       ifastswt  = 1 & !  on/off switch for testing equ. 13 (fast ion)
     , iflarmswt = 1 & !  on/off switch for testing equ. 14 (rotational effects)
     , ikinswt   = 1   !  on/off switch for testing equ. 15 (tearing nature)

!
!  chfast  = numerical multiplier for c_h (equ 16 of reference)
!  ctrho   = numerical multiplier for rho_star (equ 16 of reference)
!  ctstar  = numerical multiplier for c_star  (equ 16 of reference)
!  cwrat   = numerical multiplier depend on ratio between mode frequency 
!            and fast particle precesion freq. (c_f in equ B.8 of ref)
!  zpi     = pi
!  zmu0    = magnetic permeability
!  zepslon = small number
!
!  ifastswt  = on/off switch for testing equ. 13 (fast ion)
!  iflarmswt = on/off switch for testing equ. 14 (rotational effects)
!  ikinswt   = on/off switch for testing equ. 15 (tearing nature)

      end module porcelli_constants
!
!----------------------------------------------------------------
!
! This is a module for testing whether a sawtooth crashed had occurred.
! If a sawtooth occur module will return is_sawt=1, if not => is_sawt=0
!
       module porcelli_module
!
       implicit none
!
       contains
!
!*******
!
       subroutine porcelli_model (issawt, lprint, n_unit &
       , mzones, rminor_q_1 &
       , rminor_zone_center, beta_thermal_i, beta_thermal_e &
       , beta_total, shear, elongation &
       , rminor_zone_boundary, volume, area &
       , grad_beta_thermal, grad_beta_fast &
       , b_poloidal &
       , ion_mass, b_toroidal, r_major, denel0, ti0, efast )
!
!  Many of the variables in the argument list are profile arrays
!  which are used to prepare the scalar variables needed to
!  call "porcelli_conditions"
!
!..Output
!  -----
!  issawt  = 1 indicates that a sawtooth crash is triggered (output)
!            0 indicates that a sawtooth crash is not triggered
!
!..Input
!  -----
!  lprint  > 0 for diagnostic printout
!  n_unit = output unit number for diagnostic printout
!
!  mzones  = number of zone centers in grid
!  rminor_q_1 = minor radius of q = 1 surface [m] (input)
!        Note that the radius of the q=1 surface must be computed externally
!  rminor_zone_center(j) = minor radius profile for zone center grid [m]
!        Note that pressures are given on the zone center grid
!          j=1,mzones
!  beta_thermal_i(j) = ion thermal beta on zone center grid
!  beta_thermal_e(j) = electron thermal beta on zone center grid
!  beta_total(j)   = total pressure on zone center grid
!        Note, the total pressure includes fast ions
!  shear(j)   = magnetic shear profile on zone center grid
!    shear = r * ( d q / d r ) / q
!      where r = minor radius
!  elongation(j) = elongation of flux surfaces on zone center grid
!
!        Note that the following profiles are on the zone boundary grid
!          j=1,mzones+1
!  rminor_zone_boundary(j) = minor radius profile for zone boundary grid [m]
!  volume(j) = volume of flux sufraces on the zone boundary grid [m^3]
!  area(j)   = cross sectional area of flux sufraces [m^2]
!                on the zone boundary grid
!  grad_beta_thermal(j) = d beta_thermal / d r  [1/m] on the zone boundary grid
!  grad_beta_fast(j) = d beta_fast / d r  [1/m] on the zone boundary grid
!  b_poloidal(j) = poloidal magnetic field on zone boundary grid [Tesla]
!                  at outboard edge of midplane
!
!    Scalar variables:
!  ion_mass   = average ion mass (amu)
!  b_toroidal = toroidal magnetic field at r_major [T]
!  r_major    = major radius [m]
!  denel0  = electron density at magnetic axis 10^20 m^-3
!  ti0     = central ion temperature (keV)
!  efast   = energy of fast ions at q=1 surface (keV) (scalar)
!
!
       use porcelli_constants
!
       implicit none
!
       integer,intent(in) :: lprint, n_unit
!
       integer,intent(in) :: mzones
!
       integer,intent(out) :: issawt
!
       real, intent(in) :: rminor_q_1
!
       real, dimension(*), intent(in) :: &
         rminor_zone_center, beta_thermal_i, beta_thermal_e &
       , beta_total, shear, elongation
!
       real, dimension(*), intent(in) :: &
         rminor_zone_boundary, volume, area &
       , grad_beta_thermal, grad_beta_fast, b_poloidal
!
       real, intent(in) :: &
         ion_mass, b_toroidal, r_major, ti0, denel0, efast
!
       real :: &
         beta_i_max       & ! maximum thermal ion beta
       , beta_tot_max     & ! maximum total beta
       , pressure_at_q1   & ! total pressure at the q=1 surface
       , volume_at_q1     & ! volume of the q=1 magnetic surface
       , pressure_avg_up_to_q1 & ! average total pressure up to q=1 surface
       , elongation_at_q1 & ! elongation q=1 surface
       , b_pol_at_q1      & ! poloidal field (Tesla) at q=1 surface
                             !   (midplane outboard)
       , shear_at_q1      & ! magnetic shear at q=1 surface
       , grad_norm_press_at_q1 & ! thermal pressure scale length at q=1
       , area_at_q1      & ! cross sectional area of flux surface at q=1
       , inductance_at_q1 & ! internal inductance at q=1
       , r_minor          & ! minor radius
                             ! = rminor_zone_boundary(mzones+1)
       , cpf              & ! 
       , betapolf1          ! fast ion poloidal beta gradient parameter
!
       INTEGER :: ierror    ! error condition
!
       integer :: jz
!
!..temporary arrays
!
       real, dimension(3) :: ztemp_in, ztemp_out
!
       real, dimension(:), allocatable :: ztemp1
!
       real, dimension(:), allocatable :: ceqtmp
!
!..allocate temporary arrays
!
       ierror = 0
       allocate ( ztemp1(1:mzones), ceqtmp(1:4*mzones), STAT = ierror  )
!
       IF ( ierror > 0 ) THEN
         WRITE (n_unit,*) ' After allocate, ierror = ', ierror
       ENDIF
!
!..diagnostic output
!
       if ( lprint > 1 ) then
!
         write(n_unit,*)
         write(n_unit,101)
 101     format(' jz',t5,'rminor_zone_cent'  &
       , t24, 'beta_thermal_i', t40,'beta_thermal_e'  &
       , t56, 'beta_total', t72, 'shear', t88, 'elongation')
!
         do jz=1,mzones
           write(n_unit,102) jz, rminor_zone_center(jz)  &
           , beta_thermal_i(jz), beta_thermal_e(jz), beta_total(jz)  &
           , shear(jz), elongation(jz)
         enddo
!
 102     format(i4,1p6e16.5)
!
         write(n_unit,*)
         write(n_unit,103)
 103     format(' jz',t5,'rminor_zone_boundary'  &
       , t24, 'volume', t40,'area'  &
       , t56, 'grad_beta_thermal'  &
       , t72, 'grad_beta_fast'  &
       , t88, 'b_poloidal') 
!
         do jz=1,mzones
           write(n_unit,102) jz, rminor_zone_boundary(jz) &
           , volume(jz), area(jz) &
           , grad_beta_thermal(jz), grad_beta_fast(jz) &
           , b_poloidal(jz)
         enddo
!
         write(n_unit,*)
         write(n_unit,*) ' rminor_q_1=', rminor_q_1 &
           , ' ion_mass=', ion_mass &
           , ' b_toroidal=', b_toroidal &
           , ' r_major=', r_major
!
         write(n_unit,*) ' denel0=', denel0 &
           , ' ti0=', ti0 &
           , ' efast=', efast
!
       endif
!
!..compute, beta_i_max, the maximum thermal ion beta
!
       beta_i_max = 0.0
       beta_tot_max=0.0 
       do jz = 1,mzones
          beta_i_max = max (  beta_i_max, beta_thermal_i(jz) )
          beta_tot_max = max (  beta_tot_max, beta_total(jz) )
       enddo
!
       r_minor = rminor_zone_boundary(mzones)
!
       ztemp_in(1) = 0.0
       ztemp_in(2) = rminor_q_1
       ztemp_in(3) = rminor_zone_center(mzones-3)
!
!  Interpolate the total pressure to the q=1 surface
! 
       call cubint (rminor_zone_center, beta_total &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
        ,'abort computing total pressure at q=1 in sbrtn sawtst')
!
         pressure_at_q1 = ztemp_out(2) * b_toroidal**2 / ( 2.0 * zmu0 )
           ! total pressure nt/m**2 at q=1 surface
!
!..compute volume of the q=1 surface
!
       ztemp1 = 0.0
       do jz=1,mzones-1
         ztemp1(jz+1) = ztemp1(jz) &
          + ( volume(jz+1) - volume(jz) )
       enddo
!
       call cubint (rminor_zone_center, volume &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing volume at q=1 in sbrtn sawtst') 
!
       volume_at_q1 = ztemp_out(2)
!
! compute pressure_avg_up_to_q1 = average total pressure up to q=1 surface
!    (Nt/m^3)
!
       ztemp1 = 0.0
       do jz=1,mzones-1
         ztemp1(jz+1) = ztemp1(jz) &
          +  beta_total(jz+1) * ( volume(jz+1) - volume(jz) )
       enddo
!
       call cubint (rminor_zone_center, ztemp1 &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing sum vol*press at q=1 in sbrtn sawtst') 
!
       pressure_avg_up_to_q1 = ztemp_out(2) * b_toroidal**2 &
        / ( 2.0 * zmu0 * volume_at_q1 )
!
! compute elongation at q=1 surface
!
       call cubint (rminor_zone_center, elongation &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
         ,'abort computing elongation at q=1 in sbrtn sawtst')
!
       elongation_at_q1 = ztemp_out(2) ! elongation of q=1 surface
!
! compute b_pol_at_q1 = poloidal field (Tesla) at q=1 surface
!   (midplane outboard)
!
       call cubint (rminor_zone_center, b_poloidal &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
          ,'abort computing B-poloidal at q=1 in sbrtn sawtst')
!
       b_pol_at_q1 = ztemp_out(2)
!
! compute shear_at_q1 = magnetic shear at q=1 surface
!
       call cubint (rminor_zone_center, shear &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing shear at q=1 in sbrtn sawtst')
!
       shear_at_q1 = abs(ztemp_out(2))
!
!       if (zshearp < zepslon) zshearp = 2.0*(1.0-qsmth(2)) ! default
!
! compute 1. / pressure gradient scale length at the q=1 surface
!
       call cubint (rminor_zone_center, grad_beta_thermal &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing press. length at q=1 in sbrtn sawtst')        
!
       grad_norm_press_at_q1 = abs(ztemp_out(2) / pressure_at_q1 ) &
          * b_toroidal**2 / ( 2.0 * zmu0 )
       if ( grad_norm_press_at_q1 > 1. / zepslon) &
          grad_norm_press_at_q1 = 1. / rminor_q_1          ! default
!
! compute inductance_at_q1 = internal inductance at q=1
!
         ztemp1 = 0.0
         do jz=1,mzones-1
           ztemp1(jz+1) = ztemp1(jz) &
           + 0.5 * ( b_poloidal(jz+1)**2 + b_poloidal(jz)**2 ) &
                 * ( area(jz+1) - area(jz) )
         enddo
!
       call cubint (rminor_zone_center, ztemp1 &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing sum bpol**2da at q=1 in sbrtn sawtst') 
!
       inductance_at_q1 = ztemp_out(2)
!
       call cubint (rminor_zone_center, area &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing sum area at q=1 in sbrtn sawtst') 
!
       area_at_q1 = ztemp_out(2)
!
       inductance_at_q1 = inductance_at_q1 &
         / ( area_at_q1 * b_pol_at_q1**2 )
!
!| compute
!| \[ cpf=\frac{(5/2)\int_{0}^{1}x^{3/2}P_{\rm thermal}(x) dx}{P_{\rm i0}}\]
!
         ztemp1 = 0.0
         do jz=1,mzones-1
!
           ztemp1(jz+1) = ztemp1(jz) &
           + ( rminor_zone_boundary(jz+1) / rminor_q_1 )**1.5 &
             * ( ( beta_thermal_i(jz+1) + beta_thermal_e(jz+1) ) &
               / beta_tot_max ) &
           * ( rminor_zone_boundary(jz+1) - rminor_zone_boundary(jz) ) &
             / rminor_q_1
!
         enddo
!
          call cubint (rminor_zone_center, ztemp1 &
             ,mzones-1, 0, ceqtmp, mzones &
             ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
             ,'abort computing sum integ press at q=1 in sbrtn sawtst') 
!
!          cpf = 2.5 * ztemp_out(2)  ! correct expression
!
          cpf = 2.0 * ztemp_out(2)  ! incorrect form by Canh for testing
!
!| compute
!| \[zbetapolf1 =\frac{-2\mu_{0}}{B_{\rm p}^{2}}\int_{0}^{1}x^{3/2}
!|    \frac{dp_{\rm fast}}{dx} dx\]
!
         if ( efast < zepslon ) then
!
           betapolf1 = 0.0
!
         else
!
           ztemp1 = 0.0
           do jz=1,mzones-1
!
             ztemp1(jz+1) = ztemp1(jz) &
     &      + ( rminor_zone_boundary(jz+1) / rminor_q_1 )**1.5 &
     &      * ( rminor_zone_boundary(jz+1) - rminor_zone_boundary(jz) ) &
     &      * grad_beta_fast(jz+1) * b_toroidal**2 / ( 2.0 * zmu0 )
!
           enddo
!  
            call cubint (rminor_zone_center, ztemp1 &
     &        ,mzones-1, 0, ceqtmp, mzones &
     &        ,ztemp_in, ztemp_out, 2, 0, 0.0, 1 &
     &        ,'abort computing sum integ press at q=1 in sbrtn sawtst') 
!
            betapolf1 = - ztemp_out(2) * b_toroidal**2 / b_pol_at_q1**2
!
          endif
!
! temporary diagnostic printout
!
      write (6,*)
      write (6,*) 'cpf=',cpf,' elongation_at_q1=',elongation_at_q1
      WRITE (n_unit,*) 'ion_mass=',ion_mass
      write (6,*) 'inductance_at_q1=',inductance_at_q1
      WRITE (n_unit,*)' betapolf1=',betapolf1,' beta_i_max=',beta_i_max
      write (6,*) 'shear_at_q1=',shear_at_q1
      WRITE (n_unit,*)' r_major=',r_major &
       ,'  r_minor=', r_minor
      write (6,*) 'rminor_q_1=',rminor_q_1
      WRITE (n_unit,*)' grad_norm_press_at_q1=',grad_norm_press_at_q1
      WRITE (n_unit,*)' denel0=',denel0,' ti0=',ti0
      write (6,*) 'b_toroidal=',b_toroidal &
       ,' b_pol_at_q1=',b_pol_at_q1
      write (6,*) 'efast=',efast
      WRITE (n_unit,*)' pressure_avg_up_to_q1=',pressure_avg_up_to_q1 &
       ,' pressure_at_q1=',pressure_at_q1
      write (6,*) 'ifastswt=',ifastswt,' iflarmswt=',iflarmswt &
       ,'ikinswt=',ikinswt
      write (6,*)
!
!..deallocate temporary arrays
!
       ierror = 0
       DEALLOCATE ( ztemp1, ceqtmp, STAT=ierror )
!
       IF ( ierror > 0 ) THEN
         WRITE (n_unit,*) ' After deallocate, ierror = ', ierror
       ENDIF
!
! r_minor = minor radius (half-width) at edge of plasma
!
          call porcelli_conditions ( issawt, lprint, n_unit &
           , cpf &
           , elongation_at_q1, ion_mass, inductance_at_q1 &
           , betapolf1, beta_i_max, shear_at_q1 &
           , r_major, r_minor, rminor_q_1 &
           , grad_norm_press_at_q1, denel0 &
           , b_toroidal, b_pol_at_q1, ti0, efast &
           , pressure_avg_up_to_q1, pressure_at_q1 )
!
       return
!
       end subroutine porcelli_model
!------------------------------------------------------------------------
!
! subroutine for calculation of the porcelli model conditions
!
       subroutine porcelli_conditions ( issawt, lprint, n_unit &
       , cpf, elongation_at_q1, ion_mass &
       , inductance_at_q1, betapolf1 &
       , beta_i_max, shear_at_q1, r_major, r_minor, rminor_q_1 &
       , grad_norm_press_at_q1, denel0, b_toroidal, b_pol_at_q1, ti0 &
       , efast, pressure_avg_up_to_q1, pressure_at_q1 )
!
! All of the variables in the argument list of this routine are scalars
! that could be computed by the user or by calling sbrtn "porcelli_model".
!
!  issawt = output indicating whether a sawtooth crash should occur
!                   ( issawt=1 --> sawtooth crash,
!                     issawt=0 --> no sawtooth crash)
!
!  lprint  > 0 for diagnostic printout
!  n_unit = output unit number for diagnostic printout
!
!  cpf   = integral of x^(3/2)*p(x) p(x)=ion pressure
!                     (c_p of equ B.7 in reference) 
!  elongation_at_q1   = elongation at the q = 1 surface
!  ion_mass   = average ion mass (amu)
!  inductance_at_q1 = internal inductance at q = 1 surface
!  betapolf1   = fast ion beta parameter (equ 27 of reference)
!  beta_i_max   = peak ion toroidal beta
!  shear_at_q1   = rdq/dr at q=1 surface
!  r_major         = major radius (m)
!  r_minor      = minor radius (half-width) at edge of plasma (m)
!  rminor_q_1   = radius of q = 1 surface (m)
!  grad_norm_press_at_q1   = |dp/dr|/p at q=1 surface (m)
!    if p/|dp/dr|<<<1 take grad_norm_press_at_q1=radius of q=1 surface  
!  denel0      = electron density at magnetic axis (10^20 m-3)
!  b_toroidal  = vacuum toroidal field at r_major (testla)
!  b_pol_at_q1 = poloidal field at q = 1 surface (tesla) (outboard midplane)
!  ti0         = central ion temperature (keV)
!  efast       = energy of fast ion (keV)
!  pressure_avg_up_to_q1 = volume average pressure within q = 1 surface
!                            (pascal)
!  pressure_at_q1   = pressure at q = 1 surface (pascal)
!
       use porcelli_constants
!
       implicit none
!
       integer, intent(in) :: lprint, n_unit
!
       real, intent(in) :: cpf
       real, intent(in) :: elongation_at_q1, ion_mass &
         , inductance_at_q1 &
         , betapolf1, beta_i_max
       real, intent(in) :: shear_at_q1, r_major, r_minor, rminor_q_1 &
        , grad_norm_press_at_q1, denel0
       real, intent(in) :: b_toroidal, b_pol_at_q1, ti0, efast &
         ,pressure_avg_up_to_q1, pressure_at_q1
!
       integer, intent(out) :: issawt
!
       real :: zcel,zbetapc,zbetap1,zcmhd,zr1avg,zdiffbeta2,zdenom
       real :: zdelwbu,zdelwel,zdelwmhd,zdelwko,zdelwfast,zdelw
       real :: ztalf,zomegf,zomegdi,zrhostar,zlundnum,zgamp,zdelwc
       real :: zsub1,zsub2,zsub3,zsub4,zsub5,zsub6
!
!| 
!| \subsection{ Bussac term:}
!| Contributions from ideal MHD Bussac term is 
!| $\delta\hat{W}_{\rm Busac}=-c_{\rm MHD} \epsilon_{1}^{2}
!|   (\beta_{\rm p1}^{2}-\beta_{\rm pc}^{2})$
!| where\\
!| $\epsilon_{1}=\bar{r}_{1}/R$ (average q=1 radius/major radius)\\
!| $\beta_{\rm p1}=(8\pi/B_{\rm p1}^{2})[<p>_{1}- p(r_{1})]$
!|    (pressure: core + alpha)\\
!| $\beta_{\rm pc}=0.3[1-(5\bar{r}_{1}/3\bar{a})]$\\
!| $c_{\rm MHD}=\frac{9\pi}{s_{1}}(l_{\rm i1}-\frac{1}{2})$\\
!| $s_{1}=\bar{r}_{1}q'({\bar{r}_{1}})$ (magnetic shear parameter)\\
!| with $<>$ denoting the volume averaging within q = 1 surface,
!| $\bar{a}$ is  average plasma minor radius, and $l_{\rm i1}$ is the
!| internal inductance.
!| 
      !**********************
      !The Bussac term....delw_bu
      !**********************

       zr1avg = elongation_at_q1 * rminor_q_1 ! average q=1 radius
       zbetap1 = 8.0e-7*zpi &
                  * ( pressure_avg_up_to_q1 - pressure_at_q1 ) &
                  / b_pol_at_q1**2
 ! equ. 12 of Ref
       zbetapc = 0.3 * ( 1.0 - (5.0*zr1avg) / ( 3.0 * r_minor ) )
 ! equ. B.4 of Ref
       zcmhd = 9.0 * zpi * ( inductance_at_q1 - 0.5 )  / shear_at_q1
 !
 ! equ. B.3 of Ref
       zdiffbeta2 = zbetap1**2.0 - zbetapc**2.0
       zdelwbu=(-1.0)*zcmhd*(zr1avg/r_major)**2.0*zdiffbeta2
 ! equ. B.2 in Ref

!| 
!| \subsection{Elongation term}
!| Contributions from the elongation term is
!| $\delta\hat{W}_{\rm el}= -c_{\rm el}[\frac{\kappa_{1}-1}{2}]^{2}$
!| where\\
!| $c_{\rm el}=\frac{18\pi}{s_{1}}(l_{\rm i1}-\frac{1}{2})^{3}$\\
!| and $\kappa_{1}^{1/2}$ is the elongation of the q = 1 surface.
!|
!
      !**********************
      !The elongation term....delw_el
      !**********************

       zcel = 18.0 * zpi * ( inductance_at_q1 - 0.5 )**3.0 / shear_at_q1

  ! equ. B.6 of Ref
       zdelwel = (-1.0) * zcel * ((elongation_at_q1**2.0 - 1.0)/2.0)**2

  ! equ. B.5 of Ref
       zdelwmhd = zdelwel + zdelwbu
  ! MHD contributions
!| 
!| \subsection{Kruskal-Oberman term}
!| Contributions from the Kruskal-Oberman term is
!| $\delta\hat{W}_{KO}=0.6c_{\rm p}\epsilon_{1}^{1/2}\beta_{\rm i0}/s_{1}$\\
!| with
!| $c_{\rm p}=(5/2)\int_{0}^{1}dx x^{3/2}p_{\rm i}(x)/p_{\rm i0}$,
!|   ($x=r/r_{1}$)\\
!| where $p_{\rm i}(x)$ is the ion pressure profile, $p_{\rm i0}$
!|    is its peak value
!| and $\beta_{\rm i0}$ is the peak ion toroidal beta.  The expression
!| $\delta\hat{W}_{KO}$ applies to frequency regime
!| $\omega_{\ast\rm i}<\omega<\omega_{\rm bi}$.
!| 
!
      !**********************
      !The Kruskal-Oberman term....delw_ko
      !**********************

       zdelwko=0.6*cpf*sqrt(elongation_at_q1*rminor_q_1/r_major) &
         * beta_i_max / shear_at_q1
 ! equ. B.7 of Ref
       zdelwc = zdelwmhd+zdelwko ! core contributions
!| 
!| \subsection{Fast particle term}
!| Contributions from the fast particle term is
!| $\delta\hat{W}_{\rm fast}=c_{\rm f}\epsilon_{1}^{3/2}
!|    \beta_{\rm p\alpha}^{\ast}/s_{1}$
!| where\\
!| $\beta_{\rm p\alpha}^{\ast}=-\frac{8\pi}{B_{\rm p}^{2}(r_{1})}
!|    \int_{0}^{1}dx x^{3/2}\frac{dp}{dx}$,
!| ($x=r/r_{1}$)
!| and $c_{\rm f}$ depend depends on the ratio between the mode frequency 
!|   and the fast particle precessional frequency.

!
      !**********************
      !The fast particle term....delw_fast
      !**********************

       zdelwfast = cwrat * (elongation_at_q1*rminor_q_1/r_major)**1.5 &
         * betapolf1 / shear_at_q1
  ! equ. B.8 of Ref
       zdelw=zdelwbu+zdelwel+zdelwko+zdelwfast
  ! total contributions

      !***************
      ! conditions for sawtooth crash
      !***************

       ztalf = (0.8e-6) * (r_major / b_toroidal ) &
          * sqrt( ion_mass * denel0 )
  ! alfven time (sec)
       zomegf = 250.0 * efast &
        / ( b_toroidal * r_major * elongation_at_q1 * rminor_q_1 )
  ! fast ion prec. freq. (s^-1)
       zomegdi = 1000.0 * ti0 * grad_norm_press_at_q1 &
         /( b_toroidal * elongation_at_q1 * rminor_q_1 )
  ! ion (z=1) diamagnetic freq. (s^-1)
       zrhostar = 3.25e-3 * sqrt( ion_mass * ti0 ) &
         / ( b_toroidal * elongation_at_q1 * rminor_q_1 )
  ! ion larmor radius (z=1) to average rminor_q_1
       zlundnum = 43.7 * ti0**(1.5) &
         * elongation_at_q1 * rminor_q_1**2 / ztalf
  ! Lundquist number
       zdenom = ztalf * zlundnum**(1.0/7.0)
       zgamp  = 1.1 * shear_at_q1**(6.0/7.0) &
         * abs(zrhostar)**(4.0/7.0) / zdenom
  ! semi-collisional m =1 growth rate
       zdelwfast = cwrat &
         * ( elongation_at_q1 * rminor_q_1/r_major )**1.5 &
         * betapolf1 / max( shear_at_q1, zepslon )
  ! equ. B.8 of Ref

       zsub1 = 0.5*zomegdi*ztalf
       zsub2 = -1.0*zdelwc
       zsub3 = chfast*ztalf*zomegf
       zsub4 = -1.0*zdelw
       zsub5 = -1.0*ctrho*zrhostar
       zsub6 = 1.0*ctstar*zgamp

       issawt = 0
       IF ( ifastswt > 0) THEN
         IF ( zsub2 > zsub3 ) THEN
           issawt = 1           ! equ. 13 in Ref, condition 1
           IF ( lprint > 0 ) THEN
             WRITE (n_unit,*) " ***====> sawtooth is triggered by Eq. (13)"
           ENDIF
         ENDIF
         IF ( lprint > 0 ) THEN
           WRITE (n_unit,*) ' Porcelli sawtooth trigger condition Eq. (13)'
           WRITE (n_unit,&
&"(t5,'- delta W_core',t29,'>',t31,'c_h omega_Dh tau_A')")
           WRITE (n_unit,"(1x,3ES25.15)") zsub2,zsub3
           WRITE (n_unit,*)
         ENDIF
       ENDIF

       IF ( iflarmswt > 0) THEN
         IF ( zsub4 > zsub1 ) THEN
           issawt = 1           ! equ. 14 in Ref, condition 2
           IF ( lprint > 0 ) THEN
             WRITE (n_unit,*)" ***====> sawtooth is triggered by Eq. (14)"
           ENDIF
         ENDIF
         IF ( lprint > 0 ) THEN
           WRITE (n_unit,*) ' Porcelli sawtooth trigger condition Eq. (14)'
           WRITE (n_unit,&
&"(t5,' - delta W',t29,'>',t31,'0.5 omega_*i tau_A')")
           WRITE (n_unit,"(1x,3ES25.15)") zsub4, zsub1
           WRITE (n_unit,*)
         ENDIF
       ENDIF

       IF ( ikinswt > 0) THEN
         IF ( zsub5 < zsub4 .and. zsub4 < zsub1 &
                          .and. zomegdi < zsub6 ) THEN
           issawt = 1           ! equ. 15 on Ref, condition 3
           IF ( lprint > 0 ) THEN
             WRITE(n_unit,*) " ***====> sawtooth is triggered by Eq. (15)"
           ENDIF
         ENDIF
         IF ( lprint > 0 ) THEN
           WRITE (n_unit,*) ' Porcelli sawtooth trigger condition Eq. (15)'
           WRITE (n_unit,"(t5,' -c_rho rho',&
               &t29,'<',t31,'- delta W',&
               &t54,'<',t56,'0.5 omega_*i tau_A')")
           WRITE (n_unit,"(1x,3ES25.15)") zsub5, zsub4, zsub1 
           WRITE (n_unit,"(t2,'and omega_*i',&
               &t29,'<',t31,'c_* gamma_rho')")
           WRITE (n_unit,"(1x,3ES25.15)") zomegdi,  zsub6
           WRITE (n_unit,*)
         ENDIF
       ENDIF

       IF ( issawt == 0 .and. lprint > 0 ) then
          write(n_unit,*)&
&" ***====> no sawtooth crash triggered by Porcelli module"
       ENDIF
!
!..diagnostic printout
!
      if ( lprint > 0 ) then
        WRITE (n_unit,*)
        ! writing the various delw into file
        write(n_unit,110) zdelwbu,zdelwel,zdelwmhd,zdelwko,zdelwc, &
                          zdelwfast,zdelw,ztalf,zomegf,zomegdi,zrhostar, &
                          zr1avg,zlundnum,zgamp,zbetap1
!        write(n_unit,120) zdelwbu,zdelwel,zdelwko,zdelwfast,zdelw
!        write(n_unit,130) zsub2,zsub3,zsub4,zsub1,zsub5,zomegdi,zsub6
      endif
!
  110  format("****** ",/, &
              "delw_bu = ",1pe10.3,2x,"delw_el = ",1pe10.3,2x, &
              "delw_mhd = ",1pe10.3, /,"delw_ko = ",1pe10.3,2x, &
              "delw_c = ",1pe10.3,2x,"delw_fast = ",1pe10.3,2x,/, &
              "delw = ",1pe10.3, /,"t_alf(s) = ",1pe10.3,2x, &
              "omeg_f(s^-1) = ",1pe10.3,2x,"omeg_di(s^-1) = ",1pe10.3,/, &
              "rho_star = ",1pe10.3,2x,"r1_avg(m) = ",1pe10.3,2x, &
              "lund_num = ",1pe10.3,2x,/, &
              "semi-coll. m=1 growth rate (s^-1) = ", &
               1pe10.3,2x,"beta_p1 = ",1pe10.3,/, &
              "****************************", /)
!  120   format(2x,5(1pe10.3,2x))
!  130   format(2x,7(1pe10.3,2x)) 
!
!
        return
!
        end subroutine porcelli_conditions
!
!---------------------------------------
!
       end module porcelli_module
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{porcelli96a} F.~Pocelli, D.~Boucher, and M.~N.~Rosenbluth,
!| Plasma Physics and Controlled Fusion, {\bf 38}, 2163 (1996).
!| 
!| 
!| \end{thebibliography}
!| \end{document}

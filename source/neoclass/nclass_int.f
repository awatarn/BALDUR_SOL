!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| 
!| \documentstyle[12pt]{article}
!| \headheight 0pt \headsep 0pt
!| \topmargin 0pt  \textheight 9.0 in
!| \oddsidemargin 0pt \textwidth 6.5 in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center}
!| \large {\bf NCLASS Interface} \\
!| \normalsize  {\tt nclass_int.f} \\
!| \vspace{1pc}
!| Alexei Pankin, Glenn Bateman \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| \today
!| \end{center}
!| 
!| This module is interface to NCLASS, W. Houlberg's module that
!| computes neoclassical transport.
!| 
c@nclass_int  .../baldur/bald/nclass_int.f
c rgb 08-aug-01 changed the sign of vnneo1 to conform to BALDUR
c rgb 07-aug-01 implemented bootstrap current ajboot(jz)
c rgb 18-jul-01 replace function ftrapfn with array trapbr(jz,1)
c rgb 13-jul-01 implemented z_grbm2=delxi27b2i, z_ngrth=bdotdeltheta
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine nclass_int
!
! NCLASS options:
!
! lneocl(1)      neoclassical model 
!         = 1    to call nclass
!         = else to call trneo1
!
! lneocl(2)      version of interface to nclass
!         = 9    simplified interface for testing
!
! lneocl(3)      computation of trapped particle fraction
!         = 0    original analytic BALDUR compuation
!         = 1    Mike Hughes numerical computation
!         = else  trapped particle fraction = sqrt(1.0-|B_min/B_max|)
!
! lneocl(4)      option for output to nout (-)
!         > 0    errors only
!         = 2    errors and results
!         = else no output
! lneocl(11)      order of v moments to be solved (-)
!         = 2    u and p_q
!         = 3    u, p_q, and u2
!         = else 2
! lneocl(12)      option to include potato orbits (-)
!         = 0    off
!         = else on
!
! cfutz(19) particle pinch multiplier
!
!
! This sbrtn calls the following subprograms by Wayne Houlberg:
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! subroutine nclass ()
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        include 'nclass.f'
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! subroutine u_lu_backsub ()
!   - U_LU_BACKSUB solves the matrix equation a x = b 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        include 'u_lu_backsub.f'
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! subroutine u_lu_decomp ()
!   - LU decomposition of the matrix by W. Houlberg
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        include 'u_lu_decomp.f'
!
!
!  c_den-density cutoff below which species is ignored (/m**3)
!  p_eps-inverse aspect ratio (-)
!  p_grphi-radial electric field Phi' (V/rho)
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_q-safety factor (-)
!  r0-major radius (m)
!  a0-minor radius, scale length for rho (m)
!  e0-axial elongation (-)
!  bt0-axial toroidal field (T)
!  q0-axial safety factor (-)
!
!  mx_mz-max charge of any species (-)
!  mx_mi-max isotopes (-)
!  mx_ms-max species (-)
!
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
!
! Parameter settings
!
       include 'nclass_par.m'
!
! The following are in ../com/nclass_par.m
!
!       integer   mx_mz, mx_mi, mx_ms
!       parameter ( mx_mz=18, mx_mi=9, mx_ms=40 )
!
! mx_mz = maximum number of charge states
! mx_mi = maximum number of isotopes
! mx_ms = maximum number of densities used
!
!  Physical and conversion constants
!
       integer :: k_order, k_potato
       integer :: i_maxcharge    ! Max charge state
       integer :: m_i            ! Number of isotops
       real    :: z_aminor       ! plasma minor radius (half-width) [m]
       real    :: c_den                    
! Density cutoff below which species is ignored
       real    :: c_potb         ! kappa(0)*Bt(0)/(2*q(0)^2) (T)
       real    :: c_potl         ! q(0)*R(0) (m)
       real    :: z_b22di(kjbal) ! <|B|^2>
       real    :: z_bm22di(kjbal) ! <1/|B|^2>
       real    :: amu_i(mx_mi)   ! ion mass
       real    :: z_eb(kjbal)    ! <EB> arrays
       real    :: z_ft(kjbal)    ! trapped fraction
       real    :: z_grbm2(kjbal) ! <grad(rho)**2/B**2>
       real    :: z_fhat(kjbal)
       real    :: z_grphi(kjbal) 
! radial electric field Phi'. To be implemented later
       real    :: z_gr2phi(kjbal)
! radial electric field gradient Psi'(Phi'/Psi')'. To be implemented later
       real    :: z_fm(3, kjbal)  
! poloidal moments of the geometric factor for PS viscosity
       real    :: z_ngrth(kjbal) ! <n.grad(theta)>
       real    :: z_t(mx_mi, kjbal)      ! Temperature in keV
       real    :: z_den(mx_mi, mx_mz, kjbal)  ! Density in m^3
       real    :: z_grt(mx_mi, kjbal)    ! Isotop temperature gradient
       real    :: z_grd(mx_mi, mx_mz, kjbal)  ! Isotop density gradient
       real    :: z_grp(mx_mi, mx_mz, kjbal)
! Isotop/Chrage_State pressure gradient
       real    :: z_fex(3, mx_mi, mx_mz, kjbal)
! moments of external parallel force
!       real    :: z_rho(kjbal)   ! sqrt(toroidal flux)
!
         real    :: z_aspinv(kjbal)   ! inverse aspect ratio
         integer :: jz, jm, icharge
!
!         real    :: z_aj(kjbal, kfour)  ! Jacobian
!         real    :: z_ajint(kjbal)
! integral(Jacobian, theta, 0, 2Pi)
!
         real    :: z_rminor(kjbal)  ! r_minor
         real    :: z_epsqrt
         real    :: z_btor0    ! B_tor_0
         real    :: z_rmajor0  ! R_major_0
         real    :: z_elong0   ! elongation at mag axis
         real    :: z_q0       ! magnetic q at mag axis
!
         real    :: z_temp(kfour,kjbal)
     &       , z_temp1(kjbal)
     &       , zcav(kjbal, 3) 
     &       , z_grbm2f(kjbal, kfour)
     &       , z_grbm2r(kjbal, kfour)
!
         integer  :: ierror
         real     :: fcharge
!
! declarations from sbrtn nclass_run
!
! Local variables
!
         real, parameter      :: zpi = 3.141593
         real , parameter     :: z_j7kv = 1.6022e-16
! conversion factor for Joules to keV
!
         integer  :: m_z, im, iz, lm, lz, km, kz
         real     :: z_xineo1t, z_xineo1p
         real     :: z_z, z_dtdr
         real     :: z_zflp, z_zflt
	 real     :: z_flxt1(mx_mi,kjbal)
     &     , z_flxt2(kjbal)  ! thermal flux
	 real     :: z_flxp1(mx_mi,kjbal)
     &     , z_flxp2(kjbal)  ! particle flux
	 real     :: zxineo1(kjbal)
	 real     :: z_gfl, z_qfl
	 integer  :: lastWarning
!
! Local variables: output from NCLASS
!
         integer  :: iflag,        ! Warrning and error flag
     >               m_s,          ! Number of species
     >               im_s(mx_ms),  ! isotope number of s (-)
     >               iz_s(mx_ms)   ! charge state of s (-)
         real     :: p_bsjb,       ! <J_bs.B> (A*T/m2)
     >               p_etap,       ! parallel electrical resistivity (Ohm*m)
     >               p_exjb        
! <J_ex.B> current response to fex_iz (A*T/m2)
         real     :: calm_i(3,3,mx_mi),  ! tp eff friction matrix for i (kg/m3/s)
     >               caln_ii(3,3,mx_mi,mx_mi),! fp eff friction matrix for i1 on i2 (kg/m3/s)
     >               capm_ii(3,3,mx_mi,mx_mi),! test part (tp) friction matrix for i1 on i2 (-)
     >               capn_ii(3,3,mx_mi,mx_mi) ! field part (fp) friction matrix for i1 on i2 (-)
         real     :: bsjbp_s(mx_ms),     ! <J_bs.B> driven by unit p'/p of s (A*T*rho/m3)
     >               bsjbt_s(mx_ms),     ! <J_bs.B> driven by unit T'/T of s (A*T*rho/m3)
     >               dn_s(mx_ms),        ! diffusion coefficient (diag comp) of s (rho2/s)
     >               gfl_s(5,mx_ms),     ! radial particle flux comps of s (rho/m3/s)
     >               qfl_s(5,mx_ms),     ! radial heat flux comps of s (rho/m3/s)
     >               sqz_s(mx_ms),       ! orbit ssqueezing factor for s (-)
     >               upar_s(3,3,mx_ms),  ! parallel flow of s from force m (T*m/s)
     >               utheta_s(3,3,mx_ms),! poloidal flow of s from force m (m/s/T)
     >               vnnt_s(mx_ms),      ! convection velocity (off diag comps-p', T') of s (rho/s)
     >               vneb_s(mx_ms),      ! <E.B> particle convection velocity of s (rho/s)
     >               vnex_s(mx_ms),      ! extrenal force particle convection velocity of s (rho/s)
     >               vqnt_s(mx_ms),      ! conduction velocity (off diag comps-p',T') of s (rho2/s)
     >               vqeb_s(mx_ms),      ! <E.B> heat convection velocity of s (rho/s)
     >               vqex_s(mx_ms),      ! extrenal force heat convection velocity of s (rho/s)
     >               xi_s(mx_ms),        ! charge weighted density factor of s (-)
     >               ymu_s(3,3,mx_ms),   ! normalized viscosity for s (kg/m3/s)
     >               chip_ss(mx_ms,mx_ms),!heat cond coefficient of s2 on p'/p of s1 (rho2/s)
     >               chit_ss(mx_ms,mx_ms),!heat cond coefficient of s2 on T'/T of s1 (rho2/s)
     >               chi_s(mx_ms),       ! heat conduction coefficient (diag comp) of s (rho2/s)
     >               dp_ss(mx_ms,mx_ms), ! diffusion coefficient of s2 on p'/p of s1 (rho2/s)
     >               dt_ss(mx_ms,mx_ms)  ! diffusion coefficient of s2 on T'/T of s1 (rho2/s)
!
         character*(80) :: message             ! warning or error message (character)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!  initialize variables
!
         if (lneocl(11) == 3) then ! number of v moments to be solved
           k_order  = 3
         else
           k_order  = 2
         endif
         if (lneocl(12) == 1) then ! option to include potato orbits
           k_potato = 1
         else
           k_potato = 0
         endif
!
!..variables that are used for the interface in general
!
         z_aminor = avi(mzones,15,1) ! plasma minor radius (half-width) [m]
!
! values at magnetic axis:
!
         z_btor0     = avi(2,8,1) / avi(2,14,1)
         z_rmajor0   = avi(2,14,1)
         z_elong0    = avi(2,16,1)
         z_q0        = q(2)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!..zero out all variables:
!
      c_potb = 0.0
      c_potl = 0.0
!
      z_b22di  = 0.0
      z_bm22di = 0.0
      z_eb     = 0.0
      z_fhat   = 0.0
      z_fm     = 0.0
      z_ft     = 0.0
      z_grbm2  = 0.0
      z_grphi  = 0.0
      z_gr2phi = 0.0
      z_ngrth  = 0.0
      amu_i    = 0.0
      z_grt    = 0.0
      z_t      = 0.0
      z_den    = 0.0
      z_fex    = 0.0
      z_grp    = 0.0
!
      m_s      = 0
      im_s     = 0
      iz_s     = 0
!
      p_bsjb   = 0.0
      p_etap   = 0.0
      p_exjb   = 0.0
      calm_i   = 0.0
      caln_ii  = 0.0
      capm_ii  = 0.0
      capn_ii  = 0.0
      bsjbp_s  = 0.0
      bsjbt_s  = 0.0
      chi_s    = 0.0
      dn_s     = 0.0
      gfl_s    = 0.0
      qfl_s    = 0.0
      sqz_s    = 0.0
      upar_s   = 0.0
      utheta_s = 0.0
      vnnt_s   = 0.0
      vneb_s   = 0.0
      vnex_s   = 0.0
      vqnt_s   = 0.0
      vqeb_s   = 0.0
      vqex_s   = 0.0
      xi_s     = 0.0
      ymu_s    = 0.0
      chip_ss  = 0.0
      chit_ss  = 0.0
      dp_ss    = 0.0
      dt_ss    = 0.0
!
      iflag    = 0
!
         vineo1       = 0.
         xineo1       = 0. 
         zxineo1      = 0.        
         z_flxt1      = 0.
         z_flxt2      = 0.
         z_flxp1      = 0.
         z_flxp2      = 0.
	 dnneo1       = 0.
         vnneo1       = 0.
	 dnhhs        = 0.
	 dniis        = 0.
         dnhis        = 0.
	 dnihs        = 0.
!
! bootstrap current
!
         ajboot       = 0.
         ajtpbi       = 0.
!
	 lastWarning  = 0
!	 if (present(z_jbb)) z_jbb     = 0.
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!..set up input variables for NCLASS:
!
! Note:
! b22di and bm22di come from common /commhd/
!
!..densities, temperatures, and their gradients
!
         c_den            = 1.0e10         
! Density cutoff below which species are ignored
         amu_i            = 0.0             !
         amu_i(1)         = 5.4463e-4       ! atomic mass of electrons
         m_i              = mhyd + mimp + 1 ! number of isotops
         amu_i(2:m_i)     = aspec(1:m_i-1)  ! atomic mass numbers
         z_t              = 0.0
         z_t(1, 1:kjbal)  = tes(2, 1:kjbal) * useh         
! electron temperature [keV]
         do i = 2, m_i
           z_t(i,1:kjbal) = tis(2, 1:kjbal) * useh
! other isotops temperature = Ti in BALDUR [keV]
         end do
         i_maxcharge      = 1
         z_den            = 0.0
         z_den(1, 1, 1:kjbal) = rhoels(2, 1:kjbal) * usid
! electron density [m^3]
         do jz = 1, mzones
           do i = 2, mhyd + 1
             icharge    = 1                   ! Charge state 
             z_den(i, icharge, jz) = rhohs(i-1, 2, jz) * usid
! Density of hydrogenic particles.
           end do  
           do i = mhyd + 2, m_i
             icharge   = floor(cmean(i-mhyd-1, 2, jz)) + 1
! Charge state at zone boundaries
             fcharge   = icharge - cmean(i-mhyd-1, 2, jz)
! Fraction part of the charge state
             i_maxcharge = max(i_maxcharge, icharge)
             z_den(i, icharge, jz)   = rhois(i-mhyd-1, 2, jz)
     >                               * (1. - fcharge) * usid
! Density of impurities
! There is no distinction in charge state yet
             z_den(i, icharge-1, jz) = rhois(i-mhyd-1, 2, jz)
     >                               * fcharge * usid
! To be implemented later
           end do
         end do
         	
!         z_rho(1:kjbal)   = xbouni(1:kjbal)   ! sqrt(Toroidal flux)
!
         z_grt            = 0.0
         z_grp            = 0.0
         do i = 1, m_i
           do jz = 2, mzones
             z_grt(i, jz) = (z_t(i, jz) - z_t(i, jz-1)) /
     &                      ( (armins(jz,1) - armins(jz-1,1)) * usil )
! dT(i) /d rho
             do j = 1, i_maxcharge
               z_grp(i, j, jz) = (z_t(i, jz)*z_den(i, j, jz)
     >                           - z_t(i, jz-1)*z_den(i, j, jz-1)) /
     &                      ( (armins(jz,1) - armins(jz-1,1)) * usil )
! dp(i,j) /d rho
               z_grd(i, j, jz) = (z_den(i, j, jz) - z_den(i, j, jz-1)) /
     &                      ( (armins(jz,1) - armins(jz-1,1)) * usil )
! dn(i,j) /d rho
             end do 
           end do
         end do
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
!..crude approximations to flux surface averages for testing
!
      if ( lneocl(2) .eq. 9 ) then
!
!
         c_potb  = elong(2, 2) * z_btor0 / ( 2. * q(2) * q(2) )
! kappa(0)*Bt(0)/[2*q(0)**2]

         c_potl  = q(2) * avi(2,14,1)           ! q0*R0
!
! loop over zone boundaries
!
         do jz = 3, mzones
!
           z_rminor(jz) = avi(jz,15,1)
           z_aspinv(jz) = avi(jz,15,1) / avi(jz,14,1)
!
           z_b22di(jz)  = ( 1.0 + 0.5 * z_aspinv(jz)**2 ) * z_btor0**2
           z_bm22di(jz) = ( 1.0 + 1.5 * z_aspinv(jz)**2 ) / z_btor0**2
!
! simplified expression for z_eb(jz) = <E.B>
! vloopi(jz,1) = loop voltage at zone boundary jz at time n
! avi(jz,8,1)  = R B-toroidal   at zone boundaries
! avi(jz,10,1) = < 1. / R**2 >   at zone centers
!
           z_eb(jz) = vloopi(jz,1) * avi(jz,8,1) *
     &        avi(jz,10,1) / ( 2.0 * zpi )
!
           z_fhat(jz) = q(jz) / z_aspinv(jz)
!
           z_epsqrt    = sqrt( 1.0 - z_aspinv(jz)**2 )
!
! loop over harmonics
           do jm = 1, 3
             z_fm(jm, jz) =
     &         jm * ( ( 1.0 - z_epsqrt ) / z_aspinv(jz) )**(2.0*jm)
     &         * ( 1.0 + jm * z_epsqrt )
     &          / ( ( 1.0 - z_aspinv(jz)**2 )**1.5
     &                     * ( q(jz) * avi(jz,14,1) )**2)
! poloidal moments of geometric factor for
! PS viscosity. Approximation formula is used.
! To be replaced by actual values later
! Note:
!   avi(jz,14,1) = major radius to center of zone boundary
!   avi(jz,15,1) = minor radius (half-width) of zone boundary
!
           end do
!
           z_ft(jz) = 1.46 * sqrt ( abs ( z_aspinv(jz) ) )
! trapped fraction as computed in nclass_pt_dr.f
!
! approximate expression for < |grad rho|^2 / B^2 >
!
           z_grbm2(jz) = 1.0 / z_btor0**2
!
! z_grphi  is left at 0.0
! z_gr2phi is left at 0.0
!
           z_ngrth(jz)   = 1.0 / ( q(jz) * z_rmajor0 )
! coordinate system normalization factor
! <n.grad(theta)>. To be replaced by
! actual values later
!
         end do
! end of loop over zone boundaries

         z_fex   = 0.0  ! Moments of external parallel force
                                               ! To be replaced later
!
! 2Pi*R*Bt/(d Phi/d R)
!
         z_eb(1:2)     = z_eb(3)
         z_grbm2(1:2)  = 0.
         z_fm(1:3, 1)  = z_fm(1:3, 3)
         z_fm(1:3, 2)  = z_fm(1:3, 3)
         z_fhat(1:2)   = z_fhat(3)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! debug value of lneocl(2) = 8
!
      else if ( lneocl(2) == 8 ) then
!
!..finite aspect ratio approximations to flux surface averages
!
! Note:
!   avi(jz,8,1)  = R B-toroidal   at zone boundaries
!   avi(jz,10,1) = < 1. / R**2 >   at zone centers
!   avi(jz,14,1) = major radius to center of zone boundary
!   avi(jz,15,1) = minor radius (half-width) of zone boundary
!
         c_potb  = elong(2, 2) * z_btor0 / ( 2. * q(2) * q(2) )
! kappa(0)*Bt(0)/[2*q(0)**2]

         c_potl  = q(2) * avi(2,14,1)           ! q0*R0
!
! loop over zone boundaries
!
         do jz = 3, mzones
!
           z_rminor(jz) = avi(jz,15,1)
!
           z_aspinv(jz) = avi(jz,15,1) / avi(jz,14,1)
! inverse aspect ratio = half-width /
!  major radius to midpoint at zone bndries
!
         z_b22di(jz)  = b22di(jz)         !  <|B|^2>
         z_bm22di(jz) = bm22di(jz)        !  <1/|B|^2>
!
         z_grbm2(jz)  = delxi27b2i(jz)    !  < |grad rho|^2 / B^2 >
!
         z_ngrth(jz)  = bdotdeltheta(jz)  !  < B. grad(theta) / |B| >
!
           z_fhat(jz) = q(jz) / z_aspinv(jz)
!
! poloidal moments of geometric factor for PS viscosity.
!
           z_epsqrt    = sqrt( 1.0 - z_aspinv(jz)**2 )
! loop over harmonics
           do jm = 1, 3
             z_fm(jm, jz) =
     &         jm * ( ( 1.0 - z_epsqrt ) / z_aspinv(jz) )**(2.0*jm)
     &         * ( 1.0 + jm * z_epsqrt )
     &          / ( ( 1.0 - z_aspinv(jz)**2 )**1.5
     &                     * ( q(jz) * avi(jz,14,1) )**2)
           end do
!
!           z_ft(jz) = trapbr (jz, 1)
!
! Calculation of trapped fraction z_ft
!
! simplified expression for z_eb(jz) = <E.B>
! vloopi(jz,1) = loop voltage at zone boundary jz at time n
!
           z_eb(jz) = vloopi(jz,1) * avi(jz,8,1)
     &                * avi(jz,10,1) / ( 2.0 * zpi )
!
! approximate expression for trapped particle fraction
!
           z_ft(jz) = 1.46 * sqrt ( abs ( z_aspinv(jz) ) )
!
! approximate expression for < n . grad theta >
!
           z_ngrth(jz)   = 1.0 / ( q(jz) * z_rmajor0 )
!
! approximate expression for < |grad rho|^2 / B^2 >
!
           z_grbm2(jz) = 1.0 / z_btor0**2
!
         end do
! end of loop over zone boundaries

         z_fex   = 0.0  ! Moments of external parallel force
!
         z_eb(1:2)     = z_eb(3)
         z_grbm2(1:2)  = 0.
         z_fm(1:3, 1)  = z_fm(1:3, 3)
         z_fm(1:3, 2)  = z_fm(1:3, 3)
         z_fhat(1:2)   = z_fhat(3)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
! default value of lneocl(2)
!
      else
!
!..finite aspect ratio approximations to flux surface averages
!
! Note:
!   avi(jz,8,1)  = R B-toroidal   at zone boundaries
!   avi(jz,10,1) = < 1. / R**2 >   at zone centers
!   avi(jz,14,1) = major radius to center of zone boundary
!   avi(jz,15,1) = minor radius (half-width) of zone boundary
!
         c_potb  = elong(2, 2) * z_btor0 / ( 2. * q(2) * q(2) )
! kappa(0)*Bt(0)/[2*q(0)**2]

         c_potl  = q(2) * avi(2,14,1)           ! q0*R0
!
! loop over zone boundaries
!
         do jz = 3, mzones
!
           z_rminor(jz) = avi(jz,15,1)
!
           z_aspinv(jz) = avi(jz,15,1) / avi(jz,14,1)
! inverse aspect ratio = half-width /
!  major radius to midpoint at zone bndries
!
         z_b22di(jz)  = b22di(jz)         !  <|B|^2>
         z_bm22di(jz) = bm22di(jz)        !  <1/|B|^2>
!
         z_grbm2(jz)  = delxi27b2i(jz)    !  < |grad rho|^2 / B^2 >
!
         z_ngrth(jz)  = bdotdeltheta(jz)  !  < B. grad(theta) / |B| >
!
           z_fhat(jz) = q(jz) / z_aspinv(jz)
!
! poloidal moments of geometric factor for PS viscosity.
!
           z_epsqrt    = sqrt( 1.0 - z_aspinv(jz)**2 )
! loop over harmonics
           do jm = 1, 3
             z_fm(jm, jz) =
     &         jm * ( ( 1.0 - z_epsqrt ) / z_aspinv(jz) )**(2.0*jm)
     &         * ( 1.0 + jm * z_epsqrt )
     &          / ( ( 1.0 - z_aspinv(jz)**2 )**1.5
     &                     * ( q(jz) * avi(jz,14,1) )**2)
           end do
!
           z_ft(jz) = trapbr ( jz, 1 )
!
! Calculation of trapped fraction z_ft
!
! simplified expression for z_eb(jz) = <E.B>
! vloopi(jz,1) = loop voltage at zone boundary jz at time n
!
           z_eb(jz) = vloopi(jz,1) * avi(jz,8,1)
     &                * avi(jz,10,1) / ( 2.0 * zpi )
!
         end do
! end of loop over zone boundaries

         z_fex   = 0.0  ! Moments of external parallel force
!
         z_eb(1:2)     = z_eb(3)
         z_grbm2(1:2)  = 0.
         z_fm(1:3, 1)  = z_fm(1:3, 3)
         z_fm(1:3, 2)  = z_fm(1:3, 3)
         z_fhat(1:2)   = z_fhat(3)
!
      endif
!
!  end of if-block controlled by lneocl(2)
!
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!
         vineo1       = 0.
         xineo1       = 0. 
         zxineo1      = 0.        
         z_flxt1      = 0.
         z_flxt2      = 0.
         z_flxp1      = 0.
         z_flxp2      = 0.
	 dnneo1       = 0.
         vnneo1       = 0.
	 dnhhs        = 0.
	 dniis        = 0.
         dnhis        = 0.
	 dnihs        = 0.
	 lastWarning  = 0
!	 if (present(z_jbb)) z_jbb     = 0.
!
! diagnostic printout
!
           if (lneocl(4) == 2) then
             write(nout,*)
             write(nout,*) ' Diagnostic printout before calling NCLASS'
             write(nout,*) 'k_order = ',k_order
             write(nout,*) 'k_potato = ',k_potato
             write(nout,*) 'm_i = ',m_i
             write(nout,*) 'i_maxcharge = ',i_maxcharge
             write(nout,*) 'c_den = ',c_den
             write(nout,*) 'c_potb = ',c_potb
             write(nout,*) 'c_potl = ',c_potl
             write(nout,*) 'amu_i = ',amu_i
           end if
!
! loop over zone boundaries
!
         do jz = 3, mzones
!
! diagnostic printout
!
           if (lneocl(4) == 2) then
             write(nout, '("  -- jz = ", I4)')    jz
             write(nout,*) ' z_aspinv = ',z_aspinv(jz)
     &  ,'  rminor  = ',z_rminor(jz)
     &  ,'  rmajor0 = ',z_rmajor0
     &  ,'  elong0 = ',z_elong0
             write(nout,*)
     &   '  q(2) = ',q(2)
     &  ,'  b_tor0 = ',z_btor0
             write(nout, '("b2      = ", 1pE16.5)') z_b22di(jz)
             write(nout, '("bm2     = ", 1pE16.5)') z_bm22di(jz)
             write(nout, '("eb      = ", 1pE16.5)') z_eb(jz)
             write(nout, '("fhat    = ", 1pE16.5)') z_fhat(jz)
             write(nout, '("fm      = ", 3(1pE16.5))') z_fm(1:3, jz)
             write(nout, '("ft      = ", 1pE16.5)') z_ft(jz) 
             write(nout, '("grbm2   = ", 1pE16.5)') z_grbm2(jz) 
             write(nout, '("grphi   = ", 1pE16.5)') z_grphi(jz)
             write(nout, '("gr2phi  = ", 1pE16.5)') z_gr2phi(jz) 
             write(nout, '("ngrth   = ", 1pE16.5)') z_ngrth(jz)
             write(nout, '("amu_i   = ", 9(1pE16.5:))') amu_i(1:m_i) 
             write(nout, '("grt     = ", 9(1pE16.5:))') z_grt(1:m_i, jz) 
             write(nout, '("t       = ", 9(1pE16.5:))') z_t(1:m_i, jz) 
             write(nout, '("den1s   = ", 9(1pE16.5:))') 
     >             z_den(1, 1:7, jz)
             write(nout, '("den2s   = ", 9(1pE16.5:))') 
     >             z_den(2, 1:7, jz)
             write(nout, '("den3s   = ", 9(1pE16.5:))') 
     >             z_den(3, 1:7, jz)
             write(nout, '("den4s   = ", 9(1pE16.5:))') 
     >             z_den(4, 1:7, jz)
             write(nout, '("fex     = ", 27(1pE16.5:))') 
     >             z_fex(1:3, 1:m_i, 1, jz) 
             write(nout,*) 'grp     = ', z_grp(:,:,jz)
           end if
!
!
           call NCLASS(
     >              k_order, 
     >              k_potato, 
     >              m_i,                  ! Number of isotops
     >              i_maxcharge, 
     >              c_den, 
     >              c_potb, 
     >              c_potl, 
     >              z_b22di(jz), 
     >              z_bm22di(jz), 
     >              z_eb(jz), 
     >              z_fhat(jz), 
     >              z_fm(1:3, jz), 
     >              z_ft(jz), 
     >              z_grbm2(jz), 
     >              z_grphi(jz),
     >              z_gr2phi(jz), 
     >              z_ngrth(jz),
     >              amu_i, 
     >              z_grt(1:mx_mi, jz), 
     >              z_t(1:mx_mi, jz), 
     >              z_den(1:mx_mi, 1:mx_mz, jz),
     >              z_fex(1:3, 1:mx_mi, 1:mx_mz, jz), 
     >              z_grp(1:mx_mi, 1:mx_mz, jz),
     >              m_s,                  ! Number of species
     <              im_s,                 ! isotope number of s (-) 
     <              iz_s,                 ! charge state of s (-)  
     <              p_bsjb,               ! <J_bs.B> (A*T/m2)                           
     <              p_etap,               ! parallel electrical resistivity (Ohm*m)     
     <              p_exjb,               ! <J_ex.B> current response to fex_iz (A*T/m2)
     <              calm_i,               ! tp eff friction matrix for i (kg/m3/s)               
     <              caln_ii,              ! fp eff friction matrix for i1 on i2 (kg/m3/s)   
     <              capm_ii,              ! test part (tp) friction matrix for i1 on i2 (-) 
     <              capn_ii,              ! field part (fp) friction matrix for i1 on i2 (-)
     <              bsjbp_s,              ! <J_bs.B> driven by unit p'/p of s (A*T*rho/m3)
     <              bsjbt_s,              ! <J_bs.B> driven by unit T'/T of s (A*T*rho/m3)
     <              chi_s,                ! heat conduction coefficient (diag comp) of s (rho2/s)
     <              dn_s,                 ! diffusion coefficient (diag comp) of s (rho2/s) 
     <              gfl_s,                ! radial particle flux comps of s (rho/m3/s)      
     <              qfl_s,                ! radial heat flux comps of s (rho/m3/s)          
     <              sqz_s,                ! orbit ssqueezing factor for s (-)               
     <              upar_s,               ! parallel flow of s from force m (T*m/s)
     <              utheta_s,             ! poloidal flow of s from force m (m/s/T)
     <              vnnt_s,               ! convection velocity (off diag comps-p', T') of s (rho/s)
     <              vneb_s,               ! <E.B> particle convection velocity of s (rho/s)
     <              vnex_s,               ! extrenal force particle convection velocity of s (rho/s)
     <              vqnt_s,               ! conduction velocity (off diag comps-p',T') of s (rho2/s)
     <              vqeb_s,               ! <E.B> heat convection velocity of s (rho/s)
     <              vqex_s,               ! extrenal force heat convection velocity of s (rho/s)
     <              xi_s,                 ! charge weighted density factor of s (-)    
     <              ymu_s,
     <              chip_ss,
     <              chit_ss,
     <              dp_ss,
     <              dt_ss,
     <              iflag, 
     <              message )
!
!Check warning/error flags
!
           if (lneocl(4) > 0) then
             if (lastWarning/=iflag) then 	
               write(nout,'(1A)') message
               if (iflag > 0) stop
             end if  
             lastWarning = iflag	
           end if
!
! diagnostic printout
!
           if (lneocl(4) == 2) then
             write(nout,*) ' m_s  = ', m_s
             write(nout,*) ' im_s = ', im_s
             write(nout,*) ' iz_s = ', iz_s
           end if
!
! Compute <J_bs.B>
!           if (present(z_jbb)) then
!             do im = 2, mhyd + 1
!               z_jbb(jz)   = z_jbb(jz) + (bsjbp_s(im) + bsjbt_s(im)) /
!     >                       avi(jz, 6, 1)        
!             enddo  
!           endif
!
!|\subsection{Neoclassical effective electron/ion particle diffusivities}
!
! Note that in the BALDUR code:
! cfutz(9)  = coeff of hydroginic neoclassical transport
! cfutz(10) = coeff of impurity neoclassical transport
!
           do im = 2, m_s
             lm    = im_s(im)           ! izotop number
             lz    = abs(iz_s(im))      ! charge state
             dnneo1(lm-1, lm-1, jz) = dnneo1(lm-1,lm-1, jz)
     &          + cfutz(9) * dn_s(im)  * uisl**2
! diagonal particle diffusion coefficient
!
             vnneo1(lm-1, jz) = vnneo1(lm-1, jz) - (vnnt_s(im) 
     >                        + vneb_s(im)) * uisl * cfutz(19)
!             z_z    = z_grd(lm, lz, jz) 
!             z_z    = z_den(lm, lz, jz)         
!             z_zflp =  (dnneo1(lm-1,lm-1, jz)*z_grd(lm, lz, jz) *
!     >                 avi(jz, 5, 1) * usil-vnneo1(lm-1, jz)*        !!!
!     >                 z_den(lm, lz, jz))*usil
             do jm = 2, m_s
               km    = im_s(jm)         ! izotop number
               kz    = abs(iz_s(jm))    ! charge state
               if (int(amu_i(lm)) <= 2) then
                 if (int(amu_i(km)) <= 2) then
                   dnhhs(km-1, lm-1, jz) = dnhhs(km-1, lm-1, jz) + 
     >                cfutz(9) * ( dp_ss(jm, im) * 
     >                z_grp(km, kz, jz) / z_den(km, kz, jz) + 
     >                cfutz(9) * dt_ss(jm, im) * z_grt(km, jz)) 
     >                 / (z_t(km, jz) * z_den(lm, lz, jz) * 
     >                sign(max(abs(z_grd(km, kz, jz)), epslon), 
     >                z_grd(km, kz, jz))) 
     >                * uisl**5 
                 else
	           dnhis(km-mhyd-1, lm-1, jz) = 
     >                dnhis(km-mhyd-1, lm-1, jz)
     >                + cfutz(9) * ( dp_ss(jm, im) * 
     >                z_grp(km, kz, jz) / z_den(km, kz, jz) + 
     >                cfutz(9) * dt_ss(jm, im) * z_grt(km, jz)) 
     >                / (z_t(km, jz) * z_den(lm, lz, jz) *
     >                sign(max(abs(z_grd(km, kz, jz)), epslon), 
     >                z_grd(km, kz, jz)))
     >                * uisl**5
                 endif 
               else
                 if (int(amu_i(km)) <= 2) then
                   dnihs(km-1, lm-mhyd-1, jz) = 
     >                dnihs(km-1, lm-mhyd-1, jz) 
     >                + cfutz(10) * (dp_ss(jm, im) * 
     >                z_grp(km, kz, jz) / z_den(km, kz, jz) + 
     >                cfutz(10) * dt_ss(jm, im) * z_grt(km, jz))
     >                / (z_t(km, jz) * z_den(lm, lz, jz) * 
     >                sign(max(abs(z_grd(km, kz, jz)), epslon), 
     >                z_grd(km, kz, jz)))
     >                * uisl**5
                 else
                   dniis(km-mhyd-1, lm-mhyd-1, jz) = 
     >                dniis(km-mhyd-1, lm-mhyd-1, jz)
     >                + cfutz(10) * (dp_ss(jm, im) * 
     >                z_grp(km, kz, jz) / z_den(km, kz, jz) + 
     >                cfutz(10) * dt_ss(jm, im) * z_grt(km, jz))
     >                / (z_t(km, jz) * z_den(lm, lz, jz) * 
     >                sign(max(abs(z_grd(km, kz, jz)), epslon), 
     >                z_grd(km, kz, jz)))
     >                * uisl**5
                 endif 
               endif	
             enddo  
           enddo  
!
!|\subsection{Neoclassical effective electron/ion thermal diffusivities and convection velocities}
!|\begin{eqnarray}
!| \chi_i^{eff} = \sum_j{\chi_i^{(j)} n_j} \over \sum_j{n_j}
!|\end{eqnarray}
!|\begin{eqnarray}
!| v_{conv}^{(neo)} = \sum_j{\left(v_{q\ j}^{nT} + v_{q\ j}^{ex} + v_{q\ j}^{EB}\right) n_j} \over \sum_j{n_j}
!|\end{eqnarray}
!
         if (lneocl(4) == 2) then

             write(nout, '("chi_s   = ", 9(1pE16.5:))') 
     >             chi_s(1:m_s)
             write(nout, '("dn_s    = ", 9(1pE16.5:))') 
     >             dn_s(1:m_s)

         end if

           xeneo1(jz) = chi_s(1)
           veneo1(jz) = vqeb_s(1) + vqex_s(1) + vqnt_s(1)
           z_z        = 0.
           do im = 2, m_i
             lm         = im_s(im)                   ! izotop number
             lz         = abs(iz_s(im))              ! charge state
             z_z        = z_z + z_den(lm, lz, jz)
             xineo1(jz) = xineo1(jz) + z_den(lm, lz, jz) * chi_s(im)
             vineo1(jz) = vineo1(jz) + z_den(lm, lz, jz) * 
     >                    (vqeb_s(im) + vqex_s(im) + vqnt_s(im))
           end do
           xeneo1(jz)   = xeneo1(jz) * uisl**2 
           veneo1(jz)   = veneo1(jz) * uisl 
           xineo1(jz)   = xineo1(jz) * uisl**2 / ( z_z )
           vineo1(jz)   = vineo1(jz) * uisl / ( z_z )
!
           if (cfutz(43) > 0.) xeneo1(jz)  = min(xeneo1(jz), cfutz(43))
           if (cfutz(43) > 0.) xineo1(jz)  = min(xineo1(jz), cfutz(43))
!
             zxineo1(jz)=  zxineo1(jz) + xineo1(jz)
     >             - z_t(2, jz) * vineo1(jz)
     >           / ( sign(max(abs(z_grt(2, jz)), epslon), z_grt(2, jz)))
!
!| The effective thermal diffusivities (old definition):
!|\begin{eqnarray}
!| \chi_e = \frac{Q_e + 2.5 T_e \Gamma_e}{n_e \d T_e/\d \rho}
!| \chi_i = \sum_j{\frac{Q_i^{(j)} + 2.5 T_i \Gamma_i^{(j)}}{n_i \d T_i/\d \rho}}
!|\end{eqnarray}
!|
!|\subsection{Neoclassical transport fluxes}
!
           do im = 1, m_s
             lm    = im_s(im)                        ! izotop number
             lz    = abs(iz_s(im))                   ! charge state
             z_gfl = 0.
             z_qfl = 0.
             do iz = 1, 5
               z_gfl  = z_gfl + gfl_s(iz, im)
               z_qfl  = z_qfl + qfl_s(iz, im)
             end do  
             z_gfl           = z_gfl / avi(jz, 5, 1)
             z_qfl           = z_qfl / avi(jz, 5, 1)
             z_flxp1(im, jz) = z_gfl
             z_flxt1(im, jz) = z_qfl
             z_zflt          = 0.
             z_zflp          = 0.
             do jm = 1, m_s
               km    = im_s(jm)                    ! izotop number
               kz    = abs(iz_s(jm))               ! charge state
               z_zflp= z_zflp - 
     >                    (dp_ss(jm, im) * z_grp(km, kz, jz) /
     >                    z_den(km, kz, jz) + dt_ss(jm, im) * 
     >                    z_grt(km, jz)) / z_t(km, jz)
! alternative formulation for particle flux
               z_zflt= z_zflt - 
     >                    (chip_ss(jm, im) * z_grp(km, kz, jz) /
     >                    z_den(km, kz, jz) + chit_ss(jm, im) * 
     >                    z_grt(km, jz)) / z_t(km, jz) 
! alternative formulation for thermal flux 
             end do
             z_flxp2(jz)= z_flxp2(jz) + (z_zflp + vneb_s(im) + 
     >                    vnex_s(im)) * z_den(lm, lz, jz) /
     >                    avi(jz, 5, 1)
             z_flxt2(jz)= z_flxt2(jz) + (z_zflt + vqeb_s(im) +  
     >                    vqex_s(im)) * z_den(lm, lz, jz) * 
     >                    z_t(lm, jz) * z_j7kv / avi(jz, 5, 1)  
           end do
!
!|\subsection{Bootstrap current}
!
           ajboot(jz) = p_bsjb / avi(jz,9,1)
!
! end of loop over jz
!
         end do
!
! interpolate bootstrap current to zone centers
! Note that ajboot(jz) scales like sqrt(r) near the magnetic axis
! Hence, we interpolate ajboot(jz)*2
! This interpolation tries to preserve the sign of ajboot(jz)
!
         do jz=2,mzones-1
!
           zint = ( xbouni(jz+1) - xzoni(jz) )
     &          / ( xbouni(jz+1) - xbouni(jz) )
!
           if ( ajboot(jz) + ajboot(jz+1) .gt. 0.0 ) then
!
             ajboot(jz) =   sqrt ( zint * ajboot(jz)**2
     &                           + ( 1.0 - zint ) * ajboot(jz+1)**2 )
!
           else
!
             ajboot(jz) = - sqrt ( zint * ajboot(jz)**2
     &                           + ( 1.0 - zint ) * ajboot(jz+1)**2 )
!
           endif
!
         enddo
!
!|\subsection{Diagnostic output}
!
         if (lneocl(4) == 2) then
!
           write(nout, '(1A)') 'NCLASS diagnostic output'
!
!
           write(nout, '(A3, 3X, 8(9X, 1A8), 12(9X, 1A7, I1:))')  
     >                 ' jz', 'r', 'dTh.flux', 'dPr.flux', 'zxineo1',
     >                 'xeneo1', 'veneo1', 'xineo1', 'vineo1',  
     >                 ('dnneo1', jm, jm=1,m_s-1), 
     >                 ('vnneo1', jm, jm=1,m_s-1)
           do jz = 1, mzones 	
             write(nout, '(I3, 5X, 23(1pE16.5, 1X:))') jz, armins(jz,1),
     >                 z_flxt2(jz)-sum(z_flxt1(1:m_s,jz:jz)), 
     >                 z_flxp2(jz)-sum(z_flxp1(1:m_s,jz:jz)), 
     >                 zxineo1(jz)*uisl**2, xeneo1(jz), veneo1(jz), 
     >                 xineo1(jz), vineo1(jz), 
     >                 (dnneo1(jm,jm,jz), jm=1, m_s-1), 
     >                 (vnneo1(jm, jz), jm=1, m_s-1)
           end do
           write(nout, '(1A)')'NCLASS diagnostic output: thermal fluxes'
           write(nout, '(A3, 14X, 1A8, 8(8X, 1A7, I1:))')  ' jz', 'r', 
     >                       ('Th.flux', jm, jm=1,m_s) 
           do jz = 1, mzones 	
             write(nout, '(I3, 6X, 1pE16.5\)') jz, armins(jz,1)
             do jm = 1, m_s
               write(nout, '(1pE16.5\)') z_flxt1(jm, jz)
             enddo
             write(nout, *) ''
           end do
           write(nout,'(1A)')'NCLASS diagnostic output: particle fluxes'	
           write(nout, '(A3, 14X, 1A8, 9(8X, 1A7, I1:))')  ' jz', 'r', 
     >                       ('Pr.flux', jm, jm=1,m_s) 
           do jz = 1, mzones 	
             write(nout, '(I3, 6X, 1pE16.5\)') jz, armins(jz,1)
             do jm = 1, m_s
               write(nout, '(1pE16.5\)') z_flxp1(jm, jz)
             enddo
             write(nout, *) ''
           end do
         end if
!

        return
        end subroutine nclass_int


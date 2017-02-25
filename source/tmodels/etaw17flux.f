!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| 
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt
!| \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center}
!| {\bf {\tt etaw17flux.tex} \\
!| Toroidal Ion Temperature Gradient Mode \\
!| \vspace{1pc}
!| Glenn Bateman \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| Jan Weiland, Hans Nordman and P{\"a}r Strand\\
!| Department of Electromagnetics \\
!| Chalmers University of Technology \\
!| S-412 96 G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Jon Kinsey \\
!| General Atomics \\
!| P.O. Box 85608, San Diego, CA 92186} \\
!| \vspace{1pc}
!| \today
!| \end{center}
!| This subroutine evaluates the eigenvalue problem for $\eta_i$
!| and trapped electron modes derived by
!| Jan Weiland, H. Nordman and their group in G\"{o}teborg Sweden.
!| The equations in this routine include fast Hydrogenic ions,
!| impurities, trapped electron and finite Larmor radius effects.
!| New options include parallel ion motion, finite beta and collisional effects
!| together with an approximate treatment of $E\times B$ flow shear reduction of
!| transport.
!| 
!| The essential idea is to linearize the fluid equations
!| with magnetic drifts for a given Fourier harmonic
!| \[ u = u_0(x) + \tilde{u} \exp\left[ -i (\omega t
!|    + \vec{k} \cdot \vec{x} ) \right] \]
!| Compute the eigenvalues and eigenvectors from these equations and then
!| compute the quasi-linear heat and particle fluxes.
!| 
!| 
!| The fundamental equations used in this model are the fluid equations.
!| For ions (hydrogen or impurity), the equation of continuity is
!| \[ \Partial{n_i}{t} + \nabla \cdot ( n_i v_i ) = 0 \]
!| where
!| \[ v_i = v_E + v_{*i} + v_{Pi} + v_{\pi i}
!|    + \hat{B} v_{\parallel i} \]
!| are the fluid flows where
!| \[ v_E = E \times \hat{B} / B = \hat{B} \times \nabla \phi / B
!|  = - i k_y \phi \hat{x} / B \]
!| is the E cross B drift;
!| \[ v_{*i} = \frac{ \hat{B} \times \nabla ( n_i T_i )}{Z_i e n_i B}
!|  = i ( 1 + \eta_i ) \omega_{*i} / k_y \]
!| \[ \omega_{*i} = \frac{ - k_y T_i }{ Z_i e B L_{ni} } \]
!| is the diamagnetic drift ($ Z_i = 1 $ for hydrogen isotopes and
!| $ Z_i = -1 $ for electrons);
!| \[ L_{ni} \equiv - n_i / \hat{x} \cdot \nabla n_i \]
!| is the density gradient scale length
!| \[ \eta_i \equiv L_{ni} / L_{Ti}
!|  = \frac{n_i \hat{x} \cdot \nabla T_i}{T_i \hat{x} \cdot \nabla n_i}, \]
!| \[ v_{Pi} = \frac{ d E }{ dt } / ( B \Omega_i ) \]
!| is the polarization drift with
!| \[ \Omega_i = \frac{ Z_i e B }{ m_i}; \]
!| and
!| \[ v_{\pi i} = \frac{ \hat{B} \times \nabla \cdot \pi_i}{
!|   Z_i e n_i B } \]
!| is the drift due to the off-diagonal elements of the stress tensor.
!| 
!| The parallel ion motion $ v_{\parallel i} $ is determined by the parallel
!| ion momentum equation driven by electrostatic forces and ion pressure
!| gradient along the field lines\cite{weil92d}:
!| \[ m_i n_i \Partial{v_{\parallel i}}{t}
!|   = - e Z_i n_i \nabla_{\parallel} \phi
!|   - \nabla_{\parallel} ( n_i T_i ) \]
!| 
!| \noindent
!| The ion energy balance equation is
!| \[ \frac{3}{2} n_i \left( \Partial{ }{t} + \vec{v}_i \cdot \nabla \right) T_i
!|   + n_i T_i \nabla \cdot \vec{v}_i = - \nabla \cdot \vec{q}_{*i}
!|   = \frac{5}{2} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| where $ \vec{q}_{*i} $ is the diamagnetic ion heat flow and
!| \[ \vec{v}_{Di} = \frac{T_i}{m_i \Omega_i} \hat{B} \times
!|   \left( \frac{\nabla B}{B} + \vec{\kappa} \right) \]
!| is the drift due to $\nabla |B|$ and magnetic curvature
!| $ \vec{\kappa} = \hat{B} \cdot \nabla \hat{B} $.  Note
!| \[ \vec{k} \cdot \vec{v}_{Di} = \omega_{Di} = \frac{-2 k_y T_i}{Z_i e B R} \]
!| is the diamagnetic drift frequency.
!| 
!| The equations for trapped electron continuity and energy flow have
!| the same form except that the polarization drift and the drift
!| due to the stress tensor
!| can be neglected and the electron fluid velocity can be written
!| \[ v_e = v_E +  v_{*e}.\]
!| 
!| Let $f_t$ be the fraction of trapped electrons, $ f_t = n_{et} / n_e $,
!| where $n_{et}$ is the density of trapped electrons and
!| $n_{ef}$ and the density of free (circulating) electrons.  Then
!| \[ n_e = n_{et} + n_{ef}. \]
!| The perturbed density can be written
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_{et}}{n_e} + \frac{\tilde{n}_{ef}}{n_e}
!|  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!|        + ( 1 - f_t ) \frac{\tilde{n}_{ef}}{n_{ef}} \]
!| %  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!| %        + ( 1 - f_t ) \frac{ e \tilde{\phi} }{T_e}. \]
!| 
!| Let $ f_Z $ be the fraction of impurity ions with charge $Z$
!| relative to the electron density $ f_Z = n_Z / n_e $.
!| Charge neutrality gives
!| \[ n_e = n_H + Z n_Z + n_s \]
!| where $n_H$ is the density of the hydrogen isotopes and
!| $n_Z$ is the density of the impurity species and $n_s$ is the density of
!| superthermal hydrogenic ions.  Then
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_H}{n_e} + \frac{Z \tilde{n}_Z}{n_e}
!|  =  ( 1 - Z f_Z - f_s ) \frac{\tilde{n}_H}{n_H}
!|        + Z f_Z \frac{\tilde{n}_Z}{n_Z}  \]
!| where $ f_s \equiv n_s / n_e $.
!| For now, it is assumed that the superthermal ions
!| do not take part in the perturbation $ \tilde{n}_s = 0.$
!| 
!| A relation between the perturbed free (circulating) electron density
!| $ n_{ef} $ and the perturbed electric and magnetic potentials
!| can be obtained from the momentum equation for free electrons parallel
!| to the unperturbed magnetic field
!| \[ m_e n_{ef} \Partial{\vec{v}_{\parallel e}}{t}
!|  = -e n_{ef} ( \vec{E}_\parallel + \frac{1}{c} \vec{v} \times \vec{B} )
!|  - \nabla_\parallel p_{ef} \]
!| where
!| \[ \vec{E}_\parallel = - \nabla \phi - \frac{1}{c} \Partial{\vec{A}_\parallel}{t} \]
!| and the unperturbed perpendicular free electron velocity is the
!| diamagnetic velocity
!| \[ \vec{v}_{*ep} = - \frac{\hat{B} \times \nabla ( n_{ef} T_e )}{
!|    e n_{ef} B}. \]
!| Assume $ k_\parallel v_{the} \gg \omega $ to neglect the inertial term
!| on the left of this momentum equation.
!| 
!| Now use the fact that the electron temperature is nearly uniform along
!| the perturbed magnetic field $ B \cdot \nabla T_e \approx 0 $
!| (this follows from the free electron energy equation), by the adiabatic
!| condition
!| \[ \frac{\tilde{n}_{ef}}{n_{ef}} = \frac{e \tilde{\phi}}{T_e}. \]
!| 
!| Justification for the use of all of the above equations
!| and more complete derivations
!| can be found in the references.\cite{weil92a,nord90a}
!| 
!| Note that
!| \[ n_i \nabla \cdot \vec{v}_i = - \Partial{n_i}{t} - \vec{v}_i \cdot \vec{n}_i \]
!| from the equation of continuity can be used to eliminate
!| $ \nabla \cdot \vec{v}_i $ in the energy equation (divided by 3/2):
!| \[ n_i \left( \Partial{}{t} + \vec{v}_i \cdot \nabla \right) T_i
!| -\frac{2}{3} T_i \Partial{n_i}{t}
!|  - \frac{2}{3} T_i \vec{v}_i \cdot \nabla n_i
!|  = \frac{5}{3} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| and a similar equation for the trapped electrons.
!| The $\vec{v}_{*i}$ (and corresponding $\vec{v}_{*e}$)
!| contributions drop out because
!| \[ n_i \vec{v}_{*i} \cdot \nabla T_i
!|  - \frac{2}{3} T_i \vec{v}_{*i} \cdot \nabla n_i
!|  - \frac{5}{3} n_i \vec{v}_{*i} \cdot \nabla T_i
!|  = \frac{2}{3} \vec{v}_{*i} \cdot \nabla ( n_i T_i ) = 0 \]
!| since the diamagnetic drift
!| $ \vec{v}_{*i} = \hat{B} \times \nabla ( n_i T_i ) / ZeB n_i $
!| is perpendicular to the gradient of the pressure
!| $ \nabla ( n_i T_i ) $.
!| 
!| See the derivations in the notes by Weiland\cite{weil92a} for
!| the relations
!| \[ \nabla \cdot \delta ( {n_i v_{*i}} )
!|  = \vec{v}_{Di} \cdot \nabla \delta ( {n_i T_i} ) / T_i \]
!| (page 132 Eq. (4.20)),
!| \[ \nabla \cdot \tilde{v}_E
!|    = \frac{Z e}{T_i} \vec{v}_{Di} \cdot \nabla \tilde{\phi} \]
!| (page 132 Eq. (4.21)), and
!| \[ \nabla \cdot [ n_i ( \vec{v}_{pi} + \vec{v}_{\pi i} ) ]
!|  \approx - i n_i k_y^2 \rho_{si}^2
!|  [ \omega - \omega_{*i} ( 1 + \eta_i ) ] \frac{e \tilde{\phi}}{T_e} \]
!| (page 25 Eq. (1.28)) where
!| \[ \rho_{si}^2 = \frac{T_e}{m_i \Omega_i^2} \]
!| with $ v_{thi}^2 = 2 T_i / m_i $ and $ \Omega_i = Z e B / m_i $.
!| 
!| Using these relations in the continuity, momentum, and energy equations,
!| we obtain the following perturbed relations,
!| \[ (-\omega + \omega_{Di}) \hat{n}_i + \omega_{Di} \hat{T}_i
!|  + [ ( \omega_{Di} - \omega_{*i} ) Z_i T_e/T_i
!|  - k_y^2 \rho_{si}^2 ( \omega - \omega_{*i} ( 1 + \eta_i ) )
!|  ] \hat{\phi} + k_\parallel v_{\parallel i} = 0 \]
!| \[ - \omega m_i v_{\parallel i} + k_\parallel Z_i T_e \hat{\phi}
!|  + k_\parallel T_i ( \hat{n}_i + \hat{T}_i ) = 0 \]
!| \[ (-\omega + \frac{5}{3} \omega_{Di} ) \hat{T}_i
!|    + \frac{2}{3} \omega \hat{n}_i
!|    + \omega_{*e} \frac{L_{ne}}{L_{ni}}
!|    ( \eta_i - \frac{2}{3} ) \hat{\phi} = 0 \]
!| with corresponding equations for trapped electrons ($i=et$)
!| and impurities ($i=Z$).
!| Here $\hat{n} \equiv \tilde{n} / n$,
!| $\hat{T} \equiv \tilde{T} / T$, and
!| $\hat{\phi} \equiv e \tilde{\phi} / T_e$ are dimensionless forms
!| of the perturbation.
!| 
!| Normalizing all frequencies by
!| $\omega_{De} = \vec{k} \cdot \vec{v}_{De} = 2 k_y T_e /  e B R $,
!| and normalizing the parallel ion velocity by the speed of sound
!| in that ion,
!| $ \hat{v}_{\parallel i} \equiv v_{\parallel i} / c_{si} $
!| where $ c_{si} \equiv \sqrt{T_e / m_i} $,
!| we obtain
!| \[ ( \hat{\omega} + \frac{T_i}{Z_i T_e} ) \hat{n}_i
!|   + \frac{T_i}{Z_i T_e} \hat{T}_i
!|   + \left[ 1 - \frac{g_{ni}}{2}
!|   + k_y^2 \rho_{si}^2 \left( \hat{\omega}
!|     + \frac{T_i}{Z_i T_e} \frac{g_{ni}+g_{Ti}}{2}
!|     \right) \right] \hat{\phi}
!|   - \frac{k_\parallel c_{si}}{\omega_{De}} \hat{v}_{\parallel i}
!|   = 0 \]
!| \[ \hat{\omega} \hat{v}_{\parallel i}
!|  - \frac{k_\parallel c_{si}}{\omega_{De}} Z_i \left[ \hat{\phi}
!|    + \frac{T_i}{Z_i T_e} ( \hat{n}_i + \hat{T}_i ) \right] = 0 \]
!| \[ \left( - \hat{\omega} - \frac{5}{3} \frac{T_i}{Z_i T_e}
!|    \right) \hat{T}_i + \frac{2}{3} \hat{\omega} \hat{n}_i
!|    + \frac{1}{2} \left( g_{Ti} - \frac{2}{3} g_{ni}
!|      \right) \hat{\phi} = 0 \]
!| where $ \hat{\omega} \equiv \omega / \omega_{De} $
!| and
!| \[ g_{ni} \equiv R / L_{ni} = - R \hat{x} \cdot \nabla n_i / n_i \]
!| is the normalized gradient.
!| 
!| Now use
!| \[ f_t \hat{n}_{et} = ( 1 - Z f_Z - f_s ) \hat{n}_H
!|    + Z f_Z \hat{n}_Z - ( 1 - f_t ) \hat{\phi} \]
!| to eliminate the perturbed trapped electron density $\hat{n}_{et}$
!| in favor of the perturbed ion densities and perturbed potential.
!| The perturbed hydrogen and impurity equations remain as above
!| while the trapped electron density and energy equations become:
!| \[ (1-\hat{\omega}) ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + (1-\hat{\omega}) Z f_Z \hat{n}_Z  + f_t \hat{T}_{et}
!|  - [ 1 - f_t g_{ne} / 2 - ( 1 - f_t ) \hat{\omega} ] \hat{\phi}
!|  = 0 \]
!| \[ \frac{2}{3} \hat{\omega} ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + \frac{2}{3} \hat{\omega} Z f_Z \hat{n}_H
!|  + \left( \frac{5}{3} - \hat{\omega} \right) f_t \hat{T}_{et}
!|  + \left[ f_t \frac{g_{Te} - (2/3) g_{ne}}{2}
!|    - \frac{2}{3} \hat{\omega} ( 1 - f_t ) \right] \hat{\phi}
!|  = 0  \]
!| 
!| In all of these expressions we use the notation
!| $$ \omega_{*e} = k_y T_e / e B L_n $$
!| $$ \omega_{*i} = - k_y T_i / Z_i e B L_n
!|     = - T_i \omega_{*e} / ( Z_i T_e ) $$
!| $$ \omega_{De} = 2 k_y T_e / e B R = 2 L_n \omega_{*e} / R $$
!| $$ \omega_{Di} - 2 k_y T_i / Z_i e B R = 2 L_n \omega_{*i} / R
!|    = - 2 L_n T_i \omega_{*e} / T_e R $$
!| $$ g_{ni} = - R \hat{x} \cdot \nabla n_i / n_i $$
!| $$ \eta_i = L_{ni} / L_{T_i} = g_{Ti} / g_{ni}  $$
!| $$ \eta_e = L_{ne} / L_{T_e} = g_{Te} / g_{ne} $$
!| and
!| $$ f_t \approx \sqrt{ \frac{ 2 r/R }{ 1 + r/R } } $$
!| is the trapped electron fraction.
!| Note that here I define
!| \[ \epsilon_{ni} \equiv L_{ni} / R \]
!| while $ \epsilon_{ni} $ in the Weiland papers is defined to be
!| twice this value.
!| 
!| The basic technique is to set up the generalized eigenvalue problem
!| \[ A v = \lambda B v \]
!| where $ \lambda = \hat{\omega} + i \hat{\gamma} $ is the eigenvalue
!| and $ v $ is the corresponding eigenvector.
!| Hence the eigenvalues give the frequency and growth rates of the
!| modes while the eigenvectors give the phase of the perturbed
!| variables relative to one another.
!| 
!| In order to approximate the reduction of transport with $E\times B$
!| velocity shear the $E\times B $ shearing rate is subtracted from the
!| linear growth rates obtained as above. Currently this is implemented as
!| follows: The generalized eigenvalue problem is redefined using
!|  \[ (A - {\it i} \omega_{E\times B} B) v' = {\lambda}' B v' \]
!| and fluxes are calculated rom the new growth rates   ${\lambda}'$
!| and eigenvectors $v'$. This method do give the same results as a
!| direct reduction of the growth rates only would give but is more
!| integrated with the current framework of the model.
!| 
!| 
!| In this routine, the perturbed variables are always computed in
!| the following order:
!| \[ \hat{\phi}, \hat{T}_H, \hat{n}_H, \hat{T}_{et},
!|     \hat{n}_Z, \hat{T}_Z, \ldots \]
!| 
!| 
c@etaw17flux.f
c rap 31-jul-02 dimag function is replaced by aimag
c pis 27-may-98 replaced all calls to abortb with nerr = 2; return, statements
c pis 27-may-98 added new error message for gradients that break neautrality
c pis  7-may-98 split the routine etaw17a into two routines; this routine
c               etaw17flux and etaw17diff which defines the transport martix
c               rather than fluxes
c pis  7-may-98 removed letain(15) = 1 choice of subtracting from growth rates
c pis  5-may-98 removed some uneccessary loops and unused variables.
c pis  5-may-98 relabelled routine etaw17flux.tex and introduced tomsqz
c               instead of tomslz. QZ algorithm more stable and exact than LZ.
c pis 29-apr-98 added wexb shearing rate to argument list, currently this
c               is implemented in two different ways defined by the value of
c                   letain(15) = 0 (default) matrices redifened
c                   letain(15) = 1 shearing rate subtracted from growth rates
c                                  (this implementationwill produce err mess.)
c pis 28-apr-98 replaced nag library routines with tomslz (wrapper for lzhes
c               and lzit)
c jek 30-nov-97 added modified FLR parameter to include elongation
c               coeff in Alfven frequency is cetain(25) (default=0.0)
c rgb 17-may-97 added if ( zfns .lt. zepsmach) zgns = 0
c   multiplied zamr(11,9) by zimp
c jek 12-apr-96 added cross term to electromagnetic version with
c               parallel ion motion
c jek 25-mar-96 added parallel ion motion in strong ballooning limit
c               to nine equation model (H=S/2q). cetain(12) (default=0.0)
c               added to turn on/off
c rgb 13-feb-96 input gnein=-R(d n_e / d r)/n_e rather than gnsin
c   rearranged order of argument list to put gnein before gnhin
c rgb 12-feb-96 zgp* = zgp* + zgn*  -->  zgp* = zgt* + zgn*
c   (this should have no effect on the results)
c   corrected header for eleven eqns
c rgb 08-jan-96 diagnostic printout of frequencies, fluxes, and phases
c rgb 03-jul-95 computed zerrmax and used it for error control
c   implemented control of imatrx by letain(2)
c rgb 02-jul-95 added impurity parallel ion motion 11 equations
c rgb 01-jul-95 when computing fluxes, loop over all eigenvalues ieq
c rgb 29-jun-95 protected zalp, zalf, and kps when zgnh+zgth.lt.zepsqrt
c   corrected zamr(3,9), zamr(3,10), zbmr(3,10)
c rgb 28-jun-95 cleaned up complex matrix option
c rgb 24-jun-95 converted from eps*in to normalized gradients g*in
c   use zgnh instead of zgne in definition of k1 and k2
c   removed vef = 0.10
c   zamr(7,1) = zhalf*zgte - zone  inserted in the 7 eqn model
c   replaced 1./kpc with zanorm = cetain(11)
c   added the following equation to the 10 eqn set
c   zamr(9,8) = kps * ( zone - zfnz - zfns ) / ( zone - zft )
c   Pass nout through the argument list just after neq
c 04-apr-95 added vef, ion motion, finite beta
c 04-jan-95 added vef to argument list
c 06-dec-94 include disp9 from Weiland
c
c  THIS ROUTINE IS A MODIFICATION OF THE LINEAR PART OF ETAWN6
c  WRITTEN BY GLENN BATEMAN. INCLUDING COLLISIONS TO THE
c  TRAPPED ELECTRONS YIELDS A 7 EQUATION SYSTEM. WHEN PARALLEL
c  ION MOTION IS INCLUDED THERE ARE 8 EQUATIONS. THE 9 EQUATION
c  SYSTEM INCLUDES BOTH COLLISIONS AND PARALLEL ION MOTION.
c  WITH PARALLEL ION MOTION IN THE STRONG BALLOONING LIMIT.
c  THE 9 EQUATION SYSTEM IS MODIFIED TO INCLUDE SHEAR EFFECTS.
c  THE NEW ELECTROMAGNETIC SYSTEM INCORPORATES THE PARALLEL
c  VECTOR POTENTIAL AND ITS TIME DERIVATIVE AS NEW VARIABLES.
c  THIS LEADS TO A SYSTEM OF 10 EQUATIONS WHEN COLLISIONS AND
c  PARALLEL ION MOTION ARE INCLUDED. WITH PARALLEL ION MOTION
c  OF IMPURITIES, THE SYSTEM IS EXTENDED TO 11 EQUATIONS.
c  MOST PARAMETERS ARE TRANSFERRED THROUGH COMMON BLOCKS LIKE GRAD,
c  IMP AND ETAWN6. NOTE THE INVERSE DEFINITION OF ZTAUH AND ZTAUZ !
c
c   The variables are ordered as: e phi/Te, Th, Nh, Te, Nz, Tz, F, Vp,
c   Av, K where F is due to collisions, Vp is due to parallel ion motion,
c   Av is the vector potential and K is its time derivative.
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine etaw17flux (letain, cetain, lprintin, neq, nout
     & , gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chieff, fluxout
     & , nmodes, perform, nerr )
c
c ... the table of argument list variable names and definitions follows
c     (in latex file)
!| \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Names of variables in the argument list:
!| \begin{tabular}{lllp{3.0in}}
!| variable & status & symbol & meaning \\
!| {\tt letain(j)} & input & & Integer control variables
!|                             (see table below).\\
!| {\tt cetain(j)} & input & & Real-valued control variables
!|                             (see table below).\\
!| {\tt lprint}    & input & & Controls printout.
!|  Higher values produce more printout. \\
!| {\tt neq} & input & & number of equations \\
!| {\tt gnein} & input & $  $
!|     & $ - \hat{r} \cdot \nabla n_e / n_e $ \\
!| {\tt gnhin} & input & $  $
!|     & $ - \hat{r} \cdot \nabla n_H / n_H $ \\
!| {\tt gnzin} & input & $  $
!|     & $ - \hat{r} \cdot \nabla n_Z / n_Z $ \\
!| {\tt gtein} & input & $  $
!|     & $ - \hat{r} \cdot \nabla T_e / T_e $ \\
!| {\tt gthin} & input & $  $
!|     & $ - \hat{r} \cdot \nabla T_H / T_H $ \\
!| {\tt gtzin} & input & $  $
!|     & $ - \hat{r} \cdot \nabla T_Z / T_Z $ \\
!| {\tt tauhin} & input & $\tau_H$ & $ \tau_H = T_H / T_e $ \\
!| {\tt tauzin} & input & $\tau_Z$ & $ \tau_Z = T_Z / T_e $ \\
!| {\tt fnzin} & input & $ f_{nZ} $ & $ = n_Z / n_e $ \\
!| {\tt czin}  & input & $ Z $ & impurity charge number \\
!| {\tt azin}  & input & $ m_Z / m_H $
!|      & impurity mass to hydrogen isotope mass. \\
!|  & & & Note that ``hydrogen'' may include a deuterium or tritium mix. \\
!| {\tt fnsin}  & input & $ f_s $
!|   & $ f_s = n_s / n_e $ fraction of superthermal hydrogenic ions \\
!| {\tt betaein} & input & $\beta_e$ &
!|      $ = n_e T_e / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt betahin} & input & $\beta_H$ &
!|      $ = n_H T_H / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt betazin} & input & $\beta_Z$ &
!|      $ = n_Z T_Z / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt ftrapein}  & input & $f_{trap} $ &
!|     fraction of trapped electrons \\
!| {\tt vef}       & input & $ \nu_{th} / \omega_{De} $ &
!|      thermal collision frequency, normalized \\
!| {\tt q}         & input & $ q $ & magnetic q-value \\
!| {\tt shear}     & input & $ s $ & $ d \ln q / d \ln r $ \\
!| {\tt ekyrhoin}  & input & $ k_y \rho_s $ & normalized poloidal
!|                     wave number \\
!| {\tt ekparlin}  & input & $k_\parallel L_n$ & \\
!| {\tt wexb}      & input & $\omega_{E\times B}$& $E\times B$ shearing rate
!| (normalized with $\omega_{D_e}$) \\
!| {\tt ndim} & input & & first dimension of the 2-D array difthi
!|                and the maximum number of unstable modes allowed \\
!| {\tt omega(j)}  & output & $\omega / \omega_{De} $ &
!|      real part of the frequencies normalized by $ \omega_{De} $ \\
!| {\tt gamma(j)}  & output & $\gamma / \omega_{De} $ &
!|      growth rates normalized by $ \omega_{De} $ \\
!| {\tt difthi(i,j)}      & output & $ D \omega_{De} / k_y^2 $
!|       & diffusivity matrix normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt velthi(i,j)}      & output & $ v \omega_{De} / k_y $
!|       & convective velocities normalized by $ k_y / \omega_{De} $ \\
!| {\tt chieff(j)} & output & $ \chi_{\rm eff} \omega_{De} / k_y^2 $
!|       & effective total diffusivities
!|         for $ n_H T_H $, $ n_H $, $ n_e T_e $,
!|         $ n_Z $, $ n_Z T_Z $, \ldots
!|         normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt nmodes} & output & & number of unstable modes \\
!| 
!| \end{tabular}
!| \end{center}
!| 
!| \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Integer control variables in the argument list:
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| {\tt letain(2)} & = number of elements computed for transport matrix
!|                     (only when $> 0$) \\
!| {\tt letain(15)} & $ = 0 $ to redefine $A$ matrix with $\omega_{E\times B}$ \\
!|                  & $ = 1 $ to subtract shearing rate from growth rates directly
!| \\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| 
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Real-valued control variables in the argument list:
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| {\tt cetain(10)} & = coefficient of $ k_\parallel $ (default to 1.0) \\
!| {\tt cetain(11)} & = normalization for $ A_\parallel $ (default to 1.0) \\
!| {\tt cetain(12)} & = coefficient of $H$ (default to 0.0) \\
!| {\tt cetain(15)} & = coefficient of $ \hat{\nu} $ (default to 1.0) \\
!| {\tt cetain(20)} & = coefficient of $ \beta_{e,h,z} $ (default to 1.0) \\
!| {\tt cetain(25)} & = coefficient for elongation modified Alfven frequency (defau
!| lt to 0.0) \\
!| {\tt cetain(29)} & = radius used in printouts \\
!| {\tt cetain(30)} & = finite difference used to construct
!|                    transport matrix (see zgm(j1,jd) matrix) \\
!| {\tt cetain(32)} & = tolerance used in eigenvalue solver \\
!| 
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Effects included with different numbers of equations:
!| \begin{tabular}{lp{5.0in}}
!| neq & effects \\
!| 2 & Hydrogenic Ion Temperature Gradient mode only \\
!| 4 & Hydrogenic ITG + trapped electron modes \\
!| 5 & Hydrogenic ITG with parallel ion motion + TEM \\
!| 6 & Hydrogenic + impurity ITG + TEM \\
!| 7 & Hydrogenic + impurity ITG + TEM + collisions \\
!| 8 & Hydrogenic ITG with parallel hydrogen ion motion
!|   + impurity ITG (without parallel impurity ion motion)
!|   + TEM + collisions\\
!| 9 & Hydrogenic ITG + ipurity ITG (without any parallel ion motion)
!|   + TEM + electromagnetic (finite $\beta$) effects \\
!| 10 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG without parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!| 11 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG with parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!| 
!| \end{tabular}
!| \end{center}
!| 
!| 
!| Note, $ \omega_{*e} = k_\perp \rho_s c_s / L_n $,
!| $ \rho_s = c_s / \omega_{ci} $, $ c_s = \sqrt{ 2 T_e / m_i} $,
!| and $ \omega_{ci} = Z_i e B / m_i $,
!| $ \omega_{De} R = k_y \rho_s c_s $,
!| $v_i = \sqrt{T_i / m_i}$, and $ L_n = - n / \Partial{n}{r} $
!| in the above normalizations.
!| 
!| The diffusivity matrix $ D = {\tt difthi(j1,j2)}$
!| and convective velocity arrary $ v = {\tt velthi(j1)} $
!| are given in the following form:
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots
!|     \end{array} \right)
!|  = \nabla \cdot
!| \left( \begin{array}{llll}
!| D_{1,1} n_H & D_{1,2} T_H & D_{1,3} n_H T_H / T_e \\
!| D_{2,1} n_H / T_H & D_{2,2} & D_{2,3} n_H / T_e \\
!| D_{3,1} n_e T_e / T_H & D_{3,2} n_e T_e / n_H & D_{3,3} n_e & \vdots \\
!| D_{4,1} n_Z / T_H & D_{4,2} n_Z / n_H & D_{4,3} n_Z / T_e \\
!| D_{5,1} n_Z T_Z / T_H & D_{5,2} n_Z T_Z / n_H &
!|         D_{5,3} n_Z T_Z / T_e \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!|  \nabla
!|  \left( \begin{array}{c}  T_H \\ n_H \\  T_e \\
!|    n_Z \\  T_Z \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  + \nabla \cdot
!| \left( \begin{array}{l} {\bf v}_1 n_H T_H \\ {\bf v}_2 n_H \\
!|    {\bf v}_3 n_e T_e \\
!|    {\bf v}_4 n_Z \\ {\bf v}_5 n_Z T_Z \\ \vdots \end{array} \right) +
!|  \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right) $$
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $,
!| convective velocities are normalized by $ \omega_{De} / k_y $,
!| and all the frequencies are normalized by $ \omega_{De} $.
c
c  gnein  = - R ( d n_e / d r ) / n_e
c  gnhin  = - R ( d n_H / d r ) / n_H
c  gnzin  = - R ( d n_Z / d r ) / n_Z
c  gtein  = - R ( d T_e / d r ) / T_e
c  gthin  = - R ( d T_H / d r ) / T_H
c  gtzin  = - R ( d T_Z / d r ) / T_Z
c  tauhin = T_H / T_e
c  tauzin = T_Z / T_e
c
c  This version of etaw17flux is intended for use on workstations
c
c  Note that external subroutines need to be provided:
c
c  tomsqz   Code wrapper for ACM/TOMS routine 535 implementing the QZ algorithm
c           complex generailized eigenvalue problem ( eigenvalues and
c           eigenvectors )
c           The algorithm itself consists of the three routines
c           cqzhes, cqzval and cqzvec
c
c  Compile this routine  and routines that call it with a compiler option
c  such as -r8  to convert real to double precision when used on workstations.
c
      implicit none
c
      integer idp
      parameter ( idp = 15 )
c
      logical inital
      data inital /.true./
c
      dimension letain(32), cetain(32)
     &  , omega(*), gamma(*), fluxout(*)
     &  , chieff(*), perform(*)
c
      real cetain
     & , gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, H, ekyrhoin, ekparlin, wexb
     & , omega, gamma, chieff, fluxout, perform
c
      integer letain, lprintin, lprint, neq, ndim, nmodes, nerr, nout
     & , ieq, idim, imatrx, j1, j2, j
c
c ndim  = first dimension of the 2-D array difthi
c         and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c ieq   = number of equations
c
c imatrx = the number of elements computed along each row and
c column of the transport matrix
c
c  Note:  Normally the transport matrix is  ieq-1 by ieq-1
c    however, if there are 6 equations, compute only a 4 by 4 matrix
c    since the impurity temperature equation is not used by most
c    transport codes including BALDUR
c
      real zamr(idp,idp), zami(idp,idp), zbmr(idp,idp), zbmi(idp,idp)
      real zamrt(idp,idp),zamit(idp,idp),zbmrt(idp,idp),zbmit(idp,idp)
      real zvr(idp,idp), zvi(idp,idp), zomega(idp), zgamma(idp)
      real zalfr(idp), zalfi(idp), zbeta(idp)
c
      integer ifail
c
c ( zamr(i,j), zami(i,j) ) = matrix A
c ( zbmr(i,j), zbmi(i,j) ) = matrix B
c
c Note that the eigenvalues are
c
c   zomega(j) = zalfr(j) / zbeta(j)
c   zgamma(j) = zalfi(j) / zbeta(j)
c
c and the eigenvectors are
c
c   zevec(j) = cmplx( zvr(j), zvi(j) )
c
c Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
c
      complex  zevec(idp,idp)
c
      real  zepsmach, zepsqrt
     & , zone, ztwo, zthree, zfour, zfive, zhalf, zquarter, ztvr, zftr
     & , ztwohlf
     & , zgne, zgnh, zgnz, zgns, zgte, zgth, zgtz
     & , ztauh, ztauz, zft
     & , zimp, zfnz, zmass, zfns, zflh, zflz, zgamax
     & , zetae, zetah
c
c zepsmach = machine epsilon
c zepsqrt  = sqrt ( machine epsilon )
c zone     = 1.0
c ztwo     = 2.0
c zthree   = 3.0
c zfour    = 4.0
c zfive    = 5.0
c zhalf    = 0.5
c zquarter = 0.25
c ztwohalf = 2.5
c ztvr     = 2. / 3.
c zftr     = 5. / 3.
c
      real zflxph, zflxnh, zflxpe, zflxnz, zflxpz
     &  ,  zphsph, zphsnh, zphspe, zphsnz, zphspz
     &  ,  ztemp1, zreal, zimag
     &  ,  ztempa(idp), ztempb(idp), zerreal, zerimag, zerrmax
c
c  These are the thermal and particle fluxes and effective diffusivities
c
      real zflxm(idp), zchim(idp)
c
c  Normalized transport fluxes         Eff. Diffusivities
c       zflxm(1)   n_H T_H             zchim(1)    chi_i
c       zflxm(2)   n_H                 zchim(2)    D_i
c       zflxm(3)   n_e T_e             zchim(3)    chi_e
c       zflxm(4)   n_Z                 zchim(4)    D_q
c       zflxm(5)   n_Z T_Z             zchim(5)    chi_q
c
c..local variables added 4-jan-95
c
      real bt, bt1,  zeni, k1, k2
      real zkpsh, zkpsz, zalp, zalf, zanorm, zbetae, zrav
c
!| 
!| Definitions of some of the internal variables:
!| 
!| \renewcommand{\arraystretch}{1.4}
!| \begin{center}
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| \\
!| {\tt zgth} & $ - \frac{R}{T_H} \Partial{T_H}{r} $ \\
!| {\tt zgte} & $ - \frac{R}{T_e} \Partial{T_e}{r} $ \\
!| {\tt zgtz} & $ - \frac{R}{T_Z} \Partial{T_Z}{r} $ \\
!| {\tt zgnh} & $ - \frac{R}{n_H} \Partial{n_H}{r} $ \\
!| {\tt zgnz} & $ - \frac{R}{n_Z} \Partial{n_Z}{r} $ \\
!| {\tt zgph} & $ - \frac{R}{n_H T_H} \Partial{n_H T_H}{r} $ \\
!| {\tt zgpe} & $ - \frac{R}{n_e T_e} \Partial{n_e T_e}{r} $ \\
!| \end{tabular}  \end{center}
!| 
c
      save idim, zepsmach, zepsqrt
     &  , zone, ztwo, zthree, zfour, zfive, zhalf, zquarter, ztvr, zftr
     &  , ztwohlf, inital
c
c..initialize variables
c
      lprint = lprintin
c
      if ( nout .lt. 1  .or.  nout .gt. 99 ) nout = 6
c
      ieq = max ( 2, neq )
c
      vef = cetain(15) * vef
      imatrx = min ( ieq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )
c
      if ( inital ) then
c
        idim = idp
c
        zone   = 1.0
        ztwo   = 2.0
        zthree = 3.0
        zfour  = 4.0
        zfive  = 5.0
c
        zhalf    = zone  / ztwo
        ztwohlf  = ztwo  + zhalf
        zquarter = zone  / zfour
        ztvr     = ztwo  / zthree
        zftr     = zfive / zthree
c
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
c..Print Header
c
        write (nout,*)
        write (nout,*)
     &   'Weiland-Nordman eigenvalue equations, subroutine etaw17flux'
        write (nout,*)
c
          write (nout,*)
     &        ' Eigenvalues computed using ACM/TOMS routine 535'
c
        if ( ieq .eq. 11) then
          write (nout,*)
          write (nout,*)
     &     'Eleven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Vp_H, Av, K, and Vp_Z'
        elseif ( ieq .eq. 10) then
          write (nout,*)
          write (nout,*)
     &     'Ten eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Vp, Av, and K'
        elseif ( ieq .eq. 9) then
          write (nout,*)
          write (nout,*)
     &     'Nine eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z,
     &     F, Av, K, and Vp in strong ballooning limit'
        elseif ( ieq .eq. 8) then
          write (nout,*)
          write (nout,*)
     &     'Eight eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, and Vp'
        elseif ( ieq .eq. 7) then
          write (nout,*)
          write (nout,*)
     &     'Seven eqns e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z, and F'
        elseif ( ieq .eq. 6 ) then
          write (nout,*)
          write (nout,*)
     &     'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
               elseif ( ieq .eq. 5 ) then
          write (nout,*)
          write (nout,*)
     &     ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
        elseif ( ieq .eq. 4 ) then
          write (nout,*)
          write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
        elseif ( ieq .eq. 2 ) then
           write (nout,*)
           write (nout,*) ' Two eqns for e phi/T_e and T_H'
        else
          write (nout,*)
          write (nout,*) ieq,' = ieq in sbrtn etaw17flux'
          write(nout,*)'the value of ieq is wrong in sbrtn etaw17flux'
          nerr = 2
          return
        endif
c
        inital = .false.
c
      endif
c
c..end of initialization
c
      nerr = 0
c
c..print header
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine etaw17flux'
        write (nout,*) '(all frequencies normalized by omega_{De})'
        write (nout,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
        write (nout,108) 'gnein', gnein
        write (nout,108) 'gnhin', gnhin
        write (nout,108) 'gnzin', gnzin
        write (nout,108) 'gtein', gtein
        write (nout,108) 'gthin', gthin
        write (nout,108) 'gtzin', gtzin
        write (nout,108) 'tauhin', tauhin
        write (nout,108) 'tauzin', tauzin
        write (nout,108) 'ftrapein', ftrapein
        write (nout,108) 'fnzin', fnzin
        write (nout,108) 'czin', czin
        write (nout,108) 'azin', azin
        write (nout,108) 'fnsin', fnsin
        write (nout,108) 'betaein', betaein
        write (nout,108) 'betahin', betahin
        write (nout,108) 'betazin', betazin
        write (nout,108) 'vefin', vef
        write (nout,108) 'qin', q
        write (nout,108) 'shearin', shear
        write (nout,108) 'kappain', kappa
        write (nout,108) 'ekyrhoin', ekyrhoin
        write (nout,108) 'ekparlin', ekparlin
        write (nout,108) 'wexb', wexb
 108    format (1x,a8,' = ',1pe14.6,',')
c
      endif
c
c...copy to local variables
c
      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
c
c..check validity of input data
c
      if ( neq .lt. 2 ) then
         write(nout,*) ' neq .lt. 2 in sbrtn etaw17flux'
         nerr = 2
         return
      else if ( ndim .gt. idim ) then
         write(nout,*) ' ndim .gt. idim in sbrtn etaw17flux'
         nerr = 2
         return
      endif
c
c..initialize arrays
c
      do j1=1,ndim
        omega(j1)   = 0.0
        gamma(j1)   = 0.0
        fluxout(j1) = 0.0
        chieff(j1)  = 0.0
        perform(j1) = 0.0
      enddo
c
      do j1=1,idp
        zomega(j1)  = 0.0
        zgamma(j1)  = 0.0
        zalfr (j1)  = 0.0
        zalfi (j1)  = 0.0
        zbeta (j1)  = 0.0
        zchim(j1)   = 0.0
        zflxm(j1)   = 0.0
        do j2=1,idp
          zevec(j1,j2)  = ( 0.0, 0.0 )
          zamr (j1,j2)  = 0.0
          zami (j1,j2)  = 0.0
          zbmr (j1,j2)  = 0.0
          zbmi (j1,j2)  = 0.0
          zvr  (j1,j2)  = 0.0
          zvi  (j1,j2)  = 0.0
        enddo
      enddo
c
c..set up initial gradients
c
      zimp   = czin
      zmass  = azin
      zfnz   = fnzin * zimp
      zfns   = fnsin
c
c..Save zgns = - R ( d n_s / d r ) / n_e
c    where n_s is the fast ion density
c
      zgns = zgne - zgnh * ( zone - zfnz - zfns ) - zgnz * zfnz
c
c      if ( zfns .lt. zepsqrt .and. abs(zgns) .gt. zepsqrt) then
c        write(nout,*) ' Gradients do not fulfill neutrality condition'
c         write(nout,*) ' zgns .ne. zgne-zgnh*(1-zfnz-zfns)-zgnz*zfnz'
c         nerr = 2
c         return
c      endif
c
c
c..compute the rest of the dimensionless variables needed
c
c     ztauz = T_Z / ( Z T_e ),     zfnz  = Z n_Z / n_e
c     zimp  = Z,                   zflh  = k_y^2 \rho_{sH}^2
c     zmass = m_Z / m_H,           zflz  = k_y^2 \rho_{sZ}^2
c

      ztauh  = tauhin
      zeni = zhalf * zgne
      zbetae = max ( cetain(20) * betaein, zepsqrt )
      zimp   = max ( czin, zone )
      zmass  = max ( azin, zone )
      zfnz   = fnzin * zimp
      zft    = ftrapein
      zflh   = ekyrhoin**2
c
      if ( neq .gt. 4 ) then
        ztauz  = tauzin / czin
        zflz   = zmass * zflh / zimp**2
      else
        ztauz  = zone
        zflz   = 0.0
      endif
      zetae = zgte / zgne
      zetah = zgth / zgne
c
c..diagnostic output
c
      if ( lprint .gt. 2 ) then
        write (nout,*)
        write (nout,*) '--------------------------------------'
        write (nout,*)
        write (nout,*)
        write (nout,*) zgnh,' = zgnh'
        write (nout,*) zgne,' = zgne'
        write (nout,*) zgnz,' = zgnz'
        write (nout,*) zgns,' = zgns'
        write (nout,*) zgth,' = zgth'
        write (nout,*) zgte,' = zgte'
        write (nout,*) zgtz,' = zgtz'
        write (nout,*) zetae,' = zetae'
        write (nout,*) zetah,' = zetah'
        write (nout,*)
c
        write (nout,*) ztauh,' = ztauh'
        write (nout,*) ztauz,' = ztauz'
        write (nout,*)
        write (nout,*) zft,' = zft'
        write (nout,*) zimp,' = zimp'
        write (nout,*) zmass,' = zmass'
        write (nout,*) fnzin,' = fnz'
        write (nout,*) zfnz,' = zfnz'
        write (nout,*) zfns,' = zfns'
        write (nout,*) zflh,' = zflh'
        write (nout,*) zflz,' = zflz'
        write (nout,*)
        write (nout,*) zepsqrt,' = zepsqrt'
        write (nout,*) zepsmach,' = zepsmach'
        write (nout,109) 'gnsin', zgns
 109  format (1x,a8,' = ',1pe14.6,', computed from quasi-neutrality')

      endif

      if ( ieq .gt. 4 ) then
c
c..constants from Nilsson and Weiland, NF 34 (1994) 803  section 2
c  and from Weiland and Hirose NF 32 (1992) 151
c
c  ** Note:  zanorm = 1./kpc is used to normalize A_\parallel and K
c
        zanorm = cetain(11)
c
        bt  = 1.5
        bt1 = bt - ztwohlf
c
c  Expressions from Weiland and Hirose, NF 32 (1992) 151  Eq. (6)
c  Note that k1 rather than k1**2 appears under sqrt in zalp.
c  This is the result of a numerical fit rather than analytic solution.
c  Also, only free electrons are used from zbetae in zalf
c  Hence, zbetae -> zbetae * (zone-zft).
c
        if ( zgnh + zgth .lt. zepsqrt ) then
c
          zalp  = 0.0
          zalf  = 0.0
          zkpsh = 0.0
          zkpsz = 0.0
          zrav  = 1.0
c
        else
          k1 = zquarter * q * q * zflh *
     &         sqrt( zhalf * (zgnh + zgth) * ztauh / (zone-zft) )
          k2 = q * q * zflh * zflh * zhalf * (zgnh + zgth) *
     &         ztauh / ( zone-zft)

          zalp = zhalf * ( k1 + sqrt( k1 + shear * shear * k2) )
          zalf = zalp / (ztwo * zflh * q * q * zbetae * (zone-zft))

          zkpsh = cetain(10) * zhalf * sqrt( zalp / zflh ) / q

          zrav = zone
     &      + cetain(25) * ( zquarter * ( ztwo * shear - zone
     &      + kappa * kappa * ( shear - zone ) * ( shear - zone ) )
     &        / zalp )

          if ( zmass .gt. zepsqrt ) then
            zkpsz = zkpsh / sqrt ( zmass )
          else
            zkpsz = zkpsh
          endif
c
        endif
c
        if (lprint .gt. 5 ) then
          write (nout,*)
          write (nout,*) k1,'  = k1'
          write (nout,*) k2,'  = k2'
          write (nout,*) zalp,' = zalp'
          write (nout,*) zalf,' = zalf'
          write (nout,*) zkpsh,' = zkpsh'
          write (nout,*) zkpsz,' = zkpsz'
          write (nout,*) kappa,' = kappa'
          write (nout,*) zrav,' = rav'
          write (nout,*)
          write (nout,*) bt1,' = bt1'
          write (nout,*) zanorm,' = zanorm'
        endif
c
      endif
c
c
      if (ieq .eq. 11) then
c
c%%%%%%%%%%
c..Eleven equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel hydrogenic and impurity ion motion,
c  electromagnetic (finite beta) effects
c
c  equations for
c  e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vph, Av, K, Vpz
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &  'Eleven eqns for '
     &  ,'e phi / T_e, T_H, n_H, T_et, n_Z, T_Z, F, Vph, Av, K, Vpz'
      endif

c
c  hydrogenic density
c

      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
      zamr(5,11) = zkpsz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  hydrogenic parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
      zamr(8,9) = - zkpsh * zhalf * ( zgth + zgnh ) * ztauh
c
      zbmr(8,8) = zone
      zbmr(8,9) = zkpsh
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = (zhalf * (zgne + zgte) - zalf * zflh * zrav)*zanorm
      zamr(9,10) = - zeni * zanorm
      zamr(9,11) = zkpsz * zfnz / ( zone - zft )
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av, K = omega * Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c  impurity parallel ion motion
c
      zamr(11,1) = zimp * zkpsz
      zamr(11,5) = zimp * zkpsz * ztauz
      zamr(11,6) = zimp * zkpsz * ztauz
      zamr(11,9) = - zimp * zkpsz * zhalf * ( zgtz + zgnz ) * ztauz
c
      zbmr(11,9) = zimp * zkpsz
      zbmr(11,11) = zone
c
c
      elseif (ieq .eq. 10) then
c
c%%%%%%%%%%
c..Ten equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel ion motion, and electromagnetic effects
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Ten eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
c  hydrogenic density
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
      zamr(8,9) = - zkpsh * zhalf * ( zgth + zgnh ) * ztauh
c
      zbmr(8,8) = zone
      zbmr(8,9) = zkpsh
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = (zhalf * (zgne + zgte) - zalf * zflh * zrav)*zanorm
      zamr(9,10) = - zeni * zanorm
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av, K = omega * Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c
      elseif (ieq .eq. 9) then
c
c%%%%%%%%%
c..Nine equations with impurities, trapped electrons, FLR,
c  collisional effects, and electromagnetic effects
c  and parallel ion motion in strong ballooning limit
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout, *)
     &'Nine eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
      H = cetain(12) * zhalf * abs( shear ) / max ( q, zepsqrt )
c
c  ion continuity
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zami(1,1) = - H
      zamr(1,2) = - ztauh
      zami(1,2) = - ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = - ztauh*H
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  total electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,8) = - ( zone - zft ) * zeni * zanorm
      zami(3,8) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,9) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,9) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,9) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - ztwohlf * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(8,1) = zeni
      zamr(8,8) = ( zhalf*(zgne+zgte)-zalf*zflh*zrav ) * zanorm
      zamr(8,9) = - zeni * zanorm
c
      zbmr(8,1) = zone
      zbmr(8,8) = zanorm
      zbmr(8,9) = - zanorm
c
c   time derivative of Av
c
      zamr(9,9) = zone
c
      zbmr(9,8) = zone
c
c
      elseif (ieq .eq. 8) then
c
c%%%%%%%%%
c..Eight equations with impurities, trapped electrons, FLR,
c  collisional effects, and parallel ion motion
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Eight eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Vp'
      endif
c
c  ion continuity
c
      zamr(1,1) = -zone + zhalf*(zgnh-zflh*ztauh*(zgnh+zgth))
      zamr(1,2) = -ztauh
      zamr(1,3) = -ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
      zamr(2,2) = -ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = -ztvr
c
c  total electron density expressed through quasineutrality
c
      zamr(3,1) = zft*zeni - zone
      zami(3,1) = vef*(zone-zft)
      zamr(3,3) = zone - zfnz- zfns
      zami(3,3) = -vef*(zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
c
c  trapped electron energy
c
      zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
      zami(4,1) = vef*ztvr*(bt-ztwohlf*(zone-zft))
      zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
      zamr(4,4) = zft*zftr
      zami(4,5) = -vef*ztvr*bt1*zfnz
      zami(4,7) = -zftr*vef*zft
c
      zbmr(4,1) = (zone-zft)*ztvr
      zbmr(4,3) = - (zone - zfnz - zfns)*ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*ztvr
c
c  impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c  impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  parallel ion motion Vp = Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh*ztauh
      zamr(8,3) = zkpsh*ztauh
c
      zbmr(8,8) = zone
c
c
      elseif ( ieq .eq. 7 ) then
c
c%%%%%%%%%
c..Seven equations with impurities, trapped electrons, FLR,
c  and collisions
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here, F is defined as F = GM*e phi/T_e
c  where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (nout,*)
      write (nout,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F'
      endif
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = - zone + zft * zeni
        zami(3,1) = vef*(zone-zft)
        zamr(3,3) = zone - zfnz - zfns
        zami(3,3) = -vef*(zone - zfnz - zfns)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
        zami(4,1) = vef*ztvr*(bt-ztwohlf*(zone-zft))
        zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
        zamr(4,4) = zft * zftr
        zami(4,5) = -vef*ztvr*bt1*zfnz
        zami(4,7) = -zftr*vef*zft
c
        zbmr(4,1) = ( zone - zft ) *ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) *ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = - zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c  variable F
c
        zamr(7,1) = zhalf*zgte - zone
        zami(7,1) = vef
        zamr(7,7) = zone
        zami(7,7) = -vef
c
        zbmr(7,1) = -zone
        zbmr(7,7) = zone
c
c
      elseif ( ieq .eq. 6 ) then
c
c%%%%%%%%%
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
c
      if ( lprint .gt. 5 ) then
        write (nout,*)
        write (nout,*)
     &   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = zft*zeni - zone
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * zhalf * (zgte - ztvr*zgne)
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) * ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = -zone + zhalf * (zgnz - zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf * (zgtz - ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c
      elseif ( ieq .eq. 5 ) then
c
c%%%%%%%%%%
c..5 equations with trapped electrons, FLR effects, and parallel ion motion
c
c  equations for e phi/T_e, T_H, n_i, T_e, and Vp
c
       if ( lprint .gt. 5 ) then
         write (nout,*)
         write (nout,*)
     &    ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
       endif
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,5) = zkpsh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation.
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
c
c  Parallel ion motion Vpi/Cs
c
        zamr(5,1) = zkpsh
        zamr(5,2) = zkpsh*ztauh
        zamr(5,3) = zkpsh*ztauh
c
        zbmr(5,5) = zone
c
      elseif ( ieq .eq. 4 ) then
c
c%%%%%%%
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
       if ( lprint .gt. 5 ) then
         write (nout,*)
         write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
       endif
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
      elseif ( ieq .eq. 2 ) then
c
c%%%%%%%%
c..two equations when trapped particles and FLR effects omitted
c
c  equations for e phi/T_e and T_H
c
       if ( lprint .gt. 5 ) then
         write (nout,*)
         write (nout,*) ' Two eqns for e phi/T_e and T_H'
       endif
c
        zamr(1,1) = zhalf * zgnh - ztauh - zone
        zamr(1,2) = - ztauh
        zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(1,1) = zone
        zbmr(1,2) = 0.0
        zbmr(2,1) = - ztvr
        zbmr(2,2) = zone
c
c
      else

        write (nout,*)
        write (nout,*) ieq,' = ieq in sbrtn etaw17flux'
        write (nout,*) 'the value of ieq is wrong in sbrtn etaw17flux'
        nerr = 2
        return

      endif
c
c..find the eigenvalues and eigenvectors using ACM/TOMS routine 535
c
      ifail = -1
c
c..Update A matrix w.r.t wexb and save local copies of the matrices
c
      do j1=1,ieq
       do j2=1,ieq
        zamr(j1,j2)  = zamr(j1,j2) + abs(wexb)*zbmi(j1,j2)
        zami(j1,j2)  = zami(j1,j2) - abs(wexb)*zbmr(j1,j2)
        zamrt(j1,j2) = zamr(j1,j2)
        zamit(j1,j2) = zami(j1,j2)
        zbmrt(j1,j2) = zbmr(j1,j2)
        zbmit(j1,j2) = zbmi(j1,j2)
       enddo
      enddo

c
c..diagnostic output for defining matrices
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ieq
c
        write (nout,*)
        write (nout,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmi(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmi(j1,j2),j2=1,ieq)
        enddo
 192  format (1p10e12.4)
      endif
c
c... Solution of the generalized eigenvalue problem
c
      call tomsqz(idim,ieq,zamr,zami,zbmr,zbmi,ZALFR,ZALFI,ZBETA,
     &            ZVR,ZVI,IFAIL)
c
c... Storing the results -  eigenvalues and eigenvectors
c
      zgamax = 0.0
      do j=1,ieq
        ztemp1 = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
        zomega(j) = zalfr(j)  / ztemp1
        zgamma(j) = zalfi(j)  / ztemp1
        omega(j) = zomega(j)
        gamma(j) = zgamma(j)
        zgamax = max ( zgamax, zgamma(j) )
        do j1=1,ieq
          zevec(j1,j) = cmplx ( zvr(j1,j), zvi(j1,j))
        end do
      enddo

      nerr = ifail

c
c... If no unstable roots are found set fluxes to zero and leave
c

      if ( zgamax .lt. zepsqrt ) go to 90
c
c..check the eigenfunctions
c
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
      endif
c
c  Real and imaginary parts
c
      zerrmax = 0.0
      nmodes = 0
c
c Loop over number of possible modes
c
      do j=1,ieq
c
c Start an error check on the eigensolution comparin RHS w LHS
c
        do j1=1,ieq
          ztempa(j1) = 0.0
          ztempb(j1) = 0.0
          do j2=1,ieq
            zerreal =
     &            zamrt(j1,j2) * real (zevec(j2,j))
     &          - zamit(j1,j2) * aimag(zevec(j2,j))
     &          - zbmrt(j1,j2) * real (zevec(j2,j)) * zomega(j)
     &          + zbmrt(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * real (zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zomega(j)

            zerimag =
     &            zamrt(j1,j2) * aimag(zevec(j2,j))
     &          + zamit(j1,j2) * real (zevec(j2,j))
     &          - zbmrt(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
     &          - zbmrt(j1,j2) * real (zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          - zbmit(j1,j2) * real (zevec(j2,j)) * zomega(j)

            ztempa(j1) = ztempa(j1) + zerreal
            ztempb(j1) = ztempb(j1) + zerimag

          enddo
          zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
        enddo

        if ( lprint .gt. 12 ) then
          write (nout,*)
          write (nout,*) ' LHS - RHS for j =  ',j
          do j1=1,ieq
            write (nout,142) ztempa(j1), ztempb(j1)
          enddo
        endif
c
!| 
!| The effective diffusivities can be computed directly from the
!| eigenvectors in the following way:
!| Consider, for example, the flux of hydrogen particles produced by
!| the perturbed $ E \times B $ motion of the plasma
!| \[ \Gamma_H = \tilde{n}_H \tilde{v}_E^* + c.c.
!|  = 2 ( {\rm Re} \tilde{n}_H {\rm Im} \tilde{\phi}
!|      - {\rm Im} \tilde{n}_H {\rm Re} \tilde{\phi} ) k_y / B. \]
!| Using the quasi-linear assumption, the saturation level is
!| approximated by
!| \[ \frac{e \tilde{\phi}}{T_e}
!|  \approx \frac{1}{k_x \rho_{sH}} \frac{ \gamma }{ k_y c_{sH} }
!|  = \frac{2}{R k_x} \frac{\gamma}{\omega_{De}}. \]
!| Hence, the hydrogen particle flux is given by
!| \[ \frac{ F_H }{ n_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} \hat{n}_H {\rm Im} \hat{\phi}
!|   - {\rm Im} \hat{n}_H {\rm Re} \hat{\phi} ) / | \hat{\phi} |^2.
!|    \]
!| Correspondingly, the hydrogen heat (here, $ n_H T_H $) flux is given by
!| \[ \frac{ F_{p_H} }{ n_H T_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} ( \hat{n}_H + \hat{T}_H ) {\rm Im} \hat{\phi}
!|   - {\rm Im} ( \hat{n}_H + \hat{T}_H ) {\rm Re} \hat{\phi} )
!|     / | \hat{\phi} |^2.
!|    \]
!| 
!| The case of the electron heat flux (here, the flux of $ n_e T_e $)
!| is a little more complicated.  Note that
!| 
c
c..compute effective diffusivities directly from eigenvalues
c  assume eigenvectors are arranged in the order of
c  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
c
c... Calculate norm of solution in terms of (e\phi/T_e) **2
c
        ztemp1 =  real(zevec(1,j)) *  real(zevec(1,j))
     &         + aimag(zevec(1,j)) * aimag(zevec(1,j))
c
c... Check if current mode j is unstable
c
        if ( zgamma(j).gt.zepsqrt .and. ztemp1.gt.zepsqrt ) then

           nmodes = nmodes + 1
c
c Define fluxes : Thermal hydrogen flux
c
          zreal =  real(zevec(2,j)) +  real(zevec(3,j))
          zimag = aimag(zevec(2,j)) + aimag(zevec(3,j))
c
c...phase difference
c
          zphsph = - ( zimag * real(zevec(1,j))
     &             -   zreal * aimag(zevec(1,j)) ) / ztemp1
c
c...flux from j:th mode
c
          zflxph = 2.0 * zphsph * zgamma(j) * zgamma(j)

c
c...flux summed over all unstable modes
c
          zflxm(1) = zflxm(1) + zflxph
c
c Define hydrogen density flux - phase shift
c
          zphsnh = - ( aimag(zevec(3,j)) * real(zevec(1,j))
     &             - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c...Flux from j:th mode
c
          zflxnh = 2.0 * zphsnh * zgamma(j) * zgamma(j)
c
c...Flux summed over all unstable modes
c
          zflxm(2) = zflxm(2) + zflxnh
c
c Define thermal electron flux
c
          if ( ieq .gt. 3 ) then
            zreal =  real(zevec(4,j))
     &            + ( zone - zfnz - zfns ) * real(zevec(2,j))
     &            + zfnz * real(zevec(5,j))
            zimag = aimag(zevec(4,j))
     &            + ( zone - zfnz - zfns ) * aimag(zevec(2,j))
     &            + zfnz * aimag(zevec(5,j))
c
c  Note, the electron heat flux is reduced by the fraction of
c  trapped electrons - phase shift
c
            zphspe = - zft * ( zimag * real(zevec(1,j))
     &                       - zreal * aimag(zevec(1,j)) ) / ztemp1

c
c... Flux from the j:th mode
c
            zflxpe = 2.0 * zphspe * zgamma(j)* zgamma(j)
c
c... Flux summed over all unstable modes
c
            zflxm(3) = zflxm(3) + zflxpe
c
          endif
c
c... Impurity density flux
c
          if ( ieq .gt. 4 ) then
c
c... phase difference
c
            zphsnz = - (aimag(zevec(5,j))*real (zevec(1,j))
     &               -  real (zevec(5,j))*aimag(zevec(1,j)) ) / ztemp1

c
c... Flux from the j:th mode
c
            zflxnz =  2.0 * zphsnz * zgamma(j) * zgamma(j)
c
c... Flux summed over all unstable modes
c
            zflxm(4) = zflxm(4) + zflxnz
c
          endif
c
c Define the impurity heat flux
c
          if ( ieq .gt. 5 ) then
c
            zreal =  real(zevec(5,j)) +  real(zevec(6,j))
            zimag = aimag(zevec(5,j)) + aimag(zevec(6,j))
c
c... phase difference
c
            zphspz = - ( zimag * real (zevec(1,j))
     &               -   zreal * aimag(zevec(1,j)) ) / ztemp1
c
c... Flux from the j:th mode
c
            zflxpz = 2.0 * zphspz * zgamma(j) * zgamma(j)
c
c... Flux summed over all unstable modes
c
            zflxm(5) = zflxm(5) + zflxpz
c
          endif
c
c..header for diagnostic printout of frequencies, fluxes, and phases
c
        if ( letain(29) .gt. 0 .and. lprint .gt. 0) then
          write (nout,134)
c
 134  format (//' Diagnostic printout of frequencies, fluxes and phases'
     &  /' Note: fluxph = flux of hydrogenic thermal energy'
     &  ,' = 2.0 * gamma^2 * phaseph'
     &  /9x,' (phaseph is related to the phases of the perturbations)'
     &  /7x,'fluxnh = flux of hydrogenic ions = 2.0 * gamma^2 * phasenh'
     &  /7x,'fluxpe = flux of electron thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepe'
     &  /7x,'fluxnz = flux of impurity ions = 2.0 * gamma^2 * phasenz'
     &  /7x,'fluxpz = flux of impurity thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepz'
     &  //1x,'radius',t10,'omega',t20,'gamma'
     &  ,t30,'fluxph',t40,'phaseph',t50,'fluxnh',t60,'phasenh'
     &  ,t70,'fluxpe',t80,'phasepe',t90,'fluxnz',t100,'phasenz'
     &  ,t110,'fluxpz',t120,'phasepz  #m')
c
c..diagnostic printout of frequencies and fluxes mode by mode
c
        write (nout,135) cetain(29), zomega(j), zgamma(j)
     &    , zflxph, zphsph, zflxnh, zphsnh, zflxpe, zphspe
     &    , zflxnz, zphsnz, zflxpz, zphspz
        endif
c
 135  format (0pf7.3,1p12e10.2,' #m')
c
        endif
c
c...End of flux and diffusivity definitions loop
c
      enddo

c
c... Error diagnostics
c
      if ( lprint .gt. 0 ) then
c        if (abs(zerrmax) .gt. 100*zepsqrt) nerr = max (nerr , 1)
        if (abs(zerrmax) .gt. zepsqrt) nerr = max (nerr , 1)
        write (nout,*)
        write (nout,*) zerrmax,' = zerrmax'
        write (nout,*) nerr,
     &       ' = nerr, error in eigenvalue in etaw17flux'
      endif

 142  format (1p10e12.4)
c
c..compute effective total diffusivities
c
      zchim(1) = zflxm(1) / sign( max( abs( zgth ), zepsqrt ),  zgth )
      zchim(2) = zflxm(2) / sign( max( abs( zgnh ), zepsqrt ),  zgnh )
      zchim(3) = zflxm(3) / sign( max( abs( zgte ), zepsqrt ),  zgte )
      zchim(4) = zflxm(4) / sign( max( abs( zgnz ), zepsqrt ),  zgnz )
      zchim(5) = zflxm(5) / sign( max( abs( zgtz ), zepsqrt ),  zgtz )
c
c..save effective diffusivities and fluxes
c
      do j1=1,min(ieq-1,5)
        chieff(j1)  = zchim(j1)
        fluxout(j1) = zflxm(j1)
      enddo

c
      if ( lprint .gt. 2 ) then
c
c..print eigenvalues and eigenfunctions
c
        write (nout,121)
        do j=1,ieq
          write (nout,122) zomega(j), zgamma(j)
        enddo

 121    format (/' Solution of the eigenvalue equations'
     &   /t4,'zomega',t18,'zgamma')
 122    format (1p2e14.5,i5)
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (zchim(j1),j1=1,ieq-1)
c
      endif
c
      if ( lprint .gt. 99 ) then
c
        write (nout,*)
        write (nout,*) ' Eigenvectors zevec(j1,j2) j2->'
        do j1=1,ieq
          write (nout,124) (zevec(j1,j2),j2=1,ieq)
 124      format (2(1p12e11.3))
        enddo
c
        write (nout,*)
      endif
c
      if ( ifail .gt. 0 ) then
        write(nout,*)
     &   'ifail .gt. 0 after call tomsqz in sbrtn etaw17flux'
        nerr = 2
        return
      endif
c
c
c
c
      if ( imatrx .lt. 1 ) go to 90
c
c..diagnostic printout
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ' vector zflxm(j1)'
        do j1=1,ieq-1
          write (nout,132) zflxm(j1)
        enddo
c
      endif
c
  90  continue
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn etaw17flux'
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,imatrx)
c
      endif
c
c
      return
c
 130    format (/t3,'chi_i',t16,'dif_h',t29,'chi_e'
     &    ,t42,'dif_Z',t55,'chi_Z')
 132    format (1p11e12.4)
c
      end
!| 
!| Note that
!| $$ \omega_{De} = 2 k_\perp T_e / e B R = 2 k_\perp \rho_s c_s / R $$
!| $$ D \propto \gamma / k^2
!|   \propto \rho_s^2 \omega_{De}
!|     \frac{ \gamma / \omega_{De} }{ k^2 \rho_s^2 }. $$
!| 
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{weil92a} J. Weiland,
!| ``Low Frequency Modes Associated with Drift Motions in
!| Inhomogeneous Plasmas,''
!| report CTH--IEFT/PP-1992-17,
!| Institute for Electromagnetic Field Theory and Plasma Physics,
!| Chalmers University of Technology,
!| G\"{o}teborg, Sweden.
!| 
!| \bibitem{froj92a} M. Fr\"{o}jdh, M. Liljestr\"{o}m, H. Nordman,
!| ``Impurity effects on $\eta_i$ mode stability and transport,''
!| Nuclear Fusion {\bf 32} (1992) 419--428.
!| \bibitem{nord92a} H. Nordman and J. Weiland, ``Comments on
!| `Ion-temperature-gradient-driven
!| transport in a density modification experiment on the tokamak fusion test
!| reactor [Phys. Fluids {\bf B4} (1992) 953]' ''.
!| \bibitem{weil92b} J. Weiland,
!| ``Nonlinear effects in velocity space and drift wave
!| transport in tokamaks,'' Phys. Fluids {\bf B4} (1992) 1388--1390.
!| \bibitem{weil92c} J. Weiland and H. Nordman, ``Drift wave model for inward
!| energy transport in tokamak plasmas,'' Institute for Electromagnetic Field
!| Theory and Plasma Physics, Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.
!| \bibitem{weil92d} J.Weiland and A. Hirose, ``Electromagnetic and kinetic
!| effects on the ion temperature gradient mode,'' Nucl. Fusion {\bf 32} (1992)
!| 151--155.
!| \bibitem{nord91a} H. Nordman and J. Weiland, ``The concept of marginal
!| stability and recent experimental results from the TFTR tokamak,'' Institute
!| for Electromagnetic Field Theory and Plasma Physics, Gothenburg, Sweden, (1991)
!| CTH-IEFT/PP-1991-26 ISSN.
!| \bibitem{weil91a} J. Weiland and H. Nordman, ``Enhanced confinement regimes in
!| transport code simulations of toroidal drift wave transport,'' Nucl. Fusion
!| {\bf 31} (1991) 390--394.
!| \bibitem{nord90a} H. Nordman, J. Weiland, and A. Jarmen, ``Simulation of
!| toroidal drift mode turbulence driven by temperature gradients and electron
!| trapping,'' Nucl. Fusion {\bf 30} (1990) 983--996.
!| \bibitem{weil89a} J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, ``Diffusive
!| particle and heat pinch effects in toroidal plasmas,'' Nucl. Fusion {\bf 29}
!| (1989) 1810--1814.
!| \bibitem{nord89a} H. Nordman and J. Weiland, ``Transport due to toroidal
!| $\eta_i$ mode turbulence in tokamaks,'' Nucl. Fusion {\bf 29} (1989) 251--263.
!| \bibitem{ande88a} P. Andersson and J. Weiland, ``A fully toroidal fluid analysis
!| of the magnetohydrodynamic ballooning mode branch in tokamaks,'' Phys. Fluids
!| {\bf 31} (1988) 359--365.
!| \bibitem{jarm87a} A. Jarm\'{e}n, P. Andersson, and J. Weiland, ``Fully toroidal
!| ion temperature gradient driven drift modes,'' Nucl. Fusion {\bf 27} (1987)
!| 941--949.
!| \end{thebibliography}
!| 
!| \end{document}

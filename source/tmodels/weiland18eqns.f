!|  %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-
!| 
!| %  \documentstyle[12pt]{article}
!| \documentclass[12pt]{article}
!| \usepackage{longtable}
!| 
!|  \headheight 0pt \headsep 0pt
!|  \topmargin 0pt  \textheight 9.0in
!|  \oddsidemargin 0pt \textwidth 6.5in
!| 
!|  \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!|  \newcommand{\jacobian}{{\cal J}}
!| 
!|  \begin{document}
!| 
!|  \begin{center}
!|  {\bf {\tt weiland18eqns.f} \\
!|  Toroidal Ion Temperature Gradient Mode \\
!|  \vspace{1pc}
!|  Glenn Bateman \\
!|  Lehigh University, Physics Department \\
!|  16 Memorial Drive East, Bethlehem, PA 18015 \\
!|  \vspace{1pc}
!|  Jan Weiland and Hans Nordman \\
!|  Chalmers University of Technology \\
!|  G\"{o}teborg, Sweden \\
!|  \vspace{1pc}
!|  Jon Kinsey \\
!|  General Atomics \\
!|  P.O. Box 85608, San Diego, CA 92186} \\
!|  \vspace{1pc}
!|  \today
!|  \end{center}
!| 
!| This subroutine evaluates the eigenvalue equations for Ion Temperature
!| Gradient (ITG) and trapped electron modes (TEM) derived by Jan Weiland, 
!| H.~Nordman and their group in G\"{o}teborg Sweden. The equations in this
!| routine include fast hydrogenic ions, impurities, trapped electron and
!| finite Larmor radius effects. New options include parallel ion motion,
!| finite beta and collisional effects together with an approximate 
!| treatment of $E\times B$ flow shear reduction of transport.
!| 
!| The essential idea is to linearize the fluid equations with 
!| magnetic drifts for a given Fourier harmonic
!| \[ u = u_0(x) + \tilde{u} \exp\left[ -i (\omega t
!|    + \vec{k} \cdot \vec{x} ) \right] \]
!| Compute the eigenvalues and eigenvectors from these equations and then
!| compute the quasi-linear heat and particle fluxes.
!| Many of the steps in this derivation can be found in reference
!| \cite{weil00a}.
!| 
!| The non-linear fluid equations with their associated drift velocities
!| are given in the first section in this file.
!| The normalized and linearized are derived in the second section.
!| The matrix form of the equations are written out in the section were
!| they are implemented in the computer code below.
!| Finally, the quasi-linear form of the heat and particle fluxes are
!| given where they are computed.
!| 
!| \section{Fluid Equations}
!| 
!| The fundamental equations used in the current version of the Weiland
!| model are the fluid equations for each plasma species.
!| For hydrogenic and impurity ion densities $n_i$,
!| the equation of continuity is
!| \begin{equation}
!|   \Partial{n_i}{t} + {\mathbf \nabla} \cdot ( n_i {\mathbf v}_i ) = 0
!| \end{equation}
!| \[ {\mathbf v}_i = {\mathbf v}_E + {\mathbf v}_{*i} + {\mathbf v}_{Pi}
!|    + {\mathbf v}_{\pi i} + \hat{{\mathbf B}} v_{\parallel i} \]
!| where the fluid flows include the E cross B drift,
!| \[ {\mathbf v}_E = {\mathbf E} \times \hat{{\mathbf B}} / B
!|  = \hat{{\mathbf B}} \times {\mathbf \nabla} \phi / B; \]
!| the diamagnetic drift ($ Z_i = 1 $ for hydrogen isotopes and
!| $ Z_i = -1 $ for electrons) as well as
!| the drift due to the off-diagonal elements of the stress tensor,
!| \[ {\mathbf v}_{*i}
!|  = \frac{ \hat{{\mathbf B}} \times {\mathbf \nabla}
!|  ( n_i T_i )}{Z_i e n_i B} ;
!| \hspace{2em}
!| {\mathbf v}_{\pi i}
!|  = \frac{ \hat{{\mathbf B}} \times {\mathbf \nabla} \cdot {\mathbf \pi}_i}{
!|   Z_i e n_i B } \]
!| the polarization drift,
!| \[ v_{Pi} = \frac{ d E }{ dt } / ( B \Omega_i ); \hspace{2em}
!| \Omega_i \equiv \frac{ Z_i e B }{ m_i}; \]
!| and the ion flow along field lines
!| $ v_{\parallel i} $.
!| 
!| The parallel ion motion $ v_{\parallel i} $ is determined by the parallel
!| ion momentum equation driven by electromagnetic forces and ion pressure
!| gradient along the field lines\cite{weil92a}:
!| \begin{eqnarray}
!|   m_i n_i \Partial{v_{\parallel i}}{t}
!|  & = & - e Z_i n_i \left[ \nabla_{\parallel} \phi + \frac{1}{c} \left(
!|     \Partial{A_\parallel}{t}
!|      - ( {\mathbf v}_{*i} \times {\mathbf B}_{\perp} ) \cdot \hat{{\mathbf B}}
!|          \right) \right] \nonumber \\
!|  & &   - \nabla_{\parallel} ( n_i T_i )
!| \end{eqnarray}
!| where $A_\parallel$ is the parallel component of the vector potential and
!| $ {\mathbf B}_{\perp} = {\mathbf \nabla} \times
!|     ( A_\parallel \hat{{\mathbf B}} )
!|    \approx {\mathbf \nabla} A_\parallel \times  \hat{{\mathbf B}} $.
!| 
!| The ion energy balance equation is
!| \begin{eqnarray}
!|  \frac{3}{2} n_i \left( \Partial{ }{t}
!|   + {\mathbf v}_i \cdot {\mathbf \nabla} \right) T_i
!|   + n_i T_i {\mathbf \nabla} \cdot {\mathbf v}_i = - {\mathbf \nabla} \cdot {\mathbf q}_{*i} \nonumber \\
!|   \hspace{1em} =
!|   \frac{5}{2} n_i ( {\mathbf v}_{*i} - {\mathbf v}_{Di} ) \cdot {\mathbf \nabla} T_i
!| \end{eqnarray}
!| where $ {\mathbf q}_{*i} $ is the diamagnetic ion heat flow and
!| \[ {\mathbf v}_{Di} = \frac{T_i}{m_i \Omega_i} \hat{{\mathbf B}} \times
!|   \left( \frac{{\mathbf \nabla} B}{B}
!|   + \hat{\mathbf B} \cdot {\mathbf \nabla} \hat{\mathbf B} \right) \]
!| is the drift due to ${\mathbf \nabla} |{\mathbf B}|$ and magnetic curvature.
!| Note that the diamagnetic drift frequency is given by
!| \[ {\mathbf k} \cdot {\mathbf v}_{Di} = \omega_{Di}
!|    = \frac{-2 k_y T_i}{Z_i e B R}, \]
!| where $k_y$ is the poloidal wave number.
!| 
!| The electrons can be divided into two classes --- 
!| trapped electrons with density $n_{et}$ and fraction $f_t=n_{et}/n_e$, 
!| and free electrons with density $n_{ef}$ and fraction $ 1 - f_t $
!| \begin{equation}
!| n_e = n_{et} + n_{ef}.
!| \label{eq-trapped-free}
!| \end{equation}
!| The electron density $n_e$ is related to the density of hydrogenic ions
!| $n_H$, impurity ions $n_Z = f_Z n_e$, and superthermal hydrogenic ions
!| $n_s = f_s n_e$
!| through charge neutrality 
!| \begin{equation}
!| n_e = n_H + Z N_Z + Z_s n_s.
!| \label{eq-quasineutrality}
!| \end{equation}
!| 
!| The equations for trapped electron continuity and energy flow
!| are derived from a kinetic equation including electron-ion
!| collisions\cite{nils94a}
!| \begin{equation}
!|  \Partial{ }{t}n_{et} + {\mathbf \nabla} \cdot ( n_{et} {\mathbf v}_e )
!|   = - \hat{\nu} \omega_{De} ( \tilde{n}_{et}
!|         - \Gamma n_{et} e \tilde{\phi} / T_{e} )
!| \end{equation}
!| \begin{eqnarray}
!|  \frac{3}{2} \left( \Partial{ }{t}
!|   + {\mathbf v}_e \cdot {\mathbf \nabla} \right) (n_{et} T_{et} )
!|   + \frac{5}{2} n_{et} T_e {\mathbf \nabla} \cdot {\mathbf v}_e
!|   + {\mathbf \nabla} \cdot {\mathbf q} \nonumber \\
!|   \approx 1.5 \hat{\nu} \omega_{De} T_e ( \tilde{n}_{et}
!|         - n_{et} e \tilde{\phi} / T_e )
!| \end{eqnarray}
!| where
!| \[ \Gamma \equiv 1 +  \frac{g_{Te} }{ \omega / \omega_{De}
!|  + i \hat{\nu} - 1 } \hspace{1em} \mbox{and} \hspace{1em}
!| \hat{\nu} \equiv \frac{\nu_{ei}}{\omega_{De}} \frac{R}{r}. \]
!| The variable $ \hat{\nu} $
!| is the thermal electron-ion collision frequency $ \nu_{ei} $ normalized
!| by the electron diamagnetic drift frequency
!| $\omega_{De} = {\mathbf k} \cdot {\mathbf v}_{De} = 2 k_y T_e /  e B R $,
!| \[ {\mathbf q} = \frac{5}{2}
!|  \frac{n_{et} T_e \hat{{\mathbf B}} \times {\mathbf \nabla} T_e
!|    }{ m_e \Omega_e } \]
!| and the electron fluid velocity can be written
!| $ {\mathbf v}_e = {\mathbf v}_E +  {\mathbf v}_{*e}$.
!| 
!| A relation between the perturbed free (circulating) electron density
!| $ n_{ef} $ and the perturbed electric and magnetic potentials
!| can be obtained from the momentum equation for free electrons parallel
!| to the unperturbed magnetic field
!| \begin{eqnarray}
!|  m_e \frac{\partial v_{\parallel e}}{\partial t} & = &
!|    e \left( \hat{{\mathbf B}} \cdot {\mathbf \nabla}\phi
!|    + \frac{1}{c} \frac{\partial A_{\parallel}}{\partial t} \right)
!|    \nonumber \\
!| & & - \frac{e}{c} \left( {\mathbf v}_{*ep} \times \delta {\mathbf B}_{\perp}          \right) \cdot \hat{{\mathbf B}}
!|    - \frac{1}{n} \hat{{\mathbf e}}_{\parallel} \cdot {\mathbf \nabla} p
!|   \label{free-momentum}
!| \end{eqnarray}
!| where
!| the unperturbed perpendicular free electron velocity is the
!| diamagnetic velocity (including the full electron pressure gradient)
!| $ {\mathbf v}_{*ep} = - \hat{{\mathbf B}} \times {\mathbf \nabla}
!|    ( n_{ef} T_e ) / ( e n_{ef} B ) $.
!| The inertial term on the left of this momentum equation can be neglected
!| if $ k_\parallel v_{the} \gg \omega $.
!| 
!| The electron energy equaton is
!| \begin{equation}
!| \frac{3}{2} n_e \left( \frac{\partial}{\partial t}
!|   + {\mathbf v_e \cdot \nabla} \right) T_e
!|   + n_e T_e {\mathbf \nabla \cdot v_e}
!|   = {\mathbf - \nabla \cdot q_e}
!| \label{eq-free-electron-energy}
!| \end{equation}
!| where
!| \begin{equation}
!| {\mathbf q_e} = {\mathbf q_{*e}}
!|   - \kappa_\parallel \nabla_\parallel T_e
!| \end{equation}
!| 
!| Assuming the electron temperature is nearly uniform along
!| the perturbed magnetic field
!| $ {\mathbf B} \cdot {\mathbf \nabla} T_e \approx 0 $,
!| if follows that perturbed electron temperature is approximately
!| \begin{equation}
!|  \hat{T}_e =
!|    \frac{ \omega_{*ep} }{ k_{\parallel} c }
!|    \frac{ g_{Te} }{  g_{ne} }
!|       \hat{A}_{\parallel}
!| \end{equation}
!| where $ \hat{A}_\parallel \equiv e A_\parallel / T_e $.
!| The perturbed free electron density then follows by neglecting electron
!| inertia in Eq. (\ref{free-momentum})
!| \begin{equation}
!| \hat{n}_{ef} =
!|    \hat{\phi} - ( \omega - \omega_{*e} ) \hat{A}_\parallel
!|     / ( k_\parallel c ).
!| \end{equation}
!| 
!| In order to relate $\phi$ and $A_{\parallel}$, use the
!| electron continuity equation for the free electrons,
!| \begin{eqnarray}
!| \frac{\partial n_{ef}}{\partial t} & + &
!|   {\mathbf v}_{E}\cdot{\mathbf \nabla} n_{ef}
!|   + n_{ef}{\mathbf \nabla}\cdot{\mathbf v}_{E}  \nonumber \\
!|   & + & {\mathbf \nabla}\cdot \left( n_{ef} {\mathbf v}_{*e} \right)
!|   - \frac{1}{e}{\mathbf \nabla}\cdot
!|    \left( J_{\parallel ef}\hat{{\mathbf B}} \right)
!|    = 0
!| \end{eqnarray}
!| Only the free electrons contribute to the parallel electron current.
!| When parallel ion motion is included in the current,
!| Ampere's law implies,
!| \begin{equation}
!| J_{\parallel} = J_{\parallel ef} + J_{\parallel i} = - \frac{c}{4 \pi}
!|    \Delta A_{\parallel}
!| \end{equation}
!| Combining the last four equations,
!| eliminating $\tilde{n}_{ef}$, using the curvature relations
!| ${\mathbf \nabla} \cdot {\mathbf v}_E$ and
!|  ${\mathbf \nabla} \cdot \left( n_e {\mathbf v}_{*e} \right)$,
!| and assuming the same unperturbed temperature for the free and trapped
!| electrons, we obtain
!| \begin{equation}
!| \hat{A}_{\parallel}
!|  = \frac{k_{\parallel} c \left( \omega_{*e} - \omega \right) \hat{\phi}}
!| { \omega \left( \omega_{*e} - \omega \right)
!|    + \left( k_{\perp} \rho k_{\parallel} v_{A} \right)^2
!|    + \omega_{De} \left( \omega - \omega_{*eT} \right) }
!| \end{equation}
!| where $\omega_{*eT} = \omega_{De} ( g_{Te} + g_{ne} ) / 2 $.
!|  0.000000E+00qns (3) and (6) provide sufficient information for
!|  0eneralizing the previous Boltzmann relation for the free electrons.
!| 
!| \section{Linearized Perturbed Equations}
!| 
!| The perturbed equations are linearized and a given Fourier harmonic
!| is considered
!| \begin{equation}
!| u = u_0(x) \{ 1 + \hat{u} \exp [ - i ( \omega + {\mathbf k \cdot x} ) ] \}
!| \end{equation}
!| \[ \hat{u} \equiv \tilde{u} / u_0 \]
!| The perturbed densities and temperatures are normalized by their 
!| background values.
!| The perturbed electrostatic potential is normalized by $ T_e / e $.
!| All frequencies and growth rates are normalized by 
!| $\omega_{De} = \vec{k} \cdot \vec{v}_{De} = 2 k_y T_e /  e B R $. 
!| The parallel ion velocity is normalized by the speed
!| of sound in that ion, 
!| $ \hat{v}_{\parallel i} \equiv v_{\parallel i} / c_{si} $,
!| where $ c_{si} \equiv \sqrt{T_e / m_i} $.
!| In summary, the physical variables are normalized using
!| \begin{eqnarray*}
!| \hat{n} & \equiv & \tilde{n} / n \\
!| \hat{T} & \equiv & \tilde{T} / T \\
!| \hat{\phi} & \equiv & e \tilde{\phi} / T_e \\
!| \hat{\omega} & \equiv & \omega / \omega_{De} \\
!| \hat{v_{\parallel i}} & \equiv & v_{\parallel i} / c_{si} \\
!| \end{eqnarray*}
!| 
!| The equations for the divergence of the drifts are derived in
!| reference \cite{weil00a}
!| % The following equations for the divergence of the drifts are derived
!| 0y J.~Weiland in Chalmers University report CTH-IEFT/PP-1992-17 (1992)
!| % reference \cite{weil92b}:
!| \begin{equation}
!|  {\mathbf \nabla} \cdot \delta ( {n_i {\mathbf v}_{*i}} )
!|  = {\mathbf v}_{Di} \cdot {\mathbf \nabla} \delta ( {n_i T_i} ) / T_i 
!| \end{equation}
!| [Eq.~() on page ], 
!| % (page 132 Eq. (4.20)),
!| \begin{equation}
!|  {\mathbf \nabla} \cdot \tilde{{\mathbf v}}_E
!|    = \frac{Z e}{T_i} {\mathbf v}_{Di} \cdot {\mathbf \nabla} \tilde{\phi} 
!| \end{equation}
!| % (page 132 Eq. (4.21))
!| [Eq.~() on page ], and
!| \begin{equation}
!|  {\mathbf \nabla} \cdot [ n_i ( {\mathbf v}_{pi} + {\mathbf v}_{\pi i} ) ]
!|  \approx - i n_i k_y^2 \rho_{si}^2
!|  [ \omega - \omega_{Di} ( g_{ni}  + g_{Ti} ) / 2 ] \hat{\phi} 
!| \end{equation}
!| [Eq.~(2.49) on page 23],
!| % (page 25 Eq. (1.28))
!| where
!| $ \rho_{si}^2 = 2 T_e / ( m_i \Omega_i^2 ) $,
!| $ \hat{\phi} \equiv e \tilde{\phi} / T_e $,  and
!| \begin{equation}
!|   g_{ni} \equiv - R \hat{{\mathbf x}} \cdot {\mathbf \nabla} n_i / n_i
!|   \hspace{2em} g_{Ti} \equiv
!|  - R \hat{{\mathbf x}} \cdot {\mathbf \nabla} T_i / T_i
!|  \end{equation}
!| are the normalized gradients of ion temperature and density of
!| ion species $i$,
!| and $ \hat{{\mathbf x}} $ is a unit vector in the radial direction.
!| 
!| The perturbed ion equations for particles, momentum, and energy become
!| \begin{equation}
!| ( \hat{\omega} + \frac{T_i}{Z_i T_e} ) \hat{n}_i
!|   + \frac{T_i}{Z_i T_e} \hat{T}_i
!|   + \left[ 1 - \frac{g_{ni}}{2}
!|   + k_y^2 \rho_{si}^2 \left( \hat{\omega}
!|     + \frac{T_i}{Z_i T_e} \frac{g_{ni}+g_{Ti}}{2}
!|     \right) \right] \hat{\phi}
!|   - \frac{ c_{si}}{\omega_{De}} 
!|       {\hat{\mathbf B} \mathbf \cdot \nabla} \hat{v}_{\parallel i}
!|   = 0 
!| \end{equation}
!| \begin{equation}
!|  \hat{\omega} \hat{v}_{\parallel i}
!|  - \frac{ c_{si}}{\omega_{De}} Z_i
!|     {\hat{\mathbf B} \mathbf \cdot \nabla}  \left[ \hat{\phi}
!|    + \frac{T_i}{Z_i T_e} ( \hat{n}_i + \hat{T}_i ) \right] = 0 
!| \end{equation}
!| \begin{equation}
!|  \left( - \hat{\omega} - \frac{5}{3} \frac{T_i}{Z_i T_e}
!|    \right) \hat{T}_i + \frac{2}{3} \hat{\omega} \hat{n}_i
!|    + \frac{1}{2} \left( g_{Ti} - \frac{2}{3} g_{ni}
!|      \right) \hat{\phi} = 0 
!| \end{equation}
!| where  and
!| \begin{equation} 
!| g_{Ti} \equiv R / L_{Ti} = - R \hat{x} \cdot \nabla T_i / T_i 
!| \end{equation}
!| and
!| \begin{equation}
!| g_{ni} \equiv R / L_{ni} = - R \hat{x} \cdot \nabla n_i / n_i 
!| \end{equation}
!| is the normalized background temperature and density gradients.
!| 
!| 
!| 
!| 
!| %%::
!| 
!| Using Eq.~(\ref{eq-trapped-free}), the perturbed electron density
!| can be expressed in terms of the perturbed trapped electron density
!| $ \hat{n}_{et} $ and free electron density $ \hat{n}_{ef} $
!| \begin{equation}
!| \hat{n}_e = f_t \hat{n}_{et} + ( 1 - f_t ) \hat{n}_{ef}.
!| \label{eq-perturbed-trapped-free}
!| \end{equation}
!| Using quasineutrality, Eq.~(\ref{eq-quasineutrality}),
!| the perturbed electron density can be related to the perturbed
!| densities for hydrogenic ions $ \hat{n}_H \equiv \tilde{n}_H / n_H $
!| and impurity ions $ \hat{n}_Z \equiv \tilde{n}_Z / n_Z $
!| The electron density $n_e$ is related to the density of hydrogenic ions
!| $n_H$, impurity ions $n_Z = f_Z n_e$, and superthermal hydrogenic ions
!| $n_s = f_s n_e$
!| through charge neutrality 
!| \begin{equation}
!| \hat{n}_e = f_H \hat{n}_H + Z f_Z \hat{n}_Z
!| \label{eq-perturbed-quasineutrality}
!| \end{equation}
!| where 
!| \begin{equation}
!|  f_H \equiv n_H / n_e = 1 - Z f_Z - Z_s f_s
!| \end{equation}
!| is the fraction of hydrogenic ions
!| and it is assumed that superthermal ions do not take part in the 
!| perturbation, {\em i.e.} $ \hat{n}_s = 0.$
!| 
!| Note that quasi-neutrality imposes a constraint on the normalized
!| gradients of the background density
!| \begin{equation}
!| g_{ne} = f_H g_{ni} + Z f_Z g_{nZ} + Z_s f_s g_{ns}
!| \label{eq-gradient-quasineutrality}
!| \end{equation}
!| where the charge state of impurities $ Z $ and fast ions $ Z_s $
!| is included under the radial derivative in the definitions of the
!| normalized radial derivatives for impurities and fast ions
!| \[ g_{nZ} \equiv - R ( \partial Z n_Z / \partial r ) / ( Z n_Z ) \]
!| \[ g_{ns} \equiv - R ( \partial Z_s n_s / \partial r ) / ( Z_s n_s ). \]
!| Note that the impurity and fast ion charge state ($Z$ and $Z_s$)
!| are under the radial derivatives.
!| Actually, $ Z n_Z $ represents the sum over all the impurity 
!| densities times their charge states.
!| 
!| 
!| The basic technique is to set up the generalized eigenvalue problem
!| \[ A v = \lambda B v \]
!| where $ \lambda = \hat{\omega} + i \hat{\gamma} $ is the eigenvalue and $ v
!| $ is the corresponding eigenvector. Hence the eigenvalues give the frequency
!| and growth rates of the modes while the eigenvectors give the phase of the
!| perturbed variables relative to one another.
!| 
!| In order to approximate the reduction of transport with $E\times B$ velocity
!| shear the $E\times B $ shearing rate is subtracted from the linear growth
!| rates obtained as above. Currently this is implemented as follows: The
!| generalized eigenvalue problem is redefined using
!|  \[ (A - {\it i} \omega_{E\times B} B) v' = {\lambda}' B v' \]
!| and fluxes are calculated rom the new growth rates   ${\lambda}'$ and
!| eigenvectors $v'$. This method do give the same results as a direct
!| reduction of the growth rates only would give but is more integrated with
!| the current framework of the model.
!| 
!| 
!| The heat and particle fluxes are computed directly from the
!| eigenvectors in the following way:
!| Consider, for example, the flux of hydrogen particles produced by
!| the perturbed $ {\mathbf E} \times {\mathbf B} $ motion of the plasma
!| \[ \Gamma_H = 2 {\rm Re} ( \tilde{n}_H \tilde{v}_E^* )
!|  = 2 ( {\rm Re} \tilde{n}_H {\rm Im} \tilde{\phi}
!|      - {\rm Im} \tilde{n}_H {\rm Re} \tilde{\phi} ) k_y / B. \]
!| Using the quasi-linear assumption, the saturation level is
!| approximated by
!| \[ \hat{\phi} \equiv \frac{e \tilde{\phi}}{T_e}
!|  \approx \frac{1}{k_x \rho_{sH}} \frac{ \gamma }{ k_y c_{sH} }
!|  = \frac{2}{R k_x} \frac{\gamma}{\omega_{De}}. \]
!| Hence, the hydrogen particle flux is given by
!| \[ \frac{ \Gamma_H }{ n_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \frac{ \gamma^2 }{ \omega_{De}^2}
!|   \frac{ {\rm Re} \hat{n}_H {\rm Im} \hat{\phi}
!|        - {\rm Im} \hat{n}_H {\rm Re} \hat{\phi} }{ | \hat{\phi} |^2 }. \]
!| Correspondingly, the flux of heat flow through the hydrogen ions with
!| pressure $ n_H T_H $ is given by
!| \[ \frac{ \Gamma_{p_H} }{ n_H T_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \frac{ \gamma^2 }{ \omega_{De}^2}
!|   \frac{ {\rm Re} ( \hat{n}_H + \hat{T}_H ) {\rm Im} \hat{\phi}
!|   - {\rm Im} ( \hat{n}_H + \hat{T}_H ) {\rm Re} \hat{\phi} }{
!|     | \hat{\phi} |^2 }. \]
!| Four channels of transport are used in our simulations --- electron and
!| ion heat flux, as well as hydrogenic and impurity charged particle flux.
!| 
!| A diffusion matrix is computed by taking finite difference derivatives
!| of the fluxes with respect to the temperature and density gradients.
!| Convective velocities are computed from the difference between the total
!| fluxes and the diffusive fluxes.  There is no additional convective
!| contribution to the heat flux.
!| 
!| This model was never intended to compute the entire spectrum of the
!| turbulence with respect to poloidal mode number $k_y$.  As a working
!| assumption, the value $ k_y \rho_s = 0.316 $ was taken as a
!| representative point in the middle of the spectrum.
!| This assumption produces a gyro-Bohm transport model.
!| Based on calibrations against experimental data, the Weiland model
!| is multiplied by the coefficient $ 0.57 \kappa^{-4} $ when it is used
!| in the Multi-Mode model.
!| 
!| 
!| In this routine, the perturbed variables are always computed in the
!| following order:
!| \[ \hat{\phi}, \hat{T}_H, \hat{n}_H, \hat{T}_{et},
!|     \hat{n}_Z, \hat{T}_Z, \ldots \]
!| 
!| 
!| % \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Names of variables in the argument list:
!| \begin{tabular}{lllp{3.0in}}
!| \\ \hline
!| variable & status & symbol & meaning \\ \hline {\tt letain(j)} & input & &
!| Integer control variables
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
!| {\tt fnsin}  & input & $ Z_s f_s $
!|   & $ Z_s f_s = Z_s n_s / n_e $ fraction of superthermal hydrogenic ions \\
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
!| \\ \hline
!| variable & meaning \\ \hline {\tt letain(2)} & = number of elements computed
!| for transport matrix
!|                     (only when $> 0$) \\
!| {\tt letain(15)} & $ = 0 $ to redefine $A$ matrix with $\omega_{E\times B}$ \\
!|                  & $ = 1 $ to subtract shearing rate from growth rates directly
!| \\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| \hline \hline
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Real-valued control variables in the argument list:
!| \begin{tabular}{lp{4.0in}}
!| \\ \hline
!| variable & meaning \\ \hline
!| {\tt cetain(10)} & = coefficient of $ k_\parallel $ (default to 1.0) \\
!| {\tt cetain(11)} & = normalization for $ A_\parallel $ (default to 1.0) \\
!| {\tt cetain(12)} & = coefficient of $H$ (default to 0.0) \\
!| {\tt cetain(15)} & = coefficient of $ \hat{\nu} $ (default to 1.0) \\
!| {\tt cetain(20)} & = coefficient of $ \beta_{e,h,z} $ (default to 1.0) \\
!| {\tt cetain(25)} & = coefficient for elongation modified Alfven frequency
!| (defau
!| lt to 0.0) \\
!| {\tt cetain(29)} & = radius used in printouts \\
!| {\tt cetain(30)} & = finite difference used to construct
!|                    transport matrix (see zgm(j1,jd) matrix) \\
!| {\tt cetain(32)} & = tolerance used in eigenvalue solver \\
!| \hline \hline
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| Effects included with different numbers of equations:
!| \begin{tabular}{lp{5.0in}}
!| \\ \hline
!| neq & effects \\ \hline
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
!| \hline \hline
!| \end{tabular}
!| \end{center}
!| 
!| 
!| Note, $ \omega_{*e} = k_\perp \rho_s c_s / L_n $, $ \rho_s = c_s /
!| \omega_{ci} $, $ c_s = \sqrt{ 2 T_e / m_i} $, and $ \omega_{ci} = Z_i e B /
!| m_i $, $ \omega_{De} R = k_y \rho_s c_s $, $v_i = \sqrt{T_i / m_i}$, and $
!| L_n = - n / \Partial{n}{r} $ in the above normalizations.
!| 
!| The diffusivity matrix $ D = {\tt difthi(j1,j2)}$ and convective velocity
!| arrary $ v = {\tt velthi(j1)} $ are given in the following form:
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots \end{array} \right)
!|  = \nabla \cdot
!| \left( \begin{array}{l} 
!| n_H T_H  \\ n_H \\ n_e T_e \\ n_Z \\ n_Z T_Z \\ \vdots \end{array} \right)
!| \left( \begin{array}{llll}
!| D_{1,1} & D_{1,2} & D_{1,3} \\
!| D_{2,1} & D_{2,2} & D_{2,3} \\
!| D_{3,1} & D_{3,2} & D_{3,3} & \vdots \\
!| D_{4,1} & D_{4,2} & D_{4,3} \\
!| D_{5,1} & D_{5,2} & D_{5,3} \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!| \left( \begin{array}{l}  \nabla T_H / T_H \\
!|     \nabla n_H / n_H \\  \nabla T_e / T_e \\
!|     \nabla n_Z / n_Z \\  \nabla T_Z / T_Z \\ \vdots
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
!| Note that all the diffusivities in this routine are normalized by $
!| \omega_{De} / k_y^2 $, convective velocities are normalized by $ \omega_{De}
!| / k_y $, and all the frequencies are normalized by $ \omega_{De} $.
!| 
!| 
!| Note that
!| $$ \omega_{De} = 2 k_\perp T_e / e B R = 2 k_\perp \rho_s c_s / R $$
!| $$ D \propto \gamma / k^2
!|   \propto \rho_s^2 \omega_{De}
!|     \frac{ \gamma / \omega_{De} }{ k^2 \rho_s^2 }. $$
!| 
!| 
!| \begin{table}[h]
!| \caption{Local variables}
!| 
!| \begin{longtable}{lcl}
!| EN     & = & $ 2 / [ -R ( \partial n_e / \partial r ) / n_e ] $ \\
!| ENI    & = & $ -R ( \partial n_e / \partial r ) / n_e / 2 $ \\
!| EI     & = & $ g_{T_H} / g_{n_H} $ \\
!| EE     & = & $ g_{T_e} / g_{n_e} $ \\
!| TAU    & = & $ T_e / T_H $ \\
!| TAUI   & = & $ T_H / T_e $ \\
!| FL     & = & $ k_\theta^2 \rho_s^2 $ \\
!| RFL    & = & $ k_\theta \rho_s $ \\
!| GM     & = & $ 1 / ( 1 - f_t ) $ \\
!| BTA    & = & $ \frac{5 T_H}{T_e} ( \frac{1}{1 - f_t} + \frac{T_H}{T_e} $ \\
!| H1     & = & $ 4 q^2 k_\theta^2 \rho_s^2 $ \\
!| XT     & = & $ \frac{1}{1 + T_H / T_e} $ \\
!|  & = & $ $ \\
!|  & = & $ $ \\
!|  & = & $ $ \\
!|  & = & $ $ \\
!|  & = & $ $ \\
!| 
!| \end{longtable}
!| \end{table}
!| 
!| 
!| 
c
c@eqs.f
c
      subroutine weiland18eqns(
     &    lprintin, neq, gne, gnh, gnz, gte, gth, gtz,
     &    tauhin,tauzin,fnzin, zimp, zmass,zfns, betae, ftrapein,
     &    vef, QQ, shear, aspinv, kappa, ekyrhoin, hcn, wexb, ROT,
     &    ITL, IK, em, ITERA,
     &    zamr, zbmr, zami, zbmi,WZJ,WIMAX)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
C
      IMPLICIT none
      integer idp
      PARAMETER (idp=12)
c
      LOGICAL inital
      data inital /.true./
c
      REAL gne,gnh,gnz,gte,gth,gtz, zfns, shear, aspinv
      REAL epsnhin, epsnzin, epstein, epsthin, epstzin, tauhin, tauzin,
     &     fnzin, ftrapein, epsnin,ekyrhoin, etai, etae, bt, vef
      REAL KAPPA, GAV,ALA,SH2,hcn,He,BTA,RR
      REAL ftf,fzf,alff,WIMAX(100)
      REAL zamr(idp,idp),zami(idp,idp),zbmr(idp,idp),zbmi(idp,idp)
      REAL fft,fzft
      REAL zepsmach,zepsmin,zetah,zetaz,zetae,
     &     zepste,ztauh,ztauz,zep2nh,zep2nz,zep2ne,zimp,zfnz,
     &     zmass,zflh,zflz,em
      REAL QQ,alp,k1,k2,kps,betae,alf,kpc
      REAL WEXB, ROT

      INTEGER lprintin, neq, idim,
     &        ieq, j1, j2, j, IK
      INTEGER ITL, ITERA

      COMPLEX ALPC,ALPK,ALPHA,HQ,IU,WZJ(100)

      SAVE zepsmach
      SAVE inital, idim
C=======================================================================
C==
C DICTIONARY
C
C  INPUT
C
C    lprintin                             print flag
C    gne, gnh, gnz, gte, gth, gtz         Gradient Temperatures
C    neq                                  Number of equations
C    fnzin
C    zimp
C    zmass
C    zfns
C    betae
C    ftrapein
C    vef
C    QQ
C    shear
C    aspinv                               NEVER USED
C    kappa
C    ekyrhoin
C    hcn                                  used in 10 eqs.
C    wexb                                 NEVER USED
C    ROT
C    ITL                                  Iterations limit
C    IK
C    em                                   = 1 or 0 in 9, 10 eqs
C    ITERA                                Iterations counter
C
C  INTERMEDIATE VARIABLE NAMES
C
C    GAV                                  Used in 10 eqs
C    ALA                                  Used in 9,10 eqs
C    SH2                                  Used in 9,10 eqs (shear, kappa)
C    He                                   Used in 10 eqs.
C    BTA                                  Used in 9,10 eqs
C    RR                                   Used in 9,10 eqs
C    ftf, fzf, alff                       Used in 10 eqs
C    fft,fzft                             Used in 9,10 eqs
C    zepsmach                             machine epsilon
C    zepsmin                              = 1.e-10
C    zetah                                = etai
C    zetaz                                = epsnzin/epstzin
C    zetae                                = zep2ne/zepste
C    zepste                               = epstein
C    ztauh                                = 1./tauhin
C    ztauz                                = 1./(zimp*tauzin)
C    zep2nh                               = epsnhin
C    zep2nz                               = epsnzin
C    zep2ne                               = epsnin
C    zfnz                                 = fnzin*zimp
C    zflh                                 = ekyrhoin**2
C    zflz                                 = zmass*zflh/zimp**2
C    alp                                  Used in 8,9,10 eqs
C    k1,k2,kps                            Used in 8,9 eqs
C    alf                                  Used in 8,9,10 eqs
C    kpc                                  = 1.0 (in eqs 8,9,10)
C    idp                                  = 10
C    idim                                 = idp
C    ieq                                  = max(2, neq)
C    j1, j2, j                            Counters
C    ALPC,ALPK, ALPHA                     Used in 9, 10 eqs
C    HQ                                   Used in 9, 10 eqs
C    IU                                   = (0.,1.)
C    WZJ(100)                             Used in 9,10 eqs
C

C  OUTPUT
C
C    zamr, zbmr                           matrices - real part
C    zami, zbmi                           matrices - imaginary part
C    WIMAX(100)                           output(?)- used in 10 eqs
C    Xbest(90)
C
c
c ieq   = number of equations
c
c zamr(i,j) = matrix A
c zbmr(i,j) = matrix B
c
c    variables i=1,6: efi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
c    variables j=1,6 same as for i
C=====================================

c..initialize variables
c
      IF(inital) then
         idim     = idp
         zepsmach = 0.5
         do
           if (0.5*zepsmach+1.0 <= 1.0) exit
           zepsmach = 0.5*zepsmach
         end do
         inital = .false.
      ENDIF
c
      zepsmin = 1.e-10
c
c..variables derived from arguement list
c
      epsnin  =  2./sign(max(abs(gne),zepsmin),gne)
      epsnhin =  2./sign(max(abs(gnh),zepsmin),gnh)
      etai    = gth/sign(max(abs(gnh),zepsmin),gnh)
      etae    = gte/sign(max(abs(gne),zepsmin),gne)
      tauzin  =  1./sign(max(abs(tauzin),zepsmin),tauzin)
      epsnzin =  2./sign(max(abs(gnz), zepsmin),gnz)
      epstein =  2./sign(max(abs(gte), zepsmin),gte)
      epsthin =  2./sign(max(abs(gth), zepsmin),gth)
      epstzin =  2./sign(max(abs(gtz), zepsmin),gtz)
      tauhin  =  1./sign(max(abs(tauhin),zepsmin),tauhin)

      bt      =  1.5
c
c The following are never used
c
C     kxq     = gnz/ sign(max(abs(gne), zepsmin),gne)
C     kiq     = gnh/ sign(max(abs(gnz), zepsmin),gnz)
C     eq      = gtz/ sign(max(abs(gnz), zepsmin),gnz)
C====================
c
      IF ( lprintin .GT. 2 ) THEN
c
         write(6,*)
         write(6,*)' at beginning of weiland18eqns.f'
         write(6,*)
c
         write(6,*)gne,     ' = gne'
         write(6,*)epsnin,  ' = epsne'
         write(6,*)epsnhin, ' = epsnh'
         write(6,*)etai,    ' = etah'
         write(6,*)etae,    ' = etae'
         write(6,*)fnzin,   ' = fnz'
         write(6,*)tauzin,  ' = tauz'
         write(6,*)epsnzin, ' = epsnz'
         write(6,*)zfns,    ' = zfns'
         write(6,*)epstein, ' = epste'
         write(6,*)epsthin, ' = epsth'
         write(6,*)epstzin, ' = epstz'
         write(6,*)tauhin,  ' = tauh'
         write(6,*)zmass,   ' = zmass'
         write(6,*)ftrapein,' = ftrape'
         write(6,*)ekyrhoin,' = ekyrho'
         write(6,*)lprintin,' = lprint'
         write(6,*)bt,      ' = bt'
         write(6,*)vef,     ' = vef'
         write(6,*)betae,   ' = betae'
         write(6,*)QQ,       ' = qvalue'
         write(6,*)shear,   ' = shear'
         write(6,*)em,      ' = em'
!         write(6,*)alp,     ' = alp'
!         write(6,*)alf,     ' = alf'
!         write(6,*)kpc,     ' = kpc'
!         write(6,*)kps,     ' = kps'
!         write(6,*)GAV,     ' = GAV'
!         write(6,*)ALA,     ' = ALA'
         write(6,*)kappa,   ' = kappa'
         write(6,*)IK,      ' = IK'
         write(6,*)ITL,     ' = ITL'
         write(6,*)HCN,     ' = HCN'
!         write(6,*)He,      ' = He'
!         write(6,*)fft,     ' = fft'
!         write(6,*)fzft,    ' = fzft'
!         write(6,*)ftf,     ' = ftf'
!         write(6,*)fzf,     ' = fzf'
         write(6,*)ROT,     ' = ROT'
      ENDIF

      ieq  = max(2, neq)
      IU   = (0.,1.)
c
c..check validity of input data
c
      zepste = epstein
c
c..compute the rest of the dimensionless variables needed
c
      zetah  = etai
c  *******  NOTE THE INVERSE DEFINITION OF ZTAUH ! ******
      ztauh  = 1./tauhin
      zep2nh = epsnhin
      zep2ne = epsnin
      zetae  = zep2ne/zepste
      zflh   = ekyrhoin**2
      zfnz   = fnzin*zimp
c     zetaz  = zepsnz/zepstz
      zetaz  = epsnzin/epstzin
      ztauz  = 1./(zimp*tauzin)
      zep2nz = epsnzin
c     zep2nz = 2.0*zepsnz
      zflz   = zmass*zflh/zimp**2
c
c..diagnostic output
c
      IF(lprintin.GT.6)THEN
         write(6,*)
         write(6,*)'--------------------------------------'
         write(6,*)
         write(6,*)
         write(6,*)zetah,   ' = zetah'
         write(6,*)zetaz,   ' = zetaz'
         write(6,*)zetae,   ' = zetae'
         write(6,*)ztauh,   ' = ztauh'
         write(6,*)ztauz,   ' = ztauz'
         write(6,*)zep2nh,  ' = zep2nh'
         write(6,*)zep2nz,  ' = zep2nz'
         write(6,*)zep2ne,  ' = zep2ne'
         write(6,*)zimp,    ' = zimp'
         write(6,*)zmass,   ' = zmass'
         write(6,*)zfnz,    ' = zfnz'
         write(6,*)zflh,    ' = zflh'
         write(6,*)zflz,    ' = zflz'
         write(6,*)zepsmach,' = zepsmach'
      ENDIF
c
c..set matricies for eigenvalue equation
c
      do j1=1,idim
        do j2=1,idim
            zamr(j1,j2) = 0.0
            zami(j1,j2) = 0.0
            zbmr(j1,j2) = 0.0
            zbmi(j1,j2) = 0.0
        enddo
      enddo
C==================TEN EQUATIONS========================================
      IF(ieq.EQ.10)THEN
c
c  Ten  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta,  parallel motion of
c  impurities and edge physics
c  Equations for e phi/Te, Ti, net, nef, Tet, Tef, n_z, T_z, F, Av
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the
c  parallel magnetic vector potential.
c
         SH2   = 2.*shear-1.+(KAPPA*(shear-1.))**2
C==
C==   WZJ(IK) is the "initial guess"
C==   WIMAX(IK) is max. imag. part (FORCED to be > 0.001)
C==   WZJ(IK) stores the current value of WZ
C==
         ALA      = em*2.*QQ*QQ*betae*(1.+zetae+ztauh*(1.+zetah))/zep2nh
         ALPK     =  .5*SQRT(SH2)*SQRT(4.*QQ*QQ*zflh/(1.+ztauh)*zflh*
     &                (1.+ztauh*(1.+zetah)/(zep2nh*WZJ(IK))))
         IF(REAL(ALPK).LT.0.)ALPK=-ALPK
         ALPC     = -IU*ALPK
         ALPHA    = -IU*ABS(SQRT(SH2))*QQ*zflh
         ALPC     =  ABS(ALPHA/ALPC)*ALPC
         RR       =  2.*ABS(REAL(WZJ(IK)*ALPC))
         HQ       =  2.*ALPC/(4.*QQ*QQ*zflh)
         GAV      =  (1.+0.5*shear/RR)*EXP(-0.25/RR)
         GAV      =  GAV-0.5*ALA*(1.-EXP(-1./RR))
         IF(GAV.LT.0.001) GAV=0.001
cbate    GAV      =  GAV-(1.-1./QQ**2)/aspinv
         alp      =  0.5*RR
         IF(alp.LT.0.1) alp=0.1
         alff     =  alp/(2.D0*zflh*QQ*QQ*betae*(1.-zfnz-zfns))
         kps      =  0.5*SQRT(alp/zflh)/QQ
         kpc      =  1.D0
         He       =  ABS(SQRT(SH2))*SQRT(zflh)*SQRT(HCN*WIMAX(IK))
         ftf      =  (1.-ftrapein)/(1.-zfnz-zfns)
         fzf      =  zfnz/(1.-zfnz-zfns)
         alf      =  alff/ftf
c
      IF ( lprintin .GT. 2 ) THEN
c
        write(6,*)'=========='
c
         write(6,*)GAV,   ' = GAV'
         write(6,*)REAL(HQ),' = REAL(HQ)'
         write(6,*)em,    ' = em'
         write(6,*)betae, ' = betae'
         write(6,*)zetae,    ' = zetae'
         write(6,*)ztauh,  ' = ztauh'
         write(6,*)zep2nh,    ' = zep2nh'
         write(6,*)QQ,     ' = QQ'
         write(6,*)AIMAG(WZJ(IK)),' = WZI'
         write(6,*)ALPK,  ' = ALPK'
         write(6,*)ALPC,  ' = ALPC'
         write(6,*)ALPHA, ' = ALPHA'
         write(6,*)RR,     ' = RR'
         write(6,*)HQ,    ' = HQ'
         write(6,*)shear,     ' = shear'
         write(6,*)zflh,    ' = zflh'
         write(6,*)zimp,  ' = zimp'
         write(6,*)zflz,  ' = zflz'
         write(6,*)ztauz, ' = ztauz'
         write(6,*)zetaz, ' = zetaz'
         write(6,*)zep2nz,' = zep2nz'
c
      endif
c
c hydrogen density
c
         zamr(1,1) = -GAV+REAL(HQ)+(1.-zflh*ztauh*(1.+zetah))/zep2nh
         zami(1,1) =  AIMAG(HQ)
         zamr(1,2) =  (REAL(HQ)-GAV)*ztauh
         zami(1,2) =  ztauh*AIMAG(HQ)
         zamr(1,3) =  (REAL(HQ)-GAV)*ztauh*ftrapein/(1.-zfnz-zfns)
         zami(1,3) =  ztauh*AIMAG(HQ)*ftrapein/(1.-zfnz-zfns)
         zamr(1,4) =  (REAL(HQ)-GAV)*ztauh*ftf
         zami(1,4) =  ztauh*AIMAG(HQ)*ftf
         zamr(1,7) = -(REAL(HQ)-GAV)*ztauh*fzf
         zami(1,7) = -ztauh*AIMAG(HQ)*fzf
         zamr(1,10)= -em*ztauh*REAL(HQ)*(1.+zetah)/(kpc*zep2nh)
         zami(1,10)= -em*ztauh*AIMAG(HQ)*(1.+zetah)/(kpc*zep2nh)
c
         zbmr(1,1) =  zflh
         zbmr(1,3) =  ftrapein/(1.-zfnz-zfns)
         zbmr(1,4) =  ftf
         zbmr(1,7) = -fzf
         zbmr(1,10)=  em*REAL(HQ)/kpc
         zbmi(1,10)=  em*AIMAG(HQ)/kpc
c
c  hydrogen energy
c
         zamr(2,1) =  (zetah-2./3.)/zep2nh
         zamr(2,2) = -ztauh*5./3.
c
         zbmr(2,2) =  1.
         zbmr(2,3) = -(2./3.)*ftrapein/(1.-zfnz-zfns)
         zbmr(2,4) = -(2./3.)*ftf
         zbmr(2,7) =  (2./3.)*fzf
c
c trapped electron cont.
c
         zamr(3,1) = -1.+1./zep2ne
         zamr(3,3) =  1.
         zami(3,3) = -vef
         zamr(3,5) =  1.
         zami(3,9) =  vef
c
         zbmr(3,3) =  1.
c
c  free electon cont.
c
         zamr(4,1) =  1./zep2ne-1.
         zami(4,1) =  He
         zamr(4,4) =  1.
         zami(4,4) = -He
         zamr(4,6) =  1.
         zami(4,6) = -He
         zami(4,10)=  em*He*(1.+zetae)/(zep2ne*kpc)
c
         zbmr(4,4) =  1.
         zbmi(4,10)=  em*He/kpc
c

c  trapped electron energy
c
         zamr(5,1) =  (zetae-2./3.)/zep2ne
         zami(5,1) =  vef*2./3.*bt
         zami(5,3) =  vef*2./3.
         zamr(5,5) =  5./3.
         zami(5,9) = -5./3.*vef
c
         zbmr(5,3) = -2./3.
         zbmr(5,5) =  1.
c
c  free electron energy
c
         zamr(6,1) =  (zetae-2./3.)/zep2ne
         zamr(6,6) =  5./3.
         zami(6,6) = -1.06*He
         zami(6,10)=  em*1.06*He*zetae/(zep2ne*kpc)
c
         zbmr(6,4) = -2./3.
         zbmr(6,6) =  1.
c
c  impurity density
c
         zamr(7,1) = -GAV+zimp*REAL(HQ)/zmass+(1.-zflz*ztauz*
     &                   (1.+zetaz))/zep2nz
         zami(7,1) =  zimp*AIMAG(HQ)/zmass
         zamr(7,7) =  (REAL(HQ)*zimp/zmass-GAV)*ztauz
         zami(7,7) =  zimp*ztauz*AIMAG(HQ)/zmass
         zamr(7,8) =  (REAL(HQ)*zimp/zmass-GAV)*ztauz
         zami(7,8) =  zimp*ztauz*AIMAG(HQ)/zmass
         zamr(7,10)= -em*REAL(HQ)*zimp*ztauz*(1.+zetaz)
     &                 /(kpc*zep2nz*zmass)
         zami(7,10)= -em*AIMAG(HQ)*zimp*ztauz*(1.+zetaz)
     &                /(kpc*zep2nz*zmass)
c
         zbmr(7,1) =  zflz
         zbmr(7,7) =  1.
         zbmr(7,10)=  em*REAL(HQ)*zimp/(kpc*zmass)
         zbmi(7,10)=  em*AIMAG(HQ)*zimp/(kpc*zmass)
c
c  impurity energy
c
         zamr(8,1) =  (zetaz-2./3.)/zep2nz
         zamr(8,8) = -ztauz*5./3.
c
         zbmr(8,7) = -2./3.
         zbmr(8,8) =  1.
c
c  variable F
c
         zamr(9,1) =  zetae/zep2ne-1.
         zami(9,1) =  vef
         zamr(9,9) =  1.
         zami(9,9) = -vef
c
         zbmr(9,1) = -1.
         zbmr(9,9) =  1.
c
c  Amperes law
c
         fft  = (1.-zfnz)/(1.-ftrapein)
         fzft = zfnz/(1.-ftrapein)
c
         zamr(10,1) =  em*  REAL(HQ)*(1.+zimp*fzf/zmass)
         zami(10,1) =  em*(AIMAG(HQ)*(1.+zimp*fzf/zmass)-ftf*He)
         zamr(10,2) =  em*  REAL(HQ)*ztauh
         zami(10,2) =  em* AIMAG(HQ)*ztauh
         zamr(10,3) =  em*  REAL(HQ)*ztauh*ftrapein/(1.-zfnz-zfns)
         zami(10,3) =  em* AIMAG(HQ)*ztauh*ftrapein/(1.-zfnz-zfns)
         zamr(10,4) =  em*  REAL(HQ)*ztauh*ftf
         zami(10,4) =  em*He*ftf+em*AIMAG(HQ)*ztauh*ftf
         zami(10,6) =  em*He*ftf
         zamr(10,7) = -em* REAL(HQ)*fzf*(ztauh-zimp*ztauz/zmass)
         zami(10,7) = -em*AIMAG(HQ)*fzf*(ztauh-zimp*ztauz/zmass)
         zamr(10,8) =  em* REAL(HQ)*fzf*zimp*ztauz/zmass
         zami(10,8) =  em*AIMAG(HQ)*fzf*zimp*ztauz/zmass
         zamr(10,10)= -em*(zflh*(1.+0.25*SH2/alp)*alff+REAL(HQ)*(
     &                 (1.+zetah)*ztauh/zep2nh+fzf*zimp*ztauz*(1.+zetaz)
     &                  /(zmass*zep2nz)))/kpc
         zami(10,10)= -em*(AIMAG(HQ)*((1.+zetah)*ztauh/zep2nh+
     &                  fzf*zimp*ztauz*(1.+zetaz)/(zmass*zep2nz))+ftf*
     &                  He*(1.+zetae)/zep2ne)/kpc
c
         zbmr(10,10) = em*REAL(HQ)*(1.+zimp/zmass)/kpc
         zbmi(10,10) = em*(AIMAG(HQ)*(1.+zimp/zmass)-He*ftf)/kpc
c
        return
c
      ENDIF
C========================END=OF=TEN=EQUATIONS===========================
C==================TWO EQUATIONS========================================
      IF(ieq.EQ.2)THEN
c
c..two equations when trapped particles and FLR effects omitted
c
         zamr(1,1) =  (1.0/zep2nh)-ztauh-1.0
         zamr(1,2) = -ztauh
         zamr(2,1) =  (zetah-2./3.)/zep2nh
         zamr(2,2) = -ztauh*5./3.
c
         zbmr(1,1) =  1.0
         zbmr(1,2) =  0.0
         zbmr(2,1) = -2./3.
         zbmr(2,2) =  1.0
      ENDIF
C==================FOUR EQUATIONS=======================================
      IF(ieq.EQ.4)THEN
c
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
      IF(lprintin.GT.5)write(6,'(/,A)')
     &        ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
c
c  ion continuity
c
         zamr(1,1) =  1.0-zep2nh-zflh*ztauh*(1.0+zetah)
         zamr(1,2) = -zep2nh*ztauh
         zamr(1,3) = -zep2nh*ztauh
c
         zbmr(1,1) =  zflh*zep2nh
         zbmr(1,3) =  zep2nh
c
c  ion energy
c
         zamr(2,1) =  zetah-2./3.
         zamr(2,2) = -zep2nh*ztauh*5./3.
c
         zbmr(2,2) =  zep2nh
         zbmr(2,3) = -zep2nh*2./3.
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation. The dilution factor 1-zfnz has now
c   been added.
c
         zamr(3,1) =  ftrapein-zep2ne
         zamr(3,3) =  zep2ne*(1.-zfnz-zfns)
         zamr(3,4) =  ftrapein*zep2ne
c
         zbmr(3,1) =  (ftrapein-1.0)*zep2ne
         zbmr(3,3) =  zep2ne*(1.-zfnz-zfns)
c
c  trapped electron energy
c
         zamr(4,1) =  ftrapein*( zetae-2./3.)
         zamr(4,4) =  ftrapein*zep2ne*5./3.
c
         zbmr(4,1) =  (1.0-ftrapein)*zep2ne*2./3.
         zbmr(4,3) = -zep2ne*2./3.
         zbmr(4,4) =  ftrapein*zep2ne
      ENDIF
C==================SIX EQUATIONS========================================
      IF(ieq.EQ.6)THEN
c
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z
      IF(lprintin.GT.5)write(6,'(/,A)')
     &   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
c
c  hydrogen density
c
         zamr(1,1) = -1.0+( 1.0-zflh*ztauh*( 1.0+zetah))/zep2nh
         zamr(1,2) = -ztauh
         zamr(1,3) = -ztauh
c
         zbmr(1,1) =  zflh
         zbmr(1,3) =  1.0
c
c  hydrogen energy
c
         zamr(2,1) =  (zetah-2./3.)/zep2nh
         zamr(2,2) = -ztauh*5./3.
c
         zbmr(2,2) =  1.0
         zbmr(2,3) = -2./3.
c
c  trapped electron density
c
         zamr(3,1) = -1.0+ftrapein/zep2ne
         zamr(3,3) =  1.0-zfnz-zfns
         zamr(3,4) =  ftrapein
         zamr(3,5) =  zfnz
c
         zbmr(3,1) =  ftrapein-1.0
         zbmr(3,3) =  1.0-zfnz-zfns
         zbmr(3,5) =  zfnz
c
c  trapped electron energy
c
         zamr(4,1) =  ftrapein*( zetae-2./3.)/zep2ne
         zamr(4,4) =  ftrapein*5./3.
c
         zbmr(4,1) =  (1.0-ftrapein)*2./3.
         zbmr(4,3) = -(1.0-zfnz-zfns)*2./3.
         zbmr(4,4) =  ftrapein
         zbmr(4,5) = -zfnz*2./3.
c
c  impurity density
c
         zamr(5,1) = -1.0+( 1.0-zflz*ztauz*(1.0+zetaz))/zep2nz
         zamr(5,5) = -ztauz
         zamr(5,6) = -ztauz
c
         zbmr(5,1) =  zflz
         zbmr(5,5) =  1.0
c
c  impurity energy
c
         zamr(6,1) =  (zetaz-2./3.)/zep2nz
         zamr(6,6) = -ztauz*5./3.
c
         zbmr(6,5) = -2./3.
         zbmr(6,6) =  1.0
      ENDIF
C==================SEVEN EQUATIONS======================================
      IF(ieq.EQ.7)THEN
c
c..Seven equations with impurities, trapped electrons, parallel ion motion
c and FLR
c
c  equations for e \phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here F is defined as F = GM*e phi/T_e where GM=1+etae/(epsn*(omega-1+i*vef))
c
        IF(lprintin.GT.5)write(6,'(/,A)')
     &    'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and Vp'
c
c  hydrogen density
c
        zamr(1,1) = -1.+(1.-zflh*ztauh*(1.+zetah))/zep2nh
        zami(1,1) = -0.5*ABS(shear)/QQ
        zamr(1,2) = -ztauh
        zami(1,2) = -ztauh*0.5*ABS(shear)/QQ
        zamr(1,3) = -ztauh
        zami(1,3) = -ztauh*0.5*ABS(shear)/QQ
c
        zbmr(1,1) =  zflh
        zbmr(1,3) =  1.
c
c hydrogen energy
c
        zamr(2,1) =  (zetah-2./3.)/zep2nh
        zamr(2,2) = -ztauh*5./3.
c
        zbmr(2,2) =  1.
        zbmr(2,3) = -2./3.
c
c trapped electron density
c
        zamr(3,1) = -1.+ftrapein/zep2ne
        zamr(3,3) =  1.-zfnz-zfns
        zamr(3,4) =  ftrapein
        zamr(3,5) =  zfnz
c
        zbmr(3,1) =  ftrapein-1.
        zbmr(3,3) =  1.-zfnz-zfns
        zbmr(3,5) =  zfnz
c
c trapped electron energy
c
        zamr(4,1) =  ftrapein*(zetae-2./3.)/zep2ne
        zamr(4,4) =  ftrapein*5./3.
c
        zbmr(4,1) =  (1.-ftrapein)*2./3.
        zbmr(4,3) = -(1.-zfnz-zfns)*2./3.
        zbmr(4,4) =  ftrapein
        zbmr(4,5) = -zfnz*2./3.
c
c impurity density
c
        zamr(5,1) = -1.+( 1.-zflz*ztauz*(1.+zetaz))/zep2nz
        zami(5,1) = -zimp*0.5*ABS(shear)/QQ/zmass
        zamr(5,5) = -ztauz
        zami(5,5) = -zimp*ztauz*0.5*ABS(shear)/QQ/zmass
        zamr(5,6) = -ztauz
        zami(5,6) = -zimp*ztauz*0.5*ABS(shear)/QQ/zmass
c
        zbmr(5,1) =  zflz
        zbmr(5,5) =  1.
c
c impurity energy
c
        zamr(6,1) =  (zetaz-2./3.)/zep2nz
        zamr(6,6) = -ztauz*5./3.
c
        zbmr(6,5) = -2./3.
        zbmr(6,6) =  1.
      ENDIF
C==================EIGHT EQUATIONS======================================
      IF(ieq.EQ.8)THEN
c
c..Eight equations with impurities, trapped electrons, parallel ion
c  motion collisions and FLR
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Vp
c  Here F = omega*gamma (collisions).
c
         k1  =.25*QQ*QQ*zflh*SQRT(ABS((1.+zetah)*ztauh/
     &          ((1.-ftrapein)*zep2nh)))
         k2  =  QQ*QQ*zflh*zflh*(1.+zetah)*ztauh/((1.-ftrapein)*zep2nh)
         alp =  0.5*(k1+SQRT(k1+shear*shear*k2))
         alf =  alp/(2.D0*zflh*QQ*QQ*betae*(1.D0-ftrapein))
         kps =  0.5*SQRT(alp/zflh)/QQ
         kpc =  1.D0
c
c
c hydrogen density
c
         zamr(1,1) = -1.0+( 1.0-zflh*ztauh*( 1.0+zetah))/zep2nh
         zamr(1,2) = -ztauh
         zamr(1,3) = -ztauh
         zamr(1,8) =  kps
c
         zbmr(1,1) =  zflh
         zbmr(1,3) =  1.0
c
c  hydrogen energy
c
         zamr(2,1) =  (zetah-2./3.)/zep2nh
         zamr(2,2) = -ztauh*5./3.
c
         zbmr(2,2) =  1.0
         zbmr(2,3) = -2./3.
c
c  total electron density expressed in ion and imp densities
c
         zamr(3,1) = -1.0+ftrapein/zep2ne
         zami(3,1) =  vef*(1.-ftrapein)
         zamr(3,3) =  1.0-zfnz-zfns
         zami(3,3) = -vef*(1.-zfnz-zfns)
         zamr(3,4) =  ftrapein
         zamr(3,5) =  zfnz
         zami(3,5) = -vef*zfnz
         zami(3,7) =  vef*ftrapein
c
         zbmr(3,1) =  ftrapein-1.0
         zbmr(3,3) =  1.0-zfnz-zfns
         zbmr(3,5) =  zfnz
c
c  trapped electron energy
c
         zamr(4,1) =  ftrapein*(zetae-2./3.)/zep2ne
         zami(4,1) =  vef*(2./3.)*(bt-2.5*(1.-ftrapein))
         zami(4,3) = -vef*(2./3.)*(bt-2.5)*(1.-zfnz-zfns)
         zamr(4,4) =  ftrapein*5./3.
         zami(4,5) = -vef*(2./3.)*(bt-2.5)*zfnz
         zami(4,7) = -5./3.*vef*ftrapein
c
         zbmr(4,1) =  (1.0-ftrapein)*2./3.
         zbmr(4,3) = -(1.0-zfnz-zfns)*2./3.
         zbmr(4,4) =  ftrapein
         zbmr(4,5) = -zfnz*2./3.
c
c  impurity density
c
         zamr(5,1) = -1.0+( 1.0-zflz*ztauz*( 1.0+zetaz))/zep2nz
         zamr(5,5) = -ztauz
         zamr(5,6) = -ztauz
c
         zbmr(5,1) =  zflz
         zbmr(5,5) =  1.0
c
c  impurity energy
c
         zamr(6,1) =  (zetaz-2./3.)/zep2nz
         zamr(6,6) = -ztauz*5./3.
c
         zbmr(6,5) = -2./3.
         zbmr(6,6) =  1.0
c
c  variable F
c
         zamr(7,1) =  zetae/zep2ne-1.
         zami(7,1) =  vef
         zamr(7,7) =  1.
         zami(7,7) = -vef
c
         zbmr(7,1) = -1.
         zbmr(7,7) =  1.
c
c  Parallel ion motion Vpi/Cs
c
         zamr(8,1) =  kps
         zamr(8,2) =  kps*ztauh
         zamr(8,3) =  kps*ztauh
c
         zbmr(8,8) =  1.
      ENDIF
C==================NINE EQUATIONS=======================================
      IF(ieq.EQ.9)THEN
c
c..Nine  equations with impurities, trapped electrons, parallel ion
c  motion, collisions,  FLR , finite beta and parallel motion of impurities
c
c  Equations for e phi/Te, Ti, ni, Tet, n_z, T_z, F, Av, K
c  Here F = omega*gamma (collisions), K = omega*Av and Av is the parallel
c   magnetic vector potential.
c
         SH2  =  2.*shear-1.+(KAPPA*(shear-1.))**2
         ALPHA= -IU*ABS(SQRT(SH2))*QQ*zflh
         RR    =  2.*ABS(REAL(WZJ(IK)*ALPHA))
         ALA   =  em*QQ*QQ*betae*(1.+zetae+ztauh*(1.+zetah))/zep2nh
         ALPK  =  0.5*SQRT(SH2)*SQRT(4.*QQ*QQ*zflh/(1.+ztauh)*zflh*
     &                 (1.+ztauh*(1.+zetah)/(zep2nh*WZJ(IK))))
         IF(REAL(ALPK).LT.0.)ALPK=-ALPK
         ALPC  = -IU*ALPK
         ALPC  =  ABS(ALPHA/ALPC)*ALPC
         RR    =  2.*ABS(REAL(WZJ(IK)*ALPC))
         HQ    =  2.*ALPC/4.*QQ*QQ*zflh
  802    GAV   =  (1.+0.5*shear/RR)*EXP(-0.25/RR)
         GAV   =  GAV-0.5*ALA*(1.-EXP(-1./RR))
         IF(GAV.LT.0.001) GAV=0.001
         GAV   =  GAV-(1.-1./QQ**2)/aspinv
c
c        k1    =  .25*QQ*QQ*zflh*SQRT((1.+zetah)*ztauh/
c    &            ((1.-ftrapein)*zep2nh))
c        k2    =  QQ*QQ*zflh*zflh*(1.+zetah)*ztauh/((1.-ftrapein)*zep2nh)
c        alp   =  0.5*(k1+SQRT(k1+shear*shear*k2))
         alp   =  0.5*RR
         IF(alp.LT.0.1) alp=0.1
         alf   =  alp/(2.D0*zflh*QQ*QQ*betae*(1.D0-ftrapein))
         kps   =  0.5*SQRT(alp/zflh)/QQ
         kpc   =  1.D0
c
c hydrogen density
c
         zamr(1,1) = -GAV+REAL(HQ)+(1.0-zflh*ztauh*
     &                 (1.0+zetah))/zep2nh
         zami(1,1) =  AIMAG(HQ)
         zamr(1,2) =  (REAL(HQ)-GAV)*ztauh
         zami(1,2) =  ztauh*AIMAG(HQ)
         zamr(1,3) =  (REAL(HQ)-GAV)*ztauh
         zami(1,3) =  ztauh*AIMAG(HQ)
         zamr(1,8) = -em*ztauh*REAL(HQ)*(1.+zetah)/(kpc*zep2nh)
         zami(1,8) = -em*ztauh*AIMAG(HQ)*(1.+zetah)/(kpc*zep2nh)
         zamr(1,9) = -em*REAL(HQ)/kpc
         zami(1,9) = -em*AIMAG(HQ)/kpc
c
         zbmr(1,1) =  zflh
         zbmr(1,3) =  1.
c
c  hydrogen energy
c
         zamr(2,1) =  (zetah-2./3.)/zep2nh
         zamr(2,2) = -ztauh*5./3.
c
         zbmr(2,2) =  1.
         zbmr(2,3) = -2./3.
c
c  total electron density expressed in ion density and imp density
c
         zamr(3,1) = -1.+ftrapein/zep2ne
         zami(3,1) =  vef*(1.-ftrapein)
         zamr(3,3) =  1.-zfnz-zfns
         zami(3,3) = -vef*(1.-zfnz-zfns)
         zamr(3,4) =  ftrapein
         zamr(3,5) =  zfnz
         zami(3,5) = -vef*zfnz
         zami(3,7) =  vef*ftrapein
         zamr(3,8) = -em*(1.-ftrapein)/(kpc*zep2ne)
         zami(3,8) =  em*(1.-ftrapein)*vef/(kpc*zep2ne)
         zamr(3,9) =  em*(1.-ftrapein)*(1.+1./zep2ne)/kpc
         zami(3,9) = -em*(1.-ftrapein)*vef/kpc
c
         zbmr(3,1) =  ftrapein-1.
         zbmr(3,3) =  1.-zfnz-zfns
         zbmr(3,5) =  zfnz
         zbmr(3,9) =  em*(1.-ftrapein)/kpc
c
c  trapped electron energy
c
         zamr(4,1) =  ftrapein*(zetae-2./3.)/zep2ne
         zami(4,1) =  vef*(2./3.)*(bt-2.5*(1.-ftrapein))
         zami(4,3) = -vef*(2./3.)*(bt-2.5)*(1.-zfnz-zfns)
         zamr(4,4) =  ftrapein*5./3.
         zami(4,5) = -vef*(2./3.)*(bt-2.5)*zfnz
         zami(4,7) = -5./3.*vef*ftrapein
c
         zbmr(4,1) =  (1.-ftrapein)*(2./3.)
         zbmr(4,3) = -(1.-zfnz-zfns)*(2./3.)
         zbmr(4,4) =  ftrapein
         zbmr(4,5) = -zfnz*2./3.
c
c  impurity density
c
         zamr(5,1) = -GAV +zimp*REAL(HQ)/zmass+(1.-zflz*ztauz*
     &                    (1.+zetaz))/zep2nz
         zami(5,1) =  zimp*AIMAG(HQ)/zmass
         zamr(5,5) =  (REAL(HQ)*zimp/zmass-GAV)*ztauz
         zami(5,5) =  zimp*ztauz*AIMAG(HQ)/zmass
         zamr(5,6) =  (REAL(HQ)*zimp/zmass-GAV)*ztauz
         zami(5,6) =  zimp*ztauz*AIMAG(HQ)/zmass
         zamr(5,8) = -em*REAL(HQ)*zimp*ztauz*(1.+zetaz)/
     &                  (kpc*zep2nz*zmass)
         zami(5,8) = -em*AIMAG(HQ)*zimp*ztauz*(1.+zetaz)/
     &                  (kpc*zep2nz*zmass)
         zamr(5,9) = -em*REAL(HQ)*zimp/(kpc*zmass)
         zami(5,9) = -em*AIMAG(HQ)*zimp/(kpc*zmass)
c
         zbmr(5,1) =  zflz
         zbmr(5,5) =  1.
c
c  impurity energy
c
         zamr(6,1) =  (zetaz-2./3.)/zep2nz
         zamr(6,6) = -ztauz*5./3.
c
         zbmr(6,5) = -2./3.
         zbmr(6,6) =  1.
c
c  variable F
c
         zamr(7,1) =  zetae/zep2ne-1.
         zami(7,1) =  vef
         zamr(7,7) =  1.
         zami(7,7) = -vef
c
         zbmr(7,1) = -1.
         zbmr(7,7) =  1.
c
c
c  electromagnetic parallel vectorpotential Av = e A_par/Te
c
         fft   = (1.-zfnz)/(1.-ftrapein)
         fzft  = zfnz/(1.-ftrapein)
         zamr(8,1) =  em*kpc*(1./zep2ne+REAL(HQ)*(fft+zimp*fzft/zmass))
         zami(8,1) =  em*AIMAG(HQ)*(fft+zimp*fzft/zmass)*kpc
         zamr(8,2) =  em*REAL(HQ)*ztauh*fft*kpc
         zami(8,2) =  em*AIMAG(HQ)*ztauh*fft*kpc
         zamr(8,3) =  em*REAL(HQ)*ztauh*fft*kpc
         zami(8,3) =  em*AIMAG(HQ)*ztauh*fft*kpc
         zamr(8,5) =  em*REAL(HQ)*zimp*ztauz*fzft*kpc/zmass
         zami(8,5) =  em*AIMAG(HQ)*zimp*ztauz*fzft*kpc/zmass
         zamr(8,6) =  em*REAL(HQ)*zimp*ztauz*fzft*kpc/zmass
         zami(8,6) =  em*AIMAG(HQ)*zimp*ztauz*fzft*kpc/zmass
         zamr(8,8) =  em*((1.+zetae)/zep2ne-alf*zflh*(1.+0.25*SH2/alp))
     &               -em*REAL(HQ)*(fft*ztauh*(1.+zetah)/zep2nh
     &               +zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*zmass))
         zami(8,8) = -em*AIMAG(HQ)*(fft*ztauh*(1.+zetah)/zep2nh
     &               +zimp*fzft*ztauz*(1.+zetaz)/(zep2nz*zmass))
         zamr(8,9) = -em*(1./zep2ne+REAL(HQ)*(fft+zimp*fzft/zmass))
         zami(8,9) = -em*AIMAG(HQ)*(fft+zimp*fzft/zmass)
c
         zbmr(8,1) =  em*kpc
         zbmr(8,8) =  em
         zbmr(8,9) = -em
c
         zamr(9,9) = em
c
         zbmr(9,8) = em
      ENDIF
C=======================================================================

      RETURN
      END
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{weil00a} Jan Weiland,
!| {\em Collective Modes in Inhomogeneous Plasma, Kinetic and Advanced
!| Fluid Theory,}
!| Institue of Physics Publishing, Bristol, UK (2000)
!| 
!| \bibitem{weil92a} J. Weiland,
!| ``Low Frequency Modes Associated with Drift Motions in Inhomogeneous
!| Plasmas,'' report CTH--IEFT/PP-1992-17, Institute for Electromagnetic Field
!| Theory and Plasma Physics, Chalmers University of Technology, G\"{o}teborg,
!| Sweden.
!| 
!| \bibitem{froj92a} M. Fr\"{o}jdh, M. Liljestr\"{o}m, H. Nordman,
!| ``Impurity effects on $\eta_i$ mode stability and transport,'' Nuclear
!| Fusion {\bf 32} (1992) 419--428.
!| \bibitem{nord92a} H. Nordman and J. Weiland, ``Comments on
!| `Ion-temperature-gradient-driven transport in a density modification
!| experiment on the tokamak fusion test reactor [Phys. Fluids {\bf B4} (1992)
!| 953]' ''.
!| 
!| \bibitem{nils94a}
!| J.~Nilsson and J.~Weiland,
!| \newblock Nucl. Fusion {\bf 34}, 803 (1994).
!| 
!| \bibitem{weil92b} J. Weiland,
!| ``Nonlinear effects in velocity space and drift wave transport in
!| tokamaks,'' Phys. Fluids {\bf B4} (1992) 1388--1390.
!| \bibitem{weil92c} J. Weiland and H. Nordman, ``Drift wave model for inward
!| energy transport in tokamak plasmas,'' Institute for Electromagnetic Field
!| Theory and Plasma Physics, 
!| Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.
!| \bibitem{weil92d} J.Weiland and A. Hirose, ``Electromagnetic and kinetic
!| effects on the ion temperature gradient mode,'' 
!| Nucl. Fusion {\bf 32} (1992) 151--155.
!| \bibitem{nord91a} H. Nordman and J. Weiland, ``The concept of marginal
!| stability and recent experimental results from the TFTR tokamak,'' Institute for Electromagnetic Field Theory and Plasma Physics, 
!| Gothenburg, Sweden, (1991) CTH-IEFT/PP-1991-26 ISSN.
!| \bibitem{weil91a} J. Weiland and H. Nordman, 
!| ``Enhanced confinement regimes in
!| transport code simulations of toroidal drift wave transport,'' 
!| Nucl. Fusion {\bf 31} (1991) 390--394.
!| \bibitem{nord90a} H. Nordman, J. Weiland, and A. Jarmen, 
!| ``Simulation of toroidal drift mode turbulence driven by temperature
!| gradients and electron trapping,'' 
!| Nucl. Fusion {\bf 30} (1990) 983--996.
!| \bibitem{weil89a} J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, 
!| ``Diffusive particle and heat pinch effects in toroidal plasmas,'' 
!| Nucl. Fusion {\bf 29} (1989) 1810--1814.
!| \bibitem{nord89a} H. Nordman and J. Weiland, 
!| ``Transport due to toroidal $\eta_i$ mode turbulence in tokamaks,'' 
!| Nucl. Fusion {\bf 29} (1989) 251--263.
!| \bibitem{ande88a} P. Andersson and J. Weiland, 
!| ``A fully toroidal fluid analysis
!| of the magnetohydrodynamic ballooning mode branch in tokamaks,'' 
!| Phys. Fluids {\bf 31} (1988) 359--365.
!| \bibitem{jarm87a} A. Jarm\'{e}n, P. Andersson, and J. Weiland, 
!| ``Fully toroidal ion temperature gradient driven drift modes,'' 
!| Nucl. Fusion {\bf 27} (1987) 941--949.
!| \end{thebibliography}
!| 
!| \end{document}

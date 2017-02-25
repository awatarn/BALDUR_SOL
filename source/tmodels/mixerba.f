!| %
!| %  To extract the fortran source code, obtain the utility "xtverb" from
!|  0lenn Bateman (bateman@plasma.physics.lehigh.edu) and type:
!| 0tverb < mixerba.tex > mixerba.f
!| %
!| \documentstyle{article}    % Specifies the document style.
!| 
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| 
!| \begin{document}
!| \begin{center}
!| {\LARGE Subroutine for Computing Particle and Energy Fluxes\\ \vskip8pt
!| Using the Mixed Transport Model}\vskip1.0cm
!| Version 1.0: January 18, 1999 \\
!| Implemented by M. Erba, G. Bateman, Arnold H. Kritz\\
!|  Lehigh University
!| \end{center}
!| For questions about this routine, please contact: \\
!| Matteo Erba, Lehigh:  {\tt erba@plasma.physics.lehigh.edu}\\
!| Arnold Kritz, Lehigh: {\tt kritz@plasma.physics.lehigh.edu}
!| Glenn Bateman, Lehigh: {\tt glenn@plasma.physics.lehigh.edu}
!| 
!| 
c@mixed.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine mixerba (
     &   shear,       wexb,       rmaj,       cwexb
     & , chibohm,     zq,         zrlpe,      tekev,       tikev
     & , gradte,      zra,        zdte,       btor
     & , chii,        chie,       dh,         dimp,        ierr
     & , chbe,        chbi,       chgbe,      chgbi
     & , zcoef1,      zcoef2,     zcoef3,     zcoef4)
c
c Mixed Transport Model (by M.Erba, V.V.Parail, A.Taroni)
c
c Inputs:
c    chibohm:   Bohm diffusivity, defined as T_e/eB, in MKS units.
c    tekev:     Electron temperature in kev
c    gradte:    Electron Temperature Gradient in keV/m
c    btor:      Toroidal field in Tesla
c    rmaj:      Major Radius
c         *** all other inputs are dimensionless ***
c    zq:        The local value of q, the safety factor
c    zrlpe:      R/L_pe, where R is the major radius of the plasma,
c                    and L_pe is the local electron pressure scale length
c    zra:       normalized minor radius r/a
c    zdte:       Non-local dependence on edge electron temperature:
c               abs{[Te(0.8)-Te(1.0)]/Te(1.0)]}
c    zcoef1:    Coefficient for empirical hydrogen diffusivity
c    zcoef2:    Coefficient for empirical hydrogen diffusivity
c    zcoef3:    Coefficient for empirical impurity diffusivity
c    zcoef4:    Coefficient for empirical impurity diffusivity
c    cwexb:     Multiplier of wexb
c
c Outputs:
c    ierr:      Error code
c    chii:      The ion thermal diffusivity, in chibohm's units
c    chie:      The electron thermal diffusivity, in chibohm's units
c    chgbe:     The electron gyro-Bohm term, in chibohm's units
c    chgbi:     The ion gyro-Bohm term, in chibohm's units
c    chbe:      The electron Bohm term, in chibohm's units
c    chbi:      The ion Bohm term, in chibohm's units
c    dh:        The hydrogenic ion particle diffusivity, in [m^2/sec]
c    dimp:      The impurity ion particle diffusivity, in [m^2/sec]
c
      IMPLICIT NONE
c
c Declare variables
c NAMING CONVENTION: Dimensionless variables begin with a 'z'
c
      REAL
     &   shear,       wexb,       rmaj,       cwexb
     & , chibohm,     zq,         zrlpe,      tekev,   tikev
     & , gradte,      zra,        zdte,       btor
     & , chii,        chie,       dh,         dimp
     & , chbe,        chbi,       chgbe,      chgbi
     & , zcoef1,      zcoef2,     zcoef3,     zcoef4
     & , func,        gamma
c
      INTEGER ierr
c
c check input for validity
c
      ierr = 0
      if ((zq .lt. 0.0) .or. (zq .gt. 100.0)) then
         ierr=1
         return
      elseif (abs(zrlpe) .gt. 1000.0) then
         ierr=2
         return
      elseif ((tekev .lt. 0.01) .or. (tekev .gt. 100.0)) then
         ierr=3
         return
      elseif ((zdte .lt. 0.0) .or. (zdte .gt. 1000.0))  then
         ierr=5
         return
      elseif ((btor .lt. 0.0) .or. (btor .gt. 15.0))  then
         ierr=6
         return
      endif
!| 
!| \section{The Mixed Bohm/gyro-Bohm model}
!| 
!| The Mixed Bohm/gyro-Bohm transport model derives from an
!| originally purely Bohm-like model for electron transport
!| developed for the JET Tokamak\cite{tar94}. This preliminary model has
!| subsequently been extended to describe ion transport\cite{erb95},
!| and a gyro-Bohm term has been added in order to simulate data from
!| different machines\cite{erb98}.
!| 
!| \subsection{Bohm term}
!| 
!| The mixed model is derived using the dimensional analysis approach,
!| whereby the diffusivity in a Tokamak plasma can be written as:\\
!| \[ \chi = \chi_0 F(x_1, x_2, x_3, ...)\]
!| 
!| where $\chi_0$ is some basic transport coefficient and F is a function
!| of the plasma dimensionless parameters. We choose for $\chi_0$ the Bohm diffusivity:\\
!| \[ \chi_0 = \frac{cT_e}{eB}\]
!| 
!| The expression of the dimensionless function F is chosen according to
!| the following criteria:\
!| \begin{itemize}
!| \item{The diffusivity must be bowl-shaped, increasing towards the plasma
!| boundary}
!| \item{The functional dependencies of F must be in agreement with
!| scaling relationships of the global confinement time, reflecting
!| trends such as power degradation and linear dependence on plasma
!| current}
!| \item{The diffusivity must provide the right degree of resilience
!| of the temperature profile}
!| \end{itemize}
!| 
!| It easily shown that a very simple expression of F that satisfies
!| the above requirements is:\\
!| \[ F = aq^2/|L_{pe}^*|\]
!| 
!| where q is the safety factor and $L_{pe}^*=(dp_e/dr)^{-1}/a$, being a the
!| plasma minor radius. The resulting expression of the diffusivity
!| can be written as:\\
!| \[ \chi \propto |v_d| \Delta G\]
!| 
!| where $v_d$ is the plasma diamagnetic velocity, $\Delta=a$ and $G=q^2$,
!| so that it is clear that this model represents transport due to
!| long-wavelength turbulence.\\
!| The evidence coming up from the simulation of non-stationary
!| JET experiments \cite{erb97}(such as ELMs, cold pulses, sawteeth, etc.)
!| suggested that the above Bohm term should depend non-locally
!| on the plasma edge conditions through the temperature
!| gradient averaged over a region near the edge:\\
!| 
!| \[ <L_{T_e}^*>_{\Delta V}^{-1} = \frac{T_e(x=0.8) - T_e(x=1)}{T_e(x=1)}\]
!| 
!| where x is the normalized toroidal flux coordinate. The final
!| expression of the Bohm-like model is:\\
!| 
!| \[ \chi_{e,i}^B = \alpha_{Be,i} \frac{cT_e}{eB} L_{pe}^* q^2\]
!| 
!| where $\alpha_B$ is a parameter to be determined empirically,
!| both for ions and electrons.\\
!| 
!| \subsection{gyro-Bohm term}
!| 
!| The Bohm-like expression so derived proved to be very successful
!| in simulating JET discharges, but failed badly in smaller Tokamaks
!| such as START\cite{roa96}.\\
!| For this reason a simple gyro-Bohm-like term, also based on
!| dimensional analysis, was added:
!| 
!| \[ \chi_{e,i}^{gB} = \alpha_{gBe,i} \frac{cT_e}{eB} L_{Te}^* \rho^*\]
!| 
!| where $\rho^*$ is the normalized larmor radius:
!| 
!| \[ \rho^* =  \frac {M^{1/2}cT_e^{1/2}}{Z_ieB_t}\]
!| 
!| This expression is what can be expected from small scale
!| drift-wave turbulence. It is important to note that in large
!| Tokamaks such as JET and TFTR the gyro-Bohm term is negligible,
!| while in smaller machines, with larger values of $\rho^*$, the
!| gyro-Bohm term can play a role especially near the plasma centre.\\
!| 
!| \subsection{Final Model}
!| The resulting expressions of the diffusivities are:
!| 
!| \[ \chi_{e,i}=\chi_{Be,i}+\chi_{gBe,i}\]
!| 
!| where the Bohm and gyro-Bohm terms are defined above and the adopted values
!| of the empirical parameters are:\\
!| 
!| \[ \alpha_{Be} = 8\times10^{-5} ,\alpha_{Bi} = 2\times\alpha_{Be}\]
!| \[ \alpha_{gBe} = 3.5\times10^{-2} , \alpha_{gBi} = \alpha_{gBe}/2\]
!| 
!| The relevant coding is as follows:
!| 
c *
c * Definition of the mixed model
c *
c
c Declare the correction coeficient
c
c
c Calculate function for EXB and magnetic shear stabilization
c
         gamma = 3.0959e5 * sqrt(tikev) / (zq*rmaj)
c	 func  = (1.0 - abs(500 * wexb / gamma))
	 func  = (0.1 + shear - cwexb * abs(wexb / gamma))
c
c Calculate Bohm and gyro-Bohm terms
c

         chbe   = 8.0e-5 * zq * zq * zrlpe * zdte * chibohm
         chgbe  = 0.15811 * sqrt(tekev)*gradte/(btor* btor)
         chbi   = 2.0 * chbe
         chgbi  = chgbe / 2.0

	 if (func.lt.0) then
           chbe = 0.0
           chbi = 0.0
	 endif
c
c Now determine the actual electron and ion thermal and particle diffusivities.
c
         chie   = chbe + chgbe
         chii   = chbi + chgbi
         dh     = (zcoef1 + (zcoef2-zcoef1) * zra)
     &             * chie * chii / (chie + chii)
c
c The impurity diffusivity is not included in the mixed
c model described in Ref. [1] but is defined here using a simple
c empirical model.
c
         dimp   = zcoef3 + zcoef4 * zra * zra
c
      return
      end
!| %**********************************************************************c
!| 
!| \begin{thebibliography}{99}
!| \bibitem{tar94}
!| A. Taroni, M. Erba, E. Springmann and Tibone F.,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 36} (1994) 1629.
!| \bibitem{erb95}
!| M. Erba, V. Parail, E. Springmann and A. Taroni,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 37} (1995) 1249.
!| \bibitem{erb98}
!| M. Erba, et al.,
!| {\em Nuclear Fusion,} {\bf 38} (1998) 1013.
!| \bibitem{erb97}
!| M. Erba, et al.,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 39} (1997) 261.
!| \bibitem{roa96}
!| C.M. Roach,
!| {\em Plasma Physics and Controlled Fusion,} {\bf 38} (1996) 2187.
!| \end{thebibliography}
!| %**********************************************************************c
!| \end{document}              0.000000E+00nd of document.

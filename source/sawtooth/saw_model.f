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
!| \indent $\delta\hat{W}=\delta\hat{W}_{\rm Busac}+\delta\hat{W}_{\rm el}+\delta\hat{W}_{\rm KO}
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
!| $\omega_{\ast \rm i}(\rm sec^{-1}) = 1000T_{\rm p0}(keV)/B(\rm T)\bar{r}_{1}(m)r_{\rm p}(m)$\\
!| where $r_{\rm p}=|dp/dr|^{-1}p_{\rm i}$ is the pressure scale length.
!| 
!| \item Normalized ion larmor radius $\hat{\rho}=\rho_{\rm i}/\bar{r}_{1}$\\
!| $\hat{\rho}= 0.00323*\sqrt{T_{\rm p0}(keV)}/B(T)\bar{r}_{1}(m)$.
!| 
!| \item Semi-collisional m = 1 growth rate: $\gamma_{\rho}=1.1\hat{\rho}^{4/7}S^{-1/7}s_{1}^{6/7}\tau_{\rm A}^{-1}$\\
!| where the Lundquist number ($\tau_{\rm R}/\tau_{\rm A}$) is:\\
!| $S=43.7T_{\rm p0}^{3/2}(keV)\bar{r}_{1}^{2}/\tau_{\rm A}$.
!| 
!| \end{enumerate}
!| 
!| 
!| 
!
! cnn 5-aug-02 --->
! This is a module for testing whether a sawtooth crashed had occurred.
! If a sawtooth occur module will return is_sawt=1, if not => is_sawt=0
!
       subroutine saw_model(kt,issawt,ctrho,ctstar,cwrat,
     &               chfast,cpf,cimfast,elongr1,mui,lir1,betapf1,
     &               betaim,shearp,rmaj,rmin,rq1,rplen1,denel0,
     &               btor1,bpol11,ti0,efast,p1avg,presr1,
     &               ifastswt,iflarmswt,ikinswt,zntime)
c
!****************************
! porcelli sawtooth model reference -->
! Ref. Porcelli model: F.Porcelli, Plasm. Phys. Contr. Fus. 38, 1996 (p2163)
!
! Variables sent to module
! 
!          kt     = sawtooth model selection (k=5 or k= 6 Porcelli model)
!          issawt = output of module indicating whether a sawtooth had occurred
!                   (issawt=1 -->sawtooth crash, issawt=0 --> no sawtooth crash)
!         ctrho   = numerical multiplier for rho_star (equ 16 of reference)
!        ctstar   = numerical multiplier for c_star  (equ 16 of reference)
!         cwrat   = numerical multiplier depend on ratio between mode frequency 
!                   and fast particle precesion freq. (c_f in equ B.8 of ref)
!        chfast   = numerical multiplier for c_h (equ 16 of reference)
!           cpf   = integral of x^(3/2)*p(x) p(x)=ion pressure (c_p of equ B.7 in reference) 
!        cimfast  = charge of fast ion (not used)
!       elongr1   = elongation at the q = 1 surface
!           mui   = average ion mass (amu)
!          li_1   = internal inductance at q = 1 surface
!       betapf1   = fast ion beta parameter (equ 27 of reference)
!        betaim   = peak ion toroidal beta
!        shearp   = rdq/dr at q=1 surface
!          rmaj   = major radius (m)
!          rmin   = minor radius (m)
!           rq1   = radius of q = 1 surface (m)
!        rplen1   = p/|dp/dr| at q=1 surface (m)
!                   if p/|dp/dr|<<<1 take rplen1=radius of q=1 surface  
!        denel0   = electron density at magnetic axis (10^20 m-3)
!         btor1   = vacuum toroidal field at rmaj (testla)
!        bpol11   = poloidal field at q = 1 surface (tesla) (outboard midplane)
!           ti0   = central ion temperature (keV)
!         efast   = energy of fast ion (keV)
!         p1avg   = volume average pressure within q = 1 surface (pascal)
!        presr1   = pressure at q = 1 surface (pascal)
!        zntime   = time of step for diagnostic (sec)
!      ifastswt   = on/off switch for testing equ. 13 (fast ion)
!     iflarmswt   = on/off switch for testing equ. 14 (rotational effects)
!       ikinswt   = on/off switch for testing equ. 15 (tearing nature)
!*************
!
       use porcelli_module
!
       implicit none
       integer,intent(in)::kt,ifastswt,iflarmswt,ikinswt
       real,intent(in)::cimfast,elongr1,mui,lir1,betapf1,betaim
       real,intent(in)::shearp,rmaj,rmin,rq1,rplen1,denel0
       real,intent(in)::btor1,bpol11,ti0,efast,p1avg,presr1
       real,intent(in)::cpf,cwrat,chfast,ctstar,ctrho,zntime
       integer,intent(out)::issawt
!
       INTEGER :: iprint, i_unit  ! for diagnostic output
c
       real :: grad_norm_pres_at_q1
c
       ! writing to file for diagnostic
       call inp_rec(kt,ctrho,ctstar,cwrat,chfast,cpf,cimfast,
     &              elongr1,mui,lir1,betapf1,betaim,
     &             shearp,rmaj, rmin,rq1,rplen1,
     &             denel0,btor1,bpol11,ti0,efast,
     &             p1avg,presr1,zntime,ifastswt,iflarmswt,ikinswt)
c
       grad_norm_pres_at_q1 = 1.0 / rplen1
c
       ! call subroutine for calculations of the porcelli model
!
! temporary diagnostic printout
!
      write (6,*)
      write (6,*) 'cpf=',cpf,' elongation_at_q1=',elongr1
     & ,'ion_mass=',mui
      write (6,*) 'inductance_at_q1=',lir1
     & ,' betapolf1=',betapf1,' beta_i_max=',betaim
      write (6,*) 'shear_at_q1=',shearp,' r_major=',rmaj
     & ,'  r_minor=',rmin
      write (6,*) 'rminor_q_1=',rq1
     & ,' grad_norm_press_at_q1=',grad_norm_pres_at_q1
     & ,' denel0=',denel0
      write (6,*) 'b_toroidal=',btor1
     & ,' b_pol_at_q1=',bpol11,' ti0=',ti0
      write (6,*) 'efast=',efast
     & ,' pressure_avg_up_to_q1=',p1avg
     & ,' pressure_at_q1=',presr1
      write (6,*) 'ifastswt=',ifastswt,' iflarmswt=',iflarmswt
     & ,'ikinswt=',ikinswt
      write (6,*)
!
!..settings for diagnostic output
!
      iprint = 9
      i_unit = 6
!
       if(kt .eq. 5 .or. kt .eq. 6 )
     &    call porcelli_conditions ( issawt, iprint, i_unit, cpf,
     &                  elongr1,mui,lir1,betapf1,
     &                  betaim,shearp,rmaj, rmin,rq1,
     &                  grad_norm_pres_at_q1,
     &                  denel0,btor1,bpol11,ti0,
     &                  efast,p1avg,presr1)
c
       end subroutine saw_model
c
       ! subroutine writing module inputs into file (for diagnostic)
       subroutine inp_rec(kt,ctrho,ctstar,cwrat,chfast,cpf,cimfast,
     &              elongr1,mui,lir1,betapf1,betaim,
     &             shearp,rmaj, rmin,rq1,rplen1,
     &             denel0,btor1,bpol11,ti0,efast,
     &             p1avg,presr1,zntime,ifastswt,iflarmswt,ikinswt)
c
       implicit none
       integer,intent(in)::kt,ifastswt,iflarmswt,ikinswt
       real,intent(in)::cimfast,elongr1,mui,lir1,betapf1,betaim
       real,intent(in)::shearp,rmaj,rmin,rq1,rplen1,denel0
       real,intent(in)::btor1,bpol11,ti0,efast,p1avg,presr1
       real,intent(in)::cpf,cwrat,chfast,ctstar,ctrho,zntime
c
       write(50,100)kt,ifastswt,iflarmswt,ikinswt,zntime,ctrho,
     &          chfast,ctstar,cwrat,cpf,cimfast,
     &          elongr1,mui,lir1,betapf1,betaim,shearp,
     &          rmaj,rmin,rq1,rplen1,denel0,btor1,bpol11,
     &          ti0,efast,p1avg,presr1
c
 100   format("**************** INPUT ",/,
     &  "model = ",i3,2x,"equ. (13) = ",i3,2x,"equ. (14) = ",i3,2x,
     &  "equ. (15) = ",i3,/,"t_sec = ",f8.4, /,
     &  "ctrho = ",1pe10.3,2x,"chfast = ",1pe10.3,2x,
     &  "ctstar = ",1pe10.3, /,
     &  "cwrat = ",1pe10.3,2x,"cpf = ",1pe10.3, /,
     &  "cimfast = ",1pe10.3,2x, 
     &  "elongr1 = ",1pe10.3,2x,"mui = ",1pe10.3,2x,
     &  "lir1 = ",1pe10.3,/,"betapf1 = ",1pe10.3,2x,
     &  "betaim = ",1pe10.3,2x,"shearp = ",1pe10.3,/,
     &  "rmaj = ",1pe10.3,2x,"rmin = ",1pe10.3,2x,"rq1 = ",1pe10.3,/,
     &  "rplen1 = ",1pe10.3,2x,"denel0 = ",1pe10.3,2x,
     &  "btor1 = ",1pe10.3,/,"bpol11 = ",1pe10.3,2x,
     &  "ti0 = ",1pe10.3,2x,"efast = ",1pe10.3, /,
     &  "p1avg = ",1pe10.3,2x,"presr1 = ",1pe10.3,2x,/,
     &  "**************************************************",/,)
c
        return
c
        end subroutine inp_rec

!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{porcelli96a} F.~Pocelli, D.~Boucher, and M.~N.~Rosenbluth,
!| Plasma Physics and Controlled Fusion, {\bf 38}, 2163 (1996).
!| 
!| 
!| \end{thebibliography}
!| \end{document}

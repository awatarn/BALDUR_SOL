!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
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
!| {\bf {\tt etaw17diff.tex} \\
!| Toroidal Ion Temperature Gradient Mode  \\
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
!| {\bf This revision by P{\"a}r Strand, \today}
!| \date{}
!| \end{center}
!| This subroutine evaluates the transport matrix for $\eta_i$ and trapped
!| electron modes derived by Jan Weiland, H. Nordman and their group in
!| G\"{o}teborg Sweden.
!| The equations in this routine include fast Hydrogenic ions, impurities,
!| trapped electron and finite Larmor radius effects. New options include
!| parallel ion motion, finite beta and collisional effects together with an
!| approximate treatment of $E\times B$ flow shear reduction of transport.
!| 
!| The dimensionless diffusivity matrix $ D_{j_1,j_2} = {\tt difthi(j1,j2)}$
!| and convective velocity array $ v_{j_1} = {\tt velthi(j1)} $
!| are given as:
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
!| \vspace{0.5cm}
!| {\parindent 0cm
!| 
!|  {\bf ETAW17DIFF} calls the following routine(s)}
!| \begin{itemize}
!|    \item[{\large$\bullet$}] {\bf ETAW17FLUX} -- {\em Core routine for calculating fluxes }
!|    \begin{itemize}
!|        \item[{\normalsize $\bullet$}]
!|              {\bf TOMSQZ} -- {\em wrapper for complex QZ-algorithm routines}
!|         \begin{itemize}
!|            \item[-]{\bf CQZHES} -- {\em First step of QZ algorithm}
!|            \item[-]{\bf CQZVAL} -- {\em Second and Third step of QZ-algorithm}
!|            \item[-] {\bf CQZVEC} -- {\em Last step of QZ-algprithm}
!|        \end{itemize}
!|        \item [{\normalsize$\bullet$}]{\bf ABORTB}
!|    \end{itemize}
!| \end{itemize}
!| \newpage
c@etaw17diff.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine etaw17diff (letain, cetain, lprintin, neq, nout
     & , gnein, gnhin, gnzin, gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, difthi,velthi,chieff
     & , nmodes, perform, nerr )
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Names of variables in the argument list:}
!| \begin{tabular}{lllp{3.0in}}
!| variable & status & symbol & meaning \\
!| {\tt letain(j)} & input & & Integer control variables
!|                             (see table below).\\
!| {\tt cetain(j)} & input & & Real-valued control variables
!|                             (see table below).\\
!| {\tt lprint}    & input & & Controls printout.
!|  Higher values produce more printout. \\
!| {\tt neq} & input & & number of equations \\
!| {\tt nout} & input & & output unit for error messages \\
!| {\tt gnein} & input & $  $ & $ - R\hat{r} \cdot \nabla n_e / n_e $ \\
!| {\tt gnhin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_H / n_H $ \\
!| {\tt gnzin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_Z / n_Z $ \\
!| {\tt gtein} & input & $  $ & $ - R\hat{r} \cdot \nabla T_e / T_e $ \\
!| {\tt gthin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_H / T_H $ \\
!| {\tt gtzin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_Z / T_Z $ \\
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
!| {\tt perform(j)} &   not used &  &\\
!| {\tt nerr}       & output & & nerr $\neq 0 \rightarrow$ error\\
!| \end{tabular}
!| \end{center}
!| \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Integer control variables in the argument list:}
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| {\tt letain(2)} & = number of elements computed for transport matrix
!|                     (only when $> 0$) \\
!| {\tt letain(7)} &$ > 0 \rightarrow$ rescale transport matrix with velthi(j) = 0 \\
!| {\tt letain(9)} &$ > 0 \rightarrow$ do not produce transport matrix, only
!|                        effective diffusivities needed\\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Real-valued control variables in the argument list:}
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| {\tt cetain(10)} & = coefficient of $ k_\parallel $ (default to 1.0) \\
!| {\tt cetain(11)} & = normalization for $ A_\parallel $ (default to 1.0) \\
!| {\tt cetain(12)} & = coefficient of $H$ (default to 0.0) \\
!| {\tt cetain(15)} & = coefficient of $ \hat{\nu} $ (default to 1.0) \\
!| {\tt cetain(20)} & = coefficient of $ \beta_{e,h,z} $ (default to 1.0) \\
!| {\tt cetain(25)} & = coefficient for elongation modified Alfven frequency
!| (default to 0.0) \\
!| {\tt cetain(29)} & = radius used in printouts \\
!| {\tt cetain(30)} & = finite difference used to construct
!|                    transport matrix \\
!| \end{tabular}
!| \end{center}
!| 
!| 
c-----------------------------------------------------------------------
c
c  This version of etaw17diff is intended for use on workstations
c  Compile this routine and routines that call it with a compiler option
c  such as -r8 to convert real to double precision when used on
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: ETAW17DIFF calls the following routines
c
c    ETAW17FLUX        - Calculates fluxes and effective diffusivities
c        TOMSQZ        - Wrapper for QZ algorithm solving Ax = lambda Bx
c            CQZHES    - First step in QZ algorithm
c            CQZVAL    - Second and third step in QZ algorithm
c            CQZVEC    - Fourth step in QZ algorithm
c        ABORTB        - Fatal Error handling
c
c-----------------------------------------------------------------------

      implicit none
c
      integer idp
      parameter (idp = 15)

      logical zdiffinit
      data zdiffinit /.true./
c
      integer letain, lprintin, neq, ndim, nmodes, nerr, nout
     & ,imatrx, j1, j2, j
c
      dimension letain(*), cetain(*)
     &  , omega(*), gamma(*), chieff(*), perform(*)
     &  , difthi(ndim,*), velthi(*)
c
      real cetain
     & , gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa,  ekyrhoin, ekparlin, wexb
     & , omega, gamma, difthi, velthi, chieff, perform
c
c ndim  = first dimension of the 2-D array difthi
c         and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c imatrx= the number of elements computed along each row and
c          column of the transport matrix
c
      real  zgne, zgnh, zgnz, zgns, zgte, zgth, zgtz
      real zone, zdg, ztemp, zepsqrt,zepsmach, zdgflux(idp)
      real ztwo, zhalf, chidum(idp),           flux(idp)

      save zone, ztwo, zhalf, zepsqrt, zepsmach, zdiffinit

      if (zdiffinit) then
        zone  = 1.0
        ztwo  = 2.0
        zhalf = zone/ztwo
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
        zepsqrt = sqrt(zepsmach)
        zdiffinit = .false.
      end if

c
c..initialize arrays
c
      do j1=1,ndim
        velthi(j1) = 0.0
        do j2=1,ndim
          difthi(j1,j2) = 0.0
        end do
      enddo
c
c...Define Transport matrix sizes
c
      imatrx = min ( neq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )

c
c...copy gradients to local variables
c

      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
      zgns = zgne-zgnh*( zone-czin*fnzin-fnsin ) - zgnz*czin*fnzin
c
c... Establish fluxes
c

      call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, gnhin, gnzin, gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chieff, flux
     & , nmodes, perform, nerr )
c
c
      if (letain(9) .ne. 0) return
c
c... Define forward differencing step
c
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 0.01
c
c... Take the derivative w.r.t gth
c

      zgth = zgth + sign(zdg,zgth)
      call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, gnhin, gnzin, gtein, zgth, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chidum, zdgflux
     & , nmodes, perform, nerr )
c
      do j = 1, imatrx +1
         difthi(j,1) = (zdgflux(j) - flux(j))/sign(zdg,zgth)
      End do
      zgth = gthin

c
c... Take derivatives w.r.t gnh
c
      zgnh = zgnh + sign(zdg,zgnh)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns

      call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, gnzin, gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chidum, zdgflux
     & , nmodes, perform, nerr )
c
      do j = 1, imatrx +1
         difthi(j,2) = (zdgflux(j) - flux(j))/sign(zdg,zgnh)
      End do
      zgnh = gnhin
c
c... Take derivatives w.r.t gte
c
      zgte = zgte + sign(zdg,zgte)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns

      call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, gnhin, gnzin, zgte, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chidum, zdgflux
     & , nmodes, perform, nerr )
c
      do j = 1, imatrx +1
         difthi(j,3) = (zdgflux(j) - flux(j))/sign(zdg,zgte)
      End do
      zgte = gtein
c
c... Take derivatives w.r.t gnz
c

      if (neq .gt. 4) then

        zgnz = zgnz + sign(zdg,zgnz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, gnhin, zgnz, gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chidum, zdgflux
     & , nmodes, perform, nerr )
c
        do j = 1, imatrx + 1
           difthi(j,4) = (zdgflux(j) - flux(j))/sign(zdg,zgnz)
        End do
        zgnz = gnzin

c
c... Take derivatives w.r.t gtz
c
        zgtz = zgtz + sign(zdg,zgtz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call etaw17flux (letain, cetain, lprintin, neq, nout
     & , zgne, gnhin, gnzin, gtein, gthin, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, betahin, betazin, ftrapein
     & , vef, q, shear, kappa, ekyrhoin, ekparlin, wexb
     & , ndim, omega, gamma, chidum, zdgflux
     & , nmodes, perform, nerr )
c
        do j = 1, imatrx + 1
          difthi(j,5) = (zdgflux(j) - flux(j))/sign(zdg,zgtz)
        End do
        zgtz = gtzin
      else
        do j = 1,imatrx+1
          difthi(j,4) = 0.0
          difthi(j,5) = 0.0
        end do
      end if


c
c... Convective velocities
c
      velthi(1) = flux(1) - difthi(1,1) * gthin - difthi(1,2) * gnhin
     &                    - difthi(1,3) * gtein - difthi(1,4) * gnzin
cpis     &                    - difthi(1,5) * gtzin
c
      velthi(2) = flux(2) - difthi(2,1) * gthin - difthi(2,2) * gnhin
     &                    - difthi(2,3) * gtein - difthi(2,4) * gnzin
cpis     &                    - difthi(2,5) * gtzin

c
      velthi(3) = flux(3) - difthi(3,1) * gthin - difthi(3,2) * gnhin
     &                    - difthi(3,3) * gtein - difthi(3,4) * gnzin
cpis     &                    - difthi(3,5) * gtzin
c
      if (neq .gt. 4) then

      velthi(4) = flux(4) - difthi(4,1) * gthin - difthi(4,2) * gnhin
     &                    - difthi(4,3) * gtein - difthi(4,4) * gnzin
cpis     &                    - difthi(4,5) * gtzin


      velthi(5) = flux(5) - difthi(5,1) * gthin - difthi(5,2) * gnhin
     &                    - difthi(5,3) * gtein - difthi(5,4) * gnzin
cpis     &                    - difthi(5,5) * gtzin
      else
         velthi(4) = 0.0
         velthi(5) = 0.0
      end if
c
c... Temporary lines Needed for correspondance with etaw17a and earlier
c
      DO J = 1, imatrx+1
         DIFTHI(J,5) = 0.0
         DIFTHI(5,J) = 0.0
      END DO
      VELTHI (5) = 0.0
c
c
c... Rescale diffusivity matrix and set velthi = 0 if that is requested
c    through letain(7) > 0
c
      if (letain(7) .gt. 0) then
         do j = 1, imatrx +1
            ztemp = flux(j) - velthi(j)
            do j1 = 1, imatrx +1
                difthi(j,j1) = difthi(j,j1)*flux(j) / ztemp
            end do
            velthi(j) = 0.0
         end do
      endif
c
c... Diagnostic printout
c
      if ( lprintin .gt. 0 ) then
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn etaw17diff'
c
        write (nout,*)
        write (nout,*) ' diffusion matrix from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
        write (nout,*)
        write (nout,*) ' difthi(j1,j2)  -> j2'
        do j1=1,imatrx +1
          write (nout,132) (difthi(j1,j2),j2=1,imatrx+1)
        enddo
c
        write (nout,*)
        write (nout,*) ' convective velocities'
     &    ,' normalized by omega_{De} / k_y'
        do j1=1,imatrx+1
          write (nout,132) velthi(j1)
        enddo
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

      end
!| \end{document}
!| 
!| 
!| 
!| 

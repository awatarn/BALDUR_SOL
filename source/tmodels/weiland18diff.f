!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| \documentclass[12pt]{article}
!| \usepackage{longtable}
!| % \documentstyle{article}
!| 
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
!| {\bf {\tt weiland18diff.tex} \\
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
!| {\bf This revision by Glenn Bateman, \today}
!| \date{}
!| \end{center}
!| 
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
!| \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!| \end{array} \right) $$
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $,
!| convective velocities are normalized by $ \omega_{De} / k_y $,
!| and all the frequencies are normalized by $ \omega_{De} $.
!| \vspace{0.5cm}
!| {\parindent 0cm
!| 
!| {\bf WEILAND18DIFF} calls the following routines}
!| \begin{description} \item[WEILAND18FLUX]
!| % \item[{\large$\bullet$}] {\bf WEILAND18FLUX}
!|     -- {\em Core routine for calculating fluxes }
!| \begin{description}
!| \item[WEILAND18DISP] --{\em Solve the dispersion relation for
!|   eigenvalues and eigenfunctions}
!| \begin{description}
!| \item[WEILAND18EQNS] --{\em Computes the matrix of
!|   linearized equations}
!| \begin{description}
!| \item[TOMSQZ] -- {\em wrapper for complex QZ-algorithm routines}
!| \begin{description}
!| \item[CQZHES] -- {\em First step of QZ algorithm}
!| \item[CQZVAL] -- {\em Second and Third step of QZ-algorithm}
!| \item[CQZVEC] -- {\em Last step of QZ-algprithm}
!| \end{description}
!| \item[ABORTB]
!| \end{description}
!| \end{description}
!| \end{description}
!| \end{description}
!| 
!| % \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Names of variables in the argument list:} \\\vspace{1pc}
!| \begin{longtable}{lllp{3.0in}}
!| variable & status & symbol & meaning \\
!| {\tt letain(j)} & input & & Integer control variables
!|                           (see table below).\\
!| {\tt cetain(j)} & input & & Real-valued control variables
!|                           (see table below).\\
!| {\tt lprint}    & input & & Controls printout.
!| Higher values produce more printout. \\
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
!| & impurity mass to hydrogen isotope mass. \\
!| & & & Note that ``hydrogen'' may include a deuterium or tritium mix. \\
!| {\tt fnsin}  & input & $ f_s $
!| & $ f_s = n_s / n_e $ fraction of superthermal hydrogenic ions \\
!| {\tt betaein} & input & $\beta_e$ &
!| $ = n_e T_e / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt betahin} & input & $\beta_H$ &
!| $ = n_H T_H / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt betazin} & input & $\beta_Z$ &
!| $ = n_Z T_Z / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt ftrapein}  & input & $f_{trap} $ &
!| fraction of trapped electrons \\
!| {\tt vef}       & input & $ \nu_{th} / \omega_{De} $ &
!| thermal collision frequency, normalized \\
!| {\tt q}         & input & $ q $ & magnetic q-value \\
!| {\tt shear}     & input & $ s $ & $ d \ln q / d \ln r $ \\
!| {\tt ekyrhoin}  & input & $ k_y \rho_s $ & normalized poloidal
!|                       wave number \\
!| {\tt ekparlin}  & input & $k_\parallel L_n$ & \\
!| {\tt wexb}      & input & $\omega_{E\times B}$& $E\times B$ shearing rate
!| (normalized with $\omega_{D_e}$) \\
!| {\tt ndim} & input & & first dimension of the 2-D array difthi
!|               and the maximum number of unstable modes allowed \\
!| {\tt omega(j)}  & output & $\omega / \omega_{De} $ &
!| real part of the frequencies normalized by $ \omega_{De} $ \\
!| {\tt gamma(j)}  & output & $\gamma / \omega_{De} $ &
!| growth rates normalized by $ \omega_{De} $ \\
!| {\tt difthi(i,j)}      & output & $ D \omega_{De} / k_y^2 $
!| & diffusivity matrix normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt velthi(i,j)}      & output & $ v \omega_{De} / k_y $
!| & convective velocities normalized by $ k_y / \omega_{De} $ \\
!| {\tt chieff(j)} & output & $ \chi_{\rm eff} \omega_{De} / k_y^2 $
!| & effective total diffusivities
!| for $ n_H T_H $, $ n_H $, $ n_e T_e $,
!| $ n_Z $, $ n_Z T_Z $, \ldots
!| normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt nmodes} & output & & number of unstable modes \\
!| {\tt nerr}       & output & & nerr $\neq 0 \rightarrow$ error\\
!| \hline
!| \end{longtable}
!| \end{center}
!| 
!| % \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Integer control variables in the argument list:}
!| \begin{longtable}{lp{4.0in}}
!| variable & meaning \\
!| {\tt letain(2)} & = number of elements computed for transport matrix
!| (only when $> 0$) \\
!| {\tt letain(7)} &$ > 0
!|     \rightarrow$ rescale transport matrix with velthi(j) = 0 \\
!| {\tt letain(9)} &$ > 0 \rightarrow$ do not produce transport matrix, only
!| effective diffusivities needed\\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| \hline
!| \end{longtable}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Real-valued control variables in the argument list:}
!| \begin{longtable}{lp{4.0in}}
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
!| transport matrix \\
!| \hline
!| \end{longtable}
!| \end{center}
!| 
c@weiland18diff.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine weiland18diff (letain, cetain, lprintin, neq, nout
     & , gnein, gnhin, gnzin, gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , gradrhosqrave, gradrhoave, gradrhooutbrd
     & , hcn,       rot,       tol,       itl
     & , ndim, omega, gamma, difthi, velthi, chieff, flux
     & , nmodes, nerr, its, itera, em )
c
c-----------------------------------------------------------------------
c
c  This version of weiland18diff is intended for use on workstations
c  Compile this routine and routines that call it with a compiler option
c  such as -r8 to convert real to double precision when used on
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: WEILAND18DIFF calls the following routines
c
c  WEILAND18FLUX       - Calculates fluxes and effective diffusivities
c    WEILAND18DISP     - Computes eigenvalues and eigenvectors
c      WEILAND18EQNS   - Computes the matrix of linearized equations
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
      parameter (idp = 12)

      logical zdiffinit
      data zdiffinit /.true./
c
      integer letain(*)
c
      integer lprintin, neq, ndim, nmodes, nerr, nout
     & ,imatrx, j1, j2, j
c
      real cetain(*)
c
      real omega(*), gamma(*), chieff(*), flux(*)
     &  , difthi(ndim,*), velthi(*)
c
      real      gnein, gnhin, gnzin
     & , gtein, gthin, gtzin, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, kappa,  ekyrhoin, ekparlin, wexb
     & , em
     & , gradrhosqrave, gradrhoave, gradrhooutbrd
c
c   Metric elements:
c     Here, rho is any flux surface label that increases monotonically
c       from the magnetic axis to the edge of the plasma.
c  gradrhosqrave = <|grad rho|^2>
c  gradrhoave    = <|grad rho|>
c  gradrhooutbrd = |grad rho| at outboard edge of each flux surface
c
c ndim  = first dimension of the 2-D array difthi
c         and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c imatrx= the number of elements computed along each row and
c          column of the transport matrix
c
      integer its, itera, itl
c
      real aspinv, hydmass
c
      real  hcn, rot, tol
c
      real  zgne, zgnh, zgnz, zgns, zgte, zgth, zgtz
      real  zdg, ztemp
      real  zone, ztwo, zhalf, zepsqrt, zepsmach
      real  zomega(idp), zgamma(idp), zchidum(idp), zdgflux(idp)

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

      call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, omega, gamma, chieff, flux
     & , nmodes, nerr, its, itera, em )
c
      if ( lprintin .gt. 1 ) then
        write (6,*)
        write (6,*) ' affter calling weiland18flux'
        write (6,*) ' letain(9) = ', letain(9)
        write (6,*) '  omega   gamma'
        do j1=1,neq
          write(6,*) omega(j1), gamma(j1)
        enddo
      endif
c
      if ( letain(9) .ne. 0 ) return
c
c... Define forward differencing step
c
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 0.01
c
c... Take the derivative w.r.t gth
c
c..Note, for all the following calls to weiland18flux
c  use zomega, zgamma, zchidum, zdgflux
c  to avoid overwriting omega, gamma, chieff, flux
c
      zgth = zgth + sign(zdg,zgth)
c
      call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, zomega, zgamma, zchidum, zdgflux
     & , nmodes, nerr, its, itera, em )
c
      do j = 1, imatrx +1
         difthi(j,1) = (zdgflux(j) - flux(j))/sign(zdg,zgth)
      end do
      zgth = gthin

c
c... Take derivatives w.r.t gnh
c
      zgnh = zgnh + sign(zdg,zgnh)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns
c
      call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, zomega, zgamma, zchidum, zdgflux
     & , nmodes, nerr, its, itera, em )
c
      do j = 1, imatrx +1
         difthi(j,2) = (zdgflux(j) - flux(j))/sign(zdg,zgnh)
      end do
      zgnh = gnhin
c
c... Take derivatives w.r.t gte
c
      zgte = zgte + sign(zdg,zgte)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns

      call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, zomega, zgamma, zchidum, zdgflux
     & , nmodes, nerr, its, itera, em )
c
      do j = 1, imatrx +1
         difthi(j,3) = (zdgflux(j) - flux(j))/sign(zdg,zgte)
      end do
      zgte = gtein
c
c... Take derivatives w.r.t gnz
c

      if (neq .gt. 4) then

        zgnz = zgnz + sign(zdg,zgnz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, zomega, zgamma, zchidum, zdgflux
     & , nmodes, nerr, its, itera, em )
c
        do j = 1, imatrx + 1
           difthi(j,4) = (zdgflux(j) - flux(j))/sign(zdg,zgnz)
        end do
        zgnz = gnzin

c
c... Take derivatives w.r.t gtz
c
        zgtz = zgtz + sign(zdg,zgtz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call weiland18flux (letain, cetain, lprintin, neq, nout
     & , zgne, zgnh, zgnz, zgte, zgth, zgtz, tauhin, tauzin
     & , fnzin, czin, azin, fnsin, betaein, ftrapein
     & , vef, q, shear, aspinv, kappa, hydmass
     & , ekyrhoin, ekparlin, wexb
     & , hcn,       rot,       tol,       itl
     & , ndim, zomega, zgamma, zchidum, zdgflux
     & , nmodes, nerr, its, itera, em )
c
        do j = 1, imatrx + 1
          difthi(j,5) = (zdgflux(j) - flux(j))/sign(zdg,zgtz)
        end do
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
     &                    - difthi(1,5) * gtzin
c
      velthi(2) = flux(2) - difthi(2,1) * gthin - difthi(2,2) * gnhin
     &                    - difthi(2,3) * gtein - difthi(2,4) * gnzin
     &                    - difthi(2,5) * gtzin

c
      velthi(3) = flux(3) - difthi(3,1) * gthin - difthi(3,2) * gnhin
     &                    - difthi(3,3) * gtein - difthi(3,4) * gnzin
     &                    - difthi(3,5) * gtzin
c
      if (neq .gt. 4) then

      velthi(4) = flux(4) - difthi(4,1) * gthin - difthi(4,2) * gnhin
     &                    - difthi(4,3) * gtein - difthi(4,4) * gnzin
     &                    - difthi(4,5) * gtzin


      velthi(5) = flux(5) - difthi(5,1) * gthin - difthi(5,2) * gnhin
     &                    - difthi(5,3) * gtein - difthi(5,4) * gnzin
     &                    - difthi(5,5) * gtzin
      else
         velthi(4) = 0.0
         velthi(5) = 0.0
      end if
c
c... Temporary lines Needed for correspondance with weiland18a and earlier
c
cbate      DO J = 1, imatrx+1
cbate         DIFTHI(J,5) = 0.0
cbate         DIFTHI(5,J) = 0.0
cbate      END DO
cbate      VELTHI (5) = 0.0
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
     &    ,'  Final results from sbrtn weiland18diff'
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
!|  \end{document}
!| 
!| 
!| 
!| 

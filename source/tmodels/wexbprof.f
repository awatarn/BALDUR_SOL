!| %
!| % 29-sep-99 .../baldur/code/bald/wexbprof.tex
!| %
!| %  This is a LaTeX ASCII file.  To typeset document type: latex theory
!| %  To extract the fortran source code, obtain the utility "xtverb" from
!|  0lenn Bateman (bateman@pppl.gov) and type:
!| 0tverb <wexbprof.tex > wexbprof.f
!| %
!| \documentstyle {article}    % Specifies the document style.
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| 
!| \title{ {\tt wexbprof}: a BALDUR Subroutine \\
!|  Computing {\bf E}\times{\bf B} Flow Shear Parameters}     % title.
!| \author{
!|         Par Strand \\ Chalmers University of Technology \\
!| 	Ping Zhu \\ UT Austin }
!|                             0eclares the author's name.
!|                             0eleting the \date{} produces today's date.
!| \begin{document}            0.000000E+00nd of preamble and beginning of text.
!| \maketitle                 % Produces the title.
!| 
!| This report documents a subroutine called {\tt wexbprof}, which computes Hahm--Burrel {\bf E}\times{\bf B} shearing rate $\omega_s$.
!| 
c--------1---------2---------3---------4---------5---------6---------7-c
c@wexprof  .../baldur/code/bald/wexbprof.tex
c  pzhu 9-sep-99 removed vpol calculation part, pass vpol from argument list
c  rgb  9-oct-98 implemented zcoefpol = 0.0, ismooth1 = 5, ismooth2 = 3
c    as temporary control variables
c  me   8-oct-98 different computation of derivative zzx2
c  me   7-oct-98 smoothing applied to wexba
c  pis 30-jul-98 corrected the calc. to take into account  ghostpoints
c                in the Baldur grid
c  pis 28-jul-98 provided output to long output
c  pis 28-jul-98 added iknthe to argument list (print control)
c  pis 28-jul-98 added nprint to the argument list
c  me  28-jul-98 splitted wexba into wexbaa, wexbab, wexbac
c-----------------------------------------------------------------------
      subroutine wexbprof(nx, rho, phi, ni, ti, bpol, btor
     &         , vpol, vtor, rminor, rmajor, xbouni, wexba
     &         , iknthe, nprint)
c
      Implicit NONE
c-----------------------------------------------------------------------
      integer :: nx
      integer :: nx1
c-----------------------------------------------------------------------

      real, dimension(nx)   :: rho     ! sqrt norm. toroidal flux    [-]
      real, dimension(nx,2) :: phi     ! poloidal flux / 2 pi       [Wb]
      real, dimension(nx)   :: ni      ! main ion density         [m^-3]
      real, dimension(nx)   :: Ti      ! Ion temperature            [eV]
      real, dimension(nx)   :: Bpol    ! Poloidal magnetic field     [T]
      real, dimension(nx)   :: Btor    ! Toroidal magnetic field     [T]
      real, dimension(nx)   :: vtor    ! Toroidal rotation       [rad/s]
      real, dimension(nx)   :: vpol    ! Poroidal rotation         [m/s]
      real, dimension(nx)   :: rminor  ! Minor radius                [m]
      real, dimension(nx)   :: rmajor  ! Major radius                [m]
      real, dimension(nx)   :: xbouni  ! sq. root norm.toroidal flux [-]
      integer               :: iknthe  ! print flag (=3 -> print)
      integer               :: nprint ! print unit
c-----------------------------------------------------------------------
c Output variable
c-----------------------------------------------------------------------

      real, dimension(nx) ::   wexba   ! ExB shearing rate       [rad/s]
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------

      real,  parameter    ::   ze = 1.60219e-19  !unit charge        [C]
      real                ::   zfc, zk1
      integer, parameter  ::   intrp  = 2 ! use Akima spline
      integer, parameter  ::   kextrp = 1 ! extrapolate linearly
      integer, parameter  ::   kderiv = 1 ! calculate first derivatives
      integer, parameter  ::   kskip  = 1 ! use only 1-d arrays as input
      integer             ::   kxlast
      integer             ::   j, i , jj

c-----------------------------------------------------------------------
c Automatic arrays
c-----------------------------------------------------------------------

      real, dimension(nx)  :: zdTi          ! ion temperature gradient
      real, dimension(nx)  :: zdpi          ! ion pressure gradient
      real, dimension(nx)  :: zvpol         ! poloidal rotation
      real, dimension(nx)  :: zErad         ! radial electric field
      real, dimension(nx)  :: zErada        ! gradpi contribution to Erad
      real, dimension(nx)  :: zEradb        ! vpol contribution to Erad
      real, dimension(nx)  :: zEradc        ! vtor contribution to Erad
      real, dimension(nx)  :: zEpsr
      real, dimension(nx)  :: fpol
      real, dimension(nx)  :: zTi           ! Local ion temperature
      real, dimension(nx)  :: zpi           ! Local pressure
      real, dimension(nx)  :: zzx

      real, dimension(nx)  :: zzxa
      real, dimension(nx)  :: zzxb
      real, dimension(nx)  :: zzxc

      real, dimension(nx)  :: zzx2a
      real, dimension(nx)  :: zzx2b
      real, dimension(nx)  :: zzx2c

      real, dimension(nx)  :: zzx2

      real, dimension(nx)  :: wexbaa
      real, dimension(nx)  :: wexbab
      real, dimension(nx)  :: wexbac

      real, dimension(nx)  :: wtemp1
      real, dimension(nx)  :: wtemp2

      real, dimension(nx)  :: zr

      real                 :: pi

      real                 :: zcoefpol

      integer              :: ismooth1, ismooth2

c-----------------------------------------------------------------------
c Start of actual coding
c-----------------------------------------------------------------------
c
c... Initial setup
c
      pi =  atan2(0.,-1.)
c
      zcoefpol = 0.0
      ismooth1 = 5
      ismooth2 = 3
c
      do jj = 1, nx
        zvpol(jj)  = vpol(jj)
        zerada(jj) = 0.0
        zeradb(jj) = 0.0
        zeradc(jj) = 0.0
        zerad(jj)  = 0.0
      enddo
c
      DO jj = 2, nx
         j = jj -1
         zepsr(j) = rminor(j)/rmajor(j)           ! inverse aspect ratio
         zpi(j)  = ni(j) * Ti(j) * ze * 1.e3      ! Temporary store Pressure
         zti(j)  = Ti(j) * 1.e3                   ! and ion temp. keV -> eV
         zr (j)  = rminor(j)
         fpol(j) = phi(j,1)/2./pi*1.e-8 ! Conversion factor correct ???
      END DO


      zr(1)=-rminor(3)
      zr(2)=0.0

!ap
      kxlast = nx-1
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,zr,zpi,nx-1,zr
     &          ,kxlast,zdpi)
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,zr,zti,nx-1,zr
     &          ,kxlast,zdTi)

c
c... Loop over the grid except axis value, axis point needs exception
c    handling due to Bpol(1) = 0
c

      DO j = 1, nx-1
         zErada(j) = zdpi(j)/ze/ni(j)
         zEradb(j) = -zvpol(j) * Btor(j)
c         zEradc(j) = vtor(j+1) * zr(j) * bpol(j+1)
         zEradc(j) = vtor(j) * rmajor(j) * bpol(j)
         zErad(j) = zErada(j) + zEradb(j) +  zEradc(j)
         zzxa(j) = zErada(j)/ rmajor (j+1) / Bpol(j+1)
         zzxb(j) = zEradb(j)/ rmajor (j+1) / Bpol(j+1)
         zzxc(j) = zEradc(j)/ rmajor (j+1) / Bpol(j+1)
         zzx(j) = zzxa(j) + zzxb(j) +  zzxc(j)
      end do

c
c... Use zzx = constant towards axis => wexb = 0 on axis
c

      zzxa(1) = zzxa(2)
      zzxb(1) = zzxb(2)
      zzxc(1) = zzxc(2)
      zzx(1) = zzx(2)


c
c... Calculate the derivative w.r.t poloidal flux/ 2 pi
c

c      DO j = 2, nx-1
c         zzx2a(j) = ( zzxa(j+1)-zzxa(j) )/( fpol(j+1)-fpol(j) )
c         zzx2b(j) = ( zzxb(j+1)-zzxb(j) )/( fpol(j+1)-fpol(j) )
c         zzx2c(j) = ( zzxc(j+1)-zzxc(j) )/( fpol(j+1)-fpol(j) )
c         zzx2(j) = ( zzx(j+1)-zzx(j) )/( fpol(j+1)-fpol(j) )
c      end do

      kxlast = 0

      if (iknthe .ge. 3) then
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,fpol,zzxa,nx-1,fpol
     &           ,kxlast,zzx2a)
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,fpol,zzxb,nx-1,fpol
     &           ,kxlast,zzx2b)
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,fpol,zzxc,nx-1,fpol
     &           ,kxlast,zzx2c)
      end if
      call int1d(intrp,kextrp,kderiv,nx-1,kskip,fpol,zzx,nx-1,fpol
     &           ,kxlast,zzx2)

c
c Calculate the shearing rate
c and its three contributions (ME)
c
      DO J = 1, nx-1
         if (iknthe .ge. 3) then
         wexbaa(j) = (Rmajor(j+1)*bpol(j+1))**2/sqrt(Btor(j+1)**2
     &                 + Bpol(j+1)**2) * zzx2a(j)
         wexbab(j) = (Rmajor(j+1)*bpol(j+1))**2/sqrt(Btor(j+1)**2
     &                 + Bpol(j+1)**2) * zzx2b(j)
         wexbac(j) = (Rmajor(j+1)*bpol(j+1))**2/sqrt(Btor(j+1)**2
     &                 + Bpol(j+1)**2) * zzx2c(j)
         end if
         wexba(j)  = (Rmajor(j+1)*bpol(j+1))**2/sqrt(Btor(j+1)**2
     &                 + Bpol(j+1)**2) * zzx2(j)
      end do
c
c Smooth the shearing rate and all its terms
c

      nx1 = nx - 1
      call smooth2 ( wexbaa, 1, wtemp1, wtemp2, 1
     &    , 1, nx1, ismooth1, 0.0 )
      call smooth2 ( wexbab, 1, wtemp1, wtemp2, 1
     &    , 1, nx1, ismooth1, 0.0 )
      call smooth2 ( wexbac, 1, wtemp1, wtemp2, 1
     &    , 1, nx1, ismooth2, 0.0 )
      call smooth2 ( wexba, 1, wtemp1, wtemp2, 1
     &    , 1, nx1, ismooth1, 0.0 )
c
c Final printout
c

      if (iknthe .ge. 3) then
      write(nprint,998)
          DO j = 1, nx-1
             write(nprint,999) zr(j), wexba(j), wexbaa(j),wexbab(j)
     &                       , wexbac(j), fpol(j), xbouni(j)
          end do
      write(nprint,1000)
          Do j = 1, nx-1
             write(nprint,1001) zr(j), zErad(j), zErada(j)
     &                       , -zEradb(j), zEradc(j)
          end do
      endif
 998  FORMAT (/,3x,'radius',3x,' wexba  ',3x,'Dia.magn',3x,' vpol*Bt'
     & ,3x,' vtor*Bp', 3x, 'sqpol.fl',3x,'     rho')
 999  FORMAT (3x,0pf6.3,1x,6(1pe10.3,1x))
 1000 FORMAT (/,3x,'radius',3x,'  Er  ',3x,'gradP',3x,' vpol*Bt'
     & ,3x,' vtor*Bp')
 1001 FORMAT (3x,0pf6.3,1x,6(1pe10.3,1x))
      END
!| \end{document}

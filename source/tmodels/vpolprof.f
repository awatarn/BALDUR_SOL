!| %
!| \documentstyle {article}    % Specifies the document style.
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| 
!| \title{ {\tt ptheory}: a BALDUR Subroutine \\
!|  Interface between the BALDUR Transport Code \\
!|  and Subroutine Ntheory}     % title.
!| \author{
!|         Ping Zhu \\ UT Austin}
!|                  0eclares the author's name.
!|                  0eleting the \date{} produces today's date.
!| \begin{document}            0.000000E+00nd of preamble and beginning of text.
!| \maketitle                 % Produces the title.
!| 
!| This report documents a subroutine called {\tt vpolprof}, which
!| interfaces with different routines computing plasma poloidal velocity
!| using various neoclassical theory-based models.
!| 
!| Current included models:
!| \begin{enumerate}
!| \item Model implemented by P. Strand. ({\tt imod = 1})
!| \item Model implemented by M. Erba. ({\tt imod = 2})
!| \item NEO model implemented by P. Zhu. ({\tt imod = 3})
!| \end{enumerate}
!| 
!| The default is NEO model, which was most recently reported in 
!| Zhu~{\it etal}~\cite{zhu99}.
!| 
c@vpolprof .../baldur/bald/vpolprof.f
c rgb 06-sep-01 added maxs to the argument list and first dimensions
c   used maxs as first dimension of den, temp, and vpol
c pzhu 06-Sep-99 included the NEO model
c   See rest of changes list at end of file
c--------1---------2---------3---------4---------5---------6---------7-c
!| 
c
      subroutine vpolprof(imod, ns, maxs, ma, za, nx, den, temp, vpol,
     &  rmin, rmaj, btor, bpol, qpsi, ierr)
c
      implicit none
c
      integer imod              !model index
                                !1: Strand formula
                                !2: Erba formula
                                !3: NEO model
      integer :: ierr       	!error flag
                                !0: no error upon return
      integer :: iprint = 0	!control parameter for output
c
      integer ns                !number of particle species 
                                ! including electrons
                                !1: electron
                                !2: main ion
                                !3: impurity ion
      integer :: maxs           ! first dimension of arrays den, temp,
                                !  vpol
      real ma(ns)               !particle mass [u]
      real za(ns)               !particle charge [e]
      integer nx                !number of radial grid points
      real den(maxs,nx)         !density [m^-3]
      real temp(maxs,nx)        !temperature [eV]
      real vpol(maxs,nx)        !poloidal velocity [m/s]
      real rmin(nx)             !minor radius of plasma [m]
      real rmaj(nx)             !major radius of plasma [m]
      real btor(nx)             !toroidal B field [T]
      real bpol(nx)             !poloidal B field [T]
      real qpsi(nx)             !safty factor
c
      integer :: intrp  = 2     !use Akima spline
      integer :: kextrp = 1     !extrapolate linearly
      integer :: kderiv = 1     !calculate first derivatives
      integer :: kskip  = 1     !use only 1-d arrays as input
      integer kxlast
      integer is                !local species index
      integer jx                !local radial grid index
      real :: zcoefpol = 1.d0   !control parameter
!ap
      real, allocatable :: zrmn(:)
c
      real, allocatable :: zepsr(:)            !inverse aspect ratio
      real, allocatable :: zdti(:)             !temperature gradient
c
      real zfc, zk1
c
      real, allocatable :: xarr(:)
      real, allocatable :: xbarr(:), ztemp(:)
c
      real zrmin, zrmaj, zbtor
c
      integer j1, j2
c
      select case (imod)
c
      case(1)
c-----------------------------------------------------
c     Neoclassical expression for vpol (by P.Strand)
c-----------------------------------------------------
c
        ierr = 0
!ap	
        allocate (zepsr(nx), zdti(nx), ztemp(nx), zrmn(nx), stat = ierr)
        zrmn = rmin
	zrmn(1) = zrmn(3)
	zrmn(2) = 0.
!	
        if ( ierr .ne. 0 ) then
          write (*,*) 'unable to allocate zepsr, zdti, ztemp',
     &    ' in vpolprof.f'
          stop
        endif
c
        do jx = 1, nx-1
          zepsr(jx) = zrmn(jx)/rmaj(jx)
        end do
c
        do is = 1, ns
c
          do j1=1,nx
            ztemp(j1) = temp(is,j1)
          enddo
!ap
          kxlast = nx-1
          call int1d(intrp,kextrp,kderiv,nx-1,kskip,zrmn,ztemp,nx-1,zrmn
     &      ,kxlast,zdti)
          do jx = 1, nx-1
            zfc = 1.0 - 1.46*sqrt(zepsr(jx)) +
     &			0.46*zepsr(jx)*sqrt(zepsr(jx))
            zk1 = 0.8839 * zfc/ (0.3477 + 0.4058 * zfc)
!ap
            vpol(is,jx) = zcoefpol * zk1 * zdti(jx) * btor(jx+1) /
     &        (btor(jx+1)**2 + bpol(jx+1)**2)
          end do
        end do
!ap
        deallocate ( zepsr, zdti, ztemp, zrmn )
c
      case (2)
c-----------------------------------------------------
c     Experimental expression for vpol (by M.Erba)
c-----------------------------------------------------
c
        do j2=1,nx
          do j1=1,ns
            vpol(j1,j2) = 1.e3*temp(j1,j2)
          enddo
        enddo
c
      case (3)
c------------------------------------------------------------------
c     Neoclassical routine for vpol based on NEO model (by P.Zhu)
c------------------------------------------------------------------
c
        ierr = 0
        allocate ( xarr(nx), xbarr(nx), stat = ierr )
        if ( ierr .ne. 0 ) then
          write (*,*) 'unable to allocate xarr and xbarr in vpolprof.f'
          stop
        endif
c
        vpol = 0.0
c
        do j1=1,nx
	  xarr(j1) = rmin(j1) / rmin(nx)
        enddo
c
        zrmin = rmin(nx)
        zrmaj = rmaj(3)
        zbtor = btor(3)
c
        call neo_pol(ns, ma(1:ns), za(1:ns),
     &	  nx-2, xarr(3:nx), xarr(3:nx),
     &	  den(1:ns,1:nx-2), temp(1:ns,1:nx-2), vpol(1:ns,1:nx-2),
     &    zrmin, zrmaj, zbtor, qpsi(3:nx), ierr)
!ap{
!        call neo_pol(ns, maxs, ma, za,
!     &	  nx-2, xarr(3:nx), xarr(3:nx),
!     &	  den, temp, vpol,
!     &    zrmin, zrmaj, zbtor, qpsi(3:nx), ierr)
!}
        if(ierr.ne.0) then
          write(*,*)'vpolprof: neo_pol failed!'
          stop
        end if
	vpol(1:ns,2) = 0.e0
	vpol(1:ns,1) = -vpol(1:ns,3)
      case default
        write(*,*)'vpolprof: imod out of range!'
        ierr = 1
        stop
      end select
c
cbate      if(iprint.eq.1) then
cbate      	open(1, file='vp.d')
cbate      	do jx = 2, Nx-1
cbate      	write(1,*)xarr(jx), (vpol(is,jx),is=1,Ns)
cbate      	end do
cbate      	close(1)
cbate      end if
c
       if (allocated(xarr)) deallocate ( xarr )
       if (allocated(xbarr)) deallocate ( xbarr )
c
      return
      end
!| 
!| \end{document}

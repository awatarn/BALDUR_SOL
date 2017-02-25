!| %
!| % 15:00 27-Dec-2001 .../baldur/code/bald/ptheory.f  Bateman, Lehigh
!| %
!| %  To turn this Fortran source file into a LaTeX file, type:
!| %  python s2tex.py ptheory.f
!| %  The s2tex.py script can be obtained from
!| %  bateman@fusion.physics.lehigh.edu
!| %
!| \documentstyle {article}	% Specifies the document style.
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| 
!| \title{ {\tt ptheory}: a BALDUR Subroutine \\
!|  Interface between the BALDUR Transport Code \\
!|  and Subroutine Theory}	 % title.
!| \author{
!|         Glenn Bateman \\ Princeton Plasma Physics Laboratory}
!|                            % Declares the author's name.
!|                            % Deleting the \date{} produces today's date.
!| \begin{document}           % End of preamble and beginning of text.
!| \maketitle                 % Produces the title.
!| 
!| This report documents a subroutine called {\tt theory}, which computes plasma
!| transport coefficients using various microinstability theory-based models.
!| The default model used in subroutine {\tt theory} was developed by
!| C.~E. Singer as documented in
!| ``Theoretical Particle and Energy Flux Formulas for Tokamaks,''
!| Comments on Plasma Physics and Controlled Fusion, {\bf 11}, 165 (1988)
!| (\cite{Comments}, hereafter referred to as the Comments paper).
!| 
c@ptheory  .../baldur/code/bald/ptheory.f
c rgb 8-oct-03  changed argument list for mixed_model to use module
c rgb 2-may-02  changed zdqdrho(jz) for callglf2db
c rgb 30-apr-02 added diagnostic output for chi from TEM (zchitem, ...)
c rgb 27-dec-01 implementeded the GLF23 transport model
c rgb 19-sep-01 corrected serious error computing zimp and zmass
c   when there are two or more impurity species
c rgb 18-sep-01 set zthzmix = zthdmix after call to mixed_model
c rgb 17-sep-01 changed arrays in call to vpolprof and wexbprof
c rgb 10 Jul-01 removed if ( lthery(21) .eq. 10 ) ... call mmm95a
c rgb 17-jun-01 added veltis(jz) and veltes(jz) to neoclassical printout
c rap 29-may-01 - normalization of xineo1 is changed
c rgb 09-may-01 added growth rates and frequencies to output
c tho 23-jun-00 - add switch lthery(34) for mmm99 to use weiland18
c               - add switch lthery(40) for choosing E&M effect in weiland18
c rap 22-mar-00 - printout of diffusion matrix suppresed if lthery(39) >= 0,
c                 and zxithe/zdhthe... if lthery(39) < 0FW: µØê¡µÒºÃÒÂªÔ´«éÒÂä»àÅ
c               - coefficient cfutz(11) and cfutz(12) accounted were accounted
c                 for printout of xeneo1 and xineo1
c rap 06-mar-00 switch lthery(39) for choosing difthi/velthi or zxithe/
c               zvithe ... added
c rgb 31-dec-99 set diffusivities at magnetic axis equal to
c   diffusivities at iaxis+1 when lthery(21) > 0
c pzhu 09-dec-99 move feg to common block cbaldr; pass lreset = 1 to mmm99
c rgb 24-nov-99 take gradients wrt R_{\rm outboard} if lthery(23) = 1
c   R_{\rm outboard} 9= R_{\rm Geometric Center} + r  for major radius
c rgb 23-nov-99 write 'Sbrtn ... called from sbrtn ptheory' once
c pzhu 29-sep-99 add call to mmm99
c pzhu 9-sep-99 add call to <vpolprof> for poloidal rotation
c   routine interfacing
c tho 14-may-99 Implemented Onjun's version  mixed Bohm/gyro-Bohm model
c sep-98 Matteo Erba added the ohe and mixed Bohm/gyro-Bohm models
c rgb 13-oct-98 added sbrtn mmm98b using sda04dif
c   rearranged output
c rgb 09-oct-98 added sbrtn mmm98
c rgb 30-sep-98 implemented limit on time rate of change of ITG
c rgb 29-sep-98 if lthery(24) = 1, set zdesnsi = zdensh + zdensimp
c   and zgrdni = zgrdnh + zgrdnz
c   and zgrdpr computed from zgrdti + zgrdni + zgrdte + zgrdne
c rgb 28-sep-98 zero out diffusivities at magnetic axis
c rgb 19-sep-98 added mmm95 and removed arrays that were not used
c pis 02-jul-98 added zvrotxb for toroidal rotation
c pis 02-jul-98 added interpolation for zvrotxb using wexbint
c pis 12-jun-98 added zwexbxb for flow shear
c rgb 15-dec-96 added diagnostic output sbrtn theory argument list
c rgb 21-jan-96 changed common /cnvect/ from (55,12) to (55,9)
c rgb 11-apr-95 common blocks -> argument list
c   See rest of changes list at end of file
c--------1---------2---------3---------4---------5---------6---------7-c
!| 
!| Subroutine {\tt theory}'s calling argument,
!| {\tt nkthe}, is used to control the printout
!| (printout occurs if ${\tt nkthe} = 3$).
!| 
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
      subroutine ptheory(knthe)
c
c
cl      2.21  theoretical particle and energy fluxes
c
c
c  Input variables affecting sbrtn ptheory
c  ---------------------------------------
c
c
c  lthery(21) < 1 to call sbrtn theory
c             = 1 for sbrtn mmm95 rather than sbrtn theory
c             = 2 for sbrtn mmm98
c             = 3 for sbrtn mmm98b
c             = 4 for sbrtn mmm98c
c             = 5 for sbrtn mmm98d
c             = 6 for sbrtn ohe model
c             = 7 for sbrtn mixed_merba
c             = 8 for sbrtn mixed_model (Mixed Bohm/gyro-Bohm model)
c             = 9 for sbrtn mmm99
c             = 10 for sbrtn mmm2001
c             = 23 for sbrtn callglf2d (GLF23 model)
c
c  lthery(22) = 0 for default version of mmm9*
c             = 1 for supershot version of mmm9*
c
c  lthery(23) = 0 use major radius to geometric center of flux surfaces
c             = 1 use major radius to outboard edge of flux surfaces
c
c  lthery(24) = 1  to set desnsi = densh + densimp
c                    and grdni = grdnh + grdnz
c                    and grdpr = grdti + grdni + grdte + grdne
c
c  lthery(26) =   timestep for diagnostic output for etaw* model
c
c  lthery(27) > 0 to replace negative diffusivity with velocity
c
c  lthery(28) > 0 to smooth diffusivities lthery(28) times
c
c  lthery(29) = 0 to minimize diagnostic printout
c             > 0 larger values produce more diagnostic printout
c
c  lthery(30) = 0 take absolute value of gradient scale lengths
c             = 1 to retain sign of gradient scale lengths
c
c  lthery(31) = 1 to make r * ( gradient scale lengths ) monotonic
c                   near the magnetic axis by changing value at maxis+1
c
c  lthery(32) = 0 for the original gradient scale lengths
c             .gt. 0 for smoothing over lthery(32) orders (like lsmord)
c             .lt. 0 to smooth 1 / gradient scale lengths
c
c  lthery(33) = 0 print unsmoothed gradient scale lengths
c             .gt. 0 to print smoothed values
c
c    The following switches control mmm99, called when lthery(21) = 9
c
c  lthery(34) = 0 use sbrtn weiland14 from mmm95 for ITG mode in mmm99
c             = 1 use etaw17 from mmm98d for ITG in mmm99
c             = 2 use etaw18 from mmm99 for ITG in mmm99
c
c             = 0 use retuned GLF23 model (if lthery(21) = 23)
c             < 0 use old GLF23 model (if lthery(21) = 23)
c
c  lthery(35) = 0 use Guzdar-Drake resistive balloong mode
c                   from mmm95 in mmm99
c             = 1 use drift Alfven mode (Scott model) from MMM98d
c                   in mmm99
c
c  lthery(36) = 0 use old kinetic ballooning mode from MMM95 in mmm99
c             = 1 use new kinetic ballooning mode (by Redd) from MMM98d
c                   in mmm99
c
c  lthery(37) = 0 use Hahm-Burrel ExB shearing rate to substract
c                   from weiland model growth rates in mmm99
c             = 1 use Hamaguchi-Horton ExB shearing parameter
c                   to multiply all long-wavelength transport coeffs
c                   in mmm99
c
c  lthery(38) = 0 do not use any ETG model in mmm99
c             = 1 use Horton's ETG model in mmm99
c
cap
c  lthery(39) >=0 use zxithe, zvithe ... rather than difthi/velthi
c             < 0 use difthi/velthi rather than zxithe, zvithe ...
c
c  lthery(45) > 0 zthzmix(jz) = zthdmix(jz) for mixed_model
c
c-----------------------------------
c
c  cthery(129) 0.0 multiplier for flow shear rate wexbs
c
c  cthery(138) 0.0 multipleir for impurity transport
c
c     The following effects are turned on only when cthery(*) .gt. zepslon
c  cthery(50) 0.0  upper bound on L_{ne} / R_{major}
c  cthery(51) 0.0  upper bound on L_{ni} / R_{major}
c  cthery(52) 0.0  upper bound on L_{Te} / R_{major}
c  cthery(53) 0.0  upper bound on L_{Ti} / R_{major}
c  cthery(54) 0.0  upper bound on L_{p } / R_{major}
c
c  cthery(60) 0.0  limits local rate of change of ITG mode diffusities
c
c
c-----------------------------------
c
      use mixed_Bohm_gyro_Bohm
c
	  include 'cparm.m'
          include 'cbaldr.m'
          include 'commhd.m'
          include 'cd3he.m'
c
c..parameters for use in L-H mode determination
c..added by Boonyarit Chatthong (27/7/09)
c
          data mode /0/
c
          real zpheat, zbtorm, zrmajorm, zrminorm, zhydmass, znebar
     & , pth, totalmass
c
c..parameters for use in Kikuchi NTV model
c..added by Boonyarit Chatthong (26/9/11)
c
          real alpha_NTV, D_NTV, mu1_NTV, mu2_NTV, mu3_NTV, K1_NTV
     & , aspect_NTV, fc_NTV, g_NTV
c
      parameter ( kmatdim = 12 )
c
c  kmatdim = first dimension of matricies
c
      real zgrdne(mj), zgrdni(mj), zgrdnh(mj), zgrdnz(mj)
     & , zgrdte(mj), zgrdti(mj), zgrdpr(mj), zgrdq(mj)
c
c  Normalized gradients:
c
c  zgrdne(jz) = - R ( d n_e / d r ) / n_e
c  zgrdni(jz) = - R ( d n_i / d r ) / n_i
c     n_i = thermal ion density (sum over hydrogenic and impurity)
c  zgrdnh(jz) = - R ( d n_h / d r ) / n_h
c     n_h = thermal hydrogenic density (sum over hydrogenic species)
c  zgrdnz(jz) = - R ( d Z n_Z / d r ) / ( Z n_Z )
c     n_Z = thermal impurity density,  Z = average impurity charge
c           sumed over all impurities
c  zgrdte(jz) = - R ( d T_e / d r ) / T_e
c  zgrdti(jz) = - R ( d T_i / d r ) / T_i
c  zgrdpr(jz) = - R ( d p   / d r ) / p    for thermal pressure
c  zgrdq (jz) = R ( d q   / d r ) / q    related to magnetic shear
c
      real zsgrdne(mj), zsgrdni(mj), zsgrdnh(mj), zsgrdnz(mj)
     & , zsgrdte(mj), zsgrdti(mj), zsgrdpr(mj), zsgrdq(mj)
     & , zgfactr(mj)
c
c  zsgrd* are smoothed and preprocessed normalized gradients
c  zgfactr(jz) are factors used to convert from derivatives wrt
c    minor radius to derivatives wrt R_{\rm outboard}
c    when lthery(23) = 1
c
      real zgradrsqrave(mj), zgradrave(mj), zgradroutbrd(mj)
     &  ,  zdrdrho(mj)
c
c   Metric elements:
c     Here, r(jz) is the minor radius (half-width).
c  zgradrsqrave(jz) = <|grad r|^2>
c  zgradrave(jz)    = <|grad r|>
c  zgradroutbrd(jz) = |grad r| at outboard edge of each flux surface
c  zdrdrho(jz)      = d r / d rho
c
      integer iaxis, iedge, iseprtx, indim, iatdim, iprint
c
      real zrminor(mj), zrmajor(mj), zelong(mj), ztriang(mj)
     & , zindent(mj), zaimass(mj)
     & , zdense(mj), zdensi(mj), zdensh(mj), zdensf(mj), zdensfe(mj)
     & , zxzeff(mj), ztekev(mj), ztikev(mj), ztfkev(mj)
     & , zq(mj), zvloop(mj), zbtor(mj), zresist(mj)
     & , zwexbxb(mj), zyexbxb(mj)
     & , zvrotxb(mj), zvrotNTV(mj), zvrotNEO(mj), zvrottorque(mj)
     & , zdensimp(mj), zmassimp(mj), zavezimp(mj), zmasshyd(mj)
     & , zcharge_hyd(mj)
c
      real  zrhois(mj), zrsist(mj), ztcrit(mj)
     & , zslne(mj), zslni(mj), zslte(mj), zslti(mj), zsshr(mj)
     & , zrhohs(mj), zmlnh(mj), zlnhs(mj), zslnh(mj)
     & , zrhozs(mj), zlnzs(mj), zmlnz(mj), zslnz(mj)
     & , zlpr(mj), zslpr(mj)
c
      real
     &   zthiig(mj),   zthdig(mj),    ztheig(mj),    zthzig(mj)
     & , zthirb(mj),   zthdrb(mj),    ztherb(mj),    zthzrb(mj)
     & , zthikb(mj),   zthdkb(mj),    zthekb(mj),    zthzkb(mj)
     & , zthitem(mj),  zthdtem(mj),   zthetem(mj),   zthztem(mj)
     & , zthiitg(mj),  zthditg(mj),   ztheitg(mj),   zthzitg(mj)
     & ,                              ztheeg(mj)
     & ,                              zthetb(mj)
c
c
c  zthiig(jz) = ion thermal diffusivity from the Weiland model
c  zthdig(jz) = hydrogenic ion diffusivity from the Weiland model
c  ztheig(jz) = elelctron thermal diffusivity from the Weiland model
c  zthzig(jz) = impurity ion diffusivity from the Weiland model
c
c  zthirb(jz) = ion thermal diffusivity from resistive ballooning modes
c  zthdrb(jz) = hydrogenic ion diffusivity from resistive ballooning modes
c  ztherb(jz) = elelctron thermal diffusivity from resistive ballooning modes
c  zthzrb(jz) = impurity ion diffusivity from resistive ballooning modes
c
c  zthikb(jz) = ion thermal diffusivity from kinetic ballooning modes
c  zthdkb(jz) = hydrogenic ion diffusivity from kinetic ballooning modes
c  zthekb(jz) = elelctron thermal diffusivity from kinetic ballooning modes
c  zthzkb(jz) = impurity ion diffusivity from kinetic ballooning modes
c
c  zthitem(jz) = ion thermal diffusivity from TEM (diagnostic output)
c  zthdtem(jz) = hydrogenic ion diffusivity from TEM (diagnostic output)
c  zthetem(jz) = elelctron thermal diffusivity from TEM (diagnostic output)
c  zthztem(jz) = impurity ion diffusivity from TEM (diagnostic output)
c
c  zthiitg(jz) = ion thermal diffusivity from ITG (diagnostic output)
c  zthditg(jz) = hydrogenic ion diffusivity from ITG (diagnostic output)
c  ztheitg(jz) = elelctron thermal diffusivity from ITG (diagnostic output)
c  zthzitg(jz) = impurity ion diffusivity from ITG (diagnostic output)
c
c  ztheeg(jz) = electron thermal diffusivity from ETG modes
c  zthetb(jz) = electron thermal diffusivity from Taroni-Bohm
c
c  Arrays for time averaging when cthery(60) > 0.0
c
      real zoetai(mj), zoetae(mj), zoetad(mj), zoetaz(mj)
     &  ,  zdleti(mj), zdlete(mj), zdletd(mj), zdletz(mj)
c
      real zthigb(mj), zthegb(mj), zthibohm(mj), zthebohm(mj)
     &  ,  zthimix(mj), zthemix(mj), zthdmix(mj), zthzmix(mj)
c
      real  zgammaitg(mj), zomegaitg(mj), zgammatem(mj), zomegatem(mj)
c
c  zgammaitg(jr) = growth rate of fastest growing ITG mode
c  zomegaitg(jr) = frequency of fastest growing ITG mode ( < 0 )
c  zgammatem(jr) = growth rate of fastest growing TEM mode
c  zomegatem(jr) = frequency of fastest growing TEM mode ( > 0 )
c
      dimension zrmajm(110), znemaj(110), ztemaj(110), ztimaj(110)
     &  , zefmaj(110)
c
c   Arrays for printing out profiles along the major radius
c zrmajm(jz) major radius in meters
c znemaj(jz) electron density in m^-3
c ztemaj(jz) electron temperatur in keV
c ztimaj(jz) ion temperature in keV
c zefmaj(jz) Z_{eff}
c
c
      real ztemp1(mj), ztemp2(mj), ztemp3(mj), ztemp4(mj)
c
      real zstemp1, zstemp2, zstemp3, zstemp4, zstemp5, zstemp6
c
      real  zdhthe(mj), zvhthe(mj), zdzthe(mj), zvzthe(mj)
     &  , zxethe(mj), zxithe(mj), zweithe(mj)
     &  , zvithe(mj), zvethe(mj), vtorr(mj), boyle(mj)
c
c  zdhthe(jz)   = hydrogenic diffusivity ( m^2/sec )
c  zvhthe(jz)   = hydrogenic convective velocity ( m/sec )
c  zdzthe(jz)   = impurity diffusivity ( m^2/sec )
c  zvzthe(jz)   = impurity convective velocity ( m/sec )
c  zxethe(jz)   = electron thermal diffusivity ( m^2/sec )
c  zxithe(jz)   = ion thermal diffusivity ( m^2/sec )
c  zweithe(jz)  = anomalous electron-ion equipartition
c  zvithe(jz)   = ion thermal convective velocity (m/sec)
c  zvethe(jz)   = electron thermal convective velocity (m/sec)
c
c..control variables for Multi-Mode model mmm95
c
      integer lsuper, lreset, lmmm95(8), iswitch(8), ier
c
      real    cmmm95(29)
c
      real zgamma(kmatdim,55), zomega(kmatdim,55)
     &  , zvflux(kmatdim,55)
c
c..total diffusities for diagnostic output
c
      real zxetot(mj), zxitot(mj), zdhtot(mj), zdztot(mj)
c
c..effective convective velocities
c
        common /cnvect/ vftot(55,9), flxtot(55,9), srctot(55,9)
c
c  vftot(jz,ji)  = total effective convective velocities
c  flxtot(jz,ji) = total fluxes of particles and energies
c  srctot(jz,ji) = total sources of particles and energies
c    (these are computed in sbrtn convrt)
c
c..variables needed for GLF23 model
c
      real zbtor_axis, ztor_flux, zrmaj_axis
     &   , zsw_alpha, zsw_x_alpha
c
      real zsdifi, zschie, zschii, zsvtor, zsvpara, zsvperp
     &  , zsexchange
     &  , zscalgrdte,   zscalgrdti,   zscalgrdne,   zscalgrdni
c
      integer i_delay
c
c..variables in common blocks /comth*/
c
c                 electron/ion thermal diffusivity from:
c  thdre/i(j)   drift waves (trapped electron modes)
c  thtie/i(j)   ion temperature gradient (eta_i) modes
c  thrme/i(j)   rippling modes
c  thrbe/i(j)   resistive ballooning modes
c  thrbgb,thrbb(j)   gyro-Bohm and Bohm parts of thrbe(j)
c  thnme/i(j)   neoclassical MHD modes
c  thkbe/i(j)   kinetic ballooning modes
c  thhfe/i(j)   eta_e mode
c  thrlwe/i(j)  Rebut-Lallia-Watkins model
c
c  thdte(j)  = D_{te}
c  thdi(j)   = D_i
c
c  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c  velthi(j1,jz)    = convective velocities
c
c---------------------------------------------------------------------
c variables for vpolprof routine
c all the species indexes begin with electron
c---------------------------------------------------------------------
        integer, parameter :: maxs = 10         !maximum species
        integer, parameter :: maxx = 55         !maximum x points
        real zmasspecies(maxs)                  !mass of species
        real zcharge(maxs)                      !charge of species
        real zdensity(maxs, maxx)               !density profiles of
                                                !all species
        real ztemperature(maxs, maxx)           !temperature profiles
                                                !of all species
        real zvpolxb(maxs, maxx)                !poloidal velocity
                                                !profiles of all species
        real zvpol(maxx)                        ! hydrogenic poloidal velocity
        integer imod                            !poloidal velocity
                                                !model index
                                                !1: Strand
                                                !2: Erba
                                                !3: ZHS
c
c..arrays needed for the GLF23 model
c
        real zvtor(maxx)                        ! toroidal velocity (m/s)
        real zvpara(maxx)                       ! parallel velocity (m/s)
        real zvperp(maxx)                       ! perp velocity (m/s)
        real zdqdrho(maxx)   ! rho d q / d rho
        real zalpha(maxx)    ! 2 mu0 R q^2 ( d p / d r ) / B^2
        real zthvtor(maxx)   ! toroidal velocity diffusivity (m^2/s)
        real zthvpara(maxx)  ! parallel velocity diffusivity (m^2/s)
        real zthvperp(maxx)  ! perpendicular velocity diffusivity (m^2/s)
        real zexb_out1(maxx) ! exb shear rate output 1
        real, dimension(:,:), allocatable ::  zexb_out2
              ! exb shear rate output 2 (time delayed)
        real zexb_out3(maxx) ! parallel velocity shear rate output 
c
c---------------------------------------------------------------------
cap      real feg(4)     !multiplier for ETG mode
c
!| 
!| The full matrix form of anomalous transport has the form
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
!| 
c
c    Lengths:
c  thlni(j)  = L_{ni}
c  thlti(j)  = L_{T_i}
c  thlsh(j)  = L_s = R q / s\hat
c  thlpr(j)  = L_p
c  thlarp(j) = \rho_s
c  thrhos(j) = \rho_{\theta i}
c
c    Velocities:
c  thvthe(j) = v_{the}
c  thvthi(j) = v_{thi}
c  thsoun(j) = c_s
c  thalfv(j) = v_A
c
c     Dimensionless:
c  threti(j) = \eta_i
c  thdias(j) = resistive ballooning mode diamagnetic stabilization factor
c  thdinu(j) = \omega_e^\ast / \nu_{eff}
c  thfith(j) = f_{ith} as in eq (35)
c  thbpbc(j) = \beta^{\prime} / \beta_{c1}^{\prime}
c  thdia(j)  = \omega_e^\ast / k_\perp^2
c  thnust(j) = \nu_e^*
c  thlamb(j) = \Lambda multiplier for resistive ballooning modes
c
c  thbeta(j) = \beta
c  thetth(j) = \eta_i^{th}  threshold for \eta_i mode
c  thsrhp(j) = S = \tau_R / \tau_{hp} = r^2 \mu_0 v_A / ( \eta R_0 )
c  thvalh(jz)= Parameter determining validity of Hahm-Tang CTEM model
c
c-----------------------------------------------------------------------
!| 
!| The coding continues with the
!| OLYMPUS  ({\it cf.} \cite{BALDUR}) number (2.21),
!| and use of the OLYMPUS form for bypassing subroutines,
!| and comments on the common blocks and variables modified.
c
      data  istep /1/, itdiag /1/
c
      data               iclass /2/   , isub /21/
c
      save istep, itdiag
c
      if ( nlomt2(isub) ) then
        if (nstep .lt. 2)
     & call mesage(' *** 2.21 subroutine theory bypassed            ')
        return
      endif
c
c-----------------------------------------------------------------------
c
c
cl       common blocks and variables modified
c
c  eithes(jz), dxthes(ix,jz), vxthes(ix,jz), xethes(jz), xithes(jz),
c  weithe(jz), weiths(jz)  in /comthe/
c
c-----------------------------------------------------------------------
c
c..go to printout directly if knthe = 3
c
cbate      if ( knthe .eq. 3 ) go to 800
c
c..initialize arrays
c
      do 10 jz=1,mzones
        eithes(jz)    = 0.
        xethes(jz)    = 0.
        xithes(jz)    = 0.
        weithe(jz)    = 0.
        weiths(jz)    = 0.
        thdre(jz)     = 0.
        thdri(jz)     = 0.
        thige(jz)     = 0.
        thigi(jz)     = 0.
        thtie(jz)     = 0.
        thtii(jz)     = 0.
        thrme(jz)     = 0.
        thrmi(jz)     = 0.
        thrbgb(jz)    = 0.
        thrbb(jz)     = 0.
        thrbe(jz)     = 0.
        thnme(jz)     = 0.
        thnmi(jz)     = 0.
        thcee(jz)     = 0.
        thcei(jz)     = 0.
        thhme(jz)     = 0.
        thrbi(jz)     = 0.
        thkbe(jz)     = 0.
        thkbi(jz)     = 0.
        thhfe(jz)     = 0.
        thhfi(jz)     = 0.
        thrlwe(jz)    = 0.
        thrlwi(jz)    = 0.
        ztemp1(jz)    = 0.
        ztemp2(jz)    = 0.
        zthigb(jz)    = 0.
        zthegb(jz)    = 0.
        zthibohm(jz)  = 0.
        zthebohm(jz)  = 0.
        zthimix(jz)   = 0.
        zthemix(jz)   = 0.
        zthdmix(jz)   = 0.
        zthzmix(jz)   = 0.
        zdhthe(jz)    = 0.
        zvhthe(jz)    = 0.
        zdzthe(jz)    = 0.
        zvzthe(jz)    = 0.
        zxethe(jz)    = 0.
        zxithe(jz)    = 0.
        zweithe(jz)   = 0.
        zvithe(jz)    = 0.
        zvethe(jz)    = 0.
  10  continue
!ap{
      zgamma  = 0. 
      zxethe  = 0.
      zxithe  = 0.
      zthiig  = 0.
      zthirb  = 0.
      zthikb  = 0.
      ztherb  = 0.
      zthekb  = 0.
      ztheig  = 0.
      ztheeg  = 0.
      zthetb  = 0.
      zthdig  = 0.
      zthdrb  = 0.
      zthdkb  = 0.
      zthzig  = 0.
      zthzrb  = 0.
      zthzkb  = 0.
c
      zthitem = 0.
      zthdtem = 0.
      zthetem = 0.
      zthztem = 0.
c
      zthiitg = 0.
      zthditg = 0.
      ztheitg = 0.
      zthzitg = 0.
c
      zthvtor  = 0.
      zthvpara = 0.
      zthvperp = 0.
      zweithe  = 0.
c
      zomegaitg = 0.
      zomegatem = 0.
      zgammaitg = 0.
      zgammatem = 0.
c
      zgfactr   = 0.
      zelong    = 0.
      ztriang   = 0.
      zdhthe    = 0.
      zvhthe    = 0.
      zvzthe    = 0.
!}
c
      do 12 jz=1,mzones
        do 12 ix=1,6
          dxthes(ix,jz) = 0.
          vxthes(ix,jz) = 0.
  12  continue
c
c  converting the physical constants into local variables
c
      zcmu0=emu0
      zceps0=eps0
      zckb=cfev*10.**(fxk)
      zcme=fcme*10.**(fxme)*usim
      zcmp=fcmp*10.**(fxnucl)*usim
      zce=fce*10.**(fxe)*10.
      zcc=1./sqrt(zcmu0*zceps0)
c
      zepsqrt = sqrt ( epslon )
c
      zlgeps=log(epslon)
c
      nerr = 0
c
c
!| 
!| There is a problem with using the electron density scale length $L_{ne}$
!| in the computations below because the electron density is influenced by
!| fast beam ions which are subject to Monte-Carlo noise.
!| (Hence, $L_{ne}$ is observed to have random spatial and temporal fluctuations
!| due to this numerical artifact in the BALDUR and similar codes,
!| especially near the magnetic axis).
!| It is preferable, therefore, to compute the density scale length using
!| the thermal ion density ({\tt rhoins(2,jz)} in the BALDUR code).
!| This is computed following the same methods used in BALDUR subroutine
!| XSCALE.
!| As an alternative, {\tt zrhohs(jz)} is computed from the sum of the
!| thermal hydrogen densities in standard units
!| and {\tt zlnhs(jz)} is the gradient scale length of {\tt zrhohs(jz)}.
!| 
!| Note $ L_{n_Z} = {tt zlnzs} = Z n_Z / [ d ( Z n_Z ) / d r ] $.
!| 
c
c..compute thermal ion density scale length
c
c  zrhois(jz) = rhoins(2,jz)
c  zrhohs(jz) = hydrogen density
c  zrhozs(jz) = sum of Z_i n_i for impurities
c  rhisms(jz) = smoothed array of rhoins(2,jz)
c  slnis (jz) = 1. / ( d ln (rhoins(2,jz)) / d r )
c
c  temporarily put zrhohs in zmlnh and put zrhozs in zmlnz
c
      do jz=1,mzones
        zrhois(jz) = rhoins(2,jz)
        zmlnh(jz) = 0.0
        zmlnz(jz) = 0.0
        do js=1,mhyd
          zmlnh(jz) = zmlnh(jz) + rhohs(js,2,jz)
        enddo
        do ji=1,mimp
          zmlnz(jz) = zmlnz(jz) + cmean(ji,2,jz) * rhois(ji,2,jz)
        enddo
      enddo
c
c  Smooth profiles, filling arrays rhisms, zrhohs, and zrhozs
c
      call smooth (zrhois,rhisms,mzones,2,1,1,lcentr,dx2i,
     &             smrlow,smlwcy,lsmord)
c
      call smooth (zmlnh,zrhohs,mzones,2,1,1,lcentr,dx2i,
     &             smrlow,smlwcy,lsmord)
c
      call smooth (zmlnz,zrhozs,mzones,2,1,1,lcentr,dx2i,
     &             smrlow,smlwcy,lsmord)
!ap{
      rhisms(lcentr-1) = rhisms(lcentr)
      zrhohs(lcentr-1) = zrhohs(lcentr)
      zrhozs(lcentr-1) = zrhozs(lcentr)
!}
c
      do jz=1,mzones
        zmlnh(jz) = 0.0
        zmlnz(jz) = 0.0
      enddo
c
      zsndni = -1.
      do jz=lcentr+1,mzones
        if ( abs ( rhisms(jz) - rhisms(lcentr) ) .gt.
     &       abs ( epslon * rhisms(lcentr) ) ) then
          zsndni = ( rhisms(jz) - rhisms(lcentr) )
     &             / abs ( rhisms(jz) - rhisms(lcentr) )
          go to 24
        endif
      enddo
  24  continue
c
      zsndnh = -1.
      do jz=lcentr+1,mzones
        if ( abs ( zrhohs(jz) - zrhohs(lcentr) ) .gt.
     &       abs ( epslon * zrhohs(lcentr) ) ) then
          zsndnh = ( zrhohs(jz) - zrhohs(lcentr) )
     &             / abs ( zrhohs(jz) - zrhohs(lcentr) )
          go to 25
        endif
      enddo
  25  continue
c
      zsndnz = -1.
      do jz=lcentr+1,mzones
        if ( abs ( zrhozs(jz) - zrhozs(lcentr) ) .gt.
     &       abs ( epslon * zrhozs(lcentr) ) ) then
          zsndnz = ( zrhozs(jz) - zrhozs(lcentr) )
     &             / abs ( zrhozs(jz) - zrhozs(lcentr) )
          go to 26
        endif
      enddo
  26  continue
c
      do jz=lcentr,mzones
        slnis(jz) = -1. / ( epslon * zsndni
     &    + ( rhisms(jz) - rhisms(jz-1) )
     &    / ( ( armins(jz,2) - armins(jz-1,2) )
     &        * 0.5 * ( rhisms(jz) + rhisms(jz-1) ) ) )
        zlnhs(jz) = -1. / ( epslon * zsndnh
     &    + ( zrhohs(jz) - zrhohs(jz-1) )
     &    / ( ( armins(jz,2) - armins(jz-1,2) )
     &        * 0.5 * ( zrhohs(jz) + zrhohs(jz-1) ) ) )
        zlnzs(jz) = -1. / ( epslon * zsndnz
     &    + ( zrhozs(jz) - zrhozs(jz-1) )
     &    / ( ( armins(jz,2) - armins(jz-1,2) )
     &        * 0.5 * ( zrhozs(jz) + zrhozs(jz-1) ) ) )
      enddo
!{ap
      slnis(lcentr-1) = slnis(lcentr)
      zlnhs(lcentr-1) = zlnhs(lcentr)
      zlnzs(lcentr-1) = zlnzs(lcentr) 
!}      
c
!| 
!| Set up arrays in MKS units (keV for temperatures)
!| for argument list of sbrtn theory.
!| 
!| 
c
        do 32 jz=1,mzones
          zslne(jz) = slnes(jz) * usil
          zslni(jz) = slnis(jz) * usil
          zslnh(jz) = zlnhs(jz) * usil
          zslnz(jz) = zlnzs(jz) * usil
          zslte(jz) = sltes(jz) * usil
          zslti(jz) = sltis(jz) * usil
          zslpr(jz) = (armins(jz,1)/max(epslon,slprs(jz)))*usil
          zsshr(jz) = shear(jz)
  32    continue
c
!| If {\tt lthery(23) = 1,} use
!| \[ R_{\rm outboard} = R_{\rm Geometric Center} + r \]
!| as the major radius.  Also, take all gradients with respect to
!| $ R_{\rm outboard} $.  Hence
!| \[ d \]
c
c
      iaxis = lcentr
      iedge = mzones
      isep = mzones
      if (nadump(1) .gt. lcentr) isep = nadump(1)
c
      do jz=1,mzones
        zelong(jz)  = elong(jz,1)
        ztriang(jz) = triang(jz,1)
        zindent(jz) = dent(jz,1)
        zaimass(jz) = aimass(1,jz)
        zdense(jz)  = rhoels(1,jz)*usid
        zdensi(jz)  = rhoins(1,jz)*usid
        zdensh(jz)  = zrhohs(jz) * usid
        ztekev(jz)  = tes(1,jz)*useh
        ztikev(jz)  = tis(1,jz)*useh
        ztfkev(jz)  = 0.0
        zxzeff(jz)  = xzeff(1,jz)
       zlne = zslne(jz)
       zlni = zslni(jz)
       zlnh = zslnh(jz)
       zlnz = zslnz(jz)
       zlte = zslte(jz)
       zlti = zslti(jz)
ces   temporary numerical overflow protection, cf. statement 22 below
       zlpr(jz)=(armins(jz,1)/max(epslon,slprs(jz)))*usil
      zshear = zsshr(jz)
      zrsep  = max ( armins(isep,1), epslon ) * usil
      zrminor(jz) = max ( armins(jz,1), epslon ) * usil
      zrmajor(jz) = armajs(jz,1)*usil
      zvloop(jz)  = abs(vloopi(jz,2))
      zbtor(jz)   = bzs*usib
      zresist(jz) = eta(1,jz) * usir
c
      zdensf(jz) = ( rhobis(1,jz) + rh1fst(1,jz) + rh2fst(1,jz) )
     &  * usid
c
      zdensfe(jz) = ( rhobes(1,jz) + rh1fst(1,jz) + 2.0 * rh2fst(1,jz) )
     &  * usid
c
        znz    = 0.0
        zmass  = 0.0
        zimpz  = 0.0
c
        if ( mimp .gt. 0 ) then
          do jimp=1,mimp
            znz    = znz   + rhois(jimp,1,jz)
            zimpz  = zimpz + rhois(jimp,1,jz) * cmean(jimp,1,jz)
            zmass  = zmass + rhois(jimp,1,jz) * aspec(lhydn+jimp)
          enddo
            zimpz  = zimpz / znz
            zmass  = zmass / znz
            znz    = znz * usid
        endif
c
          zdensimp(jz) = znz
          zmassimp(jz) = max ( zmass, 1.0 )
          zavezimp(jz) = max ( zimpz, 1.0 )
c
        zmasshyd(jz) = ahmean(1,jz)
c
      enddo
c
c..set up gradient scale length arrays
c
      do jz=1,mzones
        zgrdne(jz) = zrmajor(jz)
     &    / sign(max(abs(zslne(jz)),epslon),zslne(jz))
        zgrdni(jz) = zrmajor(jz)
     &    / sign(max(abs(zslni(jz)),epslon),zslni(jz))
        zgrdnh(jz) = zrmajor(jz)
     &    / sign(max(abs(zslnh(jz)),epslon),zslnh(jz))
        zgrdnz(jz) = zrmajor(jz)
     &    / sign(max(abs(zslnz(jz)),epslon),zslnz(jz))
        zgrdte(jz) = zrmajor(jz)
     &    / sign(max(abs(zslte(jz)),epslon),zslte(jz))
        zgrdti(jz) = zrmajor(jz)
     &    / sign(max(abs(zslti(jz)),epslon),zslti(jz))
        zgrdpr(jz) = zrmajor(jz)
     &    / sign(max(abs(zslpr(jz)),epslon),zslpr(jz))
        zgrdq(jz)  = zsshr(jz) * zrmajor(jz) / zrminor(jz)
      enddo
c
!| If {\tt lthery(23) = 1,} use
!| \[ R_{\rm outboard} = R_{\rm Geometric Center} + r \]
!| as the major radius.  Also, take all gradients with respect to
!| $ R_{\rm outboard} $.  Hence
!| \[ d p / d r \rightarrow d p / d R_{\rm outboard}
!|  = ( d p / d r ) / ( d R_{\rm outboard} / d r )
!|  = ( d p / d r ) / ( 1.0 + d R_{\rm Geometric Center} / d r )  \]
c
c..if lthery(23) = 1, use
c     R_{\rm outboard} = R_{\rm Geometric Center} + r  for major radius
c     and take all gradients wrt R_{\rm outboard}
c
      if ( lthery(23) .eq. 1 ) then
        do jz=1,mzones
          zrmajor(jz) = ( armajs(jz,1) + armins(jz,1) )*usil
        enddo
c
        zgfactr(1) = 0.0
c
        do jz=2,mzones
          zgfactr(jz) =
     &      ( armajs(jz,1) + armins(jz,1) ) / ( armajs(jz,1)
     &      * max ( 0.001,  ( 1.0
     &        + ( armajs(jz,2)-armajs(jz-1,2) )
     &          / ( armins(jz,2) - armins(jz-1,2) ) ) ) )
        enddo
c
        do jz=1,mzones
          zgrdne(jz) = zgrdne(jz) * zgfactr(jz)
          zgrdni(jz) = zgrdni(jz) * zgfactr(jz)
          zgrdnh(jz) = zgrdnh(jz) * zgfactr(jz)
          zgrdnz(jz) = zgrdnz(jz) * zgfactr(jz)
          zgrdte(jz) = zgrdte(jz) * zgfactr(jz)
          zgrdti(jz) = zgrdti(jz) * zgfactr(jz)
          zgrdpr(jz) = zgrdpr(jz) * zgfactr(jz)
          zgrdq(jz)  = zgrdq(jz)  * zgfactr(jz)
        enddo
c
      endif
c
c..gradients have odd symmetry across the magnetic axis
c
      zgrdne(iaxis) = 0.0
      zgrdni(iaxis) = 0.0
      zgrdnh(iaxis) = 0.0
      zgrdnz(iaxis) = 0.0
      zgrdte(iaxis) = 0.0
      zgrdti(iaxis) = 0.0
      zgrdpr(iaxis) = 0.0
      zgrdq(iaxis)  = 0.0
c
      if ( iaxis .gt. 1 ) then
        do j=1,iaxis-1
          zgrdne(iaxis-j) = - zgrdne(iaxis+j)
          zgrdni(iaxis-j) = - zgrdni(iaxis+j)
          zgrdnh(iaxis-j) = - zgrdnh(iaxis+j)
          zgrdnz(iaxis-j) = - zgrdnz(iaxis+j)
          zgrdte(iaxis-j) = - zgrdte(iaxis+j)
          zgrdti(iaxis-j) = - zgrdti(iaxis+j)
          zgrdpr(iaxis-j) = - zgrdpr(iaxis+j)
          zgrdq(iaxis-j)  = - zgrdq(iaxis+j)
        enddo
      endif
c
c..if lthery(24) = 1, recompute zdensi, zgrdni, and zgrdpr
c
      if ( lthery(24) .eq. 1 ) then
c
        do jz=1,mzones
          zdensi(jz) = zdensh(jz) + zdensimp(jz)
        enddo
c
      endif
c
!| 
!| After setting up the variables in the argument list, call sbrtn theory.
!| 
c
 800  continue
c
      ztime  = tai * uist
      indim   = 6
      imatdim = 12
c
c Interpolate the ExB shearing rate on to the xb grid at time ztime
c
c pis 1-jun-98: Current guess of xb = xbouni, and nxb = mzones
c pis 6-jul-98: Implemented test for vrota = 0
c


      kxdim = 55 ! this has to be passed down from somewhere!!!

      zvrotxb = 0.0
      zwexbxb = 0.0
      zyexbxb = 0.0
c
c..include w_ExB calculation by Boonyarit Chatthong (27/7/09)
c  lthery(48) = 0,1   using experimental w_ExB
c  lthery(48) = 2     using calculated w_ExB by experimental vtor
c  lthery(48) = 3     using calculated w_ExB by vtor = c*Ti
c  lthery(48) = 4     using calculated w_ExB by vtor = c*vtor,0*Ti/Ti,0
c  lthery(48) = 5     using calculated w_ExB by vtor = Jtor/(e*n)
c
c  Add lthery(48) = 6 to include Kikuchi's model by Boonyarit Chatthong (25/9/11)
c  modified on 02/01/12
c  Add lthery(48) = 7 to include Neoclassical rotations (02/03/12)
c  lthery(48) = 6     using calculated vtor by kikuchi's NTV model
c  PUB modified 'lthery(48)=8' from PON (13/03/15)
c  lthery(48) = 8     using calculated vtor by NTV model (NO torque models)
c
      if ((lthery(48) .eq. 0) .or. (lthery(48) .eq. 1)) then
c
c..experimental w_ExB
c
         call wexbint (kxdim, nxwexba, xwexba, ntwexba, twexba, wexba
     & , mzones, xbouni, ztime, zwexbxb)
c
         write(nprint,*) 'Experimental w_ExB'
c
      else if (lthery(48) .ge. 2) then
c..predicted w_ExB
         if (lthery(48) .eq. 2) then
c
c..experimental vtor
c
            call wexbint (kxdim, nxwexba, xwexba, ntwexba, twexba, vrota
     &                 , mzones, xbouni, ztime, zvrotxb)
c
            write(nprint,*) 'Predicted w_ExB by experimental vtor'
c     
         else if ((lthery(48) .eq. 3) .or. (lthery(48) .eq. 4) .or.
     &   (lthery(48) .eq. 5)) then
c
c..start determining L-H mode
c..check status of plasmas
c
c
           zeheat = usip * ( geohms(mzones) + gealfs(mzones)
     &     + geauxs(mzones) + gebems(mzones) + geecrs(mzones)
     &     + geicrs(mzones))
c
           ziheat = usip * ( gialfs(mzones) + giauxs(mzones)
     &     + gibems(mzones) + giecrs(mzones) + giicrs(mzones) )
c
c..calculate the total heating power
c
           zpheat = 1.e-6 * (zeheat + ziheat)
c
c..setup the input for L-H model 
c
           zrmajorm   = max (rtmajb, 0.01)
           zrminorm   = max (rtminb, 0.01) 
           zbtorm     = max (bzs * 1.0E-4, 0.01)  
           zhydmass  = max (ahmean(1,mzones), 1.00)
c
c..calculate average electron density
c
           if (lbound(6) .eq. 2) then
              znebar = denmont * 1.e6
           else
              znebar = enlaes(mzones) * 1.e6
           endif
c
c..determine L-H transition by using power threadhold (P > P_LH)
c
c  mode = 0  for L-mode
c       = 1  for H-mode 
c
           ztime = 0.5 * ( tai + tbi ) * uiet	  
c
           if (((cbound(2)>epslon) .and. (nstep>50)) 
     &         .and. ((ztime .ge. cbound(2)))) then
c
             if (lbound(1) ==1 ) then
                call bdlhmode (
     &           lbound,     cbound,     zpheat,     zrmajorm,
     &           zrminorm,    zbtorm,      zhydmass,   znebar,
     &           pth,        mode)
             endif
c
           else   
             mode = 0
           endif
c
c..predicted vtor
c     
            if (mode .eq. 0) then
c
               do jz = 1, mzones
                  zvrotxb(jz) = 0
               end do
               write(nprint,*) 'Predicted w_ExB in L-mode'
c     
            else if ((mode .eq. 1) .and. (lthery(48) .eq. 3)) then
c
c..using predicted vtor = c*Ti
c
               do jz = 1, mzones
                  zvrotxb(jz) = cthery(140)*10000.0*ztikev(jz)
               end do
               write(nprint,*) 'Predicted w_ExB by Ti in H-mode'
c               write(nprint,'(A,E8.3)') 'ccccc is', cthery(140)
c
            else if((mode .eq. 1) .and. (lthery(48) .eq. 4)) then
c
c..using predicted vtor =  c*vtor,0*Ti/Ti,0
c
               do jz = 1, mzones
                  zvrotxb(jz) = ((-1340*(((hpowmw(1)+hpowmw(2))/2.770144
     &            )+((hpowmw(3)+hpowmw(4))/3.664556)))+29826)*ztikev(jz)
c     &             35304*((zrmajor(mzones)/(znebar/
c     &            1.0E19)*(((hpowmw(1)+hpowmw(2))/2.770144)+((hpowmw(3)
c     &            +hpowmw(4))/3.664556)))**0.61)*((eqcamp/1.0E6)**0.3)
c     &            *ztikev(jz)/ztikev(1)
               end do
               write(nprint,*) 'Predicted w_ExB by r/a in H-mode'
c               write(nprint,'(A,E8.3)') 'kkkkk is', zrmajor(mzones)
c               write(nprint,'(A,E8.3)') 'kkkkk is', ztikev(1)
c               write(nprint,'(A,E8.3)') 'kkkkk is', ztikev(mzones)
c               write(nprint,'(A,E8.3)') 'kkkkk is', zvrotxb(1)
c
            else if((mode .eq. 1) .and. (lthery(48) .eq. 5)) then
               do jz = 1, mzones
                  zvrotxb(jz) = znjtor(jz)*1.0E6/(1.602176E-19*
     &            (zdensi(jz)*zxzeff(jz)))
               end do
               write(nprint,*) 'Predicted w_ExB by Jtor/(n*e)'
c
            end if
c
         else if (lthery(48) .eq. 6) then
c..kikuchi NTV model
            aspect_NTV = zrminor(mzones)/zrmajor(1)
            fc_NTV = 1.-1.46*sqrt(aspect_NTV)+0.46*aspect_NTV*
     &      sqrt(aspect_NTV)
            g_NTV = (1.-fc_NTV)/fc_NTV
c     
            do jz = 1, mzones
               alpha_NTV = zdensimp(jz)*zcharge(3)*zcharge(3)/
     &         (zdensi(jz)*zcharge(2)*zcharge(2))
               mu1_NTV = g_NTV*(alpha_NTV+sqrt(2.)-log(1.+sqrt(2.)))
               mu2_NTV = g_NTV*(3.*alpha_NTV/2. + 4./sqrt(2.) - 
     &              5.*log(1.+sqrt(2.))/2.)
               mu3_NTV = g_NTV*(13.*alpha_NTV/4. + 39./(4.*sqrt(2.))
     &              - 25.*log(1.+sqrt(2.))/4.)
               D_NTV = mu1_NTV*(mu3_NTV+sqrt(2.)+alpha_NTV)-
     &              mu2_NTV*mu2_NTV
               K1_NTV = (sqrt(2.)+alpha_NTV)*mu2_NTV/D_NTV
               zvrotNTV(jz) = (3.54-K1_NTV)*
     &              zgrdti(jz)*ztikev(jz)*1000./(zrmajor(jz)*zcharge(2)*
     &              bpoli(jz))
               zvrotxb(jz) = zvrotNTV(jz)
c               zvrotxb(jz) = znjtor(jz)*1.0E6/(1.602176E-19*(zdensi(jz)
c     &         *zxzeff(jz)))
c               zvrotxb(jz) = znjtor(jz)*1.0E6/(1.602176E-19*(zdensimp(jz
c     &         )*12 + zdensh(jz)*2*sqrt(zmasspecies(3)/zmasspecies(2)) +
c     &         zdense(jz)*sqrt(zmasspecies(3)/zmasspecies(1))))
c               write(nprint,'(A,E8.3)') 'kekkekkekzrmajo is', zgrdti(jz)
c               write(nprint,'(A,E8.3)') 'kekkekkekk1nt is', ztikev(jz)
c               zvrotxb(jz) = znjtor(jz)*1.0E6/(1.602176E-19*zdensi(jz)
c     &         *zxzeff(jz))
c               write(nprint,'(A,E8.3)') 'kekkekkek is', zmasspecies(1)
c               write(nprint,'(A,E8.3)') 'kekkekkek is', zmasspecies(2)
c               write(nprint,'(A,E8.3)') 'kekkekkek is', zmasspecies(3)
c               write(nprint,'(A,E8.3)') 'kekkekkek is', totalmass
c             write(nprint,'(A,E8.3)') 'kekcurrent is', znjtor(jz)
c             write(nprint,'(A,E8.3)') 'kekvelocity is', zvrotxb(jz)
c             write(nprint,'(A,E8.3)') 'kekzeff is', zxzeff(jz)
c             write(nprint,'(A,E8.3)') 'kekzdense is', zdense(jz)
c             write(nprint,'(A,E8.3)') 'kekzdensi is', zdensi(jz)
            end do
c            write(nprint,'(A,E8.3)') 'kekkekkek1 is', aspect_NTV
c            write(nprint,'(A,E8.3)') 'kekkekkek2 is', fc_NTV
c            write(nprint,'(A,E8.3)') 'kekkekkek3 is', g_NTV
            write(nprint,*) 'Predicted w_ExB by Kikuchi NTV model'
c
         else if (lthery(48) .eq. 8) then
c  Rotations  model with NTV model (NO torque models)
            aspect_NTV = zrminor(mzones)/zrmajor(1)
            fc_NTV = 1.-1.46*sqrt(aspect_NTV)+0.46*aspect_NTV*
     &      sqrt(aspect_NTV)
            g_NTV = (1.-fc_NTV)/fc_NTV
c  
c
c  ajboot(j)    = < bootstrap current density / R >     [A/m**3]
c  ajtpbi(j)    = < current density due to trapped
c                 particles in the banana regime
c                 / R > at zone center j                [A/m**3]
c  zrad         = plasma half-width                     [m]
c  eta(2,j)     = resistivity at zone centers           [sec]
c  zjtor        = <J_tor>                               [MA/m**2]
c  zcdrive      = beam and RF driven current density    [MA/m**2]
c  zboot        = bootstrap current density             [MA/m**2]
c  zjtpb        = <current density due to trapped
c                 particles
c  zconv        = convective contribution to
c                 current density due to motion
c                 of toroidal flux                      [Ma/m**2]
c               = \Partial{\psi_{pol}}{\rho}
c                 \Partial{\rho}{t}
c                 \langle 1/R^2 \rangle
c                 / 2 \pi \eta \langle 1/R \rangle
c  avi(j,9,1)   = R*B_tor                               [m*tesla]
c  z7r          = <1/R>                                 [1./m]
c  avi(j,10,1)  = <1/R**2>                              [1./m**2]
c  avi(j,7,1)   = < |del xi|**2 / R**2 >                [1./m**4]
c
c  zitorr = total toroidal current                      [MA]
c  zibeam = total beam-driven current                   [MA]
c  ziboot = total bootstrap current                     [MA]
c  zijtpb = total current due to trapped particles      [MA]
c
c  Note:  if ( lneocl(1) == 1 ) then
c    ajboot and ajtpbi are computed in sbrtn nclass_int
!cap
        if ( lneocl(1) .lt. 1 ) then
          ajboot = 0.0
          ajtpbi = 0.0
          call boots  ( ajboot,ajtpbi,lcentr,ledge)
        endif
c
        zibeam = 0.0
        ziboot = 0.0
        zijtpb = 0.0
c
c  zdarea = differential cross-sectional area [m**2]
c
        do 680 j=1,mzones-1
c
          zdarea = avi(j+1,13,1) - avi(j,13,1)
c
        zrad = 0.5 * ( avi(j+1,15,1) + avi(j,15,1) )
c
c        zjb7b = 1.e-6 *  avi(j,9,1) *
c     & ( avi(j+1,2,1) * avi(j+1,3,1) * r0ref * bpoli(j+1)
c     &          * avi(j+1,7,1) / avi(j+1,8,1)
c     & - avi(j,2,1) * avi(j,3,1) * r0ref * bpoli(j)
c     &          * avi(j,7,1)/avi(j,8,1) )
c     & / ( emu0 * ( avi(j+1,12,1) - avi(j,12,1) ) )
c     & / ( 0.5  * ( avi(j+1,11,1) + avi(j,11,1) ) )
c
c        zib7bt = zib7bt + zjb7b * zdarea
c
      zjtor = 1.e-6 * r0ref *
     & ( avi(j+1,2,1) * avi(j+1,3,1) * bpoli(j+1) * avi(j+1,7,1)
     &   - avi(j,2,1) * avi(j,3,1) * bpoli(j) * avi(j,7,1) )
     &   / ( twopi * emu0 * zdarea )
c
        zcdrive = 1.e-6 * ( cjbeam * ajbs(j) * usij
     &            + cdprof(j) )
        zibeam = zibeam + zcdrive * zdarea
c
        zboot = 1.e-6 * cfutz(480) * ajboot(j)
     &          / ( 0.5 * ( avi(j+1,11,1) + avi(j,11,1) ) )
        ziboot = ziboot + zboot * zdarea
c
        zjtpb = 1.e-6 * cfutz(480) * sfutz(16) * ajtpbi(j)
     &          / ( 0.5 * ( avi(j+1,11,1) + avi(j,11,1) ) )
        zijtpb = zijtpb + zjtpb * zdarea
c
        zconv = 1.e-6 * avi(j,10,1) * r0ref * (1. / eta(2,j))
     & * 0.5 * ( bpoli(j+1) * avi(j+1,1,2) / avi(j+1,11,1)
     &         + bpoli(j)   * avi(j,1,2)   / avi(j,11,1) )
c
        z7r = 0.5 * ( avi(j+1,11,1) + avi(j,11,1) )
c        write(6,*) 'ponkris NBI'
cponkris        
c        boyle = zcdrive       
c
c        write (nprint,10658) j, zrad, eta(2,j)*usir, 1./ftrap(2,j)
c     & , zjtor, zcdrive, zboot, zjtpb, zconv, avi(j,9,1),  z7r
c     & , avi(j,7,1)
     
c        write(6,*) 'ponkris NBI', boyle

 680  continue
c
c  zitorr = total toroidal plasma current in MA
c
      zitorr = avi(mzones,2,1) * avi(mzones,3,1) * r0ref * bpoli(mzones)
     &  * avi(mzones,7,1) * 1.e-6 / (twopi * emu0)
c
c        write (nprint,10659) zitorr, zibeam, ziboot, zijtpb
cpon
c
c       if helium recycling is active, print relevant quantities
c  15.07 this is not applicable to the IRE code.
c
cbate        if (cfutz(iheflx).gt.epslon.and.(natomc.ne.3)) call impprt
c
c       if neutral-impurity influx is subject to feedback readjustment,
c       print readjustment factor
c  15.07 this is not applicable to the IRE code.
c
c
c        data    iclass /3/,     isub /3/
        data    iheflx /70/,    impneu /200/
         if(nadump(1).le.lcentr) then
          isep   = mzones
          isepm1 = ledge
          zscale = 2.0 * used
        else
          isep   = nadump(1)
          isepm1 = isep - 1
          zscale = zscale * avi(mzones,12,1) / avi(isep,12,1)
        endif
c
        if (cfutz(impneu).gt.1.1.and.(natomc.ne.3)) call fdback
c
c
cl..confinement times
c
        zvols = 2.0 * avi(isep,12,1) * uisl**3   ! plasma volume
        zsurfi = avi(isep,3,1) * avi(isep,5,1)   ! plasma surface area
c
      zlosei = 0.0  ! eloctron power loss due to radiation and ionization
      zlosii = 0.0  ! ion power loss due to ionization and charge exchange
c

        do 718 jz = lcentr, isepm1
c
        zeirs = 0.0   ! electron power lost due to impurity radiation
        if (mimp.gt.0) then
        do 704 ji = 1, mimp
          zeirs = zeirs + weirs(ji,jz)
  704   continue
      endif
c
c   les  nov-90  add synchr rad
c
        zlosei = zlosei + usip*zvols*dx2i(jz)*
     1    (weions(jz) + webrs(jz) + zeirs + wesyn(jz) + wesrs(jz) )
        zlosii = zlosii - usip*zvols*dx2i(jz)*
     1                          (wiions(jz) + wichxs(jz))
c
  718   continue
c
      zflosi = 0.0   ! loss due to heat flux
c
      do 716 jp = lelec,lion
c
          zfluxi = 0.0
        do 714 jp2 = 1, mchi
          zfluxi= zfluxi- aaaa(jp,jp2,isep)*chi(jp2,isepm1)
     &                  - bbbb(jp,jp2,isep)*chi(jp2,isep)
  714   continue
c
          zflosi = zflosi + zfluxi * zsurfi * ueit
  716 continue
c
      zloss = (zlosei + zlosii + zflosi) * uisp
c
      if (zloss .gt. epslon)
     &  ztaue = (erges(isep) + ergis(isep)) / zloss
c
c..q-cylindrical
c
      zk = elong(isep,1)
      zqcyl = ( 5. * ahalfs(isep,1)**2 * bzs*useb * ( 1. + zk**2 ) )
     &         / ( 2. * rmids(isep,1) * curnts(isep) * usei )
c
c..confinement scalings
c
      zbeta = max (epslon, cgbeta(1) * betate(isep)
     &          + cgbeta(2) * betati(isep)
     &          + cgbeta(3) * betatb(isep)
     &          + cgbeta(4) * betata(isep) ) * 100.   ! vol ave beta %
c
      zpowr = max (epslon, cgpowr(1)*geohms(isep)
     &   + cgpowr(3)*gealfs(isep) + cgpowr(4)*gialfs(isep)
     &   + cgpowr(5)*geauxs(isep) + cgpowr(6)*giauxs(isep)
     &   + cgpowr(7)*gebems(isep) + cgpowr(8)*gibems(isep)
     &   + cgpowr(9)*geecrs(isep) + cgpowr(10)*giecrs(isep)
     &   + cgpowr(11)*geicrs(isep) + cgpowr(12)*giicrs(isep) )
c
      zpowr = zpowr * usep * 1.e-6     ! convert to MW
c
c  effective mass in AMU
c   les  nov-90  d-3he mass includes helium, protons
c
      zmeff = aspec(1)
      if (mhyd .eq. 2) zmeff = ( aspec(1) + aspec(2) ) / 2.
      if ( cfutz(490) .gt. epslon )
     & zmeff=(aspec(1) + aspec(2)+aspec(3)+aspec(4)+aspec(5))/5.
c
c  Neoalcator scaling [in seconds]
c
      ztoh = 7.e-22 * enlaes(isep) * ahalfs(isep,1)
     &  * (rmids(isep,1))**2 * zqcyl
c
      ztoh = min (ztoh,1./epslon)
c
      zftoh = ztaue / ztoh  ! ratio of confinement times
c
c..Goldston scaling [Plasma Phys and Cont. Fus. 26 (1984) 87]
c  with mass correction by Kaye, Aug 1988
c  as a function of heating power
c
       ztauxg = 0.0302 * (abs(curnts(isep))*usei*1.e-3)
     &  * (ahalfs(isep,1)*0.01)**(-0.37) * zk**0.5
     &  * (rmids(isep,1)*0.01)**1.75 * zpowr**(-0.5)
     &  * (abs(zmeff))**0.5
c
      ztauxg = max ( epslon, min (ztauxg,1./epslon) )
c
      zftaug = ztaue / ztauxg  ! ratio of taue / tau_Goldston
c
c..ITER89-P scaling
c
      ziter  = 0.048 * (abs(zmeff))**0.5
     &  * ( abs(curnts(isep))*usei*1.e-3)**0.85
     &  * (rmids(isep,1)*0.01)**1.2 * (ahalfs(isep,1)*0.01)**0.3
     &  * zk**0.5 * (enlaes(isep)*1.e-14)**0.1 * (bzs*1.e-4)**0.2
     &  * zpowr**(-0.5)
c
      ziter  = max ( epslon, min ( ziter, 1./epslon ) )
c
c  ratio of tau to ITER89-P scaling
c
      zfiter = ztaue / ziter
c
c  inverse squared combination
c
cbate      ztcomb = sqrt (1. / (1./ztoh**2 + 1./ztauxg**2) )
c
cbate      zfcomb = ztaue / ztcomb  ! ratio of taue / tau-Combination
c
c      write (6,*)'Cime', tai*uist,ztaue,zftaug,ztauxg,zfiter,ziter
c     &  ,zftoh,ztoh
c 131  format (' t=',0pf10.6,' s taue =',0pf9.5
c    &  ,' = (',0pf8.3,')*',0pf9.5,' GL, = (',0pf8.3,')*',0pf9.5
c     &  ,' ITER89-P,  (',0pf8.3,')*',0pf9.5,' NA')
c
c      write (6,*) 'Cime', tai*uist,enlaes(isep),zbeta,zqcyl
c 132  format (' t=',0pf10.6
c     &  ,' s ne-bar=',1pe11.3,' cm-3, beta = ',0pf8.3,'%'
c    &  ,' qcyl=',0pf6.3)
cpon
          do jk = 1, mzones
c PUB vtorr is the rotational velocity calculated from Neoclassical model 
c Pon implemented this formula (its ref is missing)
       vtorr(jk) = 4.19e+4 * vloopi(jk,2) * zcharge(2) * ztikev(jk)**1.5/(2*zdensi(jk)*zrmajor(mzones))
       zvrotNEO(jk) = vtorr(jk)
c           
c PUB boyle is the rotational velocity from the torque models
c in which the radial current is needed in the calculation.
c However BALDUR does not have the radial current.
c The existing bolye is incorrect (need further investigation)     
       boyle(jk) = ztauxg * 1.e+25 * 1.e-6 * ( cjbeam * ajbs(jk) * usij
     &            + cdprof(jk) ) * bpoli(jk) / (2 * zdensi(jk))               
c       zvrottorque(jk) = boyle(jk)
       zvrottorque(jk) = 0
c       vtorr(jk) = zvrotNEO(jk) + boyle(jk) 
          end do

c        write(6,*) 'ion telocity', boyle(mj)     
            do jz = 1, mzones
               alpha_NTV = zdensimp(jz)*zcharge(3)*zcharge(3)/
     &         (zdensi(jz)*zcharge(2)*zcharge(2))
               mu1_NTV = g_NTV*(alpha_NTV+sqrt(2.)-log(1.+sqrt(2.)))
               mu2_NTV = g_NTV*(3.*alpha_NTV/2. + 4./sqrt(2.) - 
     &              5.*log(1.+sqrt(2.))/2.)
               mu3_NTV = g_NTV*(13.*alpha_NTV/4. + 39./(4.*sqrt(2.))
     &              - 25.*log(1.+sqrt(2.))/4.)
               D_NTV = mu1_NTV*(mu3_NTV+sqrt(2.)+alpha_NTV)-
     &              mu2_NTV*mu2_NTV
               K1_NTV = (sqrt(2.)+alpha_NTV)*mu2_NTV/D_NTV
c PUB the following 5 lines has been modified (13/03/15)
               zvrotNTV(jz) = (3.54-K1_NTV)*
     &              zgrdti(jz)*ztikev(jz)*1000./(zrmajor(jz)*zcharge(2)*
     &              bpoli(jz))
               zvrotxb(jz) = zvrottorque(jz) + zvrotNTV(jz) 
            end do
c            write(nprint,*) 'Predicted wexb  by Neoclassical 
c     &              Rotation Model'
c
         end if
c
c---------------------------------------------------------------------
c poloidal velocity piece added by pzhu (8/99)
c---------------------------------------------------------------------
!ap{
         select case (lthery(46))
         case(1)
            imod = 1            ! Strand model 
         case(2)
            imod = 2            ! Erba model
         case default
            imod = 3            ! P.Zhu model
         end select
!}
c
c           write(6,*)'imod = ',imod
c
         ispecies = 1 + mhyd + mimp
c
c        electron
c
ctontontonton     
         zmasspecies(1) = fcae*(10.0**fxae)
c         write(nprint,'(A,E8.3)') 'kekkekkek is', zmasspecies(1)
         zcharge(1) = -1.e0
         zdensity(1, 1:mzones) = zdense(1:mzones)
         ztemperature(1, 1:mzones) = ztekev(1:mzones)*1.e3
c
c        ions
c
         zmasspecies(2:ispecies) = aspec(1:mhyd+mimp)
         zcharge(2:1+mhyd) = 1.e0
         zcharge(2+mhyd:ispecies) = cmean(1:mimp, 1, mzones)
                                        !using edge charge
         zdensity(2:1+mhyd, 1:mzones) =
     &                  rhohs(1:mhyd, 1, 1:mzones)*usid
         zdensity(2+mhyd:ispecies, 1:mzones) =
     &                  rhois(1:mimp, 1, 1:mzones)*usid
         do is = 2, ispecies
         ztemperature(is, 1:mzones) = ztikev(1:mzones)*1.e3
         end do
c
         zvpolxb = 0.0
!ap
         call vpolprof (imod, ispecies, maxs, zmasspecies, zcharge,
     &          mzones, zdensity(1:maxs,1:mzones),
     &          ztemperature(1:maxs, 1:mzones),
     &          zvpolxb(1:maxs,1:mzones),
     &          zrminor, zrmajor, zbtor, bpoli, q, ierr)
c
         do jz=1,mzones
           zvpol(jz) = zvpolxb(2,jz)
         enddo
c
c---------------------------------------------------------------------
c <wexbprof> modified by pzhu (9/9/99)
c---------------------------------------------------------------------
         call wexbprof (mzones, xbouni, flpols, zdensi, ztikev
     &     , bpoli, zbtor, zvpol, zvrotxb
     &     , zrminor, zrmajor, xbouni, zwexbxb, knthe, nprint)
c---------------------------------------------------------------------
      end if
c
      iprint = 0
      if ( knthe .gt. 2 ) iprint = - max ( 1, lthery(29) )
c
c..call the transport model
c
      if ( lthery(21) .lt. 1 ) then
c
        if ( nstep .lt. 2 ) then
          write (nprint,*)
          write (nprint,* ) 'Srtn theory called from sbrtn ptheory'
          write (nprint,*)
        endif
c
        call theory( lthery, cthery
     & , iaxis, iedge, isep, imatdim
     & , zrminor, zrmajor, zelong, ztriang
     & , zdense, zdensi, zdensh, zdensimp, zdensf, zdensfe
     & , zxzeff, ztekev, ztikev, ztfkev, q, zvloop, zbtor, zresist
     & , zavezimp, zmassimp, zmasshyd, zaimass, zwexbxb
     & , zgrdne, zgrdni, zgrdnh, zgrdnz, zgrdte, zgrdti, zgrdpr, zgrdq
     & , fdr, fig, fti, frm, fkb, frb, fhf, fec, fmh, fdrint
     & , zdhthe, zvhthe, zdzthe, zvzthe, zxethe, zxithe, zweithe
     & , difthi, velthi
     & , nstep, ztime, nprint, iprint)
c
c..Some pages of long printout
c
        if ( knthe .eq. 3  .and.  lthery(29) .gt. 2 ) then
c
          zt  = tai * uist * 1000.0
          zdt = dtoldi * uist * 1000.0
c
c..Total diffusivities if more printout is desired
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c..total theory-based diffusivities
c
      write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
       if (lthery(39) >= 0) then
        write(nprint,10007)
        write(nprint,10009)
        write(nprint,10008)
c
       do jz=iaxis,mzones
         write(nprint,101) jz, zrminor(jz), zxethe(jz)
     &      , zxithe(jz), zweithe(jz), zdhthe(jz), zdzthe(jz)
       enddo
      endif ! lthery(39) >= 0
c
c.. total power interchanged
c  rgb 15-apr-95 appended factor elong(jz) to volume expression
c
      zptot=0.0
c
      do jz=iaxis,mzones
        zdvol = 2.0 * fcpi * zrmajor(jz) * fcpi
     &     * (zrminor(jz)**2 - zrminor((jz-1))**2) * zelong(jz)
c
        zptot = zptot + weithe(jz)*zdvol
      enddo
c
      write(nprint,102) zptot
c
c..totals, including neoclassical and empirical effective diffusivities
c
c..electron thermal diffusivities
c
      lpage=lpage+1
      write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write(nprint,104)
      write(nprint,105)
c
      do jz=iaxis,mzones
c
        zxeneoeff = ( cfutz(11) * xeneo1(jz) + xeneo2(jz) )*usil**2
     &    + cfutz(11) * veltes(jz) * usil * zrmajor(jz)
     &        / sign( max ( abs ( zgrdte(jz)), epslon ), zgrdte(jz) )
c
        zxetot(jz) = zxethe(jz)
     &    + zxeneoeff
     &    + xeemps(jz) * usil**2
c
        write(nprint,106) jz, zrminor(jz), zxethe(jz)
     &  , zxeneoeff
     &  , xeemps(jz)*usil**2
     &  , zxetot(jz)
      enddo
c
c
c..ion thermal diffusivities
c
      lpage=lpage+1
      write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write(nprint,107)
      write(nprint,108)
c
      do jz=iaxis,mzones
c
        zxineoeff = ( cfutz(12) * xineo1(jz) )*usil**2
     &    + cfutz(12) * veltis(jz) * usil * zrmajor(jz)
     &        / sign( max ( abs ( zgrdte(jz)), epslon ), zgrdte(jz) )
c
        zxitot(jz) = zxithe(jz)
     &    + zxineoeff
     &    + xiemps(jz) * usil**2
c
        write(nprint,106) jz, zrminor(jz), zxithe(jz)
     &    , zxineoeff
     &    , xiemps(jz)*usil**2
     &    , zxitot(jz)
      enddo
c
        endif
c
c..end of printout after sbrtn theory, beginning of other models
c
      else
c
c..preprocess gradients for the Multi-Mode model
c
c  zsgrd* are the normalized gradients for preprocessing
c
        do jz=1,mzones
          zsgrdne(jz) = zgrdne(jz)
          zsgrdni(jz) = zgrdni(jz)
          zsgrdnh(jz) = zgrdnh(jz)
          zsgrdnz(jz) = zgrdnz(jz)
          zsgrdte(jz) = zgrdte(jz)
          zsgrdti(jz) = zgrdti(jz)
          zsgrdpr(jz) = zgrdpr(jz)
          zsgrdq(jz)  = zgrdq(jz)
c--Define a local copy of normalized ExB shearing rate
c          write(6,*) 'jz = ', jz, 'zwexbxb before = ', zwexbxb(jz)
          zwexbxb(jz) = cthery(129) * zwexbxb(jz)
c          write(6,*) 'jz = ', jz, 'zwexbxb after = ', zwexbxb(jz)
        enddo
c
c  use smoothing when lthery(32) .ne. 0
c
        if ( lthery(32) .ne. 0 ) then
c
          i1     = iaxis + 1
          i2     = iedge - 1
          ismord = abs( lthery(32) )
          zlmin  = 1.e-4
          zmix   = 0.0
c
c  ismord = number of times smoothing is applied
c  zlmin  = minimum gradient scale length
c  zmix   = fraction of original array added to smoothed array
c
c  if smoothing is to be done, first divide by the minor radius
c
          do jz=1,mzones
            zsgrdne(jz) = zsgrdne(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdni(jz) = zsgrdni(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdnh(jz) = zsgrdnh(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdnz(jz) = zsgrdnz(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdte(jz) = zsgrdte(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdti(jz) = zsgrdti(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdpr(jz) = zsgrdpr(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
            zsgrdq(jz)  = zsgrdq(jz)
     &        / sign ( max( abs( zrminor(jz) ), zepsqrt ), zrminor(jz) )
          enddo
c
c  change the values at jz=maxis+1 in order to make sure zml** arrays
c    are monotonic near the magnetic axis
c    Note that this process used to be applied to zml*
c
          if ( lthery(31) .eq. 1 ) then
            zsgrdne(i1) = 2.0 * zsgrdne(i1+1) - zsgrdne(i1+2)
            zsgrdni(i1) = 2.0 * zsgrdni(i1+1) - zsgrdni(i1+2)
            zsgrdnh(i1) = 2.0 * zsgrdnh(i1+1) - zsgrdnh(i1+2)
            zsgrdnz(i1) = 2.0 * zsgrdnz(i1+1) - zsgrdnz(i1+2)
            zsgrdte(i1) = 2.0 * zsgrdte(i1+1) - zsgrdte(i1+2)
            zsgrdti(i1) = 2.0 * zsgrdti(i1+1) - zsgrdti(i1+2)
            zsgrdpr(i1) = 2.0 * zsgrdpr(i1+1) - zsgrdpr(i1+2)
            zsgrdq(i1)  = 2.0 * zsgrdq(i1+1)  - zsgrdq(i1+2)
          endif
c
c  smoothing
c
          i1 = iaxis + 1
          i2 = iedge
c
          call smooth2 ( zsgrdne, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdni, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdnh, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdnz, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdte, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdti, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdpr, 1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
          call smooth2 ( zsgrdq,  1, ztemp1, ztemp2, 1, i1, i2
     &      , ismord, zmix )
c
c  undoing the preconditioning
c
          do jz=1,mzones
            zsgrdne(jz) = zsgrdne(jz) * zrminor(jz)
            zsgrdni(jz) = zsgrdni(jz) * zrminor(jz)
            zsgrdnh(jz) = zsgrdnh(jz) * zrminor(jz)
            zsgrdnz(jz) = zsgrdnz(jz) * zrminor(jz)
            zsgrdte(jz) = zsgrdte(jz) * zrminor(jz)
            zsgrdti(jz) = zsgrdti(jz) * zrminor(jz)
            zsgrdpr(jz) = zsgrdpr(jz) * zrminor(jz)
            zsgrdq(jz)  = zsgrdq(jz) * zrminor(jz)
          enddo
c
c  end of smoothing when lthery(32) .ne. 0
c
        endif
c
c..minimum gradient
c
         if ( cthery(50) .gt. 0.5 ) then
           do jz=1,mzones
             zsgrdne(jz) = sign ( max ( abs ( zsgrdne(jz) ),
     &         1.0 / cthery(50) ), zsgrdne(jz) )
           enddo
         endif
c
         if ( cthery(51) .gt. 0.5 ) then
           do jz=1,mzones
             zsgrdni(jz) = sign ( max ( abs ( zsgrdni(jz) ),
     &         1.0 / cthery(51) ), zsgrdni(jz) )
             zsgrdnh(jz) = sign ( max ( abs ( zsgrdnh(jz) ),
     &         1.0 / cthery(51) ), zsgrdnh(jz) )
             zsgrdnz(jz) = sign ( max ( abs ( zsgrdnz(jz) ),
     &         1.0 / cthery(51) ), zsgrdnz(jz) )
           enddo
         endif
c
         if ( cthery(52) .gt. 0.5 ) then
           do jz=1,mzones
             zsgrdte(jz) = sign ( max ( abs ( zsgrdte(jz) ),
     &         1.0 / cthery(52) ), zsgrdte(jz) )
           enddo
         endif
c
         if ( cthery(53) .gt. 0.5 ) then
           do jz=1,mzones
             zsgrdti(jz) = sign ( max ( abs ( zsgrdti(jz) ),
     &         1.0 / cthery(53) ), zsgrdti(jz) )
           enddo
         endif
c
         if ( cthery(54) .gt. 0.5 ) then
           do jz=1,mzones
             zsgrdpr(jz) = sign ( max ( abs ( zsgrdpr(jz) ),
     &         1.0 / cthery(54) ), zsgrdpr(jz) )
           enddo
         endif
c
c%%%%%%%%%
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..call sbrtn mmm95
c
        if ( lthery(21) .eq. 1 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
        if ( nstep .lt. 2 ) then
          write (nprint,*)
          write (nprint,* ) 'Srtn mmm95 called from sbrtn ptheory'
          write (nprint,*)
        endif
c	
          call mmm95 (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c..call sbrtn mmm98
c
        else if ( lthery(21) .eq. 2 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* ) 'Srtn mmm98 called from sbrtn ptheory'
            write (nprint,*)
          endif
c
          call mmm98 (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c..call sbrtn mmm98b
c
        else if ( lthery(21) .eq. 3 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* ) 'Srtn mmm98b called from sbrtn ptheory'
            write (nprint,*)
          endif
c
          call mmm98b (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c..call sbrtn mmm98c
c
        else if ( lthery(21) .eq. 4 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* ) 'Srtn mmm98c called from sbrtn ptheory'
            write (nprint,*)
          endif
c
          call mmm98c (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c..call sbrtn mmm98d
c
        else if ( lthery(21) .eq. 5 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* ) 'Srtn mmm98d called from sbrtn ptheory'
            write (nprint,*)
          endif
c
          call mmm98d (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c
c..call sbrtn ohe
c
        else if ( lthery(21) .eq. 6 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'Srtn ohe_model called from sbrtn ptheory'
            write (nprint,*)
          endif
c
      call ohe_model (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1)
     & , ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1), zmasshyd(i1)
     & , zsgrdnz(i1),  zsgrdte(i1),   zsgrdq(i1),   zsgrdnh(i1)
     & , zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1)
     & , zsgrdti(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c
c..call Matteo Erba's version of the Mixed Bohm/gyro-Bohm model
c
        else if ( lthery(21) .eq. 7 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'Srtn mixed_merba called from sbrtn ptheory'
            write (nprint,*)
          endif
c
      call mixed_merba (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1)
     & , ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1), zmasshyd(i1)
     & , zsgrdnz(i1),  zsgrdte(i1),   zsgrdq(i1),   zsgrdnh(i1)
     & , zwexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1)
     & , zsgrdti(i1)
     & , zthigb(i1),   zthegb(i1),    zthibohm(i1), zthebohm(i1)
     & , zthimix(i1),   zthemix(i1),  zthdmix(i1),  zthzmix(i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb)
c
c
c..call Onjun's version of the Mixed Bohm/gyro-Bohm model
c
        else if ( lthery(21) .eq. 8 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 0
c
          lflowshear = lthery(37)
c
          zcharge_hyd(:) = 1.0
c
          zt_e_kev_edge = ztekev(iedge)
          zrminor_edge  = zrminor(iedge)
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'Srtn mixed_model called from sbrtn ptheory'
            write (nprint,*)
          endif
c
      call mixed_model (
     &   zrminor(i1),  zrmajor(i1)
     & , ztekev(i1),   ztikev(i1),    q(i1)
     & , zbtor(i1),    zmasshyd(i1),  zcharge_hyd(i1),  zwexbxb(i1)
     & , zsgrdte(i1),  zsgrdne(i1),   zsshr(i1)
     & , zt_e_kev_edge, zrminor_edge
     & , ipoints
     & , zthimix(i1),   zthemix(i1),  zthdmix(i1)
     & , zthigb(i1),   zthegb(i1),    zthibohm(i1), zthebohm(i1)
     & , ierr
     & , lflowshear=lflowshear)
c
c      call mixed_model (
c     &   zrminor(i1),  zrmajor(i1),   zelong(i1)
c     & , ztekev(i1),    ztikev(i1),    q(i1)
c     & , zbtor(i1), zmasshyd(i1)
c     & , zsgrdnz(i1),  zsgrdte(i1),   zsgrdq(i1),   zsgrdnh(i1)
c     & , zwexbxb(i1)
c     & , zsgrdne(i1),  zsgrdni(i1)
c     & , zdense(i1),   zdensh(i1),    zdensimp(i1)
c     & , zsgrdti(i1)
c     & , zthigb(i1),   zthegb(i1),    zthibohm(i1), zthebohm(i1)
c     & , zthimix(i1),   zthemix(i1),  zthdmix(i1),  zthzmix(i1)
c     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
c     & , imatdim,  ipoints,   nprint,    iprint,  ierr
c     & , lsuper,   lreset,    lmmm95,    cmmm95
c     & , fig,      frb,       fkb)
c
c Add by Boonyarit Chatthong 03/23/10
c      write(nprint,*) 'Tontime'
c      write (nprint,151) nstep, ztime
c      write(nprint,'(0pf14.6)') ztime
c      write (nprint,139)
c        do j=1,iedge
c           write(nprint,'(1p10e12.4)') zrminor(j)
c           write(nprint,'(1p10e12.4)') zsgrdti(j)
c          write (nprint,152) zrminor(j), zsgrdne(j), zsgrdni(j)
c     &      , zsgrdnh(j), zsgrdnz(j), zsgrdte(j), zsgrdti(j)
c     &      , zsgrdpr(j), zsgrdq(j)
c        enddo
c Finish Add by Boonyarit
          if ( lthery(45) .eq. 1 ) then
            do jz=1,mzones
c
               if (cthery(139) .gt. 1.0E-10) then 
                  cont_imp = cthery(139)
               else
                  cont_imp = 1.0
               endif
c               
               zthzmix(jz) = cont_imp*zthdmix(jz)
c           
            enddo
c
          elseif ( lthery(45) .eq. 2 ) then
            if (cthery(138) .gt. 1.0E-10) then 
               cont_imp = cthery(138)
            else
              cont_imp = 1.25
            endif
c
c            write(6,*)'cont_imp = ', cont_imp
c
            do jz=1,mzones
              zthzmix(jz) = abs( cont_imp*(1
     &          + zrminor(jz)**2 / zrminor(mzones)**2 ))
            enddo
          endif
c
c..pass thermal diffusivities to difthi matrix
c
          difthi(1,1,:) = zthimix(:)
          difthi(2,2,:) = zthdmix(:)
          difthi(3,3,:) = zthemix(:)
          difthi(4,4,:) = zthzmix(:)
c
c..call routine mmm99
c
        else if (lthery(21) .eq. 9) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 1
          lmmm95(1) = lthery(7)  ! Weiland ITG model weiland14 (10 eqns, no collisions)
cap       lmmm95(2) = lthery(8)  ! use effective diffusivities
          lmmm95(3) = lthery(12) ! use kappa instead of (1+\kappa^2)/2
          lmmm95(4) = lthery(27) ! replace -ve diffusivity with convective velocity
	  lmmm95(5) = 1          ! limit gradients by major radius / ion Larmor radius
c
	  cmmm95(1) = cthery(3)    ! minimum value of shear
          cmmm95(2) = cthery(8)    ! for fbeta-th in kinetic ballooning
          cmmm95(3) = cthery(12)   ! elongation scaling for Weiland model
          cmmm95(4) = cthery(14)   ! elongation scaling for RB model
          cmmm95(5) = cthery(15)   ! elongation scaling for KB mode
          cmmm95(6) = cthery(38)   ! k_y \rho_s (= 0.316 if abs(cthery(6)) < zepslon)
          cmmm95(8) = cthery(78)   ! coeff of beta_prime_1 in kinetic ballooning mode
          cmmm95(9) = cthery(86)   ! Diamagnetic stabilization in Guzdar-Drake model
          cmmm95(10) = cthery(111) ! difthi -> velthi for chi_i
          cmmm95(11) = cthery(112) ! difthi -> velthi for hydrogen
          cmmm95(12) = cthery(113) ! difthi -> velthi for chi_e
          cmmm95(13) = cthery(114) ! difthi -> velthi for impurity
          cmmm95(14) = cthery(119) ! coeff of finite beta in weiland14 = cetain(20)
          cmmm95(15) = cthery(120) ! min value of impurity charge state zimpz
          cmmm95(16) = cthery(121) ! set fast particle fraction for use in weiland14
          cmmm95(17) = cthery(123) ! coeff of k_\parallel in weiland14 = cetain(10)
          cmmm95(18) = cthery(124) ! coeff of nuhat in weiland14 = cetain(15)
          cmmm95(19) = cthery(125) ! 0.0 -> 1.0 for v_parallel in strong balloon limit = cetain(12)
          cmmm95(20) = cthery(126) ! trapping fraction used in weiland14 (when > 0.0)
                                    ! multiplies electron trapping fraction when < 0.0
          cmmm95(21) = cthery(129) ! multiplier for wexbs
          cmmm95(22) = cthery(130) ! multiplier to impurity heat flux
          cmmm95(23) = cthery(131) ! controls finite diff to construct the zgm matrix = cetain(30)
          cmmm95(27) = cthery(135) ! exponent of Hamaguchi-Horton parameter in MMM99
          cmmm95(28) = cthery(136) ! coefficient of Hamaguchi-Horton parameter in MMM99
          cmmm95(24) = 0.043        ! ???? fit coefficient of electrostatic electron thermal dif.
          cmmm95(25) = 1.88         ! ???? coefficient of critical gradient in e-st electron th diff
          cmmm95(26) = 0.082        ! ???? fit oefficient of electromagnetic electron thermal
          cmmm95(29) = 0.5         ! minimum magnetic shear used in the
                                   ! kinetic ballooning mode
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'Srtn mmm99 called from sbrtn ptheory'
            write (nprint,*)
          endif
c
        call mmm99 (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1),   zyexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & ,                              ztheeg(i1)
     & ,                              zthetb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  knthe, ierr
     & , lsuper,   lreset,    lmmm95,    cmmm95
     & , fig,      frb,       fkb,       feg
     & , lthery(34), lthery(35), lthery(36),    lthery(37)
     & , lthery(38), lthery(39), lthery(40))
c
           if( nerr .gt. 0) then
               write(nprint,*) 'Error after sbrtn mmm99 ifail = ',nerr
           endif
c
c..call routine mmm2001
c
        else if ( lthery(21) .eq. 10 ) then
c
c..set switches
c
          i1 = 3
          ipoints = iedge + 1 - i1
c
          lsuper = lthery(22)
          lreset = 1
c
c..set metric elements
c
        do jz=1,mzones
          zgradrsqrave(jz) = 1.0
          zgradrave(jz)    = 1.0
          zgradroutbrd(jz) = 1.0
        enddo
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* ) 'Srtn mmm2001 called from sbrtn ptheory'
            write (nprint,*)
          endif
c
        call mmm2001 (
     &   zrminor(i1),  zrmajor(i1),   zelong(i1),    ztriang(i1)
     & , zdense(i1),   zdensh(i1),    zdensimp(i1),  zdensfe(i1)
     & , zxzeff(i1),   ztekev(i1),    ztikev(i1),    q(i1)
     & , zbtor(i1),    zavezimp(i1),  zmassimp(i1),  zmasshyd(i1)
     & , zaimass(i1),  zwexbxb(i1),   zyexbxb(i1)
     & , zsgrdne(i1),  zsgrdni(i1),   zsgrdnh(i1),   zsgrdnz(i1)
     & , zsgrdte(i1),  zsgrdti(i1),   zsgrdq(i1)
     & , zgradrsqrave(i1), zgradrave(i1),  zgradroutbrd(i1)
     & , zthiig(i1),   zthdig(i1),    ztheig(i1),    zthzig(i1)
     & , zthirb(i1),   zthdrb(i1),    ztherb(i1),    zthzrb(i1)
     & , zthikb(i1),   zthdkb(i1),    zthekb(i1),    zthzkb(i1)
     & ,                              ztheeg(i1)
     & ,                              zthetb(i1)
     & , zthitem(i1),  zthdtem(i1),   zthetem(i1),   zthztem(i1)
     & , zthiitg(i1),  zthditg(i1),   ztheitg(i1),   zthzitg(i1)
     & , zgamma(1,i1), zomega(1,i1)
     & , difthi(1,1,i1),  velthi(1,i1), zvflux(1,i1)
     & , imatdim,  ipoints,   nprint,    iprint,  knthe, ierr
     & , lsuper,   lreset,    lthery,    cthery
     & , fig,      frb,       fkb,       feg )
c
           if( nerr .gt. 0) then
               write(nprint,*) 'Error after sbrtn mmm2001 nerr = '
     &           ,nerr
           endif
c
c
c..call routine callglf2d for the GLF23 model
c
        else if ( lthery(21) .eq. 23 ) then
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'GLF23 model called using callglf2db from sbrtn ptheory'
            write (nprint,*)
          endif
c
          i1 = 2
          ipoints = iedge - i1 + 1
c
c..compute d r / d rho at zone boundaries
c
          do jz=3,mzones
            zdrdrho(jz) = ( armins(jz,2) - armins(jz-1,2) ) * usil
     &                    / ( xzoni(jz) - xzoni(jz-1) )
          enddo
          zdrdrho(2) = zdrdrho(3)
          zdrdrho(1) = zdrdrho(3)
          zdrdrho(mzones+1) = zdrdrho(mzones)
c
c..transform zsgrdne to be - d n_e / d rho / n_e
c
c  zdqdrho(jz) = rho d q / d rho
c  zalpha(jz)  = 2 mu_0 R q^2 ( d p / d r ) / B^2
c
          do jz=1,mzones
            zdqdrho(jz) = zrminor(jz) * zsgrdq(jz) / zrmajor(jz)
            zalpha(jz)  = zsgrdpr(jz) * q(jz)**2 * betatt(jz)
          enddo
          zdqdrho(mzones+1) = zdqdrho(mzones)
          zalpha(mzones+1)  = zalpha(mzones)
c
c..set metric elements
c
c  Note:  densities must be in units of 10^{19} m^{-3}
c
c  The arrays zsgrdte ... zsgrdni are converted to logaritmic derivatives
c  with respect to rho
c
        do jz=1,mzones
          zgradrsqrave(jz) = avi(jz,6,1)
          zgradrave(jz)    = avi(jz,5,1)
          zgradroutbrd(jz) = 1.0
          ztemp1(jz)       = zdense(jz) * 1.0e-19
          ztemp2(jz)       = zdensi(jz) * 1.0e-19
          ztemp3(jz)       = zdensfe(jz) * 1.0e-19
c
          zsgrdte(jz) = zdrdrho(jz) * zsgrdte(jz) / zrmajor(jz)
          zsgrdti(jz) = zdrdrho(jz) * zsgrdti(jz) / zrmajor(jz)
          zsgrdne(jz) = zdrdrho(jz) * zsgrdne(jz) / zrmajor(jz)
          zsgrdni(jz) = zdrdrho(jz) * zsgrdni(jz) / zrmajor(jz)
        enddo
c
        zgradrsqrave(mzones+1) = zgradrsqrave(mzones)
        zgradrave(mzones+1)    = zgradrave(mzones)
        zgradroutbrd(mzones+1) = zgradroutbrd(mzones)
c
        ztekev(mzones+1)    = ztekev(mzones)
        ztikev(mzones+1)    = ztikev(mzones)
        ztemp1(mzones+1)    = ztemp1(mzones)
        ztemp2(mzones+1)    = ztemp2(mzones)
        ztemp3(mzones+1)    = ztemp3(mzones)
c
        zsgrdte(mzones+1)   = zsgrdte(mzones)
        zsgrdti(mzones+1)   = zsgrdti(mzones)
        zsgrdne(mzones+1)   = zsgrdne(mzones)
        zsgrdni(mzones+1)   = zsgrdni(mzones)
c
        zxzeff(mzones+1)    = zxzeff(mzones)
        zrminor(mzones+1)   = 2 * zrminor(mzones) - zrminor(mzones-1)
        zrmajor(mzones+1)   = zrmajor(mzones)
c
c..initialize
c
c
        allocate ( zexb_out2(1:ipoints+i1,10) )
c
        zexb_out1  = 0.
        zexb_out2  = 0.
        zexb_out3  = 0.
c
c..call the GLF23 model one zone boundary at a time
c
        do jz=i1,mzones
c
         imm = 0  ! compute full grid all at once
c
c..set switches
c
          ieigen = 1   ! use the tomsqz eigenvalue solver
          ishoot = 1   ! time-dependent code
c          igrid  = 0  ! do full grid of ipoints grid points
          i_grad = 2   ! compute gradients externally as an array
          idengrad = 2 ! for simple dilution
c
          iswitch(1) = 0  ! ion particle transport
          iswitch(2) = 1  ! electron thermal transport
          iswitch(3) = 1  ! ion thermal transport
          iswitch(4) = 0
          iswitch(5) = 0
c
c..scalars
c
          zbtor_axis = zbtor(iaxis)
                       ! toroidal field at magnetic axis (tesla)
          ztor_flux  = fltors(mzones,1) * usib * usil**2
                       ! toroidal flux at the edge of the plasma( tesla m^2)
          zrmaj_axis = zrmajor(iaxis)
                       ! major radius at magnetic axis
          zsw_alpha  = 0
          if ( cthery(129) .gt. epslon ) zsw_alpha  = 1.
                       ! include effect of flow shear stabilization
          zsw_x_alpha = 0
          if ( cthery(119).gt. epslon ) zsw_x_alpha = 1.
                       ! include alpha effec (see zalpha above)
          i_delay    = 0
                       ! time delay for flow shear
c
          zscalgrdte = 0.
          zscalgrdti = 0.
          zscalgrdne = 0.
          zscalgrdni = 0.
c
          zavezimps  = zavezimp(i1)
          zmassimps  = zmassimp(i1)
          zmasshyds  = zmasshyd(i1)
c
c..initialize
c
          zstemp1 = 0.
          zstemp2 = 0.
          zstemp3 = 0.
          zstemp4 = 0.
c
          zvtor   = 0.
          zvpara  = 0.
          zvperp  = 0.
c
          zschie  = 0.
          zschii  = 0.
          zsvtor  = 0.
          zsvpara = 0.
          zsvperp = 0.
          zsexchange = 0.
c
c..For the moment, turn off flow shear and alpha stabilization
c
          zsw_alpha  = cthery(131)
          zsw_x_alpha= cthery(132)
c
          iroot = 8
          irotstab = 1  ! to use egamma_exp
          ibt_flag = 0
c
          iglf  = 1     ! use retuned version of GLF23 model
          if ( lthery(34) .lt. 0 ) then
            iglf = 0    ! for the original GLF23 model
          endif
c
c..call the GLF23 model
c
c  Use zthikb, ... as dummy output arrays which are zeroed out later
c
          call callglf2db ( ieigen, iroot, iglf
     & , ishoot, imm, ipoints, iswitch, irotstab
     & , ztekev(i1),   ztikev(i1),   ztemp1(i1),   ztemp2(i1)
     & , ztemp3,   i_grad,   idengrad
     & , zsgrdte(i1),  zsgrdti(i1),  zsgrdne(i1),  zsgrdni(i1)
     & , zscalgrdte,   zscalgrdti,   zscalgrdne,   zscalgrdni
     & , zvrotxb(i1),  zwexbxb(i1),  zyexbxb(i1),  zvtor(i1)
     & , zvpara(i1),   zvperp(i1),   zxzeff(i1),   zbtor_axis
     & , ibt_flag
     & , xbouni(i1),   ztor_flux,    zgradrave(i1), zgradrsqrave(i1)
     & , zrminor(i1),  zrmajor(i1),  zrmaj_axis
     & , zavezimps, zmassimps, q(i1),        zdqdrho(i1)
     & , zalpha(i1),   zelong(i1),   zmasshyds
     & , zsw_alpha,    zsw_x_alpha,  i_delay,      zsdifi
     & , zschie, zschii, zsvtor, zsvpara, zsvperp, zsexchange
     & , zthdkb(i1),   zthekb(i1),   zthikb(i1)
     & , zthvtor(i1),  zthvpara(i1), zthvperp(i1), zweithe(i1)
     & , zexb_out1(i1), zexb_out2(i1,1), zexb_out3(i1)
     & , zgammaitg(i1), zomegaitg(i1)
     & , zgammatem(i1), zomegatem(i1) )
c
c..set output arrays
c
c..set impurity transport equal to ion particle transport
c
c..move negative diffusivities to convective velocities
c
c..use fig(1) to fig(4) for user control of coefficients
c
c        do jz=i1,mzones
c
          zschii = zthikb(jz)
          if ( zschii .gt. 0. ) then
            difthi(1,1,jz) = difthi(1,1,jz) + fig(3) * zschii
          else
            velthi(1,jz) = velthi(1,jz)
     &         - fig(3) * zschii * zsgrdti(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          zthiig(jz) = max ( fig(3) * zschii, 0.0 )
          zthikb(jz) = 0.
c
          zsdifi = zthdkb(jz)
          if ( zsdifi .gt. 0. ) then
            difthi(2,2,jz) = difthi(2,2,jz) + fig(1) * zsdifi
            difthi(4,4,jz) = difthi(4,4,jz) + fig(4) * zsdifi
          else
            velthi(2,jz) = velthi(2,jz)
     &         - fig(1) * zsdifi * zsgrdnh(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
             velthi(4,jz) = velthi(4,jz)
     &         - fig(4) * zsdifi * zsgrdnz(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          zthdig(jz) = max ( fig(1) * zsdifi, 0.0 )
          zthzig(jz) = max ( fig(4) * zsdifi, 0.0 )
          zthdkb(jz) = 0.
c
          zschie = zthekb(jz)
          if ( zschie .gt. 0. ) then
            difthi(3,3,jz) = difthi(3,3,jz) + fig(2) * zschie
          else
            velthi(3,jz) = velthi(3,jz)
     &         - fig(2) * zschie * zsgrdte(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          ztheig(jz) = max ( fig(2) * zschie, 0.0 )
          zthekb(jz) = 0.
c
c
         enddo
c
         deallocate ( zexb_out2 )
c
c..zero out the zthikb, ... arrays that were used as dummy output
c
         zthikb = 0.0
         zthdkb = 0.0
         zthekb = 0.0
         zthzkb = 0.0
c
c
c..call routine callglf2d for the GLF23 model
c
        else if ( lthery(21) .eq. 24 ) then
c
          if ( nstep .lt. 2 ) then
            write (nprint,*)
            write (nprint,* )
     &        'GLF23 model called using callglf2d from sbrtn ptheory'
            write (nprint,*)
          endif
c
          i1 = 3
          ipoints = iedge - i1 + 1
c
c..compute d r / d rho at zone boundaries
c
          do jz=3,mzones
            zdrdrho(jz) = ( armins(jz,2) - armins(jz-1,2) ) * usil
     &                    / ( xzoni(jz) - xzoni(jz-1) )
          enddo
          zdrdrho(2) = zdrdrho(3)
          zdrdrho(1) = zdrdrho(3)
          zdrdrho(mzones+1) = zdrdrho(mzones)
c
c..transform zsgrdne to be - d n_e / d rho / n_e
c
c  zdqdrho(jz) = rho d q / d rho
c  zalpha(jz)  = 2 mu_0 R q^2 ( d p / d r ) / B^2
c
          do jz=1,mzones
            zdqdrho(jz) = xbouni(jz) * q(jz) * zdrdrho(jz) * zsgrdq(jz)
     &                    / zrmajor(jz)
            zalpha(jz)  = zsgrdpr(jz) * q(jz)**2 * betatt(jz)
          enddo
          zdqdrho(mzones+1) = zdqdrho(mzones)
          zalpha(mzones+1)  = zalpha(mzones)
c
c..set metric elements
c
c  Note:  densities must be in units of 10^{19} m^{-3}
c
        do jz=1,mzones
          zgradrsqrave(jz) = avi(jz,6,1)
          zgradrave(jz)    = avi(jz,5,1)
          zgradroutbrd(jz) = 1.0
          ztemp1(jz)       = zdense(jz) * 1.0e-19
          ztemp2(jz)       = zdensi(jz) * 1.0e-19
          ztemp3(jz)       = zdensfe(jz) * 1.0e-19
        enddo
        zgradrsqrave(mzones+1) = zgradrsqrave(mzones)
        zgradrave(mzones+1)    = zgradrave(mzones)
        zgradroutbrd(mzones+1) = zgradroutbrd(mzones)
c
        ztekev(mzones+1)    = ztekev(mzones)
        ztikev(mzones+1)    = ztikev(mzones)
        ztemp1(mzones+1)    = ztemp1(mzones)
        ztemp2(mzones+1)    = ztemp2(mzones)
        ztemp3(mzones+1)    = ztemp3(mzones)
c
        zxzeff(mzones+1)    = zxzeff(mzones)
        zrminor(mzones+1)   = 2 * zrminor(mzones) - zrminor(mzones-1)
        zrmajor(mzones+1)   = zrmajor(mzones)
c
c..initialize
c
c
        allocate ( zexb_out2(1:55,10) )
c
        zexb_out1  = 0.
        zexb_out2  = 0.
        zexb_out3  = 0.
c
c..call the GLF23 model one zone boundary at a time
c
        do jz=i1,mzones
c
          imm = jz - 1  ! zone boundary index
c
c..set switches
c
          ieigen = 1  ! use the tomsqz eigenvalue solver
          ishoot = 0  ! time-dependent code
c          igrid  = 0  ! do full grid of ipoints grid points
          i_grad = 1  ! compute gradients externally
          idengrad = 2 ! for simple dilution
c
          iswitch(1) = 1  ! ion particle transport
          iswitch(2) = 1  ! electron thermal transport
          iswitch(3) = 1  ! ion thermal transport
          iswitch(4) = 0
          iswitch(5) = 0
c
c..scalars
c
          zbtor_axis = zbtor(iaxis)
                       ! toroidal field at magnetic axis (tesla)
          ztor_flux  = fltors(mzones,1) * usib * usil**2
                       ! toroidal flux at the edge of the plasma( tesla m^2)
          zrmaj_axis = zrmajor(iaxis)
                       ! major radius at magnetic axis
          zsw_alpha  = 0
          if ( cthery(129) .gt. epslon ) zsw_alpha  = 1
                       ! include effect of flow shear stabilization
          zsw_x_alpha = 0
          if ( cthery(119).gt. epslon ) zsw_x_alpha = 1
                       ! include alpha effec (see zalpha above)
          i_delay    = 0
                       ! time delay for flow shear
c
          zscalgrdte = zdrdrho(jz) * zsgrdte(jz) / zrmajor(jz)
          zscalgrdti = zdrdrho(jz) * zsgrdti(jz) / zrmajor(jz)
          zscalgrdne = zdrdrho(jz) * zsgrdne(jz) / zrmajor(jz)
          zscalgrdni = zdrdrho(jz) * zsgrdni(jz) / zrmajor(jz)
c
          zavezimps  = zavezimp(jz)
          zmassimps  = zmassimp(jz)
          zmasshyds  = zmasshyd(jz)
c
c..initialize
c
          zstemp1 = 0.
          zstemp2 = 0.
          zstemp3 = 0.
          zstemp4 = 0.
c
          zvtor   = 0.
          zvpara  = 0.
          zvperp  = 0.
c
          zschie  = 0.
          zschii  = 0.
          zsvtor  = 0.
          zsvpara = 0.
          zsvperp = 0.
          zsexchange = 0.
c
c..For the moment, turn off flow shear and alpha stabilization
c
          zsw_alpha  = 0.
          zsw_x_alpha= 0.
c
c..call the GLF23 model
c
c  Use zthikb, ... as dummy output arrays which are zeroed out later
c
          call callglf2d ( ieigen, ishoot, imm, mzones, iswitch
     & , ztekev,   ztikev,   ztemp1,   ztemp2
     & , ztemp3,   i_grad,       idengrad
     & , zscalgrdte,   zscalgrdti,   zscalgrdne,   zscalgrdni
     & , zvrotxb,  zwexbxb,  zyexbxb,  zvtor
     & , zvpara,   zvperp,   zxzeff,   zbtor_axis
     & , xbouni,   ztor_flux,    zgradrsqrave, zgradrave
     & , zrminor,  zrmajor,  zrmaj_axis
     & , zavezimps, zmassimps, q,        zdqdrho
     & , zalpha,   zelong,   zmasshyds
     & , zsw_alpha,    zsw_x_alpha,  i_delay,      zsdifi
     & , zschie, zschii, zsvtor, zsvpara, zsvperp, zsexchange
     & , zthdkb,   zthekb,   zthikb
     & , zthvtor,  zthvpara, zthvperp, zweithe
     & , zexb_out1, zexb_out2, zexb_out3
     & , zgammaitg, zomegaitg
     & , zgammatem, zomegatem )
c
c..set output arrays
c
c..set impurity transport equal to ion particle transport
c
c..move negative diffusivities to convective velocities
c
c..use fig(1) to fig(4) for user control of coefficients
c
          if ( zschii .gt. 0. ) then
            difthi(1,1,jz) = difthi(1,1,jz) + fig(3) * zschii
          else
            velthi(1,jz) = velthi(1,jz)
     &         - fig(3) * zschii * zsgrdti(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          zthiig(jz) = max ( fig(3) * zschii, 0.0 )
c
          if ( zsdifi .gt. 0. ) then
            difthi(2,2,jz) = difthi(2,2,jz) + fig(1) * zsdifi
            difthi(4,4,jz) = difthi(4,4,jz) + fig(4) * zsdifi
          else
            velthi(2,jz) = velthi(2,jz)
     &         - fig(1) * zsdifi * zsgrdnh(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
             velthi(4,jz) = velthi(4,jz)
     &         - fig(4) * zsdifi * zsgrdnz(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          zthdig(jz) = max ( fig(1) * zsdifi, 0.0 )
          zthzig(jz) = max ( fig(4) * zsdifi, 0.0 )
c
          if ( zschie .gt. 0. ) then
            difthi(3,3,jz) = difthi(3,3,jz) + fig(2) * zschie
          else
            velthi(3,jz) = velthi(3,jz)
     &         - fig(2) * zschie * zsgrdte(jz) * zgradrsqrave(jz)
     &           / ( zgradrave(jz) * zrmajor(jz) )
          endif
          ztheig(jz) = max ( fig(2) * zschie, 0.0 )
c
c
         enddo
c
         deallocate ( zexb_out2 )
c
c..zero out the zthikb, ... arrays that were used as dummy output
c
         zthikb = 0.0
         zthdkb = 0.0
         zthekb = 0.0
         zthzkb = 0.0
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..end of transport models
c
        endif
c
c..zero out diffusivities up to magnetic axis
c
          do jz=1,2
c
            zthigb(jz)   = 0.0
            zthegb(jz)   = 0.0
            zthibohm(jz) = 0.0
            zthebohm(jz) = 0.0
            zthimix(jz)  = 0.0
            zthemix(jz)  = 0.0
            zthdmix(jz)  = 0.0
            zthzmix(jz)  = 0.0
c
            do j1=1,4
              velthi(j1,jz) = 0.0
              do j2=1,4
                difthi(j1,j2,jz) = 0.0
              enddo
            enddo
c
          enddo
c
c..set local diffusivities as needed
c  when not using sbrtn theory
c
        if ( lthery(21) .gt. 0 ) then
cap / ...
          if ( lthery(39) .ge. 0) then
            do jz=iaxis,iedge
              zxithe(jz) = difthi(1,1,jz)
              zvithe(jz) = velthi(1,jz)
              zdhthe(jz) = difthi(2,2,jz)
              zvhthe(jz) = velthi(2,jz)
              zxethe(jz) = difthi(3,3,jz)
              zvethe(jz) = velthi(3,jz)
              zdzthe(jz) = difthi(4,4,jz)
              zvzthe(jz) = velthi(4,jz)
              weithe(jz) = 0.0
c
              velthi(1,jz) = 0.0
              velthi(3,jz) = 0.0
              difthi(1,1,jz) = 0.0
              difthi(2,2,jz) = 0.0
              difthi(3,3,jz) = 0.0
              difthi(4,4,jz) = 0.0
            enddo
          else
            zxithe = 0.
            zvithe = 0.
            zdhthe = 0.
            zvhthe = 0.
            zxethe = 0.
            zvethe = 0.
            zdzthe = 0.
            zvzthe = 0.
          endif
cap \
c
c..fill in the diffusivities at the magnetic axis
c  (needed for the non-equilibrium impurity radiation model)
c
            zxithe(iaxis) = zxithe(iaxis+1)
            zdhthe(iaxis) = zdhthe(iaxis+1)
            zxethe(iaxis) = zxethe(iaxis+1)
            zdzthe(iaxis) = zdzthe(iaxis+1)
c
c..grow rates and frequencies
c
c  zgammaitg(jr) = growth rate of fastest growing ITG mode
c  zomegaitg(jr) = frequency of fastest growing ITG mode ( < 0 )
c  zgammatem(jr) = growth rate of fastest growing TEM mode
c  zomegatem(jr) = frequency of fastest growing TEM mode ( > 0 )
c
      ieq = imatdim   !  lthery(7)   ! number of eigenvalues
c
      do jz=1,ipoints
c
        zgammaitg(jz) = 0.0
        zgammatem(jz) = 0.0
c
        do jm=1,ieq
c
          if ( zgamma(jm,jz) .gt. 1.0 ) then
c
            if ( zomega(jm,jz) .lt. 0.0 ) then
c
              if ( zgamma(jm,jz) .gt. zgammaitg(jz) ) then
                zgammaitg(jz) = zgamma(jm,jz)
                zomegaitg(jz) = zomega(jm,jz)
              endif
c
            else
c
              if ( zgamma(jm,jz) .gt. zgammatem(jz) ) then
                zgammatem(jz) = zgamma(jm,jz)
                zomegatem(jz) = zomega(jm,jz)
              endif
c
            endif
c
          endif
c
        enddo
c
      enddo
c
!| 
!| In order to limit the time rate of change of the theory-based
!| diffusivities, the old values are kept in local arrays.
!| The difference between the old and the new values is computed
!| \[ \Delta \chi_j = \chi^{N}_j - \chi^{N-1}_j \]
!| and stored in local arrays.
!| This difference is then spatially averaged
!| \[ \bar{ \Delta \chi_j }
!|     = ( \Delta \chi_{j-1} + 2 \Delta \chi_j +  \Delta \chi_{j+1} ) / 4 \]
!| Finally, the adjusted diffusivity is computed
!| \[ \chi^N_j = \chi^{N-1}_j
!|  + \Delta \chi_j /
!|  ( 1. + c_{60} |  \Delta \chi_j - \bar{ \Delta \chi_j} | ). \]
!| Here $ c_{60} = {\tt cthery(60)} $ with default value 0.0.
!| A recommended value is $ {\tt cthery(60)} = 10.0 $
!| if the diffusivities show the pattern of a numerical instability.
!| This adjustment suppresses large local changes in the diffusivities.
!| Here, this algorithm is applied only to the ITG ($\eta_i$) mode.
!| 
      if ( cthery(60) .gt. epslon .and.
     & (lthery(21) .ne. 7 .or. lthery(21) .ne. 8)) then
c
c..First, subtract the ITG diffusivites from the total diffusivities
c
        do jz=1,iedge
          zxithe(jz) = zxithe(jz) - zthiig(jz)
          zdhthe(jz) = zdhthe(jz) - zthdig(jz)
          zxethe(jz) = zxethe(jz) - ztheig(jz)
          zdzthe(jz) = zdzthe(jz) - zthzig(jz)
        enddo
c
c..Store the ITG diffusivities in zoeta* before first step
c
        if ( nstep .lt. 2 ) then
          do jz=1,iedge
            zoetai(jz) = zthiig(jz)
            zoetae(jz) = ztheig(jz)
            zoetad(jz) = zthdig(jz)
            zoetaz(jz) = zthzig(jz)
          enddo
        endif
c
c..Compute difference between diffusivities now and before
c
        do jz=1,iedge
          zdleti(jz) = zthiig(jz) - zoetai(jz)
          zdlete(jz) = ztheig(jz) - zoetae(jz)
          zdletd(jz) = zthdig(jz) - zoetad(jz)
          zdletz(jz) = zthzig(jz) - zoetaz(jz)
        enddo
c
c..Spatially average diffusivities
c
        do jz=iaxis+1,iedge-1
          ztemp1(jz)=0.25*(zdleti(jz-1)+2.0*zdleti(jz)+zdleti(jz+1))
          ztemp2(jz)=0.25*(zdlete(jz-1)+2.0*zdlete(jz)+zdlete(jz+1))
          ztemp3(jz)=0.25*(zdletd(jz-1)+2.0*zdletd(jz)+zdletd(jz+1))
          ztemp4(jz)=0.25*(zdletz(jz-1)+2.0*zdletz(jz)+zdletz(jz+1))
        enddo
c
c..compute adjusted diffusivities
c
        do jz=iaxis+1,iedge-1
          zthiig(jz) = zoetai(jz) + zdleti(jz)
     &     / ( 1. + cthery(60) * abs ( zdleti(jz) - ztemp1(jz) ) )
          ztheig(jz) = zoetae(jz) + zdlete(jz)
     &     / ( 1. + cthery(60) * abs ( zdlete(jz) - ztemp2(jz) ) )
          zthdig(jz) = zoetad(jz) + zdletd(jz)
     &     / ( 1. + cthery(60) * abs ( zdletd(jz) - ztemp3(jz) ) )
          zthzig(jz) = zoetaz(jz) + zdletz(jz)
     &     / ( 1. + cthery(60) * abs ( zdletz(jz) - ztemp4(jz) ) )
        enddo
c
c..store the new diffusivites if this step has not been repeated
c
        if (  nstep .gt. istep ) then
          istep = nstep
          do jz=1,iedge
            zoetai(jz) = zthiig(jz)
            zoetae(jz) = ztheig(jz)
            zoetad(jz) = zthdig(jz)
            zoetaz(jz) = zthzig(jz)
          enddo
        endif
c
c..add the ITG diffusivities back into the total diffusivities
c
        do jz=1,iedge
          zxithe(jz) = zxithe(jz) + zthiig(jz)
          zdhthe(jz) = zdhthe(jz) + zthdig(jz)
          zxethe(jz) = zxethe(jz) + ztheig(jz)
          zdzthe(jz) = zdzthe(jz) + zthzig(jz)
        enddo
c
      endif
c
        endif
c
c..additional output from the Multi-Mode model
c
c..Some pages of long printout
c
        if ( knthe .eq. 3  .and.  lthery(29) .gt. 2 ) then
c
          zt  = tai * uist * 1000.0
          zdt = dtoldi * uist * 1000.0
c
c..Total diffusivities if more printout is desired
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c..total theory-based diffusivities
c
       if (lthery(39) >= 0) then
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
         write(nprint,10007)
         write(nprint,10009)
         write(nprint,10008)
c
        do jz=iaxis,mzones
          write(nprint,101) jz, zrminor(jz), zxethe(jz)
     &       , zxithe(jz), zweithe(jz), zdhthe(jz), zdzthe(jz)
        enddo
       endif ! lthery(39) >= 0
c
c..totals, including neoclassical and empirical effective diffusivities
c
        do jz=1,mzones
          zxetot(jz) = zxethe(jz)
     &      + (cfutz(11) * xeneo1(jz) + xeneo2(jz)
     &      + xeemps(jz))*usil**2
          zxitot(jz) = zxithe(jz)
     &      + (cfutz(12) * xineo1(jz)
     &      + xiemps(jz))*usil**2
        enddo
c
c..ion thermal diffusivities
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write(nprint,107)
c
        if(lthery(21) .eq. 7 .or. lthery(21) .eq. 8)  then
                write(nprint,121)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz)
     &                  , zthibohm(jz), zthigb(jz), zthimix(jz)
     &                  , cfutz(12)*xineo1(jz)*usil**2
     &                  , xiemps(jz)*usil**2, zxitot(jz)
                enddo

        else
                write(nprint,122)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz)
     &                  , zthiig(jz), zthirb(jz), zthikb(jz)
     &                  , zxithe(jz), cfutz(12)*xineo1(jz)
     &                    *usil**2
     &                  , xiemps(jz)*usil**2, zxitot(jz)
                enddo
        endif
c
c..electron thermal diffusivities
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write(nprint,104)

c
        if(lthery(21) .eq. 7 .or. lthery(21) .eq. 8) then
                write(nprint,123)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz)
     &                  , zthebohm(jz), zthegb(jz), zthemix(jz)
     &                  , (cfutz(11)*xeneo1(jz)+xeneo2(jz))*usil**2
     &                  , xeemps(jz)*usil**2, zxetot(jz)
                enddo
c
        else
c
                write(nprint,124)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz)
     &                  , ztheig(jz), ztherb(jz), zthekb(jz)
     &                  , ztheeg(jz)
     &                  , zthetb(jz)
     &                  , zxethe(jz), (cfutz(11)*xeneo1(jz)
     &                    +xeneo2(jz))*usil**2
     &                  , xeemps(jz)*usil**2, zxetot(jz)
                enddo
        endif
c
c
c..hydrogenic particle diffusivities
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write(nprint,125)
c
        if(lthery(21) .eq. 7 .or. lthery(21) .eq. 8) then
                write(nprint,129)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz), zthdmix(jz)
                enddo
c
        else
c
                write(nprint,126)
c
                do jz=iaxis,mzones
                 write(nprint,120)  zrminor(jz)
     &                  , zthdig(jz), zthdrb(jz)
     &                  , zthdkb(jz), zdhthe(jz)

                enddo

        endif
c
c..impurity particle diffusivities
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write(nprint,127)
c
        if(lthery(21) .eq. 7 .or. lthery(21) .eq. 8) then
                write(nprint,131)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz), zthzmix(jz)
                enddo
c
        else
c
                write(nprint,128)
c
                do jz=iaxis,mzones
                write(nprint,120)  zrminor(jz)
     &                  , zthzig(jz), zthzrb(jz)
     &                  , zthzkb(jz), zdzthe(jz)
                enddo
        endif
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'Normalized gradients:'
c
        write (nprint,138)
 138    format ('zrminor',t15,'grdne',t27,'grdni'
     &    ,t39,'grdnh',t51,'grdnz',t63,'grdte',t75,'grdti'
     &    ,t87,'grdpr',t99,'grdq')
        do j=1,iedge
          write (nprint,152) zrminor(j), zgrdne(j), zgrdni(j)
     &      , zgrdnh(j), zgrdnz(j), zgrdte(j), zgrdti(j)
     &      , zgrdpr(j), zgrdq(j)
        enddo
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'Smoothed normalized gradients:'
c
        write (nprint,139)
 139    format ('zrminor',t15,'sgrdne',t27,'sgrdni'
     &    ,t39,'sgrdnh',t51,'sgrdnz',t63,'sgrdte',t75,'sgrdti'
     &    ,t87,'sgrdpr',t99,'sgrdq')
        do j=1,iedge
          write (nprint,152) zrminor(j), zsgrdne(j), zsgrdni(j)
     &      , zsgrdnh(j), zsgrdnz(j), zsgrdte(j), zsgrdti(j)
     &      , zsgrdpr(j), zsgrdq(j)
        enddo
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'Variables used by callglf2db:'
c
        write (nprint,149)
 149    format ('xbouni',t15,'zvrotxb',t27,'zwexbxb'
     &    ,t39,'zalpha',t51,'zvtor',t63,'zvpara',t75,'zvperp'
     &    ,t87,'zgradrsqrave',t99,' zgradrave')
        do j=1,iedge
c          if(zvrotNTV(j) .ne. zvrotNTV(j) ) then
c            zvrotNTV(j) = 0.0001
c          endif
          write (nprint,152) xbouni(j), zvrotxb(j), zwexbxb(j)
     &      , zalpha(j), zvtor(j), zvpara(j), zvperp(j)
     &      , zgradrsqrave(j), zgradrave(j)
        enddo
c
c PUB added the following line (use for writing output to jfile)
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'Toroidal velocities:'
c
        write (nprint,301)
 301    format ('xbouni',t15,'zvrotxb',t27,'zvrotNTV'
     &    ,t39,'zvrotNEO',t51,'zvrottorque',t63,'bpoli',t75,'zgrdti')
        do j=2,iedge
          if(zvrotNTV(j) .ne. zvrotNTV(j) ) then
            zvrotNTV(j) = zvrotNTV(j+1)
          endif
          if(zvrotxb(j) .ne. zvrotxb(j) ) then
            zvrotxb(j) = zvrotxb(j+1)
          endif
          write (nprint,152) xbouni(j), zvrotxb(j), zvrotNTV(j)
     &      , zvrotNEO(j), zvrottorque(j), bpoli(j), zgrdti(j)
        enddo
c
c..grow rates and frequencies
c
c  zgammaitg(jz) = growth rate of fastest growing ITG mode
c  zomegaitg(jz) = frequency of fastest growing ITG mode ( < 0 )
c  zgammatem(jz) = growth rate of fastest growing TEM mode
c  zomegatem(jz) = frequency of fastest growing TEM mode ( > 0 )
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'Growth rates and frequencies:'
c
        write (nprint,180)
  180   format ('zrminor',t15,'gamma_ITG',t27,'omega_ITG'
     &    ,t39,'gamma_TEM',t51,'omega_TEM')
c
        do jz=1,iedge
          write (nprint,152) zrminor(jz)
     &      , zgammaitg(jz), zomegaitg(jz)
     &      , zgammatem(jz), zomegatem(jz)
        enddo
c
c..trapped electron mode vs ITG mode diffusivities
c
c  zthitem(jz) = ion thermal diffusivity from TEM (diagnostic output)
c  zthdtem(jz) = hydrogenic ion diffusivity from TEM (diagnostic output)
c  zthetem(jz) = elelctron thermal diffusivity from TEM (diagnostic output)
c  zthztem(jz) = impurity ion diffusivity from TEM (diagnostic output)
c
c  zthiitg(jz) = ion thermal diffusivity from ITG (diagnostic output)
c  zthditg(jz) = hydrogenic ion diffusivity from ITG (diagnostic output)
c  ztheitg(jz) = elelctron thermal diffusivity from ITG (diagnostic output)
c  zthzitg(jz) = impurity ion diffusivity from ITG (diagnostic output)
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'ITG and TEM thermal diffusivities:'
c
        write (nprint,182)
 182    format ('zrminor',t15,'chi_i_ITG',t27,'chi_i_TEM'
     &    ,t39,'chi_e_ITG',t51,'chi_e_TEM')
c
        do jz=1,iedge
          write (nprint,152) zrminor(jz)
     &      , zthiitg(jz), zthitem(jz)
     &      , ztheitg(jz), zthetem(jz)
        enddo
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
        write (nprint,*)
        write (nprint,*) 'ITG and TEM ion diffusivities:'
c
        write (nprint,183)
 183    format ('zrminor',t15,'D_H_ITG',t27,'D_H_TEM'
     &    ,t39,'D_Z_ITG',t51,'D_Z_TEM')
c
        do jz=1,iedge
          write (nprint,152) zrminor(jz)
     &      , zthditg(jz), zthdtem(jz)
     &      , zthzitg(jz), zthztem(jz)
        enddo
c
c
        endif
c
      endif
c
c..diagnostic output as needed
c  when nstep = lthery(26)
c  or when tai*uiet .ge. tplot(j) .gt. epslon  and  lthery(26) .gt. 0
c
      idiag = 0
      if ( nstep .eq. lthery(26) ) idiag = 1
      if ( nstep .lt. -lthery(26) ) idiag = 1
      if ( itdiag .lt. 20  .and.  lthery(26) .gt. 0 ) then
        do j=itdiag,20
          if ( tplot(j) .gt. epslon
     &         .and. tai*uiet .gt. tplot(j) ) then
            idiag = 1
            itdiag = j + 1
          endif
        enddo
      endif
c
      if ( idiag .gt. 0 ) then
c
        write (nprint,150)
 150    format(/'Diagnostic output from sbrtn ptheory')
c
        write (nprint,151) nstep, ztime
 151    format ('# nstep = ',i5,' time = ',0pf14.6)
 152    format (1p10e12.4)
c
        write (nprint,154) (lthery(j),j=1,50)
 154    format ('# lthery'/,(10i5))
c
        write (nprint,155) (cthery(j),j=1,150)
 155    format ('# cthery(j)'/,(1p10e12.4))
c
        write (nprint,156) iaxis, iedge, isep, imatdim
     &    , nprint, iprint
 156    format ('# maxis = ',i5,'  medge = ',i5,'  mseprtx = ',i5
     &    ,'  matdim = ',i5,'  nprint = ',i5,'  lprint = ',i5)
c
        write (nprint,157)
 157    format ('#  xbouni',t15,'rminor',t27,'rmajor',t39,'elong'
     &    ,t51,'triang',t63,'zgfactr')
        do j=1,iedge
          write (nprint,152) xbouni(j), zrminor(j), zrmajor(j)
     &      , zelong(j), ztriang(j), zgfactr(j)
        enddo
c
        write (nprint,158)
 158    format ('#  xbouni',t15,'dense',t27,'densi'
     &    ,t39,'densh',t51,'densimp',t63,'densf',t75,'densfe')
        do j=1,iedge
          write (nprint,152) xbouni(j), zdense(j), zdensi(j)
     &      , zdensh(j), zdensimp(j), zdensf(j), zdensfe(j)
        enddo
c
        write (nprint,159)
 159    format ('#  xbouni',t15,'xzeff',t27,'tekev'
     &    ,t39,'tikev',t51,'tfkev',t63,'q',t75,'vloop'
     &    ,t87,'btor',t99,'resist')
        do j=1,iedge
          write (nprint,152) xbouni(j), zxzeff(j), ztekev(j)
     &      , ztikev(j), ztfkev(j), q(j), zvloop(j)
     &      , zbtor(j), zresist(j)
        enddo
c
        write (nprint,160)
 160    format ('#  xbouni',t15,'avezimp',t27,'amassimp'
     &    ,t39,'amasshyd',t51,'aimass')
        do j=1,iedge
          write (nprint,152) xbouni(j), zavezimp(j), zmassimp(j)
     &      , zmasshyd(j), zaimass(j)
        enddo
c
        write (nprint,161)
 161    format ('#  xbouni',t15,'grdne',t27,'grdni'
     &    ,t39,'grdnh',t51,'grdnz',t63,'grdte',t75,'grdti'
     &    ,t87,'grdpr',t99,'grdq')
        do j=1,iedge
          write (nprint,152) xbouni(j), zgrdne(j), zgrdni(j)
     &      , zgrdnh(j), zgrdnz(j), zgrdte(j), zgrdti(j)
     &      , zgrdpr(j), zgrdq(j)
        enddo
c
        write (nprint,181)
 181    format ('#  xbouni',t15,'sgrdne',t27,'sgrdni'
     &    ,t39,'sgrdnh',t51,'sgrdnz',t63,'sgrdte',t75,'sgrdti'
     &    ,t87,'sgrdpr',t99,'sgrdq')
        do j=1,iedge
          write (nprint,152) xbouni(j), zsgrdne(j), zsgrdni(j)
     &      , zsgrdnh(j), zsgrdnz(j), zsgrdte(j), zsgrdti(j)
     &      , zsgrdpr(j), zsgrdq(j)
        enddo
c
        write (nprint,162)
 162    format ('#  fdr',t15,'fig',t27,'fti'
     &    ,t39,'frm',t51,'fkb',t63,'frb',t75,'fhf'
     &    ,t87,'fec',t99,'fmh',t111,'feg')
        do j=1,5
          write (nprint,152) fdr(j), fig(j), fti(j)
     &      , frm(j), fkb(j), frb(j), fhf(j)
     &      , fec(j), fmh(j), feg(j)
        enddo
c
        write (nprint,*) '# fdrint = ',fdrint
c
       if (lthery(39) >= 0) then
         write (nprint,163)
 163     format ('#  xbouni',t15,'dhtot',t27,'vhtot'
     &     ,t39,'dztot',t51,'vztot',t63,'xetot',t75,'xitot'
     &     ,t87,'wiethe')
         do j=1,iedge
           write (nprint,152) xbouni(j), zdhthe(j), zvhthe(j)
     &       , zdzthe(j), zvzthe(j), zxethe(j), zxithe(j)
     &       , zweithe(j)
         enddo
c
       else
        write (nprint,*)
        write (nprint,*) ' Diffusion matrix:'
        write (nprint,171)

        do jr=iaxis,mzones
          write (nprint,110) zrminor(jr)
     &      , velthi(1,jr),   difthi(1,1,jr), difthi(1,2,jr)
     &      , difthi(1,3,jr), difthi(1,4,jr), difthi(1,5,jr)
        enddo
c
        write (nprint,*)
        write (nprint,172)

        do jr=iaxis,mzones
          write (nprint,110) zrminor(jr)
     &      , velthi(2,jr),   difthi(2,1,jr), difthi(2,2,jr)
     &      , difthi(2,3,jr), difthi(2,4,jr), difthi(2,5,jr)
        enddo
c
        write (nprint,*)
        write (nprint,173)

        do jr=iaxis,mzones
          write (nprint,110) zrminor(jr)
     &      , velthi(3,jr),   difthi(3,1,jr), difthi(3,2,jr)
     &      , difthi(3,3,jr), difthi(3,4,jr), difthi(3,5,jr)
        enddo
c
        write (nprint,*)
        write (nprint,174)

        do jr=iaxis,mzones
          write (nprint,110) zrminor(jr)
     &      , velthi(4,jr),   difthi(4,1,jr), difthi(4,2,jr)
     &      , difthi(4,3,jr), difthi(4,4,jr), difthi(4,5,jr)
        enddo
       endif ! lthery(39) < 0
c
      endif
c
!| 
!| Convert to cgs units, as needed.
!| 
c
      do jz=iaxis,iedge
c
        xethes(jz) = zxethe(jz)
        xithes(jz) = zxithe(jz)
c
        do ji=lhyd1,lhydn
          dxthes(ji,jz) = zdhthe(jz)
cbate          vxthes(ji,jz) = zvhthe(jz)
        enddo
c
        do ji=limp1,limpn
          dxthes(ji,jz) = zdzthe(jz)
cbate          vxthes(ji,jz) = zvzthe(jz)
        enddo
c
c  put 'weithe' in standard units to be used in sub. convrt (dsolver)
c
        weiths(jz) = weithe(jz) * uesh * uisd
c
c  rgb 15-apr-95 appended factor elong(jz) to volume expression
c
        zvolum=2*fcpi*zrmajor(jz)*fcpi*(zrminor(jz)**2-
     #       zrminor((jz-1))**2)*zelong(jz)
        eithes(jz) = weithe(jz)*zvolum
c
      enddo
c
!| 
!| \subsection{Printout}
!| 
c-----------------------------------------------------------------------
c
c  print theory's output
c
      if ( knthe .ne. 3 ) go to 990
c
cbate      entry prethprnt
c
 900  continue
c
      zt  = tai * uist * 1000.0
      zdt = dtoldi * uist * 1000.0
c
c..Profiles as a function of major radius
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      icntr = lcentr
      iedge = ledge + 1
c
      ir = 0
      do 802 jz=iedge,icntr,-1
        ir = ir + 1
        zrmajm(ir) = ( rmids(jz,2) - ahalfs(jz,2) ) * usil
        znemaj(ir) = rhoels(2,jz) * usid
        ztemaj(ir) = tes(2,jz) * useh
        ztimaj(ir) = tis(2,jz) * useh
        zefmaj(ir) = xzeff(2,jz)
 802  continue
c
      do 804 jz=icntr,iedge
        ir = ir + 1
        zrmajm(ir) = ( rmids(jz,2) + ahalfs(jz,2) ) * usil
        znemaj(ir) = rhoels(2,jz) * usid
        ztemaj(ir) = tes(2,jz) * useh
        ztimaj(ir) = tis(2,jz) * useh
        zefmaj(ir) = xzeff(2,jz)
 804  continue
c
      irmax = ir
c
      write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write (nprint,130)
c
      do 806 jz=1,irmax
        write (nprint,132) zrmajm(jz),znemaj(jz),ztemaj(jz),ztimaj(jz)
     &    ,zefmaj(jz)
 806  continue

c
c..Densities as a function of minor radius
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      if ( lthery(29) .gt. 2 ) then
c
      lpage=lpage+1
      write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write ( nprint, 133)
c
       do jz=iaxis,mzones
        write (nprint,106) jz, zrminor(jz)
     &    , zdense(jz), zdensi(jz)
     &    , zdensh(jz), zdensimp(jz)
     &    , zdensf(jz), zdensfe(jz), zxzeff(jz), zaimass(jz)
      enddo
c
      endif
c
c..print diffusivity matrix
c%%%%%%%%%%%%%%
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write (nprint,*)
      write (nprint,*) ' Convective velocities (not used):'
      write (nprint,170)
 170  format (/t2,'rminor'
     &  ,t10,'zvithe',t21,'zvhthe',t33,'zvethe'
     &  ,t45,'zvzthe')

      do jr=iaxis,mzones
        write (nprint,110) zrminor(jr)
     &    , zvithe(jr), zvhthe(jr), zvethe(jr), zvzthe(jr)
      enddo
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write (nprint,*)
      write (nprint,*) ' Diffusion matrix:'
      write (nprint,171)
 171  format (/t2,'rminor'
     &  ,t10,'vel(1,jr)',t21,'dif(1,1,jr)',t33,'dif(1,2,jr)'
     &  ,t45,'dif(1,3,jr)',t57,'dif(1,4,jr)',t69,'dif(1,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=iaxis,mzones
        write (nprint,110) zrminor(jr)
     &    , velthi(1,jr),   difthi(1,1,jr), difthi(1,2,jr)
     &    , difthi(1,3,jr), difthi(1,4,jr), difthi(1,5,jr)
      enddo
c
      write (nprint,*)
      write (nprint,172)
 172  format (/t2,'rminor'
     &  ,t10,'vel(2,jr)',t21,'dif(2,1,jr)',t33,'dif(2,2,jr)'
     &  ,t45,'dif(2,3,jr)',t57,'dif(2,4,jr)',t69,'dif(2,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=iaxis,mzones
        write (nprint,110) zrminor(jr)
     &    , velthi(2,jr),   difthi(2,1,jr), difthi(2,2,jr)
     &    , difthi(2,3,jr), difthi(2,4,jr), difthi(2,5,jr)
      enddo
c
      write (nprint,*)
      write (nprint,173)
 173  format (/t2,'rminor'
     &  ,t10,'vel(3,jr)',t21,'dif(3,1,jr)',t33,'dif(3,2,jr)'
     &  ,t45,'dif(3,3,jr)',t57,'dif(3,4,jr)',t69,'dif(3,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=iaxis,mzones
        write (nprint,110) zrminor(jr)
     &    , velthi(3,jr),   difthi(3,1,jr), difthi(3,2,jr)
     &    , difthi(3,3,jr), difthi(3,4,jr), difthi(3,5,jr)
      enddo
c
      write (nprint,*)
      write (nprint,174)
 174  format (/t2,'rminor'
     &  ,t10,'vel(4,jr)',t21,'dif(4,1,jr)',t33,'dif(4,2,jr)'
     &  ,t45,'dif(4,3,jr)',t57,'dif(4,4,jr)',t69,'dif(4,5,jr)',
     &   /t4,'[m]',t12,'[m/s]',t23,'[m^2/s]',t35,'[m^2/s]',
     &   t47,'[m^2/s]',t59,'[m^2/s]',t71,'[m^2/s]')

      do jr=iaxis,mzones
        write (nprint,110) zrminor(jr)
     &    , velthi(4,jr),   difthi(4,1,jr), difthi(4,2,jr)
     &    , difthi(4,3,jr), difthi(4,4,jr), difthi(4,5,jr)
      enddo
c
c..print effective convective velocities
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write ( nprint, 140)
c
      do jz=iaxis,mzones
c
        write (nprint,106) jz, zrminor(jz)
     &    , vftot(jz,lhyd1), vftot(jz,limp1)
     &    , vftot(jz,lion), vftot(jz,lelec)
      enddo
c
c..print fluxes
c%%%%%%%%%%%%%%
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write ( nprint, 141)
c
      do jz=iaxis,mzones
c
        write (nprint,106) jz, zrminor(jz)
     &    , flxtot(jz,lhyd1), flxtot(jz,limp1)
     &    , flxtot(jz,lion), flxtot(jz,lelec)
      enddo
c
c..print sources
c%%%%%%%%%%%%%%%
c
        lpage=lpage+1
        write(nprint,103) label1(1:48),label5(1:72),lpage,nstep,zt,zdt
c
      write ( nprint, 142)
c
      do jz=iaxis,mzones
        write (nprint,106) jz, zrminor(jz)
     &    , srctot(jz,lhyd1), srctot(jz,limp1)
     &    , srctot(jz,lion), srctot(jz,lelec)
      enddo
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
10007 format(/,10x,'transport coefficients from theory',/)
10008 format(4x,'zone',6x,'radius',9x,'chi-elc',8x,'chi-ion',
     #       8x,'intrchg',10x,'zdifh',10x,'zdifz')
10009 format(16x,'m',13x,'m*m/s',10x,'m*m/s',12x,'w',
     #       12x,'m*m/s')
c
 101  format(5x,i2,6x,0pf6.3,8x,5(1pe11.4,4x))
 102  format(/,15x,'total interchange power = ',2x,e11.4,3x,'watts')
 103  format(/2x,a48,10x,a72/
     &  2x,'-',i2,'-',2x,'*** time step',i5,' ***',14x,
     &          'time =',f12.3,2x,'millisecs.',12x,'dt =',f12.6,2x,
     &          'millisecs.')
 104  format(/,10x,'electron thermal diffusion coefficients',/,
     #       10x,43('-'))
 105  format(12x,'m',9x,'m2/s',7(8x,'m2/s'),/
     & 5x,'jz',3x,'radius',5x,'xethe',8x,'neocl',6x
     &         ,'empirc',6x,'xetot')
 106  format(5x,i2,3x,0pf6.3,4x,9(1pe10.3,2x))
 107  format(/,10x,'ion thermal diffusion coefficients',/,
     #       10x,37('-'))
 108  format(12x,'m',9x,'m2/s',7(8x,'m2/s'),/
     &  5x,'jz',3x,'radius',5x,'xithe'
     &         ,8x,'neocl',7x,'empirc',6x,'xitot')
 109  format(12x,'m',9x,'m2/s',7(8x,'m2/s'),/
     &  5x,'jz',3x,'radius',6x,'dhdr',8x,'dhig',8x,'dhti',8x
     &         ,'dhrm',8x,'dhrb',8x,'dhkb',8x,'dhnm'
     &         ,8x,'dhhf',8x,'dhtot')
c
 110  format (0p1f6.3,1p6e12.3)
c
 120  format(0pf6.3,4x,10(1pe10.3,2x))
 121  format(2x,'m',9x,'m2/s',5(8x,'m2/s'),/
     &  'radius',t12,'Xi_Bohm',t24,'Xi_gBohm',t36,'Xi_Mixed'
     &  ,t48,'Xi_Neo',t60,'Xi_Empirc',t72,'Xi_Total')
 122  format(2x,'m',9x,'m2/s',6(8x,'m2/s'),/
     &  'radius',t12,'thiig',t24,'thirb',t36,'thikb'
     &  ,t48,'xithe',t60,'neocl',t72,'empirc',t84,'xitot')
 123  format(2x,'m',9x,'m2/s',5(8x,'m2/s'),/
     &  'radius',t12,'Xe_Bohm',t24,'Xe_gBohm',t36,'Xe_Mixed'
     &  ,t48,'Xe_Neo',t60,'Xe_Empirc',t72,'Xe_Total')
 124  format(2x,'m',9x,'m2/s',6(8x,'m2/s'),/
     &  'radius',t12,'theig',t24,'therb',t36,'thekb',t48,'theeg'
     &  ,t60,'thetb'
     &  ,t72,'xethe',t84,'neocl',t96,'empirc',t108,'xetot')
 125  format(/,10x,'hydrogenic particle diffusion coefficients',/,
     #       10x,37('-'))
 126  format(2x,'m',9x,'m2/s',7(8x,'m2/s'),/
     &  'radius',t12,'thdig',t24,'thdrb',t36,'thdkb'
     &  ,t48,'dhthe',t60,'neocl',t72,'empirc',t84,'dhtot')
 127  format(/,10x,'impurity particle diffusion coefficients',/,
     #       10x,37('-'))
 128  format(2x,'m',9x,'m2/s',7(8x,'m2/s'),/
     &  'radius',t12,'tzdig',t24,'tzdrb',t36,'thzkb'
     &  ,t48,'dzthe',t60,'neocl',t72,'empirc',t84,'dztot')
 129  format(2x,'m',9x,'m2/s',(8x,'m2/s'),/
     &  'radius',t12,'X_Particle')
c
 130  format (
     & /,10x,'Profiles as a function of major radius'
     & /,t4,'rmajor(m)',t17,'ne(m^-3)',t30,'Te(keV)',t43,'Ti(keV)'
     &  ,t56,'Zeff')
 131  format(2x,'m',9x,'m2/s',(8x,'m2/s'),/
     &  'radius',t12,'X_Impuirity')
 132  format (5(2x,1pe11.4))
c
 133  format (
     & /,10x,'Densities as a function of minor radius',/
     &  ,5x,'jz',3x,'radius',8x,'ne',10x,'ni',10x,'nh',10x,'nz'
     &  ,10x,'ns',10x,'nse',8x,'zeff',9x,'mi')
c
 135  format (5x,'jz',3x,'radius',2x,'diffusivity',4x,'velthi(*)'
     &  ,'  difthi(*,1) difthi(*,2) difthi(*,3) difthi(*,4)'
     &  ,t95,'perform')
c
 140  format (/5x,'jz',3x,'radius',t22,'vftot_H'
     &  ,t34,'vftot_I',t46,'vftot_Ti',t58,'vftot_Te')
c
 141  format (/5x,'jz',3x,'radius',t22,'flxtot_H'
     &  ,t34,'flxtot_I',t46,'flxtot_Ti',t58,'flxtot_Te')
c
 142  format (/5x,'jz',3x,'radius',t22,'srctot_H',t34,'srctot_I'
     &  ,t46,'srctot_Ti',t58,'srctot_Te')
c
 990  return
      end
!| 
!| \end{document}             % End of document.
!| 

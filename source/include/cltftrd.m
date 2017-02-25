c@tftrdacl for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
c    the following cliche is from R. M. Weiland to be used to read
c  data from TFTR.  17 Oct 88
c
c=======            START OF TFTRDA.BLK        =============================

C       common /ctransps/               !SCALAR VARIABLES
C     &                 ctrmajr,        !MAJOR RADIUS (CM)
C     &                 ctrminr,        !MINOR RADIUS (CM)
C     &                 ctpcur,         !PLASMA CURRENT (AMPS)
C     &                 ctbz,           !EXT. TOROIDAL FIELD (TESLA)
C     &                 ctvsur,         !SURFACE VOLTAGE (VOLTS)
C     &                 ctbetat,        !PLASMA BETA (E+I ONLY) (POL) NOTE1
C     &                 ctbteq,         !EQUILIBRIUM BETA (TOROIDAL) NOTE1
C     &                 cttflux,        !TOROIDAL FLUX (WEBERS)
C     &                 ctpflux,        !POLOIDAL FLUX (WEBERS)
C     &                 ctli2pb,        !LAMBDA (TRANSP)
C     &                 ctul2pb,        !LAMBDA (MEASURED)
C     &                 ctlio2,         !LI/2 (TRANSP)
C     &                 ctbpeq,         !EQUILIBRIUM BETA (POLOIDAL) NOTE1
C     &                 ctq0,           !SAFETY FACTOR (q) AT EDGE
C     &                 cttime,         !TIME (SEC)
C     &                 ctrun           !TRANSP RUN #
C     &                 ctshot          !TOKAMAK SHOT #
C  NOTE1: ALL TRANSP BETAS ARE DEFINED WITH DENOMINATORS OF THE FORM
C         VOL-INT(B**2/(2U0))/VOL WHERE B=BP FOR POLOIDAL AND B=BT0 FOR
C         TOROIDAL BETAS

        common /ctransps/ 
     &                  ctrmajr,ctrminr,ctpcur,
     &                  ctbz,   ctvsur, ctbetat,
     &                  ctbteq, cttflux,ctpflux,
     &                  ctli2pb,ctul2pb,ctlio2, ctbpeq,
     &                  ctq0,   cttime, ctrun,
     &                  ctshot
        real
     &                  ctrmajr,
     &                  ctrminr,
     &                  ctpcur,
     &                  ctbz,
     &                  ctvsur,
     &                  ctbetat,
     &                  ctbteq,
     &                  cttflux,
     &                  ctpflux,
     &                  ctli2pb,
     &                  ctul2pb,
     &                  ctlio2,
     &                  ctbpeq,
     &                  ctq0,
     &                  cttime,
     &                  ctrun,
     &                  ctshot

C       THE FOLLOWING PROFILES ARE MAPPED ON ONE OF 2 EQUI-SPACED GRIDS OF 
C                 =======   SQRT( TOROIDAL FLUX ) ===========
C       SOME OF THE FOLLOWING ARRAYS ARE ON A ZONE-CENTERED GRID [ZC],
C       OTHERS ARE ON A ZONE-BOUNDARY GRID [ZB].

C       common /ctranspv/               !VECTOR (profile) VARIABLES
C     &                 ntrflx,cttrflx, !TOROIDAL FLUX (WEBERS)[ZB]
C     &                 npplas,ctpplas, !TOTAL MHD
C                                       !PRESSUREW/BEAMS(PASCALS)[ZC]
C     &                                 !{"PTOWB" not PPLAS as indicated}
C     &                                 !=PPLASMA+.5(PPERP+PPARL)+PROT
C     &                 nq,ctq,         !SAFETY FACTOR (Q)[ZB]
C     &                 nzefmd,ctzefmd, !MAGNETIC Z-EFFECTIVE[ZC]
C     &                 nzeffp,ctzeffp, !PLASMA COMPOSITION Z-EFF[ZC]
C     &                 nsshaf,ctsshaf, !SHAFRANOV SHIFT[ZB]
C     &                 nbttot,ctbttot, !TOTAL BETA TOROIDAL[ZC] -
C                                       !NOT "EQUILIBRIUM" BETA
C     &                 ncur,ctcur,     !PLASMA CUR-DENS (AMP/CM**2)[ZC]
C     &                 nne,ctne,       !ELECTRON DENSITY (#/CM**3)[ZC]
C     &                 nplflx,ctplflx, !POLOIDAL FLUX (WEBERS)[ZB]
C     &                 ngfun,ctgfun    !G-Function [R*Bt/R0/Bt0] [ZB]
C     &                 nrmnmp,ctrmnmp  !Minor Radius (Cm)  [ZB]
C     &                 nmomr,nxmomr,ctrmm(0:nmomr,1:nxmomr)
C     &                                 !R Fourier Moments [Rmm(m,x)] (Cm)
C     &                 nmomy,nxmomy,ctymm(0:nmomy,1:nxmomy)
C     &                                 !Y Fourier Moments [Rmm(m,x)] (Cm)
        common /ctranspv/               
     &                  ntrflx,cttrflx, npplas,ctpplas,
     &                  nq,ctq, nzefmd,ctzefmd,
     &                  nzeffp,ctzeffp, nsshaf,ctsshaf,
     &                  nbttot,ctbttot, ncur,ctcur,
     &                  nne,ctne, nplflx,ctplflx,
     &                  nrmzb,ctrmzb, nrmzc,ctrmzc,
     &                  ngfun,ctgfun,
     &                  nrmnmp,ctrmnmp,
     &                  nmomr,nxmomr,ctrmm,
     &                  nmomy,nxmomy,ctymm
        
        integer                 nvmax,imom
        parameter               (nvmax=20,imom=10)
        real 
     &                          cttrflx(nvmax),
     &                          ctpplas(nvmax),
     &                          ctq(nvmax),
     &                          ctzefmd(nvmax),
     &                          ctzeffp(nvmax),
     &                          ctsshaf(nvmax),
     &                          ctbttot(nvmax),
     &                          ctcur(nvmax),
     &                          ctne(nvmax),
     &                          ctplflx(nvmax),
     &                          ctrmzb(nvmax),
     &                          ctrmzc(nvmax),
     &                          ctgfun(nvmax),
     &                          ctrmnmp(nvmax),
     &                          ctrmm(0:imom,1:nvmax),
     &                          ctymm(0:imom,1:nvmax)

        integer
     &                          ntrflx,
     &                          npplas,
     &                          nq,
     &                          nzefmd,
     &                          nzeffp,
     &                          nsshaf,
     &                          nbttot,
     &                          ncur,
     &                          nne,
     &                          nplflx,
     &                          nrmzb,
     &                          nrmzc,
     &                          ngfun,
     &                          nrmnmp,
     &                          nmomr,
     &                          nxmomr,
     &                          nmomy,
     &                          nxmomy

C       PLASMA BOUNDARY - FOURIER MOMENTS DESCRIPTION (cm)
        common /ctranspb/                       
     &                          rm0,
     &                          rmm,
     &                          ymm

        real
     &                          rm0,
     &                          rmm(0:imom),
     &                          ymm(0:imom)

        integer                 nbndy



c=======            END OF TFTRDA.INC        =============================

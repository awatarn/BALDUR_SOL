c@commhd.m
c
c        COMMHD :: Equilibrium Interface Common Blocks
c        ---------------------------------------------
c     The following common blocks are available in all BALDUR routines
c  as well as in the interface, mapping, and flux surface averaging
c  routines in the equilibrium package.
c
c
      parameter  (Kjbal=mj,Knave=20,Kjflx=55,Kmhrm=7,Kixeq=65,Kjyeq=65
     &           ,kfour=32)
c
      common /comequ/
     &  leqxi , meqxi 
     &, eqrcj0(kjflx), eqrcjm(kjflx,kmhrm), eqrsjm(kjflx,kmhrm)
     &, eqycj0(kjflx), eqycjm(kjflx,kmhrm), eqysjm(kjflx,kmhrm)
c
cbate     &  leqxi , meqxi , eqxibi, eqxizi
cbate     &, eqrcj0(kjflx), eqrcjm(kjflx,kmhrm), eqrsjm(kjflx,kmhrm)
cbate     &, eqycj0(kjflx), eqycjm(kjflx,kmhrm), eqysjm(kjflx,kmhrm)
cbate     &, eqprzi, eqiotb
cbatec
cbate dimension  eqxibi(Kjflx),eqxizi(Kjflx),eqprzi(Kjflx),eqiotb(Kjflx)
c
      common /comMHD/
     >                          leqtyp , mjbal  , njav   , nnav   ,
     >                          mflxs  , mhrms  , mixeq  , mjyeq  ,
     >        r0ref  , b0ref  , pi     , twopi  , emu0   , eps0   ,
     >        torflx , qmin   , rhoedg , rbedge , eqcamp ,
     >        umib   , umid   , umie   , umii   , umij   , umil   ,
     >        umip   , umir   , umiv   ,
     >        avxiz  , avxib  , avi    , avti   , alint  ,
     >        vloopi , ajtori , dvoli  , vprimi , arms   , rgcs   ,
     >        areas  , vols   , ators  , ahalfs , rmids  , elong  ,
     >        triang , dent   , fltors , flpols , rbtors ,
     >        eqtmp  , ceqtmp , deqtmp , eeqtmp ,
     >        Requ   , Yequ   , A2dequ ,
     >        R0hr   , Rmhr   , Ymhr   , pspflx , pstflx , xsqrtf ,
     >        aprcri ,                   pspfle , pstfle ,
     >        eqerr(32)       , lsweq(32)       , cequil(32)      ,
     <        finmhd, rhoa
      dimension        avxiz(Kjbal),avxib(Kjbal),avi(Kjbal,Knave,5)
      dimension        avti(5)
      dimension        vloopi(Kjbal,2),ajtori(Kjbal,2),dvoli (Kjbal)
      dimension        vprimi(2,Kjbal),arms  (Kjbal,2),rgcs  (Kjbal,2)
      dimension        areas (Kjbal,2),vols  (Kjbal,2),ators (Kjbal,2)
      dimension        ahalfs(Kjbal,2),rmids (Kjbal,2),elong (Kjbal,2)
      dimension        triang(Kjbal,2),dent  (Kjbal,2)
      dimension        fltors(Kjbal,2),flpols(Kjbal,2),rbtors(Kjbal,2)
      dimension        eqtmp (Kjbal),  ceqtmp(Kjbal*5)
      dimension        deqtmp(Kjbal),  eeqtmp(Kjbal)
      dimension        Requ(Kixeq),Yequ(Kjyeq),A2dequ(Kixeq,Kjyeq)
      dimension        R0hr(Kjflx),Rmhr(Kmhrm,Kjflx),Ymhr(Kmhrm,Kjflx)
      dimension        pspflx(Kjflx),pstflx(Kjflx),xsqrtf(Kjflx)
      dimension        aprcri(Kjbal,2), rhoa(Kjflx)
c     ------------------------------------------------------------------
      common /comint/  cscoef(5*kjbal)
c
c  cscoef = coefficients used in IMSLMATH lib cubic spline routines
c     ------------------------------------------------------------------
      common /commht/  txtequ(10)
      character*80     txtequ
c     ------------------------------------------------------------------
      common /comeq1/
     &   rc0xbi, rcmxbi, rsmxbi, yc0xbi, ycmxbi, ysmxbi
     & , rthxbi, ythxbi, ntheta
c
      dimension
     &   rc0xbi(kjbal), rcmxbi(kmhrm,kjbal), rsmxbi(kmhrm,kjbal)
     & , yc0xbi(kjbal), ycmxbi(kmhrm,kjbal), ysmxbi(kmhrm,kjbal) 
     & , rthxbi(kfour,kjbal), ythxbi(kfour,kjbal)
c     ------------------------------------------------------------------
      common /comepr/
     &   epresm, eiotab
c
      dimension
     &   epresm(kjbal), eiotab(kjbal)
c
c**********************************************************************c
c
c  parameter list:
c  ---------------
cl       Kjbal ...BALDUR mesh                         >= njav = mjbal
cl       Knave ...dim. for # of flx.surface avergd quantities
cl       Kjflx ...dim. for # of flux surfaces                >= mflxs
cl       Kmhrm ...dim. for # of harmonics                    >= mhrms
cl       Kixeq ...x-dimension of equilibrium                 >= mixeq
cl       Kjyeq ...y-dimension of equilibrium                 >= mjyeq
c        Kfour ...number of Fourier harmonics (must be power of 2)
c
c  Variables in common block comequ:
c  --------------------------------
c
c leqxi  controls the spacing of the equilibrium grid eqxibi(j)
c                 leqxi is read in by the equilibrium namelist
c                 but, leqxi may be reset in sbrtn eqinit
c                 in order to conform to the equilibrium package
c                 being used.
c        = 0  --> equally spaced in sqrt ( toroidal flux )
c        = 1  --> equally spaced in toroidal flux
c
c meqxi  = number of equilibrium zone boundaries in grid eqxibi(j)
c
c eqxibi(j), j=1,mflxs ...normalized equilibrium flux label
c                         at equilibrium zone boundaries
c eqxizi(j)            ... same zone centered
c
c          Eqxibi(j) and eqxizi(j) define equilibrium grid boundaries
c       and centers.  They are physically equivalent
c       to the variables xbouni and xzoni in the rest of the BALDUR code.
c       They all stand for the square root of the toroidal flux normalized
c       to its value at the edge of the plasma.
c       As such, they are the fundamental independent variable against
c       which all other variables are defined.
c
c         For some equilibrium codes, eqxibi and eqxizi
c       should be equally spaced arrays.
c       For Hirshman's new equilibrium moments codes MHALF,
c       these arrays are spaced on equal intervals of toroidal flux
c       (that is, eqxibi(j)**2 are equally spaced).
c
c       Harmonic representation on the equilibrium grid:
c
c  eqrcj0(jx)      R(eqxibi(jx),theta) = eqrcj0(jx) + sum_{jm=1}^{mhrms}
c  eqrcjm(jx,jm)                 eqrcjm(jx,jm) * cos ( jm * theta )
c  eqrsjm(jx,jm)               + eqrsjm(jx,jm) * sin ( jm * theta )
c  eqycj0(jx)      Y(eqxibi(jx),theta) = eqycj0(jx) + sum_{jm=1}^{mhrms}
c  eqycjm(jx,jm)                 eqycjm(jx,jm) * cos ( jm * theta )
c  eqysjm(jx,jm)               + eqysjm(jx,jm) * sin ( jm * theta )
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  Variables in common block commhd:
c  ---------------------------------
c
c     COMMHD is the main interface common block between the equilibrium
c  package and the rest of BALDUR.  It is available in all the BALDUR
c  routines.  For more details on input variables, see the beginning
c  of file DEQBALD or see file DOCINPUT.
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  Variables in common block comeq1:
c  ---------------------------------
c
c     RC0XBI, RCMXBI, RSMXBI, YC0XBI, YCMXBI, YSMXBI represent the
c  harmonic representation of magnetic surfaces on BALDUR zone boundaries
c  computed in sbrtn EQRYTH in file DEQBALD.
c
c  R(xbouni(jx),theta) = rc0xbi(jx)
c          + sum_{jm=1}^{mhrms} [   rcmxbi(jx,jm) * cos ( jm * theta )
c                                 + rsmxbi(jx,jm) * sin ( jm * theta ) ]
c
c  Y(xbouni(jx),theta) = yc0xbi(jx)
c          + sum_{jm=1}^{mhrms} [   ycmxbi(jx,jm) * cos ( jm * theta )
c                                 + ysmxbi(jx,jm) * sin ( jm * theta ) ]
c
c  jx = 1,mjbal.
c
c     RTHXBI and YTHXBI represent magnetic surfaces on BALDUR zone
c  boundaries as a function of theta and xi:
c
c  R ( theta(jm), xbouni(jx) ) = rthxbi(jm,jx)    jm = 1 to ntheta
c  Y ( theta(jm), xbouni(jx) ) = ythxbi(jm,jx)    jx = 1 to mjbal
c
c  for theta(jm) = 2 * pi * real( jm - 1 ) / ( ntheta - 1 )
c
c
c  Variables in common block comepr
c  --------------------------------
c    Common block comepr is for plasma profiles on BALDUR grids
c    for use in the equilibrium packages.
c
c  epresm(j) = plasma pressure on baldur zone centers (avxiz(j))
c              mks units (nt/m^2)
c  eiotab(j) = iota on baldur zone boundaries (avxib(j))  dimensionless
c            = r0ref * bpoli(j) / ( b0ref * avi(j,1,1) )
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  Variables in common block comlhcd:
c  ---------------------------------
c
c    Common block comlhcd is for Lower Hybrid Current Drive and
c  Power deposition
c  Implemented by Bateman and Voitsekhovitch on 19 November 2000
c
c  xcd(jx) = normalized radial grid points, 0 < jx <= 55
c  tcd(jt) = time slices (sec), 0 < jt <= 20
c  cdrive(jx,jt) = current drive, A/m^2
c  plhcde(jx,jt) = power density deposited to electrons, W/m^3
c  cdcoef  = coefficient of cdrive (default 1.0)
c  plhcoef = coefficient of plhcde (default 1.0)
c  cdprof(jx) = profile of current drive, A/m^2
c               (computed in sbrtn solveb)
c  plhprof(jx) = profile of power density deposited to electrons, W/m^3
c  nxlhcd  = number of normalized radial grid points used
c  ntlhcd  = number of time slices used
c            (computed in sbrtn eqinit in file deqbald.f)
c
      common /comlhcd/
     &   xcd, tcd, cdrive, plhcde, cdcoef, plhcoef
     & , cdprof, plhprof
     & , nxlhcd, ntlhcd
c
      real  cdcoef, plhcoef, xcd(Kjbal), tcd(Knave)
     &    , cdrive(Kjbal,Knave), plhcde(Kjbal,Knave)
     &    , cdprof(Kjbal), plhprof(Kjbal)
   
c
      integer nxlhcd, ntlhcd
!
!--------1---------2---------3---------4---------5---------6---------7-c
!
! Common block for 2-D arrays from the equilibrium on BALDUR zone boundaries
!
cap
c bpol2di(jm,jz) = poloidal field at angle theta(jm) and zone boundary jz (T)
c btor2di(jm,jz) = toroidal field at angle theta(jm) and zone boundary jz (T)
c b22di(jz)      = <B^2>
c bm22di(jz)     = <1./B^22>
c delxi27b2i(jz) = < | del xi |^2 / B^2 >
c bdotdeltheta(jz) = < B \dot del theta / |B| >
!
      common /comb2d/ bpol2di, btor2di, b22di, bm22di
     &              , delxi27b2i, bdotdeltheta
c
      real ::         bpol2di( kjbal, kfour), btor2di( kjbal, kfour) 
     &              , b22di(kjbal), bm22di(kjbal)
     &              , delxi27b2i(kjbal), bdotdeltheta(kjbal)
!
!| {\tt delxib}( \xi, \theta ) = | \nabla \xi |
!|  = \sqrt{ ( \partial Y / \partial \theta )^2
!|        +  ( \partial R / \partial \theta )^2 } / 
!|  \left| ( \partial Y / \partial \theta ) ( \partial R / \partial \xi)
!|       - ( \partial Y / \partial \xi) ( \partial R / \partial \theta )
!|         \right|
!
!| {\tt delxthb}( \xi, \theta ) = | \nabla \theta |
!|  = \sqrt{ ( \partial Y / \partial \xi )^2
!|        +  ( \partial R / \partial \xi )^2 } / 
!|  \left| ( \partial Y / \partial \theta ) ( \partial R / \partial \xi)
!|       - ( \partial Y / \partial \xi) ( \partial R / \partial \theta )
!|         \right|
!
!| {\tt ejacob}( \xi, \theta ) = {\rm Jacobian} = R
!|  \left[ ( \partial Y / \partial \theta ) ( \partial R / \partial \xi)
!|       - ( \partial Y / \partial \xi) ( \partial R / \partial \theta )
!|         \right]
!
      common /come2d/ delxi2d, delth2d, ejacob2d
!
      real delxi2d(kjbal,kfour), delth2d(kjbal,kfour)
     &  , ejacob2d(kjbal,kfour)


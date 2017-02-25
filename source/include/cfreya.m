c@cfreya.m
CL                  C2.1     INPUT VARIABLES FOR FREYA
C     VERSION AL2     GGL/DEP     NOVEMBER 1976
       common/comloc/ yrmajr, yvols,
     r  yabeam, yiota , yangl1, yangl2, yangl3, yangl4, yphght, ypwdth,
     r  yamp  , ydvghz, ydvgvt, yhght , yhzfoc, ylengt, yvtfoc, ywdth ,
     r  clamda, cxmaxi, cxxcex, cxxel , cxximp, cxxion, yebeam, yfract,
     r  yllipt, yener , ymu   , ymub  , hrb   , yrlose, yqmean , yr    ,
     r  yrsq  , yrhob , yrhoe , yrhoi , yrhoim, yrmj  , yrmji , yrmjp,
     r  yrp   , yrpiv, yaumaj, yaumin, yte   , ytfrac, yti   , yzbeam,
     i  lybeam, mymp  , mnengy, mybeam, myzone, nyhght, nimp1 , nimp2 ,
     i  nymu  , nyprt , nyaper, nyshap, nybeam, nypart, nywdth, nyzone,
     i  mxysp, mysp, nyspec,
     l  mlyimp , nlybem
       logical
     l  mlyimp , nlybem
       dimension  yrmajr(26), yllipt(26), yvols(26),
     r   yiota(26),  yangl1(12), yangl2(12),
     r   yangl3(12), yangl4(12), yphght(12),
     r   ypwdth(12), yamp(12),  ydvghz(12), ydvgvt(12),
     r   yhght(12), yhzfoc(12), ylengt(12), yvtfoc(12),
     r   ywdth(12), clamda(12,26,2), cxmaxi(12), cxxcex(12,26,2),
     r   cxxel(12,26),  cxximp(12,26), cxxion(12,26),
     r   yebeam(12),yfract(12,3), yener(30), ymu(11),
     r   ymub(11), hrb(12,11,26,2),  yrlose(12,2), yqmean(4,26),
     r   yr(26),   yrsq(26), yrhob(26),yrhoe(26),yrhoi(26),
     r   yrhoim(4,26), yrmjp(12), yrpiv(12),
     r   yaumaj(12), yaumin(12), yte(26),
     r   ytfrac(12), yti(26), yabeam(12), yzbeam(12),
     i   lybeam(12,3), nyhght(12),nyaper(12), nyshap(12),
     i   nywdth(12), nyspec(12),
     l   nlybem(12)
c
C CLAMDA(12,26)  Effective total ionization cross section
C CXMAXI(12)     Max. CLAMDA over all JC zones.
C CXXCEX(12,26)  Charge exchange ionization cross section
C CXXEL(12,26)   Electron impact ionization cross section
C CXXIMP(12,26)  Impurity ionization cross section
C CXXION(12,26)  Ion impact ionization cross section
C HRB(12,11,26)  H(r) for energy JE, mu bin JMU, rasial zone JC.
C LYBEAM(12,3)   YENER index for beam JB, fraction JFRAC.
C MLYIMP        Fudge in FREYA to mock up impurity ionization of NB.
C MNENGY        Number of energy groups in HRB and YENER.
C MYBEAM        Max. number of beams allowed.
C MYMP          Number of impurity species.
C MYZONE        Max. number of plasma zones in HOFR.
C NIMP1         Does nothing?
C NIMP2         Does nothing?
C NLYBEM(12)     Should be .TRUE. for all beams passed to HOFR
C NYAPER(12)     Filled from NHAPER
C NYBEAM        Not set in FREDPS
C NYHGHT(12)     Filled from NHPROF
C NYMU          Set to NHMU
C NYPART        Set to NIPART
C NYPRT         Set to 5 (number of particles per injection direction)
C NYSHAP(12)     Filled from NHSHAP
C NYWDTH(12)     Filled from NHPRFV
C NYZONE        Number of JC zones, set in HOFR.
C YABEAM(12)     Atomic mass of the beam species
C YAMP(12)       Filled from HIBEAM (kA)
C YANGL1(12)     Filled from HANGLE(1
C YANGL2(12)     Filled from HANGLE(2
C YANGL3(12)     Filled from HANGLE(3
C YANGL4(12)     Filled from HANGLE(4
C YAUMAJ(12)     Optical depth to min. R
C YAUMIN(12)     Optical depth to min. r
C YDVGHZ(12)     Filled from HDIV
C YDVGVT(12)     Filled from HDIVV
C YEBEAM(12)     Filled from HEBEAM (eV)
C YENER(30)     Energies for H(r), sorted
C YFRACT(12,3)   Filled from HFRACT
C YHGHT(12)      Filled from HEIGHT
C YHZFOC(12)     Filled from HFOCL
C YIOTA(26)     Filled from 1/Q
C YLENGT(12)     Filled from HLENTH
C YLLIPT(26)    Elongation of flux surfaces on coarse FREYA grid
C YMU(11)       Filled from HMUC
C YMUB(11)      Filled from HMUB
C YPHGHT(12)     Filled from HAPERV
C YPWDTH(12)     Filled from HAPER
C YQMEAN(4,26)  Filled from CMEAN
C YR(26)        Inner boundary of zone JC (RMINS*XBOUNI)
C YRHOB(26)     Filled from RHOBIS
C YRHOE(26)     Filled from RHOELS
C YRHOI(26)     Filled from RHOHS
C YRHOIM(4,26)  Filled from RHOIS
C YRLOSE(12)     Calculated by FREYA, used to rescale HRB in FREDPS?
C YRMAJR(26)    Filled from SRMAJR.
C YRMJ          RMAJS
C YRMJI         Not set in FREDPS
C YRMJP(12)      Filled from HRMAJ
C YRP           RMINS
C YRPIV(12)      Filled from HRMIN
C YRSQ(26)      YR(JC) squared
C YTE(26)       Filled from TES
C YTFRAC(12)     Not set in FREDPS
C YTI(26)       Not set in FREDPS
C YVOLS(26)     Volume of flux surfaces on coarse FREYA grid [cm**3]
C YVTFOC(12)     Filled from HFOCLV
C YWDTH(12)      Filled from HWIDTH
C YZBEAM(12)     Charge of the beam species
c**********************************************************************c

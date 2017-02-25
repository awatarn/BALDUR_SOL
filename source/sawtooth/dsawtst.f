!| %dsawtst.f
!| %
!| %  Test to determine if sawtooth crash should occur this time-step.
!| %
!| %  To convert this FORTRAN source file into a LaTeX file, type
!| % python ./s2tex.py dsawtst.f
!| %
!| %  To typeset the resulting LaTeX file, type
!| % latex dsawtst.tex
!| % latex dsawtst.tex
!| 
!| \documentstyle {article}    % Specifies the document style.
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \title{SAWTST: Test for Sawtooth Crash in BALDUR}
!| \author{Glenn Bateman and the authors of BALDUR \\
!|         Lehigh University, Bethlehem, PA 18015}
!| 
!| \begin{document}           % end of preamble and beginning of text.
!| 
!| \maketitle                 % Produces the title.
!| 
c/ 20.09 16:50 22-dec-2004 
c BALDUR  file dsawtst.f by Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@sawtst  .../baldur/code/bald/dsawtst.tex
c  cnn 5-aug-02 added the porcelli model
c  pis 23-jun-98 changed swton to swton(1)
c  rgb 07-jul-97 reset tbsaws=-1./epslon if ztime .lt. swton
c  rgb 20.27 23-feb-92 csawth(9) multiplier for Monticello sawtooth period
c    if ( csawth(9) .lt. epslon ) csawth(9) = 1.0
c  rgb 20.09 16:50 25-aug-91 added Rogers-Zakharov trigger
c    documented with LaTeX
c    Compute sawr1 by interpolating ahalfs(jz,1) directly
c  rgb 18.78 21:00 08-nov-90 base Park-Monticello on tepost(lcentr)
c      when available and on tes(2,lcentr)*useh otherwise
c  rgb 18.66 14:00 20-sep-90 Park-Monticello sawtooth period lsawth(10)=2
c  rgb 11:00 21-dec-89 cleaned up sbrtn
c
c
c               %%%%%%%%%%%%%
c               %           %
c               %  sawtst   %
c               %           %
c               %%%%%%%%%%%%%
c
c
        subroutine sawtst (iswtst)
c
c       Test to determine if sawtooth crash should occur this time-step.
c
c       Usage:
c               CALL SAWTST (ISWTST)
c
c  ISWTST .LE. 0 ==> NO SAWTOOTH CRASH THIS STEP
c  ISWTST .GE. 1 ==>    SAWTOOTH CRASH THIS STEP
c
c  Input:
c
c  cfutz(460) or swton  = time of first sawtooth crash
c                           provided the other criteria are satisfied
c  cfutz(461) or swtoff = time sawteeth turned off
c  cfutz(462) or swperd = sawtooth period when other criteria are off
c  cfutz(463) = time at which prescribed sawtooth period changes
c                 from cfutz(462) to cfutz(464)
c  cfutz(464) = sawtooth period after cfutz(463)
c  cfutz(465) = sawtooth period given by    \tau =
c  cfutz(466)     cfutz(465) \tau_r^{cfutz(466} \tau_A^{cfutz(467)}
c  cfutz(467)     \tau_h^{1.-cfutz(466)-cfutz(466)}
c
c  cfutz(477) or swqmin = minimum q-value for sawtooth crash
c
c  lsawth(10) = 0 use sawtooth models given above
c             = 2 use Park-Monticello model
c             = 3 use Rogers-Zakharov model
c             = 4 print out diagnostic from Rogers-Zakharov model
c             = 5 use Porcelli model
c             = 6 print out diagnostic from Porcelli model
c
c  lsawth(27) = 0 to skip repeated calls to Porcelli model
c
c  csawth(9)  = multiplier for the Monticello-Park sawtooth period (1.0)
c
c  Rogers-Zakharov model:
c
c  csawth(10) coef of $s_{CD}$ critical shear for collisionless
c               tearing modes (1.0)
c  csawth(11) coef of $s_0$, critical shear for one fluid MHD mode (1.0)
c  csawth(12) coef of $s_H$, critical shear finite Larmor radius
c               stabilized ideal MHD mode (1.0)
c  csawth(13) coef of nonthermal contribution to ion pressure (1.0)
c  csawth(14) lower bound on shear for sawtooth crash (0.03)
c  csawth(15) lower bound on $\Lambda_H$ for MHD instability (0.0)
c
c
c**********************************************************************c
c
c
      use porcelli_module        
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cfokkr.m'
      include 'calpht.m'
      include 'commhd.m'
      include 'clsaw.m'
c
        common/comsaw/ tbsaws,sawte,sawne,sawdd,sawate,sawane,sawadd,
     1  sawtau,sawr1,sawri,sawrmx,sawqc,tmguir,njzq1,
     2  sawfus,sawalf,sawuth,estpol,uthavg,tauein,ohptot,
     3  tepre(55),tepost(55),teavg(55),tipre(55),tipost(55),tiavg(55),
     4  ajzpre(55),ajzpst(55),ajzavg(55),ohravg(55)
c
      dimension ztemp1(55), ztemp2(55)
c
      integer iswton, iporc
c
      real ztporck
c
c  iswton is the current index of swton(j)
c  kswmax is the maximum allowed value of iswton
c
      logical  inital, isporc
      data inital /.true./, isporc /.true./
c
      data  iclass /2/,  isub /24/, iswton /1/, iporc /0/
c
      data ztporck / -1.0e-10/
c
      save ztporck
      save iswton, iporc
c
       integer::intr1,izlsaw
       integer,save::kcnnsaw
       real,save::tpold,tpsaw
       real::zrq1,zelongr1,zdenel0,zmui,zbtor1,zefast
       real::zti0,zbpol11,zbetapf1,zntime,zrq2
       real::zbeampr,zzpr,zimfast,zshearp,ztvols
       real::zfa,zrplen1,zgradpf1,zrda
       real::zbetaim,zpresr1,zp1avg,zcpf,zlir1,zzdr
       real::zsub11,zsum,zsum11,zarea1  
       real,dimension(55)::zrdp1,zttpi,zpsum,zzfb
       real,dimension(55)::zztrapp,zztrappg,zp1avgs
       real,dimension(55)::zdvols,zbsum
       real,dimension(55,2)::zlahi
c
       real,dimension(55) :: zbetai(55), zbetae(55), zbetat(55)
     &  , zrminor1(55), zrminor2(55), zvolume(55), zarea(55)
     &  , zbpol2di(55)
c
c..Set ISWTST indicator
c  ----------
c  ISWTST .LE. 0 ==> NO SAWTOOTH CRASH THIS STEP
c  ISWTST .GE. 1 ==> SAWTOOTH CRASH THIS STEP
c
        iswtst = 0
c
c..to ensure backward compatibility
c
cpis changed swton to swton(1)
        if (cfutz(460) .gt. epslon) swton(1) = cfutz(460)
        if (cfutz(461) .gt. epslon) swtoff   = cfutz(461)
        if (cfutz(462) .gt. epslon) swperd   = cfutz(462)
        if (cfutz(477) .gt. epslon) swqmin   = cfutz(477)
c
c..initialize
c  ----------
        if(inital) then
          inital = .false.
          write (6,*) ' sbrtn SAWMIX sawtooth mixing by Bateman used'
          njzq1=0
          tmguir=0.
          tbsaws= - 1. / epslon
          sawtau=0.
          ztau=1.
cbate960613          call sawavg(2,teavg,tiavg,ajzavg,ohravg,ohptot,ztau)
        end if
c
        if(nstep.eq.0) return
c
c..test sawtooth times on and off
c
          ztime=tai*uist
        if ( ztime .lt. swton(iswton) ) then
          tbsaws = -1.0 / epslon
          return
        endif
c
        jmin = iswton
        do j=jmin, kswmax-1
          if ( ztime .ge. swton(j+1) ) then
            tbsaws = - 1.0 / epslon
            iswton = iswton + 1
          else
            go to 12
          endif
        enddo
  12    continue
c
          ztoff=1.e10
        if ( swtoff .gt. 0. ) ztoff= swtoff
        if ( ztime  .gt. ztoff ) then
          tbsaws = -1.0 / epslon
          return
        endif
c
c..iwxmin = index of the first zone with xzoni > swxmin
c
      do 16 jx=lcentr,mzones
        iwxmin = jx
        if ( xzoni(jx) .gt. swxmin ) go to 18
  16  continue
  18  continue
c
c..compute zqmin = minimum q-value for jx=iwxmin,mzones
c
        zqmin = q(iwxmin)
      do 19 jx=iwxmin+1,mzones
        zqmin = min ( zqmin, q(jx) )
  19  continue
c
c..find the outermost q=1 surface
c  do not consider q=1 surfaces with xzoni < swxmin
c
          iq1=0
        do 50 jz=iwxmin,ledge
          if( q(jz) .le. 1.0 ) iq1=jz
50      continue
          njzq1 = iq1
          intr1=iq1
c
        if(iq1.gt.lcentr) then
          zf   = 0.
          if( q(iq1) .ne. q(iq1+1) ) zf=(1.-q(iq1))/(q(iq1+1)-q(iq1))
          zxq1 = xbouni(iq1) + zf*dxzoni(iq1)
          sawr1 = ahalfs(iq1,1) + zf * (ahalfs(iq1+1,1)-ahalfs(iq1,1))
        else
          zxq1=0.
        end if
c
cbate   sawr1 = zxq1 * rmins
        sawqc = q(lcentr)
c
!| 
!| Return without resetting the sawtooth mixing trigger if there is no
!| $q=1$ magnetic surface within the plasma or if the sawteeth are not
!| turned on
!| 
c
        if ( iq1 .le. lcentr ) return
        if ( iq1 .ge. ledge  ) return
c
!| 
!| Reset the sawtooth mixing trigger and return if the minimum value of $q$
!| is less than or equal to {\tt swqmin}.
!| 
c
        if ( zqmin .le. swqmin )  then
          iswtst = 1
          return
        endif
c
c..calculate tmguir if this is the step after a sawtooth crash
c
        if( cfutz(465) .le. epslon ) go to 45
        if( njzq1 .lt. lcentr ) go to 45
          ztaur=4.*fcpi*(sawr1/(fcc*10.**fxc))**2/eta(2,njzq1)
          zrho=0.
        do 10 jh=1,lhydn
          zrho=zrho+rhohs(jh,2,lcentr)*aspec(jh)
10      continue
        if ( mimp .gt. 0 ) then
          do 20 ji=1,mimp
            ja=limp1+ji-1
            zrho=zrho+rhois(ji,2,lcentr)*aspec(ja)
  20      continue
        endif
c
          zrho=zrho*fcmp*10.**fxnucl
          ztaua=2.*fcpi*rmajs*sqrt(4.*fcpi*zrho)/bzs
          zuthe=0.
          zheat=0.
          zcool=0.
        do 40 jz=lcentr,njzq1
          zuthe=zuthe+rhoels(2,jz)*tes(2,jz)*dx2i(jz)
          zheat=zheat+dx2i(jz)*(weohms(jz)+weauxs(jz)+webems(jz)+
     1      wealfs(jz)+weecrh(jz)+weicrf(jz))
          zeirs=0.
          if ( mimp .gt. 0 ) then
            do 35 ji=1,mimp
              zeirs=zeirs+weirs(ji,jz)
35          continue
          endif
          zcool=zcool+dx2i(jz)*(webrs(jz)+weions(jz)+wesrs(jz)+zeirs)
          zei  = dx2i(jz)*cnueqs(jz)*(tes(2,jz)-tis(2,jz))
          if( zei .lt. 0.) zheat=zheat-zei
          if( zei .gt. 0.) zcool=zcool+zei
40      continue
          ztauh  = 1.5*zuthe/zheat
          zcons  = -1.
        if ( cfutz(465) .gt. 0. ) zcons=cfutz(465)
          zexp1  = 0.42
        if ( cfutz(466) .ne. 0. ) zexp1=cfutz(466)
          zexp2  = 0.14
        if ( cfutz(467) .ne. 0. ) zexp2=cfutz(467)
          zexp3  = 1.-zexp1-zexp2
          tmguir = zcons*(ztaur**zexp1)*(ztaua**zexp2)*(ztauh**zexp3)
          njzq1  = 0
45      continue
c
c..Park-Monticello sawtooth period  PPPL-2601 (1990)  when lsawth(10) = 2
c  \tau_{\mbox{sawtooth}}
c    = 0.009 R_m^2 T_{\mbox{e0 keV}}^{3/2} / Z_{\mbox{0 eff}} \;\; \mbox{sec}
c
      if ( csawth(9) .lt. epslon ) csawth(9) = 1.0
c
      if ( lsawth(10) .eq. 2 ) then
        zte0 = tes(2,lcentr)*useh
        if ( tepost(lcentr) .gt. epslon ) zte0 = tepost(lcentr)
        ztpark = csawth(9) * 0.009 * ( rmids(lcentr,2)*usil)**2
     &    * ( zte0 )**1.5 / xzeff(2,lcentr)
      else
        ztpark = 0.0
      endif
!| 
!| \subsection{Rogers-Zakharov Model}
!| 
!| The following section implements a sawtooth trigger model by Rogers
!| and Zakharov\cite{roge91a,roge91b} when ${\tt lsawth(10)} = 3$.
!| Diagnostic output from this model is presented when
!| ${\tt lsawth(10)} = 4$.
!| 
!| First, compute physical constants $m_e = {\tt zme}$, $m_p = {\tt zmp}$,
!| average ion mass $m_i = {\tt zmi}$ (including impurities,
!| $e = {\tt zee}$, and $c = {\tt zcc}$.
!| Make sure the shear and gradient scale lengths are computed
!| and establish the position of the $q=1$ and $q=2$ surfaces
!| $$ r_1 = r_{q=1}               \eqno{\tt zr1} $$
!| $$ r_2 = r_{q=2}               \eqno{\tt zr2} $$
c
      if ( lsawth(10) .eq. 3 .or. lsawth(10) .eq. 4 ) then
c
        if ( nstep .lt. 2 ) return
c
        zme = fcme * 10.**fxme
        zmp = fcmp * 10.**fxnucl
        zmi = aimass(1,iq1) * zmp
        zee = fces * 10.**fxes
        zcc = fcc * 10.**fxc
c
        call xscale
c
        iq1 = njzq1
        zr1 = sawr1
c
        zbtors = rbtors(iq1,1) / rmids(iq1,1)
c
c  find the q=2 surface
c
      do 61 jz=iq1+1,ledge
        iq2 = jz - 1
        if ( q(jz) .gt. 2.0 ) go to 62
  61  continue
  62  continue
c
      zint = 0.
      if ( q(iq2+1) .gt. q(iq2) ) zint = (2.0-q(iq2))/(q(iq2+1)-q(iq2))
      zr2  = ahalfs(iq2,1) + zint * (ahalfs(iq2+1,1)-ahalfs(iq2,1))
c
!| 
!| Interpolate variables to the $q=1$ and $q=2$ surface as needed
!| $$ j_1 = j(r_1)                    \eqno{\tt zj1} $$
!| $$ j_2 = j(r_2)                    \eqno{\tt zj2} $$
!| $$ j_1'' = j''(r_1)                \eqno{\tt zjpp1} $$
!| Note, the dimensionless form of current is given by
!| $$ j = \frac{4\pi}{c} \frac{J(r) R}{B}  \eqno{\tt zj} $$
!| $$ B_{\theta 1} = B_\theta (r_1)   \eqno{\tt zbpol1} $$
c
        ztemp1(1) = 0.0
        ztemp1(2) = zr1
        ztemp1(3) = zr2
c
      do 63 jz=1,mzones
        eqtmp(jz) = 4.0 * fcpi * ajzs(2,jz) * rmids(jz,2)
     &    / ( zcc * zbtors )
  63  continue
c
      call cubint (ahalfs(1,2),eqtmp,mzones-1,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,3,0,0.0,1
     & ,'abort computing current at r1 and r2 in sbrtn sawtst')
c
        zj1 = ztemp2(2)
        zj2 = ztemp2(3)
c
      call cubint (ahalfs(1,2),eqtmp,mzones-1,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,2,0.0,1
     & ,'abort computing d2 j / dr2 at r1 in sbrtn sawtst')
c
        zjpp1 = ztemp2(2)
c
      do 64 jz=1,mzones
        eqtmp(jz) = bpols(1,jz)
  64  continue
c
      call cubint (ahalfs(1,1),eqtmp,mzones,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,0,0.0,-1
     & ,'abort computing bpol at r1 in sbrtn sawtst')
c
        zbpol1 = ztemp2(2)
c
!| The total pressure includes the thermal pressure and
!| {\tt csawth(13)} times the nonthermal contributions
!| $$ p_1 = p(r_1)                    \eqno{\tt zpres1} $$
!| $$ \beta_1 = 8 \pi p_1 / B^2       \eqno{\tt zbeta1} $$
c
      do 65 jz=1,mzones
        eqtmp(jz) = rhoels(2,jz)*tes(2,jz) + rhoins(2,jz)*tis(2,jz)
     &   + csawth(13) * ( alphai(jz)*ealfai(jz)*uisd*uise
     &     + hebems(jz)*rhobis(2,jz)/3. )
  65  continue
c
      call cubint (ahalfs(1,2),eqtmp,mzones,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,0,0.0,1
     & ,'abort computing pressure at r1 in sbrtn sawtst')
c
        zpres1 = ztemp2(2)
        zbeta1 = 8.0 * fcpi * zpres1 / zbtors**2
c
!| $$ \beta_p = \frac{16 \pi \int_0^{r_1} [ p(r) - p(r_1) ] r \; dr
!|     }{r_1^2 B_\theta^2}        \eqno{\tt zbetap} $$
c
      do 68 jz=1,mzones
        eqtmp(jz) = ( eqtmp(jz) - zpres1 ) * ahalfs(jz,2)
  68  continue
c
      call cubint (ahalfs(1,2),eqtmp,mzones,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,3,0.0,-1
     & ,'abort computing \int (p-pr1) r in sbrtn sawtst')
c
      zbetap = 16. * fcpi * ztemp2(2) / ( zr1**2 * zbpol1**2 )
c
!| Gradient of the ion pressure
!| $$ p_i' = p_i'(r_1)                \eqno{\tt zdpi} $$
c
      do 66 jz=1,mzones
        eqtmp(jz) = rhoins(2,jz)*tis(2,jz)
     &   + csawth(13) * ( alphai(jz)*ealfai(jz)*uisd*uise
     &     + hebems(jz)*rhobis(2,jz)/3. )
  66  continue
c
      call cubint (ahalfs(1,2),eqtmp,mzones,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,1,0.0,1
     & ,'abort computing d p_i / dr at r1 in sbrtn sawtst')
c
        zdpi = ztemp2(2)
c
!| $$ b = \frac{4}{r_1^2} \int_0^{r_1} r^3
!|      \left( \frac{1}{q} - 1 \right) dr \eqno{\tt zb} $$
c
      do 67 jz=1,mzones
        eqtmp(jz) = ahalfs(jz,1)**3 * ((1/q(jz)) - 1.0)
  67  continue
c
      call cubint (ahalfs(1,1),eqtmp,mzones,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,3,0.0,-1
     &,'abort from call cubint to compute zb in sbrtn sawtst')
c
        zb = ztemp2(2) * 4.0 / zr1**4
c
!| 
!| $$ \Lambda_H = \frac{r_1^2}{R^2} \left( \frac{13}{64} b
!|     - \frac{c}{4(1-c)} \beta_p^2 \right) \eqno{\tt zlamdh} $$
!| where
!| $$ c = 1 - j_1/2 - \alpha / (1-\alpha)       \eqno{\tt zc} $$
!| $$ \alpha = \lambda^2 \; \frac{2 \lambda \lambda_2 (2-j_1)
!|    -(1-\lambda_2)(1+2\lambda_2)(j_1-j_2)}{(1-\lambda_2) \lambda_2^2
!|    (j_1-j_2) + 2 \lambda \lambda_2 (2-j_1)}  \eqno{\tt zalpha} $$
!| $$ \lambda = r_1^2 / r_2^2                   \eqno{\tt zlam} $$
!| $$ \lambda_2 = \frac{1 - j_2 - 2\lambda
!|    + j_1 \lambda}{j_1 - j_2}                 \eqno{\tt zlam2} $$
c
        zlam   = zr1**2 / zr2**2
        zlam2  = (1. - zj2 - 2.*zlam + zj1*zlam) / (zj1-zj2)
        zalpha = zlam**2 * ( 2.*zlam*zlam2*(2.-zj1)
     &    - (1.-zlam2)*(1.+2.*zlam2)*(zj1-zj2) )
     &    / ( (1.-zlam2) * zlam2**2 * (zj1-zj2)
     &      + 2. * zlam * zlam2 * (2.-zj1) )
        zc     = 1. - zj1/2. - zalpha/(1.-zalpha)
c
        zlamdh = ( zr1**2/rmids(iq1,2)**2 ) * ( 13. * zb / 64.
     &    - zc * zbetap**2 / ( 4. * ( 1. - zc ) ) )
c
!| 
!| Consider the one fluid MHD instability of the
!| $m=1$ mode at $r=r_{q=1}$:
!| If $ J''(r_1) > 0 $ or $ \Lambda_H < 0 $ then this mode is unstable
!| and we set set $ S_0 = 0.$
!| Otherwise, a lower bound on the shear is given by
!| $$ S_0 = \left(\frac{\pi \sqrt{6} \Lambda_H \sqrt{-j_1'' r_1^2}}{4}
!|          \right)^{0.4}   \eqno{\tt zs0} $$
!| 
c
        zs0 = 0.0
c
      if ( zjpp1 .lt. 0.0 .and. zlamdh .gt. 0.0 ) then
        zs0 = csawth(11) * ( fcpi * zlamdh
     &     * sqrt( - 6.0 * zjpp1 * zr1**2 ) / 4. )**0.4
      endif
c
!| 
!| Now consider the finite Larmor radius stabilized ideal mode.
!| This mode is unstable if the shear is less than
!| $$ s_H = - \Lambda_H \frac{\omega_{pi} r_1}{c}
!|  \left|  \frac{B^2}{4 p_i' R} \right| \eqno{\tt zsh} $$
!| where the ion plasma frequency is given by
!| $$ \omega_{pi} = \sqrt{4 \pi \sum_{\rm ions} n_i Z_i^2 e^2 / m_i}
!|                                     \eqno{\tt zplasi} $$
c
        zsh = 0.0
c
        zbetpi = - 4. * fcpi * zdpi * rmids(iq1,1) / zbtors**2
c
      if ( zlamdh .lt. csawth(15) .and. csawth(12) .gt. 0.0 ) then
c
          znzm = 0.0
        do 71 js=1,mhyd
          znzm = znzm + rhohs(js,1,iq1) / ( aspec(js) * zmp )
  71    continue
c
        if ( mimp .gt. 0 ) then
          do 72 js=1,mimp
            znzm = znzm + rhois(js,1,iq1) * c2mean(js,1,iq1)
     &        / ( aspec(mhyd+js) * zmp )
  72      continue
        endif
c
        zplasi = sqrt ( 4. * fcpi * zee**2 * znzm )
c
        zsh = - csawth(12) * fcpi * zlamdh * zplasi * zr1
     &        / ( zcc * max( abs(zbetpi), epslon ) )
c
      endif
c
!| 
!| Finally, consider the $m=1$ tearing mode stabilized by finite Larmor radius
!| effects (collisionless tearing modes)\cite{copp91a}.
!| the critical shear is
!| $$ S_{CD} = \sqrt{\frac{m_i}{m_e}} \left| \frac{2 \pi p_i' R}{B^2}
!|    \right| \frac{ {\tt csawth(10)} }{1 + 0.4 \sqrt{\beta m_i / m_e} }
!|                                   \eqno{\tt zscd} $$
!| where $p_i'$ ({\tt zdpi})
!| is the gradient of the thermal ion pressure at $r_1$,
!| $\beta$ is the total beta including fast ions,
!| and $m_i$ ({\tt zmi}) is the density weighted ion mass.
!| The default value of ${\tt csawth(10)}$ should be 1.0.
c
        zscd = csawth(10) * sqrt(zmi/zme)
     &    * abs( 2. * fcpi * zdpi * rmids(iq1,1) / zbtors**2 )
     &    / ( 1. + 0.4 * sqrt( zbeta1 * zmi/zme ) )
c
!| 
!| Now compare the shear at the $q=1$ surface with the critical shear.
!| The condition for a sawtooth crash is
!| $$ s_H > {\tt shear} > s_0 + s_{CD} $$
c
      call cubint (ahalfs(1,2),shear,mzones-1,0,ceqtmp,kjbal
     & ,ztemp1,ztemp2,2,0,0.0,1
     & ,'abort computing shear at q=1 in sbrtn sawtst')
c
      zshear = ztemp2(2)
c
c  diagnostic output
c
      write (6,160) nstep,tai,zs0,zscd,zsh,zshear,zlamdh
 160  format(' nstep=',i3,'  t=',0pf9.6,'  zs0=',1pe11.4
     &  ,'  zscd=',1pe11.4,'  zsh=',1pe11.4,'  shear=',1pe11.4
     &  ,'  zlamdh=',1pe11.4)
c
      write (6,161) nstep,tai,zr1,zjpp1*zr1**2,zbetpi
 161  format(' nstep=',i3,'  t=',0pf9.6,'  zr1=',1pe11.4
     &  ,'  r_1^2 d2j=',1pe11.4,'  d beta_i=',1pe11.4)
c
      if ( lsawth(10) .eq. 4 ) go to 79
c
        if ( zshear .gt. zs0 + zscd  .and. zshear .gt. csawth(14) ) then
          iswtst = 1
          write(6,*)
          write(6,*) ' Rogers-Zakharov criterion for sawtooth crash'
          write(6,*) ' shear(iq1) .gt. zs0 + zscd'
cbate     if ( tai .gt. csawth(16) ) stop
        elseif ( zshear .lt. zsh ) then
          iswtst = 1
          write(6,*)
          write(6,*) ' Rogers-Zakharov criterion for sawtooth crash'
          write(6,*) ' shear(iq1) .lt. zsh'
cbate     if ( tai .gt. csawth(16) ) stop
        endif
c
!| 
!| Return after testing the Rogers-Zakharov cirterion to close out this
!| ${\tt lsawth(10)} = 3$ option.
c
          return
      endif
c
  79  continue
c
c
!| \subsection{Porcelli Model}
!| 
!| The following section implements a sawtooth trigger model by Porcelli
!| \cite{porcelli96a} when ${\tt lsawth(10)} = 5$.
!| Here the inputs needed for the porcelli model are calculated and
!| send to the sawmodule where the actual test for whether a sawtooth
!| occured for this time step takes place.
c
      if ( lsawth(10) .eq. 5 .or. lsawth(10) .eq. 6
     &   .or. lsawth(10) .eq. 7 ) then ! porcelli model option 
         if(isporc) then  ! open file for writing diagnostic outputs
            isporc = .false.
            open (50,file='cnnsaw.out', status='unknown') 
            open (60,file='cnnsaw.dat', status='unknown') 
            open (70,file='cnntst.dat', status='unknown')
            open (80,file='cnndia.dat', status='unknown')
            open (90,file='cnntst1r.out',status='unknown')
            open (100,file='cnnper.dat', status='unknown')
            write(50,*)"Inputs and output of the Porcelli model"
            write(60,665)
            write(70,666)
            write(80,667)
            write(90,668)
            write(100,669)
 665        format("# time",8x,"delw_bu",6x,"delw_el",4x,"delw_ko",4x,
     &         "delw_fast",6x,"delw",/,
     &  "#----------------------------------------------------------")
 666         format("# time",8x,"LHS13",6x,"RHS13",6x,"LHS14",6x,
     &         "RHS14",6x,"LHS15",6x,"omeg_di",6x,"cgrowth",/,
     &  "#----------------------------------------------------------")
 667         format("# time",8x,"q-axis",6x,"shear1",4x,"lir1", /,
     &  "#----------------------------------------------------------")  
 668         format("# time",8x,"cr/no crash",6x,"period (sec)", /,
     &  "#----------------------------------------------------------") 
 669         format("# time",8x,"period (sec)", /,
     &  "#----------------------------------------------------------")    
c
             tpold = tai*uist ! sawtooth period reference time  
             tpsaw = 0.0
       endif
c
c         if(ztime .ge. ztporck) then  !check to decide whether to call the porcelli model  
c
         zrq1 = sawr1*usil    ! radius of q=1 surface (m)
         zdenel0 = rhoels(1,2)*1.0e-14 ! central electron density (10^20) m^-3
         zti0 = tis(1,2)*useh    ! central ion temperature (keV)
         zbtor1 = bzs*usib        ! vacuum toroidal field at zrmaj (tesla)
         zntime = tai*uist       ! time (sec)
c
         zbetaim=0.0
         zzpr=0.0 
         do 150 jz = 1,mzones
            if(betati(jz) > zbetaim) zbetaim = betati(jz)  ! maximum thermal ion beta
            if(totprs(jz) > zzpr)  zzpr = totprs(jz) ! maximun pressure 
 150     continue
c
         ztemp1(1) = 0.0
         ztemp1(2) = zrq1
         ztemp1(3) = ahalfs(mzones-3,2)*usil
c
         zp1avgs = 0.0
         zlahi = 0.0
         zttpi = 0.0
         zrdp1 = 0.0
         zztrapp = 0.0
         zztrappg = 0.0
         zbsum = 0.0
         zdvols = 0.0
c
         do 145 jz=1,mzones-1
            zlahi(jz,2)=ahalfs(jz,2)*usil
            zttpi(jz)=0.1*totprs(jz)   ! total pressure nt/m**2
            zrdp1(jz)=usil*armins(jz,1)/max(abs(slprs(jz)),epslon) ! thermal p/|dp/dr|  (m)
            zp1avgs(jz+1) = zp1avgs(jz)
     &           + totprs(jz+1)*(vols(jz+1,1)-vols(jz,1))
            zrda = areas(jz+1,1)-areas(jz,1)          ! da 
            zbsum(jz+1)=zbsum(jz)+ 0.5*(bpol2di(jz+1,1)**2
     &           + bpol2di(jz,1)**2)*zrda 
            zzdr = usil*(armins(jz+1,1)-armins(jz,1))/zrq1
            zfa = (usil*armins(jz+1,1)/zrq1)**1.5
            zztrapp(jz+1)=zztrapp(jz)+zfa*(thrprs(jz+1)/zzpr)*zzdr
            zgradpf1=slpts(jz+1)                           
     &           +slprs(jz+1)*thrprs(jz+1)/max(epslon,armins(jz+1,1)) ! fast_gradp=total-thermal     
            zztrappg(jz+1)=zztrappg(jz)
     &           + zfa*0.1*zgradpf1*zzdr*zrq1/usil ! mks (covert factor=0.1)
            zdvols(jz+1) = zdvols(jz)+(vols(jz+1,1)-vols(jz,1))
 145     continue   
c
c  Interpolate to the q=1 surface
c 
       call cubint (zlahi(1,2),zttpi,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &   ,'abort computing total pressure at q=1 in sbrtn sawtst')
c
         zpresr1 = ztemp2(2) ! total pressure nt/m**2 at q=1 surface
c
          call cubint (zlahi(1,2),zdvols,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing volume at q=1 in sbrtn sawtst') 
c
          ztvols = ztemp2(2)
c          
          call cubint (zlahi(1,2),zp1avgs,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing sum vol*press at q=1 in sbrtn sawtst') 
c
          zsum = ztemp2(2)
          zp1avg =  0.1*zsum/ztvols  ! average total pressure up to q=1 surface (J/m^3)
c
         call cubint (zlahi(1,2),elong,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &    ,'abort computing elongation at q=1 in sbrtn sawtst')
c
         zelongr1 = ztemp2(2) ! elongation of q=1 surface
c
         call cubint (zlahi(1,2),bpol2di,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &     ,'abort computing B-poloidal at q=1 in sbrtn sawtst')
C
         zbpol11 = ztemp2(2) ! poloidal field (Tesla) at q=1 surface (midplane outboard)
c
         call cubint (zlahi(1,2),shear,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing shear at q=1 in sbrtn sawtst')
c
         zshearp = abs(ztemp2(2))
         if (zshearp < epslon) zshearp = 2.0*(1.0-qsmth(2)) ! default
c
          call cubint (zlahi(1,2),zrdp1,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing press. length at q=1 in sbrtn sawtst')        
c
          zrplen1=abs(ztemp2(2))
          if (zrplen1 < epslon) zrplen1 = zrq1 ! default
c
          call cubint (zlahi(1,2),zbsum,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing sum bpol**2da at q=1 in sbrtn sawtst') 
c
          zsum = ztemp2(2)
c
          call cubint (zlahi(1,2),areas,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing sum area at q=1 in sbrtn sawtst') 
c
          zarea1 = ztemp2(2)
c
         if(klir1 == 1)then
           zlir1 = 0.5 + (1.0 - qsmth(2))/3.0   ! internal inductance at q=1 surface for                                                                      ! parabolic q-profile and small 1-q0
           zshearp = 6.0*(zlir1 - 0.5)
         else
            zlir1=zsum/(zarea1*zbpol11**2)
         endif
c
!| 
!| \[zbetapf1 =\frac{-2\mu_{0}}{B_{\rm p}^{2}}\int_{0}^{1}x^{3/2}\frac{dp_{\rm fast}}{dx} dx\]
!| where ($x=r/r_{1}$) 
!| \[zcpf=\frac{(5/2)\int_{0}^{1}x^{3/2}P_{\rm thermal}(x) dx}{P_{\rm i0}}\]
!| 
c
          call cubint (zlahi(1,2),zztrapp,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing sum integ press at q=1 in sbrtn sawtst') 
c
          zsum = ztemp2(2)
          zcpf = (5/2)*zsum
c 
          call cubint (zlahi(1,2),zztrappg,mzones-1,0,ceqtmp,kjbal
     &        ,ztemp1,ztemp2,2,0,0.0,1
     &        ,'abort computing sum integ press at q=1 in sbrtn sawtst') 
c
          zsum11 = ztemp2(2)
          zbetapf1 = -2.0*emu0*zsum11/zbpol11**2 ! poloidal beta of fast particle
c
         zsum = 0.0
         zsum11 = 0.0
         do 300 jz = lcentr,ledge
             zsum = zsum+dx2i(jz)*alphai(jz)*uisd
             zsum11 = zsum11+dx2i(jz)*rhobis(2,jz)
 300         continue
         zimfast = (4.0*zsum+1.0*zsum11)/(2.0*zsum+zsum11) ! effective fast ion charge
          if(zimfast < epslon) zimfast = 0.0 ! if no fast ion
c
          zsum = 0.0
          zsum11 = 0.0
c
c  aspec = atomic weight of ion species, set in file default
c
          if(mhyd == 0) goto 390
          do 380 jz = 1, mzones
             do 370 j1=1,2
                do 366 j2 = 1, mhyd ! mhyd = number of hydrogen species
                   zsum = rhohs(j2,j1,jz) + zsum
                   zsum11 = aspec(j2) * rhohs(j2,j1,jz) + zsum11
 366            continue
 370         continue
 380      continue
          zmui = zsum11 / zsum ! average ion mass (amu)
 390      continue
          if(mhyd == 0) zmui=2.0  ! assume deuterium mass
c
          zsum = 0.0
          zsum11 = 0.0
          do 400 jz = 1,mzones
             zbeampr = hebems(jz)*rhobis(2,jz)
             zsum = zsum + alphai(jz)*ealfai(jz)*uisd*uise + zbeampr
             zsum11 = zsum11 + alphai(jz)*uisd + rhobis(2,jz)
 400      continue
c
          if ( zsum11 < epslon ) then
             zefast = 0.0  ! no fast ions
             zimfast=0.0
             zbetapf1=0.0
          else
             zefast = (zsum/zsum11)*(1.0e9/1.602) ! energy of fast ion (keV)
          endif
c
          izlsaw = lsawth(10)
c
           if(kcnnsaw == 0)then ! control number of calls after sawtooth crash
c
             if ( lsawth(10) .eq. 7 ) then ! call porcelli_model with arrays
c
               zbeta_factor = 8. * fcpi / ( bzs**2 )
c
               do jz=1,mzones
                 zrminor1(jz) = ahalfs(jz,1)*usil
                 zrminor2(jz) = ahalfs(jz,2)*usil
                 zvolume(jz)  = vols(jz,1)*1.e-6
                 zarea(jz)    = areas(jz,1)*1.e-4
                 zbpol2di(jz) = bpol2di(jz,1)
                 zbetai(jz) = rhoins(2,jz)*tis(2,jz) * zbeta_factor
                 zbetae(jz) = rhoels(2,jz)*tes(2,jz) * zbeta_factor
                 zbetat(jz) = totprs(jz) * zbeta_factor
               enddo
!
               do jz=1,mzones-1
                 zrdp1(jz) = ( thpsms(jz+1)-thpsms(jz) ) * zbeta_factor
     &                     / ( ( armins(jz+1,2)-armins(jz,2) ) * usil )
                 zztrappg(jz) = zbeta_factor * ( slpts(jz)
     &             +slprs(jz+1)*thrprs(jz+1)/max(epslon,armins(jz+1,1)))
     &              / usil
               enddo
c
               iprint = 0
               iout   = 6
c
c  rgb 22 Feb 2005: zbetai --> betati
c
               call porcelli_model ( iswtst, iprint, iout
     &          , mzones, zrq1
     &          , zrminor2, betati, zbetae
     &          , zbetat, shear, elong
     &          , zrminor1, zvolume
     &          , zarea
     &          , zrdp1, zztrappg, zbpol2di
     &          , zmui, zbtor1, rmaji, zdenel0, zti0, zefast )
c
c               call porcelli_model ( iswtst, iprint, iout
c     &          , mzones, zrq1
c     &          , zrminor2, zbetai, zbetae
c     &          , zbetat, shear, elong
c     &          , zrminor1, zvolume
c     &          , zarea
c     &          , zrdp1, zztrappg, zbpol2di
c     &          , zmui, zbtor1, rmaji, zdenel0, zti0, zefast
c     &          , ifastswt, iflarmswt, ikinswt )
c
               write(6,*) 'kcnnsaw= ',kcnnsaw,' after porcelli_model'
c
             else
c
              call saw_model(izlsaw,iswtst,ctrho,ctstar,cwrat,
     &          chfast,zcpf,zimfast,zelongr1,zmui,zlir1,zbetapf1,
     &          zbetaim,zshearp,rmaji,rmini,zrq1,zrplen1,zdenel0,
     &          zbtor1,zbpol11,zti0,zefast,zp1avg,zpresr1,
     &          ifastswt,iflarmswt,ikinswt,zntime)
               if(iswtst > 0) kcnnsaw = 1
c
             endif
c
           else
               kcnnsaw = kcnnsaw + 1 
               iswtst = 0
c
             if ( lsawth(10) .eq. 7 ) then ! call porcelli_model with arrays
c
               iprint = 0
               iout   = 6
c
               if(ksawstp < 2)
     &          call porcelli_model ( iswtst, iprint, iout
     &          , mzones, zrq1
     &          , zrminor2, zbetai, zbetae
     &          , zbetat, shear, elong
     &          , zrminor1, zvolume
     &          , zarea
     &          , zrdp1, zztrappg, zbpol2di
     &          , zmui, zbtor1, rmaji, zdenel0, zti0, zefast )
c
               write(6,*) 'kcnnsaw= ',kcnnsaw,' after porcelli_model'
c
               if(kcnnsaw > ksawstp) kcnnsaw = 0
c
             else
c
               if(ksawstp < 2)call saw_model(izlsaw,iswtst,ctrho,ctstar,
     &          cwrat,chfast,zcpf,zimfast,zelongr1,zmui,zlir1,zbetapf1,
     &          zbetaim,zshearp,rmaji,rmini,zrq1,zrplen1,zdenel0,
     &          zbtor1,zbpol11,zti0,zefast,zp1avg,zpresr1,
     &          ifastswt,iflarmswt,ikinswt,zntime)
               if(kcnnsaw > ksawstp) kcnnsaw = 0
c
             endif
c
           endif
c
c        if(iswtst .eq. 1) ztporck = ztime + lsawnt*dti*uist
c
          ztnt = zntime - tpsaw
          if(ztnt > tpold) ztnt = tpold
          if(iswtst == 1) then
             tpsaw = zntime - tpold   ! sawtooth period
c             tpold = zntime  ! new sawtooth reference time
          endif
c
          write(6,*) 'zntime= ', zntime,' iswtst= ',iswtst
c
          write(80,700)zntime,qsmth(2),zshearp,zlir1
c
 700      format(f8.4,2x,3(1pe10.3,2x))
c
          iporcelli = lsawth(27)
c
          if(iswtst .gt. 0 .and. iporc .gt. iporcelli)then
             iporc = 0
             zt1 = zntime - tpold
             write(90,*)zntime,'saw crash',zt1
             write(100,*)tpold,zt1
             write(100,*)zntime,zt1
             tpold = zntime
          elseif(iswtst .gt. 0 .and. iporc .le. iporcelli)then
             iswtst = 0
             iporc = iporc + 1
             write(90,*)zntime,'no1 crash',iporc 
          else
             iswtst = 0
             iporc = iporc + 1
             write(90,*)zntime,'no2 crash',iporc
          endif   
c
          if(lsawth(10) .eq. 6 ) goto 401
c
           return ! return out of lsawth(10) = 5 option
c
       endif                    
 401   continue
c
c
!| 
!| \subsection{Sawtooth Period and Indicator}
!| 
!| Finally, compute the sawtooth period {\tt ztaus}
!| and set ${\tt iswtst} = 1$ if conditions are right for a sawtooth crash
!| to occur.
!| 
c
c  compute sawtooth period and test for sawtooth crash
c..Set indicator to produce a sawtooth crash
c
          ztaus=1.e10
        if ( cfutz(465) .gt. 0. ) ztaus = tmguir
        if ( swperd .gt.0.) ztaus= swperd
        if ( (cfutz(463) .ne. 0.) .and. (ztime .ge. cfutz(463)) )
     &    ztaus = cfutz(464)
        if ( lsawth(10) .eq. 2 ) ztaus = ztpark
        if ( ztime-tbsaws .gt. ztaus ) iswtst = 1
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{roge91a} B. Rogers and L. Zakharov,
!| ``Stability conditons for ideal and tearing $m=1$ mode for circular
!| cross-section tokamaks,''
!| private communication.
!| 
!| \bibitem{roge91b} B. Rogers and L. Zakharov,
!| ``Estimates of stability conditions for $m=1$ mode for TFTR,''
!| private communication.
!| 
!| \bibitem{copp91a} B. Coppi and P. Detragiache,
!| ``Collisionless Magnetic Reconnection,''
!| 1991 Sherwood Fusion Conference (22--24 April 1991) 3B4.
!|
!| \bibitem{porcelli96a} F. Porcelli, D. Boucher, and M. N. Rosenbluth,
!| ``Model for the sawtooth period and amplitude,''
!| Plasma Physics and Controlled Fusion, vol. 38, page 2163-2186.
!| 
!| \end{thebibliography}
!| \end{document}

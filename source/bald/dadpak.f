c  18:00 28-nov-91 .../baldur/code/adpak/dadpak.s  Bateman, PPPL
c    This cshell compiles adpak and builds yadpak.a for BALDUR
c/ 15.12 07-aug-89 /11040/bald89/wbaldn1 DADPAK, written by Hulse, PPPL
c  BALDUR file DADPAK by Stotler, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type
c cfs get /11040/bald88/wbaldn1
c end
c lib wbaldn1 ^ x dadpak ^ end
c
c**********************************************************************c
c
c   This file contains the following subprograms:
c  adset  - set up basic atomic structure data for chosen element
c  adapn  - assign the no. of electrons per shell for each ionic species
c  admore - R. More excitation energies and ionization potentials
c  adshen - shell energies computed using R. More shielding constants
c  admayr - energy levels computed using Mayer shielding constants
c  adpali - evaluates Pauli approximation to Dirac energy level equation
c  adtip  - replaces calculated ionization potentials with tabulated ones
c  adfnm  - calculates strengths of transitions from shell in to shell im
c  adeci  - calculate electron collisional ionization rate
c  adbfi  - electron collisional ionization rate by Belfast formula
c  adbfit - evaluate Belfast fitting formulas
c  adygr  - electron collisional ionization rate by Younger formula
c  fchi   - function used in Younger formula
c  adrrec - computes radiative recombination rates
c  addrec - computes dielectronic recombination rates
c  addrcc - set up dielectronic recombination rate coefficient constant
c  adecex - computes radiation rates due to collisional excitation
c  adbrem - computes bremsstrahlung radiation for each ionic species
c  expe1  - calculates exp(x) * e1(x) -> exponential integral
c  expund - calculates exp(x) without underflow error
c  adbcxr - computes recombination rate due to charge exchange 
c  adcxxs - supplies charge exchange cross sections
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine adset(knucz)
c
c     rev: 10/10/82 for new ad package
c  rap 02-may-00 control for khole added
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c
c     setup basic atomic structure data for chosen element
c     (shell populations, energies, etc.)
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adparm / laden, ladtip, ltipok, leci, leciok
c
      nucz = knucz
c
      call adapn
c
      if(laden .ne. 0 .and. laden .ne. 1) stop
      if(laden .eq. 0) call admayr
      if(laden .eq. 1) call admore
c
      if(ladtip .ne. 0 .and. ladtip .ne. 1) stop
      ltipok = 0
      if(ladtip .eq. 1) call adtip(ltipok)
c
      call adfnm
c
      return
      end
      subroutine adapn
c
c     rev: 10/9/82 for new ad package
c
c     rev: 6/15/78
c
c     input : nucz
c     output : nspc, apn, aqn, nvalnc
c
c     assign the number of electrons per shell for each
c     ionic species.  apn(jq,jn) = number of electrons in
c     shell jn for species with charge jq-1, while aqn(jq,jn) =
c     (2.*jn**2 - apn(jq,jn))/(2.*jn**2) is the fractional vacancy.
c     the shells are filled in simple consecutive ascending order.
c     also, the principal quantum number of the highest
c     occupied ('valence') shell for each species is stored
c     in nvalnc(jq).
c
c     note:  for the fully ionized species nvalnc is defaulted to 1.
c     also: nspc = nucz + 1 is set by this routine
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
c     set maximum species index nspc = nucz + 1
c
      nspc = nucz +1
c
c
      do 1 jq = 1,nspc
      icount = nucz-jq+1
c
c
      do 2 jn = 1,10
      ifull = 2*jn*jn
      ishell = min0 (icount,ifull)
      icount = icount - ishell
      if(ishell.gt.0) nvalnc(jq) = jn
      apn(jq,jn) = float(ishell)
      aqn(jq,jn) = float(ifull - ishell)/float(ifull)
c
cw    write(5,900)jq,jn,apn(jq,jn),nvalnc(jq)
cw  900 format(1x,'jq =',i3,' jn =',i3,'  apn=',f4.0,'  nvalnc=',i3)
c
    2 continue
c
c
    1 continue
c
c
      nvalnc(nspc) = 1
c
      return
      end
      subroutine admore
c
c     rev: 10/10/82 for new ad package
c
c     rev: 1 september 1981
c
c     input : nucz, nspc, apn, nvalnc
c     output : en, ein, enm, qn
c
c     excitation energies and ionization potentials are
c     tabulated for all occupied shells for each species.  energies
c     are in kev, jq = species index = ionic charge + 1, and
c     jn = shell index.
c
c     uses the r. more revised shielding constants (and new
c     calculation procedure) in lieu of the original 'xsnq'
c     treatment based on the h. mayer shielding constants.
c     (original treatment, refs. 1-4; present treatment, ref. 5).
c     certain aspects of the formalism (e.g. qn and en arrays) are
c     'historical' in nature.
c
c
c          sigma(jn,jm)   r. more shielding constants
c
c          qn(jq,jn)      effective nuclear charge at shell jn
c                         considering all electrons in the ground state
c                         ion, now computed from zsen (rather than
c                         than vice-versa).  this is *not* the
c                         same as qn in more's notation, which
c                         is the shielded charge due to inner
c                         shielding only.
c
c          en(jq,jn)      energy levels for recombination of the
c                         ground state ion.
c
c          ein(jq,jn)     ionization potential for electron in shell
c                         jn calculated from the difference in total
c                         ion energies with and without an electron
c                         removed from shell jn.
c
c          enm(jq,jn,jm)  excitation energy from shell jn up to shell
c                         jm calulated from the difference in
c                         the total ion energies of of the ground
c                         state and excited state electron populations.
c
c          zbetot(jq)     total binding energy for each ion.
c
c
c     references:
c
c     1) post et al., pppl-1352 (corrected)
c     2) lokke and grasberger, lll report ucrl-52276 (corrected)
c     3) h. mayer, 'methods of opacity calculations', los alamos
c        scientific  laboratory la-647 (1947).
c     4) xsnq listing
c     5) r. more 'atomic physics in inertial confinement fusion -
c        part i' lll report ucrl-84991 (1981) (to appear
c        in 'applied atomic collision physics', volume ii,
c        academic press)
c
c     note: there are discrepancies between these references for a few
c           of the h. mayer shielding constant values.
c
c
      dimension zsen(10), zbetot(100), zdum(10)
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
c
c     main species loop
c
      do 100 jq = 1, nspc
      jjq = jq
      ivalnc = nvalnc(jq)
c
c     compute zbetot and zsen
c
      ihole = 0
      iexcit = 0
      call adshen(jjq, ihole, iexcit, zsen(1), zbetot(jjq))
c
c     compute ein, enm, and en from differences in total
c     ion energies (rather from zsen's) in order to agree
c     with results in more's figure iii-7.  compute
c     qn from the zsen in order to have the appropriate limiting
c     integral charges at high n shells.
c
c     clear arrays
c
      do 30 jn = 1, 5
c
      ein(jq,jn) = 0.0
c
      do 20 jm = 1, 10
      enm(jq,jn,jm) = 0.0
   20 continue
   30 continue
c
c     skip fully stripped state
c
      if(jq. eq. nspc) go to 60
c
      do 50 jn = 1, ivalnc
c
      ihole = jn
      iexcit = 0
      call adshen(jjq, ihole, iexcit, zdum(1), zbeiz)
      ein(jq,jn) = zbetot(jq) - zbeiz
c
      iijm = ivalnc
      if(jn .eq. ivalnc .or. aqn(jq,ivalnc) .eq. 0.0) iijm = iijm+1
c
      do 40 jm = iijm, 10
      ihole = jn
      iexcit = jm
      call adshen(jjq, ihole, iexcit, zdum(1), zbeex)
      enm(jq,jn,jm) = zbetot(jq) - zbeex
   40 continue
   50 continue
c
c     clear arrays
c
   60 do 70 jn = 1, 10
      en(jq,jn) = 0.0
      qn(jq,jn) = 0.0
   70 continue
c
c     skip neutral state
c
      if(jq .eq. 1) go to 100
c
c     compute en and qn
c
      iijn = ivalnc
      if(aqn(jq,ivalnc) .eq. 0.0) iijn = ivalnc+1
c
      do 80 jn = iijn, 10
      ihole = 0
      iexcit = jn
      call adshen(jjq, ihole, iexcit, zdum(1), zberec)
      en(jq,jn) = zberec - zbetot(jq)
      qn(jq,jn) = float(jn) * sqrt(zsen(jn)/0.0136)
   80 continue

  100 continue
c
      return
      end
c
c---------------------------------------------------------------
c
      subroutine adshen(kq, khole, kexcit, psen, pbetot)
      external case1
c
c     rev: 10/10/82 for new ad package
c
c     rev: 1 september 1981
c
c     input: nucz, kq, khole, kexcit, apn
c     output: psen(10), pbetot
c
c     compute shell energies using new shielding constants
c     and method due to r. more (ucrl-84991).
c
c     nucz = ion nuclear charge
c     kq = species index = ionic charge + 1
c     khole, kexcit = electron removed from shell khole and/or
c                     added (excited) to shell kexcit before
c                     energies are computed (khole = 0 no electron
c                     removed, kexcit = 0 no electron added)
c     apn(100,10) = array of shell populations
c     psen = output array of shell energies (kev) (defined positive)
c     pbetot = output total ion energy (defined positive)
c
      dimension psen(10), zapn(10), zq(10), zen0(10)
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / smore / sigmor(10,10)
c
c      data ((sigmor(jn,jm),jn=1,10),jm=1,10)/
c     1 0.3125, 0.9380, 0.9840, 0.9954, 0.9970,
c     1 0.9970, 0.9990, 0.9999, 0.9999, 0.9999,
c     2 0.2345, 0.6038, 0.9040, 0.9722, 0.9979,
c     2 0.9880, 0.9900, 0.9990, 0.9990, 0.9999,
c     3 0.1093, 0.4018, 0.6800, 0.9155, 0.9796,
c     3 0.9820, 0.9860, 0.9900, 0.9920, 0.9999,
c     4 0.0622, 0.2430, 0.5150, 0.7100, 0.9200,
c     4 0.9600, 0.9750, 0.9830, 0.9860, 0.9900,
c     5 0.0399, 0.1597, 0.3527, 0.5888, 0.7320,
c     5 0.8300, 0.9000, 0.9500, 0.9700, 0.9800,
c     6 0.0277, 0.1098, 0.2455, 0.4267, 0.5764,
c     6 0.7248, 0.8300, 0.9000, 0.9500, 0.9700,
c     7 0.0204, 0.0808, 0.1811, 0.3184, 0.4592,
c     7 0.6098, 0.7374, 0.8300, 0.9000, 0.9500,
c     8 0.0156, 0.0624, 0.1392, 0.2457, 0.3711,
c     8 0.5062, 0.6355, 0.7441, 0.8300, 0.9000,
c     9 0.0123, 0.0493, 0.1102, 0.1948, 0.2994,
c     9 0.4222, 0.5444, 0.6558, 0.7553, 0.8300, 0.0100, 0.0400,
c     & 0.0900, 0.1584, 0.2450, 0.3492, 0.4655, 0.5760, 0.6723, 0.7612/
c
      do 10 jn = 1, 10
      zapn(jn) = apn(kq,jn)
   10 continue
c
c     remove electron from shell khole
cap {
      if(khole .gt. 0) then 
        if (zapn(khole) .gt. 0.0) zapn(khole) = zapn(khole) - 1.0
      endif 
cap }
c     add electron to shell kexcit
c
      if(kexcit .gt. 0) zapn(kexcit) = zapn(kexcit) + 1.0
c
c     compute nuclear charge as shielded by inner electrons only
c
      do 30 jn = 1, 10
c
      zq(jn) = float(nucz) - 0.5 * sigmor(jn,jn) * zapn(jn)
c
      if(jn .eq. 1) go to 30
      jnm1 = jn - 1
c
      do 20 jm = 1, jnm1
      zq(jn) = zq(jn) - sigmor(jn,jm) * zapn(jm)
   20 continue
   30 continue
c
c     compute potential due to outer electrons
c
      do 50 jn = 1, 10
      zen0(jn) = 0.0136 * sigmor(jn,jn) * zapn(jn)
     &           * zq(jn) / (float(jn)**2)
c
      if(jn .eq. 10) go to 50
      jnp1 = jn + 1
c
      do 40 jm = jnp1, 10
      zen0(jn) = zen0(jn) + 0.0272 * zapn(jm) * sigmor(jm,jn)
     &                      * zq(jm) / (float(jm)**2)
   40 continue
c
   50 continue
c
c     compute final shell energies and total ion binding energy,
c     both defined positive.
c
      pbetot = 0.0
c
      do 60 jn = 1, 10
      psen(jn) = 0.0136 * ((zq(jn)/float(jn))**2) - zen0(jn)
      pbetot = pbetot + 0.0136 * zapn(jn) * ((zq(jn)/float(jn))**2)
   60 continue
c
      return
      end
c******************
      BLOCK DATA case1

      common / smore / sigmor(10,10)
c
      data ((sigmor(jn,jm),jn=1,10),jm=1,10)/
     1 0.3125, 0.9380, 0.9840, 0.9954, 0.9970,
     1 0.9970, 0.9990, 0.9999, 0.9999, 0.9999,
     2 0.2345, 0.6038, 0.9040, 0.9722, 0.9979,
     2 0.9880, 0.9900, 0.9990, 0.9990, 0.9999,
     3 0.1093, 0.4018, 0.6800, 0.9155, 0.9796,
     3 0.9820, 0.9860, 0.9900, 0.9920, 0.9999,
     4 0.0622, 0.2430, 0.5150, 0.7100, 0.9200,
     4 0.9600, 0.9750, 0.9830, 0.9860, 0.9900,
     5 0.0399, 0.1597, 0.3527, 0.5888, 0.7320,
     5 0.8300, 0.9000, 0.9500, 0.9700, 0.9800,
     6 0.0277, 0.1098, 0.2455, 0.4267, 0.5764,
     6 0.7248, 0.8300, 0.9000, 0.9500, 0.9700,
     7 0.0204, 0.0808, 0.1811, 0.3184, 0.4592,
     7 0.6098, 0.7374, 0.8300, 0.9000, 0.9500,
     8 0.0156, 0.0624, 0.1392, 0.2457, 0.3711,
     8 0.5062, 0.6355, 0.7441, 0.8300, 0.9000,
     9 0.0123, 0.0493, 0.1102, 0.1948, 0.2994,
     9 0.4222, 0.5444, 0.6558, 0.7553, 0.8300, 0.0100, 0.0400,
     & 0.0900, 0.1584, 0.2450, 0.3492, 0.4655, 0.5760, 0.6723, 0.7612/

      end

      subroutine admayr
      external case2
c
c     rev: 10/10/82 for new ad package
c
c     rev: 6/15/78
c
c     input : nucz, nspc, apn
c     output : sigma, qn, en, ein, enm
c
c     computes energy levels for the ground state ionic species
c     and for species with single inner electron vacancies.
c     excitation energies and ionization potentials are thus
c     tabulated for all occupied shells for each species.  energies
c     are in kev, jq = species index = ionic charge + 1, and
c     jn = shell index.
c
c     for each configuration and shell of interest, the effective
c     nuclear charge is calculated using the h. mayer shielding
c     constants.  the energy levels are then found using this effective
c     charge in the pauli approximation to the dirac energy level
c     formula.  the shell energy is taken as the mean of the lowest
c     and highest subshell energies for that n.
c
c
c          sigma(jn,jm)   h. mayer shielding constants
c
c          qn(jq,jn)      effective nuclear charge at shell jn
c                         considering all electrons in the ground state ion.
c
c          en(jq,jn)      energy levels for ground state ion
c                         (formerly eqn(jq,jn) in previous versions)
c
c          ein(jq,jn)     ionization potential for electron in shell jn;
c                         uses shielded charge calculated with one
c                         electron removed from shell jn.
c
c          enm(jq,jn,jm)  excitation energy from shell jn up to shell jm;
c                         the difference between the jn and jm energy levels
c                         calculated with the jn shell electron removed.
c
c
c
c     references:
c
c     1) post et al., pppl-1352 (corrected)
c     2) lokke and grasberger, lll report ucrl-52276 (corrected)
c     3) h. mayer, 'methods of opacity calculations', los alamos
c        scientific  laboratory la-647 (1947).
c     4) xsnq listing
c
c     note: there are discrepancies between these references for a few
c           shielding constant values.
c
c
      dimension zsigma(10,10)
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
c
      common / smayer / sigmay(10,10)
c
c
c     start main species loop
c
c
      do 100 jq = 1 , nspc
c
      ivalnc = nvalnc(jq)
c
c
c     compute qn and en
c
      do 30 jn = 1 , 10
      zjn = float(jn)
c
      zqn = float(nucz)
c
      do 20 jm = 1 , ivalnc
      zqn = zqn - sigmay(jn,jm) * apn(jq,jm)
   20 continue
c
      zqn = zqn + apn(jq,jn) * sigmay(jn,jn) / (2. * zjn**2)
c
      qn(jq,jn) = zqn
      in = jn
      en(jq,jn) = adpali(zqn,in)
   30 continue
c
c
c     clear arrays then compute ein and enm
c
c
      do 40 jn = 1 , 5
      ein(jq,jn) = 0.
      do 40 jm = 1 , 10
      enm(jq,jn,jm) = 0.
   40 continue
c
c     skip ionization and excitation energies for fully stripped ion
c
      if(jq .eq. nspc) go to 100
c
c
      do 70 jn = 1 , ivalnc
      zjn = float(jn)
c
c     compute zqnn = effective charge at shell jn with one electron
c     removed from shell jn.
c
      zqnn = qn(jq,jn) + sigmay(jn,jn) * (1. - 1./(2.*zjn**2))
c
c     ein(jq,jn) = ionization potential for shell jn
c
      in = jn
      ein(jq,jn) = adpali(zqnn,in)
c
c     compute enm(jq,jn,jm) from zqnm, effective charge at shell
c     jm with one electron removed from shell jn.
c
      iijm = ivalnc
      if(jn .eq. ivalnc .or.
     &         aqn(jq,ivalnc) .eq. 0.) iijm = ivalnc + 1
c
      do 60 jm = iijm , 10
c
      zqnm = qn(jq,jm) + sigmay(jm,jn)
c
      im = jm
      enm(jq,jn,jm) = ein(jq,jn) - adpali(zqnm,im)
c
   60 continue
c
c
   70 continue
c
c
  100 continue
c
c
      return
      end
c
c
c******************
      BLOCK DATA case2
      common / smayer / sigmay(10,10)
c
      data ((sigmay(jn,jm),jn=1,10),jm=1,10)/
     1 .625,.938,.981,.987,.994,.997,.999,1.00,1.00,1.00,
     2 .235,.690,.893,.940,.970,.984,.990,.993,.995,1.00,
     3 .109,.397,.702,.850,.920,.955,.970,.980,.990,1.00,
     4 .0617,.235,.478,.705,.830,.900,.950,.970,.980,.990,
     5 .0398,.155,.331,.531,.720,.830,.900,.950,.970,.980,
     6 .0277,.109,.239,.400,.580,.735,.830,.900,.950,.970,
     7 .0204,.0808,.178,.310,.459,.610,.745,.830,.900,.950,
     8 .0156,.0625,.138,.243,.371,.506,.635,.750,.830,.900,
     9 .0123,.0494,.111,.194,.299,.431,.544,.656,.760,.830,
     & .010,.040,.090,.158,.245,.353,.460,.576,.670,.765/
      end
c------------------------------------------------------------------
c
c
      function adpali(pq,kn)
c
c     rev: 10/10/82 for new ad package
c
c     rev: 6/15/78
c
c     evaluate the pauli approximation to the
c     the dirac energy level equation, taking the mean of the
c     lowest and highest subshell energies to be the shell energy.
c     adpali(pq,kn) is in kev,  pq is the effective nuclear
c     charge (in units of electronic charge) and kn is the principal
c     quantum number.
c
c
      zalpha = .007298
c
c
      zkn = float(kn)
c
      adpali = .0136 * ( (pq / zkn)**2 ) *
     &     (1. + ((zalpha*pq)**2) * (0.5 - 1./(4.*zkn)) / zkn)
c
c
      return
      end
      subroutine adtip(ktipok)
c
c     substitute tabulated ionization potentials for those
c     calculated by shielding constants in setenm.  note that the
c     'number of equivalent electrons' is not correctly re-defined
c     to accompany this procedure, but the approximation is still
c     an improvement for near-neutral species.
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      dimension ztip(100)
c
      do 1 jspc = 1, 100
      ztip(jspc) = 0.0
    1 continue
c
      if(nucz .ne. 2) go to 3
c
c     helium
c
      ztip(1) = 24.58e-3
      ztip(2) = 54.4e-3
c
      go to 100
    3 if(nucz .ne. 3) go to 4
c
c     lithium
c
      ztip(1) = 5.4e-3
      ztip(2) = 75.6e-3
      ztip(3) = 122.4e-3
c
      go to 100
c
    4 if(nucz .ne. 4) go to 5
c
c     beryllium
c
      ztip(1) = 9.32e-3
      ztip(2) = 18.2e-3
      ztip(3) = 153.9e-3
      ztip(4) = 217.6e-3
c
      go to 100
c
    5 if(nucz .ne. 5) go to 6
c
c     boron
c
      ztip(1) = 8.296e-3
      ztip(2) = 25.15e-3
      ztip(3) = 37.9e-3
      ztip(4) = 259.3e-3
      ztip(5) = 340.0e-3
c
      go to 100
c
    6 if(nucz .ne. 6) go to 8
c
c     carbon
c
      ztip(1) = 11.3e-3
      ztip(2) = 24.4e-3
      ztip(3) = 47.9e-3
      ztip(4) = 64.5e-3
      ztip(5) = 0.392
      ztip(6) = 0.490
c
      go to 100
    8 if(nucz .ne. 8) go to 10
c
c     oxygen
c
      ztip(1) = 13.6e-3
      ztip(2) = 35.1e-3
      ztip(3) = 54.9e-3
      ztip(4) = 77.4e-3
      ztip(5) = 0.114
      ztip(6) = 0.138
      ztip(7) = 0.739
      ztip(8) = 0.871
c
      go to 100
   10 if(nucz .ne. 10) go to 13
c
c     neon
c
      ztip(1) = 21.6e-3
      ztip(2) = 41.0e-3
      ztip(3) = 63.5e-3
      ztip(4) = 97.1e-3
      ztip(5) = 0.126
      ztip(6) = 0.158
      ztip(7) = 0.207
      ztip(8) = 0.239
      ztip(9) = 1.196
      ztip(10) = 1.362
c
      go to 100
   13 if(nucz .ne. 13) go to 18
c
c     aluminum
c
      ztip(1) = 5.99e-3
      ztip(2) = 18.8e-3
      ztip(3) = 28.4e-3
      ztip(4) = 0.120
      ztip(5) = 0.154
      ztip(6) = 0.190
      ztip(7) = 0.241
      ztip(8) = 0.285
      ztip(9) = 0.330
      ztip(10) = 0.399
      ztip(11) = 0.442
      ztip(12) = 2.086
      ztip(13) = 2.304
c
      go to 100
   18 if(nucz .ne. 18) go to 22
c
c     argon
c
      ztip(1) = 15.8e-3
      ztip(2) = 27.6e-3
      ztip(3) = 40.7e-3
      ztip(4) = 59.8e-3
      ztip(5) = 75.0e-3
      ztip(6) = 91.0e-3
      ztip(7) = 0.124
      ztip(8) = 0.143
      ztip(9) = 0.422
      ztip(10) = 0.479
      ztip(11) = 0.539
      ztip(12) = 0.618
      ztip(13) = 0.686
      ztip(14) = 0.756
      ztip(15) = 0.855
      ztip(16) = 0.918
      ztip(17) = 4.121
      ztip(18) = 4.426
c
      go to 100
   22 if(nucz .ne. 22) go to 26
c
c     titanium
c
      ztip(1) = 6.82e-3
      ztip(2) = 13.6e-3
      ztip(3) = 27.5e-3
      ztip(4) = 43.3e-3
      ztip(5) = 99.2e-3
      ztip(6) = 0.119
      ztip(7) = 0.141
      ztip(8) = 0.169
      ztip(9) = 0.193
      ztip(10) = 0.216
      ztip(11) = 0.265
      ztip(12) = 0.292
      ztip(13) = 0.787
      ztip(14) = 0.861
      ztip(15) = 0.940
      ztip(16) = 1.042
      ztip(17) = 1.131
      ztip(18) = 1.220
      ztip(19) = 1.342
      ztip(20) = 1.425
      ztip(21) = 6.249
      ztip(22) = 6.626
c
      go to 100
   26 if(nucz .ne. 26) go to 100
c
c     iron
c
      ztip(1) = 7.87e-3
      ztip(2) = 16.2e-3
      ztip(3) = 30.7e-3
      ztip(4) = 54.8e-3
      ztip(5) = 75.0e-3
      ztip(6) = 99.0e-3
      ztip(7) = 0.125
      ztip(8) = 0.151
      ztip(9) = 0.235
      ztip(10) = 0.262
      ztip(11) = 0.290
      ztip(12) = 0.331
      ztip(13) = 0.361
      ztip(14) = 0.392
      ztip(15) = 0.457
      ztip(16) = 0.490
      ztip(17) = 1.265
      ztip(18) = 1.358
      ztip(19) = 1.456
      ztip(20) = 1.582
      ztip(21) = 1.689
      ztip(22) = 1.799
      ztip(23) = 1.950
      ztip(24) = 2.045
      ztip(25) = 8.828
      ztip(26) = 9.278
c
c     end of current tabulation; return if selected
c     element not present and set flag.
c
  100 ktipok = 0
      if(ztip(1) .eq. 0.0) return
c
      ktipok = 1
c
      do 150 jspc = 1, nucz
      ivalnc = nvalnc(jspc)
      jspcp1 = jspc + 1
c
      ein(jspc,ivalnc) = ztip(jspc)
      en(jspcp1,ivalnc) = ztip(jspc)
  150 continue
c
      return
      end
      subroutine adfnm
      external case3
c
c     rev: 10/10/82 for new ad package
c
c     rev: 8/10/79  correct zepsnn(59) (4f13) value from 3.43e-3
c                   to 3.43e-2 (error in post et al.)
c     rev: 6/15/78
c
c     input : nucz, apn, aqn, nvalnc, enm
c     output : fnm, fnn, enn
c
c     calculates absorption oscillator strengths for transitions
c     from shell in to shell im for each species jq and stores them
c     in the array fnm(jq,in,im).  f(n,m) = (hydrogenic value,
c     tabulated for low n,m calculated otherwise) x (correction
c     factor(s)) x (electron population of shell in) x (vacany
c     fraction of shell im).
c
c     for n - n (e.g. delta n equal zero) transitions we have
c     f(n,n) = ann + bnn/znuc, which represents a fit to the
c     'equivalent' single transition used to model all the physical
c     n - n transitions.  also computed is e(n,n), the 'equivalent'
c     energy, e(n,n) = epsilonnn * (zchg+1) ** alphann.  the zann, etc.
c     arrays are indexed by the total number of bound electrons
c     in their electronic configurations (see reference (1)).
c
c     note that here, as in other routines, the filling of shells
c     in simple ascending order is assumed both explicitly
c     (e.g., in the apn array contents) as well as implicity (through
c     the use of nvalnc, etc.).  also, the ann, bnn, etc.
c     coefficients are defined specifically for such states.
c
c
c     references:
c     1) post et al., pppl-1352 (with corrections)
c     2) grasberger, xsnq - tokamak memorandum #4, 7/29/76
c     3) grasberger, xsnq - tokamak memorandum #10, 1/5/77
c     4) lokke & grasberger, 'xsnq-u a non-lte emission and
c        absorption coefficient subroutine', lll report ucrl-52276
c     5) merts et al., 'the calculated power output from
c        a thin iron-seeded plasma', los alamos report la-6220-ms
c     6) xsnq fortran listing
c
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
c
      common / fnmdat / zann(60), zbnn(60), zepsnn(60), zalpnn(60),
     &          zy0(10), zy1(10), zy2(10)
c
c
c     *********************************************************
c
c
c     zero oscillator strength arrays and enn array
c
      do 1 jq = 1 , 100
      fnn(jq) = 0.
      enn(jq) = 0.
      do 1 jn = 1 , 5
      do 1 jm = 1 , 10
    1 fnm(jq,jn,jm) = 0.
c
      znuc = float ( nucz )
c
c
c     ******************************************************
c     ***** loop for all species except fully stripped *****
c     ******************************************************
c
c
      do 100 jq = 1 , nucz
c
      zchg = float(jq-1)
      ibound = nucz - jq + 1
      ivalnc = nvalnc(jq)
c
c
c     ********** calculate f(n,m) ( m > n ) **********
c
c
c     begin by evaluting hydrogenic f(n,m) equation, with
c     explicit values substituted for n < 5, m <= 5.
c
c
      do 10 jn = 1 , ivalnc
      iijm = max0 ( jn + 1 , ivalnc )
      zn = float(jn)
c
      do 10 jm = iijm , 10
      zm = float(jm)
      fnm(jq,jn,jm) = 1.96 * zn * ( ( zm/(zm*zm - zn*zn) )**3 )
   10 continue
c
c     substitute tabled values for low n , m as per ref #4
c
      fnm(jq,1,2) = 0.4161
      fnm(jq,1,3) = 0.0792
      fnm(jq,1,4) = 0.0290
      fnm(jq,1,5) = 0.0139
      fnm(jq,2,3) = 0.6370
      fnm(jq,2,4) = 0.1190
      fnm(jq,2,5) = 0.0443
      fnm(jq,3,4) = 0.8408
      fnm(jq,3,5) = 0.1499
      fnm(jq,4,5) = 1.037
c
c     modify nonzero f(n,n+1) values for n > 1 as per ref #3
c
      do 20 j = 1 , 2
      jn = ivalnc + j - 2
      if(jn .le. 1) go to 20
      jm = jn + 1
      if(aqn(jq,jm) .eq. 0.) go to 20
c
      znchg = zchg + apn(jq,jm)
      zehydr = .0136 * ((znchg+1.)**2) *
     &                   ( 1./float(jn)**2 - 1./float(jm)**2 )
      zratio = zehydr / enm(jq,jn,jm)
      if(zratio .gt. 1.0) zratio = 1.0
      zynpn = zy0(jn) + zy1(jn) * apn(jq,jn) + zy2(jn) * apn(jq,jn)**2
c
      fnm(jq,jn,jm) = fnm(jq,jn,jm) * zratio * zynpn
c
   20 continue
c
c
c     multiply by electron population of shell n and fractional
c     vacancy of shell m
c
c
      do 30 jn = 1 , 5
      do 30 jm = 1 , 10
      fnm(jq,jn,jm) = apn(jq,jn) * aqn(jq,jm) * fnm(jq,jn,jm)
   30 continue
c
c
c
c     ***** calculate n - n 'effective' oscillator strengths *****
c                   ***** and transition energies *****
c
c
c
      index = ibound
c     if more than 60 bound electrons use constants for ion with
c     32 fewer electrons
      if(index.gt.60)index = index - 32
      if(index.gt.60)index = index - 32
c
      fnn(jq) = zann(index) + zbnn(index) / znuc
c
      enn(jq) = zepsnn(index) * (zchg + 1.)**zalpnn(index)
c
c
c
c
c     check to insure all f values are positive
c
      if(fnn(jq) .ge. 0.) go to 40
      fnn(jq) = 0.
      write(5,900) jq , ivalnc , ivalnc
  900 format(1x,26h***** negative f value at ,3i4,6h *****)
   40 do 41 jn = 1 , 5
      do 41 jm = 1 , 10
      if(fnm(jq,jn,jm) .ge. 0.) go to 41
      fnm(jq,jn,jm) = 0.
      write(5,900) jq , jn , jm
   41 continue
c
c     diagnostic write statments to display f values
c
cw    write(5,901)jq,fnn(jq),enn(jq)
cw  901 format(/1x,'setfnm: jq =',i3,'  fnn =',1pe9.2,'  enn =',e9.2)
cw    write(5,902)( jn,(fnm(jq,jn,jm),jm=1,10),jn=1,5)
cw  902 format(1x,'jn =',i2,10f6.1)
c
c
c     *************** end of jq species loop ***************
c
c
  100 continue
c
      return
      end

c******************
      BLOCK DATA case3

      common / fnmdat / zann(60), zbnn(60), zepsnn(60), zalpnn(60),
     &          zy0(10), zy1(10), zy2(10)
c
c
      data (zann(j),j=1,60) / 0., 0.,
     1 -.02, -.01, -.04, .04, -.01, -.01, .01, 0.,
     2  .01, .1, .25, .17, .08, -.02, -.14, -.27,
     3 -.29, -.3, -.3, -.29, -.27, -.24, -.2, -.14,
     4 -.08, 0., .97, 1.96, 1.92, 1.89, 1.86, 1.83,
     5 1.78, 1.73, 1.41, 1.05, .674, .261, -.167,
     6 -.637, -1.14, -1.67, -2.26, -2.88, -2.9,
     7 -2.83, -2.72, -2.61, -2.45, -2.27, -2.05, -1.81,
     8 -1.55, -1.25, -.94, -.63, -.31, 0. /
c
c
      data (zbnn(j),j=1,60) / 0., 0.,
     1 2., 4.33, 4.68, 3.6, 2.85, 2.4, 1.3, 0.,
     2 10.5, 20., 16.4, 24.7, 34.1, 44.8, 56.8, 70.3,
     3 68.6, 65.8, 61.9, 56.9, 50.7, 45.5, 34.4,
     4 24.2, 12.8, 0., -25.8, -53.1, -28.1, -2.97,
     5 23., 50.3, 79.1, 109.2, 141.7, 176., 213.9,
     6 256.6, 299.6, 348.2, 400.5, 456.6, 518.1,
     7 584.0, 571.9, 550.7, 524.9, 496.8, 463.4,
     8 427., 385.3, 339.8, 291.3, 237.4, 181.3, 122.9,
     9 62.2, 0. /
c
c
      data (zepsnn(j),j=1,28) / 0., 0.,
     1 2.04e-3, 4.49e-3, 6.8e-3, 8.16e-3, 6.80e-3,
     2 1.06e-2, 1.31e-2, 0., 2.3e-3, 3.88e-3,
     3 5.71e-3, 5.44e-3, 8.16e-3, 6.8e-3, 8.3e-3,
     4 4.32e-3, 5.11e-3, 6.04e-3, 7.08e-3,
     5 8.12e-3, 9.69e-3, 1.13e-2, 1.37e-2, 1.49e-2,
     6 1.71e-2, 0. /
c
      data (zepsnn(j),j=29,60) / 1.52e-5, 1.93e-5,
     1 5.62e-5, 1.20e-4, 2.13e-4, 3.34e-4, 4.66e-4,
     2 6.66e-4, 9.10e-4, 1.25e-3, 1.69e-3, 1.98e-3,
     3 3.03e-3, 3.92e-3, 5.19e-3, 6.40e-3, 8.05e-3,
     4 9.80e-3, 1.08e-2, 1.19e-2, 1.32e-2, 1.47e-2,
     5 1.59e-2, 1.81e-2, 1.92e-2, 2.10e-2, 2.33e-2,
     6 2.57e-2, 2.82e-2, 3.10e-2, 3.43e-2, 0. /
c
c
      data (zalpnn(j),j=1,60) / 1.0, 1.0,
     1 1.0, .93, .86, .80, .87, .77, .77, 1.0,
     2 .99, .90, .83, .87, .77, .82, .79, 1.06,
     3 1.03, 1., .97, .94, .91, .88, .85,
     4 .83, .80, 1.0, 2.46, 2.4, 2.14, 1.96,
     5 1.83, 1.72, 1.64, 1.57, 1.49, 1.41, 1.34,
     6 1.30, 1.19, 1.13, 1.06, 1.01, .95, .91,
     7 .89, .87, .85, .83, .82, .80, .78, .76,
     8 .74, .72, .70, .68, .66, 1.0 /
c
c
      data (zy0(j),j=1,10) / 1.0, .511, .400, 7*.230 /
c
      data (zy1(j),j=1,10) / 0., .116, .0228, 7*.00844 /
c
      data (zy2(j),j=1,10) / 0., -.00689, .000596, 7*.000488 /
c
      end
c------------------------------------------------------------------
c
      subroutine adeci(pte, pne)
c
c     calculate the electron collisional ionization rate using
c     the formulas selected by the leci switch.  leciok returns
c     the rate formulas actually used, which may differ from
c     leci if leci was set inappropriately for the given element.
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
      common / adparm / laden, ladtip, ltipok, leci, leciok
c
c     ilmax = current maximum number of rate options
c
      ilmax = 3
c
c     check for valid leci; return leci = 1 option if not valid
c
      if(leci .lt. 1 .or. leci .gt. ilmax) go to 10
c
      go to (10, 20, 30) leci
c
c     leci = 1; use xsnq derived ionization rate calculation
c
   10 call adxeci(pte, pne)
      leciok = 1
      return
c
c     leci = 2; use belfast group light element rate fits
c
   20 call adbfi(nucz, pte, pne, rclion, cizlos, iflag)
      if(iflag .ne. 0) go to 10
      leciok = 2
      return
c
c     leci = 3; use younger's rates for available elements and shells,
c               xsnq rates otherwise on a charge state-by-charge state
c               basis
c
   30 call adxeci(pte, pne)
      call adygr(nucz, pte, pne, rclion, cizlos, iflag)
      leciok = 3
      if(iflag .ne. 0) leciok = 1
      return
c
c
c
      end
      subroutine adxeci(pte, pne)
c
c     rev: 3/30/83 adeci --> adxeci to accomodate alternative rates
c     rev: 10/10/82 for new ad package
c
c     rev: 6/15/78
c
c     computes the xsnq-derived collisional ionization rates rclion
c     at the electron temperature pte (kev), where jq=(ionic charge+1).
c     these rates have the units (sec-1).
c
c     notes:
c
c     1) electrons are assumed to fill ionic shells in simple ascending
c        order.  this is assumed both explicitly (the apn array)
c        as well as implicitly through the use of nvalnc values
c        in loop indices, etc.
c     2) through use of the function expund, exp(x) is set equal to
c        zero in the rate equations for values of x < -45.
c     3) rates for the species jq refer to transformations of
c        that species into adjacent ionization states, not the
c        rates at which that species is formed.
c
c
c     principal physics references:
c
c     1) post, jensen, tarter, grasberger, lokke, 'steady state
c        radiative cooling rates for low-density high-temperature
c        plasmas', pppl-1352 (sub. to n.f.) (with corrections)
c     2) lokke and grasberger, 'xsnq-u a non-lte emission and
c        absorption coefficient subroutine', lll report ucrl-52276
c        (with corrections).
c     3) xsnq fortran listing
c
c
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
c
      zjlkev = 1.6021e-16
c
c
      zte = pte
c
      do 100 jq = 1 , nucz
c
      zchg = float(jq - 1)
      ivalnc = nvalnc(jq)
c
c
c
c     ******************** ionization rate ********************
c
c
c     calculate electron collisional ionization rates per shell
c     and sum for all populated shells. also compute
c     ionization energy loss rate cizlos.
c
      zisum = 0.
      zesum = 0.0
c
      do 20 jn = 1 , ivalnc
      zjn = jn
      zxn = ein(jq,jn)/zte
c
      zy = log10(zxn/4.)
      zfy = .23 + .046 * zy + .1074 * zy**2
     &          - .0459 * zy**3 - .01505 * zy**4
      zgaunt = 12.18 * exp(-zjn / (zjn+5.)) * (1. + .0335*zjn)
     &         * (1. - .622/(zchg+1.) - .0745/(zchg+1.)**2)
     &         * abs(zfy)
      zitemp = apn(jq,jn) * expund(-zxn,-45.) *
     &         (1.-expund(-zxn,-45.)) * zgaunt / ein(jq,jn)**2
      zisum = zisum + zitemp
      zesum = zesum + ein(jq,jn) * zitemp
   20 continue
c
      zisum = zisum * sqrt(zte)
      zesum = zesum * sqrt(zte)
      if(zisum .lt. 1.e-26) zisum = 0.
      if(zesum .lt. 1.e-26) zesum = 0.0
      rclion(jq) = pne * 3.44e-11 * zisum
      cizlos(jq) = pne * zjlkev * 3.44e-11 * zesum
c
c
  100 continue
c
c
      rclion(nucz + 1) = 0.
      cizlos(nucz + 1) = 0.0
c
c
      return
      end
      subroutine adbfi(knucz, pte, pne, pbfi, pbfil, kflag)
      external case4 
c
c     calculate the electron collisional ionization rate coefficients
c     using fitting formulas from the belfast group, culham
c     report clm-r216 (december 1981):
c
c          atomic and molecular data for fusion, part i
c     
c          recommended cross sections and rates for electron
c          ionisation of light atoms and ions
c
c          k l bell, h b gilbody, j g hughes, a e kingston, f j smith
c          the queen's university of belfast, belfast, n ireland
c
c     knucz = element atomic number ( z=1 to z=8 only)
c     pte = electron temperature (kev)
c     pne = electron density (cm-3)
c
c     pbfi = ionization rates for the z+1 charge states (sec-1)
c     pbfil = array of ionization energy loss rates (watts/ion)
c     kflag = 0 for o.k., = 1 for invalid element
c
      common / combfi / eh,     ah(6),    aah,     bth(3),
     &                  ehe(2), ahe(6,2), aahe(2), bthe(3,2),
     &                  eli(3), ali(6,3), aali(3), btli(3,3),
     &                  ebe(4), abe(6,4), aabe(4), btbe(3,4),
     &                  eb(5),  ab(6,5),  aab(5),  btb(3,5),
     &                  ec(6),  ac(6,6),  aac(6),  btc(3,6),
     &                  en(7),  an(6,7),  aan(7),  btn(3,7),
     &                  eo(8),  ao(6,8),  aao(8),  bto(3,8)
c
c
      dimension pbfi(*), pbfil(*)
c
c----------------------------------------------------------------------
c
      inspc = knucz + 1
c
      do 3 jspc = 1, inspc
      pbfi(jspc) = 0.0
      pbfil(jspc) = 0.0
    3 continue
c
      kflag = 0
c
c     error exit if element not available
c
      if(knucz .ge. 1 .and. knucz .le. 8) go to 5
      kflag = 1
      return
c
c     branch on element to call fit evaluation routine with proper
c     constants
c
    5 go to (10, 20, 30, 40, 50, 60, 70, 80) knucz
c
   10 call adbfit(knucz, pte, pne, eh, ah, aah, bth, pbfi, pbfil)
      return
   20 call adbfit(knucz, pte, pne, ehe, ahe, aahe, bthe, pbfi, pbfil)
      return
   30 call adbfit(knucz, pte, pne, eli, ali, aali, btli, pbfi, pbfil)
      return
   40 call adbfit(knucz, pte, pne, ebe, abe, aabe, btbe, pbfi, pbfil)
      return
   50 call adbfit(knucz, pte, pne, eb, ab, aab, btb, pbfi, pbfil)
      return
   60 call adbfit(knucz, pte, pne, ec, ac, aac, btc, pbfi, pbfil)
      return
   70 call adbfit(knucz, pte, pne, en, an, aan, btn, pbfi, pbfil)
      return
   80 call adbfit(knucz, pte, pne, eo, ao, aao, bto, pbfi, pbfil)
      return
      end
c
      
c******************
      BLOCK DATA case4
c
      common / combfi / eh,     ah(6),    aah,     bth(3),
     &                  ehe(2), ahe(6,2), aahe(2), bthe(3,2),
     &                  eli(3), ali(6,3), aali(3), btli(3,3),
     &                  ebe(4), abe(6,4), aabe(4), btbe(3,4),
     &                  eb(5),  ab(6,5),  aab(5),  btb(3,5),
     &                  ec(6),  ac(6,6),  aac(6),  btc(3,6),
     &                  en(7),  an(6,7),  aan(7),  btn(3,7),
     &                  eo(8),  ao(6,8),  aao(8),  bto(3,8)
c
c
c----------------------------------------------------------------------
c
c     hydrogen
c
      data eh / 13.60 /
      data ah / 2.3742e-8, -3.6866e-9, -1.0366e-8,
     &         -3.8010e-9,  3.4159e-9,  1.6834e-9 /
      data aah / 2.4617e-8 /
      data bth / 9.5986e-8, -9.2463e-7, 3.9973e-6 /
c
c     helium
c
      data ehe / 24.6, 54.42 /
      data ahe / 1.4999e-8, 5.6656e-10, -6.0821e-9,
     1          -3.5894e-9, 1.5529e-9,   1.3207e-9,
     2           3.4356e-9, -1.6865e-9, -6.9236e-10,
     2           9.7863e-11, 1.5591e-10, 6.2236e-11 /
      data aahe / 3.1373e-8, 3.0772e-9 /
      data bthe / 4.7893e-8, -7.7359e-7, 3.7366e-6,
     2            1.1902e-8, -1.1514e-7, 5.0489e-7 /
c
c     lithium
c
      data eli / 5.39, 75.64, 122.45 /
      data ali / 9.9655e-8, -5.5941e-8, -5.5228e-8,
     1           4.0589e-8, 1.4800e-8, -1.3120e-8,
     2           3.4023e-9, -7.6588e-10, -8.6078e-10,
     2           -8.9748e-10, 4.1661e-10, 3.3188e-10,
     3           1.1786e-9, -8.7637e-10, -9.3373e-11,
     3           2.1137e-10, 1.9017e-11, -4.0679e-11 /
      data aali / 4.5456e-8, 7.3504e-9, 1.9767e-9 /
      data btli / 2.7800e-7, -1.5830e-6, 5.4650e-6,
     2            5.3459e-10, -5.6387e-8, 2.9577e-7,
     3           -1.0926e-9, 8.8700e-10, 6.0764e-9 /
c
c     beryllium
c
      data ebe / 9.32, 18.21, 153.89, 217.71 /
      data abe / 7.4206e-8, -1.5520e-8, -3.9403e-8,
     1           7.2155e-9, 1.1098e-8, -2.5501e-9,
     2           1.7136e-8, -1.3997e-8, -2.9656e-11,
     2           3.0777e-9, -1.1676e-10, -5.1930e-10,
     3           1.6039e-9, -6.4336e-10, -7.7804e-10,
     3           3.3527e-10, 2.1889e-10, -1.0600e-10,
     4           4.9587e-10, -3.6870e-10, -3.9284e-11,
     4           8.8928e-11, 8.0007e-12, -1.7114e-11 /
      data aabe / 2.1732e-7, 2.7494e-8, 2.7873e-9, 8.3163e-10 /
      data btbe / -2.1649e-7, 2.8120e-7, 5.3031e-7,
     2            -1.5900e-8, 2.5911e-8, 2.4061e-8,
     3            -2.4333e-10, -9.7973e-9, 5.2094e-8,
     4            -4.5966e-10, 3.7317e-10, 2.5564e-9 /
c
c     boron
c
      data eb / 8.3, 25.15, 37.93, 259.37, 340.22 /
      data ab / 5.8365e-8, 1.0047e-8, -3.6230e-8,
     1         -7.3448e-9, 1.0220e-8, 1.6951e-9,
     2          2.0590e-8, -9.8899e-9, -6.0949e-9,
     2          2.7762e-9, 1.6499e-9, -6.7692e-10,
     3          5.7023e-9, -4.6578e-9, -9.8688e-12,
     3          1.0242e-9, -3.8855e-11, -1.7281e-10,
     4          7.3539e-10, -2.9498e-10, -3.5672e-10,
     4          1.5372e-10, 1.0036e-10, -4.8602e-11,
     5          2.5458e-10, -1.8930e-10, -2.0169e-11,
     5          4.5657e-11, 4.1076e-12, -8.7867e-12 /
      data aab / 3.0952e-7, 4.7980e-8, 9.1493e-9, 1.2780e-9,
     &             4.2697e-10 /
      data btb / -4.9240e-7, 1.3750e-6, -2.5387e-6,
     2           -4.1361e-8, 5.5259e-8, 9.9841e-8,
     3           -5.2910e-9, 8.6226e-9, 8.0067e-9,
     4           -1.1156e-10, -4.4920e-9, 2.3885e-8,
     5           -2.3600e-10, 1.9159e-10, 1.3125e-9 /
c
c     carbon
c
      data ec / 11.26, 24.38, 47.89, 64.49, 392.08, 489.98 /
      data ac / 5.9848e-8, 1.1903e-8, -3.0140e-8,
     1          -1.3693e-8, 8.3748e-9, 4.0150e-9,
     2          2.8395e-8, -1.6698e-8, -2.3557e-9,
     2          3.2161e-10, 9.6016e-10, 5.2713e-10,
     3          9.0555e-9, -6.3206e-9, -1.3256e-9,
     3          1.7441e-9, 3.2680e-10, -3.8303e-10,
     4          2.7464e-9, -2.0070e-9, -2.3595e-11,
     4          4.2011e-10, -8.1600e-11, -3.9729e-11,
     5          3.9495e-10, -1.5842e-10, -1.9158e-10,
     5          8.2555e-11, 5.3899e-11, -2.6102e-11,
     6          1.4715e-10, -1.0941e-10, -1.1657e-11,
     6          2.6389e-11, 2.3742e-12, -5.0786e-12 /
      data aac / 3.7442e-7, 6.0150e-8, 1.4501e-8,
     &           6.0330e-9, 6.8634e-10, 2.4679e-10 /
      data btc / -6.5826e-7, 2.0521e-6, -4.4694e-6,
     2           -4.0217e-8, -2.7908e-8, 5.5499e-7,
     3           -5.3596e-9, -1.0473e-8, 1.0629e-7,
     4           -5.7710e-9, 9.9302e-9, 7.4462e-9,
     5           -5.9916e-11, -2.4124e-9, 1.2828e-8,
     6           -1.3640e-10, 1.1074e-10, 7.5862e-10 /
c
c     nitrogen
c
      data en / 14.53, 29.60, 47.45, 77.47, 97.89, 552.06, 667.03 /
      data an / 4.6209e-8, 9.2264e-9, -1.2092e-8,
     1          -2.4852e-8, 5.1361e-9, 8.3068e-9,
     2          2.4369e-8, -2.2155e-9, -1.4805e-8,
     2          -4.4218e-10, 4.5211e-9, 1.7874e-10,
     3          1.2964e-8, -8.3408e-9, -2.3684e-9,
     3          2.2485e-9, 2.6234e-10, -2.6333e-10,
     4          4.6322e-9, -3.4645e-9, -3.1014e-10,
     4          8.0576e-10, 5.7791e-11, -1.4907e-10,
     5          1.5862e-9, -9.8633e-10, 7.5130e-11,
     5          3.1005e-11, -8.1970e-12, 2.6759e-11,
     6          2.3635e-10, -9.4805e-11, -1.1465e-10,
     6          4.9404e-11, 3.2255e-11, -1.5620e-11,
     7          9.2653e-11, -6.8892e-11, -7.3402e-12,
     7          1.6616e-11, 1.4949e-12, -3.1978e-12 /
      data aan / 2.7367e-7, 4.4690e-8, 1.0243e-8,
     &           7.9743e-9, 5.7812e-9, 4.1073e-10,
     &           1.5539e-10 /
      data btn / -4.2976e-7, 9.8352e-7, -9.5745e-7,
     2           3.0430e-8, -5.2696e-7, 2.4863e-6,
     3           3.4218e-8, -2.9155e-7, 1.2324e-6,
     4           -4.9184e-9, 6.3763e-9, 1.4917e-8,
     5           -8.5163e-9, 1.8527e-8, -7.8181e-9,
     6           -3.5856e-11, -1.4437e-9, 7.6765e-9,
     7           -8.5888e-11, 6.9728e-11, 4.7767e-10 /
c
c     oxygen
c
      data eo / 13.62, 35.12, 54.93, 77.41, 113.90, 138.12, 739.32,
     &          871.39 /
      data ao / 3.3559e-8, 1.3449e-8, -6.7112e-9,
     1          -1.9976e-8, 1.6214e-9, 6.5852e-9,
     2          2.4476e-8, -5.3141e-9, -7.3316e-9,
     2          -4.4515e-9, 2.4257e-9, 1.9791e-9,
     3          1.4741e-8, -8.7905e-9, -8.6099e-10,
     3          -2.4143e-10, 1.2598e-10, 6.4901e-10,
     4          6.2130e-9, -2.5047e-9, -3.0813e-9,
     4          1.3559e-9, 8.6816e-10, -4.3189e-10,
     5           2.6145e-9, -2.0276e-9, -1.6569e-10,
     5          5.0245e-10, 3.0067e-11, -9.8231e-11,
     6          1.0099e-9, -6.5165e-10, 2.8863e-12,
     6          3.0336e-11, -1.4065e-11, 4.5106e-12,
     7          1.5258e-10, -6.1203e-11, -7.4014e-11,
     7          3.1894e-11, 2.0823e-11, -1.0084e-11,
     8          6.2090e-11, -4.6167e-11, -4.9189e-12,
     8          1.1135e-11, 1.0018e-12, -2.1430e-12 /
      data aao / 3.2721e-7, 4.9003e-8, 1.7523e-8,
     &           1.0270e-8, 4.0023e-9, 2.6228e-9,
     &           2.6516e-10, 1.0413e-10 /
      data bto / -6.7484e-7, 2.3938e-6, -6.04388e-6,
     2           2.1747e-8, -5.4612e-7, 2.7213e-6,
     3           2.8129e-8, -3.1692e-7, 1.4381e-6,
     4           4.8134e-10, -4.3936e-8, 2.1924e-7,
     5           -1.5957e-9, -7.2485e-10, 1.9724e-8,
     6           -2.6662e-9, 2.8968e-9, 1.7246e-8,
     7           -2.3148e-11, -9.3201e-10, 4.9557e-9,
     8           -5.7557e-11, 4.6727e-11, 3.2010e-10 /
      end


c#####################################################################
c@adbfit
c
c  dps 07-aug-89 15.12 Remove truncation of rates below ionization
c                potential / 10.
c*********************************************************************
c
      subroutine adbfit(knucz, pte, pne, pe, pa, paa, pbt, pbfi, pbfil)
c
c     evaluate the belfast group rate coefficient fitting formulas.
c
      dimension pe(knucz), pa(6,knucz), paa(knucz), pbt(3,knucz)
      dimension pbfi(*), pbfil(*)
c
      zjlev = 1.6021e-19
c
c     convert te to ev and zero return arrays
c
      zte = 1000.0 * pte
c
      inspc = knucz + 1
      do 10 jspc = 1, inspc
      pbfi(jspc) = 0.0
      pbfil(jspc) = 0.0
   10 continue
c
c     loop over all species except fully stripped
c
      do 100 jspc = 1, knucz
c
c..15.12 Let it calculate rates for low Te => comment out statement
c
c     rate is zero for kt < (ip/10)
c
c      if(zte .lt. pe(jspc)/10.0) go to 100
c
c     branch if high temperature fit is required
c
      if(zte .gt. 10.0 * pe(jspc)) go to 80
c
c     (ip/10) < te < (10*ip) fitting form
c
      zx = pe(jspc) / zte
      zlog = log10(1.0 / zx)
c
      zsum = pa(1,jspc)
      do 20 j = 2, 6
      zsum = zsum + pa(j,jspc) * (zlog**(j-1))
   20 continue
c
      zirc = exp(-zx) * sqrt(1.0 / zx) * zsum
      pbfi(jspc) = pne * zirc
      pbfil(jspc) = zjlev * pe(jspc) * pbfi(jspc)
      go to 100
c
c     high temperature fit for te > 10*ip
c
   80 zx = pe(jspc) / zte
c
      zirc = sqrt(zx) * (paa(jspc) * log(1.0/zx)
     &     + pbt(1,jspc) + pbt(2,jspc) * zx + pbt(3,jspc) * zx * zx)
c
      pbfi(jspc) = pne * zirc
      pbfil(jspc) = zjlev * pe(jspc) * pbfi(jspc)
c
  100 continue
c
      return
      end
      subroutine adygr(knucz, pte, pne, pygri, pygrl, kflag)
      external case5
c
c     compute direct electron ionization rates for hydrogen-like
c     through neon-like charge states of certain heavy elements
c     (presently scandium and iron).  uses formulae due to
c     s. younger, nbs.  (jqsrt 29, 61 (1983) and private
c     communications).
c
c     knucz= element (z=21 or 26 for now)
c     pte = electron temperature(kev)
c     pne = electron density (cm-3)
c
c     pygri = array of ionization rates (sec-1)
c             (only the 10 states from h-like to ne-like will be
c              overwritten by the subroutine)
c     pygrl = array of ionization energy loss rates (watts/ion)
c     kflag = 0 for o.k.;
c           = 1 for invalid element (no values are returned)
c
      common / comygr / eion(10,3), eionfe(10,3), eionsc(10,3),
     &          nss(10,3), rate(10,3),
     &          ca(3), cb(3), cc(3), cd(3)
c
      dimension pygri(*), pygrl(*)
c
c
      tev = 1000. * pte
c
      eion(1,1) = 0.0
      do 20 jss = 1, 3
      do 20 jbe = 1, 10
      if(knucz .eq. 21) eion(jbe,jss) = eionsc(jbe,jss)
      if(knucz .eq. 26) eion(jbe,jss) = eionfe(jbe,jss)
   20 continue
c
      kflag = 0
      if(eion(1,1) .ne. 0.0)go to 50
      kflag = 1
      return
c
   50 do 100 jbe = 1, 10
      ics = knucz - jbe +1
      pygri(ics) = 0.0
      pygrl(ics) = 0.0
c
      do 90 jss = 1, 3
      rate(jbe,jss) = 0.0
      if(nss(jbe,jss) .eq. 0) go to 90
c
      chi = tev / eion(jbe,jss)
      a = ca(jss) * float(nss(jbe,jss))
      b = cb(jss) * float(nss(jbe,jss))
      c = cc(jss) * float(nss(jbe,jss))
      d = cd(jss) * float(nss(jbe,jss))
      rate(jbe,jss) = fchi(chi,a,b,c,d) * (eion(jbe,jss)**(-1.5)) *
     &                2.2e-06 * sqrt(chi) * exp(-1.0/chi)
      pygri(ics) = pygri(ics) + pne * rate(jbe,jss)
      pygrl(ics) = pygrl(ics) + 1.6021e-19 * pne * eion(jbe,jss)
     &              * rate(jbe,jss)
   90 continue
  100 continue
c
      return
      end
c******************
      BLOCK DATA case5
c
      common / comygr / eion(10,3), eionfe(10,3), eionsc(10,3),
     &          nss(10,3), rate(10,3),
     &          ca(3), cb(3), cc(3), cd(3)
c
c
      data ((nss(jbe,jss),jbe=1,10),jss=1,3)/ 1,9*2,0,0,1,7*2,
     &      4*0,1,2,3,4,5,6/
c
      data (ca(j),j=1,3)/7.5e-14, 5.9e-14, 8.4e-14/
      data (cb(j), j=1,3)/ -2.605e-14, -1.635e-14, -2.68e-14/
      data (cc(j), j=1,3)/ 1.195e-14, 0.82e-14, 0.672e-14/
      data (cd(j), j=1,3)/ -6.15e-14, -3.79e-14, -5.08e-14/
c
      data ((eionsc(jbe,jss),jbe=1,10),jss=1,3)/ 6034., 5675.,
     &     5566., 5491., 5318., 5204., 5095., 4990., 4890., 4794.,
     &     0., 0., 1288., 1213., 1135., 1068., 993.3,
     &     914.0, 846.3, 769.3, 4*0., 1090., 1004., 922.4,
     &     838.6, 761.8, 687.4 /
c
      data ((eionfe(jbe,jss),jbe=1,10),jss=1,3)/ 9277., 8829.,
     &     8694., 8521., 8400., 8263., 8118., 7949., 7818., 7824.,
     &     0., 0., 2046., 1950., 1852., 1755., 1679., 1566., 1482.,
     &     1394., 4*0., 1789., 1678., 1583., 1463., 1366., 1266./
c
      end
c---------------------------------
c
      function fchi(pchi, pa, pb, pc, pd)
c
      z1 = pa + pb * (1.0 + 1.0 / pchi)
      z2 = pc - (pa + 2.0*pb + pb/pchi)/pchi
      z3 = pd / pchi
      zchisq = pchi * pchi
      zchicb = zchisq * pchi
      zalpha = (.001193 + .9764*pchi +.6604*zchisq +
     &          .02590*zchicb) / (1.0 + 1.488*pchi +
     &          .2972*zchisq + .004925*zchicb)
      zbeta = (-.0005725 + .01345*pchi +.8691*zchisq +
     &         .03404*zchicb) / (1.0 + 2.197*pchi +
     &         .2454*zchisq + .002053*zchicb)
      fchi = (3.0e13 / pchi) * (z1 + z2 * zalpha + z3 * zbeta)
      return
      end
      subroutine adrrec(pte, pne)
c
c     rev: 10/10/82 for ad package
c
c     rev: 9/9/81, 6/15/78
c
c     compute radiative recombination rates rrarec(jq) for a given
c     electron temperature pte (kev), where jq = (ionic charge + 1).
c     these rates have the units (sec-1).
c     also computed are the radiation rates associated with
c     this process, radrrc(jq).  these rates are in joules/sec.
c
c     notes:
c
c     1) electrons are assumed to fill ionic shells in simple ascending
c        order.  this is assumed both explicitly (the apn array)
c        as well as implicitly through the use of nvalnc values
c        in loop indices, etc.
c     2) through use of the function expund, exp(x) is set equal to
c        zero in the rate equations for values of x < -45.
c     3) rates for the species jq refer to transformations of
c        that species into adjacent ionization states, not the
c        rates at which that species is formed.
c     4) upon evaluation of the formulae it will be found that the
c        formalism used in xsnq/lokke and grasberger (refs 6,2)
c        in fact leads to the same equations given in
c        pppl-1352/seaton (refs 1,3) if the shielded charge for
c        level n is used for 'q'.
c
c
c     principal physics references:
c
c     1) post, jensen, tarter, grasberger, lokke, 'steady state
c        radiative cooling rates for low-density high-temperature
c        plasmas', corrected pppl-1352 (sub. to n.f.)
c     2) lokke and grasberger, 'xsnq-u a non-lte emission and
c        absorption coefficient subroutine', lll report ucrl-52276
c     3) seaton, 'radiative recombination of hydrogenic ions',
c        mnras:119,81 (1959)
c     4) cox and tucker, 'ionization equilibrium and radiative
c        cooling of a low-density plasma', ap.j.:157,1157 (1969)
c     5) grasberger, xsnq - tokamak memorandum #1, 2/11/76
c     6) xsnq fortran listing
c
c
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
      zjlkev = 1.6021e-16
c
c
c
      zte = pte
c
      rrarec(1) = 0.
      radrrc(1) = 0.
c
c     ********** loop for all non-neutral species **********
c
      do 100 jq = 2 , nspc
c
      ivalnc = nvalnc(jq)
      ivjqm1 = nvalnc(jq - 1)
c
c
c
c     calculate radiative recombination rates per shell and sum them
c     (sums (fractional vacancy of shell) * (alphan) ; alphan
c     from corrected pppl-1352).  also sum continuum and line
c     radiation rates due to radiative recombinations to each shell.
c
c
c
      zrsum = 0.
      zrsumc = 0.
      zrsuml = 0.
c
      do 31 jn = ivalnc,10
c
c     skip valance shell if full to avoid undefined en, qn
c
      if(aqn(jq,jn) .eq. 0.0) go to 31
c
      zxn = en(jq,jn)/zte
c
      zrtemp = aqn(jq,jn) * qn(jq,jn) * (zxn**1.5) * expe1(zxn)
c
c     high-lying levels correction for nlvlmx=10 is 5.13
c
      if(jn .eq. 10) zrtemp = zrtemp * 5.13
      zrsum = zrsum + zrtemp
c
      zrtmpl = enm(jq-1,ivjqm1,jn) * zrtemp
      zrsuml = zrsuml + zrtmpl
c
      zrtmpc = aqn(jq,jn) * qn(jq,jn) * sqrt(zxn) * en(jq,jn)
      zrsumc = zrsumc + zrtmpc
   31 continue
c
      rrarec(jq) = pne * 5.20e-14 * zrsum
      if(zrsumc .lt. 1.e-07) zrsumc = 0.
      zrsumc = 8.32e-30 * zrsumc
c
      if(zrsuml .lt. 1.e-7) zrsuml = 0.
      zrsuml = zjlkev * 5.20e-14 * zrsuml
      radrrc(jq) = pne * (zrsumc + zrsuml)
c
c
  100 continue
c
c
      return
      end
      subroutine addrec(pte, pne)
c7600      lcm dierec
c
c     rev: 10/10/82 for new ad package
c
c     rev: 11/11/81 add cdnn and cdnm multipliers
c
c     rev: 10/23/78 (remove 'c' from zdsum0 defining statement to fix
c                    n-n dielectronic rate bug)
c     rev: 6/15/78
c
c     computes dielectronic recombination rate rdirec(jq)
c     for a given electron temperature pte (kev) and electron
c     density pne (cm-3), where jq = (ionic charge +1).
c     these rates have the units (sec-1).  since
c     the rate coeff depend on ne, the total recombination rate
c     per ion is not strictly a linear function of ne.  however,
c     the dielectronic collisional-interruption correction is
c     only weakly dependent on the density, and thus the non-linearity
c     is generally small.
c
c     also computed are the radiation rates associated with this
c     process, raddrc(jq), which include contributions from
c     the radiative decays of both the excited and captured electrons.
c     these rates are in joules/sec.
c
c     notes:
c
c     1) electrons are assumed to fill ionic shells in simple ascending
c        order.  this is assumed both explicitly (the apn array)
c        as well as implicitly through the use of nvalnc values
c        in loop indices, etc.
c     2) through use of the function expund, exp(x) is set equal to
c        zero in the rate equations for values of x < -45.
c     3) rates for the species jq refer to transformations of
c        that species into adjacent ionization states, not the
c        rates at which that species is formed.
c
c
c     principal physics references:
c
c     1) post, jensen, tarter, grasberger, lokke, 'steady state
c        radiative cooling rates for low-density high-temperature
c        plasmas', corrected pppl-1352 (sub. to n.f.)
c     2) burgess, 'a general formula for the estimation of
c        dielectronic recombination coefficients in low-density
c        plasmas', ap.j. 141, 1588 (1965).
c     3) grasberger, xsnq - tokamak memorandum #1, 2/11/76
c     4) grasberger, xsnq - tokamak memorandum #4, 7/29/76
c     5) grasberger, xsnq - tokamak memorandum #10, 1/5/77
c     6) merts, cowan, magee, 'the calculated power output
c        from a thin iron-seeded plasma', los alamos report la-6220-ms
c     7) xsnq fortran listing
c     8) ansari, elwert and mucklich, 'on dielectronic recombination',
c        zeitschrift naturforschung 25a, 1781 (1970).
c
c
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
      common /addiel/ cdnn(100), cdnm(100), ldrmlt, yhnna, yhnnb,
     &                yhnnc, yhnma, yhnmb, yhnmc

c
      zjlkev = 1.6021e-16
c
c     call addrcc to fill cdnn and cdnm multipicative constant arrays
c
      call addrcc
c
      zte = pte
      zne = pne
c
      znuc = float(nucz)
c
      rdirec(1) = 0.
      raddrc(1) = 0.
      rdirec(nspc) = 0.
      raddrc(nspc) = 0.
c
      if(nucz .eq. 1) return
c
c
c     ********** loop for all species except neutral **********
c              ********** and fully ionized **********
c
      do 100 jq = 2 , nucz
c
c
      zchg = float(jq-1)
      ivalnc = nvalnc(jq)
c
c
c
c     ************* dielectronic recombination rate *************
c
c     (equations and variables follow pppl-1352/burgess formalism)
c
c
c
      zq = zchg
      zqp1 = zchg + 1
c
      zbq = sqrt(  (zq * (zq+1.)**5) / (zq*zq + 13.4)  )
      znt = (4.77e+18 * (zq**6) * sqrt(zte) / zne)**0.142857143
      zdnn = .0068 * znt * (1.+.001*znt) / (1.+.00002*znt*znt)
      zzz = .0015 * ( (zq+1.) * znt )**2
      zdnm = zzz / (1. + zzz)
      zzznn = 7.59e-14 * zbq * zdnn / (zte**1.5)
      zzznm = 7.59e-14 * zbq * zdnm / (zte**1.5)
      zzza = 1. + .015 * (zq**3) / ( (zq+1.)**2 )
      zzzet = ( (zq + 1.)**2 ) / (73.5 * zzza * zte)
c
c
c     compute n - n contribution to sum first
c
c
      zeijnn = enn(jq) / (.0136 * zqp1 * zqp1)
      zynn = (zq + 1.) * zeijnn
      zaynn = sqrt(zynn) / (1. + .105 * zynn + .015 * zynn**2)
c
      zdsum0 = fnn(jq) * zaynn * expund(-zzzet*zeijnn , -45.)
c
c
c     sum over n - m transitions (note that transitions from
c     inner shells to the valence shell and above are included)
c
c
      zdsum = 0.
      zdsumr = 0.
c
      do 40 jn = 1 , ivalnc
      iijm = ivalnc
      if(jn .eq. ivalnc .or.
     &          aqn(jq,ivalnc) .eq. 0.) iijm = ivalnc + 1
c
      do 40 jm = iijm , 10
c
      zeij = enm(jq,jn,jm) / (.0136 * zqp1**2)
      zy = (zq + 1.) * zeij
      zay = sqrt(zy) / (2. + .420 * zy + .06 * zy**2)
c
      zdtemp = fnm(jq,jn,jm) * zay * expund(-zzzet*zeij , -45.)
c
      zdsum = zdsum + zdtemp
      zdsumr = zdsumr + enm(jq,jn,jm) * zdtemp
   40 continue
c
      zdnnrt = zzznn * zdsum0
      zdnmrt = zzznm * zdsum
cw    write(59,909)jq,zdnnrt,zdnmrt
cw  909 format(1x,'jq =',i3,3x,'zdnnrt =',1pe12.2,'  zdnmrt =',e12.2)
      rdirec(jq) = pne * (cdnn(jq) * zdnnrt + cdnm(jq) * zdnmrt)
c
      itemp = nvalnc(jq-1)
      ziq = en(jq,itemp)
      ztemp = cdnn(jq) * (ziq + enn(jq)) * zdnnrt +
     &              cdnm(jq) * (ziq * zdnmrt + zzznm * zdsumr)
      if(ztemp .lt. 1.e-21) ztemp = 0.
      raddrc(jq) = pne * zjlkev * ztemp
c
c
c
  100 continue
c
c
      return
      end
      subroutine addrcc
c
c     rev: 10/10/82 for new ad package
c
c     setup the dielectronic recombination rate coefficient
c     multiplicative constant arrays cdnn and cdnm according
c     to the ldrmlt switch.
c     
c     ldrmlt = 0     use cdnn, cdnm arrays as given
c     ldrmlt = 1     cdnn, cdnm arrays set to 1.0
c     ldrmlt = 2     use y. hahn factors to setup cdnn, cdnm
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / addiel / cdnn(100), cdnm(100), ldrmlt,
     &                  yhnna, yhnnb, yhnnc, yhnma, yhnmb, yhnmc
c
c
      if(ldrmlt .eq. 0) return
c
      do 20 jspc = 1, 100
      cdnn(jspc) = 1.0
      cdnm(jspc) = 1.0
   20 continue
c
      if(ldrmlt .eq. 1) return
c
c
      do 100 jspc = 2, nucz
c
      znelec = float(nucz - jspc + 1)
c
      cdnn(jspc) = yhnna * exp(-yhnnb*sqrt(znelec) ) *
     &                  (float(nucz)**yhnnc)
c
      cdnm(jspc) = yhnma * exp(-yhnmb*sqrt(znelec)) *
     &                  (float(nucz)**yhnmc)
  100 continue
c
      return
      end
      subroutine adecex(pte, pne)
      external case6
c
c     rev: 10/10/82 for new ad package
c
c     rev: 6/15/78
c
c     compute radiation rates radclx(jq) (joule/sec) due to
c     collisional excitation for each ionic species jq = (ionic
c     charge + 1).  all excitations are assumed to result in
c     radiative decays.
c
c     references:
c     1) post et al. pppl-1352 (corrected)
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
      common / cexdat / zdnn(28), zgamnn(28)
c
      zjlkev = 1.6021e-16
c
c
c
      zte = pte
      znuc = float(nucz)
c
c     ******** loop for all species except fully ionized ********     
c
      do 100 jq = 1 , nucz
c
      zchg = float(jq - 1)
      ivalnc = nvalnc(jq)
      ibound = nucz - jq + 1
c
c     ***** n - n rate *****
c
      zznn = 0.
      if(fnn(jq) .eq. 0.) go to 30
c
      zxnn = enn(jq) / zte
c
      if(jq .ne. 1) go to 20
      ztemp = .06 * (sqrt(zxnn) - 2.) / (1. + zxnn)
      go to 22
c
   20 if(ibound .le. 28) go to 21
      ztemp = .7 * (1. - .8 / zchg)
      go to 22
c
   21 ztemp = zdnn(ibound) * (1. - zgamnn(ibound) / zchg)
c
   22 zgnn = ztemp + .276 * expe1(zxnn)
c
      zznn = fnn(jq) * expund(-zxnn , -45.) * zgnn
c
c     ***** sum n - m rates from all occupied shells *****
c
   30 znmsum = 0.
c
c
      do 31 jn = 1 , ivalnc
c
      zjn = float (jn)
      iijm = ivalnc
      if(jn .eq. ivalnc .or.
     &        aqn(jq,ivalnc) .eq. 0.) iijm = ivalnc + 1
c
      do 31 jm = iijm , 10
      zjm = float(jm)
c
      zxnm = enm(jq,jn,jm) / zte
      ztemp = zjm * (zjm - zjn) * (1. + (1. - 2./znuc) * zxnm) / 20.
      zgnm = .19 * (1. + .9 * (1. + ztemp) * expe1(zxnm))
c
c     note that fnm has apn and aqn factors already included
c
      zznm = fnm(jq,jn,jm) * expund(-zxnm , -45.) * zgnm
      znmsum = znmsum + zznm
   31 continue
c
c
c
      ztemp = (zznn + znmsum) / sqrt(zte)
      if(ztemp .lt. 1.e-10) ztemp = 0.
      radclx(jq) = pne * zjlkev * 4.99e-10 * ztemp
c
  100 continue
c
c     ***************end of jq species loop ***************
c
      radclx(nspc) = 0.
c
      return
      end
c******************
      BLOCK DATA case6
c
      common / cexdat / zdnn(28), zgamnn(28)
c
c
      data (zdnn(j),j=1,28) / 0., 0.,
     & 1.2, .91, .89, .92, .86, .87, .85, 0.,
     & .72, .78, .76, .78, .75, .70, .68, .65, .65,
     & 8 * .70, 0. /
c
      data (zgamnn(j),j=1,28) / 0., 0.,
     & .54, .77, .58, .57, .41, .63, .66, 0.,
     & .97, .96, .88, .86, .87, .85, .83, .81,
     & 9 * .80, 0. /
c
      end
c---------------------------------------------
      subroutine adbrem(pte, pne)
c
c     rev: 10/10/82 for new ad package
c
c     rev: 6/15/78
c
c     input : pte, nucz, nspc
c     output : radbrm
c
c     compute bremsstrahlung radiation radbrm(jq) (joule/sec)
c     for each ionic species jq = (ionic charge + 1).
c
c     references:
c     1) post et al. pppl-1352 (corrected)
c     2) karzas and latter
c     3) tucker
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10), fnn(100), enn(100)
c
      common / adelec / rclion(100), rrarec(100), rdirec(100),
     &                  cizlos(100), radrrc(100), raddrc(100),
     &                  radclx(100), radbrm(100)
c
c
      zte = pte
c
c     use fixed free-free gaunt factor zgff = 1.2
c
      zgff = 1.2
c
      do 100 jq = 1 , nspc
      zchg = float(jq - 1)
      radbrm(jq) = pne * 4.85e-31 * zchg * zchg * sqrt(zte) * zgff
  100 continue
c
      return
      end
      function expe1(px)
c
c     rev: 6/15/78
c
c     calculate exp(px) * e1(px), where the exponential integral
c     e1(x) = integral from x to infinity of (exp(-t)/t)dt.
c
c     reference: 'handbook of mathematical functions', abramowitz
c     and stegun, national bureau of standards
c
      if(px.ge.1.)go to 2
      if(px.gt.1.e-4)go to 1
c
      expe1 = exp(px) * ( - log(px) - .57721566 + px )
      return
c
c
    1 expe1 = exp(px) * ( - log(px) - .57721566 + .99999193*px
     &         - .24991055*px**2 + .05519968*px**3 - .00976004*px**4
     &         + .00107857*px**5 )
c
      return
c
c
    2 expe1 = (px + 2.334733 + .250621/px) /
     &         (px**2 + 3.330657*px + 1.681534)
c
      return
      end
c
c
c********************************************************************
c
c
      function expund(px,pchk)
c
c     rev: 6/15/78
c
c     prevent underflows in rate calculations by returning
c     0. for exp(px) if px < pchk.
c
      if( px - pchk ) 1 , 1 , 2
    1 expund = 0.
      return
c
    2 expund = exp(px)
      return
c
      end
      subroutine adbcxr(pnneut, pvneut, kncxb, kcxopt, kerr, kvunit)
c
c     rev: 10/10/82 for new ad package
c
c     rev: 7/24/79; zero radbcx before kncxb=0 test & return
c     rev: 5/10/79; add radiation due to beam charge exchange via
c                   radbcx and beqrad
c
c     compute recombination rate due to charge exchange with neutral
c     beam ions.
c
c          npsc = number of impurity species (nuclear z + 1)
c
c          pnneut = array of angular averaged neutral densities for
c                   each beam energy component
c          pvneut = array of velocities for each beam component (cm/s)
c          kncxb = number of energy components in the beam
c                 (kncxb = 0 turns off rcxrec)
c          kcxopt = 1, 2, or 3 to select cross-sections to be used
c          rcxrec = array of charge exchange recombination rates for
c                   each species, sec-1 units (e.g. rate per ion)
c          kerr = error flag from adcxxs routine (0 = o.k.)
c          radbcx = radiation rate due to cascade after
c                   charge exchange to excited level (watts / ion)
c
c
      common / adsdat / nucz, nspc, apn(100,10), aqn(100,10),
     &                  nvalnc(100), sigma(10,10), qn(100,10),
     &                  en(100,10), ein(100,5), enm(100,5,10),
     &                  fnm(100,5,10),fnn(100),enn(100)
c
      common / adneut / rcxrec(100), radbcx(100)
c
      dimension pnneut(*), pvneut(*)
      dimension zsigma(100)
c
      zjlkev = 1.602e-16
c
      inspc = nspc
c
c     zero rate arrays for accumulation of contributions from each beam
c     component
c
      do 10 jspc = 1 , 100
      rcxrec(jspc) = 0.
      radbcx(jspc) = 0.
   10 continue
c
c     kncxb = 0 or kcxopt = 0 to turn off c-x recombination
c
      if(kncxb .eq. 0 .or. kcxopt .eq. 0) return
c
c     obtain cross-sections for each velocity (energy) component
c     in the beam and compute total recombination rate.
c
      do 60 jb = 1 , kncxb
c
      call adcxxs(pvneut(jb), kvunit, zvcms, zsigma, inspc,
     &                kcxopt, kerr)
      znflux = zvcms * pnneut(jb)
c
      do 50 jspc = 2 , inspc
      rcxrec(jspc) = rcxrec(jspc) + znflux * zsigma(jspc)
   50 continue
c
   60 continue
c
      do 70 jspc = 2 , inspc
      itemp = nvalnc(jspc-1)
      radbcx(jspc) = zjlkev * en(jspc,itemp) * rcxrec(jspc)
   70 continue
c
      return
      end
      subroutine adcxxs(pvel, kvunit, pvcms, psigma, knspc,
     &                     kopt, kerr)
      external case7
c
c     rev: 10/10/82 for new ad package
c
c     supplies cross-sections for charge exchange between neutral
c     hydrogen and multi-charged ions.  the variable kopt allows
c     the selection of results from one of the various published
c     calculations.
c
c     pvel =  relative velocity, units as selected by kvunit
c     kvunit = units flag for input velocity
c                 1 = cm/sec
c                 2 = kev/amu
c                 3 = atomic units
c     psigma = array of cross-section values returned, cm**2.
c              psigma(i) = cross-section for ion with charge = i-1,
c              for 2 <= i <= knspc
c     knspc = highest species number, = max z + 1
c     kopt = switch for different results, as follows:
c            1 = olson & salop, phys. rev. a 14, 579 (1976)
c                (absorbing sphere, v <= 1.e08 cm/sec)
c            2 = grozdanov & janev, phys. rev. a 17, 880 (1978)
c                (tunneling, 2.e06 <= v <= 8.e08 cm/sec)
c            3 = olson & salop, phys. rev. a 16, 531 (1977)
c                (classical, (2 - 7)e08 cm/sec)
c     kerr = error flag. 0 = o.k., <0 = invalid input, others in text
c
c*********************************************************************
c
      dimension psigma(*)
      common / cxxdat / zfitv(4), zfitlt(4), zfitgt(4)
c---------------------------------------------------------------------
c
      kerr = 0
c
c     convert input velocity parameter
c
      if(kvunit .ne. 1) go to 2
      zvcms = pvel
      zvkevm = (zvcms / 4.39e07)**2
      zvau = zvcms / 2.19e08
      go to 5
c
    2 if(kvunit .ne. 2) go to 3
      zvkevm = pvel
      zvcms = 4.39e07 * sqrt(zvkevm)
      zvau = zvcms / 2.19e08
      go to 5
c
    3 if(kvunit .ne. 3) go to 4
      zvau = pvel
      zvcms = 2.19e08 * zvau
      zvkevm = (zvcms / 4.39e07)**2
      go to 5
c
    4 kerr = - 1
      return
c
c     return velocity in cm/s units
c
    5 pvcms = zvcms
c
c     neutral species cross-section = 0
c
      psigma(1) = 0.
c
      go to (100, 200, 300), kopt
      kerr = -2
      return
c
c--------------------------------------------------------------------
c
c     kopt = 1 ; olson & salop low - v absorbing sphere
c     nominally valid only for v < 1.0e08 cm/sec
c
  100 do 110 jspc = 3 , knspc
c
cw    write(5,198)jspc
cw  198 format(1x,/,'jspc =',i3)
c
      zchg = float(jspc - 1)
      zl = 2.648 / sqrt(zchg)
      zc = 2.864e-04 * zchg * (zchg - 1.) * zvau
      if(zc .lt. .541341133/(zl**2) ) go to 102
      kerr = 1
      return
c
c     iterative solution to equation (10)
c
c     (first guess, from linear fit to fig. 3)
c
  102 zr = 9. + 11. * zchg / 45.
      zr = - log( zc / (zr**2) ) / zl
      zfr = zr * zr * exp(- zl * zr)
      icount = 0
c
  105 if(zr .lt. 2./zl) go to 180
      if(icount .gt. 100) go to 181
      icount = icount + 1
c
cw    write(5,199)zr
cw  199 format(1x,'zr=',f10.2)
c
      zr = zr + (zc/zfr - 1.) / (2./zr - zl)
      zfr = zr * zr * exp (- zl * zr)
c
      if( abs(zfr - zc)/zc .ge. .001) go to 105
c
c     iteration complete to 0.1 0n velocity; compute sigma
c
      psigma(jspc) = 3.14159 * (0.529e-08 * zr)**2
  110 continue
c
c     fill in z = 1 cross section with z = 2 value
c
      psigma(2) = psigma(3)
      return
c
  180 kerr = 2
      return
  181 kerr = 3
      return
c
c----------------------------------------------------------------
c
c     kopt = 2; grozdanov & janev
c     results shown only for z = 10, 20, and 30, 2.e06 < v < 8.e08
c     cm/sec; all other points by interpolation or extrapolation
c
  200 do 210 jspc = 2 , knspc
      zchg = float(jspc - 1)
c
c     empirical fit to figure 3, centered at v = 2.0e08 cm/sec
c
      zslope = - (.25 + .04 * zchg) * 1.0e-14
      z0 = .11 * zchg * 1.0e-14
      psigma(jspc) = max(0., z0 + zslope * log10(zvcms/2.e08))
c
  210 continue
      return
c
c----------------------------------------------------------------------
c
c
c     kopt = 3 ; olson & salop classical trajectory.  uses a
c     power law fit to fig. 1 at 4 velocities in the range
c     2.7e08 < v < 5.e08 cm/sec.  input velocities outside this
c     range yield data for the appropriate limiting velocity.  the
c     last given data curve is for z=36, hence data for z>36
c     represents an extrapolation.
c
  300 continue
      z3v = max(2.7e08, zvcms)
      z3v = min(5.0e08, z3v)
c
      ix = int(z3v/1.e08) - 1
      if(ix .gt. 3) ix = 3
      ixp1 = ix + 1
c
c     use linear fit to z=8 curve to find x-section for z=8 at given
c     velocity
c
      zsig8 = 7.43e-15 - 1.26e-23 * z3v
c
c     interpolate to find power law exponents for this velocity; one
c     for z<8, the other for z>8.
c
      zalt = zfitlt(ix) + (zfitlt(ixp1) - zfitlt(ix)) *
     &                 ( (z3v-zfitv(ix)) / (zfitv(ixp1) - zfitv(ix)) )
      zagt = zfitgt(ix) + (zfitgt(ixp1) - zfitgt(ix)) *
     &                 ( (z3v-zfitv(ix)) / (zfitv(ixp1)-zfitv(ix)) )
c
c     fit parameters found for this velocity; compute x-sections
c
      do 315 jspc = 2 , knspc
      zdiv8 = float(jspc - 1) / 8.
      if(zdiv8 .ge. 1.) go to 310
      psigma(jspc) = zsig8 * zdiv8**zalt
      go to 315
  310 psigma(jspc) = zsig8 * zdiv8**zagt
  315 continue
c
      return
c
c---------------------------------------------------------------------
c
      end
c******************
      BLOCK DATA case7
c
      common / cxxdat / zfitv(4), zfitlt(4), zfitgt(4)
c
      data zfitv/ 2.7e08, 3.0e08, 4.0e08, 5.0e08/
      data zfitlt/1.220988, 1.277249, 1.868400, 2.575095/
      data zfitgt/1.000000, 1.050777, 1.221593, 1.668531/
c
      end

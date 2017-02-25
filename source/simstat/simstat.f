c@simstat.f  22:00 27-Feb-96  Glenn Bateman
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c    simstat.f computes statistics from temperature and density profiles
c  from BALDUR transport simulations and experimental data.
c
c
      integer kr
c
      parameter ( kr = 130 )
c
      character cdate*9, ctime*8, csimfile*32, cexpfile*32, ctemp*32
      character cline*132
c
c  cexpfile = name of the file containing experimental data
c  csimfile = name of the file containing BALDUR simulation output
c  ctemp    = temporary character string for file name
c  cline    = temporary character string for line of input
c
c..experimental data points
c
c    major radii (m), temperatures (eV) or density (m^{-3}),
c    and error bars
c
      dimension  rmteexp(kr), teexp(kr), teerr(kr), terelerr(kr)
     &         , rmtiexp(kr), tiexp(kr), tierr(kr), tirelerr(kr)
     &         , rmenexp(kr), enexp(kr), enerr(kr), enrelerr(kr)
c
c  rmteexp(j) = major radius of experimental T_e data points (m)
c  teexp(j)   = T_e electron temperature data points (keV)
c  teerr(j)   = standard deviations error bars for T_e (keV)
c  terelerr(j) = teerr(j) / ( maximum value of T_e ), ...
c
      integer  nteexp, ntiexp, nenexp, ninput, nout
c
c  nteexp = number of T_e experimental data points, ...
c
c..BALDUR simulation profiles
c
      dimension rmsim(kr), tesim(kr), tisim(kr), ensim(kr)
      integer nsim
c
c  rmsim(j)   = major radii from BALDUR simulation for all profiles
c  tesim(j)   = T_e(R) (keV)
c  tisim(j)   = T_i(R) (keV)
c  ensim(j)   = n_e(R) (m^{-3})
c  nsim       = number of elements in simulation profile arrays
c
      dimension teint(kr), tiint(kr), enint(kr)
c
c  teint(j) = simulation T_e interpolated to rmteexp(j)
c  tiint(j) = simulation T_i interpolated to rmtiexp(j)
c  enint(j) = simulation n_e interpolated to rmenexp(j)
c
      real tedev(kr), tidev(kr), endev(kr)
     & , tereldev(kr) , tireldev(kr), enreldev(kr)
c
c  tedev(j)    = deviation teexp(j) - teint(j)
c  tereldev(j) = deviation relative to maximum
c              = ( teexp(j) - teint(j) ) / temax
c
      real temax, timax, enmax
     &  , tetemp, titemp, entemp
c
c  temax = maximum value of teexp(j), ...
c
      real repsilon, epslon, epsqrt
c
c  repsilon = smallest major radius to be accepted
c  epslon   = machine epsilon
c  epsqrt   = sqrt ( epslon )
c
      epslon = 0.1
   10 continue
      epslon = epslon * 0.1
      if ( 1.0 + epslon * 0.1 .gt. 1.0 ) go to 10
c
      epsqrt = sqrt ( epslon )
c
c
      repsilon = 0.001
      nout = 7
c
c
      open (nout,file='output')
c
      call date (cdate)
      call time (ctime)
c
c..header on output
c
      call prtext (nout,'Simstat program by Glenn Bateman')
      call prtext (nout,' ')
      write (nout,*) ' epslon = ',epslon
      write (nout,*) ' epsqrt = ',epsqrt
c
c
c..read the file containing experimental data
c
      ninput = 4
c
      call expdata ( rmteexp, teexp, teerr, nteexp
     &  , rmtiexp, tiexp, tierr, ntiexp
     &  , rmenexp, enexp, enerr, nenexp
     &  , ninput )
c
c..read the file containing BALDUR simulation output
c
      ninput = 4
c
      call baldata ( rmsim, tesim, tisim, ensim, nsim, ninput )
c
c
c..interpolate simulation profiles to experimental data points
c
      ilast = 1
      call int1d (2, 0, 0, nsim, 1, rmsim, ensim
     &  , nenexp, rmenexp, ilast, enint )
c
      ilast = 1
      call int1d (2, 0, 0, nsim, 1, rmsim, tesim
     &  , nteexp, rmteexp, ilast, teint )
c
      ilast = 1
      call int1d (2, 0, 0, nsim, 1, rmsim, tisim
     &  , ntiexp, rmtiexp, ilast, tiint )
c
c
c..find maximum values, deviations, and deviations/max
c
      enmax = 0.0
      do j=1,nenexp
        enmax = max ( enmax, enexp(j) )
      enddo
c
      temax = 0.0
      do j=1,nteexp
        temax = max ( temax, teexp(j) )
      enddo
c
      timax = 0.0
      do j=1,ntiexp
        timax = max ( timax, tiexp(j) )
      enddo
c
c..find deviations and deviations/(maximum values)
c
      do j=1,nenexp
        endev(j)    = enexp(j) - enint(j)
        enreldev(j) = ( enexp(j) - enint(j) ) / enmax
        enrelerr(j) = enerr(j) / enmax
      enddo
c
      do j=1,nteexp
        tedev(j)    = teexp(j) - teint(j)
        tereldev(j) = ( teexp(j) - teint(j) ) / temax
        terelerr(j) = teerr(j) / temax
      enddo
c
      do j=1,ntiexp
        tidev(j)    = tiexp(j) - tiint(j)
        tireldev(j) = ( tiexp(j) - tiint(j) ) / timax
        tirelerr(j) = tierr(j) / timax
      enddo
c
c
c..array outputs
c
      write (7,140) nenexp, enmax
  140 format (t2,'nenexp  = ',i4,/t4,'enmax = ',1pe13.4
     &  /t4,'rmenexp',t17,'enexp',t30,'enint'
     &  ,t43,'endev',t56,'enreldev')
c
      do j=1,nenexp
        write (7,141) rmenexp(j), enexp(j), enint(j)
     &    , endev(j), enreldev(j)
      enddo
  141 format (1p5e13.4)
c
      write (7,142) nteexp, temax
  142 format (t2,'nteexp  = ',i4,/t4,'temax = ',1pe13.4
     &  /t4,'rmteexp',t17,'teexp',t30,'teint'
     &  ,t43,'tedev',t56,'tereldev')
c
      do j=1,nteexp
        write (7,141) rmteexp(j), teexp(j), teint(j)
     &    , tedev(j), tereldev(j)
      enddo
c
      write (7,144) ntiexp, timax
  144 format (t2,'ntiexp  = ',i4,/t4,'timax = ',1pe13.4
     &  /t4,'rmtiexp',t17,'tiexp',t30,'tiint'
     &  ,t43,'tidev',t56,'tireldev')
c
      do j=1,ntiexp
        write (7,141) rmtiexp(j), tiexp(j), tiint(j)
     &    , tidev(j), tireldev(j)
      enddo
c
c..compute variance, maximum values, relative RMS deviation
c
      enrelrms = 0.0
      do j=1,nenexp
        enrelrms = enrelrms + ( enreldev(j) )**2
      enddo
      enrelrms = sqrt ( enrelrms / real(nenexp) )
      enabsrms = enrelrms * enmax
      enper = 100.0 * enrelrms
c
      terelrms = 0.0
      do j=1,nteexp
        terelrms = terelrms + ( tereldev(j) )**2
      enddo
      terelrms = sqrt ( terelrms / real(nteexp) )
      teabsrms = terelrms * temax
      teper = 100.0 * terelrms
c
      tirelrms = 0.0
      do j=1,ntiexp
        tirelrms = tirelrms + ( tireldev(j) )**2
      enddo
      tirelrms = sqrt ( tirelrms / real(ntiexp) )
      tiabsrms = tirelrms * timax
      tiper = 100.0 * tirelrms
c
      call prtext (nout,' ')
      call prtext (nout,'Maximum values of experimental data:')
      call prtext (nout,' ')
      call prtrve (nout,enmax,'= enmax')
      call prtrvf (nout,temax,'= temax')
      call prtrvf (nout,timax,'= timax')
c
      call prtext (nout,' ')
      call prtext (nout,'RMS deviations with units:')
      call prtext (nout,' ')
      call prtrve (nout,enabsrms,'= enabsrms (m^-3)')
      call prtrvf (nout,teabsrms,'= teabsrms (keV)')
      call prtrvf (nout,tiabsrms,'= tiabsrms (keV)')
c
      call prtext (nout,' ')
      call prtext (nout,'RMS deviations relative to maximum, percent')
      call prtext (nout,' ')
      call prtrvf (nout,enper,'= enper %')
      call prtrvf (nout,teper,'= teper %')
      call prtrvf (nout,tiper,'= tiper %')
c
c..compute offsets and corresponding variances
c  assuming no experimental error bars
c  divide by maximum and compute percents
c
      enoff = 0.0
      do j=1,nenexp
        enoff = enoff + enint(j) - enexp(j)
      enddo
      enoff = enoff / real ( nenexp )
c
      enovr = 0.0
      do j=1,nenexp
        enovr = enovr + ( enint(j) - enexp(j) - enoff )**2
      enddo
      enovr = sqrt ( enovr / real ( nenexp ) )
c
      enpof = 100.0 * enoff / enmax
      enpov = 100.0 * enovr / enmax
c
      teoff = 0.0
      do j=1,nteexp
        teoff = teoff + teint(j) - teexp(j)
      enddo
      teoff = teoff / real ( nteexp )
c
      teovr = 0.0
      do j=1,nteexp
        teovr = teovr + ( teint(j) - teexp(j) - teoff )**2
      enddo
      teovr = sqrt ( teovr / real ( nteexp ) )
c
      tepof = 100.0 * teoff / temax
      tepov = 100.0 * teovr / temax
c
      tioff = 0.0
      do j=1,ntiexp
        tioff = tioff + tiint(j) - tiexp(j)
      enddo
      tioff = tioff / real ( ntiexp )
c
      tiovr = 0.0
      do j=1,ntiexp
        tiovr = tiovr + ( tiint(j) - tiexp(j) - tioff )**2
      enddo
      tiovr = sqrt ( tiovr / real ( ntiexp ) )
c
      tipof = 100.0 * tioff / timax
      tipov = 100.0 * tiovr / timax
c
      call prtext (nout,' ')
      call prtext (nout,'Offsets:')
      call prtext (nout,' ')
      call prtrve (nout,enoff,'= enoff')
      call prtrvf (nout,teoff,'= teoff')
      call prtrvf (nout,tioff,'= tioff')
c
      call prtext (nout,' ')
      call prtext (nout,'Variances after offsets with units:')
      call prtext (nout,' ')
      call prtrve (nout,enovr,'= enovr')
      call prtrvf (nout,teovr,'= teovr')
      call prtrvf (nout,tiovr,'= tiovr')
c
      call prtext (nout,' ')
      call prtext (nout,'Percent offset relative to maximum:')
      call prtext (nout,' ')
      call prtrvf (nout,enpof,'= enpof %')
      call prtrvf (nout,tepof,'= tepof %')
      call prtrvf (nout,tipof,'= tipof %')
c
      call prtext (nout,' ')
      call prtext (nout,'Percent variances after offsets:')
      call prtext (nout,' ')
      call prtrvf (nout,enpov,'= enpov %')
      call prtrvf (nout,tepov,'= tepov %')
      call prtrvf (nout,tipov,'= tipov %')
c
c
c..maximum likelihood offset and RMS deviations
c
      call maxlike ( nenexp, enreldev, enrelerr
     &  , epslon, enmloff, enmlrms, enchisqr, ierr )
c
      call maxlike ( nteexp, tereldev, terelerr
     &  , epslon, temloff, temlrms, techisqr, ierr )
c
      call maxlike ( ntiexp, tireldev, tirelerr
     &  , epslon, timloff, timlrms, tichisqr, ierr )
c
      call prtext (nout,' ')
      call prtext (nout
     &  ,'Estimates from the Maximum Likelihood Method:')
c
      entemp = 100.0 * enmloff
      tetemp = 100.0 * temloff
      titemp = 100.0 * timloff
c
      call prtext (nout,' ')
      call prtext (nout,'Offset relative to maximum:')
      call prtext (nout,' ')
      call prtrvf (nout,entemp,'= enmloff %')
      call prtrvf (nout,tetemp,'= temloff %')
      call prtrvf (nout,titemp,'= timloff %')
c
      entemp = 100.0 * enmlrms
      tetemp = 100.0 * temlrms
      titemp = 100.0 * timlrms
c
      call prtext (nout,' ')
      call prtext (nout,'Relative RMS deviations after offsets:')
      call prtext (nout,' ')
      call prtrvf (nout,entemp,'= enmlrms %')
      call prtrvf (nout,tetemp,'= temlrms %')
      call prtrvf (nout,titemp,'= timlrms %')
c
      call prtext (nout,' ')
      call prtext (nout,'Normalized chi square statistic:')
      call prtext (nout,' ')
      call prtrvf (nout,enchisqr,'= enchisqr')
      call prtrvf (nout,techisqr,'= techisqr')
      call prtrvf (nout,tichisqr,'= tichisqr')
c
      call prtext (nout,' ')
c
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      stop
      end

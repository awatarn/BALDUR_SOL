c profile.f  by Glenn Bateman, Lehigh University,  1 Dec 1997
c--------1---------2---------3---------4---------5---------6---------7-c
c
c    profile is a program to extract a profile from 2-D ASCII files
c  and print out a list of the profile as a function of major radius
c  or minor radius.
c  In particular, the profile program does the following:
c  1)  read ASCII 2-D U-files or concatenated U-file or RPLOT file
c        and interpolate in time to to get radial
c        profiles and R_major vs r_minor
c  2)  map the profile to be functions of R_major inboard and outboard
c  3)  write the output in ASCII column format
c
c    Namelist input:
c  filein    = name of input file
c  lrplot    = .TRUE. if input file is an ASCII RPLOT file
c              (default .FALSE.)
c  fileout   = name of output file (for profile vs major radius)
c  profname  = name of profile variable in the input file (ie, TE)
c  rmajname  = name of major radius in the input file
c  rminname  = name of minor radius in the input file
c  rnorname  = name of normalized minor radius 
c  time      = time desired for output profile
c  convrad   = conversion factor to be used for radii
c  convprof  = conversion factor to be used for profile
c  lrmajor   = .TRUE. for output as a function of major radius
c  lrminor   = .TRUE. for output as a function of minor radius
c  lrnorm    = .TRUE. for output as a function of
c                        sq-root(normalized toroid flux)
c  lintrp    = 2 for quasi-Hermite cubic interpolation
c            = 0 for linear interpolation
c
      implicit none
c
      integer kr, kd
c
      parameter (kr=200, kd=2*kr+1)
c
c..input profiles
c

      real  profin(kr), xprofin(kr)
     &    , rminin(kr), xrminin(kr)
     &    , rmajin(kr), xrmajin(kr)
     &    , rnormin(kr),  xrnormin(kr)
c
c  profin(jr)    = input profile
c  xprofin(jr)   = grid read in with profin(jr)
c  rminin(jr)    = input minor radius
c  xrminin(jr)   = grid read in with minor radius
c  rmajin(jr)    = input major radius
c  xrmajin(jr)   = grid read in with major radius
c  rnormin(jr)   = input normalized radius
c  xrnormin(jr)  = grid read in with normalized radius
c
      integer nprofin, nrminin, nrmajin, nrnormin
c
c  nprofin   = number of profin(jr) values read in
c  nrminin = number of rminin(jr) values read in
c  nrmajin = number of rmajin(jr) values read in
c  nrnormin  = number of rnormin(jr) values read in
c
c..interpolated radii
c
      real profint(kr), xprofint(kr), rminint(kr), rmajint(kr)
     &  , rnormint(kr)
c
      integer nprofint
c
c..output profiles
c
       real   profout(kd), rminout(kd), rmajout(kd), rnormout(kd)
c
c  profout(jr)  = output profile
c  rminout(jr)  = minor radius corresponding to profout(jr) values
c  rmajout(jr)  = major radius corresponding to profout(jr) values
c  rnormout(jr) = normalized radius corresponding to profout(jr) values
c
       integer nprofout
c
c  nprofout = number of values to be used in profout(jr) array
c
c      real   xc(kr), rmidxc(kr), rminxc(kr), profxc(kr), rnormxc(kr)
c     &     , xb(kr), rmidxb(kr), rminxb(kr), shafxb(kr)
c
c  xc(jr) = zone center grid sq-root(normalized toroid flux) (~r/a)
c  rmidxc(jr) = major radius to geometric centers of zone centers
c  rminxc(jr) = minor radius (halfwidth) of zone centers
c  profxc(jr) = profile at zone centers
c  xb(jr) = zone boundary grid sq-root(normalized toroid flux)
c  rmidxb(jr) = major radius to geometric centers of zone boundaries
c  rminxb(jr) = minor radius (halfwidth) of zone boundaries
c
c  profrmaj(jr) = output profile as a function of major radius
c
      integer lintrp, iextrap, ideriv, ilast
c
c  lintrp  controls method of interpolation (see options above)
c  iextrap controls extrapolation of ends
c  ideriv  controls oder of derivative
c
c..profile as a function of major radius
c
c
      character *32 filein, fileout, profname, rmajname, rminname
      character *32 xbname, xcname, rnorname
c
      integer  j, ip, lprint
      real     time, convrad, convprof
      integer  nunitin, nunitout, nxdim, nerr
c      integer  nx, nxrmaj, nxrmin, nxrnorm, nxprof
c
c  lprint controls diagnostic printout
c
      logical lrmajor, lrminor, lrnorm, lrplot
      logical lexist, lopen
c
      namelist /nin/  filein, fileout, profname, rmajname, rminname
     &  , time, convrad, convprof, lrmajor, lrminor, lrnorm, lrplot
     &  , rnorname, lprint
c
c
c  nunitin  = input unit number used to read U-files
c  nunitout = output unit number
c  nerr   = error indicator
c
c..defaults
c
      nunitin  = 10
      nunitout = 11
      nxdim = 20
      nerr  = 0
c
      rmajname = 'RMAJOR'
      rminname = 'RMINOR'
      profname = 'TE'
      rnorname = 'XZONI'
c
      time     = 0.0
      convrad  = 1.0
      convprof = 1.0
c
      lrmajor = .TRUE.
      lrminor = .FALSE.
      lrplot  = .FALSE.
      lrnorm  = .FALSE.
c
      lprint  = 0
c
c..read namelist on standard input
c
      read (5,nin,err=94,end=98)
c      
      open(unit=nunitin,  file=filein, status='OLD', action='READ', 
     >     err=90)
      open(unit=nunitout, file=fileout, status='unknown', err=90)
c
c..print header in output file
c
      write (nunitout,102) filein
 102  format ('# Profile from filein = ',(a))
      write (nunitout,103) time
 103  format ('# at time = ',0pf12.6)
c
c..Read minor and major radius as a function of zone boundaries
c
      if ( lrmajor .or. lrminor ) then
c
      inquire(file='../data/tftr_73265_2d.dat',
     & exist=lexist, opened=lopen)
        call bpread ( rminname, nunitin, nxdim, time
     &    , xrminin, nrminin, xbname, rminin, nerr )
c
        write (*,*) nrminin,' radial points for rminin ',rminname
c
        if ( lprint .gt. 0 ) then
          write (*,*) 'j ','xrminin ',rminname
          do j=1,nrminin
             write (*,902) j, xrminin(j), rminin(j)
          enddo
 902      format (i4,1p2e15.5)
        endif
c
        if ( nerr .gt. 0 ) then
          write (*,*) 'abort after bpread rminin, nerr = ',nerr
          stop
        endif
c
        call bpread ( rmajname, nunitin, nxdim, time
     &    , xrmajin, nrmajin, xbname, rmajin, nerr )
c
        write (*,*) nrmajin,' radial points for rmajin ',rmajname
c
        if ( lprint .gt. 0 ) then
          write (*,*) 'j ','xrmajin ',rmajname
          do j=1,nrmajin
             write (*,902) j, xrmajin(j), rmajin(j)
          enddo
        endif
c
        if ( nerr .gt. 0 ) then
          write (*,*) 'abort after bpread rmajin, nerr = ',nerr
          stop
        endif
c
      endif
c
c.. read normalized minor radius from RPLOT file
c
      if ( lrplot  .and.  lrnorm ) then
c
        call bpread ( rnorname, nunitin, nxdim, time
     &    , xrnormin, nrnormin, xcname, rnormin, nerr )
c
        write (*,*) nrnormin,' radial points for rnormin ',rnorname
c
      else
c
c..set rnormin = xprofin
c
        do j=1,nprofin
          xrnormin(j) = xprofin(j)
          rnormin(j)  = xprofin(j)
        enddo
c
      endif
c
c..Read profile
c
      call bpread ( profname, nunitin, nxdim, time
     &  , xprofin, nprofin, xcname, profin, nerr )
c
      write (*,*) nprofin,' radial points for profin ',profname
c
      if ( lprint .gt. 0 ) then
        write (*,*) 'j ','xprofin ',profname
        do j=1,nprofin
           write (*,902) j, xprofin(j), profin(j)
        enddo
      endif
c
      if ( nerr .gt. 0 ) then
        write (*,*) 'abort after bpread profin, nerr = ',nerr
        stop
      endif
c
      if ( nprofin .gt. kr  .or.  nprofin .lt. 2 ) then
        write (*,*)
        write (*,*) 'Abort:  problem with nprofin'
        write (*,*) 'nprofin = ',nprofin
        write (*,*) 'kr     = ',kr
        stop
      endif
c
c..profile as a function of normalized minor radius
c
      if ( lrnorm ) then
c
          write (nunitout,200) profname
 200      format ('#   xc',t17,(a))
c
        if ( lrplot ) then
c
          do j=1,nprofin
            write (nunitout,212) rnormin(j)*convrad
     &        , profin(j)*convprof
          enddo
c
        else
c
          do j=1,nprofin
            write (nunitout,212) xprofin(j)*convrad
     &        , profin(j)*convprof
          enddo
c
        endif
c
      endif
c
      if ( (.not. lrmajor) .and. (.not. lrminor) ) stop
c
c..check data consistency
c
c
c..interpolate radii to profile positions
c
      nprofint = nprofin
      do j=1,nprofint
        profint(j)  = profin(j)
        xprofint(j) = xprofin(j)
      enddo
c
      iextrap = 1
      ideriv  = 0
      ilast   = 1
      call int1d ( lintrp, iextrap, ideriv, nrminin, 1
     &  , xrminin, rminin, nprofint, xprofint, ilast, rminint )
c
      iextrap = 0
      ideriv  = 0
      ilast   = 1
      call int1d ( lintrp, iextrap, ideriv, nrmajin, 1
     &  , xrmajin, rmajin, nprofint, xprofint, ilast, rmajint )
c
      iextrap = 1
      ideriv  = 0
      ilast   = 1
      call int1d ( lintrp, iextrap, ideriv, nrnormin, 1
     &  , xrnormin, rnormin, nprofint, xprofint, ilast, rnormint )
c
c..print profiles as a function of minor radius
c
      if ( lrminor ) then
c
        write (nunitout,220) profname
 220    format ('#   rminor',t17,(a))
c
        do j=1,nprofint
          write (nunitout,212) rminint(j)*convrad, profint(j)*convprof
        enddo
c
      endif
c
c..map to major radius, if requested
c
      if ( lrmajor ) then
c
        ip = 0
c
        do j=nprofint,1,-1
c
          if ( rminint(j) .ge. 0.0 ) then
c
            ip = ip + 1
            profout(ip) = profint(j) * convprof
            rminout(ip) = rminint(j) * convrad
            rmajout(ip) = ( rmajint(j) - rminint(j) ) * convrad
            rnormout(ip) = rnormint(j)
c
          endif
c
        enddo
c
        do j=1,nprofint
c
          if ( rminint(j) .gt. 0.0 ) then
c
            ip = ip + 1
            profout(ip) = profint(j) * convprof
            rminout(ip) = rminint(j) * convrad
            rmajout(ip) = ( rmajint(j) + rminint(j) ) * convrad
            rnormout(ip) = rnormint(j)
c
          endif
c
          nprofout = ip
c
        enddo
c
c..print profiles as a function of major radius
c
        if ( lrmajor ) then
c
          write (nunitout,210) profname
 210      format ('#   rmajor',t17,(a))
c
          do j=1,nprofout
            write (nunitout,212) rmajout(j), profout(j)
          enddo
 212      format (1p8e13.4)
c
        endif
c
      endif
c
      stop
c
  90  continue
c
      write (*,*)
      write (*,*) 'Abort:  unable to open filein  ',filein
      write (*,*) 'or      unable to open fileout ',fileout
      stop
c
 94   continue
      write (*,*)
      write (*,*) 'Abort: unable to read namelist input'
      stop
c
 98   continue
      stop
c
c
c
c
c
      end

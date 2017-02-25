c list.f  by Glenn Bateman, Lehigh University, November 1996
c--------1---------2---------3---------4---------5---------6---------7-c
c
c    list is a program to extract a profile from ASCII U-file or RPLOT
c  file and print out a list of the profile as a function of time
c  at a given value of index, rho, minor or major radius.
c  In particular, the list program does the following:
c  1)  read ASCII 2-D U-files or concatenated U-file or RPLOT file
c        and interpolate in radius to get the time trace
c  2)  map the profile to be functions of R_major inboard and outboard
c  3)  write the output in ASCII column format
c
c    Namelist input in namelist /nin/:
c
c  filein    = name of input file
c  fileout   = name of output file (for list as a function of time)
c
c  nindim    = dimension of input file (1 default or 2)
c  lrplot    = .TRUE. if input file is an ASCII RPLOT file
C              (default .FALSE.)
c
c  profname  = name of profile variable in the input file (ie, TE)
c  profncdf= name of profile variable of netcdf output file
c  rmajname  = name of major radius in the input file
c  rminname  = name of minor radius in the input file
c  time      = time desired for output profile
c  convrad   = conversion factor to be used for radii
c  convprof  = conversion factor to be used for profile
c
c  lrmajor   = .TRUE. for output at a given major radius
c  lrminor   = .TRUE. for output at a given minor radius
c
c  nindex    = index of rho for time trace (default 0 to not use)
c  rho       = value of normalized independent variable
c                for time trace (default -1.0 to not use)
c  rminor    = minor radius for time trace (default -1.0 to not use)
c  rmajor    = major radius for time trace (default -1.0 to not use)
c
c  tlower    = lower bound on time for output
c  tupper    = upper bound on time for output
c  ntskip    = output only every ntskip time points
c
      use ezcdf
      implicit none
c
      integer kr, kd, kt
c
      parameter (kr=120, kd=2*kr, kt=20000)
c
c..profiles as a function of minor radius
c
      real   xc(kr), rmidxc(kr), rminxc(kr), profxc(kr)
     &     , xb(kr), rmidxb(kr), rminxb(kr), shafxb(kr)
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
      integer nindex, ntskip, ndimin
c
      real  rho, rminor, rmajor, tlower, tupper
c
c..profile as a function of major radius
c
      real  xprofile(kr), tprofile(kt), profile2d(kr,kt)
     &    , xrminor(kr),  trminor(kt),  rminor2d(kr,kt)
     &    , xrmajor(kr),  trmajor(kt),  rmajor2d(kr,kt)
     &    , profile1d(kt), time1d(kt)
c
      character *125 filein, fileout, profname, profncdf
     &  , rmajname, rminname, fileoutpath
      character *125 xbname, xcname, cradname, ctimname
     &  , crlabel, crunits
c
      real*4   time, convrad, convprof, zr0, zr2
      integer   :: zdm(3)
c
      integer  nunitint, nunitinp, nunitout, nunitcdf
     &  , nxdim, ntdim, nerr, nufile
      integer  nx, nxrmaj, nxrmin, nxprof
c
      integer  j, inx, ip, im, iprint, irmin, irmax, ierr
     &  , irad, itime, it, itmax, index, ikey
c
      logical lrmajor, lrminor, lrplot
c
      namelist /nin/  filein, fileout, profname, profncdf
     &  , rmajname, rminname
     &  , time, convrad, convprof, lrmajor, lrminor
     &  , lrplot, ndimin
     &  , nindex, ntskip, rho, rminor, rmajor, tlower, tupper
c
c
c  nunitint = input unit number used to read index RPLOT file
c  nunitinp = input unit number used to read input data files
c  nunitout = output unit number
c  nunitcdf = output netcdf unit number
c  nerr   = error indicatorprint header in output file
c
c..defaults
c
      nunitint = 9
      nunitinp = 10
      nunitout = 11
      nunitcdf = 12
      nxdim = kr
      ntdim = kt
      nerr  = 0
c
      ndimin = 1
      lrplot = .FALSE.
c
      rmajname = 'RMAJOR'
      rminname = 'RMINOR'
      profname = 'TE'
      profncdf = ' '
c
      time     = 0.0
      convrad  = 1.0
      convprof = 1.0
c
      lrmajor = .TRUE.
      lrminor = .FALSE.
c
      nindex  = 0
      index   = 0
      ntskip  = 1
      rho     = -1.0
      rminor  = -1.0
      rmajor  = -1.0
      tlower  = -1.0e10
      tupper  = 1.0e10
      irad    = 1
c
c..read namelist on standard input
c
      read (5,nin,err=94,end=98)
      if (len_trim(profncdf) == 0) profncdf = profname
c
c..does namelist data make any sense?
c
      if ( nindex .lt. 1  .and.  rho .lt. 0.0
     &  .and.   rminor .lt. 0.0  .and.  rmajor .lt. 0.0 ) then
        write (*,*) 'Abort: no radial value available for time trace'
        write (6,nin)
        stop
      endif
c      
      open(unit=nunitout, file=fileout, status='unknown', err=90)
c
      ierr = nerr
c
c..read 1-D files if ndimin = 1print header in output file
c
      if ( ndimin .eq. 1 ) then
c
        if ( lrplot ) then
c
c..read 1-D RPLOT file
c
          call r1read ( filein, nunitint, nunitinp, ntdim, profname
     &      , tprofile, itime, profile1d
     &      , ikey, crlabel, crunits, ierr )
c
          if ( ierr .gt. 0 ) then
            write (*,*) 'Abort: error from r1read, ierr = ',ierr
            stop
          endif
c
        else
c
c..read 1-D U-file
c
          open(unit=nunitinp,  file=filein, err=90)
c
          call b2read ( profname, nunitinp, ntdim, 1
     &      , tprofile, itime, ctimname
     &      , xprofile, irad, cradname
     &      , profile1d, ierr )
c
          if ( ierr .gt. 0 ) then
            write (*,*) 'Abort: error from b2read, ierr = ',ierr
            stop
          endif
c
        endif
c
c..compute time1d and profile1d arrays
c
        it = 0
        do j=1,itime
          if ( tprofile(j) .ge. tlower .and.
     &         tprofile(j) .le. tupper ) then
            it = it + 1
            time1d(it)    = tprofile(j)
            profile1d(it) = profile1d(j) * convprof
          endif
        enddo
        itmax = it
c
      else
c
        if ( lrplot ) then
c
c..read 2-D RPLOT file
c
          call r2read ( filein, nunitinp, nxdim, ntdim, profname
     &      , tprofile, itime, irad, profile2d, ierr )
c
          if ( ierr .gt. 0 ) then
            write (*,*) 'Abort: error from r2read, ierr = ',ierr
            stop
          endif
c
        else
c
c..read 2-D U-file
c
          open(unit=nunitinp,  file=filein, err=90)
c
          call b2read ( profname, nunitinp, nxdim, ntdim
     &      , xprofile, irad, cradname, tprofile, itime, ctimname
     &      , profile2d, ierr )
c
          if ( ierr .gt. 0 ) then
            write (*,*) 'Abort: error from b2read, ierr = ',ierr
            stop
          endif
c
        endif
c
c..compute time1d and profile1d arrays
c
        if ( nindex .gt. 0 ) then
          index = min ( nindex, irad )
          it = 0
          do j=1,itime
            if ( tprofile(j) .ge. tlower .and.
     &           tprofile(j) .le. tupper ) then
              it = it + 1
              time1d(it)    = tprofile(j)
              profile1d(it) = profile2d(index,j) * convprof
            endif
          enddo
        endif
        itmax = it
c
      endif
!
!..save output in netcdf file
!
      fileoutpath = fileout(1:scan(fileout, '/', .TRUE.))
      fileout     = fileout(scan(fileout, '/', .TRUE.)+1:)
      if (scan(fileout, '.', .TRUE.)>0) then
        fileout = fileout(1:scan(fileout, '.', .TRUE.)-1) //'.nc'
      else
        fileout = fileout(1:len_trim(fileout)) // '.nc'	
      endif
      call cdfopn(nunitcdf, 
     &            fileoutpath(1:len_trim(fileoutpath))//fileout, 'w')
      zdm = (/itmax, 0, 0/)
      call cdfDefVar(nunitcdf, 'TIME',   zdm,    'R8')
      call cdfDefVar(nunitcdf, profncdf,   zdm,    'R8')
      call cdfPutVar(nunitcdf, 'TIME',      time1d,    ierr)
      call cdfPutVar(nunitcdf, profncdf, profile1d, ierr)
      call cdfCls(nunitcdf)
c
c..print header in output file
c
      write (nunitout,102) filein
 102  format ('# Trace from filein = ',(a))
c
      if ( index .gt. 0 ) then
        write (nunitout,103) index
 103    format ('# at nindex = ',i5)
      endif
c
      write (nunitout,210) profname
 210  format ('#   time',t17,(a))
c
c..print trace as a function of time
c
      if ( itmax .gt. 0 ) then
        do j=1,itmax
          write ( nunitout,212) time1d(j), profile1d(j)
        enddo
      else
        write (*,*) 'Abort: no time points.  itmax = ',itmax
        stop
      endif
c
 212    format (1p8e13.4)
c
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
      end


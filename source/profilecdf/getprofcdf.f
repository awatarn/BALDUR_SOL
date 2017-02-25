c profile.f  by Glenn Bateman, Lehigh University,  1 Dec 1997
c profilecdf.f by Alexei Pankin, 1 Aug 2002
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
c  profncdf  = name of netcdf profile variable in the output file (ie, TE)
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
      include 'netcdf.inc'
c
      integer kr, kd, kt
c
      parameter (kr=200, kt=50, kd=2*kr+1)
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
     &  , rnormint(kr), profile2d(kt, kr)
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
      character *32 filein, fileout, profname, profncdf,
     &              rmajname, rminname, fileoutpath
      character *32 xbname, xcname, rnorname
c
      integer  j, ip, lprint
      real     time(kt), convrad, convprof
      integer  nunitin, nunitout, ncid, nxdim, nerr
      integer  :: zdm(3), ierr, ktime, it
c
c  lprint controls diagnostic printout
c
      logical lrmajor, lrminor, lrnorm, lrplot
!
      integer  :: TIME_dim, XI_dim             ! dimension ids
      integer  :: TIME_id,  XI_id, VAR_id      ! variable ids

      namelist /nin/  filein, fileout, profname, profncdf
     &  , rmajname, rminname
     &  , time, convrad, convprof, lrmajor, lrminor, lrnorm, lrplot
     &  , rnorname, lprint
c
c
c  nunitin  = input unit number used to read U-files
c  nunitout = output unit number
c  ncid     = netcdf id number
c  nerr   = error indicator
c
c..defaults
c
      nunitin  = 10
      nunitout = 11
      nxdim = 20
      nerr  = 0
      ierr  = 0
      xprofin = 0.
      profin  = 0.
c
      rmajname = 'RMAJOR'
      rminname = 'RMINOR'
      profname = 'TE'
      profncdf = ' '
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
      if (len_trim(profncdf) == 0) profncdf = profname
      ktime = count(time>0.)
c
      open(unit=nunitin,  file=filein, STATUS='OLD', err=90)
      open(unit=nunitout, file=fileout, status='unknown', err=90)
!
      fileoutpath = fileout(1:scan(fileout, '/', .TRUE.))
      fileout     = fileout(scan(fileout, '/', .TRUE.)+1:)
      if (scan(fileout, '.', .TRUE.)>0) then
        fileout = fileout(1:scan(fileout, '.', .TRUE.)-1) //'.nc'
      else
        fileout = fileout(1:len_trim(fileout)) // '.nc'
      endif
      ierr = nf_create(fileoutpath(1:len_trim(fileoutpath))//fileout, 
     &                 NF_CLOBBER, ncid)
      call check_err(ierr)
c
c..print header in output file
c
      write (nunitout,102) filein
 102  format ('# Profile from filein = ',(a))
      write (nunitout,103) (time(j), j=1, ktime)
 103  format ('# at time = ',500(0pf12.6, ',  ':))
c
c..Read profile
c
      do it = 1, ktime
        call bpread ( profname, nunitin, nxdim, time(it)
     &    , xprofin, nprofin, xcname, profin, nerr )
c 
!        write (*,*) nprofin,' radial points for profin ',profname
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
        do j = 1, nprofin
          profile2d(it, j) = profin(j)
        enddo
      enddo
c
c..Read minor and major radius as a function of zone boundaries
c
      if (lrmajor)  then ! major radius
c
        call bpread ( rminname, nunitin, nxdim, time(1)
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
     &     , xrmajin, rmajin, nprofin, xprofin, ilast, rmajint )
        ip = 0
        do j=nprofin,1,-1
          if ( rminin(j) .ge. 0.0 ) then
            ip = ip + 1
            profout(ip) = profin(j) * convprof
            rminout(ip) = rminint(j) * convrad
            rmajout(ip) = ( rmajint(j) - rminint(j) ) * convrad
          endif
        enddo
c
        do j=1,nprofin
          if ( rminint(j) .gt. 0.0 ) then
            ip = ip + 1
            profout(ip) = profin(j) * convprof
            rminout(ip) = rminint(j) * convrad
            rmajout(ip) = ( rmajint(j) + rminint(j) ) * convrad
          endif
        enddo
        nprofout = ip
        write (nunitout, '(1A)') '# rmajor'
!
      else if (lrminor) then ! rinor radius
        call bpread ( rmajname, nunitin, nxdim, time(1) 
     &        , xrmajin, nrmajin, xbname, rmajin, nerr )
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
        iextrap = 1
        ideriv  = 0
        ilast   = 1
        call int1d ( lintrp, iextrap, ideriv, nrminin, 1  
     &       , xrminin, rminin, nprofin, xprofin, ilast, xprofint )
c
        write (nunitout, '(1A)') '# rminor'
      else if (lrnorm .and. lrplot) then ! normalized minor radius from RPLOT file
        call bpread ( rnorname, nunitin, nxdim, time(1) 
     &    , xrnormin, nrnormin, xcname, rnormin, nerr )
c
        iextrap = 1
        ideriv  = 0
        ilast   = 1
        call int1d ( lintrp, iextrap, ideriv, nrnormin, 1 
     &       , xrnormin, rnormin, nprofin, xprofint, ilast, xprofint )
        write (*,*) nrnormin,' radial points for rnormin ',rnorname
        write (nunitout, '(1A)') '# rnorm'
      else ! 
        xprofint = xprofin
        write (nunitout, '(1A)') '# rnorm'	
      end if
!
      write (nunitout, '(1A, 500(e13.4,1A:))') 
     &	        'rnorm=', (xprofint(j)*convrad, ', ', j=1, nprofin)
      do it = 1, ktime
        write (nunitout, '(2A,i2,1A, 500(e13.4, 1A:))') 
     &	       profname(1:len_trim(profname)),'(',it,',1)=',
     &	       (profile2d(it,j)*convprof, ', ', j=1, nprofin)
      enddo
!
!...define dimensions
!      
      ierr = nf_def_dim(ncid, 'TIME', ktime, TIME_dim)
      call check_err(ierr)
      ierr = nf_def_dim(ncid, 'XI', nprofin, XI_dim)
      call check_err(ierr)
!
!..save output in netcdf file
!
      zdm = (/TIME_dim, 0, 0/)
      ierr = nf_def_var(ncid, 'TIME', NF_DOUBLE, 1, zdm(1:1), TIME_id)
      call check_err(ierr)
      zdm = (/XI_dim, 0, 0/)
      ierr = nf_def_var(ncid, 'XI', NF_DOUBLE, 1, zdm(1:1), XI_id)
      call check_err(ierr)
      zdm = (/TIME_dim, XI_dim, 0/)
      ierr = nf_def_var(ncid, profncdf, NF_DOUBLE, 2, zdm(1:2), VAR_id)
      call check_err(ierr)
!
!...leave define mode
!
      ierr = nf_enddef(ncid)
      call check_err(ierr)
!
!...store TIME
!
      ierr = nf_put_var_double(ncid, TIME_id, TIME)
      call check_err(ierr)
!
!...store XI
!      
      ierr = nf_put_var_double(ncid, XI_id, xprofin(1:nprofin)*convrad)
      call check_err(ierr)
!
!...store VARIABLE
!
      ierr = nf_put_var_double(ncid, VAR_id, 
     &                         profile2d(1:ktime,1:nprofin))
      call check_err(ierr)
       
      ierr = nf_close(ncid)
      call check_err(ierr)
c
c..map to major radius, if requested
c
 212      format (1p8e13.4)
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
!
      subroutine check_err(ierr)
      integer :: ierr
      include 'netcdf.inc'
      if (ierr .ne. NF_NOERR) then
      print *, nf_strerror(ierr)
      stop
      endif
      end


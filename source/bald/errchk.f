c--------1---------2---------3---------4---------5---------6---------7-c
c@errchk  .../baldur/code/bald/default.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 24-feb-90 removed theta check
c       dhfz 29-aug-78 add theta check
c       dhfz 29-aug-78 add density check
c       dhfz 26-july-78 change ntefit to eeteft, ntifit to eetift
c               and nbfit to eebfit
c       amck 21-jun-77 check only first nzones vals. of radius
c       amck 27-may-77 check only nzones vals. of te, ti
c       amck 26-aug-76 make rmajor, rminor arrays
c******************************************************************************
c
c
        subroutine errchk(k)
c
c       1.9     check input data for errors
c
       include 'cparm.m'
      include 'cbaldr.m'
c
        dimension
     i  iout(3)
c
c------------------------------------------------------------------------------
c
        data    iclass /1/,     isub/9/
c
c
        if(.not.nlomt1(isub)) go to 10
        call mesage(' *** 1.9 subroutine errchk bypassed ')
        return
   10   continue
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       nout, nlend (combas)
c
c------------------------------------------------------------------------------
c
c       1)      check basic code control variables
c
c       1.1)    check dtinit, dtnin, dtmax
c
        i1 = nzones + 1
        i2 = nzones + 2
c
c
        ierr = 0
        if (dtinit.gt.0.0)      go to 110
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0, 
     x	'dtinit zero or negative ')
        call error_olymp(3,dtinit,1,1,' dtinit ')
  110   continue
c
c       2)      geometry parameters
c
c       2.1)    major and minor plasma radii
c
        if (rminor(1).gt.epslon) go to 210
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0, 
     x	'plasma minor radius < or = 0 ')
        call error_olymp(3,rminor,5,mxt,' rminor ')
  210   continue
c
        if (rmajor(1)-rminor(1).gt.epslon)      go to 220
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,1,
     1          'major radius le minor radius ')
        call error_olymp(3,rminor,5,mxt,' rminor ')
        call error_olymp(3,rmajor,5,mxt,' rmajor ')
  220   continue
c
c       2.2)    tabular input for radius
c
        if (radius(1).le.0.0)   go to 230
        do 221 j = 1, nzones
        if (radius(j).le.epslon)        go to 222
  221   continue
        go to 224
  222   continue
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'radius less than or equals 0 ')
        call error_olymp(3,radius,5,nzones,' radius ')
  224   continue
        zrad = radius(1) + epslon
        do 225 j = 2, nzones
        if (zrad.ge.radius(j)) go to 226
        zrad = radius(j) + epslon
  225   continue
        go to 228
  226   continue
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'radii not strictly increasing ')
        call error_olymp(3,radius,5,nzones,' radius ')
  228   continue
  230   continue
c
c
c       3)      plasma input parameters
c
c       3.1)    check te input
c
        if (te(1).gt.0.0)       go to 310
c
c       use te0, te1, eeteft to compute te profile
c
        if (te0.gt.epslon)      go to 302
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'central elec. temp. not > 0 ')
        call error_olymp(3,te0,1,1,' te0    ')
  302   continue
c
        if (te1.gt.epslon)      go to 304
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'boundary elec. temp. not > 0 ')
        call error_olymp(3,te1,1,1,' te1   ')
  304   continue
c
        if (eeteft.gt.0)        go to 306
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'elec. temp. fitting exp. not > 0 ')
        call error_olymp(3,eeteft,2,1,' eeteft ')
  306   continue
        go to 320
c
  310   continue
c
c       tabular data for te
c
        do 312 j=1, nzones
        if (te(j).le.epslon)    go to 311
  312   continue
        go to 320
c
  311   continue
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'elec. temp. not > 0 ')
        call error_olymp(3,te,5,nzones+1,'   te   ')
  320   continue
c
c
c
c       3.2)    check ti input
c
c
        if (ti(1).gt.0.0)       go to 330
c
c       use ti0, ti1, eetift to compute ti profile
c
        if (ti0.gt.epslon)      go to 322
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'central ion temp. not > 0 ')
        call error_olymp(3,ti0,1,1,'  ti0   ')
  322   continue
c
        if (ti1.gt.epslon)      go to 324
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'boundary ion temp. not > 0 ')
        call error_olymp(3,ti1,1,1,' ti1   ')
  324   continue
c
        if (eetift.gt.0)        go to 326
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0, 
     x	'ion temp. fitting exp. not > 0 ')
        call error_olymp(3,eetift,2,1,' eetift ')
  326   continue
        go to 340
c
  330   continue
c
c       tabular data for ti
c
        do 332 j=1, nzones
        if (ti(j).le.epslon)    go to 331
  332   continue
        go to 340
c
  331   continue
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0, 'ion temp. not > 0 ')
        call error_olymp(3,ti,5,nzones+1,'   ti   ')
  340   continue
c
c       3.3)    density
c
        do 350 ii = 1,2
        if(nimp(ii).eq.0) go to 350
        zdens = denim0(ii)+denim1(ii)+denimp(ii,1)+fracti(ii)
        if(zdens.gt.epslon) go to 341
        ierr = ierr + 1
        call error_olymp(1,iclass,isub,0,
     1          'impurity density <= 0 ')
        call error_olymp(3,fracti,5,2,' fracti ')
 341    zdens = 0.
 350    continue
c
c       4.)     0 =< theta <= 1
c
cbate        if(theta.le.1.and.theta.ge.0.0) go to 360
cbate        ierr = ierr + 1
cbate        call error_olymp(1,iclass,isub,0,
cbate     1          'theta >= 1 or <=0 ')
cbate        call error_olymp(3,theta,1,1,' theta  ')
cbate 360    continue
c
c
c
c       90)     if error, quit
c
c
        if(ierr.eq.0) return
c
c
c
 9000   continue
        call error_olymp(2,ierr,2,1,'no.error')
        return
c
c" " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " "
c
c       additional checks to be added:
c
c       1)      see if iota at edge >= 1.0
c       2)      check density like te and ti
c       3)      check nrfit
c       4)      check ngas, nimp, wtgas, wtimp
c       5)      check bpoid, curent, eebfit, and tbpoid
c                       to be positive, bpoid and curent not both
c                       specified, tbpoid strictily increasing.
c       6)      check nzones <= 102
c       7)      check nrun > 0
c       8)      check fracth, fracti
c       9)      check dtmin < dtinit < dtmax
c
c" " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " " "
c
        end

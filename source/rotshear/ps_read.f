      subroutine ps_open(dev,idim,shot,info,kfile)
c
      character*4 dev, info
      integer idim, shot, kfile
      character*6 zzshot
      character*4 zzdev, zzinfo,zzdat
      character*2 zzdim
      character*128 zzpath, homedir
      character*154 filen , rfile
      character*26 file
      include 'com_logfile'
c
      zzinfo = info
      zzdev = dev
       
      write(zzshot ,'(I6)') shot
      if (idim .eq. 1) zzdim = '1d'
      if (idim .eq. 2) zzdim = '2d'

      call getenv ('TOR_DATA',zzpath)
c      zzpath='/net/elf00.elf.chalmers.se/ufs/home2/elf/elfps/DATA/' 
      call PSRBLC ( ZZPATH  , Lpath  )

      IF(ZZPATH(Lpath:lpath) .ne. '/') then
         zzpath = zzpath(1:lpath)//'/'
         lpath=lpath+1
      END IF
      if (lpath .eq. 1) then
         call getenv ('HOME', homedir)
         call psrblc (homedir, lhome)
         IF(homedir(Lhome:lhome) .ne. '/') then
            homedir = homedir(1:lhome)//'/'
            lhome=lhome+1
         END IF
         rfile = homedir(:lhome)//'.torrc'
         open  (32,file=rfile(:lhome+6),status='old',err=100)
         read(32,'(A128)') zzpath
         call PSRBLC ( ZZPATH  , Lpath  )
         IF(ZZPATH(Lpath:lpath) .ne. '/') then
            zzpath = zzpath(1:lpath)//'/'
            lpath=lpath+1
         END IF
         close (32)
      end if
      zzdat = '.dat'
      IF(ZZINFO .eq. '') THEN
         FILE = zzdev//'_'//zzshot//'_'//zzdim//zzdat
      ELSE
         FILE = zzdev//'_'//zzshot//'_'//zzdim//'_'//zzinfo//zzdat
      ENDIF
      call PSRBLC ( FILE, LENGTH1  )
      call caselc ( FILE(1:length1))

      Filen=zzpath(1:lpath)//file(1:length1)
      open (kfile,file = filen, status='old')
      write(klog,*)
      if (idim .eq. 1) then 
          write(klog,*) '1-d U-file: ', filen(1:lpath+length1)
      else
          write(klog,*) '2-d U-file: ', filen(1:lpath+length1)
      end if  
      return
 100  write(klog,*) ' Error opening data files:'
      stop
      end

      subroutine ps_wopen(dev,idim,shot,info,kfile)
c
      character*4 dev, info
      integer idim, shot, kfile
      character*6 zzshot
      character*4 zzdev, zzinfo,zzdat
      character*2 zzdim
      character*154 filen
      character*26 file
      include 'com_logfile'
c
      zzinfo = info
      zzdev = dev

      write(zzshot ,'(I6)') shot
      if (idim .eq. 1) zzdim = '1d'
      if (idim .eq. 2) zzdim = '2d'

      zzdat = '.dat'
      IF(ZZINFO .eq. '') THEN
         FILE = zzdev//'_'//zzshot//'_'//zzdim//zzdat
      ELSE
         FILE = zzdev//'_'//zzshot//'_'//zzdim//'_'//zzinfo//zzdat
      ENDIF
      call PSRBLC ( FILE, LENGTH1  )
      call caselc ( FILE(1:length1))

      Filen=file(1:length1)
      open (kfile,file = filen, status='unknown')
      write(klog,*)
      if (idim .eq. 1) then
          write(klog,*) '1-d U-file: ', filen(1:length1)
      else
          write(klog,*) '2-d U-file: ', filen(1:length1)
      end if
      end

      subroutine ps_lopen(dev,shot,info,kfile)
c
      character*4 dev, info
      integer  shot, kfile
      character*6 zzshot
      character*4 zzdev, zzinfo,zzdat
      character*154 filen
      character*26 file
c
      zzinfo = info
      zzdev = dev

      write(zzshot ,'(I6)') shot
      zzdat = '.log'
      IF(ZZINFO .eq. '') THEN
         FILE = zzdev//'_'//zzshot//zzdat
      ELSE
         FILE = zzdev//'_'//zzshot//'_'//zzinfo//zzdat
      ENDIF
      call PSRBLC ( FILE, LENGTH1  )
      call caselc ( FILE(1:length1))

      Filen=file(1:length1)
      open (kfile,file = filen, status='unknown')
      end


      subroutine ps_varnames( cvarnames, kufile, imax, it,kerr )
c
c   This routine reads concatinated 2-D ASCII U-files (Boucher files).
c
c  cvarnames = array  of independent variable namesray        (output)
c  kufile    = input unit number to read file (default 70)    (input)
c  kerr      = error indicator                            (input/output)
c

      implicit none
c
      integer kufile,  kerr
c
      integer  imaxname,  ierr
     &  , iufile, iendfile, ifdim
     &  ,  iltest, iutest, imax, it
     &  , j
c
c
c  ifdim   = 1 if a 1-D U-file is found, = 2 otherwise
c
c  iufile  = unit number for reading input file
c  iendfile = counter number of times end of input file has been reached
c
c  imaxname = maximum number of characters allowed in cvarname
c  iscalar = number of scalar labels
c  ierr    = input value of kerr

      character cline*72, cline1*72, cline2*72
      character*(*) cvarnames(imax)
      include 'com_logfile'
c
c
c..set defaults
c
      kerr     = 0
      ierr     = 0 
      imaxname = 20
      iendfile = 0
c
      it       = 0  
      ifdim    = 2
c
      cline  = ' '
      cline1 = ' '
      cline2 = ' '

c
c..open concatinated U-file
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) goto 90
c
c..read till we get all cvarnames
c
  10  continue
c
      cline1 = cline2
      cline2 = cline
c
      read (iufile,100,err=92,end=93) cline
 100  format (a)
c
c..find dependent variable line
c
      if ( index( cline, '-DEPENDENT VARIABLE LABEL-' )
     &     .gt. imaxname ) then
c
c..check to see if this is a 1-D or 2-D U-file
c
        if ( index( cline1, '-INDEPENDENT VARIABLE LABEL' )
     &       .gt. imaxname ) then
          ifdim = 2
        else
          ifdim = 1
        endif
c
c..check to see if correct dependent variable name has been found
c
        iltest = 0 
        iutest = 0
        do j=1,imaxname
           iltest = j
           if ( cline(j:j) .ne. ' ' ) go to 11
        end do
 11     iutest = iltest
        do j=iltest,imaxname
        if ( cline(j:j) .eq. ' ' ) go to 13
           iutest = j
        end do

 13     it = it + 1

         cvarnames(it) = cline(iltest:iutest)
      endif 

c
c..skip lines as needed
c
c      read (iufile,100,err=92,end=93) cline
c      read (iufile,104,err=94,end=93) kxu
c 104  format (i11)
c
c
c      if ( ifdim .eq. 1 ) then
c        kyu = 0
c      else
c        read (iufile,104,err=94,end=93) kyu
c      endif
c
c      iskip = kxu/6 + kyu/6 + (kxu*kyu)/6
c
c      do j=1,iskip
c        read (iufile,*,err=92,end=93)
c      enddo
c
c      endif
c
c..go back and read more lines
c
      go to 10

c
c
c..error conditions
c
  90  continue
c
      write (klog,*)
     & 'error reading unit number ', kufile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (klog,*)
      write (klog,*) cline
      write (klog,*)
     & 'error reading line in UFILE in sbrtn b2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  93  continue
c
c..reached the end of the input file
c  rewind input file and try again, or bug out
c
      if ( iendfile .lt. 0 ) then
        iendfile = iendfile + 1
        rewind iufile
        go to 10
      endif
c
  94  continue
c
      if ( ierr .lt. 0 ) stop
      return
c
      end



       subroutine ps_scan1d(kufile)
       Implicit None
c
       integer kufile
c

       include 'dim_readvec'
       include 'com_scan1'
       integer i,j, ierr
c
       integer nx
       real xs(kxdim), ys(kydim), parray(kydim,kxdim)
       character*10 cxname, cyname
       include 'com_logfile'
c
       write(klog,*)
       write(klog,*) 'Scanning 1-d U-file:'
       write(klog,*) '===================='
       write(klog,*)
       write(klog,*) '    Traces found:'
       call ps_varnames(cvnam1, kufile, imax1, it1, ierr)
       DO I = 1, it1
           write(klog,*) '        ',cvnam1(i)
           call b2read(cvnam1(i), kufile, kydim, kxdim
     &       ,ys,ny(i), cxname, xs, nx,cyname, parray, ierr) 
           DO J = 1, ny(I)
               Vscan(J,I) = parray(j,1)
               Yy(j,i)     = ys(j)
           END DO
       END DO
       write(klog,*)
       write(klog,*) 'Scan completed: ',it1,' Traces stored'
       write(klog,*) 
       END

       subroutine ps_scan2d(kufile)
c
       Implicit None

       integer  kufile
       include 'dim_readvec'
       include 'com_scan2'
       integer i,j,k, ierr   
c
       real xs(kxdim), ys(kydim), parray(kxdim,kydim)
       character*10 cxname, cyname
       include 'com_logfile'
c
       write(klog,*)
       write(klog,*) 'Scanning 2-d U-file:'
       write(klog,*) '===================='
       write(klog,*)
       write(klog,*) '    Traces found:'


       call ps_varnames(cvnam2, kufile,imax2, it2, ierr)  

       DO I = 1, it2
           write(klog,*) '        ',cvnam2(i)
           call b2read(cvnam2(i), kufile, kxdim, kydim
     &       ,xs,nx(i), cxname, ys, ny(i),cyname, parray, ierr) 
           DO K = 1, nx(i)
               xx2(k,i) = xs(k)
           END DO
           DO J = 1, ny(I)
               YY2(j,i)     = ys(j)
               DO K = 1, nx(i)
                  Vscan2(K,J,I) = parray(k,j)
c                   if (i.eq.3) then
c                        write(klog,*) x2(k,i),y2(j,i), Vscan2(k,j,i)
c                   end if
               END DO
           END DO
       END DO
       write(klog,*)
       write(klog,*) 'Scan completed: ',it2,'  Traces stored'                   
       write(klog,*)

       END  


      Subroutine ps_read1d(trace, ky, y, nyacts, profile, ierr)
c
      character*(*) trace
      integer       ky, nyacts, ierr
      real          y(ky), profile(ky)
c
      include 'dim_readvec'
      include 'com_scan1'
      include 'com_logfile'
c
      character*20 ctest
c
      ierr = 0
c
      call caseuc(trace)
      call psrblc(trace,l)

      DO I = 1, it1
          ctest = cvnam1(i)
          call caseuc(ctest)
          call psrblc(ctest, ltest)
          if (trace(:l) .eq. ctest(:ltest)) then
             j = i
             nyacts = ny(j)
             goto 100
          end if
      END DO

      ierr = 1
      return 

 100  Continue                      

      DO I = 1, nyacts
         y(i) = yy(i,j)
         profile(i) = Vscan (i,j)
      END DO

      Return
      END  


      Subroutine ps_read2d(trace,kx,ky,x,nxacts,y,nyacts,profile,ierr)
c
      character*(*) trace
      integer       kx, ky, nxacts, nyacts, ierr
      real          x(kx), y(ky), profile(kx, ky)
c
      include 'dim_readvec'
      include 'com_scan2'
      include 'com_logfile'
c
      character*20 ctest
c
      ierr = 0
c
      call caseuc(trace)
      call psrblc(trace,l)

      DO I = 1, it2
          ctest = cvnam2(i)
          call caseuc(ctest)
          call psrblc(ctest, ltest)
          if (trace(:l) .eq. ctest(:ltest)) then
             j = i
             nyacts = ny(j)
             nxacts = nx(j)
             goto 100
          end if
      END DO

      ierr = 1
      return 

 100  Continue                      

      DO I = 1, nxacts
        x(i) = xx2(i,j)  
      END DO
c
      DO I = 1, nyacts
         y(i) = yy2(i,j)
         DO K = 1, nxacts
            profile(k,i) = Vscan2 (k,i,j)
         END DO
      END DO

      Return
      END  

      Subroutine ps_get1d(profile, trace,  nyfix, yn, ierr)
c
      Integer       nyfix
      Character*(*) trace
c
      include 'dim_readvec'
      include 'def_interpolate'
      include 'com_filedata'
      include 'com_logfile'
c
      real profile(kydim)
c
      real y(kydim)
      real yn(nyfix)
c
      ierr = 0
      Call ps_read1d(trace, kydim,y,nyacts, profile, ierr)

      If (Ierr.NE.0) THEN
         IF( IERR.NE.1) THEN 
             write(klog,999) 'Error Reading 1-d trace',trace,ierr
         ENDIF
      ELSE
        call tlc1(nyacts,y,profile,nyfix,yn,profile,intpol,
     &            lw,w, lwmin,ierr)
      ENDIF
 999  Format( 2x,a,' ',a,' <-->','IERR= ',I5.5)

 10   Return
      END


      Subroutine ps_get2d(rho, nxp, trace, profile, nyfix, yn, ierr)
c
      Integer       nxp, nyfix
      real          rho(nxp)
      Character*(*) trace
c
      include 'dim_readvec'
      include 'com_filedata'
c
      real          profile(kxdim,*)
      real          fmax(kydim)  
      COMMON /edgedata/ fmax
c
      real data(kxdim, kydim)
      real x(kxdim) , y(kydim)
      real yn(kydim)
      include 'com_logfile'
c
      ierr = 0

      Call ps_read2d(trace, kxdim, kydim,x,nxacts,y,nyacts,data,ierr)
      DO I = 1, nyacts
          fmax(I) = data(nxacts, I)
      END DO 
      If (Ierr.NE.0) THEN
         IF( IERR.NE.1) THEN 
             write(klog,999) 'Error Reading 2-d trace',trace,ierr
         ENDIF
      ELSE
        call ps_regrid(rho, nxp, nxacts,x, nyacts,y, DATA,PROFILE)
        nyfix = nyacts
        DO I = 1, nyacts
           yn (i) = y(i)
        END DO
      ENDIF
 999  Format( 2x,a,' ',a,' <-->','IERR= ',I5.5)

 10   Return
      END
      Subroutine ps_regrid(rho, nx, nxacts, x, nyacts, y, f, p)

      Implicit None

c-------------------------------------------------------------------------
c This program is intended to regrid the ITER-EDA database files onto a
c local grid used in the transport code.
c-------------------------------------------------------------------------

      include 'dim_readvec'

      Integer nx                          ! Lenths of the passed vectors
      Real    rho(nx)                     ! Output grid rho
      Real    p(kxdim,kydim)              ! Output array p

      Integer nxacts, nyacts              ! Lengths of the read vectors
      Real    x(nxacts), y(nyacts)        ! Input grid x,y
      Real    f(kxdim, kydim)             ! Input Array f

c-------------------------------------------------------------------------
c Definitions for local help grids
c-------------------------------------------------------------------------

      Real    xn(0:kxdim), Sn(0:kxdim)
      Real    S(kxdim)
      Real    fmn(kydim), fmx(kydim), fmin, fmax, fdum
      Common /minmax_values/ fmin, fmax 
c-------------------------------------------------------------------------
c Definitions for the interpolation routines
c-------------------------------------------------------------------------

      include 'def_interpolate'
      include 'com_setaxis'
c-------------------------------------------------------------------------
c Loop variables
c-------------------------------------------------------------------------

      Integer I,J
      include 'com_logfile'

c-------------------------------------------------------------------------
c Start of Coding:
c-------------------------------------------------------------------------

      IF (X(1) .eq. 0. ) THEN              !Ordinary entry

         DO J = 1, nyacts
            DO I = 1, nxacts
               S(I) = f(I,J)
            END DO
            CALL tlc1(nxacts,x,S,nx,rho,S,intpol,lw,w,lwmin,ierr)
            call ps_minmax(nx,S, fmn(J),fmx(j))
            DO I = 1, nx
               p(I,J) = S(i)
            END DO
         END DO
         CALL ps_minmax(nyacts, fmn, fmin,fdum)
         CALL ps_minmax(nyacts, fmx, fdum,fmax) 
      ELSE IF(X(1) .ne. 0.) THEN            !Fix for X(1) > 0

         DO J = 1, nyacts
            DO I = 1, nxacts
               SN(I) = f(I,J)
               XN(I) = X(I)
            END DO
            IF (SETZERO) THEN
               SN(0) = 0.
               XN(0) = 0.
            ELSE
               SN(0) = SN(1)
               XN(0) = -X(1)
            END IF
            CALL tlc1(nxacts+1,xn,Sn,nx,rho,S,intpol,lw,w,lwmin,ierr)
            call ps_minmax(nx,S, fmn(J),fmx(j))
            DO I = 1, nx
               P(I,J) = S(I)
            END DO
         END DO
         CALL ps_minmax(nyacts, fmn, fmin,fdum)
         CALL ps_minmax(nyacts, fmx, fdum,fmax) 
      END IF
      setzero = .false.
      RETURN
      END


       Subroutine ps_linpol(t, y, p, kxdim, nxp, nyacts, f)
c
       Integer nxp, nyacts
       Real t, y(nyacts), p(kxdim, nyacts), f(nxp)
c
       Integer ilatest, il 
       COmmon /LASTREAD/ ilatest
       include 'com_logfile'
c
       IF (t .lt. y(1) ) GOTO 200
       IF (t.gt.y(nyacts)) goto 300
 10    IF (ilatest .eq.0) ilatest = 1
       IF (t.ge. y(ilatest) .and. t.lt. y(ilatest+1)) THEN
          goto 100
       ELSE
          DO I = 1, nyacts
             IF (t.ge. y(I) .and. t.lt. y(i+1)) then
                  ilatest = I
	          goto 100
             END IF
          END DO
       END IF

 100   CONTINUE
       il = ilatest
       DO I = 1, nxp
          f(i)=(p(i,il+1)-p(i,il))/(y(il+1)-y(il))*(t-y(il))+p(i,il)
       END DO
       RETURN
 200   ilatest = 1
       DO I = 1,nxp
          f(i) = p(i,ilatest)
       END DO
       RETURN
 300   ilatest = nyacts
       DO I = 1,nxp
          f(i) = p(i,ilatest)
       END DO
       RETURN
       END  
 
      subroutine ps_minmax(nxp,profile, minval, maxval)
c
      integer nxp
      real  profile(nxp)
      real  minval
      real  maxval
c
      integer i
c
      minval = profile(1)
      maxval = profile(1)
c
      Do i = 1,nxp
         minval =  min(minval,profile(i))
         maxval =  max(maxval,profile(i))
      End do
c
      RETURN
      END
      subroutine ps_file(zzpath,zzdev,zzdim,zzshot,zzinfo,zztrace,file,
     &                 lpath,lfile) 
      IMPLICIT NONE
      character*(*) ZZPATH, zzdev, zzshot, zzinfo, zztrace, file, zzdim
c

      integer lpath,lfile,  length1, length2

      FILE = zzdev//zzdim//zzshot//zzinfo//'.'

      call PSRBLC ( ZZPATH  , Lpath  )
      call PSRBLC ( FILE, LENGTH1  )
      call caselc ( FILE(1:length1))
      call PSRBLC ( ZZTRACE , length2  )
      call caseuc ( ZZTRACE (1:length2))
      FILE=FILE(1:length1)// ZZTRACE(1:length2)
      lfile=length1+length2
      return
      end






*=======================================================================
*  ##### CASELC #####
*
      SUBROUTINE CASELC(C)
      CHARACTER*(*) C !![Modified]
C
      INTEGER J,L
C
      L=LEN(C)
      J=0
00001 J=J+1
      IF(J.GT.L)RETURN
      IF((C(J:J).GE.'A').AND.(C(J:J).LE.'Z'))
     &C(J:J)=CHAR(ICHAR(C(J:J))+32)
      GOTO 00001
      END
*=======================================================================
*  ##### CASEUC #####
*
      SUBROUTINE CASEUC(C)
      CHARACTER*(*) C !![Modified]
C
      INTEGER J,L
C
      L=LEN(C)
      J=0
00001 J=J+1
      IF(J.GT.L)RETURN
      IF((C(J:J).GE.'a').AND.(C(J:J).LE.'z'))
     &C(J:J)=CHAR(ICHAR(C(J:J))-32)
      GOTO 00001
      END
C  ###MODULE###  PSRBLC.F  (1995-05-03/P.Strand)
C
C  Purpose: To move the blanks of a CHARACTER variable to the end,
C  and to compute the number of non-blank characters.
C
      SUBROUTINE PSRBLC(C,J)
      CHARACTER*(*) C !![Modified] Text string
      INTEGER       J !![Output] Number of non-blank characters found
C
      INTEGER J1,L1
      CHARACTER*1 C1
C
      L1=LEN(C)
      J=0
      DO 00001 J1=1,L1
      C1=C(J1:J1)
      IF(C1.EQ.' ')GOTO 00001
      J=J+1
      C(J:J)=C1
00001 CONTINUE
      IF(J.LT.L1)C(J+1:L1)=' '
      RETURN
      END

c
c
c ... file tlc1.f
c
c     this file contains documentation for subroutine tlc1 followed by
c     fortran code for tlc1 and additional subroutines.
c
c ... author
c
c     John C. Adams (NCAR 1994)
c
c ... subroutine tlc1(nx,x,p,mx,xx,q,intpol,lw,w,lwmin,ier)
c
c ... purpose
c
c     subroutine tlc1 interpolates the values p(i) on the grid x(i)
c     for i=1,...,nx onto q(ii) on the grid xx(ii),ii=1,...,mx.
c
c ... test program
c
c     file testlc1.f on tlcpack includes a test program for subroutine tlc1
c     and results from a Sun-4 work station and the Cray YMP
c
c ... method
c
c     linear or cubic interpolation is used (see argument intpol)
c
c ... required files
c
c     none
c
c ... requirements
c
c     x must be a strictly increasing grid and xx must be an increasing
c     grid (see ier = 4).  in addition the interval
c
c          [xx(1),xx(mx)]
c
c     must lie within the interval
c
c          [x(1),x(nx)].
c
c     extrapolation is not allowed (see ier=3).  if these intervals
c     are identical and the x and xx grids are UNIFORM then subroutine
c     tlc1u (see file tlc1u.f) should be used in place of tlc1.
c
c ... required files
c
c     none
c
c ... language
c
c     coded in portable fortran77.  innermost loops vectorize on Crays.
c
c
c *** input arguments
c
c
c ... nx
c
c     the integer dimension of the grid vector x and the dimension of p.
c     nx > 1 if intpol = 1 or nx > 3 if intpol = 3 is required.
c
c ... x
c
c     a real nx vector of strictly increasing values which defines the x
c     grid on which p is given.
c
c
c ... p
c
c     a real nx vector of values given on the x grid
c
c ... mx
c
c     the integer dimension of the grid vector xx and the dimension of q.
c     mx > 0 is required.
c
c ... xx
c
c     a real mx vector of increasing values which defines the
c     grid on which q is defined.  xx(1) < x(1) or xx(mx) > x(nx)
c     is not allowed (see ier = 3)
c
c ... intpol
c
c     an integer which sets linear or cubic
c     interpolation as follows:
c
c        intpol = 1 sets linear interpolation
c        intpol = 3 sets cubic interpolation
c
c     values other than 1 or 3 in intpol are not allowed (ier = 6).
c
c ... lw
c
c     the length of the work space w.  let
c
c          lwmin = 2*mx            if intpol(1) = 1
c          lwmin = 5*mx            if intpol(1) = 3
c
c     then lw must be greater than or equal to lwmin
c
c ... w
c
c     a real work space of length at least lw which must be provided in the
c     routine calling tlc1
c
c
c *** output arguments
c
c
c ... q
c
c     a real mx vector of values on the xx grid which are
c     interpolated from p on the x grid
c
c ... lwmin
c
c     the minimum required integer work space length for the current set of
c     input arguments (see lw).  lwmin is returned even if
c     lw < lwmin (ier = 5).
c
c ... ier
c
c     an integer error flag set as follows:
c
c     ier = 0 if no errors in input arguments are detected
c
c     ier = 1 if  mx < 1
c
c     ier = 2 if nx < 2 when intpol=1 or nx < 4 when intpol=3
c
c     ier = 3 if xx(1) < x(1) or x(nx) < xx(mx)
c
c *** to avoid this flag when end points are intended to be the
c     same but may differ slightly due to roundoff error, they
c     should be set exactly in the calling routine (e.g., if both
c     grids have the same x boundaries then xx(1)=x(1) and xx(mx)=x(nx)
c     should be set before calling tlc1)
c
c     ier = 4 if the x grid is not strictly monotonically increasing
c             or if the xx grid is not montonically increasing.  more
c             precisely if:
c
c             x(i+1) <= x(i) for some i such that 1 <= i < nx (or)
c
c             xx(ii+1) < xx(ii) for some ii such that 1 <= ii < mx
c
c     ier = 5 if lw < lwmin (insufficient work space)
c
c     ier = 6 if intpol is not equal to 1 or 3
c
c ************************************************************************
c
c     end of tlc1 documentation, fortran code follows:
c
c ************************************************************************
c
      subroutine tlc1(nx,x,p,mx,xx,q,intpol,lw,w,lwmin,ier)
      dimension x(nx),p(nx),xx(mx),q(mx),w(lw)
c
c     check arguments for errors
c
      ier = 1
c
c     check xx grid resolution
c
      if (mx .lt. 1) return
c
c     check intpol
c
      ier = 6
      if (intpol.ne.1 .and. intpol.ne.3) return
c
c     check x grid resolution
c
      ier = 2
      if (intpol.eq.1 .and. nx.lt.2) return
      if (intpol.eq.3 .and. nx.lt.4) return
c
c     check xx grid contained in x grid
c
      ier = 3
      if (xx(1).lt.x(1) .or. xx(mx).gt.x(nx)) return
c
c     check montonicity of grids
c
      do 1 i=2,nx
      if (x(i-1).ge.x(i)) then
      ier = 4
      return
      end if
    1 continue
      do 2 ii=2,mx
      if (xx(ii-1).gt.xx(ii)) then
      ier = 4
      return
      end if
    2 continue
c
c     set and check minimum work space
c
      if (intpol.eq.1) then
      lwmin = mx+mx
      else
      lwmin = 5*mx
      end if
      if (lw .lt. lwmin) then
      ier = 5
      return
      end if
c
c     arguments o.k.
c
      ier = 0

      if (intpol.eq.1) then
c
c     linear interpolation in x
c
      i1 = 1
      i2 = i1+mx
      call linmx(nx,x,mx,xx,w(i1),w(i2))
      call lint1(nx,p,mx,q,w(i1),w(i2))
      return
      else
c
c     cubic interpolation in x
c
      i1 = 1
      i2 = i1+mx
      i3 = i2+mx
      i4 = i3+mx
      i5 = i4+mx
      call cubnmx(nx,x,mx,xx,w(i1),w(i2),w(i3),w(i4),w(i5))
      call cubt1(nx,p,mx,q,w(i1),w(i2),w(i3),w(i4),w(i5))
      return
      end if
      end

      subroutine cubnmx(nx,x,mx,xx,ix,dxm,dx,dxp,dxpp)
      dimension x(nx),xx(mx),ix(mx),dxm(mx),dx(mx),dxp(mx),dxpp(mx)
      isrt = 1
      do 1 ii=1,mx
c
c     set i in [2,nx-2] closest s.t.
c     x(i-1),x(i),x(i+1),x(i+2) can interpolate xx(ii)
c
      do 2 i=isrt,nx-1
      if (x(i+1) .ge. xx(ii)) then
      ix(ii) = min0(nx-2,max0(2,i))
      isrt = ix(ii)
      go to 3
      end if
    2 continue
    3 continue
    1 continue
c
c     set cubic scale terms
c
      do 4 ii=1,mx
      i = ix(ii)
      dxm(ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/
     +        ((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
      dx(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/
     +        ((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
      dxp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2))/
     +        ((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
      dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/
     +        ((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
    4 continue
      return
      end

      subroutine cubt1(nx,p,mx,q,ix,dxm,dx,dxp,dxpp)
      dimension p(nx),q(mx)
      dimension ix(mx),dxm(mx),dx(mx),dxp(mx),dxpp(mx)
c
c     cubically interpolate p on x to q on xx
c
      do 2 ii=1,mx
      i = ix(ii)
      q(ii) = dxm(ii)*p(i-1)+dx(ii)*p(i)+dxp(ii)*p(i+1)+dxpp(ii)*p(i+2)
    2 continue
      return
      end


      subroutine linmx(nx,x,mx,xx,ix,dx)
c
c     set x grid pointers for xx grid and interpolation scale terms
c
      dimension x(nx),xx(mx)
      dimension ix(mx),dx(mx)
      isrt = 1
      do 1 ii=1,mx
c
c     find x(i) s.t. x(i) < xx(ii) <= x(i+1)
c
      do 2 i=isrt,nx-1
      if (x(i+1) .ge. xx(ii)) then
      isrt = i
      ix(ii) = i
      go to 3
      end if
    2 continue
    3 continue
    1 continue
c
c     set linear scale term
c
      do 4 ii=1,mx
      i = ix(ii)
      dx(ii) = (xx(ii)-x(i))/(x(i+1)-x(i))
    4 continue
      return
      end

      subroutine lint1(nx,p,mx,q,ix,dx)
      dimension p(nx),q(mx),ix(mx),dx(mx)
c
c     linearly interpolate p on x onto q on xx
c
      do 1 ii=1,mx
      i = ix(ii)
      q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
    1 continue
      return
      end





c@b2read.f   02 Dec 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine b2read ( cvarname, kufile, kxdim, kydim
     & , pxu, kxu, cxname, pyu, kyu, cyname, parray, kerr )
c
c   This routine reads concatinated 2-D ASCII U-files (Boucher files).
c
c  cvarname  = name of the independent variable parray        (input)
c  kufile    = input unit number to read file (default 70)    (input)
c  kxdim     = maximum number of elements allowed in x array    (input)
c              (first dimension of parray)
c  kydim     = maximum number of elements allowed in y array    (input)
c              (second dimension of parray)
c
c  pxu(jx)   = 1-D independent variable x array from file     (output)
c  kxu       = number of elements in x array from file        (output)
c  cxname    = name of independent variable x                 (output)
c  pyu(jx)   = 1-D independent variable y array from file     (output)
c  kyu       = number of elements in y array from file        (output)
c  cyname    = name of independent variable x                 (output)
c  parray(jx,jy) = 2-D array from file                        (output)
c
c  kerr      = error indicator                            (input/output)
c
c..On input, set kerr < 0 if you want the program to stop on error.
c
c..If kerr .ge. 0 on input,
c    the routine will return with the following error flag set:
c
c  kerr      = 0 for normal output
c            = 1 if the file cannot be opened or read
c            = 2 reached the end of file before finding cvarname
c            = 3 if there are too many elements in the arrays
c                (that is, if kxu > kxdim)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
      implicit none
c
      integer kxdim, kxu, kydim, kyu, kufile,  kerr
c
      real parray(kxdim,*), pxu(*), pyu(*), zu(6)
c
      integer  imaxname, iprint, ierr, ifdim
     &  , iufile, iendfile, iu, iumax, iskip
     &  , ilow1, iupp1,icpos, iltest, iutest
     &  , j, jx, jy, ju
c
c
c  iprint  = integer to control amount of diagnostic output
c
c  ifdim   = 1 if a 1-D U-file is found, = 2 otherwise
c
c  iufile  = unit number for reading input file
c  iendfile = counter number of times end of input file has been reached
c
c  ixname  = number of characters in cxname
c  iyname  = number of characters in cyname
c  imaxname = maximum number of characters allowed in cvarname
c  iscalar = number of scalar labels
c  ierr    = input value of kerr
c
c  ilow#, iupp# = lower and upper index of names in character strings
c
      character cline*72, cline1*72, cline2*72
      character cvarname*(*), cxname*(*), cyname*(*)
      include 'com_logfile'
c
c
c..set defaults
c
      ierr     = kerr
      kerr     = 0
      imaxname = 20
      iendfile = 0
c
      iprint   = 0
c
      ifdim    = 2
      kxu      = 1
      kyu      = 1
c
c..clear arrays
c
      do jx=1,kxdim
        pxu(jx)    = 0.0
      enddo
c
      do jy=1,kydim
        pyu(jy)    = 0.0
      enddo
c
      do jy=1,kydim
        do jx=1,kxdim
          parray(jx,jy) = 0.0
        enddo
      enddo
c
      cline  = ' '
      cline1 = ' '
      cline2 = ' '
c
c..find first nonblank characters in cvarname
c
      iumax = max ( 1, min ( imaxname, len ( cvarname ) ) )
      do j=1,iumax
        ilow1 = j
        if ( cvarname(j:j) .ne. ' ' ) go to 5
      enddo
c
      write (klog,*)
      write (klog,*) 'Abort: no independent variable name given '
     &    ,cvarname
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
 5    iupp1 = ilow1
      do j=ilow1,iumax
        if ( cvarname(j:j) .eq. ' ' ) go to 6
        iupp1 = j   
      enddo
 6    continue
c
c..open concatinated U-file
c
      iufile = kufile
      if ( kufile .lt. 1  .or.  kufile .gt. 99 ) goto 90
c
c
c..read till we get to the requested cvarname
c
  10  continue
c
      cline1 = cline2
      cline2 = cline
c
      read (iufile,100,err=92,end=93) cline
 100  format (a)
c
c..find dependent variable line
c
      if ( index( cline, '-DEPENDENT VARIABLE LABEL-' )
     &     .gt. imaxname ) then
c
c..check to see if this is a 1-D or 2-D U-file
c
        if ( index( cline1, '-INDEPENDENT VARIABLE LABEL' )
     &       .gt. imaxname ) then
          ifdim = 2
        else
          ifdim = 1
        endif
c
c..check to see if correct dependent variable name has been found
c
        iltest = 0 
        iutest = 0
        icpos = index( cline, cvarname(ilow1:iupp1) )
        if ( icpos .gt. 0 .and. icpos .lt. imaxname ) then 
           do j=1,imaxname
              iltest = j
              if ( cline(j:j) .ne. ' ' ) go to 11
           end do
 11        iutest = iltest
           do j=iltest,imaxname
           if ( cline(j:j) .eq. ' ' ) go to 13
              iutest = j
           end do
 13        if (iupp1-ilow1 .eq. iutest-iltest)  go to 12
        endif 

c
c..skip lines as needed
c
      read (iufile,100,err=92,end=93) cline
      read (iufile,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iufile,104,err=94,end=93) kyu
      endif
c
      iskip = kxu/6 + kyu/6 + (kxu*kyu)/6
c
      do j=1,iskip
        read (iufile,*,err=92,end=93)
      enddo
c
      endif
c
c..go back and read more lines
c
      go to 10
c
c..cvarname has been found, now set cxname and cyname
c
 12   continue
c
      iumax = len ( cline1 )
      cxname = ' '
      iu = 0
      if ( ifdim .ne. 1  .and.   iumax .gt. 0 ) then
        do j=1,iumax
          if ( cline1(j:j) .ne. ' ' ) then
            iu = iu + 1
            cxname(iu:iu) = cline1(j:j)
          else
            if ( iu .gt. 0 ) go to 21
          endif
        enddo
      endif
 21   continue
c
      iumax = len ( cline2 )
      cyname = ' '
      iu = 0
      if ( iumax .gt. 0 ) then
        do j=1,iumax
          if ( cline2(j:j) .ne. ' ' ) then
            iu = iu + 1
            cyname(iu:iu) = cline2(j:j)
          else
            if ( iu .gt. 0 ) go to 22
          endif
        enddo
      endif
 22   continue
c
c..skip line for PROC CODE
c
      read (iufile,100,err=92,end=93) cline
c
c..read number of x and y points
c
      read (iufile,104,err=94,end=93) kxu
c
      if ( ifdim .eq. 1 ) then
        kyu = 0
      else
        read (iufile,104,err=94,end=93) kyu
      endif
c
 104  format (i11)
c
      if ( iprint .gt. 8 ) then
        write (klog,*) kxu,' = kxu'
        write (klog,*) kyu,' = kyu'
      endif
c
      if ( kxu .gt. kxdim  .or.  kyu .gt. kydim ) then
        write (klog,*) 'Too many elements in 2-D file arrays'
        write (klog,*) kxu, kxdim, kyu, kydim
     &   ,' = kxu, kxdim, kyu, kydim'
        write (klog,*)
     &   ' kxu .gt. kxdim .or. kyu .gt. kydim in sbrtn b2read'
        kerr = 1
        if ( ierr .lt. 0 ) stop
        return
      endif
c
c
c..read x and y 1-D arrays
c
      iu = 7
      do jx=1,kxu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kxu + 1 - jx )
          read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pxu(jx) = zu(iu)
      enddo
c
      if ( ifdim .ne. 1 ) then
      iu = 7
      do jy=1,kyu
        iu = iu + 1
        if ( iu .gt. 6 ) then
          iumax = min ( 6, kyu + 1 - jy )
          read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
          iu = 1
        endif
        pyu(jy) = zu(iu)
      enddo
      endif
c
c
c..read 2-D profile
c
      iu = 7
      do jy=1,max(kyu,1)
        do jx=1,kxu
          iu = iu + 1
          if ( iu .gt. 6 ) then
            iumax = min ( 6, max(kyu,1)*kxu + 1 - (jy-1)*kxu - jx )
            read ( iufile, *, err=92,end=93 ) ( zu(ju),ju=1,iumax )
            iu = 1
          endif
        parray(jx,jy) = zu(iu)
c        write(klog,*) parray(jx,jy),zu(iu)
        enddo
      enddo
c
c
c..end of routine
c
      return
c
c
c..error conditions
c
  90  continue
c
      write (klog,*)
     & 'error reading unit number ', kufile
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  92  continue
c
      write (klog,*)
      write (klog,*) cline
      write (klog,*)
     & 'error reading line in UFILE in sbrtn b2read'
      kerr = 1
      if ( ierr .lt. 0 ) stop
      return
c
  93  continue
c
c..reached the end of the input file
c  rewind input file and try again, or bug out
c
      if ( iendfile .lt. 1 ) then
        iendfile = iendfile + 1
        rewind iufile
        go to 10
      endif
c
      write (klog,*)
     & 'reached the end of U-file in sbrtn b2read'
      write (klog,*) ' before finding cvarname = ',cvarname
      kerr = 2
      if ( ierr .lt. 0 ) stop
      return
c
  94  continue
c
      write (klog,*)
      write (klog,*) kxu,' = kxu'
      write (klog,*) kyu,' = kyu'
      write (klog,*)
     & 'error reading kxu or kyu in UFILE in sbrtn b2read'
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
  96  continue
c
      write (klog,*)
      write (klog,*)
     & 'error reading pxu, pyu, or parray in UFILE in sbrtn b2read'
      kerr = 3
      if ( ierr .lt. 0 ) stop
      return
c
      end

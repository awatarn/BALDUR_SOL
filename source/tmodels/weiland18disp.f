!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7--
!| 
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt          
!| \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center} 
!| {\bf {\tt weiland18disp} \\
!| Toroidal Ion Temperature Gradient Mode \\
!| \vspace{1pc}
!| Glenn Bateman \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| Jan Weiland and Hans Nordman \\
!| Chalmers University of Technology \\
!| G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Jon Kinsey \\
!| General Atomics \\
!| P.O. Box 85608, San Diego, CA 92186} \\ 
!| \vspace{1pc}
!| \today
!| \end{center}
!| 
!| 
c
c@weiland18disp.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7--
      subroutine weiland18disp(lprintin,neq,gne,gnh,gnz,gte,gth,gtz, 
     &         tauhin, tauzin,
     &         fnzin, zimp, zmass, zfns, betae, ftrapein, 
     &         vef,QQ,shear, aspinv, kappa, ekyrhoin, hcn,wexb, ROT,tol, 
     &         ITL, ndim, 
     &         omega, gamma, zamrt, zbmrt, zamit, zbmit, zvr, 
     &         zvi,zalfr,zalfi, zbeta, nerr, IK, ITS,ITERA, em)
c***********************************************************************
C
C This routine is a modification of the linear part of ETAWN6 written by
C Glenn Bateman.  The modifications consist of the inclusion of impur-
C ities in the dilution approximatioin in the system with four equations
C and the inclusion of collisions on trapped electrons in the full sys-
C tem.  This system then consists of 7 equations.  When parallel ion 
C motion is included, there are 8 equations with collisions.  Electro-
C magnetic effects are also included.
C
C***********************************************************************
c The variables are ordered as: e phi/Te, Th, Nh, Te, Nz, Tz, F, Vp,Av,
c K where F is due to collisions, Vp is due to parallel ion motion, Av 
c is the vector potential and K is its time derivative.
c---:-7--1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c      subroutine weiland18disp(neq,ZZ)
c
      IMPLICIT none
      integer, PARAMETER :: idp=12
c
      LOGICAL lmatv
      data  lmatv /.true./
c
      REAL gne, gnh, gnz, gte, gth, gtz, zfns, shear, aspinv
      REAL omega(*), gamma(*)
      REAL fnzin, ftrapein, ekyrhoin, vef, KAPPA, HCN
      REAL WIMAX(100)
      REAL zamr(idp,idp), zami(idp,idp), zbmr(idp,idp),zbmi(idp,idp),
     &     zamrt(idp,idp),zbmrt(idp,idp),zamit(idp,idp),zbmit(idp,idp),
     &     zalfr(idp),zalfi(idp),zbeta(idp),zvr(idp,idp),zvi(idp,idp),  
     &     zimp, zmass, em, QQ, betae, WEXB,ROT, TOL,
     &     zepsilon,SH2,ALPHA,
     &     tauhin, tauzin
      INTEGER lprintin, neq, idim, ndim,nerr,ieq,j1,j2,j,IK,IM
      INTEGER ifail,ITL,ITS,ITERA
      COMPLEX WZJ(100), ZZ(idp), WZ(100)

      SAVE WZJ, WZ
      SAVE WIMAX
C
c ITL = maximum allowed number of iterations (6 <= ITL <= 90)
C ITS is a convergence flag, 0-no, 1-yes
C ITC is a flag to use Muller's Method(1), (2) WZI set = 0.01 (remove)
C=======================================================================
C==
      COMPLEX Xbest(90)
      REAL    delWZ, delZ, delZmin
      INTEGER SIndx(10)
C=======================================================================
C Dictionary
C
C INPUT
C
C    INPUT for subroutine eqs only
C       lprintin                          print flag
C       gne, gnh, gnz, gte, gth, gtz      Gradient Temperatures
C       neq                               Number of equations
C       fnzin 
C       zimp
C       zmass
C       zfns
C       betae
C       ftrapein 
C       vef
C       QQ
C       shear
C       aspinv
C       kappa
C       ekyrhoin
C       hcn
C       wexb
C       ROT
C       IK
C                                         should be replaced by a test
C                                         IF neq=9 or 10 and ITERA=1
C       ITS                               Convergence flag 
C       em                                                 
C---------------------------------------------------------------------
C    ndim                                 dimension of something 
C    nerr                                 error flag
C    ITL                                  Iterations limit
C    ITERA                                Iterations limit
C    lmatv                                some sort of flag
C    idp
C    czin                                 has been replaced by zimp
C    ieq
C
C INTERMEDIATE VARIABLES
C    IM                                   index of minimum gamma          
C    j1, j2, j                            counters        
C    ifail                                error flag          
C    ZZ
C    WIMAX                                used in eqs
C    WZJ
C    SH2, ALPHA
C    x()
C    Xbest
C    delWZ, delZ
C    SIndx
C    Loopend
C----------------------------------------------------------------------
C OUTPUT
C
C    omega, gamma, zvr, zvi               computed eigenvalues
C    zamrt, zbmrt, zamit, zbmit           copies of matrices
C    zalfr, zalfi, zbeta                  output from eigvect solver
C    zbeta                                output from eigvect solver
C    tol                                  convergence tolerance   
C    
C=======================================================================
      idim = idp
      ieq  = max(2, neq)
      zepsilon = 1.e-4
C====================
c
      ITL = min ( 90, max ( 6, ITL) )
c
cbate      ITL  = 90
cbate      TOL  = 1.e-11
C====================
c
c ndim  = first dimension of the 2-D array difthi
c           and the maximum number of unstable modes allowed
c ieq   = number of equations 
c
c zamr(i,j) = matrix A
c zbmr(i,j) = matrix B
c   Note that the eigenvalues are
c omega(j) =(zalfr(j)+i zalfi(j))/zbeta(j)
c where beta(j) will be 0.0 in the case of an infinite eigenvalue
c zvr(j) = eigenvector real part
c zvi(j) = eigenvector imag part
c
c ...  STRUCTURE OF MATRIX EQUATION ...
c
c ...  omega*zbmr(i,j) = zamr(i,j)
c
c    variables i=1,6: efi/Te, dTi/Ti, dni/ni, dTe/Te, dnq/nq, dTq/Tq
c    variables j=1,6 same as for i
C=======================================================================
c
      if ( lprintin .gt. 1 ) 
     &   write(*,'(A)')'The program starts here'

cbate      OPEN(34,FILE='Vers6.out')
c
         do j = 1,ITL
            Xbest(j)=0.
         enddo
c
      IF ( lprintin .GT. 2 ) THEN
         write(6,*)
         write(6,*)
     &     'Weiland-Nordman eigenvalue equations, subroutine weiland18'
         write(6,*)'(all frequencies normalized by omega_{De})'
         write(6,*)'(all diffusivities normalized by '
     &     ,'omega_{De}/k_y^2'
      ENDIF
c
c..check validity of input data
c
      if ( lprintin .gt. 1 ) write(*,'(A,f8.3,2(A,i5))')
     &  'zimp=',zimp,' ieq = ',ieq,'neq',neq
c
      IF(neq.LT.2)call abortb(6,' neq.LT.2 in sbrtn weiland18')
      IF(zimp.LT.1.0)call abortb(6,' zimp.LT.1.0 in sbrtn weiland18')
      IF(ieq.LT.2.OR.ieq.EQ.3.OR.ieq.EQ.5.OR.ieq.GT.10)THEN 
         write(6,*)
         write(6,*) ieq,' = ieq in sbrtn disp10'
         call abortb(6,'the value of ieq is wrong in sbrtn disp10')
      ENDIF 
C
C=======================================================================
C PROGRAM LOGIC STARTS HERE
C=======================================================================
C     Needed to compute analytic expression for initial guess
C
      its   = 1
      ITERA = 0
C
      delWZ = 1.e5
      DoWhile (delWZ.GT.TOL)
         ITERA = ITERA + 1
         IF(ITERA.GT.ITL)THEN
             its = 0
             nerr = 77
             write(*,'(A)')
     &         'Excessive Iterations -- weiland18 terminated'
             RETURN
         ENDIF
         IF(ITERA.EQ.1)THEN
c
             CALL weiland18init(shear, kappa, QQ, ekyrhoin, tauhin, 
     &          gne,gnh,gth,
     &          ftrapein, IK, WZJ, WIMAX, Xbest)
c
           if ( lprintin .gt. 1 ) 
     &        write(*,'(A,2F9.3)')'Xbest( 1)=',Xbest(1)
c
             WZ(1) = Xbest(1)
         ENDIF
c
         call weiland18eqns(
     &    lprintin, neq, gne, gnh, gnz, gte, gth, gtz, 
     &    tauhin,tauzin,fnzin, zimp, zmass,zfns, betae, ftrapein, 
     &    vef, QQ, shear, aspinv, kappa, ekyrhoin, hcn, wexb, ROT,
     &    ITL, IK, em, ITERA, 
     &    zamr, zbmr, zami, zbmi,WZJ,WIMAX)
c
C         write(34,'(A)')'A-real'
C         write(34,'(10F9.3)')((zamr(j1,j2),j2=1,10),j1=1,10)
C         write(34,'(A)')'A-imag'
C         write(34,'(10F9.3)')((zami(j1,j2),j2=1,10),j1=1,10)
C         write(34,'(A)')'B-real'
C         write(34,'(10F9.3)')((zbmr(j1,j2),j2=1,10),j1=1,10)
C         write(34,'(A)')'B-Imag'
C         write(34,'(10F9.3)')((zbmr(j1,j2),j2=1,10),j1=1,10)
          
C
C        The output is the matrices zamr, zbmr, Xbest, WZJ
c
c..      find the eigenvalues using NAG14 routine f02bjf or f02gjf
c    
         ifail = 1
         do j1 = 1,ndim
            omega(j1) = 0.0
            gamma(j1) = 0.0
         enddo
c..      save copy of matrix which is over-written by sbrtn f02bjf
c
         do j2=1,ieq
           do j1=1,ieq
             zamrt(j1,j2) = zamr(j1,j2)
             zbmrt(j1,j2) = zbmr(j1,j2)
           enddo
         enddo
c
         call tomsqz(idim, ieq, zamr, zami, zbmr, zbmi, zalfr, zalfi,
     &               zbeta, zvr, zvi, ifail)
c
c..      if ifail = 0, the tomsqz routine has converged
c
         IF(ifail.GT.0)THEN
             write(*,'(A)')'In weiland18disp, tomsqz has failed'
             RETURN
         ENDIF
         do j=1,ieq
            IF(abs(zbeta(j)).LT.1.e-4)THEN
               omega(j) = zalfr(j)/zepsilon
               gamma(j) = zalfi(j)/zepsilon
            ELSE
               omega(j) = zalfr(j)/zbeta(j)
               gamma(j) = zalfi(j)/zbeta(j)
            ENDIF
            ZZ(j)=omega(j)+(0.D0,1.D0)*gamma(j)
         enddo

c
         IF ( lprintin. GT. 6 )
     &     write(*,'(2X,''wr='',G11.4,'' wi='',G11.4)')
     &                       (omega(j),gamma(j),j=1,ieq)
C==
C==      Sort the eigenvalues by gamma (used only for diagnostics)
C==
         CALL ZORT(gamma,ieq,SIndx)
c        write(*,'(20X,A,I2)')'ITERA = ',ITERA
c        write(*,'(A,12F12.7)')'Re(Z)=',(REAL(ZZ(SIndx(j))),j=1,10)
c        write(*,'(A,12F12.7)')'Im(Z)=',(AIMAG(ZZ(SIndx(j))),j=1,10)

         IF(neq.NE.9 .AND. neq.NE.10)THEN
              WZJ(IK) = ZZ(SIndx(1))
              RETURN
         ENDIF
c
      if ( lprintin .gt. 1 ) 
     &    write(*,'(I3,10F10.5)')ITERA,(gamma(SIndx(j)),j=1,ieq)
C=======================================================================
C*********** HERE THE PROCEEDURE TO FIND THE NEW WZ STARTS *************
C==
         IF(ITERA.EQ.1)THEN
            CALL ZORT(gamma,ieq,SIndx)
            WIMAX(IK) = gamma(SIndx(1))
            IF( lprintin .gt. 1 .and. WIMAX(IK) .LT. 0 ) THEN
               write(*,'(I3,2X,A)') ITERA, '+++++all gamma < 0 +++++'
C              RETURN
            ENDIF
            WZ(2)   = (ZZ(SIndx(1))+WZJ(IK))/2.
            WZJ(IK) = WZ(2)
            Xbest(2)= WZ(2)
         ELSE
C==
C==          For subsequent iterations, find the eigenvalue that
C==          is closest to the previous chosen eigenvalue.
C==
            IM = 0
            delZmin = 1000000.0
c           write(*,'(A,2F10.4)')'WZJ(IK) = ',WZJ(IK)
            do j = 1,ieq
               delZ = abs((WZJ(IK)-ZZ(j))/WZJ(IK))
               IF(delZ.LT.delZmin)THEN
                   IM      = j
                   delZmin = delZ
               ENDIF
            enddo
C           write(*,'(2(A,I2))')'Iter = ',ITERA,' IM = ',IM
            WIMAX(IK)    = gamma(IM)    

            CALL ZORT(gamma,ieq,SIndx)
            WIMAX(IK)=gamma(SIndx(1))
            IF( lprintin .gt. 1 .and. WIMAX(IK) .LT. 0 ) THEN
               write(*,'(30X,A,F10.5)')'WIMAX = ',WIMAX(IK)
               write(*,'(30X,I3,2X,A)')ITERA,'++ gamma < 0 ++'
               IF(ITERA.GT.18)RETURN
            ENDIF
            WZ(ITERA+1)   = (ZZ(IM)+WZJ(IK))/2.
C           WZ(ITERA+1)   = (ZZ(SIndx(1))+WZJ(IK))/2.
            WZJ(IK)       = WZ(ITERA+1)
            Xbest(ITERA+1)= WZJ(IK)
            CALL AITKEN(WZ,Xbest,ITERA+1)          
            delWZ=abs(Xbest(ITERA+1)-Xbest(ITERA))/abs(Xbest(ITERA))
         ENDIF  
      ENDDO
      WZJ(IK) = Xbest(ITERA+1)
c
      if ( lprintin .gt. 1 ) 
     &  write(*,'(A,I2,A,2F9.3,/)')'Xbest(',ITERA,')=',Xbest(ITERA)
c
      CALL ZORT(gamma,ieq,SIndx)
c
      if ( lprintin .gt. 1 ) 
     &  write(*,'(A,10F10.5)')'gamma=',(gamma(SIndx(j)),j=1,10)
c
C      IF(gamma(Sindx(1))-Aimag(WZJ(IK)).GT.0.05)THEN
C          redefine initial guess and start over
C      ENDIF

      RETURN     
c
C=======================================================================
      END
C  Removed Code


c..diagnostic output
c
C      IF(lprintin.GT.6)THEN
C         write(*,'(/,A)')' Matrix A-real'
C         write(*,'(10f9.4)')((zamr(j1,j2),j2=1,ieq),j1=1,ieq)
C         write(*,'(/,A)')' Matrix A-imag'
C         write(*,'(10f9.4)')((zami(j1,j2),j2=1,ieq),j1=1,ieq)
C         write(*,'(/,A)')' Matrix B-real'
C         write(*,'(10f9.4)')((zbmr(j1,j2),j2=1,ieq),j1=1,ieq)
C         write(*,'(/,A)')' Matrix B-imag'
C         write(*,'(10f9.4)')((zbmi(j1,j2),j2=1,ieq),j1=1,ieq) 
C      ENDIF
c
C      IF(lmatv.OR.lprintin.EQ.2)THEN
C         do j = 1,ieq
C            IF(gamma(j).GE.0.0001)THEN
C               write(*,'(2X,''wr='',G11.4,'' wi='',G11.4)')
C     &                    omega(j),gamma(j)
C               write(*,'('' Vr='',G11.3,'' Vi='',G11.3,'' j='',I5)')
C     &                  ((zvr(j1,j),zvi(j1,j),j1),j1=1,ieq)
C            ENDIF
C         enddo
C      ENDIF

!| \end{document}

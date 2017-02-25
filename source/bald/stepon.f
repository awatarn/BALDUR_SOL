c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine stepon
c
c       2.1     computes one timestep
c**********************************************************************c
c@stepon  .../baldur/code/bald/dsolver.f
c  les nov-90 add d3he fusion
c  rgb 18.87 12-feb-91 reread = cfutz(82) for backward compatability
c  rgb 18.73 14-oct-90 moved reread to sbrtn REREAD in file DIO
c  dps 15-may-89 15.09 add call to NONCOR: move IRE code into
c                predictor-corrector loop
c  rgb 12.70 11-oct-87 call islbal just before call neugas
c           to interface with saturated tearing mode package
c  rgb 12.58 20-aug-87 Completely revised namelist "reread"
c       added include 'clsaw.m' to sbrtn stepon
c     dps 19-mar-87 add variables to namelist reread to match 1D code
c     dps 27-jan-87 add call to XSCALE after call to GETCHI(2)
c     rgb 27-oct-85 if(lsweq(1) .gt. 9) call sprint for output to terminal
c     rgb 17-jul-85 reworked output to terminal
c     rgb 22-mar-85 added call to mhd to compute equil flux surf aves
c       drm 18-dec-84 added call to beams(4) in sawmix
c**********************************************************************c
c
c               %%%%
c               %
c               tepon          %
c               %
c               %%%%
c
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
      include 'clsaw.m'
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c
        data  iclass /2/,   isub /1/
c
c
        nstep = nstep + 1
        if(.not.nlomt2(isub)) go to 10
        call mesage('% *** 2.1 subroutine stepon bypassed *** ')
        return
   10   continue
c
c-----------------------------------------------------------------------
c
c    1.0)   output to terminal every ntty steps
c
c
        if(ntty.eq.0)go to 100
        istep=nstep-1
        if(ntty*(istep/ntty).ne.istep) go to 100
c
      if (lsweq(1) .gt. 9) then
c
      call sprint (ntychl,2)
c
      else
c
      zlined = 0.0
      do 735 jz=lcentr,ledge
        zlined = zlined + rhoels(2,jz) * dxzoni(jz)
 735  continue
c
      write (ntychl,9000) istep, tai*uset, (tbi-tai)*uset
     &  , tes(1,2)*useh, tis(1,2)*useh, rhoels(1,2), zlined
c
 9000 format (' > nstep =',i5,' time  =',0pf15.10
     & ,' sec, next dt =',0pf13.10,' sec'
     & /' te-axis =',0pf7.3,' ti-axis =',0pf7.3
     & ,' ne-axis =',1pe9.3,' ne-bar =',1pe9.3)
c
      endif
c
 100    continue
c
c..re-read the namelist input file at time REREAD [sec]
c
c  for backward compatability,  CFUTZ(82) may be used instead
c
      if ( cfutz(82) .gt. epslon ) reread = cfutz(82)
c
      if ( (tai*uist.ge.reread) .and. ((tai-dtoldi)*uist.lt.reread)
     &     .and. nstep .gt. 1 ) call redata
c
c
        call sawmix
c
        iiter = 0
c
        call resolv(1)
c
c bateman update flux surface averages or compute equilibrium
c
      call mhd
c
cend bateman 7-apr-85
c
c
c
c       "nadump(11)=2" unlocks the corrector mode in imprad(2)
        if ( cfutz(200) .gt. epslon ) nadump(11)=2
c
  200   continue
c
        call coef
c
  400   continue
c
        call solveb
c
        if (natomc.eq.3) call noncor
c
        call bounds
c
        call reduce
c
        call solve
c
        call resolv(2)
c
c
        iiter = iiter + 1
        if((nliter.or.nlrpet).and.iiter.gt.nitmax) go to 9080
        if(nliter) go to 200
        if(nlrpet)      go to 400
c
        call cmpres
        call pdrive
c
        call getchi(1)
c
c       "nadump(11)=1" allows in imprad(2) corrector values to be
c       saved and restores predictor mode
        if ( cfutz(200) .gt. epslon ) nadump(11)=1
c
        call getchi(2)
c
        call xscale
c
c   les nov-1990 d3he fusion -- update fast particle arrays
c
      if ( cfutz(490) .gt. epslon ) call fusion(0,lcentr,ledge)
c
c..test ignition
c
      if ( cfutz(490) .gt. epslon ) call igntst
c
        dtoldi = tbi - tai
        tai = tbi
        tbi = tbi + dti
c
        return
c
 9080   continue
c
        call error_olymp(1,iclass,isub,7,
     1          '? *** error *** infinite iteration loop')
        call error_olymp(2,iiter,2,1,'iteratns')
        return
        end

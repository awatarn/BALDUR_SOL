c  17:00 24-Nov-99 .../baldur/code/bald/doutput.f  Bateman, Lehigh
c--------1---------2---------3---------4---------5---------6---------7-c
c  file DOUTPUT contains modifications to BALDUR sbrtn output
c  written by Ute Schneider to produce UFILES for use with RPLOT
c  on the PPPL Vax cluser.
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c  Subprogram names in file DOUTPUT:
c
c  output - control output
c  pltrgn - parameter statements and common blocks with filenames
c           moved to crpltrgn.m
c  localp - common blocks used by this package
c           moved to crplocal.m
c  plabel - filenames and other labels
c           moved to crplabel.m
c  cnfile - moved to crpnfile.m
c  mfile  - rescale variables and then call mfwrt1
c  nfile  - computations of output variables from BALDUR
c  pltavg - finds max element and weighted averaga of input array
c  redim  - repacks a 2-D array as a 1-D array
c  rescal - rescales array
c  mfwrt1 - genergic write routine
c  mf0001 - generate output by setting up zbuf and calling mfwrt1
c  tfilsn - set labels for output of each variable
c  tfilwr - write transp plotting ft file, by McCune
c             (containing labeling and dimensioning information)
c  tfilw2 - write transp plotting indexing and labeling file
c             (this file contains all the information needed to locate
c             and plot labelled data in the corresponding "mf" and "nf"
c             files.
c  tfilhd - write header (first record)
c  
c--------1---------2---------3---------4---------5---------6---------7-c
c@output  .../baldur/code/bald/doutput.f
c  rgb 25-mar-02 commented out call hetprt
c  rgb 24-nov-99 commented out calls to tgraf and grafix
c    to remove write outs to files for03 and for04
c  rgb 21-dec-93 changed call mesage(43h...) to call mesage('...')
c  les 03-jan-91 turned off vectorization of do-loop 1211
c      for compilation under cft77, following rgb
c     changed limit of do 1211 from 50 to mxt, since tplot(20)
c     changed test in do 1211 to if(tplot(j).le.epslon) from 0.0
c  les   nov-90  added nfusn=4  ( for d-3he fusion)
c       add call to d-3he print routine dprint (DD3HFUS)
c  rgb 20.26 13:00 17-jan-92 removed call ledger
c  rgb 18.81 20-dec-90  changed do 1211 from 50 to 20
c      removed vectorization from do loop 1211
c  rgb 02:00 07-apr-90 changed snplti to snplt, tnplti to tnplt
c  rgb 23:00 06-feb-90 reconstructed from generator 
c      rx15:[bateman.bal]gbplot.com with data in bdrplt.for
c       dps 24-aug-88 15.00 add six calls to ncprnt for NC code
c       DRM 18-Feb-86 Problems with 20 entries in tplot. Fixed by
c       changing the DO loop upper limit from 50 to mxt for DO 224
c       and DO 1211, as well as adding parentheses to some IF tests
c       to make them less ambiguous. Don't really understand why
c       the runs crashed in the first place.
c       drm 10-jul-85 add calls to plabel, mfile, nfile
c       drm 4-dec-84 change all tbi back to tai
c       aes 14-jan-82 allow nfusn=3
c       aes 19-nov-81 add calls to hetprt,ncfprt,trcprt
c       aes 15-sep-81 add call to rprint on initial printout
c       aes 24-jul-81 add nllprt to prevent spurious final printouts
c       dhfz 5-apr-79 add calls to fprint
c       amck 23-jan-78 add lpage=0
c       amck 30-oct-77 add aprint calls
c       amck 13-sep-77 add iprint calls
c       amck 1-jun-77 fix setting of tnplti near #1104
c       amck 13-apr-77 add tgraf call
c       amck 8-mar-77 change nediti -> nedit
c       amck 17-feb-77 if nledge .le. 0, no call ledger
c       amck 15-feb-77 add call ledger, move zsedit, etc. to common
c       amck 17-jan-77 use ktype=0 as a flag for nprint printout to g-,hprint
c       amck 10-jan-77 eliminate extra nstep=1 grafix call
c       amck 7-jan-77 remove call grafix(2) when k=1
c       amck 7-jan-77 nploti -> nplot
c       amck 6-jan-77 no graphics if plot flags all <= 0
c       amck 6-jan-77 add calls to "grafix"
c       amck 6-aug-76 test tai (not tbi)
c       amck 23-jul-76 add k=2 hprint and gprint calls (fix bug)
c       amck 22-jul-76 add gprint, hprint calls, change initial mprint call
c       amck 21-may-76 change printout spacing
c******************************************************************************
c
c
        subroutine output(k)
c
c
cl      3.1     control the output
c
c
c
c
        include 'cparm.m'
      include 'cbaldr.m'
c
c------------------------------------------------------------------------------
c
        data    iclass/3/,      isub/1/
c
        if ( nlomt3(isub) )   then
          call mesage(' *** 3.1 subroutine output bypassed')
          return
        endif
c
c------------------------------------------------------------------------------
c
cl      common blocks and variables modified:
c
c       comout
c
c------------------------------------------------------------------------------
c
        lpage = 0
        if (k.eq.2)     go to 200
        if (k.ne.1)     go to 300
c
c
c       1.      initial and final mprint calls
c
c
  100   continue
        tnprt = tai*uiet
        if ((sedit.le.epslon) .and. (tedit(1).le.epslon) )
     1                  tnprt = epsinv
        nnprnt = nstep
        if (nedit.le.0) nnprnt = ninfin
        snprt = tai*uiet
        if(sedit.le.epslon) snprt = epsinv
c
c
        call mprint(1)
        call gprint(0,1)
        call hprint(0,1)
        call iprint(0,1)
        call ncprnt(0,1)
c
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(0,1)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(0,1)
        if(nfusn.eq.4) call dprint(0,1)
c
        call rprint(0,1)
        call mprint(2)
        call gprint(0,2)
        call hprint(0,2)
        call iprint(0,2)
c
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(0,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(0,2)
        if(nfusn.eq.4) call dprint(0,2)
c
        call sprint(nprint,2)
        call gprint(nprint,2)
        call hprint(nprint,2)
        call iprint(nprint,2)
c
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(nprint,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(nprint,2)
        if(nfusn.eq.4) call dprint(nprint,2)
c
        if (ndiary.eq.nprint) go to 104
        call sprint(ndiary,2)
        call gprint(ndiary,2)
        call hprint(ndiary,2)
        call iprint(ndiary,2)
c
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(ndiary,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(ndiary,2)
        if(nfusn.eq.4) call dprint(ndiary,2)
c
  104   continue
        go to 1000
c
c
c       2.      intermediate mprint calls
c
c
  200   continue
        nllprt = .false.
        if ( nstep .ge. nnprnt )    go to 210
        if ( tai*uiet .ge. tnprt )  go to 210
        if (nsedit.eq.0) go to 1000
        if ((nstep/nsedit) * nsedit.ne.nstep) go to 1000
        call sprint(nprint,2)
        call ncprnt(nprint,2)
        call gprint(nprint,2)
        call hprint(nprint,2)
        call iprint(nprint,2)
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(nprint,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(nprint,2)
        if(nfusn.eq.4) call dprint(nprint,2)
        go to 1000
c
  210   continue
        nllprt = .true.
        call mprint(2)
        call ncprnt(0,2)
        call ncfprt
        call trcprt
        call gprint(0,2)
cbate        call hetprt
        call hprint(0,2)
        call iprint(0,2)
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(0,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(0,2)
        if(nfusn.eq.4) call dprint(0,2)
c
        call sprint(nprint,2)
        call ncprnt(nprint,2)
        call gprint(nprint,2)
        call hprint(nprint,2)
        call iprint(nprint,2)
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(nprint,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(nprint,2)
        if(nfusn.eq.4) call dprint(nprint,2)
c
        if (ndiary.eq.nprint) go to 214
        call sprint(ndiary,2)
        call gprint(ndiary,2)
        call hprint(ndiary,2)
        call iprint(ndiary,2)
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(ndiary,2)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(ndiary,2)
        if(nfusn.eq.4) call dprint(ndiary,2)
  214   continue
c
c       reset snprt, tnprt, and nnprnt
c
c       mprint edits are generated at every timestep number
c       which is a multiple of nedit(i), and the nearest timestep
c       after a time which is a multiple of sedit(i),
c       and the first timestep after each tedit(i)
c       the array tediti is terminated by a zero-valued
c       element.
c
        if (nstep.ge.nnprnt) nnprnt = nnprnt + nedit
        if (tai*uiet .ge. snprt) snprt = snprt + sedit*ueit
        tnprt = snprt
c
        do 224 j = 1, mxt
          if (tedit(j).le.0.0)   go to 1000
          if ((tedit(j).gt.tai*uiet) .and. (tedit(j).lt.tnprt) )
     &      tnprt = tedit(j)
  224   continue
        go to 1000
c
c
cl      k not 1 or 2
c
c
  300   continue
        if (nllprt) go to 1000
        call mprint(k)
        call ncprnt(0,k)
        if (nstep.gt.0) call ncfprt
        if (nstep.gt.0) call trcprt
        call gprint(0,k)
cbate        if (nstep.gt.0) call hetprt
        call hprint(0,k)
        call iprint(0,k)
        if(nfusn.ne.2.and.nfusn.ne.4) call aprint(0,k)
        if(nfusn.ne.1.and.nfusn.ne.4) call fprint(0,k)
      if(nfusn.eq.4) call dprint(0,k)
        call ncprnt(nprint,k)
        go to 1000
c
c
c
cl      10)     graphics calls
c
c
c
 1000   continue
c
cbate        call tgraf(k)
        if ( (k.eq.2) .and. (nstep.gt.0) ) call nfile(k)
c
        if ( (nplot.le.0) .and. (splot.le.0.0) 
     &      .and. (tplot(1).le.0.0) )  go to 2000
        if (k.eq.2)     go to 1200
        if (k.ne.1) go to 1300
c
c
c       11.     initial grafix calls
c
c
 1100   continue
        nnplot = nstep + nplot
        if (nplot.le.0) nnplot = ninfin
        snplt = tai + splot
        if(splot.le.epslon) snplt = epsinv
        tnplt = snplt
c
        do 1104 jt = 1, mxt
        if (tplot(jt) .le. epslon) go to 1105
        if ( (tplot(jt).lt.tnplt) .and. (tplot(jt).gt.tai*uiet) )
     1               tnplt = tplot(jt)
 1104   continue
c
 1105   continue
c
cbate        call grafix(1)
        call plabel(label1)
        call mfile(1)
        call nfile(1)
        go to 2000
c
c
c       12.     intermediate grafix calls
c
c
 1200   continue
        if (nstep.ge.nnplot)    go to 1210
        if ( tai*uiet .ge. tnplt )      go to 1210
        go to 2000
 1210   continue
cbate        call grafix(2)
        if ( nstep .gt. 0 ) call mfile(2)
c
c       reset snplt, tnplt, and nnplot
c
c       grafix plots are generated at every timestep number
c       which is a multiple of nplot(i), and the nearest timestep
c       after a time which is a multiple of splot(i),
c       and the first timestep after each tplot(i)
c       the array tploti is terminated by a zero-valued
c       element.
c
        if (nstep.ge.nnplot) nnplot = nnplot + nplot
        if (tai.ge.snplt) snplt = snplt + splot
        tnplt = snplt
c
c dir$ novector
        do 1211 j = 1, 20
          if (tplot(j).le.0.0)   go to 2000
          if ( (tplot(j).gt.tai*uiet) .and. (tplot(j).lt.tnplt) )
     &      tnplt = tplot(j)
 1211   continue
c dir$ vector
        go to 2000
c
c
cl      13)     k not 1 or 2
c
c
 1300   continue
cbate        call grafix(k)
        go to 2000
c
c
 2000   continue
        if (nbdump.le.0.or.nledge.le.0) return
cbate        if (nstep.eq.(nstep/nbdump) * nbdump) call ledger(1)
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@mfile
c  rgb 14:40 11-nov-93 replaced call pdx(scroff,6,52) with call pdx
c rgb 20.20 30-nov-91 create... replaces with open (...)
c rgb 20.16 14-oct-91 add call to sbrtn mf0001
c
c
        subroutine mfile(icall)
c
        include 'cparm.m'
      include 'cbaldr.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'crplocal.m'
        include 'crpltrgn.m'
c
        common /cnvect/ vftot(55,9), flxtot(55,9), srctot(55,9)
c
        dimension ztemp(55)
c
        COMMON/SCRATCH/ ZBUF(999)
c
c       template for PLOTRGN to generate code for writing the
c       COMMON variables onto a direct access file
c
c       June 27, 1985  U. Schneider
c
c  zero value to fill out each record
        data zdummy/0.0/
c
c---------------------------------------------------------------------
c  all plotting output variables must be in common
c---------------------------------------------------------------------
c
c   icall=1: initialization call
c
        if(icall.ne.1) go to 10
c
c  CREATE AND OPEN TEXT FILE
c
      open (lunpl1, file=filpl1, status='unknown')
c
cbate      call create (lunpl1, filpl1, 2, -1, 9990)
c
c
 10     continue
        io=lunpl1
c
        imaxx=t1imx   ! nwr array dimension
        imaxsw=t1isw  ! nset array dimension
        if(icall.eq.1) go to 7999
c
c     *     *     *     *     *     *     *     *     *     *     *     *
c       insert local calculations for plotting output here:
c
c                prepare for rescaling the heating profiles to h(r) form
c
        call resetr(zpalfs,mxzone,0.)
        call resetr(zpbeam,mxzone,0.)
        call resetr(zpauxs,mxzone,0.)
        call resetr(ztheat,mxzone,0.)
c
        call resetr(zimrad,mxzone,0.)
c
        do 123 jz=lcentr,isepm1
        zpalfs(jz)=wealfs(jz)+wialfs(jz)
        zpbeam(jz)=webems(jz)+wibems(jz)
        zpauxs(jz)=weauxs(jz)+wiauxs(jz)
        ztheat(jz)=zpalfs(jz)+zpauxs(jz)+zpbeam(jz)+weohms(jz)+
     1  weecrh(jz)+wiecrh(jz) + weicrf(jz)+wiicrf(jz)
c
        zimrad(jz)=webrs(jz)+wesrs(jz)+weirs(1,jz)+
     1  weirs(2,jz)+weirs(3,jz)+weirs(4,jz)
123     continue
c
        izones=isep-lcentr
        ifail=0
        zdummy=0.
c
c  get rescaling factors for heating profiles
c
        call pltavg(weohms(lcentr),dx2i(lcentr),zdummy,zwapoh,
     1  izones,ifail)
        call pltavg(zpalfs(lcentr),dx2i(lcentr),zdummy,zwapal,
     1  izones,ifail)
        call pltavg(zpbeam(lcentr),dx2i(lcentr),zdummy,zwapbm,
     1  izones,ifail)
        call pltavg(zpauxs(lcentr),dx2i(lcentr),zdummy,zwapau,
     1  izones,ifail)
        call pltavg(ztheat(lcentr),dx2i(lcentr),zdummy,zwapht,
     1  izones,ifail)
c
c  get rescaling factors for temp. and dens. profiles
c
        call redim(tes,ztemp,2,2,lcentr,izones)
        call pltavg(ztemp(lcentr),dx2i(lcentr),zmxtes,zdummy,
     1  izones,ifail)
        call redim(tis,ztemp,2,2,lcentr,izones)
        call pltavg(ztemp(lcentr),dx2i(lcentr),zmxtis,zdummy,
     1  izones,ifail)
        call redim(rhoels,ztemp,2,2,lcentr,izones)
        call pltavg(ztemp(lcentr),dx2i(lcentr),zmxnes,zdummy,
     1  izones,ifail)
c
c  get rescaling factors for sawtooth profiles of te
c
        call pltavg(tepre(lcentr),dx2i(lcentr),zmxtei,zdummy,
     1  izones,ifail)
        call pltavg(tepost(lcentr),dx2i(lcentr),zmxtef,zdummy,
     1  izones,ifail)
        call pltavg(teavg(lcentr),dx2i(lcentr),zmxtea,zdummy,
     1  izones,ifail)
c
c  get rescaling factors for sawtooth profiles of ti
c
        call pltavg(tipre(lcentr),dx2i(lcentr),zmxtii,zdummy,
     1  izones,ifail)
        call pltavg(tipost(lcentr),dx2i(lcentr),zmxtif,zdummy,
     1  izones,ifail)
        call pltavg(tiavg(lcentr),dx2i(lcentr),zmxtia,zdummy,
     1  izones,ifail)
c
c  get rescaling factors for sawtooth profiles of jz
c
        call pltavg(ajzpre(lcentr),dx2i(lcentr),zmxjzi,zdummy,
     1  izones,ifail)
        call pltavg(ajzpst(lcentr),dx2i(lcentr),zmxjzf,zdummy,
     1  izones,ifail)
        call pltavg(ajzavg(lcentr),dx2i(lcentr),zmxjza,zdummy,
     1  izones,ifail)
c
c***************************
c
7999    continue
c
c***************************
c
c
c***************************
c
C**> INSERTING GENERATED CODE; TEMPLATE FILE WAS:
C**>  BALPLOT:MFILET.FOR                                             
C >>> PLTRGN GENERATED FORTRAN STARTS HERE:
      IF(ICALL.NE.1) GO TO 8000
C INCREMENT COUNTER, TEST AGAINST IMAXX PARAMETER
      INUMX=   1
      IF(INUMX.GT.IMAXX) GO TO 9000
C DEFINE NUMBER OF PTS FOR THIS AXIS
      ISIZE=NZONES                        
C DEFINE START AND FINISH ADDRESSES FOR FCN WRITES
      ISTART=LCENTR                        
      NWR(1,1,INUMX)=ISTART
      NWR(2,1,INUMX)=ISTART+ISIZE-1
      NWR(2,2,INUMX)=ISIZE
C INCREMENT COUNTER, TEST AGAINST IMAXX PARAMETER
      INUMX=   2
      IF(INUMX.GT.IMAXX) GO TO 9000
C DEFINE NUMBER OF PTS FOR THIS AXIS
      ISIZE=NZONES                        
C DEFINE START AND FINISH ADDRESSES FOR FCN WRITES
      ISTART=LCENTR+1                      
      NWR(1,1,INUMX)=ISTART
      NWR(2,1,INUMX)=ISTART+ISIZE-1
      NWR(2,2,INUMX)=ISIZE
C-----------------------------
C   NOW DEFINE THE FUNCTIONS TO BE PLOTTED (GEN. BY PLTRGN)
 8000 CONTINUE
      IF(ICALL.EQ.1) GO TO 10000
C  WRITE THE CURRENT TIME
      ZBUF(1)=
     >   UIST*TAI                                                     
      CALL MFWRT1(LUNPL1,'TIME    ',ZBUF,1)
      CALL MFIX01
c
c
c***************************
c
10000   continue
c
        return
c----
c  errors
c  branched to from generated fortran
c
 9000   continue
c
      call mesage(
     >' fatal error in MFILEX-- too many independant X ')
      call mesage(
     >' coordinates-- increase array dimension  T1IMX  ')
      call mesage(
     >' in BALDUR common and rerun  @PLOTGEN           ')
c
        stop
c
 9990  call abortb (6,'error creating file in sbrtn mfile')
       stop
c
                end
c--------1---------2---------3---------4---------5---------6---------7-c
c@mfwrt1
c rgb 18-oct-89 if label(1:6) .ne. 'scalar' do not write out label header
c   Doug McCune calls a routine to check and limit the range of output
c
c--------------------------------------------
c  mfwrt1 - generic write routine
c
c  write new records
c
c   fortran writes generated from plotting specification data
c   by program "prplot.for" which is part of automated
c   pre-processing
c
        subroutine mfwrt1(lun,label,buf,n)
c
        character*8 label
c
        real buf(*)
c
        if ( label(1:6) .ne. 'scalar' ) write(lun,101) n,label
        write(lun,102) buf(1:n)
c
        return
c  format statements
 101    format(1x,i3,1x,a8)
 102    format(1x,6(1pe12.5))
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
      SUBROUTINE MFIX01
C
C   THIS IS CALLED BY MFILE -- BREAKS UP
C   THE MF#### CALLS TO AVOID A POTENTIAL
C   COMPILER PROBLEMS (TOO MANY EXTERNALS)
C========================================
C      THIS IS PLTRGN GENERATED CODE
C========================================
C
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'cfokkr.m'
       include 'cfreya.m'
       include 'crpltrgn.m'
       include 'crplocal.m'
       include 'crpnfile.m'
C------------------/PLTRGN:  GENERATED FORTRAN ROUTINE
      IF(NSET(   2).EQ.1) CALL MF0001
C
       RETURN
       END
C------------------/PLTRGN:  GENERATED FORTRAN ROUTINE
      SUBROUTINE MF0001
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
       include 'cfokkr.m'
       include 'cfreya.m'
       include 'crpltrgn.m'
       include 'crplocal.m'
       include 'crpnfile.m'
      COMMON/SCRATCH/ ZBUF(999)
C
C========================================
C
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1001 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >XZONI(J )
 1001 CONTINUE
      CALL MFWRT1(LUNPL1,'XZONI     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1002 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >XBOUNI(J )
 1002 CONTINUE
      CALL MFWRT1(LUNPL1,'XBOUNI    ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1003 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >AHALFS(J ,2)
 1003 CONTINUE
      CALL MFWRT1(LUNPL1,'RZON      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1004 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >AHALFS(J ,1)
 1004 CONTINUE
      CALL MFWRT1(LUNPL1,'RBOUN     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1005 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RMIDS(J ,2)
 1005 CONTINUE
      CALL MFWRT1(LUNPL1,'RMAJC     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1006 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RMIDS(J ,1)
 1006 CONTINUE
      CALL MFWRT1(LUNPL1,'RMAJB     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1007 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >(RMIDS(J ,1)-RMIDS(LEDGE,1))*USIL
 1007 CONTINUE
      CALL MFWRT1(LUNPL1,'SSZB      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1008 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >(RMIDS(J ,2)-RMIDS(LEDGE,1))*USIL
 1008 CONTINUE
      CALL MFWRT1(LUNPL1,'SSZC      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1009 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >ELONG(J ,1)
 1009 CONTINUE
      CALL MFWRT1(LUNPL1,'ELONG     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1010 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >TRIANG(J ,1)
 1010 CONTINUE
      CALL MFWRT1(LUNPL1,'TRING     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1011 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >DENT(J ,1)
 1011 CONTINUE
      CALL MFWRT1(LUNPL1,'DENT      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1012 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >VOLS(J ,1)
 1012 CONTINUE
      CALL MFWRT1(LUNPL1,'VOLB      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1013 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >AVI(J ,5,1)
 1013 CONTINUE
      CALL MFWRT1(LUNPL1,'GRHO1     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1014 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >AVI(J ,6,1)
 1014 CONTINUE
      CALL MFWRT1(LUNPL1,'GRHO2     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1015 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >(VOLS(J +1,1)-VOLS(J ,1))
 1015 CONTINUE
      CALL MFWRT1(LUNPL1,'DVOL      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1016 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >(AREAS(J +1,1)-AREAS(J ,1))
 1016 CONTINUE
      CALL MFWRT1(LUNPL1,'DAREA     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1017 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >DXBOUI(J )*UISL/AVI(J ,5,1)
 1017 CONTINUE
      CALL MFWRT1(LUNPL1,'DRAV      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1018 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >AVI(J ,3,1)*AVI(J ,5,1)*UISL**2
 1018 CONTINUE
      CALL MFWRT1(LUNPL1,'SURF      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1019 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOELS(2,J )
 1019 CONTINUE
      CALL MFWRT1(LUNPL1,'NE        ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1020 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOINS(2,J )
 1020 CONTINUE
      CALL MFWRT1(LUNPL1,'NI        ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1021 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOHS(1,2,J )
 1021 CONTINUE
      CALL MFWRT1(LUNPL1,'NH1       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1022 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOHS(2,2,J )
 1022 CONTINUE
      CALL MFWRT1(LUNPL1,'NH2       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1023 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOBIS(2,J )
 1023 CONTINUE
      CALL MFWRT1(LUNPL1,'BDENS     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1024 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOIS(1,2,J )
 1024 CONTINUE
      CALL MFWRT1(LUNPL1,'NIMP1     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1025 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >RHOIS(2,2,J )
 1025 CONTINUE
      CALL MFWRT1(LUNPL1,'NIMP2     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1026 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEH*TES(2,J )
 1026 CONTINUE
      CALL MFWRT1(LUNPL1,'TE        ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1027 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEH*TIS(2,J )
 1027 CONTINUE
      CALL MFWRT1(LUNPL1,'TI        ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1028 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >XZEFF(2,J )
 1028 CONTINUE
      CALL MFWRT1(LUNPL1,'ZEFFB     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1029 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >Q(J )
 1029 CONTINUE
      CALL MFWRT1(LUNPL1,'Q         ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1030 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEJ*AJZS(2,J )
 1030 CONTINUE
      CALL MFWRT1(LUNPL1,'CUR       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1031 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEJ*CJBEAM*AJBS(J )
 1031 CONTINUE
      CALL MFWRT1(LUNPL1,'CURB      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1032 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >VLOOPI(J ,1)*UISV*UISL
 1032 CONTINUE
      CALL MFWRT1(LUNPL1,'V         ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1033 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >DETEPR(J )*USIL**2/USIT
 1033 CONTINUE
      CALL MFWRT1(LUNPL1,'XETOT     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1034 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >XETHES(J )
 1034 CONTINUE
      CALL MFWRT1(LUNPL1,'XETHE     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1035 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THDRE(J )
 1035 CONTINUE
      CALL MFWRT1(LUNPL1,'CHDRE     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1036 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THIGE(J )
 1036 CONTINUE
      CALL MFWRT1(LUNPL1,'CHIGE     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1037 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THRBE(J )
 1037 CONTINUE
      CALL MFWRT1(LUNPL1,'CHRBE     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1038 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >DITIPR(J )*USIL**2/USIT
 1038 CONTINUE
      CALL MFWRT1(LUNPL1,'XITOT     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1039 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >XITHES(J )
 1039 CONTINUE
      CALL MFWRT1(LUNPL1,'XITHE     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1040 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THDRI(J )
 1040 CONTINUE
      CALL MFWRT1(LUNPL1,'CHDRI     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1041 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THIGI(J )
 1041 CONTINUE
      CALL MFWRT1(LUNPL1,'CHIGI     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1042 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THRBI(J )
 1042 CONTINUE
      CALL MFWRT1(LUNPL1,'CHRBI     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1043 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THRETI(J )
 1043 CONTINUE
      CALL MFWRT1(LUNPL1,'ETAI      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1044 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THETTH(J )
 1044 CONTINUE
      CALL MFWRT1(LUNPL1,'ETATH     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   2
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1045 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >THFITH(J )
 1045 CONTINUE
      CALL MFWRT1(LUNPL1,'FITH      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1046 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*(WEALFS(J )+WIALFS(J ))
 1046 CONTINUE
      CALL MFWRT1(LUNPL1,'PHALF     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1047 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*(WEAUXS(J )+WIAUXS(J ))
 1047 CONTINUE
      CALL MFWRT1(LUNPL1,'PAUX      ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1048 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*(WEBEMS(J )+WIBEMS(J ))
 1048 CONTINUE
      CALL MFWRT1(LUNPL1,'PBTOT     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1049 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*(WIBEMS(J ))
 1049 CONTINUE
      CALL MFWRT1(LUNPL1,'PBI       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1050 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*(WEBEMS(J ))
 1050 CONTINUE
      CALL MFWRT1(LUNPL1,'PBE       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1051 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*WEOHMS(J )
 1051 CONTINUE
      CALL MFWRT1(LUNPL1,'POH       ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1052 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*ZTHEAT(J )
 1052 CONTINUE
      CALL MFWRT1(LUNPL1,'AHEAT     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1053 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USEP*ZIMRAD(J )
 1053 CONTINUE
      CALL MFWRT1(LUNPL1,'PRADT     ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1054 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USIP*WEIONS(J )
 1054 CONTINUE
      CALL MFWRT1(LUNPL1,'PWALLE    ',ZBUF,INUM)
C  DEFINE, WRITE PROFILE FCN AT CURRENT TIME
      INX=   1
      IOFF=NWR(1,1,INX)-1
      INUM=NWR(2,1,INX)-IOFF
      DO 1055 J=NWR(1,1,INX),NWR(2,1,INX)
      ZBUF(J-IOFF)=
     >USIP*(WIIONS(J )+WICHXS(J ))
 1055 CONTINUE
      CALL MFWRT1(LUNPL1,'PWALLI    ',ZBUF,INUM)
      RETURN
      END
c@nfile   .../baldur/code/bald/doutput
c  rgb 14-aug-96 replaced fuzz with rndeps
c  rgb 26-feb-95 change zpbabs=bpabs*usep to zpbabs=zexte*usep
c  rgb 14-oct-91 locally computed scalars passed through common cnfile
c  rgb 11-jun-89 corrected zi and zbps
c
c**************************************************************************
c
c
        subroutine nfile(icall)
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'commhd.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'crplocal.m'
        include 'crpltrgn.m'
        include 'crpnfile.m'
c
        common /cnvect/ vftot(55,9), flxtot(55,9), srctot(55,9)
c
c********************************************************************
c
c       template for PLOTRGN to generate code for writing the
c       COMMON variables onto a direct access file
c
c       June 27, 1985  U. Schneider
c
c---------------------------------------------------------------------
c  all plotting output variables must be in common
c---------------------------------------------------------------------
c
c  local parameters:
c    imaxt-- max no. fcns of time
c
        parameter (imaxt=200)
c
c  output buffer
c
        real zbuf(imaxt)
c
c  local logical variable, always true:
c
        logical iall
        data iall/.true./
c------------------------------------------------------------------------------
c
c   icall=1: initialization call
c
        if(icall.ne.1) go to 10
c
c  CREATE AND OPEN TEXT FILE
c
      open (lunpl2, file=filpl2, status='unknown')
c
cbate      call create (lunpl2, filpl2, 2, -1, 9990)
c
 10     continue
        io=lunpl2
c
        imaxx=t1imx   ! nwr array dimension
        imaxsw=t1isw  ! nset array dimension
        if(icall.eq.1) go to 7999
c
c     *     *     *     *     *     *     *     *     *     *     *     *
c
c  local calculations for BALDUR plotting output
c   profile functions which are calculated locally must be in
c     COMMON/LOCALP/  to be passed to subroutine MF002 for output
c   scalar output quantities do not need to be placed in common
c
c     *     *     *     *     *     *     *     *     *     *     *     *
c
c       code taken from mprint, called 'sixth page' there
c
c
        if(nadump(1).le.lcentr) then
        isep = mzones
        isepm1 = ledge
        else
        isep = nadump(1)
        isepm1 = isep - 1
        end if
c
        zzee = 0.0
        zzei = 0.0
        zzeb = 0.0
        zzne = 0.0
        zznion = 0.0
        zohme = 0.0
        zalfe = 0.0
        zexte = 0.0
c
c               zvols is twice the actual plasma volume; it is used to
c               normalize integrals which involve dx2i(jz)
c
      zvols = 2. * vols(mzones,1)
cbate        zvols = 4.0 * (fcpi*rmins)**2 * rmajs
c
c  plasma surface area
c
      zsurfi = avi(mzones,3,1) * avi(mzones,5,1)
cbate        zsurfi = 4.0 * fcpi**2 * rmini * rmaji
c
c
        do 701 jz = lcentr,isepm1
        zzee = zzee + chi(lelec,jz) * dx2i(jz)
        zzei = zzei + chi(lion,jz) * dx2i(jz)
        zzeb = zzeb + hebems(jz)*rhobis(2,jz) * dx2i(jz)
        zzne = zzne + rhoels(2,jz) * dx2i(jz)
        zznion = zznion + rhoins(2,jz) * dx2i(jz)
        zohme = zohme + dx2i(jz) * weohms(jz)
        zalfe = zalfe + dx2i(jz) * (wealfs(jz) + wialfs(jz))
        zexte = zexte + dx2i(jz) * (webems(jz) + wibems(jz) +
     1      weauxs(jz)+wiauxs(jz)+weecrh(jz)+wiecrh(jz)+
     2      weicrf(jz)+wiicrf(jz))
  701   continue
c
        zzee = zzee * uisd * uiee
        zzei = zzei * uisd * uiee
        zzeb = zzeb * usee
        zohme = zohme * usep
        zalfe = zalfe * usep
        zexte = zexte * usep
c
c
cl              confinement times
c
c
        ipts=4
        zpts = 1.0 / float(ipts)
        i = 1
        call resetr(zzni,mxchi,0.0)
        call resetr(ztconf,ix0,0.0)
        call resetr(zloss ,mxchi,0.0)
        call resetr(zlne  , 6,0.0)
c
        zzlne = 0.0
        do 718 jz = lcentr, isepm1
        zeirs = 0.0
        if (mimp.le.0) go to 705
        do 704 ji = 1, mimp
        zeirs = zeirs + weirs(ji,jz)
  704   continue
  705   continue
c
        zloss(lelec) = zloss(lelec) + usip*zvols*dx2i(jz)*
     1                          (weions(jz) + webrs(jz) + zeirs)
        zloss(lion) = zloss(lion) - usip*zvols*dx2i(jz)*
     1                          (wiions(jz) + wichxs(jz))
c
c
c
        z0 = zvols * uisd * dx2i(jz)
        do 712 jp = 1, mchi
        zzni(jp) = zzni(jp) + z0*chi(jp,jz)
  712   continue
c
        zzlne = zzlne + rhoels(2,jz) * dx2i(jz) * 2.
c
        if (jz.lt.isepm1.and.
     1  (xbouni(jz+1)/xbouni(isep)+rndeps).lt.(float(i)*zpts)) go to 718
c
c  half-width in external units
c
        zrconf(i) = avi(jz+1,15,1) * uiel
cbate        zrconf(i) = xbouni(jz+1) * rmins * usel
        zeloss = 0.0
        zlne(i) = zzlne
c
        do 716 jp = 1, mchi
c
        zflux = 0.0
        do 714 jp2 = 1, mchi
        zflux = zflux - aaaa(jp,jp2,jz+1)*chi(jp2,jz)
     1                - bbbb(jp,jp2,jz+1)*chi(jp2,jz+1)
  714   continue
c
c                     surface area V'(xi)*<|del xi|> internal units
        z0 = (zloss(jp) + zflux*avi(jz+1,3,1)*avi(jz+1,5,1)) * ueit
cbate        z0 = (zloss(jp) + zflux*zsurfi*xbouni(jz+1)) * ueit
        if (z0.gt.epslon) ztconf(jp,i) = zzni(jp) / z0
        if (jp.eq.lelec.or.jp.eq.lion) zeloss = zeloss + z0
  716   continue
c
        if (zeloss.gt.epslon) ztconf(mchi+1,i) =
     1  (zzni(lelec) + zzni(lion)) / zeloss
        zlne(i) = zlne(i) * ztconf(mchi+1,i) / xbouni(jz+1)**2
        i = i + 1
  718   continue
c
c
c       experimental times
c
        ztenr1 = 0.0
        if (zohme.gt.epslon)
     1          ztenr1 = (zzei + zzee) / zohme * usep*uset*uese
        ztenr2 = 0.0
        z0 = zohme + zexte + zalfe
        if (z0.gt.epslon)
     1          ztenr2 = (zzei + zzee + zzeb) / z0 * usep*uset*uese
        if (z0.gt.epslon)
     1          ztaueb = (zzei + zzee) / z0 * usep*uset*uese
c
c       average ni, ne, ti, te
c
        z01 = 2.0 * used /xbouni(isep)**2
        znibar = zznion * z01
        znebar = zzne * z01
        ztibar = zzei / zznion * ueie*uieh
        ztebar = zzee / zzne * ueie*uieh
c
c       loop voltage
c
        zv = vloopi(ledge,2)
        zvc = vloopi(lcentr,2)
c
cbate        zloop = (ajzs(2,isepm1)-cjbeam*ajbs(isepm1))*eta(2,isepm1)*usev
cbate        zv = zloop* 2 * fcpi * rmajs
cbate        zloopc = (ajzs(2,lcentr)-cjbeam*ajbs(lcentr))*eta(2,lcentr)*usev
cbate        zvc = zloopc* 2 * fcpi * rmajs
c
c       beta -- e, i, beam, total
c
        zi = avi(isep,2,1) * avi(isep,3,1) * r0ref * bpoli(isep)
        zbps=emu0*zi/(twopi*avi(isep,15,1))*uisb
        zpb = zbps**2 / (8.0*fcpi) * used*usee
cbate        zpb = bpols(1,isep)**2 / (8.0*fcpi) * used*usee
        z1 = z01 / zpb
        zbetae = 0.66666666 * zzee * z1
        zbetai = 0.66666666 * zzei * z1
        zbetab = 0.66666666 * zzeb * z1
c
        zbeta = zbetae + zbetai + zbetab
c
c               total beta (include bz)
c
        z1 = (zbps / bzs)**2
cbate        z1 = (bpols(1,isep) / bzs)**2
        ztbete = zbetae * z1
        ztbeti = zbetai * z1
        ztbetb = zbetab * z1
        ztbeta  = zbeta  * z1
c
c       compute alpha beta
c
        zalfap = 0.0
        do 730 jz = lcentr,isepm1
        zalfap = zalfap + alphai(jz) * ealfai(jz) * dx2i(jz)
  730   continue
c
        zalfap = zalfap * uisd * uise * .6666666666 * 2.
     1  * 8. * fcpi /zbps**2
cbate     1  * 8. * fcpi /(bpols(1,isep)*xbouni(isep))**2
        zalfat = zalfap * z1
        zbeta = zbeta + zalfap
        ztbeta = ztbeta + zalfat
c
c       internal inductance and lambda
c
c  ...Note: lambda here is just beta-pol + l-i/2 !
c  ...l-i = alint is computed in subroutine AVEPLA.
c
        zsum=0.
        zflux=0.
        do 732 jz=lcentr,isepm1
c  volume integral of bpols**2
        zsum=zsum+bpols(2,jz)**2*dx2i(jz)
        zflux = zflux + bpols(2,jz)*(avi(jz + 1,1,1) - avi(jz,1,1))
cbate        zflux=zflux+bpols(2,jz)*dxzoni(jz)
732     continue
        zlint=2.*zsum/bpols(1,isep)**2
        if (leqtyp.eq.0) alint=zlint
        zlambd=zbeta+alint/2.
c
c       1.e-8 converts cm**2*gauss to m**2*tesla
c
        zflux = 1.e-8 * 2. * fcpi * r0ref * zflux * uisl**2
cbate        zflux=1.e-8*2.*fcpi*rmajs*rmins*zflux
c
c       compute line averaged density
c
        zlined=0.0
        do 735 jz= lcentr , ledge
        zlined = zlined + rhoels(2,jz)
     1                     * (ahalfs(jz+1,1) - ahalfs(jz,1))
cbate        zlined = zlined + rhoels(2,jz) * dxzoni(jz)
 735    continue
        zlined = zlined / ahalfs(isep,1)
cbate        zlined = zlined / xbouni(isep)
c
  737   continue
c
c     *     *     *     *     *     *     *     *     *     *     *     *
c
c       code taken from mprint: third page-energy fluxes
c
        zohme = 0.0
        zchxe = 0.0
        zalfe = 0.0
        zrade = 0.0
        zexte = 0.0
        zeion = 0.0
c
        zthdot = 0.0
c
        do 455 jb=lcentr+1,mzones
c
        zrad = avi(jb,15,1) * uiel
cbate        zrad = xbouni(jb) * rmini * uiel
c
c
cl              compute fluxes
c
c
c  note <grad f> = (<|grad xi|**2>/<|grad xi|>) * d f / d xi,
c  as in sbrtn CONVRT (dps 09oct87)
c
        z0 = avi(jb,6,1) / (avi(jb,5,1)*dxboui(jb) * uisl)
c
cbate        z0 = 1.0 / (dxboui(jb) * rmins)
c
        zgrdte = (tes(2,jb) - tes(2,jb-1))*z0
        zgrdti = (tis(2,jb) - tis(2,jb-1))*z0
c
c
c               energy fluxes (total)
c
c
        zeeflx = 0.0
        zeiflx = 0.0
c
        do 422 jp = 1, mchi
        zeeflx = zeeflx - aaaa(lelec,jp,jb) * chi(jp,jb-1) -
     1                          bbbb(lelec,jp,jb) * chi(jp,jb)
        zeiflx = zeiflx - aaaa(lion,jp,jb) * chi(jp,jb-1) -
     1                          bbbb(lion,jp,jb) * chi(jp,jb)
  422   continue
c
        zeeflx = zeeflx * avi(jb,3,1)*avi(jb,5,1) * uiep
        zeiflx = zeiflx * avi(jb,3,1)*avi(jb,5,1) * uiep
c
cbate        zeeflx = zeeflx * zsurfi*xbouni(jb) * uiep
cbate        zeiflx = zeiflx * zsurfi*xbouni(jb) * uiep
c
c
c               conductivities
c
c
        zicond=-ditis(jb) *zgrdti*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
        zecond=-detes(jb) *zgrdte*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
        zineut=-ditins(jb)*zgrdti*avi(jb,3,1)*avi(jb,5,1)*uisl**2*usep
c
cbate        zicond = - ditis(jb)  * zgrdti * zsurfi*xbouni(jb)*uisl**2*usep
cbate        zecond = - detes(jb)  * zgrdte * zsurfi*xbouni(jb)*uisl**2*usep
cbate        zineut = - ditins(jb) * zgrdti * zsurfi*xbouni(jb)*uisl**2*usep
c
        zeconv = zeeflx - zecond
        ziconv = zeiflx - zicond - zineut
c
  430   continue
c
c
c               cumulative sources/sinks
c               (sources/sinks integrated from 0 to r(jb))
c
c
        z1ohm = 0.0
        z1chx = 0.0
        z1alf = 0.0
        z1rad = 0.0
        z1ext = 0.0
        z1eion = 0.0
c
        if(jb.le.lcentr) go to 440
        jz = jb - 1
c
        z1ohm = z1ohm + dx2i(jz) * weohms(jz)
        z1chx = z1chx - dx2i(jz)*(wiions(jz) + wichxs(jz))
        z1alf = z1alf + dx2i(jz) * (wealfs(jz) + wialfs(jz))
        z1ext = z1ext + dx2i(jz) * (webems(jz) + wibems(jz) +
     1      weauxs(jz)+wiauxs(jz)+weecrh(jz)+wiecrh(jz)+
     2      weicrf(jz)+wiicrf(jz))
        z1eion=z1eion+dx2i(jz)*cnueqs(jz)*(tes(2,jz)-tis(2,jz))
        z1rad = z1rad+dx2i(jz) * (webrs(jz)+weions(jz)+wesrs(jz))
c
        if (mimp.le.0) go to 438
c
        do 436 ji = 1, mimp
        z1rad = z1rad + dx2i(jz) * weirs(ji,jz)
  436   continue
c
  438   continue
c
        zohme = zohme + z1ohm*zvols*usep
        zchxe = zchxe + z1chx*zvols*usep
        zalfe = zalfe + z1alf*zvols*usep
        zexte = zexte + z1ext*zvols*usep
        zrade = zrade + z1rad*zvols*usep
        zeion = zeion + z1eion*zvols*usep
c
c
  440   continue
c
        ztote = zohme + zalfe + zexte - zecond - zicond - zeconv -
     1          ziconv - zineut - zchxe - zrade
        if(jz.eq.isepm1) zthdot=ztote
c
        zineut = zineut + zchxe
c
  455   continue
c
c               fusion power multiplication
c
c
c               qtech (for beam heating only)
c
c       compute the total injected beam power
        zpowr = 0.0
        do 752 jb = 1, mhbeam
        if (libeam(jb).le.0) go to 752
        ib = libeam(jb)
        do 751 jfr = 1, mxhfr
        zpowr = zpowr +
     1          yfract(ib,jfr) * hibeam(jb)*hebeam(jb) / float(jfr)
  751   continue
  752   continue
        zpowr = zpowr * 10.0**(-fxes)/fces * uesi * uesh * usep
c
c       integrate out to separatrix to get total powers
c
        zphbem=0.
        zphrf=0.
        zphalf=0.
        zpath=0.
        zpathl=0.
        zpabt=0.
        zpabtl=0.
        do 754 jz=lcentr,isepm1
        zphbem=zphbem+dx2i(jz)*(webems(jz)+wibems(jz))
        zphrf=zphrf+dx2i(jz)*(weauxs(jz)+wiauxs(jz)+
     1  weicrf(jz)+wiicrf(jz)+weecrh(jz)+wiecrh(jz))
        zphalf=zphalf+dx2i(jz)*(wealfs(jz)+wialfs(jz))
        zpath=zpath+dx2i(jz)*afuses(jz)
        zpathl=zpathl+aoloss(jz)*dx2i(jz)*afuses(jz)
        zpabt=zpabt+dx2i(jz)*halfas(jz)
        zpabtl=zpabtl+aoloss(jz)*dx2i(jz)*halfas(jz)
754     continue
        zpbinj=bpinjs*usep
        zpbabs=zexte*usep
        zpblos=bploss*usep
c
        zvols = 2.0 * vols(mzones,1)
cbate        zvols=4.*rmajs*(fcpi*rmins)**2
        zphbem=zphbem*zvols*usep
        zphrf=zphrf*zvols*usep
        zphalf=zphalf*zvols*usep
        zpath=zpath*zvols*usep*efusei*uise
        zpathl=zpathl*zvols*usep*efusei*uise
        zpabt=zpabt*zvols*usep*efusei*uise
        zpabtl=zpabtl*zvols*usep*efusei*uise
c
c       calculate loss fractions
c
        if(zpbinj.ne.0.) then
        zlbem=(zpbinj-zpbabs+zpblos)/zpbinj
        else
        zlbem=0.
        if(zpblos+zphbem.ne.0.) zlbem=zpblos/(zpblos+zphbem)
        endif
        zlrf=0.
        if(cfutz(60).lt.1.) zlrf=max(0.,cfutz(60))
        zprfin=zphrf/(1.-zlrf)
        if(zpath.ne.0.) then
        zlath=zpathl/zpath
        else
        zlath=0.
        endif
        if(zpabt.ne.0.) then
        zlabt=zpabtl/zpabt
        else
        zlabt=0.
        endif
c
c       calculate stored energy time derivatives
c
        if(zpbinj.ne.0.) then
        zebdot=zpbinj*(1.-zlbem)-zphbem
        else
        zebdot=-(zpblos+zphbem)
        endif
        zeadot=zpath*(1.-zlath)+zpabt*(1.-zlabt)-zphalf
c
c       calculate the steady state beam power
c
        if(zpbinj.ne.0.) then
        zpbnow=zpbinj
        else
        zpbnow=zpblos+zphbem
        endif
        if(zphbem.gt.0.) then
        zpbstr=(zphbem+zphrf-zthdot-zeadot+zpabt*(1.-zlabt)
     1  -zprfin*(1.-zlrf))
     1  /((1.-zlbem)*(1.+zpabt*(1.-zlabt)/zphbem))
        else
        zpbstr=-1.
        endif
c
c       calculate the steady state rf power
c
        if(zpbstr.gt.0.) then
        zprstr=zprfin
        else
        zpbstr=0.
        zprstr=(zphbem+zphrf-zthdot-zeadot+zpabt*(1.-zlabt))
     1  /(1.-zlrf)
        if(zprfin.le.0.) zprstr=0.
        endif
c
c       calculate the q values
c
        zpss=zpbstr+zprstr
        zpcomp=zphbem+zpblos+zprfin-zthdot-zeadot
        if(zpbnow+zprfin.ne.0.) then
        zqinst=5.*(zpath+zpabt)/(zpbnow+zprfin)
        else
        zqinst=0.
        endif
        zbfact=0.
        if(zphbem.ne.0.) zbfact=zpbstr*(1.-zlbem)/zphbem
        if(zpss.ne.0.) then
        zqss=5.*(zpath+zpabt*zbfact)/zpss
        else
        zqss=0.
        endif
        if(zpcomp.ne.0.) then
        zqcomp=5.*(zpath+zpabt)/zpcomp
        else
        zqcomp=0.
        endif
        zheat=zphbem+zphrf
        if(zheat.ne.0.) then
        zqdenm=(zpowr+zprfin)*(1.-zthdot/zheat)
        else
        zqdenm=0.
        endif
        if(zqdenm.ne.0.) then
        zqtech=5.*(zpath+zpabt)/zqdenm
        else
        zqtech=0.
        endif
c
        zqss=max(0.,zqss)
        zqinst=max(0.,zqinst)
        zqcomp=max(0.,zqcomp)
        zqtech=max(0.,zqtech)
        if(zheat.ne.0.) zptech=zpowr*(1.-zthdot/zheat)
c
c     *     *     *     *     *     *     *     *     *     *     *     *
c
c       integrate the toroidal currents and sum the impurity radiation
c
        zip=0.
        zib=0.
        do 1110 jz=lcentr,ledge
        zip = zip + (areas(jz+1,1)-areas(jz,1)) * ajzs(2,jz)
        zib = zib + (areas(jz+1,1)-areas(jz,1)) * ajbs(jz)
cbate        zip=zip+dx2i(jz)*ajzs(2,jz)
cbate        zib=zib+dx2i(jz)*ajbs(jz)
1110    continue
cbate        zip=2.*fcpi*rmins**2*zip
cbate        zib=2.*fcpi*rmins**2*zib
c
c  average voltage from total ohmic heating and current
c
        zva=zohme/(1000.*usei*zip)
c
c  average heating from sawtooth reconnection
c
        if(sawtau.ne.0.) then
        zpcrsh=-estpol/sawtau
        else
        zpcrsh=0.
        endif
c
c***************************
c
7999    continue
c
c
c***************************
c
c
c***************************
c
C**> INSERTING GENERATED CODE; TEMPLATE FILE WAS:
C**>  BALPLOT:NFILET.FOR                                             
C >>> PLTRGN GENERATED FORTRAN STARTS HERE:
      IF(ICALL.NE.1) GO TO 8000
C---------------------
C SET UP SWITCHING; FORTRAN GEN. BY SUBROUTINE PLSWCH
      CALL PSWITCH(IMAXSW,IER)
      IF(IER.NE.0) GO TO 9300
C---------------------
C-----------------------------
C   NOW DEFINE THE FUNCTIONS TO BE PLOTTED (GEN. BY PLTRGN)
 8000 CONTINUE
      IF(ICALL.EQ.1) GO TO 10000
C  CALCULATE A RECORD FOR THE F(T) FILE:
      IPTR=1
      ZBUF(IPTR)=
     >   UIST*TAI                                                     
C----------
      CALL NFIX01(ZBUF,IMAXT,IPTR,IER)
      IF(IER.NE.0) GO TO 9600
      WRITE(IO,101) (ZBUF(JPTR),JPTR=1,IPTR)
 101  FORMAT (1X,6(1PE12.5))
c
c
c***************************
c
10000   continue
c
        return
c----
c  errors
c  branched to from generated fortran
c
 9600   continue
        call mesage(
     >' fatal error in NFILEX-- plot output buffer     ')
        call mesage(
     >' too small-- increase local parameter IMAXT in  ')
        call mesage(
     >' the file NFILET.FOR and rerun @PLOTGEN and your')
      call mesage(
     >' BALDUR run.                                    ')
c
c
        stop
c
c
 9300   continue
        call mesage(
     >' fatal error in NFILEX-- plot switching buffer  ')
        call mesage(
     >' too small-- increase array dimension T1ISW in  ')
        call mesage(
     >' the file pltrgn.COM and rerun @PLOTGEN and your')
      call mesage(
     >' BALDUR run.                                    ')
c
c
        stop
c
 9990 call abortb (6,'error creating output file in sbrtn nfile')
      stop
c
                end
c--------1---------2---------3---------4---------5---------6---------7-c
      SUBROUTINE NFIX01(ZBUF,IMAXT,IPTR,IER)
C
C   THIS IS CALLED BY NFILE -- BREAKS UP
C   THE NF#### CALLS TO AVOID A POTENTIAL
C   COMPILER PROBLEMS (TOO MANY EXTERNALS)
C========================================
C      THIS IS PLTRGN GENERATED CODE
C========================================
C
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
       include 'cfokkr.m'
       include 'cfreya.m'
       include 'crpltrgn.m'
       include 'crplocal.m'
       include 'crpnfile.m'
C
      REAL ZBUF(IMAXT)
C===========
      IER=0
      IF(NSET(   2).EQ.1) THEN
        CALL NF0001(ZBUF,IMAXT,IPTR,IER)
        IF(IER.NE.0) RETURN
      ENDIF
C
       RETURN
       END
C------------------/PLTRGN:  GENERATED FORTRAN ROUTINE
      SUBROUTINE NF0001(ZBUF,IMAXT,IPTR,IER)
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
       include 'cfokkr.m'
       include 'cfreya.m'
       include 'crpltrgn.m'
       include 'crplocal.m'
       include 'crpnfile.m'
C
      REAL ZBUF(IMAXT)
C===========
      IER=0
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >USEI*ZIP
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >USEI*ZIB
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RMAJS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RMINS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ELONG(MZONES,1)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >TRIANG(MZONES,1)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >DENT(MZONES,1)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >VOLS(MZONES,1)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >AREAS(MZONES,1)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >0.1*BZS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOELS(2,LCENTR)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOINS(2,LCENTR)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >USEH*TES(2,LCENTR)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >USEH*TIS(2,LCENTR)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >XZEFF(2,LCENTR)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >UIST*DTI
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ADDS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >HDDS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ADDS+HDDS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >EVSINV*TES(2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >EVSINV*TIS(2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOELS(2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOHS(1,2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOHS(2,2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOIS(1,2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >RHOIS(2,2,MZONES)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTEBAR
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTIBAR
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZNEBAR
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZNIBAR
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZLINED
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZLAMBD
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZV
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZVC
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZVA
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTCONF(MCHI+1,4)
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTAUEB
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTBETE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTBETI
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTBETB
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZALFAT
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZTBETA
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZBETAE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZBETAI
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZBETAB
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZALFAP
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZBETA
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ERGTS(MZONES)*USEE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ERGES(MZONES)*USEE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ERGIS(MZONES)*USEE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZOHME
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >OHPTOT
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPCRSH
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZEXTE+ZALFE+ZOHME
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >-ZRADE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPHALF
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPATH+ZBFACT*ZPABT
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPATH
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPABT
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPHBEM
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPBSTR
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPBNOW
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPBABS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZPTECH
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZQINST
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZQSS
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZQTECH
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >ZQCOMP
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWTE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWNE
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWDD
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWR1
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWRI
      IPTR=IPTR+1
      IF(IPTR.GT.IMAXT) THEN
        IER=1
        RETURN
      ENDIF
      ZBUF(IPTR)=
     >SAWRMX
      RETURN
      END
c@plabel   .../baldur/code/bald/doutput.f
c  rgb 21-nov-96 changed RPLOT output file names
c    tfiln   'tf' // runid(1:6) --> runid(1:8) // 'TF.PLN'
c    filpl1  'mf' // runid(1:6) --> runid(1:8) // 'XF.PLN'
c    filpl2  'nf' // runid(1:6) --> runid(1:8) // 'YF.PLN'
c  rgb 07-aug-96 removed idate and call date
c  rgb 01-jul-91 fixed fatal dupplication of labels messages
c      replacing goto 99 statements as needed
c  rgb 05-oct-89 unitsb(ib)=unitsr(...) if iintb(ib)=0
c      or unitsb(ib)=unitst(...) if iintb(ib)=1, otherwise abort
c  rgb 02-oct-89 xlab, labelb, and labelr are now 1-D arrays
c  rgb 13-jun-89 added a return statement after call close
c  rgb 12-jun-89 added argument 9990 to call create for civic
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  plabel
c generate label file for baldur plotting output
c
c
c   this subroutine is called in a BALDUR run, after all namelists
c   have been read.
c   PLABEL initializes all necessary variables as I/O unit numbers,
c   file names, etc. from its given arguments.
c   Then PLABEL calls the generated subroutine TFILSN which sets up
c   the variables in cplotr.blk and constructs the mapping file for
c   BALDUR plotting output.
c
c
c   see program plotr which accesses the label file and plot data
c  files and allows interactive or batch plotting of baldur plot data
c
c
c
c
      subroutine plabel (run)
c
c
      include 'crplot.m'
      include 'crpltrgn.m'
c
c
c  integer RUN is passed; it contains the RUN-ID in positions 2:7
c
cbate960807      character*8 run, idate
      character *(*) run
      character*8 msg(6)
      data msg/6*'        '/
c
c  initialize COMMON variables
c
        lunpl1 = 23
        lunpl2 = 24
c
c..find first 8 non-blank characters in run name
c
      i1 = 0
      i2 = 0
      imax = min ( len(run), 32 )
c
      if ( imax .lt. 1 ) then
c
        runid(1:6) = 'baldur'
        filpl1 = 'baldur01XF.PLN'
        filpl2 = 'baldur01YF.PLN'
        tfiln  = 'baldur01TF.PLN'
c
      else
c
        do j=1,imax
          i1 = j
          if ( run(j:j) .ne. ' ' ) go to 12
        enddo
 12     continue
c
        i2 = i1
        do j=i1,imax
          if ( run(j:j) .eq. ' ' ) go to 14
          i2 = j
        enddo
 14     continue
c
        i2 = min ( i2, i1 + 7 )
        imax = min ( i2, i1 + 5 )
        runid = run(i1:imax)
c
        filpl1 = run(i1:i2) // 'XF.PLN'
        filpl2 = run(i1:i2) // 'YF.PLN'
        tfiln  = run(i1:i2) // 'TF.PLN'
c
      endif
c
c
cbate        runid(1:6) = run(2:7)
c
c  prepare file names for MFILE and NFILE
c
cbate         filpl1(1:2) = 'mf'
cbate         filpl1(3:8) = runid(1:6)
c
cbate         filpl2(1:2) = 'nf'
cbate         filpl2(3:8) = runid(1:6)
c
c
c  set open flag, and get date from system
c
cbate960807        call date(idate)
c
c
c  prepare file name for mapping file
c
cbate        tfiln(1:2) = 'tf'
cbate        tfiln(3:8) = runid(1:6)
        lun = lunpl1
c
        filnc='a'
c
c  create and open mapping file for TFILSN (CRAY only)
c
      open (lun, file=tfiln, status='unknown')
c
cbate        call create(lun,tfiln,2,-1,9990)
c
c  call tfilsn-- 2 passes required-- to fill label information
c
        do 1000 icall=1,2
c
        call tfilsn(icall,nzones,
     >   naxxvr,nxr,nzonex,nrecx,nfx,xlab,xndabb,
     >   naxmgp,nbal,infb,ifunb,iintb,labelb,abb,unitsb,
     >   naxfot,nft,labelt,abt,unitst,
     >   naxfxt,nfxt,labelr,abr,unitsr,itypr)
c
 1000   continue
c
c  eliminate x axes with zero fcns associated to them
c
        ixr=0
 1010   continue
c
        ixr=ixr+1
        if(ixr.gt.nxr) go to 1045
        if(nfx(ixr).gt.0) go to 1010
        ixp1=ixr+1
        if(ixr.eq.nxr) go to 1040
c
c  downshift package discriptor arrays to elimate 0 length entry
c
        do 1030 ix=ixp1,nxr
c
        ixm1=ix-1
        nfx(ixm1)=nfx(ix)
        xndabb(ixm1)=xndabb(ix)
        nrecx(ixm1)=nrecx(ix)
        nzonex(ixm1)=nzonex(ix)
c
        xlab(ixm1)=xlab(ix)
c
        do 1030 ifxt=1,nfxt
c
        if(itypr(ifxt).eq.ix) itypr(ifxt)=ixm1
c
 1030   continue
c
c------
 1040   continue
c
        nxr=nxr-1
        ixr=ixr-1
        go to 1010
c------
 1045   continue
c
c  eliminate any zero-length multigraph packages
c
        ibal=0
 1050   continue
        ibal=ibal+1
        if(ibal.gt.nbal) go to 1090
        if(infb(ibal).gt.0) go to 1050
        ibp1=ibal+1
        if(ibal.eq.nbal) go to 1080
c
c  downshift package discriptor arrays to elimate 0 length entry
c
        do 1070 ib=ibp1,nbal
c
        ibm1=ib-1
        infb(ibm1)=infb(ib)
        abb(ibm1)=abb(ib)
        iintb(ibm1)=iintb(ib)
c
        labelb(ibm1)=labelb(ib)
        unitsb(ibm1)=unitsb(ib)
c
        do 1070 i=1,15
        ifunb(i,ibm1)=ifunb(i,ib)
 1070   continue
c
c------
 1080   continue
c
        nbal=nbal-1
        ibal=ibal-1
c
        go to 1050
c------
 1090   continue
c
c  define units of multigraph packages
c
        if(nbal.gt.0) then
c
          do 1095 ib=1,nbal
c
            if ( iintb(ib) .eq. 0 ) then
              unitsb(ib)=unitsr(iabs(ifunb(1,ib)))
            elseif ( iintb(ib) .eq. 1 ) then
              unitsb(ib)=unitst(iabs(ifunb(1,ib)))
            else
              call abortb (6,
     & 'abort from sbrtn plabel after do 1095.  iintb(ib).ne. 0 or 1')
            endif
c
 1095     continue
        endif
c
c
c  check that all nontemporal independant coordinates are
c  defined; calculate mf file offsets for each fcn and indep.
c  coordinate
c
        ioff0=2
c
        do 1150 j=1,nfxt
        nrofff(j)=ioff0
        ioff0=ioff0+nrecx(itypr(j))
        jm1=j-1
        if(jm1.eq.0) go to 1140
c
        do 1130 i=1,jm1
        if(abr(i).ne.abr(j)) go to 1130
        msg(1)(1:5)=abr(i)(1:5)
        call mesage(msg)
        call mesage(
     >  '  fatal error. Abbreviation used more than once.')
        stop
 1130   continue
 1140   continue
c
c  check for duplication of abbreviations between scalar and profile
c  functions
c
        do 1145 i=1,nft
        if(abr(j).ne.abt(i)) go to 1145
        msg(1)(1:5)=abr(j)(1:5)
        call mesage(msg)
        call mesage(
     >  '  fatal error. Abbreviation used more than once.')
        stop
 1145   continue
c
        do 1150 i=1,nxr
        if(xndabb(i).eq.abr(j)) nroffx(i)=nrofff(j)
 1150   continue
c
c  check for duplication of scalar abbreviations
c
        do 1155 j=2,nft
        jm1=j-1
        do 1155 i=1,jm1
        if(abt(i).ne.abt(j)) go to 1155
        msg(1)(1:5)=abt(i)(1:5)
        call mesage(msg)
        call mesage(
     >  '  fatal error. Abbreviation used more than once.')
        stop
c
 1155   continue
c
c  calc tot number of records not. incl. time record:
        nfr=ioff0-2
c
c  check that definitions were found for all x axes
c
        do 1160 i=1,nxr
c
        if(nroffx(i).ne.0) go to 1160
        msg(1)(1:5)=xndabb(j)(1:5)
        call mesage(msg)
        call mesage(
     >  ' ? fatal error: X-coordinate undefined ')
c
        stop
 1160   continue
c
c  write the tf file
c
        nlxvar=.true.
        call tfilwr(lun)
        close(lun)
        return
c
 9990 call abortb (6,'error opening output file in sbrtn plabel')
      stop
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
C--------------------------------
C  PSWITCH -- PLOT OUTPUT LOGIC ROUTINE **GENERATED CODE**
C    DO NOT MODIFY -- SEE GENERATOR PROGRAM PLTRGN AND
C    PLTRGN SUBROUTINE PLSWCH
C
      SUBROUTINE PSWITCH(IMAXSW,IER)
C
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
       include 'cfokkr.m'
       include 'cfreya.m'
       include 'crpltrgn.m'
       include 'crplocal.m'
       include 'crpnfile.m'
C
      LOGICAL IALL
C------------------------
      IALL=.TRUE.
      ISW=1
      NSET(ISW)=1
C
C  SET NEXT SWITCH
      ISW=   2
      IF(ISW.GT.IMAXSW) GO TO 9300
      IF(
     >IALL
     >  )  NSET(ISW)=1
C==================
C  NORMAL EXIT
      IER=0
      RETURN
C  ERROR EXIT
 9300 CONTINUE
      IER=1
      RETURN
      END
c@pltavg
c
        subroutine pltavg(prin,pravg,pmax,pwavg,kel,kfail)
c
c       finds max. element and calculates the weighted average of array prin
c
c               argument index
c
c
c       prin    input array
c       pravg   weighting array for doing the weighted average
c                       for h(r) rescaling ths is dx2i(lcentr)
c       pmax    the maximum of the absolute values of all elements
c       pwavg   the weighted average of prin
c       kel     number of elements in prin, pravg, and prout
c       kfail   problem flag (reset to 0 at the start)
c
c
c
c
c
        dimension prin(*),pravg(*)
c
c
c
c       loop over the input array
c
        kfail=0
        pmax=1.
        pwavg=1.
c
        zmax=0.
        zsum=0.
        zsavg=0.
        do 100 j=1,kel
        if(abs(prin(j)).gt.zmax) zmax=abs(prin(j))
        zsum=zsum+prin(j)*pravg(j)
        zsavg=zsavg+pravg(j)
100     continue
c
        if(zmax.gt.0.) then
        pmax=zmax
        else
        go to 901
        end if
c
        if(zsavg.ne.0.) then
        zavg=zsum/zsavg
        else
        go to 902
        end if
        if(zavg.eq.0.) go to 903
c
        pwavg=zavg
c
        return
c
c       error section
c
901     continue
c               maximum value is zero
        kfail=1
        go to 999
902     continue
c               average value of weighting function pravg is zero
        kfail=2
        go to 999
903     continue
c               average value of input array prin is zero
        kfail=3
        go to 999
c
c
999     continue
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@redim
c
        subroutine redim(pin,pout,k1dim,k1,k2,kel)
c
c       repacks the two dimensional array pin into the one-d array pout
c
c               argument index
c
c       pin     input array which is to be repacked into pout
c       pout    output array
c       k1dim   dimension of first index of pin
c       k1      value of first index of pin to use in repacking
c       k2      starting value of second index of pin
c       kel     number of elements to repack
c
c
        dimension pin(k1dim,*),pout(*)
c
c
        do 100 j=k2,k2+kel-1
        pout(j)=pin(k1,j)
100     continue
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@rescal
c
        subroutine rescal(kscal,prin,pravg,prout,kel,kfail)
c
c
c               input index
c
c
c       kscal   rescaling index         1       maximum element = 1.0
c                                       2       h(r) rescaling
c       prin    input array which is to be rescaled
c       pravg   input array for the h(r) rescaling (usually dx2i(lcentr))
c       prout   output array with the rescaled elements
c       kel     number of elements in prin, pravg, and prout
c       kfail   problem flag (reset to 0 at the start)
c
c
c
c
c
        dimension prin(*),pravg(*),prout(*)
c
c
c       do first pass over the input array
c
        kfail=0
c
c
        zmax=0.
        zsum=0.
        zsavg=0.
        do 100 j=1,kel
        if(abs(prin(j)).gt.zmax) zmax=abs(prin(j))
        zsum=zsum+prin(j)*pravg(j)
        zsavg=zsavg+pravg(j)
100     continue
c
c       put the rescaled array in prout
c
c               scaled to maximum element
        if(kscal.ne.1) go to 200
        if(zmax.eq.0.) go to 901
        do 150 j=1,kel
        prout(j)=abs(prin(j))/zmax
150     continue
c
c               scaled to weighted average (as h(r))
c
200     continue
        if(kscal.ne.2) go to 300
        if(zsavg.ne.0.) then
        zavg=zsum/zsavg
        else
        go to 902
        end if
        if(zavg.eq.0.) go to 903
        do 250 j=1,kel
        prout(j)=prin(j)/zavg
250     continue
c
300     continue
c
        return
c
c       error section
c
901     continue
c               maximum value is zero
        kfail=1
        go to 999
902     continue
c               average value of weighting function pravg is zero
        kfail=2
        go to 999
903     continue
c               average value of input array prin is zero
        kfail=3
        go to 999
c
c
999     continue
c               zero the output array prout
        call resetr(prout,kel,0.0)
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@mf0002
c  rgb 09-jun-89 1-1/2-D BALDUR  modifications
c
c--------1---------2---------3---------4---------5---------6---------7-c
c@tfilsn
c  rgb 25-mar-02 changed dimension to iuntsb(*),iuntsr(*),iuntst(*)
c  rgb 02-oct-89 changed char*10 to char*30, char*5 to char*10
c    iklblx, ilablm, ilablt, ilablr changed to 1-D arrays
c    implemented iuntsb,iuntsr,iuntst as 1-D array char*16
c
c---------------------------------------------------------
c  tfils0
c  prototype module for generation of labels/ordering information
c  for plotting output
c  the >>>>>>>>>>>>>>>>>>>>>> is replaced by generated code
c  by the program PLTRGN
c  and the resulting program is named "tfilsn.for"
c
        subroutine tfilsn(icall,irecsz,
     >     imaxx,inumx,izonex,inumr,itotf,iklblx,ikabrx,
     >     imaxbl,imul,infb,ifunb,iintb,ilablm,iabbm,iuntsb,
     >     imaxft,ift,ilablt,iabt,iuntst,
     >     imxfxt,ifxt,ilablr,iabr,iuntsr,itypr)
c
c  passed quantities:
c   icall= call # (2 passes required)
c   call first with icall=1, then with icall=2, to fill all output
c   irecsz= (output) plot file record size
c    all are from common block "CPLOTR.BLK"; names changed to prevent
c    prevent possible conflict with trans.com
c  x axis information:
c  imaxx-- max number
c
        integer izonex(imaxx),inumr(imaxx),itotf(imaxx)
        character*32 iklblx(imaxx)
        character*10 ikabrx(imaxx)
c
c  multigraph information
c  imaxbl-- max number of multigraph
c
        integer infb(imaxbl),ifunb(15,imaxbl),iintb(imaxbl)
        character*32 ilablm(imaxbl)
        character*10 iabbm(imaxbl)
c
c  scalar function information
c
        character*32 ilablt(imaxft)
        character*10 iabt(imaxft)
c
c  profile function information
c
        integer itypr(imxfxt)
        character*32 ilablr(imxfxt)
        character*10 iabr(imxfxt)
c
c  units
c
        character*16 iuntsb(*),iuntsr(*),iuntst(*)
c
c  common of originator program
c
        include 'cparm.m'
      include 'cbaldr.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'crpltrgn.m'
        include 'crplocal.m'
c
        common /cnvect/ vftot(55,9), flxtot(55,9), srctot(55,9)
c----
c  the following logical variable is always .true.
c
        logical iall
        data iall/.true./
c----
c
         imaxsw=t1isw  ! nset ARRAY DIMENSION, #OF LOGICAL OUTPUT GROUPS
c
C**> INSERTING GENERATED CODE; TEMPLATE FILE WAS:
C**>  BALPLOT:TFILET.FOR                                             
C >>> PLTRGN GENERATED FORTRAN STARTS HERE:
      IF(ICALL.NE.1) GO TO 8000
C>>> FOLLOWING CODE GENERATED BY SUBROUTINE XAXSR2 (PRPLOT):
C VAX BINARY "MF" FILE RECORD SIZE:
      IRECSZ=NZONES                        
C INCREMENT COUNTER AND SET NUMBER OF ZONES
      INUMX=   1
      IF(INUMX.GT.IMAXX) GO TO 9000
      IZONEX(INUMX)=NZONES                        
C NUMBER OF RECORDS FOR THIS NUMBER OF ZONES:
      INUMR(INUMX)=1 + (IZONEX(INUMX) - 1)/IRECSZ
C CLEAR NUMBER OF FCNS COUNT FOR THIS AXIS
      ITOTF(INUMX)=0
C LABEL FOR THIS X AXIS:
      IKLBLX(INUMX)='<RADIUS>  '
      IKABRX(INUMX)='XZONI     '
C INCREMENT COUNTER AND SET NUMBER OF ZONES
      INUMX=   2
      IF(INUMX.GT.IMAXX) GO TO 9000
      IZONEX(INUMX)=NZONES                        
C NUMBER OF RECORDS FOR THIS NUMBER OF ZONES:
      INUMR(INUMX)=1 + (IZONEX(INUMX) - 1)/IRECSZ
C CLEAR NUMBER OF FCNS COUNT FOR THIS AXIS
      ITOTF(INUMX)=0
C LABEL FOR THIS X AXIS:
      IKLBLX(INUMX)='<RADIUS>  '
      IKABRX(INUMX)='XBOUNI    '
C---------------------
C SET UP SWITCHING; FORTRAN GEN. BY SUBROUTINE PLSWCH
      CALL PSWITCH(IMAXSW,IER)
      IF(IER.NE.0) GO TO 9300
C---------------------
C--------------------
C  DEFINE MULTIGRAPH LABELS AND INITIALIZE MULTIGRAPH COUNTS
C  THIS SECTION GEN. BY PLTRGN SUBROUTINE PMULTI
C
      IMUL=   1
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='DENSITIES                       '
      IABBM(IMUL)='PDENS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   2
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='IMPURITY DENSITIES              '
      IABBM(IMUL)='IDENS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   3
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='IMPURITY MEAN Z                 '
      IABBM(IMUL)='MEANZ     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   4
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='TEMPERATURES                    '
      IABBM(IMUL)='PTEMP     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   5
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='CURRENT DENSITIES               '
      IABBM(IMUL)='PCURS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   6
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SAWTOOTH ELEC TEMP              '
      IABBM(IMUL)='STTE      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   7
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SAWTOOTH ION TEMP               '
      IABBM(IMUL)='STTI      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   8
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SAWTOOTH CURR DENS              '
      IABBM(IMUL)='STJZ      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=   9
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SOURCES AND SINKS               '
      IABBM(IMUL)='SORCS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  10
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='RADIATED POWER                  '
      IABBM(IMUL)='PRADS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  11
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='THERM NEUT DENSITIES            '
      IABBM(IMUL)='DENS0     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  12
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='THERM. NEUTRAL TEMPS            '
      IABBM(IMUL)='T0        '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  13
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='THERM DIFFUSIVITY               '
      IABBM(IMUL)='KAPA      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  14
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='THEORY EL THERM DIF             '
      IABBM(IMUL)='CHTHE     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  15
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='THEORY ION THERM DIF            '
      IABBM(IMUL)='CHTHI     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  16
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='DIMENSIONLESS THEORY            '
      IABBM(IMUL)='THDML     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  17
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='TOTAL PLASMA CURRENT            '
      IABBM(IMUL)='TCURS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  18
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='CENTRAL DENSITIES               '
      IABBM(IMUL)='CDENS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  19
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='CENTRAL TEMPERATURES            '
      IABBM(IMUL)='CTEMP     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  20
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='EDGE DENSITIES                  '
      IABBM(IMUL)='NEDGE     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  21
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='EDGE TEMPERATURES               '
      IABBM(IMUL)='TEDGE     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  22
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='EDGE IMP. DENSITIES             '
      IABBM(IMUL)='NIEDG     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  23
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='AVERAGE DENSITIES               '
      IABBM(IMUL)='NAVG      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  24
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='AVERAGE TEMPERATURES            '
      IABBM(IMUL)='TAVG      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  25
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='NEUTRON EMISSION                '
      IABBM(IMUL)='XNEUT     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  26
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='POLOIDAL BETAS                  '
      IABBM(IMUL)='MBPOL     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  27
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='TOROIDAL BETAS                  '
      IABBM(IMUL)='MBTOR     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  28
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='PLASMA DIMENSIONS               '
      IABBM(IMUL)='RPLAS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  29
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SAWTOOTH AMPL.                  '
      IABBM(IMUL)='STAMP     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  30
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='SAWTOOTH RADII                  '
      IABBM(IMUL)='STRAD     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  31
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='HEATING POWERS                  '
      IABBM(IMUL)='PHEAT     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  32
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='OHMIC HEATING POWERS            '
      IABBM(IMUL)='POHS      '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  33
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='D-T ALPHA POWERS                '
      IABBM(IMUL)='PALFS     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  34
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='BEAM POWER USED IN Q            '
      IABBM(IMUL)='PBEMQ     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  35
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='D-T FUSION Q                    '
      IABBM(IMUL)='DTQ       '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  36
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='ENERGY CONFINEMENT              '
      IABBM(IMUL)='TAUES     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C
      IMUL=  37
      IF(IMUL.GT.IMAXBL) GO TO 9400
      ILABLM(IMUL)='LOOP VOLTAGES                   '
      IABBM(IMUL)='VLOOP     '
      IINTB(IMUL)=0
      INFB(IMUL)=0
C-----------------------------
C   NOW DEFINE THE FUNCTIONS TO BE PLOTTED (GEN. BY PLTRGN)
 8000 CONTINUE
      IF(ICALL.EQ.1) GO TO 10000
C  INITIALIZE FUNCTION COUNTERS
      IFT=0
      IFXT=0
        CALL       TFIX01(ICALL,IRECSZ,
     >     IMAXX,INUMX,IZONEX,INUMR,ITOTF,IKLBLX,IKABRX,
     >     IMAXBL,IMUL,INFB,IFUNB,IINTB,ILABLM,IABBM,IUNTSB,
     >     IMAXFT,IFT,ILABLT,IABT,IUNTST,
     >     IMXFXT,IFXT,ILABLR,IABR,IUNTSR,ITYPR)
C
c
c----
10000   continue
        return
c----
c  errors
c  branched to from generated fortran
c
 9300   continue
        call mesage(
     >' fatal error in TFILSN - switching buf too small')
        call mesage(
     >' try increasing dim. T1ISW  in common /rplgen/. ')
        call mesage(
     >' regenerate all plot subroutines.               ')
        stop
 9000   continue
        call mesage(
     >' fatal error in TFILSN - too many independant X ')
        call mesage(
     >' try increasing parm NAXXVR in plot common and  ')
        call mesage(
     >' regenerate all plot subroutines.               ')
        stop
 9400   continue
        call mesage(
     >' fatal error in TFILSN - too many multigraphs   ')
        call mesage(
     >' try increasing parm NAXMGP in plot common and  ')
        call mesage(
     >' regenerate all plot subroutines.               ')
        stop
        end
c--------1---------2---------3---------4---------5---------6---------7-c
        SUBROUTINE TFIX01(ICALL,IRECSZ,
     >     IMAXX,INUMX,IZONEX,INUMR,ITOTF,IKLBLX,IKABRX,
     >     IMAXBL,IMUL,INFB,IFUNB,IINTB,ILABLM,IABBM,IUNTSB,
     >     IMAXFT,IFT,ILABLT,IABT,IUNTST,
     >     IMXFXT,IFXT,ILABLR,IABR,IUNTSR,ITYPR)
C
        INTEGER IZONEX(IMAXX),INUMR(IMAXX),ITOTF(IMAXX)
        CHARACTER*32 IKLBLX(IMAXX)
        CHARACTER*10 IKABRX(IMAXX)
        INTEGER INFB(IMAXBL),IFUNB(15,IMAXBL),IINTB(IMAXBL)
        CHARACTER*32 ILABLM(IMAXBL)
        CHARACTER*10 IABBM(IMAXBL)
        CHARACTER*32 ILABLT(IMAXFT)
        CHARACTER*10 IABT(IMAXFT)
        INTEGER ITYPR(IMXFXT)
        CHARACTER*32 ILABLR(IMXFXT)
        CHARACTER*10 IABR(IMXFXT)
        CHARACTER*16 IUNTSB(*),IUNTSR(*),IUNTST(*)
        include 'cparm.m'
      include 'cbaldr.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'crpltrgn.m'
        include 'crplocal.m'
        LOGICAL IALL
C
      IF(NSET(   2).EQ.1) THEN
        CALL       TF0001(ICALL,IRECSZ,
     >     IMAXX,INUMX,IZONEX,INUMR,ITOTF,IKLBLX,IKABRX,
     >     IMAXBL,IMUL,INFB,IFUNB,IINTB,ILABLM,IABBM,IUNTSB,
     >     IMAXFT,IFT,ILABLT,IABT,IUNTST,
     >     IMXFXT,IFXT,ILABLR,IABR,IUNTSR,ITYPR)
C
      ENDIF
C
       RETURN
       END
C------------------/PLTRGN:  GENERATED FORTRAN ROUTINE
        SUBROUTINE TF0001(ICALL,IRECSZ,
     >     IMAXX,INUMX,IZONEX,INUMR,ITOTF,IKLBLX,IKABRX,
     >     IMAXBL,IMUL,INFB,IFUNB,IINTB,ILABLM,IABBM,IUNTSB,
     >     IMAXFT,IFT,ILABLT,IABT,IUNTST,
     >     IMXFXT,IFXT,ILABLR,IABR,IUNTSR,ITYPR)
C
        INTEGER IZONEX(IMAXX),INUMR(IMAXX),ITOTF(IMAXX)
        CHARACTER*32 IKLBLX(IMAXX)
        CHARACTER*10 IKABRX(IMAXX)
        INTEGER INFB(IMAXBL),IFUNB(15,IMAXBL),IINTB(IMAXBL)
        CHARACTER*32 ILABLM(IMAXBL)
        CHARACTER*10 IABBM(IMAXBL)
        CHARACTER*32 ILABLT(IMAXFT)
        CHARACTER*10 IABT(IMAXFT)
        INTEGER ITYPR(IMXFXT)
        CHARACTER*32 ILABLR(IMXFXT)
        CHARACTER*10 IABR(IMXFXT)
        CHARACTER*16 IUNTSB(*),IUNTSR(*),IUNTST(*)
        include 'cparm.m'
      include 'cbaldr.m'
        include 'cfokkr.m'
        include 'cfreya.m'
        include 'crpltrgn.m'
        include 'crplocal.m'
        LOGICAL IALL
C
C
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RADIUS                          '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='XZONI     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RADIUS                          '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='XBOUNI    '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RMINOR CENTR                    '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='RZON      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RMINOR BNDRY                    '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='RBOUN     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RMAJOR CENTR                    '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='RMAJC     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RMAJOR BNDRY                    '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='RMAJB     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='SHAF SHIFT BNDRY                '
      IUNTSR(IFXT)='M               '
      IABR(IFXT)='SSZB      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='SHAF SHIFT CENTR                '
      IUNTSR(IFXT)='M               '
      IABR(IFXT)='SSZC      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ELONGATION BNDRY                '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='ELONG     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TRIANGULARITY BNDRY             '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='TRING     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='INDENTATION BNDRY               '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='DENT      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='VOLUME BNDRY                    '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='VOLB      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='GRHO                            '
      IUNTSR(IFXT)='1/M             '
      IABR(IFXT)='GRHO1     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='GRHO2                           '
      IUNTSR(IFXT)='1/M**2          '
      IABR(IFXT)='GRHO2     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ZONE VOLUME                     '
      IUNTSR(IFXT)='CM**3           '
      IABR(IFXT)='DVOL      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ZONE CROSS-SEC. AREA            '
      IUNTSR(IFXT)='CM**2           '
      IABR(IFXT)='DAREA     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='FLUX SURF AVG <DR>              '
      IUNTSR(IFXT)='CM              '
      IABR(IFXT)='DRAV      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='FLUX SURFACE AREA               '
      IUNTSR(IFXT)='CM**2           '
      IABR(IFXT)='SURF      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ELECTRON DENSITY                '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NE        '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL ION DENSITY               '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NI        '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ION DENSITY H=1                 '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NH1       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ION DENSITY H=2                 '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NH2       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='BEAM ION DENSITY                '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='BDENS     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='IMPURITY DENSITY I=1            '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NIMP1     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='IMPURITY DENSITY I=2            '
      IUNTSR(IFXT)='N/CM**3         '
      IABR(IFXT)='NIMP2     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ELECTRON TEMPERATURE            '
      IUNTSR(IFXT)='KEV             '
      IABR(IFXT)='TE        '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   4
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ION TEMPERATURE                 '
      IUNTSR(IFXT)='KEV             '
      IABR(IFXT)='TI        '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   4
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='PLASMA ZEFF PROFILE             '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='ZEFFB     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   3
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='Q PROFILE                       '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='Q         '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL PLASMA CURRENT            '
      IUNTSR(IFXT)='KAMPS/CM2       '
      IABR(IFXT)='CUR       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   5
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='BEAM DRIVEN  CURRENT            '
      IUNTSR(IFXT)='KAMPS/CM2       '
      IABR(IFXT)='CURB      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   5
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='LOOP VOLTAGE                    '
      IUNTSR(IFXT)='VOLTS           '
      IABR(IFXT)='V         '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL CHI-E                     '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='XETOT     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  14
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='THEORY CHI-E                    '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='XETHE     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  14
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='DRIFT WAVE CHI-E                '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHDRE     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  14
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ETA_I CHI-E                     '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHIGE     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  14
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RES BALLOONING CHI-E            '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHRBE     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  14
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL CHI-I                     '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='XITOT     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  15
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='THEORY CHI-I                    '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='XITHE     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  15
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='DRIFT WAVE CHI-I                '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHDRI     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  15
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ETA_I CHI-I                     '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHIGI     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  15
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='RES BALLOONING CHI-I            '
      IUNTSR(IFXT)='M**2/SEC        '
      IABR(IFXT)='CHRBI     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  15
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ETA-I                           '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='ETAI      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  16
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ETA-I-THRESHOLD                 '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='ETATH     '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  16
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='F-ITH                           '
      IUNTSR(IFXT)='                '
      IABR(IFXT)='FITH      '
      ITYPR(IFXT)=   2
      ITOTF(   2)=ITOTF(   2)+1
      JMUL=  16
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ALPHA HEATING POWER             '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PHALF     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='AUX. HEATING POWER              '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PAUX      '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='BEAM  HEATING POWER             '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PBTOT     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='BEAM TO ION POWER               '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PBI       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='BEAM TO ELECT POWER             '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PBE       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='OHMIC HEATING POWER             '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='POH       '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL HEATING POWER             '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='AHEAT     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='TOTAL RAD. POWER                '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PRADT     '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
      JMUL=  10
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
      JMUL=   9
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFXT
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='EL. IONIZ. LOSS                 '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PWALLE    '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  LABEL F(X,T)
      IFXT=IFXT+1
      ILABLR(IFXT)='ION IONIZ.+CX LOSS              '
      IUNTSR(IFXT)='WATTS/CM3       '
      IABR(IFXT)='PWALLI    '
      ITYPR(IFXT)=   1
      ITOTF(   1)=ITOTF(   1)+1
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA CURRENT                  '
      IUNTST(IFT)='KA              '
      IABT(IFT)='TPCUR     '
      JMUL=  17
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='BEAM DRIV. CURRENT              '
      IUNTST(IFT)='KA              '
      IABT(IFT)='TBCUR     '
      JMUL=  17
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='MAJOR RADIUS                    '
      IUNTST(IFT)='CM              '
      IABT(IFT)='RMAJR     '
      JMUL=  28
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='MINOR RADIUS                    '
      IUNTST(IFT)='CM              '
      IABT(IFT)='RMINR     '
      JMUL=  28
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA ELONGATION               '
      IUNTST(IFT)='                '
      IABT(IFT)='KAPPA     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA TRIANG                   '
      IUNTST(IFT)='                '
      IABT(IFT)='DELTA     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA INDENTATION              '
      IUNTST(IFT)='                '
      IABT(IFT)='INDNT     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA VOLUME                   '
      IUNTST(IFT)='CM**3           '
      IABT(IFT)='PVOL      '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='PLASMA AREA                     '
      IUNTST(IFT)='CM**2           '
      IABT(IFT)='PAREA     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='EXTERNAL BZ AT R=0              '
      IUNTST(IFT)='TESLA           '
      IABT(IFT)='BZ        '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON DENS. @R=0             '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NE0       '
      JMUL=  18
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ION DENSITY @R=0                '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NI0       '
      JMUL=  18
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON TEMP. @R=0             '
      IUNTST(IFT)='KEV             '
      IABT(IFT)='TE0       '
      JMUL=  19
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ION TEMPERATURE @R=0            '
      IUNTST(IFT)='KEV             '
      IABT(IFT)='TI0       '
      JMUL=  19
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ZEFF AT R=0                     '
      IUNTST(IFT)='                '
      IABT(IFT)='ZEFF0     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TIME STEP                       '
      IUNTST(IFT)='SECONDS         '
      IABT(IFT)='DT        '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='THERMAL D-D NEUTRONS            '
      IUNTST(IFT)='/SEC            '
      IABT(IFT)='NEUTX     '
      JMUL=  25
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='BEAM-TARGET D-D NEUT            '
      IUNTST(IFT)='/SEC            '
      IABT(IFT)='BTNTS     '
      JMUL=  25
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL D-D NEUTRONS              '
      IUNTST(IFT)='/SEC            '
      IABT(IFT)='NEUTT     '
      JMUL=  25
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON TEMP. @EDGE            '
      IUNTST(IFT)='EV              '
      IABT(IFT)='TEEDG     '
      JMUL=  21
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ION TEMP. @EDGE                 '
      IUNTST(IFT)='EV              '
      IABT(IFT)='TIEDG     '
      JMUL=  21
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON DENS. @EDGE            '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NEEDG     '
      JMUL=  20
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='G=1 EDGE DENSITY                '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NG1EG     '
      JMUL=  20
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='G=2 EDGE DENSITY                '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NG2EG     '
      JMUL=  20
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='I=1 EDGE DENSITY                '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NI1EG     '
      JMUL=  22
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='I=2 EDGE DENSITY                '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NI2EG     '
      JMUL=  22
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='DENS. WGHTD. AVG. TE            '
      IUNTST(IFT)='KEV             '
      IABT(IFT)='TEDAV     '
      JMUL=  24
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='DENS. WGHTD. AVG. TI            '
      IUNTST(IFT)='KEV             '
      IABT(IFT)='TIDAV     '
      JMUL=  24
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='VOL. AVG. EL. DENS.             '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NEVAV     '
      JMUL=  23
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='VOL. AVG. ION. DENS.            '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='NIVAV     '
      JMUL=  23
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='LINE AVG. EL. DENS.             '
      IUNTST(IFT)='N/CM**3         '
      IABT(IFT)='DENBL     '
      JMUL=  23
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='LI/2 + BETA(P) TOTAL            '
      IUNTST(IFT)='                '
      IABT(IFT)='LI2PB     '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='CALC. SURF. VOLTAGE             '
      IUNTST(IFT)='VOLTS           '
      IABT(IFT)='VSURC     '
      JMUL=  37
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='LOOP VOLTAGE @R=0               '
      IUNTST(IFT)='VOLTS           '
      IABT(IFT)='VCNTR     '
      JMUL=  37
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='AVG. LOOP VOLTAGE               '
      IUNTST(IFT)='VOLTS           '
      IABT(IFT)='VAVG      '
      JMUL=  37
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ENERGY CONFINEMENT              '
      IUNTST(IFT)='SECONDS         '
      IABT(IFT)='TAUEA     '
      JMUL=  36
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ENERGY CONFINEMENT              '
      IUNTST(IFT)='SECONDS         '
      IABT(IFT)='TAUEB     '
      JMUL=  36
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON BETA (TOR.)            '
      IUNTST(IFT)='                '
      IABT(IFT)='BETTE     '
      JMUL=  27
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ION BETA (TOR.)                 '
      IUNTST(IFT)='                '
      IABT(IFT)='BETTI     '
      JMUL=  27
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='BEAM BETA (TOR.)                '
      IUNTST(IFT)='                '
      IABT(IFT)='BETTB     '
      JMUL=  27
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ALPHA BETA (TOR.)               '
      IUNTST(IFT)='                '
      IABT(IFT)='BETTA     '
      JMUL=  27
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL BETA (TOR.)               '
      IUNTST(IFT)='                '
      IABT(IFT)='BETTT     '
      JMUL=  27
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ELECTRON BETA (POL.)            '
      IUNTST(IFT)='                '
      IABT(IFT)='BETAE     '
      JMUL=  26
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ION BETA (POL.)                 '
      IUNTST(IFT)='                '
      IABT(IFT)='BETAI     '
      JMUL=  26
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='BEAM BETA (POL.)                '
      IUNTST(IFT)='                '
      IABT(IFT)='BETAB     '
      JMUL=  26
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ALPHA BETA (POL.)               '
      IUNTST(IFT)='                '
      IABT(IFT)='BETAA     '
      JMUL=  26
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL BETA (POL.)               '
      IUNTST(IFT)='                '
      IABT(IFT)='BETAT     '
      JMUL=  26
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL STORED ENERGY             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='WTOT      '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='WE THERMAL                      '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='WTHE      '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='WI THERMAL                      '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='WTHI      '
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL OHMIC HEATING             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='POHT      '
      JMUL=  31
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
      JMUL=  32
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ST AVG OHMIC HEATING            '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='POHAV     '
      JMUL=  32
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ST RECONNEC. HEATING            '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PCRSH     '
      JMUL=  32
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL HEATING POWER             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PHETT     '
      JMUL=  31
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='RADIATED POWER                  '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PRAD      '
      JMUL=  31
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ALPHA HEATING POWER             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PALFT     '
      JMUL=  33
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
      JMUL=  31
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='EXTRAP. ALPHA POWER             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='APSTR     '
      JMUL=  33
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='THERMO. ALPHA POWER             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='APTH      '
      JMUL=  33
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='B-TARGET ALPHA PWR.             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='APBT      '
      JMUL=  33
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL BEAM HEATING              '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='BPHTO     '
      JMUL=  34
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
      JMUL=  31
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='EXTRAP. BEAM POWER              '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='BPSTR     '
      JMUL=  34
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TOTAL INJECTED POWER            '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PINJ      '
      JMUL=  34
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='ABSORBED BEAM POWER             '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='BPABS     '
      JMUL=  34
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='QTECH DENOMINATOR               '
      IUNTST(IFT)='WATTS           '
      IABT(IFT)='PTECH     '
      JMUL=  34
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='INSTANTANEOUS Q                 '
      IUNTST(IFT)='                '
      IABT(IFT)='QINST     '
      JMUL=  35
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='STEADY STATE Q                  '
      IUNTST(IFT)='                '
      IABT(IFT)='QSS       '
      JMUL=  35
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='Q-TECH                          '
      IUNTST(IFT)='                '
      IABT(IFT)='QTECH     '
      JMUL=  35
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='Q-COMPRESSION                   '
      IUNTST(IFT)='                '
      IABT(IFT)='QCOMP     '
      JMUL=  35
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='REL. SAWTOOTH TE AMP            '
      IUNTST(IFT)='                '
      IABT(IFT)='SAWTE     '
      JMUL=  29
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='REL. SAWTOOTH NE AMP            '
      IUNTST(IFT)='                '
      IABT(IFT)='SAWNE     '
      JMUL=  29
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='REL. SAWTOOTH DD AMP            '
      IUNTST(IFT)='                '
      IABT(IFT)='SAWDD     '
      JMUL=  29
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='Q=1 RADIUS AT CRASH             '
      IUNTST(IFT)='CM              '
      IABT(IFT)='SAWR1     '
      JMUL=  30
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='TE INVERSION RADIUS             '
      IUNTST(IFT)='CM              '
      IABT(IFT)='SAWRI     '
      JMUL=  30
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
C  SCALAR FUNCTION LABEL:
      IFT=IFT+1
      ILABLT(IFT)='SAWTOOTH MIXING RAD.            '
      IUNTST(IFT)='CM              '
      IABT(IFT)='SAWMX     '
      JMUL=  30
      IINTB(JMUL)=1
      INFB(JMUL)=INFB(JMUL)+1
      IFUNB(INFB(JMUL),JMUL)= 1*IFT 
      RETURN
      END
c@tfilwr
c rgb 16-feb-92 set idum = 0 before it is used
c rgb 02-oct-89 changed xlab, labelt, labelr, labelb to 32 char var
c    formats changed from 3a10 to a30
c
c=============================================================
c  modif 28 jul 81  d. mc cune
c   a second pair of routines "tfilw2" and "tfilr2" have been
c   added to support passing added information needed for 
c   multiple x-axes of varying lengths.
c  new info-- record position offset data for fcns of time + 1
c  addl coordinate.  labeling and size data for various x-axes,
c  and x-axis abbreviations which indicate which physical variable
c  in the mf file defines the x-axis physical values and units.
c    the old routines are maintained for backwards compatibility
c  to output of transp runs pre 9-1-81... d. mc cune
c-------
c  tfilwr(lun)   write transp plotting tf file (containing
c                labeling and dimensioning information)
c                on logical unit lun
c
        subroutine tfilwr(lun)
c
c  common blocks  (plotting system) ---
        include 'crplot.m'
c
      idum = 0
c
c  create and open file
c   call new version "tfilw2"
c    if flexible x-axis features are in use
        if(nlxvar) call tfilw2(lun)
        if(nlxvar) go to 500
c  write integer descriptors
        call tfilhd(lun,nrun,runid,nshot,nzones,nft,nfr,nxr,1)
c  write plasma torus major and minor radius
        write(lun,2002) rmajor,rminor
 2002   format(2e10.4)
c  write zone center and zone bdy radial arrays
        do 10 i=1,nzones
        write(lun,2002) rzon(i),rboun(i)
 10     continue
c  write labels, units and abbreviations of functions of time
        do 20 i=1,nft
        write(lun,2003) labelt(i)(1:20),unitst(i)(1:9),abt(i),idum
 2003   format(a20,a9,'$',a5,i1)
 20     continue
c  write labels, units and zone/bdy code, and abbreviations
c  of functions of time and radius
        do 30 i=1,nfr
        write(lun,2003) labelr(i)(1:20),unitsr(i)(1:9),abr(i),itypr(i)
 30     continue
        write(lun,2005) nbal
 2005   format(i3)
        if(nbal.eq.0) go to 100
        do 50 ib=1,nbal
        write(lun,2007) labelb(ib)(1:20),unitsb(ib)(1:9)
     &    ,iintb(ib),infb(ib),abb(ib)
 2007   format(a20,a9,'$',2i5,a5)
        write(lun,2008) (ifunb(j,ib),j=1,infb(ib))
 2008   format(15i4)
 50     continue
 100    continue
c  labels and values for independant variables other than time
        do 170 i=1,nxr
        write(lun,2009) xlab(i)
 170    continue
        do 180 i=1,nxr
        write(lun,2010) (xarry(j,i),j=1,nzones)
 180    continue
 2009   format(a32)
 2010   format(5(1pe10.3))
c
c
 500    continue
        return
                end
c--------1---------2---------3---------4---------5---------6---------7-c
c@tfilw2
c rgb 16-feb-92 set idum = 0 before it is used
c
c---------------------------------------------------------------
c  tfilw2
c
c  write transp plotting indexing/labeling file
c  this file contains all information needed to locate and plot
c  labeled the data in the corresponding "mf" and "nf" plotting
c files output by transp
c
        subroutine tfilw2(lun)
c
c  common blocks--
c
        include 'crplot.m'
c
c  file is already open on unit lun
c
        ixr=-nxr
c  write integer descriptor with # of x-axes negative to identif
c  this new file format
        call tfilhd(lun,nrun,runid,nshot,nzones,nft,nfr,ixr,1)
c  write actual # of fcns of time + 1 addl coordinate (nfr is now
c  the total # of mf file records containing such data which is
c  no longer necessarily the same as the # of fcns since fcns of
c  different length are now supported)
c
      idum = 0
c
        write(lun,2000) nfxt
 2000   format(i5)
c  write plasma torus minor and major radii
        write(lun,2002) rmajor,rminor
 2002   format(2e10.4)
c  write labels, units and abbreviations of functions of time
        do 20 i=1,nft
        write(lun,2003) labelt(i)(1:20),unitst(i)(1:9),abt(i),idum
 2003   format(a20,a9,'$',a5,i1)
 20     continue
c  write labels, units, abbreviations and random-access record
c  offsets for fcns of time + 1 addl coordinate
        do 30 i=1,nfxt
        write(lun,2004) labelr(i)(1:20),unitsr(i)(1:9),abr(i),itypr(i),
     >                  nrofff(i)
 2004   format(a20,a9,'$',a5,i2,i4)
 30     continue
c  write out multigraph associations
        write(lun,2005) nbal
 2005   format(i3)
        if(nbal.eq.0) go to 100
        do 50 ib=1,nbal
        write(lun,2007) labelb(ib)(1:20),unitsb(ib)(1:9)
     &   ,iintb(ib),infb(ib),abb(ib)
 2007   format(a20,a9,'$',2i5,a5)
        write(lun,2008) (ifunb(j,ib),j=1,infb(ib))
 2008   format(15i4)
 50     continue
 100    continue
c  write out x-axis specifications
c  (quantity nxr already written out, first line)
        do 110 i=1,nxr
        write(lun,2010) nfx(i),nroffx(i),nrecx(i),nzonex(i),
     >       xlab(i)(1:10),xndabb(i)
 110    continue
 2010   format(i4,1x,i4,1x,i3,1x,i4,1x,a10,1x,a5)
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@tfilhd
c
c-----------------------------------------------------------------
c  subroutine tfilhd
c
c  read/write header (first record) * seperate routine dmc nov 1984
c  kluge stuff record to allow 6 digit shot no., but also be able
c  to read old files!
c
c  modif dmc june 1985 - run id is 6 char string - replacing run no.
c  but continue to support them old files
c
        subroutine tfilhd(lun,nrun,
     >      runid,nshot,nzones,nft,nfr,nxr,iwrite)
c
        character*10 ztemp
        character*6 runid
c
        if(iwrite.eq.1) then
          nrun=-1
          write(lun,2001) nrun,nshot,nzones,nft,nfr,nxr,runid
 2001   format(i4,i6,4i5,a6)
c
c  old format
c2001   format(6i5)
c
        else
c
c  read has "hysteresis effect"
c   use fact that **all** transp run numbers have precisely 4 digits
c   (none have 3, none have 5 digit run numbers)
c   thus distinguish old and new formats by checking if 1st char is
c   blank
c
c... june 1985: more "hysteresis":  nrun=-1 means new files where runs
c  are identified by 6 char variables.
c
          read(lun,2002) ztemp,nzones,nft,nfr,nxr,runid
 2002   format(a10,4i5,a6)
c
c  new format
          if(ztemp(1:4).eq.'  -1') then
            read(ztemp,'(4x,i6)') nshot
          else if(ztemp(1:1).eq.' ') then
c  old format
            read(ztemp,'(i5,i5)') nrun,nshot
            write(runid,'(i4,2x)') nrun
          else
c  old new format
            read(ztemp,'(i4,i6)') nrun,nshot
            write(runid,'(i4,2x)') nrun
c
          endif
          nrun=-1
c
        endif
c
        return
        end

c@sawavg  .../11040/baldur/code/bald/sawavg.f
c  rgb 30-may-96 save zteold, ztecum, ztiold, zticum, ...
c    initialize these local arrays
c  rgb 28-jan-95 dimension pteavg(1) ... -> dimension pteavg(*) ...
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
c               %%%%%%%%%%%%%
c               %
c               awavg          %
c               %
c               %%%%%%%%%%%%%
c
c
        subroutine sawavg(k,pteavg,ptiavg,pjzavg,pohavg,pohtot,ptau)
c
c
      include 'cparm.m'
      include 'cbparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
        dimension zteold(55),ztecum(55),ztiold(55),zticum(55),
     1  zjzold(55),zjzcum(55),zohold(55),zohcum(55),
     1  pteavg(*),ptiavg(*),pjzavg(*),pohavg(*)
c
      logical inital
      data inital /.true./
      save inital
c
      save zteold, ztecum,ztiold,zticum,zjzold,zjzcum,zohold,zohcum
c
c..initialize local arrays
c
      if ( inital ) then
        do jz=1,mzones
          zteold(jz)=tes(1,jz)
          ztiold(jz)=tis(1,jz)
          zjzold(jz)=ajzs(1,jz)
          zohold(jz)=eta(2,jz)*ajzs(2,jz)*(ajzs(2,jz)-cjbeam*ajbs(jz))
          ztecum(jz) = 0.0
          zticum(jz) = 0.0
          zjzcum(jz) = 0.0
          zohcum(jz) = 0.0
        enddo
        inital = .false.
      endif
c
        do 100 jz=lcentr,mzones
        zoh=eta(2,jz)*ajzs(2,jz)*(ajzs(2,jz)-cjbeam*ajbs(jz))
        if(jz.eq.mzones) zoh=0.
        if(k.eq.1) then
c
c       increment the cumulative arrays
c
        ztecum(jz)=ztecum(jz)+uist*dtoldi*0.5*(zteold(jz)+tes(1,jz))
        zticum(jz)=zticum(jz)+uist*dtoldi*0.5*(ztiold(jz)+tis(1,jz))
        zjzcum(jz)=zjzcum(jz)+uist*dtoldi*0.5*(zjzold(jz)+ajzs(1,jz))
        zohcum(jz)=zohcum(jz)+uist*dtoldi*0.5*(zohold(jz)+zoh)
        end if
c
c       save the present profiles in ???old
c
        zteold(jz)=tes(1,jz)
        ztiold(jz)=tis(1,jz)
        zjzold(jz)=ajzs(1,jz)
        zohold(jz)=zoh
100     continue
c
        if(k.ne.2) return
c
c       compute the time averages and reset the cumulative arrays
c
        if(ptau.gt.0.) then
        ztau=ptau
        else
        ztau=1.
        end if
        zsum=0.
        do jz=lcentr,mzones
          pteavg(jz)=useh*ztecum(jz)/ztau
          ztecum(jz)=0.
          ptiavg(jz)=useh*zticum(jz)/ztau
          zticum(jz)=0.
          pjzavg(jz)=usej*zjzcum(jz)/ztau
          zjzcum(jz)=0.
          pohavg(jz)=zohcum(jz)/ztau
          zohcum(jz)=0.
          zsum=zsum+dvoli(jz)*pohavg(jz)
        enddo
        pohtot=zsum*usid*usep
c
        return
        end
c--------1---------2---------3---------4---------5---------6---------7-c

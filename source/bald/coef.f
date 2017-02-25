c--------1---------2---------3---------4---------5---------6---------7-c
        subroutine coef
c
c       2.2     compute coefficients
c
c
c@coef  .../baldur/code/bald/coef.f
c  rgb 04-mar-96 call fusdrv -> fusdrv(2)
c       les nov-90 add d3he fusion, analytic aux heating; nfusn=4
c       emg 25-apr-89 add the electron-ion interchange from sub. 'theory'
c                     to the ion energy source term
c       aes 14-jan-82 allow nfusn=3
c       aes 16-sep-81 add call to icrf; add/adjust expert calls 52,53
c       aes 4-nov-80 add call to ecrh; add/adjust expert calls 51,52
c       dhfz 5-apr-79 add call to he3
c       amck 14-nov-77 add alphas(2)
c       amck 30-oct-77 heat -> heat(2)
c       amck 11-aug-77 imprad(2) -> getchi(2)
c       amck 10-aug-77 change beams(1) -> beams(2)
c       amck 15-mar-77 add (2) to imprad call
c**********************************************************************c
c
c
c       amck 8-jul-76 add "(1)" to beams call
c       amck 20-may-76 move call to trcoef, delete "(k)", skip
c               neugas through heat if repeating timestep
c       amck 15-apr-76
c
       include 'cparm.m'
      include 'cbaldr.m'
c
c-----------------------------------------------------------------------
c
        data    iclass /2/,     isub /2/
c
        if(.not.nlomt2(isub))   go to 10
        call mesage(' *** 2.2 subroutine coef bypassed')
        return
   10   continue
c
c-----------------------------------------------------------------------
c
c
        call getchi(2)
c
c
        call trcoef
c
        call islbal
c
        if(nlrpet.or.nliter) go to 600
c
        call neugas
c
        call beams(2)
c
        call heat(2)
c
        call ecrh(2)
c
        call icrf(2)
c
c   les nov-90 d-3He fusion heating
c
      if ( cfutz(490) .gt. epslon ) then
c
        call fusdrv(2)
c
      else
c
        if(nfusn.ne.2.and.nfusn.ne.4) call alphas(2)
        if(nfusn.ne.1.and.nfusn.ne.4) call he3(2)
c
      endif
c
c   les nov-90 auxiliary heating - analytic profile
c
      if ( lpaux.gt.epslon .or. lspaux.gt.epslon ) call auxanl(2)
c
  600   continue
c
        call convrt
c
        call cnvcof
c
c
        return
        end

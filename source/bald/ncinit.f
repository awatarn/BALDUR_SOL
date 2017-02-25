c--------1---------2---------3---------4---------5---------6---------7-c
c@ncinit  /11040/baldur/code/bald/dimprad.f
c  rgb 29-dec-94 moved profile initialization to sbrtn ncprof0
c  rgb 29-dec-94 moved initialization and namelist from ncdata
c  rgb 19-dec-94 created sbrtn ncinit from part of sbrtn noncor
c***********************************************************************
c
        subroutine ncinit
c
c       initialization for the non-equilibrium impurity radiation
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'comncr.m'
      include 'comadp.m'
c
c-----------------------------------------------------------------------
c
      namelist /ncodat/
     x   recflx , recscr , beta   , daqexp , vaqexp , ifneo  , nwcool ,
     x   nwline , nwioni , nwreko , nwchex , nwbrem , cermlt , cizmlt ,
     x   cxrmlt , dnnmlt , dnmmlt , yhnn   , yhnm   , ndrmlt , ncxx   ,
     x   ncxxb  , naden  , nadtip , neci
c
c-----------------------------------------------------------------------
c
c
      if ( mimp .le. 0 )    return
c
c..initialize variables for the calculation
c
      ncrept = 2
      n2 = mzones
c
      call resetr (xns,kncimp*mxzone,0.0)
      call resetr (dq,kncimp*mxzone,0.0)
      call resetr (dz,kncimp*mxzone,0.0)
      call resetr (dweirs,mximp*mzones,0.0)
c
c
c    preset input variables
c
         recflx = 1.0
         recscr = 0.0
         beta = 1.0
         ifneo = 0
         nwcool = 0
         nwline = 0
         nwioni = 0
         nwreko = 0
         nwchex = 0
         nwbrem = 0
c
c  Preset ADPAK variables
c
      ncxxb = 3
c
      do 12 jii=1,kncimp
        naden(jii) = 1
        nadtip(jii) = 1
        neci(jii) = 2
        yhnn(1,jii) = 0.2
        yhnn(2,jii) = 0.7
        yhnn(3,jii) = 1.4
        yhnm(1,jii) = 0.1
        yhnm(2,jii) = 0.55
        yhnm(3,jii) = 1.0
        ndrmlt(jii) = 1
c
        do 11 ik=1,kzmax
          dnnmlt(ik,jii) = 1.0
          dnmmlt(ik,jii) = 1.0
   11   continue
c
        cermlt(jii) = 1.0
        cizmlt(jii) = 1.0
        cxrmlt(jii) = 1.0
c
c  Other species-dependent variables
c
        daqexp(jii) = 0.0
        vaqexp(jii) = 0.0
   12 continue
c
c  Dimension of ncxx determined by number of hydrogen species
c
      do 13 jii=1,mxhyd
        ncxx(jii) = 3
   13 continue
c
c    data for diffusion, boundary-condition and sources
c
      read (nread,ncodat)
c
      write (nprint,*)
      write (nprint,*)
      write (nprint,ncodat)
c
c
c..loop over impurity species
c
      do ji=1,mimp
c
        if (nimp(ji).gt.0) then
          nkimp(ji) = nimp(ji)
        else if (nkimp(ji).eq.-4) then
          nkimp(ji) = 2
        end if
        nk = nkimp(ji)
c
c    Set neutral influx parameters
c
        call ncsorc(1,ji)
c
c    set up tables of ionisation- and recombination-rates
c
        call ncrats(1,ji)
c
c
        xntot(ji) = 0.0
        pflx(ji) = 0.0
        ploss(ji) = 0.0
        psorc(ji) = 0.0
        flout(ji) = 0.0
        flscr(ji) = 0.0
c
        do jk=0,nk
          do j =1,mzones
            xn(j,jk,ji) = 0.0
          enddo
        enddo
c
        do jk=0,nk
          do jz =1,mzones
            sa(jz,jk,ji) = 0.0
            ra(jz,jk,ji) = 0.0
          enddo
        enddo
c
      enddo
c
      return
      end

c  11:00 15-jan-92 .../baldur/code/bald/dneugas.f
c/ 16.17 18-aug-89 /11040/bald89/wbaldn1 DNEUGAS, Bateman, Stotler, PPPL
c  BALDUR  file DNEUGAS by Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c
c  To obtain this file, type
c cfs get /11040/bald89/wbaldn1
c end
c lib wbaldn1 ^ x dneugas ^ end
c
c**********************************************************************c
c
c     this file contains the following packages:
c  neugas - neutral gas sink and source package
c
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************
c                           neugas                                    *
c**********************************************************************
c
c@neugas  .../baldur/code/bald/dneugas.f
c  rap 22-mar-02 namelist degas is not used -corresponding lines 
c                are commented out
c  rgb 19-sep-01 changed 0 to 0. in resetr(znei,mzones,0.)
c  rap 23-feb-00 call ifix(...) changed to call int(...)
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 07feb00 data newname/8hbalxxxr1,... --> data newname/'balxxxr1',
c  alh  6-nov-97 called timint to interpolate intermediate and end
c                values for tcold array 
c  rgb 14-aug-96 replaced fuzz with rndeps
c  rgb 02-jun-96 save ilinit, set maxgas = 2
cbate if(tai.le.0.0.or.rplsma.le.0.0) ilinit = .false. (reinstated)
c  rgb 14-jan-92 remove common/combas/ from sbrtn mcarlo and moved
c    common/comprt/ to file cmonprt.m, common/comxs/ to file cmonxs.m,
c    and common/comspl/ to file cmonspl.m
c     rgb 18-aug-89 annotations and comments by Heifetz
c......  /014325/bald92/wbaldn1 DNEUGAS FRANDM                          |
c  dps 11-jan-89 **.** Remove argument from RANF - not needed on Cray   |
c     rgb 15-jan-85 extensive changes from 11040 .heifetz neupdx
c        preparing to interface with FRANTIC, AURORA, and DEGAS.
c     rgb 17-nov-85   1 1/2 D upgrade, use common block commhd
c        plasma surface area = <|del xi|> V'(xi)
c        changed rmins*xzoni(j) to ahalfs(j,2)
c       drm 4-dec-84 comment out write statements
c           mhr 6-sep-84 isepm1 in first do loop replaced by isep
c       mhr 16-may-84 fixed bug in summing divertor-loss terms
c       aes 03-jun-82 eliminate cfutz(103): multiply ztotal contribution
c               to zz by grecyc at do 406; remove grecyc at do 384
c       aes 19-may-82 if cfutz(103)>0., use grecyc*(outflx-fluxr) rather
c               than outflx-fluxr, to compute fluxin at do 384
c       aes 30-apr-82 remove extra factor of gfrac1 in zpuffi at do 95
c       aes 22-apr-82 fix comment heading section 1.14
c       aes 21-apr-82 make. dens. monitor work right for diverted cases
c       aes 18-apr-82 corrected ztcold bug:  don't reset after t.s. 1
c               -- also fixed comment about cfutz(198)
c       aes 04-mar-82 set nlpion=nlgpin
c       aes 02-mar-82 only reset grecyc (via cfutz(89),(102)) once, not
c               on every call to neugas
c       aes 02-mar-82 make readjustment of tcold optional via cfutz(198)
c               (this coding from fgps' expert of 20-may-81)
c       aes 29-oct-81 array dimension 52 -->55 in array zpdx
c       aes 28-oct-81 insert coding from expert after 154
c       aes 28-oct-81 test for dtoldi=0. at 75,85 (as when tinit.ne.0.)
c       aes 20-aug-81 eliminate unnecess test for t.s. repeat at do 11
c       aes 13-aug-81 use do loop, not if, to set gfrac1 below line 40;
c               eliminate dtoldi from if statements near 75 and 85;
c               elminate obsolete variable ipuff near 40
c       aes 13-aug-81 set old quantities before resetting gfract, gvsrci
c               --> fix bug.  also, skip computation of ztau if, in
c               constant-monitor mode (denmon(1).le.epslon), monitor was
c               off and has just come back on
c       aes 10-aug-81 set wmin=gwmin near 205 (allow variable min. wt.)
c       aes 7-aug-81 gfract(iti,2) now set in auxval
c       aes 27-jul-81 eliminate monitor debug printout
c       aes 24-jul-81 zfract --> gfrac1; znel1 --> gnel1; zvsorc -->
c       gvsrci; zgfdni --> gfdni.  thus, all variables stored by the
c       density monitor between t.s. are now in common.  also, zvsor2
c       --> zvsrc2.
c       aes 8-jun-81 add trap for npuff .ne. 1 or 2
c       aes 21-may-81 gfract: - = target fract., + = influx fract.
c       aes 20-may-81 new density monitor. gpers no longer used.
c               gfract(jh) now means target fraction of hyd. jh
c               with respect to total hydrogen density
c       aes 1-dec-80 fix tcoldp (near line 224); removed near 403
c       aes added tcoldp business
c       fgps 21-mar-80 modified density monitoring to depend on elec-
c                      tron levels inside plasma edge or separatrix
c       fgps 14-sep-79 corrected truncation of scrapeoff losses when
c                      the crossfield outflux is strongly negative
c       fgps 19-apr-79 added divertor terms:  (1) expanded dimension
c                      statement; (2) involved subroutine pdx(--);
c                      (3) included divertor-loss terms in outflux
c       dhfz 15-mar-79 include gfract
c       dhfz 8-mar-79 renumbered 1-100;replaced znei by rhoels
c       dhfz 6-mar-79 zdenom=dtmax
c       dhfz 9-feb-79 added density monitoring
c       amck 8-feb-78 reset fluxin before profile
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 10-jan-78 compute new fluxin after profile.
c       amck 5-jan-78 if ngsplt.le.0, no splitting surfaces.
c       amck 4-jan-78 print out always
c       amck 4-jan-78 print out ztotal, zz, zflux, comusr, if ztotal > 1
c       amck 4-jan-78 include rmins*dx2i in sigvbx term in ztotal
c       amck 3-jan-78 include sigvbx in ztotal
c       amck 3-jan-78 include zvsorc in zflux (408 loop)
c       amck 2-jan-78 finish adding vol. source
c       amck 20-dec-77 start adding volume source
c       amck 4-nov-77 use shblos, gblosi
c       amck 1-nov-77 consider shchxs as also beam ion type neutral source
c       amck 1-nov-77 use shchxs as source to influx
c       amck 21-sep-77 set ztflwi, use grecyc
c       amck 22-aug-77 use fflxoi for influx, ignore addi
c       amck 2-aug-77 reset fluxn if no influx, incl. acompi in recycl.
c       amck 17-jul-77 compute dwicxs
c       amck 14-jul-77 don't compute if no influxes
c       amck 8-jul-77 set tns right if rhons=0.0
c       amkc 6-jul-77 average flux over timestep
c       amck 5-jul-77 change name of package to "monte"
c       amck 24-jun-77 fix weions, change sp? to units 1/secs.
c       amck 23-jun-77 fix wixxx, shions loop
c       amck 23-jun-77 reset zz
c       amck 22-jun-77 jh -> jg in 408 loop
c       amck 22-jun-77 fix limpn -> limp1
c       amck 20-jun-77 more comments, fix iz1
c       amck 12-apr-77 multi-species neutral gas
c       amck 11-apr-77 (it) -> (jt) in 104 loop
c       amck 17-mar-77 fix test to skip sputtering flux at 424
c       amck 8-mar-77 fix nzspec subscript around #424
c       amck 2-mar-77 mimp -> mhyd around #100
c       amck 13-feb-77 change i -> imax  around #254
c       amck 1-feb-77 equal spaced splitting surfaces--for mike
c       amck 28-jan-77 fix sign of wichxs
c       amck 28-jan-77 choose split. surf. starting from center
c       amck 27-jan-77 set ilinit after initialization
c       amck 27-jan-77 ledge-lcentr -> mzones-lcentr around #200,
c               zero rsplit before setting up new splitting surfaces
c       amck 27-jan-77 if gflx0i=gfluxi=0, use source/sinks = 0
c       amck 25-jan-77 use ngsprf instead of nlgprf
c       amck 25-jan-77 create file
c**********************************************************************c
c
c
        subroutine neugas
c
c
cl      2.3     compute neutral gas source and sink terms
c
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'cmonte.m'
      include 'commhd.m'
c
        common/comdeg/tneutral,lbindex(mj),lneutral,lswitch,nohzs,
     &   sorcns(2,3),sorcnos(2,3),dwiions(2,25,6),dweions(25,6),
     &   rperps(2),puffings(2),
     &   dwichxs(2,25,6),recycs(2,25),recomns(2,25),recxns(2,25)
c
        dimension msg(10),msgtocte(5)
        dimension zrparls(2,mj),zrperpi(2),zrecomi(2,mj),
     &   zbchxs(2,mj),pflux(250),kplshu(mj),
     &   zrparli(2,mj),zbchxi(2,mj),zrecoms(2,mj)
c
        logical
     l  il0,            ilinit,         ilprof,         ilnega
c
        dimension
     r  znh(4),         zflux(8),       zz(8,8),        ztotal(4,8),
     r  zpuffi(4),      zsne(4,2),      zsni(4,2),      zspe(2),
     r  zspi(2),        zt0(4,2),       zn0(4,2),       zvdist(4),
     r                  znei(mj),       zflpdxi(2),
     r  zpdx(6,55),     znel2(4),
     r  zgfdi2(4),      zfrac2(4),      zgflx2(4),      zvsrc2(4),
     r  ztau(4),
     i  ipivot(8),       newname(6)
c
c----------------------------------------------------------------------
c
        namelist/nbalneut/lbindex,lneutral,nohzs,kplpuf,kplrec,kplshu,
     &   tneutral
cap        namelist/degas/sni,sne,den0,eneut,fluxr,fluxin,outflx,
cap     &   dwiions,dweions,dwichxs
cap        data newname(1) /'balxxxr1'/,
cap     &       newname(2) /'balxxxr2'/,
cap     &       newname(3) /'balxxxr3'/,
cap     &       newname(4) /'balxxxr4'/,
cap     &       newname(5) /'balxxxr5'/,
cap     &       newname(6) /'balxxxr6'/
c
c
        data    ilinit /.false./
c
      save ilinit
c
c------------------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /3/
c
c
        if (nlomt2(isub)) then
c
         if (nstep .lt. 2)
     &  call mesage(' *** 2.3 subroutine neugas bypassed ')
        return
c
         else
         if (nstep .lt. 2) then
         write (nout,9000)
         write (ntychl,9000)
 9000  format (' sbrtn nuegas modified by Bateman, 20 Jan 86')
         endif
c
      endif
c
                        call expert(iclass,isub,0)
c
c
c------------------------------------------------------------------------------
c
c
cl      common blocks and variables modified:
c
c       comusr, comneu, shions, weions, wiions, wichxs
c
c
c------------------------------------------------------------------------------
c
c
c       this subroutine calls "monte", a monte-carlo neutral gas
c       package, whenever necessary to get a gas profile.
c       at other times, it scales the profiles of neutral density and
c       source and sink terms so that the total ionization
c       equals the desired net influx.
c
c
c------------------------------------------------------------------------------
c
        denmont = 0.
        maxgas = 2
c
        if (mhyd.gt.4) go to 9010
c
c
cl      1)      compute influx, decide whether to recompute profile
c
c
c       1.10) store quantities from last timestep for density monitor
c
        if(npuff.lt.1 .or. npuff.gt.2) go to 12
        if(gftime(1).ge.1.e10) go to 12
c                                     gas puffing initialiazation
        do 11 jh = 1,mhyd
        zfrac2(jh) = gfrac1(jh)  ! fraction of hydrogen
        znel2(jh)  = gnel1(jh)   ! <ni(ih)>
        zgflx2(jh) = gfluxi(jh)  ! total influx of hydrogen ih
        zvsrc2(jh) = gvsrci(jh)  ! volume source in neugas (recomb+beamcx)
        zgfdi2(jh) = gfdni(jh)   ! density monitor flux of hydrogen
  11    continue
  12    continue
c
        call resetr(zrparls,2*55,0.0)                          ! heifetz
        call resetr(zrparli,2*55,0.0)                          ! heifetz
        call resetr(zrperpi,2,0.0)                             ! heifetz
        call resetr(rperps,2,0.0)                              ! heifetz
        call resetr(zrecomi,2*55,0.0)                          ! heifetz
        call resetr(zrecoms,2*55,0.0)                          ! heifetz
        call resetr(zbchxs,2*55,0.0)                           ! heifetz
        call resetr(zbchxi,2*55,0.0)                           ! heifetz
        call resetr(sorcns,2*3,0.0)                            ! heifetz
c                                           compute parallel ion losses
        call resetr(zflpdxi,2,0.0)                             ! heifetz
        call resetr(zpdx,6*55,0.0)                             ! heifetz
c                                               sum parallel loss terms
c                                                        <sorncs(*,1)>
c
c       1.11) recycling source
c
c isep  = index of innermost scrape-off zone
c
        isep=mzones  ! when there is no scrape-off model
          if(nadump(1).gt.lcentr) isep = nadump(1)  ! with scrape-off
        isepm1=isep-1
c
c  zsurfi = plasma surface area in internal units = <|del xi|> V'(xi)
c  zsurfs = plasma surface area in standard units (cm**2)
c  zradi  = 2 * plasma volume / plasma surface area  in internal units
c           note zradi = rmini for plasma with circular cross section
c  zrads  = 2 * plasma volume / plasma surface area  in standard units (cm)
c
c
      zsurfi = avi(mzones,5,1) * avi(mzones,3,1)
      zsurfs = zsurfi * uisl**2
      zradi  = 2. * avi(mzones,12,1) / zsurfi
      zrads  = zradi * uisl
c
      z1 = 1.0 / zsurfi
      z2 = zradi / usit
c
      zu1 = usid * zradi / usit
      zu2 = 1. / ( uist * uisl**2 )
c
c
c  sum divertor-loss terms when edge model is used
c
        if(nadump(1).le.lcentr) go to 15
c
        call pdx   ! to recompute scroff(jh,jz)  Bateman
c
        do 14 jh=1,mhyd
        do 14 jz=isep,ledge
c parallel loss rate d n / d t = n / tau_parallel
            zrparli(jh,jz)=-chi(jh,jz)*scroff(jh,jz)*dx2i(jz)*z2 ! heifetz
            zrparls(jh,jz)=zrparli(jh,jz) * zu2                ! heifetz
            sorcns(jh,1)=sorcns(jh,1)+zrparls(jh,jz)           ! heifetz
            zflpdxi(jh)=zrparli(jh,jz)+zflpdxi(jh)             ! heifetz
   14   continue
c
c  get recycling tally, perpendicular recycling
c      all in units of length^{-2} time^{-1}
c
  15    continue
c
        do 16 jh = 1, mhyd
          zrperpi(jh)=z1*max(-fflxoi(jh),0.0)                   ! heifetz
c                   ! influx rate at outer boundary / plasma surface area
c                   ! fflxoi is from SOLVE
          rperps(jh)=zrperpi(jh) * zsurfs * zsurfs                ! heifetz
          sorcns(jh,1)=sorcns(jh,1)+rperps(jh) /zsurfs            ! heifetz
c                                                              <gfluxi>
        gfluxi(jh)=grecyc*(zflpdxi(jh) + z1*max(-fflxoi(jh),0.))
c                  parallel            + perpendicular source
  16    continue
c
c       1.12) volume source -- integrate shchxs, shblos
c
        call resetr(gvsrci,maxgas,0.0)
c
        do 22 jz = lcentr, ledge
        do 20 jh = 1, mhyd
c recombination
        zrecomi(jh,jz)=recoms(jh,jz)*rhohs(jh,2,jz)*dx2i(jz)*zu1  ! heifetz
            zrecoms(jh,jz)=zrecomi(jh,jz) * zu2                   ! heifetz
            sorcns(jh,3)=sorcns(jh,3)+zrecoms(jh,jz)              ! heifetz
c charge exchange
            zbchxi(jh,jz)=shchxs(jh,jz)*dx2i(jz) * zu1            ! heifetz
            zbchxs(jh,jz)=zbchxi(jh,jz) * zu2                     ! heifetz
            sorcns(jh,3)=sorcns(jh,3)+zbchxs(jh,jz)               ! heifetz
c total volume source
        gvsrci(jh) = gvsrci(jh) + (recoms(jh,jz)*rhohs(jh,2,jz) +
     1          shchxs(jh,jz)) * dx2i(jz) * zu1
   20   continue
   22   continue
c
        do 23 jh = 1, mhyd
        gfluxi(jh) = gfluxi(jh) + gblosi(jh)
c          recycling rate  + source due to beam charge exchange
   23   continue
c--------1---------2---------3---------4---------5---------6---------7-c
c
c     At this point, the neutral particle sources have been set up
c  resulting from beams (FREYA) and recycling (PDX)
c  It remains to set up the sources from gas puffing.
c
c
c       1.13) gas puffing sources (zpuffi)
c
        call resetr(zpuffi,mhyd,0.)                               ! heifetz
        call resetr(puffings,mhyd,0.)                             ! heifetz
c
c  if there is no density monitoring, go to 100
c
        if(npuff.lt.1 .or. npuff.gt.2) go to 100
        if(gftime(1).ge.1.e10) go to 100
c
c       1.13.1) density monitor
c
c
c               initialize
c
        isep = mzones
        if(nadump(1).gt.lcentr) isep = nadump(1)
        isepm1 = isep - 1
        znel0 = 0.
        call resetr(znei,mzones,0.)
        call resetr(gflmon,mhyd,0.)
        call resetr(ztau,mhyd,0.)
        call resetr(gfrac1,mhyd,0.)
        call resetr(gfdni,mhyd,0.)
c
c
c               find timezone
c
        iti=1
        zt=tbi*uist
c
 26     continue
c
        if(zt.lt.gftime(1)) go to 100
        if(zt.ge.gftime(iti).and.zt.le.gftime(iti+1))
     1          go to 29
        if(iti.eq.19) go to 27
        iti = iti+1
        go to 26
c
 27     continue
c
        iti=20
c
 29     continue
c
        if(denmon(1).le.epslon) go to 40
c
c               fitted monitoring
c
        if(iti.ne.1.or.denmon(1).gt.2.) go to 32
        go to (50,60),npuff
c
 30     continue
c
        denmon(1)=rlined(1)
        go to 32
c
 31     continue
c
        denmon(1)=rnebar(1)
c
 32     continue
        zint=(zt-gftime(iti))/(gftime(iti+1)-gftime(iti))
        znel0=denmon(iti)+
     1          zint*(denmon(iti+1)-denmon(iti))
        denmont = znel0
c
c
        do 33 jh=1,mhyd
        gfrac1(jh)=abs(gfract(iti,jh))+
     1          zint*(abs(gfract(iti+1,jh))-abs(gfract(iti,jh)))
 33     continue
c
        go to (55,65),npuff
c
c               constant monitoring
c
 40     continue
c
        if(2*(iti/2).eq.iti)go to 100
c
        do 42 jh = 1,mhyd
        gfrac1(jh) = abs(gfract(iti,jh))
 42     continue
c
        go to (50,60),npuff
c
c               compute rlined = line avg density at
c               start of timezone iti
c  15.07 generalize to 1 1/2 D
c
 50     continue
c
        znel0=rlined(iti)
        if(rlined(iti).gt.epslon) go to 55
c
        z0 = 0.
        do 51 jz = lcentr,ledge
c
        if (versno.gt.15.06) then
          zdrs = ahalfs(jz+1,1) - ahalfs(jz,1)
        else
          zdrs = dxzoni(jz)
        end if
c
        z0 = z0+rhoels(2,jz)*zdrs
 51     continue
c
        if (versno.gt.15.06) then
          rlined(iti) = z0 / ahalfs(isep,1)
        else
          rlined(iti) = z0 / xbouni(isep)
        end if
c
        if(denmon(1).gt.epslon) go to 30
        znel0=rlined(iti)
c
c               compute present line avg dens
c  15.07 generalize to 1 1/2 D
c
 55     continue
c
        z0 = 0.
        do 56 jz = lcentr,ledge
c
        if (versno.gt.15.06) then
          zdrs = ahalfs(jz+1,1) - ahalfs(jz,1)
        else
          zdrs = dxzoni(jz)
        end if
c
        z0 = z0 + rhoels(2,jz)*zdrs
 56     continue
c
        if (versno.gt.15.06) then
          z0 = z0 / ahalfs(isep,1)
        else
          z0 = z0 / xbouni(isep)
        end if
c
        go to 70
c
c               compute rnebar = vol avg dens to the separatrix at
c               start of timezone iti
c  15.07 generalize to 1 1/2 D
c
 60     continue
c
        znel0=rnebar(iti)
        if(rnebar(iti).gt.epslon) go to 65
c
        z0 = 0.
c
        do 61 jz = lcentr,isepm1
        z0 = z0 + rhoels(2,jz)*dx2i(jz)
 61     continue
c
        if (versno.gt.15.06) then
          z0 = 2. * z0 * avi(mzones,12,1) / avi(isep,12,1)
        else
          z0 = 2. * z0 / xbouni(isep)**2
        end if
c
        rnebar(iti) = z0
c
        if(denmon(1).gt.epslon) go to 31
        znel0=rnebar(iti)
c
c               compute present vol avg density to the separatrix
c  15.07 generalize to 1 1/2 D
c
 65     continue
c
        z0 = 0.0
        do 66 jz=lcentr,isepm1
        z0 = z0 + rhoels(2,jz)*dx2i(jz)
 66     continue
c
        if (versno.gt.15.06) then
          z0 = 2. * z0 * avi(mzones,12,1) / avi(isep,12,1)
        else
          z0 = 2. * z0 / xbouni(isep)**2
        end if
c
c               compute present vol avg density *to the wall*
c               for each hyd. species and summed over all hyds.
c
 70     continue
c
        call resetr(gnel1,mhyd,0.0)
c
        z1 = 0.0
        do 72 jh = 1,mhyd
          do 71 jz = lcentr,ledge
            gnel1(jh) = gnel1(jh) + rhohs(jh,2,jz)*dx2i(jz)
 71       continue
          gnel1(jh) = gnel1(jh) * 2.
          z1 = z1 + gnel1(jh)
 72     continue
c
c               get "equivalent target density", scaled to allow
c               for impurities and for the fact that line and vol.
c               avgs. are defined in terms of the separatrix,
c               whereas gas is puffed in from the wall.  the
c               monitor must use quantities computed to the wall
c               in figuring how much gas to puff.
c               we assume that the shape of the density profile
c               and the fraction of impurities both change little
c               over one timestep.
c
        znel0 = znel0 * z1 / z0
c
c               convert densities to internal units
c
        znel0 = znel0 * usid
        z1 = z1 * usid
        call scaler(gnel1,mhyd,usid)
c
c               get (plasma volume/plasma area)
c
      z2 = avi(mzones,12,1) / ( avi(mzones,5,1) * avi(mzones,3,1) )
c
c               if gfract.lt.0., it is the target fraction for
c               each hyd. species; if .gt.0, it is the fraction
c               of the total influx comprised by each species
c
        if(sign(1.,gfract(iti,1)).gt.0.) go to 85
c
c               target fractions
c
c
c               get ztau = (1./confinement time) for each species
c               use  old values for densities, sources, dt in
c               this computation
c               if this is the first time density monitor is being
c               called, there are no old values to use, so set
c               ztau = 0.0.  likewise set ztau = 0.0 if, in constant
c               monitor mode (denmon(1).le.epslon), monitor had been off
c               and has just come back on
c
        if(tai-epslon.le.gftime(1)*usit) go to 80
        if(denmon(1).le.epslon .and.
     1                  tai-epslon.le.gftime(iti)*usit) go to 80
        if(dtoldi.le.epslon) go to 80
        do 75 jh = 1,mhyd
          ztau(jh) = ((znel2(jh)-gnel1(jh))*z2/dtoldi
     1                  +zgflx2(jh) + zvsrc2(jh) + zgfdi2(jh))
     2                  /(gnel1(jh)*z2)
  75    continue
c
c               compute density monitoring influx rate
c
  80    continue
        do 82 jh = 1,mhyd
          gfdni(jh) = (znel0*gfrac1(jh)-gnel1(jh)) * z2/dti
     1                  + znel0*gfrac1(jh)*z2*ztau(jh)
     2                  - gfluxi(jh) - gvsrci(jh)
          gfdni(jh) = max(gfdni(jh),0.0)
          gfdni(jh) = min(gfdni(jh),gflmax(iti)*uisl**2*uist)
  82    continue
        go to 92
c
c               influx fractions
c
c               get ztau for all hydrogen
c
  85    continue
        if(tai-epslon.le.gftime(1)*usit) go to 88
        if(denmon(1).le.epslon .and.
     1                  tai-epslon.le.gftime(iti)*usit) go to 88
        if(dtoldi.le.epslon) go to 88
        do 87 jh = 1,mhyd
          ztau(1) = ztau(1) + (znel2(jh)-gnel1(jh))*z2/dtoldi
     1                  + zgflx2(jh) + zvsrc2(jh) + zgfdi2(jh)
  87    continue
        ztau(1) = ztau(1)/(z1*z2)
c
c               compute total influx rate -- multiply by
c               appropriate fraction to get influx for each species
c
  88    continue
        gfdni(1) = (znel0-z1)*z2/dti
     1                  + znel0*z2*ztau(1)
        do 89 jh = 1,mhyd
          gfdni(1) = gfdni(1) - gfluxi(jh) - gvsrci(jh)
  89    continue
        gfdni(1) = max(gfdni(1),0.0)
        gfdni(1) = min(gfdni(1),gflmax(iti)*uisl**2*uist)
        z0 = gfdni(1)
c
        do 90 jh = 1,mhyd
          gfdni(jh) = z0 * gfrac1(jh)
  90    continue
c
c               convert to external units
c
  92    continue
        do 94 jh=1,mhyd
          gflmon(jh) = gfdni(jh)*ueil**2*ueit
 94     continue
c
c               add influx*time to zpuffi
c               note zpuffi is scaled by 1./dti at line 120
c
        do 95 jh = 1,mhyd
          zpuffi(jh) = gfdni(jh)*dti
 95     continue
c
c       1.13.2) fixed gas flows
 100    continue
c
c
        zt1 = tai
c
        do 118 jt = 2, mxt
          zt = gtflwi(jt)
          if (zt.gt.0.0) go to 104
c
          do 102 jh = 1, mhyd
            zpuffi(jh) = zpuffi(jh) + max( 0.0,
     1        gflowi(jh,jt-1)*(tbi - max( tai, gtflwi(jt-1) )) )
  102     continue
c
          go to 120
c
  104   continue
        if (zt.lt.tai) go to 118
        zt2 = min(zt,tbi)
c
        zint1 = 1.0
        zint2 = 1.0
        z0 = gtflwi(jt) - gtflwi(jt-1)
        if (z0.le.epslon) go to 106
        zint1 = (zt1 - gtflwi(jt-1)) / z0
        zint2 = (zt2 - gtflwi(jt-1)) / z0
  106   continue
c
        do 108 jh = 1, mhyd
          zpuffi(jh) = zpuffi(jh) + (zt2 - zt1)*0.5*((zint1 + zint2)*
     1          gflowi(jh,jt) + (2.0 - zint1 - zint2)*gflowi(jh,jt-1))
  108   continue
c
        if (zt.ge.tbi) go to 120
        zt1 = zt2
  118   continue
c
  120   continue
c
c       1.13.3) at this point, zpuffi is in (particles/area) --
c               scale it by (/time)
c
        z0 = 1.0 / dti
        call scaler(zpuffi,mhyd,z0)
c
c
c       1.14) don't bother calculating negligible sources
c
c plasma surface area = <|del xi|> V'(xi) = avi(mzones,5,1)*avi(mzones,3,1)
c
c
        ztflwi = 0.0
        z0 = (avi(mzones,5,1) * avi(mzones,3,1) * dti)
        do 134 jh = 1, mhyd
          addi(jh) = addi(jh) + zpuffi(jh) * z0 
c                  ! particles added this step
          ztflwi = ztflwi + zpuffi(jh) + gfluxi(jh) - gblosi(jh)
c                  ! recycling only
  134   continue
c
c
        ztflwi = ztflwi * rndeps
c
        do 138 jh = 1, mhyd
          if (zpuffi(jh).le.ztflwi) zpuffi(jh) = 0.0
          if (gfluxi(jh).le.ztflwi) gfluxi(jh) = 0.0
          if (gvsrci(jh).le.ztflwi) gvsrci(jh) = 0.0
  138   continue
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      1.2)    decide whether to recompute neutral gas profile
c
c
  150   continue
c
c  have we initialized the neutrals yet?
c
        if(tai.le.0.0.or.rplsma.le.0.0) ilinit = .false.
c
        ilprof = .not. ilinit
        if (ngsprf+ngprof.le.nstep) ilprof = .true. ! do a new profile
c
c  last neutral timestep
c
        il0 = .false.
        do 154 jh = 1, mhyd
          if (gfluxi(jh).gt.epslon) il0 = .true.
          if (zpuffi(jh).gt.epslon) il0 = .true.
          if (gvsrci(jh).gt.epslon) il0 = .true.
        if (zpuffi(jh).gt.epslon.and.nsrci(2,jh).le.0) ilprof = .true.
        if (gvsrci(jh).gt.epslon.and.nsrci(3,jh).le.0) ilprof = .true.
  154   continue
c
c  compression -- not used in the 1-1/2-D BALDUR code
c
        if (rplsma.le.ahalfs(ledge,2).or.nzcomp.ne.0.or.
     1         rplsma.ge.ahalfs(ledge+1,2)) ilprof = .true.
        ilprof = ilprof.and.il0
c
c
c       if cfutz(198).gt.0, tcold is continually made the minimum
c       of (1) tis(2,*) at a radial point determined by (cfutz(198)*
c       zone index of outermost non-scrapeoff zone) and (2) the
c       originally prescribed value of tcold.  there is also the mid-run
c       capability of changing the recycling coefficient to be
c       cfutz(102), where cfutz(89) prescribes the re-set time (secs).
c
c  kluge to change recycling rate or temperature of gas puffed in
c
        zt = tai*uist
        call timint (zt,ztcold,bdtime,20,tcold,1,1)
        if(cfutz(198).le.epslon) go to 157
          nzzmax=mzones
          if(nadump(1).gt.lcentr) nzzmax=nadump(1)
          jc=int(cfutz(198)*float(nzzmax-lcentr))+lcentr
          tc=tis(2,jc)*evsinv*1.e-03
          ztcold=min(ztcold,tc)
  157   continue
c
        if(cfutz(89).le.epslon) go to 159
          zt=tai*uist
          if(zt.le.cfutz(89)) go to 159
             grecyc=cfutz(102)
             cfutz(89) = 0.
  159   continue
c
                                        call expert(iclass,isub,1)
        if (ilprof) go to 200
        if (il0) go to 350     ! diferent from heifetz version neupdx
c
c
cl              if no influx, don't compute
c
c
        i1 = mzones * mxhyd
        call resetr(wichxs,mxzone,0.0)
        call resetr(dwicxs,mxzone,0.0)
        call resetr(weions,mxzone,0.0)
        call resetr(wiions,mxzone,0.0)
        call resetr(shions,i1    ,0.0)
        call resetr(tns   ,i1    ,0.0)
        call resetr(rhons ,i1    ,0.0)
        call resetr(fluxn ,maxsrc,0.0)
c
        return
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
cl      2)      set up comusr for "monte"  (AURORA package)
c
c
c..seting variables:
c  wmin     - minimum monte carlo weight
c  nsurf - number of surfaces used in sbrtn monte
c  ngases          - number of hydrogen species
c  rsurf(j)        - halfwidth of surfaces used in sbrtn monte
c  recycs
c  recomns
c  recxns
c  tein(j)         - eletron temperature
c  tiin(j)         - ion temperature (eV)
c  denein(j)       - electron density (cm-3)
c  denion(jh,j)- ion density
c  sigvbx(j)       - total n*sigma-v for beam charge-exchange
c  snvol(jh,j)     - relative volume source rate of neutrals
c  rplasma         - half-width of plasma (cm)
c  nsrci(i,jh)     - source indicies for each gas
c  nsrces          - number of influxes
c  ngasi(i)        - gas type index for each influx
c  e0in(i)         - cold neutral temperature
c  nlmono(i)       - .true. if monoenergetic source required
c  fluxin(jn)      - influx of neutrals at boundary
c  fract(jg)       - desired proportions of number of test particles
c
c  set GWMIN = \frac{\mbox{expected neutral density at center of plasma}
c    }{\mbox{expected edge neutral density}}
c  although the default is 1.e-6, try 1.e-3 to run twice as fast since 
c  the monte carlo time is linear in \log{gwmin}.
c
  200   continue
c
        if (ilinit) go to 205
          call monte(1)
c                                      read auxiliary namelist nbalneut
cheifetz        call freeus(lio)
cheifetz        call open(lio,"balneut",2,len)
cheifetz        if(len.lt.0)go to 202
cheifetz        rewind lio
cheifetz        read(lio,nbalneut)
cheifetz        call close(lio)
c
cheifetz 202    lswitch=0
cheifetz        if(lneutral.eq.0)lswitch=1
c                                                                <wmin>
          wmin = gwmin
  205   continue
c
        nsurf = min0(min0( maxsur , ngzone+1 ), mzones-lcentr+1)
c
        zfac = (float(mzones - lcentr) + 0.125) / float(nsurf - 1)
c
c
        if (zfac.le.1.0) go to 9020
c
        if (mhyd.gt.maxgas) go to 9021
        ngases = mhyd
c
        iz1 = lcentr - 1
c
c               monte boundaries are (rather arbitrarily) placed
c               every so many baldur boundaries.  the exact choice
c               of baldur boundary is not important (except for
c               accuracy of results), so long as every monte
c               surface is at the same place as some baldur boundary.
c               this makes the averaging of te, ti, etc., easier
c
c
        rsurf(1)=0.0
        call resetr(recycs,2*25,0.0)
        call resetr(recomns,2*25,0.0)
        call resetr(recxns,2*25,0.0)
c
        do 218 js = 2, nsurf
cheifetz          if(lbindex(2).le.0)ib=int(float(js-1)*zfac)+lcentr
cheifetz          if(lbindex(2).gt.0)ib=lbindex(js-1)
          ib = int(float(js-1)*zfac) + lcentr
          iz0 = iz1 + 1
          iz1 = ib - 1
c
          rsurf(js) = ahalfs(ib,1)
c
c               densities are volume-averaged, temperatures
c               particle-averaged.  iz0 is the innermost baldur zone
c               and boundary in this monte cell (cell js-1), iz1
c               the outermost zone, and ib the outermost boundary.
c
          zne = 0.0
          call resetr(znh,ngases,0.0)
          call resetr(zvdist,ngases,0.0)
          zni = 0.0
          zee = 0.0
          zei = 0.0
          zsigvb = 0.0
c
        do 208 jz = iz0, iz1
          zne = zne + rhoels(2,jz) * dx2i(jz)
          znn = 0.0
          zbloss = 0.0
c
        do 204 jh = 1, ngases
c           <recycs,recomns,recxns> [time^{-1}]
        recycs(jh,js-1)=recycs(jh,js-1)+zsurfs*zrparls(jh,jz)   ! heifetz
        recomns(jh,js-1)=recomns(jh,js-1)+zsurfs*zrecoms(jh,jz) ! heifetz
        recxns(jh,js-1)=recxns(jh,js-1)+zsurfs*zbchxs(jh,jz)    ! heifetz
c
          znh(jh) = znh(jh) + rhohs(jh,2,jz)*dx2i(jz)
          zvdist(jh) = zvdist(jh) + (recoms(jh,jz)*rhohs(jh,2,jz) +
     1          shchxs(jh,jz))*dx2i(jz)
          zbloss = zbloss + shblos(jh,jz)
          znn = znn + rhons(jh,jz)
  204   continue
c
          zni = zni + rhoins(2,jz) * dx2i(jz)
          zee = zee + tes(2,jz) * rhoels(2,jz) * dx2i(jz)
          zei = zei + tis(2,jz) * rhoins(2,jz) * dx2i(jz)
          if (znn.gt.epslon) zsigvb = zsigvb + (zbloss/znn)*dx2i(jz)
  208   continue
c
        zdx2 = 2.0 / (xbouni(ib)**2 - xbouni(iz0)**2)   !   ???
        tein(js-1) = zee * evsinv / zne
        tiin(js-1) = zei * evsinv / zni
        denein(js-1) = zne * zdx2
        sigvbx(js-1) = zsigvb * zdx2
c
        do 214 jh = 1, ngases
          denion(jh,js-1) = znh(jh) * zdx2
          snvol(jh,js-1) = zvdist(jh) * zdx2
  214   continue
c
  218   continue
        if (.not. ilprof .and. il0) go to 350                  ! heifetz
c                                         <rplsma,nlpion,nlscat,nlsput>
c
        rplsma = rmins   !   ???
        nlpion = nlgpin
        nlscat = nlgref
        nlsput = nlgspt
c                           <nsrci,npts,nvsrc,e0in,fluxin,fract,nlmono>
        zau = fcau * 10.0**fxau
        zt = tai * uist
        call reseti(nsrci,2*maxgas,0)
        call reseti(npts,maxsrc,ngpart)
        call reseti(nvsrc,maxsrc,0)
        call timint (zt,zztcold,bdtime,20,tcold,1,1)
        call resetr(e0in,maxsrc,zztcold*evsinv*uesh)
c
        call resetr(fluxin,maxsrc,1.0)
        call resetr(fract,maxgas,1.0)
        call resetl(nlmono,maxsrc,.false.)
c                                      <nsrces,nsrci,ngasi,e0in,nlmono>
        nsrces = 0
c
        do 228 jh = 1, ngases
          nsrces = nsrces + 1
          if (nsrces.gt.maxsrc) go to 9022
          nsrci(1,jh) = nsrces
          ngasi(nsrces) = jh
          if (nlglim) e0in(nsrces) = tis(2,mzones) * evsinv
          nlmono(nsrces) = .not.nlglim
c
          if (zpuffi(jh).le.0.0) go to 224
            nsrces = nsrces + 1
            if (nsrces.gt.maxsrc) go to 9022
            nsrci(2,jh) = nsrces
            ngasi(nsrces) = jh
            zt = tai * uist
            call timint (zt,ztcoldp,bdtime,20,tcoldp,1,1)
            e0in(nsrces)   = ztcoldp*evsinv*uesh
            nlmono(nsrces) = nlgmon
  224     continue
c
c               volume sources
c
        if (gvsrci(jh).le.0.0) go to 226
          nsrces = nsrces + 1
          if (nsrces.gt.maxsrc) go to 9022
          nsrci(3,jh) = nsrces
          ngasi(nsrces) = jh
          nvsrc(nsrces) = 1
          e0in(nsrces) = 0.0
          nlmono(nsrces) = .false.
  226   continue
c
        rmass(jh) = aspec(jh) * zau
  228   continue
c
c
cl              set up splitting surfaces
c
c
  250   continue
c
c
cl              ngsplt surfaces
c
c               splitting surfaces are (also rather arbitrarily) placed
c               every so many monte surfaces.
c
c
c
        call resetr(rsplit,maxrad,0.0)
        call resetr(rnu   ,maxrad,0.0)
        imax = min0( maxrad, min0( ngsplt, nsurf ) )
        if (imax.le.0) go to 255
c
        do 254 jr = 1, imax
          is = ((jr - 1) * (nsurf - 1)) / imax + 2
          if (is.ge.nsurf) go to 254
          rsplit(jr) = rsurf(is)
          rnu(jr) = 2.0
  254   continue
c
  255   continue
        ilinit = .true.
                                        call expert(iclass,isub,2)
c
c
cl      3)      get profile  from Monte Carlo computation (AURORA)
c
c
  300   continue
c
c        if(ntty.ne.0) write(ntychl,9001) nstep,ngpart
c 9001   format(x,'monte-carlo neutrals - nstep=',i4,
c     1  2h, ,'ngpart=',i5)
c
        call monte(3)
c        if(ntty.ne.0) write(ntychl,9002)
c 9002   format(x,'monte-carlo computation done')
c
        ngsprf = nstep
        gtprfi = tai
c
c--------1---------2---------3---------4---------5---------6---------7-c
c  ***  two page segement from heifetz
c  to implement, type, rp#1,#1;cheifetz ; ;
c                                                           get profile
cheifetz 300    continue
cheifetz        lneutral=1
cheifetz        if(tbi*uist.gt.tneutral)lneutral=2
c300    if(lswitch.eq.0)go to 302
c       write(ntychl,9300)nstep
c       read(ntychl,9301)lneutral,ngprof
c
cheifetz 302    if(lneutral.eq.2)go to 305
c                                                             run monte
cheifetz        if(ntty.ne.0)write(ntychl,9302)
cheifetz        call monte(3)
cheifetz        if(ntty.ne.0)write(ntychl,9304)
cheifetz        go to 340
c                                                             run degas
cheifetz 305    write(ntychl,9305)
cheifetz        call freeus(lio1)
cheifetz        call msglink(lio1,2)
c
cheifetz        do 335 jh=1,ngases
c                                                 jsrc=1: recycling
c                                                      2: puffing
c                                                      3: volume source
cheifetz        do 335 jsrc=1,3
cheifetz          if(nsrci(jsrc,jh).eq.0)go to 335
c                                     store this profile's source rates
c                                                             <sorcnos>
cheifetz          sorcnos(jh,jsrc)=sorcns(jh,jsrc)
c                         write plasma density and temperature to degas
cheifetz          call destroy("todega")
cheifetz          call freeus(lio2)
cheifetz          call create(lio2,"todega",2,100000b)
cheifetz          write(lio2,9306)ngases,nsrci(jsrc,jh)
cheifetz          do 308 js=1,nsurf-1
cheifetz 308        write(lio2,9308)js,nohzs,tein(js),js,nohzs,tiin(js),
cheifetz     &       js,nohzs,tiin(js),js,nohzs,denein(js),
cheifetz     &       js,nohzs,denion(1,js),js,nohzs,denion(2,js)
c                                                   compute ion sources
c                                                               <pflux>
cheifetz          call resetr(pflux,250,0.0)
cheifetz          zdiam=2.0*fcpi*rmajs
cheifetz          go to(310,315,320)jsrc
c                                              jsrc=1: recycling source
cheifetz 310      zrp=rperps(jh)/(float(kplpuf)*zdiam)
cheifetz          do 312 ip=1,kplpuf
cheifetz 312        pflux(ip)=zrp
cheifetz          do 314 ip=kplpuf+1,kplpuf+kplrec
cheifetz 314        pflux(ip)=recycs(jh,kplshu(ip-kplpuf))/zdiam
cheifetz          jhh=jh
cheifetz          go to 325
c                                                jsrc=2: puffing source
cheifetz 315      if(ngases.eq.2)go to 316
cheifetz          jmol=1
cheifetz          go to 318
cheifetz 316      jmol=2*jh-1
cheifetz 318      zpuf=puffings(jh)/(2.0*float(kplpuf)*zdiam)
cheifetz          do 319 ip=1,kplpuf
cheifetz 319        pflux(ip)=zpuf
cheifetz          jhh=jmol+ngases
cheifetz          go to 325
c                                                 jsrc=3: volume source
cheifetz 320      ztot=zsurfs*gvsrci(jh)/float(kplpuf)*(1.0/(uist*uisl**2))
cheifetz          ztot=ztot/zdiam
cheifetz          do 322 ip=1,kplpuf
cheifetz 322        pflux(ip)=ztot
cheifetz          jhh=jh
c                                                     write out sources
cheifetz 325      do 326 ip=1,kplpuf+kplrec
cheifetz 326        write(lio2,9326)ip,jhh,pflux(ip)
cheifetz          write(lio2,9328)
cheifetz          call close(lio2)
c                                                           start degas
cheifetz          name="xdegbal"
cheifetz          msgtocte(1)="input=ba"
cheifetz          msgtocte(2)="lxxx"
cheifetz          msglen=12
cheifetz 332      call cntelink(name,0,0)
cheifetz          call msgtoe(msgtocte,msglen,0,0,lresulta)
c                                                  check for "all done"
cheifetz          call msgfre(msg,48,48,lresultb)
cheifetz          msgtest=msg(1).and.17777777777777b
cheifetz          msgalld="   all d".and.17777777777777b
cheifetz          if(msgtest.eq.msgalld)go to 334
cheifetz          call cntename(name)
cheifetz          msgtocte(1)=8h
cheifetz          msglen=-1
cheifetz          go to 332
c                                          read source terms from degas
cheifetz 334      write(ntychl,9334)jh,jsrc
cheifetz          call close(lio1)
cheifetz          call destroy(newname(nsrci(jsrc,jh)))
cheifetz          call chgname(7hbalxxxr,newname(nsrci(jsrc,jh)))
cheifetz          call freeus(lio3)
cheifetz          call open(lio3,"tobald",2,len)
cheifetz          read(lio3,degas)
cheifetz          call close(lio3)
cheifetz 335    continue
c                                    divide sources by baldur's volumes
cheifetz        do 338 j=2,nsurf
cheifetz          zvolinv=1.0/(fcpi*(rsurf(j)**2-rsurf(j-1)**2))
cheifetz          do 337 js=1,nsrces
cheifetz            dweions(j-1,js)=dweions(j-1,js)*zvolinv
cheifetz            do 336 jh=1,ngases
cheifetz              sni(jh,j-1,js)=sni(jh,j-1,js)*zvolinv
cheifetz              sne(jh,j-1,js)=sne(jh,j-1,js)*zvolinv
cheifetz              den0(jh,j-1,js)=den0(jh,j-1,js)*zvolinv
cheifetz              dwiions(jh,j-1,js)=dwiions(jh,j-1,js)*zvolinv
cheifetz 336          dwichxs(jh,j-1,js)=dwichxs(jh,j-1,js)*zvolinv
cheifetz 337      continue
cheifetz 338    continue
c                       the source and sink terms for baldur are chosen
c                      by interpolation of the values on the monte grid
c
cheifetz 340    ngsprf=nstep
cheifetz        gtprfi=tai
c
c *** end of 2-page segment from heifetz ***
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c..map from neutral package zones to BALDUR zones
c
c
cl              compute un-normalized total ionization
cvar
c
c  Variables set in sbrtn NEUGAS
c  -----------------------------
c
c
c       in shions --
c               sne*ne is the ionization rate, in part/cc/sec
c               sni*nh is the charge-exchange rate
c
c       in weions --
c               sne*ne*eioniz is the electron energy loss due
c                       to ionizations
c
c       in wiions --
c               sne*tns*3/2*ne = spe*ne*ti is the ion energy added by
c                       the ionized ion
c               3/2*tns*shblos is the ion energy added to the plasma
c                       by neutrals charge-exchanging with beam
c                       ions.
c
c       in wichxs --
c               3/2*tis*(recoms*rhohs + shchxs)  is the ion energy loss
c                       due to recombination and charge exchange
c                       of beam neutrals with plasma
c               - (spi + spii*ti)*nh is the ion energy loss due to
c                       charge exchange
c
c               all must be scaled by the flux, i.e.,
c               if fluxn(jsrc) = the actual flux of influx jsrc / the
c                       value used to get the profile,
c               spe is  sum jsrc of fluxn(jsrc)*spe(*,jsrc)
c
c               the actual fluxes are chosen so that the net source of
c               each species of hydrogen, both from ionization,
c               and from charge-exchange, is equal to the desired influx
c               of that species which has been chosen  for particle
c               conservation.
c
c               the source and sink terms for baldur are chosen
c               by interpolation of the values on the monte grid
c
cend
c
c
  350   continue
c
cheifetz if(lneutral.eq.2)go to 400
        call resetr(ztotal,32,0.0)
        iz = 1
c
        do 378 jz = lcentr, ledge
          z2 = dx2i(jz) * zrads
          z0 = rhoels(2,jz) * z2
c
          z1 = 0.0
          do 351 jh = 1, mhyd
            z1 = z1 + rhohs(jh,2,jz)
  351     continue
          z1 = z1 * z2
c
          if (iz.ge.nsurf) go to 370
          zr = 2.0 * ahalfs(jz,2)
          if (zr.gt.(rsurf(iz) + rsurf(iz+1))) iz = iz + 1
          if (iz.ge.nsurf) go to 370
          if (iz.gt.1) go to 360
c
c               center half of first monte zone
c
          do 354 js = 1, nsrces
          do 354 jg = 1, ngases
            ztotal(jg,js) = ztotal(jg,js) + sne(jg,iz,js)*z0 +
     1                                  sni(jg,iz,js)*z1 +
     2                                  den0(jg,iz,js)*sigvbx(iz)*z2
  354     continue
c
        go to 378
c
c               between middles of center- and outer-most monte zones
c
  360   continue
          zint = (zr - rsurf(iz) - rsurf(iz-1)) /
     1                                  (rsurf(iz+1) - rsurf(iz-1))
c
          do 364 js = 1, nsrces
          do 364 jg = 1, ngases
          ztotal(jg,js) = ztotal(jg,js) +
     1          (sne(jg,iz,js)*zint + sne(jg,iz-1,js)*(1.0-zint))*z0 +
     2          (sni(jg,iz,js)*zint + sni(jg,iz-1,js)*(1.0-zint))*z1 +
     3          (den0(jg,iz,js)*sigvbx(iz)*zint +
     4                  den0(jg,iz-1,js)*sigvbx(iz-1)*(1.0-zint))*z2
  364     continue
c
        go to 378
c
c               outer half of outer-most monte zone
c
  370   continue
c
          do 374 js = 1, nsrces
          do 374 jg = 1, ngases
            ztotal(jg,js) = ztotal(jg,js) + sne(jg,iz-1,js)*z0 +
     1                                  sni(jg,iz-1,js)*z1 +
     2                                  den0(jg,iz-1,js)*sigvbx(iz-1)*z2
  374     continue
c
  378   continue
c
c
c               compute actual influx --
c               sum of all particles either absorbed by plasma
c               or escaped and not reflected (outflx-fluxr)
c               due to a given source.
c
c
        do 388 js = 1, nsrces
          fluxin(js) = 0.0
c
          do 384 jg = 1, ngases
            fluxin(js) = fluxin(js) + ztotal(jg,js) +
     1          outflx(jg,js) - fluxr(jg,js)
  384     continue
  388   continue
c
c
c               remove initial flux
c
        do 398 js = 1, nsrces
          i001 = ngasi(js)
          ztotal(i001,js) = ztotal(i001,js) - fluxin(js)
  398   continue
c
                                        call expert(iclass,isub,3)
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      4)      set source and sink terms
c
c
  400   continue
c
        i1 = mxzone*mxhyd
c
        zsputs = 0.0
        call resetr(wichxs,mxzone,0.0)
        call resetr(dwicxs,mxzone,0.0)
        call resetr(weions,mxzone,0.0)
        call resetr(wiions,mxzone,0.0)
        call resetr(shions,i1    ,0.0)
        call resetr(tns   ,i1    ,0.0)
        call resetr(rhons ,i1    ,0.0)
c
c
c               set fluxn
c
cheifetz if(lneutral.eq.2)go to 412
c                                                            lneutral=1
c
        call resetr(zflux,8,0.0)
        call resetr(zz,64,0.0)
c
        z0 = usil**2 * usit
c
        do 408 jg = 1, ngases
          i001 = nsrci(1,jg)
          i002 = nsrci(2,jg)
          i003 = nsrci(3,jg)
          zflux(i001) = gfluxi(jg) * z0
          if (i002.gt.0) zflux(i002) = zpuffi(jg) * z0
          if (i003.gt.0) zflux(i003) = gvsrci(jg) * z0
c
          is = i001
c
          do 406 js = 1, nsrces
            zz(is,js) = grecyc * ztotal(jg,js)
  406     continue
c
  408   continue
c
c..linear extrapolation of source and sink terms as a function of time
c
c
cheifetz        do 410 js=1,nsrces
cheifetz 410      zz(js,js)=zz(js,js)+fluxin(js)
c                                                              invert zz
cheifetz        call matrx1(zz,8,nsrces,ipivot,ierr)
cheifetz        if(ierr.gt.0)go to 9040
cheifetz        call matrx2(fluxn,zflux,zz,ipivot,8,nsrces,1,1,1)
cheifetz        go to 418
c                                                            lneutral=2
cheifetz 412    do 414 jg=1,ngases
cheifetz          i001=nsrci(1,jg)
cheifetz          if((i001.gt.0).and.(sorcnos(jg,1).eq.0.0))fluxn(i001)=1.0
cheifetz          if((i001.gt.0).and.(sorcnos(jg,1).gt.0.0))
cheifetz     &     fluxn(i001)=sorcns(jg,1)/sorcnos(jg,1)
cheifetz          i002=nsrci(2,jg)
cheifetz          if((i002.gt.0).and.(sorcnos(jg,2).eq.0.0))fluxn(i002)=1.0
cheifetz          if((i002.gt.0).and.(sorcnos(jg,2).gt.0.0))
cheifetz     &     fluxn(i002)=sorcns(jg,2)/sorcnos(jg,2)
cheifetz          i003=nsrci(3,jg)
cheifetz          if((i003.gt.0).and.(sorcnos(jg,3).eq.0.0))fluxn(i003)=1.0
cheifetz          if((i003.gt.0).and.(sorcnos(jg,3).gt.0.0))
cheifetz     &     fluxn(i003)=sorcns(jg,3)/sorcnos(jg,3)
cheifetz 414    continue
c                                                            sputtering
c
        do 414 js = 1, nsrces
          zz(js,js) = zz(js,js) + fluxin(js)
  414   continue
c
        call matrx1(zz,8,nsrces,ipivot,ierr)
        if (ierr.gt.0) go to 9040
        call matrx2(fluxn,zflux,zz,ipivot,8,nsrces,1,1,1)
c
c
c               sputtering
c
        do 418 js = 1, nsrces
          zsputs = zsputs + sflux(js)*fluxn(js)
  418   continue
c
c
c               calculate sources and sinks (scaled by fluxn)
c
c
        iz0 = 1
        iz1 = 1
        icells = nsurf - 1
c
        do 458 jz = lcentr, ledge
c
          zni = 0.0
          do 424 jh = 1, mhyd
            zni = zni + rhohs(jh,2,jz)
  424     continue
c
          if (iz0.ge.icells) go to 428
            zr = 2.0 * ahalfs(jz,2)
            if (zr.le.(rsurf(iz1) + rsurf(iz1+1))) go to 428
            iz1 = iz1 + 1
            iz0 = iz1 - 1
            iz1 = min0(iz1,icells)
  428     continue
c
          zint1 = (zr - rsurf(iz0) - rsurf(iz0+1)) /
     1                                  (rsurf(iz1+1) - rsurf(iz0))
          zint0 = 1.0 - zint1
c  
c
  430   continue
c
c
          zzsi = 0.0
          zzpi = 0.0
          zwbx = 0.0
          zwbxi = 0.0
c
          do 438 jg = 1, ngases
c
            zschex = 0.0
            zsi = 0.0
            ze = 0.0
            zn = 0.0
c
            do 434 jsrc = 1, nsrces
              zschex = zschex +
     1  (sni(jg,iz0,jsrc)*zint0 + sni(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
              zsi = zsi +
     1  (sne(jg,iz0,jsrc)*zint0 + sne(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
              ze = ze + (den0(jg,iz0,jsrc)*eneut(jg,iz0,jsrc)*zint0 +
     1  den0(jg,iz1,jsrc)*eneut(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
              zn = zn +
     1  (den0(jg,iz0,jsrc)*zint0 + den0(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
cheifetz              if(lneutral.eq.1)go to 434
c                                                       <wiions,wichxs>
cheifetz              wiions(jz)=wiions(jz)+(dwiions(jg,iz0,jsrc)*zint0+
cheifetz     &         dwiions(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
cheifetz              wichxs(jz)=wichxs(jz)+(dwichxs(jg,iz0,jsrc)*zint0+
cheifetz     &         dwichxs(jg,iz1,jsrc)*zint1)*fluxn(jsrc)
  434     continue
c                                                              <shions>
cheifetz            if(lneutral.eq.1)shions(jg,jz)=
cheifetz     &       zsi*rhoels(2,jz)+zschex*zni
cheifetz            if(lneutral.eq.2)shions(jg,jz)=zsi
c                                                           <rhons,tns>
c
          shions(jg,jz) = zsi*rhoels(2,jz) + zschex*zni
          rhons(jg,jz) = zn
          tns(jg,jz) = 0.0
          if (zn.gt.epslon) tns(jg,jz) = ze*evs / zn
          zzpi = zzpi + tns(jg,jz) * zsi * 1.5
          zzsi = zzsi + zsi
          zwbx = zwbx + tns(jg,jz)*shblos(jg,jz)*1.5
          zwbxi = zwbxi + recoms(jg,jz)*rhohs(jg,2,jz) + shchxs(jg,jz)
  438   continue
c                                                <weions,wiions,wichxs>
cheifetz          if(lneutral.eq.1)go to 437
cheifetz          do 436 jsrc=1,nsrces
cheifetz 436        weions(jz)=weions(jz)+(dweions(iz0,jsrc)*zint0+
cheifetz     &       dweions(iz1,jsrc)*zint1)*fluxn(jsrc)
cheifetz          go to 458
cheifetz 437      continue
c
          weions(jz) = zzsi * eioniz * uesh * rhoels(2,jz)
          wiions(jz) = zzpi * rhoels(2,jz) + zwbx
c
          zwx = 0.0
          zwxi = 0.0
c
          do 444 jsrc = 1, nsrces
            zwx = zwx +
     1  (spi(iz0,jsrc)*zint0 + spi(iz1,jsrc)*zint1)*fluxn(jsrc)
            zwxi = zwxi +
     1  (spii(iz0,jsrc)*zint0 + spii(iz1,jsrc)*zint1)*fluxn(jsrc)
  444     continue
c
          wichxs(jz) = (zwx*evs + zwxi*tis(2,jz)) * zni * 1.5
     1          - 1.5 * zwbxi * tis(2,jz)
          dwicxs(jz) = zwxi * zni * 1.5 - 1.5*zwbxi
c
c
  458   continue
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
cl      4.2)    impurity sputtering
c
c
        call resetr(gsputs,mximp,0.0)
c
        if (mimp.le.0.or..not.nlgspt) go to 470
c
        ii = 0
        do 464 ji = 1, mimp
          ii = ji
          i001 = ji + lhydn
          if (nzspec(i001).eq.26) go to 465
  464   continue
c
        go to 470
c
  465   continue
        gsputs(ii) = zsputs
c
  470   continue
c
                                        call expert(iclass,isub,4)
        return
c
c
cl      90)     error conditions
c
c
 9010   continue
        call error_olymp(0,iclass,isub,1,
     1          'neutral model cannot handle more than 4 hydrogen')
c
c
c
 9020   continue
        call error_olymp(1,iclass,isub,2,
     1          'no. of neutral zones .gt. no. of baldur zones   ')
        call error_olymp(2,zfac,1,1,'  zfac  ')
c
c
c
 9021   continue
        call error_olymp(0,iclass,isub,2,
     1          'baldur using more hyds. than monte can ')
c
c
c
 9022   continue
        call error_olymp(1,iclass,isub,2,
     1          'more influxes than monte can handle ')
        call error_olymp(3,gfluxi,5,mhyd,'recyclin')
        call error_olymp(2,zpuffi,5,mhyd,'coldflux')
c
c
c
 9040   continue
        call error_olymp(1,iclass,isub,4,
     1          'cannot solve for cold influx                    ')
        call error_olymp(3,zflux,5,8,'fluxes  ')
        call error_olymp(2,zz   ,5,64,'matrix  ')
        stop
c--------1---------2---------3---------4---------5---------6---------7-c
c                                                     format statements
cheifetz 9300   format(x,'nstep=',i3,' lneutral(1-monte,2-degas),ngprof?',
cheifetz     &   ' (i1,x,i2)')
cheifetz 9301   format(i1,x,i2)
cheifetz 9302   format(x,'monte begun')
cheifetz 9304   format(x,'monte done')
cheifetz 9305   format(x,'degas begun')
cheifetz 9306   format(x,'ngases=',i1,x,'ksource=',i1)
cheifetz 9308   format(x,'tehvt(1,',i2,',1)=',i2,'(',f8.1,')'/
cheifetz     &   x,'tihvt(1,',i2,',1,1)=',i2,'(',f8.1,')',x,
cheifetz     &   x,'tihvt(1,',i2,',1,2)=',i2,'(',f8.1,')'/
cheifetz     &   x,'denehvt(1,',i2,',1)=',i2,'(',1pg9.2,')'/
cheifetz     &   x,'denihvt(1,',i2,',1,1)=',i2,'(',1pg9.2,')',
cheifetz     &   x,'denihvt(1,',i2,',1,2)=',i2,'(',1pg9.2,')')
cheifetz 9326   format(x,'pflux(',i2,',',i2,')=',1pg9.2)
cheifetz 9328   format(x,'$end')
cheifetz 9334   format(x,'degas profile jh=',i2,x,'jsrc=',i2,' done')
c
        end
c--------1---------2---------3---------4---------5---------6---------7-c
c@gclear  .../baldur/bald/dneugas.f 
c rgb 22-nov-99 gclear generated from cleargen 
c aes 04-mar-82 add nlpion to comusr and dimension nlog24(11) 
c 
      subroutine gclear 
c 
      include 'cmonprt.m' 
      include 'cmonxs.m' 
      include 'cmonte.m' 
      include 'cmonspl.m' 
c 
c
c  cleargen 8:17 22 Nov 99
c
c sbrtn clear was generated using the cleargen program
c
      e0 = 0.0
      flxfac = 0.0
      pi = 0.0
      rprob = 0.0
      call resetr( timint, 25, 0.0 )
      vel = 0.0
      velx = 0.0
      vely = 0.0
      velz = 0.0
      weight = 0.0
      wtmin = 0.0
      wtot = 0.0
      x0 = 0.0
      y0 = 0.0
      z0 = 0.0
      msurf = 0
      ncell = 0
      nescap = 0
      ngas = 0
      ninc = 0
      njump = 0
      nsrc = 0
      nlrefl = .false.
      sigvcx = 0.0
      sigvi = 0.0
      sigvt = 0.0
      call resetr( sigvx, 2, 0.0 )
      call resetr( tablie, 25, 0.0 )
      maxtab = 0
      call resetr( den0, 300, 0.0 )
      call resetr( den0in, 12, 0.0 )
      call resetr( denein, 25, 0.0 )
      call resetr( denion, 50, 0.0 )
      call resetr( e0in, 6, 0.0 )
      call resetr( eneut, 300, 0.0 )
      call resetr( fluxin, 6, 0.0 )
      call resetr( fluxn, 6, 0.0 )
      call resetr( fluxr, 12, 0.0 )
      call resetr( fract, 2, 0.0 )
      call resetr( outflx, 12, 0.0 )
      call resetr( rmass, 2, 0.0 )
      call resetr( rnu, 10, 0.0 )
      rplsma = 0.0
      call resetr( rsplit, 10, 0.0 )
      call resetr( rsurf, 25, 0.0 )
      call resetr( sflux, 6, 0.0 )
      call resetr( sigma, 300, 0.0 )
      call resetr( sigvbx, 25, 0.0 )
      call resetr( sne, 300, 0.0 )
      call resetr( sni, 300, 0.0 )
      call resetr( snvol, 50, 0.0 )
      call resetr( spe, 150, 0.0 )
      call resetr( spi, 150, 0.0 )
      call resetr( spii, 150, 0.0 )
      call resetr( tein, 25, 0.0 )
      call resetr( tiin, 25, 0.0 )
      call resetr( veff, 12, 0.0 )
      wmin = 0.0
      maxgas = 0
      maxrad = 0
      maxsrc = 0
      maxsur = 0
      ngases = 0
      call reseti( ngasi, 6, 0 )
      call reseti( npts, 6, 0 )
      nrandm = 0
      nsrces = 0
      call reseti( nsrci, 6, 0 )
      nsurf = 0
      call reseti( nvsrc, 6, 0 )
      nlerr = .false.
      nlmod = .false.
      call resetl( nlmono, 6, .false. )
      nlpion = .false.
      nlscat = .false.
      nlsput = .false.
      call resetr( rnumb, 25, 0.0 )
      call resetr( vl, 15, 0.0 )
      call resetr( vxl, 15, 0.0 )
      call resetr( vyl, 15, 0.0 )
      call resetr( vzl, 15, 0.0 )
      call resetr( wl, 15, 0.0 )
      call resetr( xl, 15, 0.0 )
      call resetr( yl, 15, 0.0 )
      call resetr( zl, 15, 0.0 )
      maxlev = 0
      call reseti( mcell, 15, 0 )
      call reseti( msurfl, 15, 0 )
      call reseti( ngasl, 15, 0 )
      nlevel = 0
      call reseti( nodes, 15, 0 )
      nu = 0
      call resetl( nlsplt, 25, .false. )
      nlsurf = .false.

      return 
      end 
c--------1---------2---------3---------4---------5---------6---------7-c
c@monte  /11040/bald92/wbaldn1   DNEUGAS
c
c        amck 23-jan-77 move 'call splits' to after 300
c
         subroutine monte(k)
c
c p.0  organize neutral calculation
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
c
        data iclass/6/,isub/8/
c
c
         go to (100,200,300),k
c---------------------------------------------------------------------
cl              1.         initialize package
c
  100    continue
c
c     clear store
         call gclear
c     set defaults
         call setvar
c
c
         return
c
c---------------------------------------------------------------------
cl              2.         check data
c
  200    continue
c
c
c---------------------------------------------------------------------
cl              3.         monte carlo calculation
c
  300    continue
c
c     set up splitting surfaces
c
         call splits
c
c     clear accumulators
         i1=maxgas*maxsur*maxsrc
         i2=maxsur*maxsrc
         i3=maxgas*maxsrc
c
         call resetr(den0,i1,0.0)
         call resetr(eneut,i1,0.0)
         call resetr(sigma,i1,0.0)
         call resetr(sne  ,i1,0.0)
         call resetr(sni  ,i1,0.0)
         call resetr(spe  ,i2,0.0)
         call resetr(spi  ,i2,0.0)
         call resetr(spii ,i2,0.0)
c
         call resetr(outflx,i3,0.0)
         call resetr(sflux,maxsrc,0.0)
         call resetr(veff,i3,0.0)
c
c     calculate profiles
         call mcarlo
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@setvar  /11040/bald92/wbaldn1   DNEUGAS
c     aes 7-aug-81 wmin now set by baldur input variable gwmin
c     amck 16-jan-78 6 sources
c
         subroutine setvar
c
c p1.3  set default values
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonxs.m'
       include 'cmonte.m'
       include 'cmonspl.m'
c----------------------------------------------------------------------
c
        data iclass/6/,isub/9/
c
c
c     pi
         pi=3.14159
c     normalization for incoming flux
         flxfac=0.25
c     minimum acceptable weight
c        wmin is set equal to baldur input variable gwmin in
c        subroutine neugas near line 205
c     max. number of surfaces
         maxsur=25
c     size of cross-section tables
         maxtab=15
c     max. number of gasses
         maxgas=2
c     max. number of influxes
         maxsrc=6
c     number of particles
         call reseti(npts,maxsrc,500)
c     number of gasses
         ngases=1
c     number of influxes
         nsrc=1
c     atomic weight of neutrals
         call resetr(rmass,maxgas,1.6726e-24)
c     edge density and energy
         call resetr(e0in,maxsrc,3.0)
         call resetr(fluxin,maxsrc,1.0)
c     monoenergetic source
         call resetl(nlmono,maxsrc,.true.)
c     maximum number of levels for splitting
         maxlev=15
c     maximum number of splitting surfaces
         maxrad=10
c     splitting parameter
         call resetr(rnu,maxrad,2.0)
c     proportion of test particles
         call resetr(fract,maxgas,1.0)
c     influx gas type
         call reseti(ngasi,maxsrc,1)
c
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@splits  /11040/bald92/wbaldn1   DNEUGAS
c
         subroutine splits
c
c p1.8  set up splitting surfaces
c
       include 'cmonte.m'
       include 'cmonxs.m'
       include 'cmonspl.m'
c---------------------------------------------------------------------
c
        data iclass/6/,isub/10/
c
c
         call resetr(rnumb ,maxsur,0.0    )
         call resetl(nlsplt,maxsur,.false.)
c
c
         i=nsurf-1
         if(nsurf.le.2)return
         zr2=rsurf(2)*0.5
c
         do 110 js=2,i
         zr1=zr2
         zr2=(rsurf(js+1)+rsurf(js))*0.5
         do 104 jr=1,maxrad
         if(rsplit(jr).ge.zr2.or.rsplit(jr).lt.zr1)go to 104
         nlsplt(js)=.true.
         rnumb(js)=rnu(jr)
  104    continue
  110    continue
c
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@mcarlo  /11040/bald92/wbaldn1   DNEUGAS
c rgb 14-jan-92 remove common/combasOA/
c amck 27-jan-78 change subscr. exp. for icl fortran
c  amck 18-jan-78 use 1/npts instead of 1/wtot in sigma
c  amck 10-jan-78 fluxr does not include initial flux
c  amck 20-dec-77 separate snvol profile for each gas
c  amck 13-dec-77 add comments, normalize snvol, call vsourc, etc.
c
         subroutine mcarlo
c
c p2.2  monte carlo calculation
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonxs.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
       dimension zsave(4,50)
c
       data idim/200/
c
        data iclass/6/, isub/11/
c
c
c-------------------------------------------------------------------
cl              1.         set up collision cross-sections
c
         call xsect
c
cl              normalize fract to 1.0
c
         z0 = 0.0
         do 114 jg = 1, ngases
         z0 = z0 + fract(jg)
  114    continue
c
         do 118 jg = 1, ngases
         fract(jg) = fract(jg) / z0
  118    continue
c
c
c               normalize snvol to 1
c
c               source of gas jg in cell ic is
c               fluxin(nvsrc) * (area of torus / vol. of torus) *
c               snvol(jg,jc)
c               so the integral of snvol over volume
c               is the volume of the torus, for each gas.
c
c
         isurf = nsurf - 1
c
         do 130 jg = 1, ngases
         zstot = 0.0
c
         do 124 jc = 1, isurf
      zstot = zstot + snvol(jg,jc)*(rsurf(jc+1)**2 - rsurf(jc)**2) ! ???
  124    continue
c
         if (zstot.gt.0.0) zstot = rplsma**2/zstot
c
         do 128 jc = 1, isurf
         snvol(jg,jc) = snvol(jg,jc) * zstot
  128    continue
c
  130    continue
c
c
c-------------------------------------------------------------------
cl              2.         monte-carlo calculation
c
c
cl              get one set of profiles for each influx
c
c
c
        do 290 jsrc = 1, nsrces
        nsrc = jsrc
        wtot = 0.0
c
c
c               clear work area
c
        call resetr(zsave,idim,0.0)
c
c
c               follow particles
c
c
        ipts = npts(nsrc)
        do 228 jpts = 1, ipts
c
cbate        write (6,*) ' jpts = ',jpts
c
c               launch particle
c
        ngas = ngasi(nsrc)
        wtmin = wmin
        if (nvsrc(nsrc).gt.0) go to 205
c
c               from edge
c
        call locate
        go to 210
c
c               from random point in volume
c
  205   continue
        call vsourc
  210   continue
c
c
c
c
c               follow one particle
c
        call follow
c
c
c               update sum of n0(jpts)**2
c
c
        do 222 jg = 1, ngases
c
        isurf = nsurf-1
        do 221 jz = 1, isurf
        sigma(jg,jz,nsrc) = sigma(jg,jz,nsrc) +
     1                  (den0(jg,jz,nsrc) - zsave(jg,jz))**2
        zsave(jg,jz) = den0(jg,jz,nsrc)
  221   continue
  222   continue
c
c
c
  228   continue
c
c
c
cl              compute averaged quantities
c
c
c               quantities--
c               den0, sne, sni, spe, spi, spii, den0in, fluxr,
c               outflx, sflux vary linearly with the fluxes.
c
c               sigma, eneut, 1/veff are density-weighted-averaged
c
c       neutral density of gas jg in cell jc =
c               sum over jsrc of (actual flux jsrc / fluxin(jsrc)) *
c               den0(jg,jc,jsrc)
c       neutral temp. of jg in cell jc =
c               (sum jsrc (act. flx jsrc/fluxin(jsrc))*den0(jg,jc,jsrc)
c               *eneut(jg,jc,jsrc))/ neutral dens. jg in cell jc
c       ionization source of gas jg in cell jc =
c               sum jsrc (act. flx jsrc/fluxin(jsrc)) *
c               sne(jg,jc,jsrc) * ne
c       ch.ex. source of gas jg in cell jc =
c               sum jsrc (act flx jsrc/fluxin(jsrc)) *
c               sni(jg,jc,jsrc) * nh
c       ionization loss from elecs. in cell jc =
c               (ionization loss/ionization)*total ionization rate
c       ionization energy source to ions in cell jc =
c               3/2 * sum jsrc (act. flx jsrc/fluxin(jsrc)) *
c               spe(jc,jsrc) * ne * ti
c       ch.ex. loss in cell jc =
c               - 3/2 * sum jsrc (act. flx jsrc/fluxin(jsrc)) *
c               ((spi(jc,jsrc) + spii(jc,jsrc)*ti ) * nh +
c               sum jg eneut(jg,jc,jsrc)*den0(jg,jc,jsrc) * sigvbx(jc))
c
c
c               temps are in ev, dens. in part./cm**3,
c               energy loss/gain in ev/cm**3/sec
c               fluxes in part./cm**2/sec
c
c               nh is the total hydrogen density, ne the electron
c               density, and ti the local ion temperature
c
c               neutral outflux, reflected neutral flux, and
c               sputtered flux are computed from
c               outflx, fluxr, and sflux just as the neutral dens.
c               is from den0.  the edge neutral dens. is
c               computed from den0in in the same way.
c               the variation in the neutral density is density-weighted
c               -averaged from sigma just as the neut temp from eneut.
c               veff is the average of (1/v) for each particle
c               that crosses the boundary, hence 1/veff should be
c               density-weighted-averaged, not veff.
c
c
c             *             *             *             *             *
c
c
c               if volume sources are included (nvsrc(jsrc) .gt. 0 for
c               some jsrc), the charge-exchange loss includes an
c               additional term, total vol. source rate * 3/2 ti,
c               which is a loss.
c
c
c
        zpts = 1.0 / wtot
        zw = fluxin(nsrc) * 2.0*pi * rplsma * zpts   !  ???
c
        do 238 jg = 1, ngases
c
        zfac = 3.0 * 1.6021e-12 / rmass(jg)
        zpfac = 1.6021e-19
c
        isurf = nsurf - 1
        do 234 jz = 1, isurf
c
c               statistical variation
c
        zsig = 0.0
        if (den0(jg,jz,nsrc).gt.0.0) zsig = sqrt(
     1    sigma(jg,jz,nsrc)/den0(jg,jz,nsrc)**2 - 1.0/float(npts(nsrc)))
        sigma(jg,jz,nsrc) = zsig
c
c               neutral temp. and density
c
        zt = 0.0
        if (den0(jg,jz,nsrc).gt.0.0)
     1          zt = eneut(jg,jz,nsrc) / (den0(jg,jz,nsrc) * zfac)
        eneut(jg,jz,nsrc) = zt
c
        zvol = pi * (rsurf(jz+1)**2 - rsurf(jz)**2)   !   ???
        den0(jg,jz,nsrc) = den0(jg,jz,nsrc) * zw / zvol
  234   continue
c
  238   continue
c
c
c
cl              source terms (ionization and charge exchange)
c
c               note that these are at present derived just from the
c               neutral density and temperature profiles.
c
c
c
        do 258 jz = 1, isurf
c
        zni = 0.0
c
        do 244 jg = 1, ngases
        zni = zni + denion(jg,jz)
  244   continue
c
        znir = 1.0 / zni
        zner = 1.0 / denein(jz)
c
c
        do 254 jg = 1, ngases
c
        ngas = jg
        zfac = 3.0 * 1.6021e-12 / rmass(jg)
c
        iz = jz
        z = fpath(iz,sqrt(zfac*eneut(jg,jz,nsrc)))
        zdt = eneut(jg,jz,nsrc) / tiin(jz)
c
        do 248 jg2 = 1, ngases
        z0 = sigvx(jg2) * den0(jg,jz,nsrc) * znir
        sni(jg,jz,nsrc) = sni(jg,jz,nsrc) + z0
        sni(jg2,jz,nsrc) = sni(jg2,jz,nsrc) - z0
        spi(jz,nsrc) = spi(jz,nsrc) + eneut(jg,jz,nsrc)*z0
        spii(jz,nsrc) = spii(jz,nsrc) - z0
  248   continue
c
        sne(jg,jz,nsrc) = den0(jg,jz,nsrc) * sigvi*zner
        spe(jz,nsrc) = spe(jz,nsrc) + zdt * sne(jg,jz,nsrc)
  254   continue
c
  258   continue
c
c
cl              fluxes
c
c
c
        do 268 jg = 1, ngases
        den0in(jg,nsrc) = fluxin(nsrc) * zpts * veff(jg,nsrc) / flxfac
        if (veff(jg,nsrc).gt.0.0) veff(jg,nsrc) =
     1          (fluxr(jg,nsrc) + outflx(jg,nsrc)) / veff(jg,nsrc)
        fluxr(jg,nsrc) = fluxr(jg,nsrc) * fluxin(nsrc) * zpts
        outflx(jg,nsrc) = outflx(jg,nsrc) * fluxin(nsrc) * zpts
  268   continue
c
        i001 = ngasi(nsrc)
        if (nvsrc(nsrc).le.0)
     1  fluxr(i001,nsrc) = fluxr(i001,nsrc) - fluxin(nsrc)
        sflux(nsrc) = sflux(nsrc) * fluxin(nsrc) * zpts
c
  290   continue
c
c
c----------------------------------------------------------------------
cl              3.         terminate run
c
c
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@follow  /11040/bald92/wbaldn1   DNEUGAS
c       amck 10-jan-78 no splitting/russ.roulette if volume source
c
         subroutine follow
c
c 2.3  follow path of one particle
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonspl.m'
       include 'cmonxs.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
        data iclass/6/,isub/12/
c
c
c---------------------------------------------------------------------
cl              1.         prologue
c
  100    continue
c     level marker
         nlevel=0
c
c     each track starts at this point
  101    continue
         njump=0
         do 102 j=1,nsurf
  102    timint(j)=0.0
c
c
c---------------------------------------------------------------------
cl              2.         find event
c
  200    continue
c     fetch random number
         zeps=ranz()
         zlog=-log(zeps)
c     clear work variables
         zint1=0.0
         zint2=0.0
         zt1=0.0
         zt2=0.0
         zt=0.0
c
cl                  2.1      scan over segment
  210    continue
c     cell number
         ic=ncell
c     local mean free path
         zmfp=fpath(ic,vel)
c     *timestep*
         call timei(zt1)
         zdt1=zt1-zt2
c     distance travelled in current cell
         zds=vel*zdt1
c     update integral
         zint1=zint1+zds/zmfp
c     check for event
         if(zint1.ge.zlog) go to 220
c     update contribution to density and energy
         zdn=zdt1*weight
         zde0=zdn*vel*vel
         den0(ngas,ic,nsrc)=den0(ngas,ic,nsrc)+zdn
         eneut(ngas,ic,nsrc)=eneut(ngas,ic,nsrc)+zde0
c     next cell
         zt2=zt1
         zint2=zint1
         zt=zt+zdt1
c     next cell - check for escape
         ncell=ncell+ninc
         if(ncell.lt.nsurf) go to 211
c     particle has escaped
         go to 380
  211    continue
c     check if we have encountered a splitting zone
         if(nlsplt(msurf).and.nvsrc(nsrc).le.0) go to 230
c     no - continue tracking
         go to 210
c
c
cl                  2.2      point of collision
  220    continue
         zdt1=(zlog-zint2)*zmfp/vel
         zt=zt+zdt1
         x0=x0+velx*zt
         y0=y0+vely*zt
         z0=z0+velz*zt
c     correct density and temperature
         zdn1=zdt1*weight
         den0(ngas,ic,nsrc)=den0(ngas,ic,nsrc)+zdn1
         eneut(ngas,ic,nsrc)=eneut(ngas,ic,nsrc)+zdn1*vel*vel
         go to 300
c
cl                  2.3      splitting zone
  230    continue
c     check for russian roulette
         if(ninc.gt.0) go to 340
c     no - split particle
         x0=x0+velx*zt
         y0=y0+vely*zt
         z0=z0+velz*zt
         zt=0.0
c     particle on splitting surface
         nlsurf=.true.
         go to 330
c
c
c---------------------------------------------------------------------
cl              3.         process event
c
  300    continue
c     collision has ocurred - calculate a new weight
         if(weight.lt.wtmin) go to 350
c
c     decide on species
c
         z = ranz()
         do 302 jg = 1, ngases
         ngas = jg
         z = z - fract(ngas)
         if (z.le.0) go to 303
  302    continue
  303    continue
c
         weight = weight * sigvx(ngas) / (sigvt * fract(ngas))
c
cl                  3.1      ionization
  310    continue
c
cl                  3.2      charge exchange
  320    continue
c     fetch a new velocity
         call veloc(ic)
         go to 101
c
cl                  3.3      splitting
  330    continue
         nlevel=nlevel+1
         if(nlevel.le.maxlev) go to 331
c     too many levels - no more splitting
         nlevel=nlevel-1
         go to 101
  331    continue
         znu=rnumb(msurf)
         nu=int(znu)
         inu1=nu+1
         if(ranz().gt.(float(inu1)-znu)) nu=inu1
c     split into several particles with reduced weight
         weight=weight/znu
c     save location,weight and temperature at current level
         mcell(nlevel)=ncell
         msurfl(nlevel)=msurf
         xl(nlevel)=x0
         yl(nlevel)=y0
         zl(nlevel)=z0
         wl(nlevel)=weight
         vxl(nlevel)=velx
         vyl(nlevel)=vely
         vzl(nlevel)=velz
         vl(nlevel)=vel
         ngasl(nlevel) = ngas
c     number of nodes at this level
         nodes(nlevel)=nu
         go to 101
c
cl                  3.4      russian roulette
  340    continue
c     probablity of death
         znu=rnumb(msurf)
         znu1=1.0-1.0/znu
         zeps=ranz()
c     is this particle to be killed
         if(zeps.gt.znu1) go to 341
c     yes - process next node at this level
         go to 360
c
  341    weight=weight*znu
         go to 210
c
cl                  3.5      weight too small
  350    continue
c     follow particle with same weight until lost
         zeps=ranz()*sigvt
         if(zeps.le.sigvcx) go to 352
c     ionize particle
         go to 360
  352    continue
c
c     charge exchange
c
c               choose gas type
c
        do 354 jg = 1, ngases
        zeps = zeps - sigvx(jg)
        ngas = jg
        if (zeps.le.0.0) go to 356
  354   continue
  356   continue
c
c     charge exchange
         call veloc(ic)
         go to 101
c
cl                  3.6      next node at current level
  360    continue
         if(nlevel.le.0) return
c     number of remaining nodes
         inodes=nodes(nlevel)-1
         nodes(nlevel)=inodes
         if(inodes.le.0) go to 370
c     restore variables
         ncell=mcell(nlevel)
         msurf=msurfl(nlevel)
         x0=xl(nlevel)
         y0=yl(nlevel)
         z0=zl(nlevel)
         weight=wl(nlevel)
         ngas=ngasl(nlevel)
         velx=vxl(nlevel)
         vely=vyl(nlevel)
         velz=vzl(nlevel)
         vel=vl(nlevel)
         nlsurf=.true.
         go to 101
c
cl                 3.7      return to previous level
  370    continue
         nlevel=nlevel-1
         if(nlevel.le.0) return
         go to 360
c
cl                 3.8      particle escapes
  380    continue
c     point of escape
         x0=x0+velx*zt
         y0=y0+vely*zt
         z0=z0+velz*zt
c
         call escape
c     reflect particle if required
         if(nlrefl) go to 101
c
c     next node at this level
         go to 360
c
c
         end
c---------------------------------------------------------
c@escape  /11040/bald92/wbaldn1   DNEUGAS
         subroutine escape
c
c 2.6  process escaping particles
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
c
        data iclass/6/,isub/13/
c
c
cl              1.         net outflux
c
         outflx(ngas,nsrc)=outflx(ngas,nsrc)+weight
         veff(ngas,nsrc)=veff(ngas,nsrc)+(weight/vel)
         e0=0.5*rmass(ngas)*vel**2 / 1.6021e-12
c
c
c--------------------------------------------------------------------
cl              2.         sputtering
c
         if(nlsput) call sputer
c
c
c-----------------------------------------------------------
cl              3.         reflection
c
         if(nlscat) call reflec
c
c
c
         return
         end
c----------------------------------------------------
c@locate  /11040/bald92/wbaldn1   DNEUGAS
c amck 17-jan-78 wtot is total of weights after particle enters plasma
c
         subroutine locate
c
c 2.8  locate monte-carlo particle on boundary
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonspl.m'
       include 'cmonte.m'
       include 'cmonxs.m'
c---------------------------------------------------------------------
        data iclass/6/,isub/14/
c
       data zpi/3.14159/
c
c
c---------------------------------------------------------------------
cl              1.         locate position
c
c     set incoming particle weight
         weight=1.0
         ngas=ngasi(nsrc)
c     set neutral temperature
         e0=e0in(nsrc)
c     position on surface
         x0=rplsma
         y0=0.0
         z0=0.0
c     first cell
         ncell=nsurf-1
c     first surface
         msurf=nsurf
c     indicate that particle on surface
         nlsurf=.true.
c
c
c---------------------------------------------------------------------
cl              2.         maxwellian source
c
  200    continue
         if(nlmono(nsrc)) go to 300
c     select energy from maxwellian at edge temperature
  201    i=nsurf-1
         zsave=tiin(i)
         tiin(i)=e0in(nsrc)
         call veloc(i)
         tiin(i)=zsave
         if(velx.lt.0.0) velx=-velx
c     select energy from reflection model
         call reflec
         nlrefl=.false.
         go to 400
c
c
c---------------------------------------------------------------------
cl              3.         monoenergetic source
c
  300    continue
c     neutral velocity
         zvel=sqrt(e0*3.2e-12/rmass(ngas))
c
cl                  3.1      theta
  310    continue
         zeps=ranz()
         ztheta=pi*zeps
         zcthet=cos(ztheta)
         zsthet=sin(ztheta)
c
cl                  3.2      phi
  320    continue
         zeps=ranz()
         zsin2=zeps
         zcos2=1.0-zeps
         zsphi=sqrt(zsin2)
         zcphi=sqrt(zcos2)
c
cl                  3.3      velocity
  330    continue
         zvels=zvel*zsphi
         velx=-zvel*zcphi
         vely=zvels*zsthet
         velz=zvels*zcthet
         vel=zvel
c
  331    continue
c     influx
         zwt=weight/vel
         fluxr(ngas,nsrc)=fluxr(ngas,nsrc) + weight
         veff(ngas,nsrc)=veff(ngas,nsrc)+zwt
c
c     mean energy
c
         go to 400
c
c-----------------------------------------------------------------------
cl              4.         weight
c
c
  400    continue
         wtot=wtot+weight
c
         return
         end
c
c-----------------------------------------------
c@veloc  /11040/bald92/wbaldn1   DNEUGAS
         subroutine veloc(k)
c
c p2.9  fetch a new velocity
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
c
        data iclass/6/,isub/15/
c
c
         zarg=sqrt(1.6e-12*tiin(k)/rmass(ngas))
c
         velx=fgauss(zarg,0.0)
         vely=fgauss(zarg,0.0)
         velz=fgauss(zarg,0.0)
c
         vel=sqrt(velx*velx+vely*vely+velz*velz)
         e0=3.125e+11*rmass(ngas)*vel*vel
c
c
         return
         end
c
c------------------------------------------------------
c@timei  .../baldur/code/bald/dneugas.f
c rgb 16-jun-96 set zeps=0.0 when small enough
c rgb 30-may-96 added save zab, zab2, zc1, zr
c
         subroutine timei(pt)
c
c p2.10 calculate intersection times
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonspl.m'
       include 'cmonte.m'
       include 'cmonxs.m'
c---------------------------------------------------------------------
        data iclass/6/,isub/16/
c
        save zab, zab2, zc1, zr
c
cl              1.         prologue
c
        zepslon = 1.e-8
c
  100    continue
c
c     coefficients of quadratic
c
       if ( njump .eq. 0 ) then
         za=velx*velx+vely*vely
         zb=x0*velx+y0*vely
         if(zb.ge.0.0) zb=max(1.e-8,zb)
         zb=x0*velx+y0*vely
         if(zb.ge.0.0) zb=max(1.e-8,zb)
         if(zb.lt.0.0) zb=min(-1.e-8,zb)
         zb2=zb*zb
         zc1=x0*x0+y0*y0
         zab=-zb/za
         zab2=za/zb2
c
c     ensure correct result if starting on splitting surface
c
         if  ( nlsurf ) then
           zr=rsurf(msurf)
           zc1=zr*zr
         endif
c
c     test for direction
c
         ninc=-1
         if(zb.gt.0.0) ninc=1
c
      endif
c
c---------------------------------------------------------------------
cl              2.         intersection times
c
  200    continue
c
c     required surface
c
         msurf=ncell+(1+ninc)/2
         if(msurf.eq.1) go to 201
         if ( timint(msurf) .gt. 0.0 ) go to 303
         zr=rsurf(msurf)
c
c     check for root
c
         if ( abs(zc1-zr*zr) .lt. zepslon*max(abs(zc1),abs(zr*zr))) then
            zeps = 0.0
         else
            zeps = zab2*(zc1-zr*zr)
         endif
c
         if(zeps.lt.1.0) go to 202
c
c     no root - path in other direction
c
  201    continue
         ninc=-ninc
         msurf=msurf+1
         if(njump.eq.0) go to 200
         go to 303
  202    continue
c
c     reset branch marker
c
         njump=1
         nlsurf=.false.
c
c     intersection times
c
         zsqrt=1.0+sqrt(1.0-zeps)
         zt1=zab*zeps/zsqrt
         zt2=zab*zsqrt
c
c
c---------------------------------------------------------------------
cl              3.         return result
c
  300    continue
c
c     check sign of result
c
         if((zt1.gt.0.0).and.(zt2.gt.0.0)) go to 301
c
c     one root negative or zero
c
         pt=max(zt1,zt2)
c
c     one root negative or zero
c
         pt=max(zt1,zt2)
c
c
         return
c
  301    continue
c
c     both roots positive - return one and save other
c
         pt=min(zt1,zt2)
  302    timint(msurf)=max(zt1,zt2)
c
c
         return
c
c     root already known
c
  303    pt=timint(msurf)
c
c
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c
c@vsourc  /11040/bald92/wbaldn1   DNEUGAS
c       amck 9-jan-78 try again if weight.le.wmin
c       amck 4-jan-78 weight = snvol*zx, wtmin .gt. wmin**2
c       amck 2-jan-78 ncell = jc - 1
c       amck 20-dec-77 let ngas have already been set (by mcarlo)
c
        subroutine vsourc
c
c p2.11 start a monte-carlo particle at a random position
c
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonspl.m'
       include 'cmonte.m'
       include 'cmonxs.m'
c----------------------------------------------------------------------
c
c
c       snvol is assumed to have been normalized (mcarlo)
c       so that the integral over volume of
c       snvol is 1.0 for each gas species.
c
c       the particles have a uniform distribution in radius.
c       the initial weight is adjusted so
c       that for a sufficiently large number of samplings,
c       the sum of the weights of all particles produced
c       in a zone for an ion species will be proportional to
c       snvol for that zone and species.
c       the sum of the weights for all the particles produced
c       should approach the number of particles produced,
c       for sufficiently large number of samples.
c
c
c
c       if the weight of the particle produced is less
c       than "wmin", we treat it as 0, and try again.
c       nsurf is at the moment used as the number of tries,
c       for lack of a better number.
c
c
c----------------------------------------------------------------------
c
c
c               location        data iclass/6/,isub/17/
c
c
 
c
  100   continue
c
c
c
        do 124 jtries = 1, nsurf
c
        zx = ranz()
        x0 = zx*rplsma
        y0 = 0.0
        z0 = 0.0
c
c               cell number
c
        nlsurf = .false.
c
        do 105 jc = 2, nsurf
        ncell = jc - 1
        if (x0.eq.rsurf(jc)) go to 106
        if (x0.le.rsurf(jc)) go to 107
  105   continue
  106   continue
        nlsurf = .true.
  107   continue
c
c
c               weight
c
        weight = snvol(ngas,ncell) * zx
        if (weight.gt.wmin) go to 200
  124   continue
c
c
c
  200   continue
c
        wtot = wtot + weight
        wtmin = max(wmin * weight, wmin**2)
c
c
c               velocity
c
c
        call veloc(ncell)
c
c
c
        return
        end
c
c-----------------------------------------------------
c@xsect  /11040/bald92/wbaldn1   DNEUGAS
c
         subroutine xsect
c
c a.1   set up tables of cross-sections
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonxs.m'
       include 'cmonte.m'
c--------------------------------------------------------------------
        data iclass/6/,isub/18/
c
cl              1.         scan over zones
c
         isurf=nsurf-1
         do 203 j=1,isurf
c
c--------------------------------------------------------------------
cl              2.         electron impact ionization
c
c     local electron temp.
         zte=tein(j)
         zx=log10(zte)
c     cross-section (duchs)
         if(zte.gt.20.0) go to 201
         zlog=-3.054*zx-15.72*exp(-zx)+1.603*exp(-zx*zx)
         go to 202
  201    continue
         zlog=-0.5151*zx-2.563/zx-5.231
  202    continue
         tablie(j)=10.0**zlog
  203    continue
c
c
         return
         end
c
c-------------------------------------------------------
c       aes 04-mar-82 added nlpion: switches off proton ionization
c       aes 22-feb-82 added simple proton ionization model
c
         function fpath(k,pvel)
c
c a.2   calculate mean free path
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonxs.m'
       include 'cmonte.m'
c-------------------------------------------------------------------
        dimension zpicof(7)
        data zpicof/-42.03309, 3.557321, -1.045134, .3139238,
     1              -.07454475, .8459113e-2, -.3495444e-3/
        data iclass/6/,isub/19/
c
c
cl              1.         local plasma parameters
c
         zne=denein(k)
         zti=tiin(k)
c
c-------------------------------------------------------------------
cl              2.         electron ionization
c
         zigvie=tablie(k)*zne
c
c
c-------------------------------------------------------------------
cl              3.         charge exchange, proton impact ionization
c
         zigvii=0.0
         sigvcx = 0.0
c
         do 328 jg = 1, ngases
         zvrel2 = 1.6021e-12 * 8.0*zti/(rmass(jg)*pi) + pvel**2
c    5.22e-13 = 1.6726e-24 / 2 / 1.6021e-12, i.e., (cm/sec)**2
c     to ev (assuming proton mass)
         zerel = 5.22e-13 * zvrel2
         zlerel = log(zerel)
         zsig = 0.6937 * (1.0 - 0.0673*zlerel)**2 /
     1          (1.0e+14 + 0.1112 * exp(3.3*zlerel))
         sigvx(jg) = denion(jg,k) * zsig * sqrt(zvrel2)
         sigvcx = sigvcx + sigvx(jg)
c
c     proton ionization is computed just like electron ionization
c      -- fit from freeman and jones
c
         if (.not. nlpion) go to 328
         if (zerel.lt.100.) go to 328
         z10 = log(.001*zerel)
         z0 = zpicof(7)
         do 325 j=6,1,-1
         z0 = z0*z10 + zpicof(j)
  325    continue
         zigvii = zigvii +  exp(z0) * denion(jg,k) * sqrt(zvrel2)
c
  328    continue
c
c
c     total ionization
         sigvi=zigvie+zigvii
c
c     total n-sigma-v; free path
         sigvt=sigvi+sigvcx+sigvbx(k)
         fpath=pvel/sigvt
c
c         return
         end
c
c@frandm   .../baldur/code/bald/dneugas.f
c rgb 04-may-96 changed ranf to ran2
c---------------------------------------------------
         function frandm(k)
c
c a.3  fetch a random number
c--------------------------------------------------------------------
        data iclass/6/,isub/20/
c
        data idum / -1 /
        save idum
c
c
         z = ran2( idum )
         frandm=z
c
cbate      write (*,*) ' frandm = ',frandm
c
         return
         end
c@fgauss  .../baldur/code/bald/dneugas.f
c rgb 29-may-96  added save zx, zy, zsq, zsin, zcos, zt
c--------1---------2---------3---------4---------5---------6---------7-c
c
        function fgauss(zs,za)
c
c a.4 sample from a gaussian on the pdp-10
c
c
c       this function samples from a gaussian of the
c       form exp(-(x-za)**2/(2.*zs**2))
c
c       it is a variation on the box-mueller method
c
        data iclass/6/,isub/21/
c
        data i/1/
c
        save zx, zy, zsq, zsin, zcos, zt
c
cl      1. branch if this is an even numbered call
c
        if(i.eq.2) go to 400
c
c       2. compute the sine and cosine of 2*pi*ranz()
c       by a rejection technique
c
200     continue
        zx=2.*ranz()-1.
        zy=2.*ranz()-1.
        zsq=zx*zx+zy*zy
        if(zsq.gt.1.) go to 200
        zsq=1./zsq
        zsin=2.*zx*zy*zsq
        zcos=(zx*zx-zy*zy)*zsq
c
c       3.      compute rest of selection
c
        zt=sqrt(-2.*log(ranz()))
        zx=zt*zcos
        zy=zt*zsin
        fgauss=zs*zx+za
        i=2
        return
c
c       4. return to get second gaussian point
c
400     continue
        fgauss=zs*zy+za
        i=1
c
        return
        end
c--------------------------------------------------
c@sputer  /11040/bald92/wbaldn1   DNEUGAS
c amck 24-apr-78 remove zcphi reference
c amck 21-mar-78 remove random number
c
         subroutine sputer
c
c a.6  calculate sputtering
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonte.m'
c---------------------------------------------------------------------
c
       dimension zengy(9), zcoeff(9)
c
        data iclass/6/,isub/22/
       data irange/9/
c
       data zengy/70.0,300.0,600.0,1000.0,2000.0,4000.0,
     +            6000.0,8000.0,10000.0/
c
       data zcoeff/1.0e-05,1.4e-03,6.0e-03,1.0e-02,1.0e-02,
     +             6.8e-03,5.2e-03,4.0e-03,3.8e-03/
c
c-------------------------------------------------------------------
c
cl              1.         threshold
c
c     we take lower end of table as threshold - ok to within 4ev
         zth=70.0
c     test for sputtering
         if(e0.lt.zth) return
c
c
c-------------------------------------------------------------------
cl              2.         angle of incidence
c
  200    continue
         zcos=(x0*velx+y0*vely)/(rplsma*vel)
c     cut off (80 deg)
         zmin=0.17
         zfac=0.0
         if(zcos.ge.zmin) zfac=1.0/zcos
c
c
c-------------------------------------------------------------------
cl              3.         sputtering coefficient
c
  300    continue
c     interpolate to find normal incidence coefficient
         do 301 j=1,irange
         ij=j
         ze1=zengy(j)
         if(ze1.gt.e0) go to 302
  301    continue
c
  302    continue
         zc1=zcoeff(ij-1)
         zc2=zcoeff(ij)
         ze1=zengy(ij-1)
         ze2=zengy(ij)
         zalpha=(zc1*(ze2-e0)+zc2*(e0-ze1))/(ze2-ze1)
c     angular effects
         zalpha=zalpha*zfac
c
c
c--------------------------------------------------------------------
cl             4.         emission
c
  400    continue
c
c     flux into plasma
         sflux(nsrc)=sflux(nsrc)+weight*zalpha*0.5
c
c
         return
         end
c----------------------------------------------------------
c@reflec  /11040/bald92/wbaldn1   DNEUGAS
c       aes 20-aug-81 eliminate weight change at 250
c amck use zrange instead of zengy in e0 calculation
c
         subroutine reflec
c
c a.7  reflect escaping particles
c
c---------------------------------------------------------------------
       include 'cmonprt.m'
       include 'cmonspl.m'
       include 'cmonte.m'
       include 'cmonxs.m'
c---------------------------------------------------------------------
c
       dimension
     +   zrange(11), zde(11), zengy(11), zr(11), zide(11,11)
c
        data iclass/6/,isub/23/
c
c     size of tables
         data idim/11/
c     energy range
       data
     +   zrange/14.7,31.62,68.1,146.8,316.2,681.9,1468.0,
     +          3162.0,6813.0,14678.0,100000.0/
c     deltae
       data
     +   zde/14.7,16.92,36.51,78.66,169.2,365.1,786.6,1694.0,
     +       3651.0,7865.0,0.0/
c     energy
       data
     +   zengy/10.0,21.5,46.4,100.0,215.4,464.1,1000.0,2154.3,
     +         4641.3,10000.0,0.0/
c     reflection probability
       data
     +   zr/0.8,0.7,0.62,0.543,0.46,0.37,0.29,0.21,0.14,
     +      0.095,0.03/
c     i*deltae
       data
     +zide/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +     0.2,0.8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +     0.05,0.3,0.65,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +     0.025,0.075,0.35,0.55,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     +     0.02,0.03,0.1,0.4,0.45,0.0,0.0,0.0,0.0,0.0,0.0,
     +     0.005,0.015,0.035,0.12,0.45,0.375,0.0,0.0,0.0,0.0,0.0,
     +     0.003,0.008,0.018,0.05,0.145,0.475,0.301,0.0,0.0,0.0,0.0,
     +     0.003,0.004,0.009,0.024,0.065,0.196,0.47,0.23,0.0,0.0,0.0,
     +     0.001,0.002,0.004,0.011,0.032,0.09,0.26,0.43,0.17,0.0,0.0,
     +     0.001,0.002,0.003,0.005,0.014,0.048,0.142,0.29,0.36,0.135,0.,
     +     0.001,0.002,0.002,0.004,0.006,0.02,0.07,0.2,0.295,0.30,0.10/
c
c---------------------------------------------------------------------
cl              1.         prologue
c
  100    continue
c     tentatively assume no reflection
         nlrefl=.false.
c     no reflection below emin
         zmin=0.0
         zwmin=wtmin
         if((e0.le.zmin).or.(weight.lt.zwmin)) return
c
c
c---------------------------------------------------------------------
cl              2.         reflect particle
c
  200    continue
c     find energy range
         do 201 j=1,idim
         irange=j
  201    if(e0.le.zrange(j)) go to 202
  202    continue
c
cl                  2.1      test for reflection
  210    continue
c     reflection probability
         rprob=zr(irange)
c     reflect this particle
         nlrefl=.true.
         ncell=nsurf-1
         msurf=nsurf
         nlsurf=.true.
c
cl                  2.2      reflection energy
  220    continue
         zeps=ranz()
c     number of boxes in this range
         in=irange
         in1=in-1
c     find reflection energy
         za=0.0
         ze=0.0
c
         if(in1.le.0) go to 222
c
         do 221 j=1,in1
         i=j
         zdelta=zde(j)
         ze=ze+zdelta
         za=za+zide(j,irange)
  221    if(za.gt.zeps) go to 223
c
c     last box
  222    continue
         i=in
         zdelta=zengy(in)-ze
         ze=zengy(i)
         za=za+zide(i,irange)
  223    continue
c     reflection energy
         ze0=ze-(za-zeps)*zdelta/zide(i,irange)
         e0=e0*(ze0/zrange(in))
c
cl                  2.3      angle of reflection
  230    continue
         zthet=pi*ranz()
         zsthet=sin(zthet)
         zcthet=cos(zthet)
         zsphi=sqrt(ranz())
         zcphi=sqrt(1.0-zsphi*zsphi)
c
cl                  2.4      new velocity components
  240    continue
         vel=sqrt(3.2e-12*e0/rmass(nsrc))
         zvr=-vel*zcphi
         zvp=vel*zsphi*zsthet
         velx=(x0*zvr-y0*zvp)/rplsma
         vely=(y0*zvr+x0*zvp)/rplsma
         velz=vel*zsphi*zcthet
c
c
c---------------------------------------------------------------------
cl              3.         boundary conditions
c
c     update mean energy of edge neutrals
         zwt=weight/vel
         fluxr(ngas,nsrc)=fluxr(ngas,nsrc)+weight
         veff(ngas,nsrc)=veff(ngas,nsrc)+zwt*0.75
c
c
         return
         end
c--------1---------2---------3---------4---------5---------6---------7-c

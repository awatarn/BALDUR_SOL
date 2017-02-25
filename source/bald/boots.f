c--------1---------2---------3---------4---------5---------6---------7-c
c@boots  .../baldur/code/bald/dsolver.f
c     ces 28-jul-88 16.16 add flr cutoff to chang current
c     dps 18-oct-88 15.06 compute trapped particle current ajtpbi(j) from
c                   Bateman and Chang.
c     dps 07-may-87 use rhoins for ion density; switch avi(jz,2,4) to
c                   avi(jz,2,1) in zdrho; use ni*dTi/d(psi) in Jboot.
c     dps 04-mar-87 inserted with bootstrap current cf. Mike Hughes
c
         subroutine boots(pjboot,pjtpbi,kfirst,klast)
c
C. Calculate bootstrap current
c
c  ajboot(j) = < bootstrap current density / R > [Amps / m**3]
c              at BALDUR zone centers
c  ajtpbi(j) = < current density due to trapped particles in the
c              banana regime / R > [Amps / m**3]
c              at BALDUR zone centers
c              computed with an analytic formula by C. S. Chang, 1988
c
       include 'cparm.m'
       include 'cbaldr.m'
       include 'commhd.m'
       include 'comtrp.m'
c
C---------------------------------------------------------------------
       dimension   pjboot(*), pjtpbi(*)
c
c-----------------------------------------------------------------------
CL              1.         Scan Zone Centres
c
C     Conversion factor to MKS
         zfac=0.1*uish*uisd
C     Temperature will have funny units but we are interested
C     only in scale lengths.  Maybe one day PPPL will learn
C     about MKS!
C
c
c zone-independent quantities for flr cutoff of chang current
c
      zb=bzs*usib
      zckb=cfev*10.**(fxk)
      zcmp=fcmp*10.**(fxnucl)*usim
      zce=fce*10.**(fxe)*10.0
         do 220 jz=kfirst,klast
c
C     Variables in current zone
         zpe=chi(lelec,jz)*zfac
         zpi=chi(lion,jz)*zfac
         zne=rhoels(2,jz)*usid
         zni=rhoins(2,jz)*usid
         zte=zpe/zne
         zti=zpi/zni
c
CL                  1.1      Point ahead
         zpep=chi(lelec,jz+1)*zfac
         zpip=chi(lion,jz+1)*zfac
         znep=rhoels(2,jz+1)*usid
         znip=rhoins(2,jz+1)*usid
         ztep=zpep/znep
         ztip=zpip/znip
c
CL                  1.2      Point behind
         zpem=chi(lelec,jz-1)*zfac
         zpim=chi(lion,jz-1)*zfac
         znem=rhoels(2,jz-1)*usid
         znim=rhoins(2,jz-1)*usid
         ztem=zpem/znem
         ztim=zpim/znim
c
c-----------------------------------------------------------------------
CL              2.         Calculate Bootstrap Current
c
c
CL                  2.1      Gradient scale lengths
         zdrho=0.5*(avi(jz,2,1)+avi(jz+1,2,1))
         zdx=(xzoni(jz+1)-xzoni(jz-1))*zdrho
         zdxb=(xbouni(jz+1)-xbouni(jz))*zdrho
         zlp=(zpep+zpip-zpem-zpim)/(zpe*zdx)
         zlte=(ztep-ztem)/(zte*zdx)
         zlti=(ztip-ztim)/(zti*zdx)
c
CL                  2.2      Bootstrap current
         zdpsi=bpols(2,jz)*r0ref*1.0e-04
C     Transport coefficients
         zl13=rl13(jz,2)*zpe/zdpsi
         zl23=rl23(jz,2)*zpe/zdpsi
         zy=rly(jz,2)
c
C     <Jboot.Grad(phi)>
         pjboot(jz)=-zl13*(zlp-2.5*zlte+(zy-2.5)*zpi/zpe*zlti)
     +              -zl23*zlte
  220    continue
c
c..compute the trapped particle current density in the banana regime
c  ajtpbi(j) = J_{tpbi} / R  [Amps / m^3]
c
c  \right< \frac{J_{tpbi}}{R/R_o} \left>
c  = - 1.25 \frac{ (kT_e + kT_i)[ergs] n_e[cm^{-3}] \sqrt{r/R}
c      zcnvrt }{ R \overstrike{B}_{pol} [Gauss]}  [Amps / m^2]
c
c  Note: the factor zcnvrt serves to convert from the cgs units used to
c        evaluate the right hand side, to mks units for current density
c
      if ( versno .gt. 15.05 ) then
c
        zcnvrt = fcc*1.e10 * usij
c
        do 230 jz=kfirst,klast
          pjtpbi(jz) = -1.25 * (tes(2,jz) + tis(2,jz)) * rhoels(2,jz)
     &      * sqrt ( ahalfs(jz,2) / rmids(jz,2) )
     &      * zcnvrt / ( r0ref * rmids(jz,2) * bpols(2,jz) )
          zai=aimass(1,jz)
          zte=tes(1,jz)*useh
          zrmin=max(armins(jz,1),epslon)*usil
          zrmaj=armajs(jz,1)*usil
          zgyrfi=zce*zb/(zcmp*zai)
          zvthi=sqrt(2.*zckb*zti/(zcmp*zai))
          zlari=zvthi/zgyrfi
          zlari=zvthi/zgyrfi
          zrcut=sqrt(q(jz)*zlari*zrmaj)
          zlim=1.0/(1.0+sqrt(abs(sfutz(17)*zrcut/zrmin)))
            pjtpbi(jz) = zlim*pjtpbi(jz)
 230    continue
c
      endif
c
         return
         end

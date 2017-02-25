c--------1---------2---------3---------4---------5---------6---------7-c
c@ncprof0  /11040/baldur/code/bald/dimprad.f
c  rgb 31-jan-95 move call ncrats before loop over zones
c  rgb 29-dec-94 created sbrtn ncprof0 from part of sbrtn ncinit
c***********************************************************************
c
        subroutine ncprof0
c
c..initialize profiles for non-equilibrium impurity radiation package
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'cbparm.m'
      include 'comncr.m'
      include 'comadp.m'
c
      real zk(kzmax), zn(kzmax)
c
c-----------------------------------------------------------------------
c
c
      if ( mimp .le. 0 )    return
c
c
c..Solve the equilibrium ionization rate equations
c  assuming steady state with no transport, sources, or sinks
c  -S(k) n(k) + S(k-1) n(k-1) - R(k) n(k) + R(k+1) n(k+1) = 0
c  where n(k) = density of ionization state k
c  S(k) = ionization rate
c  R(k) = recombination rate
c
      zepsln = 1.e-6
c
c  loop over impurity species
c
      do ji=1,mimp
c
        nk = nkimp(ji)
c
        write (6,*) ' ji=',ji,'  nk=',nk
c
c  Compute ionization and recombination rates
c
      call ncrats (2,ji)
c
c  loop over zones
c
        do jz=2,mzones
c
          write (6,*) ' jz=',jz
c
c  Let n(k) = zk(k) n(k-1).
c  Compute zk(k) = S(k-1) / ( S(k) + R(k) - zk(k+1) R(k+1) )
c
      zk(nk) = sa(jz,nk-1,ji) / ( sa(jz,nk,ji) + ra(jz,nk,ji) )
      do jk=nk-1,1,-1
        zk(jk) = sa(jz,jk-1,ji) /
     &    max ( zepsln * sa(jz,jk-1,ji), 
     &    (sa(jz,jk,ji) + ra(jz,jk,ji) - zk(jk+1) * ra(jz,jk+1,ji)) )
      enddo
c
c..find the index where zk(jk) goes from > 1 to < 1
c
c
      iz = 1
      do jk=1,nk
        if ( zk(jk) .gt. 1.0 ) iz = jk
      enddo
c
      write (6,201) (zk(jk),jk=1,nk)
  201 format(' zk#',1p6e12.3)
c
c..compute density factors
c
      zn(iz) = 1.0
      if ( iz .lt. nk ) then
        do jk=iz+1,nk
          zn(jk) = zn(jk-1) * zk(jk)
          if ( zn(jk) .lt. zepsln ) zn(jk) = 0.0
        enddo
      endif
c
      if ( iz .gt. 1 ) then
        do jk=iz-1,1,-1
          zn(jk) = zn(jk+1) / zk(jk+1)
          if ( zn(jk) .lt. zepsln ) zn(jk) = 0.0
        enddo
      endif
c
      zsum = 0.0
      do jk=1,nk
        zsum = zsum + zn(jk)
      enddo
c
      write (6,202) (zn(jk),jk=1,nk)
  202 format(' zn#',1p6e12.3)
c
c..normalize the density of each ionization state so that the sum
c  of these densities equals the total impurity density
c
      do jk=1,nk
        xn(jz,jk,ji) = rhois(ji,2,jz) * zn(jk) / zsum
      enddo
c
c..sum over density of ionization states
c
      xns(jz,ji) = 0.0
      do jk=1,nk
        xns(jz,ji) = xns(jz,ji) + xn(jz,jk,ji)
      enddo
c
      write (6,203) jz, ji, xns(jz,ji)
  203 format (' jz=',i3,' ji=',i3,' xns(jz,ji)=',1pe13.4)
c
c  end of loop over zones
c
      enddo
c
c..even symmetry across magnetic axis
c
      do jk=1,nk
        xn(1,jk,ji) = xn(2,jk,ji)
      enddo
c
      xns(1,ji) = xns(2,ji)
c
c  Set value to be used for pedestal boundary condition. Note that
c  only the first stage pedestal value is allowed to be nonzero.
c
      do jk=0,nk
        cimped(jk,ji) = 0.0
      enddo
      cimped(1,ji) = xns(mzones,ji)
c
c..compute the mean Z and Z**2
c
      do jz=2,mzones
        zsum  = 0.0
        z2sum = 0.0
c
        do jk=1,nk
          zsum  = zsum  + real(jk) * xn(jz,jk,ji) / xns(jz,ji)
          z2sum = z2sum + real(jk)**2 * xn(jz,jk,ji) / xns(jz,ji)
        enddo
c
        cmean(ji,2,jz)  = zsum
        c2mean(ji,2,jz) = z2sum
      enddo
c
      cmean(ji,2,1)  = cmean(ji,2,2)
      c2mean(ji,2,2) = c2mean(ji,2,2)
c
c  end of loop over species
c
      enddo
c
c
      return
      end

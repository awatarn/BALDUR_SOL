c--------1---------2---------3---------4---------5---------6---------7-c
c@solve  .../baldur/code/bald/dsolver.f
c  rgb 10-aug-94 save local variables that are initializ
c       dps 24-aug-88 15.00 skip conservation check totals for impurities
c                     when using NC code
c       rgb  24-oct-86 corrected volume and surface areas for 1-1/2-D
c       fgps 2-feb-83 modified to handle more than 2 impurity species.
cdoc
c=======================================================================
c
c       -----------
c       sbrtn SOLVE   file DSOLVER
c       -----------
c
c
c      2.10    compute chi at new timestep from edge boundary condition
c               and eeee and ffff arrays.
c
c-----------------------------------------------------------------------
c
c
c      common blocks and variables modified:
c
c       chi (comsta)
c
c-----------------------------------------------------------------------
c
c
c
c       reduce has computed the coefficients eeee and ffff of the relation
c
c       eeee(j).chi(j) + ffff(j) = chi(j-1)
c
c       where the indices are zone indices.
c       we also have the boundary condition at the edge:
c
c       alpha1.(chi at r-(dr/2))
c               + (beta1 + gamma1).(chi at r+(dr/2)) = delta1
c
c       where alpha1, beta1, and gamma1 are matrices, and delta1 is a vector
c
c       we form an equation of the form:
c
c       zx.chi(ledge+1) + zy = zz.chi(ledge)
c
c       where zx and zz are matrices.
c       combining this with the relation for eeee and ffff, we have
c
c       (zz.eeee(ledge+1) - zx)chi(ledge+1) + zz.ffff(ledge+1) - zy = 0
c
c       hence:
c
c       chi(ledge+1) = - (zz.eeee(ledge+1) - zx)**-1  . (zz.ffff(ledge+1)
c               - zy)
c
c       from that, chi(j) may be found in the obvious way.
c
cend
c-----------------------------------------------------------------------
c
        subroutine solve
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'commhd.m'
        include 'cbparm.m'
c
c
        dimension
     i   ipivot(idxchi)     ,
     r   zchi0(idxchi)      , zchi1(idxchi)      , zgradz(idxchi)     ,
     r   zt(id2chi)         , zu(idxchi)         , zx(id2chi)         ,
     r   zy(idxchi)         , zz(id2chi)
c
        logical         inital
c
        data            inital /.true./
c
        save inital, ix, ix2
c
c-----------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /10/
c
c
        if (.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.10 subroutine solve bypassed')
        return
   10   continue
c
c
        if(.not.inital) go to 20
          inital=.false.
          ix=mxchi
          ix2=ix*ix
   20   continue
c
c      1)      clear local variables
c
c
        call reseti(ipivot,ix,0)
        call resetr(zt,ix2,0.0)
        call resetr(zu,ix,0.0)
        call resetr(zx,ix2,0.0)
        call resetr(zy,ix,0.0)
        call resetr(zz,ix2,0.0)
c
c
c
c      2)      difference boundary condition
c
c
        z0 = 1.0 / (dxboui(ledge+1) * rmini)
c
        do 218 j1 = 1, mchi
        do 216 j2 = 1, mchi
        i002 = j1 + mxchi*(j2-1)
        zz(i002) = - alpha1(j1,j2)
        zx(i002) = gamma1(j1,j2) + beta1(j1,j2)
  216   continue
        zy(j1) = - delta1(j1)
  218   continue
c
c
c
c      3)      compute zt, zu
c
c       zt = (zz.eeee(ledge+1) - zx)
c       zu = zy - zz.ffff(ledge+1)
c
c
        do 312 j1 = 1, mchi
        z0 = 0.0
c
        do 308 j2 = 1, mchi
        z1 = 0.0
c
        do 304 j3 = 1, mchi
        i003 = j1 + mxchi*(j3-1)
        z1 = z1 + zz(i003) * eeee(j3,j2,ledge+1)
  304   continue
c
        i002 = j1 + mxchi*(j2-1)
        zt(i002) = z1 - zx(i002)
        z0 = z0 + zz(i002) * ffff(j2,ledge+1)
  308   continue
c
        zu(j1) = zy(j1) - z0
  312   continue
c
c
c      3.2)    solve matrix equation
c
c
        call matrx1(zt,mxchi,mchi,ipivot,ierror)
        if(ierror.ne.0) go to 9030
c
        call copyr(chi(1,ledge+1),1,zchi1,1,mchi)
c
        call matrx2(chi(1,ledge+1),zu,zt,ipivot,mxchi,mchi,1,1,1)
c
        do 324 jp = 1, mchi
        zchi1(jp) = zchi1(jp)*(1.0 - thetai) + chi(jp,ledge+1)*thetai
  324   continue
c
c
c
c
c      4)      compute remaining chi, cumulative fluxes, sources,...
c
c
        call resetr(fflxii,ix,0.0)
        call resetr(fflxoi,ix,0.0)
        call resetr(fsorci,ix,0.0)
        call resetr(fsordi,ix,0.0)
c
        i1 = ledge+1
        i2 = i1 + lcentr
c
        do 499 j0 = lcentr, i1
        i0 = i2 - j0
c
c  zdvoli = differential volume between zone bndry i0 and i0-1
c
        zdvoli = avi(i0-1,4,1) * (xbouni(i0)-xbouni(i0-1))
c
        do 422 j1 = 1, mchi
        z0 = 0.0
c
        do 404 j2 = 1, mchi
        z0 = z0 + eeee(j1,j2,i0)*chi(j2,i0)
  404   continue
c
c               save old timestep value of chi
c
        z1 = chi(j1,i0-1)
c
c               compute new timestep value of chi
c
        chi(j1,i0-1) = z0 + ffff(j1,i0)
c
c               compute theta-time-centered value of chi in this zone
c
        zchi0(j1) = z1 * (1.0 - thetai) + chi(j1,i0-1)*thetai
c
c               dddd contributions to plasma
c
        fsordi(j1) = fsordi(j1) + dddd(j1,i0-1) * zdvoli
c
c               cccc contributions to plasma
c
        do 408 jp = 1, mchi
        fsorci(jp) = fsorci(jp) + cccc(jp,j1,i0-1) * zchi0(j1)*zdvoli
  408   continue
c
  422   continue
c
        if (i0.ne.lcentr + 1) go to 430
        call copyr(zchi0,1,zchi1,1,ix)
        go to 499
c
  430   continue
c
c
c               aaaa, bbbb (flux) contributions to plasma
c
c
        if (i0.ne.(ledge + 1).and.i0.ne.lcentr) go to 499
c
c  zsurfi = surface area of zone boundary i0 = V'(xi) * <|del xi|>
c
        zsurfi = avi(i0,3,1) * avi(i0,5,1)
c
c
        do 438 jp = 1, mchi
c
c               compute  flux at bound. i0 --
c               = - (aaaa.(chi-) + bbbb.(chi+))
c
 
        z1 = 0.0
c
        do 434 jp2 = 1, mchi
        z1 = z1 - (aaaa(jp,jp2,i0) * zchi0(jp2) +
     1          bbbb(jp,jp2,i0) * zchi1(jp2)) * zsurfi
  434   continue
c
        if (i0.eq.lcentr) go to 436
c
c               positive flux at outer boundary removes
c               particles/energy
c
        fflxoi(jp) = fflxoi(jp) - z1
        go to 438
c
  436   continue
c
c               positive flux at inner boundary adds
c               particles/energy
c
        fflxii(jp) = fflxii(jp) + z1
  438   continue
  499   continue
c
c
c
c
c      5)      conservation checks
c
c
  500   continue
        call resetr(totali,ix,0.0)
c
        do 508 jz = lcentr, ledge
        zdvoli = avi(jz,4,1) * (xbouni(jz+1) - xbouni(jz))
c
        do 504 jp = 1, mchi
        totali(jp) = totali(jp) + chi(jp,jz) * zdvoli
  504   continue
c
  508   continue
c
c
c       add rates*dti
c  15.00 Skip calculation for impurities when using NC code
c
c
        do 512 jp = 1, mchi
        if ((natomc.eq.3).and.(jp.ge.limp1)
     1                   .and.(jp.le.limpn)) go to 512
        aflxii(jp) = aflxii(jp) + fflxii(jp)*dti
        aflxoi(jp) = aflxoi(jp) + fflxoi(jp)*dti
        asorci(jp) = asorci(jp) + fsorci(jp)*dti
        asordi(jp) = asordi(jp) + fsordi(jp)*dti
  512   continue
c
c
        do 522 jp = 1, limpn
        if ((natomc.eq.3).and.(jp.ge.limp1)
     1                   .and.(jp.le.limpn)) go to 522
        ccons(jp) = (totali(jp) - aflxii(jp) - aflxoi(jp) - asorci(jp)
     1          - asordi(jp) - acompi(jp)) / begini(jp) - 1.0
  522   continue
c
        ccons(lelec) = (totali(lelec) + totali(lion) -
     1  aflxii(lelec) - aflxii(lion) - aflxoi(lelec) - aflxoi(lion) -
     2  asorci(lelec) - asorci(lion) - asordi(lelec) - asordi(lion) -
     3  acompi(lelec) - acompi(lion)) / (begini(lelec) + begini(lion))
     4  - 1.0
c
        ccons(lion) = ccons(lelec)
c
c
c
        return
c
c
c      90)     if boundary condition does not solve
c
c
 9030   continue
c
        call error_olymp(0,iclass,isub,3,
     1          'boundary condition problem at edge in solve')
        return
        end

c--------1---------2---------3---------4---------5---------6---------7-c
c@reduce  .../baldur/code/bald/dsolver.f
c  rgb 10-aug-94 save local variables that are initializ
c       dps  17-aug-88 14.04 use factor zdt instead of dti to multiply
c                      dddd in zs.
c       fgps 2-feb-83 modified to handle more than 2 impurity species.
cdoc
c=======================================================================
c
c       ------------
c       sbrtn REDUCE   file DSOLVER
c       ------------
c
c
c      2.9)    calculate eeee and ffff arrays for reducing equations
c               to first order, given aaaa, bbbb, cccc, and dddd
c
c
c-----------------------------------------------------------------------
c
c
c      common blocks and variables modified:
c
c       eeee and ffff (comcoe)
c
c-----------------------------------------------------------------------
c
c       the coefficients that have been computed so far are
c       aaaa, bbbb, cccc, and dddd, which are for the equation
c
c       (d/dt)chi = (1/r)(d/dr)(r * aaaa . (chi at r-(dr/2))
c                       + r*bbbb . (chi at r+(dr/2)) )
c                       + cccc . chi + dddd
c
c       where "." indicates matrix multiplication.
c       aaaa and bbbb are boundary centered matrices (remember,
c       boundary j is the *inside* boundary of zone j.)
c       cccc is a zone-centered matrix, and dddd is a zone centered vector.
c
c       this routine first converts aaaa, bbbb, cccc, and dddd to zppp,
c       zqqq, zrrr, and zs, to fit the equation
c
c       zppp(i).chi(i+1) + zqqq(i).chi(i) + zrrr(i).chi(i-1)
c               + zs(i) = 0
c       zs(i) includes the terms involving old timestep values of chi.
c       the "chi" are the new timestep values of chi, the unknows.
c
c       next, these coefficients are transformed into the coefficients
c       eeee and ffff (matrix and vector, respectively) of the
c       equation
c
c       chi(i-1) = eeee(i).chi(i) + ffff(i)
c
c       the eeee and ffff are determined by the requirement that,
c       the equation in eeee and ffff be implied by the equation in
c       zppp, zqqq, zrrr, and zs, and that it be equivalent to the
c       boundary condition at the center.
c       the first requirement implies that, given eeee and ffff at
c       zone j-1,
c
c       eeee(j) =  - (zrrr(j-1).eeee(j-1) + zqqq(j-1))**-1 .zppp(j-1)
c
c       ffff(j) =  - (zrrr(j-1).eeee(j-1) + zqqq(j-1))**-1
c                               .(zs(j-1) + zrrr(j-1).ffff(j-1))
c
c       by also making it equivalent to the boundary condition at the
c       center, eeee and ffff at the center are determined, and
c       hence may be found for all zones.
c
c       remember, the boundary condition at the magnetic axis is
c
c       (alpha0 + gamma0).(chi at r-(dr/2))
c               + beta0.(chi at r+(dr/2)) = delta0
c
c       where alpha0, beta0, and gamma0 are matrices, and delta0 is a
c       vector.
c
c       subroutine "solve" then computes the new timestep values of chi,
c       (which replace the old values)
cend
c**********************************************************************c
c
        subroutine reduce
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'commhd.m'
        include 'cbparm.m'
c
c
        dimension
     i   ipivot(idxchi)     ,
     r   zppp(id2chi)       , zqqq(id2chi)       , zrrr(id2chi)       ,
     r   zs(idxchi)         , zt(id2chi)         , zu(idxchi)
c
        logical         inital
c
        data            inital /.true./
c
        save inital, ix, ix2, i1, i2
c
c-----------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /9/
c
c
        if (.not.nlomt2(isub)) go to 10
        call mesage(' *** 2.9 subroutine reduce bypassed')
        return
   10   continue
c
c
        if(.not.inital) go to 20
          inital=.false.
          ix=mxchi
          ix2=ix*ix
          i1=mxchi*mxzone
          i2=i1*mxchi
   20   continue
c
c      1)      clear local variables
c
c
        call reseti(ipivot,ix,0)
        call resetr(zppp,ix2,0.0)
        call resetr(zqqq,ix2,0.0)
        call resetr(zrrr,ix2,0.0)
        call resetr(zs,ix,0.0)
        call resetr(zt,ix2,0.0)
        call resetr(zu,ix,0.0)
c
        if (.not.nldiff) call resetr(aaaa,i2,0.0)
        if (.not.nldiff) call resetr(bbbb,i2,0.0)
        if (.not.nlsorc) call resetr(cccc,i2,0.0)
        if (.not.nlsord) call resetr(dddd,i1,0.0)
c
c
c
c      2)      use boundary conditions at center to find
c               eeee(lcentr) and ffff(lcentr).  note: "zqqq(*)"
c               is really "-q(1)" in the equation  -q(1)*e(2)=p(1).
c
c
        do 210 j1 = 1, mchi
        do 210 j2 = 1, mchi
        i001 = j1 + mxchi*(j2-1)
        zppp(i001) = gamma0(j1,j2) + alpha0(j1,j2)
        zqqq(i001) = - beta0(j1,j2)
  210   continue
c
        i=0
        call matrx1(zppp,mxchi,mchi,ipivot,i)
        if(i.ne.0) go to 9020
c
        call matrx2(eeee(1,1,2),zqqq,zppp,ipivot,mxchi,mchi,1,mchi,1)
        call matrx2(ffff(1,2),delta0,zppp,ipivot,mxchi,mchi,1,1,1)
c
        call resetr(zppp,ix2,0.0)
        call resetr(zqqq,ix2,0.0)
c
c
c
c      3)      compute remaining  eeee's and  ffff's  by
c               computing  zppp, zqqq, zrrr, and zs at
c               each zone, and thence eeee and ffff of the
c               next zone.
c
c
        do 380 j0 = lcentr, ledge
c
c bateman  conservative finite difference 1 1/2 D transport equations
c
c  a finite difference form of the following equation is implemented:
c
c     d ( V' chi(j1) ) / dt  = V' dddd(j1)  +  sum over j2 of
c     ( d / dxi) { V' <|del xi|> [ aaaa(j1,j2) ( chi(j2) at r-dr/2 )
c                                 + bbbb(j1,j2) (chi(j2) at r+dr/2)]}
c     + cccc(j1,j2) * chi(j2)
c
c  the following local variables are used:
c
c  zdhv   = V' (step N+1, zone j0) - V' (step N+1/2, zone j0)
c  zvpnew = V' (step N+1, zone center j0)
c  zratio = V' (step N, zone j0) / V' (step N+1, zone j0)
c  zdt    = dti * V' ( N+1/2, zone j0) / V' (N+1, zone j0)
c  zdtdxi = dti / dxzoni (j0)
c  zdel   = zdtdxi * V' (N+1/2, boundary j0) * <|del xi|> (N+1/2, bound j0)
c        / [ V' (N+1, zone j0) ]
c  zdelp1 = same as zdel but at j0 + 1
c
      zdhv   = avi(j0,4,2) * 0.5 * (tbi-tai)
      zvpnew = avi(j0,4,1) + zdhv
      zratio = (avi(j0,4,1) - zdhv) / (avi(j0,4,1) + zdhv)
      zdt    = dti * avi(j0,4,1) / zvpnew
      zdtdxi = dti / (xbouni(j0+1) - xbouni(j0))
      zdel   = zdtdxi * avi(j0,3,1) * avi(j0,5,1)
     &            / (zvpnew)
      zdelp1 = zdtdxi * avi(j0+1,3,1) * avi(j0+1,5,1)
     &            / (zvpnew)
c
c  14.04 Replace dti multiplying dddd with zdt to fix bug.
c
      if (versno.gt.14.03) then
        zzdt = zdt
      else
        zzdt = dti
      end if
c
c
c      3.1)    compute zppp, zqqq, zrrr, and zs from aaaa, bbbb, cccc,
c                       and dddd
c
c
        do 312 j1 = 1, mchi
        i001 = j1 + mxchi*(j1-1)
        z7 = 0.0
c
        do 310 j2 = 1, mchi
        i002 = j1 + mxchi*(j2-1)
        z4 = zdelp1 * bbbb(j1,j2,j0+1)
        z5 = - zdel * aaaa(j1,j2,j0)
        z6 = zdt * cccc(j1,j2,j0) + zdelp1 * aaaa(j1,j2,j0+1)
     &  - zdel * bbbb(j1,j2,j0)
        z7 = z7 + chi(j2,j0+1)*z4 + chi(j2,j0)*z6 + chi(j2,j0-1)*z5
c
        zppp(i002) = z4*thetai
        zqqq(i002) = z6*thetai
        zrrr(i002) = z5*thetai
  310   continue
c
        zqqq(j1 + mxchi*(j1-1)) = zqqq(j1 + mxchi*(j1-1)) - 1.0
        zs(j1) = z7*(1.-thetai) + zratio*chi(j1,j0)
     1            + zzdt*dddd(j1,j0)
  312   continue
c
cend bateman 6-apr-85
c
c
c      3.2)    compute zt = -(zrrr(j0).eeee(j0) + zqqq(j0)) and
c               zu = (zs(j0) + zrrr(j0).ffff(j0)), so that
c               eeee(j0+1) = zt**-1.zppp(j0) and
c               ffff(j0+1) = zt**-1.zu
c
        do 332 j1 = 1, mchi
                z2 = 0.0
c
        do 328 j2 = 1, mchi
                z3 = 0.0
c
        do 324 j3 = 1, mchi
        i003 = j1 + mxchi*(j3-1)
        z3 = z3 - zrrr(i003) * eeee(j3,j2,j0)
  324   continue
c
        i002 = j1 + mxchi*(j2-1)
        zt(i002) = z3 - zqqq(i002)
        z2 = z2 + zrrr(i002)*ffff(j2,j0)
  328   continue
c
        zu(j1) = zs(j1) + z2
  332   continue
c
c
c      3.3)    compute zt**-1 . zppp(j0)  and zt**-1 . zu
c
c
        i=0
        call matrx1(zt,mxchi,mchi,ipivot,i)
        if(i.ne.0) go to 9030
c
        call matrx2(eeee(1,1,j0+1),zppp,zt,ipivot,mxchi,mchi,1,mchi,1)
        call matrx2(ffff(1,j0+1),zu,zt,ipivot,mxchi,mchi,1,1,1)
c
  380   continue
c
        return
c
c      if boundary condition has problems
c
 9020   continue
c
        call error_olymp(1,iclass,isub,2,
     1          'boundary condition problem at center')
        call error_olymp(2,z0,1,1,'  1/dr  ')
        return
c
c
c      singular transport matrix
c
c
 9030   continue
c
        call error_olymp(1,iclass,isub,3,
     1          'singular matrix in reducing transport equation')
        call error_olymp(3,j0,2,1,'  zone  ')
        call error_olymp(3,zppp,5,ix2,'p-matrix')
        call error_olymp(3,zqqq,5,ix2,'q-matrix')
        call error_olymp(2,zrrr,5,ix2,'r-matrix')
        return
        end

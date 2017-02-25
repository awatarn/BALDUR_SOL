c--------1---------2---------3---------4---------5---------6---------7-c
c@inital  .../11040/baldur/code/bald/inital.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 30-jan-95 iterate only if natomc = 3
c  rgb 29-dec-94 change call noncor to call ncinit and call nprof0
c  rgb 20-dec-94 iterate over impurity levels to get n_e to converge
c      les  nov-90  add d3he fast particles to electron density
c       dps 08-jun-89 15.09 Add cmean, c2mean settings at mzones
c       dps 15-may-89 15.09 Add initial call to NONCOR
c       fgps 9-oct-78 moved ahead the definitions of bxs, rmajs,
c                     and rmins; the latter is used in imprad
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 21-apr-77 set rhoins
c       amck 15-mar-77 call imprad -> call imprad(2)
c       amck 25-may-76 fix conversion of bpols -> bpoli
c       amck 25-may-76 add "set rmins..."
c       amck 21-may-76 "imprad(1)" -> "imprad"
c       amck 20-may-76 change -ai, -bi to -(i), add setting bpoli
c       amck 15-apr-76
c
        subroutine inital
c
c**********************************************************************c
c
cl      1.6     define physical initial conditions
c
       include 'cparm.m'
      include 'cbaldr.m'
       include 'cd3he.m'
!ap
      integer :: itmax = 1.
c
c------------------------------------------------------------------------------
c
c
        data    iclass /1/,     isub /6/
c
      save iatomc
c
        if (.not.nlomt1(isub))  go to 10
        call mesage(' *** 1.6 subroutine inital bypassed ')
        return
   10   continue
c
c
c------------------------------------------------------------------------------
c
cl      4)      set rmins, rmajs, bzs
c
        bzs = bzi * uisb
        rmajs = rmaji * uisl
        rmins = rmini * uisl
c
c
cl      1)      get impurity levels
c
c..iterate over impurity levels to get n_e to converge
c  use the equilibrium impurity radiation model on the first pass
c
      if ( natomc .eq. 3 ) call ncinit
c
c..iterate only if natomc = 3
c
      itmax = 1
      if ( natomc .eq. 3 ) itmax = 3
c
      iatomc = natomc
      if ( iatomc .gt. 2 ) natomc = 2
c
      do jk=1,itmax
c
        call imprad(2)
c
c..Initialize Impurity Rate Equations code profiles
c
        if ( natomc .eq. 3 ) call ncprof0
c
c  reset natomc after the first pass
c
        natomc = iatomc
c
c..15.09 Reset cmean, c2mean at mzones. Used to be in IMPRAD. Moved
c        to GETCHI. Duplicated here for consistency during initialization.
c
      if ( limprd(32) .eq. 0 ) then
c
c..reset cmean=1. and c2mean=1. as in original 1-D BALDUR
c
        call resetr(cmean(1,2,mzones),mimp,1.0)
        call resetr(c2mean(1,2,mzones),mimp,1.0)
c
c..linearly extrapolate cmean and c2mean to zone center outside plasma
c
      elseif ( limprd(32) .eq. 1 ) then
        zint1 = (xzoni(mzones  ) - xzoni(mzones-2))
     &        / (xzoni(mzones-1) - xzoni(mzones-2))
        zint2 = (xzoni(mzones-1) - xzoni(mzones  ))
     &        / (xzoni(mzones-1) - xzoni(mzones-2))
c
        do 201 ji=1,mimp
          cmean(ji,2,mzones) = max ( 1. ,
     &  zint1 * cmean(ji,2,mzones-1) + zint2 * cmean(ji,2,mzones-2) )
          c2mean(ji,2,mzones) = max ( 1. ,
     &  zint1 * c2mean(ji,2,mzones-1) + zint2 * c2mean(ji,2,mzones-2) )
 201  continue
c
      endif
                                        call expert(iclass,isub,1)
c
c
cl      2)      set rhoels, rhoins
c
c
        do 216 j1 = 1, mzones
        z1 = 0.0
c
        do 204 j2 = 1, mhyd
        z1 = z1 + rhohs(j2,2,j1)
  204   continue
        z0 = z1
c
c
        if (mimp.le.0) go to 210
        do 208 j2 = 1, mimp
        z0 = z0 + rhois(j2,2,j1)
        z1 = z1 + rhois(j2,2,j1) * cmean(j2,2,j1)
  208   continue
  210   continue
c
c   les  nov-90  d3he: add fast fusion particles to electrons
c
      if ( cfutz(490) .gt. epslon ) z1=z1+
     1 fnpd(j1)+fnp3(j1)+fnt(j1)+2.*(fn3(j1)+fn4t(j1)+
     2 fn43(j1)+fn42t(j1))
c
        rhoins(2,j1) = z0
        rhoels(2,j1) = z1
  216   continue
c
c..end of iteration over impurity levels
c
      enddo
c
cl      3)      set chi and bpoli from tes, tis, rhoels, rhohs, rhois,
c               bpols, and rhoins
c
c
        z0 = 1.0 / gamin1
        do 316 j1 = 1, mzones
c
        do 304 j2 = lhyd1, lhydn
        chi(j2,j1) = rhohs(j2,2,j1) * usid
  304   continue
c
        if (mimp.le.0) go to 310
        do 308 j2 = limp1, limpn
        i001 = j2 + 1 - limp1
        chi(j2,j1) = rhois(i001,2,j1) * usid
  308   continue
c
  310   continue
c
        chi(lelec,j1) = (tes(2,j1) * rhoels(2,j1)) * usih * usid
        chi(lion,j1) = (tis(2,j1) * rhoins(2,j1)) * usih * usid
c
        bpoli(j1) = bpols(1,j1) * usib
  316   continue
c
      bpoli(mzones+1) = bpols(1,mzones+1) * usib
c
                                        call expert(iclass,isub,3)
c
        return
        end

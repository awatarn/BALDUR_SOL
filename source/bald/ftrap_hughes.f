c@ftrap_hughes   .../baldur/bald/ftrap_hughes.f
c  rgb 18-jul-01 move trapf to ftrap_hughes
c  dps 11-may-87 rework interpolation onto BALDUR grids
c  dps 07-may-87 remove flux-surface-average integration weights;
c                divide Btheta by 2pi; use rbtorb, not rbtorz;
c                specify no symmetry in calls to CUBINT.
c  dps 04-mar-87 added trapf to improve trapped particle fraction
c                calculation to go with bootstrap current cf. Mike Hughes.
c
         subroutine ftrap_hughes
c
C.  Calculate trapped particle fraction implemented by Mike Hughes
c
       include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
       real :: zwork(55), zwint(55)
       real :: zbmax, zbsqav, zjacav, zlim, zdl, zint, zl
       integer :: jz, jt, jt1
c
c---------------------------------------------------------------------
       dimension ztrsq(43)
       data        int/20/
c
c---------------------------------------------------------------------
CL              1.         Scan Flux Surfaces
c
c integration weights
c
      zwint = 1.0
      zwint(1) = 0.5
      zwint(ntheta) = 0.5
c
      zwork = 0.0
c
      do 300 jz=3,mjbal
c
        zbmax  = 0.0
        zbsqav = 0.0
        zjacav = 0.0        
c
        do jt=1,ntheta
c
c  zwork(jt) = | B | on zone boundary jz at angle theta(jt)
c
          zwork(jt) = sqrt ( bpol2di(jz,jt)**2 + btor2di(jz,jt)**2 )
c
          zbsqav = zbsqav + ejacob2d(jz,jt) *
     &                     ( bpol2di(jz,jt)**2 + btor2di(jz,jt)**2 )
          zjacav = zjacav + ejacob2d(jz,jt)
c
          zbmax = max ( zbmax, zwork(jt) )
c
        enddo
c
        zbsqav = zbsqav / zjacav
c
c---------------------------------------------------------------------
CL              2.         Integrate
c
C     Upper limit of integration
c
         zlim = 1.0 / zbmax
         zdl  = zlim / float(ntheta-1)
c
         zint = 0.0
         zl   = 0.0
c
         do 202 jt=1,ntheta
c
C     Calculate surface average
c
         zden = 0.0
         do jt1=1,ntheta
           zarg = max ( 1.0e-10, (1.0-zl*zwork(jt1)) )
           zden = zden + sqrt( zarg ) * ejacob2d(jz,jt1)
         enddo
c
         zden = zden / zjacav
         zint = zint + zwint(jt) * zdl * zl / zden
  202    zl   = zl + zdl
c
c---------------------------------------------------------------------
CL              3.         Trapped Particle Fraction
c
           trapbr(jz,1) = 1.0 - 0.75 * zbsqav * zint
c
  300    continue
c
c     Magnetic axis at jz=2
c
      trapbr(2,1) = 0.0
      trapbr(1,1) = - trapbr(3,1)
c
c  Edge of plasma at mjbal
c
      trapbr(mjbal+1,1) = trapbr(mjbal,1)
c
c---------------------------------------------------------------------
CL              4.         Interpolate to Transport Grid
c
c...Interpolate f-trap ** 2 since f-trap goes like sqrt(r); can
c...then use odd symmetry.
c
         zwork = trapbr(:,1)**2
         zwork(1) = - zwork(3)
c
      do jz=1,mjbal-1
c
        zint = ( xzoni(jz) - xbouni(jz) )
     &         / ( xbouni(jz+1) - xbouni(jz) )
c
        zwork(jz) = zwork(jz) * ( 1.0 - zint ) + zwork(jz+1) * zint
c
        trapbr(jz,2) = sqrt( max( 1.e-10, zwork(jz) ) )
c
      enddo
c
         return
         end

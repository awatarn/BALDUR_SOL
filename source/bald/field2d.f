c@field2d   .../baldur/bald/field2d.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
c  Compute magnetic field components as a function of 2-D grid
c    on BALDUR zone boundaries and poloidal angle
c
c  The following 2-D arrays are in common /comb2d/ in file commhd
c
c bpol2di(jm,jz) = poloidal field
c    at BALDUR zone boundary jz and angle theta(jm) (T)
c btor2di(jm,jz) = toroidal field
c    at BALDUR zone boundary jz and angle theta(jm) (T)
c
c b22di(jz)      = <B^2>
c bm22di(jz)     = <1./B^22>
c delxi27b2i(jz) = < | del xi |^2 / B^2 >
c bdotdeltheta(jz) = < B \dot del theta / |B| >
c
      subroutine field2d
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
      integer :: jz, jm
c
      real :: zjacob  ! poloidal integral of the Jacobian
c
      bpol2di      = 0.0
      btor2di      = 0.0
      b22di        = 0.0
      bm22di       = 0.0
      delxi27b2i   = 0.0
      bdotdeltheta = 0.0
c
      do jz = 1, mzones
c
        zjacob = 0.0
c
        do jm = 1, ntheta
c
          bpol2di(jz,jm) = r0ref * bpoli(jz) * avi(jz,2,1)
     &      * delxi2d(jz,jm) / max ( rthxbi(jm,jz), epslon )
c
          btor2di(jz,jm) = rbtors(jz,1) * usil * usib
     &      / max ( rthxbi(jm,jz), epslon )
c
c..flux surface averages
c
          b22di(jz) = b22di(jz) + ejacob2d(jz,jm)
     &      * ( bpol2di(jz,jm)**2 + btor2di(jz,jm)**2 )
          bm22di(jz) = bm22di(jz) + ejacob2d(jz,jm)
     &      / ( bpol2di(jz,jm)**2 + btor2di(jz,jm)**2 )
c
          delxi27b2i(jz) = delxi27b2i(jz) + ejacob2d(jz,jm)
     &      * delxi2d(jz,jm)
     &      / ( bpol2di(jz,jm)**2 + btor2di(jz,jm)**2 )
c
          bdotdeltheta(jz) = bdotdeltheta(jz) + ejacob2d(jz,jm)
     &      *  bpol2di(jz,jm) * sqrt ( delth2d(jz,jm)
     &      / ( bpol2di(jz,jm)**2 + btor2di(jz,jm)**2 ) )
c
          zjacob = zjacob + ejacob2d(jz,jm)
c
        enddo
c
c  make sure the average Jacobian is bounded away from zero
c
        zjacob = max ( zjacob, epslon )
c
        b22di(jz) = b22di(jz) / zjacob
        bm22di(jz) = bm22di(jz) / zjacob
        delxi27b2i(jz) = delxi27b2i(jz) / zjacob
        bdotdeltheta(jz) = bdotdeltheta(jz) / zjacob
c
      enddo
c
      return
      end


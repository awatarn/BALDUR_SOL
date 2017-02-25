      subroutine weiland18init(shear, kappa, QQ, ekyrhoin, tauhin,
     &    gne, gnh,gth,
     &    ftrapein,IK, WZJ, WIMAX, Xbest)
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      IMPLICIT none
      REAL     tauhin, gne, gnh, gth,ftrapein
      REAL     shear, kappa, QQ, ekyrhoin, 
     &         zflh, FTRT, zepsmin, epsnhin,epsnin, etai, 
     &         BTA, GAV,  WIMAX(*),
     &         SH2
      INTEGER  IK
      COMPLEX  IU, ALPHA, WZJ(100), E1, EE, WZ1, WZ2, Xbest(*)
C=======================================================================
C Dictionary
C  INPUT
C     shear   =
C     kappa   =
C     QQ      =
C     ekyrhoin=
C     tauhin  =
C     gnh     =
C     ftrapein=
C     IK      =
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C INTERMEDIATE VARIABLES
C     
      IU     =  (0.,1.)     
      zflh   =  ekyrhoin**2
      SH2    =  2.*shear-1.+(KAPPA*(shear-1.))**2
      ALPHA  = -IU*ABS(SQRT(SH2))*QQ*zflh
      FTRT   =  5./3./tauhin
      zepsmin=  1.e-10
      epsnhin=  2./sign(max(abs(gnh),zepsmin),gnh)
      epsnin =  2./sign(max(abs(gne),zepsmin),gne)
      etai   =  gth/sign(max(abs(gnh),zepsmin),gnh)
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C OUTPUT
C
C     WZJ(IK)=  Current value of WZ - desired complex mode
C     Xbest(1)= Current "best" value for desired mode
C=======================================================================
      BTA    =  FTRT*(1./(1.-ftrapein)+1./tauhin)
      GAV    =  1.
      E1     =  FTRT*(1.+zflh)-(1./epsnhin)+zflh/tauhin/epsnhin*
     &              (1.+etai)+(1./(1.-ftrapein)+FTRT)*GAV
     &                +IU*0.5*ABS(SQRT(SH2))/QQ*(1.+FTRT)
      E1     =  0.5*E1/(1.+zflh)
      EE     = (1./tauhin/epsnhin/(1.-ftrapein)*(etai-2./3.)+BTA)*
     &                (GAV+IU*0.5*ABS(SQRT(SH2))/QQ)
     &                -FTRT/epsnhin*(1.-zflh/tauhin*(1.+etai))
      EE     =  EE/(1.+zflh)
      WZ1    = -E1+SQRT(E1*E1-EE)
      WZ2    = -E1-SQRT(E1*E1-EE)
      WZJ(IK)=  WZ1
      IF(AIMAG(WZ2).GT.AIMAG(WZ1))WZJ(IK)=WZ2
      Xbest(1)= WZJ(IK)
      WIMAX(IK)=  MAX(.001,.2/(epsnin*ekyrhoin))
C     write(*,'(A,I2,A,1F12.4)')'in Init WIMAX(',IK,')=',WIMAX(IK)   
C     write(*,'(A,2F12.4)')'in Init Xbest(1)=',Xbest(1)
      RETURN
      END

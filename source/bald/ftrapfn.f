c--------1---------2---------3---------4---------5---------6---------7-c
c@ftrapfn  .../baldur/bald/ftrapfn.f
c  rgb 18-jul-01 trapped particle fraction based on min-max of |B|
c    when kopt different from 0 or 1
c  rgb 18-jul-01 complete rewrite, changed from function to subroutine
c    now computes trapbr(jz,js) over the entire range of jz, js=1,2
c  rgb 11-jul-01 changed ft.f to ftrapfn.f
c  rgb 11-jul-01 removed comtrp.m and inserted cbaldr.m
c  rgb 26-jul-88  removed cliches cl1 and clintf
c  dps 04-mar-87 inserted with bootstrap current cf. Mike Hughes
c
      subroutine ftrapfn( kopt )
c
c. computation of trapped particle fraction trapbr(jz,js)
c
c KOPT = 0 for the original trapped particle fraction
c      = 1 for Mike Hughe's computation in sbrtn ftrap_hughes
c
      include 'cparm.m'
      include 'cbaldr.m'
      include 'commhd.m'
c
      integer jz, js, jm
c
c  jz = BALUDR zone number
c  js = 1 for BALDUR zone boundary
c     = 2 for BALDUR zone center
c  jm = poloidal angle index
c
      real zb2min, zb2max, zb2
c
C---------------------------------------------------------------------
CL              0.  Original BALDUR Calculation
C
      if ( kopt .eq. 0 ) then
c
        do js = 1, 2
          do jz = 1, mzones
c
           zd = abs( ahalfs( jz, js) / rmids( jz, js) )
           trapbr(jz,js) = 1.0 - ( (1.0-zd)**2 / (sqrt(1.0-zd*zd) *
     &                    (1.0+1.46*sqrt(zd)) ) )
c
          enddo
        enddo
c
        return
c
c-----------------------------------------------------------------------
CL              1.  Calculation by Mike Hughes
c
      elseif ( kopt .eq. 1 ) then
c
c..sbrtn ftrap_hughes computes trapbr(jz,js)
c
        call ftrap_hughes
c
        return
c
c-----------------------------------------------------------------------
CL              1.  Calculation by based on min-max of |B|
c
      else
c
c..Trapped particle fraction = sqrt( 1.0 - | B_min / B_max | )
c
        do jz = 3, mzones
c
          zb2min = btor2di(jz,1)**2 + bpol2di(jz,1)**2
          zb2max = zb2min
c
          do jm=2,ntheta
c
            zb2 = btor2di(jz,jm)**2 + bpol2di(jz,jm)**2
            if ( zb2 .gt. zb2max ) zb2max = zb2
            if ( zb2 .lt. zb2min ) zb2min = zb2
c
          enddo
c
          if ( zb2max .gt. zb2min)
     &       trapbr(jz,1) = 1.0 - sqrt ( zb2min / zb2max )
c
c  Note:  at this point, trapbr(jz,1) is the square of the trapped
c    particle fraction.  Interpolate that from zone boundaries to
c    zone centers before taking the square root.
c
        enddo
c
        trapbr(2,1) = 0.0
        trapbr(1,1) = - trapbr(3,1)
c
c..interpolate from zone boundaries to zone centers
c
        do jz=1,mzones-1
c
          zint =  ( xzoni(jz) - xbouni(jz) )
     &            / ( xbouni(jz+1) - xbouni(jz) )
c
          trapbr(jz,2) = trapbr(jz,1) * ( 1.0 - zint )
     &                   + trapbr(jz+1,1) * zint
c
        enddo
c
        trapbr(mzones,2) = trapbr(mzones-1,2)
c
c..Now take the sqrt to find the trapped particle fraction
c
        do jz=3,mzones
          trapbr(jz,1) = sqrt( abs( trapbr(jz,1) ) )
        enddo
        trapbr(2,1) = 0.0
        trapbr(1,1) = trapbr(3,1)
c
        do jz=2,mzones-1
          trapbr(jz,2) = sqrt( abs( trapbr(jz,2) ) )
        enddo
        trapbr(1,2) = trapbr(2,2)
c
      endif
c
      return
      end

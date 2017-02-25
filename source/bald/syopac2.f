c@syopac2  /11040/baldur/code/bald/dimprad.f
      subroutine syopac2 (nzons,tt,bb,nn)
      include 'cparm.m'
      parameter(jx=mj)
      common/alsync/  wfreq,alphae(mj),alphao(mj)
      dimension tt(mj),bb(mj)
c
c        this routine computes the direction averaged dimensionless
c        absorption coefficients in each zone from an alalytic fit.
c        this fit is valid for temperatures from about 10 to 120 kev
c        and frequencies > 3*omega-gyro.
c
c         note: original Tamor paper used denom tden=T+25/T
c           rather than (Tamor's) code version T+(25/T+4)
c
      do 200 nz=1,nzons
      wloc=wfreq/(1.76e+10*bb(nz))
      tden=1./(tt(nz)+25./(4.0+tt(nz)))
      arg1=max(0.0,.045+(wloc-2.)*tden)
      arg2=max(0.0,0.18+(wloc-1.)*tden)
      exlog=1.45-7.8*sqrt(arg1)
      ordlog=2.45-8.58*sqrt(arg2)
      alphae(nz)=10.**exlog/(wloc**2)
      alphao(nz)=10.**ordlog/(wloc**2)
  200 continue
c
      if(nn.lt.2) go to 999
      write(6,9010) wfreq
 9010 format(1h1,29h  opacity output...frequency=  ,e11.3,
     $       9h  rad/sec, /, 1x, 81(1h-), /
     $     , 5h zone
     $     ,11h  temp
     $     ,11h  bavg
     $     ,11h  alphae
     $     ,11h  alphao
     $     ,/, 1x, 81(1h-))
c
      do 300 nz=1,nzons
      write(6,9020) nz,tt(nz),bb(nz),alphae(nz),alphao(nz)
  300 continue
 9020      format (i5, 4e11.3)
c
  999 continue
      return
      end

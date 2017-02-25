c@syopac  /11040/baldur/code/bald/dimprad.f
      subroutine syopac(nzons,tt,bb,nn)
!cap
      include 'cparm.m'
      parameter ( kq=mj )
      common/alsync/  wfreq,alphae(kq),alphao(kq)
      dimension tt(kq),bb(kq)
c
c        this routine computes the direction averaged dimensionless
c        absorption coefficients in each zone from an alalytic fit.
c        this fit is valid for temperatures from about 10 to 120 kev
c        and frequencies > 3*omega-gyro.
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
 9010 format(1h1,'  opacity output...frequency=  ',e11.3,
     &  '  rad/sec', /, 1x, 81('-'), /,' zone','  temp     ',
     &  '  bavg     ','  alphae   ','  alphao   ',/, 1x, 81('-'))
c
      do 300 nz=1,nzons
      write(6,9020) nz,tt(nz),bb(nz),alphae(nz),alphao(nz)
  300 continue
 9020      format (i5, 4e11.3)
c
  999 continue
      return
      end

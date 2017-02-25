c@ofit  .../baldur/code/bald/ofit.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 29-may-96 save zfuzz, zexp10, zlog10, inital
c       amck 27-jan-78 change subscr. exp. for icl fortran
c       amck 14-jun-77 fix ic in 218 loop
c       amck 13-jun-77 add pdloss, interpolate z to 1 for te -> 0
c       amck 16-mar-77 zte -> pte
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine ofit(pte,pne,pni,ploss,pdloss,pz,pz2,kimp)
c
c
cl      p.2)    compute radiation loss, z-bar, and z squared bar
cl              for each impurity at some te.
c               otable must have been called once to initialize
c               comcor.
c
c
c---------------------------------------------------------------------
cl                  c8.1     comcor--coronal equilibrium fit coef
c     version amck 1-aug-77
       common/comcor/
     r  ocoef , orange,
     i  mocoef, mxorng
       dimension
     r   ocoef(6,20,3),      orange(2,20,3)
c
        dimension       pni(kimp),      ploss(kimp),    pdloss(kimp),
     r  pz(kimp),       pz2(kimp)
c
        logical         inital
c
        data            inital /.true./
c
        save zfuzz, zexp10, zlog10, inital
c
c------------------------------------------------------------------------------
c
c
c       ofit must be called once for each point.
c       the units are determined by the units used for the fits;
c       the calling routine must do the conversion.
c       pni, ploss, pz, pz2 must be arrays with the impurity index
c       first, and ordered the same as in the call to otable
c       pni is the impurity density, not ion density.
c
c
c------------------------------------------------------------------------------
c
cl      1)      initial values, check initialization of comcor.
c
c
  100   continue
        if (kimp.le.0) return
c
        if (mocoef.le.0) call error_olymp(0,0,0,0,
     1          'ofit called before otable--comcor not set ')
c
        if(.not.inital) go to 120
        inital=.false.
        zfuzz = 1.0001
        zexp10 = log(10.0)
        zlog10 = 1.0 / zexp10
  120   continue
c
        zlte = log(pte) * zlog10
c
c
cl      2)      impurity fit
c
c
  200   continue
c
        do 268 ji = 1, kimp
        ploss(ji) = 0.0
        pdloss(ji) = 0.0
        pz(ji) = 1.0
        pz2(ji) = 1.0
c
        do 258 jt = 1, 3
c
c
c               find correct range(s) to compute fit
c
c
        irl = 0
        iru = 0
c
        ir1 = ji * mxorng
        ir0 = ir1 + 1 - mxorng
c
        do 204 jr = ir0, ir1
        if (orange(1,jr,jt).le.0.0) go to 206
        if (orange(1,jr,jt).le.pte) irl = jr
        if (orange(2,jr,jt)*zfuzz.le.pte) go to 204
        iru = jr
        go to 206
  204   continue
c
  206   continue
c
c
        if (irl.le.0.and.iru.le.0) go to 258
        if (irl.ne.iru) go to 220
c
c
c               normal case - - te is within a range for a fit
c
c
  210   continue
c
        zq = 0.0
c
        do 214 jc = 1, mocoef
        ic = mocoef+1-jc
        zq = zq*zlte + ocoef(ic,irl,jt)
  214   continue
c
        if (jt.ne.1) go to 250
        zq = pne * pni(ji) * exp(zq*zexp10)
c
c               d(ploss) / d(te)
c
        zdq = 0.0
c
        do 218 jc = 2, mocoef
        ic = mocoef + 2 - jc
        zdq = zdq*zlte + float(ic-1)*ocoef(ic,irl,jt)
  218   continue
c
        zdq = zdq * zq / pte
        go to 250
c
c
c               special cases - - below ranges, above ranges, or
c               in a gap between ranges
c
c               iru - - range which is above te, 0 if none
c               irl - - range which is below te, 0 if none
c               if below all ranges, z and z**2 are for lowest te for
c               which we have fits, loss goes as te**2 from lowest te
c               for which we have fits on down.
c               if above all ranges, z and z**2 are for highest te for
c               which we have fits, and loss goes as te**.5
c               if in a gap between ranges, values for nearest te's
c               for fits are linearly interpolated.
c
c
  220   continue
c
        if (iru.le.0) go to 230
c
c               value for lower end of range above te
c
        zteu = orange(1,iru,jt)
        zlteu = log(zteu) * zlog10
c
        z0 = 0.0
c
        do 224 jc = 1, mocoef
        ic = mocoef+1-jc
        z0 = z0*zlteu + ocoef(ic,iru,jt)
  224   continue
c
        if (jt.eq.2) zq = (z0 - 1.0)*(pte/zteu) + 1.0
        if (jt.eq.3) zq = (z0 - 1.0)*(pte/zteu)**2 + 1.0
c
        if (jt.gt.1) go to 228
        z0 = pne * pni(ji) * exp(z0*zexp10)
        zdq = z0 * 2.0 * pte/zteu**2
        zq = zdq * 0.5 * pte
  228   continue
c
        if (irl.le.0) go to 250
  230   continue
c
c               value for upper end of range below te
c
        ztel = orange(2,irl,jt)
        zltel = log(ztel) * zlog10
c
        z1 = 0.0
c
        do 234 jc = 1, mocoef
        ic = mocoef+1-jc
        z1 = z1*zltel + ocoef(ic,irl,jt)
  234   continue
c
        zq = z1
c
        if (jt.gt.1) go to 238
        z1 = pne * pni(ji) * exp(z1*zexp10)
c
        zdq = z1 * 0.5 / sqrt(pte*ztel)
        zq = zdq * pte * 2.0
  238   continue
c
        if (iru.le.0) go to 250
c
c               interpolate upper and lower values
c
  240   continue
c
        zint = (zteu - pte) / (zteu - ztel)
        zq = z0*(1.0 - zint) + z1*zint
        if (jt.eq.1) zdq = (z0 - z1) / (zteu - ztel)
c
  250   continue
        go to (252,254,256),jt
        go to 258
c
  252   continue
        ploss(ji) = zq
        pdloss(ji) = zdq
        go to 258
c
  254   continue
        pz(ji) = zq
        go to 258
c
  256   continue
        pz2(ji) = max(zq,pz(ji)**2 * 1.01)
        go to 258
c
  258   continue
c
  268   continue
c
        return
        end

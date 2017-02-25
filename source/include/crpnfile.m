c@crpnfile.m
c
c   Variables that used to be local to sbrtn nfile must now be put
c  in a common block to pass to sbrtns 
c
      common/cnfile/  zzee, zzei, zzeb, zzne, zznion
     & , zohme, zchxe, zalfe, zrade, zexte, zeion, zthdot, zeirs
     & , ztaueb, znibr, znebar, znibar, ztibar, ztebar, zv, zvc
     & , zbetae, zbetai, zbetab, zbeta
     & , ztbete, ztbeti, ztbetb, ztbeta, zalfap, zalfat
     & , zlint, zlambd, zlined, zzlne, zflux, ztenr1, ztenr2
     & , zeeflx, zeiflx, zicond, zecond,zineut, zeconv, ziconv
     & , ztote, zpowr, zphbem, zphrf, zphalf, zpath, zpathl
     & , zpabt, zpabtl, zpbinj, zpbabs, zpblos
     & , zlbem, zlrf, zprfin, zlath, zlabt, zebdot, zeadot, zpbnow
     & , zpbstr, zpss, zpcomp, zqinst, zqss, zqcomp, zqdenm, zqtech
     & , zptech, zip, zib, zva, zpcrsh
     & , ztconf, zbfact
c

c@clear   .../baldur/bald/clear.f
c 
      subroutine clear 
c 
      include 'cparm.m'
      include 'cbaldr.m' 
      include 'cd3he.m' 
      include 'csynrd.m' 
      include 'cmg.m' 
c 
c
c  cleargen 15:30 13 Jan 99
c
c sbrtn clear was generated using the cleargen program
c
      lnumer = 0
      cnumer = 0.0
      linout = 0
      cinout = 0.0
      lstart = 0
      cstart = 0.0
      ltrnsp = 0
      ctrnsp = 0.0
      lemprc = 0
      cemprc = 0.0
      lbeams = 0
      cbeams = 0.0
      ldivrt = 0
      cdivrt = 0.0
      lauxil = 0
      cauxil = 0.0
      limprd = 0
      cimprd = 0.0
      lnugas = 0
      cnugas = 0.0
      lneocl = 0
      cneocl = 0.0
      lfsion = 0
      cfsion = 0.0
      lstabl = 0
      cstabl = 0.0
      lbound = 0
      lbound(7) = 15
      lbound(8) = 1
      cbound = 0.0
      aspec = 0.0
      bzi = 0.0
      ellipt = 0.0
      rehe = 0.0
      rehi = 0.0
      reoff = 0.0
      reon = 0.0
      revols = 0.0
      rex1 = 0.0
      rex2 = 0.0
      rey1 = 0.0
      rey2 = 0.0
      rle0 = 0.0
      rle1 = 0.0
      rleaux = 0.0
      rleprf = 0.0
      rlepwr = 0.0
      rlepwrt = 0.0
      rlfreq = 0.0
      rli0 = 0.0
      rli1 = 0.0
      rliaux = 0.0
      rliprf = 0.0
!ap      rliprw = 0.0
      rlipwrt = 0.0
      rlpara = 0.0
      rlpowr = 0.0
      rlheprf = 0.0
      rlhepwr = 0.0
      rlheaux = 0.0
      nrenp = 0
      nrldmp = 0
      nrleex = 0
      nrliex = 0
      nrlmod = 0
      nrlpar = 0
      nspec = 0
      nzspec = 0
      bpoli = 0.0
      chi = 0.0
      chn = 0.0
      aaaa = 0.0
      alpha0 = 0.0
      alpha1 = 0.0
      azzz = 0.0
      bbbb = 0.0
      beta0 = 0.0
      beta1 = 0.0
      bouna0 = 0.0
      bouna1 = 0.0
      bounb0 = 0.0
      bounb1 = 0.0
      bounc0 = 0.0
      bounc1 = 0.0
      bound = 0.0
      bpa0 = 0.0
      bpa1 = 0.0
      bpb0 = 0.0
      bpb1 = 0.0
      bpc0 = 0.0
      bpc1 = 0.0
      cccc = 0.0
      cpedst = 0.0
      dddd = 0.0
      delta0 = 0.0
      delta1 = 0.0
      eeee = 0.0
      etai = 0.0
      ffff = 0.0
      gamma0 = 0.0
      gamma1 = 0.0
      lflxng = .false.
      condei = 0.0
      condii = 0.0
      cvctei = 0.0
      cvctii = 0.0
      ctotei = 0.0
      ctotii = 0.0
      bp11 = 0.0
      bp12 = 0.0
      cl11 = 0.0
      cl12 = 0.0
      cnueqs = 0.0
      dahs = 0.0
      dais = 0.0
      dbhis = 0.0
      dbhs = 0.0
      dbihs = 0.0
      dbiis = 0.0
      dbis = 0.0
      denes = 0.0
      detepr = 0.0
      detes = 0.0
      detis = 0.0
      dinhs = 0.0
      dinis = 0.0
      ditins = 0.0
      ditipr = 0.0
      ditis = 0.0
      dnhhs = 0.0
      dnhis = 0.0
      dnhs = 0.0
      dnihs = 0.0
      dniis = 0.0
      dnis = 0.0
      dweirs = 0.0
      dwicxs = 0.0
      dxemps = 0.0
      ps11 = 0.0
      ps12 = 0.0
      recoms = 0.0
      salfs = 0.0
      shbems = 0.0
      shblos = 0.0
      shchxs = 0.0
      shfus = 0.0
      shions = 0.0
      siions = 0.0
      velnhs = 0.0
!ap      velhis = 0.0
      veltes = 0.0
      veltis = 0.0
      vewars = 0.0
      vnwars = 0.0
      vxemps = 0.0
      wealfs = 0.0
      weauxs = 0.0
      webems = 0.0
      webrs = 0.0
      weecrh = 0.0
      weicrf = 0.0
      weions = 0.0
      weirs = 0.0
      weohms = 0.0
      wesrs = 0.0
      wialfs = 0.0
      wiauxs = 0.0
      wibems = 0.0
      wichxs = 0.0
      wiecrh = 0.0
      wiicrf = 0.0
      wiions = 0.0
      xeemps = 0.0
      xiemps = 0.0
      dites = 0.0
      trdif = 0.0
      trvel = 0.0
      lkeflg = 0
      lkiflg = 0
      shfuel = 0.0
      sifuel = 0.0
      wed3fs = 0.0
      wid3fs = 0.0
      wid3fl = 0.0
      weash = 0.0
      shd3fs = 0.0
      sid3fs = 0.0
      wesyn = 0.0
      wefus = 0.0
      wifus = 0.0
      dnneo1 = 0.0
      vnneo1 = 0.0
      xeneo1 = 0.0
      xineo1 = 0.0
      xeneo2 = 0.0
      fdr = 0.0
      fig = 0.0
      fti = 0.0
      frm = 0.0
      fkb = 0.0
      frb = 0.0
      fhf = 0.0
      fmh = 0.0
      fec = 0.0
      feg = 0.0
      fdrint = 0.0
      eithes = 0.0
      dxthes = 0.0
      vxthes = 0.0
      xethes = 0.0
      xithes = 0.0
      weithe = 0.0
      weiths = 0.0
      difthi = 0.0
      velthi = 0.0
      cthery = 0.0
      sfutz = 0.0
      lthery = 0
      thdre = 0.0
      thrme = 0.0
      thrbe = 0.0
      thkbe = 0.0
      thhfe = 0.0
      thdri = 0.0
      thrmi = 0.0
      thrbi = 0.0
      thkbi = 0.0
      thhfi = 0.0
      thige = 0.0
      thigi = 0.0
      thtie = 0.0
      thtii = 0.0
      threte = 0.0
      threti = 0.0
      thdinu = 0.0
      thfith = 0.0
      thdte = 0.0
      thdi = 0.0
      thbeta = 0.0
      thlni = 0.0
      thlti = 0.0
      thdia = 0.0
      thnust = 0.0
      thrstr = 0.0
      theps = 0.0
      thlsh = 0.0
      thlpr = 0.0
      thlarp = 0.0
      thrhos = 0.0
      thdiaw = 0.0
      thlarpo = 0.0
      thvthe = 0.0
      thvthi = 0.0
      thsoun = 0.0
      thalfv = 0.0
      thbpbc = 0.0
      thetth = 0.0
      thsrhp = 0.0
      thdias = 0.0
      thlamb = 0.0
      thrlwe = 0.0
      thrlwi = 0.0
      thnme = 0.0
      thnmi = 0.0
      thhme = 0.0
      thcee = 0.0
      thcei = 0.0
      thitie = 0.0
      thitii = 0.0
      thrbgb = 0.0
      thrbb = 0.0
      thvalh = 0.0
      thdwe = 0.0
      thdwi = 0.0
      thfbth = 0.0
      thbp2 = 0.0
      thomr = 0.0
      wexba = 0.0
      vrota = 0.0
      xwexba = 0.0
      twexba = 0.0
      nxwexba = 0
      ntwexba = 0
      bedgi = 0.0
      compl = 0.0
      compp = 0.0
      dvolx = 0.0
      rcurri = 0.0
      redgi = 0.0
      rmji = 0.0
      tcompi = 0.0
      tcoef = 0.0
      coeft = 0.0
      nzcomp = 0
      nlcomp = .false.
      ajboot = 0.0
      ajtpbi = 0.0
      ahmean = 0.0
      aimass = 0.0
      aired = 0.0
      ajzs = 0.0
      bpols = 0.0
      bzs = 0.0
      c2mean = 0.0
      calph = 0.0
      ccnu = 0.0
      cdbohm = 0.0
      cdetes = 0.0
      cditis = 0.0
      cdnhis = 0.0
      cdnhs = 0.0
      cdniis = 0.0
      ceta = 0.0
      cfutz = 0.0
      cloges = 0.0
      clogis = 0.0
      cmean = 0.0
      cnuel = 0.0
      cnuhyd = 0.0
      cteles = 0.0
      ctions = 0.0
      dfutzd = 0.0
      dfutze = 0.0
      dfutzi = 0.0
      dfutzv = 0.0
      dxsemi = 0.0
      eta = 0.0
      ev50s = 0.0
      extf = 0.0
      extzef = 0.0
      ftrap = 0.0
      gspitz = 0.0
      q = 0.0
      qmhd1 = 0.0
      qmhd2 = 0.0
      qsmth = 0.0
      qstar = 0.0
      rhoels = 0.0
      rhohs = 0.0
      rhoins = 0.0
      rhois = 0.0
      rhosms = 0.0
      rmajs = 0.0
      rmins = 0.0
      rrstar = 0.0
      scroff = 0.0
      shear = 0.0
      slbps = 0.0
      slnes = 0.0
      slprs = 0.0
      slpts = 0.0
      sltes = 0.0
      sltis = 0.0
      telecs = 0.0
      tes = 0.0
      tesms = 0.0
      thpsms = 0.0
      thrprs = 0.0
      tions = 0.0
      tis = 0.0
      tisms = 0.0
      topsms = 0.0
      vxsemi = 0.0
      xesemi = 0.0
      xisemi = 0.0
      xnuel = 0.0
      xnuhyd = 0.0
      xzeff = 0.0
      memprt = 0
      tmod = 0.0
      snestr = 0.0
      snebar = 0.0
      betate = 0.0
      betati = 0.0
      betatb = 0.0
      betata = 0.0
      betatt = 0.0
      betape = 0.0
      betapi = 0.0
      betapb = 0.0
      betapa = 0.0
      betapt = 0.0
      curnts = 0.0
      erges = 0.0
      ergis = 0.0
      ergbs = 0.0
      ergas = 0.0
      ergts = 0.0
      envaes = 0.0
      envais = 0.0
      enlaes = 0.0
      enlais = 0.0
      gealfs = 0.0
      geauxs = 0.0
      gebems = 0.0
      gebrs = 0.0
      geecrs = 0.0
      geicrs = 0.0
      geions = 0.0
      geohms = 0.0
      gesrs = 0.0
      gialfs = 0.0
      giauxs = 0.0
      gibems = 0.0
      gichxs = 0.0
      giecrs = 0.0
      giicrs = 0.0
      giions = 0.0
      cgbeta = 0.0
      cgpowr = 0.0
      gcoefs = 0.0
      gxp = 0.0
      gcomb = 0.0
      armins = 0.0
      armajs = 0.0
      nrad = 0
      totprs = 0.0
      rhisms = 0.0
      slnis = 0.0
      ged3fs = 0.0
      geashs = 0.0
      gid3fs = 0.0
      gesyns = 0.0
      geirs = 0.0
      ergds = 0.0
      betatd = 0.0
      betapd = 0.0
      gblosi = 0.0
      gfdni = 0.0
      gflowi = 0.0
      gfluxi = 0.0
      gflx0i = 0.0
      gfrac1 = 0.0
      gnchek = 0.0
      gnel1 = 0.0
      gneut = 0.0
      grecyc = 0.0
      gsputs = 0.0
      gtchek = 0.0
      gtflwi = 0.0
      gtprfi = 0.0
      gvsrci = 0.0
      gwmin = 0.0
      gxdome = 0.0
      gxmaxe = 0.0
      gxmine = 0.0
      gxphi = 0.0
      gxthet = 0.0
      rhons = 0.0
      tns = 0.0
      ngpart = 0
      ngprof = 0
      ngsplt = 0
      ngsprf = 0
      ngxene = 0
      ngzone = 0
      nlglim = .false.
      nlgmon = .false.
      nlgpin = .false.
      nlgref = .false.
      nlgspt = .false.
      ajbs = 0.0
      bpabs = 0.0
      bpinjs = 0.0
      bploss = 0.0
      cjbeam = 0.0
      habeam = 0.0
      halfas = 0.0
      hangle = 0.0
      haper = 0.0
      haperv = 0.0
      hdds = 0.0
      hddtot = 0.0
      hdiv = 0.0
      hdivv = 0.0
      hebeam = 0.0
      hebems = 0.0
      height = 0.0
      hfocl = 0.0
      hfoclv = 0.0
      hfract = 0.0
      hfutz = 0.0
      hibeam = 0.0
      hlenth = 0.0
      hnchek = 0.0
      hpbeam = 0.0
      hpowmw = 0.0
      hprof = 0.0
      hprofv = 0.0
      hr = 0.0
      hrmaj = 0.0
      hrmin = 0.0
      htchek = 0.0
      htoff = 0.0
      hton = 0.0
      hwidth = 0.0
      hzbeam = 0.0
      rhobes = 0.0
      rhobis = 0.0
      lhbeam = 0
      mhbeam = 0
      mxhbem = 0
      mxhfr = 0
      nhaper = 0
      nhbeam = 0
      nhe = 0
      nhmu = 0
      nhprfv = 0
      nhprof = 0
      nhshap = 0
      nhskip = 0
      nhsrc = 0
      nipart = 0
      niprof = 0
      nizone = 0
      aalpha = 0.0
      abfuse = 0.0
      abouni = 0.0
      acons = 0.0
      adds = 0.0
      addtot = 0.0
      adifus = 0.0
      afslow = 0.0
      afuses = 0.0
      alfsrc = 0.0
      alphai = 0.0
      aoloss = 0.0
      aslows = 0.0
      atbi = 0.0
      atfuse = 0.0
      azalfa = 0.0
      ealfai = 0.0
      ebouni = 0.0
      econsi = 0.0
      efusei = 0.0
      rcwals = 0.0
      rdwals = 0.0
      ntype = 0
      delmax = 0.0
      dti = 0.0
      dtmaxi = 0.0
      dtmini = 0.0
      dtoldi = 0.0
      errmax = 0.0
      smrlow = 0.0
      smlwcy = 0.0
      tai = 0.0
      tbi = 0.0
      theta = 0.0
      thetai = 0.0
      thetap = 0.0
      tmaxi = 0.0
      tspare = 0.0
      lsmord = 0
      natomc = 0
      nbound = 0
      nfusn = 0
      nitmax = 0
      nounit = 0
      ntrans = 0
      lmpirc = .false.
      nldiff = .false.
      nlextr = .false.
      nliter = .false.
      nlrcom = .false.
      nlrpet = .false.
      nlsorc = .false.
      nlsord = .false.
      ltheor = .false.
      lholab = "  "
      evs = 0.0
      evsinv = 0.0
      fcau = 0.0
      fxau = 0.0
      gamin1 = 0.0
      ueib = 0.0
      ueid = 0.0
      ueie = 0.0
      ueih = 0.0
      ueii = 0.0
      ueij = 0.0
      ueil = 0.0
      ueim = 0.0
      ueip = 0.0
      ueir = 0.0
      ueit = 0.0
      ueiv = 0.0
      uesb = 0.0
      uesd = 0.0
      uese = 0.0
      uesh = 0.0
      uesi = 0.0
      uesj = 0.0
      uesl = 0.0
      uesm = 0.0
      uesp = 0.0
      uesr = 0.0
      uest = 0.0
      uesv = 0.0
      uieb = 0.0
      uied = 0.0
      uiee = 0.0
      uieh = 0.0
      uiei = 0.0
      uiej = 0.0
      uiel = 0.0
      uiem = 0.0
      uiep = 0.0
      uier = 0.0
      uiet = 0.0
      uiev = 0.0
      uisb = 0.0
      uisd = 0.0
      uise = 0.0
      uish = 0.0
      uisi = 0.0
      uisj = 0.0
      uisl = 0.0
      uism = 0.0
      uisp = 0.0
      uisr = 0.0
      uist = 0.0
      uisv = 0.0
      useb = 0.0
      used = 0.0
      usee = 0.0
      useh = 0.0
      usei = 0.0
      usej = 0.0
      usel = 0.0
      usem = 0.0
      usep = 0.0
      user = 0.0
      uset = 0.0
      usev = 0.0
      usib = 0.0
      usid = 0.0
      usie = 0.0
      usih = 0.0
      usii = 0.0
      usij = 0.0
      usil = 0.0
      usim = 0.0
      usip = 0.0
      usir = 0.0
      usit = 0.0
      usiv = 0.0
      uind = 0.0
      uine = 0.0
      unid = 0.0
      unie = 0.0
      usnd = 0.0
      usne = 0.0
      unsd = 0.0
      unse = 0.0
      epsinv = 0.0
      epslon = 0.0
      rndeps = 0.0
      rndup = 0.0
      ninfin = 0
      xaci1 = 0.0
      xaci2 = 0.0
      xaci3 = 0.0
      xaddi1 = 0.0
      xaddi2 = 0.0
      xaddi3 = 0.0
      xadi1 = 0.0
      xadi2 = 0.0
      xadi3 = 0.0
      xaii1 = 0.0
      xaii2 = 0.0
      xaii3 = 0.0
      xaoi1 = 0.0
      xaoi2 = 0.0
      xaoi3 = 0.0
      xbcon1 = 0.0
      xbcon2 = 0.0
      xbcon3 = 0.0
      xbpol1 = 0.0
      xbpol2 = 0.0
      xbpol3 = 0.0
      xccon1 = 0.0
      xccon2 = 0.0
      xccon3 = 0.0
      xchi1 = 0.0
      xchi2 = 0.0
      xchi3 = 0.0
      xdel = 0.0
      xdti1 = 0.0
      xdtol1 = 0.0
      xebti1 = 0.0
      xebti2 = 0.0
      xebti3 = 0.0
      xeomi1 = 0.0
      xeomi2 = 0.0
      xeomi3 = 0.0
      xepoi1 = 0.0
      xepoi2 = 0.0
      xepoi3 = 0.0
      xerr = 0.0
      xfci1 = 0.0
      xfci2 = 0.0
      xfci3 = 0.0
      xfdi1 = 0.0
      xfdi2 = 0.0
      xfdi3 = 0.0
      xfii1 = 0.0
      xfii2 = 0.0
      xfii3 = 0.0
      xfoi1 = 0.0
      xfoi2 = 0.0
      xfoi3 = 0.0
      xfutz = 0.0
      xtai1 = 0.0
      xtbi1 = 0.0
      xtot1 = 0.0
      xtot2 = 0.0
      xtot3 = 0.0
      xwbi1 = 0.0
      xwbi2 = 0.0
      xwbi3 = 0.0
      xwomi1 = 0.0
      xwomi2 = 0.0
      xwomi3 = 0.0
      xwpoi1 = 0.0
      xwpoi2 = 0.0
      xwpoi3 = 0.0
      liter = 0
      lpdel = 0
      lperr = 0
      lzdel = 0
      lzerr = 0
      acompi = 0.0
      addi = 0.0
      aflxii = 0.0
      aflxoi = 0.0
      asorci = 0.0
      asordi = 0.0
      bcons = 0.0
      begini = 0.0
      ccons = 0.0
      ebini = 0.0
      ebtoti = 0.0
      ecompi = 0.0
      eohmi = 0.0
      epoyni = 0.0
      fcompi = 0.0
      fflxii = 0.0
      fflxoi = 0.0
      fsorci = 0.0
      fsordi = 0.0
      totali = 0.0
      wbi = 0.0
      wcompi = 0.0
      wohmi = 0.0
      wpoyni = 0.0
      bint = 0.0
      dx2i = 0.0
      dx2inv = 0.0
      dxboui = 0.0
      dxzoni = 0.0
      gx = 0.0
      rmaji = 0.0
      rmini = 0.0
      xbouni = 0.0
      xzoni = 0.0
      lalpha = 0
      lcentr = 0
      ldeut = 0
      ledge = 0
      lelec = 0
      lhe3 = 0
      lhyd1 = 0
      lhydn = 0
      limp1 = 0
      limpn = 0
      lion = 0
      ltrit = 0
      mchi = 0
      mhyd = 0
      mimp = 0
      mxchi = 0
      mxhyd = 0
      mximp = 0
      mxions = 0
      mxt = 0
      mxt1 = 0
      mxzone = 0
      mzones = 0
      lprotn = 0
      cfev = 0.0
      cfevdk = 0.0
      cxev = 0.0
      cxevdk = 0.0
      fca0 = 0.0
      fcae = 0.0
      fcalfa = 0.0
      fcan = 0.0
      fcap = 0.0
      fcc = 0.0
      fcc1 = 0.0
      fcc2 = 0.0
      fce = 0.0
      fces = 0.0
      fcf = 0.0
      fcg = 0.0
      fch = 0.0
      fck = 0.0
      fclnae = 0.0
      fclnap = 0.0
      fcme = 0.0
      fcmn = 0.0
      fcmp = 0.0
      fcna = 0.0
      fcpi = 0.0
      fcr = 0.0
      fcre = 0.0
      fcrinf = 0.0
      fcsb = 0.0
      fcsigt = 0.0
      fcv0 = 0.0
      fcwien = 0.0
      fxa0 = 0.0
      fxae = 0.0
      fxalfa = 0.0
      fxc = 0.0
      fxc1 = 0.0
      fxc2 = 0.0
      fxe = 0.0
      fxes = 0.0
      fxf = 0.0
      fxg = 0.0
      fxh = 0.0
      fxk = 0.0
      fxlnae = 0.0
      fxlnap = 0.0
      fxme = 0.0
      fxmn = 0.0
      fxmp = 0.0
      fxna = 0.0
      fxnucl = 0.0
      fxr = 0.0
      fxre = 0.0
      fxrinf = 0.0
      fxsb = 0.0
      fxsigt = 0.0
      fxv0 = 0.0
      fxwien = 0.0
      bdhyde = 0.0
      bdimpe = 0.0
      bdtee = 0.0
      bdtie = 0.0
      bdtime = 0.0
      d0nmon = 0.0
      gainv = 0.0
      ftzeff = 0.0
      apresr = 0.0
      bpoid = 0.0
      bz = 0.0
      bzstar = 0.0
      cpvelc = 0.0
      cpvion = 0.0
      curent = 0.0
      denga0 = 0.0
      denga1 = 0.0
      dengas = 0.0
      denim0 = 0.0
      denim1 = 0.0
      denimp = 0.0
      denmon = 0.0
      dens = 0.0
      dens0 = 0.0
      dens1 = 0.0
      dtinit = 0.0
      dtmax = 0.0
      dtmin = 0.0
      ebfit = 0.0
      eebfit = 0.0
      eefit = 0.0
      eehfit = 0.0
      eeifit = 0.0
      eeteft = 0.0
      eetift = 0.0
      efit = 0.0
      ehfit = 0.0
      eifit = 0.0
      eioniz = 0.0
      elecd0 = 0.0
      electe = 0.0
      etamod = 0.0
      etefit = 0.0
      etifit = 0.0
      flgas = 0.0
      flimp = 0.0
      fracth = 0.0
      fracti = 0.0
      gflmax = 0.0
      gflmon = 0.0
      gfract = 0.0
      gftime = 0.0
      gpers = 0.0
      radius = 0.0
      rastar = 0.0
      rlined = 0.0
      rmajor = 0.0
      rminor = 0.0
      rnebar = 0.0
      rpa = 0.0
      rpela = 0.0
      rqs = 0.0
      tbpoid = 0.0
      tcomp = 0.0
      tcold = 0.0
      tcoldp = 0.0
      te = 0.0
      te0 = 0.0
      te1 = 0.0
      tgas = 0.0
      ti = 0.0
      ti0 = 0.0
      ti1 = 0.0
      timp = 0.0
      tinit = 0.0
      tmax = 0.0
      tpela = 0.0
      vpela = 0.0
      cfz230 = 0.0
      itvpel = 0.0
      wtgas = 0.0
      wtimp = 0.0
      frcor = 0.0
      frout = 0.0
      rpelc = 0.0
      ypa = 0.0
      nbfit = 0
      nedit = 0
      ngas = 0
      nimp = 0
      nnfit = 0
      nnhfit = 0
      nnifit = 0
      npel = 0
      npel2 = 0
      npelga = 0
      npelou = 0
      npels = 0
      npelsa = 0
      npuff = 0
      nrfit = 0
      nskip = 0
      ntefit = 0
      ntifit = 0
      ntty = 0
      ntychl = 0
      nzones = 0
      npelgc = 0
      lhgas = " "
      roo1 = 0.0
      ree1 = 0.0
      roe1 = 0.0
      reo1 = 0.0
      tacrt = 0.0
      areaw1 = 0.0
      ash4 = 0.0
      ashp = 0.0
      fploss = 0.0
      f4loss = 0.0
      ftloss = 0.0
      f3loss = 0.0
      alcmg = 0.0
      epscmg = 0.0
      alpcmg = 0.0
      gcmg = 0.0
      tohcmg = 0.0
      atim = 0.0
      pauxe = 0.0
      pauxi = 0.0
      apaux = 0.0
      apauxe = 0.0
      apauxi = 0.0
      stim = 0.0
      spd = 0.0
      spt = 0.0
      spp = 0.0
      sp3 = 0.0
      sp4 = 0.0
      spimp4 = 0.0
      spaux = 0.0
      spauxi = 0.0
      ekappa = 0.0
      edelta = 0.0
      nzons1 = 0
      iedit1 = 0
      ndbug1 = 0
      lcmg = 0
      lpaux = 0
      lspaux = 0
      lrepld = 0
      lreplt = 0
      lrepl3 = 0
      tedit = 0.0
      tplot = 0.0
      sedit = 0.0
      snplt = 0.0
      snprt = 0.0
      splot = 0.0
      tnplt = 0.0
      tnprt = 0.0
      lpage = 0
      nbdump = 0
      nbuf = 0
      nediti = 0
      nnplot = 0
      nnprnt = 0
      nplot = 0
      nrgraf = 0
      nsedit = 0
      nskpi = 0
      ntgraf = 0
      nllprt = .false.
      nlpomt = .false.
      tediti = 0.0
      tploti = 0.0
      lhspec = " "
      lhb = "  "
      lhcden = "  "
      lhcurr = "  "
      lhdenr = "  "
      lhdens = "  "
      lheden = "  "
      lhefld = "  "
      lheflx = "  "
      lhener = "  "
      lhlen = "  "
      lhnflx = "  "
      lhpden = "  "
      lhpowr = "  "
      lhtemp = "  "
      lhtime = "  "
      tppd = 0.0
      tpp3 = 0.0
      tpt = 0.0
      tp3 = 0.0
      tp4t = 0.0
      tp43 = 0.0
      tp42t = 0.0
      tepd = 0.0
      tep3 = 0.0
      tet = 0.0
      te3 = 0.0
      te4t = 0.0
      te43 = 0.0
      te42t = 0.0
      fpde = 0.0
      fp3e = 0.0
      fte = 0.0
      f3e = 0.0
      f4te = 0.0
      f43e = 0.0
      f42te = 0.0
      fddt = 0.0
      fdd3 = 0.0
      fdt = 0.0
      fdtb = 0.0
      fd3 = 0.0
      fd3b = 0.0
      ftt4 = 0.0
      fttb = 0.0
      epde = 0.0
      ep3e = 0.0
      etre = 0.0
      e3e = 0.0
      e4te = 0.0
      e43e = 0.0
      e42te = 0.0
      epdi = 0.0
      ep3i = 0.0
      etri = 0.0
      e3i = 0.0
      e4ti = 0.0
      e43i = 0.0
      e42ti = 0.0
      drnp = 0.0
      drnd = 0.0
      drnt = 0.0
      drn3 = 0.0
      drn4 = 0.0
      q1fus = 0.0
      q2fus = 0.0
      q1lfus = 0.0
      q2lfus = 0.0
      qefus = 0.0
      wntn = 0.0
      d3fast = 0.0
      fnpd0 = 0.0
      fnp30 = 0.0
      fnt0 = 0.0
      fn30 = 0.0
      fn4t0 = 0.0
      fn430 = 0.0
      fn42t0 = 0.0
      fnpd = 0.0
      fnp3 = 0.0
      fnt = 0.0
      fn3 = 0.0
      fn4t = 0.0
      fn43 = 0.0
      fn42t = 0.0
      rh1fst = 0.0
      rh2fst = 0.0
      epd0 = 0.0
      ep30 = 0.0
      et0 = 0.0
      e30 = 0.0
      e4t0 = 0.0
      e430 = 0.0
      e42t0 = 0.0
      epd = 0.0
      ep3 = 0.0
      et = 0.0
      e3 = 0.0
      e4t = 0.0
      e43 = 0.0
      e42t = 0.0
      ree = 0.0
      roo = 0.0
      reo = 0.0
      roe = 0.0
      bmin = 0.0
      bmax = 0.0
      pwall = 0.0
      ptot = 0.0
      areaw = 0.0
      taucrt = 0.0
      stemp = 0.0
      sdens = 0.0
      bavg = 0.0
      svol = 0.0
      surf = 0.0
      wsyn = 0.0
      teqp = 0.0
      deny = 0.0
      voq = 0.0
      ploss = 0.0
      nzons = 0
      nw1 = 0
      ndbug = 0
      iedit = 0
      xcmg = 0.0
      xub = 0.0
      zstara = 0.0
      dxcmg = 0.0
      vxcmg = 0.0
      xcmg0 = 0.0
      xub0 = 0.0
c
         call gclear
         call hclear
         call aclear
cbate         call iclear

      return
      end

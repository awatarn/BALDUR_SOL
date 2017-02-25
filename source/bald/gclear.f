c@gclear  .../baldur/bald/gclear.f 
c rgb 24-jan-99 gclear generated from cleargen 
c aes 04-mar-82 add nlpion to comusr and dimension nlog24(11) 
c 
      subroutine gclear 
c 
      include 'cmonprt.m' 
      include 'cmonxs.m' 
      include 'cmonte.m' 
      include 'cmonspl.m' 
c 
c
c  cleargen 15:30 13 Jan 99
c
c sbrtn clear was generated using the cleargen program
c
      e0 = 0.0
      flxfac = 0.0
      pi = 0.0
      rprob = 0.0
      timint = 0.0
      vel = 0.0
      velx = 0.0
      vely = 0.0
      velz = 0.0
      weight = 0.0
      wtmin = 0.0
      wtot = 0.0
      x0 = 0.0
      y0 = 0.0
      z0 = 0.0
      msurf = 0
      ncell = 0
      nescap = 0
      ngas = 0
      ninc = 0
      njump = 0
      nsrc = 0
      nlrefl = .false.
      sigvcx = 0.0
      sigvi = 0.0
      sigvt = 0.0
      sigvx = 0.0
      tablie = 0.0
      maxtab = 0
      den0 = 0.0
      den0in = 0.0
      denein = 0.0
      denion = 0.0
      e0in = 0.0
      eneut = 0.0
      fluxin = 0.0
      fluxn = 0.0
      fluxr = 0.0
      fract = 0.0
      outflx = 0.0
      rmass = 0.0
      rnu = 0.0
      rplsma = 0.0
      rsplit = 0.0
      rsurf = 0.0
      sflux = 0.0
      sigma = 0.0
      sigvbx = 0.0
      sne = 0.0
      sni = 0.0
      snvol = 0.0
      spe = 0.0
      spi = 0.0
      spii = 0.0
      tein = 0.0
      tiin = 0.0
      veff = 0.0
      wmin = 0.0
      maxgas = 0
      maxrad = 0
      maxsrc = 0
      maxsur = 0
      ngases = 0
      ngasi = 0
      npts = 0
      nrandm = 0
      nsrces = 0
      nsrci = 0
      nsurf = 0
      nvsrc = 0
      nlerr = .false.
      nlmod = .false.
      nlmono = .false. 
      nlpion = .false.
      nlscat = .false.
      nlsput = .false.
      rnumb = 0.0
      vl = 0.0
      vxl = 0.0
      vyl = 0.0
      vzl = 0.0
      wl = 0.0
      xl = 0.0
      yl = 0.0
      zl = 0.0
      maxlev = 0
      mcell = 0
      msurfl = 0
      ngasl = 0
      nlevel = 0
      nodes = 0
      nu = 0
      nlsplt = .false.
      nlsurf = .false.

      return
      end

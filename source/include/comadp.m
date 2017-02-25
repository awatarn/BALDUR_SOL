c@comadp.m
c                            comadp--atomic data package variablesc
c  See note in comncr regarding dimensions
c
      common /adsdat/ nucz, nspc, apn(100,10), aqn(100,10),
     &                nvalnc(100), sigma(10,10), qn(100,10),
     &                en(100,10), ein(100,5), enm(100,5,10),
     &                fnm(100,5,10), fnn(100), enn(100)
c
      common /adparm/ laden, ladtip, ltipok, leci, leciok
c
      common /addiel/ cdnn(100), cdnm(100), ldrmlt, yhnna, yhnnb,
     &                yhnnc, yhnma, yhnmb, yhnmc
c
      common /adelec/ rclion(100), rrarec(100), rdirec(100),
     &                cizlos(100), radrrc(100), raddrc(100),
     &                radclx(100), radbrm(100)
c
      common /adneut/ rcxrec(100), radbcx(100)
c
c  Note: in ncxx(2) the 2 refers to the number of hydrogenic species
c
      common /ncadp/ cermlt(2), cizmlt(2), cxrmlt(2), dnnmlt(28,2),
     1               dnmmlt(28,2), yhnn(3,2), yhnm(3,2), ndrmlt(2),
     2               ncxx(2), ncxxb, naden(2), nadtip(2), neci(2)
c
c
c  These common block contains variables needed in the NC impurity charge
c  state transport code to control and use the atomic data package. It
c  also contains variables internal to the atomic data package. The latter
c  are not documented here. Input variables are denoted by "*"; their
c  descriptions are followed by their default values.
c           
c ...In the following, a subscript i refers to impurity species i...
c  cermlt(i)    * Scales recombination reaction and radiation rates; 1.0
c  cizmlt(i)    * Scales ionization reaction and loss rates; 1.0
c  cxrmlt(i)    * Scales charge exchange reaction and radiation rates; 1.0
c  dnnmlt(k,i)  * Multiplies diel. recombination n-n rate, stage k-1; 1.0
c  dnmmlt(k,i)  * Multiplies diel. recombination n-m rate, stage k-1; 1.0
c  yhnn(l,i)    * l-th Hahn factor for diel. recomb. n-n rate; .2, .7, 1.4
c  yhnm(l,i)    * l-th Hahn factor for diel. recomb. n-m rate; .1, .55, 1.
c  ndrmlt(i)    * Controls options for calculating dn_mlt; 1 (dn_mlt=1.)
c  ncxx(ih)     * Controls CX X-section formula for hyd. species ih; 3 (OSCT)
c  ncxxb        * Same for neutral beams (not implemented); 3
c  naden(i)     * Controls energy level calculation formula; 1 (More)
c  nadtip(i)    * If = 1, adds tabulated ionization potentials; 1
c  neci(i)      * Controls ionization rate formula; 2 (Belfast)
c  cdnn(k)      ADPAK variable for dnnmlt
c  cdnm(k)        "      "      "  dnmmlt
c  ldrmlt         "      "      "  ndrmlt
c  yhnna          "      "      "  yhnn(1,i)
c  yhnnb          "      "      "  yhnn(2,i)
c  yhnnc          "      "      "  yhnn(3,i)
c  yhnma          "      "      "  yhnm(1,i)
c  yhnmb          "      "      "  yhnm(2,i)
c  yhnmc          "      "      "  yhnm(3,i)
c  laden          "      "      "  naden
c  ladtip         "      "      "  nadtip
c  leci           "      "      "  neci
c  nspc         Number of charge states + 1
c  nvalnc(k)    Principal quantum number of valence shell, charge state k-1
c  en(k,n)      Energy levels for recombination of ion, stage k-1, shell n
c ...Except as indicated, these rates are in cm**3/s and losses in W*cm**3...
c  rclion(k)    Electron collisional ionization rate of stage k-1 
c  rrarec(k)    Radiative recombination rate of stage k-1
c  rdirec(k)    Dielectron recombination rate of stage k-1
c  cizlos(k)    Electron losses coupled with rclion
c  radrrc(k)    Radiation associated with rrarec
c  raddrc(k)    Radiation associated with rdirec
c  radclx(k)    Electron losses & rad'n due to coll'al excit'n of stage k-1
c  radbrm(k)    Bremsstrahlung radiation due to stage k-1
c  rcxrec(k)    CX recombination rate due to stage k-1 in sec**(-1)
c  radbcx(k)    Radiation associated with rcxrec in W
c
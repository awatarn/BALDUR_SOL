c******************************************************************************
c 11:00 22-may-94 .../baldur/code/bald/units.f  Bateman, PPPL
c--------1---------2---------3---------4---------5---------6---------7-c
c**********************************************************************c
c@units  .../baldur/code/bald/default.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c rgb 22-may-94 added uind, uine, unid, unie, usnd, usne, unsd, unse
c   for normalization of chn(j1,jz)
c       amck 12-sep-77 add lhpden, lhdenr, lheflx, lhnflx
c       amck 13-sep-76 add uxxp
c       amck 16-aug-76 add lhefld
c       amck 11-aug-76 add u--j, u--r, u--v
c       amck 23-jul-76 add lhpowr
c--------1---------2---------3---------4---------5---------6---------7-c
c
c       -----------
c       sbrtn UNITS   file DEFAULT
c       -----------
c
c
cl      1.10    set conversion factors among external,
c               internal, and standard units
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine units
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'commhd.m'
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
        data    iclass /1/,     isub /10/
c
c
c
        if (.not.nlomt1(isub)) go to 10
        call mesage(' *** 1.10 subroutine units bypassed ')
        return
   10   continue
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      common blocks and variables modified:
c
c       comcnv, comout, commhd
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c       baldur uses three sets of units:
c
c       internal, are now mks International Standard Units
c       normalized, are used in solving the transport equations
c               and atomic physics routines,
c       external, are used for input and output,  and
c       standard, are used in computing diffusion coefficients and
c               other places where a physicist may want to "plug in"
c               his or her own physical model.
c
c       the quatities used in the code which have conversion factors are:
c       magnetic field strength, density, energy, temperature,
c       current, length, mass, and time.
c
c       in any units, the density unit is just the cubic length unit.
c       in addition, in internal units, 1 temperature = 1 energy unit/part.,
c       and 1 current unit is the current that would be flowing in a cylinder
c       one length unit in radius with a 1 magnetic field strength unit
c       poloidal field at the edge.  note that if the internal energy unit
c       is 1 ev, the internal temperature unit would be 2/3 ev.
c
c       standard units are:
c
c       magnetic field: gauss
c       density:        cm**-3
c       energy:         ergs
c       temperature:    ergs
c       current:        stat-amps
c       length:         centimeters
c       mass:           grams
c       time:           seconds
c
c       external units are:
c
c       magnetic field: kilo-gauss
c       density:        cm**-3
c       energy:         joules
c       temperature:    kev
c       current:        kilo-amperes
c       length:         centimeters
c       mass:           grams
c       time:           seconds
c
c       internal units will be referred to as:
c
c       magnetic field: ibu (internal b-units)
c       density:        idu (internal density units)
c       energy:         ieu (internal energy units)
c       temperature:    ihu ("internal heat units")
c       current:        iiu (internal intensity units)
c       length:         ilu ( etc.......)
c       mass:           imu
c       time:           itu
c
c       normalized units are:
c
c       magnetic field: tesla
c       density:        10^19 particles / m^3
c       energy:         kilo-joules
c       temperature:    kev
c       current:        kilo-amperes
c       length:         meters
c       mass:           kilograms
c       time:           seconds
c
c       analogously, abbreviations such as etu, smu, etc., will be thrown
c       in whenever necessary to make things less clear.
c
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      standard (Gaussian cgs) electricity and magnetism units
c       to fit the maxwell's equations:
c
c     D = E        ! in free space
c     H = B        ! in free space
c     curl ( H ) = ( 4 * pi / c ) * J  +  ( 1. / c ) * ( d D / d t )
c     curl ( E ) = - ( 1. / c ) * ( d B / d t )
c     div ( D ) = 4 * pi * charge density
c     div ( B ) = 0.
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
cl      internal electricity and magnetism units are chosen
c       fit the maxwell's equations:
c
c     curl ( B ) / emu0 = J + ( d E / d t ) * eps0
c     curl ( E )        = - d B / d t
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c bateman 14-may-85  constants in common block commhd
c
      pi = atan2 (0., -1.)
      twopi    = 2.0 * pi
      emu0  = pi * 4.e-7     ! volt * seconds / (amp * meter)
      eps0  = 8.854e-12      ! Coulombs / (volt * meter)
c
cend bateman
c
cl      1)      conversion from external to standard units
c
c
c       1 kilogauss = 1000 gauss
                                        uesb = 1000.0
c       1 joule = 10,000,000 ergs
                                        uese = 1.0e+07
c       1 watt = 1 joule / sec.
                                        uesp = uese
c       1 kev = 1000 ev = 1000 (#ergs/ev) ergs
                                        evs = cfev*10.0**cxev
                                        evsinv = 1.0 / evs
                                        uesh = 1000.0 * evs
c       1 kilo-ampere = c * 100 stat-amps
                                        uesi = fcc * 10.0**(2.0+fxc)
c       1 cm. = 1 cm.
                                        uesl = 1.0
c       1 gram = 1 gram
                                        uesm = 1.0
c       1 sec. = 1 sec.
                                        uest = 1.0
c       1 volt / cm. = 1 / (c * 10**-8) statvolts / cm.
                                        uesv = 10.0**(8.0-fxc) / fcc
c
c
cl      2)      conversion from internal to standard units
c
c bateman 14-may-85  change internal units to MKS ISU
c  Internal units will henceforth be MKS International Standard Units
c
c 1) magnetic field  1 tesla = 1.e4 Gauss
c 2) energy          1 Joule = 1.e7 ergs
c 3) length          1 meter = 1.e2 cm
c 4) time            1 second = 1 second
c
      uisb = 1.e4
      uise = 1.e7
      uisl = 100.
      uist = 1.
c
cend bateman
c
c               coefficient and exponent of 1 atomic unit in grams
c
        fcau = fcmp / fcap
        fxau = fxnucl
c
c
cl      3)      conversion factors dependent on previous factors
c
c
c       3.1)    density units are length units **-3.
c
        uesd = uesl**(-3)
        uisd = uisl**(-3)
c
cl      3.2)    1 internal temperature unit is 1 int. energy unit
c               per particle
c               hence, a temperature of 1 int. unit is
c               uise ergs/particle, which is  uise * (gamma - 1) ergs.
c
                gamin1  = 2.0 / 3.0
c
                uish = uise * gamin1
c
cl      3.3)    in internal units, emu0 * i = twopi * b * r in a cylinder
c               in radius with a poloidal field of b units.
c               in standard (cgs), i = (c/2) b * r
c
      uisi = uisb * uisl * 0.5 * fcc * 10.0**fxc * emu0 / twopi
c
        uisj = uisi / uisl**2
        uesj = uesi / uesl**2
c
        uisv = uisb * uisl / (fcc * uist) * 10.0**(-fxc)
c
        uisr = uisv / uisj
        uesr = uesv / uesj
c
c
cl      3.4)    1/2 mv**2 = e in internal units
c
        uism = uise * (uist/uisl)**2
c
cl      3.5)    power = e / t in internal units
c
        uisp = uise / uist
c
c
cl      4.1)    conversion from external to internal
c
        ueib = uesb / uisb
        ueid = uesd / uisd
        ueie = uese / uise
        ueih = uesh / uish
        ueii = uesi / uisi
        ueij = uesj / uisj
        ueil = uesl / uisl
        ueim = uesm / uism
        ueip = uesp / uisp
        ueir = uesr / uisr
        ueit = uest / uist
        ueiv = uesv / uisv
c
cl      4.2)    conversion from internal to external
c
        uieb = uisb / uesb
        uied = uisd / uesd
        uiee = uise / uese
        uieh = uish / uesh
        uiei = uisi / uesi
        uiej = uisj / uesj
        uiel = uisl / uesl
        uiem = uism / uesm
        uiep = uisp / uesp
        uier = uisr / uesr
        uiet = uist / uest
        uiev = uisv / uesv
c
cl      4.3)    conversion from standard to external
c
        useb = 1.0 / uesb
        used = 1.0 / uesd
        usee = 1.0 / uese
        useh = 1.0 / uesh
        usei = 1.0 / uesi
        usej = 1.0 / uesj
        usel = 1.0 / uesl
        usem = 1.0 / uesm
        usep = 1.0 / uesp
        user = 1.0 / uesr
        uset = 1.0 / uest
        usev = 1.0 / uesv
c
cl      4.4)    conversion from standard to internal
c
        usib = 1.0 / uisb
        usid = 1.0 / uisd
        usie = 1.0 / uise
        usih = 1.0 / uish
        usii = 1.0 / uisi
        usij = 1.0 / uisj
        usil = 1.0 / uisl
        usim = 1.0 / uism
        usip = 1.0 / uisp
        usir = 1.0 / uisr
        usit = 1.0 / uist
        usiv = 1.0 / uisv
c
c
cl      4.5)    conversion to normalized units (for chi only)
c
        uind = 1.e-19
        uine = 1.5 * uieh
c
        unid = 1.0 / uind
        unie = 1.0 / uine
c
        usnd = usid * uind
        usne = usie * uine
c
        unsd = 1.0 / usnd
        unse = 1.0 / usne
c
c
cl      5)      set up names of external units
c
c
        do 504 j = 1, 10
          lhb    = 'kilogauss '
          lhcden = 'kamp/sq cm'
          lhcurr = 'kiloampere'
          lhdenr = 'part/cc s '
          lhdens = 'part/cu cm'
          lheden = 'joul/cu cm'
          lhefld = 'volts / cm'
          lheflx = 'joul/cm2 s'
          lhener = '  joules  '
          lhlen  = '    cm    '
          lhnflx = 'part/cm2 s'
          lhpden = 'watt/cu cm'
          lhpowr = '  watts   '
          lhtemp = '   keV    '
          lhtime = '   secs.  '
  504   continue
c
c
c
        return
        end

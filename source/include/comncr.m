c@comncr.m /11040/bald89/wbaldn1 DCOMMON
c
cl                  c9.1     comncr--non-coronal equilibrium variables
c
c     version rhw 03/08/84
c
c  In comncr and comadp, the dimensions of the variables should match
c  the parameters set in setimp. Namely, "2" corresponds to kncimp 
c  (except for ncxx, where the dimension is governed by the number of
c  hydrogenic species). Likewise: "28" <-> kzmax.
c              
      common  /comncr/
     r     xnold(55,0:28,2), daqexp(2)       , vaqexp(2)       ,
     r     xn (55,0:28,2)  , ra (55,0:28,2)  , sa (55,0:28,2)  ,
     r     weir (55,28,2)  , wline (55,28,2) , wioni (55,28,2) ,
     r     wreko (55,28,2) , wchex (55,28,2) , wbrem (55,28,2) ,
     r     va (55,0:28,2)  , s0 (55,2)       , xns (55,2)      ,
     r     xnsold (55,2)   , da (55,0:28,2)  , dq (55,2)       ,
     r     xn0(55,2)       , xnas(55)        , aaa (55)        ,
     r     bbb (55)        , ccc (55)        , ddd (55)        ,
     r     eee (55)        , fff (55)        , pflxo (2)       ,
     r     plosso (2)      , psorco (2)      , pflx (2)        , 
     r     ploss (2)       , psorc (2)       , xntot(2)        , 
     r     xntoti (2)      , flout (2)       , flouto (2)      , 
     r     flscr (2)       , flscro (2)      , cimped(0:28,2)  ,
     r     dz(55,2)        , delts  , fcstep , xnerr  , recflx ,
     r     recscr , beta   ,
     i     nkimp(2)        , nk     , n2     , ncrept , nwcool ,
     i     nwline , nwioni , nwreko , nwchex , nwbrem , ifneo  
c
c
c  This common block contains variables used in the NC impurity
c  charge state transport code. In the following "*" denotes an
c  input variable. The default value is indicated at the end of
c  the description. All units are standard.
c
c  xn(j,k,i)     Impurity density of stage k of species i in zone j
c  xnold(j,k,i)  Same, but at previous time-step
c  xns(j,i)      Total ion density of impurity species i in zone j
c  xnsold(j,i)   Same, but at previous time-step
c  xnas(j)       Total ion density in zone j in previous substep
c  xn0(j,i)      Neutral density of impurity species i in zone j
c  ra(j,k,i)     Recombination rate, stage k -> k-1, species i, zone j
c  sa(j,k,i)     Ionization rate, stage k -> k+1, species i, zone j
c  s0(j,i)       Ionization rate of neutral impurity, species i, zone j
c  dq(j,i)       Total source of first stage of species i in zone j
c  dz(j,i)       Total source of last stage of species i in zone j
c  da(j,k,i)     Diffusion coefficient of species i, stage k, boundary j
c  va(j,k,i)     Drift velocity of species i, stage k, boundary j
c ...For the following, j indicates zone number and i indicates species...
c  wline(j,k,i)  Elec. cool. line rad'n. due to excitation of stage k-1 ions
c  wioni(j,k,i)  Electron cooling by ionization of stage k-1 ions
c  wreko(j,k,i)  Radiation due to recombination of stage k ions (no CX)
c  wchex(j,k,i)  Radiation due to charge exchange of stage k ions
c  wbrem(j,k,i)  Electron cooling radiation off stage k ions
c  weir(j,k,i)   Total electron cooling power due to stage k-1 ions
c  aaa(j)        Coefficient used in solving diffusion equation, zone j
c  bbb(j)             "        "   "    "        "        "       "   "
c  ccc(j)             "        "   "    "        "        "       "   "
c  ddd(j)             "        "   "    "        "        "       "   "
c  eee(j)             "        "   "    "        "        "       "   "
c  fff(j)             "        "   "    "        "        "       "   "
c  xntot(i)      Total number of impurity species i particles
c  xntoti(i)     Same, but at beginning of the run
c  pflx(i)       Total number of species i particles diffusing to the wall
c  pflxo(i)      Same, but at previous time-step
c  ploss(i)      Total number of species i ions lost to scrape-off
c  plosso(i)     Same, but at previous time-step
c  psorc(i)      Total source of species i particles
c  psorco(i)     Same, but at previous time-step
c  flout(i)      Number of species i ions diffusing to wall in a time-step
c  flouto(i)     Same, but at previous time-step
c  flscr(i)      Number of species i ions lost to scrap-off in a time-step
c  flscro(i)     Same, but at previous time-step
c  daqexp(i)     * Exponent for charge dependence of da, species i; 0.
c  vaqexp(i)     * Exponent for charge dependence of va, species i; 0.
c  cimped(k,i)   Pedestal boundary condition for stage k of species i
c  delts         Size of time-(sub)step
c  fcstep        Time-step factor, = 0.5 at beginning, otherwise = 1.0
c  xnerr         Greatest relative variation in a zone
c  recflx        * Recycling factor for diffusive losses to the wall; 1.0
c  recscr        * Recycling factor for scrape-off loss; 0.0
c  beta          * Multiplier for all reaction rates; 1.0
c  nk            Number of ionization stages of present impurity
c  nkimp(i)      Value of nk for impurity species i
c  n2            Number of BALDUR zones, = mzones
c  ncrept        The number of substeps in a time-step is ncrept / fcstep
c  ifneo         * If > 0, neoclassical diffusion is included; 0
c  nwcool        * Controls output for electron cooling rate; 0
c  nwline        * Controls output for line radiation; 0
c  nwioni        * Controls output for ionization losses; 0
c  nwreko        * Controls output for recombination radiation; 0
c  nwchex        * Controls output for charge exchange radiation; 0
c  nwbrem        * Controls output for bremsstrahlung radiation; 0
c
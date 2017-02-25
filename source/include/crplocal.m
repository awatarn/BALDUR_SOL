c@localp
c
c       this file contains the local common block to be used for transferring
c       local variables which are computed in mfile and nfile to the generated
c       subroutines (such as mf0001) which write into the output files
c
c       setimp contains parameters used in dimension statements
c
c       comsaw contains output from routine sawavg
c
        include 'cbparm.m'
c
        dimension zzni(idxchi),ztconf(idchp2,5),zloss(idxchi),zlne(6),
     1  zrconf(5)
        common/clocal/ zimrad(55),zpalfs(55),zpbeam(55),zpauxs(55),
     1  ztheat(55),isep,isepm1,
     2  zwapal,zwapbm,zwapau,zwapht,zwapoh,
     3  zmxtes,zmxtis,zmxnes,zmxtei,zmxtef,zmxtea,zmxtii,zmxtif,zmxtia,
     4  zmxjzi,zmxjzf,zmxjza
c
        common/comsaw/ tbsaws,sawte,sawne,sawdd,sawate,sawane,sawadd,
     1  sawtau,sawr1,sawri,sawrmx,sawqc,tmguir,njzq1,
     2  sawfus,sawalf,sawuth,estpol,uthavg,tauein,ohptot,
     3  tepre(55),tepost(55),teavg(55),tipre(55),tipost(55),tiavg(55),
     4  ajzpre(55),ajzpst(55),ajzavg(55),ohravg(55)
c
c--------1---------2---------3---------4---------5---------6---------7-c

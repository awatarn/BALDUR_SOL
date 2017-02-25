c@cbparm.m  20:00 27-nov-91
c                                     
      integer, parameter :: numimp=4
      integer, parameter :: idximp=numimp
      integer, parameter :: idxion=numimp+2
      integer, parameter :: idxchi=numimp+4
      integer, parameter :: idchp2=numimp+6
      integer, parameter :: id2chi=idxchi*idxchi
c
c  The following parameters control the size of the large local arrays
c  used in the NC code. If they are changed, the corresponding numbers
c  should be changed in cliche comncr and comadp. 
c
      integer, parameter :: kncimp=2
      integer, parameter :: kzmax=28
      integer, parameter :: knte=51
      integer, parameter :: knne=26
c
c  kncimp - maximum number of impurities in NC code
c  kzmax  - maximum Z of impurity (i.e., up to nickel)
c  knte   - number of temperature grid points used in rate interpolation
c  knne   - number of density grid points used in rate interpolation
c

c@ranz   .../baldur/code/util/ranz.f
c rgb 04-may-96 use ran1 in file random.f from Numerical Recipies
c
c  Random number function
c
      real function ranz( )
c
      data  idum / -1 /
      save idum
c
      ranz = ran2( idum )
c
cbate      write (*,*) ' ranz = ',ranz
c
      return
      end

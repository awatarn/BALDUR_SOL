c@sort1  .../baldur/code/util/sort1.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      subroutine sort1 (n, ra, rb )
c
c    Heap sort routine from Numerical Recipies (1986) p251
c
c  n     = number of elements in array to be sorted (input)
c  ra(j) = array to be sorted (input)
c  rb(j) = sorted array (output)
c
c  Note:  if ra(j) is not needed, ra and rb can share the same storage.
c
      implicit none
      integer   n, j, il, ir, i
      real      za,  ra(*), rb(*)
c
      if ( n .lt. 1 ) then
        write (*,*) 'abort: n .lt. 1 in sbrtn sort1.'
        stop
      endif
c
      do j=1,n
         rb(j) = ra(j)
      enddo
c
      if ( n .lt. 2 ) return
c
      il = n / 2 + 1
      ir = n
c
c   The index il will be decremented from its initial value down to 1
c  during the heap creation phase.  Once it reaches 1, the index ir will
c  be decremented from its initial value down to 1 during the 
c  heap selection phase.
c
 10   continue
      if ( il .gt. 1 ) then
         il = il - 1
         za = rb(il)
      else
         za = rb(ir)
         rb(ir) = rb(1)
         ir = ir - 1
         if ( ir .eq. 1 ) then
            rb(1) = za
            return
         endif
      endif
c
      i = il
      j = il + il
 20   continue
      if ( j .le. ir ) then
         if ( j .lt. ir ) then
            if ( rb(j) .lt. rb(j+1) ) j = j + 1
         endif
         if ( za .lt. rb(j) ) then
            rb(i) = rb(j)
            i = j
            j = j + 1
         else
            j = ir + 1
         endif
         go to 20
      endif
c
      rb(i) = za
c
      go to 10
c
      end

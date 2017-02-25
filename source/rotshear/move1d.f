      subroutine move1d(time, t, nt,  f, prof) 
c
      integer nt,nx,kxdim
      real time, t(nt), prof, f(nt)
c
      islice = 0
      if (nt .ne. 1) THEN
         do I = 1, nt-1
           if (time.ge.t(I) .and. time.le.t(i+1)) then
                 islice = I
                 if (t(i+1)-time .lt. time - t(i)) islice = i+1
           end if
         end do
c
         if (islice .eq.0) stop 'Time out of bounds'    
      else
          islice = 1
      end if

      prof = f(islice)

      return
      end
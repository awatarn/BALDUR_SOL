c@runtim  .../baldur/code/bald/daytim.f
c  rgb 07-feb-00 revised system calls for SGI IRIX
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
         subroutine runtim ( nout, cputime )
c
c u.12 update cpu time (secs) and print it
         real, intent(in out) ::  cputime
         real :: mclock
	   integer nout
c
cray         call second ( cputime )
         
         cputime = 0.01 * mclock() - cputime
c
         write(nout,9900) cputime
c
         return
 9900    format(5x,'cpu time used so far =',f10.2,' secs')
         end
c@daytim  .../baldur/code/bald/daytim.f
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
         subroutine daytim( nout )
c
c u.13 print date and time
c
!         character *36 :: nout
         character *10 ::  idate, itime, izone
         integer       ::  zval(8)
c
c  date is a cray fortran callable subroutine returning date.
c  clock is a cray fortran subroutine which returns the time of day.
c
c         call date ( idate )
cray         call clock ( itime )
cAIX      call clock_ ( itime )
c
         call date_and_time(idate, itime, izone, zval)
         write(nout,9900) (zval(i), i=1, 3), (zval(i), i=5, 7)
c
         return
 9900    format(5x,'date',4x,i4,2('-',i2.2),
     >      8x,'time',4x,2(i2.2,':'),i2.2)
         end

	real function mclock()
          integer :: c, cr, cm
          call system_clock(c, cr, cm)
 	  mclock=(100.*c)/cr
	end



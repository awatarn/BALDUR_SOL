c@cotrol.../baldur/code/bald/dsolver.f
c  ap  15-feb-00 changed Hollerith character representation into '...' 
c                in call mesage(...) and call error_olymp(...)
c  rgb 21-dec-93 changed call mesage(43h...) to call mesage('...')
c  rgb 14:28 11-nov-93 changed call abort to call abortb
c  rgb 20.26 16-jan-92 commented out section to resume a previous run
c   there were too many problems with common blocks redefined by
c   these routines and the restart capablility has not been used in years
c  rgb 18.83 11-jan-91 changed cliche olycom to cliche cbaldr throughout
c
         subroutine cotrol
c
c 0.3  control the run
c
c     version pa1
c
      include 'cparm.m'
      include 'cbaldr.m'
c-----------------------------------------------------------------------
c
       data iclass,isub/0,3/
       integer  :: ierr         ! error flag
       integer  :: idtmp ! fortran LUNs 
cap
         call mesage('    0.3 enter run control')
c
c-----------------------------------------------------------------------
c              1.         prologue
c
         if(nlres) go to 170
c
c                a.         new run
c
c     1.1      label the run
c
  110    call labrun
c
c     1.2      clear variables and arrays
c
 115     continue
c
 120    call clear
c
c
c     1.3      set default values
c
  130    call preset
c
c
c     1.4      define data specific to run
c
  140    call data
cap   initialization of plot output
!         select case (lplot)
!           case(1)
!             if (pgopen('/XTERM') <= 0) call 
!     >         mesage (' Device for graphical output isn''t accessible')
!           case default 
!             call mesage('Plot is disable.')
!             if (pgopen(' ') <= 0) call
!     >         mesage (' Device for graphical output isn''t accessible')
!         end select

c
c    1.5      set auxiliary values
c
  150    call auxval
c
c
c     1.6      define physical initial conditions
c
  160    call inital
c
         go to 180
c-----------------------------------------------------------------------
c
c                b.         resume a previous run
c
c                  1.7      pick up record, modify required parameters
  170    continue
c
c     label the continuation run
c
c     call labrun
c
c     clear variables and arrays
c
c        call clear
c
c     pick up record and print details
c
c         call resume
c
c     read any new data needed
c
c         call data
c
c     modify auxiliary variables as required
c
c         call auxval
c
         call abortb (6,
     & 'resume a previous run has been disabled in sbrtn cotrol')
c
c-----------------------------------------------------------------------
c
c                c.         preliminary operations
c
c                  1.8        start or restart the run
c
 180  continue
c
c bateman  initial MHD equilibrium computation
c
      call eqinit
c
cend bateman 23-feb-86
c
      call start
c
c     initial output
c
         call output(1)
         idtmp = 26  ! file id of scratch file
         if (ierr/=0) call bad_exit()
c
c-----------------------------------------------------------------------
c              2.         calculation
c
c     2.1      step on the calculation
c
  210    call stepon
c
c
c-----------------------------------------------------------------------
c              3.         output
c
c     3.1      periodic production of output
c
  310    call output(2)
         if (ierr/=0) call bad_exit()
c
c
c-----------------------------------------------------------------------
c              4.         epilogue
c
c     4.1      test for completion of run
c
  410    call tesend
c
         if(.not.nlend) go to 210
c
c     final output
c
         call output(3)
         if (ierr/=0) call bad_exit()
c
c
c     4.2      terminate the run
c
  420    call endrun
c
         return
         end

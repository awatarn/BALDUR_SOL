c@main  ../baldur/code/bald/main.f
c  rgb 24-nov-99 commented out open files for03 and for04
c  rgb 18-jan-92 changed label# to character strings
c  rgb 18-dec-91 for03 and for04 are unformatted files
c  rgb 28-nov-91 call link... converted to open (...)
c  rgb 09-jun-89 commented out the output to termimal from unit 5
c
c--------1---------2---------3---------4---------5---------6---------7-c
        program baldur
c--------:---------:---------:---------:---------:---------:---------:-c
c

      open (1, file='tidata',  status='unknown')
      open (2, file='jobxdat', status='old')
cbate      open (3, file='for03',   status='unknown',form='unformatted')
cbate      open (4, file='for04',   status='unknown',form='unformatted')
      open (5, file='short',   status='unknown')
      open (6, file='jobxlpt', status='unknown')
      open (22,file='for22',   status='unknown')
c
       call master
c
       call exit(1)
c
       stop 1
       end

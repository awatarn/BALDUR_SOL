c@labrun
c  rgb 29-may-96 move label5 from data statement to assignment statement
c--------1---------2---------3---------4---------5---------6---------7-c
c
       subroutine labrun
c
c 1.1  label the run
c                               
c---------------------------------------------------------------------
      include 'cparm.m'
      include 'cbaldr.m'
c---------------------------------------------------------------------
c
c     label5 contains version number, genealogy
c
      label5 = 'BALDPN S  21.01  1Aug02 12:00  '
c
c..VERSNO = version number
c
      versno = 20.33
c
c..read labels
c
         read (nread,9900) label1(1:48)
         read (nread,9900) label2(1:48)
         read (nread,9900) label3(1:48)
         read (nread,9900) label4(1:48)
c
c..write heading
c
         call blines(8)
         write (nout,9901)
         write (nout,9902)
         call blines(4)
c
c..write labels
c
         write (nout,9903) label5(1:8),label5(9:32)
         call blines(2)
         write (nout,*) label1(1:48)
         call blines(1)
         write (nout,*) label2(1:48)
         call blines(1)
         write (nout,*) label3(1:48)
         call blines(1)
         write (nout,*) label4(1:48)
         call blines(1)             
c
         ntychl=5
         write(ntychl,9903)(label5(1:48))
         write(ntychl,9900)(label1(1:48))
         write(ntychl,9900)(label2(1:48))
         write(ntychl,9900)(label3(1:48))
         write(ntychl,9900)(label4(1:48))
c
c
         return
c
 9900    format(a)
 9901    format(50x,'p r o g r a m   b a l d u r  ')
 9902    format(50x,27('*'))
 9903    format('   1-1/2-D BALDUR   ',a9,'   version ',1x,a24)
         end

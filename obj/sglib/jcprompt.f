	block data atcomblk
        include 'jcticcom'
        data lcrok /.false./, lechi /.false./
        data ltext /.false./
        data lmacro /-1/
        data ofile /' '/
        data cfeed /' '/
	end
        subroutine at_char(prompt,word)                         !04/05/85
c
c 07/06/98 CAL modify for f90
c              readonly,shared,   ->  action='read' via sed in Make.AIX
c              need  1 and VMS line to ship back to VMS
c              type -> at_msg
c 01/21/98 CAL modified for TRANSP
c 04/28/99 CAL general port
c
        include 'jcticcom'
        character word*(*),c*1,cn*1,prompt*(*)
        character*2 tabc 
        integer itab
        logical nil
        integer glyf
        character*130  zstr
	external atcomblk
        data itab/Z'09'/       ! TAB =HT = x09
        tabc(1:1) = char(itab)
        tabc(2:2) = '!'
        if(level.eq.0) call at__input(prompt,word)
        c=' '
        nil=.true.
        do while(nil)           !Loop til you get something
          kp=ip(level)
          if(kp.ge.ll) go to 200        !Cycle
          jp=kp                         !Find first non-blank
          do while (line(level)(jp:jp).eq.' ' .and. jp.le.ll)
            jp=jp+1
          enddo
          if(jp.le.ll) c=line(level)(jp:jp)     !Get character
          if(c.eq.'!') c=' '            !Ignore comments
          if(c.eq.' ') then
            kp=ll                       !Skip to end of line
            if(lcrok) nil=.false.       !Blank line ok
            go to 200
          endif
          ip(level)=jp                  !Start of glyph
          kp=jp+1                       !Pointer to end of glyph
          if(c.eq.''''.or.c.eq.'"') then        !quote string
            jp=kp                       !start of quote
            kp=index(line(level)(jp:ll),c)+jp-1 !Find matching quote
            if(kp.le.jp) then
              call at_msg('** need matching quote')
              c=' '
            elseif(kp.gt.jp) then
              word=line(level)(jp:kp-1)
            else
              word=' '                  !Allow null string
            endif
            kp=kp+1
            if(kp.le.ll .and. line(level)(kp:kp).eq.',') kp=kp+1
c
          elseif(c.eq.'@') then         !indirect file
            glyf=0
            do while(kp.lt.ll .and. glyf.eq.0)
              cn=line(level)(kp:kp)
              glyf=index(' !/;',cn)     !Glyph goes to blank,comment
              if(glyf.eq.0) kp=kp+1
            enddo
            ip(level)=kp                !save place
            if(kp-1.lt.jp+1) then
              call at_msg('** no indirect file after @')
              c=' '
            elseif(level.ge.8) then
              call at_msg('** more than 8 levels of indirect file')
              call at_reset             !Back to level 1
            else
              level=level+1             !Bounce up a level
              iun(level)=90+level
              open(iun(level),file=line(level-1)(jp+1:kp-1),
     >          action='read',   
     >          status='old',iostat=ios)
              if(ios.ne.0) then
                if(ios.eq.29) then
                call at_msg('** '//line(level-1)(jp+1:kp-1)//' missing')
                else
               call at_error(ios,'cant open '//line(level-1)(jp+1:kp-1))
c                   zstr = 'cant open '//line(level-1)(jp+1:kp-1)
                call at_error(ios,zstr)
                endif
                call at_reset           !Close all files
              else
                rewind iun(level)
              endif
              kp=ll             !Force reading of new line
              c=' '
            endif
c
          elseif(ltext) then            !Text string
            if(index('>/',c).ne.0) then         !check for file or exit
              call at_switch(word,line(level)(jp:),jp,ngot)
              c=word(1:1)
            endif
            word=line(level)(jp:)
            kp=jp+jc_trim(word)
c
          elseif(c.ne.',' .and. c.ne.' ') then  !normal string
            glyf=0
            do while(kp.lt.ll .and. glyf.eq.0)
              cn=line(level)(kp:kp)
              glyf=index(', !/;',cn)    !Glyph goes to blank,comma,comment
              if(glyf.eq.0) kp=kp+1
            enddo
            if(index('>/{;',c).ne.0) then       !output file or switch
              call at_switch(word,line(level)(jp:kp-1),jp,ngot)
              if(ngot.eq.0) then
                c=' '
              else
                c=word(1:1)
              endif
              if(c.eq.' ') goto 200     !skip upcase if comment
            else
              call at__upcase(word,line(level)(jp:kp-1))
            endif
            if(kp.le.ll .and. line(level)(kp:kp).eq.',') kp=kp+1
          endif
 200      ip(level)=kp
          if(c.ne.' ') nil=.false.
          if(kp.ge.ll) call at__input(prompt,word)
        enddo
c
        if(ofile.ne.' ' .and. iun(level).ne.0) then
          if(word.ne.' ') then
             write(90,'(a)',iostat=ierr)
     >		(word(1:jc_trim(word)))//tabc//prompt
c             zstr = word(1:jc_trim(word))//tabc//prompt
c             write(90,'(a)',iostat=ierr) zstr
          endif
        endif
        end
c
        subroutine at_switch(so,si,jptr,ngot)
        include 'jcticcom'
        character*(*) si,so
        integer ngot            !0=>null output, 1=>'so' ok
        character c*1,str*32
        character*80 etext
        logical setsw
        call at__upcase(str,si)         !copy to temp in case si a constant
        if(str(1:1).eq.';') str='/B'
        c=str(1:1)
        ls=jc_trim(str)
        setsw=.true.
        ngot=1          !dont change <so> unless we have something
        if(c.eq.'{' .and. str(ls:ls).eq.'}') then
          str(1:1)='/'
          call at__param(str,ls-1,so)
          return
        endif
        if(c.eq.'/' .and. ls.gt.4) then         !check for parameter
          call at__param(str,ls,so)
          return
        elseif(c.eq.'/' .and. ls.ge.2) then     !Check for switches
          c=str(2:2)
          if(c.eq.'N') then             !Strip off N or NO
            setsw=.false.
            c=str(3:3)
            if(c.eq.'O') c=str(4:4)
          endif
          is=index('ABCERSXZPH',c)
                !many switches must begin line
          if(index('CERSPH',c).ne.0 .and. jptr.gt.1) then
            ngot=0
            return
          endif
          if(is.gt.0) ngot=0            !most switches return nothing
          if(is.eq.10) then             !/H = help
            if(lhelp) then
              so=str                    !if enabled, return /H
              ngot=1
            else
              call at_msg('** help not enabled')
            endif
          elseif(is.eq.9) then          !/P = pause
            call at_pause
          elseif(is.eq.8) then          !/Z is no-op
          elseif(is.eq.7) then          !/X calls user_end
cx          call user_end
          elseif(is.eq.5) then          !/R = rewind
            if(level.gt.1) rewind iun(level)
          elseif(is.eq.4) then          !/E = echo user input
            lechi=setsw
            ngot=0
          elseif(is.eq.3) then          !/C = treat CR same as blank
            lcrok=setsw
          elseif(is.eq.2) then          !/B returns '/B'/
            so=str
            ngot=1
          elseif(is.eq.1) then          !/A aborts task
            call good_exit  !   
          else
            call at_msg('** invalid switch: '//str)
          endif
c
        elseif(c.eq.'/') then
          so=c
          ngot=1
        elseif(c.eq.'>') then
c* note: line is written after parameter is translated.
          ofile=si(2:)
          if(ofile.eq.' ') then
            close(90,iostat=ierr)
          else
            open(90,file=ofile,
     >        status='new',iostat=ierr)
            if(ierr.ne.0) then
            call at_error(ierr,'** cant open output file')
            ofile=' '
            endif
          endif
        endif
        end
        subroutine at__param(str,ns,so)
        character str*(*),so*(*)
        parameter (mtext=128)
        character texts*(mtext),values*(mtext)
        logical lstor
        character*130 zstr
!?put nt,nv in common or have entry to re-init
        integer nt, nv
        save nt,nv
	data nt/0/,nv/0/
        ls=ns           !local copy, length of string
        lq=index(str,'=')
        lt=ls           !length of text part
        if(lq.gt.0) lt=lq-1
        kt=1            !index in texts (/param)
        kv=1            !index in values (=value)
        lstor=.false.   !search for param in texts
        do while (kt.lt.nt .and. .not.lstor)
          knt=index(texts(kt:),char(0))+kt-1
          knv=index(values(kv:),char(0))+kv-1
cx        if(knt.le.kt) go to 200       !leave
          if(str(2:lt).eq.texts(kt:knt-1)) then
            lstor=.true.        !found it
          else
            kt=knt+1
            kv=knv+1
          endif
        enddo
        so=' '
        if(lstor .and. lq.eq.0) then    !/param returns value
          so=values(kv:knv-1)
        elseif(lstor) then              !redefine param
          values(kv:knv-1)=str(lq+1:ls)
        elseif(lq.gt.0) then            !not there, new param
          if(nt+lq-1.gt.mtext) then
            call at_msg('[too many parameters]')
            return
          endif
          nt=nt+1
          texts(nt:)=str(2:lt)  !store param text
          nt=nt+lt-1
          texts(nt:nt)=char(0)          !null terminate
          nv=nv+1
          values(nv:)=str(lq+1:ls)      !store param value
          nv=nv+ls-lq
          values(nv:nv)=char(0)
        else            !not there on get
	  call at_msg(' ** '//str(2:ls)//' undefined')
c           zstr = ' ** '//str(2:ls)//' undefined'
c          call at_msg(zstr)
        endif
        end
c
        subroutine at__input(prompt,word)       !Input a line           !!04/03/85
        character prompt*(*),word*(*)
        include 'jcticcom'              !8 lines, 8 pointers
        integer klass
        logical search
        integer idummy
	data klass /-1/
        if(level.le.1) then
          nw=jc_trim(word)              !Measure size of default word
          lp=jc_trim(prompt)            !Returns len>=1
        endif
        search=.true.
        do while (search)               !Loop until you get results
c
          ip(max(1,level))=1            !Point at beginning of line
          if(level.eq.0) then           !At start, input from command buffer
            call lib_get_foreign(line(1),narg,idumy)  !   
            search=narg.eq.0  !   
            iun(1)=-1                   !not a macro
            level=1
c                                       !Next try the input stream
          elseif(level.eq.1 .and. iun(1).ne.0) then
c                !lib_get_input doesnt do Fortran carriage control
            if(lp+nw.gt.72) then
               write(*,'(1x,a)') prompt(:lp)
               nww=min(nw,74)
               write(*,'(1x,a)') '['//word(:nww)//']? '
            else
               write(*,'(1x,a)') prompt(:lp)//'['//word(:nw)//']? '
            endif
            if(cfeed.ne.' ') then
              line(1)=cfeed
              ios=0
            else
              read(*,'(a)',iostat=ios) line(1)
            endif
            if(ios.eq.-1) line(1)='/A'  !abort on EOF
            search=ios.gt.0             !Try again on error
c
          else          !Input from an indirect file
            if(iun(level).ne.0) then    !Read if not macro
              read(iun(level),'(a)',iostat=ios) line(level)
              if(ios.ne.0) close(iun(level),iostat=ierr)
              if(ios.gt.0)              !error msg unless eof
     >          call at_error(ios,'error reading indirect file')
            else                        !macros are one line
              ios=-1
            endif
            if(ios.ne.0) then
              level=level-1             !Return to prior line
              lmacro=-1
            endif
            search=.false.
c?          if(nc.eq.0) search=.true.
            if(level.eq.0) search=.true.
          endif
          nc=jc_trim(line(level))
c         if(lechi) type *, line(level)(1:nc)
          if(lechi) call at_msg(line(level)(1:nc))
        enddo
c
        if(ip(level).eq.1) then
          do i=1,ll                     !Replace tabs,etc, with blanks
            if(line(level)(i:i).lt.' ') line(level)(i:i)=' '
          enddo
        endif
        end
c
        subroutine at_reset
        include 'jcticcom'
        do while (level.gt.1)
          if(iun(level).gt.0) close(iun(level))
          level=level-1
        enddo
        level=1
        iun(1)=-1
        line(1)=' '
        ip(1)=1
        ltext=.false.
        end
c
c...Force feed the tic input buffer
        subroutine at_text(prompt,word) !Get entire line        !06/20/88
        include 'jcticcom'
        character word*(*),prompt*(*)
        ltext=.true.
        call at_char(prompt,word)
        ltext=.false.
        end
c
c...return next glyph from input stream
        subroutine at_macro(mline)
        character*(*) mline
        include 'jcticcom'
        lmacro=level            !level to return to
        if(level.lt.8) level=level+1
        iun(level)=0
        line(level)=mline
        do i=1,ll               !replace tabs,etc, with blanks
          if(line(level)(i:i).lt.' ') line(level)(i:i)=' '
        enddo
        ip(level)=1
        end
c
        subroutine at_mend
        include 'jcticcom'
        if(lmacro.ge.0) then
cx        close(iun(level),iostat=ierr)
          level=lmacro
          lmacro=-1
        endif
        end
c
        subroutine at__upcase(stro,stri)
        character*(*) stro,stri
        integer quo
	integer hex5f
	data hex5f /Z'5f'/
        ls=len(stro)
        quo=0
        i=1
        do while (i.le.min(len(stri),ls))
          ich=ichar(stri(i:i))
          if(quo.ne.0) then
            if(ich.eq.quo) quo=0
          elseif(ich.eq.ichar('"') .or. ich.eq.ichar('''')) then
            quo=ich
          elseif(ich.ge.ichar('a') .and. ich.le.ichar('z')) then
C            ich=iand(ich,'5f'x)         !mask out '40'x
            ich=iand(ich,hex5f)         !mask out '40'x
          endif
          stro(i:i)=char(ich)
          i=i+1
        enddo
        if(i.le.ls) stro(i:)=' '
        end
c
        subroutine user_end             !dummy user_end         03/30/88
        call good_exit                    !   
        end
        subroutine at_stkln(nrem)       !Return length of "stack" line  !89
        include 'jcticcom'
        nrem=0
        if(level.eq.0) return           !before first input
        kp=ip(level)
        if(kp.ge.ll) return             !pointer past end of line
        jp=kp                   !skip over any blanks
        do while (line(level)(jp:jp).eq.' ' .and. jp.le.ll)
          jp=jp+1
        enddo
        if(line(level)(jp:jp).eq.'!') jp=ll+1   !ignore comments
        if(jp.le.ll) then
          nrem=jc_trim(line(level)(jp:))        !non-blank chars remaining
        endif
        end
        subroutine at_force(str)
        character*(*) str
        include 'jcticcom'
        cfeed=str
        end
        subroutine at_wipe
        include 'jcticcom'
        if(level.eq.0) return
        ip(level)=ll            !force next line
        end
        subroutine lib_get_foreign(line,narg,ll)  !   
        character line*(*)  !   
        line=' '  !   
        call get_arg_count(narg)  !   
        la=0  !   
        do i=1,narg  !   
          call get_arg(i,line(la+1:))  !   
          la=jc_len(line)+1  !   
        enddo  !   
        end  !   

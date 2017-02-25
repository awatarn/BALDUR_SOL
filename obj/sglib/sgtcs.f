c  sgtcs.for	sg terminal control routines (vax)
c
C  Dec-95       Ed Lazarus changes for labeling and post-script
C               enhancements. Put in repository by T.B. Terpstra
c
c  mar-82	change initt and finitt to service terminal types.
c		add termck to obtain the terminal type.
c		move iowait,tinput,toutpt to agiod.mar
c
c  aug-83	revise xycnvt, movea, drawa, pointa, dasha, agcsca and
c		 ivclip to do all clipping in screen units.
c
c  nov-84	update wincot to avoid integer overflow. [hht]
c
c  feb-85	change finitt to put mac terminal into alpha mode.  [hht]
c
c  jun-85	reorder subroutines; replace termck with ag_inio.  [met]
c
c  nov-86	add calls to aggsend and tclos_disk to finitt.
c		change initt not to call newpag if ibaud=0.
c		change finitt not to call movabs if args = 0,0.
c		add agsalf and agsgrf routines.  [met]
c
c  oct-87	update vecmod,dshmod,movabs,drwabs for 4105.
c		restore pntmod routine (for marker lines).
c		add tminit,linclr,cmdout,intout,intray from new plot-10.
c		add agmark,agtxrot,agtxclr,agpnclr.
c		translate terminal_4105 in initt and call tminit if yes.
c		handle 4105 screen switch; do all switches by calling by
c		 agsalf and agsgrf.
c		do 4105 text output in anstr; toutst,toutpt still 4014 only.
c		add agseet to access k4105,kfactr,kmaxsc. [met]
c
c  jul-88	call erase before reset in initt (no extra blank page on
c		 QMS printer); more initialization in reset.
c		change screen switch for (Versaterm) mac 4105 -> VT100.
c		add support for graph-on (terminal_type GO).
c		add iaghres,iaglres. [met]
c
c  nov-88	move check for log of zero/negative to wincot (from rescal).
c		convert source to lower case.
c		test bytes rather than bits in ivclip (faster).
c		fix length of comstr in agsgrf for 4025.
c		non-mac 4105 -> VT100 switch supports Intecolor.
c		replace recovx with alfmod in scursr.
c		purge typeahead (agtpurg) on first call to scursr after initt.
c		mac VT100 -> 4014 switch overrides 4105 in VTPRO. [met]
c
c  sep-89	on log of zero/neg, set to data min instead of 0 (log(1)).
c		always purge typeahead on entry to scursr.
c		enable broadcast trapping in agsgrf, disable in agsalf and
c		 type out any buffered messages. [met]
c
c  jan-90	put ksizef in common so returned ok in agseet.
c		tminit just set flag; initt do all initialization.
c		in 4105 mode, init text precision, line width.
c		add kscrn to common to avoid screen switch when already there.
c		flush tty buffer in agsgrf so user sees switch immediately.
c		save ksres in common so only do translate once.
c		move reset into initt; for 4014, call alfmod first (US),
c		 then erase, then chrsiz.
c		add XTERM, VS3100 terminal types. [met]
c
c  jan-91	add klcolr,ktcolr to common; always reset colors in initt.
c		init bypass cancel, error threshold, panel color in INIT.
c		redefine kterm: 2=sw dashes/pts, 3=hw; hw points in 4105.
c		remove kin,kcm (factors wrong in any case).
c		add MODGRF, DW4125 terminal types; remove ADM3M.
c		add agsmodl (see model). [met]
c
c  aug-91	add TK4107 terminal type; convert some subs to entries.
c		remove vwindo, recovx. [met]
c
c  feb-92	change size of csize 3/4 to fit new Versaterm.
c		only translate terminal_4105 if tminit not called.
c		revise agidds per MAC version; move into dshmod.
c		only set kkmode=4 for hardware dash, not software dash.
c		convert byte variables to integer; more subs to entries.
c		convert ainst from byte to character string.
c		add kmovef=-1 so move, then draw to same location = point.
c		comment out reset. [met]
c  sep-93	go to alfmod from vecmod before going to pntmod (some
c		  emulators think vecmod->pntmod is special point plot). [met]
c  aug-97	replace "abs(xi)" in wincot to speed execution.
c		revert to bit test in ivclip (faster on RISC).  [met]
c
C  jul-98       Ed Lazarus changes (from Dec. 1995) for labeling and post-script
C               enhancements. Put in repository by T.B. Terpstra
c               These changes are for use in Unix for TEK2PS stream file printing.
c               If the changes are not there, there is a font error after the first 
c               of printing.
c
c----------------------------------------------------------------/sgtcsblk
	block data sgtcsblk
	include 'sgtcs_inc'
	data k4105 /-1/
	data kfactr /4/	!in case initt not linked
	end
c----------------------------------------------------------------/initt
c * initialization routine
c
	subroutine initt(ibaud)
	include 'sgcopr_inc'
	include 'sgtcs_inc'
	character scres*4
	integer   icolmp(32)
	integer   istat
	integer   isgood
	external  sgtcsblk
	data icolmp /		!color map (?tek version has triplets)
     >	 0,000,00,000, 1,00,100,000,	!black, white
     >	 2,120,50,100, 3,240,50,100,	!red, green
     >	 4,000,50,100, 5,300,50,100,	!blue, cyan
     >	 6,060,50,100, 7,180,50,100/	!magenta, yellow
c
c * kkmode  0=alfmod (text), 1=vecmod (move,draw), 2=pntmod, 4=dshmod
c * kmovef  1=need a move (after text), 0=just did a move, -1=did a draw
c * kterm   <3=software dash lines/points, 3=hardware dashes/points
c
	kterm=2			!software dash lines
	kfactr=4		!=4 for 1024 screen size; =1 for 4096
c * set screen resolution if terminal_res is defined
	call usr_trnlog(' ','TERMINAL_RES',scres,lsres,isgood)
	if (isgood .eq. 1) then
	   if (lsres .eq. 4) then
	      read(scres(1:lsres),'(bn,i4)',iostat=ios) ksres
	   else if (lsres .eq. 3) then
	      read(scres(1:lsres),'(bn,i3)',iostat=ios) ksres
	   else
	      read(scres(1:lsres),'(bn,i5)',iostat=ios) ksres
	   end if
	   if(ksres.eq.4096) then
	      kfactr=1
	   endif
	endif
	kmxscr=4095/kfactr
c * set hardware dash line if terminal_dash is defined
	call usr_trnlog(' ','TERMINAL_DASH',scres,lsres,isgood)
	if (isgood .eq. 1) then
	   if(scres(1:lsres).eq.'YES') kterm=3
	endif
c * set 4105 mode if terminal_4105 is defined
	if(k4105.lt.0) then
	  call usr_trnlog(' ','TERMINAL_4105',scres,lsres,isgood)
	  k4105=0
	  if (isgood .eq. 1) then
	     if(scres(1:lsres).eq.'YES') k4105=1
	  endif
	endif
c * open output device(s)
	call ag_inio(' ',istat)	!get device type, etc.
	if(istat .eq. 0) call good_exit 
cx	call tmark_disk(0)	!back up over header
	call agsgrf		!switch to graph screen, output us (alfmod)
c
	kbeamx=-1		!initialize xycnvt
	kbeamy=-1
	do ii=1,4
	  kpchar(ii)=-1
	enddo
c * output 4105 header
	if(k4105.ne.0) then
	  call cmdout('KC')	!cancel (mode,queues)
c set error threshold to suppress most messages
	  call cmdout('KT')	!set error threshold level (no 4510)
	  call intout(3)	!for DW4125
c set bypass cancel for DW4125 GIN report
	  call cmdout('NU')	!set bypass cancel character
	  call intout(0)	!null = disable bypass
c set view window to default
	  kbeamx=0
	  kbeamy=0
	  call cmdout('RW')	!set view window
	  call xycnvt		!first corner = 0,0
	  call xycnvt		!second corner
c set a default color map (4105, 8 colors)
	  call cmdout('TG')	!set surface color map
	  call intout(1)	!surface number
	  call intray(32,icolmp)
c set graphtext writing mode to overstrike
cx	  call cmdout('MG')	!set graphtext writing mode (no 4510)
cx	  call intout(1)
c page clear resets text/line params in 4510 so do it now
	  if(ibaud.gt.0) call erase
c set graphtext rotation to 0
	  call agtxrot(0)
c set line width
	  call cmdout('MW')
	  call intout(0)	!1 pixel
c set graphtext precision to stroke
	  call cmdout('MQ')
	  call intout(2)	!1=string, 2=stroke
c set marker type to dot
	  call cmdout('MM')	!set marker type
	  call intout(0)
c set line index to white (black if on white screen)
	  klcolr=-1
	  call linclr(1)	!assumes black screen
c set text index to white
	  ktcolr=-1
	  call agtxclr(1)	!set text color
c set panel index to white
	  kpcolr=-256
	  call agpnclr(-256)	!set panel color
	else		!k4105.eq.0
c
	  if(ibaud.gt.0) call erase
	endif
c
cx	call dshmod(0)		!clear hardware dash line
	call chrsiz(1)		!set character size
	call tmark_disk		!end of header
c
	kdasht=0		!dash parameters
	kclip=0			!set default to clip
c
	tvsmx=1.		!transform factors
	tvsmy=1.
	tmintx=0.
	tminty=0.
c
	kminsx=0
	kmaxsx=0		!4095/kfactr
	kminsy=0
	kmaxsy=0		!3120/kfactr
	tminvx=0.
	tmaxvx=kmaxsx
	tminvy=0.
	tmaxvy=kmaxsy
c
	kxtype=1		!set default linear, linear
	kytype=1
c
	return
	end
c%
c----------------------------------------------------------------/tminit
c * set 4105 flag
c
	subroutine tminit(iterm)
	include 'sgtcs_inc'
	if (iterm.eq.1) then
	  k4105=1
	else
	  k4105=0
	endif
	return
	end
c%
c----------------------------------------------------------------/finitt
c * finish plotting
c
	subroutine finitt(ix,iy)
c
	call agsalf		!switch screen to alpha, dump buffer
	return
	end
c%
c----------------------------------------------------------------/reset
c * initialize tcs common
c
	subroutine reset
cx	include 'sgtcs_inc'
cxc
cx	kbeamx=-1		!initialize common
cx	kbeamy=-1
cx	do ii=1,4
cx	  kpchar(ii)=-1
cx	enddo
cxc
cx	call alfmod
cx	kdasht=0
cx	kclip=0			!set default to clip
cxc
cx	call chrsiz(1)		!set character size
cxc
cx	tvsmx=1.		!transform factors
cx	tvsmy=1.
cx	tmintx=0.
cx	tminty=0.
cxc
cx	kminsx=0
cx	kmaxsx=0	!4095/kfactr
cx	kminsy=0
cx	kmaxsy=0	!3120/kfactr
cx	tminvx=0.
cx	tmaxvx=kmaxsx
cx	tminvy=0.
cx	tmaxvy=kmaxsy
cxc
cx	kxtype=1		!set default linear, linear
cx	kytype=1
	return
	end
c%
c----------------------------------------------------------------/agsalf
c * switch to alpha screen for various terminal types
c
	subroutine agsalf
c
	include 'sgtcs_inc'
	character comstr*16
	character*120 msg
c
	if(kscrn.eq.0) return	!already on alpha
	call agfoff		!turn off file output
	if(kkmode.ne.0) call alfmod
c..device dependent actions
	call toutpt(24)		!reset v550
	if(k4105.ne.0) then	!4105 sequences
	  if(model.eq.'VS3100' .or. model.eq.'MAC') then
c				!leave visibility off
	  else
	    call cmdout('LV')	!set dialog area visibility
	    call intout(1)	!on
	  endif
	  if(model.eq.'MAC') then
	    call cmdout('LI')	!set dialog area index
	    call intout(1)	!so VT100 will become front window
	    call intout(0)
	    call intout(1)	!opaque
	    call cmdout('%!')	!select code
	    call intout(1)	!ansi
	  elseif(model.eq.'V550') then	!also IBM emulator
	    call cmdout('%!')	!select code
	    call intout(2)	!vt100
	  elseif(model.eq.'TK4107') then
	    call cmdout('%!')	!select code
	    call intout(0)	!dialog
	  else			!DW4125, ColorTrend
	    call cmdout('%!')	!select code
	    call intout(1)	!1=ansi; 2=edit (VT100)
	  endif
	else			!4014 sequences
	  if(model.eq.'MAC') then
	    call toutpt(24)	!extra can for versaterm
	  elseif(model.eq.'GO') then
	    call toutpt(0)	!null for graph-on
	  endif
	endif
c
c  these terminals dont support 4105 so always do 4014 sequence
c
	if (model.eq.'TK4025') then
	  comstr = '  '//comchr//'MON 34 H K'//char(13)	!return
	  do ii=1,len(comstr)
	    call toutpt(ichar(comstr(ii:ii)))
	  enddo
	elseif(model.eq.'VT240') then
	  comstr = char(27)//'[?38l'	!esc [ ? 3 8 l
	  do ii=1,6
	    call toutpt(ichar(comstr(ii:ii)))
	  enddo
	elseif(model.eq.'XTERM') then
	  comstr = char(27)//char(3)	!esc etx
	  do ii=1,2
	    call toutpt(ichar(comstr(ii:ii)))
	  enddo
	endif
	call agfonn
c
	call tsend		!dump the output buffer
c..type out any broadcast messages and deassign terminal (drop mailbox).
	lmsg=1
	do while(lmsg.gt.0)
	  call agtgbm(msg,lmsg)
	  if(lmsg.gt.0) write(*,'(1x,a)') msg(1:lmsg)
	enddo
	call tclos_term
c
	kscrn = 0		!on alpha screen
	return
	end
c%
c----------------------------------------------------------------/ixitblk
	block data ixitblk
	integer ixit_blok(4)
	integer ixit_cond
	common /ixit/ ixit_blok,ixit_cond
	data ixit_blok /0,0,1,0/
	end 
c----------------------------------------------------------------/agsgrf
c * switch to graphics screen for various terminal types
c
	subroutine agsgrf
c
	include 'sgtcs_inc'
	character*24 comstr
cdec$	psect	/ixit/ noshr
	integer usr_dclexh  
	integer ixit_blok(4)	
	integer ixit_cond
	common /ixit/ ixit_blok,ixit_cond
	external ixitblk
c
	if(kscrn.eq.1) return	!already on graph screen
	call agfoff		!turn off file output
c..declare exit handler
	if(ixit_blok(2).eq.0) then
c	  ixit_blok(2) = loc(agexit)   !   
c	  ixit_blok(4) = loc(ixit_cond)  !   
	  igood = usr_dclexh(ixit_blok)
	endif
c..device dependent actions
	if(k4105.ne.0) then	!4105 sequences
	  call cmdout('%!')	!select code
	  call intout(0)	!tek
	  if(model.eq.'VS3100') then
	  else
	    call cmdout('LV')	!set dialog area visibility
	    call intout(0)	!off
	  endif
	else			!4014 sequences
	  if(model.eq.'MAC') then
	    call cmdout('%!')	!force mac from 4105 to 4014
	    call intout(8)
	  endif
	endif
c
c  these terminals dont support 4105 so always do 4014 sequence
c
	if (model.eq.'TK4025') then
	  comstr = comchr//'WOR 33 H'//comchr//'GRA 1,35'//comchr//
     >		'SHR'//char(13)	!return
	  do ii=1,len(comstr)
	    call toutpt(ichar(comstr(ii:ii)))
	  enddo
c
	elseif(model.eq.'VT240' .or. model.eq.'XTERM') then
	  comstr = char(27)//'[?38h'
	  do ii=1,6
	    call toutpt(ichar(comstr(ii:ii)))
	  enddo
	endif
c
	call tsend_term		!dump terminal buffer
	call agfonn
c
	call alfmod		!if(kkmode.eq.1)
	kscrn = 1		!on graph screen
	return
	end
c%
c----------------------------------------------------------------/agexit
c * called by VMS on forced exit
c
	subroutine agexit(icond)
c
	call agfonn		!in case program was in upause
	call tback_term		!clear terminal buffer
	call agsalf		!back to alpha, release broadcast mbx
	call tsend_disk
	call tclos_disk		!close/delete disk file
	return
	end
c%
c======routines that do alpha output====================
c
c----------------------------------------------------------------/alfmod
c * set terminal to alphanumeric mode
c
	subroutine alfmod
	include 'sgtcs_inc'
c
	kkmode=0		!a/n mode
	call toutpt(31)		!us sets a/n mode
	return
	end
c%
c----------------------------------------------------------------/anmode
c * set terminal to alphanumeric mode and flush buffer
c
	subroutine anmode
	include 'sgtcs_inc'
c
	kkmode=0		!a/n mode
	call toutpt(31)		!us sets a/n mode
	call tsend_term		!flush tty output buffer
	return
	end
c%
c----------------------------------------------------------------/anstr
c * output a string stored as integers
c
	subroutine anstr(nchar,str) !   
	character*(*) str	!   
	include 'sgtcs_inc'
cx This routine should be converted to accept a character string argument.
cx However the Andrew Holland routines used by John Schivell expect an
cx integer array.
c
	call alfmod
	if(k4105.ne.0) then
	  call cmdout('LT')
	  call intout(nchar)
	endif
	do j=1,nchar
          call toutpt(ichar(str(j:j))) !   
	enddo
	kbeamx=kbeamx+(nchar*khorsz+kfactr/2)/kfactr
	return
	end
c%
c----------------------------------------------------------------/toutst
c * output text 1 char/word
c
	subroutine toutst(len,iade)
	include 'sgtcs_inc'
	integer iade(*)
c
	do i=1,len
	  call toutpt(iade(i))
	enddo
	return
	end
c%
c----------------------------------------------------------------/erase
c * erase page
c
	subroutine erase
	include 'sgtcs_inc'
c
	call toutpt(27)		!escape
	call toutpt(12)		!form feed
	call iowait(9)		!wait .9 seconds
c
	do ii=1,4
	  kpchar(ii)=-1
	enddo
	kmovef=1
	return
	end
c%
c----------------------------------------------------------------/newpag
c * erase page and move
c
	subroutine newpag_old
c
	call erase
	call movabs(0,767)		!top/left
	return
	end
c%
c----------------------------------------------------------------/newpag
c * erase page and move
c
	subroutine newpag
        character*12 usernm,pname*31,id*72,formats*72
C CAL 08/06/98 Cant do datebuf=fdate() on IBM
C        character*24 datebuf,fdate, one*2
        character*26 datebuf, one*2
        data init/1/
        if(init.gt.0)then
CCC      call get_arg(0,pname)
        call get_arg(0,id)
        l1=0
        l2=0
        do i=len(id),1,-1
        if(l2.eq.0.and.id(i:i).ne.' ')l2=i
        if(l1.eq.0.and.id(i:i).eq.'/')l1=i
        end do
        l2=min0(l2,len(pname))
        pname=id(l1+1:l2)
c         call getlog(usernm)
         call sget_user(usernm)
c         call fdate(datebuf)
         call C24DATE(datebuf)
         formats='(1x,a'
         write(one(1:2),fmt='(i2.2)')lenstr(pname)
         formats=formats(1:lenstr(formats))//one(1:2)
         formats=formats(1:lenstr(formats))//',14h generated by ,'
         write(one(1:2),fmt='(i2.2)')lenstr(usernm)
         formats=formats(1:lenstr(formats))//'a'//one(1:2)
         formats=formats(1:lenstr(formats))//',4h on ,a24,1h$)'
         write(id,fmt=formats)pname,usernm,datebuf
         init=0
        endif
        call csize(i,j)
	call erase
        call chrsiz(4)
        call agtmhs(0,0,id)
        if(j.eq.22)then
         call chrsiz(1)
        elseif(j.eq.20)then
         call chrsiz(2)
        elseif(j.eq.13)then
         call chrsiz(3)
        elseif(j.eq.12)then
         call chrsiz(4)
        else
         call chrsiz(1)
        endif
	call movabs(0,767)		!top/left
	return
	end
c%
        integer function lenstr(c)
        character*(*) c
        do i=len(c),1,-1
        lenstr=i
        if(c(i:i).ne.' ')return
        end do
        return
        end
c%
c----------------------------------------------------------------/bell
c * ring bell at terminal
c
	subroutine bell
c
	call toutpt(7)
	return
	end
c%
c----------------------------------------------------------------/ainst
c * read input string 4 chars/word
c
	subroutine ainst(nchar,string)
	character string*(*)
	data maxlen/72/
c
	len=min(nchar,maxlen)
	do ii=1,len
	  if(ii.gt.1) call agtpurg(.false.)
	  call tinput(kchar)
	  string(ii:ii)=char(kchar)
	enddo
	return
	end
c%
c----------------------------------------------------------------/tinstr
c * read input characters one/word
c
	subroutine tinstr(nchar,iade)
	integer iade(*)
c
	do i=1,nchar
	  call tinput(kchar)
	  iade(i)=kchar
	enddo
	return
	end
c%
c----------------------------------------------------------------/scursr
c * get cursor position in screen units
c
	subroutine scursr(ichr,ix,iy)
	include 'sgtcs_inc'
	character instr*8
c	integer*2 ini(6)
	integer ini(6)
	integer   hex1f
	data hex1f /z'1f'/
c
c * output (esc) (sub) to turn on cursor
	call agfoff		!turn off file output
	call toutpt(27)		!escape
	call toutpt(26)		!sub
	call tsend
c * flush any type-ahead
cx	call agtpurg(.true.)
c * determine how many chars to input
	nchar=5
cx	if(model.eq.'GO') nchar=5	!graph-on
	if(model.eq.'TK4025') nchar=5
	if(model.eq.'XTERM') nchar=5
	if(model.eq.'DW4125') nchar=6
	if(model.eq.'MAC')  nchar=6
	if(model.eq.'MODGRF') nchar=6
	if(model.eq.'TK4107') nchar=6	!??
	if(model.eq.'VS3100') nchar=6
	if(model.eq.'VT240') nchar=6	!if setup is GIN term=CR
	if(model.eq.'XTC') nchar=6
	if(model.eq.'V550') nchar=7
c
	call ainst(nchar,instr)
c * restore the terminal status
	call agfonn
	call alfmod
	ichr=ichar(instr(1:1))
c * convert coordinates to integer
	do ii=2,5
	  ini(ii)=iand(ichar(instr(ii:ii)),hex1f)
	enddo
c * decode screen co-ordinates
	ix=int(ini(2))*32+ini(3)
	iy=int(ini(4))*32+ini(5)
c * apply screen scale factor
	ix=ix*(4/kfactr)
	iy=iy*(4/kfactr)
	return
	end
c%
c----------------------------------------------------------------/dcursr
c * read graphics cursor in screen coordinates
c
	subroutine dcursr(ichar,ix,iy)
c
	call scursr(ichar,ix,iy)
	return
	end
c%
c----------------------------------------------------------------/vcursr
c * get cursor position in virtual units
c
	subroutine vcursr(ichar,x,y)
c
	call scursr(ichar,ixa,iya)
	call revcot(ixa,iya,x,y)
	return
	end
c%
c----------------------------------------------------------------/hdcopy
c * make a hardcopy of the screen
c
	subroutine hdcopy
c
	call agfoff		!turn off file output
	call toutpt(27)		!escape
	call toutpt(23)		!etb
	call tsend_term		!flush output buffer
	call iowait(180)	!wait 18 seconds
	call agfonn
	return
	end
c%
c======routines that do graphic output====================
c
c----------------------------------------------------------------/vecmod
c * set terminal to graphics mode
c
	subroutine vecmod
	include 'sgtcs_inc'
c
	if(kkmode.ne.1) then
c * output (us) to enter a/n mode and reset for vector mode
	  if(kkmode.ne.0) call alfmod
	  do ii=1,4
	    kpchar(ii)=-1
	  enddo
	  kkmode=1	  	!vector mode
	endif
c * output (gs) to enter vector mode
	call toutpt(29)		!send gs
	kmovef=1		!need a move
	return
	end
c%
c----------------------------------------------------------------/movea
c * move beam (virtual, absolute)
c
	subroutine movea(x,y)
	include 'sgtcs_inc'
c
	trealx=x		!update unclipped coords
	trealy=y
	go to 100
c
c----------------------------------------------------------------/mover
c * move beam (virtual, relative)
c
	entry mover(dx,dy)
c
	trealx=trealx+dx
	trealy=trealy+dy
 100	continue		!convert to screen coord
	call wincot(trealx,trealy,krealx,krealy)
	if(kclip.eq.0) then	!quit if outside
	  if(ivclip(krealx,krealy).ne.0) return
	endif
	call movabs(krealx,krealy)	!move
	return
	end
c%
c----------------------------------------------------------------/movabs
c * move beam (screen, absolute)
c
	subroutine movabs(ix,iy)
	include 'sgtcs_inc'
c
	ibx=ix
	iby=iy
	go to 100
c
c----------------------------------------------------------------/movrel
c * move beam (screen, relative)
c
	entry movrel(idx,idy)
c
	ibx=kbeamx+idx
	iby=kbeamy+idy
 100	continue
	if(kdasht.ne.0) call dshmod(0)	!clear dash line type
	if(kmovef.le.0) then	!last output move or draw
	  if(kbeamx.eq.ibx .and. kbeamy.eq.iby) return !already here
	endif
	call vecmod
	kbeamx=ibx		!update coordinates
	kbeamy=iby
	call xycnvt		! and send them
	return
	end
c%
c----------------------------------------------------------------/drawa
c * draw line (virtual, relative)
c
	subroutine drawa(x,y)
	include 'sgtcs_inc'
c
	trealx=x		!update unclipped coords
	trealy=y
	go to 100
c
c----------------------------------------------------------------/drawr
c * draw line (virtual, relative)
c
	entry drawr(dx,dy)
c
	trealx=trealx+dx
	trealy=trealy+dy
 100	continue
	itx=krealx		!get previous point
	ity=krealy
	call wincot(trealx,trealy,ix,iy)	!convert this point
	krealx=ix		!save unclipped coords
	krealy=iy
	if(kclip.eq.0) then	!if clip option
	  iclast=ivclip(itx,ity)	!check previous point
	  call agcsca(itx,ity,ix,iy,lout)	!move either/both
	  if(lout.ne.0) return		!quit if both outside
c * if last was outside, move to clipped coordinates
	  if(iclast.ne.0) call movabs(itx,ity)
	endif
c * draw to new coordinates
	call drwabs(ix,iy)
	return
	end
c%
c----------------------------------------------------------------/drwabs
c * draw line (screen, absolute)
c
	subroutine drwabs(ix,iy)
	include 'sgtcs_inc'
c
	ibx=ix
	iby=iy
	go to 100
c
c----------------------------------------------------------------/drwrel
c * draw line (screen, relative)
c
	entry drwrel(idx,idy)
c
	ibx=kbeamx+idx
	iby=kbeamy+idy
 100	continue
	if(kkmode.ne.1) call vecmod
	if(kmovef.eq.-1) then	!last output draw
	  if(kbeamx.eq.ibx .and. kbeamy.eq.iby) return !already here
	endif
	if(kmovef.eq.1) call xycnvt
	kbeamx=ibx		!update coordinates
	kbeamy=iby
	call xycnvt		!send coordinates
	kmovef=-1		!last output draw
	return
	end
c%
c----------------------------------------------------------------/pntmod
c * set terminal to point mode
c
	subroutine pntmod
	include 'sgtcs_inc'
c
	if(kkmode.eq.4) call alfmod
	if(kkmode.eq.1) call alfmod
	do ii=1,4
	  kpchar(ii)=-1
	enddo
	kkmode=2		!point mode
	if(kterm.ge.3 .or. k4105.gt.0) then
c * output (fs) to enter point mode
	  call toutpt(28)	!set hardware point
cx	  kpoins=-1
	  kmovef=1
	endif
c
	return
	end
c%
c----------------------------------------------------------------/pointa
c * put point (virtual, absolute)
c
	subroutine pointa(x,y)
	include 'sgtcs_inc'
c
	trealx=x		!update unclipped coords
	trealy=y
	go to 100
c
c----------------------------------------------------------------/pointr
c * put point (virtual, relative)
c
	entry pointr(dx,dy)
c
	trealx=trealx+dx
	trealy=trealy+dy
 100	continue
	call wincot(trealx,trealy,krealx,krealy)
	if(kclip.eq.0) then
	  if(ivclip(krealx,krealy).ne.0) return
	endif
 	call pntabs(krealx,krealy)
	return
	end
c%
c----------------------------------------------------------------/pntabs
c * put point (screen, absolute)
c
	subroutine pntabs(ix,iy)
	include 'sgtcs_inc'
c
	ibx=ix
	iby=iy
	go to 100
c
c----------------------------------------------------------------/pntrel
c * put point (screen, relative)
c
	entry pntrel(idx,idy)
c
	ibx=kbeamx+idx
	iby=kbeamy+idy
 100	continue
	if(kkmode.ne.2) call pntmod
	if(kmovef.le.0) then	!last output move or draw
	  if(kbeamx.eq.ibx .and. kbeamy.eq.iby) return !already here
	endif
	kbeamx=ibx
	kbeamy=iby
	if(kterm.lt.3 .and. k4105.le.0) then
	  call toutpt(29)	!move to point (gs)
	  call xycnvt
	endif
cx	do ii=1,4
cx	  kpchar(ii)=-1
cx	enddo
	call xycnvt		!draw to same to make a point
	kmovef=-1
	return
	end
c%
c----------------------------------------------------------------/dshmod
c * set terminal to dash mode
c
	subroutine dshmod(ldash)
	include 'sgtcs_inc'
	integer itable(4)
	data itable/ 5,10,25,50/
c
	if(kkmode.eq.2) call alfmod	!point mode -> alpha
	if(ldash.eq.0) then		!if no dash
	  if(kdasht.ne.0) then
	  call toutpt(27)		!clear hardware dash
	  call toutpt(96)
	  kdasht=0
	  endif
cx	  call vecmod
	  return
	endif
	if(kkmode.eq.0) then		!if a/n mode
	  call toutpt(29)		!vector mode
	  call xycnvt			!move to this spot
	endif
c
	if(kterm.ge.3 .and. ldash.gt.0 .and. ldash.le.7) then
	  kkmode=4			!dash mode
	  call toutpt(27)		!set hardware dash
	  kdasht=ldash
	  call toutpt(96+kdasht)	!dash code
	  kdashs=-1		!flag for hardware dash
	  kmovef=1
	  return
	endif
c
	if(ldash.ne.kdasht) then
c * determine dash type, then decode
	  idash = ldash
	  if(idash.lt.0) idash=1
	  kdasht=idash
	  if(idash.ge.11) go to 500
	  kmdash(0)=1		!draw
	  kmdash(1)=0		!move
	  kmdash(2)=1		!draw
	  kmdash(3)=0		!move
	  go to (100,200,300,400,100,100,100,100,1000,100),idash
c * dot
  100	  idash=1
	  kdash(0)=1		!dot
	  kdash(1)=10		!space
	  kdash(2)=1		!dot
	  kdash(3)=10		!space
	  go to 1000
c * dot dash
  200	  kdash(0)=1		!dot
	  kdash(1)=10		!space
	  kdash(2)=15		!dash
	  kdash(3)=10		!space
	  go to 1000
c * short dash
  300	  kdash(0)=15		!dash
	  kdash(1)=15		!space
	  kdash(2)=15		!dash
	  kdash(3)=15		!space
	  go to 1000
c * long dash
  400	  kdash(0)=35		!dash
	  kdash(1)=10		!space
	  kdash(2)=35		!dash
	  kdash(3)=10		!space
	  go to 1000
c * decode integer dash specification string
  500	  irem=idash		!process digits from
	  do j=3,0,-1		! right to left (low to high)
	    if(irem.eq.0) then	!duplicate earlier pattern
	      kdash(j)=kdash(j+2)
	      kmdash(j)=kmdash(j+2)
	    else
	      isub=mod(irem,10)
	      irem=irem/10
	      kdash(j)=itable((isub+1)/2)	!dash length
	      kmdash(j)=iand(isub,1)	!even=move, odd=draw
	    endif
	  enddo
c
 1000	  continue
	  kdashs=3		!init dash step and
	  tdashr=0.		! remaining dash length
	endif
	return
	end
c%
c----------------------------------------------------------------/dasha
c * draw dash line (virtual, absolute)
c
	subroutine dasha(x,y,ldash)
	include 'sgtcs_inc'
c
 	trealx=x		!update unclipped coords
	trealy=y
	go to 100
c
c----------------------------------------------------------------/dashr
c * draw dash line (virtual, relative)
c
	entry dashr(dx,dy,ldash)
c
	trealx=trealx+dx
	trealy=trealy+dy
 100	continue
	itx=krealx		!get previous point
	ity=krealy
	call wincot(trealx,trealy,ix,iy)	!convert this point
	krealx=ix		!save unclipped coords
	krealy=iy
	if(kclip.eq.0) then	!if clip option
	  iclast=ivclip(itx,ity)	!check previous point
	  call agcsca(itx,ity,ix,iy,lout)	!move either/both
	  if(lout.ne.0) return		!quit if both outside
c * if last was outside, move to clipped coordinates
	  if(iclast.ne.0) call movabs(itx,ity)
	endif
c * draw to new coordinates
 	call dshabs(ix,iy,ldash)
	return
	end
c%
c----------------------------------------------------------------/dshabs
c * draw dash line (screen, absolute)
c
	subroutine dshabs(ix,iy,ldash)
	include 'sgtcs_inc'
c
	ibx=ix
	iby=iy
	go to 100
c
c----------------------------------------------------------------/dshrel
c * draw dash line (screen, relative)
c
	entry dshrel(idx,idy,ldash)
c
	ibx=kbeamx+idx
	iby=kbeamy+idy
 100	continue
	if(ldash.eq.0) then	!if no dash pattern
	  call drwabs(ibx,iby)	!draw solid line
	  return
	endif
	if(ldash.ne.kdasht) call dshmod(ldash)	!init dashed vector
	if(kdashs.lt.0) then	!if hardware dash
	  kkmode=1
	  call drwabs(ibx,iby)	!draw line
	  kkmode=4
	  return
	endif
	if(ldash.eq.9) then	!alternate move,draw
	  if(kdashs.eq.0) then
	    kdashs=1		!set to draw next time
	    call movabs(ibx,iby)	!move
	  else
	    kdashs=0		!set to move next time
	    call drwabs(ibx,iby)	!draw
	  endif
	  return
	endif
c dashed vector
	dx=ibx-kbeamx
	dy=iby-kbeamy
	d0=sqrt(dx*dx+dy*dy)		!distance
	kkmode=1
	d=d0
	do while(d.ne.0.)		!loop until distance exhausted
c
	  if(tdashr.eq.0.) then		!get next
	    kdashs=mod(kdashs+1,4)	! move or draw and its
	    tdashr=kdash(kdashs)	! dash length from dash spec
	  endif
c
	  if(d.ge.tdashr) then		!remaining
	    d=d-tdashr			! distance to point
	    tdashr=0.			! use total dash
	  else
	    tdashr=tdashr-d		!unused dash fragment
	    d=0.			! with no distance left
	  endif
c
	  dod0=d/d0			!with fraction
	  jx=ibx-ifix(dod0*dx)		! of distance being drawn
	  jy=iby-ifix(dod0*dy)		! update coordinates
	  if(kmdash(kdashs).eq.0) then
	    call toutpt(29)		!move
	    kbeamy=jy
	    kbeamx=jx
	    call xycnvt
	  else
	    call drwabs(jx,jy)		!draw
	  endif
	enddo
	return
	end
c%
c----------------------------------------------------------------/linclr
c * set the line color
c
	subroutine linclr(index)
	include 'sgtcs_inc'
	if(index.eq.klcolr) return
	klcolr=index
	if(k4105.ne.0 .and. index.ge.0) then
	  call cmdout('ML')
	  call intout(index)
	endif
	return
	end
c%
c----------------------------------------------------------------/agpnclr
c * set the panel color (-7 to 0 are solid colors)
c
	subroutine agpnclr(index)
	include 'sgtcs_inc'
	if(index.eq.kpcolr) return
	kpcolr=index
	if(k4105.ne.0 .and. index.ge.-255) then
	  call cmdout('MP')
	  call intout(index)
	endif
	return
	end
c%
c----------------------------------------------------------------/agtxclr
c * set the text color
c
	subroutine agtxclr(index)
	include 'sgtcs_inc'
	if(index.eq.ktcolr) return
	ktcolr=index
	if(k4105.ne.0 .and. index.ge.0) then
	  call cmdout('MT')
	  call intout(index)
	endif
	return
	end
c%
c----------------------------------------------------------------/agtxrot
c * set the text rotation
c
	subroutine agtxrot(iangl)
	include 'sgtcs_inc'
	if(k4105.ne.0) then
	call cmdout('MR')
	call intout(iangl)
	call intout(0)
	endif
	return
	end
c%
c----------------------------------------------------------------/agmark
c * output a marker
c
	subroutine agmark(lmark)
c
	include 'sgtcs_inc'
	if(k4105.eq.0) return
	call cmdout('MM')
	call intout(lmark)
	return
	end
c%
c======graphic support and clipping routines==================
c
c----------------------------------------------------------------/vbgpnl
c * begin panel with virtual coordinates
c
	subroutine vbgpnl(x,y,ibound)
	include 'sgtcs_inc'
	call cmdout('LP')
	call wincot(x,y,kbeamx,kbeamy)
	call xycnvt
	call intout(ibound)
	return
	end
c
c----------------------------------------------------------------/begpnl
c * begin panel with screen coordinates
c
	subroutine begpnl(ix,iy,ibound)
	include 'sgtcs_inc'
	kbeamx=ix
	kbeamy=iy
c
c * begin panel at screen location
	entry agbgpnl(ibound)
	call cmdout('LP')
	call xycnvt
	call intout(ibound)
	return
	end
c
c----------------------------------------------------------------/endpnl
c * end panel
c
	subroutine endpnl
	call cmdout('LE')
	return
	end
c
c----------------------------------------------------------------/cmdout
c * output a 4105 command (esc,char,char)
c
	subroutine cmdout(cmd)
	character cmd*2
	call toutpt(27)
	call toutpt(ichar(cmd(1:1)))
	call toutpt(ichar(cmd(2:2)))
	return
	end
c
c----------------------------------------------------------------/intray
c * output an array of 4105 command integers
c
	subroutine intray(len,iray)
	integer iray(*)
	call intout(len)
	if (len.le.0) return
	do i=1,len
	  call intout(iray(i))
	enddo
	return
	end
c%
c----------------------------------------------------------------/intout
c * output an integer as part of a 4105 command
c
	subroutine intout(int)
c * find the absolute value
	jint=abs(int)
c * compute the two hi-i is and the lo-i
	jhi1=jint/1024+64
	jhi2=mod(jint/16,64)+64
	jloi=mod(jint,16)+32
	if (int.ge.0) jloi=jloi+16
c * see if hi-i is needed
	if (jhi1.ne.64) go to 10
	if (jhi2.ne.64) go to 30
	go to 50
 10	call toutpt(jhi1)
c * insert second hi-i
 30	call toutpt(jhi2)
c * insert lo-i
 50	call toutpt(jloi)
c * send the string
	return
	end
c%
c----------------------------------------------------------------/iaghres
c * convert 1024 coordinate to 4096 resolution
c
	function iaghres(ixy)
	include 'sgtcs_inc'
c
	iaghres=ixy*4/kfactr
	return
	end
c%
c----------------------------------------------------------------/iaglres
c * convert 4096 resolution to 1024 coordinates
c
	function iaglres(ixy)
	include 'sgtcs_inc'
c
	iaglres=ixy*kfactr/4
	return
	end
c%
c----------------------------------------------------------------/agcsca
c * determine if need to clip
c
	subroutine agcsca(ix1,iy1,ix2,iy2,lout)
	include 'sgtcs_inc'
	logical swap
	integer ic1,ic2
c	parameter (isleft ='01'x)	!bit mask to move left
c	parameter (isright='02'x)	! move right
c	parameter (isbelow='04'x)	! move down
c	parameter (isabove='08'x)	! move up
	parameter (isleft =1)	!bit mask to move left
	parameter (isright=2)	! move right
	parameter (isbelow=4)	! move down
	parameter (isabove=8)	! move up
c  Cohen/Sutherland clipping algorithm.  Newman and Sproull,
c  Principles of Interactive Computer Graphics.  pp 123-124.
c
	swap=.false.
	do while (.true.)
	  ic2=ivclip(ix2,iy2)		!compute clip code for
	  ic1=ivclip(ix1,iy1)		! each endpoint
c  line entirely inside
	  if(ic1.eq.0 .and. ic2.eq.0) go to 1000
c  line entirely outside
	  if(iand(ic1,ic2).ne.0) then
	    lout=1
	    return
	  endif
c
	  if(ic1.eq.0) then
	    swap=.not.swap
	    ic1=ic2
	    ic2=0
	    iz=ix1
	    ix1=ix2
	    ix2=iz
	    iz=iy1
	    iy1=iy2
	    iy2=iz
	  endif
c
	  if(iand(ic1,isleft ).ne.0) then	!push toward left edge
	    iy1=iy1+(iy2-iy1)*(kminsx-ix1)/(ix2-ix1)
	    ix1=kminsx
	  elseif(iand(ic1,isright).ne.0) then	!push toward right edge
	    iy1=iy1+(iy2-iy1)*(kmaxsx-ix1)/(ix2-ix1)
	    ix1=kmaxsx
	  elseif(iand(ic1,isbelow).ne.0) then	!push toward bottom edge
	    ix1=ix1+(ix2-ix1)*(kminsy-iy1)/(iy2-iy1)
	    iy1=kminsy
	  elseif(iand(ic1,isabove).ne.0) then	!push toward top edge
	    ix1=ix1+(ix2-ix1)*(kmaxsy-iy1)/(iy2-iy1)
	    iy1=kmaxsy
	  endif
	enddo
c
 1000	if(swap) then
	  iz=ix1
	  ix1=ix2
	  ix2=iz
	  iz=iy1
	  iy1=iy2
	  iy2=iz
	endif
	lout=0
	return
	end
c%
c----------------------------------------------------------------/ivclip
c  inout=ivclip(ix,iy)
c
	function ivclip(ix,iy)
	include 'sgtcs_inc'
	integer ivcl
	parameter (isleft =1)	!bit mask to move left
	parameter (isright=2)	! move right
	parameter (isbelow=4)	! move down
	parameter (isabove=8)	! move up
c
	ivcl=0			!inside
	if(ix.lt.kminsx) then
	  ivcl=ior(ivcl,isleft )	!left
	elseif(ix.gt.kmaxsx) then
	  ivcl=ior(ivcl,isright)	!right
	endif
	if(iy.lt.kminsy) then
	  ivcl=ior(ivcl,isbelow)	!below
	elseif(iy.gt.kmaxsy) then
	  ivcl=ior(ivcl,isabove)	!above
	endif
	ivclip=ivcl
	return
	end
c%
c----------------------------------------------------------------/wincot
c * transform from virtual window to screen coordinates
c
	subroutine wincot(x,y,ix,iy)
	include 'sgtcs_inc'
c
	xx=x
	yy=y
	if(kxtype.eq.2) then
	  if(xx.le.0.) then
	    xx = log10(tminvx)	!avoid logzerneg problems
	  else
	    xx = log10(xx)
	  endif
	endif
	if(kytype.eq.2) then
	  if(yy.le.0.) then
	    yy = log10(tminvy)	!avoid logzerneg problems
	  else
	    yy = log10(yy)
	  endif
	endif
c transform to screen coordinates
	xi = xx*tvsmx+tvsbx
	yi = yy*tvsmy+tvsby
cx	if(abs(xi).gt.3.e4) xi=sign(3.e4,xi)
cx	if(abs(yi).gt.3.e4) yi=sign(3.e4,yi)
	xi = min(xi,+3.e4)		!avoid float->int error
	xi = max(xi,-3.e4)
	yi = min(yi,+3.e4)
	yi = max(yi,-3.e4)
	ix = xi
	iy = yi
	return
	end
c%
c----------------------------------------------------------------/revcot
c * transform from screen to virtual window coordinates
c
	subroutine revcot(ix,iy,x,y)
	include 'sgtcs_inc'
c
	x=float(ix-kminsx)/tvsmx+tmintx
	y=float(iy-kminsy)/tvsmy+tminty
	if(kxtype.eq.2) x=10**x
	if(kytype.eq.2) y=10**y
	return
	end
c%
c----------------------------------------------------------------/xycnvt
c * convert x/y coordinates to tektronix beam commands
c
	subroutine xycnvt
	include 'sgtcs_inc'
c  see 4010 maintenance manual (6-33) graphic vector description
	integer hix,hiy,lsb,lox,loy
	integer hex1f, hex03
	data hex1f,hex03 /z'1f', z'03'/
c
c * clip off screen data
	kbeamx=min(kbeamx,kmxscr)
	kbeamx=max(kbeamx,0)
	kbeamy=min(kbeamy,kmxscr)	!?780
	kbeamy=max(kbeamy,0)
c * decode coordinates for reduced data transmission
	kx=kbeamx*kfactr
	ky=kbeamy*kfactr
	hiy=iand(ky/128,hex1f)+32
	lsb=iand(ky,hex03)*4+iand(kx,hex03)+96
	loy=iand(ky/4,hex1f)+96
	hix=iand(kx/128,hex1f)+32
	lox=iand(kx/4,hex1f)+64
c * send coordinates
	if(hiy.ne.kpchar(1)) then
	  kpchar(1)=hiy
	  call toutpt(hiy)
	endif
	if(lsb.ne.kpchar(2)) then	!.and. kterm.gt.1
	  kpchar(2)=lsb
	  call toutpt(lsb)
	  kpchar(3)=-1
	endif
	if(loy.ne.kpchar(3) .or. hix.ne.kpchar(4)) then
	  kpchar(3)=loy
	  call toutpt(loy)
	endif
	if(hix.ne.kpchar(4)) then
	  kpchar(4)=hix
	  call toutpt(hix)
	endif
	call toutpt(lox)
	kmovef=0
	return
	end
c%
c======routines that get/set values in common===============
c
c----------------------------------------------------------------/term
c * set (tektronix) terminal type and screen resolution
c
	subroutine term(iterm,iscal)
	include 'sgtcs_inc'
c * set term type up, but not down
	if(iterm.ge.kterm) kterm=iterm
c * default iscal = 1024
	if(iscal.gt.0) kfactr=4	!1024 screen
c * set screen resolution unless override by terminal_res
	if(iscal.gt.1024 .and. ksres.ne.1024) kfactr=1
	kmxscr=4095/kfactr
cx	call reset
	return
	end
c%
c----------------------------------------------------------------/dwindo
c * set the virtual window
c
	subroutine dwindo(xmin,xmax,ymin,ymax)
	include 'sgtcs_inc'
c
	tminvx=xmin
	tmaxvx=xmax
	tminvy=ymin
	tmaxvy=ymax
	go to 100
c%
c----------------------------------------------------------------/swindo
c * set the screen window
c
	entry swindo(minx,lx,miny,ly)
c
	kminsx=minx
	kmaxsx=minx+lx
	kminsy=miny
	kmaxsy=miny+ly
	go to 100
c%
c----------------------------------------------------------------/twindo
c * set the screen window
c
	entry twindo(minx,maxx,miny,maxy)
c
	kminsx=minx
	kmaxsx=maxx
	kminsy=miny
	kmaxsy=maxy
	!go to 100
c%
c----------------------------------------------------------------/rescal
c * set virtual to screen windown transformation
c
	entry rescal
c
 100	continue
	tmintx=tminvx
	tmaxtx=tmaxvx
	tminty=tminvy
	tmaxty=tmaxvy
c * if log mode, transform user space to log
	if(kxtype.eq.2) then	!log x
	  if(tmintx.gt.0.) tmintx=log10(tmintx)
	  if(tmaxtx.gt.0.) tmaxtx=log10(tmaxtx)
	endif
	if(kytype.eq.2) then	!log y
	  if(tminty.gt.0.) tminty=log10(tminty)
	  if(tmaxty.gt.0.) tmaxty=log10(tmaxty)
	endif
c
	dtx=tmaxtx-tmintx
	if(dtx.eq.0.) dtx=1.
	dty=tmaxty-tminty
	if(dty.eq.0.) dty=1.
c  set user to screen transformation
	tvsmx=float(kmaxsx-kminsx)/dtx
	tvsbx=float(kminsx)-tvsmx*tmintx
	tvsmy=float(kmaxsy-kminsy)/dty
	tvsby=float(kminsy)-tvsmy*tminty
	return
	end
c%
c----------------------------------------------------------------/csize
c * return character size in raster units from common
c
	subroutine csize(ihorz,ivert)
	include 'sgtcs_inc'
c
	if(kfactr.eq.0) then
	  stop 'STOP in sglib csize: initt not called'
	endif
	ihorz=(khorsz+kfactr/2)/kfactr
	ivert=(kversz+kfactr/2)/kfactr
	return
	end
c%
c----------------------------------------------------------------/aglnlg
c * communicate between agii and tcs common blocks
c
	subroutine aglnlg(ixtype,iytype)
	include 'sgtcs_inc'
c
	kxtype=ixtype
	kytype=iytype
	return
	end
c%
c----------------------------------------------------------------/setclp
c * set clip parameter in common
c
	subroutine setclp(ivalue)
	include 'sgtcs_inc'
c
	kclip=ivalue
	return
c
c----------------------------------------------------------------/agclip
c
	entry agclip
c
	kclip=0		!clipping on
	return
c
c----------------------------------------------------------------/agnclp
c
	entry agnclp
c
	kclip=1		!clipping off
	return
	end
c%
c----------------------------------------------------------------/seetrm
c * return terminal parameters from common
c
	subroutine seetrm(ibaud,iterm,icsize,maxscr)
	include 'sgtcs_inc'
c
	ibaud=0		!baud no longer stored
	iterm=kterm
	icsize=ksizef
	maxscr=kmxscr
	return
	end
c%
c----------------------------------------------------------------/agseet
c * return terminal parameters from common
c
	subroutine agseet(i4105,ifactr,icsize,maxscr)
	include 'sgtcs_inc'
c 04/15/99 CAL: moved to block data
c	data kfactr /4/	!in case initt not linked
c
	i4105=k4105
	ifactr=4/kfactr
	icsize=ksizef
	maxscr=kmxscr
	return
	end
c%
c----------------------------------------------------------------/agsmodl
c * return terminal model from common
c
	subroutine agsmodl(cmodel)
	character cmodel*(*)
	include 'sgtcs_inc'
c
	cmodel=model
	return
	end
c%
c----------------------------------------------------------------/seeloc
c * return screen position from common
c
	subroutine seeloc(ix,iy)
	include 'sgtcs_inc'
c
	ix=kbeamx	!*kfactr/4
	iy=kbeamy	!*kfactr/4
	return
	end
c%
c----------------------------------------------------------------/seetw
c * return screen window from common
c
	subroutine seetw(minx,maxx,miny,maxy)
	include 'sgtcs_inc'
c
	minx=kminsx	!*kfactr/4
	maxx=kmaxsx	!*kfactr/4
	miny=kminsy	!*kfactr/4
	maxy=kmaxsy	!*kfactr/4
	return
	end
c%
c----------------------------------------------------------------/seedw
c * return virtual window from common
c
	subroutine seedw(xmin,xmax,ymin,ymax)
	include 'sgtcs_inc'
c
	xmin=tminvx
	xmax=tmaxvx
	ymin=tminvy
	ymax=tmaxvy
	return
	end
c%
c----------------------------------------------------------------/linwdt
c * return width of character string in raster units
c
	function linwdt(numchr)
	include 'sgtcs_inc'
c
	linwdt=(khorsz*numchr+kfactr/2)/kfactr
	return
	end
c%
c----------------------------------------------------------------/linhgt
c * return height of vertical string in raster units
c
	function linhgt(numlin)
	include 'sgtcs_inc'
c
	linhgt=(kversz*numlin+kfactr/2)/kfactr
	return
	end
c%
c----------------------------------------------------------------/chrsiz
c * set character size
c
	subroutine chrsiz(k)
	include 'sgtcs_inc'
	integer ichrtb(2,4,2)
	data ichrtb /
     >	  56,88, 51,82, 34,53, 31,48,	!4014 sizes
     >	  51,88, 51,88, 31,51, 31,51 /	!4105 sizes
c
	ksizef=k
	ksizef=max(ksizef,1)
	ksizef=min(ksizef,4)
	if(k4105.ne.0) then
	  khorsz=ichrtb(1,ksizef,2)
	  kversz=ichrtb(2,ksizef,2)
c * only height is used by VT, 4510; but set spacing for good measure
	  ispc=16
	  iwid=.69*khorsz	!default is 40 for 41xx, 39 for 4510
	  ihgt=.69*kversz	!default is 60 for 41xx, 59 for 4510
	  call cmdout('MC')
	  call intout(iwid)	!width
	  call intout(ihgt)	!height
	  call intout(ispc)	!spacing
	else
	  khorsz=ichrtb(1,ksizef,1)
	  kversz=ichrtb(2,ksizef,1)
	  icode=55+ksizef
	  call toutpt(27)
	  call toutpt(icode)
	endif
	return
	end
c%
c----------------------------------------------------------------/icclip
c * return clip parameter from common
c
	function icclip(idum)
	include 'sgtcs_inc'
c
	icclip=kclip
	return
	end
c%
c%

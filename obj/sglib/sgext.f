c  sgext.for - set/examine bppcom common
c
c
c  may-79	create module [jane murphy]
c
c  mar-86	revise slimx, slimy for 4096 resolution [met]
c
c  mar-87	add agyband, hbarst, vbarst; reorder routines. [met]
c
c  jun-88	modify slim/dlim/xlen/ylen for 4096; f77 changes.
c		add agccurv, agcgrid, agctext.  [met]
c
c  jan-90	make y routines entries of x routines to speed linking. [met]
c
c  oct-90	add agsteps, agnpts; error in agbasx if > one time base.
c		allow upper or lower case ascii argument in place.
c		set minor tick form on call to xfrm/yfrm. [met]
c
c  feb-92	reorder and make many more subroutines into entries.
c		remove x/ylabnd (set flag that was always overwritten). [met]
c
c  may-00       if CRAY use "numarg()" - similar to VMS. [cal]
c-----------------------------------------------------------/dlimx/dlimy
c * set data limits for the x/y axis
c
	subroutine dlimx(xmin,xmax)
c
	include 'sgcopr_inc'
c
	call comset(ibasex(11),xmin)	!cxdmin, cxdmax
	call comset(ibasex(12),xmax)
	return
c
	entry dlimy(ymin,ymax)
	call comset(ibasey(11),ymin)	!cydmin, cydmax
	call comset(ibasey(12),ymax)
	return
c%
c-----------------------------------------------------------/slimx/slimy
c * set screen limits for the x/y axis
c
	entry slimx(ixmin,ixmax)
	call agseet(i4105,kfac,ic,is)
cx should check 0<ixmin<4096, 0<ixmax<4096
	call comset(ibasex(13),float(ixmin*kfac))	!cxsmin
	call comset(ibasex(14),float(ixmax*kfac))	!cxsmax
	return
c
	entry slimy(iymin,iymax)
	call agseet(i4105,kfac,ic,is)
	call comset(ibasey(13),float(iymin*kfac))	!cysmin
	call comset(ibasey(14),float(iymax*kfac))	!cysmax
	return
c%
c----------------------------------------------------------/dinitx/dinity
c * reinitialize axis variables
c
	entry dinitx
	nbase=ibasex(0)
	call comset(nbase+11,0.)	!data limits
	call comset(nbase+12,0.)
	call comset(nbase+5,0.)		!label width
	call comset(nbase+17,0.)	!number of ticks
	call comset(nbase+10,0.)	!number of decimal places
	call comset(nbase+18,0.)	!label expenent
	call comset(nbase+8,0.)		!number of minor tick marks
	return
c
	entry dinity
	nbase=ibasey(0)
	call comset(nbase+11,0.)	!data limits
	call comset(nbase+12,0.)
	call comset(nbase+5,0.)		!label width
	call comset(nbase+17,0.)	!number of ticks
	call comset(nbase+10,0.)	!number of decimal places
	call comset(nbase+18,0.)	!label exponent
	call comset(nbase+8,0.)		!number of minor tick marks
	return
	end
c%
c------------------------------------------------------------/xfrm/yfrm
c * set form of axis major tick mark (short ticks or grid lines)
c
	subroutine xfrm(ivalue)
	common /bppcom/ comget(80)
	integer mform(0:6)
	data mform /0,1,2,3,4,2,3/
	if(ivalue.ge.0 .and. ivalue.le.6) then
	  call comset(ibasex(7),float(ivalue))		!cxfrm
	  call comset(ibasex(9),float(mform(ivalue)))	!cxmfrm
	endif
	return
c
	entry yfrm(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.6) then
	  call comset(ibasey(7),float(ivalue))		!cyfrm
	  call comset(ibasey(9),float(mform(ivalue)))	!cymfrm
	endif
	return
c%
c-----------------------------------------------------------/xmfrm/ymfrm
c * set form of axis minor tick marks
c
	entry xmfrm(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.6) then
	  call comset(ibasex(9),float(ivalue))	!cxmfrm
	endif
	return
c
	entry ymfrm(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.6) then
	  call comset(ibasey(9),float(ivalue))	!cymfrm
	endif
	return
c%
c------------------------------------------------------------/xlab/ylab
c * set which tick marks to label (0=none, 1=all major, 2=ends only)
c
	entry xlab(ivalue)
	call comset(ibasex(3),float(ivalue))	!cxlab
	return
c
	entry ylab(ivalue)
	call comset(ibasey(3),float(ivalue))	!cylab
	return
c%
cxc----------------------------------------------------------/xlabnd/ylabnd
cxc * set flag to label end ticks only
cxc
cx	entry xlabnd(ivalue)
cx	call comset(ibasex(3),2.)	!was ibasex(28)=float(ivalue)
cx	return
cxc
cx	entry ylabnd(ivalue)
cx	call comset(ibasey(3),2.)	!was ibasey(28)=float(ivalue)
cx	return
cxc%
c------------------------------------------------------------/xlen/ylen
c * set length of major tick marks
c
	entry xlen(ivalue)
	call agseet(i4105,kfac,ic,is)
	if(ivalue.ge.0) then
	  call comset(ibasex(6),float(ivalue*kfac))	!cxlen
	endif
	return
c
	entry ylen (ivalue)
	call agseet(i4105,kfac,ic,is)
	if(ivalue.ge.0) then
	  call comset(ibasey(6),float(ivalue*kfac))	!cylen
	endif
	return
c%
c------------------------------------------------------------/xloc/yloc
c * set axis location
c
	entry xloc(ivalue)
	call agseet(i4105,kfac,ic,is)
	call comset(ibasex(2),float(ivalue*kfac))	!cxloc
	return
c
	entry yloc(ivalue)
	call agseet(i4105,kfac,ic,is)
	call comset(ibasey(2),float(ivalue*kfac))	!cyloc
	return
c%
c----------------------------------------------------------/xloctp/ylocrt
c * set axis location relative to the top/right edge
c
	entry xloctp(irast)
	max = ibasey(14)
	min = ibasey(13)
	iloc = abs(comget(max)-comget(min))
	call agseet(i4105,kfac,ic,is)
	call comset(ibasex(2),float(iloc+irast*kfac))
	return
c
	entry ylocrt(irast)
	max = ibasex(14)
	min = ibasex(13)
	iloc = abs(comget(max)-comget(min))
	call agseet(i4105,kfac,ic,is)
	call comset(ibasey(2),float(iloc+irast*kfac))
	return
c%
c-----------------------------------------------------------/xneat/yneat
c * set axis "neatness" (0=fit data exactly, 1=extend to next minor tick,
c *			2=to next major tick)
c
	entry xneat(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.2) then
	  call comset(ibasex(0),float(ivalue))	!cxneat
	endif
	return
c
	entry yneat(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.2) then
	  call comset(ibasey(0),float(ivalue))	!cyneat
	endif
	return
c%
c-----------------------------------------------------------/xtype/ytype
c * set type of data (0=linear, 1=log)
c
	entry xtype(ivalue)
	if(ivalue.ge.1 .and. ivalue.le.2) then
	  call comset(ibasex(15),float(ivalue))	!cxtype
	endif
	return
c
	entry ytype(ivalue)
	if(ivalue.ge.1 .and. ivalue.le.2) then
	  call comset(ibasey(15),float(ivalue))	!cytype
	endif
	return
	end
c%
c-----------------------------------------------------------/xtics/ytics
c * set number of tick mark intervals
c
	subroutine xtics(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.19) then
	  call comset(ibasex(5),float(ivalue))	!cxtics
	endif
	return
c
	entry ytics(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.19) then
	  call comset(ibasey(5),float(ivalue))	!cytics
	endif
	return
c%
c-----------------------------------------------------------/xmtcs/ymtcs
c * set number of minor tick mark intervals
c
	entry xmtcs(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.19) then
	  call comset(ibasex(8),float(ivalue))	!cxmtcs
	endif
	return
c
	entry ymtcs(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.19) then
	  call comset(ibasey(8),float(ivalue))	!cymtcs
	endif
	return
c%
c-----------------------------------------------------------/xwdth/ywdth
c * set max number of characters in tick mark labels
c  if ivalue < 5, sets the max characters in the tick label
c  if ivalue >=5, sets the min value of the tick remote exponent
c
	entry xwdth(ivalue)
	if(ivalue.ge.1 .and. ivalue.le.8) then
	  call comset(ibasex(16),float(ivalue))	!cxwdth
	endif
	return
c
	entry ywdth(ivalue)
	if(ivalue.ge.1 .and. ivalue.le.8) then
	  call comset(ibasey(16),float(ivalue))	!cywdth
	endif
	return
c%
c-----------------------------------------------------------/xzero/yzero
c * set flag to force axis to include origin
c
	entry xzero(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.1) then
	  call comset(ibasex(1),float(ivalue))	!cxzero
	endif
	return
c
	entry yzero(ivalue)
	if(ivalue.ge.0 .and. ivalue.le.1) then
	  call comset(ibasey(1),float(ivalue))	!cyzero
	endif
	return
	end
c%
c------------------------------------------------------------/line
c * set type of line to draw
c
	subroutine line(ivalue)
	call comset(ibasec(0),float(ivalue))
	return
c%
c------------------------------------------------------------/npts/agnpts
c * set number of points in the data array
c
	entry npts(ivalue)
	call comset(ibasec(4),float(ivalue))	!npty
	call comset(ibasec(14),float(ivalue))	!nptx
	return
c
	entry agnpts(ivalx,ivaly)
	call comset(ibasec(4),float(ivaly))	!npty
	call comset(ibasec(14),float(ivalx))	!nptx
	return
c%
c-----------------------------------------------------------/sizes
c * set the size of the data point symbol
c
	entry sizes(value)
	call comset(ibasec(7),value)	!sizes
	return
c%
c-----------------------------------------------------------/stepl
c * set increment between data points for line connections
c
	entry stepl(ivalue)
	call comset(ibasec(5),float(ivalue))	!cstepl
	return
c%
c-----------------------------------------------------------/steps/agsteps
c * set increment between data point symbols (and offset to first symbol)
c
	entry steps(ivalue)
	call comset(ibasec(2),float(ivalue))	!csteps	call comset(ibasec(6),0.)	!cfrsts
	return
c
	entry agsteps(ivals,ivali)
	call comset(ibasec(2),float(ivals))	!csteps
	call comset(ibasec(6),float(ivali))	!cfrsts
	return
c%
c-----------------------------------------------------------/symbl
c * set type of data point symbol to draw
c
	entry symbl(ivalue)
	call comset(ibasec(1),float(ivalue))	!symbl
	return
	end
c%
c----------------------------------------------------------/hbarst/vbarst
c * set line type to horizontal/vertical bar
c
	subroutine hbarst(ishade,iwbar,idbar)
	common /bppcom/ comget(80)
	nbase=ibasec(0)
	call comset(nbase,-3.)
	if(ishade.ge.0 .and. ishade.le.15) then
	  call comset(nbase+1,float(ishade))
	endif
	call comset(nbase+8,float(iwbar))
	call comset(nbase+7,float(idbar))
c * alter tick mark specification to prevent lines thru bars
	nbase=ibasey(0)
	ticfrm=comget(nbase+7)
	if(ticfrm.eq.5.) ticfrm=2.
	if(ticfrm.eq.6.) ticfrm=1.
	call comset(nbase+7,ticfrm)
	return
c
	entry vbarst(ishade,iwbar,idbar)
	nbase=ibasec(0)
	call comset(nbase,-2.)
	if(ishade.ge.0 .and. ishade.le.15) then
	  call comset(nbase+1,float(ishade))
	endif
	call comset(nbase+8,float(iwbar))
	call comset(nbase+7,float(idbar))
c * alter tick mark specification to prevent lines thru bars
	nbase=ibasex(0)
	ticfrm=comget(nbase+7)
	if(ticfrm.eq.5.) ticfrm=2.
	if(ticfrm.eq.6.) ticfrm=1.
	call comset(nbase+7,ticfrm)
	return
	end
c%
c-----------------------------------------------------------/place
c * set the screen window from ascii location
c
	subroutine place(code)
	character code*3
	character*3 litary(13),ucode
	real screen(4,13)
c * table of ascii names for locations
	data litary/'STD','UPH','LOH','UL4','UR4','LL4','LR4','UL6',
     >  'UC6','UR6','LL6','LC6','LR6'/
c * table of co-ordinates for locations (xmin,xmax,ymin,ymax)
	data screen/150.,900.,125.,700.,
     1		 150.,850.,525.,700.,
     2		 150.,850.,150.,325.,
     3		 150.,450.,525.,700.,
     4		 650.,950.,525.,700.,
     5		 150.,450.,150.,325.,
     6		 650.,950.,150.,325.,
     7		 150.,325.,525.,700.,
     8		 475.,650.,525.,700.,
     9		 800.,975.,525.,700.,
     1		 150.,325.,150.,325.,
     2		 475.,650.,150.,325.,
     3		 800.,975.,150.,325./
c * convert to upper case
	call iag_upcase(ucode,code)
c * search for a literal match
	do i=1,13
	  if(ucode.eq.litary(i)) go to 200
	enddo
c * see if it is an integer index
	i=ichar(code(1:1))
c * default to full screen
	if(i.le.0 .or. i.gt.13) i=1
c
200	call agseet(i4105,kfac,ic,is)
c * set the screen co-ordinates into common
	call comset(ibasex(13),screen(1,i)*kfac)
	call comset(ibasex(14),screen(2,i)*kfac)
	call comset(ibasey(13),screen(3,i)*kfac)
	call comset(ibasey(14),screen(4,i)*kfac)
	return
	end
c%
c---------------------------------------------------------/agccurv
c * set the curve color for 4105 output
c
	subroutine agccurv(ivalue)
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) call comset(ibasec(10),float(ivalue))
	return
c
c---------------------------------------------------------/agcgrid
c * set the grid color for 4105 output
c
	entry agcgrid(ivalue)
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) call comset(ibasec(11),float(ivalue))
	return
c
c---------------------------------------------------------/agctext
c * set the text color for 4105 output
c
	entry agctext(ivalue)
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) call comset(ibasec(12),float(ivalue))
	return
	end
c%
c----------------------------------------------------------/agbasx
c * set the number of bases in the data array
c
	subroutine agbasx(ivalue)
	if(ivalue.gt.1) then
	  call agsalf
c	  type *, 'agbasx: multiple bases no longer allowed'
	  call at_msg('agbasx: multiple bases no longer allowed')
	  call agsgrf
	endif
	call comset(ibasec(9),float(ivalue))	!agbasx
	return
c%
c-----------------------------------------------------------/infin
c * specify largest acceptable floating point number
c
	entry infin(value)
	if(value.gt.0.) call comset(ibasec(3),value)
	return
c%
c---------------------------------------------------------/agyband
c * set max variation for fast cplot "straight" line as percent of window
c
	entry agyband(value)
	if(value.ge.0 .and. value.le.20.) then
	  call comset(ibasey(21),value)	!cyband
	endif
	return
	end
c%
cxc-----------------------------------------------------------/sizel
cxc * set width of line
cxc
cx	entry sizel(value)
cx	call comset(ibasec(8),value)	!sizel
cx	return
cx	end
cxc
c%
c  sgtxt.for - text routines
c
c  jan-82	written by jane murphy
c
c  mar-82	change agths/agtvs/agtlng so string can be char or integer.
c		(jm)
c
c  apr-84	move remote exponent to end of x-axis label. (jam)
c
c  jun-85	add agvcur. (met)
c
c  apr-86	modify agtlng to get length from char string descriptor
c		if a character variable. (met)
c
c  nov-86	add local string astrin and center argument to agtlng to
c		drop leading blanks if centering. (met)
c
c  oct-87	add agtmrs/agtrs; set text color in 4105 mode. (met)
c		write vertical strings at 90 rotation in 4105 mode.
c		check for log(0) in agtcvs/agtmvs.
c
c  jun-88	fix 4096 resolution (my*kfac) for agtchs/agtmhs/agtcgt. (met)
c
c  may-90	fix bug - agtcgt not centering string.
c		move 4105 vertical strings to the left.
c
c  jan-91	scale offsets based on character size.
c		revise logic in agtchs,agtcvs to match label.
c
c  may-91	dont drop trailing blanks if string is $ terminated. [met]
c
c--------------------------------------------------------/agtchs
c * center string on x axis
c
	subroutine agtchsl(string,length) !   
	character string*(*)
	common /bppcom/ comget(80)
	character astrin*80
	logical noout(0:6)	!no part of tick outside axis
	data noout /.true.,.true.,.false.,.true.,.false.,.false.,.true./
c
	nchars=length		!   
c  copy to astrin and get number of characters, centered
	call agtlng(string,.true.,astrin,nchars)
c  calculate center position
	nbasex=ibasex(0)
	near=comget(nbasex+13)
	nfar=comget(nbasex+14)
	maxis=(nfar-near)/2+near
c  half string length
	call csize(ihorz,ivert)
	ihalf=(nchars*ihorz)/2
	mx=maxis-ihalf
c if only end ticks are labeled, push unit label between tick labels
	lend=comget(nbasex+3)		!label ends only
	iloc=comget(nbasex+2)		!top/bottom
	lmaj=comget(nbasex+6)		!tick length
	lform=comget(nbasex+7)		!tick form
	if(noout(lform)) lmaj=7
	iexp=comget(nbasex+18)		!remote exponent
	if(iloc.le.0) then		!x bottom
	  loff=-ivert*(1.1+1.4)		!was *2.6
	  if(lend.eq.2) then		!if tick labels only at ends
	    loff=-ivert*(1.1+.5)	!push up between labels
	  elseif(iexp.ne.0) then	!if remote exponent
	    loff=loff-ivert/3		!move down to level of exponent
	    mx=mx-(4*ihorz/2)		!shorten space by exp length
	  endif
	  lmaj=-lmaj
	else				!x top
	  loff=ivert*(.45+1.5)		!was *2.0
	endif
	ilbpos=comget(ibasey(13))
	my=ilbpos+iloc+lmaj+loff
c  keep label on screen
	call agseet(i4105,kfac,ic,maxscr)
	my=min(max(my,0),maxscr*3/4)
c  move and output horizontal string
	call movabs(mx,my)
	go to 1000
c
c---------------------------------------------------------/agtmhs
c * move and output horizontal string
c
	entry agtmhsl(ix,iy,string,length) !   
	nchars=length  !   
	go to 200
c
	entry agxmhs(ix,iy,string)
	nchars=0
 200	continue
	call agseet(i4105,kfac,ic,is)
	call movabs(ix*kfac,iy*kfac)
	go to 400
c
c--------------------------------------------------------/agtcgt
c * output centered graph title
c
	entry agtcgtl(iy,string,length)  !   
	nchars=length  !   
	call agtlng(string,.true.,astrin,nchars)
c  calculate center position
	nbasex=ibasex(0)
	near=comget(nbasex+13)
	nfar=comget(nbasex+14)
	maxis=(nfar-near)/2+near
c  half string length
	ihalf=(nchars*linwdt(1))/2
	mx=maxis-ihalf
	call agseet(i4105,kfac,ic,is)
	my=iy*kfac
	if(my.eq.0) my=ifix(comget(ibasey(14)))+15*kfac
c  move and output horizontal string
	call movabs(mx,my)
	go to 1000
c
c---------------------------------------------------------/agths
c * output horizontal string
c
	entry agthsl(string,length)  !   
	nchars=length  !   
	go to 400
c
	entry agxhs(string)
	nchars=0
 400	continue
c  copy to astrin and get length in nchars, not centered
	call agtlng(string,.false.,astrin,nchars)
c
 1000	continue
c  set text color
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) then
	  itext=comget(ibasec(12))
	  call agtxclr(itext)	!set text color
	  call cmdout('LT')
	  call intout(nchars)
	else
	  call alfmod		! enter ascii mode
	endif
c  output string
	do ii=1,nchars
	  ich=ichar(astrin(ii:ii))
	  call toutpt(ich)
	enddo
c  update beam position (for next character)
	call movrel(linwdt(nchars),0)
	return
	end
c%
c--------------------------------------------------------/agtcvs
c * center string on y-axis
c
	subroutine agtcvsl(string,length)  !   
	character string*(*)
	common /bppcom/ comget(80)
	logical noout(0:6)	!no part of tick outside axis
	data noout /.true.,.true.,.false.,.true.,.false.,.false.,.true./
	character astrin*80
c
	nchars=length  !   
	call agtlng(string,.true.,astrin,nchars)
c  calculate center position
	nbasey=ibasey(0)
	near=comget(nbasey+13)
	nfar=comget(nbasey+14)
	maxis=(nfar-near)/2+near
	call agseet(i4105,if,ic,is)
c  half character length
	call csize(ihorz,ivert)
	if(i4105.ne.0) then
	  ihalf=(nchars*ihorz)/2
	else
	  ihalf=(nchars*ivert)/2
	endif
	my=(maxis+ihalf)+4
	iwidth=comget(nbasey+17)	!tick label width
c  compute x position
	iloc=comget(nbasey+2)		!top/bottom
	lend=comget(nbasey+3)		!label ends only
	lmaj=comget(nbasey+6)		!tick length
	lform=comget(nbasey+7)		!tick form
	if(noout(lform)) lmaj=7
	if(iloc.le.0) then		!bottom
	  iwidth=-iwidth
	  ioff=-ihorz*2.6
	  if(lend.eq.2) ioff=-ihorz	!was (iwidth=-1)
	  lmaj=-lmaj
	else				!top
	  ioff=ihorz*2.0
	endif
	ilbpos=comget(ibasex(13))
	mx=ilbpos+iloc+lmaj+iwidth*ihorz+ioff
	mx=max(mx,0)			!keep label on screen
	if(i4105.ne.0) mx=max(mx,ihorz/2)
	call agseet(i4105,if,ic,maxscr)
	mx=min(mx,maxscr-ihorz)
c  move and output vertical string
	call movabs(mx,my)
	go to 1000
c
c---------------------------------------------------------/agtmvs
c * move and output vertical string
c
	entry agtmvsl(ix,iy,string,length)  !   
	nchars=length		!   
	go to 200
c
	entry agxmvs(ix,iy,string)
	nchars=0
 200	continue
	call agseet(i4105,kfac,ic,maxscr)
	call movabs(ix*kfac,iy*kfac)
	go to 400
c
c---------------------------------------------------------/agtvs
c * output vertical string
c
	entry agtvsl(string,length)  !   
	nchars=length		!   
	go to 400
c
	entry agxvs(string)
	nchars=0
 400	continue
	call agtlng(string,.false.,astrin,nchars)
c
 1000	continue
	ivert=linhgt(1)		! get character height
c see if can draw label vertically
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) then
c at top of string - move down length of string
	  call movrel(0,-linwdt(nchars))
	  call movrel(ivert,0)
	  itext=comget(ibasec(12))
	  call agtxclr(itext)	!set text color
	  call agtxrot(90)	!set rotation
	  call cmdout('LT')	!declare text string
	  call intout(nchars)
	endif
c loop over characters
	  do ii=1,nchars
	    if(i4105.eq.0) then
	      call movrel(0,-ivert)	! move down 1 character height
	      call toutpt(31)	  	! ascii mode
	    endif
	    ich = ichar(astrin(ii:ii))
	    call toutpt(ich)
	  enddo
c  update beam position (for next character)
	  if(i4105.eq.0) then
	    call movrel(0,-ivert)
cx	    call alfmod
	  endif
	  if(i4105.ne.0) then
	    call agtxrot(0)		!clear rotation
	  endif
cx	endif
	return
	end
c%
c---------------------------------------------------------/agtmrs
c * move and output rotated string
c
	subroutine agtmrsl(ix,iy,irot,string,length)  !   
	character string*(*)
	character astrin*80
c
	nchars=length  !   
	call agseet(i4105,kfac,ic,maxscr)
	call movabs(ix*kfac,iy*kfac)
	go to 400
c
c---------------------------------------------------------/agtrs
c * output rotated string
c
	entry agtrsl(irot,string,length)  !   
	nchars=length  !   
c
 400	continue
	call agtlng(string,.false.,astrin,nchars)
c
c see if 4105
	call agseet(i4105,if,ic,is)
	if(i4105.ne.0) then
	  itext=comget(ibasec(12))
	  call agtxclr(itext)	!set text color
	  call agtxrot(irot)	!set rotation
	  call cmdout('LT')	!declare text string
	  call intout(nchars)
	endif
c
c loop over characters
	do ii=1,nchars
	  ich = ichar(astrin(ii:ii))
	  call toutpt(ich)
	enddo
	if(i4105.ne.0) then
	  call agtxrot(0)	!clear rotation
	endif
	return
	end
c%
c---------------------------------------------------------/agtlng
c * determine length of string
c
	subroutine agtlng(string,center,astrin,nchars)
	character string*(*)
	logical center		!true if string is to be centered
	character astrin*80
	character chr*1
c
	nloop=nchars		!user specified length
	if(nloop.eq.0) nloop=len(string)
	if(nloop.le.0) return	!error if still nothing
cx	  call lib$signal
	nloop=min(nloop,80)
c  count the characters
	nc=0
	do i=1,nloop
	  chr=string(i:i)
	  if(ichar(chr).eq.0) go to 2000
	  if(chr.eq.'$' .and. nchars.eq.0) then
c				!if ends in $, drop $,
	    nchars=nc-1		! but not trailing blanks
	    go to 2000
	  endif
	  if(center .and. nc.eq.0 .and. chr.eq.' ') then
c				!drop leading spaces if center
	  else
	    nc=nc+1
	    astrin(nc:nc) = chr
	  endif
	enddo
 2000	continue
c		!drop trailing spaces if center or no user length
	if(center .or. nchars.eq.0) then
	  do while (nc.gt.0 .and. astrin(nc:nc).eq.' ')
	    nc=nc-1
	  enddo
	endif
	nchars = nc
	return
	end
c%
c---------------------------------------------------------/agvcur
c * get cursor location and ascii char
c
	subroutine agvcur(uch,x,y)
	character*1 uch
c
	call scursr(ich,ixa,iya)
	call revcot(ixa,iya,x,y)
	uch = char(ich)
	return
	end
c%
c  sgerr.for - error bar routines
c
c
c  feb-82	create new routine names. (jam)
c
c  oct-85	clip negative log error bars. (jane murphy)
c
c  jul-88	set line color for y bars.  (met)
c		remove obsolete entry names.
c
c  nov-88	get common variables only once for ageya/agexa. (met)
c		set line color for x bars.
c		handle unclipped log data - get type before checking iclip.
c
c  aug-90	make common noshr. (met)
c
c  nov-90	get data values same as in cplot (allow all data formats)
c		  in ageya/agexa. (met)
c
c  feb-91	use char width as width of tick instead of constant. [met]
c
c  feb-92	add agexsr/ageysr, single error bars of unequal lengths. [met]
c
c----------------------------------------------------------------/ageya
c * draw an array of y error bars
c
	subroutine ageya(x,y,e)
c
	dimension x(*), y(*), e(*)
cdec$   psect  /sgerr/ noshr
	common /sgerr/ kxtype,kytype,xmin,xmax,ymin,ymax,iclip
c
	isym=comget(ibasec(1))
	istepl=comget(ibasec(5))
	icurv=comget(ibasec(10))
	call linclr(icurv)	!if(icurv.ne.1)
	call agkeys(x,y,ifastx,ifasty,limit)
	if(ifastx.eq.-1) then
	  xfirs=x(3)
	  xdelt=x(4)
	elseif(ifastx.eq.-2) then
	  xfirs=x(2)
	  xdelt=x(3)
	endif
	if(limit.eq.0) return
c
	iclip=icclip(0)		!get params once for entire array
	if(iclip.eq.0) then
	  kxtype=comget(ibasex(15))	!if clipping on, need to check x
	  xmin=comget(ibasex(26))
	  xmax=comget(ibasex(27))
	  if(kxtype.eq.2) then
	    xmin=10**xmin
	    xmax=10**xmax
	  endif
	endif
	kytype=comget(ibasey(15))	!always need to check y
	ymin=comget(ibasey(26))
	ymax=comget(ibasey(27))
	if(kytype.eq.2) then
	  ymin=10**ymin
	  ymax=10**ymax
	endif
	do i=1,limit,istepl
	  if(ifastx.gt.0) then
	    xpoint=x(ifastx)
	    ifastx=ifastx+istepl
	  else	!if(ifastx.lt.0) then
	    xpoint=xfirs+xdelt*float(i-1)
	  endif
	  if(ifasty.gt.0) then
	    ypoint=y(ifasty)
	    ebar=e(ifasty)
	    ifasty=ifasty+istepl
	  else	!if(ifasty.lt.0) then
	    ypoint=y(3)+y(4)*float(i-1)
	    ebar=e(3)+e(4)*float(i-1)
	  endif
	  call agey1(xpoint,ypoint,ebar)
	  if(isym.ne.0) call bsyms(xpoint,ypoint,isym)
	enddo
c
	return
	end
c%
c----------------------------------------------------------------/agexa
c * draw an array of x error bars
c
	subroutine agexa(x,y,e)
c
	dimension x(*), y(*), e(*)
cdec$   psect  /sgerr/ noshr
	common /sgerr/ kxtype,kytype,xmin,xmax,ymin,ymax,iclip
c
	isym=comget(ibasec(1))
	istepl=comget(ibasec(5))
	icurv=comget(ibasec(10))
	call linclr(icurv)	!if(icurv.ne.1)
	call agkeys(x,y,ifastx,ifasty,limit)
	if(ifastx.eq.-1) then
	  xfirs=x(3)
	  xdelt=x(4)
	elseif(ifastx.eq.-2) then
	  xfirs=x(2)
	  xdelt=x(3)
	endif
	if(limit.eq.0) return
c
	iclip=icclip(0)		!get params once for entire array
	kxtype=comget(ibasex(15))	!always need to check x
	xmin=comget(ibasex(26))
	xmax=comget(ibasex(27))
	if(kxtype.eq.2) then
	  xmin=10**xmin
	  xmax=10**xmax
	endif
	if(iclip.eq.0) then
	  kytype=comget(ibasey(15))	!if clipping on, need to check y
	  ymin=comget(ibasey(26))
	  ymax=comget(ibasey(27))
	  if(kytype.eq.2) then
	    ymin=10**ymin
	    ymax=10**ymax
	  endif
	endif
c
	do i=1,limit,istepl
	  if(ifastx.gt.0) then
	    xpoint=x(ifastx)
	    ifastx=ifastx+istepl
	  else	!if(ifastx.lt.0) then
	    xpoint=xfirs+xdelt*float(i-1)
	  endif
	  if(ifasty.gt.0) then
	    ypoint=y(ifasty)
	    ebar=e(ifasty)
	    ifasty=ifasty+istepl
	  else	!if(ifasty.lt.0) then
	    ypoint=y(3)+y(4)*float(i-1)
	    ebar=e(3)+e(4)*float(i-1)
	  endif
	  call agex1(xpoint,ypoint,ebar)
	  if(isym.ne.0) call bsyms(xpoint,ypoint,isym)
	enddo
c
	return
	end
c%
c----------------------------------------------------------------/ageys
c * draw a single y error bar
c
	subroutine ageys(x,y,e)
c
cdec$   psect  /sgerr/ noshr
	common /sgerr/ kxtype,kytype,xmin,xmax,ymin,ymax,iclip
c
	em=e
	ep=e
	go to 100
c----------------------------------------------------------------/ageysr
c * draw single y error bar with unequal lengths
c
	entry ageysr(x,y,elo,ehi)
c
	em=elo
	ep=ehi
c
 100	continue
	iclip=icclip(0)
	if(iclip.eq.0) then
	  kxtype=comget(ibasex(15))
	  xmin=comget(ibasex(26))
	  xmax=comget(ibasex(27))
	  if(kxtype.eq.2) then
	    xmin=10**xmin
	    xmax=10**xmax
	  endif
	endif
	kytype=comget(ibasey(15))
	ymin=comget(ibasey(26))
	ymax=comget(ibasey(27))
	if(kytype.eq.2) then
	  ymin=10**ymin
	  ymax=10**ymax
	endif
	ym=y-em
	yp=y+ep
	go to 200
c
	entry agey1(x,y,e)
c
	ym=y-e
	yp=y+e
c
 200	continue
c
	if(iclip.eq.0) then
	  if(x.lt.xmin .or. x.gt.xmax) return	!x out of range
	  if(ym.le.ymin .and. yp.le.ymin) return	!off bottom
	  if(ym.ge.ymax .and. yp.ge.ymax) return	!off top
	endif
	ltic=linwdt(1)
c
	if(yp.le.ymax .or. iclip.ne.0) then
	  call movea(x,yp)		!cross tick at top
	  call movrel(-ltic,0)
	  call drwrel(ltic*2,0)
	endif
	call movea(x,yp)		!draw vertical line
	call drawa(x,ym)		!movea/drawa do clipping
c
	if(ym.ge.ymin .or. iclip.ne.0) then
	  call movrel(-ltic,0)		!cross tick at bottom
	  call drwrel(ltic*2,0)
	  call movea(x,y)
	endif
c
	return
	end
c%
c----------------------------------------------------------------/agexs
c * draw a single x error bar
c
	subroutine agexs(x,y,e)
c
cdec$   psect  /sgerr/ noshr
	common /sgerr/ kxtype,kytype,xmin,xmax,ymin,ymax,iclip
c
	em=e
	ep=e
	go to 100
c----------------------------------------------------------------/agexs
c * draw single x error bar with unequal lengths
c
	entry agexsr(x,y,elo,ehi)
c
	em=elo
	ep=ehi
c
 100	continue
	iclip=icclip(0)
	kxtype=comget(ibasex(15))
	xmin=comget(ibasex(26))
	xmax=comget(ibasex(27))
	if(kxtype.eq.2) then
	  xmin=10**xmin
	  xmax=10**xmax
	endif
	if(iclip.eq.0) then
	  kytype=comget(ibasey(15))
	  ymin=comget(ibasey(26))
	  ymax=comget(ibasey(27))
	  if(kytype.eq.2) then
	    ymin=10**ymin
	    ymax=10**ymax
	  endif
	endif
	xm=x-em
	xp=x+ep
	go to 200
c
	entry agex1(x,y,e)
c
	xm=x-e
	xp=x+e
c
 200	continue
	if(iclip.ne.0) then
	  if(y.lt.ymin .or. y.gt.ymax) return	!y out of range
	  if(xm.le.xmin .and. xp.le.xmin) return	!off left
	  if(xm.ge.xmax .and. xp.ge.xmax) return	!off right
	endif
c
	ltic=linwdt(1)
	if(xp.le.xmax .or. iclip.ne.0) then
	  call movea(xp,y)		!cross tick at right
	  call movrel(0,-ltic)
	  call drwrel(0,ltic*2)
	endif
	call movea(xp,y)		!draw horizontal line
	call drawa(xm,y)		!movea/drawa do clipping
c
	if(xm.ge.xmin .or. iclip.ne.0) then
	  call movrel(0,-ltic)		!cross tick at bottom
	  call drwrel(0,ltic*2)
	  call movea(x,y)
	endif
c
	return
	end

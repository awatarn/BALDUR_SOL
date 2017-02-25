c file: sggrf.for	sglib source (vax)
c		advanced graphing routines except those in sgcom.for
c
c
c  may-79	scientific advanced graphing ii release 3.3
c		without calendar plotting [jane murphy]
c
c  sep-80	rewrite graph routines to support irrational scaling
c
c  mar-82	change aglint, loptim and setwin to correct scaling for
c		second scale. reduce number of digits in labels.
c
c  aug-83	correct ivclip call in bsyms from das copy [jam]
c
c  jun-84	modify remlab to put x axis remote exponent at right end
c		of axis label [jam]
c
c  feb-86	remove eform (e format) option from numset, replace fonly
c		(f format), iform (i format) with fortran write.
c		remove xtent stuff from bppcom common.
c		change x/ydmin/max adjustment if y data range is zero to
c		a percent of y if y is large. [met]
c
c  mar-86	make label placement a function of linhgt for 4096 screen
c		check binitt called in check [met]
c
c  nov-86	speed up teksym by replacing most symbols with line
c		drawings; make x smaller and recenter triangle/del.
c		replace cplot with Tom Gibneys fast plot version.
c		switch sgaxi common to arrays; replace loptim args with
c		 x/y index.
c		replace ibasec calls with integer constant ibase;
c		replace comget calls with comget common. [met]
c
c  mar-87	restore bar/filbox routines to draw bar charts.
c		shrink width/grid by using common x/y code. [met]
c
c  oct-87	fix aglogt for irrational (raw data) endpoints.
c		shrink label to use common x/y code.
c		call linclr from grid/cplot/teksym using new variables
c		 in ibasec common (ccurv,cgrid,ctext,cfill).
c		add marker line (line style -5) to cplot.
c		use fast cplot algorithm for point lines; use same code
c		 for non-std and long std line type.
c		add agdspl - draw grid and axes but no lines.
c		make all tests for infinity use abs(value).
c		replace roundd/roundu from new plot10 code. [met]
c
c  jun-88	change symout to call agths instead of toutpt to get color.
c		adjust for 4096 resolution in bar.
c		call agtxclr in grid to get right color if no initt (only
c		 binitt (may not have fixed problem). [met]
c
c  sep-89	dont set neg values to finfin in cplot; let wincot clip
c		 them. [met]
c
c  jan-90	convert to lower case.
c		dont call term in binitt to set resolution; done by initt.
c		move lwidth into width; fform/iform into numset.
c		close X symbol so filled symbol looks like bowtie.
c		change aglogt to use dmax when user dmax>0, dmin<=0.
c		ignore zeros when finding min/max of log data. [met]
c
c  aug-90	fix bug where no tick labels if data constant < -10.
c		fix bug where top tick of log axis usually missing.
c		move tset into grid; simplify label.
c		do "fast" data access in cplot for short form data. [met]
c
c  oct-90	allow non-standard Y with short form X in mnmx/rgchek/cplot.
c		create solid symbols in teksym using panel fill in 4105 mode
c		 and repeatedly smaller symbols in 4014 mode.
c		add symbol types 12-15; use nint in teksym.
c		add Doug McCunes size adjustment factors for each symbol.
c		allow arbitrary start point for symbols.
c		add agkeys; "fast" data access in cplot for all data.
c		more cleanup in grid and label.
c		error in cplot if more than 1 X base. [met]
c
c  apr-91	move notate into juster; remove notate.
c		convert numset and juster to character variables.
c		rewrite roundd/roundu to fix problem with gran=.1.
c		remove umnmx. [met]
c
c  nov-91	check for dmax/dspan>1.e5 to avoid loop in roundd/u.
c
c  feb-92	add agmnmx with new args; restore mnmx to original.
c		misc fixes from MAC version.
c		output message and return if invalid data type in cplot.
c		remove keyset,datget,uline,upoint. [met]
c
c  feb-93	change texp calculation in aglint to avoid 1.**(-2)=
c		1.0000001e-2, which makes roundd inaccurate.  [met]
c
c  oct-93	fix rmin/rmax in aglint after round; texp calc yields
c		1.**(-2)=1.0000001e-2 and 1./10.**2=9.9999998E-03.
c		convert ibyte to word in juster.  [met]
c
c  aug-97	use x/ywdth value as mininum value of tick remote exponent
c		as well as minimum tick label width; remove users.  [met]
c
c  nov-2002     dmc: added comments to aglint and aglogt; took log10(dmin)
c               and log10(dmax) out of aglogt call
c%
c----------------------------------------------------------------/binitt
c * set initial values in ag common
c
	subroutine binitt
c
	include 'sgcopr_inc'
c
	common /bppcom/ cline,csymbl,csteps,
     > cinfin,cnpts,cstepl,cfrsts,csizes,csizel,cnbasx,
     > ccurv,cgrid,ctext,cfill,cnptx,cunused,
     > cxneat,cxzero,cxloc,cxlab,cxden,cxtics,
     > cxlen,cxfrm,cxmtcs,cxmfrm,cxdec,
     > cxdmin,cxdmax,cxsmin,cxsmax,cxtype,
     > cxtrex,cxwdth,cxepon,cxstep,cxstag,cxetyp,
     > cxstrt,cxend,cxmstt,cxmend,cxamin,cxamax,cxlabn,
     > cyneat,cyzero,cyloc,cylab,cyden,cytics,
     > cylen,cyfrm,cymtcs,cymfrm,cydec,
     > cydmin,cydmax,cysmin,cysmax,cytype,
     > cytrex,cywdth,cyepon,cystep,cystag,cyband,
     > cystrt,cyend,cymstt,cymend,cyamin,cyamax,cylabn,
     > unused(6)
c *
c * ibasec+0 cline	type of line used in plots
c * ibasec+1 csymbl	symbol used for point plots
c * ibasec+2 csteps	increment between symbols
c * ibasec+3 cinfin	infinity
c * ibasec+4 cnpts	# of y points for non standard arrays
c * ibasec+5 cstepl	increment between points in lines
c * ibasec+6 cfrsts	first point to put a symbol
c * ibasec+7 csizes	symbol size
c * ibasec+8 csizel	line size (bars only)
c * ibasec+9 cnbasx	# of multiple x bases
c * ibasec+10 ccurv	color of all curves
c * ibasec+11 cgrid	color of axes/ticks
c * ibasec+12 ctext	color of text
c * ibasec+13 cfill	color of filled panels
c * ibasec+14 cnptx	# of x points for non-standard arrays
c * ---ibase = common vector pointer
c * nbase+0    c(x/y)neat  neat tick flag
c * nbase+1    c(x/y)zero  zero suppresion flag
c * nbase+2    c(x/y)loc   location of tick marks
c * nbase+3    c(x/y)lab   type of labels
c * nbase+4    c(x/y)den   density of ticks - not used
c * nbase+5    c(x/y)tics  actual number of ticks
c * nbase+6    c(x/y)len   length of tick marks
c * nbase+7    c(x/y)frm   form of tick marks
c * nbase+8    c(x/y)mtcs  number of minor ticks
c * nbase+9    c(x/y)mfrm  form of minor ticks
c * nbase+10   c(x/y)dec   number of decimal places
c * nbase+11   c(x/y)dmin  data minimum
c * nbase+12   c(x/y)dmax  data maximum
c * nbase+13   c(x/y)smin  screen minimum
c * nbase+14   c(x/y)smax  screen maximum
c * nbase+15   c(x/y)type  type of data (lin/log)
c * nbase+16   c(x/y)lsig  least significant digit - not used
c * nbase+16   c(x/y)twid  user set tick label width - not used
c * nbase+16   c(x/y)trex  minimum tick remote exponent
c * nbase+17   c(x/y)wdth  tick label width
c * nbase+18   c(x/y)epon  exponent
c * nbase+19   c(x/y)step  number of label steps
c * nbase+20   c(x/y)stag  staggered tick labels - not used
c * nbase+21   c(x/y)etyp  exponent type format - not used
c * nbase+21   c(  y)band  allowable jitter in y data for fast cplot
c * nbase+22   c(x/y)strt  start location for tick marks - not used
c * nbase+23   c(x/y)end   ending location for tick marks - not used
c * nbase+24   c(x/y)mstt  start location for minor tick marks - not used
c * nbase+25   c(x/y)mend  ending location for minor tick marks - not used
c * nbase+26   c(x/y)amin  calculated data minimum
c * nbase+27   c(x/y)amax  calculated data maximum
c * nbase+28   c(x/y)labn  label end major ticks only
c * ---nbase = x or y parm vector pointer
	call agseet(i4105,kfac,ic,is)	!get screen factor
c
	cline=0.
	csymbl=0.
	csteps=1.
	cinfin=1.e36		!oct92, oct87=1.e32, orig=1.e38
	cnpts=0.
	cnptx=0.
	cnbasx=0.
	cstepl=1.
	csizes=1.
	csizel=1.
	cfrsts=0.
	ccurv=1.
	cgrid=1.
	ctext=1.
	cfill=-1.
c
	iv=linhgt(1)
	cxneat=2.
	cxzero=1.
	cxloc=0.
	cxlab=1.
	cxlabn=0.
	cxtics=0.
	cxlen=iv
	cxfrm=5.
	cxmtcs=0.
	cxmfrm=2.
	cxdec=0.
	cxdmin=0.
	cxdmax=0.
	cxsmin=150*kfac
	cxsmax=900*kfac
	cxtype=1.
	cxtrex=0.
	cxwdth=0.
	cxepon=0.
	cxstep=1.
	cxstrt=0.
	cxend=0.
	cxmstt=0.
	cxmend=0.
	cxamin=0.
	cxamax=0.
c
	cyneat=2.
	cyzero=1.
	cyloc=0.
	cylab=1.
	cylabn=0.
	cytics=0.
	cylen=iv
	cyfrm=5.
	cymtcs=0.
	cymfrm=2.
	cydec=0.
	cydmin=0.
	cydmax=0.
	cysmin=125*kfac
	cysmax=700*kfac
	cytype=1.
	cytrex=0.
	cyband=0.
	cyepon=0.
	cystep=1.
	cywdth=0.
	cystrt=0.
	cyend=0.
	cymstt=0.
	cymend=0.
	cyamin=0.
	cyamax=0.
	call setclp(0)		!set clipping on
	return
	end
c%
c----------------------------------------------------------------/check
c * set flags for grid and label routines
c
	subroutine check(x,y)
	dimension x(*), y(*)
cdec$	psect  /sgaxi/ noshr
	common /sgaxi/ k1(2),k2(2),nmod(2),gran(2),ibeg(2),dtic(2)
c
	nbasex=ibasex(0)
	nbasey=ibasey(0)
c
c set data limits
	call rgchek(nbasex,x)
	call rgchek(nbasey,y)
c
c determine parameters for ticks and tick labels
	call loptim(nbasex,1)
	call loptim(nbasey,2)
c
c set the virtual to screen window
	call setwin
c
	return
	end
c%
c----------------------------------------------------------------/dsplay
c * draw axes, ticks, and tick labels; call cplot to draw data curve
c
	subroutine dsplay(x,y)
	dimension x(*), y(*)
	common /bppcom/ comget(80)
c * draw axes and ticks
	call grid
c * label the axes if label type > 0
	call label(ibasex(0))
	call label(ibasey(0))
c * draw curve
	call cplot(x,y)
c
	return
	end
c----------------------------------------------------------------/agdspl
c * draw axes, ticks, and tick labels
c
	subroutine agdspl(x,y)
	dimension x(*), y(*)
	common /bppcom/ comget(80)
c
c * draw axes
	call grid
c * label the axes
	call label(ibasex(0))
	call label(ibasey(0))
c
	return
	end
c%
c----------------------------------------------------------------/rgchek
c * determine the minimum and maximum for the data array
c
	subroutine rgchek(nbase,array)
	real array(*)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
c * is this an incremental x array for non-standard short x
	nxbase=comget(ibasec+9)
	finfin=comget(ibasec+3)
	amin=comget(nbase+11)
	amax=comget(nbase+12)
c * linear or logarithmic
	itype=comget(nbase+15)
cc	select code itype(1,2)
	if(itype.eq.1) then
c * linear check
	  if(amin.ne.amax) go to 400
	  dmin=finfin
	  if(comget(nbase+1).ne.1.) dmin=0.	!zero include
	  dmax=-dmin
c * find the min + max of the array
	  if(nxbase.gt.0 .and. nbase.eq.ibasex(0)) then
	    call agbmnx(array,dmin,dmax,nxbase)
	  else
	    if(nbase.eq.ibasex(0)) then
	    npts=comget(ibasec+14)
	    else
	    npts=comget(ibasec+4)
	    endif
	    call agmnmx(array,dmin,dmax,npts,itype)
	  endif
c * if min = max, expand data range to prevent one dimensional graph
	  if(dmin.eq.dmax) then
	  if(dmin.gt.10.) then	!(abs(dmin).gt.10.
	    dmin=dmin*0.95	!*(1.-sign(.05,dmin))
	    dmax=dmax*1.05
	  elseif(dmin.lt.-10.) then
	    dmin=dmin*1.05
	    dmax=dmax*0.95
	  else
	    dmin=dmin-.5
	    dmax=dmax+.5
	  endif
	  endif
	else
c * logarithmic check
	  if(amin.gt.0. .and. amax.gt.0.) then
	    if(amin.ne.amax) go to 400
c * min and max both > 0, but =; so, expand range.
	    dmax=amax
	    dmin=dmax/10.
	  else
c * search array for min and max
	    dmin=finfin
	    dmax=-finfin
	    if(nxbase.gt.0 .and. nbase.eq.ibasex(0)) then
	      call agbmnx(array,dmin,dmax,nxbase)
	    else
	      if(nbase.eq.ibasex(0)) then
	      npts=comget(ibasec+14)
	      else
	      npts=comget(ibasec+4)
	      endif
	      call agmnmx(array,dmin,dmax,npts,itype)
	    endif
c * check that both min and max are > 0.
	    if(dmax.le.0.) dmax=10.
	    if(dmin.le.0.) dmin=dmax/10.
	  endif
	endif
c * set data range in common
	call comset(nbase+11,dmin)
	call comset(nbase+12,dmax)
c
400	return
	end
c%
c----------------------------------------------------------------/loptim
c * determine endpoints and # of ticks for each axis
c
	subroutine loptim(nbase,n)
	common /bppcom/ comget(80)
cdec$	psect  /sgaxi/ noshr
	common /sgaxi/ k1(2),k2(2),nmod(2),gran(2),ibeg(2),dtic(2)
c
c axis tick selection
	dmin=comget(nbase+11)		!data limits
	dmax=comget(nbase+12)
	ismin=comget(nbase+13)		!screen limits
	ismax=comget(nbase+14)
	itype=comget(nbase+15)		!axis type
	iround=nint(comget(nbase+0))	!tick rounding option
	call agseet(i4105,kfac,ic,is)
	dslen=ismax-ismin
	nmjtic=comget(nbase+5)
	if(nmjtic.gt.0) then
c
c user tick selection
	  nmntic=comget(nbase+8)
	  nmntic=max(nmntic,1)
	  if(itype.eq.2) then
	    dmin=log10(dmin)
	    dmax=log10(dmax)
	  endif
	  gran(n)=(dmax-dmin)/nmjtic/nmntic
	  nmod(n)=nmntic
	  k1(n)=dmin/gran(n) * 1.01
	  k2(n)=dmax/gran(n) * 1.01
	  iround=1
c
	elseif(itype.eq.1) then
c
c linear tick selection
	  slngth=abs(dslen)/kfac
	  call aglint(dmin,dmax,slngth,iround,gran(n),nmod(n),
     >		k1(n),k2(n))
	else
c
c logarithmic tick selection
	  slngth=abs(dslen)/kfac
          dmin=log10(dmin)
          dmax=log10(dmax)
	  call aglogt(dmin,dmax,slngth,iround,gran(n),nmod(n),
     >		k1(n),k2(n))
	endif
c
 	call comset(nbase+26,dmin)	!data limits
	call comset(nbase+27,dmax)	! for labeling
c
c set screen position of first tick and tick spacing
	if(iround.gt.0) then		!extend to ticks
	  if(itype.eq.2 .and. nmjtic.eq.0) then
	    dtic(n)=dslen/(dmax-dmin)	!log
	  else
	    dtic(n)=dslen/(k2(n)-k1(n))	!linear
	  endif
	  ibeg(n)=ismin
	else				!extend to data
c
	  tvsm=dslen/(dmax-dmin)	!user to screen transformation
	  tvsb=float(ismin)-tvsm*dmin
	  abeg=roundu(dmin,gran(n))	!if irrational endpoints
	  aend=roundd(dmax,gran(n))	! find position of first tick
	  if(aend.eq.abeg) aend=abeg+1.
	  isbeg=abeg*tvsm+tvsb+0.5
	  isend=aend*tvsm+tvsb+0.5
  	  if(itype.eq.2) then	!log
	    dtic(n)=float(isend-isbeg)/(aend-abeg)
	  else			!linear
	    dtic(n)=float(isend-isbeg)/(k2(n)-k1(n))
	  endif
	  ibeg(n)=isbeg
	endif
c
	call width(nbase,n)	!?put here
	return
	end
c%
c----------------------------------------------------------------/aglint
c * select linear tick parameters
c
	subroutine aglint(dmin,dmax,slngth,iround,gran,nmod,k1,k2)
c
c arguments:
c
        real, intent(inout) :: dmin,dmax ! min & max of data
c               ...may be extended to match min/max tic marks
        real, intent(in) :: slngth      ! axis length in "rasters"
        integer, intent(in) :: iround   ! neatness control
c               ...iround=0 -- use data min & max, all tic marks inside
c                              this range
c               ...iround=1 -- axes limits at minor tic marks
c               ...iround=2 -- axes limits at major tic marks
c
        real, intent(out) :: gran       ! minor tic mark granularity
        integer, intent(out) :: nmod    ! #tic marks / major tic mark
        integer, intent(out) :: k1,k2   ! tics at k1*gran ... k2*gran
c
c  ticks will be drawn from 0.0 to 1.2 by 0.1 with major
c  ticks at 0.0, 0.5, 1.0.
c
	real spans(4), tgran(3,3)
	integer nmodt(3,3)
	character espan*16
c
	data (spans(i),i=1,4)	 /10., 20., 50.,100./	!data span
	data (tgran(i,1),i=1,3)  / 1.,  1.,  2./	!tick granularity
	data (nmodt(i,1),i=1,3)  / 5,   5,   5/		!major tick modulus
	data (tgran(i,2),i=1,3)	 / 2.,  2.,  5./
	data (nmodt(i,2),i=1,3)	 / 5,   5,   2/
	data (tgran(i,3),i=1,3)	 / 5.,  5., 10./
	data (nmodt(i,3),i=1,3)	 / 2,   2,   5/
c
c compute data span
	dspan=dmax-dmin
	if(dspan.eq.0 .or. abs(dmax).gt.abs(1.e5*dspan)) then
	  if(abs(dmin).gt.10.) then
	    dmin=dmin*0.95
	    dmax=dmax*1.05
	  else
	    dmin=dmin-.5
	    dmax=dmax+.5
	  endif
	  dspan=dmax-dmin
	endif
	write(espan,100) dspan
  100	format(2pe12.4)
	read(espan,101) fmant,nexp
  101	format(1x,f7.4,1x,i3)
c
c given span select tick option column
	j=1
	do while (fmant.ge.spans(j+1) .and. j.lt.3)
	  j=j+1
	enddo
c
c find suitable number of ticks based on screen length
	texp=10.**nexp
	!texp=1./10.**(-nexp)	!is also inaccurate for nexp=-2
	do nden=1,3
	  gran=tgran(j,nden)*texp	!tick granularity
	  nmod=nmodt(j,nden)
c
	  if(iround.eq.0) then		!0=raw data
	    rmin=roundu(dmin,gran)		!round in to minor ticks
	    rmax=roundd(dmax,gran)		! for irrational axes	
	  elseif(iround.eq.1) then	!1=minor ticks
	    rmin=roundd(dmin,gran)		!round out to
	    rmax=roundu(dmax,gran)		! minor ticks
	  else				!3=major
	    rmin=roundd(dmin,gran*nmod)		!round out to
	    rmax=roundu(dmax,gran*nmod)		! major ticks
	  endif
	  if(iround.gt.0) then
	    rmin=min(rmin,dmin)		!ensure end points are included
	    rmax=max(rmax,dmax)
	  endif
c
	  k1=rmin/gran+sign(0.1,rmin)	!set tick counts s.t. k*gran
	  k2=rmax/gran+sign(0.1,rmax)	! equals tick value
c
	  if(slngth/(k2-k1).gt.20) go to 500
	enddo
c
  500	if(iround.gt.0) then
	  dmin=rmin
	  dmax=rmax
	endif
	return
	end
c%
c----------------------------------------------------------------/aglogt
c * select log tick parameters
c
	subroutine aglogt(dmin,dmax,slngth,iround,gran,nmod,k1,k2)
c
c  nmod=2, 4, or 5	minor and major at decade intervals
c  nmod=3		minor at 2 and 5, major at decades
c  nmod=9		minor at 2 to 9, major at decades
c arguments:
c
        real, intent(inout) :: dmin,dmax ! min & max of log10(data)
c               ...may be extended to match min/max tic marks
        real, intent(in) :: slngth      ! axis length in "rasters"
        integer, intent(in) :: iround   ! neatness control
c               ...iround=0 -- use data min & max, all tic marks inside
c                              this range -- but show a decade.
c               ...iround=1 or 2 -- axes limits at decades.
c
        real, intent(out) :: gran       ! minor tic mark granularity, log(data)
        integer, intent(out) :: nmod    ! #tic marks / major tic mark
        integer, intent(out) :: k1,k2   ! tics at k1*gran ... k2*gran
c
c          actual data values at 10**(k1*gran) to 10**(k2*gran)
c
	real tablog(9)
	data tablog /	0.,
     >		0.30103, 0.47712, 0.60206, 0.69897,
     >		0.77815, 0.84510, 0.90309, 0.95424/
c
        amin=dmin
        amax=dmax
	if(iround.eq.0) then		!0=raw data
	  rmin=roundu(amin,1.)		!round in for
	  rmax=roundd(amax,1.)		! irrational endpoints
c
c				(sign(1.,amin).eq.sign(1.,amax))
	  if((int(amin).eq.int(amax)) .and. (amin*amax.gt.0.)) then
	    tmin=roundd(amin,1.)	!must cross a decade
	    tmax=roundu(amax,1.)
	    if(abs(amin-tmin).le.abs(amax-tmax)) then
	      rmin=tmin			!move one end to decade
	      amin=tmin
	    else
	      rmax=tmax
	      amax=tmax
	    endif
	  endif
	  beg=abs(amin-roundd(amin,1.))
c
	else				!1=minor ticks, 2=major
cx	  amin=amin	!+0.000001
cx	  amax=amax	!-0.000001
	  rmin=roundd(amin,1.)		!round out to decades
	  rmax=roundu(amax,1.)
	  amin=rmin
	  amax=rmax
	  beg=0.
	endif
	gran=1.				!tick granularity
	itint=slngth/(amax-amin)	! screen interval of one decade
c
	if(itint.lt.60) then
	  k1=rmin			!minors will
	  k2=rmax			! be at decades
	  nmod=2			!tick modulus
	  if(itint.lt.30) nmod=4
	  if(itint.lt.15) nmod=5
c
	elseif(itint.lt.120) then	!1,2,5 ticks
	  k1=0
	  if(beg.gt.0) then		!.gt.tablog(1)
	  if(beg.le.tablog(5)) k1=2	!k1 = first tick
	  if(beg.le.tablog(2)) k1=1
	  endif
	  k2=int(rmax-rmin)*3		!3 ticks per decade
	  if(k1.ne.0) k2=k2+3-k1
	  end=amax-rmax
	  if(end.ge.tablog(2)) k2=k2+1
	  if(end.ge.tablog(5)) k2=k2+1
	  k2=k1+k2			!k2 = last tick
	  nmod=3
c
	else				!1 to 9 ticks
	  k1=0
	  do while(k1.lt.9 .and. beg.gt.tablog(k1+1))
	    k1=k1+1			!first minor tick
	  enddo
	  ntic=int(rmax-rmin)*9		!9 ticks per decade
	  if(k1.ne.0) ntic=ntic+9-k1
	  end=amax-rmax
	  k2=0
	  do while(k2.lt.8 .and. end.ge.tablog(k2+2))
	    k2=k2+1			!minors past last decade
	  enddo
	  k2=ntic+k2+k1			!last (minor) tick
	  nmod=9
	endif
c
	dmin=amin
	dmax=amax
c
	return
	end
c%
c----------------------------------------------------------------/width
c * set label width and remote exponent
c
	subroutine width(nbase,n)
	common /bppcom/ comget(80)
cdec$	psect  /sgaxi/ noshr
	common /sgaxi/ k1(2),k2(2),nmod(2),gran(2),ibeg(2),dtic(2)
	character lit*16
	logical logm
c
	labend=comget(nbase+3)		!which ticks to label
	dmin=comget(nbase+26)
	dmax=comget(nbase+27)
	logm=nmod(n).eq.3 .or. nmod(n).eq.9
	if(.not.logm) then
	  rmaj1=roundu(dmin,gran(n)*nmod(n))	!insure at least
	  rmaj2=roundd(dmax,gran(n)*nmod(n))	! 2 labels
	  if(abs((rmaj1-rmaj2)/gran(n)).lt.1e-3) labend=2	!end only
	endif
	call comset(nbase+28,float(labend))	!work flag vs user request
c
	if(labend.eq.2 .or. logm) then
	  ticint=1000.			!label end ticks only
	  granu=gran(n)			! may be minor
	  if(mod(k1(n),nmod(n)).eq.0 .and. mod(k2(n),nmod(n)).eq.0)
     >		granu=gran(n)*nmod(n)
	else
	  ticint=dtic(n)*nmod(n)	!label major ticks
	  granu=gran(n)*nmod(n)
	endif
c
101	format(1pe12.5)
102	format(9x,i3)
c * calc label width
	rmax=max(abs(dmin),abs(dmax))	!order of larger magnitude
	itype=comget(nbase+15)
	if(itype.eq.1) then	!linear
	  write(lit,101) granu
	  read(lit,102) igexp		!exponent of granularity
	  lsig=4			!least significant digit
	  do while (lsig.gt.0 .and. lit(lsig+3:lsig+3).eq.'0')
	    lsig=lsig-1			!check digits 4,7
	  enddo
	  ndec=max(0,lsig-igexp)	!number of decimal places
	  iwid=ndec+1
	  if(ndec.gt.0) iwid=iwid+1	!allow for point
	  !?irexp=log10(abs(rmax)+.00005)
	  write(lit,101) rmax		!exponent of larger limit
	  read(lit,102) irexp		!number of leading digits
	  if(irexp.ge.0) iwid=iwid+irexp
	  if(dmin.lt.0.) iwid=iwid+1	!allow for sign
	  iexp=0
	  mrexp=comget(nbase+16)
	  if(iwid.gt.max(5,mrexp)) then
	    iexp=irexp			!remote exponent
	    ndec=irexp-igexp+lsig
	    iwid=ndec+1
	    if(ndec.gt.0) iwid=iwid+1	!point
	    if(dmin.lt.0.) iwid=iwid+1	!sign
	  endif
	  mwidth=comget(nbase+16)
	  if (mwidth.ge.5) mwidth=0	!param also width for remote exp
	  if(mwidth.ne.0) then		!user specified width
	    ndec=mwidth-(iwid-ndec)	!decimal places to fill label
	    ndec=max(ndec,0)
	    iwid=mwidth
	  endif
c
	else			!log (itype.eq.2)
 	  iexp=0
	  ndec=0
	  iwid=2			!10
 	  if(rmax.ge.1.) iwid=3		!10 and 1 digit
 	  if(rmax.ge.10.) iwid=4	!10 and 2 digits
	  if(dmin.lt.0.) iwid=iwid+1
	endif
c set label width in common
 	call comset(nbase+17,float(iwid))
	call comset(nbase+10,float(ndec))
	call comset(nbase+18,float(iexp))
c set label spacing
	call agseet(i4105,kfac,ic,is)
	if(n.eq.1) then
	  chspce=iwid*kfac*14	!?ihorz
	else
	  chspce=34		!?ivert*?
	endif
c
	istep=1
	do while(chspce.ge.0.7*ticint*istep .and. istep.le.3)
	  istep=istep+1
	enddo
	call comset(nbase+19,float(istep))
c
	return
	end
c%
c----------------------------------------------------------------/agmnmx
c * determine minimum and maximum values in an array
c
	subroutine agmnmx(array,amin,amax,npta,itype)
	real array(*)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
	logical logs
c
	logs=itype.ge.2
	npts=npta
	go to 100
c
	entry mnmx(array,amin,amax)
c
	logs=.false.
	npts=comget(ibasec+4)
 100	continue
c * determine form of data array
	if(npts.gt.0) then	!non-standard long form
	  nstart=1
	else			!standard long/short form
	  npts=nint(array(1))+1
	  nstart=2
	  if(array(1).lt.0.) then
	      xmax=array(3)+(array(2)-1.)*array(4)
	      amin=min(array(3),xmax,amin)
	      amax=max(array(3),xmax,amax)
	    go to 600
	  endif
	endif
c * data value eq to machine infinity considered missing data
c * ignore values <=0 if log scale
	finfin=comget(ibasec+3)
	do i=nstart,npts
	  if(logs .and. array(i).le.0) go to 500	!(2/90)
	  if(abs(array(i)).lt.finfin) then
	    amin=min(array(i),amin)
	    amax=max(array(i),amax)
	  endif
 500	enddo
 600	return
	end
c%
c----------------------------------------------------------------/setwin
c * set the proper virtual space window transformation
c
	subroutine setwin
	common /bppcom/ comget(80)
	nbase=ibasex(0)
	ixmin=comget(nbase+13)
	ixmax=comget(nbase+14)
	dxmin=comget(nbase+26)
	dxmax=comget(nbase+27)
	itypex=comget(nbase+15)
	if(itypex.eq.2) then
	  dxmin=10.**dxmin
	  dxmax=10.**dxmax
	endif
	nbase=ibasey(0)
	iymin=comget(nbase+13)
	iymax=comget(nbase+14)
	dymin=comget(nbase+26)
	dymax=comget(nbase+27)
	itypey=comget(nbase+15)
	if(itypey.eq.2) then
	  dymin=10.**dymin
	  dymax=10.**dymax
	endif
c * set linear or log transformation parameters in tcs common
	call aglnlg(itypex,itypey)
c * set the screen window
	call twindo(ixmin,ixmax,iymin,iymax)
	call dwindo(dxmin,dxmax,dymin,dymax)
	return
	end
c%
c----------------------------------------------------------------/grid
c * draw axis lines, tick marks and grid lines
c
	subroutine grid
	common /bppcom/ comget(80)
cdec$	psect  /sgaxi/ noshr
	common /sgaxi/ k1(2),k2(2),nmod(2),gran(2),ibeg(2),dtic(2)
	integer ibasec
	parameter (ibasec=1)
	logical logm
	real tablog(9)
	data tablog /0, 0.30103, 0.47712, 0.60206, 0.69897,
     >	  0.77815, 0.84510, 0.90309, 0.95424/
c
	igrid=comget(ibasec+11)		!set line/text color - 10/87
	call linclr(igrid)
	call agtxclr(igrid)
	do n=1,2			!x,y
	  if(n.eq.1) nbase=ibasex(0)
	  if(n.eq.2) nbase=ibasey(0)
	  lform=comget(nbase+7)		!major axis form
	  if(lform.eq.0) go to 300	!no axis
c * get the min max of opposite axis
	  mbase=iother(nbase)
	  ipmin=comget(mbase+13)	!other axis screen min
	  ipmax=comget(mbase+14)
	  near=min(ipmin,ipmax)
	  nfar=max(ipmin,ipmax)
	  iaxpos=near+int(comget(nbase+2))	!axis position
c * draw axis lines
	  ismin=comget(nbase+13)	!screen limits
	  ismax=comget(nbase+14)
	  if(n.eq.1) then
	    call movabs(ismin,iaxpos)	!draw axis line
	    call drwabs(ismax,iaxpos)
	  else
	    call movabs(iaxpos,ismin)	!draw axis line
	    call drwabs(iaxpos,ismax)
	  endif
	  if(lform.eq.1) go to 300	!no ticks
c
	  mjbeg=0
	  mjend=0
	  lmaj=comget(nbase+6)		!length of major ticks
	  if(iaxpos.ge.(nfar+near)/2) then
	    lmaj=-lmaj
	    nfar=near
	  endif
c * find end points of tick marks
	  call tset2(iaxpos,nfar,lmaj,lform,mjbeg,mjend)
	  mform=comget(nbase+9)
	  mnbeg=0
	  mnend=0
	  if(mform.gt.1) then
	    lmin=lmaj/2
c * find end points of minor ticks
	    call tset2(iaxpos,nfar,lmin,mform,mnbeg,mnend)
	  endif
c * draw ticks
	  ntics=0
	  logm=nmod(n).eq.3 .or. nmod(n).eq.9
	  do k=k1(n),k2(n)
	    i=mod(k,nmod(n))
	    if(i.eq.0 .or. (.not.logm)) then
c * linear, or logarithmic at major ticks
	      ixy=ibeg(n)+ifix(ntics*dtic(n))
	      ntics=ntics+1
	    else
c * logarithmic at minor ticks
	      if(nmod(n).eq.9) then
		toff=tablog(i+1)
	      elseif(i.eq.1) then
		toff=tablog(2)
	      else
		toff=tablog(5)
	      endif
	      ixy=ibeg(n)+ifix((ntics-1+toff)*dtic(n))
	    endif
	    if(n.eq.1) then		!x axis
	      if(i.eq.0) then		!major tick
		call movabs(ixy,mjbeg)
		call drwabs(ixy,mjend)
	      elseif(mform.gt.1) then	!minor tick
		call movabs(ixy,mnbeg)
		call drwabs(ixy,mnend)
	      endif
	    else			!y axis
	      if(i.eq.0) then		!major tick
		call movabs(mjbeg,ixy)
		call drwabs(mjend,ixy)
	      elseif(mform.gt.1) then	!minor tick
		call movabs(mnbeg,ixy)
		call drwabs(mnend,ixy)
	      endif
	    endif
	  enddo
c
 300	enddo
	return
	end
c%
c----------------------------------------------------------------/tset2
c * utility routine for grid
c
	subroutine tset2(newloc,nfar,nlen,nfrm,kstart,kend)
	kstart=newloc-nlen
	kend=newloc+nlen
	if(nfrm.eq.3 .or. nfrm.eq.6) kstart=newloc
	if(nfrm.eq.2) kend=newloc
	if(nfrm.eq.5 .or. nfrm.eq.6) kend = nfar
cx	call agseet(i4105,if,ic,maxs)
cx	kstart=max(0,min(maxs,kstart))
cx	kend=max(0,min(maxs,kend))
	return
	end
c%
c----------------------------------------------------------------/label
c * display tick mark labels
c
	subroutine label(nbase)
	common /bppcom/ comget(80)
cdec$	psect  /sgaxi/ noshr
	common /sgaxi/ k1(2),k2(2),nmod(2),gran(2),ibeg(2),dtic(2)
	character clabl*12
	logical logm,putlab
	logical noout(0:6)	!no part of tick outside axis
	data noout /.true.,.true.,.false.,.true.,.false.,.false.,.true./
c
	if(nbase.eq.ibasex(0)) then
	  n=1
	else
	  n=2
	endif
	if(comget(nbase+3).eq.0.) return
	mbase=iother(nbase)
	near=comget(mbase+13)
	nfar=comget(mbase+14)
	iaxoff=comget(nbase+2)
	ilbpos=min(near,nfar)+iaxoff	!axis position
c
	lform=comget(nbase+7)		!major tick form
	lmaj=comget(nbase+6)		!major tick length
	if(noout(lform)) lmaj=7		!?*kfac
	if(ilbpos.lt.(nfar+near)/2) lmaj=-lmaj
	itype=comget(nbase+15)  	!linear/log
	call csize(ihorz,ivert)
	ih2=ihorz/2
	if(n.eq.1) then			!x axis
	  if(lmaj.lt.0) then
	    loff=-ivert*1.1		!x bottom (was 24)
	    if(itype.eq.2) loff=loff-4	!x bottom log
	  else
	    loff=ivert*.46		!x top (10)
	  endif
	else				!y axis
	  if(lmaj.lt.0) then
	    loff=-ih2			!y left (7, was 8)
	  else
	    loff=ih2			!y right (7, was 8)
	  endif
	endif
c * adjust label position by string size and tick length
	ilbpos=ilbpos+lmaj+loff
c
	logm=nmod(n).eq.3 .or. nmod(n).eq.9	!log with minor ticks
	labend=comget(nbase+28)		!label major or end ticks
	lstep=comget(nbase+19)		!label skip
	if(n.eq.2) lstep=1		!ignore for y axis
	if(logm) lstep=1
	nmjtic=comget(nbase+5)		!?user set ticks
	ntics=0
c
	dmin=comget(nbase+26)		!value of first tick
	if(logm) dmin=roundu(dmin,1.)
	iexp=comget(nbase+18)		!remote exponent
	fac=10.**(-iexp)
	iwidth=comget(nbase+17)		!label width
	nlab=0
	do k=k1(n),k2(n)
	  i=mod(k,nmod(n))
	  if(logm) valab=dmin+ntics
c
	  putlab=.false.
	  if(labend.eq.1) then		!label all major
	    putlab=i.eq.0 .and. mod(nlab,lstep).eq.0
	  elseif(labend.eq.2) then	!label ends only
	    if(logm) then
	      putlab=i.eq.0 .and. ((k-k1(n)).lt.nmod(n) .or.
     >		(k2(n)-k).lt.nmod(n))
	    else
	      putlab=k.eq.k1(n) .or. k.eq.k2(n)
	    endif
	  endif
c
	  if(putlab) then
	    ixy=ibeg(n)+ifix(ntics*dtic(n))
	    if(n.eq.2) ixy=ixy-ih2	!?ivert/3 (was 7)
	    if(.not.logm) then
	      if(nmjtic.ne.0) then
	        valab=(dmin+ntics*gran(n))*fac
	      else
	        valab=k*gran(n)*fac
	      endif
	    endif
	    len=iwidth
	    call numset(valab,nbase,clabl,len)
	    if(n.eq.1) then		!x axis
	      call movabs(ixy,ilbpos)
	      call juster(clabl,len,0)
	    else			!y axis
	      call movabs(ilbpos,ixy)
	      call juster(clabl,len,-lmaj)
	    endif
	  endif
	  if(i.eq.0 .or. .not.logm) ntics=ntics+1
	  if(i.eq.0) nlab=nlab+1
 110	enddo
c * put remote exponent
	if(iexp.ne.0) then
	  ismax=comget(nbase+14)
	  if(n.eq.2) then
	    irx=ilbpos
	    iry=ismax+ivert		!(22)
	  else
	    irx=ismax			!move to end of axis (apr84)
	    if(lmaj.lt.0) then
	      iry=ilbpos-ivert*1.8	!x bottom (40)
	    else
	      iry=ilbpos+ivert*1.5	!x top (30)
	    endif
	  endif
	  call remlab(nbase,iaxoff,idum,irx,iry)
	endif
c
	return
	end
c%
c----------------------------------------------------------------/cplot
c * plot a curve
c
	subroutine cplot(x,y)
	real x(*),y(*)
	common /bppcom/ comget(80)
	real yy(4)		! array of y values to be plotted
	integer iyy(4)		! corresponding y-raster positions
	integer iydelt		! max variation in y raster position
	logical none
	integer ibasec
	parameter (ibasec=1)
c
c * get params from common
	line=comget(ibasec)		!line style
	isym=comget(ibasec+1)		!symbol type
	if (line.eq.-1 .and. isym.eq.0)  return		! nothing to do
	isteps=comget(ibasec+2)		!step between symbols
	istepl=comget(ibasec+5)		!step between points in line
cx	if(istepl.eq.0) then
cx	  stop 'STOP in sglib cplot: binitt not called'
cx	endif
	finfin=comget(ibasec+3)
cx	nxbase=comget(ibasec+9)
c * set the curve color for 4105
	icurv=comget(ibasec+10)
	call linclr(icurv)
c * get percent wiggle in y; convert to screen units
	nbasey=ibasey(0)
	yband=comget(nbasey+21)		!max variation in %
	smin=comget(nbasey+13)
	smax=comget(nbasey+14)
	iydelt=(smax-smin)*yband*.01
c * set fast plot line type
	if (icclip(0).eq.0) then
	  if (line.eq.0) then
	    iopt = 1			!solid line, clipped
	  else if (line.gt.0) then
	    iopt = 2			! dashed line, clipped
	  endif
	else
	  if (line.eq.0) then
	    iopt = 3			!solid line, not clipped
	  else if (line.gt.0) then
	    iopt = 4			! dashed line, not clipped
	  endif
	endif
c * set up fast data get
	call agkeys(x,y,ifastx,ifasty,limit)
	if(limit.eq.0) return
	if(ifastx.eq.0 .or. ifasty.eq.0) then
C	  type *, 'CPLOT: invalid data form'
	  call at_msg('CPLOT: invalid data form')
	  return
	endif
	if(ifastx.eq.-1) then
	  xfirs=x(3)
	  xdelt=x(4)
	else	!if(ifastx.eq.-2) then
	  xfirs=x(2)
	  xdelt=x(3)
	endif
	lines=1				!solid
	if(line.gt.0)  lines=3		!dashed
	if(line.eq.-1) lines=2		!symbols only
	if(line.eq.-2 .or. line.eq.-3) then
	  lines=5			!bar
	  isym=0			!comget used for shade code
	endif
	if(line.eq.-4) lines=4	!points
	if(line.eq.-5) then
	  lines=4			!marker
	  call agmark(isym)		!set marker type
	  isym=0
	  if(isteps.gt.0) istepl=istepl*isteps	!put point only at symbols
	  isteps=0
	endif
	if(line.lt.-10) lines=6		!user drawn
c * a regular symbol is plotted unless no symbol or bar chart specified
	i=1
	none=.true.
c * get first non-infinite point
	do while (i.le.limit .and. none)
	  if(ifastx.gt.0) then
	    xpoint=x(ifastx)
	    ifastx=ifastx+istepl
cx	  elseif(ifastx.eq.0) then
cx	    xpoint=agbdgt(x,nxbase,i)
	  else	!if(ifastx.lt.0) then
	    xpoint=xfirs+xdelt*float(i-1)
	  endif
	  if(ifasty.gt.0) then
	    ypoint=y(ifasty)
	    ifasty=ifasty+istepl
cx	  elseif(ifasty.eq.0) then
cx	    call uline
	  else	!if(ifasty.lt.0) then
	    ypoint=y(3)+y(4)*float(i-1)
	  endif
	  i=i+istepl
	  if(abs(xpoint).lt.finfin .and. abs(ypoint).lt.finfin)
     >	    none=.false.
	enddo
c * move to first data location
	call movea(xpoint,ypoint)
c * test if can condense y values that have same x screen coordinate
	if(line.lt.0 .or. isym.ne.0) go to 400
	call wincot(xpoint,ypoint,ix,iy)
	ixa = ix		!initialize raster parameters
	ixd = 0			!x direction
	iydb4 = iy		!end of previous raster
	iyy(1) = iy		!start of raster
	iyd = iy		!last y for this raster
	iylo = iy		!low y at this x
	iyhi = iy		!high y at this x
	xa = xpoint
	yy(1) = ypoint
	yd = ypoint
	ylo = ypoint
	yhi = ypoint
	do while(i.le.limit)
	  if(ifastx.gt.0) then
	    xpoint=x(ifastx)
	    ifastx=ifastx+istepl
cx	  elseif(ifastx.eq.0) then
cx	    xpoint=agbdgt(x,nxbase,i)
	  else	!if(ifastx.lt.0) then
	    xpoint=xfirs+xdelt*float(i-1)
	  endif
	  if(ifasty.gt.0) then
	    ypoint=y(ifasty)
	    ifasty=ifasty+istepl
cx	  elseif(ifasty.eq.0) then
cx	    call uline
	  else	!if(ifasty.lt.0) then
	    ypoint=y(3)+y(4)*float(i-1)
	  endif
	  i=i+istepl
	  if (abs(xpoint).ge.finfin .or. abs(ypoint).ge.finfin) goto 380
	  call wincot(xpoint,ypoint,ix,iy)
c				!if new x value, examine iy rasters
	  if(ix.ne.ixa) then
	    jnum = 1		! number of points to plot
	    ixr = (ix-ixa)*ixd	!see if x values change direction
	    ixd = sign(1,ix-ixa)
	    if (iyhi-iylo.gt.iydelt) then
	      if (iyy(1).gt.iyd) then	!end lower than start
		if(iyhi.ne.iyy(1)) then
		jnum = jnum+1	!go high first
		iyy(jnum) = iyhi
		yy(jnum) = yhi
		endif
		jnum = jnum+1	!then low
		iyy(jnum) = iylo
		yy(jnum) = ylo
	      else			!end higher than start
		if(iylo.ne.iyy(1)) then
		jnum = jnum+1	!go low first
		iyy(jnum) = iylo
		yy(jnum) = ylo
		endif
		jnum = jnum+1	!then high
		iyy(jnum) = iyhi
		yy(jnum) = yhi
	      endif
	      if (iyd.ne.iyy(jnum)) then
		jnum = jnum+1
		yy(jnum) = yd
		iyy(jnum) = iyd
	      endif
	    endif
c see if need to plot
	    if (jnum.gt.1 .or. iy.ne.iydb4 .or. iyy(1).ne.iydb4 .or.
     >		ixr.lt.0) then
	      do j=1,jnum
	        go to (341,342,343,344),iopt
341		call drawa(xa,yy(j))	!solid line, clipped
		go to 350
342		call dasha(xa,yy(j),line) !dashed line, clipped
		go to 350
343		call drwabs(ixa,iyy(j))	!solid line, not clipped
		go to 350
344		call dshabs(ixa,iyy(j),line) !dashed line, not clipped
350	      enddo
	      iydb4 = iyd
	    endif
c
	    ixa = ix	    	!initialize new x-raster values
	    iyy(1) = iy
	    iylo = iy
	    iyhi = iy
	    iyd = iy
	    xa = xpoint
	    yy(1) = ypoint
	    ylo = ypoint
	    yhi = ypoint
	    yd = ypoint
c same x screen coordinate as last point
	  else		!(ix.eq.ixa)
	    iyd = iy		!exit point for this raster
	    yd = ypoint
	    if (iy.lt.iylo) then
	      iylo = iy
	      ylo = ypoint
	    else if (iy.gt.iyhi) then
	      iyhi = iy
	      yhi = ypoint
	    endif
	  endif
 380	  enddo
c  plot last raster for this plot segment
	  jnum = 1		!number of points to plot
	  if (iyy(1).gt.iyd) then	!end lower than start
	    if(iyhi.ne.iyy(1)) then
	    jnum = jnum+1	!go high first
	    iyy(jnum) = iyhi
	    yy(jnum) = yhi
	    endif
	    jnum = jnum+1	!then low
	    iyy(jnum) = iylo
	    yy(jnum) = ylo
	  else	    		!end higher than start
	    if(iylo.ne.iyy(1)) then
	    jnum = jnum+1	!go low first
	    iyy(jnum) = iylo
	    yy(jnum) = ylo
	    endif
	    jnum = jnum+1	!then high
	    iyy(jnum) = iyhi
	    yy(jnum) = yhi
	  endif
	  if (iyd.ne.iyy(jnum)) then
	    jnum = jnum+1
	    yy(jnum) = yd
	    iyy(jnum) = iyd
	  endif
	  do j=1,jnum
	    if(iopt.le.2) then
	      call dasha(xa,yy(j),line)		!dash/solid clipped
	    elseif(iopt.le.4) then
	      call dshabs(ixa,iyy(j),line)	!unclipped line
	    endif
	  enddo
	go to 900
c ** original cplot code for symbols and special line types
400	continue
	if(lines.eq.4) call pointa(xpoint,ypoint)
	if(lines.eq.5) call bar(xpoint,ypoint,line)
cx	if(lines.eq.6) call uline(xpoint,ypoint,1)
	ibets=comget(ibasec+6)
	ibets=ibets-1		!could start at 0 or 1
	if(isym.ne.0 .and. ibets.le.0) then
	  call bsyms(xpoint,ypoint,isym)
	  ibets=isteps
	endif
c * save line assignment in event of "missing data"
	linsav=lines
	do while(i.le.limit)
	  if(ifastx.gt.0) then
	    xpoint=x(ifastx)
	    ifastx=ifastx+istepl
cx	  elseif(ifastx.eq.0) then
cx	    xpoint=agbdgt(x,nxbase,i)
	  else
	    xpoint=xfirs+xdelt*float(i-1)
	  endif
	  if(ifasty.gt.0) then
	    ypoint=y(ifasty)
	    ifasty=ifasty+istepl
cx	  elseif(ifasty.eq.0) then
cx	    call uline
	  else
	    ypoint=y(3)+y(4)*float(i-1)
	  endif
	  i=i+istepl
c * allow plots to skip missing data points
	  if(abs(xpoint).ge.finfin .or. abs(ypoint).ge.finfin) then
c * change lines assignment so next data point is a move
	    if(lines.ne.6) lines=2
	    go to 800
	  endif
	  go to (510,520,530,540,550,560), lines
510	  call drawa(xpoint,ypoint)	!draw solid line
	  go to 600
520	  call movea(xpoint,ypoint)	!move to next point
c * restore line assignment after "missing data"
	  lines=linsav
	  go to 600
530	  call dasha(xpoint,ypoint,line)	!draw dashed line
	  go to 600
540	  call pointa(xpoint,ypoint)	!plot point
	  go to 600
550	  call bar(xpoint,ypoint,line)	!draw a bar
	  go to 800
560	  continue
cx	  call uline(xpoint,ypoint,i)	!user written line
c * see if time to put symbol
600	  if(isym.ne.0) then
	    ibets=ibets-1
	    if(ibets.gt.0) go to 800
c * draw symbol at point
	    call bsyms(xpoint,ypoint,isym)
	    ibets=isteps
	  endif
 800	enddo		!end of 400 outer loop
	if(lines.eq.4) call agmark(0)	!clear marker type
 900	continue
	return
	end
c%
c----------------------------------------------------------------/agkeys
c * set index/limit based on type of data
c
	subroutine agkeys(x,y,ifastx,ifasty,limit)
	real x(*),y(*)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
c
	nptx=comget(ibasec+14)		!number of x points
	nxbase=comget(ibasec+9)
c * determine type of line to be used in this plot
	if(nptx.gt.0) then
	  if(nxbase.eq.0) then
	    ifastx=1			!non-standard
	  else
	    ifastx=-2			!basex short form
	  endif
	  limitx=nptx
	else
	  keyx=x(1)
	  if(keyx.gt.0) then
	    ifastx=2			!standard long form
	    limitx=keyx
	  elseif(keyx.eq.-1) then
	    ifastx=-1			!standard short form
	    limitx=x(2)
	  else			!unknown - error
	    ifastx=0
	    limitx=0
cx	    call agseter
	  endif
	endif
c
	npty=comget(ibasec+4)		!number of y points
	if(npty.gt.0) then
	  ifasty=1			!non-standard
	  limity=npty
	else
	  keyy=y(1)
	  if(keyy.gt.0) then
	    ifasty=2			!standard long form
	    limity=keyy
	  elseif(keyy.eq.-1) then
	    ifasty=-1			!standard short form
	    limity=y(2)
	  else			!keyy.eq.4 = user line
	    ifasty=0
	    limity=npty
	  endif
	endif
	limit=min(limitx,limity)
	end
c----------------------------------------------------------------/agbdgt
c * get the ith x given the x base array and the # of bases in x
c
	function agbdgt(xbase,nxbase,i)
c
	dimension xbase(3,nxbase)
	ib=i
	do nb=1,nxbase
	  nt=xbase(1,nb)
	  if(ib.le.nt) then
c x=t0+dt*(i-1) of the nbth base
 	    agbdgt=xbase(2,nb)+xbase(3,nb)*(ib-1)
	    return
	  endif
	  ib=ib-nt
	enddo
c i is outside the defined x bases
	agbdgt=0
	return
	end
c%
c----------------------------------------------------------------/agbmnx
c * find min and max of an incremental x base array
c
	subroutine agbmnx(array,amin,amax,nxbase)
c
	dimension array(3,nxbase)
	amin=min(amin,array(2,1))
	nlim=array(1,nxbase)-1.
	amax=max(amax,array(2,nxbase)+array(3,nxbase)*nlim)
	return
	end
c%
c----------------------------------------------------------------/roundd
c * round down routine
c
	function roundd(value,fint)
c * determine the number of tick intervals between 0 and value
	test=value/fint
	if(test.gt.0) test=test+1.e-4	!?1./(2**3)
	round=float(ifix(test))*fint
	do while ((round-value).gt.fint*1.e-4)
	  round=round-fint
	enddo
	roundd=round
	return
	end
c%
c----------------------------------------------------------------/roundu
c * round up routine
c
	function roundu(value,fint)
c * determine the number of tick intervals between 0 and value
	test=value/fint
	if(test.lt.0) test=test-1.e-4	!?1./(2**3)
	round=float(ifix(test))*fint
	do while ((value-round).gt.fint*1.e-4)
	  round=round+fint
	enddo
	roundu=round
	return
	end
c%
c----------------------------------------------------------------/frame
c * draw a box around the screen window
c
	subroutine frame
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
	igrid=comget(ibasec+11)
	call linclr(igrid)
	ix=ibasex(0)
	iy=ibasey(0)
	minx=comget(ix+13)
	maxx=comget(ix+14)
	miny=comget(iy+13)
	maxy=comget(iy+14)
	call movabs(maxx,miny)
	call drwabs(maxx,maxy)
	call drwabs(minx,maxy)
	call drwabs(minx,miny)
	call drwabs(maxx,miny)
	return
	end
c%
c----------------------------------------------------------------/bsyms
c * draw symbols for data points
c
	subroutine bsyms(x,y,isym)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
c * is clipping desired?
	if(icclip(0).eq.0) then
c * yes, is x,y in range?
	  call wincot(x,y,ix,iy)
	  if(ivclip(ix,iy).ne.0) go to 500
	endif
	factr=comget(ibasec+7)
	if(isym.ge.-15) then	!was .ge.0 (Nov/90)
c * call general symbol routine
	  call symout(isym,factr)
	else
c * user symbols
cx	  call users(x,y,isym)
	endif
	call movea(x,y)
500	return
	end
c%
c----------------------------------------------------------------/symout
c * put out symbols at the current location
c
	subroutine symout(isym,factr)
	character csym*1
c * get current location
	call seeloc(ix,iy)
	if(isym.gt.127) then
	  call softek(isym)
	elseif(isym.ge.33) then
c * put ascii character
	  call csize(ihorz,ivert)
	  ihorz=float(ihorz)*.3572
	  ivert=float(ivert)*.3182
	  call movrel(-ihorz,-ivert)
	  csym(1:1)=char(isym)
	  call agths(csym,1)
	elseif(isym.le.15) then		!was 11 (Nov/90)
	  call teksym(isym,factr)
	endif
c * return beam to original position
400	call movabs(ix,iy)
	return
	end
c%
c----------------------------------------------------------------/teksym
c * draw standard (plot-10) symbol
c
	subroutine teksym(isym,amult)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
	logical lpanel,lmfill
	real    sadj(15)	!adjustment for proportional sizes
	logical solid(15)	!true if can be filled
	data sadj /1.0,1.25,1.125,1.125,1.25,1.125,1.25,1.25,1.25,1.25,
     >	  1.125,1.0,1.0,1.125,1.125/
	data solid /.true.,.true.,.true.,.true.,.true.,.true.,.false.,
     >	  .true.,.false.,.false.,.true.,.false.,.false.,.false.,.false./
c
	isymb=abs(isym)
	if(isymb.gt.15 .or. isymb.eq.0) return	!was 11 (Nov/90)
	call agseet(i4105,kfac,ic,is)	!kfac = 1 for 1024, 4 for 4096
	rod=8.*amult*kfac
	adjust=sadj(isymb)
	rod=rod*adjust
	ihalf=nint(rod)	!half height/width
	ifull=ihalf+ihalf
c * set line/panel color
	if(i4105.ne.0) then
	  icurv=comget(ibasec+10)
	  call linclr(icurv)
	  call agpnclr(-icurv)
	endif
	lpanel=i4105.ne.0 .and. isym.lt.0 .and. solid(isymb)
	lmfill=i4105.eq.0 .and. isym.lt.0 .and. solid(isymb)
	if(lmfill) call seeloc(ix,iy)
	go to (310,320,330,340,350,360,370,380,390,400,410,420,
     >	  310,340,360),isymb
c * draw circle
 310	ic30=nint(rod*.8660254)		!cos(30)=sin(60)
	is30=nint(rod*.5)		!sin(30)=cos(60)
	call movrel(ihalf,0)		!move right
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(ic30-ihalf,is30)
	call drwrel(is30-ic30,ic30-is30)
	call drwrel(-is30,ihalf-ic30)
	call drwrel(-is30,ic30-ihalf)
	call drwrel(-ic30+is30,is30-ic30)
	call drwrel(-ihalf+ic30,-is30)
	call drwrel(-ic30+ihalf,-is30)
	call drwrel(-is30+ic30,-ic30+is30)
	call drwrel(is30,-ihalf+ic30)
	call drwrel(is30,-ic30+ihalf)
	call drwrel(ic30-is30,-is30+ic30)
	call drwrel(ihalf-ic30,0+is30)
	if(isymb.eq.13) go to 382	!circle with plus
	if(lmfill.and.ihalf.gt.1) then
	  rod=rod-.7071068
	  ihalf=nint(rod)
	  call movabs(ix,iy)
	  go to 310
	endif
	go to 600
c * draw an x
 320	ic45=nint(rod*.75)		!bit more than cos(45)
	it45=ic45+ic45
	call movrel( ic45, ic45)	!move up,right
	if(lpanel) call agbgpnl(1)	!draw panel border
 322	continue
	call drwrel(-it45,-it45)
	call movrel(0,it45)
	call drwrel(it45,-it45)
	call movrel(0,it45)
	go to 600
c * draw a triangle
 330	ic30=nint(rod*.8660254)		!cos(30)
	is30=nint(rod*.5)		!sin(30)
	call movrel(0,ihalf)		!move up
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(-ic30,-is30-ihalf)
	call drwrel(ic30+ic30,0)
	call drwrel(-ic30,ihalf+is30)
	if(lmfill.and.ihalf.gt.1) then
	  rod=rod-1.
	  ihalf=nint(rod)
	  call movabs(ix,iy)
	  go to 330
	endif
	go to 600
c * draw a square
 340	ic45=nint(rod*.7071068)	!cos/sin(45)
	it45=ic45+ic45
	call movrel( ic45, ic45)	!move up,right
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(-it45, 0)
	call drwrel( 0,   -it45)
	call drwrel( it45, 0)
	call drwrel( 0,    it45)
	if(isymb.eq.14) go to 322	!square with x
	if(lmfill.and.ihalf.gt.1) then
	  rod=rod-1.
	  ihalf=nint(rod)
	  call movabs(ix,iy)
	  go to 340
	endif
	go to 600
c * draw a star
 350	ic234=nint(rod*.5878)		!(-)cos(234)
	is234=nint(rod*.8090)		!(-)sin(234)
	ic378=nint(rod*.9511)		!cos(378)
	is378=nint(rod*.3090)		!sin(378)
	call movrel( 0,           ihalf)
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(-ic234,      -is234-ihalf)
	call drwrel( ic378+ic234, is378+is234)
	call drwrel(-ic378-ic378, 0)
	call drwrel( ic234+ic378,-is234-is378)
	call drwrel(-ic234,       ihalf+is234)
	go to 600
c * draw a diamond
 360	call movrel(ihalf,0)		!move right
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(-ihalf, ihalf)
	call drwrel(-ihalf,-ihalf)
	call drwrel( ihalf,-ihalf)
	call drwrel( ihalf, ihalf)
	if(isymb.eq.15) go to 382	!diamond with plus
	if(lmfill.and.ihalf.gt.1) then
	  ihalf=ihalf-1
	  call movabs(ix,iy)
	  go to 360
	endif
	go to 600
c * draw a vertical bar
 370	call movrel(0, ihalf)
	call drwrel(0,-ifull)
	go to 600
c * draw a plus
 380	call movrel(ihalf,0)		!move right
	if(lpanel) call agbgpnl(1)	!draw panel border
 382	continue
	call drwrel(-ifull, 0)
	call movrel( ihalf, ihalf)
	call drwrel( 0    ,-ifull)
	call movrel( ihalf, ihalf)
	go to 600
c * draw up arrow
 390	ihead=nint(kfac*amult*adjust)
	call drwrel(-2*ihead,-6*ihead)
	call drwrel( 4*ihead, 0)
	call drwrel(-2*ihead, 6*ihead)
	call drwrel( 0      ,-ifull)
	go to 600
c * draw down arrow
 400	ihead=nint(kfac*amult*adjust)
	call drwrel(-2*ihead, 6*ihead)
	call drwrel( 4*ihead, 0)
	call drwrel(-2*ihead,-6*ihead)
	call drwrel( 0      , ifull)
	go to 600
c * draw a del
 410	ic30=nint(rod*.8660254)		!cos(30)
	is30=nint(rod*.5)		!sin(30)
	call movrel(0,-ihalf)		!move down
	if(lpanel) call agbgpnl(1)	!draw panel border
	call drwrel(ic30,is30+ihalf)
	call drwrel(-ic30-ic30,0)
	call drwrel(ic30,-ihalf-is30)
	if(lmfill.and.ihalf.gt.1) then
	  rod=rod-1.
	  ihalf=nint(rod)
	  call movabs(ix,iy)
	  go to 410
	endif
	go to 600
c * draw MAC (VT-Pro) circle
 420	if(i4105.eq.0) go to 310	!normal circle if not 4105
	call cmdout('DA')
	iff=1
	if(isym.lt.0) iff=3
	call intout(iff)		!1=frame, 2=fill, 3=frame+fill
	call intout(0)
	call intout(360)		!from 0 to 360
	ihalf=nint(rod*4.)
	call intout(ihalf)		!x/y radius (?diam)
	call intout(ihalf)
	call vecmod			!end esc sequence (for XPLOT)
600	if(lpanel) call endpnl
	return
	end
c%
c----------------------------------------------------------------/bar
c * set up barchart params and call filbox to draw bar
c
	subroutine bar(x,y,line)
	common /bppcom/ comget(80)
	integer ibasec
	parameter (ibasec=1)
c
	call agseet(i4105,kfac,ic,is)
	if(line.ne.0) then		!need to set up bar
	  isymb=comget(ibasec+1)	!type of shading is symbol type
	  ihalf=comget(ibasec+8)*.5	!width of bar is symbol size
	  if(ihalf.lt.2) ihalf=20*kfac	!default is 40 pixels wise
	  lspace=comget(ibasec+7)	!spacing of line is line size
	  if(lspace.le.1) lspace=20*kfac
	  ivrhor=abs(line)		!bars vertical unless line is 3
	  nbase=ibasex(0)
	  minx=comget(nbase+13)
	  maxx=comget(nbase+14)
	  nbase=ibasey(0)
	  miny=comget(nbase+13)
	  maxy=comget(nbase+14)
	  call wincot(0.,0.,ix,iy)
	  ibegx=ix
	  ibegy=iy
	  it=max(minx,maxx)
	  minx=min(minx,maxx)
	  maxx=it
	  it=max(miny,maxy)
	  miny=min(miny,maxy)
	  maxy=it
	endif
c * find point in screen space
	call wincot(x,y,ix,iy)
	if(ivrhor.eq.2) then
c * vertical bars
	  iyl=min(ibegy,iy)
	  iyh=max(ibegy,iy)
	  ixl=min(ix-ihalf,ix+ihalf)
	  ixh=max(ix-ihalf,ix+ihalf)
	else		!ivrhor.eq.3
c * horizontal bars
	  iyl=min(iy-ihalf,iy+ihalf)
	  iyh=max(iy-ihalf,iy+ihalf)
	  ixl=min(ibegx,ix)
	  ixh=max(ibegx,ix)
	endif
c * adjust sides of bar to be within window
	ixl=max(ixl,minx)
	ixh=min(ixh,maxx)
	iyl=max(iyl,miny)
	iyh=min(iyh,maxy)
c * if bar is too small, return without drawing
	if(ixh-ixl.ge.2 .and. iyh-iyl.ge.2) then
	  call filbox(ixl,iyl,ixh,iyh,isymb,lspace)
	endif
	return
	end
c
c----------------------------------------------------------------/filbox
c * draw a box and fill it with shade lines if desired
c
	subroutine filbox(minx,miny,maxx,maxy,isymb,lspace)
c
	logical hor,ver,dup,dwn
	iminx=min(minx,maxx)
	iminy=min(miny,maxy)
	imaxx=max(minx,maxx)
	imaxy=max(miny,maxy)
c * draw box outline
	call movabs(iminx,iminy)
	call drwabs(imaxx,iminy)
	call drwabs(imaxx,imaxy)
	call drwabs(iminx,imaxy)
	call drwabs(iminx,iminy)
	if(iminx.eq.imaxx .or. iminy.eq.imaxy) return
	if(isymb.gt.15 .or. isymb.lt.0) go to 800
c * decode bar type (bit map)
	hor=btest(isymb,0)	!horizontal lines
	ver=btest(isymb,1)	!vertical lines
	dup=btest(isymb,2)	!diagonal up lines (45 deg)
	dwn=btest(isymb,3)	!diagonal down lines (-45 deg)
c * draw horizontal lines
	if(hor) then
	  il=iminy+lspace
	  do while(il.lt.imaxy)
	    call movabs(iminx,il)
	    call drwabs(imaxx,il)
	    il=il+lspace
	  enddo
	endif
c * draw vertical lines
	if(ver) then
	  il=iminx+lspace
	  do while(il.lt.imaxx)
	    call movabs(il,iminy)
	    call drwabs(il,imaxy)
	    il=il+lspace
	  enddo
	endif
c
	if(.not.dup .and. .not.dwn) go to 800
c * set window/scaling to limits of shading area
	ximin=iminx
	ximax=imaxx
	yimin=iminy
	yimax=imaxy
	call twindo(iminx,imaxx,iminy,imaxy)
	call dwindo(ximin,ximax,yimin,yimax)
	dely=ximax-ximin
	space=lspace
c * diagonal lines up or down
	if(dup) then
	  y=yimin-dely+space
	else	!dwn
	  y=yimin+space
	  yimax=yimax+dely
	  dely=-dely
	endif
	do while(y.lt.yimax)
	  call movea(ximin,y)
	  y2=y+dely
	  call drawa(ximax,y2)
	  y=y+space
	enddo
c * restore window/scaling to variables in common
	call setwin		!twindo/dwindo
800	return
	end
c----------------------------------------------------------------/remlab
c * produce off axis labels for linear scale factors
c
	subroutine remlab(nbase,iloc,labtyp,irx,iry)
	parameter (llabl=8)
	character clabl*(llabl)
	common /bppcom/ comget(80)
	iexp=comget(nbase+18)
	if(iexp.eq.0) return
c * set justification according to axis and side
	iposit=-sign(1,iloc-1)		!y axis
	if(nbase.eq.ibasex(0)) iposit=1	!x axis
	call expout(nbase,iexp,clabl,llabl)
c * put out label string
	len=llabl
	call movabs(irx,iry)
	call juster(clabl,len,iposit)
	return
	end
c%
c----------------------------------------------------------------/expout
c * create exponential label
c
	subroutine expout(nbase,iexp,clabl,nchars)
	character clabl*(*)
	common /bppcom/ comget(80)
	integer jup
	data jup/-1/
	clabl=' '
	log=comget(nbase+15)
	nexp=iabs(iexp)
	ic=nchars
	if(nexp.eq.0) then	!10**0
	  clabl(ic:ic)='1'
	  go to 800
	endif
	if(iexp.eq.1) go to 380
c * x10 to exponent processing
	if(nexp.ge.100) then
	  clabl(ic-1:ic)='**'
	  ic=ic-2
	else			!exponents between 1 and 99
	  if(nexp.ge.10) then
	    clabl(ic:ic)=char(mod(nexp,10)+ichar('0'))
	    nexp=nexp/10
	    ic=ic-1
	  endif
	  clabl(ic:ic)=char(nexp+ichar('0'))
	  ic=ic-1
	endif
	if(iexp.lt.0) then	!negative exponent
	  clabl(ic:ic)='-'
	  ic=ic-1
	endif
	clabl(ic:ic)=char(jup)	!up shift character
	ic=ic-1
c * times 10
380	clabl(ic-1:ic)='10'
	ic=ic-2
	if(log.ne.2) then	!x if not log
	  clabl(ic:ic)='X'
	  ic=ic-1
	endif
 800	continue
	return
	end
c%
c----------------------------------------------------------------/numset
c * reformat floating point for axis labeling
c
	subroutine numset(fnum,nbase,clabl,iwidth)
	character clabl*12
	common /bppcom/ comget(80)
	character cform*8 !c@u
	character stars*12
	data stars /'************'/
	idec=comget(nbase+10)
	itype=comget(nbase+15)
	if(itype.eq.2) then	!log scale
	  iexp=fnum+sign(0.00005,fnum)
	  iwidth=iwidth+1	!increase iwidth to allow for up shift
	  call expout(nbase,iexp,clabl,iwidth)
	else			!linear
	  cform=' ' !c@u
	  if(idec.eq.0) then
	    write(cform,'(2h(i,i2.2,1h))') iwidth !   
	    write(clabl,cform,iostat=ierr) nint(fnum) !   
	  else
	    write(cform,'(2h(f,i2.2,1h.,i1,1h))') iwidth,idec !   
	    write(clabl,cform,iostat=ierr) fnum !   
	  endif
	  if(ierr.ne.0) clabl=stars
	endif
	return
	end
c%
c----------------------------------------------------------------/juster
c * output a string justified right(-1), center(0), or left(1)
c
cx	options /nounderflow
	subroutine juster(clabl,lenchr,iposit)
	character clabl*(*)
        integer jup, jdn, ich, hex00
   	integer ilabl(12)	!c@v
cx If anstr is ever converted to accept a character string argument
cx instead of an integer array, copying to ilabl can be omitted.
	data jup/-1/,jdn/-2/
        data hex00/z'ffffff00'/
c * bypass fill characters
	istart=1
	do while(istart.le.lenchr .and. clabl(istart:istart).eq.' ')
	  istart=istart+1
	enddo
	lentru=lenchr-istart+1
c * count special (invisible) characters
	i=istart
	do while(i.le.lenchr)
           ich = ichar(clabl(i:i))
	  ibyte=ior(ich,hex00)
	  if(ibyte.eq.jup .or. ibyte.eq.jdn) lentru=lentru-1
   	  ilabl(i)=ibyte	!c@v
	  i=i+1
	enddo
c * adjust offset according to position code
	if(iposit.lt.0) then		!left justify
	  ioff=0
	elseif(iposit.eq.0) then	!center
	  ioff=-linwdt(lentru*7-2)/14
	else				!right justify
	  ioff=-linwdt(lentru*7-2)/7
	endif
c * move to output position
	if(ioff.ne.0) call movrel(ioff,0)
	iv2=linhgt(1)/2
	do i=istart,lenchr
c					!extend sign bit
	  ibyte=ior(ichar(clabl(i:i)),hex00)
	  if(ibyte.eq.jup .or. ibyte.eq.jdn) then
c * when special character found put out preceeding string
	    call anstr(i-istart,clabl(istart:)) !   
c * if jump up move for superscript
	    if(ibyte.eq.jup) call movrel(0,iv2)
c * if jump down move for subscript
	    if(ibyte.eq.jdn) call movrel(0,-iv2)
	    istart=i+1
	  endif
	enddo
	call anstr(lenchr-istart+1,clabl(istart:)) !   
	return
	end
c%
c----------------------------------------------------------------/softek
c * dummy routine for user created symbol
c
	subroutine softek(isym)
	return
	end
c%
cxc----------------------------------------------------------------/users
cxc * dummy routine for user created symbol
cxc
cx	subroutine users(x,y,i)
cx	return
cx	end
cxc%
c----------------------------------------------------------------/ibasec
c * return pointer to parameter in common portion of bppcom
c
	function ibasec(ioff)
	common /bppcom/ comget(80)
	ibasec=1+ioff	!cvectr+ioff
	return
c%
c----------------------------------------------------------------/ibasex
c * return pointer to parameter in x vector portion of common
c
	entry ibasex(ioff)
	ibasex=17+ioff	!xvectr+ioff
	return
c%
c----------------------------------------------------------------/ibasey
c * return pointer to parameter in y vector portion of common
c
	entry ibasey(ioff)
	ibasey=46+ioff	!yvectr+ioff
	return
	end
c%
c----------------------------------------------------------------/bppblk
        block data bppblk
        real comget
	common /bppcom/ comget(80)
	data comget(13) /1./	!ccurv (in case no binitt)
        end
c----------------------------------------------------------------/iother
c * find the other base
c
	function iother(nbase)
	nx=ibasex(0)
	ny=ibasey(0)
	if(nbase.eq.nx) iother=ny
	if(nbase.eq.ny) iother=nx
	return
	end
c%
c----------------------------------------------------------------/comget
c * read value from common block
c
	real function comget(item)
        real table
	common /bppcom/ table(80)
	comget=table(item)
	return
	end
c%
c----------------------------------------------------------------/comset
c * write value into common block
c
	subroutine comset(item,value)
        implicit none
        integer item
        real    table, value
	common /bppcom/ table(80)
	if(item.ge.1 .and. item.lt.80) table(item)=value
	return
	end

c  sgnco.for	new contour routines derived from ncar library
c
c  oct-89	created module from conrec,clgen,stline,drline
c
c  feb-93	fix calculation of dashed vs solid lines in agculev.
c		modify agcxget/agcyget for non-VMS compilsers. [met]
c
c  jan-98 CAL   Adapt to TRANSP
c  Apr-99 CAL   make portable
C
C  Sep-00 CAL   pointers / entry need rework for f90
c----------------------------------------------------------------/agclini
c
	subroutine agclini(liru,lplu,lamu)
c
	USE sgnco_mod
	implicit none
	integer i,   j
	integer liru,lplu,lamu
	integer nklr,ifill
	integer nums,istep
	integer numd
	integer iklrs(*),isyms(*),idshs(*)
c
	integer ufxy  !
c
	kfill=0
	nsym=0		!no symbols
	do i=1,nlm
	  kmap(i)=1	!init all colors to black
	  kdsh(i)=0	!all solid lines
	  ksym(i)=0	!no symbols
	enddo
	jufxy=0
c
	nir=max(0,liru)
	npl=max(0,lplu)
	lam=80000
	if(lamu.gt.0) lam=lamu
	return
c
	entry agcolrs(nklr,iklrs,ifill)
c
c  nklr=0  to turn off color lines and color fill
c       > 0  is number of entries in the color sequence
c  iklrs      contains the sequence of color codes
c  ifill > 0 to turn on color fill
c
	if(nklr.gt.0) then
	  kmap(1)=0	!base color is white
	  do i=2,nlm
	    j=mod(i-2,nklr)+1
	    kmap(i)=iklrs(j)
	  enddo
	endif
	kfill=ifill
	return
c
	entry agcsymb(nums,isyms,istep)
c
c  isymb	is symbol code
c  istep	is step between symbols
c
	nsym=nums
	if(nums.gt.0) then
	  do i=1,nlm
	    j=mod(i-1,nums)+1
	    ksym(i)=isyms(j)
	  enddo
	  call agsteps(istep,istep/2)
	endif
	return
c
	entry agcdash(numd,idshs)
c
c  idshs	is array of dash codes
c
	if(numd.gt.0) then
	  do i=1,nlm
	    j=mod(i-1,numd)+1
	    kdsh(i)=idshs(j)
	  enddo
	endif
	return
c
	entry agctran(ufxy)
c
	jufxy=1                 !  Dummy for 1 code
	return
	end
c%
c----------------------------------------------------------------/agclins
c
	subroutine agclins(zdata,mxd,x,nx,y,ny,cl,ncl)
	USE sgnco_mod
	implicit none
	integer mxd, nx, ny, ncl
	real, dimension(:), target :: x
	real, dimension(:), target :: y
	real, dimension(:,:)       :: zdata
	real, dimension(:)         :: cl
	integer, dimension(:), allocatable :: ir
	integer ince
cdec$	psect  /conrer/ noshr
	common /conrer/ ince
	integer iai(nlm),iag(nlm)
	real xfram(5),yfram(5)
	integer lvir,lvpl
	save lvir,lvpl
	character ctemp*12
	logical igood,lib_get_vm
	integer  jx,   jy,  lc, lrsize, i, ios
	integer  i4105,kfac,ic
	integer  kolor
	real     xdir, ydir
	real     contr
	integer  nfloat_size
	external agncof
c
c arguments
c	mxd	first (x) dimension of zdata
c	nx	number of data values to be contoured in the x-direction
c	ny	number of data values to be contoured in the y-direction
c
	ince=0		!no error
	xdir=sign(1.,x(2)-x(1))
	jx=3
	do while (ince.eq.0 .and. jx.le.nx)
	  if(sign(1.,x(jx)-x(jx-1)).ne.xdir) then
	    call ageset('agclins: x vector not monotonic',5,1)
	    go to 800
	  endif
	  jx=jx+1
	enddo
	ydir=sign(1.,y(2)-y(1))
	jy=3
	do while (ince.eq.0 .and. jy.le.ny)
	  if(sign(1.,y(jy)-y(jy-1)).ne.ydir) then
	    call ageset('agclins: y vector not monotonic',5,1)
	    go to 800
	  endif
	  jy=jy+1
	enddo
c
	if(lvir.eq.0) then
	  call usr_trnlog(' ','SG_IR_BUFFER',ctemp,lc,igood)
	  read(ctemp,'(bn,i12)',iostat=ios) lvir
	  if(lvir.eq.0) lvir=-1		!no translation
	  call usr_trnlog(' ','SG_PL_BUFFER',ctemp,lc,igood)
	  read(ctemp,'(bn,i12)',iostat=ios) lvpl
	  if(lvpl.eq.0) lvpl=-1		!no translation
	endif
c
	if(nir.eq.0) then
	  nir=max(nx,ny)*4
	  nir=max(nir,nx*ny/16)
	  nir=max(nir,1024)
	  if(lvir.gt.0) nir=lvir
	endif
	igood = .true.
	allocate (ir(nir), stat=i)
	if (i .ne. 0) igood=.false.
c
c set x,y scaling
	xorg=x(1)
	xrang=1.
	if(nx.gt.1) xrang=xrang/float(nx-1)
	xscal=(x(nx)-x(1))
	nxs=nx
	yorg=y(1)
	yrang=1.
	if(ny.gt.1) yrang=yrang/float(ny-1)
	yscal=(y(ny)-y(1))
	nys=ny
	lax => x
	lay => y
cx	call xneat(0)
cx	call yneat(0)
	call dlimx(x(1),x(nx))
	call dlimy(y(1),y(ny))
c
	call agseet(i4105,kfac,ic,is)
	if (i4105.eq.0) kfill=0
	if(kfill.ne.0 .or. nsym.ne.0) then
	  if(npl.eq.0) then
	    npl=max(nx,ny)*4
	    npl=max(npl,512)
	    if(lvpl.gt.0) npl=lvpl
	  endif
	  igood=.true.
	  allocate (pl1(npl), stat=i)
	  allocate (pl2(npl), stat=i)
	  if (i .ne. 0) igood=.false.
	  if (.not.igood) then
	    call ageset('agclins: no virtual memory',6,1)
	    go to 800
	  endif
	endif
	if (kfill.ne.0) then
	   igood=.true.
	   allocate (am(lam), stat=i)
	   if (i .ne. 0) igood=.false.
	  if (.not.igood) then
	    call ageset('agclins: no virtual memory',6,1)
	    go to 800
	  endif
c
	  call arinam(am,lam)
	  call ageget(1,ince,0)
	  if(ince.gt.0) go to 800
	  xfram(1)=x(1)	!define frame counterclockwise
	  xfram(2)=x(1)
	  xfram(3)=x(nx)
	  xfram(4)=x(nx)
	  xfram(5)=xfram(1)
	  yfram(1)=y(1)
	  yfram(2)=y(ny)
	  yfram(3)=y(ny)
	  yfram(4)=y(1)
	  yfram(5)=yfram(1)
c		!outside (left) is -1, inside (right) is area 0
	  call aredam(am,xfram,yfram,5,1,-1,0)
	else
	endif
c
	ince=0		!no error
c
c find major and minor lines - removed
c
cx	nml=0
cx	if(ilag.ne.0) then
cx	  call reord(cl,ncl,rwork,nml,nulbll+1)
cx	endif
c
c process each level
c
cx	fpart=.5
	mir=0		!high water mark for ir
	mpl=0		!high water mark for pl
	do i=1,ncl
	  contr=cl(i)
	  if (i4105.ne.0) then
	    kolor=max(kmap(i),1)
	    call linclr(kolor)
	    call agccurv(kolor)
	  endif
	  call symbl(ksym(i))
	  call line(kdsh(i))
	  igr=i		!group
	  iar=0		!area
c
c draw all lines at this level
c
	  call agcstln(zdata,mxd,nx,ny,contr,ir)
	  if(ince.gt.0) go to 800
	enddo
c
c  find relative minimums and maximums and mark values - removed
c
cx	if(nhi.eq.0) then
cx	  call minmax(zdata,mxd,nx,ny,isizem,ash,ioffdt)
cx	elseif(nhi.gt.0) then
cx	  call minmax(zdata,mxd,nx,ny,isizem,-ash,ioffdt)
cx	endif
c
	if (kfill.ne.0) then
	  call arpram(am,1,1,1)
	  call ageget(1,ince,0)
	  if(ince.gt.0) go to 800
	  npu=0
	  call arscam(am,pl1,pl2,npl,iai,iag,nlm,agncof)
	endif
	if (i4105.ne.0) then
	  call linclr(1)	!black
	endif
c
 800	continue
	deallocate (ir)
	if (allocated(pl1)) then
	  deallocate(pl1)
	  deallocate(pl2)
	endif
	  if (allocated(am)) deallocate(am)
	return
	end
c----------------------------------------------------------------/agclgen
c
	subroutine agclgen(z,mxd,nx,ny,cclo,chi,cinc,cl,ncl,icnst)
	implicit none
cx	save
	integer   mxd,nx,ny,ncl,icnst
	real      cclo,chi,cinc
	integer   nlm, nla
	parameter (nlm=40)
	parameter (nla=16)
	real cl(nlm), z(mxd,ny)
	real      clo, glo, ha, fanc, crat, p, cc
	integer   i,   j,  k,    kk
c
c agclgen puts the values of the contour levels in cl.
c variable names match those in conrec, with the following additions.
c	flo	value of lowest contour level; if flo=fhi=0., a value
c		rounded up from the minimum zdata is generated
c	fhi	value of highest contour level
c	finc	>0 => increment between contour levels
c		=0 => "nice" value (between 10 and 30) generated
c		<0 => abs(finc)=number of levels
c         ncl     -number of contour levels put in cl.
c         icnst   -flag to tell conrec if a constant field was detected.
c                 .icnst=0 means non-constant field.
c                 .icnst non-zero means constant field.
c
c to produce non-uniform contour level spacing, replace the code in this
c routine with code to produce whatever spacing is desired.
c
	icnst=0
	clo=cclo
	glo=clo
	ha=chi
	fanc=cinc
	crat=nla
	if(ha.lt.glo) then
	  glo=ha
	  ha=clo
	elseif(ha.eq.glo) then
	  if (glo.ne.0.) then
	    cl(1)=glo
	    ncl=1
	    return
	  endif
	  glo=z(1,1)
	  ha=z(1,1)
	  do j=1,ny
	    do i=1,nx
	      glo=min(z(i,j),glo)
	      ha=max(z(i,j),ha)
	    enddo
	  enddo
	  if (glo.ge.ha) then
	    icnst=1
	    ncl=1
	    cclo=glo
	    return
	  endif
	endif
	if (fanc.lt.0) crat=max(1.,-fanc)
	if (fanc.le.0) then
	  fanc=(ha-glo)/crat
	  p=10.**(ifix(alog10(fanc)+5000.)-5000)
	  fanc=max(aint(fanc/p),1.)*p
	endif
c
	if (chi.eq.clo) then
	  glo=aint(glo/fanc)*fanc
	  ha=aint(ha/fanc)*fanc*(1.+sign(1.e-6,ha))
	endif
	do k=1,nlm
	  cc=glo+float(k-1)*fanc
	  if (cc.gt.ha) go to 118
	  kk=k
	  cl(k)=cc
	enddo
  118	ncl=kk
	cclo=cl(1)
	chi=cl(ncl)
	cinc=fanc
	return
	end
c----------------------------------------------------------------/agcstln
c
	subroutine agcstln(z,ll,mm,nn,conv,ir)
	USE sgnco_mod
	implicit none
	integer ll, mm, nn
	real    z(ll,nn)
	integer ir(*)
	real    conv
c
c this routine finds the beginnings of all contour lines at level conv.
c first the edges are searched for lines intersecting the edge (open
c lines) then the interior is searched for lines which do not intersect
c the edge (closed lines).  beginnings are stored in ir to prevent
c retracing of lines.  if ir is filled, an error is output.
c
cdec$	psect  /conrer/ noshr
	integer ince
	common /conrer/ ince
	integer iybits
	parameter (iybits=16)	!?12
	logical drawh
	integer i,j,l,m,n
	integer ixypak, ixx, iyy
	integer jp1, ixy, ip1
	integer inn
	ixypak(ixx,iyy)=ishft(ixx,iybits)+iyy	!ishift
c
	l=ll
	m=mm
	n=nn
	cv=conv
	kir=0
	iss=.true.	!lines intersect frame
	drawh=n.le.m	!draw horizontally if wider than higher
	inir=7
	if(drawh) inir=1
c
	do jp1=2,n		!bottom edge, left to right
 	  if (z(1,jp1).lt.cv .and. z(1,jp1-1).ge.cv) then
	    ix=1
	    iy=jp1-1
	    is=3
	    call agcdrln(z,l,m,n,ir)
	  endif
	  if(ince.gt.0) return
	enddo
	do ip1=2,m		!right edge, bottom to top
	  if (z(ip1,n).lt.cv .and. z(ip1-1,n).ge.cv) then
	    ix=ip1-1
	    iy=n
	    is=5
	    call agcdrln(z,l,m,n,ir)
	  endif
	  if(ince.gt.0) return
	enddo
	do jp1=n,2,-1		!top edge, right to left
	  if (z(m,jp1-1).lt.cv .and. z(m,jp1).ge.cv) then
	    ix=m
	    iy=jp1
	    is=7
	    call agcdrln(z,l,m,n,ir)
	  endif
	  if(ince.gt.0) return
	enddo
	do ip1=m,2,-1		!left edge, top to bottom
	  if (z(ip1-1,1).lt.cv .and. z(ip1,1).ge.cv) then
	    ix=ip1
	    iy=1
	    is=1
	    call agcdrln(z,l,m,n,ir)
	  endif
	  if(ince.gt.0) return
	enddo
	iss=.false.	!lines will close inside
	if(drawh) then
	do j=n-1,2,-1		!vertical
	  do ip1=m,2,-1		!horizontal
	  if (z(ip1-1,j).lt.cv .and. z(ip1,j).ge.cv) then
	    ixy=ixypak(ip1,j)
	    call agcirck(ir,kir,ixy,inn)
	    if(inn.gt.0) go to 107		!cycle
	    call agcirst(ir,kir,nir,ixy,ince)
	    if(ince.gt.0) return
	    mir=max(mir,kir)
	    ix=ip1
	    iy=j
	    is=1
	    call agcdrln(z,l,m,n,ir)
	    if(ince.gt.0) return
	  endif
 107	  enddo
	enddo
	else	!(.not.drawh)
	do i=m-1,2,-1		!horizontal
	  do jp1=n,2,-1		!vertical
	  if (z(i,jp1-1).lt.cv .and. z(i,jp1).ge.cv) then
	    ixy=ixypak(i,jp1)
	    call agcirck(ir,kir,ixy,inn)
	    if(inn.gt.0) go to 109		!cycle
	    call agcirst(ir,kir,nir,ixy,ince)
	    if(ince.gt.0) return
	    mir=max(mir,kir)
	    ix=i
	    iy=jp1
	    is=7
	    call agcdrln(z,l,m,n,ir)
	    if(ince.gt.0) return
	  endif
 109	  enddo
	enddo
	endif
	return
	end
c----------------------------------------------------------------/agcdrln
c
	subroutine agcdrln(z,l,mm,nn,ir)
	USE sgnco_mod
	USE agcGet_mod
	implicit none
	integer  l,mm,nn
	real z(l,nn)
	integer ir(*)
c
c this routine traces a contour line when given the beginning by agcstln.
c transformations can be added by deleting the statement functions for
c fx and fy in agcdrln and minmax and adding external functions.
c x=1. at z(1,j), x=float(m) at z(m,j). x takes on non-integer values.
c y=1. at z(i,1), y=float(n) at z(i,n). y takes on non-integer values.
c
cdec$	psect  /conrer/ noshr
	integer ince, iybits
	common /conrer/ ince
	integer inx(8),iny(8)
	data inx /-1,-1,0,1,1,1,0,-1/
	data iny /0,1,1,1,0,-1,-1,-1/
	logical insect,closed
c
	integer ixypak,ixx,iyy 
	integer m,n,ix0,iy0,idx,idy, ix2,iy2,is0
	real    y,x, fxx, fyy
	parameter (iybits=16)	!?12
	real    c, p1, p2
	ixypak(ixx,iyy)=ishft(ixx,iybits)+iyy	!ishift
	c(p1,p2)=(p1-cv)/(p1-p2)
c
	m=mm
	n=nn
	ix0=ix
	iy0=iy
	is0=is
	idx=inx(is)
	idy=iny(is)
	if (idx.ne.0) then
	  y=iy
	  x=c(z(ix,iy),z(ix+idx,iy))*float(idx)+float(ix)
	else
	  x=ix
	  y=c(z(ix,iy),z(ix,iy+idy)) *float(idy)+float(iy)
	endif
	if(linx) then
	  fxx=xorg + (x-1.)*xrang*xscal
	else
	  call agcxget(x,fxx,lax,nxs)
	endif
	if(liny) then
	  fyy=yorg + (y-1.)*yrang*yscal
	else
	  call agcyget(y,fyy,lay,nys)
	endif
	if(jufxy.ne.0) then
	  x=fxx
	  y=fyy
	  call agcutrn(uufxy,x,y,fxx,fyy)
	endif
	call agcfrst (fxx,fyy)
	closed=.false.
	insect=.false.
	do while (.not.closed .and. .not.insect)
	  is=is+1
	  if (is.gt.8) is=is-8
	  idx=inx(is)
	  idy=iny(is)
	  ix2=ix+idx
	  iy2=iy+idy
	  if (iss) then		!if frame
	    insect=ix2.gt.m .or. iy2.gt.n .or. ix2.lt.1 .or. iy2.lt.1
	    if (insect) go to 120	!leave
	  endif
	  if (cv.le.z(ix2,iy2)) then
	    is=is+4
	    ix=ix2
	    iy=iy2
	    go to 106	!cycle
	  endif
	  if (is/2*2.eq.is) go to 106	!cycle
	  if (idx.ne.0) then
	    y=iy
	    x=c(z(ix,iy),z(ix+idx,iy))*float(idx)+float(ix)
	  else
	    x=ix
	    y=c(z(ix,iy),z(ix,iy+idy))*float(idy)+float(iy)
	  endif
	  if(linx) then
	    fxx=xorg + (x-1.)*xrang*xscal
	  else
	  call agcxget(x,fxx,lax,nxs)
	  endif
	  if(liny) then
	    fyy=yorg + (y-1.)*yrang*yscal
	  else
	     call agcyget(y,fyy,lay,nys)
	  endif
	  if(jufxy.ne.0) then
	    x=fxx
	    y=fyy
	  call agcutrn(uufxy,x,y,fxx,fyy)
	  endif
	  call agcvect (fxx,fyy)
	  if(ince.gt.0) return
	  if (is.eq.inir) then
	    call agcirst(ir,kir,nir,ixypak(ix,iy),ince)
	    if(ince.gt.0) return
	    mir=max(mir,kir)
	  endif
	  if (.not.iss) then
	    closed=ix.eq.ix0 .and. iy.eq.iy0 .and. is.eq.is0
cx	    if (closed) go to 120	!leave
	  endif
 106	enddo
c
c end of line
c
 120	continue
	iar=iar+1		!area number
	if (closed.or.insect) then
	  ilft=0
	  irgt=iar
	else
	  ilft=iar
	  irgt=0
	endif
	call agclast
	return
	end
c
c----------------------------------------------------------------/agcxget
	subroutine agcxget(xind,fxx,x,nxs)
	implicit none
	real     xind, fxx
	integer  nxs
	real, dimension(:), pointer :: x
	integer ixx,ixn
c
	ixx=xind		!1<=xind<=nxs
	ixn=ixx+1
	if(ixn.le.1) then
	  fxx=x(1)
	elseif(ixx.ge.nxs) then
	  fxx=x(nxs)
	else
	  fxx=x(ixx)+(xind-float(ixx))*(x(ixn)-x(ixx))
	endif
	return
	end
c
	subroutine agcyget(yind,fyy,y,nys)
	implicit none
	real     yind, fyy
	integer  nys
	real, dimension(:), pointer :: y
	integer iyy,iyn
c
	iyy=yind		!1<=yind<=nys
	iyn=iyy+1
	if(iyn.le.1) then
	  fyy=y(1)
	elseif(iyy.ge.nys) then
	  fyy=y(nys)
	else
	  fyy=y(iyy)+(yind-float(iyy))*(y(iyn)-y(iyy))
	endif
	return
	end
c
c----------------------------------------------------------------/agcfrst
	subroutine agcfrst(x,y)
c
	USE sgnco_mod
	implicit none
cdec$	psect  /conrer/ noshr
	integer ince
	common /conrer/ ince
	real    comget
	common /bppcom/ comget(80)
	integer ibasec
	data ibasec /1/
	real    x,   y
	integer kpl
	integer line, i, in
	save kpl, line
c
	if (kfill.ne.0 .or. nsym.ne.0) then
	  kpl=0
c	  !? store ix,iy instead of x,y
	  kpl=kpl+1
	  pl1(kpl)=x
	  pl2(kpl)=y
	  mpl=max(mpl,kpl)
	else
	  line=comget(ibasec)
	  call movea(x,y)
	  kpl=0
	endif
	return
c----------------------------------------------------------------/agcvect
c
	entry agcvect(x,y)
	if(kpl.gt.0) then
	  if(kpl.ge.npl) then
	    call ageset('agcdrln: pl buffer overflow',4,1)
	    ince=4	!call ageget(1,ince,0)
	    return
	  endif
	  kpl=kpl+1
	  pl1(kpl)=x
	  pl2(kpl)=y
	  mpl=max(mpl,kpl)
	else
	  if(line.le.0) then
	  call drawa(x,y)
	  else
	  call dasha(x,y,line)
	  endif
	endif
	return
c----------------------------------------------------------------/agclast
c
	entry agclast
	if(kpl.gt.0) then
	  call npts(kpl)
	  call cplot(pl1,pl2)
	  if(kfill.ne.0) then
	  call aredam(am,pl1,pl2,kpl,igr,ilft,irgt)
	  call ageget(1,ince,0)
	  endif
	endif
	return
	end
c
	subroutine agcirst(ir,kir,nir,ixy,ince)
	implicit none
	integer ir(*)
	character cnir*40
	integer kir,nir,ixy,ince, ios, in, i
	if (kir.ge.nir) then
	  write(cnir,'("agclins: ir buffer overflow",i12)',
     >          iostat=ios) nir
	  call ageset(cnir,2,1)
	  ince=2
	  return
	endif
	kir=kir+1
	ir(kir)=ixy
	return
c
	entry agcirck(ir,kir,ixy,in)
	do i=1,kir
	  if(ir(i).eq.ixy) then
	    in=1
	    return
	  endif
	enddo
	in=0
	return
	end
c
	subroutine agcplst(pl1,pl2,kpl,x,y)
	implicit none
	real    pl1(*),pl2(*)
	real    x,     y
	integer kpl
	kpl=kpl+1
	pl1(kpl)=x
	pl2(kpl)=y
	end
c
	subroutine agncof(xcs,ycs,ncs,iai,iag,nai)
	USE sgnco_mod
	implicit none
	integer ncs
	real xcs(*),ycs(*)
	integer iai(*),iag(*),nai
	integer i, kolor, itop
	integer minx,maxx,miny,maxy,ixt,iyt
c
	itop=0
	do i=1,nai
	  if(iai(i).gt.0) itop=max(itop,iag(i))
	enddo
	if(itop.eq.0) return
	kolor=kmap(itop)
	if(kolor.eq.0) return
	if(npu.eq.0) call seetw(minx,maxx,miny,maxy)
	ixt=nint(xcs(1)*float(maxx-minx)+float(minx))
	iyt=nint(ycs(1)*float(maxy-miny)+float(miny))
	call agpnclr(-kolor)
	call begpnl(ixt,iyt,0)
	do i=2,ncs
	  ixt=nint(xcs(i)*float(maxx-minx)+float(minx))
	  iyt=nint(ycs(i)*float(maxy-miny)+float(miny))
	  call movabs(ixt,iyt)
	enddo
	call endpnl
	npu=max(npu,ncs)
	return
	end
c%
	subroutine agcutrn(ufxy,x,y,fxx,fyy)
	USE  agcu_mod
	implicit none
	real x,y,fxx,fyy
	integer  ufxy           !   
	call agcuxy(x,y,fxx,fyy) !    -c@s
	end
c%
c----------------------------------------------------------------/rcontr
	subroutine rcontr(k1,c,k2,a,maxa,xarray,imin,imax,istep,
     >	  yarray,jmin,jmax,jstep)
	use sgint_mod
	implicit none
	integer nlm,maxa,k1,k2,imin,imax,istep,jmin,jmax,jstep
	real a(maxa,*),c(*),xarray(*),yarray(*)
	parameter (nlm=40)
	real cl(nlm)
	real, dimension(:),   allocatable, target :: xar, yar
	real, dimension(:,:), allocatable  :: zar
	integer :: istat
	logical igood,lib_get_vm
	integer  nfloat_size
	integer nx,ny,nxy,i,ii,i9,j,jj,j9,lrsize,ncl
c  assumes that a(i,j)=fcn(xarray(i),yarray(j))
	nx=(imax-imin+istep)/istep
	ny=(jmax-jmin+jstep)/jstep
	nxy=nx*ny
	allocate(zar(nx,ny),stat=istat)
	allocate(xar(nx),stat=istat)
	allocate(yar(ny), stat=istat)
	ii=0
	do i=imin,imax,istep
	  ii=ii+1
	  xar(ii)=xarray(i)
	enddo
	jj=0
	do j=jmin,jmax,jstep
	  jj=jj+1
	  yar(jj)=yarray(j)
	enddo
	ii=0
	j9=0
	do j=jmin,jmax,jstep
	   j9=j9+1
	   i9=0
	  do i=imin,imax,istep
	     i9=i9+1
	     zar(i9,j9)=a(i,j)
	  enddo
	enddo
c
	call agclini(0,0,0)
	call agculev(k1,c,k2,zar,nxy,cl,ncl)
	call agclins(zar,nx,xar,nx,yar,ny,cl,ncl)
c
	deallocate(xar)
	deallocate(yar)
	deallocate(zar)
	end
c%
c----------------------------------------------------------------/contur
	subroutine contur(k1,c,k2,a,maxa,xarray,imin,imax,istep,
     >	  yarray,jmin,jmax,jstep)
	USE sgInt_mod
	USE agcu_mod
	implicit none
	integer nlm,k1,k2,maxa,imin,imax,istep,jmin,jmax,jstep
	real a(maxa,*),c(*)
	parameter (nlm=40)
	real cl(nlm)
	real xarray(maxa,jmax),yarray(maxa,jmax)
	real, dimension(:),   allocatable :: xar, yar
	real, dimension(:,:), allocatable :: zar
	logical igood,lib_get_vm
	integer nx,ny,nxy,i,ii,i9,j,j9,lrsize,ncl
	integer istat
	integer nfloat_size
c  assumes that a(i,j)=fcn(xarray(i,j),yarray(i,j))
	nx=(imax-imin+istep)/istep
	ny=(jmax-jmin+jstep)/jstep
	nxy=nx*ny
	allocate(zar(nx,ny),stat=istat)
	allocate(xar(nx),stat=istat)
	allocate(yar(ny),stat=istat)
	call agclxy(xar,nx,imin)
	call agclxy(yar,ny,jmin)
	ii=0
	j9=0
	do j=jmin,jmax,jstep
	   i9=0
	   j9=j9+1
	  do i=imin,imax,istep
	     i9=i9+1
	     zar(i9,j9) = a(i,j)
	  enddo
	enddo
	call agcuset(xarray,yarray,maxa,nx,ny)
c
	call agclini(0,0,0)
	call agculev(k1,c,k2,zar,nxy,cl,ncl)
	call agctran(0)
	call agclins(zar,nx,xar,nx,yar,ny,cl,ncl)
c
	deallocate(xar)
	deallocate(yar)
	deallocate(zar)
	return
	end
c%
	subroutine agclxy(xar,nx,imin)
	real xar(*)
	do i=1,nx
	  xar(i)=i+imin-1	!?*istep
	enddo
	end
c%
	subroutine agcuset(xarray,yarray,maxa,nx,ny)
	USE sgInt_mod
	implicit none
	integer :: maxa,nx,ny
	real, dimension(maxa,*), target :: xarray
	real, dimension(maxa,*), target :: yarray
	real, dimension(:,:), allocatable, save :: arrx
	real, dimension(:,:), allocatable, save :: arry
	integer :: i,j,maxy
	integer nxx,nyy,maxaa
	real xraw,yraw,xdat,ydat
	save nxx,nyy,maxaa
	nxx=nx
	nyy=ny
	maxaa=maxa
	if (allocated (arrx)) deallocate(arrx)
	if (allocated (arry)) deallocate(arry)
	maxy=ny
	allocate(arrx(maxa,maxy))
	allocate(arry(maxa,maxy))
	do i=1,maxa
	   do j=1,maxy
	      arrx(i,j) = xarray(i,j)
	      arry(i,j) = yarray(i,j)
	   end do
	end do
	return
c
	entry agcuxy(xraw,yraw,xdat,ydat)
	call agcuxyi(arrx,arry,maxaa,nxx,nyy,xraw,yraw,xdat,ydat)
	return
	end
c%
	subroutine agcuxyi(xarray,yarray,maxa,nx,ny,xraw,yraw,xdat,ydat)
	implicit none
	integer maxa,nx,ny 
	real, dimension(maxa,*) :: xarray
	real, dimension(maxa,*) :: yarray
	real xraw,yraw,xdat,ydat
	integer ix,iy,ixn,iyn
	real    xfrac,yfrac,xdat1,xdat2,ydat1, ydat2
	ix=max(int(xraw),1)
	iy=max(int(yraw),1)
	ixn=min(ix+1,nx)
	iyn=min(iy+1,ny)
	xfrac=xraw-float(ix)
	yfrac=yraw-float(iy)
	xdat1=xarray(ix,iy)*(1.-xfrac)+xarray(ixn,iy)*(xfrac)
	xdat2=xarray(ix,iyn)*(1.-xfrac)+xarray(ixn,iyn)*(xfrac)
	xdat=xdat1*(1.-yfrac)+xdat2*yfrac
	ydat1=yarray(ix,iy)*(1.-yfrac)+yarray(ix,iyn)*(yfrac)
	ydat2=yarray(ixn,iy)*(1.-yfrac)+yarray(ixn,iyn)*(yfrac)
	ydat=ydat1*(1.-xfrac)+ydat2*xfrac
	end
c%
	subroutine agculev(k1,c,k2,zar,nz,cl,ncl)
	implicit none
	integer k1,k2,nz,ncl
	real c(*),zar(*),cl(*)
	integer nlm
	parameter (nlm=40)
	integer idash(nlm)
	integer i, nsolid, ndashd
	real    cmin,cmax,dff
	if(k1.eq.0) then	!user specified levels
	  call agmnmx(zar,cmin,cmax,nz,1)
	  ncl=0
	  nsolid=nint((cmax-c(1))/c(2))
	  nsolid=min(nsolid,nlm-1)
	  ndashd=nint((c(1)-cmin)/c(2))
	  ndashd=min(ndashd,nlm-1-nsolid)
	  do i=ndashd,1,-1
	    ncl=ncl+1
	    cl(ncl)=c(1)-(i)*c(2)
	    idash(ncl)=1
	  enddo
	  do i=0,nsolid
	    ncl=ncl+1
	    cl(ncl)=c(1)+(i)*c(2)
	    idash(ncl)=0
	  enddo
	  call agcdash(ncl,idash)
	elseif(k1.gt.0) then
	  ncl=min(k1,nlm)
	  do i=1,ncl
	    cl(i)=c(i)
	  enddo
	  if(k2.gt.0) then
	    do i=1,k2-1
	      idash(i)=1
	    enddo
	    do i=k2,ncl
	      idash(i)=0
	    enddo
	    call agcdash(ncl,idash)
	  endif
	else		!k1.lt.0
	  ncl=min(-k1,40)
	  dff=(c(2)-c(1))/(ncl-1)
	  do i=1,ncl
	    cl(i)=c(1)+(i-1)*dff
	  enddo
cx	  cl(ncl)=c(2)
	  do i=2,ncl
	    c(i)=cl(i)
	  enddo
	endif
	end
c%
c%
c----------------------------------------------------------------/agraph
c  draw graph and map virtual space for contour plot
c
	subroutine agraph(minx,maxx,miny,maxy,xmin,xmax,ymin,ymax,
     >	  xlabl,ylabl)
	implicit none
	integer minx,maxx,miny,maxy
	real    xmin,xmax,ymin,ymax
	character xlabl*(*), ylabl*(*)
	real x(2), y(2)
c
	call binitt
	call xfrm(2)		!short ticks
	call yfrm(2)
	call xneat(0)		!frame data exactly
	call yneat(0)
	call slimx(minx,maxx)	!screen limits
	call slimy(miny,maxy)
	call dlimx(xmin,xmax)	!data limits
	call dlimy(ymin,ymax)
	call npts(0)
	call check(x,y)
	call agdspl(x,y)	!draw left/bottom axes, label ticks
	call xloctp(0)		!setup for top/right axes
	call ylocrt(0)
	call grid		!put ticks on top/right, no tick labels
	call xloc(0)		!back to bottom/left for agt routines
	call yloc(0)
	call agtchs(xlabl,0)	!put user labels
	call agtcvs(ylabl,0)
	return
	end
c%
c----------------------------------------------------------------/agrcfr
c  draw grid with ticks which just encloses a contour plot
c
	subroutine agrcfr(x,imin,imax,istep,y,jmin,jmax,jstep,igrid)
	implicit none
	real    x(*), y(*)
	integer imin,imax,istep,jmin,jmax,jstep,igrid
	real    xx(2),yy(2)
	integer irem
	call xfrm(2)		!short ticks, outside
	call yfrm(2)
	call xneat(0)		!frame data exactly
	call yneat(0)
	call npts(2)
	xx(1)=x(imin)		!first x value
	irem=mod((imax-imin),istep)
	xx(2)=x(imax-irem)
	yy(1)=y(jmin)
	irem=mod((jmax-jmin),jstep)
	yy(2)=y(jmax-irem)
	call check(xx,yy)
	call agdspl(xx,yy)	!draw grid&ticks, label ticks
	if(igrid.eq.2) then
	  call xloctp(0)	!set up for top/right grid
	  call ylocrt(0)
	  call grid		!draw grid&ticks, no tick labels
	  call xloc(0)		!back to bottom/left for agt routines
	  call yloc(0)
	else
	  call frame		!just draw box
	endif
	return
	end
c----------------------------------------------------------------/agcfram
c  draw grid with ticks which just encloses a contour plot
	subroutine agcfram(x,nx,y,ny,igrid)
	implicit none
	real x(*),y(*)
	real xx(2),yy(2)
	integer nx,ny,igrid 
	call xfrm(2)		!short ticks, outside
	call yfrm(2)
	call xneat(0)		!frame data exactly
	call yneat(0)
	call npts(2)
	xx(1)=x(1)		!first x value
	xx(2)=x(nx)
	yy(1)=y(1)
	yy(2)=y(ny)
	call check(xx,yy)
	call agdspl(xx,yy)	!draw grid, ticks, labels
	if(igrid.eq.2) then
	  call xloctp(0)	!set up for right/top grid
	  call ylocrt(0)
	  call grid		!draw grid/ticks (no tick labels)
	  call xloc(0)		!back to bottom/left for agt routines
	  call yloc(0)
	else
	  call frame		!just draw box
	endif
	return
	end

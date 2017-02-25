c@zeff.f  Program to compute impurity densities given Z_eff, etc
c  by Glenn Bateman, Lehigh University, 29 December 1996
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c  Given in namelist zin:
c  zeff(i)   = 1 + sum_j Z_j^2 n_Zj / n_e, at breakpoint time i
c  scalezeff(i) = scale factor for zeff(i)
c  nimpset   = index of impurity density to determine from Z_eff
c  zimp(i,j) = charge state of j-th impurity, j=1,N_Z
c  dene(i)   = electron density
c  scalene(i) = scale factor for dene(i) (default all 1.0)
c  rimp(i,j) = denz(j) / dene, ratio of impurity to electron density
c              (for j .ne. nimpset)
c  rhyd(i,j) = denh(j) / denh(1), ratio of jth to 1st hydrogen density
c  cdenh     = output name of denh (default bdhyde)
c  cdenz     = output name of denz (default bdimpe)
c
c  Compute:
c  nhyd      = number of hydrogen ions
c  nimp      = number of impurities
c  denz(i,j) = impurity densities, j=1,nimp
c  denh(i,j) = hydrogenic densities, j=1,nhyd
c
c  Input data given in namelist zin in the standard input file
c  Output to standard output in a format suitable for insertion
c  directly into the BALDUR namelist input file
c
c  ntmax = number of breakpoint times, is determined by
c    zeff(i) > 1.0 and dene(i) > 1.0 for 1 .le. i .le. ntmax
c  nhyd = number of hydrogen ions, is determined by
c    rhyd(i,j) > epslon for 1 .le. j .le. nhyd
c    for at least one value of 1 .le. i .le. ntmax
c  nimp = number of impurities, is determined by
c    zimp(i,j) > 1.0 for all i with 1 .le. i .le. ntmax
c    1 .le. j .le. nimp
c
c Note, with 1 impurity
c n_H / n_e = ( Z - Z_eff ) / ( Z - 1 )
c n_Z / n_e = ( Z_eff - 1 ) / [ Z ( Z - 1 ) ]
c n_Z / n_H = ( Z_eff - 1 ) / [ Z ( Z - Z_eff ) ]

c With 2 impurities
c Z_eff - 1 - Z_2 ( Z_2 - 1 ) n_Z2 / n_e = Z_1 ( Z_1 - 1 ) n_Z1 / n_e
c n_H / n_e = 1 - Z_1 n_Z1 / n_e - Z_2 n_Z2 / n_e
c Z_eff -1 = Z_1 ( Z_1 - Z_eff ) n_Z1 / n_H
c              + Z_2 ( Z_2 - Z_eff ) n_Z1 / n_H
c
c  Routines used:  ctrim
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
      implicit none
c
      integer  ktdim, ksdim, nimp, nhyd, nimpset
c
      parameter ( ktdim = 30, ksdim = 4 )
c
      real  zeff(ktdim), scalezeff(ktdim), zimp(ktdim, ksdim)
     &    , dene(ktdim), scalene(ktdim)
     &    , rimp(ktdim, ksdim), rhyd(ktdim, ksdim)
     &    , denz(ktdim, ksdim), denh(ktdim, ksdim)
     &    , zdene(ktdim)
     &    , epslon, ztemp, ztemp1, ztemp2, zdiff, zmax, zdens
c
      character *8 cdenh, cdenz, ctemp*15
c
      integer  ntmax, ifirst, ilast, itemp, j, jt, js
c
      namelist /zin/  nimpset
     &  , zeff, scalezeff, zimp, dene, scalene, rimp, rhyd
     &  , cdenh, cdenz
c
c..defaults
c
      nimp = 1
      nhyd = 1
      nimpset = 1
      ntmax = 0
c
      epslon = 1.e-10
c
      do jt=1,ktdim
        zeff(jt) = 0.0
        scalezeff(jt) = 1.0
        dene(jt) = 0.0
        scalene(jt) = 1.0
        do js=1,ksdim
          zimp(jt,js) = 0.0
          rimp(jt,js) = 0.0
          rhyd(jt,js) = 0.0
          denz(jt,js) = 0.0
          denh(jt,js) = 0.0
        enddo
      enddo
c
      cdenh = 'bdhyde'
      cdenz = 'bdimpe'
      ctemp = ' '
c
c..Input
c
      open ( unit=8, file='temp', status='unknown', err=90 )
c
      call stripx ( 5, 8, 6 )
c
      read ( 8, zin, err=92, end=92 )
c
c..scale electron density
c
      do jt=1,ktdim
        zeff(jt)  = scalezeff(jt) * zeff(jt)
        zdene(jt) = scalene(jt) * dene(jt)
      enddo
c
c..find number of timesteps ntmax
c  from number of zeff(jt) > 1.0 and zdene(jt) > 1.0
c
      ntmax = 0
      do jt=1,ktdim
        if ( zeff(jt) .lt. 1.0 ) go to 2
        if ( zdene(jt) .lt. 1.0 ) go to 2
        ntmax = jt
      enddo
   2  continue
c
      if ( ntmax .lt. 1 ) then
        write (*,*) ntmax,' = ntmax .lt. 1 ### Abort'
        stop
      endif
c
c..find the number of hydrogen species
c
      nhyd = 1
      do js=2,ksdim
        do jt=1,ntmax
          if ( rhyd(jt,js) .gt. epslon ) then
            nhyd = js
            go to 3
          endif
        enddo
        go to 4
   3    continue
      enddo
   4  continue
c
c..find the number of impurity species
c
      nimp = 0
      do js=1,ksdim
        ztemp  = 1.0
        ztemp1 = 0.0
        do jt=1,ntmax
          ztemp  = min ( ztemp,  zimp(jt,js) )
          ztemp1 = max ( ztemp1, rimp(jt,js) )
        enddo
        if ( ztemp .lt. 0.999 ) go to 6
        if ( js .ne. nimpset .and. ztemp1 .lt. epslon ) go to 6
        nimp = js
      enddo
   6  continue
c
c..check for input errors
c
      if ( nimpset .gt. nimp .or. nimp .lt. 1 ) then
        write (*,*) 'Abort: ',nimpset,' = nimpset .gt. nimp .or.'
        write (*,*) 'Abort: ',nimp,' = nimp .lt. 1'
        stop
      endif
c
      if ( nhyd .lt. 1 ) then
        write (*,*) 'Abort: ',nhyd,' = nhyd .lt. 1'
        stop
      endif
c
c..compute output
c
      do jt=1,ntmax
        ztemp = zeff(jt) - 1.0
        do js=1,nimp
c
          if ( zimp(jt,js) .lt. 1.01 ) then
            write (*,*) 'Abort: ',zimp(jt,js)
     &        ,' = zimp(',jt,js,') .lt. 1.01'
            stop
          endif
c
          if ( js .ne. nimpset ) then
            ztemp = ztemp
     &        - zimp(jt,js) * ( zimp(jt,js) - 1.0 ) * rimp(jt,js)
            denz(jt,js) = rimp(jt,js) * zdene(jt)
          endif
        enddo
c
        rimp(jt,nimpset) = ztemp
     &    / ( zimp(jt,nimpset) * ( zimp(jt,nimpset) - 1.0 ) )
c
        denz(jt,nimpset) = ztemp * zdene(jt)
     &    / ( zimp(jt,nimpset) * ( zimp(jt,nimpset) - 1.0 ) )
c
c..compute ztemp = sum_js denh(jt,js) / zdene(jt), js=1,nhyd
c  ztemp2 = sum_js rhyd(jt,js)
c
        ztemp = 1.0
c
        do js=1,nimp
          ztemp = ztemp - zimp(jt,js) * rimp(jt,js)
        enddo
c
        rhyd(jt,1) = 1.0
        ztemp2 = 1.0
c
        if ( nhyd .gt. 1 ) then
          do js=2,nhyd
            ztemp2 = ztemp2 + rhyd(jt,js)
          enddo
        endif
c
        denh(jt,1) = ztemp * zdene(jt) / ztemp2
c
        if ( nhyd .gt. 1 ) then
          do js=2,nhyd
            denh(jt,js) = denh(jt,1) * rhyd(jt,js)
          enddo
        endif
c
      enddo
          
c
c..output suitable for inclusion in BALDUR input
c
      write (*,100)
 100  format (' ! densities computed by the zeff code, given:')
c
      write (*,101) ' zeff(1) = ', (zeff(jt),jt=1,ntmax)
 101  format (' !',a,5(1pe10.4,', '),
     &  /(' !',11x,5(1pe10.4,', ')))
c
      write (*,101) ' zdene(1) = ', (zdene(jt),jt=1,ntmax)
c
c
c
      call ctrim ( cdenh, ifirst, ilast )
      ctemp = ' ' // cdenh(ifirst:ilast) // '(1,'
      itemp = ilast - ifirst + 5
      do js=1,nhyd
        write (*,110) ctemp(1:itemp), js, (denh(jt,js),jt=1,ntmax)
 110    format ( a,i1,') = ',5(1pe10.4,', '),
     &    /(11x,5(1pe10.4,', ')))
      enddo
c
      call ctrim ( cdenz, ifirst, ilast )
      ctemp = ' ' // cdenz(ifirst:ilast) // '(1,'
      itemp = ilast - ifirst + 5
      do js=1,nimp
        write (*,110) ctemp(1:itemp), js, (denz(jt,js),jt=1,ntmax)
      enddo
c
c
c..check correctness of results
c
c  zdens = computed electron density
c  ztemp = computed Z_eff
c
      do jt=1,ntmax
        ztemp = 0.0
        do js=1,nhyd
          zdens = zdens + denh(jt,js)
        enddo
        do js=1,nimp
          zdens = zdens + zimp(jt,js) * denz(jt,js)
          ztemp = ztemp + zimp(jt,js)**2 * denz(jt,js)
        enddo
        if ( abs ( dene(jt) - zdens ) .lt. dene(jt) * epslon ) then
          write (*,*) 'Error: computed zdens = ',zdens
          write (*,*) '             dene(jt) = ',dene(jt)
          write (*,*) '               for jt = ',jt
        endif
        ztemp = 1.0 + ztemp / zdens
        if ( abs ( ztemp - zeff(jt) ) .lt. epslon ) then
          write (*,*) 'Error: computed Z_eff = ',ztemp
          write (*,*) '             zeff(jt) = ',zeff(jt)
          write (*,*) '               for jt = ',jt
        endif
      enddo
c
c
      stop
c
  90  continue
c
      write (*,*)
      write (*,*) 'Abort:  unable to open temp file  '
      write (*,*) 'or      unable to open input file '
      stop
c
   92 continue
      write (*,*)
      write (*,*) 'Abort: unable to read namelist input'
      stop
c
   94 continue
      write (*,*)
      write (*,*) 'Abort: error before end of input list'
      stop
c
 98   continue
      stop
      end
c@ctrim.f   29 December 1996   Glenn Bateman, Lehigh University
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine ctrim ( string, nfirst, nlast)
c
c     This routine finds the first non-blank string of characters
c
c  string  = input character string
c  nfirst  = first non-blank character in string
c  nlast   = last non-blank character in string
c
      implicit none
c
      character string*(*)
c
      integer nfirst, nlast
c
      integer ichar1, ichar2, imax, j
c
c  imax  = number of characters in the input string
c
      imax = len ( string )
c
      nfirst = 0
      nlast  = 0
c
      if ( imax .lt. 1 ) return
c
c..find the first non-blank character in string
c
      do j=1,imax
        if ( string(j:j) .ne. ' ' ) then
          nfirst = j
          go to 10
        endif
      enddo
  10  continue
c
      if ( nfirst .lt. 1 ) return
c
c..find the last non-blank character in string
c
      nlast = nfirst
      do j=nfirst,imax
        if ( string(j:j) .eq. ' ' ) go to 20
        nlast = j
      enddo
  20  continue
c
      return
      end
c--------1---------2---------3---------4---------5---------6---------7-c
c@stripx   /baldur/code/bald/stripx.f
c  rgb 25-aug-96 clean the data set to conform to namelist requirements
c  rgb 16-aug-96 insert commas as needed in namelist data
c  rgb 06-oct-94 changed line length from 80 to 132 characters
c
      subroutine stripx (nin,ntrans,nout)
c
c  This sbrtn reads input data file logical unit number nin
c  one line at a time up to 80 characters long,
c  prints out the input verbatum to output file logical unit number nout
c  then searches for the first appearance of the character ! on each line
c  and prints out everything to the left of ! to output unit number ntrans
c
c  In this way, everything to the right of ! is treated as a comment
c
      parameter (kc = 132)
c
      character line * 132
c
      logical lcomma, lequil, lspace, lnlist
c
c..lcomma  a comma has been encountered
c  lequil  an "=" has been encountered
c  lspace  the last character was a space
c  lnlist  this line is in a namelist
c
c..iline = line number
c
      iline = 0
      lnlist = .false.
c
  10  read (nin,100,end=20) line
 100  format (a132)
c
      iline = iline + 1
c
c..find number of characters before spaces
c
      ilength = 0
      do j=1,kc
cap      
        if (( line(j:j) .ne. ' ' ) .and. ( line(j:j) .ne. char(9) )) 
     >	   ilength = j
      enddo
c
c..echo input data on output file
c
      if ( ilength .gt. 0 ) then
        write (nout,110) line(1:ilength)
      else
        write (nout,*)
      endif
c
c..ixlen = number of nontrivial characters before "!"
c  ieq   = position of "=" in line
c
      ixlen = 0
      ieq   = 1
      do j=1,kc
        if ( line(j:j) .eq. '!' ) go to 14
cap	
        if (( line(j:j) .ne. ' ' ) .and. ( line(j:j) .ne. char(9) )) 
     >	    ixlen = j
        if ( line(j:j) .eq. '=' ) ieq = j
      enddo
  14  continue
c
c..clean the data set to conform to namelist requirements
c
      i1 = 2
      i2 = min( max(ixlen+1,i1), kc )
c
      if ( index(line(1:ixlen+1),'$end') .gt. 1
     &   .or. index(line(1:ixlen+1),'/') .gt. 1 ) then
        lnlist = .false.
        go to 19
      elseif ( ( index(line(1:ixlen+1),'$') .gt. 1 )
     &    .or. ( index(line(1:ixlen+1),'&') .gt. 1 )
     &    ) then
        lnlist = .true.
        go to 19
      endif
c
c
      if ( lnlist 
     &  .and. ( index(line(1:ixlen+1),'$') .lt. 1 )
     &  .and. ( index(line(1:ixlen+1),'&') .lt. 1 )
     &  ) then
c
      lequil = .false.
      lcomma = .false.
      lspace = .true.
c
      do 18 j=i2,i1,-1
cap
        select case (line(j-1:j-1))
	  case (' ', char(9)) 
	    lspace = .true.
	    cycle
	  case ('=')
	    lequil = .true.
	    cycle
	  case (',')
	    lcomma = .true.
	    cycle	            	    
	end select
c
        if ( lequil .or. lcomma ) then
          lequil = .false.
          lcomma = .false.
          lspace = .false.
          go to 18
        endif
c
        if ( lspace ) then
          line(j:j) = ','
          lspace = .false.
          if ( j .gt. ixlen ) ixlen = j
          go to 18
        endif
c
  18  continue
c
      endif
c
c--------1---------2---------3---------4---------5---------6---------7-c
c
c
c..echo input data up to "!" on temporary file
c
  19  continue
c
      if ( ixlen .gt. 0 ) then
        write (ntrans,110) line(1:ixlen)
      else
cbate        write (ntrans,*)
      endif
c
 110  format (a)
c
c
      go to 10
c
  20  continue
      rewind ntrans
c
      return
      end

c--------1---------2---------3---------4---------5---------6---------7-c
c@cnvcof  .../baldur/code/bald/dsolver.f
c  rgb 10-jan-96 redefine i2=mxchi*mxchi before do 216
c  rgb 15-aug-94 compute flxtot and srctot as well as vftot
c  rgb 10-aug-94 save local variables that are initialized
c  rgb 06-may-90 18.34 compute ctotei(jb) and ctotii(jb) power fluxes
c  rgb 02-may-90 cleaned up style
c      removed condition on tai for initialization
c  fgps 24-feb-83 modified to handle more than 2 impurity species.
cdoc
c=======================================================================
c       ------------
c       sbrtn CNVCOF   file DSOLVER
c       ------------
c
c      2.16    convert aaaa, bbbb, cccc, dddd and eta to internal units
c
c
c      common blocks and variables modified:
c
c       aaaa, bbbb, cccc, dddd, eta (comcoe)
c
c-----------------------------------------------------------------------
c
c
c       aaaa, bbbb, cccc, and dddd are multiplied, element by element,
c       by zconva, zconva, zconvc, and zconvd, respectively
c       eta is multiplied by zconve
c
c       aaaa(i,j) has units  l / t * u(i)/u(j), as does bbbb(i,j)
c       cccc(i,j) has units  u(i)/u(j)  / t
c       dddd(i)   has units  u(i) / t
c       eta       has units  l**2 / t
c       where u(i) is the unit of the ith parameter in chi
c
cend
c**********************************************************************c
c
        subroutine cnvcof
c
        include 'cparm.m'
        include 'cbaldr.m'
        include 'cbparm.m'
        include 'commhd.m'
c
c
        common /cnvect/ vftot(55,9), flxtot(55,9), srctot(55,9)
c
c  vftot(jz,ji)  = total effective convective velocities
c  flxtot(jz,ji) = total fluxes of particles and energies
c  srctot(jz,ji) = total sources of particles and energies
c
        dimension
     r   zconva(id2chi)     , zconvc(id2chi)     , zconvd(idxchi)
c
        logical
     l  ilinit
c
        data    ilinit /.false./
c
        save  i2, zconva, zconvc, zconvd, ilinit
c
c-----------------------------------------------------------------------
c
c
        data    iclass /2/,     isub /16/
c
c
      if ( nlomt2(isub) )  then
        call mesage(' *** 2.16 subroutine cnvcof bypassed')
        return
      endif
c
c
c-----------------------------------------------------------------------
c
        if ( ilinit ) go to 200
c
c
c      1)      initialize zconva, zconvc, zconvd, zconve
c
c
        if (mxchi.gt.12) go to 9010
        i2 = mxchi*mxchi
        call resetr (zconvc,i2,uist)
c
        do 104 j1 = lhyd1, limpn
          i001 = lelec + mxchi*(j1-1)
          i002 = lion  + mxchi*(j1-1)
          i003 = j1 + mxchi*(lelec-1)
          i004 = j1 + mxchi*(lion -1)
          zconvc(i001) = uist * usie
          zconvc(i002) = uist * usie
          zconvc(i003) = uist * uise
          zconvc(i004) = uist * uise
          zconvd(j1) = uist * usid
  104   continue
c
        zconvd(lelec) = uist * usid * usie
        zconvd(lion)  = uist * usid * usie
c
        do 108 j1 = 1, i2
          zconva(j1) = zconvc(j1) * usil
  108   continue
c
c eta(internal) = emu0 * c**2 / (4 pi) * (units conversion) * eta(standard)
c eta(internal) = usir * eta(standard)
c
c       since in internal units
c
c       db/dt = (d/dr) eta d/(r dr) r*b / emu0
c
c       but in standard units
c
c       db/dt = c**2/(4 pi) (d/dr) eta d/(r dr) r*b
c
c
        ilinit = .true.
  200 continue
c
c
c      2)      convert aaaa, bbbb, cccc, and dddd
c
cbate      i2 = mxchi*mxchi
c
      do jz = 1, mzones
c
        do j = 1, mxchi
          do j2 = 1, mxchi
            j3 = j + (j2-1)*mxchi
            aaaa(j,j2,jz) = aaaa(j,j2,jz) * zconva(j3)
            bbbb(j,j2,jz) = bbbb(j,j2,jz) * zconva(j3)
            cccc(j,j2,jz) = cccc(j,j2,jz) * zconvc(j3)
          enddo
        enddo
c
        do j = 1, mxchi
          dddd(j,jz)   = dddd(j,jz)   * zconvd(j)
        enddo
c
        etai(jz) = etai(jz) * usir
c
      enddo
c
c..compute total electron and ion power fluxes
c
      do jz = 2, mzones
c
        zsurfi = avi(jz,3,1) * avi(jz,5,1)
        ctotei(jz) = 0.0
        ctotii(jz) = 0.0
c
        do js=1,mchi
          ctotei(jz) = ctotei(jz) -
     &               ( aaaa(lelec,js,jz) * chi(js,jz-1)
     &               + bbbb(lelec,js,jz) * chi(js,jz) ) * zsurfi
          ctotii(jz) = ctotii(jz) -
     &               ( aaaa(lion ,js,jz) * chi(js,jz-1)
     &               + bbbb(lion ,js,jz) * chi(js,jz) ) * zsurfi
        enddo
c
      enddo
c
c..compute the total effective convective velocities
c
      do jz=lcentr,mzones
        do j1=1,mchi
          vftot(jz,j1)  = 0.0
          flxtot(jz,j1) = 0.0
          srctot(jz,j1) = dddd(j1,jz)
          do j2=1,mchi
            vftot(jz,j1) = vftot(jz,j1) - 2.0
     &        * ( aaaa(j1,j2,jz) * chi(j2,jz-1)
     &          + bbbb(j1,j2,jz) * chi(j2,jz) )
     &        / ( chi(j1,jz-1) + chi(j1,jz) )
            flxtot(jz,j1) = flxtot(jz,j1)
     &        - ( aaaa(j1,j2,jz) * chi(j2,jz-1)
     &          + bbbb(j1,j2,jz) * chi(j2,jz) )
            srctot(jz,j1) = srctot(jz,j1)
     &        + cccc(j1,j2,jz) * chi(j2,jz)
          enddo
        enddo
      enddo
c
c
        return
c
c
c       90)     errors
c
c
 9010   continue
        call error_olymp(0,iclass,isub,1,
     &         'mxchi .gt. 12. redimension internal arrays')
        return
        end

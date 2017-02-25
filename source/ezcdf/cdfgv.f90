! automatic conversion to free f90 compatible form
! free.pl cdfgv.for
! linewidth: 72
! file names: cdfgv.for
!
!23456789-123456789-123456789-123456789-123456789-123456789-123456789-23
      SUBROUTINE cdfgv(ncid,varnam,varid,dimlens,sizes,type,status)
!
!     Get Variable id, etc.
!     02/11/99 C.Ludescher
! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
!
      implicit none
      include "netcdf.inc"
! Input
      integer                 :: ncid
      character*(*)           :: varnam
      character*1             :: type
      integer, dimension(:)   :: sizes
! Output
      integer, dimension(:)   :: dimlens
      integer                 :: varid, status
! Local
      integer                 :: i, vartyp, ndims, atts, rank
      integer, dimension(3)   :: dimids
      character*(nf_max_name) :: name
 
      status = nf_inq_varid(ncid,varnam,varid)
      call handle_err(status,varnam,'cdfgv','nf_inq_varid')
      if (status .ne. 0) return
      status = nf_inq_var(ncid,varid,name,vartyp,ndims,dimids,atts)
      call handle_err(status,varnam,'cdfgv','nf_inq_var')
      if (status .ne. 0) return
! Verify input dimension is correct
      rank = size(sizes)
      status = 1
      if (ndims .eq. 3 .and. rank .ne. 3) then
         print "('% cdfgv: --E-- The variable ',a,                      &
     &         ' is 3 dimensional')",varnam
         return
      else if (ndims .eq. 2 .and. rank .ne. 2) then
         print "('% cdfgv: --E-- The variable ',a,                      &
     &         ' is 2 dimensional')",varnam
         return
      endif
      if (ndims .eq. 0 .and. sizes(1) .ne. 0) then
         print "('% cdfgv: --E-- The variable ',a,                      &
     &           ' is a Scalar')",varnam
         return
      else if (ndims .eq. 1 .and. rank .ne. 1 ) then
         print "('% cdfgv: --E-- The variable ',a,                      &
     &           ' is 1 dimensional')",                                 &
     &         varnam
         return
      endif
! Verify data type is matching
      select case (type)
      case ('i')
         if(vartyp .ne. nf_int) then
            print "('% cdfgv: --E-- ',a,' is of type Integer !')",      &
     &             varnam
            return
         end if
      case ('d')
         if(vartyp .ne. nf_double) then
            print "('% cdfgv: --E-- ',a,' is of type REAL*8 !')",       &
     &            varnam
            return
         end if
      case ('r')
         if(vartyp .ne. nf_real) then
            print "('% cdfgv: --E-- ',a,' is of type default REAL !')",       &
     &            varnam
            return
         end if
      case ('c')
         if(vartyp .ne. nf_char) then
            print "('% cdfgv: --E-- ',a,' is of type Character !')",         &
     &            varnam
            return
         end if
      end select
      status = 0
      do i=1,ndims
         dimlens(i) = 0
         status = nf_inq_dim(ncid,dimids(i),name,dimlens(i))
         call handle_err(status,varnam,'cdfgv','nf_inq_dim')
      end do
! Check array size is big enough
      select case (ndims)
      case (1)
         if(dimlens(1) .gt. sizes(1) ) then
            print "('% cdfgv: --W-- Output array size =',I6,/           &
     &           '                is smaller than ',                    &
     &           a,' size =',I6/,                                       &
     &           '                output will be truncated !')",        &
     &         sizes(1),varnam,dimlens(1)
      endif
      case (2)
      if(dimlens(1) .gt. sizes(1) .or. dimlens(2) .gt. sizes(2)) then
         print "('% cdfgv: --W-- Output array size =',I6,' *',I6,/      &
     &           '                is smaller than ',                    &
     &           a,' size =',I6,' *',I6/,                               &
     &           '                output will be truncated !')",        &
     &         sizes(1),sizes(2),varnam,dimlens(1),dimlens(2)
      endif
      case (3)
      if(dimlens(1) .gt. sizes(1) .or. dimlens(2) .gt. sizes(2)         &
     &   .or. dimlens(3) .gt. sizes(3)) then
         print "('% cdfgv: --W-- Output array size =',                  &
     &          I5,' *',I5,' *',I5/,                                    &
     &           '                is smaller than ',                    &
     &           a,' size =',I5,' *',I5,' *',I5/,                       &
     &           '                output will be truncated !')",        &
     &         sizes(1),sizes(2),sizes(3),varnam,                       &
     &         dimlens(1),dimlens(2),dimlens(3)
      endif
      end select
      end

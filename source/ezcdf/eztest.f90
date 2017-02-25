! Sample to Write / Read netCDF dataset
! C. Ludescher/A. Pletzer Tue Apr  4 10:11:33 EDT 2000
      USE ezcdf
      IMPLICIT NONE
      integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
      integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
! Variables to be written
      character*25                   :: title =                         &
     &                                  "Test of EZcdf Interface"
      character*8, dimension(2)      :: comment =                       &
     &                                  (/"Written ","02/15/99"/)
      real (KIND=r8)                 :: scalar = 999.99_r8
      integer,        dimension(5)   :: ival = (/1,2,3,4,5/)
      real (KIND=r8), dimension(3,4) :: dval =                          &
     &       reshape ( (/11.1_r8, 12.2_r8, 13.3_r8,                     &
     &                   21.1_r8, 22.2_r8, 23.3_r8,                     &
     &                   31.1_r8, 32.2_r8, 33.3_r8,                     &
     &                   41.1_r8, 42.2_r8, 43.3_r8/),                   &
     &                  (/3,4/))
      real (KIND=r4), dimension(3,4,2) :: fval =            &
     &       reshape ( (/                                   &
     &                   1.1, 1.2, 1.3,                     &
     &                   2.1, 2.2, 2.3,                     &
     &                   3.1, 3.2, 3.3,                     &
     &                   4.1, 4.2, 4.3,                     &
     &                   5.1, 5.2, 5.3,                     &
     &                   6.1, 6.2, 6.3,                     &
     &                   7.1, 7.2, 7.3,                     &
     &                   8.1, 8.2, 8.3                      &
     &                  /),                                 &
     &                  (/3,4,2/))
      complex(KIND=r8) :: c16val_0,  c16val_1(5), c16val_2(3,4), c16val_3(3,4,2)
      complex(KIND=r4) :: c8val_0,  c8val_1(5), c8val_2(3,4), c8val_3(3,4,2)

! to define Variables
      integer, dimension(3) :: dimlens = (/1,1,1/)
      integer               :: ncid
      integer               :: ier
!
!----------------------------------------
! Create File
      call ezcdf_open(ncid,'EZtest.cdf','w')
      if (ncid .eq. 0) then
         print *,'Error creating file'
         stop
      end if
! Define Variables
! Title
      dimlens(1)=len(title)
      call cdfDefVar(ncid,'Title',dimlens,'CHAR',ier)  ! ier = optional
      if (ier .ne. 0) stop
! Comment
      dimlens(1)=len(comment(1))
      dimlens(2)=size(comment)
      call cdfDefVar(ncid,'Comment',dimlens,'CHAR')
! Scalar
      dimlens(1)=1
      dimlens(2)=1
      call cdfDefVar(ncid,'Scalar',dimlens,'R8')
! 1D-INT
      dimlens(1)=size(ival)
      call cdfDefVar(ncid,'1D-INT',dimlens,'INT')
! 2D-R8
      dimlens(1)=3
      dimlens(2)=4
      call cdfDefVar(ncid,'2D-R8',dimlens,'R8')
! 3D-R4
      dimlens(1)=3
      dimlens(2)=4
      dimlens(3)=2
      call cdfDefVar(ncid,'3D-R4',dimlens,'R4')
! Complex
      dimlens = 0
      call cdfDefVar(ncid,'Scalar_C16',dimlens,'C16')
      call cdfDefVar(ncid,'Scalar_C8',dimlens,'C8')
      dimlens = (/ 5, 1, 1 /)
      call cdfDefVar(ncid,'1D_C16',dimlens,'C16')
      call cdfDefVar(ncid,'1D_C8',dimlens,'C8')
      dimlens = (/ 3, 4, 1 /)
      call cdfDefVar(ncid,'2D_C16',dimlens,'C16')
      call cdfDefVar(ncid,'2D_C8',dimlens,'C8')
      dimlens = (/ 3, 4, 2 /)
      call cdfDefVar(ncid,'3D_C16',dimlens,'C16')
      call cdfDefVar(ncid,'3D_C8',dimlens,'C8')

      c16val_0 = scalar * (1._r8, 2._r8)
      c16val_1 = ival * (1._r8, 2._r8)
      c16val_2 = dval * (1._r8, 2._r8)
      c16val_3 = fval * (1._r8, 2._r8)
      c8val_0 = scalar * (3._r4, 4._r4)
      c8val_1 = ival * (3._r4, 4._r4)
      c8val_2 = dval * (3._r4, 4._r4)
      c8val_3 = fval * (3._r4, 4._r4)

!
! Write Variables
      call cdfPutVar(ncid,'Title',title,ier)   ! ier = optional
      if (ier .ne. 0) stop
      call cdfPutVar(ncid,'Comment',comment)
      call cdfPutVar(ncid,'Scalar',scalar)
      call cdfPutVar(ncid,'2D-R8',dval)
      call cdfPutVar(ncid,'3D-R4',fval)
      call cdfPutVar(ncid,'1D-INT',ival)
      call cdfPutVar(ncid,'Scalar_C16',c16val_0)
      call cdfPutVar(ncid,'1D_C16',c16val_1)
      call cdfPutVar(ncid,'2D_C16',c16val_2)
      call cdfPutVar(ncid,'3D_C16',c16val_3)
      call cdfPutVar(ncid,'Scalar_C8',C8val_0)
      call cdfPutVar(ncid,'1D_C8',C8val_1)
      call cdfPutVar(ncid,'2D_C8',C8val_2)
      call cdfPutVar(ncid,'3D_C8',C8val_3)

      call ezcdf_close(ncid ,ier)
!
! Now read back
      call SampleRead
      end
!========================================================================
      SUBROUTINE SampleRead
      USE ezcdf
      implicit none
      integer, PARAMETER :: r8 = SELECTED_REAL_KIND(12,100)
      integer, PARAMETER :: r4 = SELECTED_REAL_KIND(6,37)
! netCDF data sets:
      real (KIND=r8), dimension(:,:), allocatable :: dval
      real (KIND=r4), dimension(:,:,:), allocatable :: fval
      integer,        dimension(:),   allocatable :: ival
      real (KIND=r8)                              :: scalar
      character*25                                :: title
      character*8, dimension(2)                   :: comment
      complex(KIND=r8) :: c16val_0
      complex(KIND=r8), allocatable :: c16val_1(:), c16val_2(:,:), &
           & c16val_3(:,:,:)
      complex(KIND=r4) :: c8val_0
      complex(KIND=r4), allocatable :: c8val_1(:), c8val_2(:,:), &
           & c8val_3(:,:,:)
! Local
      integer                 :: ncid
      integer                 :: ier
      integer, dimension(3)   :: dimlens = (/1,1,1/)
      character*4             :: type
!
      integer                 :: i, j, k, ner
!----------------------------------------------------------------

      ner = 0
      print "(/'Read again:'/)"
! Open File
      call ezcdf_open(ncid,'EZtest.cdf','r')
! Inquire dimensions and type of Variable
! 1D-INT
      call cdfInqVar(ncid,'1D-INT',dimlens,type,ier)      ! ier = optional
      ALLOCATE( ival(dimlens(1)), STAT=ier )
! 2D-R8
      call cdfInqVar(ncid,'2D-R8',dimlens,type)
      ALLOCATE( dval(dimlens(1),dimlens(2)), STAT=ier )
! 2D-R4
      call cdfInqVar(ncid,'3D-R4',dimlens,type)
      ALLOCATE( fval(dimlens(1),dimlens(2),dimlens(3)), STAT=ier )
! Comment
      call cdfInqVar(ncid,'Comment',dimlens,type,ier)
! Complex
      call cdfInqVar(ncid,'Scalar_C16',dimlens,type,ier)
      print "('Scalar_C16 dims=',3I4, 2x, a)", dimlens, type
      call cdfInqVar(ncid,'1D_C16',dimlens,type,ier)
      print "('1D_C16 dims=',3I4, 2x, a)", dimlens, type
      allocate(c16val_1(dimlens(1)))
      call cdfInqVar(ncid,'1D_C8',dimlens,type,ier)
      allocate(c8val_1(dimlens(1)))
      call cdfInqVar(ncid,'2D_C16',dimlens,type,ier)
      print "('2D_C16 dims=',3I4, 2x, a)", dimlens, type
      allocate(c16val_2(dimlens(1), dimlens(2)))
      call cdfInqVar(ncid,'2D_C8',dimlens,type,ier)
      allocate(c8val_2(dimlens(1), dimlens(2)))
      call cdfInqVar(ncid,'3D_C16',dimlens,type,ier)
      print "('3D_C16 dims=',3I4, 2x, a)", dimlens, type
      allocate(c16val_3(dimlens(1), dimlens(2), dimlens(3)))
      call cdfInqVar(ncid,'3D_C8',dimlens,type,ier)
      allocate(c8val_3(dimlens(1), dimlens(2), dimlens(3)))
!
! Get data
      call cdfGetVar(ncid,'Title',title,ier)             ! ier = optional
      print *,'Title:   ',title
      call cdfGetVar(ncid,'Comment',comment)
      print *,'Comment: ',comment
      call cdfGetVar(ncid,'Scalar',scalar)
      print *,'Scalar:  ',scalar
      call cdfGetVar(ncid,'1D-INT',ival)
      print *,'1d-INT:  ',ival
      call cdfGetVar(ncid,'2D-R8',dval)
      call cdfGetVar(ncid,'3D-R4',fval)
      call cdfGetVar(ncid,'Scalar_C16',c16val_0)
      call cdfGetVar(ncid,'1D_C16',c16val_1)
      call cdfGetVar(ncid,'2D_C16',c16val_2)
      call cdfGetVar(ncid,'3D_C16',c16val_3)
      call cdfGetVar(ncid,'Scalar_C8',C8val_0)
      call cdfGetVar(ncid,'1D_C8',C8val_1)
      call cdfGetVar(ncid,'2D_C8',C8val_2)
      call cdfGetVar(ncid,'3D_C8',C8val_3)
      print *,'2D-R8:'
      do i=1,4
         print *,dval(:,i)
      end do
      print *,'3D-R4:'
      do i=1,2
         print *,(fval(:,j,i),j=1,4)
      end do

      if(c16val_0 /= scalar * (1._r8, 2._r8)) ner = ner + 1

      do i = 1, size(ival)
         if(c16val_1(i) /= ival(i) * (1._r8, 2._r8)) ner = ner + 1
      enddo

      do j = 1, size(dval, 2)
         do i = 1, size(dval, 1)
            if(c16val_2(i,j) /= dval(i,j) * (1._r8, 2._r8)) ner = ner + 1
         enddo
      enddo

      do k = 1, size(fval, 3)
         do j = 1, size(fval, 2)
            do i = 1, size(fval, 1)
               if(c16val_3(i,j,k) /= fval(i,j,k) * (1._r8, 2._r8)) ner = ner + 1
            enddo
         enddo
      enddo
      print *,ner,' errors occurred while reading in complex data'
      
    end SUBROUTINE SampleRead

! Alexei Pankin, Lehigh Univ., 29 Jul 2002
! interface module to EzCdf module
!   - save output netcdf file
module netcdf_int_mod
  !
  use ezcdf
  !
  !IMPLICIT NONE
  save
  private
  integer                    :: ilen          ! length of record in a scratch file
  real, allocatable          :: ztime(:)      ! times
  type zt2d100vars
    character(23) :: cvar                     ! variable name to be saved
    real, pointer :: pvar(:)                  ! pointer to 1d array
    real          :: mvar                     ! multiplier for the variable
  end type
  type zt2d110vars
    character(23) :: cvar                     ! variable name to be saved
    real, pointer :: pvar(:,:)                ! pointer to 1d array
    real          :: mvar                     ! multiplier for the variable
  end type
  type zt2d111vars
    character(23) :: cvar                     ! variable name to be saved
    real, pointer :: pvar(:,:,:)              ! pointer to 1d array
    real          :: mvar                     ! multiplier for the variable
  end type
  type(zt2d100vars), allocatable:: z2dvars1(:)  
  type(zt2d110vars), allocatable:: z2dvars2(:)  
  type(zt2d111vars), allocatable:: z2dvars3(:)  
  integer                       :: istep         ! =number of records/number of variables
  integer                       :: in2dvars(3)   ! number of 2d vars to be saved
  integer                       :: ilaststep = 0 ! 
  integer                       :: ilastrec  = 0 ! 
  public                        :: netcdf_output
contains
  !.....................................................
  ! Subroutine netcdf_output to open temporary scratch
  ! file, save data, and convert the data into the
  ! netcdf file
  !.....................................................
  subroutine netcdf_output( iluntmp, iluncdf, iaction, ier )
   !
   include 'cparm.m'
   include 'cbaldr.m'
   include 'ctrplot.m'
   !
   integer, intent(in)       :: iluntmp    ! file id of scratch file 
   integer, intent(in)       :: iluncdf    ! file id of netcdf file 
   integer, intent(in)       :: iaction    ! 1 to open temp file, 2 to save data, 
                                           ! 3 to copy data to netcdf file and 
                                           ! close temp file
   integer, intent(out)      :: ier        ! error indicator
   !
   integer                   :: zdm(3)
   integer                   :: i1, i2, in123
   real, allocatable         :: z2dtmp(:,:)
   character*32              :: zcrun
   !
   ier   = 0
   zcrun = trim(adjustl(label1))
   select case (iaction)
   case(1)
      inquire(IOLENGTH=ilen) te(1:nzones+2)
      open(iluntmp, STATUS='SCRATCH', ACCESS='DIRECT', RECL=ilen, IOSTAT=ier)
      in2dvars  = (/2, 2, 0/) 
      in123     = sum(in2dvars)
      !in2dvars = (/2, 2, 2/) 
      allocate(z2dvars1(in2dvars(1)), z2dvars2(in2dvars(2)), z2dvars3(in2dvars(3)), ztime(ntime))
      z2dvars1(1) = zt2d100vars('TE',      tesms(1:nzones+2),   useh)
      z2dvars1(2) = zt2d100vars('TI',      tisms(1:nzones+2),   useh)
      z2dvars2(1) = zt2d110vars('NE',      rhoels,              1.)
      z2dvars2(2) = zt2d110vars('NI',      rhoins,              1.)
      !z2dvars3(1) = zt2d111vars('NH',      rhohs,               1.)
      !z2dvars3(2) = zt2d111vars('NX',      rhois,               1.)
   case(2)
      if (ilaststep>0) then
        if (atbi*uist < ztime(ilaststep)+splot) return
      endif 
      ilaststep = ilaststep + 1
      ztime(ilaststep) = atbi * uist
      do i1 = 1, in2dvars(1)
        ilastrec = ilastrec + 1
        write(iluntmp, REC=ilastrec) z2dvars1(i1)%pvar*z2dvars1(i1)%mvar
      enddo
      !
      do i1 = 1, in2dvars(2)
        ilastrec = ilastrec + 1
        write(iluntmp, REC=ilastrec) z2dvars2(i1)%pvar(1:1, 1:nzones+2)*z2dvars2(i1)%mvar
      enddo
      !
      !do i1 = 1, in2dvars(3)
      !  ilastrec = ilastrec + 1
      !  write(iluntmp, REC=ilastrec) z2dvars3(i1)%pvar(1:1, 1:nzones+2)*z2dvars3(i1)%mvar
      !enddo
   case(3)
      call cdfopn(iluncdf, zcrun(1:scan(zcrun, ' ')-1)//'.nc', 'w')
      zdm = (/ilaststep, 0, 0/)
      call cdfDefVar(iluncdf, 'TIME',   zdm,    'R8')
      zdm = (/nzones+2, 0, 0/)
      call cdfDefVar(iluncdf, 'XBOUNI', zdm,    'R8')
      zdm = (/ilaststep, nzones+2, 0/)
      do i1 = 1, in2dvars(1)
        call cdfDefVar(iluncdf, z2dvars1(i1)%cvar, zdm, 'R8')
      enddo
      do i1 = 1, in2dvars(2)
        call cdfDefVar(iluncdf, z2dvars2(i1)%cvar, zdm, 'R8')
      enddo
      call cdfPutVar(iluncdf, 'TIME',     ztime, ier)
      call cdfPutVar(iluncdf, 'XBOUNI',   xbouni, ier)
      !
      allocate(z2dtmp(ilaststep, nzones+2))
      do i1 = 1, in2dvars(1)
        do i2 = 1, ilaststep
          read(iluntmp, REC=(i2-1)*in123+i1) z2dtmp(i2, 1:nzones+2)
        end do
        call cdfPutVar(iluncdf, z2dvars1(i1)%cvar, z2dtmp, ier)
      end do
      !
      do i1 = in2dvars(1)+1, in123
        do i2 = 1, ilaststep
          read(iluntmp, REC=(i2-1)*in123+i1) z2dtmp(i2, 1:nzones+2)
        end do
        call cdfPutVar(iluncdf, z2dvars2(i1-in2dvars(1))%cvar, z2dtmp, ier)
      end do	
      call cdfCls(iluncdf)       
      deallocate(z2dvars1, z2dvars2, z2dvars3, ztime, z2dtmp)
      close(iluntmp)
   case default
      call bad_exit()
   end select
  end subroutine netcdf_output
  !
end module netcdf_int_mod

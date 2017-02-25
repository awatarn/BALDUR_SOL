! Alexei Pankin, Lehigh Univ., 30 Jul 2002
! Module to read inout netcdf profile and to interpolate
! to the radial code grid
module ncprofile_mod
  ! implicit none
  public
  save
  !
  type ncprofile
    real, pointer :: x(:,:)
    real, pointer :: t(:)
  end type ncprofile
  !
  type(ncprofile)            :: qprofile   ! q-profile
  type(ncprofile)            :: cdprofile  ! current drive profile
  ! private variables
  real, allocatable, private :: zrho(:)
  real, allocatable, private :: zvalue(:,:)
  integer, private           :: idm1, idm2
  !
  real, allocatable, target, private :: zqt(:), zqx(:,:) 
  real, allocatable, target, private :: zcdt(:), zcdx(:,:) 
  !
  private :: ncprofile_interpolate
  !
contains
  !-------------------------------------------------
  ! read netcdf file and save interpolated profile
  ! in module public section
  !-------------------------------------------------
  subroutine ncprofile_read(zvarname, zfilename, ierr)
   !
   use ezcdf
   include 'cparm.m'
   include 'cbaldr.m'
   !
   character*(*), intent(in) :: zvarname  ! name of the variable
   character*(*), intent(in) :: zfilename ! name of the netcdf
   integer, intent(out)      :: ierr      ! output error indicator
   character*12              :: zvartype
   ! local variable
   integer                   :: ilun = 25      ! file unit number
   integer                   :: idm(3)
   integer, parameter        :: icubspl = 0, & ! interpolate with cubic spline
                                ihermit = 1    ! interpolate with quasi Hermite
   integer, parameter        :: ieven   = 1, & ! even function
                                iodd    = -1,& ! odd function
                                inosym  = 0    ! no symmetry
   ! Set default values
   ! ---------------------------
   ierr = 0
   ! Open netcdf file
   ! ---------------------------
   call ezcdf_open(ilun, zfilename,'r', ierr)
   if (ierr/=0) call errmsg_exit('ncprofile_read: cannot open netcdf file '// zfilename);
   ! Retrieving of array dimensions
   ! ---------------------------
   call cdfInqVar(ilun, 'XI',   idm, zvartype)
   idm2 = idm(1);  allocate(zrho(idm2))
   call cdfGetVar(ilun, 'XI',   zrho, ierr)
   if (ierr/=0) call errmsg_exit('ncprofile_read: cannot read from netcdf file '// zfilename); 
   call cdfInqVar(ilun, 'TIME', idm, zvartype)
   idm1 = idm(1);  allocate(zvalue(idm1, idm2))
   select case(zvarname)
     case ('Q')
       allocate(zqt(idm1)) 
       qprofile%t=>zqt(1:idm1)
       call cdfGetVar(ilun, 'TIME', qprofile%t, ierr)
       if (ierr/=0) call errmsg_exit('ncprofile_read: cannot read from netcdf file '// zfilename); 
       call cdfInqVar(ilun, 'Q',   idm, zvartype)
       call cdfGetVar(ilun, 'Q',   zvalue, ierr)
       if (ierr/=0) call errmsg_exit('ncprofile_read: cannot read from netcdf file '// zfilename); 
       allocate(zqx(idm1,mzones+1))
       qprofile%x=>zqx(1:idm1,1:mzones+1)
       call ncprofile_interpolate(ihermit, ieven, 2, qprofile%x)
     case ('CDRIVE')
       allocate(zcdt(idm1))
       cdprofile%t=>zcdt(1:idm1)
       call cdfGetVar(ilun, 'TIME', cdprofile%t, ierr)
       if (ierr/=0) call errmsg_exit('ncprofile_read: cannot read from netcdf file '// zfilename); 
       call cdfInqVar(ilun, 'CDRIVE',   idm, zvartype)
       call cdfGetVar(ilun, 'CDRIVE',   zvalue, ierr)
       if (ierr/=0) call errmsg_exit('ncprofile_read: cannot read from netcdf file '// zfilename); 
       allocate(zcdx(idm1,mzones+1))
       cdprofile%x=>zcdx(1:idm1,1:mzones+1)
       call ncprofile_interpolate(ihermit, ieven, 1, cdprofile%x)
     case default 
       call errmsg_exit('ncerror_read: variable ' // zvarname // ' is not implemented yet.')
   end select
   deallocate(zrho, zvalue)
   ! Close netcdf file
   ! ---------------------------
   call ezcdf_close(ilun, ierr)
  end subroutine ncprofile_read
  !
  !-------------------------------------------------
  ! interpolate profile to BALDUR code radial grid
  !-------------------------------------------------
  subroutine ncprofile_interpolate(iinterp, isym, izone, zresult)
   !
   include 'cparm.m'
   include 'cbaldr.m'
   !
   integer, intent(in) :: iinterp    ! interpolation parameter
   integer, intent(in) :: isym       ! type of the function symmetry
   integer, intent(in) :: izone      ! =1 for interpolation to zone centers
                                     ! =2 for interpolation to zone boundaries
   real, intent(out)   :: zresult(idm1, *) ! interpolated values
   ! local variable
   real, allocatable   :: zc(:)      ! coeffcients for interp.
   integer             :: it
   !
   allocate(zc(5*idm2))
   do it = 1, idm1
     select case(izone)
       case(1) ! interpolation to zone centers
         call cubint(zrho, zvalue(it:it, 1:idm2), idm2, iinterp, zc, 5*idm2, &
              xzoni, zresult(it:it, 1:mzones-1), mzones-1, 0, 0., &
              isym, 'interpolation of the netcdf profile to the radial grid.')
       case(2) ! interpolation to zone boundaries
         call cubint(zrho, zvalue(it:it, 1:idm2), idm2, iinterp, zc, 5*idm2, &
              xbouni, zresult(it:it, 1:mzones+1), mzones+1, 0, 0., &
              isym, 'interpolation of the netcdf profile to the radial grid.')
       case default
         call errmsg_exit('ncprofile_interpolcate: error izones should be 1 or 2.')
       end select
   enddo
   deallocate(zc)
  end subroutine ncprofile_interpolate
  !
  !-------------------------------------------------
  ! time interpolation 
  !-------------------------------------------------
  subroutine ncprofile_tinterpolate(ztime, zncp, zx)
    real, intent(in)            :: ztime  ! time for interpolation
    type(ncprofile), intent(in) :: zncp   ! ncprofile structure 
    real, intent(out)           :: zx(*)  ! interploated value
    ! local variables
    integer                     :: is, it, j
    !
    if (associated(zncp%t)) then 
      it = size(zncp%t)
      is = size(zncp%x)/it
      do j = 1, is
        call timint(ztime, zx(j), zncp%t, it, transpose(zncp%x), j, is)
      enddo
    endif
  end subroutine ncprofile_tinterpolate
end module ncprofile_mod

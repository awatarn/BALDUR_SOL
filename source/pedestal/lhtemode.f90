subroutine lhtemode(zbt, zr, zhmass, zslne, zzeff, zte, lbound, iedge,  mode)
  IMPLICIT NONE
  !
  include 'cparm.m'
  save imode
  !
  ! Te_cr=0.45*(zbt*zbt*zzeff(ij)/sqrt(zr*za))**(1./3.)*sqrt(zslne(ij))
  !
  ! input variables
  !
  real, intent(in)       :: zbt            ! vacuum Bt at geometric center
  real, intent(in)       :: zr             ! geometric center
  real, intent(in)       :: zhmass         ! hydrogenic isotope mass AMU at flux surface
  real, intent(in)       :: zslne(1:iedge) ! electron density
  real, intent(in)       :: zzeff(1:iedge) ! Zeff local to flux surface at peak value of -grad(ne)/ne
  real, intent(in)       :: zte(1:iedge)   ! electron temperature
  integer, intent(in)    :: iedge          ! array dimension
  integer, intent(in)    :: lbound(*)      ! paramters: 
                                           !   lbound(7) - number of grid points from the edge
					   !   lbound(8) - number of consequent time steps to set mode=1
  !
  ! output variables
  !
  integer, intent(inout) :: mode     ! = 0 for L mode
                                     ! = 1 for H mode
  !
  ! local variables
  !
  integer                :: ij(1)    ! index
  integer                :: i, i1    ! index
  real                   :: te_crit  ! critical temperature
  real                   :: zdslne1, zdslne2
  real, allocatable      :: zslnemin(:)
  integer, allocatable   :: imin(:)
  integer                :: imode = 0  ! counter
  !
  ! if (mode == 1) return  ! no H->L transition 
  mode = 0
  !
  ! Finding multiple minimums
  !
  allocate(zslnemin(10), imin(10))
  i1 = 1
  imin(i1) = iedge-1
  zslnemin(i1) = abs(zslne(iedge-1))
  zdslne1 = abs(zslne(iedge-2)) - abs(zslne(iedge-1))
  do i = iedge-1,iedge-lbound(7),-1
    zdslne2 = abs(zslne(i-1))-abs(zslne(i))
    if ((zdslne1<0) .and. (zdslne2>0)) then
      i1 = i1+1
      imin(i1) = i
      zslnemin(i1) = abs(zslne(i))
    endif
    zdslne1 = zdslne2
  enddo
  deallocate(zslnemin, imin)
  !
  ij  = minloc(abs(zslne(iedge-lbound(7):iedge-1)))
  i   = ij(1) + iedge - lbound(7) - 1
  !
  te_crit = 0.045*(zbt*zbt*zzeff(i)/sqrt(zr*zhmass))**(1./3.)* &
            sqrt(abs(zslne(i)))
  !
  if (te_crit < zte(i)) then
    imode = imode + 1
    if (imode>lbound(8)) mode = 1     ! H mode
  else
    mode  = 0     ! L mode
    imode = 0     ! reset the counter
  endif
  !
  !write(6,'(1A, i2, 1A, e13.5, 1A, i2, 1A, e13.5, 1A)') &
  !        'Mode = ', mode, '  (t_cr=', te_crit, ' < te(', i, ')=', zte(i), ')'
  !
end subroutine lhtemode

subroutine ezcdf_close(ncid, ier)
  include "netcdf.inc"
  INTEGER, INTENT(in) ::  ncid
  integer, optional,         intent(out) :: ier
  INTEGER status
  status = nf_close(ncid)
  call handle_err(status,' ','cdfcls','nf_close')
  if(PRESENT (ier)) ier = status
end subroutine ezcdf_close

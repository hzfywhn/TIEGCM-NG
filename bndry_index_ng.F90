subroutine bndry_index_ng(domain,nlon,nlat,lon0,lon1,lat0,lat1)
! find actual boundaries in subdomains

  implicit none

  integer,dimension(4),intent(in) :: domain
  integer,intent(in) :: nlon,nlat
  integer,intent(out) :: lon0,lon1,lat0,lat1

  if (domain(1) == 1) then
    lon0 = -1
  else
    lon0 = domain(1)
  endif

  if (domain(2) == nlon) then
    lon1 = nlon+2
  else
    lon1 = domain(2)
  endif

  if (domain(3) == 1) then
    lat0 = -1
  else
    lat0 = domain(3)
  endif

  if (domain(4) == nlat) then
    lat1 = nlat+2
  else
    lat1 = domain(4)
  endif

end subroutine bndry_index_ng

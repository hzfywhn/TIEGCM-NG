module nudge_ng_module
! force the model fields near lower boundary with external fields

  use params_module, only: mx_ng
  implicit none

! if external fields contain lb, then the model lower boundary will use external fields
  character(len=*), dimension(*), parameter :: lb = (/'TN', 'UN', 'VN', 'Z '/)
  integer, parameter :: nlb = size(lb)

  logical :: wrap ! whether the external field has a full longitude cycle
  integer :: nlon, nlat, nlev, ntime, nfile, nf4d, ifile, itime, ncid
  real :: y0, y1, z0, z1
  integer, dimension(nlb) :: lb_idx
  logical, dimension(:), allocatable :: no_fill
  integer, dimension(:), allocatable :: time, nt, t0, t1, f4d_idx, varid
  real, dimension(:), allocatable :: lon, fill_value
! note that lb_idx and f4d_idx are referenced upon different fields
! lb_idx is the index of lb fields in the external fields
! f4d_idx is the index of external fields in the model fields

  type fields
    integer :: maxlev, latbeg, latend, lonbeg, lonend, offbeg, offend, lonbeg1, lonend1
    real, dimension(:,:), allocatable :: vert_weight, hori_weight ! vertical and horizontal relaxation factors
    real, dimension(:,:,:,:), allocatable :: lbc
    real, dimension(:,:,:,:,:), allocatable :: f4d
  end type fields

  type(fields), dimension(0: mx_ng) :: flds

  contains
!-----------------------------------------------------------------------
  subroutine check
! check the compliance of the external field

    use params_module, only: mxhvols, zibot
    use input_module, only: nudge_ncpre, nudge_ncfile, nudge_ncpost, nudge_flds, nudge_sponge, nudge_tshift
    use fields_module, only: f4d
    use char_module, only: find_index
    use netcdf, only: nf90_open, nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, &
      nf90_get_var, nf90_inq_var_fill, nf90_close, nf90_nowrite, nf90_nofill

    integer, parameter :: maxnt = 1440
    real, parameter :: eps = 1e-12
    integer :: stat, dimid_lon, dimid_lat, dimid_lev, dimid_time, &
      varid_lon, varid_lat, varid_lev, varid_time, ifld, nofill
    real :: dlon
    real(kind=4) :: fillvalue
    integer, dimension(mxhvols, maxnt) :: t
    real, dimension(1) :: tmp
    external :: shutdown

    stat = nf90_open(trim(nudge_ncpre)//trim(nudge_ncfile(1))//trim(nudge_ncpost), nf90_nowrite, ncid)

    stat = nf90_inq_dimid(ncid, 'lon', dimid_lon)
    stat = nf90_inquire_dimension(ncid, dimid_lon, len=nlon)
    allocate(lon(nlon))
    stat = nf90_inq_varid(ncid, 'lon', varid_lon)
    stat = nf90_get_var(ncid, varid_lon, lon)

! if the external field covers the full longitude cycle,
! then model fields at lon=1,2 and lon=nlonp4-1,nlonp4 will be obtained from the external field
    dlon = (lon(nlon)-lon(1)) / (nlon-1)
    if (abs(lon(1)+360-lon(nlon)-dlon) < eps*dlon) then
      wrap = .true.
    else
      wrap = .false.
    endif

    stat = nf90_inq_dimid(ncid, 'lat', dimid_lat)
    stat = nf90_inquire_dimension(ncid, dimid_lat, len=nlat)
    stat = nf90_inq_varid(ncid, 'lat', varid_lat)
    stat = nf90_get_var(ncid, varid_lat, tmp, count=(/1/))
    y0 = tmp(1)
    stat = nf90_get_var(ncid, varid_lat, tmp, start=(/nlat/), count=(/1/))
    y1 = tmp(1)

    stat = nf90_inq_dimid(ncid, 'lev', dimid_lev)
    stat = nf90_inquire_dimension(ncid, dimid_lev, len=nlev)
    stat = nf90_inq_varid(ncid, 'lev', varid_lev)
    stat = nf90_get_var(ncid, varid_lev, tmp, count=(/1/))
    z0 = tmp(1)
    stat = nf90_get_var(ncid, varid_lev, tmp, start=(/nlev/), count=(/1/))
    z1 = tmp(1)

    nf4d = count(len_trim(nudge_flds) > 0)
    allocate(f4d_idx(nf4d))
    allocate(varid(nf4d))
    allocate(no_fill(nf4d))
    allocate(fill_value(nf4d))
    do ifld = 1, nf4d
      stat = nf90_inq_varid(ncid, trim(nudge_flds(ifld)), varid(ifld))
      stat = nf90_inq_var_fill(ncid, varid(ifld), nofill, fillvalue)
      if (nofill == nf90_nofill) then
        no_fill(ifld) = .true.
      else
        no_fill(ifld) = .false.
      endif
      fill_value(ifld) = fillvalue
      f4d_idx(ifld) = find_index(nudge_flds(ifld), f4d%short_name)
      if (f4d_idx(ifld) == 0) call shutdown(trim(nudge_flds(ifld))//' is not a valid model field')
    enddo

    stat = nf90_close(ncid)

    if (z0>zibot .or. z1<zibot) call shutdown('nudge_ncfile dimension lev must include model lbc')
    if (nudge_sponge(1)*2>=max(lon(nlon)-lon(1), y1-y0) .or. nudge_sponge(2)>z1-z0) &
      call shutdown('nudge_sponge cannot exceed nudge_ncfile dimension range')

! if there are multiple external data files, the model will read the correct data file based on time
    nfile = count(len_trim(nudge_ncfile) > 0)
    allocate(nt(nfile))
    do ifile = 1, nfile
      stat = nf90_open(trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost), nf90_nowrite, ncid)
      stat = nf90_inq_dimid(ncid, 'time', dimid_time)
      stat = nf90_inquire_dimension(ncid, dimid_time, len=nt(ifile))
      stat = nf90_inq_varid(ncid, 'time', varid_time)
      stat = nf90_get_var(ncid, varid_time, t(ifile, 1: nt(ifile)))
      stat = nf90_close(ncid)
    enddo

    allocate(t0(nfile))
    allocate(t1(nfile))
    t0(1) = 1
    t1(1) = t0(1) + nt(1) - 1
    do ifile = 2, nfile
      t0(ifile) = t1(ifile-1) + 1
      t1(ifile) = t0(ifile) + nt(ifile) - 1
    enddo

    ntime = t1(nfile)
    if (ntime <= 1) call shutdown('nudge_ncfile must have at least 2 time points')

    allocate(time(ntime))
    do ifile = 1, nfile
      time(t0(ifile): t1(ifile)) = t(ifile, 1: nt(ifile))
    enddo
    time = time + nudge_tshift

    do ifld = 1, nlb
      lb_idx(ifld) = find_index(lb(ifld), nudge_flds)
    enddo

  end subroutine check
!-----------------------------------------------------------------------
  subroutine init
! synchronize external field paramters among processes and init subdomain grids

    use params_module, only: n_ng, nlevp1, nlevp1_ng, zmbot, zibot, &
      zpmid, zpint, glon0, glat, zpmid_ng, zpint_ng, glon_ng, glat_ng, ispval
    use cons_module, only: dtr
    use input_module, only: nudge_ncpre, nudge_ncfile, nudge_ncpost, &
      nudge_flds, nudge_sponge, nudge_delta, nudge_power
    use fields_module, only: f4d
    use fields_ng_module, only: ng_flds=>flds
    use mpi_module, only: mytid, lon0, lon1, lat0, lat1
    use mpi_f08, only: mpi_request, mpi_bcast, mpi_ibcast, mpi_waitall, &
      mpi_integer, mpi_logical, mpi_real8, mpi_comm_world, mpi_statuses_ignore
    use netcdf, only: nf90_open, nf90_inq_varid, nf90_nowrite

    logical :: intersect
    integer :: stat, if4d, i_ng, nk, k, j0, j1, lat, i0, i1, i, shift, ls, rs
    real :: lb, rb, tb, bb, slon, latpart1, latpart2, lonpart
    integer, dimension(1) :: idx
    real, dimension(max(nlevp1, maxval(nlevp1_ng))) :: zpm, zpi, zout
    type(mpi_request), dimension(3) :: request
    logical, dimension(:), allocatable :: buffer_l
    integer, dimension(:), allocatable :: buffer_i
    real, dimension(:), allocatable :: buffer_r
    real, dimension(min(lat0, minval(ng_flds%latd0)): max(lat1, maxval(ng_flds%latd1))) :: model_lat, yc
    real, dimension(min(lon0, minval(ng_flds%lond0)): max(lon1, maxval(ng_flds%lond1))) :: model_lon, xc
    real, dimension(min(lon0, minval(ng_flds%lond0)): max(lon1, maxval(ng_flds%lond1)), &
      min(lat0, minval(ng_flds%latd0)): max(lat1, maxval(ng_flds%latd1))) :: dist

    allocate(buffer_i(4))

    if (mytid == 0) buffer_i = (/nlon, ntime, nfile, nf4d/)

    call mpi_bcast(buffer_i, 4, mpi_integer, 0, mpi_comm_world)

    if (mytid /= 0) then
      nlon = buffer_i(1)
      ntime = buffer_i(2)
      nfile = buffer_i(3)
      nf4d = buffer_i(4)
    endif

    deallocate(buffer_i)
    allocate(buffer_l(nf4d+1))
    allocate(buffer_i(nlb+ntime+nfile*3+nf4d+2))
    allocate(buffer_r(nlon+nf4d+4))

    if (mytid == 0) then
      buffer_l(1: nf4d) = no_fill
      buffer_l(nf4d+1) = wrap

      buffer_i(1: nlb) = lb_idx
      buffer_i(nlb+1: nlb+ntime) = time
      buffer_i(nlb+ntime+1: nlb+ntime+nfile) = nt
      buffer_i(nlb+ntime+nfile+1: nlb+ntime+nfile*2) = t0
      buffer_i(nlb+ntime+nfile*2+1: nlb+ntime+nfile*3) = t1
      buffer_i(nlb+ntime+nfile*3+1: nlb+ntime+nfile*3+nf4d) = f4d_idx
      buffer_i(nlb+ntime+nfile*3+nf4d+1: nlb+ntime+nfile*3+nf4d+2) = (/nlat, nlev/)

      buffer_r(1: nlon) = lon
      buffer_r(nlon+1: nlon+nf4d) = fill_value
      buffer_r(nlon+nf4d+1: nlon+nf4d+4) = (/y0, y1, z0, z1/)
    endif

    call mpi_ibcast(buffer_l, nf4d+1, mpi_logical, 0, mpi_comm_world, request(1))
    call mpi_ibcast(buffer_i, nlb+ntime+nfile*3+nf4d+2, mpi_integer, 0, mpi_comm_world, request(2))
    call mpi_ibcast(buffer_r, nlon+nf4d+4, mpi_real8, 0, mpi_comm_world, request(3))
    call mpi_waitall(3, request, mpi_statuses_ignore)

    if (mytid /= 0) then
      allocate(no_fill(nf4d))
      allocate(time(ntime))
      allocate(nt(nfile))
      allocate(t0(nfile))
      allocate(t1(nfile))
      allocate(f4d_idx(nf4d))
      allocate(lon(nlon))
      allocate(fill_value(nf4d))
      allocate(varid(nf4d))

      no_fill = buffer_l(1: nf4d)
      wrap = buffer_l(nf4d+1)

      lb_idx = buffer_i(1: nlb)
      time = buffer_i(nlb+1: nlb+ntime)
      nt = buffer_i(nlb+ntime+1: nlb+ntime+nfile)
      t0 = buffer_i(nlb+ntime+nfile+1: nlb+ntime+nfile*2)
      t1 = buffer_i(nlb+ntime+nfile*2+1: nlb+ntime+nfile*3)
      f4d_idx = buffer_i(nlb+ntime+nfile*3+1: nlb+ntime+nfile*3+nf4d)
      nlat = buffer_i(nlb+ntime+nfile*3+nf4d+1)
      nlev = buffer_i(nlb+ntime+nfile*3+nf4d+2)

      lon = buffer_r(1: nlon)
      fill_value = buffer_r(nlon+1: nlon+nf4d)
      y0 = buffer_r(nlon+nf4d+1)
      y1 = buffer_r(nlon+nf4d+2)
      z0 = buffer_r(nlon+nf4d+3)
      z1 = buffer_r(nlon+nf4d+4)
    endif

! a sponge layer is imposed to allow smooth transition
    lb = lon(1) + nudge_sponge(1)
    rb = lon(nlon) - nudge_sponge(1)
    tb = y0 + nudge_sponge(1)
    bb = y1 - nudge_sponge(1)

    flds%lonbeg = ispval
    flds%lonend = ispval
    flds%offbeg = ispval
    flds%offend = ispval
    flds%lonbeg1 = ispval
    flds%lonend1 = ispval

    do i_ng = 0, n_ng
      if (i_ng == 0) then
        nk = nlevp1
        zpm(1: nk) = zpmid(1: nk)
        zpi(1: nk) = zpint(1: nk)
        j0 = lat0
        j1 = lat1
        model_lat(j0: j1) = glat(j0: j1)
        i0 = lon0
        i1 = lon1
        model_lon(i0: i1) = glon0(i0: i1)
      else
        nk = nlevp1_ng(i_ng)
        zpm(1: nk) = zpmid_ng(i_ng, 1: nk)
        zpi(1: nk) = zpint_ng(i_ng, 1: nk)
        j0 = ng_flds(i_ng)%latd0
        j1 = ng_flds(i_ng)%latd1
        model_lat(j0: j1) = glat_ng(i_ng, j0: j1)
        i0 = ng_flds(i_ng)%lond0
        i1 = ng_flds(i_ng)%lond1
        model_lon(i0: i1) = glon_ng(i_ng, i0: i1)
      endif

! smooth vertical transition from external fields to model fields, exponential decay
      do k = 1, nk
        if (zpm(k)>zmbot+nudge_sponge(2) .or. zpi(k)>zibot+nudge_sponge(2)) exit
      enddo
      flds(i_ng)%maxlev = k - 1
      allocate(flds(i_ng)%vert_weight(flds(i_ng)%maxlev, nf4d))
      do if4d = 1, nf4d
        if (trim(f4d(f4d_idx(if4d))%vcoord) == 'midpoints') then
          zout(1: flds(i_ng)%maxlev) = zpm(1: flds(i_ng)%maxlev)
        else
          zout(1: flds(i_ng)%maxlev) = zpi(1: flds(i_ng)%maxlev)
        endif
        flds(i_ng)%vert_weight(:, if4d) = exp(-((zout(1: flds(i_ng)%maxlev) - zibot) / nudge_delta(2))**nudge_power(2))
      enddo

      do lat = j0, j1
        if (model_lat(lat) >= y0) exit
      enddo
      flds(i_ng)%latbeg = lat

      do lat = j1, j0, -1
        if (model_lat(lat) <= y1) exit
      enddo
      flds(i_ng)%latend = lat

      if (wrap) then
! if the external field covers the full longitude cycle, then the model subdomain is fully embedded
! move the edge of the external domain to cover the model subdomain
        slon = model_lon((i0 + i1) / 2) - 180
        do shift = -1, 1
          if (lon(1)+shift*360<=slon .and. slon<=lon(nlon)+shift*360) exit
        enddo
        idx = minloc(abs(lon + shift*360 - slon))
        flds(i_ng)%lonbeg = idx(1)
        flds(i_ng)%offbeg = shift
      else

! find the suitable 360 degree wrap of the external domain to intersect the model subdomain
! the following part finds the left edge of the model subdomain
        intersect = .false.
        left: do i = i0, i1
! the shift loop can be omitted if the external domain is also from -180 to 180 (shift=0)
! if the external domain is not from -180 to 180, the 360 degree wrap will be only one of -1,0,1
          do shift = -1, 1
            if ((lon(1)+shift*360<=model_lon(i) .and. model_lon(i)<=lon(nlon)+shift*360)) then
              intersect = .true.
              exit left
            endif
          enddo
        enddo left
        if (intersect) then
          flds(i_ng)%lonbeg = i
          flds(i_ng)%offbeg = shift
        endif

! similar for the right edge of the mode subdomain
        intersect = .false.
        right: do i = i1, i0, -1
          do shift = -1, 1
            if (lon(1)+shift*360<=model_lon(i) .and. model_lon(i)<=lon(nlon)+shift*360) then
              intersect = .true.
              exit right
            endif
          enddo
        enddo right
        if (intersect) then
          flds(i_ng)%lonend = i
          flds(i_ng)%offend = shift
        endif

! after the 360 degree wrap of the external domain, there will be three possibilities, discussed as follows:
! 1. lonbeg<lonend && offbeg==offend, left edge and right edge doesn't cross the 180 degree boundary,
! the overlapping region is (model_lon(lonbeg), model_lon(lonend)) and (lon(1)+offbeg*360, lon(nlon)+offend*360)
! 2. lonbeg<lonend && offbeg+1==offend, left edge doesn't cross the 180 degree boundary but right edge does,
! the overlapping region is splitted into two:
!   - 1. (model_lon(lonbeg), model_lon(lonend1)) vs (lon(1)+offbeg*360, lon(nlon)+offbeg*360)
!   - 2. (model_lon(lonbeg1), model_lon(lonend)) vs (lon(1)+offend*360, lon(nlon)+offend*360)
! the additional parameters lonbeg1,lonend1 of this situation are calculated below
! 3. lonbeg>lonend && offbeg==offend+1, left edge crosses the 180 degree boundary but right edge doesn't,
! left edge and right edge are swapped, the overlapping region also consists of two parts
!   - 1. (model_lon(lonbeg), model_lon(i1)) vs (lon(1)+offbeg*360, lon(nlon)+offbeg*360)
!   - 2. (model_lon(i0), model_lon(lonend)) vs (lon(1)+offend*360, lon(nlon)+offend*360)

        if (flds(i_ng)%lonbeg<flds(i_ng)%lonend .and. flds(i_ng)%offbeg+1==flds(i_ng)%offend) then
          do i = flds(i_ng)%lonbeg, flds(i_ng)%lonend
            if (model_lon(i) > lon(nlon)+flds(i_ng)%offbeg*360) exit
          enddo
          flds(i_ng)%lonend1 = i - 1

          do i = flds(i_ng)%lonend, flds(i_ng)%lonbeg, -1
            if (model_lon(i) < lon(1)+flds(i_ng)%offend*360) exit
          enddo
          flds(i_ng)%lonbeg1 = i + 1
        endif
      endif

      if (flds(i_ng)%latbeg < flds(i_ng)%latend) then
        allocate(flds(i_ng)%hori_weight(i0: i1, flds(i_ng)%latbeg: flds(i_ng)%latend))
        allocate(flds(i_ng)%lbc(i0: i1, flds(i_ng)%latbeg: flds(i_ng)%latend, 2, nlb))
        allocate(flds(i_ng)%f4d(flds(i_ng)%maxlev, i0: i1, flds(i_ng)%latbeg: flds(i_ng)%latend, 2, nf4d))

        do lat = flds(i_ng)%latbeg, flds(i_ng)%latend
          if (model_lat(lat) < tb) yc(lat) = tb
          if (tb<=model_lat(lat) .and. model_lat(lat)<=bb) yc(lat) = model_lat(lat)
          if (model_lat(lat) > bb) yc(lat) = bb
        enddo

        if (wrap) then
          xc(i0: i1) = model_lon(i0: i1)
        else
          do i = i0, i1
            do ls = -1, 1
              if (abs(lb+ls*360 - model_lon(i)) <= 180) exit
            enddo
            do rs = -1, 1
              if (abs(rb+rs*360 - model_lon(i)) <= 180) exit
            enddo
            if (model_lon(i) < lb+ls*360) xc(i) = lb + ls*360
            if (lb+ls*360<=model_lon(i) .and. model_lon(i)<=rb+rs*360) xc(i) = model_lon(i)
            if (model_lon(i) > rb+rs*360) xc(i) = rb + rs*360
          enddo
        endif

! if the point is inside the inner square, the distance is zero
! if the point is outside the inner square, the distance is the great-circle distance to the nearest point
        do lat = flds(i_ng)%latbeg, flds(i_ng)%latend
          latpart1 = sin((yc(lat) - model_lat(lat)) * dtr / 2)**2
          latpart2 = 1 - latpart1 - sin((yc(lat) + model_lat(lat)) * dtr / 2)**2
          do i = i0, i1
            lonpart = sin((xc(i) - model_lon(i)) * dtr / 2)**2
            dist(i, lat) = 2 * asin(sqrt(latpart1 + latpart2*lonpart))
          enddo
        enddo

! the relaxation function is an exponential function of the great-circle distance
        flds(i_ng)%hori_weight = exp(-(dist(i0: i1, flds(i_ng)%latbeg: flds(i_ng)%latend) / (nudge_delta(1)*dtr))**nudge_power(1))
      endif
    enddo

    ifile = 1
    itime = 1

    stat = nf90_open(trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost), nf90_nowrite, ncid)
    do if4d = 1, nf4d
      stat = nf90_inq_varid(ncid, trim(nudge_flds(if4d)), varid(if4d))
    enddo

    call read_data(itime, 1)
    call read_data(itime+1, 2)

    deallocate(buffer_l)
    deallocate(buffer_i)
    deallocate(buffer_r)

  end subroutine init
!-----------------------------------------------------------------------
  subroutine update(modelsec)
! move the cursor itime to include modelsec inside the interval [time(itime), time(itime+1)]

    use params_module, only: n_ng

    integer, intent(in) :: modelsec

    integer :: increment, i_ng

    increment = 0
    do while (itime <= ntime-2)
      if (time(itime)<=modelsec .and. modelsec<=time(itime+1)) exit
      itime = itime + 1
      increment = increment + 1
    enddo

    if (increment >= 1) then
      if (increment == 1) then
! if itime is moved forward by 1, the previous time is then filled the the current time
        do i_ng = 0, n_ng
          if (flds(i_ng)%latbeg < flds(i_ng)%latend) then
            flds(i_ng)%lbc(:, :, 1, :) = flds(i_ng)%lbc(:, :, 2, :)
            flds(i_ng)%f4d(:, :, :, 1, :) = flds(i_ng)%f4d(:, :, :, 2, :)
          endif
        enddo
      else
! only at the first step will both steps be re-calculated
        call read_data(itime, 1)
      endif
      call read_data(itime+1, 2)
    endif

  end subroutine update
!-----------------------------------------------------------------------
  subroutine read_data(itime, it)

    use params_module, only: n_ng, zibot, zpmid, zpint, glon0, glat, glon_ng, glat_ng
    use input_module, only: nudge_ncpre, nudge_ncfile, nudge_ncpost, &
      nudge_lbc, nudge_f4d, nudge_pert, nudge_flds, nudge_level
    use fields_module, only: f4d
    use fields_ng_module, only: ng_flds=>flds, zlog
    use interp_module, only: interp3d
    use char_module, only: ismember
    use mpi_module, only: lon0, lon1
    use netcdf, only: nf90_close, nf90_open, nf90_inq_varid, nf90_get_var, nf90_nowrite

    integer, intent(in) :: itime, it

    logical :: flag
    integer :: stat, ifld, i_ng, latbeg, latend, lonbeg, lonend, offbeg, offend, lonbeg1, lonend1, i0, i1, ilev
    real, dimension(1) :: zlbc
    real, dimension(maxval(flds%maxlev)) :: zout
    real, dimension(min(lon0, minval(ng_flds%lond0)): max(lon1, maxval(ng_flds%lond1))) :: model_lon
    real, dimension(minval(flds%latbeg): maxval(flds%latend)) :: model_lat
    real(kind=4), dimension(nlon, nlat, nlev, nf4d) :: f4d0
    real, dimension(nlev, nlon+1, nlat, nf4d) :: ncf4d
    real, dimension(1, min(lon0, minval(ng_flds%lond0)): max(lon1, maxval(ng_flds%lond1)), &
      minval(flds%latbeg): maxval(flds%latend)) :: lbc0

! open a new file to read if the current modelsec is over the last time index of the file
    if (itime > t1(ifile)) then
      stat = nf90_close(ncid)

      ifile = ifile + 1
      stat = nf90_open(trim(nudge_ncpre)//trim(nudge_ncfile(ifile))//trim(nudge_ncpost), nf90_nowrite, ncid)
      do ifld = 1, nf4d
        stat = nf90_inq_varid(ncid, trim(nudge_flds(ifld)), varid(ifld))
      enddo
    endif

    do ifld = 1, nf4d
      stat = nf90_get_var(ncid, varid(ifld), f4d0(:, :, :, ifld), &
        start=(/1, 1, 1, itime-t0(ifile)+1/), count=(/nlon, nlat, nlev, 1/))
    enddo

    zlbc(1) = zibot

! interpolate based on the subdomain intersection discussed in init

    do i_ng = 0, n_ng
      latbeg = flds(i_ng)%latbeg
      latend = flds(i_ng)%latend
      lonbeg = flds(i_ng)%lonbeg
      lonend = flds(i_ng)%lonend
      offbeg = flds(i_ng)%offbeg
      offend = flds(i_ng)%offend
      lonbeg1 = flds(i_ng)%lonbeg1
      lonend1 = flds(i_ng)%lonend1

      if (latbeg<latend .and. nudge_level(i_ng)) then
        if (i_ng == 0) then
          i0 = lon0
          i1 = lon1
          model_lon(i0: i1) = glon0(i0: i1)
          model_lat(latbeg: latend) = glat(latbeg: latend)
        else
          i0 = ng_flds(i_ng)%lond0
          i1 = ng_flds(i_ng)%lond1
          model_lon(i0: i1) = glon_ng(i_ng, i0: i1)
          model_lat(latbeg: latend) = glat_ng(i_ng, latbeg: latend)
        endif

        if (nudge_lbc) then
          if (wrap) then
            do ilev = 1, nlev
              ncf4d(ilev, 1: nlon-lonbeg, :, :) = f4d0(lonbeg+1: nlon, :, ilev, :)
              ncf4d(ilev, nlon-lonbeg+1: nlon, :, :) = f4d0(1: lonbeg, :, ilev, :)
            enddo
            ncf4d(:, nlon+1, :, :) = ncf4d(:, 1, :, :)
          else
            do ilev = 1, nlev
              ncf4d(ilev, 1: nlon, :, :) = f4d0(:, :, ilev, :)
            enddo
          endif

          do ifld = 1, nlb
            if (lb_idx(ifld) /= 0) then
              if (nudge_pert) then
                flag = .false.
              else
                flag = ismember(lb(ifld), zlog)
              endif

              if (wrap) then
                lbc0(:, i0: i1, latbeg: latend) = interp3d( &
                  zlbc, model_lon(i0: i1), model_lat(latbeg: latend), &
                  z0, z1, lon(lonbeg+1)+offbeg*360, lon(lonbeg+1)+offbeg*360+360, y0, y1, &
                  ncf4d(:, :, :, lb_idx(ifld)), flag)
                flds(i_ng)%lbc(:, :, it, ifld) = lbc0(1, i0: i1, latbeg: latend)
              else
                if (lonbeg<lonend .and. offbeg==offend) then
                  lbc0(:, lonbeg: lonend, latbeg: latend) = interp3d( &
                    zlbc, model_lon(lonbeg: lonend), model_lat(latbeg: latend), &
                    z0, z1, lon(1)+offbeg*360, lon(nlon)+offend*360, y0, y1, &
                    ncf4d(:, 1: nlon, :, lb_idx(ifld)), flag)
                  flds(i_ng)%lbc(lonbeg: lonend, :, it, ifld) = lbc0(1, lonbeg: lonend, latbeg: latend)
                endif

                if (lonbeg<lonend .and. offbeg+1==offend) then
                  lbc0(:, lonbeg: lonend1, latbeg: latend) = interp3d( &
                    zlbc, model_lon(lonbeg: lonend1), model_lat(latbeg: latend), &
                    z0, z1, lon(1)+offbeg*360, lon(nlon)+offbeg*360, y0, y1, &
                    ncf4d(:, 1: nlon, :, lb_idx(ifld)), flag)
                  flds(i_ng)%lbc(lonbeg: lonend1, :, it, ifld) = lbc0(1, lonbeg: lonend1, latbeg: latend)
                  lbc0(:, lonbeg1: lonend, latbeg: latend) = interp3d( &
                    zlbc, model_lon(lonbeg1: lonend), model_lat(latbeg: latend), &
                    z0, z1, lon(1)+offend*360, lon(nlon)+offend*360, y0, y1, &
                    ncf4d(:, 1: nlon, :, lb_idx(ifld)), flag)
                  flds(i_ng)%lbc(lonbeg1: lonend, :, it, ifld) = lbc0(1, lonbeg1: lonend, latbeg: latend)
                endif

                if (lonbeg>lonend .and. offbeg==offend+1) then
                  lbc0(:, lonbeg: i1, latbeg: latend) = interp3d( &
                    zlbc, model_lon(lonbeg: i1), model_lat(latbeg: latend), &
                    z0, z1, lon(1)+offbeg*360, lon(nlon)+offbeg*360, y0, y1, &
                    ncf4d(:, 1: nlon, :, lb_idx(ifld)), flag)
                  flds(i_ng)%lbc(lonbeg: i1, :, it, ifld) = lbc0(1, lonbeg: i1, latbeg: latend)
                  lbc0(:, i0: lonend, latbeg: latend) = interp3d( &
                    zlbc, model_lon(i0: lonend), model_lat(latbeg: latend), &
                    z0, z1, lon(1)+offend*360, lon(nlon)+offend*360, y0, y1, &
                    ncf4d(:, 1: nlon, :, lb_idx(ifld)), flag)
                  flds(i_ng)%lbc(i0: lonend, :, it, ifld) = lbc0(1, i0: lonend, latbeg: latend)
                endif
              endif
            endif
          enddo
        endif

        if (nudge_f4d) then
          if (wrap) then
            do ilev = 1, nlev
              ncf4d(ilev, 1: nlon-lonbeg, :, :) = f4d0(lonbeg+1: nlon, :, ilev, :)
              ncf4d(ilev, nlon-lonbeg+1: nlon, :, :) = f4d0(1: lonbeg, :, ilev, :)
            enddo
            ncf4d(:, nlon+1, :, :) = ncf4d(:, 1, :, :)
          else
            do ilev = 1, nlev
              ncf4d(ilev, 1: nlon, :, :) = f4d0(:, :, ilev, :)
            enddo
          endif

          do ifld = 1, nf4d
            if (nudge_pert) then
              flag = .false.
            else
              flag = ismember(f4d(f4d_idx(ifld))%short_name, zlog)
            endif
            if (trim(f4d(f4d_idx(ifld))%vcoord) == 'midpoints') then
              zout(1: flds(i_ng)%maxlev) = zpmid(1: flds(i_ng)%maxlev)
            else
              zout(1: flds(i_ng)%maxlev) = zpint(1: flds(i_ng)%maxlev)
            endif

            if (wrap) then
              flds(i_ng)%f4d(:, :, :, it, ifld) = interp3d( &
                zout(1: flds(i_ng)%maxlev), model_lon(i0: i1), model_lat(latbeg: latend), &
                z0, z1, lon(lonbeg+1)+offbeg*360, lon(lonbeg+1)+offbeg*360+360, y0, y1, &
                ncf4d(:, :, :, ifld), flag)
            else
              if (lonbeg<lonend .and. offbeg==offend) &
                flds(i_ng)%f4d(:, lonbeg: lonend, :, it, ifld) = interp3d( &
                zout(1: flds(i_ng)%maxlev), model_lon(lonbeg: lonend), model_lat(latbeg: latend), &
                z0, z1, lon(1)+offbeg*360, lon(nlon)+offend*360, y0, y1, &
                ncf4d(:, 1: nlon, :, ifld), flag)

              if (lonbeg<lonend .and. offbeg+1==offend) then
                flds(i_ng)%f4d(:, lonbeg: lonend1, :, it, ifld) = interp3d( &
                  zout(1: flds(i_ng)%maxlev), model_lon(lonbeg: lonend1), model_lat(latbeg: latend), &
                  z0, z1, lon(1)+offbeg*360, lon(nlon)+offbeg*360, y0, y1, &
                  ncf4d(:, 1: nlon, :, ifld), flag)
                flds(i_ng)%f4d(:, lonbeg1: lonend, :, it, ifld) = interp3d( &
                  zout(1: flds(i_ng)%maxlev), model_lon(lonbeg1: lonend), model_lat(latbeg: latend), &
                  z0, z1, lon(1)+offend*360, lon(nlon)+offend*360, y0, y1, &
                  ncf4d(:, 1: nlon, :, ifld), flag)
              endif

              if (lonbeg>lonend .and. offbeg==offend+1) then
                flds(i_ng)%f4d(:, lonbeg: i1, :, it, ifld) = interp3d( &
                  zout(1: flds(i_ng)%maxlev), model_lon(lonbeg: i1), model_lat(latbeg: latend), &
                  z0, z1, lon(1)+offbeg*360, lon(nlon)+offbeg*360, y0, y1, &
                  ncf4d(:, 1: nlon, :, ifld), flag)
                flds(i_ng)%f4d(:, i0: lonend, :, it, ifld) = interp3d( &
                  zout(1: flds(i_ng)%maxlev), model_lon(i0: lonend), model_lat(latbeg: latend), &
                  z0, z1, lon(1)+offend*360, lon(nlon)+offend*360, y0, y1, &
                  ncf4d(:, 1: nlon, :, ifld), flag)
              endif
            endif
          enddo
        endif
      endif
    enddo

  end subroutine read_data
!-----------------------------------------------------------------------
  subroutine finalize

    use netcdf, only: nf90_close

    integer :: stat

    stat = nf90_close(ncid)

  end subroutine finalize
!-----------------------------------------------------------------------
end module nudge_ng_module

subroutine init_fields_ng
! initialize nested grid fields by interpolating global fields

  use params_module,only: n_ng,nlevp1,nlonp1,nlonp2,nlonp4,nlat, &
    nlevp1_ng,nlon_ng,zibot,zitop,glon0,glat,zpint_ng,glon_ng,glat_ng
  use fields_module,only: tlbc,ulbc,vlbc,tlbc_nm,ulbc_nm,vlbc_nm,f4d,nf4d,itp
  use fields_ng_module,only: flds,itp_ng=>itp,itc_ng=>itc,nmap,fmap,ubfill,zlog
  use interp_module,only: interp3d,interp2d
  use char_module,only: ismember
  use mpi_module,only: mxlon,mxlat,ntask,tasks,lon0,lon1,lat0,lat1
  use mpi_f08,only: mpi_allgather,mpi_real8,mpi_comm_world

  integer :: nx,ny,n,ifld,cnt3d,cnt2d,itask,i0,i1,j0,j1,sidx,i_ng
  real :: x0,x1,mid_lon
  integer,dimension(1) :: idx
  real,dimension(nlevp1,mxlon,mxlat,nmap) :: sendbuf3d
  real,dimension(nlevp1,mxlon,mxlat,nmap,0:ntask-1) :: recvbuf3d
  real,dimension(nlevp1,nlonp4,nlat,nmap) :: full3d,tmp3d
  real,dimension(mxlon,mxlat,6) :: sendbuf2d
  real,dimension(mxlon,mxlat,6,0:ntask-1) :: recvbuf2d
  real,dimension(nlonp4,nlat,6) :: full2d,tmp2d
  real,dimension(minval(flds%lond0):maxval(flds%lond1),minval(flds%latd0):maxval(flds%latd1),6) :: tmplbc

! initialize sendbuf and recvbuf to 0
  nx = lon1-lon0+1
  ny = lat1-lat0+1

  sendbuf3d = 0.
  n = 1
  do ifld = 1,nf4d
    if (ismember(f4d(ifld)%short_name,fmap)) then
      sendbuf3d(:,1:nx,1:ny,n) = f4d(ifld)%data(:,lon0:lon1,lat0:lat1,itp)
! fill upper boundary with valid numbers (extrapolation)
      if (ismember(f4d(ifld)%short_name,ubfill)) then
        if (ismember(f4d(ifld)%short_name,zlog)) then
          sendbuf3d(nlevp1,:,:,n) = sendbuf3d(nlevp1-1,:,:,n)**2/sendbuf3d(nlevp1-2,:,:,n)
        else
          sendbuf3d(nlevp1,:,:,n) = sendbuf3d(nlevp1-1,:,:,n)*2-sendbuf3d(nlevp1-2,:,:,n)
        endif
      endif
      n = n+1
    endif
  enddo

  sendbuf2d = 0.
  sendbuf2d(1:nx,1:ny,1) = tlbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,2) = ulbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,3) = vlbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,4) = tlbc_nm(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,5) = ulbc_nm(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,6) = vlbc_nm(lon0:lon1,lat0:lat1)

  recvbuf3d = 0.
  recvbuf2d = 0.

  cnt3d = nlevp1*mxlon*mxlat*nmap
  cnt2d = mxlon*mxlat*6

  call mpi_allgather(sendbuf3d,cnt3d,mpi_real8,recvbuf3d,cnt3d,mpi_real8,mpi_comm_world)
  call mpi_allgather(sendbuf2d,cnt2d,mpi_real8,recvbuf2d,cnt2d,mpi_real8,mpi_comm_world)

! reconstruct global fields from received global subdomains
  do itask = 0,ntask-1
    i0 = tasks(itask)%lon0
    i1 = tasks(itask)%lon1
    j0 = tasks(itask)%lat0
    j1 = tasks(itask)%lat1
    full3d(:,i0:i1,j0:j1,:) = recvbuf3d(:,1:i1-i0+1,1:j1-j0+1,:,itask)
    full2d(i0:i1,j0:j1,:) = recvbuf2d(1:i1-i0+1,1:j1-j0+1,:,itask)
  enddo

  x0 = glon0(1)
  x1 = glon0(nlonp4)
  nx = nlonp4

! move global range to cover nested grid range
  if (.not. (glon0(1)<=glon_ng(1,-1) .and. glon_ng(1,nlon_ng(1)+2)<=glon0(nlonp4))) then
    mid_lon = (glon_ng(1,-1)+glon_ng(1,nlon_ng(1)+2))/2
    if (glon0(1)<=mid_lon+180 .and. mid_lon+180<=glon0(nlonp4)) then
      idx = minloc(abs(glon0-(mid_lon+180)))
      sidx = idx(1)
      x0 = glon0(sidx)-360
      x1 = glon0(sidx)
    else
      idx = minloc(abs(glon0-(mid_lon-180)))
      sidx = idx(1)
      x0 = glon0(sidx)
      x1 = glon0(sidx)+360
    endif
    tmp3d = full3d
    full3d(:,1:nlonp2+1-sidx,:,:) = tmp3d(:,sidx:nlonp2,:,:)
    full3d(:,nlonp4-sidx:nlonp1,:,:) = tmp3d(:,3:sidx,:,:)
    tmp2d = full2d
    full2d(1:nlonp2+1-sidx,:,:) = tmp2d(sidx:nlonp2,:,:)
    full2d(nlonp4-sidx:nlonp1,:,:) = tmp2d(3:sidx,:,:)
    nx = nlonp1
  endif

! interpolate to nested grid subdomain
  do i_ng = 1,n_ng
    i0 = flds(i_ng)%lond0
    i1 = flds(i_ng)%lond1
    j0 = flds(i_ng)%latd0
    j1 = flds(i_ng)%latd1

    n = 1
    do ifld = 1,nf4d
      if (ismember(f4d(ifld)%short_name,fmap)) then
        flds(i_ng)%f4d(ifld)%data(:,:,:,itp_ng(i_ng)) = &
          interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)), &
          glon_ng(i_ng,i0:i1),glat_ng(i_ng,j0:j1), &
          zibot,zitop,x0,x1,glat(1),glat(nlat), &
          full3d(:,1:nx,:,n),ismember(f4d(ifld)%short_name,zlog))
        flds(i_ng)%f4d(ifld)%data(:,:,:,itc_ng(i_ng)) = flds(i_ng)%f4d(ifld)%data(:,:,:,itp_ng(i_ng))
        n = n+1
      else
        flds(i_ng)%f4d(ifld)%data(:,:,:,itp_ng(i_ng)) = 0.
        flds(i_ng)%f4d(ifld)%data(:,:,:,itc_ng(i_ng)) = 0.
      endif
    enddo

    do ifld = 1,6
      tmplbc(i0:i1,j0:j1,ifld) = interp2d( &
        glon_ng(i_ng,i0:i1),glat_ng(i_ng,j0:j1), &
        x0,x1,glat(1),glat(nlat),full2d(1:nx,:,ifld))
    enddo
    flds(i_ng)%tlbc = tmplbc(i0:i1,j0:j1,1)
    flds(i_ng)%ulbc = tmplbc(i0:i1,j0:j1,2)
    flds(i_ng)%vlbc = tmplbc(i0:i1,j0:j1,3)
    flds(i_ng)%tlbc_nm = tmplbc(i0:i1,j0:j1,4)
    flds(i_ng)%ulbc_nm = tmplbc(i0:i1,j0:j1,5)
    flds(i_ng)%vlbc_nm = tmplbc(i0:i1,j0:j1,6)
  enddo

end subroutine init_fields_ng

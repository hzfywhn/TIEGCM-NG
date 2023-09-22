subroutine init_fields_ng
! initialize nested grid fields by interpolating global fields
! please see level_outer_ng and levels_inner_ng for details
! this subroutine is called once per model run

  use params_module,only: n_ng,nlevp1,nlon,nlonp1,nlonp2,nlonp4,nlat, &
    nlevp1_ng,nlon_ng,nlat_ng,zpint,glon0,glat,zpint_ng,glon_ng,glat_ng
  use fields_module,only: f4d,nf4d,itp
  use pdynamo_module,only: ex,ey,ez
  use lbc,only: t_lbc,u_lbc,v_lbc,z_lbc,flx_he
  use amie_module,only: ekvg,efxg
  use fields_ng_module,only: nf3din,nf2din,flds,itp_ng=>itp,itc_ng=>itc,ubfill,zlog,bndry
  use interp_module,only: interp3d,interp2d
  use char_module,only: ismember
  use mpi_module,only: mxlon,mxlat,ntask,tasks,lon0,lon1,lat0,lat1,TIEGCM_WORLD
  use mpi
  implicit none

  integer :: nx,ny,ifld,cnt4d,cnt3d,cnt2d,itask,i0,i1,j0,j1,sidx,i_ng,i,lat,n,ierror
  real :: mid_lon
  integer,dimension(1) :: idx
  real,dimension(1) :: bnd
  real,dimension(nlonp4) :: glon0_r
  real,dimension(nlevp1,mxlon,mxlat,nf4d) :: sendbuf4d
  real,dimension(nlevp1,mxlon,mxlat,nf4d,0:ntask-1) :: recvbuf4d
  real,dimension(nlevp1,nlonp4,nlat,nf4d) :: full4d,tmp4d
  real,dimension(nlevp1,mxlon,mxlat,nf3din) :: sendbuf3d
  real,dimension(nlevp1,mxlon,mxlat,nf3din,0:ntask-1) :: recvbuf3d
  real,dimension(nlevp1,nlonp4,nlat,nf3din) :: full3d,tmp3d
  real,dimension(mxlon,mxlat,5) :: sendbuf2d
  real,dimension(mxlon,mxlat,5,0:ntask-1) :: recvbuf2d
  real,dimension(nlonp4,nlat,nf2din) :: full2d,tmp2d
  real,dimension(maxval(nlevp1_ng),1,minval(flds%latd0):maxval(flds%latd1)) :: bndx
  real,dimension(maxval(nlevp1_ng),minval(flds%lond0):maxval(flds%lond1),1) :: bndy
  external :: shutdown

! initialize sendbuf and recvbuf to 0 for 4d/3d/2d fields
  nx = lon1-lon0+1
  ny = lat1-lat0+1

  sendbuf4d = 0.
  do ifld = 1,nf4d
    sendbuf4d(:,1:nx,1:ny,ifld) = f4d(ifld)%data(:,lon0:lon1,lat0:lat1,itp)
! fill upper boundary with valid numbers (extrapolation)
    if (ismember(f4d(ifld)%short_name,ubfill)) then
      if (ismember(f4d(ifld)%short_name,zlog)) then
        sendbuf4d(nlevp1,:,:,ifld) = sendbuf4d(nlevp1-1,:,:,ifld)**2/sendbuf4d(nlevp1-2,:,:,ifld)
      else
        sendbuf4d(nlevp1,:,:,ifld) = sendbuf4d(nlevp1-1,:,:,ifld)*2-sendbuf4d(nlevp1-2,:,:,ifld)
      endif
    endif
  enddo
  recvbuf4d = 0.
  cnt4d = nlevp1*mxlon*mxlat*nf4d

  sendbuf3d = 0.
  sendbuf3d(:,1:nx,1:ny,1) = ex(:,lon0:lon1,:)
  sendbuf3d(:,1:nx,1:ny,2) = ey(:,lon0:lon1,:)
  sendbuf3d(:,1:nx,1:ny,3) = ez(:,lon0:lon1,:)
  recvbuf3d = 0.
  cnt3d = nlevp1*mxlon*mxlat*nf3din

  sendbuf2d = 0.
  sendbuf2d(1:nx,1:ny,1) = t_lbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,2) = u_lbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,3) = v_lbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,4) = z_lbc(lon0:lon1,lat0:lat1)
  sendbuf2d(1:nx,1:ny,5) = flx_he(lon0:lon1,lat0:lat1)
  recvbuf2d = 0.
  cnt2d = mxlon*mxlat*5

! collective communication is used only in initialization
  call mpi_allgather(sendbuf4d,cnt4d,mpi_real8,recvbuf4d,cnt4d,mpi_real8,TIEGCM_WORLD,ierror)
  if (ierror /= mpi_success) call shutdown('failed to gather 4d fields to each process')

  call mpi_allgather(sendbuf3d,cnt3d,mpi_real8,recvbuf3d,cnt3d,mpi_real8,TIEGCM_WORLD,ierror)
  if (ierror /= mpi_success) call shutdown('failed to gather 3d fields to each process')

  call mpi_allgather(sendbuf2d,cnt2d,mpi_real8,recvbuf2d,cnt2d,mpi_real8,TIEGCM_WORLD,ierror)
  if (ierror /= mpi_success) call shutdown('failed to gather 2d fields to each process')

! reconstruct global fields from received global subdomains
  do itask = 0,ntask-1
    i0 = tasks(itask)%lon0
    i1 = tasks(itask)%lon1
    j0 = tasks(itask)%lat0
    j1 = tasks(itask)%lat1
    full4d(:,i0:i1,j0:j1,:) = recvbuf4d(:,1:i1-i0+1,1:j1-j0+1,:,itask)
    full3d(:,i0:i1,j0:j1,:) = recvbuf3d(:,1:i1-i0+1,1:j1-j0+1,:,itask)
    full2d(i0:i1,j0:j1,1:5) = recvbuf2d(1:i1-i0+1,1:j1-j0+1,:,itask)
  enddo

  full2d(3:nlonp2,:,6) = ekvg(1:nlon,1:nlat)
  full2d(1:2,:,6) = ekvg(nlon-1:nlon,1:nlat)
  full2d(nlonp2+1:nlonp4,:,6) = ekvg(1:2,1:nlat)
  full2d(3:nlonp2,:,7) = efxg(1:nlon,1:nlat)
  full2d(1:2,:,7) = efxg(nlon-1:nlon,1:nlat)
  full2d(nlonp2+1:nlonp4,:,7) = efxg(1:2,1:nlat)

  glon0_r = glon0
  nx = nlonp4

! move global range to cover nested grid range
  if (.not. (glon0(1)<=glon_ng(1,-1) .and. glon_ng(1,nlon_ng(1)+2)<=glon0(nlonp4))) then
    mid_lon = (glon_ng(1,-1)+glon_ng(1,nlon_ng(1)+2))/2
    if (glon0(1)<=mid_lon+180 .and. mid_lon+180<=glon0(nlonp4)) then
      idx = minloc(abs(glon0-(mid_lon+180)))
      sidx = idx(1)
      glon0_r(1:nlonp2+1-sidx) = glon0(sidx:nlonp2)-360
      glon0_r(nlonp4-sidx:nlonp1) = glon0(3:sidx)
    else
      idx = minloc(abs(glon0-(mid_lon-180)))
      sidx = idx(1)
      glon0_r(1:nlonp2+1-sidx) = glon0(sidx:nlonp2)
      glon0_r(nlonp4-sidx:nlonp1) = glon0(3:sidx)+360
    endif
    tmp4d = full4d
    full4d(:,1:nlonp2+1-sidx,:,:) = tmp4d(:,sidx:nlonp2,:,:)
    full4d(:,nlonp4-sidx:nlonp1,:,:) = tmp4d(:,3:sidx,:,:)
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

! initialize 4d/3d/2d fields
    do ifld = 1,nf4d
      flds(i_ng)%f4d(ifld)%data(:,:,:,itp_ng(i_ng)) = &
        interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)),glon_ng(i_ng,i0:i1),glat_ng(i_ng,j0:j1), &
        zpint,glon0_r(1:nx),glat,full4d(:,1:nx,:,ifld),ismember(f4d(ifld)%short_name,zlog))
      flds(i_ng)%f4d(ifld)%data(:,:,:,itc_ng(i_ng)) = flds(i_ng)%f4d(ifld)%data(:,:,:,itp_ng(i_ng))
    enddo

    do ifld = 1,nf3din
      flds(i_ng)%f3d_save(:,:,:,0,ifld) = &
        interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)),glon_ng(i_ng,i0:i1),glat_ng(i_ng,j0:j1), &
        zpint,glon0_r(1:nx),glat,full3d(:,1:nx,:,ifld))
    enddo

    do ifld = 1,nf2din
      flds(i_ng)%f2d_save(:,:,0,ifld) = &
        interp2d(glon_ng(i_ng,i0:i1),glat_ng(i_ng,j0:j1),glon0_r(1:nx),glat,full2d(1:nx,:,ifld))
    enddo

! initialize lateral boundaries
    do i = 1,2
      if (flds(i_ng)%is_bndry(i)) then
        if (i == 1) then
          bnd(1) = glon_ng(i_ng,-1)
        else
          bnd(1) = glon_ng(i_ng,nlon_ng(i_ng)+2)
        endif
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,bndry)) then
            bndx(1:nlevp1_ng(i_ng),:,j0:j1) = &
              interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)),bnd,glat_ng(i_ng,j0:j1), &
              zpint,glon0_r(1:nx),glat,full4d(:,1:nx,:,ifld),ismember(f4d(ifld)%short_name,zlog))
            flds(i_ng)%lon_b(:,i,:,0,n) = bndx(1:nlevp1_ng(i_ng),1,j0:j1)
            n = n+1
          endif
        enddo
! now should have n==nbnd
      endif
    enddo

    do lat = 3,4
      if (flds(i_ng)%is_bndry(lat)) then
        if (lat == 3) then
          bnd(1) = glat_ng(i_ng,-1)
        else
          bnd(1) = glat_ng(i_ng,nlat_ng(i_ng)+2)
        endif
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,bndry)) then
            bndy(1:nlevp1_ng(i_ng),i0:i1,:) = &
              interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)),glon_ng(i_ng,i0:i1),bnd, &
              zpint,glon0_r(1:nx),glat,full4d(:,1:nx,:,ifld),ismember(f4d(ifld)%short_name,zlog))
            flds(i_ng)%lat_b(:,:,lat,0,n) = bndy(1:nlevp1_ng(i_ng),i0:i1,1)
            n = n+1
          endif
        enddo
      endif
    enddo

! initialize other essential fields
    do ifld = 1,nf3din
      flds(i_ng)%f3din(ifld)%data = flds(i_ng)%f3d_save(:,:,:,0,ifld)
    enddo
    do ifld = 1,nf2din
      flds(i_ng)%f2din(ifld)%data = flds(i_ng)%f2d_save(:,:,0,ifld)
    enddo
    flds(i_ng)%tlbc = flds(i_ng)%t_lbc
    flds(i_ng)%ulbc = flds(i_ng)%u_lbc
    flds(i_ng)%vlbc = flds(i_ng)%v_lbc
  enddo

end subroutine init_fields_ng

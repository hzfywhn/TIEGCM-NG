module levels_inner_ng_module
! data transfer between different levels of nested grid subdomains, see level_outer_ng for descriptions

  implicit none

  integer,dimension(:,:),allocatable :: send_outer,recv_outer,send_inner,recv_inner
  integer,dimension(:,:,:),allocatable :: send_bndry,recv_bndry

  contains
!-----------------------------------------------------------------------
  subroutine init_mapping

    use params_module,only: n_ng,nlon_ng,nlat_ng,glon_ng,glat_ng
    use fields_ng_module,only: flds,domain
    use mpi_module,only: ntask
    use mpi_f08,only: mpi_proc_null

    integer :: i_ng,itask
    real :: inner_lb,inner_rb,inner_tb,inner_bb, &
      inner_left,inner_right,inner_top,inner_bottom, &
      outer_left,outer_right,outer_top,outer_bottom, &
      cur_inner_left,cur_inner_right,cur_inner_top,cur_inner_bottom, &
      cur_outer_left,cur_outer_right,cur_outer_top,cur_outer_bottom

    allocate(send_outer(n_ng,0:ntask-1))
    allocate(recv_outer(n_ng,0:ntask-1))
    allocate(send_inner(n_ng,0:ntask-1))
    allocate(recv_inner(n_ng,0:ntask-1))
    allocate(send_bndry(n_ng,4,0:ntask-1))
    allocate(recv_bndry(n_ng,4,0:ntask-1))

    send_outer = mpi_proc_null
    recv_outer = mpi_proc_null
    send_inner = mpi_proc_null
    recv_inner = mpi_proc_null
    send_bndry = mpi_proc_null
    recv_bndry = mpi_proc_null

! i_ng==1 is handled in level_outer_ng, skip
! longitude wrap is handled in level_outer_ng, no need to wrap here
    do i_ng = 2,n_ng
      cur_outer_left = glon_ng(i_ng-1,flds(i_ng-1)%lond0)
      cur_outer_right = glon_ng(i_ng-1,flds(i_ng-1)%lond1)
      cur_outer_top = glat_ng(i_ng-1,flds(i_ng-1)%latd0)
      cur_outer_bottom = glat_ng(i_ng-1,flds(i_ng-1)%latd1)

      cur_inner_left = glon_ng(i_ng,flds(i_ng)%lond0)
      cur_inner_right = glon_ng(i_ng,flds(i_ng)%lond1)
      cur_inner_top = glat_ng(i_ng,flds(i_ng)%latd0)
      cur_inner_bottom = glat_ng(i_ng,flds(i_ng)%latd1)

      do itask = 0,ntask-1
        outer_left = glon_ng(i_ng-1,domain(i_ng-1,1,itask)-2)
        outer_right = glon_ng(i_ng-1,domain(i_ng-1,2,itask)+2)
        outer_top = glat_ng(i_ng-1,domain(i_ng-1,3,itask)-2)
        outer_bottom = glat_ng(i_ng-1,domain(i_ng-1,4,itask)+2)

        inner_left = glon_ng(i_ng,domain(i_ng,1,itask)-2)
        inner_right = glon_ng(i_ng,domain(i_ng,2,itask)+2)
        inner_top = glat_ng(i_ng,domain(i_ng,3,itask)-2)
        inner_bottom = glat_ng(i_ng,domain(i_ng,4,itask)+2)

        if (cur_outer_left<=inner_right .and. inner_left<=cur_outer_right .and. &
          cur_outer_top<=inner_bottom .and. inner_top<=cur_outer_bottom) &
          send_outer(i_ng,itask) = itask

        if (outer_left<=cur_inner_right .and. cur_inner_left<=outer_right .and. &
          outer_top<=cur_inner_bottom .and. cur_inner_top<=outer_bottom) &
          recv_outer(i_ng,itask) = itask

        if (cur_inner_left<=outer_right .and. outer_left<=cur_inner_right .and. &
          cur_inner_top<=outer_bottom .and. outer_top<=cur_inner_bottom) &
          send_inner(i_ng,itask) = itask

        if (inner_left<=cur_outer_right .and. cur_outer_left<=inner_right .and. &
          inner_top<=cur_outer_bottom .and. cur_outer_top<=inner_bottom) &
          recv_inner(i_ng,itask) = itask
      enddo

      inner_lb = glon_ng(i_ng,-1)
      inner_rb = glon_ng(i_ng,nlon_ng(i_ng)+2)
      inner_tb = glat_ng(i_ng,-1)
      inner_bb = glat_ng(i_ng,nlat_ng(i_ng)+2)

      do itask = 0,ntask-1
        outer_left = glon_ng(i_ng-1,domain(i_ng-1,1,itask)-2)
        outer_right = glon_ng(i_ng-1,domain(i_ng-1,2,itask)+2)
        outer_top = glat_ng(i_ng-1,domain(i_ng-1,3,itask)-2)
        outer_bottom = glat_ng(i_ng-1,domain(i_ng-1,4,itask)+2)

        if (cur_outer_left<=inner_lb .and. inner_lb<=cur_outer_right .and. &
          cur_outer_top<=inner_bb .and. inner_tb<=cur_outer_bottom .and. &
          domain(i_ng,1,itask)==1) send_bndry(i_ng,1,itask) = itask

        if (cur_outer_left<=inner_rb .and. inner_rb<=cur_outer_right .and. &
          cur_outer_top<=inner_bb .and. inner_tb<=cur_outer_bottom .and. &
          domain(i_ng,2,itask)==nlon_ng(i_ng)) send_bndry(i_ng,2,itask) = itask

        if (cur_outer_left<=inner_rb .and. inner_lb<=cur_outer_right .and. &
          cur_outer_top<=inner_tb .and. inner_tb<=cur_outer_bottom .and. &
          domain(i_ng,3,itask)==1) send_bndry(i_ng,3,itask) = itask

        if (cur_outer_left<=inner_rb .and. inner_lb<=cur_outer_right .and. &
          cur_outer_top<=inner_bb .and. inner_bb<=cur_outer_bottom .and. &
          domain(i_ng,4,itask)==nlat_ng(i_ng)) send_bndry(i_ng,4,itask) = itask

        if (flds(i_ng)%is_bndry(1) .and. &
          outer_left<=inner_lb .and. inner_lb<=outer_right .and. &
          outer_top<=inner_bb .and. inner_tb<=outer_bottom) &
          recv_bndry(i_ng,1,itask) = itask

        if (flds(i_ng)%is_bndry(2) .and. &
          outer_left<=inner_rb .and. inner_rb<=outer_right .and. &
          outer_top<=inner_bb .and. inner_tb<=outer_bottom) &
          recv_bndry(i_ng,2,itask) = itask

        if (flds(i_ng)%is_bndry(3) .and. &
          outer_left<=inner_rb .and. inner_lb<=outer_right .and. &
          outer_top<=inner_tb .and. inner_tb<=outer_bottom) &
          recv_bndry(i_ng,3,itask) = itask

        if (flds(i_ng)%is_bndry(4) .and. &
          outer_left<=inner_rb .and. inner_lb<=outer_right .and. &
          outer_top<=inner_bb .and. inner_bb<=outer_bottom) &
          recv_bndry(i_ng,4,itask) = itask
      enddo
    enddo

  end subroutine init_mapping
!-----------------------------------------------------------------------
  subroutine map_in(i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,zpint_ng,glon_ng,glat_ng,zibot,zitop
    use input_module,only: nstep_ng
    use fields_ng_module,only: flds,nf3din,nf2din,maxlon,maxlat,domain
    use interp_module,only: interp3d,interp2d
    use mpi_module,only: ntask
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng

    integer :: lon0,lon1,lat0,lat1,cnt3d,cnt2d,itask,ifld
    real,dimension(nlevp1_ng(i_ng-1),maxlon(i_ng-1),maxlat(i_ng-1),nf3din) :: sendbuf3d
    real,dimension(nlevp1_ng(i_ng-1),maxlon(i_ng-1),maxlat(i_ng-1),nf3din,0:ntask-1) :: recvbuf3d
    real,dimension(nlevp1_ng(i_ng-1),-1:nlon_ng(i_ng-1)+2,-1:nlat_ng(i_ng-1)+2,nf3din) :: full3d
    real,dimension(maxlon(i_ng-1),maxlat(i_ng-1),nf2din) :: sendbuf2d
    real,dimension(maxlon(i_ng-1),maxlat(i_ng-1),nf2din,0:ntask-1) :: recvbuf2d
    real,dimension(-1:nlon_ng(i_ng-1)+2,-1:nlat_ng(i_ng-1)+2,nf2din) :: full2d
    type(mpi_request),dimension(0:ntask*4-1) :: request
    external :: bndry_index_ng

    call bndry_index_ng((/flds(i_ng-1)%lon0,flds(i_ng-1)%lon1,flds(i_ng-1)%lat0,flds(i_ng-1)%lat1/), &
      nlon_ng(i_ng-1),nlat_ng(i_ng-1),lon0,lon1,lat0,lat1)

    sendbuf3d = 0.
    do ifld = 1,nf3din
      sendbuf3d(:,1:lon1-lon0+1,1:lat1-lat0+1,ifld) = flds(i_ng-1)%f3din(ifld)%data(:,lon0:lon1,lat0:lat1)
    enddo

    sendbuf2d = 0.
    do ifld = 1,nf2din
      sendbuf2d(1:lon1-lon0+1,1:lat1-lat0+1,ifld) = flds(i_ng-1)%f2din(ifld)%data(lon0:lon1,lat0:lat1)
    enddo

    recvbuf3d = 0.
    recvbuf2d = 0.

    cnt3d = nlevp1_ng(i_ng-1)*maxlon(i_ng-1)*maxlat(i_ng-1)*nf3din
    cnt2d = maxlon(i_ng-1)*maxlat(i_ng-1)*nf2din

    do itask = 0,ntask-1
      call mpi_isend(sendbuf3d,cnt3d,mpi_real8,send_outer(i_ng,itask), &
        1,mpi_comm_world,request(itask))
      call mpi_isend(sendbuf2d,cnt2d,mpi_real8,send_outer(i_ng,itask), &
        2,mpi_comm_world,request(itask+ntask))
    enddo

    do itask = 0,ntask-1
      call mpi_irecv(recvbuf3d(:,:,:,:,itask),cnt3d,mpi_real8,recv_outer(i_ng,itask), &
        1,mpi_comm_world,request(itask+ntask*2))
      call mpi_irecv(recvbuf2d(:,:,:,itask),cnt2d,mpi_real8,recv_outer(i_ng,itask), &
        2,mpi_comm_world,request(itask+ntask*3))
    enddo

    call mpi_waitall(ntask*4,request,mpi_statuses_ignore)

    do itask = 0,ntask-1
      call bndry_index_ng(domain(i_ng-1,:,itask),nlon_ng(i_ng-1),nlat_ng(i_ng-1),lon0,lon1,lat0,lat1)
      full3d(:,lon0:lon1,lat0:lat1,:) = recvbuf3d(:,1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
      full2d(lon0:lon1,lat0:lat1,:) = recvbuf2d(1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
    enddo

    do ifld = 1,nf3din
      flds(i_ng)%f3d_save(:,:,:,nstep_ng(i_ng),ifld) = interp3d( &
        zpint_ng(i_ng,1:nlevp1_ng(i_ng)), &
        glon_ng(i_ng,flds(i_ng)%lond0:flds(i_ng)%lond1), &
        glat_ng(i_ng,flds(i_ng)%latd0:flds(i_ng)%latd1), &
        zibot,zitop, &
        glon_ng(i_ng-1,-1),glon_ng(i_ng-1,nlon_ng(i_ng-1)+2), &
        glat_ng(i_ng-1,-1),glat_ng(i_ng-1,nlat_ng(i_ng-1)+2), &
        full3d(:,:,:,ifld))
    enddo

    do ifld = 1,nf2din
      flds(i_ng)%f2d_save(:,:,nstep_ng(i_ng),ifld) = interp2d( &
        glon_ng(i_ng,flds(i_ng)%lond0:flds(i_ng)%lond1), &
        glat_ng(i_ng,flds(i_ng)%latd0:flds(i_ng)%latd1), &
        glon_ng(i_ng-1,-1),glon_ng(i_ng-1,nlon_ng(i_ng-1)+2), &
        glat_ng(i_ng-1,-1),glat_ng(i_ng-1,nlat_ng(i_ng-1)+2), &
        full2d(:,:,ifld))
    enddo

  end subroutine map_in
!-----------------------------------------------------------------------
  subroutine map_out(i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,zpint_ng,glon_ng,glat_ng,zibot,zitop
    use fields_module,only: f4d,nf4d
    use fields_ng_module,only: flds,maxlon,maxlat,domain,itc,nmap,fmap,zlog
    use interp_module,only: interp3d
    use char_module,only: ismember
    use mpi_module,only: ntask
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng

    integer :: lon0,lon1,lat0,lat1,n,ifld,cnt,itask,lonbeg,lonend,latbeg,latend
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng),nmap) :: sendbuf
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng),nmap,0:ntask-1) :: recvbuf
    real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,nmap) :: full
    type(mpi_request),dimension(0:ntask*2-1) :: request
    external :: bndry_index_ng

    call bndry_index_ng((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
      nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    n = 1
    do ifld = 1,nf4d
      if (ismember(f4d(ifld)%short_name,fmap)) then
        sendbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,n) = flds(i_ng)%f4d(ifld)%data(:,lon0:lon1,lat0:lat1,itc(i_ng))
        n = n+1
      endif
    enddo
    recvbuf = 0.
    cnt = nlevp1_ng(i_ng)*maxlon(i_ng)*maxlat(i_ng)*nmap

    do itask = 0,ntask-1
      call mpi_isend(sendbuf,cnt,mpi_real8,send_inner(i_ng,itask), &
        0,mpi_comm_world,request(itask))
    enddo

    do itask = 0,ntask-1
      call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_inner(i_ng,itask), &
        0,mpi_comm_world,request(itask+ntask))
    enddo

    call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

    do lonbeg = flds(i_ng-1)%lond0,flds(i_ng-1)%lond1
      if (glon_ng(i_ng-1,lonbeg) > glon_ng(i_ng,-1)) exit
    enddo
    do lonend = flds(i_ng-1)%lond1,flds(i_ng-1)%lond0,-1
      if (glon_ng(i_ng-1,lonend) < glon_ng(i_ng,nlon_ng(i_ng)+2)) exit
    enddo
    do latbeg = flds(i_ng-1)%latd0,flds(i_ng-1)%latd1
      if (glat_ng(i_ng-1,latbeg) > glat_ng(i_ng,-1)) exit
    enddo
    do latend = flds(i_ng-1)%latd1,flds(i_ng-1)%latd0,-1
      if (glat_ng(i_ng-1,latend) < glat_ng(i_ng,nlat_ng(i_ng)+2)) exit
    enddo
    if (lonbeg>=lonend .or. latbeg>=latend) return

    do itask = 0,ntask-1
      call bndry_index_ng(domain(i_ng,:,itask),nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)
      full(:,lon0:lon1,lat0:lat1,:) = recvbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
    enddo

    n = 1
    do ifld = 1,nf4d
      if (ismember(f4d(ifld)%short_name,fmap)) then
        flds(i_ng-1)%f4d(ifld)%data(:,lonbeg:lonend,latbeg:latend,itc(i_ng-1)) = &
          interp3d(zpint_ng(i_ng-1,1:nlevp1_ng(i_ng-1)), &
          glon_ng(i_ng-1,lonbeg:lonend),glat_ng(i_ng-1,latbeg:latend), &
          zibot,zitop, &
          glon_ng(i_ng,-1),glon_ng(i_ng,nlon_ng(i_ng)+2), &
          glat_ng(i_ng,-1),glat_ng(i_ng,nlat_ng(i_ng)+2), &
          full(:,:,:,n),ismember(f4d(ifld)%short_name,zlog))
        n = n+1
      endif
    enddo

  end subroutine map_out
!-----------------------------------------------------------------------
  subroutine set_bndry(i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,zpint_ng,glon_ng,glat_ng,zibot,zitop
    use input_module,only: nstep_ng
    use fields_module,only: f4d,nf4d
    use fields_ng_module,only: flds,maxlon,maxlat,domain,itc,nbnd,bndry,zlog
    use char_module,only: ismember
    use interp_module,only: interp3d
    use mpi_module,only: ntask
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng

    integer :: lon0,lon1,lat0,lat1,n,if4d,cnt,itask,i,lat,ibnd
    real,dimension(1) :: bnd
    real,dimension(nlevp1_ng(i_ng-1),maxlon(i_ng-1),maxlat(i_ng-1),nbnd) :: sendbuf
    real,dimension(nlevp1_ng(i_ng-1),maxlon(i_ng-1),maxlat(i_ng-1),nbnd,0:ntask-1) :: recvbuf
    real,dimension(nlevp1_ng(i_ng-1),-1:nlon_ng(i_ng-1)+2,-1:nlat_ng(i_ng-1)+2,nbnd) :: full
    real,dimension(nlevp1_ng(i_ng),1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: bndx
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,1) :: bndy
    type(mpi_request),dimension(0:ntask*2-1) :: request
    external :: bndry_index_ng

    call bndry_index_ng((/flds(i_ng-1)%lon0,flds(i_ng-1)%lon1,flds(i_ng-1)%lat0,flds(i_ng-1)%lat1/), &
      nlon_ng(i_ng-1),nlat_ng(i_ng-1),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    n = 1
    do if4d = 1,nf4d
      if (ismember(f4d(if4d)%short_name,bndry)) then
        sendbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,n) = flds(i_ng-1)%f4d(if4d)%data(:,lon0:lon1,lat0:lat1,itc(i_ng-1))
        n = n+1
      endif
    enddo
    recvbuf = 0.
    cnt = nlevp1_ng(i_ng-1)*maxlon(i_ng-1)*maxlat(i_ng-1)*nbnd

    do i = 1,2
      do itask = 0,ntask-1
        call mpi_isend(sendbuf,cnt,mpi_real8,send_bndry(i_ng,i,itask), &
          0,mpi_comm_world,request(itask))
      enddo

      do itask = 0,ntask-1
        call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_bndry(i_ng,i,itask), &
          0,mpi_comm_world,request(itask+ntask))
      enddo

      call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

      if (flds(i_ng)%is_bndry(i)) then
        do itask = 0,ntask-1
          call bndry_index_ng(domain(i_ng-1,:,itask),nlon_ng(i_ng-1),nlat_ng(i_ng-1),lon0,lon1,lat0,lat1)
          full(:,lon0:lon1,lat0:lat1,:) = recvbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
        enddo

        if (i == 1) then
          bnd(1) = glon_ng(i_ng,-1)
        else
          bnd(1) = glon_ng(i_ng,nlon_ng(i_ng)+2)
        endif

        do ibnd = 1,nbnd
          bndx = interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)), &
            bnd,glat_ng(i_ng,flds(i_ng)%latd0:flds(i_ng)%latd1), &
            zibot,zitop, &
            glon_ng(i_ng-1,-1),glon_ng(i_ng-1,nlon_ng(i_ng-1)+2), &
            glat_ng(i_ng-1,-1),glat_ng(i_ng-1,nlat_ng(i_ng-1)+2), &
            full(:,:,:,ibnd),ismember(bndry(ibnd),zlog))
          flds(i_ng)%lon_b(:,i,:,nstep_ng(i_ng),ibnd) = bndx(:,1,:)
        enddo
      endif
    enddo

    do lat = 3,4
      do itask = 0,ntask-1
        call mpi_isend(sendbuf,cnt,mpi_real8,send_bndry(i_ng,lat,itask), &
          0,mpi_comm_world,request(itask))
      enddo

      do itask = 0,ntask-1
        call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_bndry(i_ng,lat,itask), &
          0,mpi_comm_world,request(itask+ntask))
      enddo

      call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

      if (flds(i_ng)%is_bndry(lat)) then
        do itask = 0,ntask-1
          call bndry_index_ng(domain(i_ng-1,:,itask),nlon_ng(i_ng-1),nlat_ng(i_ng-1),lon0,lon1,lat0,lat1)
          full(:,lon0:lon1,lat0:lat1,:) = recvbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
        enddo

        if (lat == 3) then
          bnd(1) = glat_ng(i_ng,-1)
        else
          bnd(1) = glat_ng(i_ng,nlat_ng(i_ng)+2)
        endif

        do ibnd = 1,nbnd
          bndy = interp3d(zpint_ng(i_ng,1:nlevp1_ng(i_ng)), &
            glon_ng(i_ng,flds(i_ng)%lond0:flds(i_ng)%lond1),bnd, &
            zibot,zitop, &
            glon_ng(i_ng-1,-1),glon_ng(i_ng-1,nlon_ng(i_ng-1)+2), &
            glat_ng(i_ng-1,-1),glat_ng(i_ng-1,nlat_ng(i_ng-1)+2), &
            full(:,:,:,ibnd),ismember(bndry(ibnd),zlog))
          flds(i_ng)%lat_b(:,:,lat,nstep_ng(i_ng),ibnd) = bndy(:,:,1)
        enddo
      endif
    enddo

  end subroutine set_bndry
!-----------------------------------------------------------------------
end module levels_inner_ng_module

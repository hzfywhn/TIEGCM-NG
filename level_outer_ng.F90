module level_outer_ng_module
! data transfer between global and nested grid subdomains
! use individual send/recv instead of collective operations to avoid unnecessary data transfers
! nested grid longitude range can be different from global domain, situations are discussed seperately

  implicit none

  integer,dimension(:),allocatable :: send_glb,recv_glb,send_ng,recv_ng
  integer,dimension(:,:),allocatable :: send_bndry,recv_bndry

  contains
!-----------------------------------------------------------------------
  subroutine init_mapping

    use params_module,only: nlonp4,nlon_ng,nlat_ng,glon0,glat,glon_ng,glat_ng,dlon,dlat
    use fields_ng_module,only: flds,domain
    use mpi_module,only: ntask,tasks,lon0,lon1,lat0,lat1
    use mpi_f08,only: mpi_proc_null

    integer :: itask
    real :: ng_lb,ng_rb,ng_tb,ng_bb, &
      ng_left,ng_right,ng_top,ng_bottom, &
      glb_left,glb_right,glb_top,glb_bottom, &
      cur_ng_left,cur_ng_right,cur_ng_top,cur_ng_bottom, &
      cur_glb_left,cur_glb_right,cur_glb_top,cur_glb_bottom

    allocate(send_glb(0:ntask-1))
    allocate(recv_glb(0:ntask-1))
    allocate(send_ng(0:ntask-1))
    allocate(recv_ng(0:ntask-1))
    allocate(send_bndry(4,0:ntask-1))
    allocate(recv_bndry(4,0:ntask-1))

    send_glb = mpi_proc_null
    recv_glb = mpi_proc_null
    send_ng = mpi_proc_null
    recv_ng = mpi_proc_null
    send_bndry = mpi_proc_null
    recv_bndry = mpi_proc_null

! global subdomains count from lon0/lon1,lat0/lat1
! nested grid subdomains count from lond0/lond1,latd0/latd1
! expand the global subdomain coverage to allow overlapping across processes
! in case nested grid subdomain boundaries lie between global subdomain boundaries
    cur_glb_left = glon0(lon0)-dlon*2
    cur_glb_right = glon0(lon1)+dlon*2
    cur_glb_top = glat(lat0)-dlat*2
    cur_glb_bottom = glat(lat1)+dlat*2

    cur_ng_left = glon_ng(1,flds(1)%lond0)
    cur_ng_right = glon_ng(1,flds(1)%lond1)
    cur_ng_top = glat_ng(1,flds(1)%latd0)
    cur_ng_bottom = glat_ng(1,flds(1)%latd1)

    do itask = 0,ntask-1
      glb_left = glon0(tasks(itask)%lon0)-dlon*2
      glb_right = glon0(tasks(itask)%lon1)+dlon*2
      glb_top = glat(tasks(itask)%lat0)-dlat*2
      glb_bottom = glat(tasks(itask)%lat1)+dlat*2

      ng_left = glon_ng(1,domain(1,1,itask)-2)
      ng_right = glon_ng(1,domain(1,2,itask)+2)
      ng_top = glat_ng(1,domain(1,3,itask)-2)
      ng_bottom = glat_ng(1,domain(1,4,itask)+2)

! If nested grid domain has a different range from global domain (0<lon<360 vs -180<lon<180),
! it is necessary to move the range before comparison. There are three situations:
! 1. ng_right <= glon0(nlonp4), no need to move
! 2. ng_left >= glon0(nlonp4), move the whole nested grid range to glon0(1)<=lon<glon0(nlonp4)
! 3. ng_left < glon0(nlonp4) < ng_right, compare lon<glon0(nlonp4) and lon>glon0(nlonp4) separately

! all processes: find processes whose nested grid subdomain intersect the current process's global subdomain
! used to send global subdomains (init_f4d/map_in), should be the same as recv_ng
      if (((ng_right<=glon0(nlonp4) .and. cur_glb_left<=ng_right .and. ng_left<=cur_glb_right) .or. &
        (ng_left>=glon0(nlonp4) .and. cur_glb_left<=ng_right-360 .and. ng_left-360<=cur_glb_right) .or. &
        (ng_left<glon0(nlonp4) .and. ng_right>glon0(nlonp4) .and. &
        (cur_glb_left<=ng_right-360 .or. ng_left<=cur_glb_right))) .and. &
        cur_glb_top<=ng_bottom .and. ng_top<=cur_glb_bottom) send_glb(itask) = itask

! all processes: find processes whose global subdomain intersect the current process's nested grid subdomain
! used to receive global subdomains (init_f4d/map_in), should be the same as send_ng
      if (((cur_ng_right<=glon0(nlonp4) .and. glb_left<=cur_ng_right .and. cur_ng_left<=glb_right) .or. &
        (cur_ng_left>=glon0(nlonp4) .and. glb_left<=cur_ng_right-360 .and. cur_ng_left-360<=glb_right) .or. &
        (cur_ng_left<glon0(nlonp4) .and. cur_ng_right>glon0(nlonp4) .and. &
        (glb_left<=cur_ng_right-360 .or. cur_ng_left<=glb_right))) .and. &
        glb_top<=cur_ng_bottom .and. cur_ng_top<=glb_bottom) recv_glb(itask) = itask

! all processes: find processes whose global subdomain intersect the current process's nested grid subdomain
! used to send nested grid subdomains (map_out), should be the same as recv_glb
      if (((cur_ng_right<=glon0(nlonp4) .and. cur_ng_left<=glb_right .and. glb_left<=cur_ng_right) .or. &
        (cur_ng_left>=glon0(nlonp4) .and. cur_ng_left-360<=glb_right .and. glb_left<=cur_ng_right-360) .or. &
        (cur_ng_left<glon0(nlonp4) .and. cur_ng_right>glon0(nlonp4) .and. &
        (cur_ng_left<=glb_right .or. glb_left<=cur_ng_right-360))) .and. &
        cur_ng_top<=glb_bottom .and. glb_top<=cur_ng_bottom) send_ng(itask) = itask

! all processes: find processes whose nested grid subdomain intersect the current process's global subdomain
! used to receive nested grid subdomains (map_out), should be the same as send_glb
      if (((ng_right<=glon0(nlonp4) .and. ng_left<=cur_glb_right .and. cur_glb_left<=ng_right) .or. &
        (ng_left>=glon0(nlonp4) .and. ng_left-360<=cur_glb_right .and. cur_glb_left<=ng_right-360) .or. &
        (ng_left<glon0(nlonp4) .and. ng_right>glon0(nlonp4) .and. &
        (ng_left<=cur_glb_right .or. cur_glb_left<=ng_right-360))) .and. &
        ng_top<=cur_glb_bottom .and. cur_glb_top<=ng_bottom) recv_ng(itask) = itask
    enddo

! nested grid boundaries are handled in a similar way
    ng_lb = glon_ng(1,-1)
    ng_rb = glon_ng(1,nlon_ng(1)+2)
    ng_tb = glat_ng(1,-1)
    ng_bb = glat_ng(1,nlat_ng(1)+2)

    do itask = 0,ntask-1
      glb_left = glon0(tasks(itask)%lon0)-dlon*2
      glb_right = glon0(tasks(itask)%lon1)+dlon*2
      glb_top = glat(tasks(itask)%lat0)-dlat*2
      glb_bottom = glat(tasks(itask)%lat1)+dlat*2

! processes whose global subdomain intersect nested grid boundary: find processes with nested grid boundary
! used to send global subdomains (set_bndry)
      if (((ng_lb<=glon0(nlonp4) .and. cur_glb_left<=ng_lb .and. ng_lb<=cur_glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. cur_glb_left<=ng_lb-360 .and. ng_lb-360<=cur_glb_right)) .and. &
        cur_glb_top<=ng_bb .and. ng_tb<=cur_glb_bottom .and. domain(1,1,itask)==1) &
        send_bndry(1,itask) = itask

      if (((ng_rb<=glon0(nlonp4) .and. cur_glb_left<=ng_rb .and. ng_rb<=cur_glb_right) .or. &
        (ng_rb>=glon0(nlonp4) .and. cur_glb_left<=ng_rb-360 .and. ng_rb-360<=cur_glb_right)) .and. &
        cur_glb_top<=ng_bb .and. ng_tb<=cur_glb_bottom .and. domain(1,2,itask)==nlon_ng(1)) &
        send_bndry(2,itask) = itask

      if (((ng_rb<=glon0(nlonp4) .and. cur_glb_left<=ng_rb .and. ng_lb<=cur_glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. cur_glb_left<=ng_rb-360 .and. ng_lb-360<=cur_glb_right) .or. &
        (ng_lb<glon0(nlonp4) .and. ng_rb>glon0(nlonp4) .and. (cur_glb_left<=ng_rb-360 .or. ng_lb<=cur_glb_right))) .and. &
        cur_glb_top<=ng_tb .and. ng_tb<=cur_glb_bottom .and. domain(1,3,itask)==1) &
        send_bndry(3,itask) = itask

      if (((ng_rb<=glon0(nlonp4) .and. cur_glb_left<=ng_rb .and. ng_lb<=cur_glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. cur_glb_left<=ng_rb-360 .and. ng_lb-360<=cur_glb_right) .or. &
        (ng_lb<glon0(nlonp4) .and. ng_rb>glon0(nlonp4) .and. (cur_glb_left<=ng_rb-360 .or. ng_lb<=cur_glb_right))) .and. &
        cur_glb_top<=ng_bb .and. ng_bb<=cur_glb_bottom .and. domain(1,4,itask)==nlat_ng(1)) &
        send_bndry(4,itask) = itask

! processes with nested grid boundary: find processes whose global subdomain intersect nested grid boundary
! used to receive global subdomains (set_bndry)
      if (flds(1)%is_bndry(1) .and. &
        ((ng_lb<=glon0(nlonp4) .and. glb_left<=ng_lb .and. ng_lb<=glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. glb_left<=ng_lb-360 .and. ng_lb-360<=glb_right)) .and. &
        glb_top<=ng_bb .and. ng_tb<=glb_bottom) recv_bndry(1,itask) = itask

      if (flds(1)%is_bndry(2) .and. &
        ((ng_rb<=glon0(nlonp4) .and. glb_left<=ng_rb .and. ng_rb<=glb_right) .or. &
        (ng_rb>=glon0(nlonp4) .and. glb_left<=ng_rb-360 .and. ng_rb-360<=glb_right)) .and. &
        glb_top<=ng_bb .and. ng_tb<=glb_bottom) recv_bndry(2,itask) = itask

      if (flds(1)%is_bndry(3) .and. &
        ((ng_rb<=glon0(nlonp4) .and. glb_left<=ng_rb .and. ng_lb<=glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. glb_left<=ng_rb-360 .and. ng_lb-360<=glb_right) .or. &
        (ng_lb<glon0(nlonp4) .and. ng_rb>glon0(nlonp4) .and. (glb_left<=ng_rb-360 .or. ng_lb<=glb_right))) .and. &
        glb_top<=ng_tb .and. ng_tb<=glb_bottom) recv_bndry(3,itask) = itask

      if (flds(1)%is_bndry(4) .and. &
        ((ng_rb<=glon0(nlonp4) .and. glb_left<=ng_rb .and. ng_lb<=glb_right) .or. &
        (ng_lb>=glon0(nlonp4) .and. glb_left<=ng_rb-360 .and. ng_lb-360<=glb_right) .or. &
        (ng_lb<glon0(nlonp4) .and. ng_rb>glon0(nlonp4) .and. (glb_left<=ng_rb-360 .or. ng_lb<=glb_right))) .and. &
        glb_top<=ng_bb .and. ng_bb<=glb_bottom) recv_bndry(4,itask) = itask
    enddo

  end subroutine init_mapping
!-----------------------------------------------------------------------
  subroutine map_in
! essential 2d/3d fields are mapped at every time

    use params_module,only: nlevp1,nlon,nlonp1,nlonp2,nlonp4,nlat,nlevp1_ng, &
      zibot,zitop,glon0,glat,zpint_ng,glon_ng,glat_ng
    use input_module,only: nstep_ng
    use pdynamo_module,only: ex,ey,ez
    use lbc,only: t_lbc,u_lbc,v_lbc,z_lbc,flx_he
    use amie_module,only: ekvg,efxg
    use fields_ng_module,only: nf3din,nf2din,flds
    use interp_module,only: interp3d,interp2d
    use mpi_module,only: mxlon,mxlat,ntask,tasks,lon0,lon1,lat0,lat1
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer :: nx,ny,cnt3d,cnt2d,itask,i0,i1,j0,j1,sidx,ifld
    real :: x0,x1,mid_lon
    integer,dimension(1) :: idx
    real,dimension(nlevp1,mxlon,mxlat,nf3din) :: sendbuf3d
    real,dimension(nlevp1,mxlon,mxlat,nf3din,0:ntask-1) :: recvbuf3d
    real,dimension(nlevp1,nlonp4,nlat,nf3din) :: full3d,tmp3d
    real,dimension(mxlon,mxlat,5) :: sendbuf2d
    real,dimension(mxlon,mxlat,5,0:ntask-1) :: recvbuf2d
    real,dimension(nlonp4,nlat,nf2din) :: full2d,tmp2d
    type(mpi_request),dimension(0:ntask*4-1) :: request

! initialize sendbuf and recvbuf to 0
    nx = lon1-lon0+1
    ny = lat1-lat0+1

    sendbuf3d = 0.
    sendbuf3d(:,1:nx,1:ny,1) = ex(:,lon0:lon1,:)
    sendbuf3d(:,1:nx,1:ny,2) = ey(:,lon0:lon1,:)
    sendbuf3d(:,1:nx,1:ny,3) = ez(:,lon0:lon1,:)

    sendbuf2d = 0.
    sendbuf2d(1:nx,1:ny,1) = t_lbc(lon0:lon1,lat0:lat1)
    sendbuf2d(1:nx,1:ny,2) = u_lbc(lon0:lon1,lat0:lat1)
    sendbuf2d(1:nx,1:ny,3) = v_lbc(lon0:lon1,lat0:lat1)
    sendbuf2d(1:nx,1:ny,4) = z_lbc(lon0:lon1,lat0:lat1)
    sendbuf2d(1:nx,1:ny,5) = flx_he(lon0:lon1,lat0:lat1)

    recvbuf3d = 0.
    recvbuf2d = 0.

    cnt3d = nlevp1*mxlon*mxlat*nf3din
    cnt2d = mxlon*mxlat*5

! send to processes whose nested grid subdomain intersect the current process's global subdomain
    do itask = 0,ntask-1
      call mpi_isend(sendbuf3d,cnt3d,mpi_real8,send_glb(itask), &
        1,mpi_comm_world,request(itask))
      call mpi_isend(sendbuf2d,cnt2d,mpi_real8,send_glb(itask), &
        2,mpi_comm_world,request(itask+ntask))
    enddo

! receive from processes whose global subdomain intersect the current process's nested grid subdomain
    do itask = 0,ntask-1
      call mpi_irecv(recvbuf3d(:,:,:,:,itask),cnt3d,mpi_real8,recv_glb(itask), &
        1,mpi_comm_world,request(itask+ntask*2))
      call mpi_irecv(recvbuf2d(:,:,:,itask),cnt2d,mpi_real8,recv_glb(itask), &
        2,mpi_comm_world,request(itask+ntask*3))
    enddo

! wait until data transfer complete
    call mpi_waitall(ntask*4,request,mpi_statuses_ignore)

! reconstruct global fields from received global subdomains
    do itask = 0,ntask-1
      i0 = tasks(itask)%lon0
      i1 = tasks(itask)%lon1
      j0 = tasks(itask)%lat0
      j1 = tasks(itask)%lat1
      full3d(:,i0:i1,j0:j1,:) = recvbuf3d(:,1:i1-i0+1,1:j1-j0+1,:,itask)
      full2d(i0:i1,j0:j1,1:5) = recvbuf2d(1:i1-i0+1,1:j1-j0+1,:,itask)
    enddo

    full2d(3:nlonp2,:,6) = ekvg(1:nlon,1:nlat)
    full2d(1:2,:,6) = ekvg(nlon-1:nlon,1:nlat)
    full2d(nlonp2+1:nlonp4,:,6) = ekvg(1:2,1:nlat)
    full2d(3:nlonp2,:,7) = efxg(1:nlon,1:nlat)
    full2d(1:2,:,7) = efxg(nlon-1:nlon,1:nlat)
    full2d(nlonp2+1:nlonp4,:,7) = efxg(1:2,1:nlat)

    x0 = glon0(1)
    x1 = glon0(nlonp4)
    nx = nlonp4

! move global range to cover nested grid range
    if (.not. (glon0(1)<=glon_ng(1,flds(1)%lond0) .and. glon_ng(1,flds(1)%lond1)<=glon0(nlonp4))) then
      mid_lon = (glon_ng(1,flds(1)%lond0)+glon_ng(1,flds(1)%lond1))/2
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
    do ifld = 1,nf3din
      flds(1)%f3d_save(:,:,:,nstep_ng(1),ifld) = interp3d(zpint_ng(1,1:nlevp1_ng(1)), &
        glon_ng(1,flds(1)%lond0:flds(1)%lond1), &
        glat_ng(1,flds(1)%latd0:flds(1)%latd1), &
        zibot,zitop,x0,x1,glat(1),glat(nlat),full3d(:,1:nx,:,ifld))
    enddo

    do ifld = 1,nf2din
      flds(1)%f2d_save(:,:,nstep_ng(1),ifld) = interp2d( &
        glon_ng(1,flds(1)%lond0:flds(1)%lond1), &
        glat_ng(1,flds(1)%latd0:flds(1)%latd1), &
        x0,x1,glat(1),glat(nlat),full2d(1:nx,:,ifld))
    enddo

  end subroutine map_in
!-----------------------------------------------------------------------
  subroutine map_out

    use params_module,only: nlevp1,nlonp4,nlevp1_ng,nlon_ng,nlat_ng, &
      zpint,glon0,glat,zibot,zitop,glon_ng,glat_ng
    use fields_module,only: f4d,nf4d,itc
    use fields_ng_module,only: flds,maxlon,maxlat,domain,itc_ng=>itc,nmap,fmap,zlog
    use interp_module,only: interp3d
    use char_module,only: ismember
    use mpi_module,only: ntask,lon0,lon1,lat0,lat1
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer :: lond0,lond1,latd0,latd1,n,ifld,cnt,itask,latbeg,latend,lonbeg,lonend
    real :: ng_lb,ng_rb,ng_tb,ng_bb
    real,dimension(nlevp1_ng(1),maxlon(1),maxlat(1),nmap) :: sendbuf
    real,dimension(nlevp1_ng(1),maxlon(1),maxlat(1),nmap,0:ntask-1) :: recvbuf
    real,dimension(nlevp1_ng(1),-1:nlon_ng(1)+2,-1:nlat_ng(1)+2,nmap) :: full
    type(mpi_request),dimension(0:ntask*2-1) :: request
    external :: bndry_index_ng

    call bndry_index_ng((/flds(1)%lon0,flds(1)%lon1,flds(1)%lat0,flds(1)%lat1/), &
      nlon_ng(1),nlat_ng(1),lond0,lond1,latd0,latd1)

    sendbuf = 0.
    n = 1
    do ifld = 1,nf4d
      if (ismember(f4d(ifld)%short_name,fmap)) then
        sendbuf(:,1:lond1-lond0+1,1:latd1-latd0+1,n) = &
          flds(1)%f4d(ifld)%data(:,lond0:lond1,latd0:latd1,itc_ng(1))
        n = n+1
      endif
    enddo
    recvbuf = 0.
    cnt = nlevp1_ng(1)*maxlon(1)*maxlat(1)*nmap

! send to processes whose global subdomain intersect the current process's nested grid subdomain
    do itask = 0,ntask-1
      call mpi_isend(sendbuf,cnt,mpi_real8,send_ng(itask), &
        0,mpi_comm_world,request(itask))
    enddo

! receive from processes whose nested grid subdomain intersect the current process's global subdomain
    do itask = 0,ntask-1
      call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_ng(itask), &
        0,mpi_comm_world,request(itask+ntask))
    enddo

    call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

! reconstruct nested grid fields from received nested grid subdomains
    do itask = 0,ntask-1
      call bndry_index_ng(domain(1,:,itask),nlon_ng(1),nlat_ng(1),lond0,lond1,latd0,latd1)
      full(:,lond0:lond1,latd0:latd1,:) = recvbuf(:,1:lond1-lond0+1,1:latd1-latd0+1,:,itask)
    enddo

! find global subdomain falling into the nested grid full domain
    ng_tb = glat_ng(1,-1)
    ng_bb = glat_ng(1,nlat_ng(1)+2)
    do latbeg = lat0,lat1
      if (glat(latbeg) > ng_tb) exit
    enddo
    do latend = lat1,lat0,-1
      if (glat(latend) < ng_bb) exit
    enddo
    if (latbeg >= latend) return

! note that there are also three longitudinal situations that need to be considered separately
    ng_lb = glon_ng(1,-1)
    ng_rb = glon_ng(1,nlon_ng(1)+2)

    if (ng_rb <= glon0(nlonp4)) then
      do lonbeg = lon0,lon1
        if (glon0(lonbeg) > ng_lb) exit
      enddo
      do lonend = lon1,lon0,-1
        if (glon0(lonend) < ng_rb) exit
      enddo
      if (lonbeg < lonend) then
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,fmap)) then
            f4d(ifld)%data(1:nlevp1-1,lonbeg:lonend,latbeg:latend,itc) = &
              interp3d(zpint(1:nlevp1-1),glon0(lonbeg:lonend),glat(latbeg:latend), &
              zibot,zitop,ng_lb,ng_rb,ng_tb,ng_bb, &
              full(:,:,:,n),ismember(f4d(ifld)%short_name,zlog))
            n = n+1
          endif
        enddo
      endif
    endif

    if (ng_lb >= glon0(nlonp4)) then
      do lonbeg = lon0,lon1
        if (glon0(lonbeg) > ng_lb-360) exit
      enddo
      do lonend = lon1,lon0,-1
        if (glon0(lonend) < ng_rb-360) exit
      enddo
      if (lonbeg < lonend) then
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,fmap)) then
            f4d(ifld)%data(1:nlevp1-1,lonbeg:lonend,latbeg:latend,itc) = &
              interp3d(zpint(1:nlevp1-1),glon0(lonbeg:lonend),glat(latbeg:latend), &
              zibot,zitop,ng_lb-360,ng_rb-360,ng_tb,ng_bb, &
              full(:,:,:,n),ismember(f4d(ifld)%short_name,zlog))
            n = n+1
          endif
        enddo
      endif
    endif

    if (ng_lb<glon0(nlonp4) .and. ng_rb>glon0(nlonp4)) then
      do lonbeg = lon0,lon1
        if (glon0(lonbeg) > ng_lb) exit
      enddo
      if (lonbeg < lon1) then
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,fmap)) then
            f4d(ifld)%data(1:nlevp1-1,lonbeg:lon1,latbeg:latend,itc) = &
              interp3d(zpint(1:nlevp1-1),glon0(lonbeg:lon1),glat(latbeg:latend), &
              zibot,zitop,ng_lb,ng_rb,ng_tb,ng_bb, &
              full(:,:,:,n),ismember(f4d(ifld)%short_name,zlog))
            n = n+1
          endif
        enddo
      endif

      do lonend = lon1,lon0,-1
        if (glon0(lonend) < ng_rb-360) exit
      enddo
      if (lon0 < lonend) then
        n = 1
        do ifld = 1,nf4d
          if (ismember(f4d(ifld)%short_name,fmap)) then
            f4d(ifld)%data(1:nlevp1-1,lon0:lonend,latbeg:latend,itc) = &
              interp3d(zpint(1:nlevp1-1),glon0(lon0:lonend),glat(latbeg:latend), &
              zibot,zitop,ng_lb-360,ng_rb-360,ng_tb,ng_bb, &
              full(:,:,:,n),ismember(f4d(ifld)%short_name,zlog))
            n = n+1
          endif
        enddo
      endif
    endif

  end subroutine map_out
!-----------------------------------------------------------------------
  subroutine set_bndry

    use params_module,only: nlevp1,nlonp1,nlonp2,nlonp4,nlat,nlevp1_ng,nlon_ng,nlat_ng, &
      zibot,zitop,glon0,glat,zpint_ng,glon_ng,glat_ng
    use input_module,only: nstep_ng
    use fields_module,only: f4d,nf4d,itc
    use fields_ng_module,only: flds,nbnd,bndry,ubfill,zlog
    use interp_module,only: interp3d
    use char_module,only: ismember
    use mpi_module,only: mxlon,mxlat,ntask,tasks,lon0,lon1,lat0,lat1
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall, &
      mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer :: n,if4d,cnt,itask,i0,i1,j0,j1,i,lat,ibnd,nx,sidx
    real :: x0,x1,mid_lon
    integer,dimension(1) :: idx
    real,dimension(1) :: bnd
    real,dimension(nlevp1,mxlon,mxlat,nbnd) :: sendbuf
    real,dimension(nlevp1,mxlon,mxlat,nbnd,0:ntask-1) :: recvbuf
    real,dimension(nlevp1,nlonp4,nlat,nbnd) :: full,tmp
    real,dimension(nlevp1_ng(1),1,flds(1)%latd0:flds(1)%latd1) :: bndx
    real,dimension(nlevp1_ng(1),flds(1)%lond0:flds(1)%lond1,1) :: bndy
    type(mpi_request),dimension(0:ntask*2-1) :: request

    sendbuf = 0.
    n = 1
    do if4d = 1,nf4d
! only necessary fields are transferred
      if (ismember(f4d(if4d)%short_name,bndry)) then
        sendbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,n) = f4d(if4d)%data(:,lon0:lon1,lat0:lat1,itc)
        if (ismember(f4d(if4d)%short_name,ubfill)) then
          if (ismember(f4d(if4d)%short_name,zlog)) then
            sendbuf(nlevp1,:,:,n) = sendbuf(nlevp1-1,:,:,n)**2/sendbuf(nlevp1-2,:,:,n)
          else
            sendbuf(nlevp1,:,:,n) = sendbuf(nlevp1-1,:,:,n)*2-sendbuf(nlevp1-2,:,:,n)
          endif
        endif
        n = n+1
      endif
    enddo
    recvbuf = 0.
    cnt = nlevp1*mxlon*mxlat*nbnd

! set up lon boundaries
    do i = 1,2
! processes whose global subdomain intersect nestes grid boundary:
! send the global subdomain to processes with nested grid boundary
! otherwise do nothing
      do itask = 0,ntask-1
        call mpi_isend(sendbuf,cnt,mpi_real8,send_bndry(i,itask), &
          0,mpi_comm_world,request(itask))
      enddo

! processes with nested grid boundary:
! receive the global subdomain from processes whose global subdomain intersect nestes grid boundary
! otherwise do nothing
      do itask = 0,ntask-1
        call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_bndry(i,itask), &
          0,mpi_comm_world,request(itask+ntask))
      enddo

      call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

      if (flds(1)%is_bndry(i)) then
        do itask = 0,ntask-1
          i0 = tasks(itask)%lon0
          i1 = tasks(itask)%lon1
          j0 = tasks(itask)%lat0
          j1 = tasks(itask)%lat1
          full(:,i0:i1,j0:j1,:) = recvbuf(:,1:i1-i0+1,1:j1-j0+1,:,itask)
        enddo

        if (i == 1) then
          bnd(1) = glon_ng(1,-1)
        else
          bnd(1) = glon_ng(1,nlon_ng(1)+2)
        endif

        if (bnd(1) >= glon0(nlonp4)) bnd(1) = bnd(1)-360

! interpolate to nested grid boundaries
        do ibnd = 1,nbnd
          bndx = interp3d(zpint_ng(1,1:nlevp1_ng(1)), &
            bnd,glat_ng(1,flds(1)%latd0:flds(1)%latd1), &
            zibot,zitop,glon0(1),glon0(nlonp4),glat(1),glat(nlat), &
            full(:,:,:,ibnd),ismember(bndry(ibnd),zlog))
          flds(1)%lon_b(:,i,:,nstep_ng(1),ibnd) = bndx(:,1,:)
        enddo
      endif
    enddo

! repeat for lat boundaries
    do lat = 3,4
      do itask = 0,ntask-1
        call mpi_isend(sendbuf,cnt,mpi_real8,send_bndry(lat,itask), &
          0,mpi_comm_world,request(itask))
      enddo

      do itask = 0,ntask-1
        call mpi_irecv(recvbuf(:,:,:,:,itask),cnt,mpi_real8,recv_bndry(lat,itask), &
          0,mpi_comm_world,request(itask+ntask))
      enddo

      call mpi_waitall(ntask*2,request,mpi_statuses_ignore)

      if (flds(1)%is_bndry(lat)) then
        do itask = 0,ntask-1
          i0 = tasks(itask)%lon0
          i1 = tasks(itask)%lon1
          j0 = tasks(itask)%lat0
          j1 = tasks(itask)%lat1
          full(:,i0:i1,j0:j1,:) = recvbuf(:,1:i1-i0+1,1:j1-j0+1,:,itask)
        enddo

        if (lat == 3) then
          bnd(1) = glat_ng(1,-1)
        else
          bnd(1) = glat_ng(1,nlat_ng(1)+2)
        endif

        x0 = glon0(1)
        x1 = glon0(nlonp4)
        nx = nlonp4

        if (.not. (glon0(1)<=glon_ng(1,flds(1)%lond0) .and. glon_ng(1,flds(1)%lond1)<=glon0(nlonp4))) then
          mid_lon = (glon_ng(1,flds(1)%lond0)+glon_ng(1,flds(1)%lond1))/2
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
          tmp = full
          full(:,1:nlonp2+1-sidx,:,:) = tmp(:,sidx:nlonp2,:,:)
          full(:,nlonp4-sidx:nlonp1,:,:) = tmp(:,3:sidx,:,:)
          nx = nlonp1
        endif

        do ibnd = 1,nbnd
          bndy = interp3d(zpint_ng(1,1:nlevp1_ng(1)), &
            glon_ng(1,flds(1)%lond0:flds(1)%lond1),bnd, &
            zibot,zitop,x0,x1,glat(1),glat(nlat), &
            full(:,1:nx,:,ibnd),ismember(bndry(ibnd),zlog))
          flds(1)%lat_b(:,:,lat,nstep_ng(1),ibnd) = bndy(:,:,1)
        enddo
      endif
    enddo

  end subroutine set_bndry
!-----------------------------------------------------------------------
end module level_outer_ng_module

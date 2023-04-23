module sync_ng_module

  use params_module,only: mx_ng
  implicit none

  integer,dimension(mx_ng) :: left,right,above,below

  contains
!-----------------------------------------------------------------------
  subroutine init_sync

    use params_module,only: n_ng
    use fields_ng_module,only: flds,domain
    use mpi_module,only: ntask
    use mpi_f08,only: mpi_proc_null

    integer :: i_ng,lon0,lon1,lat0,lat1,itask

    left = mpi_proc_null
    right = mpi_proc_null
    above = mpi_proc_null
    below = mpi_proc_null

    do i_ng = 1,n_ng
      lon0 = flds(i_ng)%lon0
      lon1 = flds(i_ng)%lon1
      lat0 = flds(i_ng)%lat0
      lat1 = flds(i_ng)%lat1

      do itask = 0,ntask-1
        if (domain(i_ng,2,itask)==lon0-1 .and. domain(i_ng,3,itask)==lat0 .and. domain(i_ng,4,itask)==lat1) &
          left(i_ng) = itask
        if (domain(i_ng,1,itask)==lon1+1 .and. domain(i_ng,3,itask)==lat0 .and. domain(i_ng,4,itask)==lat1) &
          right(i_ng) = itask
        if (domain(i_ng,1,itask)==lon0 .and. domain(i_ng,2,itask)==lon1 .and. domain(i_ng,4,itask)==lat0-1) &
          above(i_ng) = itask
        if (domain(i_ng,1,itask)==lon0 .and. domain(i_ng,2,itask)==lon1 .and. domain(i_ng,3,itask)==lat1+1) &
          below(i_ng) = itask
      enddo
    enddo

  end subroutine init_sync
!-----------------------------------------------------------------------
  subroutine sync_var2d_lon(var,i_ng)

    use fields_ng_module,only: flds,maxlat
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall,mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: lon0,lon1,ny,cnt
    real,dimension(2,maxlat(i_ng)) :: send_left,send_right,recv_left,recv_right
    type(mpi_request),dimension(4) :: request

    lon0 = flds(i_ng)%lon0
    lon1 = flds(i_ng)%lon1
    ny = flds(i_ng)%latd1-flds(i_ng)%latd0+1

! load to work array
    send_left = 0.
    send_right = 0.
    send_left(:,1:ny) = var(lon0:lon0+1,:)
    send_right(:,1:ny) = var(lon1-1:lon1,:)
    recv_left = 0.
    recv_right = 0.

! sync in longitude
    cnt = 2*maxlat(i_ng)

    call mpi_isend(send_left,cnt,mpi_real8,left(i_ng),1,mpi_comm_world,request(1))
    call mpi_isend(send_right,cnt,mpi_real8,right(i_ng),2,mpi_comm_world,request(2))

    call mpi_irecv(recv_right,cnt,mpi_real8,right(i_ng),1,mpi_comm_world,request(3))
    call mpi_irecv(recv_left,cnt,mpi_real8,left(i_ng),2,mpi_comm_world,request(4))

! wait for sync complete
    call mpi_waitall(4,request,mpi_statuses_ignore)

! unpack to model fields
    if (.not. flds(i_ng)%is_bndry(1)) var(lon0-2:lon0-1,:) = recv_left(:,1:ny)
    if (.not. flds(i_ng)%is_bndry(2)) var(lon1+1:lon1+2,:) = recv_right(:,1:ny)

  end subroutine sync_var2d_lon
!-----------------------------------------------------------------------
  subroutine sync_var3d_lon(var,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds,maxlat
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall,mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: lon0,lon1,ny,cnt
    real,dimension(nlevp1_ng(i_ng),2,maxlat(i_ng)) :: send_left,send_right,recv_left,recv_right
    type(mpi_request),dimension(4) :: request

    lon0 = flds(i_ng)%lon0
    lon1 = flds(i_ng)%lon1
    ny = flds(i_ng)%latd1-flds(i_ng)%latd0+1

    send_left = 0.
    send_right = 0.
    send_left(:,:,1:ny) = var(:,lon0:lon0+1,:)
    send_right(:,:,1:ny) = var(:,lon1-1:lon1,:)
    recv_left = 0.
    recv_right = 0.

    cnt = nlevp1_ng(i_ng)*2*maxlat(i_ng)

    call mpi_isend(send_left,cnt,mpi_real8,left(i_ng),1,mpi_comm_world,request(1))
    call mpi_isend(send_right,cnt,mpi_real8,right(i_ng),2,mpi_comm_world,request(2))

    call mpi_irecv(recv_right,cnt,mpi_real8,right(i_ng),1,mpi_comm_world,request(3))
    call mpi_irecv(recv_left,cnt,mpi_real8,left(i_ng),2,mpi_comm_world,request(4))

    call mpi_waitall(4,request,mpi_statuses_ignore)

    if (.not. flds(i_ng)%is_bndry(1)) var(:,lon0-2:lon0-1,:) = recv_left(:,:,1:ny)
    if (.not. flds(i_ng)%is_bndry(2)) var(:,lon1+1:lon1+2,:) = recv_right(:,:,1:ny)

  end subroutine sync_var3d_lon
!-----------------------------------------------------------------------
  subroutine sync_var2d_lat(var,i_ng)

    use fields_ng_module,only: flds,maxlon
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall,mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: nx,lat0,lat1,cnt
    real,dimension(maxlon(i_ng),2) :: send_above,send_below,recv_above,recv_below
    type(mpi_request),dimension(4) :: request

    nx = flds(i_ng)%lond1-flds(i_ng)%lond0+1
    lat0 = flds(i_ng)%lat0
    lat1 = flds(i_ng)%lat1

! load to work array
    send_above = 0.
    send_below = 0.
    send_above(1:nx,:) = var(:,lat0:lat0+1)
    send_below(1:nx,:) = var(:,lat1-1:lat1)
    recv_above = 0.
    recv_below = 0.

! sync in latitude
    cnt = maxlon(i_ng)*2

    call mpi_isend(send_above,cnt,mpi_real8,above(i_ng),1,mpi_comm_world,request(1))
    call mpi_isend(send_below,cnt,mpi_real8,below(i_ng),2,mpi_comm_world,request(2))

    call mpi_irecv(recv_below,cnt,mpi_real8,below(i_ng),1,mpi_comm_world,request(3))
    call mpi_irecv(recv_above,cnt,mpi_real8,above(i_ng),2,mpi_comm_world,request(4))

! wait for sync complete
    call mpi_waitall(4,request,mpi_statuses_ignore)

! unpack to model fields
    if (.not. flds(i_ng)%is_bndry(3)) var(:,lat0-2:lat0-1) = recv_above(1:nx,:)
    if (.not. flds(i_ng)%is_bndry(4)) var(:,lat1+1:lat1+2) = recv_below(1:nx,:)

  end subroutine sync_var2d_lat
!-----------------------------------------------------------------------
  subroutine sync_var3d_lat(var,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds,maxlon
    use mpi_f08,only: mpi_request,mpi_isend,mpi_irecv,mpi_waitall,mpi_real8,mpi_comm_world,mpi_statuses_ignore

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: nx,lat0,lat1,cnt
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),2) :: send_above,send_below,recv_above,recv_below
    type(mpi_request),dimension(4) :: request

    nx = flds(i_ng)%lond1-flds(i_ng)%lond0+1
    lat0 = flds(i_ng)%lat0
    lat1 = flds(i_ng)%lat1

    send_above = 0.
    send_below = 0.
    send_above(:,1:nx,:) = var(:,:,lat0:lat0+1)
    send_below(:,1:nx,:) = var(:,:,lat1-1:lat1)
    recv_above = 0.
    recv_below = 0.

    cnt = nlevp1_ng(i_ng)*maxlon(i_ng)*2

    call mpi_isend(send_above,cnt,mpi_real8,above(i_ng),1,mpi_comm_world,request(1))
    call mpi_isend(send_below,cnt,mpi_real8,below(i_ng),2,mpi_comm_world,request(2))

    call mpi_irecv(recv_below,cnt,mpi_real8,below(i_ng),1,mpi_comm_world,request(3))
    call mpi_irecv(recv_above,cnt,mpi_real8,above(i_ng),2,mpi_comm_world,request(4))

    call mpi_waitall(4,request,mpi_statuses_ignore)

    if (.not. flds(i_ng)%is_bndry(3)) var(:,:,lat0-2:lat0-1) = recv_above(:,1:nx,:)
    if (.not. flds(i_ng)%is_bndry(4)) var(:,:,lat1+1:lat1+2) = recv_below(:,1:nx,:)

  end subroutine sync_var3d_lat
!-----------------------------------------------------------------------
end module sync_ng_module

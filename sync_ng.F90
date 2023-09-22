module sync_ng_module
! sync ghost points, this module is called whenever a procedure violating the consistency is executed
! some variables need not sync, but I am enforcing it for easier bug detecting

  use params_module,only: mx_ng
  implicit none

  integer,dimension(mx_ng) :: left,right,above,below

  contains
!-----------------------------------------------------------------------
  subroutine init_sync

    use params_module,only: n_ng
    use fields_ng_module,only: flds,domain
    use mpi_module,only: ntask
    use mpi

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

! the following assignment will happen exactly once
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
    use mpi_module,only: TIEGCM_WORLD
    use mpi

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: lon0,lon1,ny,cnt,ierror
    real,dimension(2,maxlat(i_ng)) :: send_left,send_right,recv_left,recv_right
    integer,dimension(4) :: request
    external :: shutdown

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

    call mpi_isend(send_left,cnt,mpi_real8,left(i_ng),1,TIEGCM_WORLD,request(1),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 2d fields to left task')

    call mpi_isend(send_right,cnt,mpi_real8,right(i_ng),2,TIEGCM_WORLD,request(2),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 2d fields to right task')

    call mpi_irecv(recv_right,cnt,mpi_real8,right(i_ng),1,TIEGCM_WORLD,request(3),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 2d fields from right task')

    call mpi_irecv(recv_left,cnt,mpi_real8,left(i_ng),2,TIEGCM_WORLD,request(4),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 2d fields from left task')

! wait for sync complete
    call mpi_waitall(4,request,mpi_statuses_ignore,ierror)
    if (ierror /= mpi_success) call shutdown('failed to wait for all mpi tasks to complete')

! unpack to model fields
    if (.not. flds(i_ng)%is_bndry(1)) var(lon0-2:lon0-1,:) = recv_left(:,1:ny)
    if (.not. flds(i_ng)%is_bndry(2)) var(lon1+1:lon1+2,:) = recv_right(:,1:ny)

  end subroutine sync_var2d_lon
!-----------------------------------------------------------------------
  subroutine sync_var3d_lon(var,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds,maxlat
    use mpi_module,only: TIEGCM_WORLD
    use mpi

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: lon0,lon1,ny,cnt,ierror
    real,dimension(nlevp1_ng(i_ng),2,maxlat(i_ng)) :: send_left,send_right,recv_left,recv_right
    integer,dimension(4) :: request
    external :: shutdown

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

    call mpi_isend(send_left,cnt,mpi_real8,left(i_ng),1,TIEGCM_WORLD,request(1),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 3d fields to left task')

    call mpi_isend(send_right,cnt,mpi_real8,right(i_ng),2,TIEGCM_WORLD,request(2),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 3d fields to right task')

    call mpi_irecv(recv_right,cnt,mpi_real8,right(i_ng),1,TIEGCM_WORLD,request(3),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 3d fields from right task')

    call mpi_irecv(recv_left,cnt,mpi_real8,left(i_ng),2,TIEGCM_WORLD,request(4),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 3d fields from left task')

    call mpi_waitall(4,request,mpi_statuses_ignore,ierror)
    if (ierror /= mpi_success) call shutdown('failed to wait for all mpi tasks to complete')

    if (.not. flds(i_ng)%is_bndry(1)) var(:,lon0-2:lon0-1,:) = recv_left(:,:,1:ny)
    if (.not. flds(i_ng)%is_bndry(2)) var(:,lon1+1:lon1+2,:) = recv_right(:,:,1:ny)

  end subroutine sync_var3d_lon
!-----------------------------------------------------------------------
  subroutine sync_var2d_lat(var,i_ng)

    use fields_ng_module,only: flds,maxlon
    use mpi_module,only: TIEGCM_WORLD
    use mpi

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: nx,lat0,lat1,cnt,ierror
    real,dimension(maxlon(i_ng),2) :: send_above,send_below,recv_above,recv_below
    integer,dimension(4) :: request
    external :: shutdown

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

    call mpi_isend(send_above,cnt,mpi_real8,above(i_ng),1,TIEGCM_WORLD,request(1),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 2d fields to above task')

    call mpi_isend(send_below,cnt,mpi_real8,below(i_ng),2,TIEGCM_WORLD,request(2),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 2d fields to below task')

    call mpi_irecv(recv_below,cnt,mpi_real8,below(i_ng),1,TIEGCM_WORLD,request(3),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 2d fields from below task')

    call mpi_irecv(recv_above,cnt,mpi_real8,above(i_ng),2,TIEGCM_WORLD,request(4),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 2d fields from above task')

! wait for sync complete
    call mpi_waitall(4,request,mpi_statuses_ignore,ierror)
    if (ierror /= mpi_success) call shutdown('failed to wait for all mpi tasks to complete')

! unpack to model fields
    if (.not. flds(i_ng)%is_bndry(3)) var(:,lat0-2:lat0-1) = recv_above(1:nx,:)
    if (.not. flds(i_ng)%is_bndry(4)) var(:,lat1+1:lat1+2) = recv_below(1:nx,:)

  end subroutine sync_var2d_lat
!-----------------------------------------------------------------------
  subroutine sync_var3d_lat(var,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds,maxlon
    use mpi_module,only: TIEGCM_WORLD
    use mpi

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: var

    integer :: nx,lat0,lat1,cnt,ierror
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),2) :: send_above,send_below,recv_above,recv_below
    integer,dimension(4) :: request
    external :: shutdown

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

    call mpi_isend(send_above,cnt,mpi_real8,above(i_ng),1,TIEGCM_WORLD,request(1),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 3d fields to above task')

    call mpi_isend(send_below,cnt,mpi_real8,below(i_ng),2,TIEGCM_WORLD,request(2),ierror)
    if (ierror /= mpi_success) call shutdown('failed to send 3d fields to below task')

    call mpi_irecv(recv_below,cnt,mpi_real8,below(i_ng),1,TIEGCM_WORLD,request(3),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 3d fields from below task')

    call mpi_irecv(recv_above,cnt,mpi_real8,above(i_ng),2,TIEGCM_WORLD,request(4),ierror)
    if (ierror /= mpi_success) call shutdown('failed to receive 3d fields from above task')

    call mpi_waitall(4,request,mpi_statuses_ignore,ierror)
    if (ierror /= mpi_success) call shutdown('failed to wait for all mpi tasks to complete')

    if (.not. flds(i_ng)%is_bndry(3)) var(:,:,lat0-2:lat0-1) = recv_above(:,1:nx,:)
    if (.not. flds(i_ng)%is_bndry(4)) var(:,:,lat1+1:lat1+2) = recv_below(:,1:nx,:)

  end subroutine sync_var3d_lat
!-----------------------------------------------------------------------
end module sync_ng_module

module gather2root_ng_module
! wrapper of mpi_gather to collect nested grid fields to a designated root process
! currently only used in hdif

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine gather2root_var2d(var,full,root,i_ng)

    use params_module,only: nlon_ng,nlat_ng
    use fields_ng_module,only: flds,maxlon,maxlat,domain
    use mpi_module,only: ntask,mytid,TIEGCM_WORLD
    use mpi

    integer,intent(in) :: root,i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    real,dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2),intent(out) :: full

    integer :: lon0,lon1,lat0,lat1,cnt,itask,ierror
    real,dimension(maxlon(i_ng),maxlat(i_ng)) :: sendbuf
    real,dimension(maxlon(i_ng),maxlat(i_ng),0:ntask-1) :: recvbuf
    external :: bndry_index_ng,shutdown

    call bndry_index_ng((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
      nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    sendbuf(1:lon1-lon0+1,1:lat1-lat0+1) = var(lon0:lon1,lat0:lat1)
    recvbuf = 0.
    cnt = maxlon(i_ng)*maxlat(i_ng)

    call mpi_gather(sendbuf,cnt,mpi_real8,recvbuf,cnt,mpi_real8,root,TIEGCM_WORLD,ierror)
    if (ierror /= mpi_success) call shutdown('failed to gather one 2d field from other processes')

    if (mytid == root) then
      do itask = 0,ntask-1
        call bndry_index_ng(domain(i_ng,:,itask),nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)
        full(lon0:lon1,lat0:lat1) = recvbuf(1:lon1-lon0+1,1:lat1-lat0+1,itask)
      enddo
    endif

  end subroutine gather2root_var2d
!-----------------------------------------------------------------------
  subroutine gather2root_vars2d(vars,full,nfld,root,i_ng)

    use params_module,only: nlon_ng,nlat_ng
    use fields_ng_module,only: flds,maxlon,maxlat,domain
    use mpi_module,only: ntask,mytid,TIEGCM_WORLD
    use mpi

    integer,intent(in) :: nfld,root,i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,nfld),intent(in) :: vars
    real,dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,nfld),intent(out) :: full

    integer :: lon0,lon1,lat0,lat1,cnt,itask,ierror
    real,dimension(maxlon(i_ng),maxlat(i_ng),nfld) :: sendbuf
    real,dimension(maxlon(i_ng),maxlat(i_ng),nfld,0:ntask-1) :: recvbuf
    external :: bndry_index_ng,shutdown

    call bndry_index_ng((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
      nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    sendbuf(1:lon1-lon0+1,1:lat1-lat0+1,:) = vars(lon0:lon1,lat0:lat1,:)
    recvbuf = 0.
    cnt = maxlon(i_ng)*maxlat(i_ng)*nfld

    call mpi_gather(sendbuf,cnt,mpi_real8,recvbuf,cnt,mpi_real8,root,TIEGCM_WORLD,ierror)
    if (ierror /= mpi_success) call shutdown('failed to gather one 3d field from other processes')

    if (mytid == root) then
      do itask = 0,ntask-1
        call bndry_index_ng(domain(i_ng,:,itask),nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)
        full(lon0:lon1,lat0:lat1,:) = recvbuf(1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
      enddo
    endif

  end subroutine gather2root_vars2d
!-----------------------------------------------------------------------
  subroutine gather2root_var3d(var,full,root,i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
    use fields_ng_module,only: flds,maxlon,maxlat,domain
    use mpi_module,only: ntask,mytid,TIEGCM_WORLD
    use mpi

    integer,intent(in) :: root,i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2),intent(out) :: full

    integer :: lon0,lon1,lat0,lat1,cnt,itask,ierror
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng)) :: sendbuf
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng),0:ntask-1) :: recvbuf
    external :: bndry_index_ng,shutdown

    call bndry_index_ng((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
      nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    sendbuf(:,1:lon1-lon0+1,1:lat1-lat0+1) = var(:,lon0:lon1,lat0:lat1)
    recvbuf = 0.
    cnt = nlevp1_ng(i_ng)*maxlon(i_ng)*maxlat(i_ng)

    call mpi_gather(sendbuf,cnt,mpi_real8,recvbuf,cnt,mpi_real8,root,TIEGCM_WORLD,ierror)
    if (ierror /= mpi_success) call shutdown('failed to gather 2d fields from other processes')

    if (mytid == root) then
      do itask = 0,ntask-1
        call bndry_index_ng(domain(i_ng,:,itask),nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)
        full(:,lon0:lon1,lat0:lat1) = recvbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,itask)
      enddo
    endif

  end subroutine gather2root_var3d
!-----------------------------------------------------------------------
  subroutine gather2root_vars3d(vars,full,nfld,root,i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
    use fields_ng_module,only: flds,maxlon,maxlat,domain
    use mpi_module,only: ntask,mytid,TIEGCM_WORLD
    use mpi

    integer,intent(in) :: nfld,root,i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,nfld),intent(in) :: vars
    real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,nfld),intent(out) :: full

    integer :: lon0,lon1,lat0,lat1,cnt,itask,ierror
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng),nfld) :: sendbuf
    real,dimension(nlevp1_ng(i_ng),maxlon(i_ng),maxlat(i_ng),nfld,0:ntask-1) :: recvbuf
    external :: bndry_index_ng,shutdown

    call bndry_index_ng((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
      nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)

    sendbuf = 0.
    sendbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,:) = vars(:,lon0:lon1,lat0:lat1,:)
    recvbuf = 0.
    cnt = nlevp1_ng(i_ng)*maxlon(i_ng)*maxlat(i_ng)*nfld

    call mpi_gather(sendbuf,cnt,mpi_real8,recvbuf,cnt,mpi_real8,root,TIEGCM_WORLD,ierror)
    if (ierror /= mpi_success) call shutdown('failed to gather 3d fields from other processes')

    if (mytid == root) then
      do itask = 0,ntask-1
        call bndry_index_ng(domain(i_ng,:,itask),nlon_ng(i_ng),nlat_ng(i_ng),lon0,lon1,lat0,lat1)
        full(:,lon0:lon1,lat0:lat1,:) = recvbuf(:,1:lon1-lon0+1,1:lat1-lat0+1,:,itask)
      enddo
    endif

  end subroutine gather2root_vars3d
!-----------------------------------------------------------------------
end module gather2root_ng_module

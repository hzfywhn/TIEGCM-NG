subroutine swdot_ng(w,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use fields_ng_module,only: flds,itp,dlamda,dphi,dz,expzmid,bndry
  use char_module,only: find_index
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: w

  integer :: nk,latd0,latd1,k,lat,idx_omega
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: racs
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: un,vc,w_divrg,du,dvc
  external :: div_ng

  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vc = flds(i_ng)%vc(:,:,:,itp(i_ng))

  racs = flds(i_ng)%racs

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry

! divrg
  call div_ng(un,du,1,i_ng)
  call div_ng(vc,dvc,-1,i_ng)
  w_divrg = du/dlamda(i_ng)+dvc/dphi(i_ng)
  do lat = latd0,latd1
    w_divrg(:,:,lat) = w_divrg(:,:,lat)*racs(lat)
  enddo

  w(nk,:,:) = w_divrg(nk-1,:,:)

  do k = nk-1,1,-1
    w(k,:,:) = expzmid(i_ng)*(expzmid(i_ng)*w(k+1,:,:)+dz(i_ng)*w_divrg(k,:,:))
  enddo

  idx_omega = find_index('OMEGA',bndry)
  if (is_bndry(1)) w(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_omega)
  if (is_bndry(2)) w(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_omega)
  if (is_bndry(3)) w(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_omega)
  if (is_bndry(4)) w(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_omega)

end subroutine swdot_ng

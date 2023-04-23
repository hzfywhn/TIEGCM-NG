subroutine advec_ng(f,hadvec,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use cons_module,only: re_inv
  use fields_ng_module,only: flds,itp,dlamda,dphi
  use dffm_ng_module,only: df_3d,fm_3d
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: hadvec

  integer :: latd0,latd1,nx,ny,lat
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: racs
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    un,vn,ubarl1,ubarl2,dfdx1,dfdx2,wk1,vbarp1,vbarp2,dfdy1,dfdy2,wk2

  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))
  racs = flds(i_ng)%racs

  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  nx = nlon_ng(i_ng)
  ny = nlat_ng(i_ng)
  is_bndry = flds(i_ng)%is_bndry

  call fm_3d(un,ubarl1,1,i_ng)
  call fm_3d(un,ubarl2,2,i_ng)
  call df_3d(f,dfdx1,1,i_ng)
  call df_3d(f,dfdx2,2,i_ng)

  wk1 = dfdx1*2.*2./3.*ubarl1-dfdx2*4./12.*ubarl2
  if (is_bndry(1)) wk1(:,-1:0,:) = dfdx1(:,-1:0,:)*ubarl1(:,-1:0,:)
  if (is_bndry(2)) wk1(:,nx+1:nx+2,:) = dfdx1(:,nx+1:nx+2,:)*ubarl1(:,nx+1:nx+2,:)
  do lat = latd0,latd1
    wk1(:,:,lat) = wk1(:,:,lat)*racs(lat)
  enddo

  call fm_3d(vn,vbarp1,-1,i_ng)
  call fm_3d(vn,vbarp2,-2,i_ng)
  call df_3d(f,dfdy1,-1,i_ng)
  call df_3d(f,dfdy2,-2,i_ng)

  wk2 = dfdy1*2.*2./3.*vbarp1-dfdy2*4./12.*vbarp2
  if (is_bndry(3)) wk2(:,:,-1:0) = dfdy1(:,:,-1:0)*vbarp1(:,:,-1:0)
  if (is_bndry(4)) wk2(:,:,ny+1:ny+2) = dfdy1(:,:,ny+1:ny+2)*vbarp1(:,:,ny+1:ny+2)

  hadvec = wk1/dlamda(i_ng)+wk2*re_inv/dphi(i_ng)

end subroutine advec_ng

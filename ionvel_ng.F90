subroutine ionvel_ng(ui,vi,wi,Etot,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: re
  use fields_ng_module,only: flds,itc
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: ui,vi,wi,Etot

  integer :: nk,k
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: xb,yb,zb,bmod
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: z,ex,ey,ez,eex,eey,eez
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,2,2) :: rjac

  z = flds(i_ng)%z(:,:,:,itc(i_ng))

  ex = flds(i_ng)%ex
  ey = flds(i_ng)%ey
  ez = flds(i_ng)%ez
  rjac = flds(i_ng)%rjac
  xb = flds(i_ng)%xb
  yb = flds(i_ng)%yb
  zb = flds(i_ng)%zb
  bmod = flds(i_ng)%bmod

  nk = nlevp1_ng(i_ng)

  do k = 1,nk
    eex(k,:,:) = rjac(:,:,1,1)*ex(k,:,:)+rjac(:,:,2,1)*ey(k,:,:)
    eey(k,:,:) = rjac(:,:,1,2)*ex(k,:,:)+rjac(:,:,2,2)*ey(k,:,:)
  enddo
  eex = eex/(re+z)
  eey = eey/(re+z)

  do k = 2,nk-1
    eez(k,:,:) = ez(k,:,:)/(z(k+1,:,:)-z(k-1,:,:))
  enddo
  eez(1,:,:) = 2.*eez(2,:,:)-eez(3,:,:)
  eez(nk,:,:) = 2.*eez(nk-1,:,:)-eez(nk-2,:,:)

  do k = 1,nk
    ui(k,:,:) = -(eey(k,:,:)*zb+eez(k,:,:)*xb)*1.e6/bmod**2
    vi(k,:,:) = (eez(k,:,:)*yb+eex(k,:,:)*zb)*1.e6/bmod**2
    wi(k,:,:) = (eex(k,:,:)*xb-eey(k,:,:)*yb)*1.e6/bmod**2
  enddo

  ui = ui*100.
  vi = vi*100.
  wi = wi*100.

  Etot = sqrt(eex**2+eey**2+eez**2)*1E2

  call addfld(ui,'UI_ExB',i_ng)
  call addfld(vi,'VI_ExB',i_ng)
  call addfld(wi,'WI_ExB',i_ng)

end subroutine ionvel_ng

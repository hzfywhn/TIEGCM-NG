subroutine comp_ar_ng(ar_out,arnm_out,istep,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmass_ar
  use ar_module,only: ar_glbm
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: ar_out,arnm_out

  real,parameter :: alfa_ar = 0.17
  real,dimension(4),parameter :: phi_ar = (/1.042,1.509,0.732,1.176/)
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: ar_ubc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: ar_lbc
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: ar,ar_nm,ar_prod,ar_loss
  external :: minor_ng

  ar = flds(i_ng)%ar(:,:,:,itp(i_ng))
  ar_nm = flds(i_ng)%ar_nm(:,:,:,itp(i_ng))

  ar_lbc(:,:,1) = 0.
  ar_lbc(:,:,2) = 1.
  ar_lbc(:,:,3) = -sqrt(ar_glbm(1)*ar_glbm(2))
  ar_ubc = 0.

  ar_prod = 0.
  ar_loss = 0.

  call minor_ng(ar,ar_nm,ar_out,arnm_out,ar_loss,ar_prod,ar_lbc,ar_ubc,rmass_ar,phi_ar,alfa_ar,'AR',istep,i_ng)

end subroutine comp_ar_ng

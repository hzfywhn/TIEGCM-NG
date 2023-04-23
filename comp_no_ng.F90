subroutine comp_no_ng(no_out,no_nm1_out,istep,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmass_no,rmassinv_o2,rmassinv_n4s,rmassinv_n2d,rmassinv_o1,rmassinv_n2
  use chemrates_module,only: beta2,beta6,rk5
  use no_module,only: phi_no,alfa_no
  use fields_ng_module,only: flds,itp,itc
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    no_out,no_nm1_out

  real,parameter :: nob = 4.e6
  integer :: nk,k
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tlbc,no_ubc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: no_lbc
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,n2,xnmbar,xnmbari,n4s,n2d,o2p,no,no_nm1,beta1,beta3,beta8,beta9,beta9n,beta17,no_prod,no_loss,betai
  external :: minor_ng

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  xnmbar = flds(i_ng)%xnmbarm
  xnmbari = flds(i_ng)%xnmbari
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itc(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  no_nm1 = flds(i_ng)%no_nm(:,:,:,itp(i_ng))

  beta1 = flds(i_ng)%beta1
  beta3 = flds(i_ng)%beta3
  beta8 = flds(i_ng)%beta8
  beta9 = flds(i_ng)%beta9
  beta9n = flds(i_ng)%beta9n
  beta17 = flds(i_ng)%beta17
  tlbc = flds(i_ng)%tlbc

  nk = nlevp1_ng(i_ng)

  no_ubc = 0.
  no_lbc(:,:,1) = 0.
  no_lbc(:,:,2) = 1.
  no_lbc(:,:,3) = -nob*rmass_no/xnmbari(1,:,:)

  no_prod = xnmbar**2*o2*rmassinv_o2*(beta1*n4s*rmassinv_n4s+beta2*n2d*rmassinv_n2d)+ &
    xnmbar**3*beta17*o1*rmassinv_o1*n2*rmassinv_n2*n4s*rmassinv_n4s

  do k = 1,nk-1
    betai(k,:,:) = .5*(beta8(k,:,:)+beta9(k,:,:)+beta9n(k,:,:)+beta8(k+1,:,:)+beta9(k+1,:,:)+beta9n(k+1,:,:))
  enddo
  betai(nk,:,:) = 1.5*(beta8(nk,:,:)+beta9(nk,:,:)+beta9n(nk,:,:))- &
    .5*(beta8(nk-1,:,:)+beta9(nk-1,:,:)+beta9n(nk-1,:,:))

  no_loss = -xnmbar*(beta3*n4s*rmassinv_n4s+beta6*n2d*rmassinv_n2d)-rk5*o2p-betai

  call minor_ng(no,no_nm1,no_out,no_nm1_out,no_loss,no_prod,no_lbc,no_ubc,rmass_no,phi_no,alfa_no,'NO',istep,i_ng)

end subroutine comp_no_ng

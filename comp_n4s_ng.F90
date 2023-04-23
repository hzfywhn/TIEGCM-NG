subroutine comp_n4s_ng(n4s_out,n4s_nm1_out,istep,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: brn2d,rmass_n4s,rmassinv_n2d,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no
  use chemrates_module,only: beta4,beta7,rk4,rk6,rk8
  use n4s_module,only: phi_n4s,alfa_n4s
  use fields_ng_module,only: flds,itp,itc
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    n4s_out,n4s_nm1_out

  integer :: nk,k
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tlbc,n4s_ubc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: n4s_lbc
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,n2,xnmbarm,xnmbari,no,n2d,ne,o2p,op,n2p,nplus,nop,n4s,n4s_nm1, &
    qtef,beta1,beta3,beta5,beta8,beta17,ra1,ra3,rk2,n4s_prod,n4s_loss,qtefi,beta8i,nei1,nei2
  external :: minor_ng

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  xnmbarm = flds(i_ng)%xnmbarm
  xnmbari = flds(i_ng)%xnmbari
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itc(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  n2p = flds(i_ng)%n2p
  nplus = flds(i_ng)%nplus
  nop = flds(i_ng)%nop
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  n4s_nm1 = flds(i_ng)%n4s_nm(:,:,:,itp(i_ng))

  qtef = flds(i_ng)%qtef
  beta1 = flds(i_ng)%beta1
  beta3 = flds(i_ng)%beta3
  beta5 = flds(i_ng)%beta5
  beta8 = flds(i_ng)%beta8
  beta17 = flds(i_ng)%beta17
  ra1 = flds(i_ng)%ra1
  ra3 = flds(i_ng)%ra3
  rk2 = flds(i_ng)%rk2
  tlbc = flds(i_ng)%tlbc

  nk = nlevp1_ng(i_ng)

  n4s_lbc(:,:,1) = 0.
  n4s_lbc(:,:,2) = 1.
  n4s_lbc(:,:,3) = -rmass_n4s/xnmbari(1,:,:)*(qtef(1,:,:)*(1.-brn2d)/xnmbari(1,:,:)+ &
    n2d(1,:,:)*rmassinv_n2d*(beta4*xnmbari(1,:,:)*o1(1,:,:)*rmassinv_o1+beta5(1,:,:)*ne(1,:,:)+beta7)+ &
    beta8(1,:,:)*no(1,:,:)*rmassinv_no)/ &
    (beta1(1,:,:)*o2(1,:,:)*rmassinv_o2+beta3(1,:,:)*no(1,:,:)*rmassinv_no+ &
    xnmbari(1,:,:)*beta17(1,:,:)*o1(1,:,:)*rmassinv_o1*n2(1,:,:)*rmassinv_n2)

  n4s_ubc = 0.

  do k = 1,nk-1
    qtefi(k,:,:) = .5*(qtef(k,:,:)+qtef(k+1,:,:))
    beta8i(k,:,:) = .5*(beta8(k,:,:)+beta8(k+1,:,:))
    nei1(k,:,:) = .5*(ne(k,:,:)+ne(k+1,:,:))
    nei2(k,:,:) = sqrt(ne(k,:,:)*ne(k+1,:,:))
  enddo
  qtef(nk,:,:) = 1.5*qtef(nk,:,:)-.5*qtef(nk-1,:,:)
  beta8i(nk,:,:) = 1.5*beta8(nk,:,:)-.5*beta8(nk-1,:,:)
  nei1(nk,:,:) = 1.5*ne(nk,:,:)-.5*ne(nk-1,:,:)
  nei2(nk,:,:) = sqrt(ne(nk,:,:)**3/ne(nk-1,:,:))

  n4s_prod = qtefi*(1.-brn2d)+ &
    xnmbarm*(n2d*rmassinv_n2d*(xnmbarm*beta4*o1*rmassinv_o1+beta5*nei1+beta7)+beta8i*no*rmassinv_no)+ &
    xnmbarm*(rk2*op*n2*rmassinv_n2+rk6*nplus*o2*rmassinv_o2+rk8*nplus*o1*rmassinv_o1)+ &
    nei2*(ra1*nop*0.15+ra3*n2p*1.1)

  n4s_loss = -xnmbarm*(beta1*o2*rmassinv_o2+beta3*no*rmassinv_no+xnmbarm*beta17*o1*rmassinv_o1*n2*rmassinv_n2)-rk4*o2p

  call minor_ng(n4s,n4s_nm1,n4s_out,n4s_nm1_out,n4s_loss,n4s_prod,n4s_lbc,n4s_ubc,rmass_n4s,phi_n4s,alfa_n4s,'N4S',istep,i_ng)

end subroutine comp_n4s_ng

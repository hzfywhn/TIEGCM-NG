subroutine comp_n2d_ng(n2d,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: brn2d,rmassinv_o2,rmassinv_o1,rmassinv_no,rmass_n2d
  use chemrates_module,only: beta2,beta4,beta6,beta7,rk10
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: n2d

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,no,ne,op,n2p,nop,xnmbar,rk3,ra1,ra3,beta5,qtef,n2d_prod,n2d_loss,qtefi,nei

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  n2p = flds(i_ng)%n2p
  nop = flds(i_ng)%nop
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))

  rk3 = flds(i_ng)%rk3
  ra1 = flds(i_ng)%ra1
  ra3 = flds(i_ng)%ra3
  beta5 = flds(i_ng)%beta5
  qtef = flds(i_ng)%qtef

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    qtefi(k,:,:) = .5*(qtef(k,:,:)+qtef(k+1,:,:))
    nei(k,:,:) = sqrt(ne(k,:,:)*ne(k+1,:,:))
  enddo
  qtefi(nk,:,:) = 1.5*qtef(nk,:,:)-.5*qtef(nk-1,:,:)
  nei(nk,:,:) = sqrt(ne(nk,:,:)**3/ne(nk-1,:,:))

  n2d_prod = qtefi*brn2d+rk3*n2p*xnmbar*o1*rmassinv_o1+(ra1*nop*0.85+ra3*n2p*0.9)*nei
  n2d_loss = xnmbar*(beta2*o2*rmassinv_o2+beta4*o1*rmassinv_o1+beta6*no*rmassinv_no)+beta7+beta5*nei+rk10*op
  n2d = rmass_n2d*n2d_prod/(n2d_loss*xnmbar)

end subroutine comp_n2d_ng

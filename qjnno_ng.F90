subroutine qjnno_ng(qtotal,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: evergs,avo,rmassinv_n4s,rmassinv_no,rmassinv_o2,rmassinv_o1,rmassinv_n2d
  use chemrates_module,only: beta2,beta4,beta6
  use fields_ng_module,only: flds,itp,itc
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: qtotal

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,ne,no,n4s,n2d,xnmbar,beta1,beta3,beta5,deltaq,nei,dq

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itc(i_ng))
  xnmbar = flds(i_ng)%xnmbarm

  beta1 = flds(i_ng)%beta1
  beta3 = flds(i_ng)%beta3
  beta5 = flds(i_ng)%beta5

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    nei(k,:,:) = .5*(ne(k,:,:)+ne(k+1,:,:))
  enddo
  nei(nk,:,:) = 1.5*ne(nk,:,:)-.5*ne(nk-1,:,:)

  deltaq = evergs*avo*xnmbar*(n4s*rmassinv_n4s*(beta1*o2*rmassinv_o2*1.4+beta3*no*rmassinv_no*2.68)+ &
    n2d*rmassinv_n2d*(beta2*o2*rmassinv_o2*1.84+beta4*o1*rmassinv_o1*2.38+beta5*nei*2.38/xnmbar+beta6*no*rmassinv_no*5.63))

  do k = 2,nk-1
    dq(k,:,:) = deltaq(k-1,:,:)*deltaq(k,:,:)
  enddo
  dq = sqrt(max(dq,0.))
  dq(1,:,:) = 1.5*deltaq(1,:,:)-0.5*deltaq(2,:,:)
  dq(nk,:,:) = 1.5*deltaq(nk-1,:,:)-0.5*deltaq(nk-2,:,:)
  qtotal = qtotal+dq

end subroutine qjnno_ng

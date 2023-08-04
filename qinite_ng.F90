subroutine qinite_ng(qo2p,qop,qn2p,qnp,qnop,qtef,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no
  use fields_ng_module,only: fb,b,flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: &
    qo2p,qop,qn2p,qnp,qnop,qtef

  real,dimension(3),parameter :: al = (/1.5E7,1.5E6,1.5E6/)
  real,dimension(3,3),parameter :: &
    sa = reshape((/1.6E-18,0.,0.,22.0E-18,10.24E-18,23.11E-18,16.0E-18,8.40E-18,11.61E-18/),(/3,3/)), &
    si = reshape((/1.0E-18,0.,0.,22.0E-18,10.24E-18,23.11E-18,16.0E-18,8.40E-18,11.61E-18/),(/3,3/))
  integer :: nk,k,l
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,he,n2,no,xnmbari,vo2,vo1,vn2,beta9n,o2i,o1i,hei,n2i,tau,qbo2,qbo1,qbn2,noi

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  xnmbari = flds(i_ng)%xnmbari(:,:,:,itp(i_ng))
  vo2 = flds(i_ng)%vo2
  vo1 = flds(i_ng)%vo1
  vn2 = flds(i_ng)%vn2

  beta9n = flds(i_ng)%beta9n

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    o2i(k+1,:,:) = 0.5*(o2(k,:,:)+o2(k+1,:,:))
    o1i(k+1,:,:) = 0.5*(o1(k,:,:)+o1(k+1,:,:))
    n2i(k+1,:,:) = 0.5*(n2(k,:,:)+n2(k+1,:,:))
  enddo

  o2i(1,:,:) = .5*((b(1,1)+1.)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:)+fb(1))
  o1i(1,:,:) = .5*(b(2,1)*o2(1,:,:)+(b(2,2)+1.)*o1(1,:,:)+b(2,3)*he(1,:,:)+fb(2))
  hei(1,:,:) = .5*(b(3,1)*o2(1,:,:)+b(3,2)*o1(1,:,:)+(b(3,3)+1.)*he(1,:,:)+fb(3))
  n2i(1,:,:) = 1.-o2i(1,:,:)-o1i(1,:,:)-hei(1,:,:)

  qbo2 = 0.
  qbo1 = 0.
  qbn2 = 0.
  do l = 1,3
    tau = exp(-(sa(1,l)*vo2+sa(2,l)*vo1+sa(3,l)*vn2))
    qbo2 = qbo2+al(l)*si(1,l)*o2i*tau*rmassinv_o2
    qbo1 = qbo1+al(l)*si(2,l)*o1i*tau*rmassinv_o1
    qbn2 = qbn2+al(l)*si(3,l)*n2i*tau*rmassinv_n2
  enddo

  qo2p = qo2p+0.67*qbo2*xnmbari
  qop = qop+(0.33*qbo2+qbo1)*xnmbari
  qn2p = qn2p+0.86*qbn2*xnmbari
  qnp = qnp+0.14*qbn2*xnmbari
  qtef = qtef+1.57*0.86*qbn2*xnmbari

  do k = 2,nk
    noi(k,:,:) = .5*(no(k,:,:)+no(k-1,:,:))
  enddo
  noi(1,:,:) = no(1,:,:)

  qnop = qnop+beta9n*noi*xnmbari*rmassinv_no

end subroutine qinite_ng

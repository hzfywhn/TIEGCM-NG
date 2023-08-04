subroutine elden_ng(nplus,n2p,nop,o2p,electrons,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmassinv_o2,rmassinv_n2d,rmassinv_o1,rmassinv_n2,rmassinv_no,rmassinv_n4s
  use chemrates_module,only: rk4,rk5,rk6,rk7,rk8,rk9,rk10,rk16,rk23,rk26
  use fields_ng_module,only: flds,itp,itc
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    nplus,n2p,nop,o2p,electrons

  real,parameter :: err = 1.e-300
  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    xnmbar,op,op_upd,o2,o1,n2,n2d,no,n4s,xiop2p,xiop2d,rk1,rk2,rk3,beta9,ra1,ra2,ra3,qnp,qnop,qo2p,qn2p, &
    a0,a1,a2,a3,a4,a,b,c,d,e,fg,h,root,w1,w2,w3,qnpi,qnopi,beta9i,qo2pi,qn2pi

  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  op_upd = flds(i_ng)%op(:,:,:,itc(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  n2d = flds(i_ng)%n2d(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  xiop2p = flds(i_ng)%xiop2p
  xiop2d = flds(i_ng)%xiop2d

  rk1 = flds(i_ng)%rk1
  rk2 = flds(i_ng)%rk2
  rk3 = flds(i_ng)%rk3
  beta9 = flds(i_ng)%beta9
  ra1 = flds(i_ng)%ra1
  ra2 = flds(i_ng)%ra2
  ra3 = flds(i_ng)%ra3
  qnp = flds(i_ng)%qnp
  qnop = flds(i_ng)%qnop
  qo2p = flds(i_ng)%qo2p
  qn2p = flds(i_ng)%qn2p

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    qnpi(k,:,:) = .5*(qnp(k,:,:)+qnp(k+1,:,:))
    qnopi(k,:,:) = .5*(qnop(k,:,:)+qnop(k+1,:,:))
    beta9i(k,:,:) = .5*(beta9(k,:,:)+beta9(k+1,:,:))
    qo2pi(k,:,:) = .5*(qo2p(k,:,:)+qo2p(k+1,:,:))
    qn2pi(k,:,:) = .5*(qn2p(k,:,:)+qn2p(k+1,:,:))
  enddo
  qnpi(nk,:,:) = 1.5*qnp(nk,:,:)-.5*qnp(nk-1,:,:)
  qnopi(nk,:,:) = 1.5*qnop(nk,:,:)-.5*qnop(nk-1,:,:)
  beta9i(nk,:,:) = 1.5*beta9(nk,:,:)-.5*beta9(nk-1,:,:)
  qo2pi(nk,:,:) = 1.5*qo2p(nk,:,:)-.5*qo2p(nk-1,:,:)
  qn2pi(nk,:,:) = 1.5*qn2p(nk,:,:)-.5*qn2p(nk-1,:,:)

  nplus = (qnpi+rk10*op*n2d*xnmbar*rmassinv_n2d)/(xnmbar*((rk6+rk7)*o2*rmassinv_o2+rk8*o1*rmassinv_o1))

  a = qnopi+xnmbar*(rk2*op_upd*n2*rmassinv_n2+rk7*nplus*o2*rmassinv_o2+beta9i*no*rmassinv_no)
  b = qo2pi+xnmbar*(rk1*op_upd+rk6*nplus)*o2*rmassinv_o2+rk26*xiop2d*o2*rmassinv_o2
  c = xnmbar*(rk4*n4s*rmassinv_n4s+rk5*no*rmassinv_no)
  d = qn2pi+(rk16*xiop2p+rk23*xiop2d)*n2*rmassinv_n2
  e = xnmbar*(rk3*o1*rmassinv_o1+rk9*o2*rmassinv_o2)
  fg = op_upd+nplus
  h = xnmbar*rk9*o2*rmassinv_o2

  a0 = -e*c*(a+b+d)
  a1 = -(ra1*(e*(c*fg+b)+d*(c+h))+ra2*(e*(a+d)-h*d)+ra3*c*(a+b))/4.
  a2 = (ra1*(e*c-(ra2*e+ra3*c)*fg-ra2*d-ra3*b)-ra2*ra3*a)/6.
  a3 = ra1*(ra2*e+ra3*c-ra2*ra3*fg)/4.
  a4 = ra1*ra2*ra3

! vquart
  w1 = -(a4*a0-4.*a3*a1+3.*a2**2)/12.
  w2 = (a4*(a2*a0-a1**2)-a3*(a3*a0-a1*a2)+a2*(a3*a1-a2**2))/4.
  root = -2.*real((.5*(w2+sqrt(cmplx(w2**2+4.*w1**3+err)))+err)**(1./3.))
  w1 = sqrt(max(a4*root+a3**2-a4*a2+err,0.))
  w2 = sqrt((2.*root+a2)**2-a4*a0+err)
  w3 = 2.*a3*root+a3*a2-a4*a1+err
  w1 = sign(w1,w2*w3)
  w3 = w1-a3
  root = (w3+sqrt(w3**2-a4*(a2+2.*root-w2)))/a4

  root = max(root,1.0)
  n2p = d/(e+ra3*root)
  o2p = (b+h*d/(e+ra3*root))/(c+ra2*root)
  nop = (a+d*(e-h)/(e+ra3*root)+c*(b+h*d/(e+ra3*root))/(c+ra2*root))/(ra1*root)

  do k = 1,nk-2
    electrons(k+1,:,:) = sqrt(root(k,:,:)*root(k+1,:,:))
  enddo
  electrons(1,:,:) = sqrt(root(1,:,:)**3/root(2,:,:))
  electrons(nk,:,:) = sqrt(root(nk-1,:,:)**3/root(nk-2,:,:))

end subroutine elden_ng

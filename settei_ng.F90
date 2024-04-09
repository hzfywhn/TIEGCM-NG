subroutine settei_ng(te_out,ti_out,qtotal,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: pi,rtd,evergs,rmassinv_o2,rmassinv_o1, &
    rmassinv_he,rmassinv_n2,avo,rmass_o1,rmassinv_n4s,rmassinv_no
  use input_module,only: f107,et,electron_heating
  use lbc,only: fb,b
  use fields_ng_module,only: flds,itp,dipmin,dz
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    te_out,ti_out,qtotal

  real,parameter :: fpolar = -3.0e+9, del = 1.e-6, alam = 0.0069, ad = 0.0091, sd = 2.3e-11
  integer :: nk,k
  real :: f107te
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    t_lbc,chi,rlatm,dipmag,qteaur,a,fed,fen,fe,sindipmag
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,he,n2,ne,te,ti,op,o2p,nplus,n2p,nop,xnmbar,xnmbari,scht,schti,qji_ti,qop2p,qop2d,qo2p,qop,qn2p,qnp,qnop, &
    te_int,tn_int,o2n,o1n,hen,n2n,root_te,root_tn,root_ne,tek0,h_mid,h_int,p_coef,q_coef,r_coef,rhs, &
    qtot,qe,q_eni,coll_en2v,loss_en2v,loss_eo2,loss_eo1d,loss_eo1,loss_ehe,loss_en2,loss_en,loss_xen, &
    loss_ei,loss_in,tek0_h_int,tek0_h_int_1,q_eni_i,Q1,Q2
  external :: trsolv_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  nplus = flds(i_ng)%nplus
  n2p = flds(i_ng)%n2p
  nop = flds(i_ng)%nop
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  xnmbari = flds(i_ng)%xnmbari(:,:,:,itp(i_ng))
  scht = flds(i_ng)%scht(:,:,:,itp(i_ng))
  schti = flds(i_ng)%schti(:,:,:,itp(i_ng))
  qji_ti = flds(i_ng)%qji_ti

  t_lbc = flds(i_ng)%t_lbc
  chi = flds(i_ng)%chi
  rlatm = flds(i_ng)%rlatm
  dipmag = flds(i_ng)%dipmag
  qteaur = flds(i_ng)%qteaur
  qop2p = flds(i_ng)%qop2p
  qop2d = flds(i_ng)%qop2d
  qo2p = flds(i_ng)%qo2p
  qop = flds(i_ng)%qop
  qn2p = flds(i_ng)%qn2p
  qnp = flds(i_ng)%qnp
  qnop = flds(i_ng)%qnop
  Q1 = flds(i_ng)%Q1
  Q2 = flds(i_ng)%Q2

  nk = nlevp1_ng(i_ng)

  f107te = min(f107,235.)

  where (abs(rlatm) >= pi/4.5) a = 1.
! where (abs(rlatm) < pi/4.5) a = .5*(1.+sin(pi*(abs(rlatm)-pi/9.)/(pi/4.5)))
  where (abs(rlatm)<pi/4.5 .and. abs(rlatm)>pi/18) &
    a = .5*(1.+cos(abs(rlatm)*6.-pi*4./3.))
  where (abs(rlatm) <= pi/18) a = 0.

! fed = (-5.0e+7*f107te*a-4.0e+7*f107te)*1.2
  fed = -9.0e+7*f107te*a
  fen = fed/2.
  fed = fed+qteaur
  fen = fen+qteaur

  where (chi >= .5*pi)
    fe = fen
  elsewhere
    fe = fed
  endwhere
  where ((chi*rtd-80.)*(chi*rtd-100.) >= 0.)
    fe = fe*evergs
  elsewhere
    fe = (.5*(fed+fen)+.5*(fed-fen)*cos(pi*(chi*rtd-80.)/20.))*evergs
  endwhere
  where (abs(rlatm) >= pi/3.) fe = fe+fpolar*evergs

  do k = 2,nk-1
    te_int(k,:,:) = .5*(te(k,:,:)+te(k-1,:,:))
    o2n(k,:,:) = .5*(o2(k,:,:)+o2(k-1,:,:))
    o1n(k,:,:) = .5*(o1(k,:,:)+o1(k-1,:,:))
    hen(k,:,:) = .5*(he(k,:,:)+he(k-1,:,:))
    n2n(k,:,:) = .5*(n2(k,:,:)+n2(k-1,:,:))
    tn_int(k,:,:) = .5*(tn(k,:,:)+tn(k-1,:,:))
  enddo

  te_int(1,:,:) = max(1.5*te(1,:,:)-.5*te(2,:,:),t_lbc)
  o2n(1,:,:) = .5*(fb(1)+(b(1,1)+1.)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:))
  o1n(1,:,:) = .5*(fb(1)+b(1,1)*o2(1,:,:)+(b(1,2)+1.)*o1(1,:,:)+b(1,3)*he(1,:,:))
  hen(1,:,:) = .5*(fb(1)+b(1,1)*o2(1,:,:)+b(1,2)*o1(1,:,:)+(b(1,3)+1.)*he(1,:,:))
  n2n(1,:,:) = max(1.-o2n(1,:,:)-o1n(1,:,:)-hen(1,:,:),0.)
  tn_int(1,:,:) = t_lbc

  te_int(nk,:,:) = 1.5*te(nk-1,:,:)-.5*te(nk-2,:,:)
  o2n(nk,:,:) = .5*(o2(nk,:,:)+o2(nk-1,:,:))
  o1n(nk,:,:) = .5*(o1(nk,:,:)+o1(nk-1,:,:))
  hen(nk,:,:) = .5*(he(nk,:,:)+he(nk-1,:,:))
  n2n(nk,:,:) = .5*(n2(nk,:,:)+n2(nk-1,:,:))
  tn_int(nk,:,:) = 1.5*tn(nk-1,:,:)-.5*tn(nk-2,:,:)

  n2n = max(n2n,0.)

  o2n = xnmbari*o2n*rmassinv_o2
  o1n = xnmbari*o1n*rmassinv_o1
  hen = xnmbari*hen*rmassinv_he
  n2n = xnmbari*n2n*rmassinv_n2

  root_te = sqrt(te_int)
  tek0 = 7.5e5/(1.+3.22e4*te_int**2/ne* &
    ((2.20e-16 + 7.92e-18*root_te)*o2n+ &
    1.10e-16*(1.+5.7e-4*te_int)*o1n+ &
    5.60e-16*hen+ &
    (2.82e-17 - 3.41e-21*te_int)*root_te*n2n))* &
    evergs

  h_mid = scht
  h_int = schti

  sindipmag = max(sin(max(abs(dipmag),dipmin(i_ng)))**2,.10)

  p_coef = 2./7./(h_mid*dz(i_ng)**2)
  do k = 1,nk
    p_coef(k,:,:) = p_coef(k,:,:)*sindipmag
  enddo

  tek0_h_int = tek0/h_int
  do k = 1,nk-1
    tek0_h_int_1(k,:,:) = tek0_h_int(k+1,:,:)
  enddo
  tek0_h_int_1(nk,:,:) = 2.*tek0_h_int(nk,:,:)-tek0_h_int(nk-1,:,:)

  r_coef = p_coef*tek0_h_int_1
  p_coef = p_coef*tek0_h_int
  q_coef = -(p_coef+r_coef)
  rhs = 0.

  q_coef(1,:,:) = q_coef(1,:,:)-p_coef(1,:,:)
  rhs(1,:,:) = rhs(1,:,:)-2.*p_coef(1,:,:)*tn_int(1,:,:)**3.5
  p_coef(1,:,:) = 0.

  q_coef(nk-1,:,:) = q_coef(nk-1,:,:)+r_coef(nk-1,:,:)
  rhs(nk-1,:,:) = rhs(nk-1,:,:)+r_coef(nk-1,:,:)*dz(i_ng)*3.5*h_int(nk,:,:)*fe/tek0(nk,:,:)
  r_coef(nk-1,:,:) = 0.

  qtot = max(qo2p+qop+qn2p+qnop+qnp+qop2d+qop2p,1.e-20)

  do k = 1,nk-1
    qtot(k,:,:) = sqrt(qtot(k,:,:)*qtot(k+1,:,:))
  enddo
  qtot(nk,:,:) = 0.

  do k = 1,nk-1
    root_ne(k,:,:) = ne(k,:,:)*ne(k+1,:,:)
  enddo
  root_ne(nk,:,:) = ne(nk,:,:)**3/ne(nk-1,:,:)
  root_ne = sqrt(max(root_ne,1.e4))

  o2n = xnmbar*o2*rmassinv_o2
  o1n = xnmbar*o1*rmassinv_o1
  hen = xnmbar*he*rmassinv_he
  n2n = xnmbar*n2*rmassinv_n2

  if (electron_heating == 6) then
    qe = log(root_ne/(o2n+n2n+o1n))
    qe = exp(((((((-1.249e-5*qe-5.755e-4)*qe-9.346e-3)*qe-5.900e-2)*qe-4.392e-2)*qe+1.056)*qe+5.342))
  else
    qe = log(root_ne/(o2n+n2n+0.1*o1n))
    qe = exp(-((((0.001996*qe+0.08034)*qe+1.166)*qe+6.941)*qe+12.75))
  endif

  rhs = rhs-qe*qtot*evergs

  if (et) rhs = rhs-(Q1+Q2)*10.0

  root_te = sqrt(te)

  where (te >= 1000.)
    coll_en2v = 2.e-7*exp(-4605.2/te)
  elsewhere
    coll_en2v = 5.71e-8*exp(-3352.6/te)
  endwhere
  where (te > 2000.) coll_en2v = 2.53e-6*root_te*exp(-17620./te)

  where (abs(te-tn) < del)
    loss_en2v = 3200./tn**2
  elsewhere
    loss_en2v = 1./(te-tn)*(1.-exp(-3200.*(te-tn)/(te*tn)))
  endwhere
  loss_en2v = 1.3e-4*loss_en2v*coll_en2v

  loss_eo2 = o2n*(1.21e-18*(1.+3.6e-2*root_te)*root_te+6.9e-14/root_te+3.125e-21*te**2)

! where (abs(te-tn) < del)
!   loss_eo1d = 22713./tn**2
! elsewhere
!   loss_eo1d = 1./(te-tn)*(1.-exp(-22713.*(te-tn)/(te*tn)))
! endwhere
! loss_eo1d = 1.57e-12*loss_eo1d* &
!   exp((2.4e4+0.3*(te-1500.)-1.947e-5*(te-1500.)*(te-4000.))*(te-3000.)/(3000.*te))
  loss_eo1d = 0.

  loss_eo1 = o1n*(7.9e-19*(1.+5.7e-4*te)*root_te+3.4e-12*(1.-7.e-5*te)/tn*(150./te+0.4))

  loss_ehe = hen*2.46e-17*root_te

  loss_en2 = n2n*(1.77E-19*(1.-1.21E-4*te)*te+2.9e-14/root_te+loss_en2v)

  loss_en = loss_eo2+loss_eo1+loss_ehe+loss_en2

  if (et) then
    where (te>500.0 .and. Q2>0.0) loss_en = loss_en*exp(-7.54E-4*(te-500.0))
  endif

  loss_xen = (loss_en+o1n*(1.-alam/(ad+sd*n2n))*loss_eo1d)*root_ne*evergs

  loss_en = (loss_en+o1n*loss_eo1d)*root_ne*evergs

  loss_ei = 3.2e-8*root_ne/(root_te*te)*15.*evergs*rmass_o1* &
    (op*rmassinv_o1+o2p*rmassinv_o2+nplus*rmassinv_n4s+n2p*rmassinv_n2+nop*rmassinv_no)

  root_tn = sqrt(2*tn)

  loss_in = 1e-14*evergs* &
    (op*(5.8*o2n+0.21*o1n*root_tn+2.8*hen+6.6*n2n)+ &
    o2p*(0.14*o2n*root_tn+4.36*o1n+1.63*hen+5.81*n2n)+ &
    nplus*(5.84*o2n+5.84*o1n+3.05*hen+6.56*n2n)+ &
    n2p*(5.54*o2n+4.65*o1n+1.82*hen+0.27*n2n*root_tn)+ &
    nop*(5.45*o2n+4.5*o1n+1.72*hen+5.92*n2n))

  q_coef = q_coef-(loss_en+loss_ei)/te**2.5

  rhs = rhs-loss_en*tn-loss_ei*ti

  q_eni = loss_ei*max(te-ti,0.)
  q_eni = (loss_xen*(te-tn)+q_eni)*avo/xnmbar

  do k = 1,nk-2
    q_eni_i(k+1,:,:) = .5*(q_eni(k,:,:)+q_eni(k+1,:,:))
  enddo
  q_eni_i(1,:,:) = 1.5*q_eni(1,:,:)-0.5*q_eni(2,:,:)
  q_eni_i(nk,:,:) = 1.5*q_eni(nk-1,:,:)-0.5*q_eni(nk-2,:,:)
  qtotal = qtotal+q_eni_i

  call trsolv_ng(p_coef,q_coef,r_coef,rhs,te_out,nk-1,i_ng)

  te_out = te_out**(2./7.)
  te_out(nk,:,:) = te_out(nk-1,:,:)**2/te_out(nk-2,:,:)
  te_out = max(te_out,tn)

  ti_out = max((qji_ti*xnmbar/avo+loss_ei*te_out+loss_in*tn)/(loss_ei+loss_in),tn)

end subroutine settei_ng

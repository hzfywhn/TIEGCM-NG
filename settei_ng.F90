subroutine settei_ng(te_out,ti_out,qtotal,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: pi,rtd,evergs,rmassinv_o2,rmassinv_o1,rmassinv_n2,gask,grav,avo
  use input_module,only: f107
  use fields_ng_module,only: flds,itp,dipmin,dz
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    te_out,ti_out,qtotal

  real,parameter :: fpolar = -3.0e+9, del = 1.e-6, alam = 0.0069, ad = 0.0091, sd = 2.3e-11
  integer :: nk,k
  real :: f107te
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    chi,rlatm,dipmag,qteaur,a,fed,fen,fe,sindipmag
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,n2,ne,te,ti,op,o2p,nop,mbar,barm,xnmbar,xnmbari,qji_ti,qop2p,qop2d,qo2p,qop,qn2p,qnp,qnop, &
    te_int,tn_int,o2n,o1n,n2n,root_te,root_tn,root_ne,tek0,h_mid,h_int,p_coef,q_coef,r_coef,rhs, &
    qtot,qe,q_eni,coll_en2v,loss_en2v,loss_en2,loss_eo2,loss_eo1d,loss_eo1,loss_xen,loss_en,loss_ei,loss_in, &
    tek0_h_int,tek0_h_int_1,q_eni_i
  external :: trsolv_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  nop = flds(i_ng)%nop
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  barm = flds(i_ng)%barm(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar
  xnmbari = flds(i_ng)%xnmbari
  qji_ti = flds(i_ng)%qji_ti

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

  nk = nlevp1_ng(i_ng)

  f107te = min(f107,235.)

  where (abs(rlatm) >= pi/4.5) a = 1.
  where (abs(rlatm)<pi/4.5 .and. abs(rlatm)<=pi/18) a = 0.
  where (abs(rlatm) > pi/18) a = .5*(1.+sin(pi*(abs(rlatm)-pi/4.5)/(pi/6.)))

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
    n2n(k,:,:) = .5*(n2(k,:,:)+n2(k-1,:,:))
    tn_int(k,:,:) = .5*(tn(k,:,:)+tn(k-1,:,:))
  enddo

  te_int(1,:,:) = 1.5*te(1,:,:)-.5*te(2,:,:)
  o2n(1,:,:) = 1.5*o2(1,:,:)-.5*o2(2,:,:)
  o1n(1,:,:) = 1.5*o1(1,:,:)-.5*o1(2,:,:)
  n2n(1,:,:) = 1.5*n2(1,:,:)-.5*n2(2,:,:)
  tn_int(1,:,:) = 1.5*tn(1,:,:)-.5*tn(2,:,:)

  te_int(nk,:,:) = 1.5*te(nk-1,:,:)-.5*te(nk-2,:,:)
  o2n(nk,:,:) = 1.5*o2(nk-1,:,:)-.5*o2(nk-2,:,:)
  o1n(nk,:,:) = 1.5*o1(nk-1,:,:)-.5*o1(nk-2,:,:)
  n2n(nk,:,:) = 1.5*n2(nk-1,:,:)-.5*n2(nk-2,:,:)
  tn_int(nk,:,:) = 1.5*tn(nk-1,:,:)-.5*tn(nk-2,:,:)

  n2n = max(n2n,0.)

  o2n = xnmbari*o2n*rmassinv_o2
  o1n = xnmbari*o1n*rmassinv_o1
  n2n = xnmbari*n2n*rmassinv_n2

  root_te = sqrt(te_int)
  tek0 = 7.5e5/(1.+3.22e4*te_int**2/ne*(root_te*(2.82e-17-3.41e-21*te_int)*n2n+ &
    (2.20e-16+7.92e-18*root_te)*o2n+1.10e-16*(1.+5.7e-4*te_int)*o1n))*evergs

  h_mid = gask*tn/(mbar*grav)
  h_int = gask*tn_int/(barm*grav)

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

  do k = 1,nk-1
    root_ne(k,:,:) = ne(k,:,:)*ne(k+1,:,:)
  enddo
  root_ne(nk,:,:) = ne(nk,:,:)**3/ne(nk-1,:,:)
  root_ne = sqrt(max(root_ne,1.e4))

  o2n = xnmbar*o2*rmassinv_o2
  o1n = xnmbar*o1*rmassinv_o1
  n2n = xnmbar*n2*rmassinv_n2

!  qe = log(root_ne/(o2n+n2n+0.1*o1n))
!  qe = exp(-((((0.001996*qe+0.08034)*qe+1.166)*qe+6.941)*qe+12.75))
  qe = log(root_ne/(o2n+n2n+o1n))
  qe = exp(-((((((0.00001249*qe+0.0005755)*qe+0.009346)*qe+0.059)*qe+0.04392)*qe-1.056)*qe-5.342))

  rhs = rhs-qe*qtot*evergs
  root_te = sqrt(te)

  where (te >= 1000.)
    coll_en2v = 2.e-7*exp(-4605.2/te)
  elsewhere
    coll_en2v = 5.71e-8*exp(-3352.6/te)
  endwhere
  where (te > 2000.) coll_en2v = 2.53e-6*root_te*exp(-17620./te)

  loss_en2v = 3200.*(1./te-1./tn)
  loss_en2v = sign(abs(loss_en2v)+del,loss_en2v)
  loss_en2v = -3200./(te*tn)*(1.-exp(loss_en2v))/loss_en2v
  loss_en2v = 1.3e-4*loss_en2v*coll_en2v

  loss_en2 = n2n*(1.77E-19*(1.-1.21E-4*te)*te+2.9e-14/root_te+loss_en2v)
  loss_en = loss_en2

  loss_eo2 = o2n*(1.21e-18*(1.+3.6e-2*root_te)*root_te+6.9e-14/root_te+3.125e-21*te**2)
  loss_en = loss_en+loss_eo2

! loss_eo1d = 22713.*(1./te-1./tn)
! loss_eo1d = sign(abs(loss_eo1d)+del,loss_eo1d)
! loss_eo1d = 22713./(te*tn)*(1.-exp(loss_eo1d))/loss_eo1d
  loss_eo1d = 0.

  loss_eo1 = o1n*(7.9e-19*(1.+5.7e-4*te)*root_te+3.4e-12*(1.-7.e-5*te)/tn*(150./te+0.4))

  loss_en = loss_en+loss_eo1

  loss_xen = (loss_en+o1n*(1.-alam/(ad+sd*n2n))*loss_eo1d)*root_ne*evergs

  loss_en = (loss_en+o1n*loss_eo1d)*root_ne*evergs

  loss_ei = 3.2e-8*root_ne/(root_te*te)*15.*(op+0.5*o2p+0.53*nop)*evergs

  root_tn = sqrt(tn)

  loss_in = ((6.6e-14*n2n+5.8e-14*o2n+0.21e-14*o1n*sqrt(2.)*root_tn)*op+ &
    (5.45e-14*o2n+5.9e-14*n2n+4.5e-14*o1n)*nop+ &
    (5.8e-14*n2n+4.4e-14*o1n+0.14e-14*o2n*root_tn)*o2p)*evergs

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

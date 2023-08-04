subroutine qjion_ng(qtotal,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: avo,evergs,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no,rmassinv_n2d,rmassinv_n4s
  use chemrates_module,only: rk4,rk5,rk6,rk7,rk8,rk9,rk10,rk16,rk17,rk18,rk21,rk22,rk23,rk24,rk26,rk27
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: qtotal

  real,parameter :: aureff = 0.05
  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o2,o1,n2,o2p,op,n4s,n2d,no,ne,xnmbar,xnmbari,n2p,nplus,nop,xiop2p,xiop2d, &
    rk1,rk2,rk3,ra1,ra2,ra3,rk19,rk20,rk25,qop2p,qop2d,qo2p,qop,qn2p,qnp,qnop,qtot,qphoto,qic,nei

  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  xnmbari = flds(i_ng)%xnmbari(:,:,:,itp(i_ng))
  n2p = flds(i_ng)%n2p
  nplus = flds(i_ng)%nplus
  nop = flds(i_ng)%nop
  xiop2p = flds(i_ng)%xiop2p
  xiop2d = flds(i_ng)%xiop2d

  rk1 = flds(i_ng)%rk1
  rk2 = flds(i_ng)%rk2
  rk3 = flds(i_ng)%rk3
  ra1 = flds(i_ng)%ra1
  ra2 = flds(i_ng)%ra2
  ra3 = flds(i_ng)%ra3
  rk19 = flds(i_ng)%rk19
  rk20 = flds(i_ng)%rk20
  rk25 = flds(i_ng)%rk25
  qop2p = flds(i_ng)%qop2p
  qop2d = flds(i_ng)%qop2d
  qo2p = flds(i_ng)%qo2p
  qop = flds(i_ng)%qop
  qn2p = flds(i_ng)%qn2p
  qnp = flds(i_ng)%qnp
  qnop = flds(i_ng)%qnop

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    nei(k,:,:) = .5*(ne(k,:,:)+ne(k+1,:,:))
  enddo
  nei(nk,:,:) = 1.5*ne(nk,:,:)-.5*ne(nk-1,:,:)

  qtot = qo2p+qop+qn2p+qnop+qnp+qop2d+qop2p

  qphoto = qtot*aureff*35.*avo*evergs/xnmbari

  qic = (avo*(o2*rmassinv_o2*(rk1*op*1.555+(rk6*2.486+rk7*6.699)*nplus+rk9*n2p*3.52)+ &
    op*(rk2*n2*rmassinv_n2*1.0888+rk10*n2d*rmassinv_n2d*1.45)+ &
    o1*rmassinv_o1*(rk3*n2p*0.70+rk8*nplus*0.98)+ &
    o2p*(rk4*n4s*rmassinv_n4s*4.21+rk5*no*rmassinv_no*2.813))+ &
    nei*(ra1*nop*0.854+ra2*o2p*5.2755+ra3*n2p*3.678)/xnmbar)*evergs+ &
    (avo*(((rk16*3.02+rk17*0.7)*n2*rmassinv_n2+rk18*o1*rmassinv_o1*5.0)*xiop2p+ &
    (rk23*n2*rmassinv_n2*1.33+rk24*o1*rmassinv_o1*3.31+rk26*4.87*o2*rmassinv_o2)*xiop2d)+ &
    (nei*((rk19*5.0+rk20*1.69)*xiop2p+rk25*3.31*xiop2d)-(rk21*5.02+rk22*1.69)*xiop2p-rk27*3.33*xiop2d)/xnmbar)*evergs
  qic = max(qic,1.e-30)

  do k = 1,nk-2
    qtot(k+1,:,:) = sqrt(qic(k,:,:)*qic(k+1,:,:))
  enddo
  qtot(1,:,:) = sqrt(qic(1,:,:)**3/qic(2,:,:))
  qtot(nk,:,:) = sqrt(qic(nk-1,:,:)**3/qic(nk-2,:,:))

  qtot = qphoto+qtot

  qtotal = qtotal+qtot

end subroutine qjion_ng

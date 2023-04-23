subroutine comp_o2o_ng(fs,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmassinv_o2,rmassinv_o1,rmass_no,rmass_n4s,rmass_n2d,rmass_o2,rmass_o1
  use chemrates_module,only: rk4,rk5,rk6,rk7,rk8,rk9,rk10,beta2,beta6
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,0:3),intent(out) :: fs

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    o1,mbar,xnmbar,op,no,n4s,n2d,o2p,ne,n2p,nplus,nop,rj,qo2p,qop,ra1,ra2,beta3,beta8,rk1,beta1, &
    rkm12,rk3,rji,qo2pi,qopi,pox1,pox2,lox1,lox2,lox3,po21,po22,po23,lo21,lo22,beta8i,nei

  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  n4s = flds(i_ng)%n4s(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  n2p = flds(i_ng)%n2p
  nplus = flds(i_ng)%nplus
  nop = flds(i_ng)%nop

  rj = flds(i_ng)%rj
  qo2p = flds(i_ng)%qo2p
  qop = flds(i_ng)%qop
  ra1 = flds(i_ng)%ra1
  ra2 = flds(i_ng)%ra2
  beta3 = flds(i_ng)%beta3
  beta8 = flds(i_ng)%beta8
  rk1 = flds(i_ng)%rk1
  beta1 = flds(i_ng)%beta1
  rkm12 = flds(i_ng)%rkm12
  rk3 = flds(i_ng)%rk3

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    rji(k,:,:) = 0.5*(rj(k,:,:)+rj(k+1,:,:))
    qo2pi(k,:,:) = 0.5*(qo2p(k,:,:)+qo2p(k+1,:,:))
    qopi(k,:,:) = 0.5*(qop(k,:,:)+qop(k+1,:,:))
    beta8i(k,:,:) = 0.5*(beta8(k,:,:)+beta8(k+1,:,:))
    nei(k,:,:) = sqrt(ne(k,:,:)*ne(k+1,:,:))
  enddo
  rji(nk,:,:) = 1.5*rj(nk,:,:)-.5*rj(nk-1,:,:)
  qo2pi(nk,:,:) = 1.5*qo2p(nk,:,:)-.5*qo2p(nk-1,:,:)
  qopi(nk,:,:) = 1.5*qop(nk,:,:)-.5*qop(nk-1,:,:)
  beta8i(nk,:,:) = 1.5*beta8(nk,:,:)-.5*beta8(nk-1,:,:)
  nei(nk,:,:) = sqrt(ne(nk,:,:)**3/ne(nk-1,:,:))

  pox1 = xnmbar**2*(beta3*n4s/rmass_n4s*no/rmass_no+beta6*n2d/rmass_n2d*no/rmass_no)+ &
    beta8i*no/rmass_no*xnmbar+xnmbar*(rk4*o2p*n4s/rmass_n4s+rk10*op*n2d/rmass_n2d)+ &
    (ra1*nop+2.*ra2*o2p)*nei

  pox2 = xnmbar*(beta1*n4s/rmass_n4s+beta2*n2d/rmass_n2d)+rk1*op+rk7*nplus+2.*rji

  lox1 = 2.*rkm12*xnmbar/mbar
  lox2 = rk3*n2p+rk8*nplus
  lox3 = qopi

  po21 = rkm12*xnmbar/mbar
  po22 = 0.
  po23 = rk5*no/rmass_no*o2p*xnmbar

  lo21 = xnmbar*(beta1*n4s/rmass_n4s+beta2*n2d/rmass_n2d)+rk1*op+(rk6+rk7)*nplus+rk9*n2p+rji
  lo22 = qo2pi

  fs(:,:,:,1,1) = -lo21
  fs(:,:,:,1,2) = xnmbar*po21*o1*rmassinv_o1**2*rmass_o2
  fs(:,:,:,2,1) = pox2*rmass_o1*rmassinv_o2
  fs(:,:,:,2,2) = -lox2-lox1*o1*rmassinv_o1*xnmbar
  fs(:,:,:,1,0) = (po23-lo22)*rmass_o2/xnmbar
  fs(:,:,:,2,0) = (pox1-lox3)*rmass_o1/xnmbar
  fs(:,:,:,1:2,3) = 0.
  fs(:,:,:,3,:) = 0.

end subroutine comp_o2o_ng

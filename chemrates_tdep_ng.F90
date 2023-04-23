subroutine chemrates_tdep_ng(rk1,rk2,rk3,ra1,ra2,ra3,beta1,beta3, &
    beta5,beta8,beta9,beta9n,beta17,rk19,rk20,rk25,rkm12,i_ng)

  use params_module,only: nlevp1_ng
  use input_module,only: f107
  use init_module,only: sfeps
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) ::  &
    rk1,rk2,rk3,ra1,ra2,ra3,beta1,beta3,beta5,beta8,beta9,beta9n,beta17,rk19,rk20,rk25,rkm12

  integer :: nk
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,te,ti,fno2,fnvo2,ti1,ti2,ti3,etvib

  nk = nlevp1_ng(i_ng)

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  fno2 = flds(i_ng)%sco2
  fnvo2 = flds(i_ng)%vo2

  ti1 = (0.667*ti+0.333*tn)/300.
  ti2 = (0.6363*ti+0.3637*tn)/300.
  ti3 = .5*(ti+tn)/300.

  rk1 = 1.6E-11*ti1**(-0.52)+5.5E-11*exp(-22.85/ti1)
  rk2 = (8.6E-14*ti2-5.92E-13)*ti2+1.533E-12
  etvib = exp(-3353./tn)
  rk2 = (((((270.*etvib+220.)*etvib+85.)*etvib+38.)*etvib+1.)*rk2*etvib+rk2)*(1.-etvib)
  where (300.*ti3 >= 1500.)
    rk3 = 5.2E-11*ti3**0.2
  elsewhere
    rk3 = 1.4E-10*ti3**(-0.44)
  endwhere
  rk19 = 4.0E-8*sqrt(300./te)
  rk20 = 1.5E-7*sqrt(300./te)
  rk25 = 6.6E-8*sqrt(300./te)
  rkm12 = 9.59E-34*exp(480./tn)

  ra1 = 4.2E-7*(300./te)**0.85
  where (te >= 1200.)
    ra2 = 1.6E-7*(300./te)**0.55
  elsewhere
    ra2 = 2.7E-7*(300./te)**0.7
  endwhere
  ra3 = 1.8E-7*(300./te)**0.39

  beta1 = 1.5E-11*exp(-3600./tn)
  beta3 = 3.4e-11*sqrt(tn/300.)
  beta5 = 3.6E-10*sqrt(te/300.)
  beta8 = 4.5E-6*(1.+0.11*(f107-65.)/165.)*exp(-1.E-8*fno2**0.38)*sfeps
  beta9 = 2.91E11*(1.+0.2*(f107-65.)/100.)*2.E-18*exp(-8.E-21*fno2)*sfeps
  beta9n = 5.E9*(1.+0.2*(f107-65.)/100.)*2.E-18*exp(-8.E-21*fnvo2)*sfeps
  beta17 = 1.E-32*sqrt(300./te)

  beta8(nk,:,:) = 4.5E-6*(1.+0.11*(f107-65.)/165.)*exp(-1.E-8*fno2(nk,:,:)**0.38)
  beta9(nk,:,:) = 2.91E11*(1.+0.2*(f107-65.)/100.)*2.E-18*exp(-8.E-21*fno2(nk,:,:))
  beta9n(nk,:,:) = 5.E9*(1.+0.2*(f107-65.)/100.)*2.E-18*exp(-8.E-21*fnvo2(nk,:,:))

end subroutine chemrates_tdep_ng

subroutine cpktkm_ng(fcp,fkt,fkm,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2,gask
  use fields_ng_module,only: flds,itp,t0
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: fcp,fkt,fkm

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,n2,he,mbar,po2,po1,pn2,phe,tni

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))

  nk = nlevp1_ng(i_ng)

  po2 = mbar*o2*rmassinv_o2
  po1 = mbar*o1*rmassinv_o1
  phe = mbar*he*rmassinv_he
  pn2 = mbar*n2*rmassinv_n2
  fcp = gask*.5*(po2*7./32.+po1*5./16.+phe*5./4.+pn2*7./28.)
  fkm = po2*4.03+po1*3.9+phe*3.84+pn2*3.42
  fkt = po2*56.+po1*75.9+phe*299.+pn2*56.

  tni = tn+t0
  tni(nk,:,:) = tn(nk-1,:,:)+t0

  fkm = fkm*tni**0.69*1.e-6
  fkt = fkt*tni**0.69

  do k = nk,2,-1
    fcp(k,:,:) = .5*(fcp(k,:,:)+fcp(k-1,:,:))
    fkm(k,:,:) = .5*(fkm(k,:,:)+fkm(k-1,:,:))
    fkt(k,:,:) = .5*(fkt(k,:,:)+fkt(k-1,:,:))
  enddo
  fcp(1,:,:) = 2.*fcp(1,:,:)-fcp(2,:,:)
  fkm(1,:,:) = 2.*fkm(1,:,:)-fkm(2,:,:)
  fkt(1,:,:) = 2.*fkt(1,:,:)-fkt(2,:,:)

end subroutine cpktkm_ng

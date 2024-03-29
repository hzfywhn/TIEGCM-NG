subroutine addiag_ng(vc,mbar,barm,xnmbar,xnmbari,scht,schti,z,zg,n2,i_ng)

  use params_module,only: nlevp1_ng,glat_ng
  use cons_module,only: dzgrav,gask,grav,p0,boltz, &
    rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2
  use fields_ng_module,only: flds,itp,dz,expz,expzmid_inv
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    vc,mbar,barm,xnmbar,xnmbari,scht,schti,z,zg,n2

  real,parameter :: dgtr = 1.74533E-2
  integer :: nk,latd0,latd1,k,lat
  real :: g0,r0,c2
  real,dimension(nlevp1_ng(i_ng)) :: expzi
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cs
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1) :: g
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tlbc,z_lbc,barm1
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,he,vn,tni,w1

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))

  cs = flds(i_ng)%cs
  tlbc = flds(i_ng)%tlbc
  z_lbc = flds(i_ng)%z_lbc

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1

  n2 = 1.-o2-o1-he
  call addfld(n2,'N2',i_ng)

  do lat = latd0,latd1
    vc(:,:,lat) = cs(lat)*vn(:,:,lat)
  enddo

  mbar = 1./(o2*rmassinv_o2+o1*rmassinv_o1+he*rmassinv_he+n2*rmassinv_n2)

  xnmbar = p0*mbar/(boltz*tn)
  do k = 1,nk
    xnmbar(k,:,:) = xnmbar(k,:,:)*expz(i_ng,k)
  enddo

  barm1 = 1.5*mbar(1,:,:)-0.5*mbar(2,:,:)
  do k = nk,2,-1
    barm(k,:,:) = 0.5*(mbar(k,:,:)+mbar(k-1,:,:))
  enddo
  barm(1,:,:) = barm1

  tni(1,:,:) = tlbc
  do k = 2,nk-1
    tni(k,:,:) = .5*(tn(k-1,:,:)+tn(k,:,:))
  enddo
  tni(nk,:,:) = tn(nk-1,:,:)

  expzi = expzmid_inv(i_ng)*expz(i_ng,1:nk)

  xnmbari = p0*barm/(boltz*tni)
  do k = 1,nk
    xnmbari(k,:,:) = xnmbari(k,:,:)*expzi(k)
  enddo

  scht = gask*tn/(mbar*grav)
  schti = gask*tni/(barm*grav)

  z(1,:,:) = z_lbc
  w1 = dz(i_ng)/dzgrav*tn/mbar
  do k = 1,nk-1
    z(k+1,:,:) = w1(k,:,:)+z(k,:,:)
  enddo

! calczg
  zg(1,:,:) = z(1,:,:)
  do lat = latd0,latd1
    c2 = cos(2.*dgtr*glat_ng(i_ng,lat))
    g0 = 980.616*(1.-.0026373*c2)
    r0 = 2.*g0/(3.085462e-6+2.27e-9*c2)
    g(1,:) = g0*(r0/(r0+0.5*(z(1,:,lat)+z(2,:,lat))))**2
    do k = 2,nk-1
      zg(k,:,lat) = zg(k-1,:,lat)+dz(i_ng)*scht(k-1,:,lat)*grav/g(k-1,:)
      g(k,:) = g0*(r0/(r0+0.5*(zg(k,:,lat)+z(k+1,:,lat))))**2
    enddo
  enddo
  zg(nk,:,:) = 2.0*zg(nk-1,:,:)-zg(nk-2,:,:)

  call addfld(zg,'ZG',i_ng)

end subroutine addiag_ng

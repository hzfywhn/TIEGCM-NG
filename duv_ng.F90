subroutine duv_ng(un_upd,unm_upd,vn_upd,vnm_upd,ulbc,vlbc,ulbc_nm,vlbc_nm,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use init_module,only: iday
  use cons_module,only: p0,gask,grav,re,dtsmooth_div2,dtsmooth,dzgrav,re_inv
  use fields_ng_module,only: flds,itp,itc,shapiro,dtx2inv,dlamda,dphi,dz,expz,xmue,bndry
  use char_module,only: find_index
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    un_upd,unm_upd,vn_upd,vnm_upd
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: ulbc,vlbc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: ulbc_nm,vlbc_nm

  real,parameter :: wt = 0.225
  integer :: nk,latd0,latd1,k,lat,idx_un,idx_vn
  logical,dimension(4) :: is_bndry
  real,dimension(nlevp1_ng(i_ng)) :: rtxmue
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cor,tanphi,racs
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: u_lbc,v_lbc,tlbc_nm
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,tn_upd,tn_nm,un,vn,un_nm,vn_nm,w_upd,mbar,barm,z,hdu,hdv,ui,vi,lxx,lyy,lxy,lyx,km,Fe,Fn, &
    g,dwdz,tni,advec_un,advec_vn,zl,zp,unm_smooth,vnm_smooth,eddyvisc,ss_un,ss_vn, &
    ztmp,tbar,dztbar,lxxi,lyyi,lxyi,lyxi,uii,vii,uni,vni,wi,g_1,qq_a
  real,dimension(2,nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: ss,xx,yy
  real,dimension(2,2,nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: pp,qq,rr,beta,gamma
  external :: advec_ng,div_ng,smooth_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  tn_upd = flds(i_ng)%tn(:,:,:,itc(i_ng))
  tn_nm = flds(i_ng)%tn_nm(:,:,:,itp(i_ng))
  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))
  un_nm = flds(i_ng)%un_nm(:,:,:,itp(i_ng))
  vn_nm = flds(i_ng)%vn_nm(:,:,:,itp(i_ng))
  w_upd = flds(i_ng)%w(:,:,:,itc(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  barm = flds(i_ng)%barm(:,:,:,itp(i_ng))
  z = flds(i_ng)%z(:,:,:,itp(i_ng))
  hdu = flds(i_ng)%hdu
  hdv = flds(i_ng)%hdv
  ui = flds(i_ng)%ui
  vi = flds(i_ng)%vi
  lxx = flds(i_ng)%lxx
  lyy = flds(i_ng)%lyy
  lxy = flds(i_ng)%lxy
  lyx = flds(i_ng)%lyx
  km = flds(i_ng)%km
  Fe = flds(i_ng)%Fe
  Fn = flds(i_ng)%Fn

  u_lbc = flds(i_ng)%u_lbc
  v_lbc = flds(i_ng)%v_lbc
  cor = flds(i_ng)%cor
  tanphi = flds(i_ng)%tanphi
  tlbc_nm = flds(i_ng)%tlbc_nm
  racs = flds(i_ng)%racs

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry

  call advec_ng(un,advec_un,i_ng)
  call advec_ng(vn,advec_vn,i_ng)

! glp
  tbar = (tn*(-2.)+tn_nm+tn_upd)*wt+tn

  dztbar = dz(i_ng)/dzgrav*tbar/mbar

  ztmp(1,:,:) = z(1,:,:)
  do k = 1,nk-1
    ztmp(k+1,:,:) = dztbar(k,:,:)+ztmp(k,:,:)
  enddo

! dldp
  call div_ng(ztmp,zl,1,i_ng)
  call div_ng(z,zp,-1,i_ng)

  do k = 1,nk-1
    zl(k,:,:) = (zl(k,:,:)+zl(k+1,:,:))*.5
    zp(k,:,:) = (zp(k,:,:)+zp(k+1,:,:))*.5
  enddo

  do lat = latd0,latd1
    zl(:,:,lat) = zl(:,:,lat)*racs(lat)
  enddo
  zl = zl*grav/dlamda(i_ng)
  zp = zp*grav*re_inv/dphi(i_ng)

  advec_un = advec_un+zl-hdu
  advec_vn = advec_vn+zp-hdv

  call smooth_ng(un_nm,unm_smooth,shapiro(i_ng),i_ng)
  call smooth_ng(vn_nm,vnm_smooth,shapiro(i_ng),i_ng)

  advec_un = advec_un-dtx2inv(i_ng)*unm_smooth
  advec_vn = advec_vn-dtx2inv(i_ng)*vnm_smooth

  do k = 1,nk-1
    lxxi(k,:,:) = .5*(lxx(k,:,:)+lxx(k+1,:,:))
    lyyi(k,:,:) = .5*(lyy(k,:,:)+lyy(k+1,:,:))
    lxyi(k,:,:) = .5*(lxy(k,:,:)+lxy(k+1,:,:))
    lyxi(k,:,:) = .5*(lyx(k,:,:)+lyx(k+1,:,:))
    uii(k,:,:) = .5*(ui(k,:,:)+ui(k+1,:,:))
    vii(k,:,:) = .5*(vi(k,:,:)+vi(k+1,:,:))
    uni(k,:,:) = .5*(un(k,:,:)+un(k+1,:,:))
    vni(k,:,:) = .5*(vn(k,:,:)+vn(k+1,:,:))
    wi(k,:,:) = .5*(w_upd(k,:,:)+w_upd(k+1,:,:))
  enddo
  lxxi(nk,:,:) = 1.5*lxx(nk,:,:)-.5*lxx(nk-1,:,:)
  lyyi(nk,:,:) = 1.5*lyy(nk,:,:)-.5*lyy(nk-1,:,:)
  lxyi(nk,:,:) = 1.5*lxy(nk,:,:)-.5*lxy(nk-1,:,:)
  lyxi(nk,:,:) = 1.5*lyx(nk,:,:)-.5*lyx(nk-1,:,:)
  uii(nk,:,:) = 1.5*ui(nk,:,:)-.5*ui(nk-1,:,:)
  vii(nk,:,:) = 1.5*vi(nk,:,:)-.5*vi(nk-1,:,:)
  uni(nk,:,:) = 1.5*un(nk,:,:)-.5*un(nk-1,:,:)
  vni(nk,:,:) = 1.5*vn(nk,:,:)-.5*vn(nk-1,:,:)
  wi(nk,:,:) = 1.5*w_upd(nk,:,:)-.5*w_upd(nk-1,:,:)

  ss(1,:,:,:) = lxxi*uii+lxyi*vii-advec_un+Fe
  ss(2,:,:,:) = lyyi*vii-lyxi*uii-advec_vn+Fn
  do k = 1,nk
    ss(:,k,:,:) = ss(:,k,:,:)*expz(i_ng,k)
  enddo
  ss_un = ss(1,:,:,:)
  ss_vn = ss(2,:,:,:)

  do k = 2,nk-1
    rtxmue(k) = sqrt(xmue(i_ng,k-1,iday)*xmue(i_ng,k,iday))
  enddo
  rtxmue(1) = sqrt(xmue(i_ng,1,iday)**3/xmue(i_ng,2,iday))
  rtxmue(nk) = sqrt(xmue(i_ng,nk-1,iday)**3/xmue(i_ng,nk-2,iday))

  do k = 2,nk-1
    tni(k,:,:) = .5*(tn(k-1,:,:)+tn(k,:,:))
  enddo
  tni(1,:,:) = tlbc_nm
  tni(nk,:,:) = tn(nk-1,:,:)

  eddyvisc = p0*gask*tn/(mbar*grav**2)
  do k = 1,nk-1
    eddyvisc(k,:,:) = eddyvisc(k,:,:)*expz(i_ng,k)*xmue(i_ng,k,iday)
  enddo
  eddyvisc(nk,:,:) = sqrt(eddyvisc(nk-1,:,:)**3/eddyvisc(nk-2,:,:))

  do k = nk,2,-1
    eddyvisc(k,:,:) = sqrt(eddyvisc(k-1,:,:)*eddyvisc(k,:,:))
  enddo
  eddyvisc(1,:,:) = eddyvisc(1,:,:)**2/eddyvisc(2,:,:)

  g = grav**2*(km+eddyvisc)*barm/(p0*gask*tni*dz(i_ng)**2)

  dwdz = wi/(2.*dz(i_ng))
  do k = 1,nk
    dwdz(k,:,:) = dwdz(k,:,:)*expz(i_ng,k)
  enddo

  do k = 1,nk-1
    g_1(k,:,:) = g(k+1,:,:)
  enddo
  g_1(nk,:,:) = 2.*g(nk,:,:)-g(nk-1,:,:)

  pp(1,1,:,:,:) = -(g+dwdz)
  rr(1,1,:,:,:) = -(g_1-dwdz)
  pp(2,2,:,:,:) = pp(1,1,:,:,:)
  rr(2,2,:,:,:) = rr(1,1,:,:,:)
  pp(1,2,:,:,:) = 0.
  rr(1,2,:,:,:) = 0.
  pp(2,1,:,:,:) = 0.
  rr(2,1,:,:,:) = 0.

  qq(1,1,:,:,:) = g+g_1
  qq(2,2,:,:,:) = qq(1,1,:,:,:)

  qq_a = dtx2inv(i_ng)+lxxi
  do k = 1,nk
    qq_a(k,:,:) = qq_a(k,:,:)*expz(i_ng,k)
  enddo
  qq(1,1,:,:,:) = qq(1,1,:,:,:)+qq_a

  qq_a = dtx2inv(i_ng)+lyyi
  do k = 1,nk
    qq_a(k,:,:) = qq_a(k,:,:)*expz(i_ng,k)
  enddo
  qq(2,2,:,:,:) = qq(2,2,:,:,:)+qq_a

  do lat = latd0,latd1
    qq(1,2,:,:,lat) = cor(lat)+un(:,:,lat)/re*tanphi(lat)
  enddo
  qq(2,1,:,:,:) = qq(1,2,:,:,:)

  qq_a = qq(1,2,:,:,:)-lxyi
  do k = 1,nk
    qq_a(k,:,:) = qq_a(k,:,:)*expz(i_ng,k)
  enddo
  qq(1,2,:,:,:) = -qq_a

  qq_a = qq(2,1,:,:,:)-lyxi
  do k = 1,nk
    qq_a(k,:,:) = qq_a(k,:,:)*expz(i_ng,k)
  enddo
  qq(2,1,:,:,:) = qq_a

  qq(:,:,1,:,:) = qq(:,:,1,:,:)-pp(:,:,1,:,:)
  ss(1,1,:,:) = ss(1,1,:,:)-2.*(pp(1,1,1,:,:)*u_lbc+pp(1,2,1,:,:)*v_lbc)
  ss(2,1,:,:) = ss(2,1,:,:)-2.*(pp(2,1,1,:,:)*u_lbc+pp(2,2,1,:,:)*v_lbc)
  ss_un(1,:,:) = ss(1,1,:,:)
  ss_vn(1,:,:) = ss(2,1,:,:)
  pp(:,:,1,:,:) = 0.

  qq(:,:,nk-1,:,:) = qq(:,:,nk-1,:,:)+rr(:,:,nk-1,:,:)
  rr(:,:,nk-1,:,:) = 0.

! blktri
  yy(1,1,:,:) = qq(1,1,1,:,:)*qq(2,2,1,:,:)-qq(1,2,1,:,:)*qq(2,1,1,:,:)
  beta(1,1,1,:,:) = qq(2,2,1,:,:)/yy(1,1,:,:)
  beta(1,2,1,:,:) = -qq(1,2,1,:,:)/yy(1,1,:,:)
  beta(2,1,1,:,:) = -qq(2,1,1,:,:)/yy(1,1,:,:)
  beta(2,2,1,:,:) = qq(1,1,1,:,:)/yy(1,1,:,:)
  yy(1,1,:,:) = beta(1,1,1,:,:)*ss(1,1,:,:)+beta(1,2,1,:,:)*ss(2,1,:,:)
  yy(2,1,:,:) = beta(2,1,1,:,:)*ss(1,1,:,:)+beta(2,2,1,:,:)*ss(2,1,:,:)
  do k = 2,nk-1
    gamma(1,1,k-1,:,:) = beta(1,1,k-1,:,:)*rr(1,1,k-1,:,:)+beta(1,2,k-1,:,:)*rr(2,1,k-1,:,:)
    gamma(1,2,k-1,:,:) = beta(1,1,k-1,:,:)*rr(1,2,k-1,:,:)+beta(1,2,k-1,:,:)*rr(2,2,k-1,:,:)
    gamma(2,1,k-1,:,:) = beta(2,1,k-1,:,:)*rr(1,1,k-1,:,:)+beta(2,2,k-1,:,:)*rr(2,1,k-1,:,:)
    gamma(2,2,k-1,:,:) = beta(2,1,k-1,:,:)*rr(1,2,k-1,:,:)+beta(2,2,k-1,:,:)*rr(2,2,k-1,:,:)
    gamma(1,1,k,:,:) = qq(1,1,k,:,:)-pp(1,1,k,:,:)*gamma(1,1,k-1,:,:)-pp(1,2,k,:,:)*gamma(2,1,k-1,:,:)
    gamma(1,2,k,:,:) = qq(1,2,k,:,:)-pp(1,1,k,:,:)*gamma(1,2,k-1,:,:)-pp(1,2,k,:,:)*gamma(2,2,k-1,:,:)
    gamma(2,1,k,:,:) = qq(2,1,k,:,:)-pp(2,1,k,:,:)*gamma(1,1,k-1,:,:)-pp(2,2,k,:,:)*gamma(2,1,k-1,:,:)
    gamma(2,2,k,:,:) = qq(2,2,k,:,:)-pp(2,1,k,:,:)*gamma(1,2,k-1,:,:)-pp(2,2,k,:,:)*gamma(2,2,k-1,:,:)
    yy(1,k,:,:) = gamma(1,1,k,:,:)*gamma(2,2,k,:,:)-gamma(1,2,k,:,:)*gamma(2,1,k,:,:)
    beta(1,1,k,:,:) = gamma(2,2,k,:,:)/yy(1,k,:,:)
    beta(1,2,k,:,:) = -gamma(1,2,k,:,:)/yy(1,k,:,:)
    beta(2,1,k,:,:) = -gamma(2,1,k,:,:)/yy(1,k,:,:)
    beta(2,2,k,:,:) = gamma(1,1,k,:,:)/yy(1,k,:,:)
    xx(1,k,:,:) = ss(1,k,:,:)-pp(1,1,k,:,:)*yy(1,k-1,:,:)-pp(1,2,k,:,:)*yy(2,k-1,:,:)
    xx(2,k,:,:) = ss(2,k,:,:)-pp(2,1,k,:,:)*yy(1,k-1,:,:)-pp(2,2,k,:,:)*yy(2,k-1,:,:)
    yy(1,k,:,:) = beta(1,1,k,:,:)*xx(1,k,:,:)+beta(1,2,k,:,:)*xx(2,k,:,:)
    yy(2,k,:,:) = beta(2,1,k,:,:)*xx(1,k,:,:)+beta(2,2,k,:,:)*xx(2,k,:,:)
  enddo
  xx(:,nk-1,:,:) = yy(:,nk-1,:,:)
  do k = nk-2,1,-1
    xx(1,k,:,:) = yy(1,k,:,:)-gamma(1,1,k,:,:)*xx(1,k+1,:,:)-gamma(1,2,k,:,:)*xx(2,k+1,:,:)
    xx(2,k,:,:) = yy(2,k,:,:)-gamma(2,1,k,:,:)*xx(1,k+1,:,:)-gamma(2,2,k,:,:)*xx(2,k+1,:,:)
  enddo

  xx(:,nk,:,:) = 2*xx(:,nk-1,:,:)-xx(:,nk-2,:,:)

  un_upd = xx(1,:,:,:)
  vn_upd = xx(2,:,:,:)

  idx_un = find_index('UN',bndry)
  idx_vn = find_index('VN',bndry)
  if (is_bndry(1)) then
    un_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_un)
    vn_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_vn)
  endif
  if (is_bndry(2)) then
    un_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_un)
    vn_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_vn)
  endif
  if (is_bndry(3)) then
    un_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_un)
    vn_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_vn)
  endif
  if (is_bndry(4)) then
    un_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_un)
    vn_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_vn)
  endif

  ulbc_nm = ulbc
  vlbc_nm = vlbc
  ulbc = u_lbc
  vlbc = v_lbc

  unm_upd = dtsmooth_div2*(un_nm+un_upd)+dtsmooth*un
  vnm_upd = dtsmooth_div2*(vn_nm+vn_upd)+dtsmooth*vn

end subroutine duv_ng

subroutine dt_ng(tn_upd,tn_nm_upd,tlbc,tlbc_nm,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use init_module,only: iday
  use cons_module,only: rmassinv_o1,tsurplus,p0,avo,grav,gask,dtsmooth,dtsmooth_div2
  use fields_ng_module,only: flds,itp,itc,shapiro,dtx2inv,dz,expzmid_inv,expz,dift,bndry
  use char_module,only: find_index
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    tn_upd,tn_nm_upd
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: tlbc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: tlbc_nm

  integer :: nk,k,idx_tn
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: t_lbc,ulbc,vlbc
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,tn_nm,un,vn,o1,mbar,xnmbar,scht,schti,cp,kt,km,hdt,qji_tn,cool_imp,cool_exp,w_upd, &
    qtotal,rkm12,cptn,qm,total_heat,dudz,dvdz,g,f,h,rho,tni,p,q,r,rhs,qpart,solarq,recombq, &
    conductq,coolq,advecq,adiaq,tnsmooth,advec_tn,advec,cpi,kmi,wi,total_heat_a,gpart,gpart_a,p_1,r_0,qpart_a
  external :: advec_ng,smooth_ng,trsolv_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  tn_nm = flds(i_ng)%tn_nm(:,:,:,itp(i_ng))
  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  scht = flds(i_ng)%scht(:,:,:,itp(i_ng))
  schti = flds(i_ng)%schti(:,:,:,itp(i_ng))
  cp = flds(i_ng)%cp
  kt = flds(i_ng)%kt
  km = flds(i_ng)%km
  hdt = flds(i_ng)%hdt
  qji_tn = flds(i_ng)%qji_tn
  cool_imp = flds(i_ng)%cool_implicit
  cool_exp = flds(i_ng)%cool_explicit
  w_upd = flds(i_ng)%w(:,:,:,itc(i_ng))

  t_lbc = flds(i_ng)%t_lbc
  qtotal = flds(i_ng)%qtotal
  rkm12 = flds(i_ng)%rkm12
  ulbc = flds(i_ng)%ulbc
  vlbc = flds(i_ng)%vlbc

  nk = nlevp1_ng(i_ng)
  is_bndry = flds(i_ng)%is_bndry

  call advec_ng(tn,advec_tn,i_ng)

! advecv
  advec(1,:,:) = (tn(1,:,:)-t_lbc)*w_upd(1,:,:)*2.
  advec(nk,:,:) = 0.
  do k = 2,nk-1
    advec(k,:,:) = (tn(k,:,:)-tn(k-1,:,:))*w_upd(k,:,:)
  enddo
  do k = 1,nk-1
    advec(k,:,:) = (advec(k,:,:)+advec(k+1,:,:))*.5
  enddo
  advec_tn = advec_tn+advec/dz(i_ng)

  call smooth_ng(tn_nm,tnsmooth,shapiro(i_ng),i_ng)

  do k = 1,nk-1
    cpi(k,:,:) = .5*(cp(k,:,:)+cp(k+1,:,:))
    kmi(k,:,:) = .5*(km(k,:,:)+km(k+1,:,:))
    wi(k,:,:) = .5*(w_upd(k,:,:)+w_upd(k+1,:,:))
  enddo
  cpi(nk,:,:) = 1.5*cp(nk,:,:)-.5*cp(nk-1,:,:)
  kmi(nk,:,:) = 1.5*km(nk,:,:)-.5*km(nk-1,:,:)
  wi(nk,:,:) = 1.5*w_upd(nk,:,:)-.5*w_upd(nk-1,:,:)

  cptn = cpi*(advec_tn-dtx2inv(i_ng)*tnsmooth)
  do k = 1,nk
    cptn(k,:,:) = cptn(k,:,:)*expz(i_ng,k)
  enddo

  do k = 2,nk-2
    dudz(k,:,:) = (un(k+1,:,:)-un(k-1,:,:))/(2.*dz(i_ng))
    dvdz(k,:,:) = (vn(k+1,:,:)-vn(k-1,:,:))/(2.*dz(i_ng))
  enddo
  dudz(1,:,:) = (un(1,:,:)+1./3.*un(2,:,:)-4./3.*ulbc)/dz(i_ng)
  dvdz(1,:,:) = (vn(1,:,:)+1./3.*vn(2,:,:)-4./3.*vlbc)/dz(i_ng)
  dudz(nk-1,:,:) = dudz(nk-2,:,:)/3.
  dvdz(nk-1,:,:) = dvdz(nk-2,:,:)/3.
  dudz(nk,:,:) = dudz(nk-1,:,:)/3.
  dvdz(nk,:,:) = dvdz(nk-1,:,:)/3.

  do k = 1,nk-1
    solarq(k,:,:) = .5*(qtotal(k,:,:)+qtotal(k+1,:,:))
  enddo
  solarq(nk,:,:) = 1.5*qtotal(nk,:,:)-.5*qtotal(nk-1,:,:)

  recombq = tsurplus*rkm12*(xnmbar*o1*rmassinv_o1)**2*avo/mbar

  qm = grav*kmi/(p0*scht)*(dudz**2+dvdz**2)
  do k = 1,nk
    qm(k,:,:) = qm(k,:,:)/expz(i_ng,k)
  enddo

  total_heat = solarq+hdt+recombq+qji_tn+qm

  call addfld(solarq,'DT_SOLAR',i_ng)
  call addfld(hdt,'DT_HORDIF',i_ng)
  call addfld(recombq,'DT_RECOMB',i_ng)
  call addfld(qji_tn,'DT_JOULE',i_ng)
  call addfld(qm,'DT_MOLDIF',i_ng)

  do k = 1,nk
    total_heat_a(k,:,:) = total_heat(k,:,:)*expz(i_ng,k)
  enddo
  cptn = cptn-total_heat_a

  do k = 2,nk-1
    tni(k,:,:) = .5*(tn(k-1,:,:)+tn(k,:,:))
  enddo
  tni(1,:,:) = tlbc
  tni(nk,:,:) = tn(nk-1,:,:)

  h = schti

  rho = p0*expzmid_inv(i_ng)/(h*grav)
  do k = 1,nk
    rho(k,:,:) = rho(k,:,:)*expz(i_ng,k)
  enddo

  call addfld(rho,'DEN',i_ng)

  gpart = h**2*rho*cp
  do k = 1,nk
    gpart_a(k,:,:) = gpart(k,:,:)*dift(i_ng,k,iday)
  enddo
  g = grav*(kt+gpart_a)/(p0*h*dz(i_ng)**2)

  f = grav**2*h**2*rho/(tni*p0*2.*dz(i_ng))
  do k = 1,nk
    f(k,:,:) = f(k,:,:)*dift(i_ng,k,iday)
  enddo

  p = g-f

  r_0 = g+f
  do k = 1,nk-1
    p_1(k,:,:) = p(k+1,:,:)
    r(k,:,:) = r_0(k+1,:,:)
  enddo
  p_1(nk,:,:) = 2.*p(nk,:,:)-p(nk-1,:,:)
  r(nk,:,:) = 2.*r_0(nk,:,:)-r_0(nk-1,:,:)

  q = -r_0-p_1
  rhs = cptn

  q(nk-1,:,:) = -r_0(nk-1,:,:)
  q(nk,:,:) = -r_0(nk,:,:)
  r(nk-1,:,:) = 0.
  r(nk,:,:) = 0.

  qpart = cpi*(dtx2inv(i_ng)+cool_imp)+wi*gask/mbar
  rhs = rhs+cool_exp
  do k = 1,nk
    qpart_a(k,:,:) = qpart(k,:,:)*expz(i_ng,k)
  enddo
  q = q-qpart_a

  q(1,:,:) = q(1,:,:)-p(1,:,:)
  rhs(1,:,:) = rhs(1,:,:)-2.*p(1,:,:)*t_lbc
  p(1,:,:) = 0.

  call trsolv_ng(p,q,r,rhs,tn_upd,nk-1,i_ng)

  tn_upd(nk,:,:) = tn_upd(nk-1,:,:)**2/tn_upd(nk-2,:,:)

  do k = 2,nk-1
    conductq(k,:,:) = (tn_upd(k-1,:,:)*(-f(k,:,:)+g(k,:,:))+ &
      tn_upd(k,:,:)*(f(k+1,:,:)-f(k,:,:)-g(k+1,:,:)-g(k,:,:))+ &
      tn_upd(k+1,:,:)*(f(k+1,:,:)+g(k+1,:,:)))/expz(i_ng,k)
  enddo
  coolq = cpi*cool_imp*tn_upd
  do k = 1,nk
    coolq(k,:,:) = -cool_exp(k,:,:)/expz(i_ng,k)-coolq(k,:,:)
  enddo
  advecq = -cpi*advec_tn
  adiaq = -wi*gask/mbar*tn_upd
  call addfld(conductq,'DT_CONDUCT',i_ng)
  call addfld(coolq,'DT_COOL',i_ng)
  call addfld(advecq,'DT_ADVEC',i_ng)
  call addfld(adiaq,'DT_ADIA',i_ng)

! a hook to apply Dirichlet lateral boundary conditions, also see duv,comp,oplus,swdot
  idx_tn = find_index('TN',bndry)
  if (is_bndry(1)) tn_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_tn)
  if (is_bndry(2)) tn_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_tn)
  if (is_bndry(3)) tn_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_tn)
  if (is_bndry(4)) tn_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_tn)

  tn_nm_upd = dtsmooth_div2*(tn_nm+tn_upd)+dtsmooth*tn

  tlbc_nm = tlbc
  tlbc = t_lbc

  tn_upd = max(tn_upd,100.)

end subroutine dt_ng

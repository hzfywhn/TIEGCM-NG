subroutine oplus_ng(op,optm1,opout,optm1out,xiop2p,xiop2d,Fe,Fn,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,zpmid_ng,glat_ng
  use cons_module,only: rmass_op,gask,grav,re,rmassinv_o2,rmassinv_o1, &
    rmassinv_n2,rmassinv_n2d,dtsmooth,dtsmooth_div2,pi,rtd,rmassinv_he
  use chemrates_module,only: rk10,rk16,rk17,rk18,rk21,rk22,rk23,rk24,rk26,rk27
  use input_module,only: enforce_opfloor,colfac,opdiffcap,nstep_sub
  use fields_ng_module,only: flds,itp,shapiro,dtx2inv,dlamda,dphi,dlev,dz
  use dffm_ng_module,only: df_2d,df_3d
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: op,optm1
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    opout,optm1out,xiop2p,xiop2d,Fe,Fn

  real,parameter :: explic = 1., opmin = 3000., phid = 2.0e8, phin = -2.0e8, ppolar = 0.
  integer :: nk,latd0,latd1,k,lat
  real :: gmr,opfloor,mgr
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cs
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    bx,by,bz,bmod2,dipmag,sndec,csdec,chi,rlatm,opflux,dvb,ubca,ubcb,a,fed,fen,cs_by,dbxdx,dcs_bydy,sncsdip
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,te,ti,o2,o1,n2,he,n2d,ne,u,v,w,mbar,ui,vi,wi,xnmbar,qop2p,qop2d,qop,rk1,rk2,rk19,rk20,rk25, &
    bdzdvb_op,explicit,hdz,tphdz1,tphdz0,djint,divbz,hdzmbz,hdzpbz,p_coeff,q_coeff,r_coeff,bdotu, &
    op_loss,tp1,hj,bvel,diffj,tp,tr,bdotdh_op,bdotdh_opj,bdotdh_diff,dj,optm1_smooth,vni,wd, &
    wni,uii,vii,wii,qop2pi,qop2di,qopi,hdzi,tp_op,op_bmod2,dbveldx,dbveldy,dop_bmod2dx,dop_bmod2dy, &
    exp1,exp2,tp1_0,tp1_1,dj_bz,ddj_bzdx,ddj_bzdy,divbz1,divbz2,djint_tphdz0,djint_tphdz0_1, &
    djint_tphdz1,djint_tphdz1_1,bz_bdotu,bz_bdotu_1,dvb_bz
  logical,external :: isclose
  external :: smooth_ng,trsolv_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2d = flds(i_ng)%n2d(:,:,:,itp(i_ng))
  ne = flds(i_ng)%ne(:,:,:,itp(i_ng))
  u = flds(i_ng)%un(:,:,:,itp(i_ng))
  v = flds(i_ng)%vn(:,:,:,itp(i_ng))
  w = flds(i_ng)%w(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  ui = flds(i_ng)%ui
  vi = flds(i_ng)%vi
  wi = flds(i_ng)%wi
  xnmbar = flds(i_ng)%xnmbar

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry
  cs = flds(i_ng)%cs
  qop2p = flds(i_ng)%qop2p
  qop2d = flds(i_ng)%qop2d
  qop = flds(i_ng)%qop
  rk1 = flds(i_ng)%rk1
  rk2 = flds(i_ng)%rk2
  rk19 = flds(i_ng)%rk19
  rk20 = flds(i_ng)%rk20
  rk25 = flds(i_ng)%rk25
  bx = flds(i_ng)%bx
  by = flds(i_ng)%by
  bz = flds(i_ng)%bz
  bmod2 = flds(i_ng)%bmod2
  dipmag = flds(i_ng)%dipmag
  sndec = flds(i_ng)%sndec
  csdec = flds(i_ng)%csdec
  chi = flds(i_ng)%chi
  rlatm = flds(i_ng)%rlatm

! oplus_flux
  where (abs(rlatm) >= pi/24.)
    a = 1.
  elsewhere
    a = .5*(1.+sin(pi*(abs(rlatm)-pi/48.)/(pi/24.)))
  endwhere
  where (a < 0.05) a = 0.05
  fed = phid*a
  fen = phin*a
  where (chi >= 0.5*pi)
    opflux = fen
  elsewhere
    opflux = fed
  endwhere
  where ((chi*rtd-80.)*(chi*rtd-100.) < 0.) &
    opflux = .5*(fed+fen)+.5*(fed-fen)*cos(pi*(chi*rtd-80.)/20.)
  where (abs(rlatm) >= pi/3.) opflux = opflux+ppolar

! divb
  do lat = latd0,latd1
    cs_by(:,lat) = cs(lat)*by(:,lat)
  enddo
  call df_2d(bx,dbxdx,1,i_ng)
  call df_2d(cs_by,dcs_bydy,-1,i_ng)

  dvb = dbxdx/dlamda(i_ng)+dcs_bydy/dphi(i_ng)
  do lat = latd0,latd1
    dvb(:,lat) = dvb(:,lat)/cs(lat)
  enddo
  dvb = (dvb+2.*bz)/re

  tr = 0.5*(tn+ti)

  do k = 1,nk-1
    wni(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
    uii(k,:,:) = .5*(ui(k,:,:)+ui(k+1,:,:))
    vii(k,:,:) = .5*(vi(k,:,:)+vi(k+1,:,:))
    wii(k,:,:) = .5*(wi(k,:,:)+wi(k+1,:,:))
    qop2pi(k,:,:) = .5*(qop2p(k,:,:)+qop2p(k+1,:,:))
    qop2di(k,:,:) = .5*(qop2d(k,:,:)+qop2d(k+1,:,:))
    qopi(k,:,:) = .5*(qop(k,:,:)+qop(k+1,:,:))
  enddo
  wni(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)
  uii(nk,:,:) = 1.5*ui(nk,:,:)-.5*ui(nk-1,:,:)
  vii(nk,:,:) = 1.5*vi(nk,:,:)-.5*vi(nk-1,:,:)
  wii(nk,:,:) = 1.5*wi(nk,:,:)-.5*wi(nk-1,:,:)
  qop2pi(nk,:,:) = 1.5*qop2p(nk,:,:)-.5*qop2p(nk-1,:,:)
  qop2di(nk,:,:) = 1.5*qop2d(nk,:,:)-.5*qop2d(nk-1,:,:)
  qopi(nk,:,:) = 1.5*qop(nk,:,:)-.5*qop(nk-1,:,:)

! rrk
  vni = o1*rmassinv_o1*sqrt(tr)*(1.-0.064*log10(tr))**2*colfac+ &
    18.6*n2*rmassinv_n2+18.1*o2*rmassinv_o2+3.6*he*rmassinv_he
  dj = 1.42E17/(xnmbar*vni)
  vni = 16*3.53E-11*vni
  if (isclose(opdiffcap,0.)) dj = min(dj,opdiffcap)

  tp = te+ti

  hj = gask*tn/(mbar*grav)

  do k = 1,nk
    bdotu(k,:,:) = bx*u(k,:,:)+by*v(k,:,:)+hj(k,:,:)*bz*wni(k,:,:)
  enddo

  bvel = bdotu*op

! diffus
  mgr = rmass_op*grav/gask
  tp_op = tp*op
  do k = 1,nk-2
    diffj(k+1,:,:) = (tp_op(k+2,:,:)-tp_op(k,:,:))/2.
  enddo
  diffj(nk-1,:,:) = tp_op(nk-1,:,:)-tp_op(nk-2,:,:)
  diffj(nk,:,:) = tp_op(nk,:,:)-tp_op(nk-1,:,:)
  diffj(1,:,:) = tp_op(2,:,:)-tp_op(1,:,:)
  diffj = 1./(hj*dlev(i_ng))*diffj+mgr*op

  tp = tp_op

  call smooth_ng(optm1,optm1_smooth,shapiro(i_ng)/nstep_sub,i_ng)

  wd = vni*diffj*dj
  sncsdip = sin(dipmag)*cos(dipmag)
  do k = 1,nk
    Fe(k,:,:) = wd(k,:,:)*sncsdip*sndec
    Fn(k,:,:) = wd(k,:,:)*sncsdip*csdec
  enddo

  call bdotdh(diffj,bdotdh_op)
  do k = 1,nk
    bdotdh_op(k,:,:) = dj(k,:,:)*bz*bdotdh_op(k,:,:)
  enddo

  call bdotdh(tp,bdotdh_opj)
  bdotdh_opj = bdotdh_opj*dj

  call bdotdh(bdotdh_opj,bdotdh_diff)

! bdzdvb
  do k = 2,nk-2
    bdzdvb_op(k,:,:) = bz/(2.*hj(k,:,:)*dz(i_ng))*(bdotdh_opj(k+1,:,:)-bdotdh_opj(k-1,:,:))+ &
      dvb*bdotdh_opj(k,:,:)
  enddo
  bdzdvb_op(nk-1,:,:) = bz/(hj(nk-1,:,:)*dz(i_ng))*(bdotdh_opj(nk-1,:,:)-bdotdh_opj(nk-2,:,:))+ &
    dvb*bdotdh_opj(nk-1,:,:)
  bdzdvb_op(nk,:,:) = bz/(hj(nk,:,:)*dz(i_ng))*(bdotdh_opj(nk,:,:)-bdotdh_opj(nk-1,:,:))+ &
    dvb*bdotdh_opj(nk,:,:)
  bdzdvb_op(1,:,:) = bz/(hj(1,:,:)*dz(i_ng))*(bdotdh_opj(2,:,:)-bdotdh_opj(1,:,:))+ &
    dvb*bdotdh_opj(1,:,:)

  explicit = -explic*(bdzdvb_op+bdotdh_diff+bdotdh_op)

  do k = 1,nk
    op_bmod2(k,:,:) = op(k,:,:)/bmod2**2
  enddo
  call df_3d(bvel,dbveldx,1,i_ng)
  call df_3d(op_bmod2,dop_bmod2dx,1,i_ng)
  call df_3d(bvel,dbveldy,-1,i_ng)
  call df_3d(op_bmod2,dop_bmod2dy,-1,i_ng)

  do k = 1,nk
    exp1(k,:,:) = bx*dbveldx(k,:,:)+uii(k,:,:)*bmod2**2*dop_bmod2dx(k,:,:)
    exp2(k,:,:) = by*dbveldy(k,:,:)+vii(k,:,:)*bmod2**2*dop_bmod2dy(k,:,:)
  enddo
  do lat = latd0,latd1
    exp1(:,:,lat) = exp1(:,:,lat)/cs(lat)
  enddo

  explicit = explicit+1./re*(exp1/dlamda(i_ng)+exp2/dphi(i_ng))

  dvb = dvb/bz
  hdz = 1./(hj*dz(i_ng))
  tp1 = 0.5*(ti+te)

  do k = 1,nk-2
    hdzi(k+1,:,:) = 0.5*(hdz(k,:,:)+hdz(k+1,:,:))
  enddo
  hdzi(1,:,:) = 1.5*hdz(1,:,:)-0.5*hdz(2,:,:)
  hdzi(nk,:,:) = 1.5*hdz(nk-1,:,:)-0.5*hdz(nk-2,:,:)

  gmr = grav*rmass_op/(2.*gask)

  tp1_0 = tp1
  tp1_0(nk,:,:) = 2.*tp1(nk-1,:,:)-tp1(nk-2,:,:)
  tphdz1 = 2.*tp1_0*hdzi+gmr

  do k = 1,nk-1
    tp1_1(k+1,:,:) = tp1(k,:,:)
  enddo
  tp1_1(1,:,:) = 2.*tp1(1,:,:)-tp1(2,:,:)
  tphdz0 = 2.*tp1_1*hdzi-gmr

  do k = 1,nk-2
    djint(k+1,:,:) = 0.5*(dj(k,:,:)+dj(k+1,:,:))
  enddo
  djint(1,:,:) = 1.5*dj(1,:,:)-0.5*dj(2,:,:)
  djint(nk,:,:) = 1.5*dj(nk-1,:,:)-0.5*dj(nk-2,:,:)

  do k = 1,nk
    dj_bz(k,:,:) = dj(k,:,:)*bz
  enddo
  call df_3d(dj_bz,ddj_bzdx,1,i_ng)
  call df_3d(dj_bz,ddj_bzdy,-1,i_ng)

  do k = 1,nk
    divbz1(k,:,:) = bx*ddj_bzdx(k,:,:)
    divbz2(k,:,:) = by*ddj_bzdy(k,:,:)
  enddo
  do lat = latd0,latd1
    divbz1(:,:,lat) = divbz1(:,:,lat)/cs(lat)
  enddo

  divbz = 1./(re*dj)*(divbz1/dlamda(i_ng)+divbz2/dphi(i_ng))
  do k = 1,nk
    divbz(k,:,:) = dvb+divbz(k,:,:)/bz**2
  enddo

  hdzmbz = hdz-0.5*divbz
  hdzpbz = hdz+0.5*divbz
  do k = 1,nk
    hdzmbz(k,:,:) = hdzmbz(k,:,:)*bz**2
    hdzpbz(k,:,:) = hdzpbz(k,:,:)*bz**2
  enddo

  explicit = explicit-optm1_smooth*dtx2inv(i_ng)*nstep_sub

  djint_tphdz0 = djint*tphdz0
  do k = 1,nk-1
    djint_tphdz0_1(k,:,:) = djint_tphdz0(k+1,:,:)
  enddo
  djint_tphdz0_1(nk,:,:) = 2*djint_tphdz0(nk,:,:)-djint_tphdz0(nk-1,:,:)

  djint_tphdz1 = djint*tphdz1
  do k = 1,nk-1
    djint_tphdz1_1(k,:,:) = djint_tphdz1(k+1,:,:)
  enddo
  djint_tphdz1_1(nk,:,:) = 2*djint_tphdz1(nk,:,:)-djint_tphdz1(nk-1,:,:)

  p_coeff = hdzmbz*djint_tphdz0
  q_coeff = -(hdzpbz*djint_tphdz0_1+hdzmbz*djint_tphdz1)
  r_coeff = hdzpbz*djint_tphdz1_1

  do k = 1,nk
    bz_bdotu(k,:,:) = bz*bdotu(k,:,:)
  enddo

  do k = 1,nk-1
    bz_bdotu_1(k+1,:,:) = bz_bdotu(k,:,:)
  enddo
  bz_bdotu_1(1,:,:) = 2.*bz_bdotu(1,:,:)-bz_bdotu(2,:,:)
  p_coeff = p_coeff+(bz_bdotu_1+wii)*0.5*hdz

  q_coeff = q_coeff-wii*6./re

  do k = 1,nk-2
    bz_bdotu_1(k,:,:) = bz_bdotu(k+1,:,:)
  enddo
  bz_bdotu_1(nk-1,:,:) = 2.*bz_bdotu(nk-1,:,:)-bz_bdotu(nk-2,:,:)
  bz_bdotu_1(nk,:,:) = 2.*bz_bdotu(nk,:,:)-bz_bdotu(nk-1,:,:)
  r_coeff = r_coeff-(bz_bdotu_1+wii)*0.5*hdz

  do k = 1,nk
    dvb_bz(k,:,:) = dvb*bz
  enddo
  q_coeff = q_coeff-bdotu*dvb_bz-dtx2inv(i_ng)*nstep_sub

  ubcb = -bz**2*djint(nk,:,:)*tphdz0(nk,:,:)
  ubca = -bz**2*djint(nk,:,:)*tphdz1(nk,:,:)

  q_coeff(nk-1,:,:) = q_coeff(nk-1,:,:)+ubcb/ubca*r_coeff(nk-1,:,:)

  explicit(nk-1,:,:) = explicit(nk-1,:,:)-opflux*r_coeff(nk-1,:,:)/ubca
  r_coeff(nk-1,:,:) = 0.

  xiop2p = qop2pi/(xnmbar*((rk16+rk17)*n2*rmassinv_n2+rk18*o1*rmassinv_o1)+(rk19+rk20)*ne+rk21+rk22)

  xiop2d = (qop2di+(rk20*ne+rk22)*xiop2p)/ &
    (xnmbar*(rk23*n2*rmassinv_n2+rk24*o1*rmassinv_o1+rk26*o2*rmassinv_o2)+rk25*ne+rk27)

  op_loss = xnmbar*(rk1*o2*rmassinv_o2+rk2*n2*rmassinv_n2+rk10*n2d*rmassinv_n2d)

  q_coeff = q_coeff-op_loss

  explicit = explicit-qopi-(rk19*ne+rk21)*xiop2p-(rk25*ne+rk27)*xiop2d-(rk18*xiop2p+rk24*xiop2d)*o1*rmassinv_o1*xnmbar

  q_coeff(1,:,:) = q_coeff(1,:,:)-p_coeff(1,:,:)
  explicit(1,:,:) = explicit(1,:,:)-2.*p_coeff(1,:,:)*qop(1,:,:)/(1.5*op_loss(1,:,:)-0.5*op_loss(2,:,:))
  p_coeff(1,:,:) = 0.

  call trsolv_ng(p_coeff,q_coeff,r_coeff,explicit,opout,nk-1,i_ng)

  opout(nk,:,:) = opout(nk-1,:,:)**2/opout(nk-2,:,:)

  if (is_bndry(1)) opout(:,-1,:) = flds(i_ng)%op_lon_b(:,1,:,istep)
  if (is_bndry(2)) opout(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%op_lon_b(:,2,:,istep)
  if (is_bndry(3)) opout(:,:,-1) = flds(i_ng)%op_lat_b(:,:,3,istep)
  if (is_bndry(4)) opout(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%op_lat_b(:,:,4,istep)

  opout = max(opout,1.e-5)
  optm1out = max(dtsmooth*op+dtsmooth_div2*(optm1+opout),1.e-5)

  if (enforce_opfloor > 0) then
    do k = 1,nk-1
      do lat = latd0,latd1
        opfloor = opmin*exp(-(glat_ng(i_ng,lat)/90.0)**2/0.3)*exp(-((zpmid_ng(i_ng,k)-4.25)/zpmid_ng(i_ng,nk))**2/0.1)
        opout(k,:,lat) = max(opout(k,:,lat),opfloor)
      enddo
    enddo
  endif

  contains
!-----------------------------------------------------------------------
  subroutine bdotdh(phij,ans)

    real,dimension(nk,flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: phij
    real,dimension(nk,flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: ans

    integer :: k,lat
    real,dimension(nk,flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: dphijdx,dphijdy,ans1,ans2

    call df_3d(phij,dphijdx,1,i_ng)
    call df_3d(phij,dphijdy,-1,i_ng)

    do k = 1,nk
      ans1(k,:,:) = bx*dphijdx(k,:,:)
      ans2(k,:,:) = by*dphijdy(k,:,:)
    enddo
    do lat = latd0,latd1
      ans1(:,:,lat) = ans1(:,:,lat)/cs(lat)
    enddo

    ans = 1./re*(ans1/dlamda(i_ng)+ans2/dphi(i_ng))

  end subroutine bdotdh
!-----------------------------------------------------------------------
end subroutine oplus_ng

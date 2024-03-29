subroutine oplus_ng(op,optm1,opout,optm1out,xiop2p,xiop2d,Fe,Fn,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,zpint_ng
  use cons_module,only: rmass_op,gask,grav,re,rmassinv_o2,rmassinv_o1, &
    rmassinv_n2,rmassinv_n2d,dtsmooth,dtsmooth_div2,pi,rtd,rmassinv_he
  use chemrates_module,only: rk10,rk16,rk17,rk18,rk21,rk22,rk23,rk24,rk26,rk27
  use input_module,only: colfac,opdiffcap,nstep_sub
  use fields_ng_module,only: flds,itp,shapiro,dtx2inv,dlamda,dphi,dlev,dz
  use dffm_ng_module,only: df_2d,df_3d
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: op,optm1
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    opout,optm1out,xiop2p,xiop2d,Fe,Fn

  real,parameter :: explic = 1., phid = 2.0e8, phin = -2.0e8, ppolar = 0.
  integer :: nk,latd0,latd1,k,lat
  real :: gmr
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cs
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    bx,by,bz,bmod2,dipmag,sndec,csdec,chi,rlatm,opflux,dvb,ubca,ubcb,a,fed,fen,cs_by,dbxdx,dcs_bydy,sncsdip
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,te,ti,o2,o1,n2,he,n2d,ne,u,v,w,ui,vi,wi,xnmbar,scht,qop2p,qop2d,qop,rk1,rk2,rk19,rk20,rk25, &
    bdzdvb_op,explicit,hdz,tphdz1,tphdz0,djint,divbz,hdzmbz,hdzpbz,p_coeff,q_coeff,r_coeff,tp1,wd, &
    bdotdh_djbz,op_prod,op_loss_out,dopdt,diffsum,driftsum,windsum,hj,bdotu,bvel,diffj,tp,tr, &
    bdotdh_op,bdotdh_opj,bdotdh_diff,dj,optm1_smooth2,op_loss,vni,djbz,vdotn_h,bdotdh_bvel,diffexp, &
    diffp,diffq,diffr,driftp,driftq,driftr,windp,windq,windr,wni,uii,vii,wii,qop2pi,qop2di,qopi,hdzi, &
    tp_op,dtp_opdz,dbdotdh_opjdz,op_bmod2,dop_bmod2dx,dop_bmod2dy,vdotn_x,vdotn_y
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
  ui = flds(i_ng)%ui
  vi = flds(i_ng)%vi
  wi = flds(i_ng)%wi
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  scht = flds(i_ng)%scht(:,:,:,itp(i_ng))

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
  where (abs(rlatm) >= pi/4.5) a = 1.
  where (abs(rlatm)<pi/4.5 .and. abs(rlatm)>pi/18) &
    a = .5*(1.+cos(abs(rlatm)*6.-pi*4./3.))
  where (abs(rlatm) <= pi/18) a = 0.

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
  vni = 18.1*o2*rmassinv_o2+ &
    o1*rmassinv_o1*sqrt(tr)*(1.-0.064*log10(tr))**2*colfac+ &
    3.6*he*rmassinv_he+18.6*n2*rmassinv_n2
  dj = 1.42E17/(xnmbar*vni)
  vni = 16*3.53E-11*vni
  if (.not. isclose(opdiffcap,0.)) then
    do k = 1,nk
      dj(k,:,:) = min(dj(k,:,:),opdiffcap/(1+2**(zpint_ng(i_ng,k)-6)))
    enddo
  endif

  tp = te+ti

  hj = scht

  do k = 1,nk
    bdotu(k,:,:) = bx*u(k,:,:)+by*v(k,:,:)+hj(k,:,:)*bz*wni(k,:,:)
  enddo

  bvel = bdotu*op

! diffus
  tp_op = tp*op
  do k = 1,nk-2
    dtp_opdz(k+1,:,:) = (tp_op(k+2,:,:)-tp_op(k,:,:))/2.
  enddo
  dtp_opdz(nk-1,:,:) = tp_op(nk-1,:,:)-tp_op(nk-2,:,:)
  dtp_opdz(nk,:,:) = tp_op(nk,:,:)-tp_op(nk-1,:,:)
  dtp_opdz(1,:,:) = tp_op(2,:,:)-tp_op(1,:,:)
  diffj = 1./(hj*dlev(i_ng))*dtp_opdz+rmass_op*grav/gask*op

  tp = tp_op

  call smooth_ng(optm1,optm1_smooth2,shapiro(i_ng)/nstep_sub,i_ng)

  wd = vni*diffj*dj
  sncsdip = sin(dipmag)*cos(dipmag)
  do k = 1,nk
    Fe(k,:,:) = wd(k,:,:)*sncsdip*sndec
    Fn(k,:,:) = wd(k,:,:)*sncsdip*csdec
  enddo

  do k = 1,nk
    djbz(k,:,:) = dj(k,:,:)*bz
  enddo

  call bdotdh(diffj,bdotdh_op)
  bdotdh_op = djbz*bdotdh_op

  call bdotdh(tp,bdotdh_opj)
  bdotdh_opj = bdotdh_opj*dj

  call bdotdh(bdotdh_opj,bdotdh_diff)

! bdzdvb
  do k = 2,nk-2
    dbdotdh_opjdz(k,:,:) = (bdotdh_opj(k+1,:,:)-bdotdh_opj(k-1,:,:))/2.
  enddo
  dbdotdh_opjdz(nk-1,:,:) = bdotdh_opj(nk-1,:,:)-bdotdh_opj(nk-2,:,:)
  dbdotdh_opjdz(nk,:,:) = bdotdh_opj(nk,:,:)-bdotdh_opj(nk-1,:,:)
  dbdotdh_opjdz(1,:,:) = bdotdh_opj(2,:,:)-bdotdh_opj(1,:,:)
  do k = 1,nk
    bdzdvb_op(k,:,:) = bz/(hj(k,:,:)*dz(i_ng))*dbdotdh_opjdz(k,:,:)+dvb*bdotdh_opj(k,:,:)
  enddo

  diffexp = bdzdvb_op+bdotdh_diff+bdotdh_op
  explicit = -explic*diffexp

  call bdotdh(bvel,bdotdh_bvel)

  do k = 1,nk
    op_bmod2(k,:,:) = op(k,:,:)/bmod2**2
  enddo
  call df_3d(op_bmod2,dop_bmod2dx,1,i_ng)
  call df_3d(op_bmod2,dop_bmod2dy,-1,i_ng)
  vdotn_x = uii*dop_bmod2dx
  vdotn_y = vii*dop_bmod2dy
  do lat = latd0,latd1
    vdotn_x(:,:,lat) = vdotn_x(:,:,lat)/cs(lat)
  enddo
  vdotn_h = 1./re*(vdotn_x/dlamda(i_ng)+vdotn_y/dphi(i_ng))
  do k = 1,nk
    vdotn_h(k,:,:) = vdotn_h(k,:,:)*bmod2**2
  enddo

  explicit = explicit+vdotn_h+bdotdh_bvel

  dvb = dvb/bz
  hdz = 1./(hj*dz(i_ng))
  tp1 = 0.5*(ti+te)

  do k = 1,nk-2
    hdzi(k+1,:,:) = 0.5*(hdz(k,:,:)+hdz(k+1,:,:))
  enddo
  hdzi(1,:,:) = 1.5*hdz(1,:,:)-0.5*hdz(2,:,:)
  hdzi(nk,:,:) = 1.5*hdz(nk-1,:,:)-0.5*hdz(nk-2,:,:)

  gmr = grav*rmass_op/(2.*gask)

  tphdz1 = 2.*tp1*hdzi+gmr
  do k = 1,nk-1
    tphdz0(k+1,:,:) = 2.*tp1(k,:,:)*hdzi(k+1,:,:)-gmr
  enddo
  tphdz1(nk,:,:) = 2.*(2.*tp1(nk-1,:,:)-tp1(nk-2,:,:))*hdzi(nk,:,:)+gmr
  tphdz0(1,:,:) = 2.*(2.*tp1(1,:,:)-tp1(2,:,:))*hdzi(1,:,:)-gmr

  do k = 1,nk-2
    djint(k+1,:,:) = 0.5*(dj(k,:,:)+dj(k+1,:,:))
  enddo
  djint(1,:,:) = 1.5*dj(1,:,:)-0.5*dj(2,:,:)
  djint(nk,:,:) = 1.5*dj(nk-1,:,:)-0.5*dj(nk-2,:,:)

  call bdotdh(djbz,bdotdh_djbz)
  do k = 1,nk
    divbz(k,:,:) = dvb+bdotdh_djbz(k,:,:)/(dj(k,:,:)*bz**2)
  enddo

  hdzmbz = hdz-0.5*divbz
  hdzpbz = hdz+0.5*divbz
  do k = 1,nk
    hdzmbz(k,:,:) = hdzmbz(k,:,:)*bz**2
    hdzpbz(k,:,:) = hdzpbz(k,:,:)*bz**2
  enddo

  explicit = explicit-optm1_smooth2*dtx2inv(i_ng)*nstep_sub

  diffp = hdzmbz*djint*tphdz0
  do k = 1,nk-1
    diffq(k,:,:) = -(hdzpbz(k,:,:)*djint(k+1,:,:)*tphdz0(k+1,:,:)+hdzmbz(k,:,:)*djint(k,:,:)*tphdz1(k,:,:))
    diffr(k,:,:) = hdzpbz(k,:,:)*djint(k+1,:,:)*tphdz1(k+1,:,:)
  enddo
  diffq(nk,:,:) = -(hdzpbz(nk,:,:)*(2.*djint(nk,:,:)*tphdz0(nk,:,:)-djint(nk-1,:,:)*tphdz0(nk-1,:,:))+ &
    hdzmbz(nk,:,:)*djint(nk,:,:)*tphdz1(nk,:,:))
  diffr(nk,:,:) = hdzpbz(nk,:,:)*(2.*djint(nk,:,:)*tphdz1(nk,:,:)-djint(nk-1,:,:)*tphdz1(nk-1,:,:))

  driftp = wii*0.5*hdz
  driftq = -wii*6./re
  driftr = -wii*0.5*hdz

  do k = 1,nk-2
    windp(k+1,:,:) = bz*bdotu(k,:,:)*0.5*hdz(k+1,:,:)
    windq(k,:,:) = -bdotu(k,:,:)*dvb*bz
    windr(k,:,:) = -bz*bdotu(k+1,:,:)*0.5*hdz(k,:,:)
  enddo
  windp(1,:,:) = bz*(2.*bdotu(1,:,:)-bdotu(2,:,:))*0.5*hdz(1,:,:)
  windq(nk-1,:,:) = -bdotu(nk-1,:,:)*dvb*bz
  windr(nk-1,:,:) = bz*(2.*bdotu(nk-1,:,:)-bdotu(nk-2,:,:))*0.5*hdz(nk-1,:,:)
  windq(nk,:,:) = -bdotu(nk,:,:)*dvb*bz
  windr(nk,:,:) = bz*(2.*bdotu(nk,:,:)-bdotu(nk-1,:,:))*0.5*hdz(nk,:,:)

  p_coeff = diffp+driftp+windp
  q_coeff = diffq+driftq+windq-dtx2inv(i_ng)*nstep_sub
  r_coeff = diffr+driftr+windr

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

  op_prod = qopi+(rk19*ne+rk21)*xiop2p+(rk25*ne+rk27)*xiop2d+(rk18*xiop2p+rk24*xiop2d)*o1*rmassinv_o1*xnmbar
  call addfld(op_prod,'OP_PROD',i_ng)

  explicit = explicit-op_prod

  q_coeff(1,:,:) = q_coeff(1,:,:)-p_coeff(1,:,:)
  explicit(1,:,:) = explicit(1,:,:)-2.*p_coeff(1,:,:)*qop(1,:,:)/(1.5*op_loss(1,:,:)-0.5*op_loss(2,:,:))
  p_coeff(1,:,:) = 0.

  call trsolv_ng(p_coeff,q_coeff,r_coeff,explicit,opout,nk-1,i_ng)

  opout(nk,:,:) = opout(nk-1,:,:)**2/opout(nk-2,:,:)

  op_loss_out = op_loss*opout
  dopdt = (opout-optm1_smooth2)*dtx2inv(i_ng)*nstep_sub
  do k = 2,nk-1
    diffsum(k,:,:) = diffexp(k,:,:)+diffp(k,:,:)*opout(k-1,:,:)+diffq(k,:,:)*opout(k,:,:)+diffr(k,:,:)*opout(k+1,:,:)
    driftsum(k,:,:) = -vdotn_h(k,:,:)+driftp(k,:,:)*opout(k-1,:,:)+driftq(k,:,:)*opout(k,:,:)+driftr(k,:,:)*opout(k+1,:,:)
    windsum(k,:,:) = -bdotdh_bvel(k,:,:)+windp(k,:,:)*opout(k-1,:,:)+windq(k,:,:)*opout(k,:,:)+windr(k,:,:)*opout(k+1,:,:)
  enddo
  diffsum(1,:,:) = 2.*diffsum(2,:,:)-diffsum(3,:,:)
  driftsum(1,:,:) = 2.*driftsum(2,:,:)-driftsum(3,:,:)
  windsum(1,:,:) = 2.*windsum(2,:,:)-windsum(3,:,:)
  diffsum(nk,:,:) = 2.*diffsum(nk-1,:,:)-diffsum(nk-2,:,:)
  driftsum(nk,:,:) = 2.*driftsum(nk-1,:,:)-driftsum(nk-2,:,:)
  windsum(nk,:,:) = 2.*windsum(nk-1,:,:)-windsum(nk-2,:,:)
  call addfld(-op_loss_out,'OP_LOSS',i_ng)
  call addfld(dopdt,'DOPDT',i_ng)
  call addfld(diffsum,'OP_DIFF',i_ng)
  call addfld(driftsum,'OP_DRIFT',i_ng)
  call addfld(windsum,'OP_WIND',i_ng)

  if (is_bndry(1)) opout(:,-1,:) = flds(i_ng)%op_lon_b(:,1,:,istep)
  if (is_bndry(2)) opout(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%op_lon_b(:,2,:,istep)
  if (is_bndry(3)) opout(:,:,-1) = flds(i_ng)%op_lat_b(:,:,3,istep)
  if (is_bndry(4)) opout(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%op_lat_b(:,:,4,istep)

  optm1out = dtsmooth*op+dtsmooth_div2*(optm1+opout)

  opout = max(opout,0.)
  optm1out = max(optm1out,0.)

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

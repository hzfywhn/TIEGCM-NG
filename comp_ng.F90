subroutine comp_ng(o2_upd,o2nm_upd,o1_upd,o1nm_upd,he_upd,henm_upd,flx_he,istep,i_ng)
! changed the 3d layout to lev,lon,lat (previously lon,lev) in accordance with other subroutines

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,glat_ng,zpmid_ng
  use init_module,only: iday
  use cons_module,only: pi,dtr,rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2, &
    rmass_o2,rmass_o1,rmass_he,dtsmooth,dtsmooth_div2,difhor,grav,p0
  use input_module,only: calc_helium
  use fields_ng_module,only: hor,flds,itp,shapiro,dtx2inv,dz,expzmid,expzmid_inv,expz,difk,bndry
  use lbc,only: b,fb,pshelb
  use matutil_module,only: matinv3
  use char_module,only: find_index
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    o2_upd,o2nm_upd,o1_upd,o1nm_upd,he_upd,henm_upd
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: flx_he

  integer,parameter :: io2 = 1, io1 = 2, ihe = 3
  real,parameter :: tau = 1.86e+3, thdiffalpha = -0.38, t00 = 273., small = 1.e-9
  real,dimension(3),parameter :: ss = (/1.710,1.749,1.718/)
  real,dimension(3,3),parameter :: delta = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
  real,dimension(3,4),parameter :: phi = reshape((/0.,0.673,0.270,1.35,0.,0.404,2.16,1.616,0.,1.11,0.769,0.322/),(/3,4/))
  integer :: nk,lond0,lond1,latd0,latd1,k,i,lat,isp,km,kp,ktmp,m,n,ktop,idx_o2,idx_o1,idx_he
  real :: zmtop,n2top,n2top1,dn2dz
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: dfactor,rlat
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tlbc,wks1,wks3,wks4,embar0,alpha23,alpha22,alpha33,alpha32,flx00
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: fk,wkv1,ps0
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,2) :: ep
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3) :: pk,qk,rk,wkm1,alpha
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3,2) :: ak
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o2_nm,o1,o1_nm,he,he_nm,w,mbar,hdo2,hdo1,hdhe,xnmbar,n2,n2nm,normalize, &
    o2nm_smooth,o1nm_smooth,henm_smooth,o2_advec,o1_advec,he_advec, &
    o2i,o1i,hei,embari,tni,dembardz,dtndz,wi,eddyp,eddyq,eddyr
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: &
    zz,upd,diff_fac,moldif,eddydif,veradv
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3) :: &
    gama,coefp,coefq,coefr
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,0:3) :: fs
  external :: advec_ng,smooth_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o2_nm = flds(i_ng)%o2_nm(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  o1_nm = flds(i_ng)%o1_nm(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  he_nm = flds(i_ng)%he_nm(:,:,:,itp(i_ng))
  w = flds(i_ng)%w(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  hdo2 = flds(i_ng)%hdo2
  hdo1 = flds(i_ng)%hdo1
  hdhe = flds(i_ng)%hdhe
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))

  fs = flds(i_ng)%fs
  tlbc = flds(i_ng)%tlbc

  nk = nlevp1_ng(i_ng)
  lond0 = flds(i_ng)%lond0
  lond1 = flds(i_ng)%lond1
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry

! flx_he is from interpolation
  if (calc_helium == 0) flx_he = 0.

  ps0(:,:,io2) = b(1,1)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:)+fb(1)
  ps0(:,:,io1) = b(2,1)*o2(1,:,:)+b(2,2)*o1(1,:,:)+b(2,3)*he(1,:,:)+fb(2)
  ps0(:,:,ihe) = b(3,1)*o2(1,:,:)+b(3,2)*o1(1,:,:)+b(3,3)*he(1,:,:)+fb(3)
  embar0 = 1./(ps0(:,:,io2)*rmassinv_o2+ps0(:,:,io1)*rmassinv_o1+ps0(:,:,ihe)*rmassinv_he+ &
    (1.-ps0(:,:,io2)-ps0(:,:,io1)-ps0(:,:,ihe))*rmassinv_n2)

  o2i(1,:,:) = .5*(ps0(:,:,io2)+o2(1,:,:))
  o1i(1,:,:) = .5*(ps0(:,:,io1)+o1(1,:,:))
  hei(1,:,:) = .5*(ps0(:,:,ihe)+he(1,:,:))
  embari(1,:,:) = .5*(embar0+mbar(1,:,:))
  tni(1,:,:) = tlbc
  dembardz(1,:,:) = mbar(1,:,:)-embar0
  dtndz(1,:,:) = (tn(1,:,:)-tlbc)*2
  do k = 2,nk
    o2i(k,:,:) = .5*(o2(k-1,:,:)+o2(k,:,:))
    o1i(k,:,:) = .5*(o1(k-1,:,:)+o1(k,:,:))
    hei(k,:,:) = .5*(he(k-1,:,:)+he(k,:,:))
    embari(k,:,:) = .5*(mbar(k-1,:,:)+mbar(k,:,:))
    tni(k,:,:) = .5*(tn(k-1,:,:)+tn(k,:,:))
    dembardz(k,:,:) = mbar(k,:,:)-mbar(k-1,:,:)
    dtndz(k,:,:) = tn(k,:,:)-tn(k-1,:,:)
  enddo
  dembardz = dembardz/dz(i_ng)
  dtndz = dtndz/dz(i_ng)

  do k = 1,nk-1
    wi(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
  enddo
  wi(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)

! advecl
  call advec_ng(o2,o2_advec,i_ng)
  call advec_ng(o1,o1_advec,i_ng)
  call advec_ng(he,he_advec,i_ng)

  call addfld(o2_advec,'O2_HORADV',i_ng)
  call addfld(o1_advec,'O1_HORADV',i_ng)
  call addfld(he_advec,'HE_HORADV',i_ng)
  call addfld(hdo2,'O2_HORDIF',i_ng)
  call addfld(hdo1,'O1_HORDIF',i_ng)
  call addfld(hdhe,'HE_HORDIF',i_ng)

  call smooth_ng(o2_nm,o2nm_smooth,shapiro(i_ng),i_ng)
  call smooth_ng(o1_nm,o1nm_smooth,shapiro(i_ng),i_ng)
  call smooth_ng(he_nm,henm_smooth,shapiro(i_ng),i_ng)

  if (difhor > 0) then
    rlat = glat_ng(i_ng,latd0:latd1)*dtr
    where (abs(rlat) >= pi/4.5)
      dfactor = hor+1.
    elsewhere
      dfactor = hor+.5*(1.+sin(pi*(abs(rlat)-pi/9.)/(pi/4.5)))
    endwhere
  else
    dfactor = 1.
  endif

  do n = 1,3
    diff_fac(:,:,:,n) = (tni/t00)**(1.75-ss(n))
  enddo

  wks4 = dembardz(1,:,:)/(embari(1,:,:)*2.)

  km = 1
  kp = 2

  ep(:,:,io2,kp) = 1.-1./embari(1,:,:)*(rmass_o2+dembardz(1,:,:))
  ep(:,:,io1,kp) = 1.-1./embari(1,:,:)*(rmass_o1+dembardz(1,:,:))
  ep(:,:,ihe,kp) = 1.-1./embari(1,:,:)*(rmass_he+dembardz(1,:,:))- &
    thdiffalpha*dtndz(1,:,:)/tni(1,:,:)
  zz(1,:,:,:) = 0.

  alpha(:,:,io2,1) = -(phi(io2,4)+(phi(io2,io1)-phi(io2,4))*o1i(1,:,:)+ &
    (diff_fac(1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*hei(1,:,:))
  alpha(:,:,io1,2) = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(1,:,:)+ &
    (diff_fac(1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(1,:,:))
  alpha(:,:,ihe,3) = -(diff_fac(1,:,:,3)*phi(ihe,4)+ &
    (diff_fac(1,:,:,io2)*phi(ihe,io2)-diff_fac(1,:,:,3)*phi(ihe,4))*o2i(1,:,:)+ &
    (diff_fac(1,:,:,io1)*phi(ihe,io1)-diff_fac(1,:,:,3)*phi(ihe,4))*o1i(1,:,:))
  alpha(:,:,io2,2) = (phi(io2,io1)-phi(io2,4))*o2i(1,:,:)
  alpha(:,:,io2,3) = (diff_fac(1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*o2i(1,:,:)
  alpha(:,:,io1,1) = (phi(io1,io2)-phi(io1,4))*o1i(1,:,:)
  alpha(:,:,io1,3) = (diff_fac(1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(1,:,:)
  alpha(:,:,ihe,1) = (diff_fac(1,:,:,io2)*phi(ihe,io2)-diff_fac(1,:,:,3)*phi(ihe,4))*hei(1,:,:)
  alpha(:,:,ihe,2) = (diff_fac(1,:,:,io1)*phi(ihe,io1)-diff_fac(1,:,:,3)*phi(ihe,4))*hei(1,:,:)

  do lat = latd0,latd1
    do i = lond0,lond1
      ak(i,lat,:,:,kp) = matinv3(alpha(i,lat,:,:))
    enddo
  enddo

  wks1 = embari(1,:,:)*rmassinv_n2*(t00/tlbc)**0.25/tau

  do m = 1,3
    do isp = io2,ihe
      ak(:,:,isp,m,kp) = ak(:,:,isp,m,kp)*wks1
    enddo
  enddo
  gama(1,:,:,:,:) = 0.

  km = 1
  kp = 2
  do k = 1,nk-1
    ktmp = km
    km = kp
    kp = ktmp

    ep(:,:,io2,kp) = 1.-1./embari(k+1,:,:)*(rmass_o2+dembardz(k+1,:,:))
    ep(:,:,io1,kp) = 1.-1./embari(k+1,:,:)*(rmass_o1+dembardz(k+1,:,:))
    ep(:,:,ihe,kp) = 1.-1./embari(k+1,:,:)*(rmass_he+dembardz(k+1,:,:))
    if (k /= nk-1) ep(:,:,ihe,kp) = ep(:,:,ihe,kp)-thdiffalpha*dtndz(k+1,:,:)/tni(k+1,:,:)

    alpha(:,:,io2,1) = -(phi(io2,4)+(phi(io2,io1)-phi(io2,4))*o1i(k+1,:,:)+ &
      (diff_fac(k+1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*hei(k+1,:,:))
    alpha(:,:,io1,2) = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(k+1,:,:)+ &
      (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(k+1,:,:))
    alpha(:,:,ihe,3) = -(diff_fac(k+1,:,:,3)*phi(ihe,4)+ &
      (diff_fac(k+1,:,:,io2)*phi(ihe,io2)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o2i(k+1,:,:)+ &
      (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o1i(k+1,:,:))
    alpha(:,:,io2,2) = (phi(io2,io1)-phi(io2,4))*o2i(k+1,:,:)
    alpha(:,:,io2,3) = (diff_fac(k+1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*o2i(k+1,:,:)
    alpha(:,:,io1,1) = (phi(io1,io2)-phi(io1,4))*o1i(k+1,:,:)
    alpha(:,:,io1,3) = (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(k+1,:,:)
    alpha(:,:,ihe,1) = (diff_fac(k+1,:,:,io2)*phi(ihe,io2)-diff_fac(k+1,:,:,3)*phi(ihe,4))*hei(k+1,:,:)
    alpha(:,:,ihe,2) = (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*hei(k+1,:,:)

    do lat = latd0,latd1
      do i = lond0,lond1
        ak(i,lat,:,:,kp) = matinv3(alpha(i,lat,:,:))
      enddo
    enddo

    wks1 = embari(k+1,:,:)*rmassinv_n2*(t00/tni(k+1,:,:))**0.25/tau

    wks3 = wks4
    wks4 = dembardz(k+1,:,:)/(embari(k+1,:,:)*2.)

    eddyp(k,:,:) = expzmid_inv(i_ng)*difk(i_ng,k,iday)*(1./dz(i_ng)-wks3)
    eddyq(k,:,:) = &
      expzmid    (i_ng)*difk(i_ng,k+1,iday)*(1./dz(i_ng)-wks4)+ &
      expzmid_inv(i_ng)*difk(i_ng,k  ,iday)*(1./dz(i_ng)+wks3)
    eddyr(k,:,:) = expzmid(i_ng)*difk(i_ng,k+1,iday)*(1./dz(i_ng)+wks4)

    do m = 1,3
      do isp = io2,ihe
        ak(:,:,isp,m,kp) = ak(:,:,isp,m,kp)*wks1
      enddo
    enddo

    do isp = io2,ihe
      coefp(k,:,:,isp,:) = -ak(:,:,isp,:,km)*(1./dz(i_ng)+ep(:,:,:,km)/2.)
      coefq(k,:,:,isp,:) = &
        ak(:,:,isp,:,km)*(1./dz(i_ng)-ep(:,:,:,km)/2.)+ &
        ak(:,:,isp,:,kp)*(1./dz(i_ng)+ep(:,:,:,kp)/2.)
      coefr(k,:,:,isp,:) = -ak(:,:,isp,:,kp)*(1./dz(i_ng)-ep(:,:,:,kp)/2.)
    enddo

    do m = 1,3
      do isp = io2,ihe
        do lat = latd0,latd1
          pk(:,lat,isp,m) = -coefp(k,:,lat,isp,m)/dz(i_ng)- &
            expz(i_ng,k)*delta(isp,m)/dz(i_ng)* &
            (eddyp(k,:,lat)*dfactor(lat)+.5*wi(k,:,lat))
          rk(:,lat,isp,m) = -coefr(k,:,lat,isp,m)/dz(i_ng)- &
            expz(i_ng,k)*delta(isp,m)/dz(i_ng)* &
            (eddyr(k,:,lat)*dfactor(lat)-.5*wi(k,:,lat))
          qk(:,lat,isp,m) = -coefq(k,:,lat,isp,m)/dz(i_ng)+ &
            expz(i_ng,k)*delta(isp,m)*(eddyq(k,:,lat)*dfactor(lat)/dz(i_ng)+dtx2inv(i_ng))- &
            expz(i_ng,k)*fs(k,:,lat,isp,m)
        enddo
      enddo
    enddo

    fk(:,:,io2) = expz(i_ng,k)*(o2nm_smooth(k,:,:)*dtx2inv(i_ng)-(o2_advec(k,:,:)-fs(k,:,:,io2,0))+hdo2(k,:,:))
    fk(:,:,io1) = expz(i_ng,k)*(o1nm_smooth(k,:,:)*dtx2inv(i_ng)-(o1_advec(k,:,:)-fs(k,:,:,io1,0))+hdo1(k,:,:))
    fk(:,:,ihe) = expz(i_ng,k)*(henm_smooth(k,:,:)*dtx2inv(i_ng)-(he_advec(k,:,:)-fs(k,:,:,ihe,0))+hdhe(k,:,:))

    if (k == 1) then
      do m = 1,3
        do n = 1,3
          do isp = io2,ihe
            qk(:,:,isp,m) = qk(:,:,isp,m)+pk(:,:,isp,n)*b(n,m)
          enddo
        enddo
      enddo
      do m = 1,3
        do isp = io2,ihe
          fk(:,:,isp) = fk(:,:,isp)-pk(:,:,isp,m)*fb(m)
        enddo
      enddo
      pk = 0.
    elseif (k == nk-1) then
      do m = 1,3
        do isp = io2,ihe
          qk(:,:,isp,m) = qk(:,:,isp,m)+(1.+.5*ep(:,:,m,kp)*dz(i_ng))/(1.-.5*ep(:,:,m,kp)*dz(i_ng))*rk(:,:,isp,m)
        enddo
      enddo
      alpha22 = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(k+1,:,:)+ &
        (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(k+1,:,:))
      alpha33 = -(diff_fac(k+1,:,:,3)*phi(ihe,4)+ &
        (diff_fac(k+1,:,:,io2)*phi(ihe,io2)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o2i(k+1,:,:)+ &
        (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o1i(k+1,:,:))
      alpha23 = (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(k+1,:,:)
      alpha32 = (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*hei(k+1,:,:)
      flx00 = embari(k,:,:)*rmassinv_n2*p0*(t00/tn(k,:,:))**0.25/(tau*grav)
      do isp = io2,ihe
        fk(:,:,isp) = fk(:,:,isp)- &
          rk(:,:,isp,2)*(alpha23-alpha22)*flx_he/(flx00*(1./dz(i_ng)-0.5*ep(:,:,2,kp)))- &
          rk(:,:,isp,3)*(alpha33-alpha32)*flx_he/(flx00*(1./dz(i_ng)-0.5*ep(:,:,3,kp)))
      enddo
      rk = 0.
    endif

    do m = 1,3
      do n = 1,3
        do isp = io2,ihe
          qk(:,:,isp,m) = qk(:,:,isp,m)-pk(:,:,isp,n)*gama(k,:,:,n,m)
        enddo
      enddo
    enddo

    do lat = latd0,latd1
      do i = lond0,lond1
        wkm1(i,lat,:,:) = matinv3(qk(i,lat,:,:))
      enddo
    enddo

    wkv1 = fk
    gama(k+1,:,:,:,:) = 0.

    do m = 1,3
      do isp = io2,ihe
        wkv1(:,:,isp) = wkv1(:,:,isp)-pk(:,:,isp,m)*zz(k,:,:,m)
      enddo
      do n = 1,3
        do isp = io2,ihe
          gama(k+1,:,:,isp,m) = gama(k+1,:,:,isp,m)+wkm1(:,:,isp,n)*rk(:,:,n,m)
        enddo
      enddo
    enddo

    zz(k+1,:,:,:) = 0.
    do m = 1,3
      do isp = io2,ihe
        zz(k+1,:,:,isp) = zz(k+1,:,:,isp)+wkm1(:,:,isp,m)*wkv1(:,:,m)
      enddo
    enddo
  enddo

  upd(nk,:,:,:) = 0.

  do k = nk-1,1,-1
    upd(k,:,:,:) = zz(k+1,:,:,:)
    do isp = io2,ihe
      do m = 1,3
        upd(k,:,:,isp) = upd(k,:,:,isp)-gama(k+1,:,:,isp,m)*upd(k+1,:,:,m)
      enddo
    enddo
  enddo

  do isp = io2,ihe
    upd(nk,:,:,isp) = (1.+.5*ep(:,:,isp,kp)*dz(i_ng))/(1.-.5*ep(:,:,isp,kp)*dz(i_ng))*upd(nk-1,:,:,isp)
  enddo
  upd(nk,:,:,io1) = upd(nk,:,:,io1)+(alpha23-alpha22)*flx_he/(flx00*(1./dz(i_ng)-0.5*ep(:,:,io1,kp)))
  upd(nk,:,:,ihe) = upd(nk,:,:,ihe)+(alpha33-alpha32)*flx_he/(flx00*(1./dz(i_ng)-0.5*ep(:,:,ihe,kp)))

  moldif = 0.
  eddydif = 0.
  veradv = 0.
  do k = 2,nk-1
    do isp = io2,ihe
      do m = 1,3
        moldif(k,:,:,isp) = moldif(k,:,:,isp)+ &
          coefp(k,:,:,isp,m)*upd(k-1,:,:,m)+ &
          coefq(k,:,:,isp,m)*upd(k  ,:,:,m)+ &
          coefr(k,:,:,isp,m)*upd(k+1,:,:,m)
      enddo
      eddydif(k,:,:,isp) = &
        eddyp(k,:,:)*upd(k-1,:,:,isp)+ &
        eddyq(k,:,:)*upd(k  ,:,:,isp)+ &
        eddyr(k,:,:)*upd(k+1,:,:,isp)
      veradv(k,:,:,isp) = wi(k,:,:)*(upd(k-1,:,:,isp)-upd(k+1,:,:,isp))
    enddo
    moldif(k,:,:,:) = moldif(k,:,:,:)/dz(i_ng)/expz(i_ng,k)
    do lat = latd0,latd1
      eddydif(k,:,lat,:) = eddydif(k,:,lat,:)*dfactor(lat)/expz(i_ng,k)
    enddo
    veradv(k,:,:,:) = veradv(k,:,:,:)/(2*dz(i_ng))
  enddo
  call addfld(moldif(:,:,:,io2),'O2_MOLDIF',i_ng)
  call addfld(moldif(:,:,:,io1),'O1_MOLDIF',i_ng)
  call addfld(moldif(:,:,:,ihe),'HE_MOLDIF',i_ng)
  call addfld(eddydif(:,:,:,io2),'O2_EDDYDIF',i_ng)
  call addfld(eddydif(:,:,:,io1),'O1_EDDYDIF',i_ng)
  call addfld(eddydif(:,:,:,ihe),'HE_EDDYDIF',i_ng)
  call addfld(veradv(:,:,:,io2),'O2_VERADV',i_ng)
  call addfld(veradv(:,:,:,io1),'O1_VERADV',i_ng)
  call addfld(veradv(:,:,:,ihe),'HE_VERADV',i_ng)

  o2_upd = upd(:,:,:,io2)
  o1_upd = upd(:,:,:,io1)
  he_upd = upd(:,:,:,ihe)

  if (calc_helium == 0) then
    he_upd = pshelb
    henm_upd = pshelb
  endif

  idx_o2 = find_index('O2',bndry)
  idx_o1 = find_index('O1',bndry)
  idx_he = find_index('HE',bndry)
  if (is_bndry(1)) then
    o2_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_o2)
    o1_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_o1)
    he_upd(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_he)
  endif
  if (is_bndry(2)) then
    o2_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_o2)
    o1_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_o1)
    he_upd(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_he)
  endif
  if (is_bndry(3)) then
    o2_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_o2)
    o1_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_o1)
    he_upd(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_he)
  endif
  if (is_bndry(4)) then
    o2_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_o2)
    o1_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_o1)
    he_upd(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_he)
  endif

  o2nm_upd = dtsmooth_div2*(o2_nm+o2_upd)+dtsmooth*o2
  o1nm_upd = dtsmooth_div2*(o1_nm+o1_upd)+dtsmooth*o1
  henm_upd = dtsmooth_div2*(he_nm+he_upd)+dtsmooth*he

  n2 = 1.-o2_upd-o1_upd-he_upd
  n2nm = 1.-o2nm_upd-o1nm_upd-henm_upd

  do lat = latd0,latd1
    do i = lond0,lond1
      do ktop = 1,nk
        if (n2(ktop,i,lat) < 0.01) exit
      enddo
      if (ktop < nk) then
        zmtop = zpmid_ng(i_ng,ktop)
        n2top = log(n2(ktop,i,lat)*xnmbar(ktop,i,lat))
        n2top1 = log(n2(ktop-1,i,lat)*xnmbar(ktop-1,i,lat))
        dn2dz = (n2top-n2top1)/dz(i_ng)
        do k = ktop+1,nk
          n2(k,i,lat) = exp(n2top+dn2dz*(zpmid_ng(i_ng,k)-zmtop))/xnmbar(k,i,lat)
        enddo
      endif
    enddo
  enddo

  o2_upd = max(o2_upd,small)
  o1_upd = max(o1_upd,small)
  he_upd = max(he_upd,small)
  n2 = max(n2,small)

  o2nm_upd = max(o2nm_upd,small)
  o1nm_upd = max(o1nm_upd,small)
  henm_upd = max(henm_upd,small)
  n2nm = max(n2nm,small)

  normalize = o2_upd+o1_upd+he_upd+n2
  o2_upd = o2_upd/normalize
  o1_upd = o1_upd/normalize
  he_upd = he_upd/normalize

  normalize = o2nm_upd+o1nm_upd+henm_upd+n2nm
  o2nm_upd = o2nm_upd/normalize
  o1nm_upd = o1nm_upd/normalize
  henm_upd = henm_upd/normalize

end subroutine comp_ng

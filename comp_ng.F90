subroutine comp_ng(n2,o2_upd,o2nm_upd,o1_upd,o1nm_upd,he_upd,henm_upd,flx_he,istep,i_ng)
! changed the 3d layout to lev,lon,lat (previously lon,lev) in accordance with other subroutines

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,glat_ng
  use init_module,only: iday
  use cons_module,only: pi,dtr,rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2, &
    rmass_o2,rmass_o1,rmass_he,dtsmooth,dtsmooth_div2,difhor,grav,p0
  use input_module,only: calc_helium
  use fields_ng_module,only: hor,b,fb,flds,itp,shapiro,dtx2inv,dz,expzmid,expzmid_inv,expz,difk,bndry
  use char_module,only: find_index
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: n2
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    o2_upd,o2nm_upd,o1_upd,o1nm_upd,he_upd,henm_upd
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: flx_he

  integer,parameter :: io2 = 1, io1 = 2, ihe = 3
  real,parameter :: tau = 1.86e+3, thdiffalpha = -0.38, t00 = 273., small = 1.e-6, small_he = 1.e-9
  real,dimension(3),parameter :: ss = (/1.710,1.749,1.718/)
  real,dimension(3,3),parameter :: delta = reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
  real,dimension(3,4),parameter :: phi = reshape((/0.,0.673,0.270,1.35,0.,0.404,2.16,1.616,0.,1.11,0.769,0.322/),(/3,4/))
  integer :: nk,latd0,latd1,k,lat,isp,km,kp,ktmp,m,n,idx_o2,idx_o1,idx_he
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: dfactor,rlat
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tlbc,wks1,wks3,wks4,embar0,detalpha,alpha23,alpha22,alpha33,alpha32,flx00
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: fk,wkv1,ps0
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,2) :: ep
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3) :: pk,qk,rk,wkm1,alpha
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3,2) :: ak
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o2_nm,o1,o1_nm,he,he_nm,w,mbar,hdo2,hdo1,hdhe,normalize, &
    o2nm_smooth,o1nm_smooth,henm_smooth,o2_advec,o1_advec,he_advec,wi
  real,dimension(0:nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: o2i,o1i,hei,embari,tni
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3) :: zz,upd,diff_fac
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3,3) :: gama
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

  fs = flds(i_ng)%fs
  tlbc = flds(i_ng)%tlbc

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry

! flx_he is from interpolation
  if (calc_helium /= 1) then
    he_upd = 0.
    henm_upd = 0.
    flx_he = 0.
  endif

  ps0(:,:,io2) = b(1,1)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:)+fb(1)
  ps0(:,:,io1) = b(2,1)*o2(1,:,:)+b(2,2)*o1(1,:,:)+b(2,3)*he(1,:,:)+fb(2)
  ps0(:,:,ihe) = b(3,1)*o2(1,:,:)+b(3,2)*o1(1,:,:)+b(3,3)*he(1,:,:)+fb(3)
  embar0 = 1./(ps0(:,:,io2)*rmassinv_o2+ps0(:,:,io1)*rmassinv_o1+ps0(:,:,ihe)*rmassinv_he+ &
    (1.-ps0(:,:,io2)-ps0(:,:,io1)-ps0(:,:,ihe))*rmassinv_n2)

  o2i(0,:,:) = .5*(ps0(:,:,io2)+o2(1,:,:))
  o1i(0,:,:) = .5*(ps0(:,:,io1)+o1(1,:,:))
  hei(0,:,:) = .5*(ps0(:,:,ihe)+he(1,:,:))
  embari(0,:,:) = .5*(embar0+mbar(1,:,:))
  tni(0,:,:) = .5*(tlbc+tn(1,:,:))
  do k = 1,nk-1
    o2i(k,:,:) = .5*(o2(k,:,:)+o2(k+1,:,:))
    o1i(k,:,:) = .5*(o1(k,:,:)+o1(k+1,:,:))
    hei(k,:,:) = .5*(he(k,:,:)+he(k+1,:,:))
    embari(k,:,:) = .5*(mbar(k,:,:)+mbar(k+1,:,:))
    tni(k,:,:) = .5*(tn(k,:,:)+tn(k+1,:,:))
    wi(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
  enddo
  o2i(nk,:,:) = 1.5*o2(nk,:,:)-.5*o2(nk-1,:,:)
  o1i(nk,:,:) = 1.5*o1(nk,:,:)-.5*o1(nk-1,:,:)
  hei(nk,:,:) = 1.5*he(nk,:,:)-.5*he(nk-1,:,:)
  embari(nk,:,:) = 1.5*mbar(nk,:,:)-.5*mbar(nk-1,:,:)
  tni(nk,:,:) = 1.5*tn(nk,:,:)-.5*tn(nk-1,:,:)
  wi(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)

  call advec_ng(o2,o2_advec,i_ng)
  call advec_ng(o1,o1_advec,i_ng)
  call advec_ng(he,he_advec,i_ng)

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
    diff_fac(nk,:,:,n) = (tn(nk-1,:,:)/t00)**(1.75-ss(n))
    do k = 1,nk-1
      diff_fac(k,:,:,n) = (tni(k-1,:,:)/t00)**(1.75-ss(n))
    enddo
  enddo

  wks4 = (mbar(1,:,:)-embar0)/(dz(i_ng)*embari(0,:,:)*2.)

  km = 1
  kp = 2

  ep(:,:,io2,kp) = 1.-1./embari(0,:,:)*(rmass_o2+(mbar(1,:,:)-embar0)/dz(i_ng))
  ep(:,:,io1,kp) = 1.-1./embari(0,:,:)*(rmass_o1+(mbar(1,:,:)-embar0)/dz(i_ng))
  ep(:,:,ihe,kp) = 1.-1./embari(0,:,:)*(rmass_he+(mbar(1,:,:)-embar0)/dz(i_ng))- &
    thdiffalpha*(tn(1,:,:)-tlbc)/(dz(i_ng)*tni(0,:,:))
  zz(1,:,:,:) = 0.

  alpha(:,:,io2,1) = -(phi(io2,4)+(phi(io2,io1)-phi(io2,4))*o1i(0,:,:)+ &
    (diff_fac(1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*hei(0,:,:))
  alpha(:,:,io1,2) = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(0,:,:)+ &
    (diff_fac(1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(0,:,:))
  alpha(:,:,ihe,3) = -(diff_fac(1,:,:,3)*phi(ihe,4)+ &
    (diff_fac(1,:,:,io2)*phi(ihe,io2)-diff_fac(1,:,:,3)*phi(ihe,4))*o2i(0,:,:)+ &
    (diff_fac(1,:,:,io1)*phi(ihe,io1)-diff_fac(1,:,:,3)*phi(ihe,4))*o1i(0,:,:))
  alpha(:,:,io2,2) = (phi(io2,io1)-phi(io2,4))*o2i(0,:,:)
  alpha(:,:,io2,3) = (diff_fac(1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*o2i(0,:,:)
  alpha(:,:,io1,1) = (phi(io1,io2)-phi(io1,4))*o1i(0,:,:)
  alpha(:,:,io1,3) = (diff_fac(1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(0,:,:)
  alpha(:,:,ihe,1) = (diff_fac(1,:,:,io2)*phi(ihe,io2)-diff_fac(1,:,:,3)*phi(ihe,4))*hei(0,:,:)
  alpha(:,:,ihe,2) = (diff_fac(1,:,:,io1)*phi(ihe,io1)-diff_fac(1,:,:,3)*phi(ihe,4))*hei(0,:,:)

  detalpha = alpha(:,:,1,1)*(alpha(:,:,2,2)*alpha(:,:,3,3)-alpha(:,:,2,3)*alpha(:,:,3,2))+ &
    alpha(:,:,1,2)*(alpha(:,:,2,3)*alpha(:,:,3,1)-alpha(:,:,2,1)*alpha(:,:,3,3))+ &
    alpha(:,:,1,3)*(alpha(:,:,2,1)*alpha(:,:,3,2)-alpha(:,:,2,2)*alpha(:,:,3,1))

  ak(:,:,1,1,kp) = alpha(:,:,2,2)*alpha(:,:,3,3)-alpha(:,:,2,3)*alpha(:,:,3,2)
  ak(:,:,2,2,kp) = alpha(:,:,1,1)*alpha(:,:,3,3)-alpha(:,:,1,3)*alpha(:,:,3,1)
  ak(:,:,3,3,kp) = alpha(:,:,1,1)*alpha(:,:,2,2)-alpha(:,:,1,2)*alpha(:,:,2,1)
  ak(:,:,1,2,kp) = alpha(:,:,1,3)*alpha(:,:,3,2)-alpha(:,:,1,2)*alpha(:,:,3,3)
  ak(:,:,1,3,kp) = alpha(:,:,1,2)*alpha(:,:,2,3)-alpha(:,:,1,3)*alpha(:,:,2,2)
  ak(:,:,2,3,kp) = alpha(:,:,1,3)*alpha(:,:,2,1)-alpha(:,:,1,1)*alpha(:,:,2,3)
  ak(:,:,2,1,kp) = alpha(:,:,2,3)*alpha(:,:,3,1)-alpha(:,:,2,1)*alpha(:,:,3,3)
  ak(:,:,3,1,kp) = alpha(:,:,2,1)*alpha(:,:,3,2)-alpha(:,:,2,2)*alpha(:,:,3,1)
  ak(:,:,3,2,kp) = alpha(:,:,3,1)*alpha(:,:,1,2)-alpha(:,:,1,1)*alpha(:,:,3,2)

  wks1 = embari(0,:,:)*rmassinv_n2*(t00/tlbc)**0.25/(tau*detalpha)

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

    ep(:,:,io2,kp) = 1.-1./embari(k,:,:)*(rmass_o2+(mbar(k+1,:,:)-mbar(k,:,:))/dz(i_ng))
    ep(:,:,io1,kp) = 1.-1./embari(k,:,:)*(rmass_o1+(mbar(k+1,:,:)-mbar(k,:,:))/dz(i_ng))
    ep(:,:,ihe,kp) = 1.-1./embari(k,:,:)*(rmass_he+(mbar(k+1,:,:)-mbar(k,:,:))/dz(i_ng))
    if (k /= nk-1) ep(:,:,ihe,kp) = ep(:,:,ihe,kp)-thdiffalpha*(tn(k+1,:,:)-tn(k,:,:))/(dz(i_ng)*tni(k,:,:))

    alpha(:,:,io2,1) = -(phi(io2,4)+(phi(io2,io1)-phi(io2,4))*o1i(k,:,:)+ &
      (diff_fac(k+1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*hei(k,:,:))
    alpha(:,:,io1,2) = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(k,:,:)+ &
      (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(k,:,:))
    alpha(:,:,ihe,3) = -(diff_fac(k+1,:,:,3)*phi(ihe,4)+ &
      (diff_fac(k+1,:,:,io2)*phi(ihe,io2)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o2i(k,:,:)+ &
      (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*o1i(k,:,:))
    alpha(:,:,io2,2) = (phi(io2,io1)-phi(io2,4))*o2i(k,:,:)
    alpha(:,:,io2,3) = (diff_fac(k+1,:,:,io2)*phi(io2,ihe)-phi(io2,4))*o2i(k,:,:)
    alpha(:,:,io1,1) = (phi(io1,io2)-phi(io1,4))*o1i(k,:,:)
    alpha(:,:,io1,3) = (diff_fac(k+1,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(k,:,:)
    alpha(:,:,ihe,1) = (diff_fac(k+1,:,:,io2)*phi(ihe,io2)-diff_fac(k+1,:,:,3)*phi(ihe,4))*hei(k,:,:)
    alpha(:,:,ihe,2) = (diff_fac(k+1,:,:,io1)*phi(ihe,io1)-diff_fac(k+1,:,:,3)*phi(ihe,4))*hei(k,:,:)

    detalpha = alpha(:,:,1,1)*(alpha(:,:,2,2)*alpha(:,:,3,3)-alpha(:,:,2,3)*alpha(:,:,3,2))+ &
      alpha(:,:,1,2)*(alpha(:,:,2,3)*alpha(:,:,3,1)-alpha(:,:,2,1)*alpha(:,:,3,3))+ &
      alpha(:,:,1,3)*(alpha(:,:,2,1)*alpha(:,:,3,2)-alpha(:,:,2,2)*alpha(:,:,3,1))

    ak(:,:,1,1,kp) = alpha(:,:,2,2)*alpha(:,:,3,3)-alpha(:,:,2,3)*alpha(:,:,3,2)
    ak(:,:,2,2,kp) = alpha(:,:,1,1)*alpha(:,:,3,3)-alpha(:,:,1,3)*alpha(:,:,3,1)
    ak(:,:,3,3,kp) = alpha(:,:,1,1)*alpha(:,:,2,2)-alpha(:,:,1,2)*alpha(:,:,2,1)
    ak(:,:,1,2,kp) = alpha(:,:,1,3)*alpha(:,:,3,2)-alpha(:,:,1,2)*alpha(:,:,3,3)
    ak(:,:,1,3,kp) = alpha(:,:,1,2)*alpha(:,:,2,3)-alpha(:,:,1,3)*alpha(:,:,2,2)
    ak(:,:,2,3,kp) = alpha(:,:,1,3)*alpha(:,:,2,1)-alpha(:,:,1,1)*alpha(:,:,2,3)
    ak(:,:,2,1,kp) = alpha(:,:,2,3)*alpha(:,:,3,1)-alpha(:,:,2,1)*alpha(:,:,3,3)
    ak(:,:,3,1,kp) = alpha(:,:,2,1)*alpha(:,:,3,2)-alpha(:,:,2,2)*alpha(:,:,3,1)
    ak(:,:,3,2,kp) = alpha(:,:,3,1)*alpha(:,:,1,2)-alpha(:,:,1,1)*alpha(:,:,3,2)

    if (k == nk-1) then
      wks1 = embari(nk-1,:,:)*rmassinv_n2*(t00/tn(nk-1,:,:))**0.25/(tau*detalpha)
    else
      wks1 = embari(k,:,:)*rmassinv_n2*(t00/tni(k,:,:))**0.25/(tau*detalpha)
    endif

    wks3 = wks4
    wks4 = (mbar(k+1,:,:)-mbar(k,:,:))/(dz(i_ng)*embari(k,:,:)*2.)

    do m = 1,3
      do isp = io2,ihe
        ak(:,:,isp,m,kp) = ak(:,:,isp,m,kp)*wks1
        do lat = latd0,latd1
          pk(:,lat,isp,m) = (ak(:,lat,isp,m,km)*(1./dz(i_ng)+ep(:,lat,m,km)/2.)- &
            expz(i_ng,k)*(expzmid_inv(i_ng)*difk(i_ng,k,iday)*dfactor(lat)*(1./dz(i_ng)-wks3(:,lat))+.5*wi(k,:,lat))* &
            delta(isp,m))/dz(i_ng)
          rk(:,lat,isp,m) = (ak(:,lat,isp,m,kp)*(1./dz(i_ng)-ep(:,lat,m,kp)/2.)- &
            expz(i_ng,k)*(expzmid(i_ng)*difk(i_ng,k+1,iday)*dfactor(lat)*(1./dz(i_ng)+wks4(:,lat))-.5*wi(k,:,lat))* &
            delta(isp,m))/dz(i_ng)
          qk(:,lat,isp,m) = -(ak(:,lat,isp,m,km)*(1./dz(i_ng)-ep(:,lat,m,km)/2.)+ &
            ak(:,lat,isp,m,kp)*(1./dz(i_ng)+ep(:,lat,m,kp)/2.))/dz(i_ng)+ &
            expz(i_ng,k)*(((expzmid(i_ng)*difk(i_ng,k+1,iday)*(1./dz(i_ng)-wks4(:,lat))+ &
            expzmid_inv(i_ng)*difk(i_ng,k,iday)*(1./dz(i_ng)+wks3(:,lat)))* &
            dfactor(lat)/dz(i_ng)+dtx2inv(i_ng))*delta(isp,m)-fs(k,:,lat,isp,m))
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
      alpha22 = -(phi(io1,4)+(phi(io1,io2)-phi(io1,4))*o2i(nk-1,:,:)+ &
        (diff_fac(nk,:,:,io1)*phi(io1,ihe)-phi(io1,4))*hei(nk-1,:,:))
      alpha33 = -(diff_fac(nk,:,:,3)*phi(ihe,4)+ &
        (diff_fac(nk,:,:,io2)*phi(ihe,io2)-diff_fac(nk,:,:,3)*phi(ihe,4))*o2i(nk-1,:,:)+ &
        (diff_fac(nk,:,:,io1)*phi(ihe,io1)-diff_fac(nk,:,:,3)*phi(ihe,4))*o1i(nk-1,:,:))
      alpha23 = (diff_fac(nk,:,:,io1)*phi(io1,ihe)-phi(io1,4))*o1i(nk-1,:,:)
      alpha32 = (diff_fac(nk,:,:,io1)*phi(ihe,io1)-diff_fac(nk,:,:,3)*phi(ihe,4))*hei(nk-1,:,:)
      flx00 = embari(nk-1,:,:)*rmassinv_n2*p0*(t00/tn(nk-1,:,:))**0.25/(tau*grav)
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

    wks1 = qk(:,:,1,1)*(qk(:,:,2,2)*qk(:,:,3,3)-qk(:,:,2,3)*qk(:,:,3,2))+ &
      qk(:,:,1,2)*(qk(:,:,2,3)*qk(:,:,3,1)-qk(:,:,2,1)*qk(:,:,3,3))+ &
      qk(:,:,1,3)*(qk(:,:,2,1)*qk(:,:,3,2)-qk(:,:,2,2)*qk(:,:,3,1))
    wkm1(:,:,io2,1) = (qk(:,:,2,2)*qk(:,:,3,3)-qk(:,:,2,3)*qk(:,:,3,2))/wks1
    wkm1(:,:,io2,2) = (qk(:,:,1,3)*qk(:,:,3,2)-qk(:,:,1,2)*qk(:,:,3,3))/wks1
    wkm1(:,:,io2,3) = (qk(:,:,1,2)*qk(:,:,2,3)-qk(:,:,1,3)*qk(:,:,2,2))/wks1
    wkm1(:,:,io1,1) = (qk(:,:,2,3)*qk(:,:,3,1)-qk(:,:,2,1)*qk(:,:,3,3))/wks1
    wkm1(:,:,io1,2) = (qk(:,:,1,1)*qk(:,:,3,3)-qk(:,:,1,3)*qk(:,:,3,1))/wks1
    wkm1(:,:,io1,3) = (qk(:,:,1,3)*qk(:,:,2,1)-qk(:,:,1,1)*qk(:,:,2,3))/wks1
    wkm1(:,:,ihe,1) = (qk(:,:,2,1)*qk(:,:,3,2)-qk(:,:,2,2)*qk(:,:,3,1))/wks1
    wkm1(:,:,ihe,2) = (qk(:,:,1,2)*qk(:,:,3,1)-qk(:,:,1,1)*qk(:,:,3,2))/wks1
    wkm1(:,:,ihe,3) = (qk(:,:,1,1)*qk(:,:,2,2)-qk(:,:,1,2)*qk(:,:,2,1))/wks1

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

  o2_upd = upd(:,:,:,io2)
  o1_upd = upd(:,:,:,io1)
  he_upd = upd(:,:,:,ihe)

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

  o2_upd = max(o2_upd,small)
  o1_upd = max(o1_upd,small)
  he_upd = max(he_upd,small_he)

  o2nm_upd = max(o2nm_upd,small)
  o1nm_upd = max(o1nm_upd,small)
  henm_upd = max(henm_upd,small_he)

  normalize = o2_upd+o1_upd+he_upd
  where (1.-small-normalize < 0.)
    o2_upd = o2_upd*(1.-small)/normalize
    o1_upd = o1_upd*(1.-small)/normalize
    he_upd = he_upd*(1.-small)/normalize
  endwhere

  normalize = o2nm_upd+o1nm_upd+henm_upd
  where (1.-small-normalize < 0.)
    o2nm_upd = o2nm_upd*(1.-small)/normalize
    o1nm_upd = o1nm_upd*(1.-small)/normalize
    henm_upd = henm_upd*(1.-small)/normalize
  endwhere

  n2 = 1.-o2_upd-o1_upd-he_upd

end subroutine comp_ng

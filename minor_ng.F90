subroutine minor_ng(fcomp,fcomp_tm1,fcomp_out,fcomp_tm1_out,sloss,sprod,flbc,fubc,rmx,phix,alfax,name,istep,i_ng)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng,glat_ng
  use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,p0,rmass_o2,rmass_o1,rmass_n2, &
    dtr,pi,grav,avo,dtsmooth,dtsmooth_div2,difhor,rmassinv_he
  use init_module,only: iday
  use fields_ng_module,only: hor,b,fb,flds,itp,itc,shapiro,dtx2inv,dzp,expzmid_inv,expzmid,expz,difk,bndry
  use char_module,only: find_index
  implicit none

  integer,intent(in) :: istep,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: &
    fcomp,fcomp_tm1,sloss,sprod
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    fcomp_out,fcomp_tm1_out
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,3),intent(in) :: flbc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: fubc
  real,intent(in) :: rmx,alfax
  real,dimension(3),intent(in) :: phix
  character(len=*),intent(in) :: name

  real,parameter :: small = 1.e-12, tau = 1.86e+3, t00 = 273.
  real,dimension(2,3),parameter :: phi = reshape((/0.,1.35,1.11,0.673,0.,0.769/),(/2,3/),order=(/2,1/))
  real,parameter :: salfa12 = phi(1,2)-phi(1,3), salfa21 = phi(2,1)-phi(2,3)
  integer :: nk,latd0,latd1,k,lat,idx_minor
  real :: salfax1,salfax2
  logical,dimension(4) :: is_bndry
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: rlat,dfactor
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tlbc,xmbari,bo2,bo1,bhe
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,he,mbar,barm,xnmbar,w,hadvec,do2dz,do1dz,pso2,pso1,dmdz,xmbar_k,tni,s0prod, &
    alfa11,alfa12,alfa21,alfa22,ex,ax,thdiff,p_coef,q_coef,r_coef,f_rhs,ftm1_smooth,wi,dtnidz
  external :: smooth_ng,advec_ng,trsolv_ng

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  barm = flds(i_ng)%barm(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))

  w = flds(i_ng)%w(:,:,:,itc(i_ng))
  tlbc = flds(i_ng)%tlbc

  nk = nlevp1_ng(i_ng)
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  is_bndry = flds(i_ng)%is_bndry

  salfax1 = phix(1)-phix(3)
  salfax2 = phix(2)-phix(3)

  call smooth_ng(fcomp_tm1,ftm1_smooth,shapiro(i_ng),i_ng)

  call advec_ng(fcomp,hadvec,i_ng)

  bo2 = b(1,1)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:)+fb(1)
  bo1 = b(2,1)*o2(1,:,:)+b(2,2)*o1(1,:,:)+b(2,3)*he(1,:,:)+fb(2)
  bhe = b(3,1)*o2(1,:,:)+b(3,2)*o1(1,:,:)+b(3,3)*he(1,:,:)+fb(3)
  xmbari = 1./(bo2*rmassinv_o2+bo1*rmassinv_o1+(1.-bo2-bo1-bhe)*rmassinv_n2+bhe*rmassinv_he)

  xmbar_k = barm
  xmbar_k(1,:,:) = .5*(xmbari+mbar(1,:,:))

  dmdz(1,:,:) = (mbar(1,:,:)-xmbari)/dzp(i_ng)
  pso1(1,:,:) = .5*(bo1+o1(1,:,:))
  pso2(1,:,:) = .5*(bo2+o2(1,:,:))
  do1dz(1,:,:) = (o1(1,:,:)-bo1)/dzp(i_ng)
  do2dz(1,:,:) = (o2(1,:,:)-bo2)/dzp(i_ng)

  do k = 2,nk
    dmdz(k,:,:) = (mbar(k,:,:)-mbar(k-1,:,:))/dzp(i_ng)
    pso1(k,:,:) = .5*(o1(k,:,:)+o1(k-1,:,:))
    pso2(k,:,:) = .5*(o2(k,:,:)+o2(k-1,:,:))
    do1dz(k,:,:) = (o1(k,:,:)-o1(k-1,:,:))/dzp(i_ng)
    do2dz(k,:,:) = (o2(k,:,:)-o2(k-1,:,:))/dzp(i_ng)
  enddo

  tni(1,:,:) = tlbc
  tni(nk,:,:) = tn(nk-1,:,:)
  do k = 2,nk-1
    tni(k,:,:) = .5*(tn(k,:,:)+tn(k-1,:,:))
  enddo

  s0prod = sprod*rmx/xnmbar

  alfa11 = -(phi(1,3)+salfa12*pso1)
  alfa12 = salfa12*pso2
  alfa21 = salfa21*pso1
  alfa22 = -(phi(2,3)+salfa21*pso2)

  ex = ((salfax1*alfa22-salfax2*alfa21)*(do2dz-(1.-(rmass_o2+dmdz)/xmbar_k)*pso2)+ &
    (salfax2*alfa11-salfax1*alfa12)*(do1dz-(1.-(rmass_o1+dmdz)/xmbar_k)*pso1))/ &
    (alfa11*alfa22-alfa12*alfa21)+1.-(rmx+dmdz)/xmbar_k

  dmdz = dmdz/xmbar_k

  ax = -xmbar_k/(tau*rmass_n2)*(t00/tni)**0.25/(phix(3)+salfax1*pso2+salfax2*pso1)

  do k = 2,nk-1
    dtnidz(k,:,:) = (tni(k+1,:,:)-tni(k-1,:,:))/2.
  enddo
  dtnidz(1,:,:) = tni(2,:,:)-tni(1,:,:)
  dtnidz(nk,:,:) = tni(nk,:,:)-tni(nk-1,:,:)
  thdiff = ex-alfax*dtnidz/(dzp(i_ng)*tni)

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

  do k = 1,nk-1
    wi(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
  enddo
  wi(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)

  do k = 1,nk-1
    do lat = latd0,latd1
      p_coef(k,:,lat) = ax(k,:,lat)/dzp(i_ng)*(1./dzp(i_ng)+.5*thdiff(k,:,lat))- &
        expz(i_ng,k)*(expzmid_inv(i_ng)*difk(i_ng,k,iday)*dfactor(lat)*(1./dzp(i_ng)-.5*dmdz(k,:,lat))+ &
        .5*wi(k,:,lat))/dzp(i_ng)
      r_coef(k,:,lat) = ax(k+1,:,lat)/dzp(i_ng)*(1./dzp(i_ng)-.5*thdiff(k+1,:,lat))- &
        expz(i_ng,k)*(expzmid(i_ng)*difk(i_ng,k+1,iday)*dfactor(lat)*(1./dzp(i_ng)+.5*dmdz(k+1,:,lat))- &
        .5*wi(k,:,lat))/dzp(i_ng)
      q_coef(k,:,lat) = -(ax(k,:,lat)/dzp(i_ng)*(1./dzp(i_ng)-.5*thdiff(k,:,lat))+ &
        ax(k+1,:,lat)/dzp(i_ng)*(1./dzp(i_ng)+.5*thdiff(k+1,:,lat)))+ &
        expz(i_ng,k)*((expzmid_inv(i_ng)*difk(i_ng,k,iday)*(1./dzp(i_ng)+.5*dmdz(k,:,lat))+ &
        expzmid(i_ng)*difk(i_ng,k+1,iday)*(1./dzp(i_ng)-.5*dmdz(k+1,:,lat)))* &
        dfactor(lat)/dzp(i_ng)-sloss(k,:,lat)+dtx2inv(i_ng))
    enddo
  enddo

  f_rhs = ftm1_smooth*dtx2inv(i_ng)-hadvec+s0prod
  do k = 1,nk
    f_rhs(k,:,:) = f_rhs(k,:,:)*expz(i_ng,k)
  enddo

  q_coef(1,:,:) = q_coef(1,:,:)+p_coef(1,:,:)*(flbc(:,:,1)+.5*flbc(:,:,2)*dzp(i_ng))/(flbc(:,:,1)-.5*flbc(:,:,2)*dzp(i_ng))
  f_rhs(1,:,:) = f_rhs(1,:,:)-p_coef(1,:,:)*flbc(:,:,3)*dzp(i_ng)/(flbc(:,:,1)-.5*flbc(:,:,2)*dzp(i_ng))
  p_coef(1,:,:) = 0.

  p_coef(nk,:,:) = 1.+.5*dzp(i_ng)*thdiff(nk,:,:)
  q_coef(nk,:,:) = p_coef(nk,:,:)-2.
  r_coef(nk,:,:) = 0.
  f_rhs(nk,:,:) = -grav*rmx*fubc*dzp(i_ng)/(p0*ax(nk,:,:)*avo)

  call trsolv_ng(p_coef,q_coef,r_coef,f_rhs,fcomp_out,nk,i_ng)

  idx_minor = find_index(name,bndry)
  if (is_bndry(1)) fcomp_out(:,-1,:) = flds(i_ng)%lon_b(:,1,:,istep,idx_minor)
  if (is_bndry(2)) fcomp_out(:,nlon_ng(i_ng)+2,:) = flds(i_ng)%lon_b(:,2,:,istep,idx_minor)
  if (is_bndry(3)) fcomp_out(:,:,-1) = flds(i_ng)%lat_b(:,:,3,istep,idx_minor)
  if (is_bndry(4)) fcomp_out(:,:,nlat_ng(i_ng)+2) = flds(i_ng)%lat_b(:,:,4,istep,idx_minor)

  fcomp_tm1_out = dtsmooth*fcomp+dtsmooth_div2*(fcomp_tm1+fcomp_out)

  fcomp_out = max(fcomp_out,small)
  fcomp_tm1_out = max(fcomp_tm1_out,small)

end subroutine minor_ng

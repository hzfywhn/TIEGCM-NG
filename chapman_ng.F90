subroutine chapman_ng(vo2,vo1,vn2,sco2,sco1,scn2,chi,i_ng)

  use params_module,only: nlevp1_ng,glon_ng,glat_ng
  use init_module,only: secs,sin_sundec,cos_sundec
  use cons_module,only: pi,re,rmass_o2,rmass_o1,rmass_n2,dtr
  use fields_ng_module,only: flds,itp,itc
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    vo2,vo1,vn2,sco2,sco1,scn2
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: chi

  integer :: nk,lond0,lond1,latd0,latd1,k,lat
  real :: rlat
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1) :: slt
  integer,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: idn
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tlbc,sin_chi,cos_chi,rt_sinchi
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: z,tn,o2,o1,n2,rp,ti

  z = flds(i_ng)%z(:,:,:,itc(i_ng))
  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2

  nk = nlevp1_ng(i_ng)
  lond0 = flds(i_ng)%lond0
  lond1 = flds(i_ng)%lond1
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  tlbc = flds(i_ng)%tlbc

  slt = modulo(secs/3600.+glon_ng(i_ng,lond0:lond1)/15.,24.)
  do lat = latd0,latd1
    rlat = glat_ng(i_ng,lat)*dtr
    chi(:,lat) = acos(sin_sundec*sin(rlat)+cos_sundec*cos(rlat)*cos(pi*(slt-12.)/12.))
  enddo
  sin_chi = sin(chi)
  cos_chi = cos(chi)
  rt_sinchi = sqrt(sin_chi)

  idn = 1
  where (chi > 1.8326) idn = 0

  rp = z+re

  ti(1,:,:) = tlbc
  do k = 2,nk-1
    ti(k,:,:) = .5*(tn(k-1,:,:)+tn(k,:,:))
  enddo
  ti(nk,:,:) = tn(nk-1,:,:)

  call line_integ(o2,rmass_o2,vo2,sco2)
  call line_integ(o1,rmass_o1,vo1,sco1)
  call line_integ(n2,rmass_n2,vn2,scn2)

  contains
!-----------------------------------------------------------------------
  subroutine line_integ(f,fmass,v,s)

    use cons_module,only: grav,p0,avo,gask
    use fields_ng_module,only: dz,expz

    real,intent(in) :: fmass
    real,dimension(nk,lond0:lond1,latd0:latd1),intent(in) :: f
    real,dimension(nk,lond0:lond1,latd0:latd1),intent(out) :: v,s

    real,parameter :: big = 1.e80, logbig = log(big)
    integer :: k
    real :: rtpi
    real,dimension(lond0:lond1,latd0:latd1) :: exparg,r2ig
    real,dimension(nk,lond0:lond1,latd0:latd1) :: barm,rtrp,yp

    barm = flds(i_ng)%barm(:,:,:,itp(i_ng))

    v = max(f,0.)
    v(nk,:,:) = avo*p0*expz(i_ng,nk-1)*exp(-.5*dz(i_ng))/(fmass**2*grav)* &
      .5*(v(nk-1,:,:)+v(nk,:,:))*barm(nk,:,:)

    do k = nk-1,1,-1
      v(k,:,:) = v(k+1,:,:)+avo*p0*expz(i_ng,k)/(fmass*grav)*dz(i_ng)*v(k,:,:)
    enddo

    rtpi = sqrt(pi)
    rtrp = sqrt(rp*fmass*grav/(2.*gask)/ti)
    do k = 1,nk
      yp(k,:,:) = rtrp(k,:,:)*abs(cos_chi)
    enddo
    where (yp >= 8.)
      yp = v*rtpi*rtrp*0.56498823/(0.06651874+yp)
    elsewhere
      yp = v*rtpi*rtrp*(1.0606963+0.5564383*yp)/((yp+1.7245609)*yp+1.0619896)
    endwhere

    do k = 1,nk
      exparg = rp(k,:,:)*(1.-sin_chi)*grav*fmass/gask/ti(k,:,:)
      where (idn==1 .and. exparg<logbig)
        r2ig = 2.*v(k,:,:)*exp(exparg)*rtpi*rt_sinchi*rtrp(k,:,:)
      elsewhere
        r2ig = big
      endwhere

      where (cos_chi >= 0.)
        s(k,:,:) = big
      elsewhere
        s(k,:,:) = r2ig-yp(k,:,:)
      endwhere
      where (rp(k,:,:)*sin_chi < re) s(k,:,:) = big
      where (cos_chi >= 0.) s(k,:,:) = yp(k,:,:)
    enddo

  end subroutine line_integ
!-----------------------------------------------------------------------
end subroutine chapman_ng

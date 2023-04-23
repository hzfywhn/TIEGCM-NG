subroutine hdif12_ng(fkldt,fkldu,fkldv,fkldo2,fkldo1,fkldhe,fnrh,fkmh,i_ng)
! changed the domains of fnrh,fkmh to lat0:lat1 (previously lat0-1:lat1-1)

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use cons_module,only: re_inv
  use fields_ng_module,only: flds,itp,dlamda,dphi,t0
  use gather2root_ng_module,only: gather2root_var3d
  use mpi_f08,only: mpi_bcast,mpi_real8,mpi_comm_world
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    fkldt,fkldu,fkldv,fkldo2,fkldo1,fkldhe,fnrh,fkmh

  real,parameter :: cp2 = 0.2, con3 = 2.*cp2**2
  integer :: nk,lond0,lond1,latd0,latd1,k,i,lat,lonbeg,latend
  real :: con1,con2,n
  real,dimension(nlevp1_ng(i_ng)) :: dup,du,dvp,dv,ump,um,vmp,vm,delt,dels
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cs
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn_nm,un_nm,vn_nm,o2_nm,o1_nm,he_nm,mbar,avkmh,rhokmh
  real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2) :: full
  external :: lsqdsq_ng

  tn_nm = flds(i_ng)%tn_nm(:,:,:,itp(i_ng))
  un_nm = flds(i_ng)%un_nm(:,:,:,itp(i_ng))
  vn_nm = flds(i_ng)%vn_nm(:,:,:,itp(i_ng))
  o2_nm = flds(i_ng)%o2_nm(:,:,:,itp(i_ng))
  o1_nm = flds(i_ng)%o1_nm(:,:,:,itp(i_ng))
  he_nm = flds(i_ng)%he_nm(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))

  nk = nlevp1_ng(i_ng)
  lond0 = flds(i_ng)%lond0
  lond1 = flds(i_ng)%lond1
  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  cs = flds(i_ng)%cs

  con1 = re_inv*.5/dlamda(i_ng)
  do lat = latd0,latd1-1
    con2 = re_inv/(dphi(i_ng)*(cs(lat)+cs(lat+1)))
    do i = lond0,lond1-1
      dup = (un_nm(:,i+1,lat+1)-un_nm(:,i,lat+1))/cs(lat+1)
      du = (un_nm(:,i+1,lat)-un_nm(:,i,lat))/cs(lat)
      vmp = (vn_nm(:,i+1,lat+1)+vn_nm(:,i,lat+1))*cs(lat+1)
      vm = (vn_nm(:,i+1,lat)+vn_nm(:,i,lat))*cs(lat)
      dvp = (vn_nm(:,i+1,lat+1)-vn_nm(:,i,lat+1))/cs(lat+1)
      dv = (vn_nm(:,i+1,lat)-vn_nm(:,i,lat))/cs(lat)
      ump = (un_nm(:,i+1,lat+1)+un_nm(:,i,lat+1))*cs(lat+1)
      um = (un_nm(:,i+1,lat)+un_nm(:,i,lat))*cs(lat)
      delt = con1*(dup+du)-con2*(vmp-vm)
      dels = con1*(dvp+dv)+con2*(ump-um)
      fkmh(:,i,lat) = con3*sqrt(dels**2+delt**2)
    enddo
    fkmh(:,lond1,lat) = 2.*fkmh(:,lond1-1,lat)-fkmh(:,lond1-2,lat)
  enddo
  fkmh(:,:,latd1) = 2.*fkmh(:,:,latd1-1)-fkmh(:,:,latd1-2)

  call gather2root_var3d(fkmh,full,i_ng)
  call mpi_bcast(full,nk*(nlon_ng(i_ng)+4)*(nlat_ng(i_ng)+4),mpi_real8,0,mpi_comm_world)

  fkmh = full(:,lond0:lond1,latd0:latd1)

  fnrh = mbar/(tn_nm+t0)

  do lat = latd0,latd1
    latend = min(lat+1,nlat_ng(i_ng)+2)
    do i = lond0,lond1
      lonbeg = max(i-1,-1)
      n = (i-lonbeg+1)*(latend-lat+1)
      do k = 1,nk
        avkmh(k,i,lat) = sum(full(k,lonbeg:i,lat:latend))/n
      enddo
    enddo
  enddo

  rhokmh = avkmh*fnrh

  call lsqdsq_ng(un_nm,avkmh,i_ng)
  fkldu = avkmh*rhokmh

  call lsqdsq_ng(vn_nm,avkmh,i_ng)
  fkldv = avkmh*rhokmh

  call lsqdsq_ng(tn_nm,avkmh,i_ng)
  fkldt = avkmh*rhokmh

  call lsqdsq_ng(o2_nm,avkmh,i_ng)
  fkldo2 = avkmh*rhokmh

  call lsqdsq_ng(o1_nm,avkmh,i_ng)
  fkldo1 = avkmh*rhokmh

  call lsqdsq_ng(he_nm,avkmh,i_ng)
  fkldhe = avkmh*rhokmh

end subroutine hdif12_ng
!-----------------------------------------------------------------------
subroutine hdif3_ng(hdt,hdu,hdv,hdo2,hdo1,hdhe,i_ng)

  use params_module,only: nlevp1_ng
  use fields_ng_module,only: flds
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    hdt,hdu,hdv,hdo2,hdo1,hdhe

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    cp,kldt,kldu,kldv,kldo2,kldo1,kldhe,fnrh,fnrh_inv,hdout,cpi
  external :: lsqdsq_ng

  cp = flds(i_ng)%cp
  kldt = flds(i_ng)%kldt
  kldu = flds(i_ng)%kldu
  kldv = flds(i_ng)%kldv
  kldo2 = flds(i_ng)%kldo2
  kldo1 = flds(i_ng)%kldo1
  kldhe = flds(i_ng)%kldhe

  fnrh = flds(i_ng)%fnrh

  nk = nlevp1_ng(i_ng)
 
  fnrh_inv = -1./fnrh

  call lsqdsq_ng(kldu,hdout,i_ng)
  hdu = hdout*fnrh_inv

  call lsqdsq_ng(kldv,hdout,i_ng)
  hdv = hdout*fnrh_inv

  call lsqdsq_ng(kldt,hdout,i_ng)
  do k = 1,nk-1
    cpi(k,:,:) = .5*(cp(k,:,:)+cp(k+1,:,:))
  enddo
  hdt = hdout*fnrh_inv*cpi

  call lsqdsq_ng(kldo2,hdout,i_ng)
  hdo2 = hdout*fnrh_inv

  call lsqdsq_ng(kldo1,hdout,i_ng)
  hdo1 = hdout*fnrh_inv

  call lsqdsq_ng(kldhe,hdout,i_ng)
  hdhe = hdout*fnrh_inv

end subroutine hdif3_ng
!-----------------------------------------------------------------------
subroutine lsqdsq_ng(fj,fout,i_ng)

  use params_module,only: nlevp1_ng
  use fields_ng_module,only: flds
  use dffm_ng_module,only: fm_3d
  use sync_ng_module,only: sync_var3d_lat
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: fj
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: fout

  integer :: latd0,latd1,lat
  real,dimension(flds(i_ng)%latd0:flds(i_ng)%latd1) :: cs
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    dpsi,dpsip,dpsim,fjm1,dlambda

  cs = flds(i_ng)%cs

  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1

  call fm_3d(fj,fjm1,1,i_ng)
  dlambda = fj*(-2.)+fjm1*2.

  do lat = latd0,latd1-1
    dpsip(:,:,lat) = (fj(:,:,lat+1)-fj(:,:,lat))*0.5*(cs(lat+1)+cs(lat))
  enddo
  dpsip(:,:,latd1) = (fj(:,:,latd1)-fj(:,:,latd1-1))*cs(latd1)
  call sync_var3d_lat(dpsip,i_ng)

  dpsim(:,:,latd0) = (fj(:,:,latd0+1)-fj(:,:,latd0))*cs(latd0)
  do lat = latd0+1,latd1
    dpsim(:,:,lat) = (fj(:,:,lat)-fj(:,:,lat-1))*0.5*(cs(lat)+cs(lat-1))
  enddo
  call sync_var3d_lat(dpsim,i_ng)

  dpsi = dpsip-dpsim
  do lat = latd0,latd1
    dpsi(:,:,lat) = dpsi(:,:,lat)/cs(lat)
  enddo

  fout = dlambda+dpsi

end subroutine lsqdsq_ng

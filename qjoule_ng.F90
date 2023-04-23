subroutine qjoule_ti_ng(qji_ti,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: grav,gask,rmass_o2,rmass_o1,rmass_n2,rmass_n4s,rmass_no
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: qji_ti

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    un,vn,w,ui,vi,wi,lam1,op,nplus,n2p,nop,o2p,mbar,tn,vert_vel,uii,vii,wii,m_i,mfac,wni,lam1i

  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))
  w = flds(i_ng)%w(:,:,:,itp(i_ng))
  ui = flds(i_ng)%ui
  vi = flds(i_ng)%vi
  wi = flds(i_ng)%wi
  lam1 = flds(i_ng)%lam1
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  nplus = flds(i_ng)%nplus
  n2p = flds(i_ng)%n2p
  nop = flds(i_ng)%nop
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    wni(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
    uii(k,:,:) = .5*(ui(k,:,:)+ui(k+1,:,:))
    vii(k,:,:) = .5*(vi(k,:,:)+vi(k+1,:,:))
    wii(k,:,:) = .5*(wi(k,:,:)+wi(k+1,:,:))
    lam1i(k,:,:) = .5*(lam1(k,:,:)+lam1(k+1,:,:))
  enddo
  wni(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)
  uii(nk,:,:) = 1.5*ui(nk,:,:)-.5*ui(nk-1,:,:)
  vii(nk,:,:) = 1.5*vi(nk,:,:)-.5*vi(nk-1,:,:)
  wii(nk,:,:) = 1.5*wi(nk,:,:)-.5*wi(nk-1,:,:)
  lam1i(nk,:,:) = 1.5*lam1(nk,:,:)-.5*lam1(nk-1,:,:)

  m_i = (op*rmass_o1+o2p*rmass_o2+nplus*rmass_n4s+n2p*rmass_n2+nop*rmass_no)/(op+o2p+nplus+n2p+nop)
  mfac = mbar/(m_i+mbar)
  vert_vel = wni*gask*tn/(mbar*grav)
  qji_ti = mfac*lam1i*((uii-un)**2+(vii-vn)**2+(wii-vert_vel)**2)

end subroutine qjoule_ti_ng
!-----------------------------------------------------------------------
subroutine qjoule_tn_ng(qji_tn,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: grav,gask
  use input_module,only: joulefac
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: qji_tn

  integer :: nk,k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    un,vn,w,ui,vi,wi,lam1,mbar,tn,vel_zonal,vel_merid,vel_vert,scheight,wni,uii,vii,wii,lam1i

  un = flds(i_ng)%un(:,:,:,itp(i_ng))
  vn = flds(i_ng)%vn(:,:,:,itp(i_ng))
  w = flds(i_ng)%w(:,:,:,itp(i_ng))
  ui = flds(i_ng)%ui
  vi = flds(i_ng)%vi
  wi = flds(i_ng)%wi
  lam1 = flds(i_ng)%lam1
  mbar = flds(i_ng)%mbar(:,:,:,itp(i_ng))
  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))

  nk = nlevp1_ng(i_ng)

  do k = 1,nk-1
    wni(k,:,:) = .5*(w(k,:,:)+w(k+1,:,:))
    uii(k,:,:) = .5*(ui(k,:,:)+ui(k+1,:,:))
    vii(k,:,:) = .5*(vi(k,:,:)+vi(k+1,:,:))
    wii(k,:,:) = .5*(wi(k,:,:)+wi(k+1,:,:))
    lam1i(k,:,:) = .5*(lam1(k,:,:)+lam1(k+1,:,:))
  enddo
  wni(nk,:,:) = 1.5*w(nk,:,:)-.5*w(nk-1,:,:)
  uii(nk,:,:) = 1.5*ui(nk,:,:)-.5*ui(nk-1,:,:)
  vii(nk,:,:) = 1.5*vi(nk,:,:)-.5*vi(nk-1,:,:)
  wii(nk,:,:) = 1.5*wi(nk,:,:)-.5*wi(nk-1,:,:)
  lam1i(nk,:,:) = 1.5*lam1(nk,:,:)-.5*lam1(nk-1,:,:)

  scheight = gask*tn/(mbar*grav)
  vel_zonal = uii-un
  vel_merid = vii-vn
  vel_vert = wii-scheight*wni
  qji_tn = lam1i*(vel_zonal**2+vel_merid**2+vel_vert**2)*joulefac

end subroutine qjoule_tn_ng

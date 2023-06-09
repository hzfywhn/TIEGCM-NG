subroutine lamdas_ng(lxx,lyy,lxy,lyx,lamda1,ped_out,hall_out,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: rmassinv_o2,rmassinv_o1,rmassinv_n2,avo
  use input_module,only: colfac
  use fields_ng_module,only: flds,itp,dipmin
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    lxx,lyy,lxy,lyx,lamda1,ped_out,hall_out

  real,parameter :: qe = 1.602e-19, qeomeo10 = 1.7588028E7, qeoNao10 = 9.6489E3, rmass_nop = 30., &
    rmassinv_nop = 1./rmass_nop, rnu_op_o2 = 6.64E-10, rnu_nop_o2 = 4.27E-10, rnu_o2p_o = 2.31E-10, &
    rnu_nop_o = 2.44E-10, rnu_o2p_n2 = 4.13E-10, rnu_op_n2 = 6.82E-10, rnu_nop_n2 = 4.34E-10
  integer :: nk,k
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    bmod2,sn2dec,csdec,sndec,dipmag,qe_fac,sindip,cosdip,cos2dip,sin2dip,cos2dec, &
    omega_o2p,omega_op,omega_nop,omega_o2p_inv,omega_op_inv,omega_nop_inv,omega_e,omega_e_inv,dip
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,xnmbar,o2,o1,n2,ti,te,o2p,op,nop,tnti,o2_cm3,o1_cm3,n2_cm3,sigma_ped,sigma_hall, &
    ne,lamda2,lamda1tmp,lamda2tmp,lxxnorot,lyynorot,lxynorot,lyxnorot, &
    rnu_o2p_o2,rnu_op_o,rnu_o2p,rnu_op,rnu_nop,rnu_ne,sqrt_te

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  nop = flds(i_ng)%nop

  bmod2 = flds(i_ng)%bmod2
  sn2dec = flds(i_ng)%sn2dec
  csdec = flds(i_ng)%csdec
  sndec = flds(i_ng)%sndec
  dipmag = flds(i_ng)%dipmag

  nk = nlevp1_ng(i_ng)

  qe_fac = qe*1.e10/bmod2
  omega_op = qeoNao10*bmod2*rmassinv_o1
  omega_o2p = qeoNao10*bmod2*rmassinv_o2
  omega_nop = qeoNao10*bmod2*rmassinv_nop
  omega_op_inv = 1./omega_op
  omega_o2p_inv = 1./omega_o2p
  omega_nop_inv = 1./omega_nop
  omega_e = qeomeo10*bmod2
  omega_e_inv = 1./omega_e

  dip = sign(max(abs(dipmag),dipmin(i_ng)),dipmag)
  sindip = sin(dip)
  cosdip = cos(dip)
  cos2dip = cosdip**2
  sin2dip = sindip**2
  cos2dec = csdec**2

  tnti = 0.5*(ti+tn)
  rnu_o2p_o2 = 2.59E-11*sqrt(tnti)*(1.-0.073*log10(tnti))**2
  rnu_op_o = 3.67e-11*sqrt(tnti)*(1.-0.064*log10(tnti))**2*colfac

  o2_cm3 = o2*xnmbar*rmassinv_o2
  o1_cm3 = o1*xnmbar*rmassinv_o1
  n2_cm3 = n2*rmassinv_n2*xnmbar

  rnu_o2p = rnu_o2p_o2*o2_cm3+rnu_o2p_o*o1_cm3+rnu_o2p_n2*n2_cm3
  rnu_op = rnu_op_o2*o2_cm3+rnu_op_o*o1_cm3+rnu_op_n2*n2_cm3
  rnu_nop = rnu_nop_o2*o2_cm3+rnu_nop_o*o1_cm3+rnu_nop_n2*n2_cm3
  sqrt_te = sqrt(te)
  rnu_ne = 2.33e-11*n2_cm3*te*(1.-1.21e-4*te)+ &
    1.82e-10*o2_cm3*sqrt_te*(1.+3.60e-2*sqrt_te)+ &
    8.90e-11*o1_cm3*sqrt_te*(1.+5.70e-4*te)
  do k = 1,nk
    rnu_o2p(k,:,:) = rnu_o2p(k,:,:)*omega_o2p_inv
    rnu_op(k,:,:) = rnu_op(k,:,:)*omega_op_inv
    rnu_nop(k,:,:) = rnu_nop(k,:,:)*omega_nop_inv
    rnu_ne(k,:,:) = rnu_ne(k,:,:)*omega_e_inv
  enddo
  rnu_ne = rnu_ne*4.

  ne = op+o2p+nop

  sigma_ped = op*rnu_op/(1.+rnu_op**2)+o2p*rnu_o2p/(1.+rnu_o2p**2)+nop*rnu_nop/(1.+rnu_nop**2)+ne*rnu_ne/(1.+rnu_ne**2)
  sigma_hall = ne/(1.+rnu_ne**2)-op/(1.+rnu_op**2)-o2p/(1.+rnu_o2p**2)-nop/(1.+rnu_nop**2)
  do k = 1,nk
    sigma_ped(k,:,:) = sigma_ped(k,:,:)*qe_fac
    sigma_hall(k,:,:) = sigma_hall(k,:,:)*qe_fac
  enddo

  lamda1tmp = sigma_ped*avo/(1.e3*xnmbar)
  lamda2tmp = sigma_hall*avo/(1.e3*xnmbar)
  do k = 1,nk
    lamda1tmp(k,:,:) = lamda1tmp(k,:,:)*(1.e-4*bmod2)**2
    lamda2tmp(k,:,:) = lamda2tmp(k,:,:)*(1.e-4*bmod2)**2
  enddo

  do k = 1,nk-2
    lamda1(k+1,:,:) = sqrt(lamda1tmp(k,:,:)*lamda1tmp(k+1,:,:))
    lamda2(k+1,:,:) = sqrt(lamda2tmp(k,:,:)*lamda2tmp(k+1,:,:))
  enddo
  lamda1(1,:,:) = sqrt(lamda1tmp(1,:,:)**3/lamda1tmp(2,:,:))
  lamda2(1,:,:) = sqrt(lamda2tmp(1,:,:)**3/lamda2tmp(2,:,:))
  lamda1(nk,:,:) = sqrt(lamda1tmp(nk-1,:,:)**3/lamda1tmp(nk-2,:,:))
  lamda2(nk,:,:) = sqrt(lamda2tmp(nk-1,:,:)**3/lamda2tmp(nk-2,:,:))

  lxxnorot = lamda1
  do k = 1,nk
    lyynorot(k,:,:) = lamda1(k,:,:)*sin2dip
    lxynorot(k,:,:) = lamda2(k,:,:)*sindip
  enddo
  lyxnorot = lxynorot

  do k = 1,nk
    lxx(k,:,:) = lxxnorot(k,:,:)*cos2dec+lyynorot(k,:,:)*sn2dec
    lyy(k,:,:) = lyynorot(k,:,:)*cos2dec+lxxnorot(k,:,:)*sn2dec
    lyx(k,:,:) = lxynorot(k,:,:)-(lyynorot(k,:,:)-lxxnorot(k,:,:))*sndec*csdec
    lxy(k,:,:) = lxynorot(k,:,:)+(lyynorot(k,:,:)-lxxnorot(k,:,:))*sndec*csdec
  enddo

  ped_out = sigma_ped
  hall_out = sigma_hall

end subroutine lamdas_ng

subroutine lamdas_ng(lxx,lyy,lxy,lyx,lamda1,ped_out,hall_out,Q1,Q2,i_ng)

  use params_module,only: nlevp1_ng
  use cons_module,only: avo,rtd, &
    rmass_o2,rmass_o1,rmass_he,rmass_n2,rmass_no,rmass_n4s, &
    rmassinv_o2,rmassinv_o1,rmassinv_he,rmassinv_n2,rmassinv_no,rmassinv_n4s
  use input_module,only: colfac
  use fields_ng_module,only: flds,itp,dipmin
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    lxx,lyy,lxy,lyx,lamda1,ped_out,hall_out,Q1,Q2

  real,parameter :: qe = 1.602e-19, qeomeo10 = 1.7588028E7, qeoNao10 = 9.6489E3, &
    rnu_op_o2 = 6.64E-10, rnu_np_o2 = 7.25E-10, rnu_n2p_o2 = 4.49E-10, rnu_nop_o2 = 4.27E-10, &
    rnu_o2p_o = 2.31E-10, rnu_np_o = 4.42E-10, rnu_n2p_o = 2.58E-10, rnu_nop_o = 2.44E-10, &
    rnu_o2p_he = 0.70E-10, rnu_op_he = 1.32E-10, rnu_np_he = 1.49E-10, rnu_n2p_he = 0.79E-10, rnu_nop_he = 0.74E-10, &
    rnu_o2p_n2 = 4.13E-10, rnu_op_n2 = 6.82E-10, rnu_np_n2 = 7.47E-10, rnu_nop_n2 = 4.34E-10, &
    Me = 9.109E-31, Mp = 1.6726E-27, Kb = 1.38E-23
  integer :: nk,k
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    bmod2,sn2dec,csdec,sndec,dipmag,qe_fac,sindip,cosdip,cos2dip,sin2dip,cos2dec,dip,rlatm, &
    omega_o2p,omega_op,omega_np,omega_n2p,omega_nop,omega_e, &
    omega_o2p_inv,omega_op_inv,omega_np_inv,omega_n2p_inv,omega_nop_inv,omega_e_inv
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,xnmbar,o2,o1,he,n2,ti,te,o2p,op,nplus,n2p,nop,tnti,o2_cm3,o1_cm3,he_cm3,n2_cm3,sigma_ped,sigma_hall, &
    ne,lamda2,lamda1tmp,lamda2tmp,lxxnorot,lyynorot,lxynorot,lyxnorot, &
    rnu_o2p_o2,rnu_op_o,rnu_n2p_n2,rnu_o2p,rnu_op,rnu_np,rnu_n2p,rnu_nop,rnu_ne,sqrt_te,Etot,Ki,Mi,E1

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  ti = flds(i_ng)%ti(:,:,:,itp(i_ng))
  te = flds(i_ng)%te(:,:,:,itp(i_ng))
  o2p = flds(i_ng)%o2p(:,:,:,itp(i_ng))
  op = flds(i_ng)%op(:,:,:,itp(i_ng))
  nplus = flds(i_ng)%nplus
  n2p = flds(i_ng)%n2p
  nop = flds(i_ng)%nop
  Etot = flds(i_ng)%Etot

  bmod2 = flds(i_ng)%bmod2
  sn2dec = flds(i_ng)%sn2dec
  csdec = flds(i_ng)%csdec
  sndec = flds(i_ng)%sndec
  dipmag = flds(i_ng)%dipmag
  rlatm = flds(i_ng)%rlatm

  nk = nlevp1_ng(i_ng)

  qe_fac = qe*1.e10/bmod2
  omega_op = qeoNao10*bmod2*rmassinv_o1
  omega_o2p = qeoNao10*bmod2*rmassinv_o2
  omega_np = qeoNao10*bmod2*rmassinv_n4s
  omega_n2p = qeoNao10*bmod2*rmassinv_n2
  omega_nop = qeoNao10*bmod2*rmassinv_no
  omega_op_inv = 1./omega_op
  omega_o2p_inv = 1./omega_o2p
  omega_np_inv = 1./omega_np
  omega_n2p_inv = 1./omega_n2p
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
  rnu_n2p_n2 = 5.14E-11*sqrt(tnti)*(1.-0.069*log10(tnti))**2

  o2_cm3 = o2*xnmbar*rmassinv_o2
  o1_cm3 = o1*xnmbar*rmassinv_o1
  he_cm3 = he*xnmbar*rmassinv_he
  n2_cm3 = n2*xnmbar*rmassinv_n2

  rnu_o2p = rnu_o2p_o2*o2_cm3+rnu_o2p_o*o1_cm3+rnu_o2p_he*he_cm3+rnu_o2p_n2*n2_cm3
  rnu_op = rnu_op_o2*o2_cm3+rnu_op_o*o1_cm3+rnu_op_he*he_cm3+rnu_op_n2*n2_cm3
  rnu_np = rnu_np_o2*o2_cm3+rnu_np_o*o1_cm3+rnu_np_he*he_cm3+rnu_np_n2*n2_cm3
  rnu_n2p = rnu_n2p_o2*o2_cm3+rnu_n2p_o*o1_cm3+rnu_n2p_he*he_cm3+rnu_n2p_n2*n2_cm3
  rnu_nop = rnu_nop_o2*o2_cm3+rnu_nop_o*o1_cm3+rnu_nop_he*he_cm3+rnu_nop_n2*n2_cm3
  sqrt_te = sqrt(te)
  rnu_ne = 1.82e-10*o2_cm3*sqrt_te*(1.+3.60e-2*sqrt_te)+ &
    8.90e-11*o1_cm3*sqrt_te*(1.+5.70e-4*te)+ &
    4.60e-10*he_cm3*sqrt_te+ &
    2.33e-11*n2_cm3*te*(1.-1.21e-4*te)
  do k = 1,nk
    rnu_o2p(k,:,:) = rnu_o2p(k,:,:)*omega_o2p_inv
    rnu_op(k,:,:) = rnu_op(k,:,:)*omega_op_inv
    rnu_np(k,:,:) = rnu_np(k,:,:)*omega_np_inv
    rnu_n2p(k,:,:) = rnu_n2p(k,:,:)*omega_n2p_inv
    rnu_nop(k,:,:) = rnu_nop(k,:,:)*omega_nop_inv
    rnu_ne(k,:,:) = rnu_ne(k,:,:)*omega_e_inv
  enddo
  rnu_ne = rnu_ne*4.

  ne = op+o2p+nplus+n2p+nop

  sigma_ped = op*rnu_op/(1.+rnu_op**2)+ &
    o2p*rnu_o2p/(1.+rnu_o2p**2)+ &
    nplus*rnu_np/(1.+rnu_np**2)+ &
    n2p*rnu_n2p/(1.+rnu_n2p**2)+ &
    nop*rnu_nop/(1.+rnu_nop**2)+ &
    ne*rnu_ne/(1.+rnu_ne**2)
  sigma_hall = ne/(1.+rnu_ne**2)- &
    op/(1.+rnu_op**2)- &
    o2p/(1.+rnu_o2p**2)- &
    nplus/(1.+rnu_np**2)- &
    n2p/(1.+rnu_n2p**2)- &
    nop/(1.+rnu_nop**2)
  do k = 1,nk
    sigma_ped(k,:,:) = sigma_ped(k,:,:)*qe_fac
    sigma_hall(k,:,:) = sigma_hall(k,:,:)*qe_fac
  enddo
  sigma_hall = max(sigma_hall,1e-20)

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

  call addfld(lamda1,'LAMDA_PED',i_ng)
  call addfld(lamda2,'LAMDA_HAL',i_ng)

  do k = 1,nk
    lxx(k,:,:) = lxxnorot(k,:,:)*cos2dec+lyynorot(k,:,:)*sn2dec
    lyy(k,:,:) = lyynorot(k,:,:)*cos2dec+lxxnorot(k,:,:)*sn2dec
    lyx(k,:,:) = lxynorot(k,:,:)-(lyynorot(k,:,:)-lxxnorot(k,:,:))*sndec*csdec
    lxy(k,:,:) = lxynorot(k,:,:)+(lyynorot(k,:,:)-lxxnorot(k,:,:))*sndec*csdec
  enddo

  ped_out = sigma_ped
  hall_out = sigma_hall

  call addfld(ped_out,'SIGMA_PED',i_ng)
  call addfld(hall_out,'SIGMA_HAL',i_ng)

  Ki = 1/ne*(op/rnu_op+o2p/rnu_o2p+nplus/rnu_np+n2p/rnu_n2p+nop/rnu_nop)
  Mi = 1/ne*(op*rmass_o1+o2p*rmass_o2+nplus*rmass_n4s+n2p*rmass_n2+nop*rmass_no)

  E1 = (1.0+rnu_ne/Ki)*sqrt(Kb*(1.0+Ki**2)/(1.0-Ki**2)*(te+ti)/(Mi*Mp))
  do k = 1,nk
    E1(k,:,:) = E1(k,:,:)*bmod2*1E-4
  enddo

  Q1 = 0.0
  Q2 = 0.0
  do k = 1,nk
    where (abs(rlatm)*rtd>50.0 .and. Ki(k,:,:)<1.0 .and. Etot(k,:,:)>E1(k,:,:))
      Q1(k,:,:) = Me*ne(k,:,:)*1E6*Etot(k,:,:)**2/ &
        (bmod2*1E-4)**2*rnu_ne(k,:,:)/omega_e_inv
      Q2(k,:,:) = Mi(k,:,:)*Mp*1E6* &
        Ki(k,:,:)**2*(Etot(k,:,:)-E1(k,:,:))**2/(1.0+Ki(k,:,:)**2)* &
        (Etot(k,:,:)/E1(k,:,:)*(1.0+rnu_ne(k,:,:)/Ki(k,:,:))-1.0)/ &
        (bmod2*1E-4)**2* &
        (op(k,:,:)*rnu_op(k,:,:)/omega_op_inv+ &
        o2p(k,:,:)*rnu_o2p(k,:,:)/omega_o2p_inv+ &
        nplus(k,:,:)*rnu_np(k,:,:)/omega_np_inv+ &
        n2p(k,:,:)*rnu_n2p(k,:,:)/omega_n2p_inv+ &
        nop(k,:,:)*rnu_nop(k,:,:)/omega_nop_inv)
    endwhere
  enddo

end subroutine lamdas_ng

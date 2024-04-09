subroutine dynamics_ng(istep,i_ng)
! input parameters are passed through module use, output parameters are passed by subroutine calls

  use params_module,only: nlevp1_ng
  use input_module,only: nstep_sub,aurora
  use fields_ng_module,only: flds,itp,itc
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: istep,i_ng

  integer :: itp_sub,itc_sub,istep_sub,itmp
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,2) :: op_i,op_nm_i
  external :: cpktkm_ng,swdot_ng,ionvel_ng,chapman_ng,chemrates_tdep_ng,qrj_ng,qinite_ng,aurora_ng, &
    oplus_ng,elden_ng,lamdas_ng,comp_n2d_ng,comp_n4s_ng,comp_no_ng,comp_ar_ng,qjnno_ng,qjion_ng, &
    qjoule_ti_ng,settei_ng,newton_ng,hdif3_ng,qjoule_tn_ng,dt_ng,duv_ng,comp_o2o_ng,comp_ng

  call cpktkm_ng( &
    flds(i_ng)%cp, &
    flds(i_ng)%kt, &
    flds(i_ng)%km, &
    i_ng)

  call swdot_ng( &
    flds(i_ng)%w(:,:,:,itc(i_ng)), &
    istep,i_ng)

  call ionvel_ng( &
    flds(i_ng)%ui, &
    flds(i_ng)%vi, &
    flds(i_ng)%wi, &
    flds(i_ng)%Etot, &
    i_ng)

  call chapman_ng( &
    flds(i_ng)%vo2, &
    flds(i_ng)%vo1, &
    flds(i_ng)%vn2, &
    flds(i_ng)%sco2, &
    flds(i_ng)%sco1, &
    flds(i_ng)%scn2, &
    flds(i_ng)%chi, &
    i_ng)

  call chemrates_tdep_ng( &
    flds(i_ng)%rk1, &
    flds(i_ng)%rk2, &
    flds(i_ng)%rk3, &
    flds(i_ng)%ra1, &
    flds(i_ng)%ra2, &
    flds(i_ng)%ra3, &
    flds(i_ng)%beta1, &
    flds(i_ng)%beta3, &
    flds(i_ng)%beta5, &
    flds(i_ng)%beta8, &
    flds(i_ng)%beta9, &
    flds(i_ng)%beta9n, &
    flds(i_ng)%beta17, &
    flds(i_ng)%rk19, &
    flds(i_ng)%rk20, &
    flds(i_ng)%rk25, &
    flds(i_ng)%rkm12, &
    i_ng)

  call qrj_ng( &
    flds(i_ng)%rj, &
    flds(i_ng)%qtef, &
    flds(i_ng)%qtotal, &
    flds(i_ng)%qop2p, &
    flds(i_ng)%qop2d, &
    flds(i_ng)%qo2p, &
    flds(i_ng)%qop, &
    flds(i_ng)%qn2p, &
    flds(i_ng)%qnp, &
    flds(i_ng)%qnop, &
    i_ng)

  call qinite_ng( &
    flds(i_ng)%qo2p, &
    flds(i_ng)%qop, &
    flds(i_ng)%qn2p, &
    flds(i_ng)%qnp, &
    flds(i_ng)%qnop, &
    flds(i_ng)%qtef, &
    i_ng)

  if (aurora > 0) &
    call aurora_ng( &
    flds(i_ng)%qteaur, &
    flds(i_ng)%qo2p, &
    flds(i_ng)%qop, &
    flds(i_ng)%qn2p, &
    flds(i_ng)%qnp, &
    flds(i_ng)%qtef, &
    flds(i_ng)%ui, &
    flds(i_ng)%vi, &
    i_ng)

  itp_sub = itp(i_ng)
  itc_sub = itc(i_ng)
  op_i = flds(i_ng)%op
  op_nm_i = flds(i_ng)%op_nm
  do istep_sub = 1,nstep_sub
    call oplus_ng( &
      op_i(:,:,:,itp_sub), &
      op_nm_i(:,:,:,itp_sub), &
      op_i(:,:,:,itc_sub), &
      op_nm_i(:,:,:,itc_sub), &
      flds(i_ng)%xiop2p, &
      flds(i_ng)%xiop2d, &
      flds(i_ng)%Fe, &
      flds(i_ng)%Fn, &
      (istep-1)*nstep_sub+istep_sub,i_ng)
    if (istep_sub /= nstep_sub) then
      itmp = itp_sub
      itp_sub = itc_sub
      itc_sub = itmp
    endif
  enddo
  flds(i_ng)%op(:,:,:,itc(i_ng)) = op_i(:,:,:,itc_sub)
  flds(i_ng)%op_nm(:,:,:,itc(i_ng)) = op_nm_i(:,:,:,itc_sub)

  call elden_ng( &
    flds(i_ng)%nplus, &
    flds(i_ng)%n2p, &
    flds(i_ng)%nop, &
    flds(i_ng)%o2p(:,:,:,itc(i_ng)), &
    flds(i_ng)%ne(:,:,:,itc(i_ng)), &
    i_ng)

  call addfld(flds(i_ng)%nplus,'NP_diag',i_ng)
  call addfld(flds(i_ng)%n2p,'N2P_diag',i_ng)
  call addfld(flds(i_ng)%nop,'NOP_diag',i_ng)

  call lamdas_ng( &
    flds(i_ng)%lxx, &
    flds(i_ng)%lyy, &
    flds(i_ng)%lxy, &
    flds(i_ng)%lyx, &
    flds(i_ng)%lam1, &
    flds(i_ng)%ped, &
    flds(i_ng)%hall, &
    flds(i_ng)%Q1, &
    flds(i_ng)%Q2, &
    i_ng)

  call comp_n2d_ng( &
    flds(i_ng)%n2d(:,:,:,itc(i_ng)), &
    flds(i_ng)%n2d_nm(:,:,:,itc(i_ng)), &
    istep,i_ng)

  call comp_n4s_ng( &
    flds(i_ng)%n4s(:,:,:,itc(i_ng)), &
    flds(i_ng)%n4s_nm(:,:,:,itc(i_ng)), &
    istep,i_ng)

  call comp_no_ng( &
    flds(i_ng)%no(:,:,:,itc(i_ng)), &
    flds(i_ng)%no_nm(:,:,:,itc(i_ng)), &
    istep,i_ng)

  call comp_ar_ng( &
    flds(i_ng)%ar(:,:,:,itc(i_ng)), &
    flds(i_ng)%ar_nm(:,:,:,itc(i_ng)), &
    istep,i_ng)

  call qjnno_ng( &
    flds(i_ng)%qtotal, &
    i_ng)

  call qjion_ng( &
    flds(i_ng)%qtotal, &
    i_ng)

  call qjoule_ti_ng( &
    flds(i_ng)%qji_ti, &
    i_ng)

  call settei_ng( &
    flds(i_ng)%te(:,:,:,itc(i_ng)), &
    flds(i_ng)%ti(:,:,:,itc(i_ng)), &
    flds(i_ng)%qtotal, &
    i_ng)

  call newton_ng( &
    flds(i_ng)%cool_implicit, &
    flds(i_ng)%cool_explicit, &
    i_ng)

  call hdif3_ng( &
    flds(i_ng)%hdt, &
    flds(i_ng)%hdu, &
    flds(i_ng)%hdv, &
    flds(i_ng)%hdo2, &
    flds(i_ng)%hdo1, &
    flds(i_ng)%hdhe, &
    i_ng)

  call qjoule_tn_ng( &
    flds(i_ng)%qji_tn, &
    i_ng)

  call dt_ng( &
    flds(i_ng)%tn(:,:,:,itc(i_ng)), &
    flds(i_ng)%tn_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%tlbc, &
    flds(i_ng)%tlbc_nm, &
    istep,i_ng)

  call duv_ng( &
    flds(i_ng)%un(:,:,:,itc(i_ng)), &
    flds(i_ng)%un_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%vn(:,:,:,itc(i_ng)), &
    flds(i_ng)%vn_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%ulbc, &
    flds(i_ng)%vlbc, &
    flds(i_ng)%ulbc_nm, &
    flds(i_ng)%vlbc_nm, &
    istep,i_ng)

  call comp_o2o_ng( &
    flds(i_ng)%fs, &
    i_ng)

  call comp_ng( &
    flds(i_ng)%o2(:,:,:,itc(i_ng)), &
    flds(i_ng)%o2_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%o1(:,:,:,itc(i_ng)), &
    flds(i_ng)%o1_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%he(:,:,:,itc(i_ng)), &
    flds(i_ng)%he_nm(:,:,:,itc(i_ng)), &
    flds(i_ng)%flx_he, &
    istep,i_ng)

end subroutine dynamics_ng

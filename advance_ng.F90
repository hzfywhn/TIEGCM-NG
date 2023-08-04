recursive subroutine advance_ng(i_ng)
! instead of linear procedure in advance, advance_ng has a stack-like working procedure
! enter low levels first before high levels, but finish high levels before low levels
! a recursive subroutine is suitable for such stack-like procedures

  use params_module,only: n_ng,nlevp1_ng
  use input_module,only: nstep_ng,nstep_sub,nudge_level,nudge_lbc,nudge_f4d,nudge_alpha
  use fields_ng_module,only: flds,nf3din,nf2din,itp,itc,modeltime,step
  use level_outer_ng_module,only: map_in_1=>map_in,map_out_1=>map_out,set_bndry_1=>set_bndry
  use levels_inner_ng_module,only: map_in,map_out,set_bndry
  use nudge_ng_module,only: nudge_fields=>flds,nlb,nf4d,time,itime,wrap,lb_idx,f4d_idx
  use output_ng_module,only: addfld
  implicit none

  integer,intent(in) :: i_ng

  integer :: istep,ifld,nk,k,latbeg,latend,lonbeg,lonend,offbeg,offend,lonbeg1,lonend1,itmp
  real :: delta,fac1,fac2
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: omega,omegai,wn
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tmp_lbc
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,nudge_fields(i_ng)%latbeg:nudge_fields(i_ng)%latend) :: wt,ext_fld
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,nudge_fields(i_ng)%latbeg:nudge_fields(i_ng)%latend,nlb) :: nclbc
  real,dimension(nudge_fields(i_ng)%maxlev,flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: tmp_f4d
  real,dimension(nudge_fields(i_ng)%maxlev,flds(i_ng)%lond0:flds(i_ng)%lond1, &
    nudge_fields(i_ng)%latbeg:nudge_fields(i_ng)%latend,nf4d) :: ncf4d
  external :: interp_fields_ng,addiag_ng,hdif12_ng,dynamics_ng

  nk = nlevp1_ng(i_ng)

! return beyond inner most level
  if (i_ng >= n_ng+1) return

! set up essential 2d/3d fields and boundaries
  if (i_ng == 1) then
    call map_in_1
    call set_bndry_1
  else
    call map_in(i_ng)
    call set_bndry(i_ng)
  endif

! index 0 is either from initialization or the previous iteration
  call interp_fields_ng(i_ng)

  if (nudge_lbc .or. nudge_f4d) then
    latbeg = nudge_fields(i_ng)%latbeg
    latend = nudge_fields(i_ng)%latend
    lonbeg = nudge_fields(i_ng)%lonbeg
    lonend = nudge_fields(i_ng)%lonend
    offbeg = nudge_fields(i_ng)%offbeg
    offend = nudge_fields(i_ng)%offend
    lonbeg1 = nudge_fields(i_ng)%lonbeg1
    lonend1 = nudge_fields(i_ng)%lonend1
  endif

  do istep = 1,nstep_ng(i_ng)
    modeltime(i_ng) = modeltime(i_ng)+step(i_ng)

! each sub-cycle starts with replacing lower boundary with external fields
    if (nudge_lbc) then
      if (nudge_level(i_ng) .and. latbeg<latend .and. &
        time(itime)<=modeltime(i_ng) .and. modeltime(i_ng)<=time(itime+1)) then

        delta = time(itime+1)-time(itime)
        fac1 = (time(itime+1)-modeltime(i_ng))/delta
        fac2 = (modeltime(i_ng)-time(itime))/delta
        nclbc = fac1*nudge_fields(i_ng)%lbc(:,:,1,:)+fac2*nudge_fields(i_ng)%lbc(:,:,2,:)

        wt = nudge_alpha*nudge_fields(i_ng)%hori_weight

        do ifld = 1,nlb
          if (lb_idx(ifld) /= 0) then
            tmp_lbc = flds(i_ng)%f2d_save(:,:,istep,ifld)
            ext_fld = wt*nclbc(:,:,ifld)+(1-wt)*tmp_lbc(:,latbeg:latend)

            if (wrap) then
              tmp_lbc(:,latbeg:latend) = ext_fld
            else
              if (lonbeg<lonend .and. offbeg==offend) &
                tmp_lbc(lonbeg:lonend,latbeg:latend) = ext_fld(lonbeg:lonend,:)

              if (lonbeg<lonend .and. offbeg+1==offend) then
                tmp_lbc(lonbeg:lonend1,latbeg:latend) = ext_fld(lonbeg:lonend1,:)
                tmp_lbc(lonbeg1:lonend,latbeg:latend) = ext_fld(lonbeg1:lonend,:)
              endif

              if (lonbeg>lonend .and. offbeg==offend+1) then
                tmp_lbc(lonbeg:flds(i_ng)%lond1,latbeg:latend) = ext_fld(lonbeg:flds(i_ng)%lond1,:)
                tmp_lbc(flds(i_ng)%lond0:lonend,latbeg:latend) = ext_fld(flds(i_ng)%lond0:lonend,:)
              endif
            endif

            flds(i_ng)%f2d_save(:,:,istep,ifld) = tmp_lbc
          endif
        enddo
      endif
    endif

! essential 2d/3d fields (lower boundary conditions and electric fields) are stored in f2d_save/f3d_save
! before each sub-cycle, they will replace the outdated f2din/f3din fields in the previous sub-cycle
    do ifld = 1,nf3din
      flds(i_ng)%f3din(ifld)%data = flds(i_ng)%f3d_save(:,:,:,istep,ifld)
    enddo
    do ifld = 1,nf2din
      flds(i_ng)%f2din(ifld)%data = flds(i_ng)%f2d_save(:,:,istep,ifld)
    enddo

    call addfld(flds(i_ng)%t_lbc,'TLBC',i_ng)
    call addfld(flds(i_ng)%u_lbc,'ULBC',i_ng)
    call addfld(flds(i_ng)%v_lbc,'VLBC',i_ng)
    call addfld(flds(i_ng)%z_lbc,'ZLBC',i_ng)

    call addiag_ng( &
      flds(i_ng)%vc(:,:,:,itp(i_ng)), &
      flds(i_ng)%mbar(:,:,:,itp(i_ng)), &
      flds(i_ng)%barm(:,:,:,itp(i_ng)), &
      flds(i_ng)%xnmbar(:,:,:,itp(i_ng)), &
      flds(i_ng)%xnmbari(:,:,:,itp(i_ng)), &
      flds(i_ng)%scht(:,:,:,itp(i_ng)), &
      flds(i_ng)%schti(:,:,:,itp(i_ng)), &
      flds(i_ng)%z(:,:,:,itp(i_ng)), &
      flds(i_ng)%zg, &
      flds(i_ng)%n2, &
      i_ng)

    flds(i_ng)%z(:,:,:,itc(i_ng)) = flds(i_ng)%z(:,:,:,itp(i_ng))

    call hdif12_ng( &
      flds(i_ng)%kldt, &
      flds(i_ng)%kldu, &
      flds(i_ng)%kldv, &
      flds(i_ng)%kldo2, &
      flds(i_ng)%kldo1, &
      flds(i_ng)%kldhe, &
      flds(i_ng)%fnrh, &
      flds(i_ng)%fkmh, &
      i_ng)

    call dynamics_ng(istep,i_ng)

! calc_wn
    omega = flds(i_ng)%w(:,:,:,itp(i_ng))
    do k = 1,nk-1
      omegai(k,:,:) = 0.5*(omega(k,:,:)+omega(k+1,:,:))
    enddo
    omegai(nk,:,:) = 1.5*omega(nk,:,:)-0.5*omega(nk-1,:,:)
    wn = omegai*flds(i_ng)%scht(:,:,:,itp(i_ng))
    call addfld(wn,'WN',i_ng)

    if (nudge_f4d) then
      if (nudge_level(i_ng) .and. latbeg<latend .and. &
        time(itime)<=modeltime(i_ng) .and. modeltime(i_ng)<=time(itime+1)) then

        delta = time(itime+1)-time(itime)
        fac1 = (time(itime+1)-modeltime(i_ng))/delta
        fac2 = (modeltime(i_ng)-time(itime))/delta
        ncf4d = fac1*nudge_fields(i_ng)%f4d(:,:,:,1,:)+fac2*nudge_fields(i_ng)%f4d(:,:,:,2,:)

        do ifld = 1,nf4d
          tmp_f4d = flds(i_ng)%f4d(f4d_idx(ifld))%data(1:nudge_fields(i_ng)%maxlev,:,:,itc(i_ng))

          do k = 1,nudge_fields(i_ng)%maxlev
            wt = nudge_alpha*nudge_fields(i_ng)%hori_weight*nudge_fields(i_ng)%vert_weight(k,ifld)
            ext_fld = wt*ncf4d(k,:,:,ifld)+(1-wt)*tmp_f4d(k,:,latbeg:latend)

            if (wrap) then
              tmp_f4d(k,:,latbeg:latend) = ext_fld
            else
              if (lonbeg<lonend .and. offbeg==offend) &
                tmp_f4d(k,lonbeg:lonend,latbeg:latend) = ext_fld(lonbeg:lonend,:)

              if (lonbeg<lonend .and. offbeg+1==offend) then
                tmp_f4d(k,lonbeg:lonend1,latbeg:latend) = ext_fld(lonbeg:lonend1,:)
                tmp_f4d(k,lonbeg1:lonend,latbeg:latend) = ext_fld(lonbeg1:lonend,:)
              endif

              if (lonbeg>lonend .and. offbeg==offend+1) then
                tmp_f4d(k,lonbeg:flds(i_ng)%lond1,latbeg:latend) = ext_fld(lonbeg:flds(i_ng)%lond1,:)
                tmp_f4d(k,flds(i_ng)%lond0:lonend,latbeg:latend) = ext_fld(flds(i_ng)%lond0:lonend,:)
              endif
            endif
          enddo

          flds(i_ng)%f4d(f4d_idx(ifld))%data(1:nudge_fields(i_ng)%maxlev,:,:,itc(i_ng)) = tmp_f4d
        enddo
      endif
    endif

! go to the next level before swapping itp and itc
    call advance_ng(i_ng+1)

! fields are updated with the next level after advance_ng, then swap itp and itc
    itmp = itp(i_ng)
    itp(i_ng) = itc(i_ng)
    itc(i_ng) = itmp
  enddo

! save the last step for interpolation in the next iteration
  flds(i_ng)%f3d_save(:,:,:,0,:) = flds(i_ng)%f3d_save(:,:,:,nstep_ng(i_ng),:)
  flds(i_ng)%f2d_save(:,:,0,:) = flds(i_ng)%f2d_save(:,:,nstep_ng(i_ng),:)
  if (flds(i_ng)%is_bndry(1) .or. flds(i_ng)%is_bndry(2)) &
    flds(i_ng)%lon_b(:,:,:,0,:) = flds(i_ng)%lon_b(:,:,:,nstep_ng(i_ng),:)
  if (flds(i_ng)%is_bndry(3) .or. flds(i_ng)%is_bndry(4)) &
    flds(i_ng)%lat_b(:,:,:,0,:) = flds(i_ng)%lat_b(:,:,:,nstep_ng(i_ng),:)

! note that only 4d fields at time index itc are mapped out
  if (i_ng == 1) then
    call map_out_1
  else
    call map_out(i_ng)
  endif

end subroutine advance_ng

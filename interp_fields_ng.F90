subroutine interp_fields_ng(i_ng)
! interpolate essential 2d/3d fields and boundaries in time

  use params_module,only: nlevp1_ng
  use input_module,only: nstep_ng,nstep_sub
  use fields_ng_module,only: flds,nf3din,nf2din,nbnd,bndry
  use interp_module,only: interp1d
  implicit none

  integer,intent(in) :: i_ng

  integer :: t1,t1_sub,k,i,lat,it,ifld
  real,dimension(2) :: fsave
  real,dimension(1:nstep_ng(i_ng)-1) :: t
  real,dimension(0:nstep_ng(i_ng)*nstep_sub) :: t_sub

  t1 = nstep_ng(i_ng)
  t1_sub = nstep_ng(i_ng)*nstep_sub
  t = (/(it,it=1,t1-1)/)
  t_sub = (/(it,it=0,t1_sub)/)

! index 0 is from initialization or previous iteration
! index nstep is from boundary mapping from global fields
! in case of round-off errors, index 0 and nstep are excluded from interpolation
  do lat = flds(i_ng)%latd0,flds(i_ng)%latd1
    do i = flds(i_ng)%lond0,flds(i_ng)%lond1
      do ifld = 1,nf3din
        do k = 1,nlevp1_ng(i_ng)
          fsave = (/flds(i_ng)%f3d_save(k,i,lat,0,ifld),flds(i_ng)%f3d_save(k,i,lat,t1,ifld)/)
          flds(i_ng)%f3d_save(k,i,lat,1:t1-1,ifld) = interp1d(t,(/0.,real(t1)/),fsave)
        enddo
      enddo
      do ifld = 1,nf2din
        fsave = (/flds(i_ng)%f2d_save(i,lat,0,ifld),flds(i_ng)%f2d_save(i,lat,t1,ifld)/)
        flds(i_ng)%f2d_save(i,lat,1:t1-1,ifld) = interp1d(t,(/0.,real(t1)/),fsave)
      enddo
    enddo
  enddo

  do ifld = 1,nbnd
    do k = 1,nlevp1_ng(i_ng)
      do i = 1,2
        if (flds(i_ng)%is_bndry(i)) then
          do lat = flds(i_ng)%latd0,flds(i_ng)%latd1
            fsave = (/flds(i_ng)%lon_b(k,i,lat,0,ifld),flds(i_ng)%lon_b(k,i,lat,t1,ifld)/)
            flds(i_ng)%lon_b(k,i,lat,1:t1-1,ifld) = interp1d(t,(/0.,real(t1)/),fsave)
            if (trim(bndry(ifld)) == 'OP') &
              flds(i_ng)%op_lon_b(k,i,lat,:) = interp1d(t_sub,(/0.,real(t1_sub)/),fsave)
          enddo
        endif
      enddo

      do lat = 3,4
        if (flds(i_ng)%is_bndry(lat)) then
          do i = flds(i_ng)%lond0,flds(i_ng)%lond1
            fsave = (/flds(i_ng)%lat_b(k,i,lat,0,ifld),flds(i_ng)%lat_b(k,i,lat,t1,ifld)/)
            flds(i_ng)%lat_b(k,i,lat,1:t1-1,ifld) = interp1d(t,(/0.,real(t1)/),fsave)
            if (trim(bndry(ifld)) == 'OP') &
              flds(i_ng)%op_lat_b(k,i,lat,:) = interp1d(t_sub,(/0.,real(t1_sub)/),fsave)
          enddo
        endif
      enddo
    enddo
  enddo

end subroutine interp_fields_ng

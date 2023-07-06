module output_ng_module
! a simplified output module for nested grid fields
! the layout of the output fields is (lon,lat,lev,time) compared to the internal field (lev,lon,lat,time)

  use params_module,only: mx_ng
  use fields_module,only: shortname_len
  implicit none

  logical :: defined

  integer,dimension(mx_ng) :: ncid,dim_lev,dim_ilev,dim_lon,dim_lat,dim_time,varid_time,irec

! this number can be adjusted if there are a lot of output fields
  integer,parameter :: maxvarout = 30

! all levels have the same set of output fields
  integer :: nf2d,nf3d
  character(len=shortname_len),dimension(maxvarout) :: f2d_name,f3d_name

! fields at different levels have different dimensions,
! and they may correspond to different IDs in the output file,
! therefore they are saved separately at each level
  type fields_out
    integer,dimension(maxvarout) :: f2d_id,f3d_id
    real,dimension(:,:,:),allocatable :: f2d
    real,dimension(:,:,:,:),allocatable :: f3d
  end type fields_out

! output buffer
  type(fields_out),dimension(mx_ng) :: fout

  interface addfld
    module procedure addfld_2d,addfld_3d
  end interface addfld

  contains
!-----------------------------------------------------------------------
  subroutine init

    use params_module,only: nlevp1_ng,n_ng,nlon_ng,nlat_ng,zpmid_ng,zpint_ng,glon_ng,glat_ng
    use input_module,only: start_year,start_day,fileout_ng,varout_ng
    use fields_module,only: nf4d,f4d
    use fields_ng_module,only: flds
    use char_module,only: ismember
    use mpi_module,only: mytid
    use netcdf,only: nf90_netcdf4,nf90_unlimited,nf90_double, &
      nf90_create,nf90_def_dim,nf90_def_var,nf90_enddef,nf90_put_var,nf90_put_att

    integer :: month,day,if4d,i_ng,stat,varid_lon,varid_lat,varid_lev,varid_ilev
    character(len=80) :: units
    external :: to_month_day

    irec = 1
    call to_month_day(start_year,start_day,month,day)
    write(units,"('seconds since ',i4,2('-',i2.2),' 00:00:00')") start_year,month,day

    nf2d = 0
    nf3d = 0

! copy 4d fields to the output buffer
    do if4d = 1,nf4d
      if (ismember(f4d(if4d)%short_name,varout_ng)) then
        nf3d = nf3d+1
        f3d_name(nf3d) = trim(f4d(if4d)%short_name)
      endif
    enddo

    do i_ng = 1,n_ng
      allocate(fout(i_ng)%f2d(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,maxvarout))
      allocate(fout(i_ng)%f3d(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,maxvarout))
    enddo

! fields are not defined in the output file until the first call of output
    defined = .false.

! root only: define coordinates at each level
    if (mytid == 0) then
      do i_ng = 1,n_ng
        stat = nf90_create(trim(fileout_ng(i_ng)),nf90_netcdf4,ncid(i_ng))

        stat = nf90_def_dim(ncid(i_ng),'lon',nlon_ng(i_ng)+4,dim_lon(i_ng))
        stat = nf90_def_dim(ncid(i_ng),'lat',nlat_ng(i_ng)+4,dim_lat(i_ng))
        stat = nf90_def_dim(ncid(i_ng),'lev',nlevp1_ng(i_ng),dim_lev(i_ng))
        stat = nf90_def_dim(ncid(i_ng),'ilev',nlevp1_ng(i_ng),dim_ilev(i_ng))
        stat = nf90_def_dim(ncid(i_ng),'time',nf90_unlimited,dim_time(i_ng))

        stat = nf90_def_var(ncid(i_ng),'lon',nf90_double,dim_lon(i_ng),varid_lon)
        stat = nf90_def_var(ncid(i_ng),'lat',nf90_double,dim_lat(i_ng),varid_lat)
        stat = nf90_def_var(ncid(i_ng),'lev',nf90_double,dim_lev(i_ng),varid_lev)
        stat = nf90_def_var(ncid(i_ng),'ilev',nf90_double,dim_ilev(i_ng),varid_ilev)
        stat = nf90_def_var(ncid(i_ng),'time',nf90_double,dim_time(i_ng),varid_time(i_ng))

        stat = nf90_enddef(ncid(i_ng))

        stat = nf90_put_var(ncid(i_ng),varid_lon,glon_ng(i_ng,-1:nlon_ng(i_ng)+2))
        stat = nf90_put_var(ncid(i_ng),varid_lat,glat_ng(i_ng,-1:nlat_ng(i_ng)+2))

        stat = nf90_put_var(ncid(i_ng),varid_lev,zpmid_ng(i_ng,1:nlevp1_ng(i_ng)))
        stat = nf90_put_var(ncid(i_ng),varid_ilev,zpint_ng(i_ng,1:nlevp1_ng(i_ng)))

        stat = nf90_put_att(ncid(i_ng),varid_lon,'long_name','geographic longitude (-west, +east)')
        stat = nf90_put_att(ncid(i_ng),varid_lon,'units','degrees_east')
        stat = nf90_put_att(ncid(i_ng),varid_lat,'long_name','geographic latitude (-south, +north)')
        stat = nf90_put_att(ncid(i_ng),varid_lat,'units','degrees_north')
        stat = nf90_put_att(ncid(i_ng),varid_lev,'long_name','midpoint levels')
        stat = nf90_put_att(ncid(i_ng),varid_lev,'short_name','ln(p0/p)')
        stat = nf90_put_att(ncid(i_ng),varid_lev,'units','')
        stat = nf90_put_att(ncid(i_ng),varid_ilev,'long_name','interface levels')
        stat = nf90_put_att(ncid(i_ng),varid_ilev,'short_name','ln(p0/p)')
        stat = nf90_put_att(ncid(i_ng),varid_ilev,'units','')
        stat = nf90_put_att(ncid(i_ng),varid_time(i_ng),'long_name','time')
        stat = nf90_put_att(ncid(i_ng),varid_time(i_ng),'units',trim(units))
      enddo
    endif

  end subroutine init
!-----------------------------------------------------------------------
  subroutine output

    use params_module,only: n_ng,nlevp1_ng,nlon_ng,nlat_ng
    use fields_module,only: f4d
    use fields_ng_module,only: flds,itc,modeltime
    use char_module,only: find_index
    use mpi_module,only: mytid
    use gather2root_ng_module,only: gather2root_vars2d,gather2root_vars3d
    use netcdf,only: nf90_float,nf90_redef,nf90_def_var,nf90_enddef,nf90_put_att,nf90_put_var,nf90_sync
    use mpi_f08,only: mpi_barrier,mpi_comm_world

    integer :: i_ng,stat,if2d,if3d,if4d,ilev,dim_v
    real,dimension(-1:maxval(nlon_ng)+2,-1:maxval(nlat_ng)+2,nf2d) :: full2d
    real,dimension(maxval(nlevp1_ng),-1:maxval(nlon_ng)+2,-1:maxval(nlat_ng)+2,nf3d) :: full3d
    real(kind=4),dimension(-1:maxval(nlon_ng)+2,-1:maxval(nlat_ng)+2) :: fout2d
    real(kind=4),dimension(-1:maxval(nlon_ng)+2,-1:maxval(nlat_ng)+2,maxval(nlevp1_ng)) :: fout3d

! define output fields in the output files, this finishes initialization
    if (.not. defined) then
      do i_ng = 1,n_ng
        stat = nf90_redef(ncid(i_ng))

        do if2d = 1,nf2d
          stat = nf90_def_var(ncid(i_ng),trim(f2d_name(if2d)),nf90_float, &
            (/dim_lon(i_ng),dim_lat(i_ng),dim_time(i_ng)/),fout(i_ng)%f2d_id(if2d))
        enddo

        do if3d = 1,nf3d
          dim_v = dim_lev(i_ng)
          if4d = find_index(f3d_name(if3d),f4d%short_name)
          if (if4d /= 0) then
            if (trim(f4d(if4d)%vcoord) == 'interfaces') dim_v = dim_ilev(i_ng)
          endif
          stat = nf90_def_var(ncid(i_ng),trim(f3d_name(if3d)),nf90_float, &
            (/dim_lon(i_ng),dim_lat(i_ng),dim_v,dim_time(i_ng)/),fout(i_ng)%f3d_id(if3d))
        enddo

        stat = nf90_enddef(ncid(i_ng))

        do if3d = 1,nf3d
          if4d = find_index(f3d_name(if3d),f4d%short_name)
          if (if4d /= 0) then
            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'long_name',f4d(if4d)%long_name)
            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'short_name',f4d(if4d)%short_name)
            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'units',f4d(if4d)%units)
          endif
        enddo
      enddo

      defined = .true.
    endif

! fields through addfld have been updated, leaving 4d fields not updated
    do i_ng = 1,n_ng
      do if3d = 1,nf3d
        if4d = find_index(f3d_name(if3d),f4d%short_name)
        if (if4d /= 0) fout(i_ng)%f3d(:,:,:,if3d) = flds(i_ng)%f4d(if4d)%data(:,:,:,itc(i_ng))
      enddo
    enddo

    do i_ng = 1,n_ng
      call gather2root_vars2d(fout(i_ng)%f2d(:,:,1:nf2d),full2d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,:),nf2d,i_ng)
      call gather2root_vars3d(fout(i_ng)%f3d(:,:,:,1:nf3d), &
        full3d(1:nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,:),nf3d,i_ng)

! root only: write to files
      if (mytid == 0) then
        stat = nf90_put_var(ncid(i_ng),varid_time(i_ng),modeltime(i_ng),start=(/irec(i_ng)/))

        do if2d = 1,nf2d
          fout2d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2) = full2d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,if2d)
          stat = nf90_put_var(ncid(i_ng),fout(i_ng)%f2d_id(if2d),fout2d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2), &
            start=(/1,1,irec(i_ng)/),count=(/nlon_ng(i_ng)+4,nlat_ng(i_ng)+4,1/))
        enddo

        do if3d = 1,nf3d
          do ilev = 1,nlevp1_ng(i_ng)
            fout3d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,ilev) = full3d(ilev,-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,if3d)
          enddo
          stat = nf90_put_var(ncid(i_ng),fout(i_ng)%f3d_id(if3d),fout3d(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,1:nlevp1_ng(i_ng)), &
            start=(/1,1,1,irec(i_ng)/),count=(/nlon_ng(i_ng)+4,nlat_ng(i_ng)+4,nlevp1_ng(i_ng),1/))
        enddo

! flush output buffer in case the program quited abnormally
        stat = nf90_sync(ncid(i_ng))

        irec(i_ng) = irec(i_ng)+1
      endif

      call mpi_barrier(mpi_comm_world)
    enddo

  end subroutine output
!-----------------------------------------------------------------------
  subroutine addfld_2d(var,varname,i_ng)

    use input_module,only: varout_ng
    use fields_ng_module,only: flds
    use char_module,only: ismember,find_index

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    character(len=*),intent(in) :: varname

    integer :: if2d
    external :: shutdown

    if (ismember(varname,varout_ng)) then
      if2d = find_index(varname,f2d_name)
      if (if2d == 0) then
        nf2d = nf2d+1
        if (nf2d > maxvarout) call shutdown('2d output fields exceed upper limit')
        f2d_name(nf2d) = varname
        if2d = nf2d
      endif
      fout(i_ng)%f2d(:,:,if2d) = var
    endif

  end subroutine addfld_2d
!-----------------------------------------------------------------------
  subroutine addfld_3d(var,varname,i_ng)

    use params_module,only: nlevp1_ng
    use input_module,only: varout_ng
    use fields_ng_module,only: flds
    use char_module,only: ismember,find_index

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    character(len=*),intent(in) :: varname

    integer :: if3d
    external :: shutdown

    if (ismember(varname,varout_ng)) then
      if3d = find_index(varname,f3d_name)
      if (if3d == 0) then
        nf3d = nf3d+1
        if (nf3d > maxvarout) call shutdown('3d output fields exceed upper limit')
        f3d_name(nf3d) = varname
        if3d = nf3d
      endif
      fout(i_ng)%f3d(:,:,:,if3d) = var
    endif

  end subroutine addfld_3d
!-----------------------------------------------------------------------
  subroutine finalize

    use params_module,only: n_ng
    use mpi_module,only: mytid
    use netcdf,only: nf90_close

    integer :: i_ng,stat

    if (mytid == 0) then
      do i_ng = 1,n_ng
        stat = nf90_close(ncid(i_ng))
      enddo
    endif

  end subroutine finalize
!-----------------------------------------------------------------------
end module output_ng_module

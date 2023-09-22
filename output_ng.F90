module output_ng_module
! a simplified output module for nested grid fields, parallel write
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

! f2d is purely from addfld
    real,dimension(:,:,:),allocatable :: f2d

! f3d also includes internal f4d fields
    real,dimension(:,:,:,:),allocatable :: f3d
  end type fields_out

! output buffer
  type(fields_out),dimension(mx_ng) :: fout

! netcdf operations (define/write) are done in output, addfld only put the array into buffer
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
    use mpi_module,only: TIEGCM_WORLD
    use mpi
    use netcdf,only: nf90_create,nf90_def_dim,nf90_def_var,nf90_put_att,nf90_var_par_access,nf90_enddef, &
      nf90_put_var,nf90_strerror,nf90_netcdf4,nf90_mpiio,nf90_unlimited,nf90_double,nf90_collective,nf90_noerr

    integer :: month,day,if4d,i_ng,stat,varid_lon,varid_lat,varid_lev,varid_ilev
    character(len=80) :: units
    external :: to_month_day,shutdown

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
      allocate(fout(i_ng)%f2d(flds(i_ng)%lon0:flds(i_ng)%lon1,flds(i_ng)%lat0:flds(i_ng)%lat1,maxvarout))
      allocate(fout(i_ng)%f3d(nlevp1_ng(i_ng),flds(i_ng)%lon0:flds(i_ng)%lon1,flds(i_ng)%lat0:flds(i_ng)%lat1,maxvarout))
    enddo

! fields are not defined in the output file until the first call of output
    defined = .false.

! define coordinates at each level
    do i_ng = 1,n_ng
      stat = nf90_create(trim(fileout_ng(i_ng)),ior(nf90_netcdf4,nf90_mpiio),ncid(i_ng),comm=TIEGCM_WORLD,info=mpi_info_null)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at creating '//trim(fileout_ng(i_ng)))

      stat = nf90_def_dim(ncid(i_ng),'lon',nlon_ng(i_ng),dim_lon(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining dimension lon')

      stat = nf90_def_dim(ncid(i_ng),'lat',nlat_ng(i_ng),dim_lat(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining dimension lat')

      stat = nf90_def_dim(ncid(i_ng),'lev',nlevp1_ng(i_ng),dim_lev(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining dimension lev')

      stat = nf90_def_dim(ncid(i_ng),'ilev',nlevp1_ng(i_ng),dim_ilev(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining dimension ilev')

      stat = nf90_def_dim(ncid(i_ng),'time',nf90_unlimited,dim_time(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining dimension time')

      stat = nf90_def_var(ncid(i_ng),'lon',nf90_double,dim_lon(i_ng),varid_lon)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable lon')

      stat = nf90_put_att(ncid(i_ng),varid_lon,'long_name','geographic longitude (-west, +east)')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lon.long_name')

      stat = nf90_put_att(ncid(i_ng),varid_lon,'units','degrees_east')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lon.units')

      stat = nf90_def_var(ncid(i_ng),'lat',nf90_double,dim_lat(i_ng),varid_lat)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable lat')

      stat = nf90_put_att(ncid(i_ng),varid_lat,'long_name','geographic latitude (-south, +north)')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lat.long_name')

      stat = nf90_put_att(ncid(i_ng),varid_lat,'units','degrees_north')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lat.units')

      stat = nf90_def_var(ncid(i_ng),'lev',nf90_double,dim_lev(i_ng),varid_lev)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable lev')

      stat = nf90_put_att(ncid(i_ng),varid_lev,'long_name','midpoint levels')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lev.long_name')

      stat = nf90_put_att(ncid(i_ng),varid_lev,'short_name','ln(p0/p)')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lev.short_name')

      stat = nf90_put_att(ncid(i_ng),varid_lev,'units','')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute lev.units')

      stat = nf90_def_var(ncid(i_ng),'ilev',nf90_double,dim_ilev(i_ng),varid_ilev)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable ilev')

      stat = nf90_put_att(ncid(i_ng),varid_ilev,'long_name','interface levels')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute ilev.long_name')

      stat = nf90_put_att(ncid(i_ng),varid_ilev,'short_name','ln(p0/p)')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute ilev.short_name')

      stat = nf90_put_att(ncid(i_ng),varid_ilev,'units','')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute ilev.units')

      stat = nf90_def_var(ncid(i_ng),'time',nf90_double,dim_time(i_ng),varid_time(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable time')

! variable with unlimited dimension requires collective operation
      stat = nf90_var_par_access(ncid(i_ng),varid_time(i_ng),nf90_collective)
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at setting variable time collective')

      stat = nf90_put_att(ncid(i_ng),varid_time(i_ng),'long_name','time')
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute time.long_name')

      stat = nf90_put_att(ncid(i_ng),varid_time(i_ng),'units',trim(units))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting attribute time.units')

      stat = nf90_enddef(ncid(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at exiting def mode')

      stat = nf90_put_var(ncid(i_ng),varid_lon,glon_ng(i_ng,1:nlon_ng(i_ng)))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable lon')

      stat = nf90_put_var(ncid(i_ng),varid_lat,glat_ng(i_ng,1:nlat_ng(i_ng)))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable lat')

      stat = nf90_put_var(ncid(i_ng),varid_lev,zpmid_ng(i_ng,1:nlevp1_ng(i_ng)))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable lev')

      stat = nf90_put_var(ncid(i_ng),varid_ilev,zpint_ng(i_ng,1:nlevp1_ng(i_ng)))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable ilev')
    enddo

  end subroutine init
!-----------------------------------------------------------------------
  subroutine output

    use params_module,only: n_ng,nlevp1_ng
    use input_module,only: fileout_ng
    use fields_module,only: f4d
    use fields_ng_module,only: flds,itc,modeltime
    use char_module,only: find_index
    use netcdf,only: nf90_redef,nf90_def_var,nf90_var_par_access,nf90_put_att,nf90_enddef, &
      nf90_put_var,nf90_sync,nf90_strerror,nf90_float,nf90_collective,nf90_noerr

    integer :: i_ng,stat,if2d,if3d,if4d,ilev,dim_v,lon0,lon1,lat0,lat1
    real(kind=4),dimension(minval(flds%lon0):maxval(flds%lon1),minval(flds%lat0):maxval(flds%lat1)) :: fout2d
    real(kind=4),dimension(minval(flds%lon0):maxval(flds%lon1),minval(flds%lat0):maxval(flds%lat1),maxval(nlevp1_ng)) :: fout3d
    external :: shutdown

! define output fields in the output files, this finishes initialization
    if (.not. defined) then
      do i_ng = 1,n_ng
        stat = nf90_redef(ncid(i_ng))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at entering def mode')

        do if2d = 1,nf2d
          stat = nf90_def_var(ncid(i_ng),trim(f2d_name(if2d)),nf90_float, &
            (/dim_lon(i_ng),dim_lat(i_ng),dim_time(i_ng)/),fout(i_ng)%f2d_id(if2d))
          if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable '//trim(f2d_name(if2d)))

! variable with unlimited dimension requires collective operation
          stat = nf90_var_par_access(ncid(i_ng),fout(i_ng)%f2d_id(if2d),nf90_collective)
          if (stat /= nf90_noerr) &
            call shutdown(trim(nf90_strerror(stat))//' at setting variable '//trim(f2d_name(if2d))//' collective')
        enddo

        do if3d = 1,nf3d
! set default vertical coordinate to be midpoints
          dim_v = dim_lev(i_ng)

! if it is an internal f4d field, vertical coordinate can be retrieved from fields module
          if4d = find_index(f3d_name(if3d),f4d%short_name)
          if (if4d /= 0) then
            if (trim(f4d(if4d)%vcoord) == 'interfaces') dim_v = dim_ilev(i_ng)
          endif

          stat = nf90_def_var(ncid(i_ng),trim(f3d_name(if3d)),nf90_float, &
            (/dim_lon(i_ng),dim_lat(i_ng),dim_v,dim_time(i_ng)/),fout(i_ng)%f3d_id(if3d))
          if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at defining variable '//trim(f3d_name(if3d)))

! variable with unlimited dimension requires collective operation
          stat = nf90_var_par_access(ncid(i_ng),fout(i_ng)%f3d_id(if3d),nf90_collective)
          if (stat /= nf90_noerr) &
            call shutdown(trim(nf90_strerror(stat))//' at setting variable '//trim(f3d_name(if2d))//' collective')

! if it is an internal f4d field, there will be additional associated information in fields module
          if (if4d /= 0) then
            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'long_name',f4d(if4d)%long_name)
            if (stat /= nf90_noerr) &
              call shutdown(trim(nf90_strerror(stat))//' at putting attribute '//trim(f3d_name(if3d))//'.long_name')

            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'short_name',f4d(if4d)%short_name)
            if (stat /= nf90_noerr) &
              call shutdown(trim(nf90_strerror(stat))//' at putting attribute '//trim(f3d_name(if3d))//'.short_name')

            stat = nf90_put_att(ncid(i_ng),fout(i_ng)%f3d_id(if3d),'units',f4d(if4d)%units)
            if (stat /= nf90_noerr) &
              call shutdown(trim(nf90_strerror(stat))//' at putting attribute '//trim(f3d_name(if3d))//'.units')
          endif
        enddo

        stat = nf90_enddef(ncid(i_ng))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at exiting def mode')
      enddo

      defined = .true.
    endif

! write to files
    do i_ng = 1,n_ng
      lon0 = flds(i_ng)%lon0
      lon1 = flds(i_ng)%lon1
      lat0 = flds(i_ng)%lat0
      lat1 = flds(i_ng)%lat1

      stat = nf90_put_var(ncid(i_ng),varid_time(i_ng),modeltime(i_ng),start=(/irec(i_ng)/))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable time')

      do if2d = 1,nf2d
! real*8 (nf90_double) -> real*4 (nf90_float)
        fout2d(lon0:lon1,lat0:lat1) = fout(i_ng)%f2d(:,:,if2d)

        stat = nf90_put_var(ncid(i_ng),fout(i_ng)%f2d_id(if2d),fout2d(lon0:lon1,lat0:lat1), &
          start=(/lon0,lat0,irec(i_ng)/),count=(/lon1-lon0+1,lat1-lat0+1,1/))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable '//trim(f2d_name(if2d)))
      enddo

      do if3d = 1,nf3d
! fields through addfld have been updated, leaving internal 4d fields not updated, so update them now
        if4d = find_index(f3d_name(if3d),f4d%short_name)
        if (if4d /= 0) fout(i_ng)%f3d(:,:,:,if3d) = flds(i_ng)%f4d(if4d)%data(:,lon0:lon1,lat0:lat1,itc(i_ng))

! real*8 (nf90_double) -> real*4 (nf90_float)
        do ilev = 1,nlevp1_ng(i_ng)
          fout3d(lon0:lon1,lat0:lat1,ilev) = fout(i_ng)%f3d(ilev,:,:,if3d)
        enddo

        stat = nf90_put_var(ncid(i_ng),fout(i_ng)%f3d_id(if3d),fout3d(lon0:lon1,lat0:lat1,1:nlevp1_ng(i_ng)), &
          start=(/lon0,lat0,1,irec(i_ng)/),count=(/lon1-lon0+1,lat1-lat0+1,nlevp1_ng(i_ng),1/))
        if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at putting variable '//trim(f3d_name(if3d)))
      enddo

! flush output buffer in case the program quited abnormally
      stat = nf90_sync(ncid(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at syncing '//trim(fileout_ng(i_ng)))

      irec(i_ng) = irec(i_ng)+1
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
      fout(i_ng)%f2d(:,:,if2d) = var(flds(i_ng)%lon0:flds(i_ng)%lon1,flds(i_ng)%lat0:flds(i_ng)%lat1)
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
      fout(i_ng)%f3d(:,:,:,if3d) = var(:,flds(i_ng)%lon0:flds(i_ng)%lon1,flds(i_ng)%lat0:flds(i_ng)%lat1)
    endif

  end subroutine addfld_3d
!-----------------------------------------------------------------------
  subroutine finalize

    use params_module,only: n_ng
    use input_module,only: fileout_ng
    use netcdf,only: nf90_close,nf90_strerror,nf90_noerr

    integer :: i_ng,stat
    external :: shutdown

    do i_ng = 1,n_ng
      stat = nf90_close(ncid(i_ng))
      if (stat /= nf90_noerr) call shutdown(trim(nf90_strerror(stat))//' at closing '//trim(fileout_ng(i_ng)))
    enddo

  end subroutine finalize
!-----------------------------------------------------------------------
end module output_ng_module

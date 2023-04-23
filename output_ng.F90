module output_ng_module
! a simplified output module for nested grid fields
! output fields remain the same layout as the internal field (lev,lon,lat,time)

  use params_module,only: mx_ng
  use fields_module,only: nf4d,shortname_len
  implicit none

  integer,parameter :: max_addfld = 30

  integer :: n_varout
  logical,dimension(nf4d) :: varout_flag
  integer,dimension(mx_ng) :: ncid,dim_lev,dim_lon,dim_lat,dim_time,varid_time,irec,nvar2d,nvar3d
  integer,dimension(mx_ng,nf4d) :: varout_id,var2d_id,var3d_id
  character(len=shortname_len),dimension(mx_ng,max_addfld) :: var2d_name,var3d_name

  interface addfld
    module procedure addfld_2d,addfld_3d
  end interface addfld

  contains
!-----------------------------------------------------------------------
  subroutine init

    use params_module,only: nlevp1_ng,n_ng,nlon_ng,nlat_ng,zpmid_ng,zpint_ng,glon_ng,glat_ng
    use input_module,only: fileout_ng,varout_ng,start_year,start_day,dlev_ng
    use fields_module,only: f4d
    use char_module,only: ismember
    use mpi_module,only: mytid
    use netcdf,only: nf90_netcdf4,nf90_unlimited,nf90_double,nf90_float, &
      nf90_create,nf90_def_dim,nf90_def_var,nf90_enddef,nf90_put_var,nf90_put_att

    integer :: if4d,n,imon,day,i_ng,stat,dim_ilev,dim_v,varid_lev,varid_ilev,varid_lon,varid_lat
    integer,dimension(12) :: month
    character(len=80) :: units

    n_varout = 0
    do if4d = 1,nf4d
      if (ismember(f4d(if4d)%short_name,varout_ng)) then
        varout_flag(if4d) = .true.
        n_varout = n_varout+1
      else
        varout_flag(if4d) = .false.
      endif
    enddo

    if (mytid /= 0) return

    month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    if ((mod(start_year,4)==0 .and. mod(start_year,100)/=0) .or. mod(start_year,400)==0) month(2) = 29

    imon = 1
    day = start_day
    do while (day > month(imon))
      day = day-month(imon)
      imon = imon+1
    enddo

    write(units,"('seconds since ',i4,'-',i2.2,'-',i2.2,' 00:00:00')") start_year,imon,day

    irec = 1
    nvar2d = 0
    nvar3d = 0
    var2d_id = 0
    var3d_id = 0
    var2d_name = ''
    var3d_name = ''

    do i_ng = 1,n_ng
      stat = nf90_create(trim(fileout_ng(i_ng)),nf90_netcdf4,ncid(i_ng))

      stat = nf90_def_dim(ncid(i_ng),'lon',nlon_ng(i_ng)+4,dim_lon(i_ng))
      stat = nf90_def_dim(ncid(i_ng),'lat',nlat_ng(i_ng)+4,dim_lat(i_ng))
      stat = nf90_def_dim(ncid(i_ng),'lev',nlevp1_ng(i_ng),dim_lev(i_ng))
      stat = nf90_def_dim(ncid(i_ng),'ilev',nlevp1_ng(i_ng),dim_ilev)
      stat = nf90_def_dim(ncid(i_ng),'time',nf90_unlimited,dim_time(i_ng))

      stat = nf90_def_var(ncid(i_ng),'lon',nf90_double,dim_lon(i_ng),varid_lon)
      stat = nf90_def_var(ncid(i_ng),'lat',nf90_double,dim_lat(i_ng),varid_lat)
      stat = nf90_def_var(ncid(i_ng),'lev',nf90_double,dim_lev(i_ng),varid_lev)
      stat = nf90_def_var(ncid(i_ng),'ilev',nf90_double,dim_ilev,varid_ilev)
      stat = nf90_def_var(ncid(i_ng),'time',nf90_double,dim_time(i_ng),varid_time(i_ng))

      n = 1
      do if4d = 1,nf4d
        if (varout_flag(if4d)) then
          if (trim(f4d(if4d)%vcoord) == 'midpoints') then
            dim_v = dim_lev(i_ng)
          else
            dim_v = dim_ilev
          endif
          stat = nf90_def_var(ncid(i_ng),trim(f4d(if4d)%short_name),nf90_float, &
            (/dim_lon(i_ng),dim_lat(i_ng),dim_v,dim_time(i_ng)/),varout_id(i_ng,n))
          n = n+1
        endif
      enddo

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

      n = 1
      do if4d = 1,nf4d
        if (varout_flag(if4d)) then
          stat = nf90_put_att(ncid(i_ng),varout_id(i_ng,n),'long_name',f4d(if4d)%long_name)
          stat = nf90_put_att(ncid(i_ng),varout_id(i_ng,n),'short_name',f4d(if4d)%short_name)
          stat = nf90_put_att(ncid(i_ng),varout_id(i_ng,n),'units',f4d(if4d)%units)
          n = n+1
        endif
      enddo
    enddo

  end subroutine init
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
  subroutine output(i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
    use input_module,only: start_day
    use fields_ng_module,only: flds,itc,modeltime
    use mpi_module,only: mytid
    use gather2root_ng_module,only: gather2root_vars3d
    use netcdf,only: nf90_put_var,nf90_sync

    integer,intent(in) :: i_ng

    integer :: if4d,n,ilev,stat
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1,n_varout) :: vars
    real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,n_varout) :: full
    real(kind=4),dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,nlevp1_ng(i_ng)) :: fout

    n = 1
    do if4d = 1,nf4d
      if (varout_flag(if4d)) then
        vars(:,:,:,n) = flds(i_ng)%f4d(if4d)%data(:,:,:,itc(i_ng))
        n = n+1
      endif
    enddo
    call gather2root_vars3d(vars,full,n_varout,i_ng)

    if (mytid == 0) then
      do if4d = 1,n_varout
        do ilev = 1,nlevp1_ng(i_ng)
          fout(:,:,ilev) = full(ilev,:,:,if4d)
        enddo
        stat = nf90_put_var(ncid(i_ng),varout_id(i_ng,if4d),fout, &
          start=(/1,1,1,irec(i_ng)/),count=(/nlon_ng(i_ng)+4,nlat_ng(i_ng)+4,nlevp1_ng(i_ng),1/))
      enddo

      stat = nf90_put_var(ncid(i_ng),varid_time(i_ng),modeltime(i_ng),start=(/irec(i_ng)/))
      irec(i_ng) = irec(i_ng)+1

      stat = nf90_sync(ncid(i_ng))
    endif

  end subroutine output
!-----------------------------------------------------------------------
  subroutine addfld_2d(var,varname,i_ng)

    use params_module,only: nlon_ng,nlat_ng
    use fields_ng_module,only: flds
    use char_module,only: ismember,find_index
    use mpi_module,only: mytid
    use gather2root_ng_module,only: gather2root_var2d
    use netcdf,only: nf90_double,nf90_redef,nf90_def_var,nf90_enddef,nf90_put_var

    integer,intent(in) :: i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    character(len=*),intent(in) :: varname

    integer :: if2d,stat
    real,dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2) :: full
    real(kind=4),dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2) :: fout

    call gather2root_var2d(var,full,i_ng)

    if (mytid == 0) then
      if (ismember(varname,var2d_name(i_ng,1:nvar2d(i_ng)))) then
        if2d = find_index(varname,var2d_name(i_ng,1:nvar2d(i_ng)))
      else
        if2d = nvar2d(i_ng)+1
        nvar2d(i_ng) = if2d
        var2d_name(i_ng,if2d) = trim(varname)
        stat = nf90_redef(ncid(i_ng))
        stat = nf90_def_var(ncid(i_ng),trim(varname),nf90_double, &
          (/dim_lon(i_ng),dim_lat(i_ng),dim_time(i_ng)/),var2d_id(i_ng,if2d))
        stat = nf90_enddef(ncid(i_ng))
      endif
      fout = full
      stat = nf90_put_var(ncid(i_ng),var2d_id(i_ng,if2d),fout, &
        start=(/1,1,irec(i_ng)/),count=(/nlon_ng(i_ng)+4,nlat_ng(i_ng)+4,1/))
    endif

  end subroutine addfld_2d
!-----------------------------------------------------------------------
  subroutine addfld_3d(var,varname,i_ng)

    use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
    use fields_ng_module,only: flds
    use char_module,only: ismember,find_index
    use mpi_module,only: mytid
    use gather2root_ng_module,only: gather2root_var3d
    use netcdf,only: nf90_double,nf90_redef,nf90_def_var,nf90_enddef,nf90_put_var

    integer,intent(in) :: i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: var
    character(len=*),intent(in) :: varname

    integer :: if3d,ilev,stat
    real,dimension(nlevp1_ng(i_ng),-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2) :: full
    real(kind=4),dimension(-1:nlon_ng(i_ng)+2,-1:nlat_ng(i_ng)+2,nlevp1_ng(i_ng)) :: fout

    call gather2root_var3d(var,full,i_ng)

    if (mytid == 0) then
      if (ismember(varname,var3d_name(i_ng,1:nvar3d(i_ng)))) then
        if3d = find_index(varname,var3d_name(i_ng,1:nvar3d(i_ng)))
      else
        if3d = nvar3d(i_ng)+1
        nvar3d(i_ng) = if3d
        var3d_name(i_ng,if3d) = trim(varname)
        stat = nf90_redef(ncid(i_ng))
        stat = nf90_def_var(ncid(i_ng),trim(varname),nf90_double, &
          (/dim_lon(i_ng),dim_lat(i_ng),dim_lev(i_ng),dim_time(i_ng)/),var3d_id(i_ng,if3d))
        stat = nf90_enddef(ncid(i_ng))
      endif
      do ilev = 1,nlevp1_ng(i_ng)
        fout(:,:,ilev) = full(ilev,:,:)
      enddo
      stat = nf90_put_var(ncid(i_ng),var3d_id(i_ng,if3d),fout, &
        start=(/1,1,1,irec(i_ng)/),count=(/nlon_ng(i_ng)+4,nlat_ng(i_ng)+4,nlevp1_ng(i_ng),1/))
    endif

  end subroutine addfld_3d
!-----------------------------------------------------------------------
end module output_ng_module

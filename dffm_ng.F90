module dffm_ng_module
! calculate finite difference and centered mean

  implicit none

  contains
!-----------------------------------------------------------------------
  subroutine df_2d(f,df,order,i_ng)
! two point centered derivative in lon (order=1,2) or lat (order=-1,-2)
! normalized by dividing the spacing between two points
! edges use one-sided derivatives

    use fields_ng_module,only: flds
    use sync_ng_module,only: sync_var2d_lon,sync_var2d_lat

    integer,intent(in) :: order,i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: df

    integer :: lond0,lond1,latd0,latd1,i,lat
    external :: shutdown

    lond0 = flds(i_ng)%lond0
    lond1 = flds(i_ng)%lond1
    latd0 = flds(i_ng)%latd0
    latd1 = flds(i_ng)%latd1

    select case (order)
      case (0)
        df = 0.

      case (1,2)
        df(lond0,:) = f(lond0+1,:)-f(lond0,:)
        df(lond1,:) = f(lond1,:)-f(lond1-1,:)
        if (order == 2) then
          df(lond0+1,:) = (f(lond0+2,:)-f(lond0,:))/2.
          df(lond1-1,:) = (f(lond1,:)-f(lond1-2,:))/2.
        endif
        do i = lond0+order,lond1-order
          df(i,:) = (f(i+order,:)-f(i-order,:))/(order*2)
        enddo
        call sync_var2d_lon(df,i_ng)

      case (-1,-2)
        df(:,latd0) = f(:,latd0+1)-f(:,latd0)
        df(:,latd1) = f(:,latd1)-f(:,latd1-1)
        if (order == -2) then
          df(:,latd0+1) = (f(:,latd0+2)-f(:,latd0))/2.
          df(:,latd1-1) = (f(:,latd1)-f(:,latd1-2))/2.
        endif
        do lat = latd0-order,latd1+order
          df(:,lat) = (f(:,lat-order)-f(:,lat+order))/(-order*2)
        enddo
        call sync_var2d_lat(df,i_ng)

      case default
        call shutdown('Illegal order in df_2d')
    endselect

  end subroutine df_2d
!-----------------------------------------------------------------------
  subroutine df_3d(f,df,order,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds
    use sync_ng_module,only: sync_var3d_lon,sync_var3d_lat

    integer,intent(in) :: order,i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: df

    integer :: lond0,lond1,latd0,latd1,i,lat
    external :: shutdown

    lond0 = flds(i_ng)%lond0
    lond1 = flds(i_ng)%lond1
    latd0 = flds(i_ng)%latd0
    latd1 = flds(i_ng)%latd1

    select case (order)
      case (0)
        df = 0.

      case (1,2)
        df(:,lond0,:) = f(:,lond0+1,:)-f(:,lond0,:)
        df(:,lond1,:) = f(:,lond1,:)-f(:,lond1-1,:)
        if (order == 2) then
          df(:,lond0+1,:) = (f(:,lond0+2,:)-f(:,lond0,:))/2.
          df(:,lond1-1,:) = (f(:,lond1,:)-f(:,lond1-2,:))/2.
        endif
        do i = lond0+order,lond1-order
          df(:,i,:) = (f(:,i+order,:)-f(:,i-order,:))/(order*2)
        enddo
        call sync_var3d_lon(df,i_ng)

      case (-1,-2)
        df(:,:,latd0) = f(:,:,latd0+1)-f(:,:,latd0)
        df(:,:,latd1) = f(:,:,latd1)-f(:,:,latd1-1)
        if (order == -2) then
          df(:,:,latd0+1) = (f(:,:,latd0+2)-f(:,:,latd0))/2.
          df(:,:,latd1-1) = (f(:,:,latd1)-f(:,:,latd1-2))/2.
        endif
        do lat = latd0-order,latd1+order
          df(:,:,lat) = (f(:,:,lat-order)-f(:,:,lat+order))/(-order*2)
        enddo
        call sync_var3d_lat(df,i_ng)

      case default
        call shutdown('Illegal order in df_3d')
    endselect

  end subroutine df_3d
!-----------------------------------------------------------------------
  subroutine fm_2d(f,fm,order,i_ng)
! two point centered mean in lon (order=1,2) or lat (order=-1,-2)
! edges use the input itself

    use fields_ng_module,only: flds
    use sync_ng_module,only: sync_var2d_lon,sync_var2d_lat

    integer,intent(in) :: order,i_ng
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: fm

    integer :: lond0,lond1,latd0,latd1,i,lat
    external :: shutdown

    lond0 = flds(i_ng)%lond0
    lond1 = flds(i_ng)%lond1
    latd0 = flds(i_ng)%latd0
    latd1 = flds(i_ng)%latd1

    select case (order)
      case (0)
        fm = f

      case (1,2)
        fm(lond0,:) = f(lond0,:)
        fm(lond1,:) = f(lond1,:)
        if (order == 2) then
          fm(lond0+1,:) = (f(lond0+2,:)+f(lond0,:))/2.
          fm(lond1-1,:) = (f(lond1,:)+f(lond1-2,:))/2.
        endif
        do i = lond0+order,lond1-order
          fm(i,:) = (f(i+order,:)+f(i-order,:))/2.
        enddo
        call sync_var2d_lon(fm,i_ng)

      case (-1,-2)
        fm(:,latd0) = f(:,latd0)
        fm(:,latd1) = f(:,latd1)
        if (order == -2) then
          fm(:,latd0+1) = (f(:,latd0+2)+f(:,latd0))/2.
          fm(:,latd1-1) = (f(:,latd1)+f(:,latd1-2))/2.
        endif
        do lat = latd0-order,latd1+order
          fm(:,lat) = (f(:,lat-order)+f(:,lat+order))/2.
        enddo
        call sync_var2d_lat(fm,i_ng)

      case default
        call shutdown('Illegal order in fm_2d')
    endselect

  end subroutine fm_2d
!-----------------------------------------------------------------------
  subroutine fm_3d(f,fm,order,i_ng)

    use params_module,only: nlevp1_ng
    use fields_ng_module,only: flds
    use sync_ng_module,only: sync_var3d_lon,sync_var3d_lat

    integer,intent(in) :: order,i_ng
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: fm

    integer :: lond0,lond1,latd0,latd1,i,lat
    external :: shutdown

    lond0 = flds(i_ng)%lond0
    lond1 = flds(i_ng)%lond1
    latd0 = flds(i_ng)%latd0
    latd1 = flds(i_ng)%latd1

    select case (order)
      case (0)
        fm = f

      case (1,2)
        fm(:,lond0,:) = f(:,lond0,:)
        fm(:,lond1,:) = f(:,lond1,:)
        if (order == 2) then
          fm(:,lond0+1,:) = (f(:,lond0+2,:)+f(:,lond0,:))/2.
          fm(:,lond1-1,:) = (f(:,lond1,:)+f(:,lond1-2,:))/2.
        endif
        do i = lond0+order,lond1-order
          fm(:,i,:) = (f(:,i+order,:)+f(:,i-order,:))/2.
        enddo
        call sync_var3d_lon(fm,i_ng)

      case (-1,-2)
        fm(:,:,latd0) = f(:,:,latd0)
        fm(:,:,latd1) = f(:,:,latd1)
        if (order == -2) then
          fm(:,:,latd0+1) = (f(:,:,latd0+2)+f(:,:,latd0))/2.
          fm(:,:,latd1-1) = (f(:,:,latd1)+f(:,:,latd1-2))/2.
        endif
        do lat = latd0-order,latd1+order
          fm(:,:,lat) = (f(:,:,lat-order)+f(:,:,lat+order))/2.
        enddo
        call sync_var3d_lat(fm,i_ng)

      case default
        call shutdown('Illegal order in fm_3d')
    endselect

  end subroutine fm_3d
!-----------------------------------------------------------------------
end module dffm_ng_module

subroutine div_ng(f,df,dir,i_ng)
! horizontal divergence in lon (dir=1) or lat (dir=-1)
! normalized by dividing the spacing between two points
! edges use lower order derivatives

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use fields_ng_module,only: flds
  use dffm_ng_module,only: df_3d
  implicit none

  integer,intent(in) :: dir,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: f
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: df

  integer :: nx,ny
  logical,dimension(4) :: is_bndry
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: dfdx1,dfdx2,dfdy1,dfdy2
  external :: shutdown

  is_bndry = flds(i_ng)%is_bndry

  select case (dir)
    case (1)
      nx = nlon_ng(i_ng)
      call df_3d(f,dfdx1,1,i_ng)
      call df_3d(f,dfdx2,2,i_ng)
      df = dfdx1*2.*2./3.-dfdx2*4./12.
      if (is_bndry(1)) df(:,-1:0,:) = dfdx1(:,-1:0,:)
      if (is_bndry(2)) df(:,nx+1:nx+2,:) = dfdx1(:,nx+1:nx+2,:)

    case (-1)
      ny = nlat_ng(i_ng)
      call df_3d(f,dfdy1,-1,i_ng)
      call df_3d(f,dfdy2,-2,i_ng)
      df = dfdy1*2.*2./3.-dfdy2*4./12.
      if (is_bndry(3)) df(:,:,-1:0) = dfdy1(:,:,-1:0)
      if (is_bndry(4)) df(:,:,ny+1:ny+2) = dfdy1(:,:,ny+1:ny+2)

    case default
      call shutdown('Illegal dir in div_ng')
  endselect

end subroutine div_ng

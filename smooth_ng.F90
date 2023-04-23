subroutine smooth_ng(fin,fout,shapiro,i_ng)
! 5 point stencil at middle points, 3 point stencil at boundaries

  use params_module,only: nlevp1_ng,nlon_ng,nlat_ng
  use fields_ng_module,only: flds
  use dffm_ng_module,only: fm_3d
  implicit none

  integer,intent(in) :: i_ng
  real,intent(in) :: shapiro
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: fin
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: fout

  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: ftmp,fm1,fm2

  call fm_3d(fin,fm1,-1,i_ng)
  call fm_3d(fin,fm2,-2,i_ng)
  ftmp = fin-shapiro*(2.*fm2-8.*fm1+6.*fin)

  if (flds(i_ng)%is_bndry(3)) then
    ftmp(:,:,0) = fin(:,:,0)-shapiro*(-8.*fm1(:,:,0)+8.*fin(:,:,0))
    ftmp(:,:,-1) = fin(:,:,-1)
  endif
  if (flds(i_ng)%is_bndry(4)) then
    ftmp(:,:,nlat_ng(i_ng)+1) = fin(:,:,nlat_ng(i_ng)+1)-shapiro*(-8.*fm1(:,:,nlat_ng(i_ng)+1)+8.*fin(:,:,nlat_ng(i_ng)+1))
    ftmp(:,:,nlat_ng(i_ng)+2) = fin(:,:,nlat_ng(i_ng)+2)
  endif

  call fm_3d(ftmp,fm1,1,i_ng)
  call fm_3d(ftmp,fm2,2,i_ng)
  fout = ftmp-shapiro*(2.*fm2-8.*fm1+6.*ftmp)

  if (flds(i_ng)%is_bndry(1)) then
    fout(:,0,:) = ftmp(:,0,:)-shapiro*(-8.*fm1(:,0,:)+8.*ftmp(:,0,:))
    fout(:,-1,:) = ftmp(:,-1,:)
  endif
  if (flds(i_ng)%is_bndry(2)) then
    fout(:,nlon_ng(i_ng)+1,:) = ftmp(:,nlon_ng(i_ng)+1,:)-shapiro*(-8.*fm1(:,nlon_ng(i_ng)+1,:)+8.*ftmp(:,nlon_ng(i_ng)+1,:))
    fout(:,nlon_ng(i_ng)+2,:) = ftmp(:,nlon_ng(i_ng)+2,:)
  endif

end subroutine smooth_ng

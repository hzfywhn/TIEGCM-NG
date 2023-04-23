subroutine trsolv_ng(a,b,c,f,x,k2,i_ng)

  use params_module,only: nlevp1_ng
  use fields_ng_module,only: flds
  implicit none

  integer,intent(in) :: k2,i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(in) :: a,b,c,f
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: x

  integer :: k
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: w1,w2,w3

  w1(1,:,:) = b(1,:,:)
  do k = 2,k2
    w2(k-1,:,:) = c(k-1,:,:)/w1(k-1,:,:)
    w1(k,:,:) = b(k,:,:)-a(k,:,:)*w2(k-1,:,:)
  enddo
  w3(1,:,:) = f(1,:,:)/w1(1,:,:)
  do k = 2,k2
    w3(k,:,:) = (f(k,:,:)-a(k,:,:)*w3(k-1,:,:))/w1(k,:,:)
  enddo
  x(k2,:,:) = w3(k2,:,:)
  do k = k2-1,1,-1
    x(k,:,:) = w3(k,:,:)-w2(k,:,:)*x(k+1,:,:)
  enddo

end subroutine trsolv_ng

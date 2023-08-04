subroutine newton_ng(cool_implicit,cool_explicit,i_ng)

  use params_module,only: nlevp1_ng,zpint_ng,zibot
  use cons_module,only: rmass_co2,p0,boltz,avo,rmassinv_o2,rmassinv_o1,rmassinv_n2,rmassinv_no
  use init_module,only: iday,iyear,secs
  use fields_ng_module,only: flds,itp,dz,expzmid_inv,expz
  use interp_module,only: interp1d
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    cool_implicit,cool_explicit

  real,dimension(3),parameter :: an = (/0.835E-18,0.6,0.2/), bn = (/228.,228.,325./)
  integer :: nk,k
  real :: co2u,uttime,yfrac
  real,dimension(nlevp1_ng(i_ng)) :: xfac,xfaci
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    tn,o2,o1,n2,no,barm,xnmbar,cp,xmco2,nco2,aco2,bco2,ano,co2_cool,no_cool,fac1,fac2,cpi,fac1_a
  real,dimension(29),parameter :: zpint_5 = (/(zibot+(k-1)*0.5, k=1,29)/), &
    xfac_5 = (/0.1000E-01, 0.1000E-01, 0.1000E-01, 0.5000E-01, 0.1000E+00, &
      0.2000E+00, 0.4000E+00, 0.5500E+00, 0.7000E+00, 0.7500E+00, &
      (0.8000E+00, k=1,19)/)

  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  barm = flds(i_ng)%barm(:,:,:,itp(i_ng))
  xnmbar = flds(i_ng)%xnmbar(:,:,:,itp(i_ng))
  cp = flds(i_ng)%cp

  nk = nlevp1_ng(i_ng)

  uttime = secs/3600.
  if (mod(iyear,4) == 0) then
    yfrac = iyear+(iday-1+uttime/24.)/366.
  else
    yfrac = iyear+(iday-1+uttime/24.)/365.
  endif
  if (yfrac > 1953.0) then
    co2u = 248.0+30.0*exp(0.014*(yfrac-1900.))
  else
    co2u = 280.0+0.30*(yfrac-1850.0)
  endif
  co2u = max(co2u,280.0)*1.e-6

  xmco2(1,:,:) = .25*rmass_co2*(1./barm(1,:,:)+1./barm(2,:,:))*dz(i_ng)
  do k = 2,nk
    xmco2(k,:,:) = xmco2(k-1,:,:)+rmass_co2/barm(k,:,:)*dz(i_ng)
  enddo

  nco2 = co2u*p0*expz(i_ng,1)*expzmid_inv(i_ng)/(boltz*tn)*exp(-xmco2)

  aco2 = 2.5e-15*(1.+0.03*max(tn-200.,0.))

  where (tn < 260.)
    bco2 = 1.56e-12
  elsewhere (tn>=260. .and. tn<=300.)
    bco2 = (2.6-0.004*tn)*1.0e-12
  elsewhere
    bco2 = 1.4e-12
  endwhere

  do k = 1,nk-1
    cpi(k,:,:) = .5*(cp(k,:,:)+cp(k+1,:,:))
  enddo
  cpi(nk,:,:) = 1.5*cp(nk,:,:)-.5*cp(nk-1,:,:)

  co2_cool = 2.65e-13*nco2*exp(-960./tn)*avo*((o2*rmassinv_o2+n2*rmassinv_n2)*aco2+o1*rmassinv_o1*bco2)

  ano = xnmbar*(4.2e-11*o1*rmassinv_o1+2.4e-14*o2*rmassinv_o2)

  no_cool = 4.956e-12*avo*no*rmassinv_no*ano/(ano+13.3)*exp(-2700./tn)

  cool_implicit = (2700.*no_cool+960.*co2_cool)/(cpi*tn**2)

  cool_explicit = (1.-2700./tn)*no_cool+(1.-960./tn)*co2_cool
  do k = 1,nk
    cool_explicit(k,:,:) = cool_explicit(k,:,:)*expz(i_ng,k)
  enddo

  co2_cool(nk,:,:) = co2_cool(nk-1,:,:)
  no_cool(nk,:,:) = no_cool(nk-1,:,:)

! newto3p
  xfac = interp1d(zpint_ng(i_ng,1:nk),zpint_5,xfac_5)

  do k = 1,nk-1
    xfaci(k) = .5*(xfac(k)+xfac(k+1))
  enddo
  xfaci(nk) = 1.5*xfac(nk)-.5*xfac(nk-1)

  fac1 = an(1)*avo*rmassinv_o1*o1*exp(-bn(1)/tn)
  do k = 1,nk
    fac1(k,:,:) = fac1(k,:,:)*xfaci(k)
  enddo
  fac2 = 1.+an(2)*exp(-bn(2)/tn)+an(3)*exp(-bn(3)/tn)
  fac1 = fac1/fac2
  fac2 = fac1/fac2/tn**2*(bn(1)+(bn(1)-bn(3))*an(3)*exp(-bn(3)/tn))

  cool_implicit = cool_implicit+fac2/cpi

  fac1 = fac1-tn*fac2
  do k = 1,nk
    fac1_a(k,:,:) = fac1(k,:,:)*expz(i_ng,k)
  enddo
  cool_explicit = cool_explicit+fac1_a

end subroutine newton_ng

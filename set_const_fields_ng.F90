subroutine set_const_fields_ng
! set up constant fields, called once per run

  implicit none

  call set_const
  call set_magfield

  contains
!-----------------------------------------------------------------------
  subroutine set_const
! grid related parameters

    use params_module,only: n_ng,nlevp1_ng,zibot,glat_ng
    use cons_module,only: default_step,smooth_fac,dtr,re,difk,dift,xmue
    use input_module,only: step,dlev_ng,dlat_ng,dlon_ng,nstep_ng
    use lbc,only: fb,b
    use fields_ng_module,only: flds,fb_ng=>fb,b_ng=>b, &
      step_ng=>step,shapiro,dtx2inv, &
      dlev,dphi,dlamda,dz,dzp,expzmid,expzmid_inv,expz, &
      difk_ng=>difk,dift_ng=>dift,xmue_ng=>xmue

    real,parameter :: omega = 7.292E-5
    integer :: i_ng,k,latd0,latd1
    real :: expdz

    fb_ng = fb(1,:)
    b_ng = b(1,:,:)

    step_ng(1) = step/nstep_ng(1)
    do i_ng = 2,n_ng
      step_ng(i_ng) = step_ng(i_ng-1)/nstep_ng(i_ng)
    enddo
    shapiro = step_ng/default_step*smooth_fac
    dtx2inv = 1./(2.*step_ng)

    dlev = 0.
    dphi = 0.
    dlamda = 0.
    dz = 0.
    dzp = 0.
    expzmid = 0.
    expzmid_inv = 0.
    expz = 0.
    difk_ng = 0.
    dift_ng = 0.
    xmue_ng = 0.

    do i_ng = 1,n_ng
      dlev(i_ng) = dlev_ng(i_ng)
      dphi(i_ng) = dlat_ng(i_ng)*dtr
      dlamda(i_ng) = dlon_ng(i_ng)*dtr
      dz(i_ng) = dlev_ng(i_ng)
      dzp(i_ng) = dlev_ng(i_ng)
      expzmid(i_ng) = exp(-.5*dlev_ng(i_ng))
      expzmid_inv(i_ng) = 1./expzmid(i_ng)

      expdz = exp(-dlev_ng(i_ng))
      difk_ng(i_ng,1,:) = difk(1,:)
      dift_ng(i_ng,1,:) = dift(1,:)
      xmue_ng(i_ng,1,:) = xmue(1,:)
      expz(i_ng,1) = exp(-zibot-.5*dlev_ng(i_ng))
      do k = 2,nlevp1_ng(i_ng)
        expz(i_ng,k) = expz(i_ng,k-1)*expdz
        difk_ng(i_ng,k,:) = difk_ng(i_ng,k-1,:)*expdz
        dift_ng(i_ng,k,:) = dift_ng(i_ng,k-1,:)*expdz
        if (k == nlevp1_ng(i_ng)) then
          xmue_ng(i_ng,k,:) = difk_ng(i_ng,k-1,:)
        else
          xmue_ng(i_ng,k,:) = difk_ng(i_ng,k,:)
        endif
      enddo

      latd0 = flds(i_ng)%latd0
      latd1 = flds(i_ng)%latd1
      flds(i_ng)%cs = cos(glat_ng(i_ng,latd0:latd1)*dtr)
      flds(i_ng)%tanphi = tan(glat_ng(i_ng,latd0:latd1)*dtr)
      flds(i_ng)%cor = 2.*omega*sin(glat_ng(i_ng,latd0:latd1)*dtr)
      flds(i_ng)%racs = 1./(re*flds(i_ng)%cs)
    enddo

  end subroutine set_const
!-----------------------------------------------------------------------
  subroutine set_magfield
! magnetic field related parameters

    use params_module,only: n_ng,glon_ng,glat_ng
    use cons_module,only: hs,dtr
    use input_module,only: dlat_ng
    use apex,only: geolon,nglon,apex_mall
    use fields_ng_module,only: flds,dipmin

    integer :: i_ng,i,lat,ist
    real :: alt,bmag,alon,xlatm,vmp,w,d,be3,sim,xlatqd,f,si,x
    real,dimension(2) :: f1,f2
    real,dimension(3) :: b,bhat,d1,d2,d3,e1,e2,e3

    alt = hs*1.e-5

    do i_ng = 1,n_ng
      do lat = flds(i_ng)%latd0,flds(i_ng)%latd1
        do i = flds(i_ng)%lond0,flds(i_ng)%lond1
          if (glon_ng(i_ng,i) <= geolon(nglon)) then
            x = glon_ng(i_ng,i)
          else
            x = glon_ng(i_ng,i)-360
          endif
          call apex_mall(glat_ng(i_ng,lat),x,alt,alt,b,bhat, &
            bmag,si,alon,xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3,xlatqd,f,f1,f2,ist)
          flds(i_ng)%alatm(i,lat) = xlatm*dtr
          flds(i_ng)%alonm(i,lat) = alon*dtr
          flds(i_ng)%xb(i,lat) = b(2)*1.e-5
          flds(i_ng)%yb(i,lat) = b(1)*1.e-5
          flds(i_ng)%zb(i,lat) = -b(3)*1.e-5
          flds(i_ng)%bmod(i,lat) = bmag*1.e-5
          flds(i_ng)%rjac(i,lat,1,1) = f2(2)
          flds(i_ng)%rjac(i,lat,1,2) = -f2(1)
          flds(i_ng)%rjac(i,lat,2,1) = -f1(2)
          flds(i_ng)%rjac(i,lat,2,2) = f1(1)
        enddo
      enddo

      flds(i_ng)%rlatm = flds(i_ng)%alatm
      flds(i_ng)%rlonm = flds(i_ng)%alonm
      flds(i_ng)%dipmag = atan(flds(i_ng)%zb/sqrt(flds(i_ng)%xb**2+flds(i_ng)%yb**2))
      flds(i_ng)%decmag = atan2(flds(i_ng)%yb,flds(i_ng)%xb)
      flds(i_ng)%sndec = sin(flds(i_ng)%decmag)
      flds(i_ng)%csdec = cos(flds(i_ng)%decmag)
      flds(i_ng)%sn2dec = flds(i_ng)%sndec**2
      flds(i_ng)%bx = flds(i_ng)%yb/flds(i_ng)%bmod
      flds(i_ng)%by = flds(i_ng)%xb/flds(i_ng)%bmod
      flds(i_ng)%bz = -flds(i_ng)%zb/flds(i_ng)%bmod
      flds(i_ng)%bmod2 = flds(i_ng)%bmod

      dipmin(i_ng) = sin(dlat_ng(i_ng)*2*dtr)
      where (abs(flds(i_ng)%bz) < dipmin(i_ng))
        flds(i_ng)%bx = flds(i_ng)%bx*sqrt(1.-dipmin(i_ng)**2)/sqrt(1.-flds(i_ng)%bz**2)
        flds(i_ng)%by = flds(i_ng)%by*sqrt(1.-dipmin(i_ng)**2)/sqrt(1.-flds(i_ng)%bz**2)
        flds(i_ng)%bz = sign(dipmin(i_ng),flds(i_ng)%bz)
      endwhere
    enddo

  end subroutine set_magfield
!-----------------------------------------------------------------------
end subroutine set_const_fields_ng

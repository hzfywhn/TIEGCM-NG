subroutine aurora_ng(qteaur,qo2p,qop,qn2p,qnp,qtef,i_ng)

  use params_module,only: nlevp1_ng,glat_ng
  use fields_ng_module,only: flds
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: qteaur
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(inout) :: &
    qo2p,qop,qn2p,qnp,qtef

  real,parameter :: aurlat = 30.
  integer :: latd0,latd1,latbeg,latend

  latd0 = flds(i_ng)%latd0
  latd1 = flds(i_ng)%latd1
  qteaur = 0.

! relations between glat(latd0), glat(latd1) and +/- aurlat
! should be included in the following 5 situations
  if (glat_ng(i_ng,latd0)>=-aurlat .and. glat_ng(i_ng,latd1)<=aurlat) return

  if (glat_ng(i_ng,latd0)>=aurlat .or. glat_ng(i_ng,latd1)<=-aurlat) then
    call aurora(latd0,latd1)
    return
  endif

  if (glat_ng(i_ng,latd0)>=-aurlat .and. glat_ng(i_ng,latd0)<=aurlat .and. glat_ng(i_ng,latd1)>=aurlat) then
    do latbeg = latd0,latd1
      if (glat_ng(i_ng,latbeg) >= aurlat) exit
    enddo
    call aurora(latbeg,latd1)
    return
  endif

  if (glat_ng(i_ng,latd0)<=-aurlat .and. glat_ng(i_ng,latd1)>=-aurlat .and. glat_ng(i_ng,latd1)<=aurlat) then
    do latend = latd1,latd0,-1
      if (glat_ng(i_ng,latend) <= -aurlat) exit
    enddo
    call aurora(latd0,latend)
    return
  endif

  if (glat_ng(i_ng,latd0)<=-aurlat .and. glat_ng(i_ng,latd1)>=aurlat) then
    do latend = latd1,latd0,-1
      if (glat_ng(i_ng,latend) <= -aurlat) exit
    enddo
    call aurora(latd0,latend)
    do latbeg = latd0,latd1
      if (glat_ng(i_ng,latbeg) >= aurlat) exit
    enddo
    call aurora(latbeg,latd1)
    return
  endif

  contains
!-----------------------------------------------------------------------
  subroutine aurora(latbeg,latend)

    use cons_module,only: pi,p0,grav,gask,rmass_o2,rmass_n2,rmassinv_o2,rmassinv_o1,rmassinv_n2,avo
    use input_module,only: iamie
    use magfield_module,only: sunlons
    use aurora_module,only: fac_p2e,alfad,fd,alfac,fc,add_sproton,alfa_sp, &
      flx_sp,add_helectron,alfa30,e30,theta0,offa,dskofa,phid,rrad,rradp, &
      alfa0,ralfa,ralfa2,rrote,rroth,h0,rh,e0,e20,ree,re2,alfa20
    use fields_ng_module,only: itp,expz

    integer,intent(in) :: latbeg,latend

    logical,parameter :: use_cion = .false.
    real,parameter :: s5 = .08726646, s20 = .34906585, pi_cusp = 3.1415927, s10 = 0.174532925
    integer :: nk,lond0,lond1,k,ihem
    real :: ofda,cosofa,sinofa,aslona,p0ez
    real,dimension(5) :: qia
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend) :: &
      rlatm,rlonm,ekvg,efxg,dlat_aur,dlon_aur,colat,sinlat,coslat,coslon,sinlon,alon, &
      alfa,alfa2,alfa3,flux,flux2,flux3,drizl,cusp,eflux,eflux2,eflux3, &
      coslamda,halfwidth,dtheta,clat,falfa1,falfa2,fcusp,fdrizl,falfa_sp,falfa3
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend) :: &
      tn,o2,o1,n2,mbar,barm,xnmbarm,xalfa1,xalfa2,xcusp,xdrizl,xalfa_sp,xalfa3, &
      flux1_ion,flux2_ion,cusp_ion,drizl_ion,alfa1_ion,alfa2_ion,alfa3_ion,alfasp_bion, &
      barm_t,qsum,denom,p0ez_mbar,tk_mbar,qo2p_aur,qop_aur,qn2p_aur,qaurora, &
      falfa3_alfa3_ion,falfa_sp_alfasp_bion,qo2p_aur_a,qop_aur_a,qn2p_aur_a

    tn = flds(i_ng)%tn(:,:,latbeg:latend,itp(i_ng))
    o2 = flds(i_ng)%o2(:,:,latbeg:latend,itp(i_ng))
    o1 = flds(i_ng)%o1(:,:,latbeg:latend,itp(i_ng))
    n2 = flds(i_ng)%n2(:,:,latbeg:latend)
    mbar = flds(i_ng)%mbar(:,:,latbeg:latend,itp(i_ng))
    barm = flds(i_ng)%barm(:,:,latbeg:latend,itp(i_ng))
    xnmbarm = flds(i_ng)%xnmbarm(:,:,latbeg:latend)

    rlatm = flds(i_ng)%rlatm(:,latbeg:latend)
    rlonm = flds(i_ng)%rlonm(:,latbeg:latend)
    ekvg = flds(i_ng)%ekvg(:,latbeg:latend)
    efxg = flds(i_ng)%efxg(:,latbeg:latend)

    nk = nlevp1_ng(i_ng)
    lond0 = flds(i_ng)%lond0
    lond1 = flds(i_ng)%lond1

    dlat_aur = rlatm
    dlon_aur = rlonm-sunlons(1)
    ihem = int(dlat_aur((lond0+lond1)/2,(latbeg+latend)/2)*2./pi+2.)

    ofda = sqrt(offa(ihem)**2+dskofa(ihem)**2)
    cosofa = cos(ofda)
    sinofa = sin(ofda)
    aslona = asin(dskofa(ihem)/ofda)
    sinlat = sin(abs(dlat_aur))
    coslat = cos(dlat_aur)
    sinlon = sin(dlon_aur+aslona)
    coslon = cos(dlon_aur+aslona)
    colat = acos(cosofa*sinlat-sinofa*coslat*coslon)
    alon = mod(atan2(sinlon*coslat,sinlat*sinofa+cosofa*coslat*coslon)-aslona+pi,pi*2.)

! aurora_cusp
    cusp = (exp(-((theta0(ihem)-colat)/s5)**2)+exp(-((pi_cusp-theta0(ihem)-colat)/s5)**2))* &
      exp(-(atan2(sin(alon-phid(ihem)),cos(alon-phid(ihem)))/s20)**2)

! aurora_heat
    coslamda = cos(atan2(sin(alon-rrote),cos(alon-rrote)))
    halfwidth = h0*(1.-rh*cos(atan2(sin(alon-rroth),cos(alon-rroth))))
    dtheta = colat-rrad(ihem)

    alfa = alfa0*(1.-ralfa*coslamda)
    eflux = e0*(1.-ree*coslamda)*exp(-(dtheta/halfwidth)**2)
    flux = eflux/(alfa*fac_p2e)

    drizl = exp(-((dtheta+abs(dtheta))/(2.*h0))**2)

    alfa2 = alfa20*(1.-ralfa2*coslamda)
    eflux2 = e20*(1.-re2*coslamda)*exp(-(dtheta/halfwidth)**2)
    flux2 = eflux2/(alfa2*fac_p2e)

    alfa3 = alfa30
    flux3 = e30*exp(-(dtheta/halfwidth)**2)/1.602e-6
    eflux3 = flux3*fac_p2e*alfa3

    qteaur(:,latbeg:latend) = -7.e+8*exp(-(dtheta/halfwidth)**2)

    if (iamie > 0) then
      alfa = max(ekvg,1.)/2.
      alfa2 = alfa20
      clat = acos(sin(abs(dlat_aur)))
      flux = max(efxg/(alfa*fac_p2e),1.e-20)
      drizl = exp(-((clat-rradp(ihem)+abs(clat-rradp(ihem)))/s10)**2)
      flux2 = e20*(1.-re2*coslamda)*exp(-((clat-rrad(ihem))/halfwidth)**2)/(alfa2*fac_p2e)
    endif

! aurora_ions
    qia = 0.

    do k = 1,nk
      p0ez = (p0*expz(i_ng,k)/(grav*4.e-6))**0.606
      xalfa1(k,:,:) = p0ez/alfa
      xalfa2(k,:,:) = p0ez/alfa2
      xalfa3(k,:,:) = p0ez/alfa3
      xcusp(k,:,:) = p0ez/alfac
      xdrizl(k,:,:) = p0ez/alfad
    enddo

    p0ez_mbar = xnmbarm/avo
    tk_mbar = gask*tn/(grav*barm)*p0ez_mbar
    xalfa_sp = (tk_mbar/0.00271)**0.58140/alfa_sp

    flux1_ion = 1.e-20
    flux2_ion = 1.e-20
    alfa1_ion = 1.e-20
    alfa2_ion = 1.e-20
    cusp_ion = 1.e-20
    drizl_ion = 1.e-20

    if (.not. use_cion) then
      call aion(xalfa1,alfa1_ion,latbeg,latend)
      call aion(xalfa2,alfa2_ion,latbeg,latend)
      call aion(xcusp,cusp_ion,latbeg,latend)
      call aion(xdrizl,drizl_ion,latbeg,latend)
    else
      call cion(xalfa1,alfa1_ion,alfa,latbeg,latend)
      call cion(xalfa2,alfa2_ion,alfa,latbeg,latend)
      call cion(xcusp,cusp_ion,alfa,latbeg,latend)
      call cion(xdrizl,drizl_ion,alfa,latbeg,latend)
    endif
    if (flx_sp > 1.e-19) call bion(xalfa_sp,alfasp_bion,latbeg,latend)
    if (e30 > 1.e-19) call aion(xalfa3,alfa3_ion,latbeg,latend)

    falfa1 = alfa*flux
    falfa2 = alfa2*flux2
    fcusp = cusp*alfac*fc
    fdrizl = drizl*alfad*fd
    falfa3 = alfa3*flux3

    barm_t = grav*mbar/(35.e-3*gask*tn)
    do k = 1,nk
      falfa_sp = drizl*p0ez_mbar(k,:,:)*flx_sp*1.e6/(tk_mbar(k,:,:)*35.)
      qsum(k,:,:) = falfa1*alfa1_ion(k,:,:)+falfa2*alfa2_ion(k,:,:)+fcusp*cusp_ion(k,:,:)+fdrizl*drizl_ion(k,:,:)
      qaurora(k,:,:) = barm_t(k,:,:)*(falfa1*alfa1_ion(k,:,:)+falfa2*alfa2_ion(k,:,:))
      falfa3_alfa3_ion(k,:,:) = falfa3*alfa3_ion(k,:,:)
      falfa_sp_alfasp_bion(k,:,:) = falfa_sp*alfasp_bion(k,:,:)
    enddo
    if (add_helectron) qsum = qsum+falfa3_alfa3_ion
    if (add_sproton) qsum = qsum+falfa_sp_alfasp_bion
    qsum = qsum*barm_t

    denom = 0.92*n2*rmassinv_n2+1.5*o2*rmassinv_o2+0.56*o1*rmassinv_o1
    qo2p_aur = qsum*o2/(rmass_o2*denom)+qia(2)
    qop_aur = qsum*(0.5*o2*rmassinv_o2+0.56*o1*rmassinv_o1)/denom+qia(3)
    qn2p_aur = qsum*0.7*n2/(rmass_n2*denom)+qia(1)

    do k = 2,nk-1
      qo2p_aur_a(k,:,:) = sqrt(qo2p_aur(k,:,:)*qo2p_aur(k-1,:,:))
      qop_aur_a(k,:,:) = sqrt(qop_aur(k,:,:)*qop_aur(k-1,:,:))
      qn2p_aur_a(k,:,:) = sqrt(qn2p_aur(k,:,:)*qn2p_aur(k-1,:,:))
    enddo

    qo2p_aur_a(1,:,:) = 1.5*qo2p_aur(1,:,:)-0.5*qo2p_aur(2,:,:)
    qop_aur_a(1,:,:) = 1.5*qop_aur(1,:,:)-0.5*qop_aur(2,:,:)
    qn2p_aur_a(1,:,:) = 1.5*qn2p_aur(1,:,:)-0.5*qn2p_aur(2,:,:)

    qo2p_aur_a(nk,:,:) = 1.5*qo2p_aur(nk-1,:,:)-0.5*qo2p_aur(nk-2,:,:)
    qop_aur_a(nk,:,:) = 1.5*qop_aur(nk-1,:,:)-0.5*qop_aur(nk-2,:,:)
    qn2p_aur_a(nk,:,:) = 1.5*qn2p_aur(nk-1,:,:)-0.5*qn2p_aur(nk-2,:,:)

    qo2p(:,:,latbeg:latend) = qo2p(:,:,latbeg:latend)+qo2p_aur_a
    qop(:,:,latbeg:latend) = qop(:,:,latbeg:latend)+qop_aur_a
    qn2p(:,:,latbeg:latend) = qn2p(:,:,latbeg:latend)+qn2p_aur_a
    qnp(:,:,latbeg:latend) = qnp(:,:,latbeg:latend)+.22/.7*qn2p_aur_a
    qtef(:,:,latbeg:latend) = qtef(:,:,latbeg:latend)+1.57*qn2p_aur_a

  end subroutine aurora
!-----------------------------------------------------------------------
  subroutine aion(si,so,latbeg,latend)

    integer,intent(in) :: latbeg,latend
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(in) :: si
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(out) :: so

    real,dimension(8),parameter :: cc = &
      (/3.2333134511131,2.5658873458085,2.2540957232641,0.72971983372673, &
        1.1069072431948,1.7134937681128,1.8835442312993,0.86472135072090/)
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend) :: xlog

    xlog = log(si)
    so = cc(1)*exp(cc(2)*xlog-cc(3)*exp(cc(4)*xlog))+cc(5)*exp(cc(6)*xlog-cc(7)*exp(cc(8)*xlog))

  end subroutine aion
!-----------------------------------------------------------------------
  subroutine bion(si,so,latbeg,latend)

    integer,intent(in) :: latbeg,latend
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(in) :: si
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(out) :: so

    real,dimension(8),parameter :: cc = (/0.12718,4.9119,1.8429,0.99336,0.52472,1.5565,0.85732,1.4116/)
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend) :: xlog

    xlog = log(si)
    so = cc(1)*exp(cc(2)*xlog-cc(3)*exp(cc(4)*xlog))+cc(5)*exp(cc(6)*xlog-cc(7)*exp(cc(8)*xlog))

  end subroutine bion
!-----------------------------------------------------------------------
  subroutine cion(si,so,alpha,latbeg,latend)

    integer,intent(in) :: latbeg,latend
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(in) :: si
    real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(out) :: so
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend),intent(in) :: alpha

    real,dimension(4,8),parameter :: P = reshape(       &
      (/3.49979e-1,-6.18200e-2,-4.08124e-2, 1.65414e-2, &
        5.85425e-1,-5.00793e-2, 5.69309e-2,-4.02491e-3, &
        1.69692e-1,-2.58981e-2, 1.96822e-2, 1.20505e-3, &
       -1.22271e-1,-1.15532e-2, 5.37951e-6, 1.20189e-3, &
        1.57018,    2.87896e-1,-4.14857e-1, 5.18158e-2, &
        8.83195e-1, 4.31402e-2,-8.33599e-2, 1.02515e-2, &
        1.90953,   -4.74704e-2,-1.80200e-1, 2.46652e-2, &
       -1.29566,   -2.10952e-1, 2.73106e-1,-2.92752e-2/),(/4,8/))
    integer :: i,k
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend) :: logE,logY
    real,dimension(flds(i_ng)%lond0:flds(i_ng)%lond1,latbeg:latend,8) :: cc

    logE = log(alpha)
    do i = 1,8
      cc(:,:,i) = exp(P(1,i)+P(2,i)*logE+P(3,i)*logE**2.+P(4,i)*logE**3.)
    enddo

    do k = 1,nlevp1_ng(i_ng)
      where (alpha >= 0.1)
        logY = log(si(k,:,:))
        so(k,:,:) = cc(:,:,1)*exp(cc(:,:,2)*logY-cc(:,:,3)*exp(cc(:,:,4)*logY))+ &
          cc(:,:,5)*exp(cc(:,:,6)*logY-cc(:,:,7)*exp(cc(:,:,8)*logY))
      elsewhere
        so(k,:,:) = 0.
      endwhere
    enddo

  end subroutine cion
!-----------------------------------------------------------------------
end subroutine aurora_ng

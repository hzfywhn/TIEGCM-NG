subroutine qrj_ng(rj,qtef,qtotal,qop2p,qop2d,qo2p,qop,qn2p,qnp,qnop,i_ng)

  use params_module,only: nlevp1_ng
  use input_module,only: f107
  use init_module,only: sfeps
  use cons_module,only: avo,rmassinv_n4s,rmassinv_no,rmassinv_o2,rmassinv_o1,rmassinv_n2
  use qrj_module,only: lmax,l1,sigeuv,rlmeuv,feuv,fsrc,sigsrc,rlmsrc,sigin4s,quench,BPhotonI,BElectronI, &
    brop2pPh,brop2dPh,brop4sPh,bro2DPh,brn2DPh,bro2DIPh,brn2DIPh,brop2pEl,brop2dEl,brop4sEl,bro2DIEl,brn2DIEl,brn2DEl,bro2DEl
  use lbc,only: fb,b
  use fields_ng_module,only: flds,itp
  implicit none

  integer,intent(in) :: i_ng
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1),intent(out) :: &
    rj,qtef,qtotal,qop2p,qop2d,qo2p,qop,qn2p,qnp,qnop

  real,parameter :: do2 = 8.203E-12, do22 = 1.1407E-11, aband = 0.143, bband = 9.64E8, cband = 9.03E-19, &
    e3 = 0.33, hc = 1.9845E-16, euveff = 0.05
  integer :: nk,k,l
  real,dimension(l1) :: rlmsrcinv
  real,dimension(lmax) :: rlmeuvinv
  real,dimension(nlevp1_ng(i_ng),flds(i_ng)%lond0:flds(i_ng)%lond1,flds(i_ng)%latd0:flds(i_ng)%latd1) :: &
    sco2,sco1,scn2,tn,no,o2,o1,he,n2,xnmbari,beta9,o2i,o1i,hei,n2i,n4si,quenchfac,sigchap,p3f, &
    absorp_o,absorp_o2,absorp_n2,ioniz_o,ioniz_o2,ioniz_n2,di_o2,di_n2,mn_o2,mn_o1,mn_n2,mn_n,sum1,sum2,noi

  sco2 = flds(i_ng)%sco2
  sco1 = flds(i_ng)%sco1
  scn2 = flds(i_ng)%scn2
  tn = flds(i_ng)%tn(:,:,:,itp(i_ng))
  no = flds(i_ng)%no(:,:,:,itp(i_ng))
  o2 = flds(i_ng)%o2(:,:,:,itp(i_ng))
  o1 = flds(i_ng)%o1(:,:,:,itp(i_ng))
  he = flds(i_ng)%he(:,:,:,itp(i_ng))
  n2 = flds(i_ng)%n2
  xnmbari = flds(i_ng)%xnmbari(:,:,:,itp(i_ng))

  beta9 = flds(i_ng)%beta9

  nk = nlevp1_ng(i_ng)

  rlmeuvinv = 1./rlmeuv
  rlmsrcinv = 1./rlmsrc

  do k = 1,nk-1
    o2i(k+1,:,:) = 0.5*(o2(k,:,:)+o2(k+1,:,:))
    o1i(k+1,:,:) = 0.5*(o1(k,:,:)+o1(k+1,:,:))
    hei(k+1,:,:) = 0.5*(he(k,:,:)+he(k+1,:,:))
    n2i(k+1,:,:) = 0.5*(n2(k,:,:)+n2(k+1,:,:))
  enddo
  n4si = 0.

  o2i(1,:,:) = .5*((b(1,1)+1.)*o2(1,:,:)+b(1,2)*o1(1,:,:)+b(1,3)*he(1,:,:)+fb(1))
  o1i(1,:,:) = .5*(b(2,1)*o2(1,:,:)+(b(2,2)+1.)*o1(1,:,:)+b(2,3)*he(1,:,:)+fb(2))
  hei(1,:,:) = .5*(b(3,1)*o2(1,:,:)+b(3,2)*o1(1,:,:)+(b(3,3)+1.)*he(1,:,:)+fb(3))

  n2i(1,:,:) = 1.-o2i(1,:,:)-o1i(1,:,:)-hei(1,:,:)

  mn_o2 = o2i*rmassinv_o2
  mn_o1 = o1i*rmassinv_o1
  mn_n2 = n2i*rmassinv_n2
  mn_n = n4si*rmassinv_n4s

  rj = 0.
  qtef = 0.
  qtotal = 0.
  qop2p = 0.
  qop2d = 0.
  qo2p = 0.
  qop = 0.
  qn2p = 0.
  qnp = 0.
  qnop = 0.
  di_o2 = 0.
  di_n2 = 0.

  do l = l1+1,lmax
    sum1 = feuv(l)*exp(-(sigeuv(1,l)*sco2+sigeuv(2,l)*sco1+sigeuv(3,l)*scn2))
    sum2 = sigeuv(1,l)*mn_o2+sigeuv(2,l)*mn_o1+sigeuv(3,l)*mn_n2

    absorp_o2 = sum1*sigeuv(1,l)
    absorp_o = sum1*sigeuv(2,l)
    absorp_n2 = sum1*sigeuv(3,l)
    ioniz_o2 = absorp_o2*BPhotonI(1,l)
    ioniz_o = absorp_o*BPhotonI(2,l)
    ioniz_n2 = absorp_n2*BPhotonI(3,l)

    di_o2 = di_o2+absorp_o2*bro2DIPh(l)+ioniz_o2*bro2DIEl(l)
    di_n2 = di_n2+absorp_n2*brn2DIPh(l)+ioniz_n2*brn2DIEl(l)
    qnp = qnp+sigin4s(l)*sum1
    qn2p = qn2p+ioniz_n2+ioniz_n2*BElectronI(3,l)
    qo2p = qo2p+ioniz_o2+ioniz_o2*BElectronI(1,l)
    qop2p = qop2p+absorp_o*brop2pPh(l)+ioniz_o*brop2pEl(l)
    qop2d = qop2d+absorp_o*brop2dPh(l)+ioniz_o*brop2dEl(l)
    qop = qop+absorp_o*brop4sPh(l)+ioniz_o*brop4sEl(l)

    rj = rj+absorp_o2*bro2DPh(l)+ioniz_o2*bro2DEl(l)
    qtef = qtef+absorp_n2*brn2DPh(l)+ioniz_n2*brn2DEl(l)
    qtotal = qtotal+hc*rlmeuvinv(l)*sum1*sum2
  enddo

  qtotal = qtotal*euveff*avo

  mn_o2 = mn_o2*xnmbari
  mn_n2 = mn_n2*xnmbari
  mn_o1 = mn_o1*xnmbari
  mn_n = mn_n*xnmbari
  di_o2 = di_o2*mn_o2
  di_n2 = di_n2*mn_n2
  qo2p = qo2p*mn_o2-di_o2
  qn2p = qn2p*mn_n2-di_n2
  qnp = qnp*mn_n+di_n2
  qop = qop*mn_o1+0.54*di_o2
  qop2p = qop2p*mn_o1+0.22*di_o2
  qop2d = qop2d*mn_o1+0.24*di_o2
  qtef = 2.*qtef*mn_n2+di_n2

  do k = 2,nk
    noi(k,:,:) = .5*(no(k,:,:)+no(k-1,:,:))
  enddo
  noi(1,:,:) = no(1,:,:)

  qnop = qnop+beta9*noi*xnmbari*rmassinv_no

  quenchfac = xnmbari*(quench(1)*n2i*rmassinv_n2+quench(2)*o2i*rmassinv_o2)
  quenchfac = quench(3)*quenchfac/(quench(4)+quenchfac)

  sum1 = 0.
  do l = 1,l1
    sigchap = sigsrc(l)*fsrc(l)*exp(-sigsrc(l)*sco2)
    sum1 = sum1+sigchap*(hc*rlmsrcinv(l)-do22+quenchfac)
    rj = rj+sigchap
  enddo

  qtotal = qtotal+sum1*avo*o2i*rmassinv_o2

  where (sco2 >= 1.e18)
    p3f = 1./(aband*sco2+bband*sqrt(sco2))*(1.+0.11*(f107-65.)/165.)*sfeps
  elsewhere
    p3f = cband*(1.+0.11*(f107-65.)/165.)*sfeps
  endwhere
  qtotal = qtotal+p3f*avo*o2i*rmassinv_o2*e3
  rj = rj+p3f/do2

end subroutine qrj_ng

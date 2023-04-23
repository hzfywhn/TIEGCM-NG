module fields_ng_module

  use params_module,only: mx_ng
  use fields_module,only: nf4d,fields_4d,fields_3d,fields_2d
  implicit none

  integer,parameter :: nf3din = 3, nf2din = 7

  type fields
    integer :: lon0,lon1,lat0,lat1,lond0,lond1,latd0,latd1
    logical,dimension(4) :: is_bndry

! const fields per process
    real,dimension(:),allocatable :: racs,cs,tanphi,cor
    real,dimension(:,:),allocatable :: &
      alatm,alonm,xb,yb,zb,bmod,bx,by,bz,bmod2,rlatm,rlonm,dipmag,decmag,sndec,csdec,sn2dec
    real,dimension(:,:,:,:),allocatable :: rjac

! essential 2d/3d/4d fields, use pointer to associate with f2d/f3d/f4d
    real,dimension(:,:),allocatable :: tlbc,ulbc,vlbc,tlbc_nm,ulbc_nm,vlbc_nm
    real,dimension(:,:),pointer :: t_lbc,u_lbc,v_lbc,z_lbc,flx_he,ekvg,efxg
    real,dimension(:,:,:),pointer :: ex,ey,ez
    real,dimension(:,:,:,:),pointer :: &
      tn,    un,    vn,    o2,    o1,    n4s,   no,    he,    ar,    op,    n2d,   &
      ti,    te,    ne,    o2p,   w,     z,     poten, tn_nm, un_nm, vn_nm, o2_nm, &
      o1_nm, n4s_nm,no_nm, op_nm, he_nm, ar_nm, barm,  vc,    mbar
    type(fields_2d),dimension(nf2din) :: f2din
    type(fields_3d),dimension(nf3din) :: f3din
    type(fields_4d),dimension(nf4d) :: f4d

! essential 2d/3d fields at previous and current times
    real,dimension(:,:,:,:),allocatable :: f2d_save
    real,dimension(:,:,:,:,:),allocatable :: f3d_save

! boundaries (essential)
    real,dimension(:,:,:,:,:),allocatable :: lon_b,lat_b
    real,dimension(:,:,:,:),allocatable :: op_lon_b,op_lat_b

! auxiliary fields
    real,dimension(:,:),allocatable :: qteaur,chi
    real,dimension(:,:,:),allocatable :: &
      kldt,   kldu,   kldv,   kldo2,  kldo1,  kldhe,  xnmbar, xnmbari,xnmbarm, &
      cp,     kt,     km,     ui,     vi,     wi,     vo2,    vo1,    vn2,     &
      sco2,   sco1,   scn2,   xiop2p, xiop2d, nplus,  n2p,    nop,    lxx,     &
      lyy,    lxy,    lyx,    qji_ti, qji_tn, cool_implicit,  cool_explicit,   &
      hdt,    hdu,    hdv,    hdo2,   hdo1,   hdhe,   ped,    hall,   lam1,    &
      zg,     n2,     wn,     Fe,     Fn,     fnrh,   fkmh,   rk1,    rk2,     &
      rk3,    ra1,    ra2,    ra3,    beta1,  beta3,  beta5,  beta8,  beta9,   &
      beta9n, beta17, rk19,   rk20,   rk25,   rkm12,  rj,     qtef,   qtotal,  &
      qop2p,  qop2d,  qo2p,   qop,    qn2p,   qnp,    qnop
    real,dimension(:,:,:,:,:),allocatable :: fs
  end type fields

  type(fields),dimension(mx_ng) :: flds

! level unrelated fields
  real,parameter :: hor = .25
  real,parameter :: t0 = 0. ! t0 is not currently used, not sure whether it will be used in the future
  real,dimension(3) :: fb
  real,dimension(3,3) :: b ! fb,b are constant along lon, simplify to lon independent constants

! common fields across processes
  integer,dimension(mx_ng) :: itp,itc,maxlon,maxlat
  integer,dimension(:,:,:),allocatable :: domain
  real,dimension(mx_ng) :: shapiro,dtx2inv,dlev,dphi,dlamda,dz,dzp,expzmid,expzmid_inv,dipmin,modeltime,step
  real,dimension(:,:),allocatable :: expz
  real,dimension(:,:,:),allocatable :: difk,dift,xmue

  character(len=*),dimension(*),parameter :: fmap = &
    (/'TN    ','UN    ','VN    ','O2    ','O1    ','N4S   ', &
      'NO    ','HE    ','AR    ','OP    ','N2D   ','TI    ', &
      'TE    ','NE    ','O2P   ','OMEGA ','Z     ','TN_NM ', &
      'UN_NM ','VN_NM ','O2_NM ','O1_NM ','N4S_NM','NO_NM ', &
      'OP_NM ','HE_NM ','AR_NM '/)
  character(len=*),dimension(*),parameter :: bndry = &
    (/'TN   ','UN   ','VN   ','O2   ','O1   ','N4S  ', &
      'NO   ','HE   ','AR   ','OP   ','OMEGA'/)
  integer,parameter :: nmap = size(fmap), nbnd = size(bndry)

  character(len=*),dimension(*),parameter :: ubfill = &
    (/'TN   ','UN   ','VN   ','OP   ','N2D  ','TI   ', &
      'TE   ','O2P  ','TN_NM','UN_NM','VN_NM','OP_NM'/)
  character(len=*),dimension(*),parameter :: zlog = &
    (/'TN    ','O2    ','O1    ','N4S   ','NO    ','HE    ', &
      'AR    ','OP    ','N2D   ','TI    ','TE    ','NE    ', &
      'O2P   ','Z     ','TN_NM ','O2_NM ','O1_NM ','N4S_NM', &
      'NO_NM ','OP_NM ','HE_NM ','AR_NM '/)

  contains
!-----------------------------------------------------------------------
  subroutine init_ng

    use params_module,only: n_ng,nlon_ng,nlat_ng
    use input_module,only: mkntask
    use mpi_module,only: ntask,mytid,distribute_1d
    use mpi_f08,only: mpi_allgather,mpi_integer,mpi_comm_world

    integer :: i_ng,ntaski,ntaskj,ierror

    allocate(domain(n_ng,4,0:ntask-1))

    do i_ng = 1,n_ng
      call mkntask(ntask,nlon_ng(i_ng),nlat_ng(i_ng),ntaski,ntaskj,ierror)
      call distribute_1d(1,nlon_ng(i_ng),ntaski,mod(mytid,ntaski),flds(i_ng)%lon0,flds(i_ng)%lon1)
      call distribute_1d(1,nlat_ng(i_ng),ntaskj,mytid/ntaski,flds(i_ng)%lat0,flds(i_ng)%lat1)

      call mpi_allgather((/flds(i_ng)%lon0,flds(i_ng)%lon1,flds(i_ng)%lat0,flds(i_ng)%lat1/), &
        4,mpi_integer,domain(i_ng,:,:),4,mpi_integer,mpi_comm_world)
      maxlon(i_ng) = maxval(domain(i_ng,2,:)-domain(i_ng,1,:))+5
      maxlat(i_ng) = maxval(domain(i_ng,4,:)-domain(i_ng,3,:))+5

      flds(i_ng)%lond0 = flds(i_ng)%lon0-2
      flds(i_ng)%lond1 = flds(i_ng)%lon1+2
      flds(i_ng)%latd0 = flds(i_ng)%lat0-2
      flds(i_ng)%latd1 = flds(i_ng)%lat1+2

      flds(i_ng)%is_bndry = .false.
      if (flds(i_ng)%lon0 == 1) flds(i_ng)%is_bndry(1) = .true.
      if (flds(i_ng)%lon1 == nlon_ng(i_ng)) flds(i_ng)%is_bndry(2) = .true.
      if (flds(i_ng)%lat0 == 1) flds(i_ng)%is_bndry(3) = .true.
      if (flds(i_ng)%lat1 == nlat_ng(i_ng)) flds(i_ng)%is_bndry(4) = .true.
    enddo

    itp = 1
    itc = 2

  end subroutine init_ng
!-----------------------------------------------------------------------
  subroutine alloc_fields

    use params_module,only: n_ng,nlevp1_ng
    use cons_module,only: ndays
    use input_module,only: nstep_ng,nstep_sub

    integer :: i_ng,nk,lond0,lond1,latd0,latd1

    allocate(expz(n_ng,maxval(nlevp1_ng)))
    allocate(difk(n_ng,maxval(nlevp1_ng),ndays))
    allocate(dift(n_ng,maxval(nlevp1_ng),ndays))
    allocate(xmue(n_ng,maxval(nlevp1_ng),ndays))

    do i_ng = 1,n_ng
      nk = nlevp1_ng(i_ng)
      lond0 = flds(i_ng)%lond0
      lond1 = flds(i_ng)%lond1
      latd0 = flds(i_ng)%latd0
      latd1 = flds(i_ng)%latd1

! constant fields
      allocate(flds(i_ng)%racs(latd0:latd1))
      allocate(flds(i_ng)%cs(latd0:latd1))
      allocate(flds(i_ng)%tanphi(latd0:latd1))
      allocate(flds(i_ng)%cor(latd0:latd1))

      allocate(flds(i_ng)%alatm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%alonm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xb(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%yb(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%zb(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%bmod(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%bx(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%by(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%bz(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%bmod2(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rlatm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rlonm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%dipmag(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%decmag(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%sndec(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%csdec(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%sn2dec(lond0:lond1,latd0:latd1))

      allocate(flds(i_ng)%rjac(lond0:lond1,latd0:latd1,2,2))

! essential fields
      allocate(flds(i_ng)%tlbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ulbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vlbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%tlbc_nm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ulbc_nm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vlbc_nm(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%t_lbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%u_lbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%v_lbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%z_lbc(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%flx_he(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ekvg(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%efxg(lond0:lond1,latd0:latd1))

      flds(i_ng)%f2din(1)%data => flds(i_ng)%t_lbc
      flds(i_ng)%f2din(2)%data => flds(i_ng)%u_lbc
      flds(i_ng)%f2din(3)%data => flds(i_ng)%v_lbc
      flds(i_ng)%f2din(4)%data => flds(i_ng)%z_lbc
      flds(i_ng)%f2din(5)%data => flds(i_ng)%flx_he
      flds(i_ng)%f2din(6)%data => flds(i_ng)%ekvg
      flds(i_ng)%f2din(7)%data => flds(i_ng)%efxg

      allocate(flds(i_ng)%ex(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ey(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ez(nk,lond0:lond1,latd0:latd1))

      flds(i_ng)%f3din(1)%data => flds(i_ng)%ex
      flds(i_ng)%f3din(2)%data => flds(i_ng)%ey
      flds(i_ng)%f3din(3)%data => flds(i_ng)%ez

      allocate(flds(i_ng)%tn(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%un(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%vn(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%o2(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%o1(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%n4s(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%no(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%he(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%ar(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%op(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%n2d(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%ti(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%te(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%ne(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%o2p(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%w(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%z(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%poten(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%tn_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%un_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%vn_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%o2_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%o1_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%n4s_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%no_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%op_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%he_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%ar_nm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%barm(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%vc(nk,lond0:lond1,latd0:latd1,2))
      allocate(flds(i_ng)%mbar(nk,lond0:lond1,latd0:latd1,2))

      flds(i_ng)%f4d(1)%data => flds(i_ng)%tn
      flds(i_ng)%f4d(2)%data => flds(i_ng)%un
      flds(i_ng)%f4d(3)%data => flds(i_ng)%vn
      flds(i_ng)%f4d(4)%data => flds(i_ng)%o2
      flds(i_ng)%f4d(5)%data => flds(i_ng)%o1
      flds(i_ng)%f4d(6)%data => flds(i_ng)%n4s
      flds(i_ng)%f4d(7)%data => flds(i_ng)%no
      flds(i_ng)%f4d(8)%data => flds(i_ng)%he
      flds(i_ng)%f4d(9)%data => flds(i_ng)%ar
      flds(i_ng)%f4d(10)%data => flds(i_ng)%op
      flds(i_ng)%f4d(11)%data => flds(i_ng)%n2d
      flds(i_ng)%f4d(12)%data => flds(i_ng)%ti
      flds(i_ng)%f4d(13)%data => flds(i_ng)%te
      flds(i_ng)%f4d(14)%data => flds(i_ng)%ne
      flds(i_ng)%f4d(15)%data => flds(i_ng)%o2p
      flds(i_ng)%f4d(16)%data => flds(i_ng)%w
      flds(i_ng)%f4d(17)%data => flds(i_ng)%z
      flds(i_ng)%f4d(18)%data => flds(i_ng)%poten
      flds(i_ng)%f4d(19)%data => flds(i_ng)%tn_nm
      flds(i_ng)%f4d(20)%data => flds(i_ng)%un_nm
      flds(i_ng)%f4d(21)%data => flds(i_ng)%vn_nm
      flds(i_ng)%f4d(22)%data => flds(i_ng)%o2_nm
      flds(i_ng)%f4d(23)%data => flds(i_ng)%o1_nm
      flds(i_ng)%f4d(24)%data => flds(i_ng)%n4s_nm
      flds(i_ng)%f4d(25)%data => flds(i_ng)%no_nm
      flds(i_ng)%f4d(26)%data => flds(i_ng)%op_nm
      flds(i_ng)%f4d(27)%data => flds(i_ng)%he_nm
      flds(i_ng)%f4d(28)%data => flds(i_ng)%ar_nm
      flds(i_ng)%f4d(29)%data => flds(i_ng)%barm
      flds(i_ng)%f4d(30)%data => flds(i_ng)%vc
      flds(i_ng)%f4d(31)%data => flds(i_ng)%mbar

      allocate(flds(i_ng)%f2d_save(lond0:lond1,latd0:latd1,0:nstep_ng(i_ng),nf2din))
      allocate(flds(i_ng)%f3d_save(nk,lond0:lond1,latd0:latd1,0:nstep_ng(i_ng),nf3din))

      if (flds(i_ng)%is_bndry(1) .or. flds(i_ng)%is_bndry(2)) then
        allocate(flds(i_ng)%lon_b(nk,1:2,latd0:latd1,0:nstep_ng(i_ng),nbnd))
        allocate(flds(i_ng)%op_lon_b(nk,1:2,latd0:latd1,0:nstep_ng(i_ng)*nstep_sub))
      endif
      if (flds(i_ng)%is_bndry(3) .or. flds(i_ng)%is_bndry(4)) then
        allocate(flds(i_ng)%lat_b(nk,lond0:lond1,3:4,0:nstep_ng(i_ng),nbnd))
        allocate(flds(i_ng)%op_lat_b(nk,lond0:lond1,3:4,0:nstep_ng(i_ng)*nstep_sub))
      endif

! auxiliary fields
      allocate(flds(i_ng)%qteaur(lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%chi(lond0:lond1,latd0:latd1))

      allocate(flds(i_ng)%kldt(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kldu(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kldv(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kldo2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kldo1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kldhe(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xnmbar(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xnmbari(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xnmbarm(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%cp(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%kt(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%km(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ui(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vi(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%wi(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vo2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vo1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%vn2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%sco2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%sco1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%scn2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xiop2p(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%xiop2d(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%nplus(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%n2p(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%nop(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%lxx(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%lyy(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%lxy(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%lyx(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qji_ti(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qji_tn(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%cool_implicit(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%cool_explicit(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdt(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdu(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdv(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdo2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdo1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hdhe(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ped(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%hall(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%lam1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%zg(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%n2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%wn(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%Fe(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%Fn(nk,lond0:lond1,latd0:latd1))

      allocate(flds(i_ng)%fnrh(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%fkmh(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk3(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ra1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ra2(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%ra3(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta1(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta3(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta5(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta8(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta9(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta9n(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%beta17(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk19(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk20(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rk25(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rkm12(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%rj(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qtef(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qtotal(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qop2p(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qop2d(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qo2p(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qop(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qn2p(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qnp(nk,lond0:lond1,latd0:latd1))
      allocate(flds(i_ng)%qnop(nk,lond0:lond1,latd0:latd1))

      allocate(flds(i_ng)%fs(nk,lond0:lond1,latd0:latd1,3,0:3))
    enddo

  end subroutine alloc_fields
!-----------------------------------------------------------------------
end module fields_ng_module

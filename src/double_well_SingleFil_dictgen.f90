program double_well_test
  use mpi
  implicit none
  include 'param.inc'

  integer :: ierr, rank, procs
  integer, parameter :: np = 1000
  integer, parameter :: SwitchZupdate = 0
  !Energy Parameters
  real(8), parameter :: KB_T = 0.0138 * 310                !pN*nm Boltzmann Constact x Temparature
  real(8), parameter :: E_ATP = 22.5*KB_T                  !pN*nm Energy of ATP hydrolysis
  real(8), parameter :: c_pre = 8.0, c_pos = 8.0           !pN/nm Curvature of the two wells
  real(8), parameter :: E_pre = 0.7*E_ATP, E_pos = 0.0     !pN*nm Assuming 70% of E_ATP is used for powerstroke
  real(8), parameter :: x_pre = 0.0, x_pos = 8.5           !nm    Levearm end position
  !org  real(8), parameter :: delta = 4.0*KB_T                   !pN*nm Barrier relaxation
  real(8), parameter :: delta = 2.0*KB_T                   !pN*nm Barrier relaxation
  real(8), parameter :: omega_stiff = 1.0, c_minus = 2.5   !Unitless Stiffness coefficient
  real(8), parameter :: k_spring = 2.8                     !pN/nm Spring
  real(8), parameter :: x_S0 = 0.0                     !nm   Spring energy = 0.5*k_spring*(x_S+x)
  real(8) :: x_barrier
  !Friction
  real(8), parameter :: fric_x = 80.0 !pN*ns/nm  friction for leverarm rotation
  real(8), parameter :: fric_d = 80.0 !pN*ns/nm  friction for x_S during detachment
  real(8), parameter :: x_min_detach = -10.0
  !Transition
  real(8), parameter :: a_trans = 500, d_trans = 5000, g_trans=100  !1/s Transition rate constants
  !Sarcomere (per one active filament(AF))
  real(8) :: z                                             !nm Contraction distance
  real(8) :: FzPerAF                                       !pN Contraction Force per one AF
  integer, parameter :: nM = 20                            !Number of myosins per one AF
  real(8), parameter :: gamma_sarco = 1.d-5                !pN*s/nm Viscosity per one AF
  real(8), parameter :: kZ = 0.5                           !pN/nm   Spring constant per one AF
  !Time step
  !dt = 0.5
!  integer, parameter :: time=2000
!  integer, parameter :: time=500
! integer, parameter :: time=20000
!  integer, parameter :: time=25000
!   integer, parameter :: time=100
!  integer, parameter :: time=50000
!C in param.inc  integer, parameter :: time=100000
  real(8), parameter :: dt = 0.5   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 200000, nt_out = 10000
  !integer, parameter :: nt_in = 20000, nt_out = 500000  !Fine record
  !integer, parameter :: nt_in = 200 
  integer, parameter :: nt_in = int(time/dt+0.1)
!  integer, parameter :: nt_in = 1000
  ! integer, parameter :: nt_in = 2000
  !dt = 1.0
  !  real(8), parameter :: dt = 1.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 100000, nt_out = 10000
  !dt = 2.0
  !  real(8), parameter :: dt = 2.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 50000, nt_out = 10000
  !dt = 10.0
  !  real(8), parameter :: dt =  10.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 10000, nt_out = 10000
  integer :: it_out, it_in
  real(8) :: total_time
  !Random Force
  integer, parameter :: NP_RandomForce = 10000000
  real(8) :: RandomForceArray(NP_RandomForce)
  !Statistical Data
  integer :: npre, npos, ndet
  real(8) :: NattachPerMol, Ndeatch_prePerMol, Ndeatch_posPerMol

  integer :: state(np)
  integer :: attach_count(np), detach_pre_count(np), detach_pos_count(np)
  real(8) :: x_L(np), force(np), force_random(np)
  real(8) :: x_S(np)
  real(8) :: phi(np), dphi(np), d2phi(np), stiff_ps(np)
  integer :: RndForceSeedArray(np)  ! random force seeds
  integer :: RndStateSeedArray(np)  ! state transition seeds

  !C in param.inc integer, parameter :: dict_2d_xlowb = -7
  !C in param.inc integer, parameter :: dict_2d_xupb = 14
  !C in param.inc integer, parameter :: dict_2d_xsz = 100   ! 読み込むファイル数 (適宜変更)      
  !C in param.inc integer, parameter :: dict_2d_slowb = -19
  !C in param.inc integer, parameter :: dict_2d_supb = 23
  !C in param.inc integer, parameter :: dict_2d_shift_sz = 200   ! 読み込むファイル数 (適宜変更)
  !C in param.inc integer, parameter :: dict_1d_xlowb = -6
  !C in param.inc integer, parameter :: dict_1d_xupb = 6
  !C in param.inc integer, parameter :: dict_1d_xsz = 100
  integer, parameter :: dict_psz=1000

  !C in param.inc integer, parameter :: dictp_2d_xlowb = -3
  !C in param.inc integer, parameter :: dictp_2d_xupb = 10
  !integer, parameter :: dictp_2d_xsz = 100   ! 読み込むファイル数 (適宜変更)      
  !C in param.inc integer, parameter :: dictp_2d_xsz = 200   ! 読み込むファイル数 (適宜変更)        
  !C in param.inc integer, parameter :: dictp_2d_slowb = -3
  !C in param.inc integer, parameter :: dictp_2d_supb = 3
  !integer, parameter :: dictp_2d_slowb = -2
  !integer, parameter :: dictp_2d_supb = 2
!  integer, parameter :: dictp_2d_shift_sz = 300   ! 読み込むファイル数 (適宜変更)
  !C in param.inc integer, parameter :: dictp_2d_shift_sz = 30   ! 読み込むファイル数 (適宜変更)

!  real(8) :: dict_2d(dict_psz, dict_2d_xsz, dict_2d_shift_sz)  
!  real(8) :: dict_1d(dict_psz, dict_1d_xsz)
!  real(8) :: dictp_2d(dict_psz, dictp_2d_xsz, dictp_2d_shift_sz)

  real(8), allocatable :: dict_2d(:, :, :)  
  real(8), allocatable :: dict_1d(:, :)
  real(8), allocatable :: dictp_2d(:, :, :)

  real(8),allocatable :: local_2d(:,:,:)
  real(8),allocatable :: local_1d(:,:)
  integer, allocatable :: sendcounts(:), displs(:)
  integer :: sz, local_sz, j_begin, j_end
  integer :: lowerdimension_sz

  character(len=150) :: filename
  real(8) :: dx,ds
  integer :: in_x, in_s,i

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, procs, ierr)
  write(*,*) "rank=",rank

  if (rank == 0) then
    allocate(dict_2d(dict_psz, dict_2d_xsz, dict_2d_shift_sz) )
    allocate(dict_1d(dict_psz, dict_1d_xsz) )
    allocate(dictp_2d(dict_psz, dictp_2d_xsz, dictp_2d_shift_sz) )
  end if

  write(*,*) "dt, dt in dict, num steps:",time, dt, nt_in
  call draw_powerstroke_potential()
  call initRandomForce()

  sz = dict_2d_shift_sz
  lowerdimension_sz = (dict_psz*dict_2d_xsz)

  call set_begin_end(j_begin, j_end, local_sz, procs, sz, rank)

  ! ローカル配列（nx × local_ny）を各ランクで確保
  allocate(local_2d(dict_psz, dict_2d_xsz, local_sz))
  allocate(sendcounts(procs), displs(procs))

  ! ルートが全体配列を用意して初期化（例：A(i,j)=i+j）
  if (rank == 0) then
     ! Scatterv/Gatherv 用の sendcounts, displs（要素数・要素オフセット）
     call set_counts_dspls(sendcounts, displs, procs, sz, lowerdimension_sz)
  end if

  dx=(dict_2d_xupb-dict_2d_xlowb)/real((dict_2d_xsz-1),kind=8)
  ds=(dict_2d_supb-dict_2d_slowb)/real((dict_2d_shift_sz-1),kind=8)  
  do in_s = j_begin-1, j_end-1
    do in_x = 0, dict_2d_xsz-1
      call initialize(in_x,in_s,dict_2d_xsz,0)
      state(:)=1
      x_L(:)      = dict_2d_xlowb+in_x*dx
      x_S(:)= dict_2d_slowb+in_s*ds
      do it_in = 1, nt_in
        call simu_force_random()
        call simu_force_powerstroke()
        call simu_update()
!       call simu_state_trans()
      end do
  
      !$omp parallel do default(none) &
      !$omp private(i) &
      !$omp shared(local_2d,in_x,in_s,x_L,j_begin) 
      do i=1, np
        local_2d(i,in_x+1,in_s+1-j_begin+1)=x_L(i)
      end do
    end do
    write(*,*) "dict_2d:",in_s+1," in ",dict_2d_shift_sz    
  end do

  call MPI_Gatherv(local_2d, lowerdimension_sz*local_sz, MPI_DOUBLE_PRECISION, &
                    dict_2d, sendcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  deallocate(local_2d,sendcounts,displs)



  sz = dict_1d_xsz
  lowerdimension_sz = (dict_psz)

  call set_begin_end(j_begin, j_end, local_sz, procs, sz, rank)

  ! ローカル配列（nx × local_ny）を各ランクで確保
  allocate(local_1d(dict_psz, local_sz))
  allocate(sendcounts(procs), displs(procs))

  ! ルートが全体配列を用意して初期化（例：A(i,j)=i+j）
  if (rank == 0) then
     ! Scatterv/Gatherv 用の sendcounts, displs（要素数・要素オフセット）
     call set_counts_dspls(sendcounts, displs, procs, sz, lowerdimension_sz)
  end if

  dx = (dict_1d_xupb-dict_1d_xlowb)/real((dict_1d_xsz-1),kind=8)  
  state(:)=0
  do in_x = j_begin-1, j_end-1
    call initialize(in_x,dict_2d_shift_sz,dict_2d_xsz,0)
    x_S(:)= dict_1d_xlowb+in_x*dx
    do it_in = 1, nt_in
      call simu_force_random()
      call simu_force_powerstroke()
      call simu_update()
!       call simu_state_trans()
    end do

    !$omp parallel do default(none) &
    !$omp private(i) &
    !$omp shared(local_1d,in_x,x_S,j_begin) 
    do i=1, np
      local_1d(i,in_x+1-j_begin+1)=x_S(i)
    end do
    if(mod(in_x,10)==0) then
      write(*,*) "dict_1d:",in_x+1,"in",dict_1d_xsz
    end if
  end do

  call MPI_Gatherv(local_1d, lowerdimension_sz*local_sz, MPI_DOUBLE_PRECISION, &
                    dict_1d, sendcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  deallocate(local_1d,sendcounts,displs)


  dx=(dictp_2d_xupb-dictp_2d_xlowb)/real((dictp_2d_xsz-1),kind=8)
  ds=(dictp_2d_supb-dictp_2d_slowb)/real((dictp_2d_shift_sz-1),kind=8)  

  sz = dictp_2d_shift_sz !C ??bug??
  lowerdimension_sz = (dict_psz*dictp_2d_xsz)
  call set_begin_end(j_begin, j_end, local_sz, procs, sz, rank)

  ! ローカル配列（nx × local_ny）を各ランクで確保
  allocate(local_2d(dict_psz, dictp_2d_xsz, local_sz))
  allocate(sendcounts(procs), displs(procs))

  ! ルートが全体配列を用意して初期化（例：A(i,j)=i+j）
  if (rank == 0) then
     ! Scatterv/Gatherv 用の sendcounts, displs（要素数・要素オフセット）
     call set_counts_dspls(sendcounts, displs, procs, sz, lowerdimension_sz)
  end if

  do in_s = j_begin-1, j_end-1
    do in_x = 0, dictp_2d_xsz-1
      call initialize(in_x,in_s,dictp_2d_xsz,2*np*(dict_1d_xsz+dict_2d_xsz*dict_2d_shift_sz))
      state(:)=1
      x_L(:)      = dictp_2d_xlowb+in_x*dx
      x_S(:)= dictp_2d_slowb+in_s*ds
      do it_in = 1, nt_in
        call simu_force_random()
        call simu_force_powerstroke()
        call simu_update()
!       call simu_state_trans()
      end do

    !$omp parallel do default(none) &
    !$omp private(i) &
    !$omp shared(local_2d,in_x,in_s,x_L,j_begin) 
      do i=1, np
        local_2d(i,in_x+1,in_s+1-j_begin+1)=x_L(i)
      end do
    end do 
    write(*,*) "dictp_2d:",in_s+1," in ",dictp_2d_shift_sz    
  end do

  call MPI_Gatherv(local_2d, lowerdimension_sz*local_sz, MPI_DOUBLE_PRECISION, &
                  dictp_2d, sendcounts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
  deallocate(local_2d,sendcounts,displs)

  if(rank== 0) then
    write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',&
      & time,'_',dict_2d_xlowb, '_',dict_2d_xupb,'_',dict_2d_xsz,'_',dict_2d_slowb,'_',dict_2d_supb,&
      & '_',dict_2d_shift_sz,'.csv'  
    write(*,*) filename

    open(20,file=filename, form="unformatted")
    write(20) dict_2d(:,:,:)
    close(20)

    write(filename, '(A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',time, &
      & '_',dict_1d_xlowb,'_',dict_1d_xupb,'_',dict_1d_xsz,'.csv'  
    write(*,*) filename

    open(20,file=filename, form="unformatted")
    write(20) dict_1d(:,:)
    close(20)
  
    write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',time,'_',dictp_2d_xlowb,&
      & '_',dictp_2d_xupb,'_',dictp_2d_xsz,'_',dictp_2d_slowb,'_',dictp_2d_supb,'_', &
      & dictp_2d_shift_sz,'.csv'  
    write(*,*) filename
    open(20,file=filename, form="unformatted")
    write(20) dictp_2d(:,:,:)
    close(20)

    deallocate(dict_2d, dict_1d, dictp_2d)
  end if 

  call MPI_Finalize(ierr)
  
contains

  subroutine set_begin_end(j_begin, j_end, local_sz, procs, sz, rank)
    implicit none
    integer, intent(out) :: j_begin, j_end, local_sz
    integer, intent(in) :: procs, sz, rank
    integer :: extra
    local_sz = sz / procs
    extra    = mod(sz, procs)
    if (rank < extra) then
      local_sz = local_sz + 1
      j_begin  = rank * local_sz + 1
    else
      j_begin  = extra*(local_sz+1) + (rank-extra)*local_sz + 1
    end if
    j_end = j_begin + local_sz - 1
  end subroutine set_begin_end

  subroutine set_counts_dspls(sendcounts, displs, procs, sz, lowerdimension_sz)
    implicit none
    integer :: i, extra
    integer :: sendcounts(:), displs(:)
    integer, intent(in) :: lowerdimension_sz, procs, sz

    extra = mod(sz, procs)
    do i = 0, procs-1
        integer :: locsz_i, jb_i
        locsz_i = sz / procs
        if (i < extra) locsz_i = locsz_i + 1
        if (i < extra) then
           jb_i = i*locsz_i + 1
        else
           jb_i = extra*( (sz/procs)+1 ) + (i-extra)*(sz/procs) + 1
        end if
        sendcounts(i+1) = lowerdimension_sz * locsz_i
        displs(i+1)     = lowerdimension_sz * (jb_i - 1)       ! 先頭からの要素オフセット
     end do
  end subroutine set_counts_dspls


  subroutine draw_powerstroke_potential()
    implicit none
    integer, parameter :: nx = 1000
    real(8), parameter :: x_min = -5.0, x_max = 15.0
    integer :: i, ic
    real(8) :: omega, x, phi, dphi, d2phi
    real(8) :: phi_save(0:nx)

!    open(10,file="PowerStrokePot.csv")
!    write(10,*) "#x, phi, dphi/dx, d^2phi/dx^2"
    do i = 0, nx
      omega = dble(i)/dble(nx)
      x = (1.d0-omega)*x_min + omega*x_max
      call powerstroke_potential(x, phi, dphi, d2phi)
!      write(10,'(4f12.4)') x, phi, dphi, d2phi
      phi_save(i) = phi
    end do
!    close(10)

    ic = 0
    do i = 1, nx-1
      omega = dble(i)/dble(nx)
      x = (1.d0-omega)*x_min + omega*x_max
      if (phi_save(i) >= max(phi_save(i-1),phi_save(i+1))) then
        x_barrier = x
        write(6,*) "x_barrier=", x_barrier
        ic = ic + 1
      end if
    end do
    if (ic .ne. 1) then
      write(6,*) "Error in draw_powerstroke_potential ic=", ic
      stop
    end if
  end subroutine draw_powerstroke_potential

  subroutine initialize(in_x,in_s,xsz,shift)
    implicit none
    integer :: i
    integer, intent(in) :: in_x, in_s, xsz, shift
    z = 0.d0
    x_L = 0.d0
    x_S = x_S0
    state = 0
    attach_count = 0
    detach_pre_count = 0
    detach_pos_count = 0
    do i = 1, np
      RndForceSeedArray(i) = i     + 2*np*(in_x+xsz*in_s)+shift
      RndStateSeedArray(i) = np+i  + 2*np*(in_x+xsz*in_s)+shift
    end do
  end subroutine initialize

  subroutine initRandomForce()
    implicit none
    integer, parameter :: NX = 10000
    real(8), parameter :: PI = 3.14159265358979323846264338 !< Circular constant (Pi)
    real(8) :: sigma
    real(8) :: w(NX)

    real(8) :: x, X_MIN, X_MAX, dx, wSum, rnd, CF, w0, w1, EPS
    integer :: ix, ip, iw0, iw1
    real(8), parameter :: Fsigma = 5.0
    integer :: i0, nmod0, nmod1
    integer :: i, j, k, ic


    EPS = 0.5d0/NP_RandomForce
    sigma = 1.d0
    CF = 1.d0 / (sqrt(2.d0 * pi) * sigma)
    X_MIN = -Fsigma * sigma
    X_MAX =  Fsigma * sigma
    dx = (X_MAX - X_MIN) / dble(NX)
    wSum = 0.d0
    do ix = 1, NX
      x = X_MIN + dx * (dble(ix - 1) + 0.5d0)
      w(ix) = CF * exp(-0.5d0 * (x / sigma) ** 2)
      wSum = wSum + w(ix)
    end do
    w(1:NX) = w(1:NX) / wSum

    w0 = 0.d0
    iw0 = 0
    do ix = 1, NX
      x = X_MIN + dx * (dble(ix - 1) + 0.5d0)
      w1 = w0 + w(ix)
      iw1 = min(int(NP_RandomForce * (w1 + EPS)), NP_RandomForce)
      do ip = iw0 + 1, iw1
        RandomForceArray(ip) = x
      end do
      iw0 = iw1
      w0 = w1
    end do

  end subroutine InitRandomForce

  subroutine simu_state_trans()
    implicit none
    real(8), parameter :: t_scale = 1.d-9   !ns to s
    integer :: i
    real(8) :: unifrd, rnd

    !$omp parallel do default(none) &
    !$omp private(i,rnd)&
    !$omp shared(state,x_L,RndStateSeedArray,x_barrier,x_S)&
    !$omp shared(attach_count,detach_pre_count,detach_pos_count)
    do i = 1, np
      rnd = unifrd(RndStateSeedArray(i))
      if (state(i) == 0) then
        if (rnd <= t_scale*a_trans*dt)  then
          state(i) = 1
          attach_count(i) = attach_count(i) + 1
        end if
      else
        if (x_L(i) < x_barrier) then
          if (rnd <= t_scale*d_trans*dt) then
            state(i) = 0
            x_L(i) = 0.d0
            x_S(i) = 0.d0
            detach_pre_count(i) = detach_pre_count(i)+1
          end if
        else
          if (rnd <= t_scale*g_trans*dt .or. x_L(i)+x_S(i) <= x_min_detach) then
            state(i) = 0
            x_L(i) = 0.d0
            x_S(i) = 0.d0
            detach_pos_count(i) = detach_pos_count(i)+1
          end if
        end if
      end if
    end do

  end subroutine simu_state_trans

  subroutine simu_force_random()
    implicit none

    integer :: i, ir
    real(8) :: rnd
    real(8) :: unifrd

    !$omp parallel do default(none) &
    !$omp private(i,rnd,ir)&
    !$omp shared(force_random,RndForceSeedArray,RandomForceArray)
    do i = 1, np
      rnd = unifrd(RndForceSeedArray(i))
      ir = max(int(rnd * NP_RandomForce), 1)
      force_random(i) = RandomForceArray(ir)
    end do
    
  end subroutine simu_force_random

  subroutine simu_force_powerstroke()
    implicit none
    integer :: i

    !$omp parallel do default(none)&
    !$omp private(i)&
    !$omp shared(x_L,phi,dphi,d2phi,stiff_ps,force,x_S,state)
    do i = 1, np
      if (state(i) == 1) then
        call powerstroke_potential_bind(x_L(i), x_S(i), phi(i), dphi(i), d2phi(i))
        force(i) = -dphi(i)
        if ( d2phi(i) > 0.d0) then
          stiff_ps(i)  = omega_stiff*d2phi(i)
        else
          stiff_ps(i)  = c_minus*dabs(d2phi(i))
        end if
      else
        force(i) = -k_spring*x_S(i)
        stiff_ps(i) = k_spring
      end if
    end do
    
  end subroutine simu_force_powerstroke

  subroutine powerstroke_potential(x, phi, dphi, d2phi)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(out) :: phi, dphi, d2phi
    
    real(8) :: A, dA, d2A, B, dB, d2B, DS, dDS, d2DS
    
    A = 0.5*c_pre*(x-x_pre)**2 + E_pre
    dA = c_pre*(x-x_pre)
    d2A = c_pre
    B = 0.5*c_pos*(x-x_pos)**2 + E_pos
    dB = c_pos*(x-x_pos)
    d2B = c_pos
    DS = sqrt((A - B )**2 + 2.0*delta**2)
      
    dDS = (A-B)*(dA-dB)/DS
    d2DS = (d2A-d2B)*(A-B)/DS + (dA-dB)**2/DS - (dA-dB)*(A-B)*dDS/DS**2

    phi = 0.5*( A + B - DS)
    dphi = 0.5*(dA + dB - dDS)
    d2phi = 0.5*(d2A + d2B - d2DS)
  end subroutine powerstroke_potential

  subroutine powerstroke_potential_bind(x_L, x_S, phi, dphi, d2phi)
    implicit none
    real(8), intent(in) :: x_L, x_S
    real(8), intent(out) :: phi, dphi, d2phi
    
    real(8) :: A, dA, d2A, B, dB, d2B, DS, dDS, d2DS
    
    A = 0.5*c_pre*(x_L-x_pre)**2 + E_pre
    dA = c_pre*(x_L-x_pre)
    d2A = c_pre
    B = 0.5*c_pos*(x_L-x_pos)**2 + E_pos
    dB = c_pos*(x_L-x_pos)
    d2B = c_pos
    DS = sqrt((A - B )**2 + 2.0*delta**2)
      
    dDS = (A-B)*(dA-dB)/DS
    d2DS = (d2A-d2B)*(A-B)/DS + (dA-dB)**2/DS - (dA-dB)*(A-B)*dDS/DS**2

    phi = 0.5*( A + B - DS) + 0.5*k_spring*(x_S+x_L)**2
    dphi = 0.5*(dA + dB - dDS) + k_spring*(x_S+x_L)
    d2phi = 0.5*(d2A + d2B - d2DS) + k_spring
  end subroutine powerstroke_potential_bind

  subroutine simu_update()
    implicit none
    integer :: i
    real(8) :: coef, fric_all, vel, fz, dz
    real(8) :: time_scale = 1.d-9   !ns to s
    integer :: na
    real(8) :: naPerAF

    fz = 0.d0
    na = 0
    !$omp parallel do default(none) &
    !$omp private(i,coef,vel) &
    !$omp shared(force,force_random,x_L,state,x_S) 
    do i = 1, np
      if (state(i) == 1) then
!C        fz = fz + k_spring*(x_S(i) + x(i))
!C        na = na + 1
        coef = sqrt(2.0*fric_x*KB_T/dt)
        vel = (1.d0/fric_x)*(force(i) + coef*force_random(i))
        x_L(i) = x_L(i) + dt*vel
      else
        coef = sqrt(2.0*fric_d*KB_T/dt)
        vel = (1.d0/fric_d)*(force(i) + coef*force_random(i))
        x_S(i) = x_S(i) + dt*vel
      end if
    end do
  end subroutine simu_update

  subroutine  get_StatisticalData()
    implicit none
    integer :: i
    integer :: n_attach, n_detach_pre, n_detach_pos
    
    npre = 0
    npos = 0
    ndet = 0
    n_attach = 0
    n_detach_pre = 0
    n_detach_pos = 0
    do i = 1, np
      n_attach = n_attach + attach_count(i)
      n_detach_pre = n_detach_pre + detach_pre_count(i)
      n_detach_pos = n_detach_pos + detach_pos_count(i)
      if (state(i) == 1) then
        if (x_L(i) <= x_barrier) then
          npre = npre + 1
        else
          npos = npos + 1
        end if
      else
        ndet = ndet + 1
      end if
    end do
    NattachPerMol = dble(n_attach)/dble(np)
    Ndeatch_prePerMol = dble(n_detach_pre)/dble(np)
    Ndeatch_posPerMol = dble(n_detach_pos)/dble(np)
  end subroutine get_StatisticalData

end program double_well_test

real(8) function unifrd(rndSeed)
  implicit none
  integer, intent(inout) :: rndSeed
  integer, parameter :: LAMBDA = 843314861, C = 453816693, T30 = 2**30
  real(8), parameter :: MYU = 2.d0**31
  rndSeed = LAMBDA * rndSeed + C
  if (rndSeed < 0) rndSeed = (rndSeed + T30) + T30
  unifrd = dble(rndSeed) / MYU
end function unifrd

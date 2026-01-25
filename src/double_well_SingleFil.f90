program double_well_test
  implicit none
  integer, parameter :: np = 40320
    integer, parameter :: SwitchZupdate = 0
  !Energy Parameters
  real(8), parameter :: KB_T = 0.0138 * 310                !pN*nm Boltzmann Constact x Temparature
  real(8), parameter :: E_ATP = 22.5*KB_T                  !pN*nm Energy of ATP hydrolysis
  real(8), parameter :: c_pre = 8.0, c_pos = 8.0           !pN/nm Curvature of the two wells
  real(8), parameter :: E_pre = 0.7*E_ATP, E_pos = 0.0     !pN*nm Assuming 70% of E_ATP is used for powerstroke
  real(8), parameter :: x_pre = 0.0, x_pos = 8.5           !nm    Levearm end position
  real(8), parameter :: delta = 2.0*KB_T                   !pN*nm Barrier relaxation
  real(8), parameter :: omega_stiff = 1.0, c_minus = 2.5   !Unitless Stiffness coefficient
  real(8), parameter :: k_spring = 2.8                     !pN/nm Spring
  real(8), parameter :: x_S0 = 0.0                         !nm   Spring energy = 0.5*k_spring*(x_S+x_L)
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
  integer, parameter :: nM = 80                            !Number of myosins per one AF
  real(8), parameter :: gamma_sarco = 1.d-5                !pN*s/nm Viscosity per one AF
  real(8), parameter :: kZ = 8.0                           !pN/nm   Spring constant per one AF

  !Time step
  real(8), parameter :: dt = 0.5   !ns less than fric_x/max(c_pre,c_pos)
  integer, parameter :: nt_in = 100000/dt, nt_out = 200  !Fine record
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
  

  call draw_powerstroke_potential()
  call initialize()
  call initRandomForce()

  do it_out = 1, nt_out
    total_time = it_out*dt*nt_in
    do it_in = 1, nt_in
      call simu_force_random()
      call simu_force_powerstroke()
      call simu_update()
      call simu_state_trans()
    end do
    call get_StatisticalData()
    write(6,'("it_out,time[s],npre,npos,Trans(3),z,F=",i8,e13.5,3i8,5e13.5)') it_out,1.d-9*total_time,npre,npos,ndet,&
      NattachPerMol,Ndeatch_prePerMol,Ndeatch_posPerMol,z,FzPerAF
    if (it_out == 1) then
      open(20,file="trans.csv")
      open(21,file="data.bin",form="unformatted")
     else
       open(20,file="trans.csv",position="append")
       open(21,file="data.bin",position="append",form="unformatted")
     end if
     write(20,'(i8,10e13.5)') it_out,1.d-9*total_time,dble(npre)/dble(np),dble(npos)/dble(np),dble(ndet)/dble(np),&
       NattachPerMol,Ndeatch_prePerMol,Ndeatch_posPerMol,z,FzPerAF
     close(20)
     write(21) it_out
     write(21) z
     write(21) state
     write(21) x_L
     write(21) x_S
     close(21)
     
  end do
contains
  subroutine draw_powerstroke_potential()
    implicit none
    integer, parameter :: nx = 1000
    real(8), parameter :: x_min = -5.0, x_max = 15.0
    integer :: i, ic
    real(8) :: omega, x, phi, dphi, d2phi
    real(8) :: phi_save(0:nx)

    open(10,file="PowerStrokePot.csv")
    write(10,*) "#x, phi, dphi/dx, d^2phi/dx^2"
    do i = 0, nx
      omega = dble(i)/dble(nx)
      x = (1.d0-omega)*x_min + omega*x_max
      call powerstroke_potential(x, phi, dphi, d2phi)
      write(10,'(4f12.4)') x, phi, dphi, d2phi
      phi_save(i) = phi
    end do
    close(10)

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

  subroutine initialize()
    implicit none
    integer :: i
    z = 0.d0
    x_L = 0.d0
    x_S = x_S0
    state = 0
    attach_count = 0
    detach_pre_count = 0
    detach_pos_count = 0
    do i = 1, np
!C      RndForceSeedArray(i) = i
!C      RndStateSeedArray(i) = np+i
      RndForceSeedArray(i) = i+20
      RndStateSeedArray(i) = np+i+20
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
    !$omp shared(force,force_random,x_L,state,x_S) &
    !$omp reduction(+:fz,na)
    do i = 1, np
      if (state(i) == 1) then
        fz = fz + k_spring*(x_S(i) + x_L(i))
        na = na + 1
        coef = sqrt(2.0*fric_x*KB_T/dt)
        vel = (1.d0/fric_x)*(force(i) + coef*force_random(i))
        x_L(i) = x_L(i) + dt*vel
      else
        coef = sqrt(2.0*fric_d*KB_T/dt)
        vel = (1.d0/fric_d)*(force(i) + coef*force_random(i))
        x_S(i) = x_S(i) + dt*vel
      end if
    end do
    
    FzPerAF = fz/(dble(np)/dble(nM))
    naPerAF = dble(na)/(dble(np)/dble(nM))
    if (SwitchZupdate == 2) then
      ! Solve gamma_sarco*(z(t+dt)-z(t))/(dt*time_scale) = FzPerAF - naPerAF*k_spring*dz - kZ*z(t+dt)
      ! i.e. gamma_sarco*dz/(dt*time_scale) = FzPerAF - naPerAF*k_spring*dz - kZ*(z(t)+dz)
      ! i.e. [gamma_sarco/(dt*time_scale) + naPerAF*k_spring + kZ]*dz = FzPerAF - kZ*z(t)
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + naPerAF*k_spring + kZ)
    else if (SwitchZupdate == 1) then
      ! Solve gamma_sarco*(z(t+dt)-z(t))/(dt*time_scale) = FzPerAF - kZ*z(t+dt)
      ! i.e. gamma_sarco*dz/(dt*time_scale) = FzPerAF - kZ*(z(t)+dz)
      ! i.e. [gamma_sarco/(dt*time_scale)+kZ]*dz = FzPerAF - kZ*z(t)
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + kZ)
    else
      ! Solve gamma_sarco*(z(t+dt)-z(t))/(dt*time_scale) = FzPerAF - kZ*z(t)
      ! i.e. gamma_sarco*dz/(dt*time_scale) = FzPerAF - kZ*z(t)
      ! i.e. [gamma_sarco/(dt*time_scale)]*dz = FzPerAF - kZ*z(t)
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale))      
    end if
    z = z + dz

    !$omp parallel do default(none) &
    !$omp private(i) &
    !$omp shared(state,x_S,dz)
    do i = 1, np
      if (state(i) == 1) x_S(i) = x_S(i) - dz
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

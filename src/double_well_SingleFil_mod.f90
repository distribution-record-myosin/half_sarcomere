program double_well_test
  implicit none
  include 'param.inc'
!C  integer, parameter :: np = 800
!C  integer, parameter :: np = 5040
!C  integer, parameter :: np = 10080
!C  integer, parameter :: np = 20160
  integer, parameter :: np = 40320
  logical, parameter :: SI_Switch = .false.
  integer, parameter :: SwitchZupdate = 2
  !integer, parameter :: SwitchZupdate = 21
!   integer, parameter :: SampleSwitch = 31
!C in parac.inc  integer, parameter :: SampleSwitch = 3
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
  real(8), parameter :: x_shift0 = 0.0                     !nm   Spring energy = 0.5*k_spring*(x_shift+x)
  !C 1/2 kspring xpos**2 == E_ATP　ATP加水分解のエネルギー

  real(8) :: x_barrier
  !Friction
  real(8), parameter :: fric_x = 80.0 !pN*ns/nm  friction for leverarm rotation
  real(8), parameter :: fric_d = 80.0 !pN*ns/nm  friction for x_shift during detachment
  real(8), parameter :: x_min_detach = -10.0
  !Transition
  real(8), parameter :: a_trans = 500, d_trans = 5000, g_trans=100  !1/s Transition rate constants
  !Sarcomere (per one active filament(AF))
  real(8) :: z                                             !nm Contraction distance
  real(8) :: FzPerAF                                       !pN Contraction Force per one AF
!C  integer, parameter :: nM = 20                            !Number of myosins per one AF
  integer, parameter :: nM = 80                            !Number of myosins per one AF
  real(8), parameter :: gamma_sarco = 1.d-5                !pN*s/nm Viscosity per one AF
!C  real(8), parameter :: kZ = 0.5                           !pN/nm   Spring constant per one AF
  real(8), parameter :: kZ = 8.0                           !pN/nm   Spring constant per one AF
  !Time step
  real(8), parameter :: dict_dt = 0.5
  !real(8), parameter :: dt = 0.5   !ns less than fric_x/max(c_pre,c_pos)
  !integer, parameter :: nt_in = 200000, nt_out = 10000
  !integer, parameter :: nt_in = 20000, nt_out = 500000  !Fine record
  !dt = 1.0
  !  real(8), parameter :: dt = 1.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 100000, nt_out = 10000
  !dt = 2.0
  !  real(8), parameter :: dt = 2.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 50000, nt_out = 10000
  !dt = 10.0
  !  real(8), parameter :: dt =  10.0   !ns less than fric_x/max(c_pre,c_pos)
  !  integer, parameter :: nt_in = 10000, nt_out = 10000
  
!  real(8), parameter :: dt =  100.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 100, nt_out = 100000
!  integer, parameter :: nt_in = 1000, nt_out = 10000
  
  !real(8), parameter :: dt =  500.0   !ns less than fric_x/max(c_pre,c_pos)
  !integer, parameter :: nt_in = 200, nt_out = 10000  

  real(8) :: dt   !ns less than fric_x/max(c_pre,c_pos)
  integer :: nt_in, nt_out

 !real(8), parameter :: dt =  1000.0   !ns less than fric_x/max(c_pre,c_pos)
 !integer, parameter :: nt_in = 100, nt_out = 10000
 !integer, parameter :: nt_in = 100, nt_out = 5000
 
 !real(8), parameter :: dt =  700.0   !ns less than fric_x/max(c_pre,c_pos)
 !integer, parameter :: nt_in = 100, nt_out = 10000

!  real(8), parameter :: dt =  2000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 50, nt_out = 10000  

!  real(8), parameter :: dt =  5000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 20, nt_out = 10000  

!  real(8), parameter :: dt =  20000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 5, nt_out = 10000  

!  real(8), parameter :: dt =  25000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 4, nt_out = 10000  
!
!  real(8), parameter :: dt =  50000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 2, nt_out = 10000  

!  real(8), parameter :: dt =  100000.0   !ns less than fric_x/max(c_pre,c_pos)
!  integer, parameter :: nt_in = 1, nt_out = 10000  


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
  real(8) :: x(np), force(np), force_random(np)
  real(8) :: x_shift(np)
  real(8) :: phi(np), dphi(np), d2phi(np), stiff_ps(np)
  integer :: RndForceSeedArray(np)  ! random force seeds
  integer :: RndStateSeedArray(np)  ! state transition seeds
  
  integer, parameter :: dict_psz = 1000 ! 各ファイルに含まれる数値の数
  ! integer, parameter :: dict_2d_xlowb = -7
  ! integer, parameter :: dict_2d_xupb = 13
  ! integer, parameter :: dict_2d_xsz = 30   ! 読み込むファイル数 (適宜変更)      
  ! integer, parameter :: dict_2d_slowb = -19
  ! integer, parameter :: dict_2d_supb = 20
  ! integer, parameter :: dict_2d_shift_sz = 40   ! 読み込むファイル数 (適宜変更)
  !C in parac.inc integer, parameter :: dict_2d_xlowb = -7
  !C in parac.inc integer, parameter :: dict_2d_xupb = 14
  !C in parac.inc integer, parameter :: dict_2d_xsz = 100   ! 読み込むファイル数 (適宜変更)      
  !C in parac.inc integer, parameter :: dict_2d_slowb = -19
  !C in parac.inc integer, parameter :: dict_2d_supb = 23
  !C in parac.inc integer, parameter :: dict_2d_shift_sz = 200   ! 読み込むファイル数 (適宜変更)
  !C in parac.inc integer, parameter :: dict_1d_xlowb = -6
  !C in parac.inc integer, parameter :: dict_1d_xupb = 6
  !C in parac.inc integer, parameter :: dict_1d_xsz = 100

  !C in parac.inc integer, parameter :: dictp_2d_xlowb = -3
  !C in parac.inc integer, parameter :: dictp_2d_xupb = 10
  !C in parac.inc integer, parameter :: dictp_2d_xsz = 200   ! 読み込むファイル数 (適宜変更)      
  !C in parac.inc integer, parameter :: dictp_2d_slowb = -3
  !C in parac.inc integer, parameter :: dictp_2d_supb = 3
!C  integer, parameter :: dictp_2d_slowb = -2
!C  integer, parameter :: dictp_2d_supb = 2
!C  integer, parameter :: dictp_2d_shift_sz = 200   ! 読み込むファイル数 (適宜変更)
!C  integer, parameter :: dictp_2d_shift_sz = 300   ! 読み込むファイル数 (適宜変更)
  !C in parac.inc integer, parameter :: dictp_2d_shift_sz = 30   ! 読み込むファイル数 (適宜変更)


  real(8) :: dict_2d(dict_psz, dict_2d_xsz, dict_2d_shift_sz)
!  real(8) :: dict_2d_band(dict_2d_xsz, dict_2d_shift_sz)
  real(8) :: dict_2d_ave(dict_2d_xsz, dict_2d_shift_sz)

  real(8) :: dict_1d(dict_psz, dict_1d_xsz)
!  real(8) :: dict_1d_band(dict_1d_xsz)  
  real(8) :: dict_1d_ave(dict_1d_xsz)  

  real(8) :: dictp_2d(dict_psz, dictp_2d_xsz, dictp_2d_shift_sz)
! 

  !C-- SwitchZupdate == 4
  real(8) :: zdot_p, zdot, gamma

  dt = real(time,kind=8)
  nt_in  = 100000/time !C output every 100000 ns
  nt_out = 10000.0*(time/1000.0) !C 10^9 ns

  if(SwitchZupdate==4) then
    zdot_p = 0
    gamma = 0.9
  end if
  
  !--
  ! integer :: index1, index2, index3
  ! real(8) :: v

  ! write(*,*) "hoehoe"
  call read_files(int(dt))
!  write(*,*) maxval(dict_2d_band), "hoho"


  call draw_powerstroke_potential()
  call initialize()
  call initRandomForce()


  do it_out = 1, nt_out
!C  do it_out = 1, 52819 !C ---
    total_time = it_out*dt*nt_in
    do it_in = 1, nt_in
      call simu_force_random()
      call simu_force_powerstroke()
      call simu_update()
      call simu_state_trans()
!C      if(it_out>51319) write(22,'(3e13.5)') z, maxval(x_shift),minval(x_shift)!C --- itout: 51319...52819の出力
    end do
    call get_StatisticalData()
    write(6,'("it_out,time[s],npre,npos,Trans(3),z,F=",i8,e13.5,3i8,5e13.5)') it_out,1.d-9*total_time,npre,npos,ndet,&
      NattachPerMol,Ndeatch_prePerMol,Ndeatch_posPerMol,z,FzPerAF
    if (it_out == 1) then
      open(20,file="trans.csv")
      open(21,file="data.bin",form="unformatted")
!C      open(22,file="debug.csv")!C --- 
    else
      open(20,file="trans.csv",position="append")
      open(21,file="data.bin",position="append",form="unformatted")
    end if
    write(20,'(i8,10e13.5)') it_out,1.d-9*total_time,dble(npre)/dble(np),dble(npos)/dble(np),dble(ndet)/dble(np),&
       NattachPerMol,Ndeatch_prePerMol,Ndeatch_posPerMol,z,FzPerAF
    close(20)
    if(mod(it_out,50)==0) then
      write(21) it_out 
      write(21) z 
      write(21) state 
      write(21) x 
      write(21) x_shift 
    end if 
   close(21)
  end do
!C  close(22)!C ---   
  
contains
  
  subroutine sample(v, vector, ave,i)
    implicit none
    real(8), intent(inout) :: v
    real(8), dimension(:)  :: vector
    real(8), intent(in) :: ave
    real(8) :: rnd, m
    integer, intent(in) :: i
    integer :: index,ir
    real(8) :: unifrd
    rnd = unifrd(RndForceSeedArray(i))
    index = max(min(int(rnd*dict_psz+1),dict_psz),1)
    m =vector(index)
!    rnd = unifrd(RndForceSeedArray(1))
!    ir = max(int(rnd * NP_RandomForce), 1)
!    v = m+haba*RandomForceArray(ir)
    v = v + m
  end subroutine sample

  subroutine sample_w_shift(v, vector, xlowb, xupb, dx, index, dict_ave, i)
    implicit none
    real(8), intent(inout) :: v
    real(8), dimension(:)  :: vector, dict_ave
    real(8), intent(in) :: xlowb, xupb, dx
    real(8) :: hosei, deltax
    real(8) :: rnd, m
    integer :: index,ir,dindex,ind
    integer, intent(in) :: i
    real(8) :: unifrd

    if (v < xlowb .or. v >= xupb) then
      hosei = 0
    else
      !C ind is set as (ind-1)th sampling point < x(i) < ind th sampling point
      if (xlowb + dx*(index-1) > v) then
        ind = index
      else
        ind = index+1
      end if
      deltax = v-(xlowb + dx*(ind-2))
!C    hosei = (1-deltax/dx)*dict_2d_ave(ind-1) + deltax/dx*(dict_2d_ave(ind)+dx) - deltax - dict_2d_ave(ind)
      hosei = (1-deltax/dx)*dict_ave(ind-1) + (deltax/dx-1)*dict_ave(ind) 
    end if

    rnd = unifrd(RndForceSeedArray(i))
    dindex = max(min(int(rnd*dict_psz+1),dict_psz),1)
    m =vector(dindex)
!    rnd = unifrd(RndForceSeedArray(1))
!    ir = max(int(rnd * NP_RandomForce), 1)
!    v = m+haba*RandomForceArray(ir)
    v = v + m + hosei

  end subroutine sample_w_shift

  subroutine sample_2d_w_shift(x, s, vector, xlowb, xupb, dx, index1, slowb, supb, ds, index2, dict_ave, i)
    implicit none
    real(8), intent(inout)  :: x
    real(8), intent(in)     :: s
    real(8), dimension(:)   :: vector
    real(8), dimension(:,:) :: dict_ave
    real(8), intent(in) :: xlowb, xupb, dx
    real(8), intent(in) :: slowb, supb, ds    
    real(8) :: hosei, deltax, deltas, save1, save2
    real(8) :: rnd, m, w_0_0,w_0_1,w_1_0,w_1_1
    integer, intent(in) :: index1, index2,i
    integer :: index,ir,dindex,xind,sind
    real(8) :: unifrd

    if(s>=slowb .and. s<supb) then
      if (slowb + ds*(index2-1) > s) then
        sind = index2
      else
        sind = index2+1
      end if
      deltas = s-(slowb + ds*(sind-2))      
    end if

    if(x>=xlowb .and. x<xupb) then
      if (xlowb + dx*(index1-1) > x) then
        xind = index1
      else
        xind = index1+1
      end if
      deltax = x-(xlowb + dx*(xind-2))      
    end if

    if((s < slowb .or. s >= supb) .and. (x < xlowb .or. x >= xupb)) then !C points out of range in both directions
      hosei = 0
    else if(s >= slowb .and. s < supb .and. x >= xlowb .and. x<xupb) then !C points in the range of both directions
      save1 = (1-deltas/ds)*dict_ave(xind-1,sind-1) + deltas/ds*dict_ave(xind-1, sind)
      save2 = (1-deltas/ds)*dict_ave(xind,  sind-1) + deltas/ds*dict_ave(xind,   sind)      

!C    hosei = (1-deltax/dx)*dict_2d_ave(ind-1) + deltax/dx*(dict_2d_ave(ind)+dx) - deltax - dict_2d_ave(ind)
!C    hosei = (1-deltax/dx)*save1 +  deltax/dx*(save2+dx) - deltax - dict_ave(index1,index2)
      hosei = (1-deltax/dx)*save1 +  deltax/dx*save2 - dict_ave(index1,index2)
    else if(x >= xlowb .and. x<xupb) then !C points out of range in x_shift
      save1 = dict_ave(xind-1,index2)
      save2 = dict_ave(xind,  index2)
      hosei = (1-deltax/dx)*save1 + deltax/dx*save2 - dict_ave(index1, index2)
    else                            !C points out of range in x
      save1 = (1-deltas/ds)*dict_ave(index1,sind-1) + deltas/ds*dict_ave(index1, sind)
      hosei = save1 - dict_ave(index1,index2)
    end if

    rnd = unifrd(RndForceSeedArray(i))
    dindex = max(min(int(rnd*dict_psz+1),dict_psz),1)
    m =vector(dindex)
    x = x + m + hosei

  end subroutine sample_2d_w_shift

  subroutine sample_2d(x, s, dict_2d, xlowb, xupb, dx, index1, slowb, supb, ds, &
    & index2, i, xsz, shift_sz)
    implicit none
    real(8), intent(inout)  :: x
    real(8), intent(in)     :: s
    real(8), dimension(:,:,:) :: dict_2d
    real(8), intent(in) :: xlowb, xupb, dx
    real(8), intent(in) :: slowb, supb, ds    
    real(8) :: hosei, deltax, deltas, save1, save2
    real(8) :: rnd, m, w_0_0,w_0_1,w_1_0,w_1_1
    integer, intent(in) :: index1, index2,i,xsz,shift_sz
    integer :: index,ir,dindex,xind,sind, dict_index1, dict_index2
    real(8) :: unifrd

    if(s>slowb .and. s<supb) then
      if (slowb + ds*(index2-1) >= s) then
        sind = index2
      else
        sind = index2+1
      end if
      deltas = s-(slowb + ds*(sind-2))      
    end if

    if(x>xlowb .and. x<xupb) then
      if (xlowb + dx*(index1-1) >= x) then
        xind = index1
      else
        xind = index1+1
      end if
      deltax = x-(xlowb + dx*(xind-2))      
    end if

    if((s <= slowb .or. s >= supb) .and. (x <= xlowb .or. x >= xupb)) then !C points out of range in both directions
      dict_index1=index1
      dict_index2=index2
    else if(s > slowb .and. s < supb .and. x > xlowb .and. x<xupb) then !C points in the range of both directions
      w_0_0=(1-deltax/dx)*(1-deltas/ds)
      w_0_1=(1-deltax/dx)*(deltas/ds)
      w_1_0=(deltax/dx)*(1-deltas/ds)
      w_1_1=(deltax/dx)*(deltas/ds)

      rnd = unifrd(RndForceSeedArray(i))
      if (rnd <= w_0_0) then
        dict_index1=xind-1
        dict_index2=sind-1
      else if (rnd <= w_0_0+w_1_0) then
        dict_index1=xind
        dict_index2=sind-1
      else if (rnd <= w_0_0+w_1_0+w_0_1) then
        dict_index1=xind-1
        dict_index2=sind
      else
        dict_index1=xind
        dict_index2=sind
      end if
    else if(x > xlowb .and. x<xupb) then !C points out of range in x_shift
      w_0_0 = (1-deltax/dx)
      rnd = unifrd(RndForceSeedArray(i))
      dict_index2=index2
      if (rnd <= w_0_0) then
        dict_index1=xind-1
      else 
        dict_index1=xind
      end if
    else                            !C points out of range in x
      w_0_0 = (1-deltas/ds)
      dict_index1=index1      
      rnd = unifrd(RndForceSeedArray(i))

      if (rnd <= w_0_0) then
        dict_index2=sind-1
      else 
        dict_index2=sind
      end if
    end if

    rnd = unifrd(RndForceSeedArray(i))
    dindex     = max(min(int(rnd*dict_psz+1), dict_psz),1)

    !C debug---    
    if(dict_index1<=0 .or. dict_index1>xsz) then
      write(*,*) "dict_index1,x,xlowb,xupb,dx,xsz", dict_index1, x,xlowb,xupb,dx,xsz
    end if
    if(dict_index2<=0 .or. dict_index2>shift_sz) then
      write(*,*) "dict_index2,x,xlowb,xupb,dx,ssz", dict_index2, s,slowb,supb,ds,shift_sz
    end if
    !C debug---        

    dict_index1= max(min(dict_index1,xsz),1)
    dict_index2= max(min(dict_index2,shift_sz),1)    
    m =dict_2d(dindex,dict_index1, dict_index2)
    x = x + m
  end subroutine sample_2d


  subroutine sample_1d(x, dict_1d, xlowb, xupb, dx, index, i)
    implicit none
    real(8), intent(inout)  :: x
    real(8), dimension(:,:) :: dict_1d
    real(8), intent(in) :: xlowb, xupb, dx
    real(8) :: deltax,w_0_0,w_0_1,w_1_0,w_1_1
    real(8) :: rnd, m
    integer, intent(in) :: index,i
    integer :: ir,dindex,xind,sind, dict_index
    real(8) :: unifrd

    if(x>xlowb .and. x<xupb) then
      if (xlowb + dx*(index-1) >= x) then
        xind = index
      else
        xind = index+1
      end if
      deltax = x-(xlowb + dx*(xind-2))      
    end if

    if(x <= xlowb .or. x >= xupb) then
      dict_index=index
    else
      w_0_0=(1-deltax/dx)
      rnd = unifrd(RndForceSeedArray(i))
      if (rnd <= w_0_0) then
        dict_index = xind-1
      else
        dict_index = xind
      end if
    end if

    rnd = unifrd(RndForceSeedArray(i))
    dindex     = max(min(int(rnd*dict_psz+1),dict_psz),1)

    !C debug---    
    if(dict_index<=0 .or. dict_index>dict_1d_xsz) then
      write(*,*) "dict_index,x,xlowb,xupb,dx,dict_1d_xsz", dict_index, x,xlowb,xupb,dx,dict_1d_xsz
    end if
    !C debug---        
    dict_index = max(min(dict_index, dict_1d_xsz),1)
    m =dict_1d(dindex,dict_index)
    x = x + m
  end subroutine sample_1d

  subroutine index_determine(index, x, lowb, upb, xsz)
    implicit none
    integer, intent(out) :: index
    integer, intent(in) :: xsz
    real(8), intent(in) :: x, lowb, upb
    real(8) :: step

        ! ステップ幅を計算
    step = (upb - lowb) / dble(xsz - 1)

        ! x の範囲チェック
    if (x <= lowb) then
        index = 1
    else if (x >= upb) then
        index = xsz
    else
        ! 直接 index を計算
        index = int((x - lowb) / step + 1.5)
    end if
  end subroutine index_determine

  subroutine read_files(time)
    implicit none
    integer, intent(in) :: time
    integer :: i, j, k, ios
    character(len=100) :: filename
    real(8) :: dx, v
    
    ! すべてのファイルをループして読み込む
    write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',&
      & time,'_',dict_2d_xlowb, '_',dict_2d_xupb,'_',dict_2d_xsz,'_',dict_2d_slowb,'_',dict_2d_supb,&
      & '_',dict_2d_shift_sz,'.csv'  
    write(*,*) filename
    open(unit=10, file=filename, status='old', action='read', iostat=ios, form="unformatted")
    if (ios /= 0) then
      print *, "Error opening file:", trim(filename)
      stop
    end if

    read(10) dict_2d(:,:,:)
    close(10)

    write(filename, '(A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',time, &
      & '_',dict_1d_xlowb,'_',dict_1d_xupb,'_',dict_1d_xsz,'.csv'  
    write(*,*) filename
    open(unit=10, file=filename, status='old', action='read', iostat=ios, form="unformatted")
    if (ios /= 0) then
      print *, "Error opening file:", trim(filename)
      stop
    end if
    read(10) dict_1d(:,:)
    close(10)

    write(filename, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'dist_records/dict1000_0.5_',time,'_',dictp_2d_xlowb,&
      & '_',dictp_2d_xupb,'_',dictp_2d_xsz,'_',dictp_2d_slowb,'_',dictp_2d_supb,'_', &
      & dictp_2d_shift_sz,'.csv'  
    write(*,*) filename
    open(unit=10, file=filename, status='old', action='read', iostat=ios, form="unformatted")
    if (ios /= 0) then
      print *, "Error opening file:", trim(filename)
      stop
    end if

    read(10) dictp_2d(:,:,:)
    close(10)

    dx = (dict_2d_xupb - dict_2d_xlowb) / dble(dict_2d_xsz-1)
    do i = 0, dict_2d_shift_sz-1    
      do j = 0, dict_2d_xsz-1
        !C--
        dict_2d(:, j+1, i+1) = dict_2d(:, j+1, i+1) - (dict_2d_xlowb+dx*j)
        v=0.0
        do k=1, dict_psz
          v = v+dict_2d(k,j+1,i+1)
        end do
        dict_2d_ave(j+1,i+1)=v/dble(dict_psz)
        !C--        
      end do
    end do

    dx = (dict_1d_xupb - dict_1d_xlowb) / dble(dict_1d_xsz-1)
    ! すべてのファイルをループして読み込む
    do j = 0, dict_1d_xsz-1
      !C--
      dict_1d(:, j+1) = dict_1d(:, j+1) - (dict_1d_xlowb + dx*j)
      v=0.0
      do k=1, dict_psz
        v = v+dict_1d(k,j+1)
      end do
      dict_1d_ave(j+1)=v/dble(dict_psz)
      !C--        
    end do


    dx = (dictp_2d_xupb - dictp_2d_xlowb) / dble(dictp_2d_xsz-1)
    do i = 0, dictp_2d_shift_sz-1    
      do j = 0, dictp_2d_xsz-1
        !C--
        dictp_2d(:, j+1, i+1) = dictp_2d(:, j+1, i+1) - (dictp_2d_xlowb+dx*j)
!C        v=0.0
!C        do k=1, dict_psz
!C          v = v+dictp_2d(k,j+1,i+1)
!C        end do
!C        dict_2d_ave(j+1,i+1)=v/dble(dict_psz)
      end do
    end do


    print *, "Data successfully loaded."
  end subroutine read_files

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
    x = 0.d0
    x_shift = x_shift0
    state = 0
    attach_count = 0
    detach_pre_count = 0
    detach_pos_count = 0
    do i = 1, np
      RndForceSeedArray(i) = i
      RndStateSeedArray(i) = np+i
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
!$omp shared(state,x,RndStateSeedArray,x_barrier,x_shift)&
!$omp shared(attach_count,detach_pre_count,detach_pos_count,dt)
    do i = 1, np
      rnd = unifrd(RndStateSeedArray(i))
      if (state(i) == 0) then
        if (rnd <= t_scale*a_trans*dt)  then
          state(i) = 1
          attach_count(i) = attach_count(i) + 1
        end if
      else
        if (x(i) < x_barrier) then
          if (rnd <= t_scale*d_trans*dt) then
            state(i) = 0
            x(i) = 0.d0
            x_shift(i) = 0.d0
            detach_pre_count(i) = detach_pre_count(i)+1
          end if
        else
          if (rnd <= t_scale*g_trans*dt .or. x(i)+x_shift(i) <= x_min_detach) then
            state(i) = 0
            x(i) = 0.d0
            x_shift(i) = 0.d0
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

  subroutine simu_force_random_i(i)
    implicit none

    integer :: i, ir
    real(8) :: rnd
    real(8) :: unifrd

    rnd = unifrd(RndForceSeedArray(i))
    ir = max(int(rnd * NP_RandomForce), 1)
    force_random(i) = RandomForceArray(ir)
    
  end subroutine simu_force_random_i

  subroutine simu_force_powerstroke()
    implicit none
    integer :: i

!$omp parallel do default(none)&
!$omp private(i)&
!$omp shared(x,phi,dphi,d2phi,stiff_ps,force,x_shift,state)
    do i = 1, np
      if (state(i) == 1) then
        call powerstroke_potential_bind(x(i), x_shift(i), phi(i), dphi(i), d2phi(i))
        force(i) = -dphi(i)
        if ( d2phi(i) > 0.d0) then
          stiff_ps(i)  = omega_stiff*d2phi(i)
        else
          stiff_ps(i)  = c_minus*dabs(d2phi(i))
        end if
      else
        force(i) = -k_spring*x_shift(i)
        stiff_ps(i) = k_spring
      end if
    end do
    
  end subroutine simu_force_powerstroke

  subroutine simu_force_powerstroke_i(i)
    implicit none
    integer :: i

    if (state(i) == 1) then
      call powerstroke_potential_bind(x(i), x_shift(i), phi(i), dphi(i), d2phi(i))
      force(i) = -dphi(i)
      if ( d2phi(i) > 0.d0) then
        stiff_ps(i)  = omega_stiff*d2phi(i)
      else
        stiff_ps(i)  = c_minus*dabs(d2phi(i))
      end if
    else
      force(i) = -k_spring*x_shift(i)
      stiff_ps(i) = k_spring
    end if
    
  end subroutine simu_force_powerstroke_i

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

  subroutine powerstroke_potential_bind(x, x_shift, phi, dphi, d2phi)
    implicit none
    real(8), intent(in) :: x, x_shift
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

    phi = 0.5*( A + B - DS) + 0.5*k_spring*(x_shift+x)**2
    dphi = 0.5*(dA + dB - dDS) + k_spring*(x_shift+x)
    d2phi = 0.5*(d2A + d2B - d2DS) + k_spring
  end subroutine powerstroke_potential_bind

  subroutine simu_update()
    implicit none
    integer :: it_dt,i
    real(8) :: coef, fric_all, vel, fz, dz
    real(8) :: time_scale = 1.d-9   !ns to s
    integer :: na
    real(8) :: naPerAF

    integer :: index1, index2, index3
    real(8) :: hosei, dx_1d, dx_2d, ds_2d,  dx_p2d, ds_p2d

    fz = 0.d0
    na = 0
    dx_2d=   (dict_2d_xupb-dict_2d_xlowb)/dble(dict_2d_xsz-1)
    ds_2d=   (dict_2d_supb-dict_2d_slowb)/dble(dict_2d_shift_sz-1)
    dx_1d=   (dict_1d_xupb-dict_1d_xlowb)/dble(dict_1d_xsz-1)
    dx_p2d=   (dictp_2d_xupb-dictp_2d_xlowb)/dble(dictp_2d_xsz-1)
    ds_p2d=   (dictp_2d_supb-dictp_2d_slowb)/dble(dictp_2d_shift_sz-1)

!$omp parallel do default(none) &
!$omp private(i,coef,vel,index1,index2,index3,it_dt) &
!$omp shared(force,force_random,x,state,x_shift,dict_2d,dict_2d_ave) &
!$omp shared(dx_2d,ds_2d,dict_1d,dict_1d_ave,dx_1d) &
!$omp shared(dx_p2d,ds_p2d,dictp_2d,dt) &
!$omp reduction(+:fz,na)
    do i = 1, np
      if (state(i) == 1) then
!C          fz = fz + k_spring*(x_shift(i) + x(i))
!C          na = na + 1
        coef = sqrt(2.0*fric_x*KB_T/dt)
!          vel = (1.d0/fric_x)*(force(i) + coef*force_random(i))
!          x(i) = x(i) + dt*vel
        call index_determine(index1, x(i),       real(dict_2d_xlowb,kind=8), real(dict_2d_xupb,kind=8), dict_2d_xsz)
        call index_determine(index2, x_shift(i), real(dict_2d_slowb,kind=8), real(dict_2d_supb,kind=8), dict_2d_shift_sz)  

        if(SampleSwitch == 1) then !C nearest distribution record is used
            call sample(x(i), dict_2d(:,index1,index2), dict_2d_ave(index1, index2),i)
        else if(SampleSwitch == 2) then !C nearest record is shifted to match average of neighboring records before sampling
            call sample_2d_w_shift(x(i), x_shift(i), dict_2d(:,index1,index2),         &
              &      dble(dict_2d_xlowb), dble(dict_2d_xupb), dx_2d, index1,             &
              &      dble(dict_2d_slowb), dble(dict_2d_supb), ds_2d, index2, dict_2d_ave, i) 
        else if(SampleSwitch == 3) then !C weighted sampling with probabilistic record switching
            call sample_2d(x(i), x_shift(i), dict_2d,                          &
              & dble(dict_2d_xlowb), dble(dict_2d_xupb), dx_2d, index1,      &
              & dble(dict_2d_slowb), dble(dict_2d_supb), ds_2d, index2,i,    &
              & dict_2d_xsz, dict_2d_shift_sz)
        else if(SampleSwitch == 31) then !C in addition, precise record is added.
          if(x(i) >= dictp_2d_xlowb .and. x(i) <= dictp_2d_xupb .and. &
              & x_shift(i) >= dictp_2d_slowb .and. x_shift(i) <= dictp_2d_supb) then
            call index_determine(index1, x(i),       real(dictp_2d_xlowb,kind=8), &
              &   real(dictp_2d_xupb,kind=8), dictp_2d_xsz)
            call index_determine(index2, x_shift(i), real(dictp_2d_slowb,kind=8), &
              &   real(dictp_2d_supb,kind=8), dictp_2d_shift_sz)  

            call sample_2d(x(i), x_shift(i), dictp_2d,                          &
                & dble(dictp_2d_xlowb), dble(dictp_2d_xupb), dx_p2d, index1,       &
                & dble(dictp_2d_slowb), dble(dictp_2d_supb), ds_p2d, index2,i,     &
                & dictp_2d_xsz, dictp_2d_shift_sz)
          else
            call sample_2d(x(i), x_shift(i), dict_2d,                          &
                & dble(dict_2d_xlowb), dble(dict_2d_xupb), dx_2d, index1,          &
                & dble(dict_2d_slowb), dble(dict_2d_supb), ds_2d, index2,i,       &
                & dict_2d_xsz, dict_2d_shift_sz)
          endif
        else
          write(*,*) "Unknown Switch:",SampleSwitch
          stop
        end if

        fz = fz + k_spring*(x_shift(i) + x(i))
        na = na + 1

      else
        coef = sqrt(2.0*fric_d*KB_T/dt)
!          vel = (1.d0/fric_d)*(force(i) + coef*force_random(i))
!          x_shift(i) = x_shift(i) + dt*vel
        call index_determine(index3, x_shift(i), real(dict_1d_xlowb,kind=8), real(dict_1d_xupb,kind=8), dict_2d_xsz)
        if(SampleSwitch == 1) then
            call sample(x_shift(i),dict_1d(:,index3),dict_1d_ave(index3), i)
        else if(SampleSwitch == 2) then 
            call sample_w_shift(x_shift(i),dict_1d(:,index3), dble(dict_1d_xlowb), dble(dict_1d_xupb), &
            &       dx_1d, index3, dict_1d_ave(:),i)
        else if(SampleSwitch == 3 .or. SampleSwitch==31) then
            call sample_1d(x_shift(i), dict_1d, dble(dict_1d_xlowb), dble(dict_1d_xupb), &
              &      dx_1d, index3, i)
        else 
            write(*,*) "Unknown Switch"
            stop
        end if
      end if
    end do
      
    FzPerAF = fz/(dble(np)/dble(nM))
    naPerAF = dble(na)/(dble(np)/dble(nM))
    if(SwitchZupdate == 4) then
      zdot = (FzPerAF-naPerAF*k_spring*dt*(1-gamma)*zdot_p-kZ*(z+dt*(1-gamma)*zdot_p))/ &
      & ((gamma_sarco/time_scale)+naPerAF*k_spring*dt*gamma+kZ*dt*gamma)
      dz = dt*(gamma * zdot + (1-gamma) * zdot_p)
      zdot_p=zdot
    else if(SwitchZupdate == 3 ) then
      ! average between explicit and implicit methods
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + 0.5*naPerAF*k_spring + 0.5*kZ)
    else if (SwitchZupdate == 2) then
      ! Solve gamma_sarco*(z(t+dt)-z(t))/(dt*time_scale) = FzPerAF - naPerAF*k_spring*dz - kZ*z(t+dt)
      ! i.e. gamma_sarco*dz/(dt*time_scale) = FzPerAF - naPerAF*k_spring*dz - kZ*(z(t)+dz)
      ! i.e. [gamma_sarco/(dt*time_scale) + naPerAF*k_spring + kZ]*dz = FzPerAF - kZ*z(t)
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + naPerAF*k_spring + kZ)
    else if (SwitchZupdate == 21) then !C modified implicit method
!      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + naPerAF*k_spring + kZ)
      dz = (FzPerAF - kZ*z)/(gamma_sarco/(dt*time_scale) + naPerAF*(k_spring*8.0/(8.0+k_spring)) + kZ)
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
!$omp shared(state,x_shift,dz)
    do i = 1, np
      if (state(i) == 1) x_shift(i) = x_shift(i) - dz
    end do

!$omp parallel do default(none) &
!$omp private(i) &
!$omp shared(x,x_shift)
    do i=1, np
      if (x(i) < dict_2d_xlowb) then
        x(i)=dict_2d_xlowb
      else if (x(i) > dict_2d_xupb) then
        x(i)= dict_2d_xupb
      end if
      if (x_shift(i) < dict_2d_slowb) then
        x_shift(i)=dict_2d_slowb
      else if (x_shift(i) > dict_2d_supb) then
        x_shift(i)=dict_2d_supb
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
        if (x(i) <= x_barrier) then
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

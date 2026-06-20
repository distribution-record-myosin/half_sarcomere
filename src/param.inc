  !C time[ns] for one step of distribution record
  integer, parameter :: time=1000
  logical, parameter :: ConstantVelocitySwitch = .true.
!  logical, parameter :: ConstantVelocitySwitch = .false.
  real(8) :: ConstantVelocity  !100-10000 nm/s

!C simulation length: nt_out x 0.1 ms
  integer, parameter :: nt_out = 1000

!C potential modification
  integer, parameter :: DELTA_SCALE = 2         !C delta=DELTA_SCALE x KB_T 2,4,8

  !C record size specification
  integer, parameter :: dict_psz=15000

  !C 1d record
  integer, parameter :: dict_1d_xlowb = -6
  integer, parameter :: dict_1d_xupb = 6
  integer, parameter :: dict_1d_xsz = 100
  
  !C 2d record
  integer, parameter :: dict_2d_xlowb = -7
  integer, parameter :: dict_2d_xupb = 14
  integer, parameter :: dict_2d_xsz = 100   
  integer, parameter :: dict_2d_slowb = -19
  integer, parameter :: dict_2d_supb = 23
  integer, parameter :: dict_2d_shift_sz = 200 

  !C 2d precise record
  integer, parameter :: dictp_2d_xlowb = -3
  integer, parameter :: dictp_2d_xupb = 10
  integer, parameter :: dictp_2d_xsz = 200   
  integer, parameter :: dictp_2d_slowb = -3
  integer, parameter :: dictp_2d_supb = 3
  integer, parameter :: dictp_2d_shift_sz = 30   

  !C distribution record based method: 3 or 31
  !C 3: only use 1d and 2d records, 31:will use precise record in addition to 1d and 2d records.
  integer, parameter :: SampleSwitch = 3

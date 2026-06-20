program get_dist
  implicit none
  integer, parameter :: np = 40320
!C  integer, parameter :: np = 20160
  real(8), parameter :: x_min = -7.0, x_max = 14.0
  real(8), parameter :: x_shift_min = -19.0, x_shift_max = 23.0
  real(8) :: z, x(np), x_shift(np)
  integer :: state(np)
  integer :: it, it_out, ip

  write(6,*) "input it_out"
  read(5,*) it_out


!  open(11,file="fdic_0.5_1000/Data_2_31_np_40320_-3_3_300.bin",form="unformatted")
!  open(11,file="fdic_0.5_1000/Data_2_3_np_40320_-3_3_300.bin",form="unformatted")
  open(11,file="data.bin",form="unformatted")
!  open(11,file="fdic_0.5_1000/bData_2_31_np_40320_-3_3_300.bin",form="unformatted")
!  open(11,file="data.bin",form="unformatted")
100 read(11) it
  read(11) z
  read(11) state
  read(11) x
  read(11) x_shift
  do ip = 1, np
    if (state(ip) == 1) then
      if (x(ip) < x_min .or. x(ip) > x_max) then
        write(30,*) it, x(ip), x_shift(ip)
        write(6,*) it, x(ip), x_shift(ip)
      end if
      if (x_shift(ip) < x_shift_min .or. x_shift(ip) > x_shift_max) then
        write(30,*) it, x(ip), x_shift(ip)
        write(6,*) it, x(ip), x_shift(ip)
      end if
    end if
  end do
  if (it .ne. it_out) goto 100
  close(11)

  open(20,file="dist.csv")
  do ip = 1, np
    if (state(ip) == 1) then
      write(20,'(i8,f12.4,1x,f12.4)') ip, x(ip), x_shift(ip)
    end if
  end do
  close(20)

end program get_dist

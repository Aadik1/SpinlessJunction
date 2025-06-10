program DJDV
  implicit none
  integer :: i, N
  real*8, allocatable, dimension(:) :: Voltage, Current, dIdV
  real*8 :: dV

  call FileLen(N, 'Volt_Current_2.dat')

  allocate(Voltage(N), Current(N), dIdV(N))

  open(1, file='Volt_Current_2.dat', status='old')
  do i = 1, N
     read(1, *) Voltage(i), Current(i)
  end do
  close(1)

  !...Forward difference at first point
  dIdV(1) = (Current(2) - Current(1)) / (Voltage(2) - Voltage(1))
  
  !...Central difference for interior points
  do i = 2, N - 1
     dV = Voltage(i+1) - Voltage(i-1)
     dIdV(i) = (Current(i+1) - Current(i-1)) / dV
  end do
  
  !...Backward difference at last point
  dIdV(N) = (Current(N) - Current(N-1)) / (Voltage(N) - Voltage(N-1))

  
  open(12, file='dI_dV.dat', status='unknown')
  do i = 1, N
     write(12, *) Voltage(i), dIdV(i)
  end do
  close(12)
  
  deallocate(Voltage, Current, dIdV) 
end program DJDV

subroutine Filelen(N, filename)
  implicit none
  integer :: N, iostat
  character(len=300) :: line
  character(*), intent(in) :: filename
  
  open(unit=11, file=filename, status='old', iostat=iostat)
  if (iostat /= 0) then
     write(*,*) 'Error reading file', iostat
     stop
  end if
  N = 0
  do
     read(11, '(A)', iostat=iostat) line
     if (iostat /= 0) exit
     N = N + 1
  end do
  close(11) !...backspace instead of close? 
  
end subroutine Filelen

program StypeJunction_Spinless
  use DefineHamiltonian
  use GreensFunctions
  implicit none
  real*8 :: V1, Current, total_time
  integer ::  k, i, start_tick, end_tick, rate, max_count
  character(len=30) :: vfn
  logical :: first

!.......................Reads all the bond parameters ,length of the molecule etc. 
  
 !....creates runtime datasheet 
  open(22, file='runtime_datasheet.dat', status='unknown')  
  call SYSTEM_CLOCK(COUNT_RATE=rate, COUNT_MAX=max_count)
  if (rate .eq. 0) then
     print *, "Error: SYSTEM_CLOCK not supported or rate is 0."
     stop
  end if
  call SYSTEM_CLOCK(COUNT=start_tick)

  write(22,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>','Runtime Data Sheet','>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

  !.........................Deifnes Hamiltonian
  call input()
  
  allocate(H(Natoms,Natoms))
  allocate(GammaL(Natoms, Natoms)); allocate(GammaR(Natoms, Natoms))
  
  call Central_Hamiltonian(H, Natoms, 4.d0, -2.d0)
  
  GammaL = (0.d0, 0.d0); GammaL(1,1) = 1.d0
  GammaR = (0.d0, 0.d0); GammaR(Natoms, Natoms) = 1.d0
  
  allocate(Hub(Natoms))
  open(11, file='Hubbard_info.dat', status='unknown')
  do i = 1, Natoms
     Hub(i) = Gamma0/10.d0 !...read(11, *) Hub(i)
  end do
  close(11)

  print *, 'HUBBARD:', Hub

!.......................Finds the eigenvalues of the Central Hamiltonian and uses
!                       it to define the w-grid
  
  allocate(Eigenvec(Natoms, Natoms))  
  allocate(Ev(Natoms))

  Eigenvec = H
  call complex_eigen_symm_martix(Eigenvec, Natoms, Ev)

  w_init = Ev(1)-dw
  w_fin = Ev(Natoms)+up
  N_of_w = (w_fin - w_init)/delta
  write(*,*) 'N_of_w:', N_of_w, 'w_fin:', w_fin, 'w_init:', w_init
  allocate(omega(N_of_w))

  do i = 1, N_of_w 
     omega(i) = w_init + i*delta
  end do
  
  call PrintFunctions()
  deallocate(Eigenvec, Ev)

!.......................Allocates Greens Functions for the integrals within the Sigmas
!                       across all omegas and organises them into user defined type GF0
!                       and GFf in 'mod_Gfunctions.f90'

  allocate(GF0%r(Natoms, Natoms, N_of_w)) ;  allocate(GF0%a(Natoms, Natoms, N_of_w))
  allocate(GF0%L(Natoms, Natoms, N_of_w)) ;  allocate(GF0%G(Natoms, Natoms, N_of_w))

  allocate(GFf%L(Natoms, Natoms, N_of_w)) ;  allocate(GFf%G(Natoms, Natoms, N_of_w))
  allocate(GFf%R(Natoms, Natoms, N_of_w)) ;  allocate(GFf%A(Natoms, Natoms, N_of_w))

  allocate(G_nil(Natoms))
  
!.......................Allocates all the full self-energies and full Greens functions
!                       needed in the current

!  allocate(SigmaL(Natoms, Natoms)) ; allocate(SigmaR(Natoms, Natoms))
  
!...............calculate GR and GA for all voltages on the omega grid
  allocate(work1(Natoms, Natoms)) ; allocate(work2(Natoms, Natoms)) ; allocate(work3(Natoms, Natoms))
  
!.......................Calculates and plots Voltage vs Current curve  

  open(3, file='Print.dat', status='unknown')
  write(vfn,'(i0)') order
  open(30, file='Volt_Current_'//trim(vfn)//'.dat', status='unknown')  
  !  print *, 'Pre-voltage'
  first=.true.
  do k = 0, Volt_range
     V1 = 0.d0!V + k*0.05

     call SCF_GFs(V1,first)
     
     GF0%r=GFf%r ; GF0%a=GFf%a ; GF0%L=GFf%L ; GF0%G=GFf%G
     if(first) first=.false.
     
     write(30, *) V1, Current(V1)
     print *, 'Progress:', k/(Volt_range*0.01), '%', Current(V1)
     STOP
  end do
  
  close(30) 
  close(3)
  
  call SYSTEM_CLOCK(COUNT=end_tick)
  total_time = real(end_tick - start_tick)/real(rate)
  total_time = total_time/60.d0

  write(22,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
  write(22,*) 'Total Runtime:', total_time, 'mins'
  write(22,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
  close(22)
  deallocate(GFf%L, GFf%G,GFf%R, GFf%A)
  deallocate(GF0%r, GF0%a, GF0%L, GF0%G) 
  deallocate(work1, work2, work3)
  deallocate(H, Hub, omega)
  !deallocate(SigmaL, SigmaR)
  deallocate(GammaL, GammaR, G_nil)
end program StypeJunction_Spinless

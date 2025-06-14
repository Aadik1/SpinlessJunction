program StypeJunction_one
  use DefineHamiltonian
  use GreensFunctions
  implicit none
  real*8 :: V1, Current, start_time, end_time
  integer ::  k, i
  character(len=1) :: Spin_orbit
  character(len=30) :: vfn

!.......................Reads all the bond parameters ,length of the molecule etc. 
  
  !.........................Deifnes Hamiltonian
  open(22, file='runtime_datasheet.dat', status='unknown')  
  call CPU_TIME(start_time)

  write(22,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>','Runtime Data Sheet','>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

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

!.......................Finds the eigenvalues of the Central Hamiltonian and uses
!                       it to define the w-grid
  
  allocate(Eigenvec(Natoms, Natoms))  
  allocate(Ev(Natoms))

  Eigenvec = H
  call eigen_symm_matrix(Eigenvec, Natoms, Ev)

  w_init = Ev(1)-dw
  w_fin = Ev(Natoms)+up
  N_of_w = (w_fin - w_init)/delta
  
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

  !open(3, file='Print.dat', status='replace')
  write(vfn,'(i0)') order
  open(30, file='Volt_Current_'//trim(vfn)//'.dat', status='unknown')  
!  print *, 'Pre-voltage' 
  do k = 0, Volt_range
     V1 = V + k*0.05
  call SCF_GFs(V1)
     write(30, *) V1, Current(V1)
     print *, 'Progress:', k/(Volt_range*0.01), '%', Current(V1)
  end do

  close(30) 
  ! close(3)

  call CPU_TIME(end_time)

  write(22,*) 'Total Runtime:', (start_time-end_time), 'mins'
  close(22)
  deallocate(GFf%L, GFf%G,GFf%R, GFf%A)
  deallocate(GF0%r, GF0%a, GF0%L, GF0%G) 
  deallocate(work1, work2, work3)
  deallocate(H, Hub, omega)
  !deallocate(SigmaL, SigmaR)
  deallocate(GammaL, GammaR, G_nil)
end program StypeJunction_one



subroutine input()
  use GreensFunctions 
  implicit none
  character :: ver
  
  open(2, file='input.dat', status='old')
  !...order--> 0 : non-interacting, 1: first order interactions, 2: second order interactions 
  read(2, *) T, V, mu, Volt_range, U_int
  write(*,*) 'T:', T, 'V:', V, 'mu:', mu, 'Volt_range:', Volt_range, U_int
  read(2,*) order, method
  read(2,*) Natoms
  write(*,*) 'Order:', order, 'Natoms:', Natoms
  read(2,*) dw,up,delta
  write(*,*) 'dw:', dw, 'up:', up,'delta:', delta
  read(2,*) pullay
  write(*,*) 'pulay:', pullay
  read(2,'(a)') ver
  verb=.false. ; if(ver.eq.'Y') verb=.true.
  close(2)
  beta = 1.d0/(kb * T)
  
end subroutine input
  
subroutine PrintFunctions()
  use GreensFunctions
  implicit none
  integer :: i, j
  complex*16 :: diff

  open(12, file='info_Hamiltonian.dat', status='replace')
  
  write(12,'(/a)') '... Hamiltonian:'
  do i = 1, Natoms
     write(12, '(10(a,2f10.5,a))') (' [',H(i,j),'] ', j=1,Natoms)
  end do

!.... checking if Hermitian  

  write(12,'(/a)') '... Checking the Hamiltonian is Hermitian:'
  do i=1,Natoms
     do j=i,Natoms
        diff=H(i,j)-conjg(H(j,i))
        write(12,*) i,j,diff
     end do
  end do
  
  write(12,'(/a)') '... Eigenvalues:'
  do i = 1, Natoms
     write(12, *) i,'.', Ev(i)
  end do
  
  write(12,'(/a)') '... Eigenvectors:'
  do j = 1, Natoms
     write(12, '(i3,10(a,2f10.5,a))') j,(' [',Eigenvec(i,j),'] ', i = 1, Natoms)
  end do

  close(12)
end subroutine PrintFunctions


real*8 function trans(iw, Volt) !....square bracket terms of Eq. (2) in CHE
  use GreensFunctions
  implicit none
  integer :: iw,i,j
  real*8 :: Volt, w
  complex*16 :: trace1

  w = omega(iw)
!  work3 = (0.d0, 0.d0)

  work1 = GF0%L(:,:,iw)
  work2 = GF0%G(:,:,iw)

  work3 = im*matmul(GammaL, (fermi_dist(w, Volt)-1.d0)*work1 - fermi_dist(w, Volt)*work2)
  
  call trace_of_A(work3, Natoms, trace1)
  trans = real(trace1)/(2.d0*pi)

end function trans

real*8 function Current(Volt)
  use GreensFunctions
  implicit none
  real*8 :: Volt, J_L, trans, trace !, w_in, w_f
  integer :: iw !, iw_in, iw_f
  
!......useful for calculating the non-interacting case, finds the integration window
  !call int_window(w_in, w_f, Volt) 
  !iw_in = w_grid(w_in)
  !iw_f = w_grid(w_f)

  J_L = 0.d0; trace = 0.d0
  do iw = 1, N_of_w
     trace = trans(iw, Volt)
     J_L = J_L + trace
  end do
  Current = J_L*(delta/hbar)
end function Current

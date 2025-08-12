subroutine input()
  use GreensFunctions 
  implicit none
  character :: ver
  
  open(2, file='input.dat', status='old')
  !...order--> 0 : non-interacting, 1: first order interactions, 2: second order interactions 
  read(2, *) T, V, Vf, delv, mu, U_int
 
  read(2,*) order
  read(2,*) Natoms
  read(2,*) dw,up,delta
  read(2,*) pullay
  read(2,'(a)') ver

  verb=.false. ; if(ver.eq.'Y') verb=.true.
  read(2,*) restart
  close(2)
  beta = 1.d0/(kb * T)
  if(restart) then
     write(*,*) '.... RESTART from the previous GFs ....'
  else
     write(*,*) '.... NEW start ....'
  end if

  Volt_range = (Vf - V)/delv
  
  write(*,*) 'T:', T, 'V:', V, 'mu:', mu, 'Volt_range:', Volt_range, 'U_int:', U_int
  write(*,*) 'Order:', order, 'Natoms:', Natoms
  write(*,*) 'dw:', dw, 'up:', up,'delta:', delta, 'delv:', delv
  write(*,*) 'pulay:', pullay
  
end subroutine input

subroutine print_3matrix(iteration,X,name,Natoms,N_of_w)
  integer :: iteration,n,iw,j,i,Natoms,N_of_w
  character :: fn*3,name*3
  complex*16 :: X(Natoms,Natoms,N_of_w)

  write(*,*) 'print_3 => ',iteration,Natoms,name,N_of_w
  
  n=Natoms ; if(n.gt.10) n=10
  write(fn,'(i0)') iteration
  open(7,file=name//'_'//trim(fn)//'.dat',status='unknown')
  do iw=1,N_of_w
     write(7,'(/i3,x,f10.5)') iw
     do i=1,n
        write(7,*) i,j,(X(i,j,iw),j=1,n)
     end do
  end do
  close(7)
  write(*,*) '... written '//name//'_'//trim(fn)//'.dat'
end subroutine print_3matrix

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



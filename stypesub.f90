subroutine input()
  use GreensFunctions 
  implicit none
  character :: ver
  character(len=256) :: line, key
  integer :: ios, eq_pos
  
  open(2, file='input.dat', status='old')
  !...order--> 0 : non-interacting, 1: first order interactions, 2: second order interactions 

  do
    read(2,'(A)', iostat=ios) line
    if (ios .ne. 0) exit   ! end of file
    if (trim(line) .eq. '' .or. line(1:1) .eq. '#') cycle  ! skip blank or comment

    eq_pos = index(line, "=")
    if (eq_pos .eq. 0) cycle  ! skip malformed lines

    key = adjustl(trim(line(:eq_pos-1)))

    select case (trim(key))
    case("T");            read(line(eq_pos+1:),*) T
    case("Natoms");       read(line(eq_pos+1:),*) Natoms
    case("V");            read(line(eq_pos+1:),*) V
    case("Vf");           read(line(eq_pos+1:),*) Vf
    case("delv");         read(line(eq_pos+1:),*) delv
       
    case("mu");           read(line(eq_pos+1:),*) mu
    case("order");        read(line(eq_pos+1:),*) order
    case("U_int");        read(line(eq_pos+1:),*) U_int
    case("dw");           read(line(eq_pos+1:),*) dw
    case("up");           read(line(eq_pos+1:),*) up
    case("delta");        read(line(eq_pos+1:),*) delta

    case("pulay");        read(line(eq_pos+1:),*) pullay
    case("ver");          read(line(eq_pos+1:),*) ver
    case("restart");      read(line(eq_pos+1:),*) restart

    case("E_CC");         read(line(eq_pos+1:),*) E_CC 
    case("t_hop");        read(line(eq_pos+1:),*) t_hop
    end select
  end do

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



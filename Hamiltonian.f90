module DefineHamiltonian
  implicit none
  integer :: Natoms, N_of_w, Volt_range, order, method

  real*8, dimension(:,:) :: C(3,3)
  real*8, dimension(3) :: Rij, w0
  real*8 :: T, mu, beta, V, Vf, delta, w_init, w_fin,up,dw,delv
  real*8 :: E_CC, t_hop, lamb, U_int

  real*8, parameter :: hbar = 1.d0 !6.582119569e-16 !ev s
  real*8, parameter :: kb = 8.6173303e-5 !ev/K
  real*8, parameter :: pi = 4.d0*atan(1.d0)
  real*8, parameter :: el = 1.d0 !.1.602176634e-19 !C
  
  complex*16, parameter :: im = (0.d0, 1.d0)
  complex*16, dimension(:,:) :: Pauli_x(2,2), Pauli_y(2,2), Pauli_z(2,2)
  complex*16, allocatable, dimension(:,:) :: H

  logical :: restart
contains 
 
  subroutine Central_Hamiltonian(HCC, n, ECC, thop)
    implicit none 
    integer :: i, j, n
    real*8 :: ECC, thop
    complex*16,  dimension(n,n) :: HCC
    
    HCC = 0.d0
    do i = 1, n
       HCC(i,i) = ECC
       do j = 1, n
          if(i .eq. j + 1 .or. j .eq. i + 1) then
             HCC(i,j) = thop
          end if
       end do
    end do
    
  end subroutine Central_Hamiltonian

  subroutine Island_Hamiltonian(H,N)
    implicit none
    integer :: i, j, N, N_tip, N_surface, connect
    complex*16, dimension(N,N) :: H
    real*8 :: t_tip, E_tip, t0, t_surface, E_surface

    open(12, file='Island_Hamiltonian.dat', status='old')

    read(12,*) N_tip, N_surface
    read(12,*) E_tip, E_surface
    read(12,*) t_tip, t0, t_surface
    read(12,*) connect

    close(12)

    H = (0.d0, 0.d0)
    !....Triangle tip
    do i = 1, N_tip
       do j = 1, N_tip
          if (i .eq. j) then 
             H(i,j) = E_tip
          else
             H(i,j) = t_tip
          end if
       end do
    end do

    !....bridge atom (tip apex)
    H(N_tip, connect) = t0; H(connect, N_tip) = t0

    !....Surface
    do i = N_tip+1, N_surface+N_tip
       do j = N_tip+1, N_surface+N_tip

          if (i .eq. j) then
             H(i,j) = E_surface
          else if (i .eq. j + 1 .or. j .eq. i + 1) then
             H(i,j) = t_surface
          end if
       end do
    end do
    
  end subroutine Island_Hamiltonian

  
end module DefineHamiltonian

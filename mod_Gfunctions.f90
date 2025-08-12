module GreensFunctions
  use DefineHamiltonian
  use OMP_LIB
  implicit none
  integer :: INFO

  real*8, parameter :: epsilon = 1e-9
  real*8, allocatable, dimension(:) ::  Hub, Ev
  real*8, allocatable, dimension(:) :: omega
  real*8 :: pullay
  
  complex*16, allocatable, dimension(:,:) :: GammaL, GammaR, Eigenvec
  complex*16, allocatable, dimension(:) :: G_nil
  complex*16, allocatable, dimension(:,:) ::  SigmaL, SigmaR
  complex*16, allocatable, dimension(:,:) :: work1, work2, work3, work4
  
  complex*16  :: prodr, proda, prodL, prodG
  complex*16  :: IntL, IntG, OmCon, OmCon1, Om_Con_Inf

  logical :: verb
  
  type :: GF
     complex*16, allocatable, dimension(:,:,:) :: r, a, L, G
  end type GF
  type(GF) :: GF0
  
  type :: GF_full
     complex*16, allocatable, dimension(:,:,:) :: R, A, L, G
  end type GF_full
  type(GF_full) :: GFf

CONTAINS 
  
   !.....fermi-dirac distribution
  real*8 function fermi_dist(w, V)
    implicit none
    real*8 ::  arg, w, V

    arg = (w - mu + V/hbar)*beta
    fermi_dist =  1.d0/(exp(hbar*arg) +1.d0)
    
  end function fermi_dist

!====================================================
!========== Self consistency field calculations =====
!====================================================

!.............need G0f%L to calculate GL_of_0
!.............Calculates GFs at every omega and Voltage simultaneously 
subroutine SCF_GFs(Volt,first)
  implicit none
  integer :: iw, iteration,i, Vname,wheel,it
  real*8 :: Volt, err, diff, st, et,current
  character(len=30) :: fn, fn1
  logical :: first

 1 iteration = 0

  if (first) then 
     call G0_R_A()
     call G0_L_G(Volt)
  end if

  print *, '>>>>>>>>>>VOLTAGE:', Volt

  DO
     iteration = iteration + 1
     write(*,*) '.... ITERATION = ',iteration,' ....'

!.......real variable interactions turns off the Interaction component of the sigmas 
!.................full Gr and Ga, Eq. (5) and (6)

     call GL_of_0()

     !$OMP PARALLEL DO &
     !$OMP& PRIVATE(iw, INFO)
     
     do iw = 1, N_of_w
        call G_full(iw, Volt)
     end do
     !$OMP END PARALLEL DO

        GF0%r = pullay*GFf%r + (1.0d0-pullay)*GF0%r
        GF0%l = pullay*GFf%l + (1.0d0-pullay)*GF0%l
!...... do the advanced and greater components
     
     do iw=1,N_of_w
        work1=GF0%r(:,:,iw) 
        call Hermitian_Conjg(work1, Natoms, work2)
        GF0%a(:,:,iw)=work2
     end do
     GF0%G = GF0%L + GF0%R - GF0%A

!..... calculation of the error     

     err=0.0d0
     do iw = 1, N_of_w
        do i=1,Natoms
           diff=2.d0*hbar*(AIMAG(GFf%R(i,i,iw))-AIMAG(GF0%R(i,i,iw)))
           err=err +diff*diff
        end do
     end do
     write(*,*) 'err = ',sqrt(err)

!____________ useful if one would like to check the convergence     
!     write(13,*) iteration, current(Volt)
     
     if (sqrt(err) .lt. epsilon .or. order .eq. 0) then
        write(*,*)'... REACHED REQUIRED ACCURACY ...'
        exit
     end if
     
  END DO
  
end subroutine SCF_GFs
  
!=====================================================
!================== Full GFs =========================
!===================================================== 

  subroutine G_full(iw, Volt)
!... Full Greens function, leaves Retarded and Advanced in the work arrays, application
!    of Eq. (16) and (17), but with the full Sigmas, Eq. (3), (7) and (8) in CHE
   implicit none
   integer :: i, j, iw
   real*8 :: Volt, w 
   complex*16 :: Omr, SigG, SigL
   complex*16, allocatable, dimension(:,:) ::  SigmaL, SigmaR, SigmaG,  work_1, work_2
   
   allocate(SigmaL(Natoms, Natoms), SigmaG(Natoms, Natoms), SigmaR(Natoms, Natoms))
   allocate(work_1(Natoms, Natoms), work_2(Natoms, Natoms))
!............full SigmaR due to interactions Eq. (7)
   SigmaR = (0.d0, 0.d0); SigmaL = (0.d0, 0.d0)
      
   if (order .eq. 1) then
      do i = 1, Natoms
         SigmaR(i,i) = SigmaR(i,i) - im*hbar*Hub(i)*G_nil(i) ! LK <== must be minus!
      end do
  
   else if (order .eq. 2) then 
      do i = 1, Natoms
         do j = 1, Natoms
            Omr = Omega_R(i,j,iw)
            if (i .eq. j) then 
               SigmaR(i,j) =  - im*hbar*Hub(i)*G_nil(i) + (Hub(j)*Hub(i)*OmR)*hbar**2
            else 
               SigmaR(i,j) = (Hub(j)*Hub(i)*OmR)*hbar**2
            end if
         end do
      end do
      
!..............full SigmaL and SigmaG, Eq. (3) and (4)     
!.....Interaction contribution of both Sigmas     
      do i = 1, Natoms
         do j = 1, Natoms
            call Omega_int_SigL_SigG(i,j, iw, SigL, SigG)  
!            SigmaG(i,j) =  Hub(i)*Hub(j)*SigG*hbar**2 
            SigmaL(i,j) =  Hub(i)*Hub(j)*SigL*hbar**2
         end do
      end do
   end if

!..real variable interactions turns off the Interaction component of the sigmas 
!............full Gr and Ga, Eq. (5) and (6)
   
  w = omega(iw)
  work_1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR 
  do i = 1 , Natoms
     work_1(i,i) = work_1(i,i) + hbar*(w+im*0.01)
  end do
  
  call Inverse_complex(Natoms, work_1, info)
  call Hermitian_Conjg(work_1, Natoms, work_2)
  
  GFf%R(:,:,iw) = work_1; GFf%A(:,:,iw) = work_2

!.....Embedding contribution of both Sigmas

  SigmaL =  im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)/hbar + SigmaL 
!  SigmaG =  im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)/hbar + SigmaG

  !.............full GL and GG, Eq. (16) and (17) 
  
  GFf%L(:,:,iw) = matmul(matmul(work_1, SigmaL), work_2) 
!  GFf%G(:,:,iw) = matmul(matmul(work_1, SigmaG), work_2)

  deallocate(SigmaL,SigmaG,SigmaR); deallocate(work_1, work_2)
end subroutine G_full

!=====================================================
!======== Non-interacting GFs ========================
!=====================================================  

  subroutine G0_R_A()
!............non-interacting Greens functions: GR and Ga,  Eq. (5) and Eq. (6) in 'Current_Hubbard_Equations' document (CHE)
    implicit none
    integer :: j, i
    real*8 :: w
   
    do j = 1, N_of_w
       work1 = -H + 0.5d0*(im/hbar) * (GammaL + GammaR) !LK <========= must be +
       w = omega(j)
       do i = 1 , Natoms
          work1(i,i) = work1(i,i) + hbar*(w+im*0.01)
       end do
       
       call Inverse_complex(Natoms, work1, info)
       call Hermitian_Conjg(work1, Natoms, work2)
       
       GF0%r(:,:,j) = work1 ; GF0%a(:,:,j) = work2
    end do  
    
  end subroutine G0_R_A
  
  subroutine G0_L_G(Volt)
!............non-interacting Greens functions: G> and G< for all omega on the grid, Eq. (16) and (17) in CHE
    implicit none
    real*8 :: Volt, w
    integer :: j, i
    
    work1 = (0.d0, 0.d0) ; work2 =(0.d0, 0.d0) ; work3 = (0.d0, 0.d0); work4 = (0.d0, 0.d0)
    do j = 1 , N_of_w
       w = omega(j)
       work1 = GF0%r(:,:,j) 
       work2 = GF0%a(:,:,j)

       work3 = matmul(matmul(work1, im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)), work2) 
!      work4 = matmul(matmul(work1, im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)), work2)
      GF0%L(:,:,j) = work3
!       GF0%G(:,:,j) = work4
    end do
    GF0%G =  GF0%L + GF0%R - GF0%A
    
  end subroutine G0_L_G 

!=====================================================
!========Calcualtions needed for full GFs=============
!===================================================== 
  
!......................Calculation of Omega terms for the self-energies, Eq. (9) in CHE
  complex*16 function Omega_R(i, j, iw)
    implicit none
    integer :: i, j, iw, k_1, k_2, k_3
    complex*16 :: Omr
    real*8 :: pp
    
    Omr = (0.d0, 0.d0)
    do k_1 = 1, N_of_w
       do k_2 = 1, N_of_w
          
          k_3 = iw- k_1 +k_2
          if (k_3 .ge. 1 .and. k_3 .le. N_of_w) then
             Omr = Omr + GF0%r(i,j,k_1)*GF0%L(j,i,k_2)*GF0%G(i,j,k_3)  &
                  + GF0%L(i,j,k_1)*GF0%L(j,i,k_2)*GF0%r(i,j,k_3) &
                  + GF0%L(i,j,k_1)*GF0%a(j,i,k_2)*GF0%L(i,j,k_3)
          end if

       end do
    end do

    pp = delta/(2.d0*pi) 
    Omega_R = Omr*pp*pp
  end function Omega_R
  
  subroutine Omega_int_SigL_SigG(i,j, iw, SigL, SigG) !... interaction contributions of Eq. (3) and (4) in CHE
    implicit none
    integer :: i, j, k1, k2, k3, iw
    complex*16 :: SigL, SigG
    real*8 :: pp

    
    pp = delta/(2.d0*pi) 
    SigL=(0.0d0, 0.0d0); SigG = (0.d0, 0.d0) 
    do k1=1,N_of_w
       do k2=1,N_of_w        
          k3=iw-k1+k2
          if(k3 .ge. 1 .and. k3 .le. N_of_w) then
             SigL = SigL + GF0%L(i,j,k1)*GF0%G(j,i,k2)*GF0%L(i,j,k3) 
!             SigG = SigG + GF0%G(i,j,k1)*GF0%L(j,i,k2)*GF0%G(i,j,k3) 
          end if
       end do
    end do
    
    SigL = SigL*pp*pp !; SigG = SigG*pp*pp
  end subroutine Omega_int_SigL_SigG

!......Calculates lesser Greens function at time = 0
  subroutine GL_of_0()        !LK <======= a slight change
    !.....Lesser Greens function at time = 0
    !.....ei - eigenvalues of the Hamiltonian 
    implicit none
    integer :: i, k1
    real*8 :: pp
    complex*16 :: s
    
    pp=delta/(2.d0*pi)
    do i = 1, Natoms 
       s =(0.d0, 0.d0)
       do k1 = 1, N_of_w
          s = s + GF0%L(i,i,k1)
       end do
       G_nil(i)=s*pp
    end do

  end subroutine GL_of_0
  
!====================================================
!========== calcualtion of the current ==============
!====================================================

subroutine trans(iw, Volt, trns) 
  implicit none
  integer :: iw,i,j
  real*8 :: Volt, w,trns
  complex*16 :: trace1
  complex*16, allocatable, dimension(:,:) :: work_1, work_2,work_3

  allocate(work_1(Natoms, Natoms),work_2(Natoms, Natoms),work_3(Natoms, Natoms))
  
  w = omega(iw)
  work_3 = (0.d0, 0.d0) ; work_1 = GF0%L(:,:,iw) ; work_2 = GF0%G(:,:,iw)

  work_3 = im*matmul(GammaL,(fermi_dist(w, Volt)-1.d0)*work_1 - fermi_dist(w, Volt)*work_2)
  
  call trace_of_A(work_3, Natoms, trace1)
  trns = real(trace1)/(2.d0*pi)

  deallocate(work_1,work_2,work_3)
  
end subroutine trans

!====================================================
!========== Printing Routines ================= =====
!====================================================

subroutine print_sf(iteration)
  integer :: iteration,n,iw,j
  character :: fn*3 
  
  n=Natoms ; if(n.gt.10) n=10
  write(fn,'(i0)') iteration
  open(7,file='sf_'//trim(fn)//'.dat',status='unknown')
  do iw=1,N_of_w
     write(7,'(f10.5,x,10(e12.5,x))') omega(iw),(- 2.d0*hbar*AIMAG(GF0%R(j,j,iw)),j=1,n)
  end do
  close(7)
  write(*,*) '... written sf_'//trim(fn)//'.dat'
end subroutine print_sf

!====================================================
!========== restart routines  =======================
!====================================================

subroutine save_GFs()
  implicit none
  open(1,file='Greens_functions.dat',form='unformatted',status='unknown')
  write(1) GF0%r,GF0%a,GF0%l,GF0%g
  close(1)
end subroutine save_GFs

subroutine read_saved_GFs()
  implicit none
  open(1,file='Greens_functions.dat',form='unformatted',status='old')
  read(1) GF0%r,GF0%a,GF0%l,GF0%g
  close(1)
end subroutine read_saved_GFs

end module GreensFunctions

module GreensFunctions
  use DefineHamiltonian
  use OMP_LIB
  implicit none
  integer :: INFO

  real*8, parameter :: epsilon = 1e-6
  real*8, allocatable, dimension(:) ::  Hub, Ev
  real*8, allocatable, dimension(:) :: omega
  real*8 :: pullay
  
  complex*16, allocatable, dimension(:,:) :: GammaL, GammaR, Eigenvec
  complex*16, allocatable, dimension(:) :: G_nil
!  complex*16, allocatable, dimension(:,:) ::  SigmaL, SigmaG, SigmaR
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
  integer :: iw, iteration,i, Vname
  real*8 :: Volt, err, diff, st, et
  character(len=30) :: fn, fn1
  logical :: first

  iteration = 0

  write(22,*) '........SCF Calculations at Voltage:', Volt, '..........'

  if (first) then 
     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call G0_R_A()
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'G0_R_A runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
     
     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call G0_L_G(Volt)
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'G0_L_G runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
     
  end if

  print *, '>>>>>>>>>>VOLTAGE:', Volt

  if(verb) then
     call print_3matrix(0,GF0%R,'GFR')
     call print_3matrix(0,GF0%L,'GFL')
  end if
  
! LK........ printing the spectral function at 0 iteration (embedding only) 
!            and calculating spec(1,iw)

  call print_sf(iteration)
  
  Vname = abs(Volt)
  write(fn1,'(i0)') Vname
  if (Volt .ge. 0) then 
     open(17,file='err_V_'//trim(fn1)//'.dat',status='unknown')
  else
     open(17,file='err_V_n'//trim(fn1)//'.dat',status='unknown')
  end if

  DO
     iteration = iteration + 1
     write(*,*) '.... ITERATION = ',iteration,' ....'

!.......real variable interactions turns off the Interaction component of the sigmas 
!.................full Gr and Ga, Eq. (5) and (6)

     st =0.d0; et= 0.d0
     call CPU_TIME(st)
     call GL_of_0()
     call CPU_TIME(et)
     write(22,'(A,F10.8,A,A,A,I4)') 'GL_of_0 runtime:', (et-st), 'seconds', '   ','Iteration:', iteration
     
     st = 0.d0; et = 0.d0
     st = OMP_GET_WTIME()

     !$OMP PARALLEL DO &
     !$OMP& PRIVATE(iw, INFO)
     
     do iw = 1, N_of_w
        call G_full(iw, Volt)
     end do
     !$OMP END PARALLEL DO
     
     et = OMP_GET_WTIME()
     write(22,'(A,F10.8,A,A,A,I4)') 'G_full w-loop runtime:', (et-st), 'seconds', '   ', 'Iteration:', iteration
     
     
     if(verb) then
        call print_3matrix(iteration,GFf%R,'GFR')
        call print_3matrix(iteration,GFf%L,'GFL')
     end if
     
!..... calculation of the error     

     err=0.0d0
     !$OMP PARALLEL DO PRIVATE(iw, i, diff) REDUCTION(+:err)
     do iw = 1, N_of_w
        do i=1,Natoms
           diff=2.d0*hbar*(AIMAG(GFf%R(i,i,iw))-AIMAG(GF0%R(i,i,iw)))
           err=err +diff*diff
        end do
     end do
     !$OMP END PARALLEL DO
     write(*,*) 'err = ',sqrt(err)
     write(17,*) iteration, sqrt(err)
 
     !$OMP CRITICAL     
     GF0%R = pullay*GFf%R + (1.0d0-pullay)*GF0%R
     GF0%A = pullay*GFf%A + (1.0d0-pullay)*GF0%A
     GF0%L = pullay*GFf%L + (1.0d0-pullay)*GF0%L
     GF0%G = pullay*GFf%G + (1.0d0-pullay)*GF0%G
     !$OMP END CRITICAL
     
     call print_sf(iteration)  
     if (sqrt(err) .lt. epsilon .or. order .eq. 0) then
        write(*,*)'... REACHED REQUIRED ACCURACY ...'
        exit
     end if
  END DO
  close(17)
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
            !SigmaG(i,j) =  Hub(i)*Hub(j)*SigG*hbar**2 
            SigmaL(i,j) =  Hub(i)*Hub(j)*SigL*hbar**2
         end do
      end do
   end if

   if(verb) then
      write(3,*) iw,' SigmaR'
      do i=1,Natoms
         write(3,*) i, (SigmaR(i,j),j=1,Natoms)
      end do
      write(3,*) iw,' SigmaL'
      do i=1,Natoms
         write(3,*) i, (SigmaL(i,j),j=1,Natoms)
      end do
      write(3,*) iw,' SigmaG'
   end if
!..real variable interactions turns off the Interaction component of the sigmas 
!............full Gr and Ga, Eq. (5) and (6)
   
  w = omega(iw)
  work_1 = -H + 0.5d0*(im/hbar)*(GammaL + GammaR) - SigmaR ! LK <== must be + for emb and minus for interaction sigma
  do i = 1 , Natoms
     work_1(i,i) = work_1(i,i) + hbar*(w+im*0.01)
  end do
  
  call Inverse_complex(Natoms, work_1, info)
  call Hermitian_Conjg(work_1, Natoms, work_2)
  
  GFf%R(:,:,iw) = work_1; GFf%A(:,:,iw) = work_2

!.....Embedding contribution of both Sigmas

  SigmaL =  im*(fermi_dist(w, Volt)*GammaL + fermi_dist(w, 0.d0)*GammaR)/hbar + SigmaL 
  SigmaG =  im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)/hbar + SigmaG

  !.............full GL and GG, Eq. (16) and (17)
  
  GFf%L(:,:,iw) = matmul(matmul(work_1, SigmaL), work_2) !.. GL = Gr * SigmaL * Ga
  GFf%G(:,:,iw) = matmul(matmul(work_1, SigmaG), work_2) !.. GG = Gr * SigmaG * Ga
  !GFf%G(:,:,iw) = GFf%L(:,:,iw) + GFf%R(:,:,iw) - GFf%A(:,:,iw)

  deallocate(SigmaL,SigmaG,SigmaR); deallocate(work_1, work_2)
end subroutine G_full


!=====================================================
!======== Non-interacting GFs ========================
!=====================================================  

  
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
  
  subroutine G0_R_A()
    !............non-interacting Greens functions: GR and Ga,  Eq. (5) and Eq. (6) in 'Current_Hubbard_Equations' document (CHE)
    implicit none
    integer :: j, i
    real*8 :: w
   
    do j = 1, N_of_w
       work1 = -H + 0.5d0*(im/hbar) * (GammaL + GammaR) !LK <========= must be +
       w = omega(j)
       do i = 1 , Natoms
          work1(i,i) = work1(i,i) + hbar*(w+ im*0.01)
       end do
       
       call Inverse_complex(Natoms, work1, info)
       call Hermitian_Conjg(work1, Natoms, work2)
       
       GF0%r(:,:,j) = work1
       GF0%a(:,:,j) = work2
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
       work4 = matmul(matmul(work1, im*((fermi_dist(w, Volt)-1.d0)*GammaL + (fermi_dist(w, 0.d0)-1.d0)*GammaR)), work2)
       GF0%L(:,:,j) = work3
       GF0%G(:,:,j) = work4
       !GF0%G(:,:,j) =  GF0%L(:,:,j) + GF0%R(:,:,j) - GF0%A(:,:,j)
    end do
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
             Omr = Omr + GF0%r(i,j,k_1)*GF0%L(j,i,k_2)*GF0%L(i,j,k_3)  + GF0%L(i,j,k_1)*GF0%a(j,i,k_2)*GF0%L(i,j,k_3) &
                  + GF0%L(i,j,k_1)*GF0%G(j,i,k_2)*GF0%a(i,j,k_3)
          end if

       end do
    end do

    pp = delta/(2.d0*pi) !LK <== error in brackets
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
             SigG = SigG + GF0%G(i,j,k1)*GF0%L(j,i,k2)*GF0%G(i,j,k3) 
          end if
       end do
    end do
    
    SigL = SigL*pp*pp ; SigG = SigG*pp*pp
  end subroutine Omega_int_SigL_SigG

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

subroutine print_3matrix(iteration,X,name)
  integer :: iteration,n,iw,j,i
  character :: fn*3,name*3
  complex*16 :: X(Natoms,Natoms,N_of_w)
   
  n=Natoms ; if(n.gt.10) n=10
  write(fn,'(i0)') iteration
  open(7,file=name//'_'//trim(fn)//'.dat',status='unknown')
  do iw=1,N_of_w
     write(7,'(/i3,x,f10.5)') iw,omega(iw)
     do i=1,n
        write(7,*) i,j,(X(i,j,iw),j=1,n)
     end do
  end do
  close(7)
  write(*,*) '... written '//name//'_'//trim(fn)//'.dat'
end subroutine print_3matrix

end module GreensFunctions

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

  !........ the wheel memory for Pulay
  integer :: iP,size
  type :: GX
     complex*16, allocatable, dimension(:,:,:,:) :: r_in,r_out,l_out
  end type GX
  type(GX) :: GP
  real*8,allocatable,dimension(:,:) :: Ov,Bm,Rhs
  real*8,allocatable,dimension(:) :: C_coeff
  integer,dimension(:),allocatable :: IPIV


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

  iteration = 0

  if (first) then 
     call G0_R_A()
     call G0_L_G(Volt)
  end if

  print *, '>>>>>>>>>>VOLTAGE:', Volt

!.... initialise the wheel: this is position on the wheel for the curent GFs
!                           to be placed
  GP%r_in=(0.0d0,0.0d0) ; GP%r_out=(0.0d0,0.0d0) ; Ov=0.0d0
  
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

!
!==========================............... one step Pulay.........===========
!
     if(method.eq.1) then

        GF0%r = pullay*GFf%r + (1.0d0-pullay)*GF0%r
        GF0%l = pullay*GFf%l + (1.0d0-pullay)*GF0%l

!
!==========================............. many steps Pulay.........===========
!
     else if(method.eq.2) then
        
!....... place the GFf%r_in/out into the wheel
     
        IF(iteration.eq.1) THEN
        
           GP%r_in(:,:,:,1)=GF0%r
           GF0%r = pullay*GFf%r + (1.0d0-pullay)*GF0%r
           GF0%l = pullay*GFf%l + (1.0d0-pullay)*GF0%l
           GP%r_out(:,:,:,1)=GF0%r ; GP%l_out(:,:,:,1)=GF0%l
           
           size=1
           call overlap(1,1)
           
        ELSE IF(iteration.eq.2) THEN
           
           GP%r_in(:,:,:,2)=GFf%r
           GF0%r = pullay*GFf%r + (1.0d0-pullay)*GF0%r
           GF0%l = pullay*GFf%l + (1.0d0-pullay)*GF0%l
           GP%r_out(:,:,:,2)=GF0%r ; GP%l_out(:,:,:,2)=GF0%l
           
           size=2
           call overlap(2,2) ; call overlap(2,1) ; Ov(1,2)=Ov(2,1)
           
        ELSE IF(iteration.le.iP) THEN

!_________ get c-coefficients from the previous Ov matrix        

           call C_coefficients()

!_________ update GF0%r and GF0%l (out) by mixing with the previos iterations

           GP%r_in(:,:,:,iteration)=GFf%r
           GF0%r=pullay*GFf%r ; GF0%l=pullay*GFf%l
           do it=1,iteration-1
              GF0%r=GF0%r+(1-pullay)*C_coeff(it)*GP%r_out(:,:,:,it)
              GF0%l=GF0%l+(1-pullay)*C_coeff(it)*GP%l_out(:,:,:,it)
           end do
           GP%r_out(:,:,:,iteration)=GF0%r ; GP%l_out(:,:,:,iteration)=GF0%l
        
!_________ updting the ov matrix for the next iteration       

           size=iteration
           call overlap(iteration,iteration) 
           do i=1,iteration-1
              call overlap(iteration,i) ; Ov(i,iteration)=Ov(iteration,i)
           end do
        
        ELSE IF(iteration.gt.iP) THEN
           
           call C_coefficients()

           wheel=mod(iteration,iP) ; if(wheel.eq.0) wheel=iP
        
!_________ update GF0%r and GF0%l (out) by mixing with the previos iterations

           GP%r_in(:,:,:,wheel)=GFf%r
           GF0%r=pullay*GFf%r ; GF0%l=pullay*GFf%l
           do it=1,iP
              GF0%r=GF0%r+(1-pullay)*C_coeff(it)*GP%r_out(:,:,:,it)
              GF0%l=GF0%l+(1-pullay)*C_coeff(it)*GP%l_out(:,:,:,it)
           end do
           GP%r_out(:,:,:,wheel)=GF0%r ; GP%l_out(:,:,:,wheel)=GF0%l
        
!_________ updating the ov matrix for the next iteration       

           size=iP
           call overlap(wheel,wheel)
           do i=1,iP
              if(i.ne.wheel) then
                 call overlap(wheel,i) ; Ov(i,wheel)=Ov(wheel,i)
              end if
           end do
        
        END IF

     end if
     
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
  
!  write(*,*) 'exiting the convergence loop'
!  call print_3matrix(100,GF0%g,'GFg',Natoms,N_of_w)

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
            SigmaG(i,j) =  Hub(i)*Hub(j)*SigG*hbar**2 
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
 ! GFf%G(:,:,iw) = matmul(matmul(work_1, SigmaG), work_2) !.. GG = Gr * SigmaG * Ga
  GFf%G(:,:,iw) = GFf%L(:,:,iw) + GFf%R(:,:,iw) - GFf%A(:,:,iw)

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


!====================================================
!========== Pulay's Routines ========================
!====================================================

subroutine overlap(i,j)
  integer :: i,j,i1,iw
  real*8 :: ci,cj,oo,pp

  oo=(0.0d0,0.0d0)
  do i1=1,Natoms
     do iw=1,N_of_w
        ci=aimag(GP%r_in(i1,i1,iw,i)-GP%r_out(i1,i1,iw,i))
        cj=aimag(GP%r_in(i1,i1,iw,j)-GP%r_out(i1,i1,iw,j))
        oo = oo + ci*cj
     end do
  end do
  pp = delta/(2.d0*pi) 
  Ov(i,j)=oo*pp
  
end subroutine overlap

subroutine C_coefficients()
  integer :: i,j,N,INFO
  real*8 :: s(10)

  N=size+1
  allocate(Bm(N,N))

  Bm(1:size,1:size)=Ov ; Bm(N,N)=0.0d0
  do i=1,size
     Bm(i,N)=-1.0d0 ;  Bm(N,i)=-1.0d0
  end do
  Rhs(1:size,1)=0.0d0 ; Rhs(N,1)=-1.0d0

!...... calling linear system of eqs: Bm * X = Rhs, 
!       where on output X is in Rhs and Bm destroyed

  call dgesv(N,1,Bm,N,IPIV,Rhs,N,INFO)

!.... solution

  C_coeff(1:size)=Rhs(1:size,1)
  deallocate(Bm)

end subroutine C_coefficients

end module GreensFunctions

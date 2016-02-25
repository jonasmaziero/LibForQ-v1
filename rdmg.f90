!###################################################################################################################################
!                                             Random density matrix generators - RDMG
!###################################################################################################################################
subroutine rdm_std(d, rdm) ! Generates a random density matrix using the standard method: rdm = \sum_j p_j U|c_j><c_j|U^†
! Ref: Maziero, J. (2015). Random sampling of quantum states: A survey of methods. Braz. J. Phys. 45, 575.
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: ru(:,:)  ! Random unitary matrix
real(8), allocatable :: rpv(:)  ! Random probability vector
integer :: j, k, l  ! Auxiliary variable for counters

allocate ( ru(1:d,1:d), rpv(1:d) ) ;   call rpvg(d, rpv) ;   call rug(d, ru)  ! Allocate memory for and get these random variables

rdm = (0.d0,0.d0)  ! Generates the rdm
do j = 1, d ;   do k = 1, d  ;  do l = 1, d ;   rdm(j,k) = rdm(j,k) + rpv(l)*ru(j,l)*conjg(ru(k,l)) ;   enddo ;   enddo ;   enddo

deallocate ( ru, rpv )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_ginibre(d, rdm) ! Generates a random density matrix normalizing G*G^†, with G being a Ginibre matrix
! Ref:  \.{Z}yczkowski, K., and Sommers, H.-J. (2001). Induced measures in the space of mixed quantum states. 
!       J. Phys. A: Math. Gen. 34, 7111.
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: G(:,:), GGd(:,:)  ! For the Ginibre matrix and its product with its adjoint
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm function

allocate ( G(1:d,1:d), GGd(1:d,1:d) ) ;   call ginibre(d, d, G) ;   call  matmul_AAd(d, d, G, GGd) 
rdm = GGd/((norm_hs(d, d, G))**2.d0)  ! Defines the density matrix
deallocate ( G, GGd )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_bures(d, rdm) ! Generates a random density matrix normalizing (id+U)G*G^†(id+U^†), with G being a Ginibre matrix and U a random unitary
! Ref: Al Osipov, V., Sommers, H.-J., and \.{Z}yczkowski, K. (2010). Random Bures mixed states and the distribution of their purity. 
!      J. Phys. A: Math. Theor. 43, 055302.
implicit none
integer :: d  ! Dimension of the matrices
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: G(:,:)  ! For the Ginibre matrix
complex(8), allocatable :: U(:,:)  ! For the random unitary matrix
complex(8), allocatable :: id(:,:)  ! For the indentity matrix
complex(8), allocatable :: A(:,:), AAd(:,:)  ! Auxiliary matrices
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm function

allocate ( G(1:d,1:d), U(1:d,1:d), id(1:d,1:d), A(1:d,1:d), AAd(1:d,1:d) )
call ginibre(d, d, G) ;   call rug(d, U) ;   call identity_c(d, id) ;   A = matmul((id+U),G) ;   call  matmul_AAd(d, d, A, AAd)
rdm = AAd/((norm_hs(d, d, A))**2.d0)  ! Defines the density matrix
deallocate ( G, U, id, A, AAd )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_ptrace(d, rdm) ! Generates a random density matrix via partial tracing over a random state vector
! Ref: Mej\'ia, J., Zapata, C., and Botero, A. (2015). The difference between two random mixed quantum states: Exact and asymptotic 
!      spectral analysis. arXiv:1511.07278.
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: rsv(:)  ! For the random state vector
complex(8), allocatable :: proj(:,:)  ! For the projector

allocate ( rsv(1:d*d), proj(1:d*d,1:d*d) )
call rsvg(d*d, rsv) ;   call projector(rsv, d*d, proj) ;   call partial_trace_a(d, d, proj, rdm)
deallocate ( rsv, proj )

end
!###################################################################################################################################
!                                                   Calling subroutines - RDMG
!###################################################################################################################################
subroutine rdmg(d, rdm)  ! Calls the choosed random density matrix generator
use meths
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! The random density matrix

     if ( opt_rdmg == "std" ) then ;   call rdm_std(d, rdm)
else if ( opt_rdmg == "ginibre" ) then ;   call rdm_ginibre(d, rdm)
else if ( opt_rdmg == "bures" ) then ;   call rdm_bures(d, rdm)
else if ( opt_rdmg == "ptrace" ) then ;   call rdm_ptrace(d, rdm)
     endif

end
!###################################################################################################################################
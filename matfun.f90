!###################################################################################################################################
!                                                    Auxiliary Matrices
!###################################################################################################################################
subroutine ginibre(m, n, G)  ! Returns a m x n complex matrix from the Ginibre ensemble
implicit none
integer :: m, n  ! No. of rows and columns of G
complex(8) :: G(1:m,1:n)  ! The Ginibre matrix
real(8), allocatable :: grn(:)  ! Vector of gaussianily distributed random numbers
integer :: j, k  ! Auxiliary variables for counters

allocate( grn(1:2*m) )
do j = 1, n
  call rng_gauss(2*m, grn) ;   forall( k = 1:m ) G(k,j) = grn(k) + (0.d0,1.d0)*grn(m+k)  ! Generates G column by column
enddo
deallocate( grn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine identity_c(d, identity)  ! Returns the complex dxd identity matrix
implicit none
integer :: d  ! Dimension of the identity matrix
complex(8) :: identity(1:d,1:d)  ! Identity matrix
integer :: j, k  ! Auxiliary variable for counters

forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = (0.d0,0.d0) ;   forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = (1.d0,0.d0)
  
end
!------------------------------------------------------------------------------------------------------------------------------------
module pauli_group  ! Defines the three Pauli's matrices (sigma_1,sigma_2,sigma_3) and the identity (sigma_0)
implicit none

! identity matrix
complex(kind=8), parameter :: sigma_0(1:2,1:2) = reshape( (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),(1.d0,0.d0) /) , (/2,2/) )
! Pauli's matrix sigma_x    
complex(kind=8), parameter :: sigma_1(1:2,1:2) = reshape( (/ (0.d0,0.d0),(1.d0,0.d0),(1.d0,0.d0),(0.d0,0.d0) /) , (/2,2/) )
! Pauli's matrix sigma_y    
complex(kind=8), parameter :: sigma_2(1:2,1:2) = reshape( (/ (0.d0,0.d0),(0.d0,1.d0),-(0.d0,1.d0),(0.d0,0.d0) /) , (/2,2/) )
! Pauli's matrix sigma_z   
complex(kind=8), parameter :: sigma_3(1:2,1:2) = reshape( (/ (1.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0),-(1.d0,0.d0) /) , (/2,2/) )   

end module
!------------------------------------------------------------------------------------------------------------------------------------
module bell_basis  ! Defines the Bell basis states
implicit none

complex(8), parameter :: psi_p(1:4) = (/ 0.d0, 1.d0/sqrt(2.d0), 1.d0/sqrt(2.d0), 0.d0 /)
complex(8), parameter :: psi_m(1:4) = (/ 0.d0, 1.d0/sqrt(2.d0), -1.d0/sqrt(2.d0), 0.d0 /)
complex(8), parameter :: phi_p(1:4) = (/ 1.d0/sqrt(2.d0), 0.d0, 0.d0, 1.d0/sqrt(2.d0) /)
complex(8), parameter :: phi_m(1:4) = (/ 1.d0/sqrt(2.d0), 0.d0, 0.d0, -1.d0/sqrt(2.d0) /)  

end module
!###################################################################################################################################
!                                                      Matrix functions
!###################################################################################################################################
subroutine rand_perm(d, rperm)  ! Returns a random permutation of {1,2,...,d-1,d}
! Used e.g. in the normalization and trigonometric methods for rpvg
implicit none
integer :: d  ! Dimension of the random permutation vector
integer :: rperm(1:d)  ! Random permutation vector
integer :: j  ! Counter for do (auxiliary variable)
integer, allocatable :: counter(:)  ! Counter for the no. of times a component of rand_perm is randomly choosed (auxiliary variable)
integer :: ind  ! Identify the component of rand_perm choosed (auxiliary variable)
real(8) :: rn(1)  ! Random numbers

allocate( counter(1:d) )  ! Allocate memory for the counter for the componets of rand_perm
 counter = 0
 
do j = 1, d
  do
    call rng(1,rn)  ! Returns one random number using the method choosed via option op_rng 
    if ( rn(1) <= 0.d0 ) rn(1) = 1.d-10  ! These two lines are a precaution for the determination of ind (avoid rn = 0 and rn = 1) 
    if ( rn(1) >= 1.d0 ) rn(1) = 1.d0 - 1.d-10
    ind = aint( dble(d)*rn(1) + 1.d0 )  ! aint(x) returns the largest integer smaller than x (i.e., ind is in [1,d])
    counter(ind) = counter(ind) + 1  ! No. of times the value ind appeared in rperm
    if ( counter(ind) >= 2 ) cycle  ! cycle returns the pointer to the top of the do
    if ( counter(ind) == 1 ) then
      rperm(j) = ind ;   exit  ! exit the do going to the next j, i.e., the next component of rand_perm
    endif
  enddo
enddo
 
deallocate( counter )  ! Deallocate the memory used with the counter

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine projector(vec, d, proj)  ! Returns a projector on the provided vector 
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: vec(1:d)  ! Vector we want the projector on
complex(8) :: proj(1:d,1:d)  ! Projector on vec
integer :: j, k  ! Auxiliary variables for counters

forall ( j=1:d, k=1:d, j <= k ) proj(j,k) = vec(j)*conjg(vec(k)) ;   forall ( j=1:d, k=1:d, j > k ) proj(j,k) = conjg(proj(k,j))

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kron_prod_pauli_mat(ord_pm, n, kp_pauli_mat)  ! Returns the Kronecker's product of n Pauli matrices. 
! The especification of the order of the matrices is given as input within the vector vec_ord_pm, whose dimension is n
use pauli_group
implicit none
integer, intent(in):: n
integer, intent(in):: ord_pm(1:n)
complex(kind=8), intent(out) :: kp_pauli_mat(1:2**n,1:2**n)
integer, parameter :: d2 = 2
complex(kind=8) :: M2(1:d2,1:d2)
complex(kind=8), allocatable :: M1(:,:), M1_kp_M2(:,:)
integer :: d1, i

do i = 1, n - 1
       if (ord_pm(i+1) == 0) then ;   M2 = sigma_0 ;   else if (ord_pm(i+1) == 1) then ;   M2 = sigma_1
  else if (ord_pm(i+1) == 2) then ;   M2 = sigma_2 ;   else if (ord_pm(i+1) == 3) then ;   M2 = sigma_3
       end if
  d1 = 2**i
  if( i == 1 ) then               ! Used to initiate the sequence of tensor products
    allocate(M1(1:d1,1:d1))
         if (ord_pm(i) == 0) then ;   M1 = sigma_0 ;   else if (ord_pm(i) == 1) then ;   M1 = sigma_1
    else if (ord_pm(i) == 2) then ;   M1 = sigma_2 ;   else if (ord_pm(i) == 3) then ;   M1 = sigma_3
         end if
  end if
  allocate(M1_kp_M2(1:d1*d2,1:d1*d2)) ;   call kronecker_product_C(M1, d1, d1, M2, d2, d2, M1_kp_M2)  
  deallocate(M1) ;   allocate(M1(1:d1*d2,1:d1*d2)) ;   M1 = M1_kp_M2 ;   deallocate(M1_kp_M2)
end do

 kp_pauli_mat = M1 ;   deallocate(M1)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_id_bra(m, n, ket, kp)  ! Returns the tensor product of the mxm identity with the bra (adjoint) of the ket
implicit none
integer :: m  ! Dimension of the identity (given as input)
complex(8) :: id(1:m,1:m)  ! Identify
integer :: n  ! Dimension of the vector
complex(8) :: ket(1:n), bra(1,1:n)  ! Vector and its adjoint (we use a 2d array for bra for practical issues)
complex(8) :: kp(1:m,1:m*n)  ! The Kronecker product
integer :: j  ! Auxiliary variable for counters

forall(j=1:n) bra(1,j) = conjg(ket(j))
kp = 0.d0 ;   forall(j=1:m) kp(j,((j-1)*n+1):j*n) = bra(1,1:n)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_id_ket(m, n, ket, kp)  ! Returns the tensor product of the mxm identity with a ket
implicit none
integer :: m  ! Dimension of the identity (given as input)
complex(8) :: id(1:m,1:m)  ! Identify
integer :: n  ! Dimension of the vector
complex(8) :: ket(1:n), keta(1,1:n)  ! Vector and its array version
complex(8) :: kp(1:m*n,1:m)  ! The Kronecker product
integer :: j  ! Auxiliary variable for counters

forall(j=1:n) keta(j,1) = ket(j)
kp = 0.d0 ;   forall(j=1:m) kp(((j-1)*n+1):j*n , j) = keta(1:n,1)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_c(M1, nr1, nc1, M2, nr2, nc2, M1_kp_M2)  ! Returns the tensor product of two general complex matrices
implicit none
integer :: nr1, nc1, nr2, nc2  ! Number of rows and columns of the two matrices
complex(8) :: M1(1:nr1,1:nc1), M2(1:nr2,1:nc2)  ! Matrices to take the tensor product of
complex(8) :: M1_kp_M2(1:nr1*nr2,1:nc1*nc2)  ! Matrix containing the tensor product of M1 and M2
integer :: i, j  ! Auxiliary variables for counters

M1_kp_M2 = (0.d0,0.d0) ;   forall ( i = 1:nr1 , j = 1:nc1 ) M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine gram_schmidt(m, n, A, B)  ! Returns, in the coloumns of B, an orthonormal basis obtained from linearly independent vectors 
! given as imput in the columns of A
! Ref: Golub, G. H., and Van Loan, C. F. (2013). Matrix Computations (4th ed.). Baltimore: The Johns Hopkins University Press.
implicit none
integer :: m, n  ! Dimensions of the matrices A and B
complex(8) :: A(1:m,1:n)  ! Matrix with the linearly independent vectors
complex(8) :: B(1:m,1:n)  ! Matrix with the orthonormal basis
real(8) :: norm  ! Vector norm function
complex(8) :: inner_prod  ! Function for the inner product of two complex vectors
integer :: j, k  ! Auxiliary variables for counters

B(:,1) = A(:,1)/norm(m, A(:,1))
do j = 2, n
  B(:,j) = A(:,j)
  do k = 1, j-1
    B(:,j) = B(:,j) - inner_prod(m, B(:,k), A(:,j))*B(:,k)
  enddo
  B(:,j) = B(:,j)/norm(m, B(:,j))
enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine gram_schmidt_modified(m, n, A, B)  ! Returns, in the coloumns of B, an orthonormal basis obtained from linearly 
! independent vectors given as imput in the columns of A
! Ref: Golub, G. H., and Van Loan, C. F. (2013). Matrix Computations (4th ed.). Baltimore: The Johns Hopkins University Press.
implicit none
integer :: m, n  ! Dimensions of the matrices A and B
complex(8) :: A(1:m,1:n)  ! Matrix with the linearly independent vectors
complex(8) :: B(1:m,1:n)  ! Matrix with the orthonormal basis
real(8) :: norm  ! Vector norm function
complex(8) :: inner_prod  ! Function for the inner product of two complex vectors
integer :: j, k  ! Auxiliary variables for counters

do j = 1, n ;   B(:,j) = A(:,j)/norm(m, A(:,j))
  if ( j < n ) then ;   do k = j+1, n ;   A(:,k) = A(:,k) - inner_prod(m, B(:,j), A(:,k))*B(:,j) ;   enddo ;   endif
enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine adjoint(m, n, A, Ad)  ! Returns the adjoint (conjugate transpose) of a m by n complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n)  ! Matrix we want to compute the adjoint of
complex(8) :: Ad(1:n,1:m)  ! Adjoint of the matrix A
integer :: j, k  ! Auxiliary variables for counters

forall( j = 1:m, k = 1:n ) Ad(k,j) = conjg(A(j,k))

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine outer_product(d, psi, phi, op)  ! Returns the outer product of two vectors
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! The vectors
complex(8) :: op(1:d,1:d)  ! The outer product matrix
integer :: j, k  ! Auxiliary variables for counters

forall (j=1:d, k=1:d) op(j,k) = psi(j)*conjg(phi(k))

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine matmul_AAd(m, n, A, AAd)  ! Returns the product of a matrix A by its adjoint
implicit none
integer :: m, n  ! Dimensions of the matrix
complex(8) :: A(1:m,1:n)  ! The Matrix
complex(8) :: AAd(1:m,1:m)  ! Product of A by its adjoint
integer :: j, k, l  ! Auxiliary variables for counters

AAd = 0.d0
do j = 1, m ;   do k = 1, m
  do l = 1, n ;   AAd(j,k) = AAd(j,k) + A(j,l)*conjg(A(k,l)) ;   enddo
enddo ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine mm(m, n, o, A, B, AB) ! Returns the product of the matrices A and B. But only multiplies elements greater than 1.d-15,
! to avoid error propragation.
implicit none
integer :: m, n, o ! Dimensions of the matrices
complex(8) :: A(1:m,1:n), B(1:n,1:o), AB(1:m,1:o)  ! The Matrices
integer :: j, k, l  ! Auxiliary variables for counters

AB = (0.d0,0.d0)
do j = 1, m ;   do k = 1, o
  do l = 1, n ;   if ( (abs(A(j,l)) > 1.d-15) .and. (abs(B(l,k)) > 1.d-15) )  AB(j,k) = AB(j,k) + A(j,l)*B(l,k) ;   enddo
enddo ;   enddo

end
!###################################################################################################################################
!                                       Trace, partial trace, and partial transpose
!###################################################################################################################################
real(8) function trace_he(d, A)  ! Returns the trace of a Hermitian matrix A
implicit none
integer :: d  ! Dimension of the matrix
complex(8) :: A(1:d,1:d)  ! Matrix whose trace is computed
integer :: j  ! Auxiliary variable for counters

trace_he = 0.d0 ;   do j = 1, d ;   trace_he = trace_he + dble(A(j,j)) ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_b(da, db, rho, rho_a)  ! Returns the reduced state of sub-system A
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_a(1:da,1:da)  ! Bipartite and reduced density matrices, respectively
integer :: j, k, l  ! Auxiliary variables for counters

rho_a = (0.d0,0.d0)  
do j = 1, da ;   do k = 1, da ;   do l = 1, db ;   rho_a(j,k) = rho_a(j,k) + rho((j-1)*db+l,(k-1)*db+l) ;   enddo; enddo; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_a(da, db, rho, rho_b)  ! Returns the reduced state of sub-system B
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_b(1:db,1:db)  ! Bipartite and reduced density matrices, respectively
integer :: j, k, l  ! Auxiliary variable for counters

rho_b = (0.d0,0.d0)
do j = 1, db ;   do k = 1, db ;   do l = 1, da ;   rho_b(j,k) = rho_b(j,k) + rho((j-1)*da+l,(k-1)*da+l) ;   enddo; enddo; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_a(da, db, rho, rho_ta)  ! Returns its partial transpose with relation to system A
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_ta(1:da*db,1:da*db)  ! Bipartite original and transposed states
integer :: j, k  ! Auxiliary variable for counters

do j = 1, db ;   do k = 1, db  ! Each pair (j,k) corresponds to a block sub-matrix
  rho_ta(((j-1)*da+1):j*da , ((k-1)*da+1):k*da) = transpose( rho(((j-1)*da+1):j*da , ((k-1)*da+1):k*da) )
enddo ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_b(da, db, rho, rho_tb)  ! Returns the partial transpose with relation to system B
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_tb(1:da*db,1:da*db)  ! Bipartite original and transposed states
integer :: j, k  ! Auxiliary variable for counters

do j = 1, da ;   do k = 1, da  ! Each pair (j,k) corresponds to a block sub-matrix
  rho_tb(((j-1)*db+1):j*db , ((k-1)*db+1):k*db) = transpose( rho(((j-1)*db+1):j*db , ((k-1)*db+1):k*db) )
enddo ;   enddo

end
!###################################################################################################################################
!                                                     Scalar functions
!###################################################################################################################################
real(8) function stokes_parameter(rho, ord_pm, n) !  Returns the Stokes' parameter Tr(rho (sigma_j1 otimes cdots sigma_jn),
! given a n-qubit density matrix and a sequence of n Pauli matrices
implicit none
integer :: n  ! No. of qubits
integer :: ord_pm(n)  ! Vector with the indexes for the order of the Pauli matrices
complex(8) :: rho(1:2**n,1:2**n)  ! Density matrix
complex(8), allocatable :: kp_pauli_mat(:,:) ! Matrix for the Kronecker product of nqb Pauli matrices 
integer :: j, k  ! Auxiliary variables for counters

allocate ( kp_pauli_mat(1:2**n,1:2**n) )

call kron_prod_pauli_mat(ord_pm, n, kp_pauli_mat)
stokes_parameter = 0.d0
do j = 1, 2**n ;   do k = 1, 2**n ;   stokes_parameter = stokes_parameter + dble(kp_pauli_mat(j,k)*rho(k,j)) ;   end do ;   end do
 
deallocate ( kp_pauli_mat )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine stokes_parameters_2qb(rho, ma, mb, corr)  ! Computes the Stokes parameters for a two-qubit density matrix
implicit none
complex(8) :: rho(1:4,1:4)  ! Density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
integer :: ord_pm(1:2)  ! For the order of the Pauli matrix, to be sent to the subroutine which computes the Stokes parameters
real(8) :: stokes_parameter  ! For the Stokes parameter function
integer :: j, k  ! Auxiliary variables for counters

ord_pm(2) = 0 ;   do j = 1, 3;   ord_pm(1) = j ;   ma(j) = stokes_parameter(rho, ord_pm, 2) ;   enddo
ord_pm(1) = 0 ;   do j = 1, 3;   ord_pm(2) = j ;   mb(j) = stokes_parameter(rho, ord_pm, 2) ;   enddo
do j = 1, 3 ;   do k = 1, 3 ;   ord_pm(1) = j ;   ord_pm(2) = k ;   corr(j,k) = stokes_parameter(rho, ord_pm, 2) ;   enddo ; enddo
     
end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function log2(x)  ! Returns the base two log of x
implicit none
real(8):: x

log2 = dlog(x)/dlog(2.d0)

end
!###################################################################################################################################
!                                                        Entropies
!###################################################################################################################################
real(8) function purity(d, rho)  ! Returns the purity of a density matrix: P = Tr(rho*rho)
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! Density matrix
real(8) :: trace_he  ! Trace function

purity = trace_he(d, matmul(rho,rho))

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function neumann(d, rho) ! Returns the von Neumann entropy of a density matrix
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! Density matrix
complex(8) :: A(1:d,1:d)  ! Auxiliary matrices for LAPACK calling
real(8) :: W(1:d)  ! Vector of eigenvalues coming from LAPACK eigensolver
real(8) :: shannon  ! Variables for the shannon and von Neumann entropies

A = rho ;   call lapack_zheevd('N', d, A, W) ;   neumann = shannon(d, W)

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function shannon(d, pv)  ! Returns the Shannon entropy of a probability vector
implicit none
integer :: d  ! Dimension of the probability vector
real(8) :: pv(1:d)  ! probability vector
real(8) :: log2  ! Base two log function
integer :: j  ! Auxiliary variable for counters

shannon = 0.d0
do j = 1, d ;   if( (pv(j) >= 1.d-16) .and. (pv(j) <= (1.d0-1.d-16)) ) shannon = shannon - pv(j)*log2(pv(j)) ;   end do

end
!###################################################################################################################################
!                                         Inner products, Norms, fidelities, and Distances
!###################################################################################################################################
complex(8) function inner_prod(d, psi, phi)  ! Returns the Euclidean inner product of two complex vectors psi and phi
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! Complex vectors whose inner product is to be computed
integer :: j  ! Auxiliary variable

inner_prod = (0.d0,0.d0) ;   do j = 1, d ;   inner_prod = inner_prod + conjg(psi(j))*phi(j) ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function inner_prod_r(d, psi, phi)  ! Returns the Euclidean inner product of two real vectors psi and phi
implicit none
integer :: d  ! Dimension of the vectors
real(8) :: psi(1:d), phi(1:d)  ! Real vectors whose inner product is to be computed
integer :: j  ! Auxiliary variable

inner_prod_r = 0.d0 ;   do j = 1, d ;   inner_prod_r = inner_prod_r + psi(j)*phi(j) ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm(d, psi)  ! Returns the Euclidean norm of the complex vector psi
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: psi(1:d)  ! Complex vector whose norm is to be computed
complex(8) :: inner_prod  ! Inner product function to be used

norm = sqrt(dble(inner_prod(d, psi, psi)))

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_r(d, psi)  ! Returns the Euclidean norm of the real vector psi
implicit none
integer :: d  ! Dimension of the vector
real(8) :: psi(1:d)  ! Real vector whose norm is to be computed
real(8) :: inner_prod_r  ! Inner product function to be used

norm_r = sqrt( inner_prod_r(d, psi, psi) )

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function fidelity_pp(d, psi, phi)  ! Returns the fidelity between two pure states
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! Complex vectors whose fidelity is to be computed
complex(8) :: inner_prod, ip  ! Inner product function

ip = inner_prod(d, psi, phi) ;   fidelity_pp = (dreal(ip))**2.d0 + (dimag(ip))**2.d0

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function fidelity_pm(d, psi, rho)  ! Returns the fidelity between a pure and a mixed state
implicit none
integer :: d  ! Dimension of the vector and matrix
complex(8) :: psi(1:d)  ! Pure state vector
complex(8) :: rho(1:d,1:d)  ! Density matrix
complex(8) :: fid  ! Auxiliary variable
integer :: j, k  ! Auxiliary variables for counters

fid = 0.d0
do j = 1, d ;   do k = 1, d ;   fid = fid + conjg(psi(j))*rho(j,k)*psi(k) ;   enddo ;   enddo
fidelity_pm = dble(fid)

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function fidelity_mm(d, rho, zeta)  ! Returns the fidelity between two mixed states
implicit none
integer :: d  ! Dimension of the density matrices
complex(8) :: rho(1:d,1:d), zeta(1:d,1:d)  ! Density matrices
real(8) :: r(1:d), z(1:d)  ! Vectors for the eigenvalues of rho and zeta, respectively. 
complex(8) :: A(1:d,1:d)  ! Auxiliary matrix. We use r also for the eigenvalues if A=sqrt(rho)*zeta*sqrt(rho)
integer :: j, k, l  ! Auxiliary variables for counters
complex(8) :: op(1:d,1:d)  ! For the outer product
complex(8) ::  inner_prod  ! For the inner product function

call lapack_zheevd('V', d, rho, r) ;   call lapack_zheevd('V', d, zeta, z)

A = 0.d0
do j = 1, d ;   do k = 1, d ;   do l = 1, d
  call outer_product(d, rho(:,j), rho(:,l), op)
  A = A + sqrt(r(j)*r(l))*z(k)*inner_prod(d, rho(:,j), zeta(:,k))*inner_prod(d, zeta(:,k), rho(:,l))*op
enddo ;   enddo ;   enddo

call lapack_zheevd('N', d, A, r) ;   fidelity_mm = 0.d0 ;   do j = 1, d ;   fidelity_mm = fidelity_mm + sqrt(r(j)) ;   enddo
fidelity_mm = (fidelity_mm)**2.d0

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_hs(m, n, A)  ! Returns the Hilbert-Schmidt norm (or 2-norm) of a complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n)  ! Matrix of which want to compute the 2-norm
complex(8), allocatable :: Ad(:,:)  ! Adjoint of A
real(8) :: trace_he  ! Function for the trace of an Hermitian matrix

allocate( Ad(1:n,1:m) ) ;   call adjoint(m, n, A, Ad) ;   norm_hs = dsqrt( trace_he(m, matmul(A,Ad)) ) ;   deallocate( Ad )

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_tr(d, A) ! Returns the trace norm (or 1-norm) of an hermitian matrix A
! ||A||_1 = \sum_j |a_j|, where a_j are the eigenvalues of A 
implicit none
integer :: d  ! Dimension of the matrix A
complex(8) :: A(1:d,1:d)  ! Matrix of which want to compute the 1-norm
real(8) :: egv(1:d)  ! Eigenvalues of A
integer :: j  ! Auxiliary variable for counters

call lapack_zheevd('N', d, A, egv) ;   norm_tr = 0.d0 ;   do j = 1, d ;   norm_tr = norm_tr + abs(egv(j)) ; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_tr_ge(d, A) ! Returns the trace norm (or 1-norm) of a general normal matrix A
! ||A||_1 = \sum_j |a_j|, where a_j are the eigenvalues of A 
implicit none
integer :: d  ! Dimension of the matrix A
complex(8) :: A(1:d,1:d)  ! Matrix of which want to compute the 1-norm
complex(8) :: egv(1:d)  ! Eigenvalues of A
integer :: j  ! Auxiliary variable for counters

call lapack_zgeev('N', d, A, egv) ;   norm_tr_ge = 0.d0 ;   do j = 1, d ;   norm_tr_ge = norm_tr_ge + abs(egv(j)) ; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_l1(m, n, A)  ! Returns the l1-norm of a complex mxn complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n) ! The matrix
integer :: j, k

norm_l1 = 0.d0 ;   do j=1,m ;   do k=1,n ;   if ( j /= k ) norm_l1 = norm_l1 + abs(A(j,k)) ;   enddo ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function d_hs(d, rho, zeta)  ! Returns the Hilbert-Schmidt (2-norm) distance between two density matrices
implicit none
integer :: d  ! Dimension of the density matrices
complex(8) :: rho(1:d,1:d), zeta(1:d,1:d)  ! Density matrices
real(8) :: norm_hs  ! Function for the Hilbert-Schmidt norm

d_hs = norm_hs(d, d, rho - zeta)

end
!###################################################################################################################################
!                                                          LAPACK callers
!###################################################################################################################################
subroutine lapack_zgeev(JOBVR, N, A, Wc)  ! Calls LAPACK's eigensolver for general complex matrices
! ZGEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will 
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray 
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
!character(1) :: JOBVL  !  JOBVL is CHARACTER*1; = 'N': left eigenvectors of A are not computed; = 'V': left eigenvectors of are computed.
character(1) :: JOBVR  !  JOBVR is CHARACTER*1; = 'N': right eigenvectors of A are not computed;  = 'V': right eigenvectors of A are computed.
character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA = N  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(1:N,1:N)  ! A is COMPLEX*16 array, dimension (LDA,N); On entry, the N-by-N matrix A. On exit, A has been overwritten.
complex(8) :: Wc(1:N)  ! Wc is COMPLEX*16 array, dimension (N). Wc contains the computed eigenvalues.
!integer :: LDVL = N  !LDVL is INTEGER. The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
complex(8) :: VL(1:N,1:N)  ! VL is COMPLEX*16 array, dimension (LDVL,N); If JOBVL = 'V', the left eigenvectors u(j) are stored one
                        ! after another in the columns of VL, in the same order as their eigenvalues.
                        ! If JOBVL = 'N', VL is not referenced. u(j) = VL(:,j), the j-th column of VL.
!integer :: LDVR = N  ! LDVR is INTEGER. The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
complex(8) :: VR(1:N,1:N)  ! VR is COMPLEX*16 array, dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one
                                     ! after another in the columns of VR, in the same order their eigenvalues.
                                     ! If JOBVR = 'N', VR is not referenced. v(j) = VR(:,j), the j-th column of VR.
                                     
!integer :: LWORK = 2*N  ! LWORK is INTEGER; The dimension of the array WORK.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.
                  ! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                  ! this value as the first entry of the WORK array, and no error related to LWORK is issued by XERBLA.                                     
complex(8) :: WORK(1:2*N)  ! WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
real(8) :: RWORK(1:2*N)  ! RWORK is DOUBLE PRECISION array, dimension (2*N)                    
integer :: INFO   ! INFO is INTEGER
                  ! = 0:  successful exit
                  ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                  ! > 0:  if INFO = i, the QR algorithm failed to compute all the
                  ! eigenvalues, and no eigenvectors have been computed; elements and i+1:N of W contain eigenvalues which have converged.
                
call zgeev ('N', JOBVR, N, A, N, Wc, VL, N, VR, N, WORK, 2*N, RWORK, INFO)
!call zgeev (JOBVL, JOBVR, N, A, LDA, Wc, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zheevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for Hermitian complex matrices
! ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will 
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray 
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(N,N)  ! A is COMPLEX array, dimension (LDA, N). On entry, the Hermitian matrix A.  If UPLO = 'U', the
                      ! leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A.  
                      ! If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower triangular part of the 
                      ! matrix A. On exit, if JOBZ = 'V', then if INFO = 0, A contains the orthonormal eigenvectors of the 
                      ! matrix A. If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') or the upper triangle 
                      ! (if UPLO='U') of A, including the diagonal, is destroyed.
real(8) :: W(N)  ! W is REAL array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
          !The length of the array WORK.
          !If N <= 1,                LWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
          !If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.

          !If LWORK = -1, then a workspace query is assumed; the routine
          !only calculates the optimal sizes of the WORK, RWORK and
          !IWORK arrays, returns these values as the first entries of
          !the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
complex(8), allocatable :: WORK(:)  !WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LRWORK   ! LRWORK is INTEGER
          !The dimension of the array RWORK.
          !If N <= 1,                LRWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
          !If JOBZ  = 'V' and N > 1, LRWORK must be at least
           !              1 + 5*N + 2*N**2.

          !If LRWORK = -1, then a workspace query is assumed; the
          !routine only calculates the optimal sizes of the WORK, RWORK
          !and IWORK arrays, returns these values as the first entries
          !of the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
real(8), allocatable :: RWORK(:)    ! RWORK is DOUBLE PRECISION array,
                                        ! dimension (LRWORK)
         ! On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
integer :: LIWORK   ! LIWORK is INTEGER
          !The dimension of the array IWORK.
          !If N <= 1,                LIWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
          !If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.

          !If LIWORK = -1, then a workspace query is assumed; the
          !routine only calculates the optimal sizes of the WORK, RWORK
          !and IWORK arrays, returns these values as the first entries
          !of the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
          !On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
                  ! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements 
                  ! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the 
                  ! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns 
                  ! INFO/(N+1) through mod(INFO,N+1).
                  
     if (JOBZ == 'N') then ;   LWORK = N + 1 ;   LRWORK = N ;   LIWORK = 1
else if (JOBZ == 'V') then ;   LWORK = 2*N + N**2 ;   LRWORK = 1 + 5*N + 2*N**2 ;   LIWORK = 3 + 5*N
     endif
allocate( WORK(1:LWORK), RWORK(1:LRWORK), IWORK(1:LIWORK) )

call zheevd(JOBZ,'U',N,A,N,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
!call zheevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
  
deallocate( WORK, RWORK, IWORK )
  
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_dsyevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for symmetric real matrices
! DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
! real symmetric matrix A. If eigenvectors are desired, it uses a
! divide and conquer algorithm.!

! The divide and conquer algorithm makes very mild assumptions about
! floating point arithmetic. It will work on machines with a guard
! digit in add/subtract, or on those binary machines without guard
! digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
! Cray-2. It could conceivably fail on hexadecimal or decimal machines
! without guard digits, but we know of none.

! Because of large use of BLAS of level 3, DSYEVD needs N**2 more
! workspace than DSYEVX.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
real(8) :: A(1:N,1:N)  ! A is DOUBLE PRECISION array, dimension (LDA, N)
          !On entry, the symmetric matrix A.  If UPLO = 'U', the
          !leading N-by-N upper triangular part of A contains the
          !upper triangular part of the matrix A.  If UPLO = 'L',
          !the leading N-by-N lower triangular part of A contains
          !the lower triangular part of the matrix A.
          !On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          !orthonormal eigenvectors of the matrix A.
          !If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
          !or the upper triangle (if UPLO='U') of A, including the
          !diagonal, is destroyed.
real(8) :: W(1:N)  ! W is DOUBLE PRECISION array, dimension (N)
          !If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
          !The dimension of the array WORK.
          !If N <= 1,               LWORK must be at least 1.
          !If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
          !If JOBZ = 'V' and N > 1, LWORK must be at least
          !                                      1 + 6*N + 2*N**2.

          !If LWORK = -1, then a workspace query is assumed; the routine
          !only calculates the optimal sizes of the WORK and IWORK
          !arrays, returns these values as the first entries of the WORK
          !and IWORK arrays, and no error message related to LWORK or
          !LIWORK is issued by XERBLA.
real(8), allocatable :: WORK(:)  !WORK is DOUBLE PRECISION array,
           !                              dimension (LWORK)
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LIWORK   ! LIWORK is INTEGER
          !The dimension of the array IWORK.
          !If N <= 1,                LIWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
          !If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.!

         ! If LIWORK = -1, then a workspace query is assumed; the
         ! routine only calculates the optimal sizes of the WORK and
         ! IWORK arrays, returns these values as the first entries of
         ! the WORK and IWORK arrays, and no error message related to
         ! LWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
          !On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
                  ! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements 
                  ! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the 
                  ! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns 
                  ! INFO/(N+1) through mod(INFO,N+1).
                  
     if (JOBZ == 'N') then ;   LWORK = 2*N+1 ;   LIWORK = 1 
else if (JOBZ == 'V') then ;   LWORK =  1 + 6*N + 2*N**2 ;   LIWORK = 3 + 5*N
     endif
allocate( WORK(1:LWORK), IWORK(1:LIWORK) )

call dsyevd(JOBZ,'U',N,A,N,W,WORK,LWORK,IWORK,LIWORK,INFO)
!call dsyevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO)

deallocate(WORK, IWORK)
    
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zgeqrfp(N, Z, Q, R)  ! Calls LAPACK QR factorization 
! For a complex square matrix A, it returns the orthogonalized Q and upper triangular R matrices
implicit none
!ZGEQRFP computes a QR factorization of a complex M-by-N matrix A:
! A = Q * R. The diagonal entries of R are real and nonnegative.

!integer :: M !is INTEGER
          !The number of rows of the matrix A.  M >= 0.
integer :: N !is INTEGER
          !The number of columns of the matrix A.  N >= 0.
!integer :: LDA !is INTEGER
!          !The leading dimension of the array A.  LDA >= max(1,M).
complex(8) :: A(1:N,1:N) !is COMPLEX*16 array, dimension (LDA,N)
          !On entry, the M-by-N matrix A.
          !On exit, the elements on and above the diagonal of the array
          !contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          !upper triangular if m >= n); the elements below the diagonal,
          !with the array TAU, represent the orthogonal matrix Q as a
          !product of min(m,n) elementary reflectors (see Further
          !Details).
complex(8) :: TAU(1:N) !is COMPLEX*16 array, dimension (min(M,N))
          !The scalar factors of the elementary reflectors
          !Further Details
          !The matrix Q is represented as a product of elementary reflectors
          !  Q = H(1) H(2) . . . H(k), where k = min(m,n).
          !Each H(i) has the form
          !  H(i) = I - tau * v * v**H 
          !where tau is a real scalar, and v is a real vector with
          !v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
          !and tau in TAU(i).
!integer :: LWORK !is INTEGER
          !The dimension of the array WORK. The dimension can be divided into three parts.
          ! 1) The part for the triangular factor T. If the very last T is not bigger 
          !    than any of the rest, then this part is NB x ceiling(K/NB), otherwise, 
          !   NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T              
          !2) The part for the very last T when T is bigger than any of the rest T. 
          !   The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB,
          !   where K = min(M,N), NX is calculated by
           !        NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
          ! 3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB)
          ! So LWORK = part1 + part2 + part3
          ! If LWORK = -1, then a workspace query is assumed; the routine
          ! only calculates the optimal size of the WORK array, returns
          ! this value as the first entry of the WORK array, and no error
          ! message related to LWORK is issued by XERBLA.
complex(8) :: WORK(1:3*N*N) !is COMPLEX*16 array, dimension (MAX(1,LWORK))
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: INFO !is INTEGER
          != 0:  successful exit
          !< 0:  if INFO = -i, the i-th argument had an illegal value
integer :: j, k  ! Auxiliary variables
complex(8) :: identity(1:N,1:N), proj(1:N,1:N), v(1:N)
complex(8) :: Z(1:N,1:N), Q(1:N,1:N), R(1:N,1:N)  ! Input and output matrices

A = Z
call zgeqrfp (N, N, A, N, TAU, WORK, 3*N*N, INFO)
!call zgeqrfp (M, N, A, LDA, TAU, WORK, LWORK, INFO)

R = 0.d0 ;   forall (j = 1:N, k = 1:N , k >= j ) R(j,k) = A(j,k)

call identity_c(N, identity) ;   Q = identity
do j = 1, N
  if (j > 1) v(1:j-1) = 0.d0 ;   v(j) = 1.d0 ;   v(j+1:N) = A(j+1:N,j) ;   call projector(v, N, proj)
  Q = matmul(Q,(identity-TAU(j)*proj)) 
enddo

end
!###################################################################################################################################

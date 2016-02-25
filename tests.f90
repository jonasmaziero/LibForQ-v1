!###################################################################################################################################
!                                               Tests - LAPACK eigensolvers
!###################################################################################################################################
subroutine test_lapack()  ! Some simple tests on LAPACK eigensolvers
!use qsort_mod  ! Used to sort the eigenvalues obtained with zgeev
implicit none
integer, parameter :: nqb = 8  ! No. of qubits
integer :: ord_pm(1:nqb)  ! Vector for the order of the Pauli matrices
integer, parameter :: d = 2**nqb  ! Dimension of the matrices and vectors
complex(8) :: kp(1:d,1:d), A(1:d,1:d)  ! Matrix for the Kronecker product and the one to be diagonalized
real(8) :: W(1:d)  ! Vector of real eigenvalues
complex(8) :: Wc(1:d)  ! Vector of complex eigenvalues
real(8) :: rn(1:nqb)  ! For the random numbers use to choose ord_pm
integer :: j  ! Auxiliary variable for counters
real(8) :: t1, t2  ! Times

call cpu_time(t1)

write(*,*) "## Performing some tests for the LAPACK eigensolvers"

call rng(nqb, rn)  ! Gets the random numbers
do j = 1, nqb  ! Choose the indexes for the order of the Pauli matrices
     if ( (rn(j) >= 0.d0) .and. (rn(j) < 0.25d0) ) then; ord_pm(j) = 0
else if ( (rn(j) >= 0.25d0) .and. (rn(j) < 0.5d0) ) then; ord_pm(j) = 1
else if ( (rn(j) >= 0.5d0) .and. (rn(j) < 0.75d0) ) then; ord_pm(j) = 2
else if ( (rn(j) >= 0.75d0) .and. (rn(j) <= 1.d0) ) then; ord_pm(j) = 3
     endif
enddo
write(*,*) "ord_pm",ord_pm

call kron_prod_pauli_mat(ord_pm, nqb, kp)  ! Gets the tensor product of the choosed Pauli matrices
A = kp ;   call lapack_zheevd('N', d, A, W)  ! Calls LAPACK zheevd, which computes the eigenvalues of Hermitian matrices
write(*,*) "zheevd", sum(W)  ! The sum should be zero

A = kp ;   call lapack_zgeev('N', d, A, Wc)  ! Calls LAPACK zgeev, which computes the eigenvalues of general matrices
W = real(Wc) ;   call qsort(W,d)  ! We define A again because it could be change by zheevd
write(*,*) "zgeev", sum(W), ". The sum should be zero."

call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!###################################################################################################################################
!                                                       Tests - RNG
!###################################################################################################################################
subroutine rng_tests()  ! Some simple tests for the random number generators
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
use meths
implicit none
integer, parameter :: ni = 100  ! No. of intervals we divide the domain of the random numbers
integer, parameter :: ns = 10**6  ! No. of samples for the averages and relative frequencies (must be even for gaussian rn)
real(8) :: rn(1:ns), rnu(1:ns), rnga(1:ns), rne(1:ns)  ! Vectors of random numbers
integer :: ct(1:ni) = 0, ctu(1:ni) = 0, ctg(1:ni) = 0, cte(1:ni) = 0  ! Counter for the probability distribution estimation
real(8), parameter :: xmin = -5.d0, xmax = 5.d0  ! Domain for the random numbers
real(8), parameter :: delta = (xmax - xmin)/dble(ni)  ! Step in the above domain
real(8) :: x1, x2, x3  ! For the 3D scatter plot
integer :: j, l  ! Auxiliary variables for counters
real(8) :: t1, t2  ! Times

call cpu_time(t1)

write(*,*) "## Performing some tests for the random number generators"

opt_rng = "mt" ;   call rng_init()  ! Sets the method to be used in the generators ans Initializes it


write(*,*) "# Computing the probability distributions"
open(unit = 11, file = 'rn_pd.dat', status = 'unknown')
! Gnuplot commands to see the results
! plot 'rn_pd.dat' u 1:2 w p t 'u in [0,1]', '' u 1:3 w p t 'u in [a,b]', '' u 1:4 w p t 'gaussian', '' u 1:5 w p t 'exponential'
call rng(ns,rn) ;   call rng_unif(ns, xmin, xmax, rnu) ;    call rng_gauss(ns, rnga) ;   call rng_exp(ns, rne)
do j = 1, ns ;   do l = 1, ni-1
  if ( (rn(j) >= (xmin + dble(l)*delta)) .and. (rn(j) < (xmin + dble(l+1)*delta)) ) ct(l) = ct(l) + 1
  if ( (rnu(j) >= (xmin + dble(l)*delta)) .and. (rnu(j) < (xmin + dble(l+1)*delta)) ) ctu(l) = ctu(l) + 1
  if ( (rnga(j) >= (xmin + dble(l)*delta)) .and. (rnga(j) < (xmin + dble(l+1)*delta)) ) ctg(l) = ctg(l) + 1
  if ( (rne(j) >= (xmin + dble(l)*delta)) .and. (rne(j) < (xmin + dble(l+1)*delta)) ) cte(l) = cte(l) + 1
enddo ;   enddo
do l = 1, ni-1
  write(11,*) (xmin +(dble(l)+0.5d0)*delta),dble(ct(l))/dble(ns),dble(ctu(l))/dble(ns),dble(ctg(l))/dble(ns),dble(cte(l))/dble(ns)
enddo
close(11)


write(*,*) "# Generating data for the 2d and 3d scatter plots"
open(unit = 12, file = 'rn_2d.dat', status = 'unknown')
open(unit = 13, file = 'rn_3d.dat', status = 'unknown')
! Gnuplot commands to see the results
! plot 'rn_2d.dat' 
! splot 'rn_3d.dat'
call rng(10**4, rn) ;   do j = 1, 10**4, 2 ;   write(12,*) rn(j), rn(j+1) ;   enddo  ! For the 2D scatter plot
close(12)
do j = 1, 10**4, 3  ! For the 3D scatter plot
  x1 = rn(j)/dsqrt(rn(j)**2.d0 + rn(j+1)**2.d0 + rn(j+2)**2.d0)
  x2 = rn(j+1)/dsqrt(rn(j)**2.d0 + rn(j+1)**2.d0 + rn(j+2)**2.d0)
  x3 = rn(j+2)/dsqrt(rn(j)**2.d0 + rn(j+1)**2.d0 + rn(j+2)**2.d0)
  write(13,*) x1, x2, x3
enddo
close(13)


write(*,*) "# Computing the average, moments, and correlations"
call rng_corr() 

call cpu_time(t2) ;   write(*,*) 'Time:',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!--------------------
subroutine rng_corr()  
! Computes some averages, moments, and correlations and write them on a file
implicit none
integer :: ns ! No. of samples
integer :: j, l  ! Aux variables for counters
real(8) :: average  ! Average of the random numbers
real(8) :: corr12, corr13, corr14, corr15  ! Correlations between subsequent random numbers
real(8) :: mu1, mu2, mu3, mu4, mu5  ! Moments of the distribution of random numbers
real(8), allocatable :: rn(:)  ! Vectors of random numbers
open(unit = 14, file = 'rn_corr.dat', status = 'unknown')
! Gnuplot commands to see the results
! plot 'rn_corr_gnu.dat' u 1:2 w lp t 'avg'
! plot 'rn_corr.dat' u 1:3 w lp t 'mu1', '' u 1:4 w lp t 'mu2', '' u 1:5 w lp t 'mu3', '' u 1:6 w lp t 'mu4', '' u 1:7 w lp t 'mu5'
! plot 'rn_corr.dat' u 1:8 w lp t 'c12', '' u 1:9 w lp t 'c13', '' u 1:10 w lp t 'c14', '' u 1:11 w lp t 'c15'

do l = 3, 25
  ns = 2**l
  allocate( rn(ns) )
  average = 0.d0 ;   corr12 = 0.d0 ;   corr13 = 0.d0 ;   corr14 = 0.d0 ;   corr15 = 0.d0
  mu1 = 0.d0 ;   mu2 = 0.d0 ;   mu3 = 0.d0 ;   mu4 = 0.d0 ;   mu5 = 0.d0
  call rng(ns, rn)
  do j = 1, ns 
    average = average + rn(j) ! average & correlation functions
    mu1 = mu1 + rn(j)**1.d0; mu2 = mu2 + rn(j)**2.d0; mu3 = mu3 + rn(j)**3.d0; mu4 = mu4 + rn(j)**4.d0; mu5 = mu5 + rn(j)**5.d0
    if ( j >= 2 ) corr12 = corr12 + rn(j)*rn(j-1); if ( j >= 3 ) corr13 = corr13 + rn(j)*rn(j-2)
    if ( j >= 4 ) corr14 = corr14 + rn(j)*rn(j-3); if ( j >= 5 ) corr15 = corr15 + rn(j)*rn(j-4)
  enddo
  average = average/dble(ns) ! average
  mu1 = dabs( (mu1/ns) - (1.d0/(1.d0+1.d0)) ) ;   mu2 = dabs( (mu2/ns) - (1.d0/(1.d0+2.d0)) )
  mu3 = dabs( (mu3/ns) - (1.d0/(1.d0+3.d0)) ) ;   mu4 = dabs( (mu4/ns) - (1.d0/(1.d0+4.d0)) )
  mu5 = dabs( (mu5/ns) - (1.d0/(1.d0+5.d0)) )
  corr12 = corr12/(dble(ns)-1.d0) - average**2.d0 ;   corr13 = corr13/(dble(ns)-2.d0) - average**2.d0
  corr14 = corr14/(dble(ns)-3.d0) - average**2.d0 ;   corr15 = corr15/(dble(ns)-4.d0) - average**2.d0
  write(14,*) log(dble(ns)), average-0.5, mu1, mu2, mu3, mu4, mu5, corr12, corr13, corr14, corr15
  deallocate( rn )
enddo
close(14)

end
!###################################################################################################################################
!                                                       Tests - RPVG
!###################################################################################################################################
subroutine rpvg_tests() ! Some simple tests for the random probability vector generators
use meths
implicit none
integer, parameter :: d = 3  ! Dimension of the random probability vectors
integer, parameter :: ns = 10**6  ! No. of samples for the averages and probability distributions
integer, parameter :: ni = 100  ! No. of intervals for the domain of the RPV components
real(8) :: delta = 1.d0/dble(ni) ! Step for the domain of the RPV components
real(8) :: rpv(1:d), avg_rpv(1:d)  ! Random probability vector and its average
real(8) :: t1, t2  ! Times
integer :: j, k, l  ! Auxiliary variables for counters
integer :: ct(1:ni,1:d) = 0  ! Counter for the probability densities

call cpu_time(t1)

write(*,*) "## Performing some tests for the random probability vector generators"

! Sets the methods to be used in the generators
opt_rng = "mt" ;   opt_rpvg = "zhsl"


write(*,*) "# Generating a sampling for 2d scatter plots"
open(unit = 11, file = 'rpv_2d_norm_d3.dat', status = 'unknown') ! Gnuplot commands to see the results: plot 'rpv_2d.dat'
do j = 1, 5*10**3 ;   call rpvg(d, rpv) ;  write(11,*) (rpv(k) , k=1,2) ;   enddo
!do j = 1, 5*10**3 ;   call rpvg(d, rpv) ;  write(11,*) ((1.d0 - rpv(k)) , k=1,2) ;   enddo  ! To put the points above the 'diagonal'
close(11)


write(*,*) "# Computing the probability distribution for the components of the RPVs"
open(unit = 12, file = 'rpv_pd_norm_d3.dat', status = 'unknown') ! Gnuplot commands to see the results: plot 'rpv_pd.dat'
avg_rpv = 0.d0
do j = 1, ns
  call rpvg(d, rpv) ;   avg_rpv = avg_rpv + rpv
  do k = 1, d ;   if (rpv(k) == 1.d0) rpv(k) = 1.d0 - 1.d-10
    do l = 1, ni ;   if ( (rpv(k) >= (dble(l)-1.d0)*delta) .and. (rpv(k) < dble(l)*delta)) ct(l,k) = ct(l,k) + 1 ;   enddo
  enddo
enddo
! Writes  on the screen the average of the (up to) first 5 components of the RPV
if ( d <= 5 ) write(*,*) 'avg_rpv = ', avg_rpv/dble(ns)
do l = 1, ni
  write(12,*) (dble(l)-0.5)*delta, dble(ct(l,1))/dble(ns)  ! Writes the probability density for the first component of the RPV on a file
  !write(12,*) (dble(l)-0.5)*delta, (dble(ct(l,k))/dble(ns), k=1,d)  ! Writes the probability density for all components of the RPV on a file
enddo
close(12)

call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!###################################################################################################################################
!                                                       Tests - RUG
!###################################################################################################################################
subroutine rug_tests()  ! Some simple tests for the random unitary generators
! Ref: \.{Z}yczkowski, K., and Ku\'{s}, M. (1994). Random unitary matrices, J. Phys. A: Math. Gen. 27, 4235.
use meths
!use qsort_mod
implicit none
integer, parameter :: d = 20  ! Dimension of the unitary matrices
integer, parameter :: ns = 10**4  ! No. of samples for estimating the probability distributions
integer, parameter :: ni = 100  ! No. of intervals we divide the domain under analysis
complex(8), allocatable :: ru(:,:)  ! Random unitary matrix
complex(8), allocatable :: egv(:)  ! Vector of eigenvalues of ru
real(8) :: delta  ! Step in the domain of the eigenvalues and of its spacings
real(8), allocatable :: egp(:)  ! Phases in the eigevalues of the random unitaries, and its sorted version
real(8), allocatable :: sp(:)  ! Spacings of the eigenvalues
integer :: j, k, l, m, n  ! Auxiliary variables for counters
integer, allocatable :: ct(:), cts(:)  ! Counters for the probability distributions
real(8) :: asp  ! Average eigenphase spacing for the random matrices eigenvalues
real(8), parameter :: pi = 4.d0*atan(1.d0)
complex(8), allocatable :: rud(:,:)  ! For the adjoint of ru
real(8) :: norm_l1, trace_he  ! For the l1-norm and trace of UU^†
real(8) :: ip_jj, ip_jk  ! For the sum of the inner products for j=j or j\=k columns of U
complex(8) :: inner_prod  ! Variable for the inner product fucntion
real(8), allocatable :: pv(:)  ! Probability vector
real(8) :: avg_ent, shannon, exact, log2  ! For entropies
real(8) :: t1, t2  ! For the computation time

call cpu_time(t1)

write(*,*) "## Performing some tests for the random unitary generators"

! Sets the methods to be used in the generators
opt_rng = "mt" ;   opt_rug = "hurwitz"


allocate( ru(1:d,1:d), egv(1:d), egp(1:d), sp(1:d-1), ct(1:ni), cts(1:ni), rud(1:d,1:d) )


write(*,*) "# Testing the orthonormality of the column vectors of U"
call rug(d, ru) ;   ip_jj = 0.d0 ;   ip_jk = 0.d0 
do j = 1, d ;   ip_jj = ip_jj + abs(inner_prod(d, ru(:,j), ru(:,j))) ;   enddo
do j = 1, d ;   do k = 1, d ;   if ( j /= k ) ip_jk = ip_jk + abs(inner_prod(d, ru(:,j), ru(:,k))) ;   enddo ;   enddo
write(*,*) 'd:',d,', sum |<v_j|v_j>|:', ip_jj, ', sum |<v_j|v_k>|:', ip_jk


write(*,*) "# Computing the L1-norm and trace of UU^†"
call rug(d, ru) ;   call adjoint(d, d, ru, rud)
write(*,*) 'l1-norm of UU^†:', norm_l1(d, d, matmul(ru,rud)), ',     d:', d, ',   Tr(UU^†):',trace_he(d, matmul(ru,rud))


write(*,*) "# Generating data to plot the circle of eigenvalues"
open(unit = 11, file = 'ru_egval_circle_hur.dat', status = 'unknown')  ! Gnuplot commands to see the results: plot 'ru_egval_circle.dat'
j = -d+1
do
  j = j + d ;   if ( j > 10**3 ) exit ;   call rug(d, ru) ;   call lapack_zgeev('N', d, ru, egv) 
  do k = 1, d ;   write(11,*) dble(egv(k)), dimag(egv(k)) ;   enddo
enddo
close(11)


write(*,*) "# Computing the probability density for the eigevalues of the unitaries and for their spacings"
open(unit = 12, file = 'ru_egval_dist.dat', status = 'unknown')
! Gnuplot commands to see the results: plot [0:2*pi][0:] 'ru_egval_dist.dat', '' u 1:3
delta = (2.d0*pi)/dble(ni)  ! Step for the eigenphases
 ct = 0 ;   cts = 0  ! Initializes the counters
do j = 1, ns
  call rug(d, ru) ;   call lapack_zgeev('N', d, ru, egv)  ! Gets the random unitary and its eigenvalues
  forall(m=1:d) egp(m) = dble(-(0.d0,1.d0)*log(egv(m))) + pi ;   call qsort(egp, d)  ! Puts the eigenphases in [0,2*pi] and sort them
  forall(m=1:d-1) sp(m) = egp(m+1)-egp(m) ;   asp = sum(sp)/dble(d-1) ;   sp = sp/asp
  do k = 1, d  ! Counter for the different eigenvalues and spacings
  do l = 1, ni  ! Counter for finding where a certain eigenvalues or spacing is in the whole range
    if ( (egp(k) >= ((dble(l)-1.d0)*delta)) .and. (egp(k) < (dble(l)*delta)) ) then  ! Density of the eigenvalues
      ct(l) = ct(l) + 1 
    endif
    if ( (k<d) .and. (sp(k) >= ((dble(l)-1.d0)*delta)) .and. (sp(k) < (dble(l)*delta)) ) then  ! Density of the eigenvalues spacings                                
       cts(l) = cts(l) + 1                                                                                       
    endif                                                                                                                   
enddo ;   enddo ;   enddo   
write(*,*) 'Sum egval counting:',dble(sum(ct))/dble(ns*d), ', Sum egval spacing counting:', dble(sum(cts))/dble(ns*(d-1))  ! Verify normalizations
do l = 1, ni ;   write(12,*) (dble(l)-0.5)*delta, dble(ct(l))/dble(ns*d), dble(cts(l))/dble(ns*(d-1)) ;   enddo  ! Write the densities on a file
close(12)

deallocate( ru, egv, egp, sp, ct, cts, rud )

call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!###################################################################################################################################
!                                                       Tests - RSVG
!###################################################################################################################################
subroutine rsvg_tests() ! Some simple tests for the random state vector generator
use meths
implicit none
integer :: ns = 10**3  ! No. of samples for the averages
integer :: d  ! Dimension of the state vector
real(8) :: avg_ent, neumann, trace_he  ! Average and von Neumann entropies, and the trace function
complex(8), allocatable :: rsv(:), rsv2(:)  ! Random state vectors
complex(8), allocatable :: rho(:,:), rhoa(:,:)  ! Density matrix and reduced density matrix
complex(8), allocatable :: proj(:,:), identity(:,:) ! Projector and identity matrix
real(8) :: norm  ! Norm of a vector
integer :: j, k, l  ! Auxiliary variables for counters
integer :: N  ! For the mixture, for comparison with the identity
real(8) :: t1, t2  ! Times
real(8) :: d_hs, avg_d_hs  ! Hilbert-Schmidt distance and its average
real(8) :: fidelity_pp, avgf  ! Fidelity and its average

call cpu_time(t1)

write(*,*) "## Performing some tests for the random state vector generators"

! Sets the methods to be used in the generators
opt_rng = "mt" ;   opt_rpvg = "zhsl" ;   opt_rug = "gso" ;   opt_rsvg = "std"


write(*,*) "# Computing the NORM of some RSVs"
do j = 1, 5 ;   d = 2**j 
 allocate( rsv(1:d) ) ;   call rsvg(d,rsv) ;  write(*,*) "No. qubits", j, "norm", norm(d, rsv) ;   deallocate(rsv)
enddo

  
write(*,*) "# Computing the average fidelity between random state vectors"
open(unit = 11, file = 'rsv_fid.dat', status = 'unknown')
! Gnuplot commands to see the results: plot 1/x w l t 'exact', 'rsv_fid.dat' w p t 'numerics'
do  d = 2, 50, 2
  allocate( rsv(1:d), rsv2(1:d) ) ;   avgf = 0.d0
  do j = 1, ns ;   call rsvg(d, rsv) ;   call rsvg(d, rsv2) ;   avgf = avgf + fidelity_pp(d, rsv, rsv2) ;   enddo
  write(11,*) d, avgf/dble(ns) ;   write(*,*) 'd:',d, ', avg fidelity:', avgf/dble(ns)
  deallocate( rsv, rsv2)
enddo
close(11)


write(*,*) "# Computing the average entanglement entropy of random state vectors"
open(unit = 12, file = 'rsv_EE.dat', status = 'unknown') 
! Gnuplot commands to see the results: plot 'rsv_entanglement.dat' w lp
do d = 2, 25, 2  ! Changes of dimension of the RSV, which is = d*d
  allocate( rsv(1:d*d), rho(1:d*d,1:d*d), rhoa(1:d,1:d) ) ;   avg_ent = 0.d0
  do j = 1, ns
    call rsvg(d*d, rsv) ;   call projector(rsv, d*d, rho) ;   call partial_trace_b(d, d, rho, rhoa)
    avg_ent = avg_ent + neumann(d, rhoa)
  enddo
  avg_ent = avg_ent/dble(ns) ;   write(12,*) d, avg_ent ;   write(*,*) 'd:',d, ', avg entaglement:', avg_ent/dble(ns)
  deallocate( rsv, rho, rhoa )
enddo
close(12)


write(*,*) "# Computing the distance from a mixture of RSVs to the identity"
open(unit = 13, file = 'rsv_d_id.dat', status = 'unknown')  ! For the distance of a mixture of RSVs to the identity
! Gnuplot commands to see the results: plot 'rsv_d_id.dat' w lp
d = 2**2  ! You can change the dimension to see what changes 
allocate( rsv(1:d), proj(1:d,1:d), rho(1:d,1:d) , identity(1:d,1:d) ) ;   call identity_c(d, identity)
do N = 1, 10**2, 10**1  ! Changes the no. of RSVs used in the mixture
  avg_d_hs = 0.d0
  do k = 1, ns
    rho = 0.d0
    do j = 1, N ;   call rsvg(d,rsv) ;   call projector(rsv, d, proj) ;  rho = rho + proj ;   enddo
    avg_d_hs = avg_d_hs + d_hs(d, rho/real(N), identity/real(d))
  enddo
  write(13,*) N, avg_d_hs/real(ns)
enddo
do N = 10**2, 10**3, 2*10**2  ! Same as the previous do, but for higher values of N
  avg_d_hs = 0.d0
  do k = 1, ns
    rho = 0.d0
    do j = 1, N ;   call rsvg(d,rsv) ;   call projector(rsv, d, proj) ;  rho = rho + proj ;   enddo
    avg_d_hs = avg_d_hs + d_hs(d, rho/real(N), identity/real(d))
  enddo
  write(13,*) N, avg_d_hs/real(ns)
enddo
close(13)

call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!###################################################################################################################################
!                                                       Tests - RDMG
!###################################################################################################################################
subroutine rdmg_tests() ! Some simple tests for the random density matrix generator
use meths
implicit none
integer :: ns = 10**4  ! No. of samples for the averages
integer :: d, da, db  ! For the dimension(s) of the density operator, or state vector
complex(8), allocatable :: rsv(:), rdm(:,:), rdm_pt(:,:)  ! For the random state vector, density matrix, and its partial transpose
real(8), allocatable :: egv(:)  ! For the eigenvalues of rdm_pt
integer :: ct  ! For the density of positive partial transpose RDMs
integer :: j, k, l  ! Auxiliary variables for counters
real(8) :: t1, t2  ! Times
real(8) :: fidelity_pm, avgf  ! Fidelity function and its average
real(8) :: trace_he  ! Trace function for Hermitian matrices
real(8) :: coh_re, coh_l1n, avg_coh  ! For the quantum coherence functions and for its average
real(8) :: log2  ! For the base log function

call cpu_time(t1)

write(*,*) "## Performing some tests for the random density matrix generators"

! Sets the methods to be used in the generators
opt_rng = "mt" ;   opt_rpvg = "zhsl" ;   opt_rug = "gso" ;   opt_rsvg = "std" ;   opt_rdmg = "std"


write(*,*) "# Computing the trace and verifying positivity for some RDMs"
do j = 1, 5 ;   d = 2**j 
 allocate( rdm(1:d,1:d), egv(1:d) ) ;   call rdmg(d, rdm)
 call lapack_zheevd('N', d, rdm, egv) ;   ct = 0 ;   do k = 1, d ;   if ( egv(k) < 0.d0 ) ct = ct + 1 ;   enddo
 write(*,*) 'Tr(rho)=', trace_he(d, rdm), ",   No. of negative eigenvalues:", ct ;   deallocate(rdm, egv)
enddo


write(*,*) '# Computing the volume of PPT (separable) states'
open(unit = 13, file = 'rdm_ppt.dat', status = 'unknown')
! Gnuplot commands to see the results: plot [4:][0:] 1.8*exp(-0.26*x) w l t 'exact', 'rdm_ppt.dat' w p t 'numerics'
da = 2
do  d = 4, 30, 2 ;   db = d/da 
  allocate( rdm(1:d,1:d), rdm_pt(1:d,1:d), egv(1:d) )
  ct = 0
  do j = 1, ns
    call rdmg(d, rdm) ;   call partial_transpose_b(da, db, rdm, rdm_pt) ;   call lapack_zheevd('N', d, rdm_pt, egv)
    l = 0 ;   do k = 1, d ;   if ( egv(k) < 0.d0 ) l = l + 1 ;   enddo ! Verify if the RDM is PPT
    if ( l == 0 ) ct = ct + 1
  enddo
  write(13,*) d, dble(ct)/dble(ns) ;   write(*,*) d, dble(ct)/dble(ns)
  deallocate( rdm, rdm_pt, egv )
enddo
close(13)


write(*,*) '# Computing the average quantum coherence'
open(unit = 11, file = 'rdm_coh.dat', status = 'unknown')
! Gnuplot commands to see the results: 
do  d = 2, 30, 2 ;   allocate( rdm(1:d,1:d) ) ;   avg_coh = 0.d0
  do j = 1, ns
    call rdmg(d, rdm)
    !avg_coh = avg_coh + coh_l1n(d, rdm)
    avg_coh = avg_coh + coh_re(d, rdm)
  enddo
  !write(11,*) d, avg_coh/(dble(ns)*log2(dble(d))) ;   write(*,*) d, avg_coh/(dble(ns)*log2(dble(d)))
  write(11,*) d, avg_coh/dble(ns) ;   write(*,*) d, avg_coh/dble(ns)
  deallocate( rdm)
enddo
close(11)


write(*,*) '# Computing the average fidelity <F(psi,rho)>'
open(unit = 12, file = 'rdm_fid_ptr.dat', status = 'unknown')
! Gnuplot commands to see the results: plot 1/x w l t 'exact', 'rdm_fid.dat' w p t 'numerics'
do  d = 2, 30, 2 ;   allocate( rsv(1:d), rdm(1:d,1:d) ) ;   avgf = 0.d0
  do j = 1, ns
    call rsvg(d, rsv) ;   call rdmg(d, rdm) ;   avgf = avgf + fidelity_pm(d, rsv, rdm)
  enddo
  write(12,*) d, avgf/dble(ns) ;   write(*,*) d, avgf/dble(ns)
  deallocate( rsv, rdm)
enddo
close(12)


call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this subroutine

end
!###################################################################################################################################
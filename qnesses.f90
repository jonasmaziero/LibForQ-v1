!###################################################################################################################################
!                                    Tests quantumnesses measures using the two-qubit Werner state
!###################################################################################################################################
subroutine werner_2qb()  ! Tests for quantumnesses quantifiers using a two-qubit Werner's state
! Gnuplot command to see the results
! plot 'werner.dat' u 1:2 t 'Conc', '' u 1:3 t 'EoF', '' u 1:4 t 'Neg', '' u 1:5 t 'logNeg', '' u 1:6 t '2*Disc_hs', 
!      '' u 1:7 t 'Disc_oz_bds', '' u 1:8 t 'max(0,CHSH-2)', '' u 1:9 t 'coherence_l1n', '' u 1:10 t 'coherence_re'
use bell_basis
implicit none
complex(8) :: rho(1:4,1:4)  ! The 2-qubit Werner density matrix
complex(8) :: identity(1:4,1:4) ! The 4x4 identiy matrix
complex(8) :: proj(1:4,1:4)  ! For the projector on a Bell state 
real(8) :: w, dw  ! For the weight appearing in the Werner states
real(8) :: concurrence_2qb, EoF_2qb, negativity, log_negativity  ! Entanglement quantifiers
real(8) :: discord_hs_2qb, discord_oz_bds  ! Discord quantifiers
real(8) :: CHSH_2qb  ! Non-locality
real(8) :: coh_l1n, coh_re  ! Coherence quantifiers
open(unit = 11, file = 'werner.dat', status = 'unknown')

call identity_c(4, identity) ;   call projector(psi_m, 4, proj)

dw = 1.d0/100.d0 ;   w = 0.d0 - dw  ! Sets the width of the step in w
do
  w = w + dw ;   if( w > 1.d0 ) exit
  rho = (1.d0-w)*(identity/4.d0) + w*proj
  write(11,*) w, concurrence_2qb(rho), EoF_2qb(rho), negativity(2,2,'a',rho), log_negativity(2,2,'a',rho), & 
                 2.d0*discord_hs_2qb('a',rho), discord_oz_bds(rho), & 
                 max(CHSH_2qb(rho)-2.d0,0.d0), coh_l1n(4, rho), coh_re(4, rho)
enddo

end
!###################################################################################################################################
!                                                       Entanglement
!###################################################################################################################################
real(8) function concurrence_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
use pauli_group
implicit none
complex(8) :: rho(1:2**2,1:2**2)  ! Density matrix we want to compute the concurrence
complex(8) :: R(1:2**2,1:2**2), rho_tilde(1:2**2,1:2**2), s2_kp_s2(1:2**2,1:2**2)  ! Auxiliary matrices
complex(8) :: egv(1:2**2) ! Eigenvalues of R = rho*rho^tilde
real(8) :: egv_max  ! The greater eigenvalue of R

 call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, s2_kp_s2) ;   rho_tilde = matmul( matmul(s2_kp_s2,conjg(rho)) , s2_kp_s2 )
 R = matmul(rho,rho_tilde) ;   call lapack_zgeev('N', 4, R, egv)

 egv_max = max( real(egv(1)), real(egv(2)), real(egv(3)), real(egv(4)) )
 concurrence_2qb = max( 0.d0, (2.d0*sqrt(egv_max)-sqrt(real(egv(1)))-sqrt(real(egv(2)))-sqrt(real(egv(3)))-sqrt(real(egv(4)))))
     
end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function EoF_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(1:4,1:4)  ! Density matrix we want to compute the concurrence
real(8) :: concurrence_2qb  ! For the concurrence function
real(8) :: pv(1:2), shannon  ! Probability vector and Shannon's entropy

pv(1) = (1.d0 + sqrt(1.d0 - concurrence_2qb(rho)**2.d0))/2.d0 ;   pv(2) = 1.d0 - pv(1) ;   EoF_2qb = shannon(2, pv)

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function negativity(da, db, ssys, rho)  ! Returns the entanglement negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
implicit none
character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
real(8) :: norm_tr  ! For the trace norm function

if (ssys == 'a') then
  call partial_transpose_a(da, db, rho, rho_pt)
else if (ssys == 'b') then
  call partial_transpose_b(da, db, rho, rho_pt)
endif

negativity = norm_tr(da*db, rho_pt) - 1.d0
     
end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function log_negativity(da, db, ssys, rho)  ! Returns the entanglement negativity of a bipartite system
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
implicit none
character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: neg ! Negativity
real(8) :: norm_tr  ! For the trace norm function
real(8) :: log2  ! For the log base two

if (ssys == 'a') then
  call partial_transpose_a(da, db, rho, rho_pt)
else if (ssys == 'b') then
  call partial_transpose_b(da, db, rho, rho_pt)
endif

log_negativity = log2( norm_tr(da*db, rho_pt) )
     
end
!###################################################################################################################################
!                                                        Coherence
!###################################################################################################################################
real(8) function coh_l1n(d, rho)  ! Returns the l1-norm quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifyingcoherence, Phys. Rev. Lett. 113, 140401 (2014).
implicit none
integer :: d ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
real(8) :: norm_l1  ! For the l1-norm function

 coh_l1n = norm_l1(d, d, rho)

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function coh_re(d, rho)  ! Returns the relative entropy of quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifyingcoherence, Phys. Rev. Lett. 113, 140401 (2014).
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
real(8) :: shannon, neumann  ! For the Shannon and von Neumman entropy functions
real(8) :: pv(1:d)  ! Auxiliary probability vector
integer :: j  ! Auxiliary variable for counters

forall(j=1:d) pv(j) = real(rho(j,j)) ;   coh_re = shannon(d, pv) - neumann(d, rho)

end
!###################################################################################################################################
!                                                        Discord
!###################################################################################################################################
real(8) function discord_hs_2qb(ssys, rho)  ! Returns the 2-norm (or Hilbert-Schmidt) discord for 2-qubits states
! Ref: B. Dakić, V. Vedral e Č. Brukner, Necessary and sufficient condition for nonzero quantum discord, Phys. Rev. Lett. 105, 190502 (2010).
implicit none
character(1) :: ssys  ! Tells if sub-system a or b is classical (in the minimization)
complex(8) :: rho(1:4,1:4)  ! The density matrix
integer :: m, n  ! Auxiliary variables for counters
real(8) :: norm_r, norm_hs  ! For the vector norm and Hilbert-Schmidt norm functions
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: mK(1:3,1:3)  ! Auxiliary matrix
real(8) :: W(1:3)  ! For the eigenvalues of K
real(8) :: trace_he  ! For the trace function

call stokes_parameters_2qb(rho, ma, mb, corr)

if (ssys == 'a' ) then
  do m = 1, 3 ;   do n = 1, 3 ;   mK(m,n) = ma(m)*ma(n) ;   enddo ;   enddo ;   mK = mK + matmul(corr,transpose(corr))
  call lapack_dsyevd('N', 3, mK, W)
  discord_hs_2qb = 0.25d0*( (norm_r(3, ma))**2.d0 + (norm_hs(3, 3, corr))**2.d0 - maxval(W) )
else if ( ssys == 'b' ) then
  do m = 1, 3 ;   do n = 1, 3 ;   mK(m,n) = mb(m)*mb(n) ;   enddo ;   enddo ;   mK = mK + matmul(corr,transpose(corr))
  call lapack_dsyevd('N', 3, mK, W)
  discord_hs_2qb = 0.25d0*( (norm_r(3, mb))**2.d0 + (norm_hs(3, 3, corr))**2.d0 - maxval(W) )
endif

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_oz_bds(rho)  ! Returns the Ollivier-Zurek discord for 2-qubits Bell-diagonal states
! Ref: S. Luo, Quantum discord for two-qubit systems, Phys. Rev. A 77, 042303 (2008).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: l00, l01, l10, l11, MI, c, CC ! For mutual information and classical correlation, and related auxiliary variables
real(8) :: W(1:3)  ! For the eigenvalues of K
real(8) :: log2  ! For the log base two function

call stokes_parameters_2qb(rho, ma, mb, corr)  ! Get the Stokes parameters

! Eigenvalues of the Bell-diagonal two-qubit density matrix
l00 = ( 1.d0 + corr(1,1) - corr(2,2) + corr(3,3) )/4.d0 ;   l01 = ( 1.d0 + corr(1,1) + corr(2,2) - corr(3,3) )/4.d0
l10 = ( 1.d0 - corr(1,1) + corr(2,2) + corr(3,3) )/4.d0 ;   l11 = ( 1.d0 - corr(1,1) - corr(2,2) - corr(3,3) )/4.d0
! Mutual information (total correlation)
MI = l00*log2(4.d0*l00) + l01*log2(4.d0*l01) + l10*log2(4.d0*l10) + l11*log2(4.d0*l11)

! Classical correlation
 c = max( abs(corr(1,1)), abs(corr(2,2)), abs(corr(3,3)) ) ;   CC = 0.5d0*(1.d0-c)*log2(1.d0-c) + 0.5d0*(1.d0+c)*log2(1.d0+c)
 
! Ollivier-Zurek discord for Bell-diagonal states
discord_oz_bds = MI - CC

end
!###################################################################################################################################
!                                                      Non-Locality
!###################################################################################################################################
real(8) function CHSH_2qb(rho)  ! Returns the CHSH parameter (it max)
! Ref: R. Horodecki, P. Horodecki e M. Horodecki, “Violating Bell inequality by mixed spin-1 states: Necessary and sufficient 
!      condition”, Phys. Lett. A 200, 3402(1995).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: mK(1:3,1:3)  ! Auxiliary matrix
real(8) :: W(1:3)  ! For the eigenvalues of K

call stokes_parameters_2qb(rho, ma, mb, corr) ;   mK = matmul(transpose(corr),corr) ;   call lapack_dsyevd('N', 3, mK, W)
 CHSH_2qb = 2.d0*sqrt( W(2) + W(3) ) 

end
!###################################################################################################################################
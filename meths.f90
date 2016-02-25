!###################################################################################################################################
module meths
! This module is utilized to share the options for the methods to be used with the several related sub-programs. The conventions 
! are shown below. All these options can be changed at any time or place by using this module (with, use meths) in your sub-program. 
! But if you change the random number generator, you should use "call rng_init()" to initialize the new one (if necessary).
implicit none
character(10) :: opt_rng = "mt"
character(10) :: opt_rpvg = "zhsl"
character(10) :: opt_rug = "gso"
character(10) :: opt_rsvg = "std" 
character(10) :: opt_rdmg = "std"

! Below follows the conventions for the possible choices of methods

!! opt_rng      Selects the random number generator to be used
! = mt          (Mersenne Twister pseudo-rng)
! = gnu         (Gnu's pseudo-rng)
! = netlib      (Netlib's pseudo-rng)
! = ...         (...)

!! opt_rpvg     Selects the method to be used when generationg random probability vectors
! = zhsl        (zyczkowski-Horodecki-Sanpera-Lewenstein method)
! = norm        (normalization method)
! = trig        (trigonometric method)
! = iid         (iid method)
! = devroye     (Devroye's method)
! = kraemer     (Kraemer's method)
! = ...         (...)

!! opt_rug      Selects the method for unitary matrix generation
! = hhr         (Uses the QR factorization with Householder reflections)
! = gso         (Uses the QR factorization with modified Gram-Schmidt orthogonalization)
! = hurwitz     (Uses Hurwitz parametrization)
! = ...         (...)

!! opt_rsvg     Selects the method to be used in the generation of random pure quantum states (quantum state vectors)
! = "std"       (standard method |psi> = \sum_j \sqrt(p_j)\exp(i*theta_j)|c_j>)
! = "gauss"     (standard method |psi> = \sum_j c_j|c_j> with gaussianily distributed complex coefficients)
! = "ru"        (uses the first column of a randum unitary matrix)
! = ...         (...)

!! opt_rdmg     Selects the method to be used in the generation of random density matrices
! = "std"       (Standard method, rho = sum_j p_j U|c_j><c_j|U^†, where |c_j> is the computational basis)
! = "ginibre"   (Normalize G*G^†, with G being a Ginibre matrix)
! = "bures"     (Normalize (id+U)G*G^†(id+U^†), with G being a Ginibre matrix, U a random unitary, and id the identity matrix)
! = "ptrace"    (Takes the partial trace of a random state vector)
! = ...         (...)

end module
!###################################################################################################################################
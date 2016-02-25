!###################################################################################################################################
!                                          Random probability vector generators
!###################################################################################################################################
subroutine rpv_zhsl(d, rpv) ! Generates a random probability vector of dimension d using the Zyczkowski-Horodecki-Sanpera-Lewenstein method
! Zyczkowski, K., Horodecki, P., Sanpera, A., and Lewenstein, M. (1998). Volume of the set of separable states, Phys. Rev. A 58, 883. 
implicit none
integer :: d   ! Dimension of the probability vector
real(8) :: rpv(1:d)  ! Random probability vector 
real(8), allocatable :: rn(:)  ! Vector of random numbers
real(8) :: norm ! Normalization for the prv (auxiliary variable)
integer :: j  ! Auxiliary variable, counter for do

allocate( rn(1:d-1) ) ! Allocate memory for the vector of random numbers
call rng(d-1, rn) ! Call the choosed random number generator
 
rpv(1) = 1.d0 - (rn(1))**(1.d0/(dble(d-1))) ;  norm = rpv(1) ! These lines implement the ZHSL method
if ( d >= 3 ) then
  do j = 2, d-1 ;   rpv(j) = (1.d0 - norm)*(1.d0 - (rn(j))**(1.d0/(dble(d-j)))) ;   norm = norm + rpv(j) ;   enddo
endif
rpv(d) = 1.d0 - norm

deallocate( rn ) ! Deallocate the memory used with the vector of random numbers

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_norm(d, rpv)  ! Generates a random-unbiased probability vector of dimension d using the normalization method & shuffling
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d   ! Dimension of probability vector
real(8) :: rpv(1:d)  ! Random-unbiased probability vector
real(8), allocatable :: rn(:)  ! Vector of random numbers
real(8) :: norm  ! Normalization (auxiliary variable)
integer :: j  ! Counter for do (auxiliary variable)
integer, allocatable :: rperm(:)  ! vector for the random permutation of {1,2,...,d} 

allocate( rn(1:d-1) )  ! Allocate memory for the vector of random numbers
call rng(d-1, rn)  ! Call the choosed random number generator

allocate( rperm(1:d) )  ! Allocate memory for the vector of random permutation of {1,2,...,d}
call rand_perm(d, rperm)  ! Gets the random permutation 

rpv(rperm(1)) = rn(1) ;  norm = rpv(rperm(1))  ! These lines implement the normilization method with shuffling
if ( d >= 3 ) then
  do j = 2, d-1
    rpv(rperm(j)) = (1.d0 - norm)*rn(j) ;  norm = norm + rpv(rperm(j))
  enddo
endif
rpv(rperm(d)) = 1.d0 - norm

deallocate( rn, rperm ) ! Deallocate the memory used with these two variables

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_trig_b(d, rpv) ! Generates a random probability vector of dimension d using the trigonometric method
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d   ! Dimension of the random probability vector
real(8) :: rpv(1:d)  ! Random probability vector
real(8), allocatable :: rn(:)  ! Vector of random numbers
real(8), allocatable :: vec_theta(:)  ! Vector of angles used in this method
real(8) :: prod_cos  ! Store the product of the squared cossines
integer :: j, k  ! Auxiliary variables

allocate( rn(1:d-1), vec_theta(1:d-1) ) ! Allocate memory for the vector of random numbers and for the vector of angles
call rng(d-1, rn) ! Call the choosed random number generator

prod_cos = 1.d0  ! These lines implement the trigonometric method
do k = d-1, 1, -1
  vec_theta(k) = acos(sqrt(rn(d-k)))
  rpv(k+1) = ((sin(vec_theta(k)))**2.d0)*prod_cos
  prod_cos = prod_cos*((cos(vec_theta(k)))**2.d0)
enddo
rpv(1) = prod_cos

deallocate( rn, vec_theta )  ! Liberate the memory used by these variables

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_trig(d, rpv) ! Generates a random-unbiased probability vector of dimension d using the trigonometric method & shuffling
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d   ! Dimension of the random(-unbiased) probability vector
real(8), allocatable :: rpv_(:)  ! Random probability vector
real(8) :: rpv(1:d)  ! Random-unbiased probability vector 
integer :: j  ! Auxiliary variable
integer, allocatable :: rperm(:)  ! Vector for the random permutation of {1,2,...,d} 

allocate( rpv_(1:d) )  ! Allocates memory for the random probability vector
call rpv_trig_b(d, rpv_)  ! Gets the random probability vector via the trigonometric method

allocate( rperm(1:d) )  ! Allocates memory for the random permutation of {1,2,...,d}
call rand_perm(d, rperm)  ! Gets the random permutation 

forall( j = 1:d )  ! Shuffles the components of the random probability vector to avoid biasing
  rpv(j) = rpv_(rperm(j))
end forall

deallocate( rpv_, rperm ) ! Deallocates the memory used by these variables

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_iid(d, rpv) ! This subroutine generates and random probability vector of dimension d using the iid method
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d  ! Dimension of the random probability vector
real(8) :: rpv(1:d)  ! Random probability vector
real(8), allocatable :: rn(:)  ! Vector of random numbers

allocate( rn(1:d) ) ;   call rng(d, rn)  ! Allocate memory for and get the uniformly distributed random numbers
rpv = rn/sum(rn)   ! Divides each of the d independent random numbers by their sum
deallocate( rn )  ! Frees the memory used with this variable

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_devroye(d, rpv) ! Generates a random probability vector of dimension d using Devroye's method
! Ref: Devroye, L. (1986). Non-Uniform Random Variate Generation. New York: Springer.
implicit none
integer :: d  ! Dimension of the random probability vector
real(8) :: rpv(1:d)  ! Random probability vector
real(8), allocatable :: ern(:)  ! Vector of exponentially distributed random numbers

allocate( ern(1:d) ) ;  call rng_exp(d, ern)  ! Allocate memory for and get the exponentially distributed random numbers
rpv = ern/sum(ern)  ! This line implements Devroye's method
deallocate( ern )  ! Frees the memory used by this variable

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rpv_kraemer(d, rpv) ! Generates a random probability vector of dimension d using Kraemer's method 
! (in contrast to the others, this method uses sorting)
! Ref: Kraemer, H. (1999). Post on MathForum on December 20. Topic: Sampling uniformly from the n-simplex.
!use qsort_mod  ! Uses the modules which implements the Quicksort algorithm
implicit none
integer :: d  ! Dimension of the random probability vector
real(8) :: rpv(1:d)  ! Random probability vector 
real(8), allocatable :: rn(:)  ! Vector of random numbers
integer :: j  ! Auxiliary variable for counters

allocate( rn(1:d-1) ) ;   call rng(d-1, rn)  ! Allocates memory for and gets the vector of random numbers
call qsort(rn, d-1)  ! Sort the random numbers in non-decreasing order
rpv(1) = rn(1) ;  rpv(d) = 1.d0 - rn(d-1)  ! These two lines implement the Kraemer's method
forall ( j = 2:(d-1) ) rpv(j) = rn(j) - rn(j-1)

deallocate( rn ) ! Frees the memory used by this variable

end
!###################################################################################################################################
!                                              Calling subroutine (for the RPVG)
!###################################################################################################################################
subroutine rpvg(d, rpv) ! Calls the choosed random probability vector generator
use meths
implicit none
integer :: d
real(8) :: rpv(1:d)

     if ( opt_rpvg == "zhsl" ) then ;   call rpv_zhsl(d, rpv)
else if ( opt_rpvg == "norm" ) then ;   call rpv_norm(d, rpv)
else if ( opt_rpvg == "trig" ) then ;   call rpv_trig(d, rpv)
else if ( opt_rpvg == "iid" ) then ;   call rpv_iid(d, rpv)
else if ( opt_rpvg == "devroye" ) then ;   call rpv_devroye(d, rpv)
else if ( opt_rpvg == "kraemer" ) then ;   call rpv_kraemer(d, rpv)
     endif

end
!###################################################################################################################################
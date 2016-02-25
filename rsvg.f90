!###################################################################################################################################
!                                                     State vector generators - RSVG
!###################################################################################################################################
subroutine rsv_std(d, rsv)  ! Generates a random state vector (pure quantum state) using the standard method
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! Random state vector
real(8), allocatable :: rpv(:)  ! Random probability vector
real(8), allocatable :: rn(:)  ! Vector of random numbers
integer :: j  ! Auxiliary variable for counters
real(8), parameter :: pi = 4.d0*atan(1.d0)

allocate( rpv(1:d), rn(1:d) ) ;  call rpvg(d, rpv) ;   call rng(d, rn)  ! Allocate memory for and get rpv and rn
forall( j = 1 : d ) rsv(j) = sqrt(rpv(j))*exp( (0.d0,1.d0)*(rn(j)*2.d0*pi) )  ! The random phases theta_j = rn*2*pi are in [0,2*pi] 
deallocate( rpv, rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rsv_gauss(d, rsv) ! Generates a random state vector (pure quantum state) using the standard method with gaussianily distributed complex coefficients
! Ref: Maziero, J. (2015). Random sampling of quantum states: A survey of methods. Braz. J. Phys. 45, 575.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! (output) Random state vector
real(8), allocatable :: grn(:)  ! Random vector of gaussianily distributed random numbers
real(8) :: norm  ! Auxiliary variable for normalization
integer :: j  ! Auxiliary variable for counters

allocate( grn(1:2*d) ) ;   call rng_gauss(2*d, grn) ;   forall ( j = 1:d ) rsv(j) = grn(j) + (0.d0,1.d0)*grn(j+d)
rsv = rsv/norm(d, rsv)  ! Normalize the vector

deallocate( grn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rsv_ru(d, rsv) ! Generates a random quantum state vector using the first column of a random unitary matrix
! Ref: Zyczkowski, K. (1999). Volume of the set of separable states. II, Phys. Rev. A 60, 3496.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! (output) Random state vector
complex(8), allocatable :: ru(:,:)  ! Random unitary matrix
integer :: j  ! Auxiliary variable for counters

allocate( ru(1:d,1:d) ) ;  call rug(d, ru) ;   rsv(:) = ru(:,1) ;   deallocate( ru )

end
!###################################################################################################################################
!                                                        Calling subroutines - RQS
!###################################################################################################################################
subroutine rsvg(d, rsv)  ! Calls the choosed random state vector generator
use meths
implicit none
integer :: d
complex(8) :: rsv(1:d)

     if ( opt_rsvg == "std" ) then ;   call rsv_std(d, rsv)
else if ( opt_rsvg == "gauss" ) then ;   call rsv_gauss(d, rsv)
else if ( opt_rsvg == "ru" ) then ;   call rsv_ru(d, rsv)
     endif

end
!###################################################################################################################################
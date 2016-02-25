!###################################################################################################################################
!                                                   Random unitary matrix generators
!###################################################################################################################################
subroutine ru_gram_schmidt(d, ru)  ! Returns a dxd random unitary matrix by applying the Gram-Schmidt procedure to a random gaussian complex matrix 
! Ref: Diaconis, P. (2005). What is ... a random matrix?, Notices of the AMS 52, 1348.
implicit none
integer :: d  ! Dimension of the rando unitary
complex(8) :: ru(1:d,1:d)  ! Random unitary matrix
complex(8), allocatable :: A(:,:)  ! Auxiliary Ginibre matrix

allocate( A(1:d,1:d) ) ;   call ginibre(d, d, A) ;   call gram_schmidt_modified(d, d, A, ru) ;   deallocate( A )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine ru_householder(d, ru) ! Returns a dxd random unitary matrix using one of LAPACK's QR factorization subroutine (uses Hoseholder reflections)
! Ref: Mezzadri, F. (2007). How to generate random matrices from the classical compact groups, Notices of the AMS 54, 592.
implicit none
integer :: d  ! Dimension of the random unitary
complex(8) :: ru(1:d,1:d)  ! The random unitary
complex(8), allocatable :: Z(:,:), Q(:,:), R(:,:)  ! For the LAPACK QR factorization (Z = Q*R)

allocate( Z(1:d,1:d), Q(1:d,1:d), R(1:d,1:d) )
call ginibre(d, d, Z) ;   call lapack_zgeqrfp(d, Z, Q, R) ;   ru = Q
deallocate( Z, Q, R )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine ru_hurwitz(d, ru) ! Returns a dxd random unitary matrix from the Circular Unitary Ensemble (CUE) using Hurwitz's parametrization
! Ref: Zyczkowski, K., and Kus, M. (1994). Random unitary matrices, J. Phys. A: Math. Gen. 27, 4235.
implicit none
integer :: d  ! Dimension of the unitary matrix
complex(8) :: ru(1:d,1:d)  ! Random unitary matrix
complex(8), allocatable :: id(:,:), Ei(:,:), Ejk(:,:)  ! Identity and other auxiliary matrices
integer :: i, j, k  ! Auxiliary variables for the counters
real(8) :: alpha, psi, phi, chi, xi  ! Phases
real(8) :: rn(1:1)  ! Random number
real(8), parameter :: pi = 4.d0*datan(1.d0)

allocate( id(1:d,1:d), Ei(1:d,1:d), Ejk(1:d,1:d) ) ;  call identity_c(d, id)

do i = 1, d-1   ! Loop for the product of Ei: E_1*E_2*...*E_(d-1)
  k = i + 1  ! k is fixed, given i   
  do j = i, 1, -1  ! Loop for the products among the Ejk matrices, needed for obtaining Ei      
    Ejk = id
    call rng(1,rn) ;   xi = rn(1) ;   phi = dasin(xi**(1.d0/(2.d0*dble(j))))
    call rng(1,rn) ;   psi = 2.d0*pi*rn(1) 
    Ejk(j,j) = dcos(phi)*(dcos(psi) + (0.d0,1.d0)*dsin(psi)) ;   Ejk(k,k) = conjg(Ejk(j,j))
    if (j == 1) then
      call rng(1,rn) ;   chi = 2.d0*pi*rn(1)
      Ejk(j,k) = dsin(phi)*(dcos(chi) + (0.d0,1.d0)*dsin(chi)) ;   Ejk(k,j) = - conjg(Ejk(j,k))
    else if (j /= 1) then
      Ejk(j,k) = sin(phi) ;   Ejk(k,j) = - Ejk(j,k)
    endif
    if ( j == i ) then ;   Ei = Ejk ;   else if ( j /= i ) then ;  Ei(1:k,1:k) = matmul(Ei(1:k,1:k), Ejk(1:k,1:k)) ;  endif
  enddo
  if ( i == 1 ) then ;   ru = Ei ;   else if ( i /= 1 ) then ;  ru(1:k,1:k) = matmul(ru(1:k,1:k), Ei(1:k,1:k)) ;  endif
enddo

call rng(1,rn) ;   alpha = 2.d0*pi*rn(1)
do j = 1, d ;   do k = 1, d 
 if ( abs(ru(j,k)) > 1.d-15 ) ru(j,k) = (dcos(alpha) + (0.d0,1.d0)*dsin(alpha))*ru(j,k)
enddo ;   enddo

deallocate( id, Ei, Ejk)

end
!###################################################################################################################################
!                                                     Calling subroutines for RU
!###################################################################################################################################
subroutine rug(d, ru)  ! Calls the choosed random unitary matrix generator
use meths
implicit none
integer :: d  ! Dimension of the Random Unitary Matrix
complex(8) :: ru(1:d,1:d)  ! Random Unitary Matrix

     if ( opt_rug == "gso" ) then ;   call ru_gram_schmidt(d, ru)
else if ( opt_rug == "hhr" ) then ;   call ru_householder(d, ru)
else if ( opt_rug == "hurwitz" ) then ;   call ru_hurwitz(d, ru)
     endif

end
!###################################################################################################################################
!###################################################################################################################################
program main  ! You can change by your main program

! Initializes the random number generator. Call it again if you change the RNG afterwards. 
! Depending on the RNG, it is advisable doing this also if you generate a very large amount of RNs.
call rng_init()  

!call test_lapack() ;   call werner_2qb()  ! Calls the tests for the LAPACK eigensolvers

!! Uncomment the tests you want to perform and see the description in the correspondent subroutine
!call rng_tests()  ! Calls the tests for the random number generators
!call rpvg_tests()  ! Calls the tests for the random probability vectors generators
!call rug_tests()  ! Calls the tests for the random unitary generators
!call rsvg_tests()  ! Calls the tests for the random state vector generators
call rdmg_tests()  ! Calls the tests for the random density matrix generators

end
!###################################################################################################################################
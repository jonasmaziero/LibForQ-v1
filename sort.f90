!###################################################################################################################################
!                                                         Sorting
!###################################################################################################################################
!module qsort_mod
! Given nA real numbers in the vector A, it returns them sorted in non-decreasing order
! Adapted from the Rosetta Code Project: http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
! Remark: We modified it to allow choosing the random number generator
! Remark: qsort is used in Kraemer's method for the random probability vector generator.
!implicit none
! contains
recursive subroutine qsort(A, nA)
implicit none
! DUMMY ARGUMENTS
integer, intent(in) :: nA  ! Dimension of the list
real(8) :: A(1:nA)  ! List of real numbers to be sorted
 
! LOCAL VARIABLES
integer :: left, right
real(8) :: random(1)
real(8) :: pivot
real(8) :: temp
integer :: marker
 
    if (nA > 1) then
 
        call rng(1,random)
        pivot = A(int(random(1)*real(nA-1))+1)   ! random pivot (not best performance, but avoids worst-case)
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right) > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left) < pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call QSort(A(:marker-1), marker-1)
        call QSort(A(marker:), nA-marker+1)
 
    end if
 
end subroutine qsort 
!end module qsort_mod
!###################################################################################################################################
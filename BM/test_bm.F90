!Author: Rahul Kashyap (rahulkashyap@iitb.ac.in)
!Date: 2025-08-06
!Description: Fortran program to test the Blandford-Mckee solution for a blast wave.
!This program uses the module defined in blandford_mckee_module.f90
program test_bm
    use blandford_mckee_module
    implicit none
    integer, parameter :: n = 100
    real(8) :: r_array(n), p(n), g(n), nprime(n)
    real(8) :: t
    integer :: i
    character(len=*), parameter :: filename = "bm_solution.txt"
    integer :: unit

    ! Define time and radius array
    t = 0.01d0
    do i = 1, n
         r_array(i) = 0.0001d0 + (1.0d0 - 0.0d0) * real(i-1)/(n-1)
    end do

    call solve_blandford_mckee(r_array, n, t, p, g, nprime)

    open(newunit=unit, file=filename, status='replace', action='write')
    do i = 1, n
         write(unit, '(4ES20.10)') r_array(i), p(i), g(i), nprime(i)
    end do
    close(unit)
end program test_bm
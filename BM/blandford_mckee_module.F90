!Author: Rahul Kashyap (rahulkashyap@iitb.ac.in)
!Date: 2025-08-06
!Description: This module is used in test_bm.F90
module blandford_mckee_module
  implicit none
  private
  public :: solve_blandford_mckee

contains

  subroutine solve_blandford_mckee(r_array, n_points, t, pressure, lorentz, density)
    implicit none
    ! Input
    integer, intent(in) :: n_points
    real(8), intent(in) :: r_array(n_points)
    real(8), intent(in) :: t

    ! Output
    real(8), intent(out) :: pressure(n_points)
    real(8), intent(out) :: lorentz(n_points)
    real(8), intent(out) :: density(n_points)

    ! Parameters
    real(8), parameter :: Gamma = 100.0d0
    real(8), parameter :: m = 3.0d0
    real(8), parameter :: w1 = 1.0d0
    real(8), parameter :: n1 = 1.0d0
    real(8), parameter :: chi_min = 1.0d0
    real(8), parameter :: chi_max = 1.0d0 + 2.0d0*(m + 1.0d0)*5.0d0
    integer, parameter :: N_chi = 1000

    ! Internal arrays
    real(8) :: chi_vals(N_chi), f_vals(N_chi), g_vals(N_chi), h_vals(N_chi)
    real(8) :: dchi, chi
    integer :: i

    ! Temporary variables
    real(8) :: k1(3), k2(3), k3(3), k4(3), y(3), chi_mid
    real(8) :: xi, chi_interp, f, g, h
    real(8) :: r_target

    ! Initialize chi array
    dchi = (chi_max - chi_min) / real(N_chi - 1)
    chi_vals(1) = chi_min
    f_vals(1) = 1.0d0
    g_vals(1) = 1.0d0
    h_vals(1) = 1.0d0

    ! Integrate using Runge-Kutta 4
    y = (/ f_vals(1), g_vals(1), h_vals(1) /)
    do i = 2, N_chi
      chi = chi_vals(i-1)
      call bm_rhs(chi, y, k1)
      call bm_rhs(chi + 0.5d0*dchi, y + 0.5d0*dchi*k1, k2)
      call bm_rhs(chi + 0.5d0*dchi, y + 0.5d0*dchi*k2, k3)
      call bm_rhs(chi + dchi, y + dchi*k3, k4)

      y = y + (dchi/6.0d0)*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
      chi_vals(i) = chi_vals(i-1) + dchi
      f_vals(i) = y(1)
      g_vals(i) = y(2)
      h_vals(i) = y(3)
    end do

    ! For each radius, compute corresponding chi and interpolate f, g, h
    do i = 1, n_points
      r_target = r_array(i)
      xi = Gamma**2 * (1.0d0 - r_target / t)
      chi_interp = 1.0d0 + 2.0d0*(m + 1.0d0)*xi

      ! Simple linear interpolation in chi_vals
      call interpolate_chi(chi_vals, f_vals, N_chi, chi_interp, f)
      call interpolate_chi(chi_vals, g_vals, N_chi, chi_interp, g)
      call interpolate_chi(chi_vals, h_vals, N_chi, chi_interp, h)

      pressure(i) = (2.0d0/3.0d0) * w1 * Gamma**2 * f
      lorentz(i)  = sqrt(0.5d0 * Gamma**2 * g)
      density(i)  = 2.0d0 * n1 * Gamma**2 * h
    end do

  end subroutine solve_blandford_mckee

  !----------------------------------------------------------

  subroutine bm_rhs(chi, y, dydchi)
    implicit none
    real(8), intent(in) :: chi, y(3)
    real(8), intent(out) :: dydchi(3)
    real(8), parameter :: m = 3.0d0
    real(8) :: f, g, h, denom, denom_h
    real(8) :: dlnf, dlng, dlnh

    f = y(1)
    g = y(2)
    h = y(3)

    denom = (m + 1.0d0)*(4.0d0 - 8.0d0*g*chi + g**2 * chi**2)
    dlnf = (8.0d0*(m - 1.0d0) - (m - 4.0d0)*g*chi) / denom
    dlng = ((7.0d0*m - 4.0d0) - (m + 2.0d0)*g*chi) / denom

    denom_h = denom * (2.0d0 - g*chi)
    dlnh = (2.0d0*(9.0d0*m - 8.0d0) - 2.0d0*(5.0d0*m - 6.0d0)*g*chi + (m - 2.0d0)*g**2 * chi**2) / denom_h

    dydchi(1) = g*f * dlnf
    dydchi(2) = g*g * dlng
    dydchi(3) = g*h * dlnh

  end subroutine bm_rhs

  !----------------------------------------------------------

  !----------------------------------------------------------------------
  ! interpolate_chi:
  !   Performs linear interpolation to estimate y_interp at x using arrays
  !   x_vals and y_vals of length n.
  !   Parameters:
  !     x_vals(n)   - array of x values (must be sorted ascending)
  !     y_vals(n)   - array of y values corresponding to x_vals
  !     n           - number of points in arrays
  !     x           - target x value for interpolation
  !     y_interp    - interpolated y value at x
  !----------------------------------------------------------------------
  subroutine interpolate_chi(x_vals, y_vals, n, x, y_interp)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x_vals(n), y_vals(n), x
    real(8), intent(out) :: y_interp
    integer :: i

    ! Clamp to bounds
    if (x <= x_vals(1)) then
      y_interp = y_vals(1)
    else if (x >= x_vals(n)) then
      y_interp = y_vals(n)
    else
      do i = 1, n - 1
        if (x >= x_vals(i) .and. x < x_vals(i+1)) then
          if (x_vals(i+1) == x_vals(i)) then
            y_interp = y_vals(i)
          else
            y_interp = y_vals(i) + (x - x_vals(i)) * (y_vals(i+1) - y_vals(i)) / (x_vals(i+1) - x_vals(i))
          end if
          exit
        end if
      end do
    end if
  end subroutine interpolate_chi

end module blandford_mckee_module
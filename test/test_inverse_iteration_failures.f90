! Analytical tests for mathematical and numerical inverse-iteration limits.
!
! These cases distinguish recoverable library errors from inherent properties
! of shifted inverse iteration. Degenerate eigenspaces may return any retained
! direction, absent starting-vector components cannot be created, equidistant
! shifts can oscillate, close clusters converge slowly, singular shifts cannot
! be solved, and a non-positive overlap norm violates the S-positive-definite
! precondition.
program test_inverse_iteration_failures
  use qrlinalg
  use test_support, only: check, finish_suite
  implicit none

  integer :: failures

  failures = 0
  call test_degenerate_eigenspaces(failures)
  call test_equidistant_shift_oscillation(failures)
  call test_missing_target_component(failures)
  call test_slow_cluster(failures)
  call test_precision_scale_singular_shift(failures)
  call test_invalid_overlap_and_tiny_start(failures)
  call finish_suite('qrlinalg inverse-iteration failure regimes', failures)

contains

  ! An exactly repeated eigenvalue determines an eigenspace rather than one
  ! eigenvector. Inverse iteration must retain the starting vector's projection
  ! into that eigenspace, suppress components belonging to other eigenvalues,
  ! and return a valid eigenpair for both real and complex arithmetic.
  subroutine test_degenerate_eigenspaces(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 3
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp), parameter :: shift = 2.5_wp, exact_lambda = 2.0_wp
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    real(wp) :: real_expected(n), real_residual_vector(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n), complex_expected(n)
    complex(wp) :: complex_residual_vector(n)
    real(wp) :: lambda, rel_acc, tolerance, eigen_tolerance
    real(wp) :: direction_overlap, residual
    integer :: info, num_iter

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    eigen_tolerance = 100000.0_wp * epsilon(1.0_wp)
    real_s = 0.0_wp
    real_h = 0.0_wp
    real_s(1,1) = 1.0_wp
    real_s(2,2) = 1.0_wp
    real_s(3,3) = 1.0_wp
    real_h(1,1) = exact_lambda
    real_h(2,2) = exact_lambda
    real_h(3,3) = 5.0_wp
    real_v = [1.0_wp, 2.0_wp, 1.0_wp]
    real_expected = [1.0_wp, 2.0_wp, 0.0_wp] / sqrt(5.0_wp)

    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, shift, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, 200, 1, &
                          rel_acc, num_iter, info)
    real_residual_vector = matmul(real_h, real_x) - &
                           lambda * matmul(real_s, real_x)
    residual = sqrt(sum(real_residual_vector**2))
    direction_overlap = abs(dot_product(real_expected, real_x))
    call check(info == QR_SUCCESS .and. num_iter < 200, &
               'real repeated eigenvalue converges to its eigenspace', failures)
    call check(abs(lambda - exact_lambda) <= eigen_tolerance .and. &
               residual <= eigen_tolerance, &
               'real degenerate result is an analytical eigenpair', failures)
    call check(abs(direction_overlap - 1.0_wp) <= eigen_tolerance, &
               'real degeneracy retains the initial eigenspace projection', &
               failures)

    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_v = [cmplx(1.0_wp, 1.0_wp, kind=wp), &
                 cmplx(2.0_wp, -1.0_wp, kind=wp), &
                 cmplx(1.0_wp, 0.5_wp, kind=wp)]
    complex_expected = [complex_v(1), complex_v(2), &
                        cmplx(0.0_wp, 0.0_wp, kind=wp)]
    complex_expected = complex_expected / sqrt(sum(abs(complex_expected)**2))

    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, shift, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, 200, 1, rel_acc, num_iter, info)
    complex_residual_vector = matmul(complex_h, complex_x) - &
      cmplx(lambda, 0.0_wp, kind=wp) * matmul(complex_s, complex_x)
    residual = sqrt(sum(abs(complex_residual_vector)**2))
    direction_overlap = abs(dot_product(complex_expected, complex_x))
    call check(info == QR_SUCCESS .and. num_iter < 200, &
               'complex repeated eigenvalue converges to its eigenspace', &
               failures)
    call check(abs(lambda - exact_lambda) <= eigen_tolerance .and. &
               residual <= eigen_tolerance, &
               'complex degenerate result is an analytical eigenpair', failures)
    call check(abs(direction_overlap - 1.0_wp) <= eigen_tolerance, &
               'complex degeneracy retains its projected direction up to phase', &
               failures)
  end subroutine test_degenerate_eigenspaces

  ! A shift exactly midway between two eigenvalues gives inverse factors with
  ! equal magnitude and opposite sign. An equal-weight start alternates between
  ! two directions and cannot satisfy the directional convergence test.
  subroutine test_equidistant_shift_oscillation(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2, max_iter = 8
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp), parameter :: shift = 2.0_wp
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n)
    real(wp) :: lambda, rel_acc, tolerance, norm_of_x
    integer :: info, num_iter

    tolerance = 100.0_wp * epsilon(1.0_wp)
    real_h = 0.0_wp
    real_s = 0.0_wp
    real_h(1,1) = 1.0_wp
    real_h(2,2) = 3.0_wp
    real_s(1,1) = 1.0_wp
    real_s(2,2) = 1.0_wp
    real_v = [1.0_wp, 1.0_wp]
    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, shift, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, &
                          max_iter, 1, rel_acc, num_iter, info)
    norm_of_x = sqrt(sum(real_x**2))
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == max_iter, &
               'real midpoint shift reports persistent oscillation', failures)
    call check(rel_acc > 0.5_wp .and. &
               abs(lambda - shift) <= 1000.0_wp * epsilon(1.0_wp) .and. &
               abs(norm_of_x - 1.0_wp) <= 1000.0_wp * epsilon(1.0_wp), &
               'real oscillation returns a finite normalized approximation', &
               failures)

    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = [cmplx(1.0_wp, 1.0_wp, kind=wp), &
                 cmplx(1.0_wp, -1.0_wp, kind=wp)]
    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, shift, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, max_iter, 1, rel_acc, num_iter, info)
    norm_of_x = sqrt(sum(abs(complex_x)**2))
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == max_iter, &
               'complex midpoint shift reports persistent oscillation', failures)
    call check(rel_acc > 0.5_wp .and. &
               abs(lambda - shift) <= 1000.0_wp * epsilon(1.0_wp) .and. &
               abs(norm_of_x - 1.0_wp) <= 1000.0_wp * epsilon(1.0_wp), &
               'complex oscillation returns a finite normalized approximation', &
               failures)
  end subroutine test_equidistant_shift_oscillation

  ! Inverse iteration cannot create a component in an eigenvector that is
  ! exactly absent from the start. The solver should therefore converge
  ! successfully to the invariant subspace supplied by the caller, even when
  ! another eigenvalue is closer to the shift.
  subroutine test_missing_target_component(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp), parameter :: shift = 1.25_wp
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n)
    real(wp) :: lambda, rel_acc, tolerance
    integer :: info, num_iter

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    real_h = 0.0_wp
    real_s = 0.0_wp
    real_h(1,1) = 1.0_wp
    real_h(2,2) = 4.0_wp
    real_s(1,1) = 1.0_wp
    real_s(2,2) = 1.0_wp
    real_v = [0.0_wp, 1.0_wp]
    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, shift, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, 20, 1, &
                          rel_acc, num_iter, info)
    call check(info == QR_SUCCESS .and. abs(lambda - 4.0_wp) <= tolerance .and. &
               abs(abs(real_x(2)) - 1.0_wp) <= tolerance .and. &
               abs(real_x(1)) <= tolerance, &
               'real start confined to wrong eigenspace stays there', failures)

    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = [cmplx(0.0_wp, 0.0_wp, kind=wp), &
                 cmplx(0.0_wp, 1.0_wp, kind=wp)]
    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, shift, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, 20, 1, rel_acc, num_iter, info)
    call check(info == QR_SUCCESS .and. abs(lambda - 4.0_wp) <= tolerance .and. &
               abs(abs(complex_x(2)) - 1.0_wp) <= tolerance .and. &
               abs(complex_x(1)) <= tolerance, &
               'complex start confined to wrong eigenspace stays there', &
               failures)
  end subroutine test_missing_target_component

  ! Two eigenvalues separated by sqrt(epsilon) are distinguishable in storage
  ! but converge extremely slowly when the shift is far from both. A short
  ! iteration limit must return nonconvergence rather than a false success.
  subroutine test_slow_cluster(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2, max_iter = 8
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n)
    real(wp) :: cluster_gap, lambda, rel_acc, tolerance, norm_of_x
    integer :: info, num_iter

    cluster_gap = sqrt(epsilon(1.0_wp))
    tolerance = 100.0_wp * epsilon(1.0_wp)
    real_h = 0.0_wp
    real_s = 0.0_wp
    real_h(1,1) = 1.0_wp
    real_h(2,2) = 1.0_wp + cluster_gap
    real_s(1,1) = 1.0_wp
    real_s(2,2) = 1.0_wp
    real_v = [1.0_wp, 1.0_wp]
    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, 0.0_wp, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, &
                          max_iter, 1, rel_acc, num_iter, info)
    norm_of_x = sqrt(sum(real_x**2))
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == max_iter .and. &
               rel_acc > tolerance, &
               'real close cluster does not report premature convergence', &
               failures)
    call check(lambda >= 1.0_wp .and. lambda <= 1.0_wp + cluster_gap .and. &
               abs(norm_of_x - 1.0_wp) <= 1000.0_wp * epsilon(1.0_wp), &
               'real clustered result remains finite and normalized', failures)

    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = [cmplx(1.0_wp, 1.0_wp, kind=wp), &
                 cmplx(1.0_wp, -1.0_wp, kind=wp)]
    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.0_wp, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, max_iter, 1, rel_acc, num_iter, info)
    norm_of_x = sqrt(sum(abs(complex_x)**2))
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == max_iter .and. &
               rel_acc > tolerance, &
               'complex close cluster does not report premature convergence', &
               failures)
    call check(lambda >= 1.0_wp .and. lambda <= 1.0_wp + cluster_gap .and. &
               abs(norm_of_x - 1.0_wp) <= 1000.0_wp * epsilon(1.0_wp), &
               'complex clustered result remains finite and normalized', &
               failures)
  end subroutine test_slow_cluster

  ! xGEQRF accepts rank-deficient and extremely ill-conditioned matrices.
  ! Place the shift one epsilon from an eigenvalue while another diagonal sets
  ! an O(1) factor scale. The explicit precision-scaled R test must reject the
  ! triangular solve before any iteration is attempted.
  subroutine test_precision_scale_singular_shift(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n)
    real(wp) :: shift, lambda, rel_acc
    integer :: info, num_iter

    shift = 1.0_wp + epsilon(1.0_wp)
    real_h = 0.0_wp
    real_s = 0.0_wp
    real_h(1,1) = 1.0_wp
    real_h(2,2) = 4.0_wp
    real_s(1,1) = 1.0_wp
    real_s(2,2) = 1.0_wp
    real_v = [1.0_wp, 1.0_wp]
    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, shift, info)
    call real_state%solve(real_s, real_v, real_x, lambda, &
                          100.0_wp * epsilon(1.0_wp), 20, 1, &
                          rel_acc, num_iter, info)
    call check(info == QR_ERR_SINGULAR .and. num_iter == 0 .and. &
               all(abs(real_x) <= 0.0_wp), &
               'real precision-scale singular shift is rejected before solve', &
               failures)

    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = cmplx(real_v, real_v, kind=wp)
    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, shift, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             100.0_wp * epsilon(1.0_wp), 20, 1, &
                             rel_acc, num_iter, info)
    call check(info == QR_ERR_SINGULAR .and. num_iter == 0 .and. &
               all(abs(complex_x) <= 0.0_wp), &
               'complex precision-scale singular shift is rejected before solve', &
               failures)
  end subroutine test_precision_scale_singular_shift

  ! S must be positive definite and the starting norm must exceed the library's
  ! tiny threshold. An indefinite overlap can produce x^H*S*x<=0 after a valid
  ! iteration; a roundoff-scale start is rejected before BLAS or factor use.
  subroutine test_invalid_overlap_and_tiny_start(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_real_state) :: real_state
    type(qr_complex_state) :: complex_state
    real(wp) :: real_h(n,n), real_s(n,n), real_v(n), real_x(n)
    complex(wp) :: complex_h(n,n), complex_s(n,n)
    complex(wp) :: complex_v(n), complex_x(n)
    real(wp) :: lambda, rel_acc, tolerance, tiny_component
    integer :: info, num_iter

    tolerance = 100.0_wp * epsilon(1.0_wp)
    real_h = 0.0_wp
    real_s = 0.0_wp
    real_h(1,1) = 2.0_wp
    real_h(2,2) = 4.0_wp
    real_s(1,1) = 1.0_wp
    real_s(2,2) = -1.0_wp
    real_v = [0.0_wp, 1.0_wp]
    call real_state%initialize(n, info)
    call real_state%factorize_fresh(real_h, real_s, 0.0_wp, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, 20, 1, &
                          rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 1 .and. &
               abs(lambda) <= 0.0_wp .and. all(abs(real_x) < huge(1.0_wp)), &
               'real non-positive overlap norm is reported after iteration', &
               failures)

    complex_h = cmplx(real_h, 0.0_wp, kind=wp)
    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = [cmplx(0.0_wp, 0.0_wp, kind=wp), &
                 cmplx(0.0_wp, 1.0_wp, kind=wp)]
    call complex_state%initialize(n, info)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.0_wp, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, 20, 1, rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 1 .and. &
               abs(lambda) <= 0.0_wp .and. &
               all(abs(complex_x) < huge(1.0_wp)), &
               'complex non-positive overlap norm is reported after iteration', &
               failures)

    ! The squared norm is deliberately below tiny(1.0_wp), exercising the
    ! numerical threshold rather than the exact-zero shortcut.
    tiny_component = 0.5_wp * sqrt(tiny(1.0_wp))
    real_s(2,2) = 1.0_wp
    real_v = [tiny_component, 0.0_wp]
    call real_state%factorize_fresh(real_h, real_s, 0.0_wp, info)
    call real_state%solve(real_s, real_v, real_x, lambda, tolerance, 20, 1, &
                          rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 0 .and. &
               all(abs(real_x) <= 0.0_wp), &
               'real roundoff-scale starting vector is rejected', failures)

    complex_s = cmplx(real_s, 0.0_wp, kind=wp)
    complex_v = [cmplx(tiny_component, 0.0_wp, kind=wp), &
                 cmplx(0.0_wp, 0.0_wp, kind=wp)]
    call complex_state%factorize_fresh(complex_h, complex_s, 0.0_wp, info)
    call complex_state%solve(complex_s, complex_v, complex_x, lambda, &
                             tolerance, 20, 1, rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 0 .and. &
               all(abs(complex_x) <= 0.0_wp), &
               'complex roundoff-scale starting vector is rejected', failures)
  end subroutine test_invalid_overlap_and_tiny_start

end program test_inverse_iteration_failures

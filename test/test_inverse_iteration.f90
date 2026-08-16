! Analytical inverse-iteration and solve-error tests for qrlinalg.
!
! The generalized eigenproblems have exact, noncommuting H and S matrices so
! the tests exercise the complete QR-backed GSEPIIS and GHEPIIS paths.
program test_inverse_iteration
  use qrlinalg
  use test_support, only: check, finish_suite, real_frobenius, &
                          complex_frobenius
  implicit none

  integer :: failures

  failures = 0
  call test_real_inverse_iteration(failures)
  call test_complex_inverse_iteration(failures)
  call test_solve_error_paths(failures)
  call finish_suite('qrlinalg inverse iteration', failures)

contains

  ! Test the QR-backed GSEPIIS path on a noncommuting generalized problem with
  ! an exact eigendecomposition. For
  !
  !   X = [1 1; 0 1],  S=X^(-T)*X^(-1),  H=X^(-T)*diag(2,6)*X^(-1),
  !
  ! columns of X are S-orthonormal generalized eigenvectors with eigenvalues
  ! 2 and 6. A shift of 2.25 is closest to x1=(1,0), which already has unit
  ! S norm. H and S do not commute, preventing the test from collapsing into
  ! an ordinary eigenproblem in a shared orthogonal basis.
  subroutine test_real_inverse_iteration(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_real_state) :: state
    real(wp), parameter :: shift = 2.25_wp, exact_lambda = 2.0_wp
    real(wp) :: h(n,n), s(n,n), s_before(n,n), v(n), v_before(n), x(n)
    real(wp) :: lambda, rel_acc, tolerance, eigen_tolerance
    real(wp) :: converged_lambda, converged_rel_acc
    real(wp) :: residual, s_norm, euclidean_norm
    integer :: info, num_iter, converged_num_iter

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    eigen_tolerance = 10000.0_wp * epsilon(1.0_wp)
    s = reshape([1.0_wp, -1.0_wp, -1.0_wp, 2.0_wp], [n,n])
    h = reshape([2.0_wp, -2.0_wp, -2.0_wp, 8.0_wp], [n,n])
    v = [1.0_wp, 0.5_wp]
    s_before = s
    v_before = v

    call state%initialize(n, info)
    call check(info == QR_SUCCESS, &
               'real inverse-iteration state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, &
               'real inverse-iteration matrix factorizes', failures)

    call state%solve(s, v, x, lambda, tolerance, 100, 0, &
                     rel_acc, num_iter, info)
    residual = sqrt(sum((matmul(h, x) - lambda * matmul(s, x))**2)) / &
               ((real_frobenius(h) + abs(lambda) * real_frobenius(s)) * &
                sqrt(sum(x * x)))
    s_norm = dot_product(x, matmul(s, x))
    converged_lambda = lambda
    converged_rel_acc = rel_acc
    converged_num_iter = num_iter
    call check(info == QR_SUCCESS, &
               'real inverse iteration converges', failures)
    call check(num_iter > 0 .and. num_iter < 100 .and. rel_acc <= tolerance, &
               'real inverse iteration reports converged iteration data', &
               failures)
    call check(abs(lambda - exact_lambda) <= eigen_tolerance, &
               'real inverse iteration finds the analytical eigenvalue', &
               failures)
    call check(maxval(abs(abs(x) - [1.0_wp, 0.0_wp])) <= &
               eigen_tolerance, &
               'real S-normalized vector matches analytical magnitudes', &
               failures)
    call check(abs(s_norm - 1.0_wp) <= eigen_tolerance, &
               'real norm mode zero produces unit S norm', failures)
    call check(residual <= eigen_tolerance, &
               'real analytical generalized residual is small', failures)
    call check(all(abs(s - s_before) <= 0.0_wp) .and. &
               all(abs(v - v_before) <= 0.0_wp), &
               'real solve preserves S and the initial vector', failures)

    ! Exercise the remaining normalization modes using the same factors.
    call state%solve(s, v, x, lambda, tolerance, 100, 1, &
                     rel_acc, num_iter, info)
    euclidean_norm = sqrt(sum(x * x))
    call check(info == QR_SUCCESS .and. &
               abs(euclidean_norm - 1.0_wp) <= eigen_tolerance, &
               'real norm mode one produces unit Euclidean norm', failures)

    call state%solve(s, v, x, lambda, tolerance, 100, 2, &
                     rel_acc, num_iter, info)
    call check(info == QR_SUCCESS .and. &
               abs(maxval(abs(x)) - 1.0_wp) <= eigen_tolerance, &
               'real other norm modes retain largest-component scaling', &
               failures)

    ! GSEPIIS returns its best approximation with a distinct status when the
    ! iteration limit is exhausted. The QR-backed API preserves that contract.
    call state%solve(s, v, x, lambda, epsilon(1.0_wp), 1, 1, &
                     rel_acc, num_iter, info)
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == 1, &
               'real solve reports an exhausted iteration limit', failures)
    call check(abs(sqrt(sum(x * x)) - 1.0_wp) <= eigen_tolerance, &
               'nonconverged real solve still normalizes its best vector', &
               failures)

    ! A negative tolerance cannot accept the first apparently accurate step;
    ! it waits for the direction error to turn upward. On this monotonically
    ! converging case, a two-iteration limit must therefore be exhausted.
    call state%solve(s, v, x, lambda, -1.0_wp, 2, 1, &
                     rel_acc, num_iter, info)
    call check(info == QR_ERR_NO_CONVERGENCE .and. num_iter == 2, &
               'negative tolerance requires a confirming iteration', failures)

    write(*,'(a,i0,3(a,es12.4),a,i0)') '  real solve wp=', wp, &
      ' lambda=', converged_lambda, ' residual=', residual, &
      ' rel_acc=', converged_rel_acc, ' iterations=', converged_num_iter
  end subroutine test_real_inverse_iteration

  ! Hermitian noncommuting counterpart. For X=[1 i; 0 1],
  !
  !   S=X^(-H)*X^(-1),       H=X^(-H)*diag(1,5)*X^(-1).
  !
  ! The columns x1=(1,0) and x2=(i,1) are S-orthonormal generalized
  ! eigenvectors with eigenvalues 1 and 5. This exercises conjugate
  ! transposes, Hermitian S multiplication, phase-insensitive comparison, and
  ! all complex normalization modes without reducing to commuting matrices.
  subroutine test_complex_inverse_iteration(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_complex_state) :: state
    real(wp), parameter :: shift = 1.25_wp, exact_lambda = 1.0_wp
    complex(wp) :: h(n,n), s(n,n), s_before(n,n)
    complex(wp) :: v(n), v_before(n), x(n), residual_vector(n)
    real(wp) :: lambda, rel_acc, tolerance, eigen_tolerance
    real(wp) :: converged_lambda, converged_rel_acc
    real(wp) :: residual, s_norm, euclidean_norm, component_scale
    integer :: info, num_iter, converged_num_iter

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    eigen_tolerance = 10000.0_wp * epsilon(1.0_wp)
    s(1,1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    s(1,2) = cmplx(0.0_wp, -1.0_wp, kind=wp)
    s(2,1) = conjg(s(1,2))
    s(2,2) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    h(1,1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    h(1,2) = cmplx(0.0_wp, -1.0_wp, kind=wp)
    h(2,1) = conjg(h(1,2))
    h(2,2) = cmplx(6.0_wp, 0.0_wp, kind=wp)
    v = [cmplx(1.0_wp, 0.0_wp, kind=wp), &
         cmplx(1.0_wp, 0.0_wp, kind=wp)]
    s_before = s
    v_before = v

    call state%initialize(n, info)
    call check(info == QR_SUCCESS, &
               'complex inverse-iteration state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, &
               'complex inverse-iteration matrix factorizes', failures)

    call state%solve(s, v, x, lambda, tolerance, 100, 0, &
                     rel_acc, num_iter, info)
    residual_vector = matmul(h, x) - &
                      cmplx(lambda, 0.0_wp, kind=wp) * matmul(s, x)
    residual = sqrt(sum(abs(residual_vector)**2)) / &
               ((complex_frobenius(h) + abs(lambda) * complex_frobenius(s)) * &
                sqrt(sum(abs(x)**2)))
    s_norm = real(dot_product(x, matmul(s, x)), wp)
    converged_lambda = lambda
    converged_rel_acc = rel_acc
    converged_num_iter = num_iter
    call check(info == QR_SUCCESS, &
               'complex inverse iteration converges', failures)
    call check(num_iter > 0 .and. num_iter < 100 .and. rel_acc <= tolerance, &
               'complex inverse iteration reports converged iteration data', &
               failures)
    call check(abs(lambda - exact_lambda) <= eigen_tolerance, &
               'complex inverse iteration finds the analytical eigenvalue', &
               failures)
    call check(maxval(abs(abs(x) - [1.0_wp, 0.0_wp])) <= &
               eigen_tolerance, &
               'complex S-normalized vector matches analytical magnitudes', &
               failures)
    call check(abs(s_norm - 1.0_wp) <= eigen_tolerance, &
               'complex norm mode zero produces unit S norm', failures)
    call check(residual <= eigen_tolerance, &
               'complex analytical generalized residual is small', failures)
    call check(all(abs(s - s_before) <= 0.0_wp) .and. &
               all(abs(v - v_before) <= 0.0_wp), &
               'complex solve preserves S and the initial vector', failures)

    call state%solve(s, v, x, lambda, tolerance, 100, 1, &
                     rel_acc, num_iter, info)
    euclidean_norm = sqrt(sum(abs(x)**2))
    call check(info == QR_SUCCESS .and. &
               abs(euclidean_norm - 1.0_wp) <= eigen_tolerance, &
               'complex norm mode one produces unit Euclidean norm', failures)

    call state%solve(s, v, x, lambda, tolerance, 100, 2, &
                     rel_acc, num_iter, info)
    component_scale = max(maxval(abs(real(x, wp))), maxval(abs(aimag(x))))
    call check(info == QR_SUCCESS .and. &
               abs(component_scale - 1.0_wp) <= eigen_tolerance, &
               'complex other norm modes retain GHEPIIS component scaling', &
               failures)

    write(*,'(a,i0,3(a,es12.4),a,i0)') '  complex solve wp=', wp, &
      ' lambda=', converged_lambda, ' residual=', residual, &
      ' rel_acc=', converged_rel_acc, ' iterations=', converged_num_iter
  end subroutine test_complex_inverse_iteration

  ! Check recoverable solve failures that do not have analytical eigenpairs.
  subroutine test_solve_error_paths(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2
    type(qr_real_state) :: empty_state, singular_state, valid_state
    real(wp) :: h(n,n), s(n,n), v(n), zero_v(n), x(n)
    real(wp) :: lambda, rel_acc
    integer :: info, num_iter

    s = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [n,n])
    h = reshape([2.0_wp, 0.0_wp, 0.0_wp, 4.0_wp], [n,n])
    v = [1.0_wp, 1.0_wp]
    zero_v = 0.0_wp

    call empty_state%solve(s, v, x, lambda, &
                           100.0_wp * epsilon(1.0_wp), 20, 1, &
                           rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 0, &
               'solve rejects an unfactorized state', failures)

    call singular_state%initialize(n, info)
    call singular_state%factorize_fresh(h, s, 2.0_wp, info)
    call check(info == QR_SUCCESS, &
               'exactly shifted singular matrix still has QR storage', failures)
    call singular_state%solve(s, v, x, lambda, &
                              100.0_wp * epsilon(1.0_wp), 20, 1, &
                              rel_acc, num_iter, info)
    call check(info == QR_ERR_SINGULAR .and. num_iter == 0, &
               'solve detects a singular R factor before iteration', failures)

    call valid_state%initialize(n, info)
    call valid_state%factorize_fresh(h, s, 2.25_wp, info)
    call valid_state%solve(s, zero_v, x, lambda, &
                           100.0_wp * epsilon(1.0_wp), 20, 1, &
                           rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 0, &
               'solve rejects a zero starting vector', failures)
    call valid_state%solve(s, v, x, lambda, &
                           100.0_wp * epsilon(1.0_wp), 0, 1, &
                           rel_acc, num_iter, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. num_iter == 0, &
               'solve rejects a nonpositive iteration limit', failures)
  end subroutine test_solve_error_paths

end program test_inverse_iteration

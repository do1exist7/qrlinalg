! Numerical and state-contract tests for the qrlinalg version 0.1 methods.
!
! This source is compiled only by `make test`, against a library built with
! QRLINALG_TESTING. That flag exposes otherwise-private state components so
! the tests can inspect Q and R without adding copy-out accessors to the
! production API. No test-only numerical implementation is substituted.
program test_qrlinalg
  use qrlinalg
  implicit none

  integer :: failures

  failures = 0
  call test_invalid_initialization(failures)
  call test_real_fresh_factorization(failures)
  call test_complex_fresh_factorization(failures)
  call test_real_inverse_iteration(failures)
  call test_complex_inverse_iteration(failures)
  call test_solve_error_paths(failures)

  if (failures /= 0) then
    write(*,'(a,i0,a,i0)') 'FAIL: wp=', wp, ', failed checks=', failures
    error stop 1
  end if
  write(*,'(a,i0)') &
    'PASS: qrlinalg initialization, fresh QR, and inverse iteration, wp=', wp

contains

  ! Record a failed condition without aborting immediately, allowing one test
  ! run to report every violated invariant.
  subroutine check(condition, message, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write(*,'(a)') '  failed: ' // message
    end if
  end subroutine check

  ! Verify rejection of invalid initialization and use of an uninitialized
  ! state. These are recoverable library errors and must never error-stop.
  subroutine test_invalid_initialization(failures)
    integer, intent(inout) :: failures
    type(qr_real_state) :: state
    real(wp) :: h(1,1), s(1,1)
    integer :: info

    h = 1.0_wp
    s = 1.0_wp
    call state%factorize_fresh(h, s, 0.0_wp, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'factorization rejects an uninitialized state', failures)

    call state%initialize(0, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'initialize rejects zero capacity', failures)
    call check(state%capacity == 0 .and. state%n == 0, &
               'failed initialization leaves zero dimensions', failures)
    call check(.not. state%valid, &
               'failed initialization leaves factors invalid', failures)
    call check(.not. allocated(state%q) .and. .not. allocated(state%r), &
               'failed initialization owns no matrix storage', failures)
  end subroutine test_invalid_initialization

  ! Test against a symmetric matrix whose columns are analytically orthogonal:
  !
  !             [ 3   4 ]
  !       M  =  [       ],       ||M(:,j)||_2 = 5.
  !             [ 4  -3 ]
  !
  ! Therefore any valid QR factorization has |diag(R)|=(5,5), zero R(1,2),
  ! and |Q(:,j)|=|M(:,j)|/5, independent of LAPACK's allowed column signs.
  subroutine test_real_fresh_factorization(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2, capacity = 5
    type(qr_real_state) :: state
    real(wp), parameter :: shift = 0.5_wp
    real(wp) :: h(n,n), h_before(n,n), s(n,n), s_before(n,n)
    real(wp) :: m(n,n), identity(n,n), expected_abs_q(n,n)
    real(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: bad_h(2,3), bad_s(2,2)
    real(wp) :: residual, orthogonality, q_error, r_error, tolerance
    integer :: info

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    m = reshape([3.0_wp, 4.0_wp, 4.0_wp, -3.0_wp], [n,n])
    s = reshape([2.0_wp, 0.25_wp, 0.25_wp, 1.5_wp], [n,n])
    h = m + shift * s
    h_before = h
    s_before = s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'real initialize succeeds', failures)
    call check(state%capacity == capacity .and. state%n == 0, &
               'real initialize records capacity and zero active order', failures)
    call check(.not. state%valid, &
               'real initialize does not claim valid factors', failures)
    call check(all(shape(state%q) == [capacity, capacity]) .and. &
               all(shape(state%r) == [capacity, capacity]), &
               'real Q and R use full capacity storage', failures)
    call check(size(state%tau) == capacity, &
               'real tau has capacity elements', failures)
    call check(size(state%factor_work) >= capacity, &
               'real LAPACK workspace meets its minimum size', failures)
    call check(size(state%update_work) == 4 * capacity .and. &
               size(state%solve_work) == 2 * capacity, &
               'real update and solve workspaces have planned sizes', failures)

    ! Fresh factorization must reset only the since-fresh count. The lifetime
    ! structural count represents successful updates across rebuilds.
    state%structural_updates = 7
    state%updates_since_fresh = 3
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'real fresh factorization succeeds', failures)
    call check(state%valid .and. state%n == n, &
               'real factorization commits validity and active order', failures)
    call check(abs(state%shift - shift) <= 0.0_wp, &
               'real factorization records its shift', failures)
    call check(state%structural_updates == 7, &
               'real factorization preserves lifetime update count', failures)
    call check(state%updates_since_fresh == 0, &
               'real factorization resets since-fresh update count', failures)
    call check(all(abs(h - h_before) <= 0.0_wp) .and. &
               all(abs(s - s_before) <= 0.0_wp), &
               'real factorization does not modify caller matrices', failures)

    identity = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [n,n])
    expected_abs_q = reshape([0.6_wp, 0.8_wp, 0.8_wp, 0.6_wp], [n,n])
    residual = real_frobenius(m - matmul(state%q(1:n,1:n), &
                                        state%r(1:n,1:n))) / &
               real_frobenius(m)
    orthogonality = real_frobenius(identity - &
                      matmul(transpose(state%q(1:n,1:n)), &
                             state%q(1:n,1:n))) / sqrt(real(n,wp))
    q_error = maxval(abs(abs(state%q(1:n,1:n)) - expected_abs_q))
    r_error = max(abs(abs(state%r(1,1)) - 5.0_wp), &
                  abs(abs(state%r(2,2)) - 5.0_wp), &
                  abs(state%r(1,2))) / 5.0_wp

    call check(residual <= tolerance, &
               'real Q*R reconstructs the analytical shifted matrix', failures)
    call check(orthogonality <= tolerance, &
               'real Q is orthogonal', failures)
    call check(abs(state%r(2,1)) <= 0.0_wp, &
               'real R is explicitly upper triangular', failures)
    call check(q_error <= tolerance, &
               'real |Q| matches the analytical sign-invariant answer', failures)
    call check(r_error <= tolerance, &
               'real R matches analytical column norms and orthogonality', failures)

    ! Invalid arguments are checked before buffers are touched, so an existing
    ! valid factorization must survive an invalid subsequent request.
    q_before = state%q
    r_before = state%r
    bad_h = 0.0_wp
    bad_s = 0.0_wp
    call state%factorize_fresh(bad_h, bad_s, shift, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real factorization rejects nonsquare H', failures)
    call check(state%valid .and. state%n == n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp), &
               'invalid real request preserves existing factors', failures)

    write(*,'(a,i0,2(a,es12.4))') '  real wp=', wp, &
      ' residual=', residual, ' orthogonality=', orthogonality
  end subroutine test_real_fresh_factorization

  ! Test a Hermitian analytical counterpart. With z=3+4i,
  !
  !             [ 3       z ]
  !       M  =  [           ],    ||M(:,j)||_2 = sqrt(34),
  !             [ conjg(z) -3 ]
  !
  ! and its two columns are exactly orthogonal. Comparing magnitudes removes
  ! the arbitrary complex phase of each QR column.
  subroutine test_complex_fresh_factorization(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 2, capacity = 5
    type(qr_complex_state) :: state
    real(wp), parameter :: shift = -0.25_wp
    real(wp) :: column_norm, residual, unitarity, q_error, r_error
    real(wp) :: tolerance
    real(wp) :: expected_abs_q(n,n)
    complex(wp) :: h(n,n), h_before(n,n), s(n,n), s_before(n,n)
    complex(wp) :: m(n,n), identity(n,n), z
    integer :: info

    tolerance = 1000.0_wp * epsilon(1.0_wp)
    z = cmplx(3.0_wp, 4.0_wp, kind=wp)
    m(1,1) = cmplx(3.0_wp, 0.0_wp, kind=wp)
    m(1,2) = z
    m(2,1) = conjg(z)
    m(2,2) = cmplx(-3.0_wp, 0.0_wp, kind=wp)
    s(1,1) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    s(1,2) = cmplx(0.1_wp, 0.2_wp, kind=wp)
    s(2,1) = conjg(s(1,2))
    s(2,2) = cmplx(1.5_wp, 0.0_wp, kind=wp)
    h = m + cmplx(shift, 0.0_wp, kind=wp) * s
    h_before = h
    s_before = s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'complex initialize succeeds', failures)
    call check(state%capacity == capacity .and. state%n == 0 .and. &
               .not. state%valid, &
               'complex initialize records empty-state metadata', failures)
    call check(all(shape(state%q) == [capacity, capacity]) .and. &
               all(shape(state%r) == [capacity, capacity]), &
               'complex Q and R use full capacity storage', failures)
    call check(size(state%tau) == capacity .and. &
               size(state%factor_work) >= capacity, &
               'complex factorization arrays meet LAPACK requirements', failures)
    call check(size(state%update_work) == 4 * capacity .and. &
               size(state%solve_work) == 2 * capacity .and. &
               size(state%real_work) == capacity, &
               'complex update, solve, and real workspaces have planned sizes', &
               failures)

    state%structural_updates = 9
    state%updates_since_fresh = 4
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, &
               'complex fresh factorization succeeds', failures)
    call check(state%valid .and. state%n == n .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'complex factorization commits metadata', failures)
    call check(state%structural_updates == 9 .and. &
               state%updates_since_fresh == 0, &
               'complex factorization preserves/resets update counters', failures)
    call check(all(abs(h - h_before) <= 0.0_wp) .and. &
               all(abs(s - s_before) <= 0.0_wp), &
               'complex factorization does not modify caller matrices', failures)

    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    identity(1,1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    identity(2,2) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    column_norm = sqrt(34.0_wp)
    expected_abs_q = reshape([3.0_wp, 5.0_wp, 5.0_wp, 3.0_wp], [n,n]) / &
                     column_norm
    residual = complex_frobenius(m - matmul(state%q(1:n,1:n), &
                                           state%r(1:n,1:n))) / &
               complex_frobenius(m)
    unitarity = complex_frobenius(identity - &
                  matmul(conjg(transpose(state%q(1:n,1:n))), &
                         state%q(1:n,1:n))) / sqrt(real(n,wp))
    q_error = maxval(abs(abs(state%q(1:n,1:n)) - expected_abs_q))
    r_error = max(abs(abs(state%r(1,1)) - column_norm), &
                  abs(abs(state%r(2,2)) - column_norm), &
                  abs(state%r(1,2))) / column_norm

    call check(residual <= tolerance, &
               'complex Q*R reconstructs the analytical shifted matrix', failures)
    call check(unitarity <= tolerance, &
               'complex Q is unitary', failures)
    call check(abs(state%r(2,1)) <= 0.0_wp, &
               'complex R is explicitly upper triangular', failures)
    call check(q_error <= tolerance, &
               'complex |Q| matches the analytical phase-invariant answer', &
               failures)
    call check(r_error <= tolerance, &
               'complex R matches analytical column norms and orthogonality', &
               failures)

    write(*,'(a,i0,2(a,es12.4))') '  complex wp=', wp, &
      ' residual=', residual, ' unitarity=', unitarity
  end subroutine test_complex_fresh_factorization

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

    call empty_state%solve(s, v, x, lambda, 100.0_wp * epsilon(1.0_wp), &
                           20, 1, rel_acc, num_iter, info)
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

  ! Frobenius norm of a real matrix, implemented with kind-correct intrinsics
  ! so the reference check works unchanged at wp=8, 10, and 16.
  function real_frobenius(a) result(norm)
    real(wp), intent(in) :: a(:,:)
    real(wp) :: norm

    norm = sqrt(sum(a * a))
  end function real_frobenius

  ! Frobenius norm of a complex matrix. ABS returns real(wp), so no fixed
  ! double BLAS or LAPACK is introduced into the analytical reference path.
  function complex_frobenius(a) result(norm)
    complex(wp), intent(in) :: a(:,:)
    real(wp) :: norm

    norm = sqrt(sum(abs(a)**2))
  end function complex_frobenius

end program test_qrlinalg

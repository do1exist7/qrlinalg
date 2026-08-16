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

  if (failures /= 0) then
    write(*,'(a,i0,a,i0)') 'FAIL: wp=', wp, ', failed checks=', failures
    error stop 1
  end if
  write(*,'(a,i0)') 'PASS: qrlinalg initialization and fresh QR, wp=', wp

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

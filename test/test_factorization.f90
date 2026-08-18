! Analytical fresh-factorization tests for real and complex qrlinalg states.
!
! The test matrices have orthogonal columns, which provides known QR factor
! magnitudes while allowing the signs or phases selected by LAPACK to vary.
program test_factorization
  use qrlinalg
  use test_support, only: check, finish_suite, real_frobenius, &
                          complex_frobenius
  implicit none

  integer :: failures

  failures = 0
  call test_real_fresh_factorization(failures)
  call test_complex_fresh_factorization(failures)
  call finish_suite('qrlinalg fresh factorization', failures)

contains

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

    !The lower triangle defines each symmetric input. Deliberately inconsistent
    !upper entries verify that fresh factorization never references them.
    h(1,2) = 1234.0_wp
    s(1,2) = -5678.0_wp
    h_before = h
    s_before = s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, &
               'real factorization state initializes', failures)

    ! Fresh factorization resets only the since-fresh count. The lifetime
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
    call check(orthogonality <= tolerance, 'real Q is orthogonal', failures)
    call check(abs(state%r(2,1)) <= 0.0_wp, &
               'real R is explicitly upper triangular', failures)
    call check(q_error <= tolerance, &
               'real |Q| matches the analytical sign-invariant answer', failures)
    call check(r_error <= tolerance, &
               'real R matches analytical column norms and orthogonality', failures)

    ! Argument validation precedes all writes to factor storage. An existing
    ! valid factorization must therefore survive a rejected subsequent call.
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
    complex(wp) :: invalid_diagonal_h(n,n), invalid_diagonal_s(n,n)
    complex(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    complex(wp) :: bad_h(2,3), bad_s(2,2)
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

    !Only the lower Hermitian triangle is part of the input contract. Poison
    !the upper entries with values unrelated to their lower conjugates.
    h(1,2) = cmplx(1234.0_wp, -4321.0_wp, kind=wp)
    s(1,2) = cmplx(-5678.0_wp, 8765.0_wp, kind=wp)
    h_before = h
    s_before = s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, &
               'complex factorization state initializes', failures)

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
    call check(unitarity <= tolerance, 'complex Q is unitary', failures)
    call check(abs(state%r(2,1)) <= 0.0_wp, &
               'complex R is explicitly upper triangular', failures)
    call check(q_error <= tolerance, &
               'complex |Q| matches the analytical phase-invariant answer', &
               failures)
    call check(r_error <= tolerance, &
               'complex R matches analytical column norms and orthogonality', &
               failures)

    q_before = state%q
    r_before = state%r
    bad_h = cmplx(0.0_wp, 0.0_wp, kind=wp)
    bad_s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call state%factorize_fresh(bad_h, bad_s, shift, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex factorization rejects nonsquare H', failures)
    call check(state%valid .and. state%n == n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp), &
               'invalid complex request preserves existing factors', failures)

    !A Hermitian diagonal is real. Reject a physical matrix that violates this
    !condition before overwriting the previously valid factorization.
    invalid_diagonal_h = h
    invalid_diagonal_h(1,1) = cmplx(real(h(1,1), wp), 1.0_wp, kind=wp)
    call state%factorize_fresh(invalid_diagonal_h, s, shift, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex factorization rejects a non-real H diagonal', failures)
    call check(state%valid .and. state%n == n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp), &
               'non-real complex diagonal preserves existing factors', failures)

    invalid_diagonal_s = s
    invalid_diagonal_s(2,2) = cmplx(real(s(2,2), wp), -1.0_wp, kind=wp)
    call state%factorize_fresh(h, invalid_diagonal_s, shift, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex factorization rejects a non-real S diagonal', failures)
    call check(state%valid .and. state%n == n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp), &
               'non-real overlap diagonal preserves existing factors', failures)

    write(*,'(a,i0,2(a,es12.4))') '  complex wp=', wp, &
      ' residual=', residual, ' unitarity=', unitarity
  end subroutine test_complex_fresh_factorization

end program test_factorization

! Analytical tests for symmetric and Hermitian row-and-column replacement.
!
! Each test forms the updated shifted matrix independently from the defining
! rank-two expression. The stored factors are judged by reconstruction and
! orthogonality/unitarity, so the tests do not depend on LAPACK's arbitrary QR
! column signs or phases.
program test_replacement
  use qrlinalg
  use test_support, only: check, finish_suite, real_frobenius, &
                          complex_frobenius
  implicit none

  integer :: failures

  failures = 0
  call test_real_replacement(failures)
  call test_complex_replacement(failures)
  call finish_suite('qrlinalg symmetric replacement', failures)

contains

  ! Verify the real update
  !
  !       M_new = M + d*e_i^T + e_i*(d-d(i)*e_i)^T,
  !       d     = delta_h - shift*delta_s,
  !
  ! including preservation of caller vectors and state on rejected requests.
  subroutine test_real_replacement(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 3, capacity = 5, idx = 2
    type(qr_real_state) :: state, empty_state
    real(wp), parameter :: shift = 0.25_wp
    real(wp) :: h(n,n), s(n,n), shifted(n,n), expected(n,n)
    real(wp) :: delta_h(n), delta_s(n), delta_h_before(n), delta_s_before(n)
    real(wp) :: identity(n,n), q_before(capacity,capacity)
    real(wp) :: r_before(capacity,capacity), short_vector(n-1)
    real(wp) :: change(n), residual, orthogonality, tolerance
    integer :: info, j

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    shifted = reshape([4.0_wp, 1.0_wp, -0.5_wp, &
                       1.0_wp, 3.0_wp,  0.75_wp, &
                      -0.5_wp, 0.75_wp, 5.0_wp], [n,n])
    s = reshape([2.0_wp, 0.2_wp,  0.1_wp, &
                 0.2_wp, 1.5_wp, -0.1_wp, &
                 0.1_wp, -0.1_wp, 1.2_wp], [n,n])
    h = shifted + shift * s
    delta_h = [0.3_wp, -0.4_wp, 0.7_wp]
    delta_s = [0.1_wp, 0.2_wp, -0.05_wp]
    delta_h_before = delta_h
    delta_s_before = delta_s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'real replacement state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'real replacement base factorizes', failures)
    state%structural_updates = 4
    state%updates_since_fresh = 2

    change = delta_h - shift * delta_s
    expected = shifted
    expected(:,idx) = expected(:,idx) + change
    do j = 1, n
      if (j /= idx) expected(idx,j) = expected(idx,j) + change(j)
    end do

    call state%replace_symmetric(idx, delta_h, delta_s, info)
    call check(info == QR_SUCCESS, 'real symmetric replacement succeeds', failures)
    call check(state%valid .and. state%n == n .and. &
               state%capacity == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'real replacement preserves factorization metadata', failures)
    call check(state%structural_updates == 5 .and. &
               state%updates_since_fresh == 3, &
               'real replacement increments each counter exactly once', failures)
    call check(all(abs(delta_h - delta_h_before) <= 0.0_wp) .and. &
               all(abs(delta_s - delta_s_before) <= 0.0_wp), &
               'real replacement preserves caller vectors', failures)

    identity = 0.0_wp
    do j = 1, n
      identity(j,j) = 1.0_wp
    end do
    residual = real_frobenius(expected - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(real_frobenius(expected), tiny(1.0_wp))
    orthogonality = real_frobenius(identity - &
      matmul(transpose(state%q(1:n,1:n)), state%q(1:n,1:n))) / &
      sqrt(real(n,wp))
    call check(residual <= tolerance, &
               'real updated Q*R reconstructs the analytical matrix', failures)
    call check(orthogonality <= tolerance, &
               'real updated Q remains orthogonal', failures)
    call check(abs(state%r(2,1)) <= tolerance .and. &
               abs(state%r(3,1)) <= tolerance .and. &
               abs(state%r(3,2)) <= tolerance, &
               'real updated R remains upper triangular', failures)

    ! Every validation failure must precede the destructive qr1up calls.
    q_before = state%q
    r_before = state%r
    call state%replace_symmetric(0, delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real replacement rejects an invalid index', failures)
    short_vector = 0.0_wp
    call state%replace_symmetric(idx, short_vector, short_vector, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real replacement rejects incorrect vector extents', failures)
    call check(all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 5 .and. &
               state%updates_since_fresh == 3, &
               'rejected real replacements preserve the complete state', failures)
    call empty_state%replace_symmetric(1, delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real replacement rejects an uninitialized state', failures)

    write(*,'(a,i0,2(a,es12.4))') '  real replacement wp=', wp, &
      ' residual=', residual, ' orthogonality=', orthogonality
  end subroutine test_real_replacement

  ! Verify the Hermitian update. The independently formed reference row uses
  ! conjugated off-diagonal changes, while its diagonal is changed only once.
  ! A non-real diagonal change is also checked as a pre-update error.
  subroutine test_complex_replacement(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 3, capacity = 5, idx = 2
    type(qr_complex_state) :: state
    real(wp), parameter :: shift = -0.4_wp
    complex(wp) :: h(n,n), s(n,n), shifted(n,n), expected(n,n)
    complex(wp) :: delta_h(n), delta_s(n), delta_h_before(n), delta_s_before(n)
    complex(wp) :: identity(n,n), q_before(capacity,capacity)
    complex(wp) :: r_before(capacity,capacity), bad_delta_h(n), change(n)
    real(wp) :: residual, unitarity, hermitian_error, tolerance
    integer :: info, j

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    shifted(1,1) = cmplx(4.0_wp, 0.0_wp, kind=wp)
    shifted(1,2) = cmplx(0.5_wp, 0.3_wp, kind=wp)
    shifted(1,3) = cmplx(-0.2_wp, 0.4_wp, kind=wp)
    shifted(2,1) = conjg(shifted(1,2))
    shifted(2,2) = cmplx(3.5_wp, 0.0_wp, kind=wp)
    shifted(2,3) = cmplx(0.6_wp, -0.25_wp, kind=wp)
    shifted(3,1) = conjg(shifted(1,3))
    shifted(3,2) = conjg(shifted(2,3))
    shifted(3,3) = cmplx(5.0_wp, 0.0_wp, kind=wp)
    s(1,1) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    s(1,2) = cmplx(0.1_wp, -0.05_wp, kind=wp)
    s(1,3) = cmplx(0.0_wp, 0.08_wp, kind=wp)
    s(2,1) = conjg(s(1,2))
    s(2,2) = cmplx(1.6_wp, 0.0_wp, kind=wp)
    s(2,3) = cmplx(-0.12_wp, 0.03_wp, kind=wp)
    s(3,1) = conjg(s(1,3))
    s(3,2) = conjg(s(2,3))
    s(3,3) = cmplx(1.3_wp, 0.0_wp, kind=wp)
    h = shifted + cmplx(shift, 0.0_wp, kind=wp) * s
    delta_h = [cmplx(0.4_wp, -0.2_wp, kind=wp), &
               cmplx(0.75_wp, 0.0_wp, kind=wp), &
               cmplx(-0.3_wp, 0.5_wp, kind=wp)]
    delta_s = [cmplx(0.1_wp, 0.05_wp, kind=wp), &
               cmplx(-0.2_wp, 0.0_wp, kind=wp), &
               cmplx(0.2_wp, -0.1_wp, kind=wp)]
    delta_h_before = delta_h
    delta_s_before = delta_s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, &
               'complex replacement state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'complex replacement base factorizes', failures)
    state%structural_updates = 6
    state%updates_since_fresh = 1

    change = delta_h - cmplx(shift, 0.0_wp, kind=wp) * delta_s
    expected = shifted
    expected(:,idx) = expected(:,idx) + change
    do j = 1, n
      if (j /= idx) expected(idx,j) = expected(idx,j) + conjg(change(j))
    end do

    call state%replace_symmetric(idx, delta_h, delta_s, info)
    call check(info == QR_SUCCESS, &
               'complex Hermitian replacement succeeds', failures)
    call check(state%valid .and. state%n == n .and. &
               state%capacity == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'complex replacement preserves factorization metadata', failures)
    call check(state%structural_updates == 7 .and. &
               state%updates_since_fresh == 2, &
               'complex replacement increments each counter exactly once', failures)
    call check(all(abs(delta_h - delta_h_before) <= 0.0_wp) .and. &
               all(abs(delta_s - delta_s_before) <= 0.0_wp), &
               'complex replacement preserves caller vectors', failures)

    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do j = 1, n
      identity(j,j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do
    residual = complex_frobenius(expected - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(complex_frobenius(expected), tiny(1.0_wp))
    unitarity = complex_frobenius(identity - &
      matmul(conjg(transpose(state%q(1:n,1:n))), state%q(1:n,1:n))) / &
      sqrt(real(n,wp))
    hermitian_error = maxval(abs(expected - conjg(transpose(expected))))
    call check(residual <= tolerance, &
               'complex updated Q*R reconstructs the analytical matrix', failures)
    call check(unitarity <= tolerance, &
               'complex updated Q remains unitary', failures)
    call check(hermitian_error <= tolerance, &
               'complex analytical replacement remains Hermitian', failures)
    call check(abs(state%r(2,1)) <= tolerance .and. &
               abs(state%r(3,1)) <= tolerance .and. &
               abs(state%r(3,2)) <= tolerance, &
               'complex updated R remains upper triangular', failures)

    q_before = state%q
    r_before = state%r
    bad_delta_h = delta_h
    bad_delta_h(idx) = bad_delta_h(idx) + &
      cmplx(0.0_wp, 0.01_wp, kind=wp)
    call state%replace_symmetric(idx, bad_delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex replacement rejects a non-real diagonal change', failures)
    call check(all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 7 .and. &
               state%updates_since_fresh == 2, &
               'rejected complex replacement preserves the complete state', &
               failures)

    write(*,'(a,i0,2(a,es12.4))') '  complex replacement wp=', wp, &
      ' residual=', residual, ' unitarity=', unitarity
  end subroutine test_complex_replacement

end program test_replacement

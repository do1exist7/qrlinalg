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
    integer, parameter :: n = 3, capacity = 5, num_updates = 100
    type(qr_real_state) :: state, empty_state
    real(wp), parameter :: shift = 0.25_wp
    real(wp) :: h(n,n), s(n,n), shifted(n,n), expected(n,n)
    real(wp) :: delta_h(n), delta_s(n), delta_h_before(n), delta_s_before(n)
    real(wp) :: identity(n,n), q_before(capacity,capacity)
    real(wp) :: r_before(capacity,capacity), short_vector(n-1)
    real(wp) :: change(n), residual, orthogonality, lower_triangle_error
    real(wp) :: max_residual, max_orthogonality, max_lower_triangle_error
    real(wp) :: tolerance
    integer :: info, j, update_idx, update_number
    logical :: all_updates_succeeded, caller_vectors_preserved
    logical :: all_reconstructions_accurate, all_factors_orthogonal
    logical :: all_factors_triangular

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    shifted = reshape([4.0_wp, 1.0_wp, -0.5_wp, &
                       1.0_wp, 3.0_wp,  0.75_wp, &
                      -0.5_wp, 0.75_wp, 5.0_wp], [n,n])
    s = reshape([2.0_wp, 0.2_wp,  0.1_wp, &
                 0.2_wp, 1.5_wp, -0.1_wp, &
                 0.1_wp, -0.1_wp, 1.2_wp], [n,n])
    h = shifted + shift * s
    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'real replacement state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'real replacement base factorizes', failures)
    state%structural_updates = 4
    state%updates_since_fresh = 2

    expected = shifted
    identity = 0.0_wp
    do j = 1, n
      identity(j,j) = 1.0_wp
    end do
    all_updates_succeeded = .true.
    caller_vectors_preserved = .true.
    max_residual = 0.0_wp
    max_orthogonality = 0.0_wp
    max_lower_triangle_error = 0.0_wp
    all_reconstructions_accurate = .true.
    all_factors_orthogonal = .true.
    all_factors_triangular = .true.

    ! Cycle through every index with deterministic changes of both H and S.
    ! The reference matrix is accumulated directly from the mathematical
    ! symmetric update, independently of the two qr1up calls under test.
    do update_number = 1, num_updates
      update_idx = mod(update_number - 1, n) + 1
      do j = 1, n
        delta_h(j) = 0.002_wp * &
          real(mod(update_number * (j + 1), 11) - 5, wp)
        delta_s(j) = 0.001_wp * &
          real(mod((update_number + 2) * (2 * j + 1), 13) - 6, wp)
      end do
      delta_h_before = delta_h
      delta_s_before = delta_s
      change = delta_h - shift * delta_s
      expected(:,update_idx) = expected(:,update_idx) + change
      do j = 1, n
        if (j /= update_idx) then
          expected(update_idx,j) = expected(update_idx,j) + change(j)
        end if
      end do

      call state%replace_symmetric(update_idx, delta_h, delta_s, info)
      all_updates_succeeded = all_updates_succeeded .and. info == QR_SUCCESS
      caller_vectors_preserved = caller_vectors_preserved .and. &
        all(abs(delta_h - delta_h_before) <= 0.0_wp) .and. &
        all(abs(delta_s - delta_s_before) <= 0.0_wp)

      ! Inspect every intermediate factorization. A later replacement must not
      ! be able to conceal an earlier reconstruction or orthogonality failure.
      residual = real_frobenius(expected - &
        matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
        max(real_frobenius(expected), tiny(1.0_wp))
      orthogonality = real_frobenius(identity - &
        matmul(transpose(state%q(1:n,1:n)), state%q(1:n,1:n))) / &
        sqrt(real(n,wp))
      lower_triangle_error = max(abs(state%r(2,1)), &
        abs(state%r(3,1)), abs(state%r(3,2)))
      max_residual = max(max_residual, residual)
      max_orthogonality = max(max_orthogonality, orthogonality)
      max_lower_triangle_error = max(max_lower_triangle_error, &
                                     lower_triangle_error)
      all_reconstructions_accurate = all_reconstructions_accurate .and. &
                                     residual <= tolerance
      all_factors_orthogonal = all_factors_orthogonal .and. &
                               orthogonality <= tolerance
      all_factors_triangular = all_factors_triangular .and. &
                               lower_triangle_error <= tolerance
    end do

    call check(all_updates_succeeded, &
               'all 100 real symmetric replacements succeed', failures)
    call check(state%valid .and. state%n == n .and. &
               state%capacity == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'real replacement preserves factorization metadata', failures)
    call check(state%structural_updates == 4 + num_updates .and. &
               state%updates_since_fresh == 2 + num_updates, &
               '100 real replacements increment each counter 100 times', failures)
    call check(caller_vectors_preserved, &
               'all real replacements preserve caller vectors', failures)

    call check(all_reconstructions_accurate, &
               'every real updated Q*R reconstructs its analytical matrix', &
               failures)
    call check(all_factors_orthogonal, &
               'real updated Q remains orthogonal after every update', failures)
    call check(all_factors_triangular, &
               'real updated R remains upper triangular after every update', &
               failures)

    ! Every validation failure must precede the destructive qr1up calls.
    q_before = state%q
    r_before = state%r
    call state%replace_symmetric(0, delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real replacement rejects an invalid index', failures)
    short_vector = 0.0_wp
    call state%replace_symmetric(2, short_vector, short_vector, info)
    call check(info == QR_ERR_DIMENSION_MISMATCH, &
               'real replacement rejects incorrect vector extents', failures)
    call check(all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 4 + num_updates .and. &
               state%updates_since_fresh == 2 + num_updates, &
               'rejected real replacements preserve the complete state', failures)
    call empty_state%replace_symmetric(1, delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_STATE, &
               'real replacement rejects an uninitialized state', failures)

    write(*,'(a,i0,2(a,es12.4))') '  real replacement wp=', wp, &
      ' max residual=', max_residual, &
      ' max orthogonality=', max_orthogonality
  end subroutine test_real_replacement

  ! Verify the Hermitian update. The independently formed reference row uses
  ! conjugated off-diagonal changes, while its diagonal is changed only once.
  ! A non-real diagonal change is also checked as a pre-update error.
  subroutine test_complex_replacement(failures)
    integer, intent(inout) :: failures
    integer, parameter :: n = 3, capacity = 5, num_updates = 100
    type(qr_complex_state) :: state, empty_state
    real(wp), parameter :: shift = -0.4_wp
    complex(wp) :: h(n,n), s(n,n), shifted(n,n), expected(n,n)
    complex(wp) :: delta_h(n), delta_s(n), delta_h_before(n), delta_s_before(n)
    complex(wp) :: identity(n,n), q_before(capacity,capacity)
    complex(wp) :: r_before(capacity,capacity), bad_delta_h(n), change(n)
    complex(wp) :: short_vector(n-1)
    real(wp) :: residual, unitarity, lower_triangle_error, hermitian_error
    real(wp) :: max_residual, max_unitarity, max_lower_triangle_error
    real(wp) :: tolerance
    integer :: info, j, update_idx, update_number
    logical :: all_updates_succeeded, caller_vectors_preserved
    logical :: all_reconstructions_accurate, all_factors_unitary
    logical :: all_factors_triangular

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
    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, &
               'complex replacement state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'complex replacement base factorizes', failures)
    state%structural_updates = 6
    state%updates_since_fresh = 1

    expected = shifted
    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do j = 1, n
      identity(j,j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do
    all_updates_succeeded = .true.
    caller_vectors_preserved = .true.
    max_residual = 0.0_wp
    max_unitarity = 0.0_wp
    max_lower_triangle_error = 0.0_wp
    all_reconstructions_accurate = .true.
    all_factors_unitary = .true.
    all_factors_triangular = .true.

    ! Exercise every Hermitian row and column repeatedly. The diagonal entries
    ! are made exactly real; off-diagonal entries carry changing real and
    ! imaginary parts so both complex qrupdate paths remain active.
    do update_number = 1, num_updates
      update_idx = mod(update_number - 1, n) + 1
      do j = 1, n
        delta_h(j) = cmplx( &
          0.0015_wp * real(mod(update_number * (j + 2), 13) - 6, wp), &
          0.001_wp * real(mod((update_number + 1) * (2 * j + 1), 11) - 5, wp), &
          kind=wp)
        delta_s(j) = cmplx( &
          0.0008_wp * real(mod((update_number + 3) * (j + 1), 9) - 4, wp), &
          0.0006_wp * real(mod((update_number + 2) * (j + 3), 7) - 3, wp), &
          kind=wp)
      end do
      delta_h(update_idx) = &
        cmplx(real(delta_h(update_idx), wp), 0.0_wp, kind=wp)
      delta_s(update_idx) = &
        cmplx(real(delta_s(update_idx), wp), 0.0_wp, kind=wp)
      delta_h_before = delta_h
      delta_s_before = delta_s
      change = delta_h - cmplx(shift, 0.0_wp, kind=wp) * delta_s
      expected(:,update_idx) = expected(:,update_idx) + change
      do j = 1, n
        if (j /= update_idx) then
          expected(update_idx,j) = expected(update_idx,j) + conjg(change(j))
        end if
      end do

      call state%replace_symmetric(update_idx, delta_h, delta_s, info)
      all_updates_succeeded = all_updates_succeeded .and. info == QR_SUCCESS
      caller_vectors_preserved = caller_vectors_preserved .and. &
        all(abs(delta_h - delta_h_before) <= 0.0_wp) .and. &
        all(abs(delta_s - delta_s_before) <= 0.0_wp)

      residual = complex_frobenius(expected - &
        matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
        max(complex_frobenius(expected), tiny(1.0_wp))
      unitarity = complex_frobenius(identity - &
        matmul(conjg(transpose(state%q(1:n,1:n))), &
               state%q(1:n,1:n))) / sqrt(real(n,wp))
      lower_triangle_error = max(abs(state%r(2,1)), &
        abs(state%r(3,1)), abs(state%r(3,2)))
      max_residual = max(max_residual, residual)
      max_unitarity = max(max_unitarity, unitarity)
      max_lower_triangle_error = max(max_lower_triangle_error, &
                                     lower_triangle_error)
      all_reconstructions_accurate = all_reconstructions_accurate .and. &
                                     residual <= tolerance
      all_factors_unitary = all_factors_unitary .and. unitarity <= tolerance
      all_factors_triangular = all_factors_triangular .and. &
                               lower_triangle_error <= tolerance
    end do

    call check(all_updates_succeeded, &
               'all 100 complex Hermitian replacements succeed', failures)
    call check(state%valid .and. state%n == n .and. &
               state%capacity == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'complex replacement preserves factorization metadata', failures)
    call check(state%structural_updates == 6 + num_updates .and. &
               state%updates_since_fresh == 1 + num_updates, &
               '100 complex replacements increment each counter 100 times', &
               failures)
    call check(caller_vectors_preserved, &
               'all complex replacements preserve caller vectors', failures)

    hermitian_error = maxval(abs(expected - conjg(transpose(expected))))
    call check(all_reconstructions_accurate, &
               'every complex updated Q*R reconstructs its analytical matrix', &
               failures)
    call check(all_factors_unitary, &
               'complex updated Q remains unitary after every update', failures)
    call check(hermitian_error <= tolerance, &
               'complex analytical replacement remains Hermitian', failures)
    call check(all_factors_triangular, &
               'complex updated R remains triangular after every update', &
               failures)

    q_before = state%q
    r_before = state%r
    bad_delta_h = delta_h
    bad_delta_h(update_idx) = bad_delta_h(update_idx) + &
      cmplx(0.0_wp, 0.01_wp, kind=wp)
    call state%replace_symmetric(update_idx, bad_delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex replacement rejects a non-real diagonal change', failures)
    call check(all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 6 + num_updates .and. &
               state%updates_since_fresh == 1 + num_updates, &
               'rejected complex replacement preserves the complete state', &
               failures)

    short_vector = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call state%replace_symmetric(update_idx, short_vector, short_vector, info)
    call check(info == QR_ERR_DIMENSION_MISMATCH, &
               'complex replacement distinguishes incorrect extents', failures)
    call empty_state%replace_symmetric(update_idx, delta_h, delta_s, info)
    call check(info == QR_ERR_INVALID_STATE, &
               'complex replacement distinguishes an invalid state', failures)

    write(*,'(a,i0,2(a,es12.4))') '  complex replacement wp=', wp, &
      ' max residual=', max_residual, ' max unitarity=', max_unitarity
  end subroutine test_complex_replacement

end program test_replacement

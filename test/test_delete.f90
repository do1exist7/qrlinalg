! Analytical tests for symmetric and Hermitian principal-submatrix deletion.
!
! Real and complex states are first expanded and reduced in an append/delete
! round trip. They then delete middle, first, and last indices successively
! until order one. The reference matrices are reduced directly by removing the
! requested rows and columns, independently of qrdec and qrder.
program test_delete
  use qrlinalg
  use test_support, only: check, finish_suite, real_frobenius, &
                          complex_frobenius
  implicit none

  integer :: failures

  failures = 0
  call test_real_delete(failures)
  call test_complex_delete(failures)
  call finish_suite('qrlinalg symmetric deletion', failures)

contains

  ! Exercise deletion for a real symmetric shifted matrix, including
  ! composition with append and state preservation on every rejected request.
  subroutine test_real_delete(failures)
    integer, intent(inout) :: failures
    integer, parameter :: initial_n = 5, capacity = 6, num_deletes = 4
    integer, parameter :: delete_indices(num_deletes) = [3, 1, 3, 2]
    type(qr_real_state) :: state, empty_state
    real(wp), parameter :: shift = 0.2_wp
    real(wp) :: h(initial_n,initial_n), s(initial_n,initial_n)
    real(wp) :: expected(capacity,capacity), shifted_column(capacity)
    real(wp) :: h_column(capacity), s_column(capacity)
    real(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: residual, orthogonality, lower_error, tolerance
    real(wp) :: max_residual, max_orthogonality
    integer :: delete_number, i, idx, info, j, new_n, old_n
    integer :: successful_updates
    logical :: factors_accurate, inactive_storage_cleared

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    expected = 0.0_wp
    do j = 1, initial_n
      expected(j,j) = 4.0_wp + real(j,wp)
      do i = j + 1, initial_n
        expected(i,j) = 0.03_wp * real(i + 2 * j,wp)
        expected(j,i) = expected(i,j)
      end do
    end do
    do j = 1, initial_n
      do i = 1, initial_n
        if (i == j) then
          s(i,j) = 1.0_wp + 0.1_wp * real(i,wp)
        else
          s(i,j) = 0.005_wp * real(i + j,wp)
        end if
      end do
    end do
    h = expected(1:initial_n,1:initial_n) + shift * s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'real deletion state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'real deletion base factorizes', failures)
    state%structural_updates = 4
    state%updates_since_fresh = 2
    successful_updates = 0

    ! Append a sixth row and column and immediately delete it. This checks that
    ! the two structural operations compose without requiring fresh factors.
    do i = 1, initial_n
      shifted_column(i) = 0.025_wp * real(2 * i - 3,wp)
      s_column(i) = 0.004_wp * real(i - 2,wp)
    end do
    shifted_column(capacity) = 8.0_wp
    s_column(capacity) = 1.4_wp
    h_column = shifted_column + shift * s_column
    expected(1:capacity,capacity) = shifted_column
    expected(capacity,1:initial_n) = shifted_column(1:initial_n)
    call state%append_symmetric(h_column, s_column, info)
    call check(info == QR_SUCCESS .and. state%n == capacity, &
               'real append/delete round trip appends successfully', failures)
    successful_updates = successful_updates + 1
    old_n = state%n
    call state%delete_symmetric(capacity, info)
    call remove_real_index(expected, old_n, capacity)
    successful_updates = successful_updates + 1
    call real_factor_errors(state, expected, initial_n, residual, &
                            orthogonality, lower_error)
    call check(info == QR_SUCCESS .and. state%n == initial_n .and. &
               residual <= tolerance .and. orthogonality <= tolerance .and. &
               lower_error <= tolerance, &
               'real append/delete round trip restores the original matrix', &
               failures)

    ! Invalid indices must be rejected before qrdec uses R's final column as
    ! scratch storage.
    q_before = state%q
    r_before = state%r
    call state%delete_symmetric(0, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real deletion rejects index zero', failures)
    call state%delete_symmetric(state%n + 1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 4 + successful_updates .and. &
               state%updates_since_fresh == 2 + successful_updates, &
               'invalid real deletion preserves the complete state', failures)

    factors_accurate = .true.
    inactive_storage_cleared = .true.
    max_residual = residual
    max_orthogonality = orthogonality
    do delete_number = 1, num_deletes
      old_n = state%n
      idx = delete_indices(delete_number)
      new_n = old_n - 1
      call state%delete_symmetric(idx, info)
      call remove_real_index(expected, old_n, idx)
      successful_updates = successful_updates + 1
      call real_factor_errors(state, expected, new_n, residual, &
                              orthogonality, lower_error)
      factors_accurate = factors_accurate .and. info == QR_SUCCESS .and. &
        state%n == new_n .and. state%valid .and. &
        residual <= tolerance .and. orthogonality <= tolerance .and. &
        lower_error <= tolerance .and. &
        state%structural_updates == 4 + successful_updates .and. &
        state%updates_since_fresh == 2 + successful_updates
      inactive_storage_cleared = inactive_storage_cleared .and. &
        all(abs(state%q(old_n,:)) <= 0.0_wp) .and. &
        all(abs(state%q(:,old_n)) <= 0.0_wp) .and. &
        all(abs(state%r(old_n,:)) <= 0.0_wp) .and. &
        all(abs(state%r(:,old_n)) <= 0.0_wp)
      max_residual = max(max_residual, residual)
      max_orthogonality = max(max_orthogonality, orthogonality)
    end do
    call check(factors_accurate, &
               'every real deletion preserves accurate factors and counters', &
               failures)
    call check(inactive_storage_cleared, &
               'real deletion clears every newly inactive factor row/column', &
               failures)
    call check(state%n == 1 .and. state%valid .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'real deletion retains a valid order-one factorization', failures)

    q_before = state%q
    r_before = state%r
    call state%delete_symmetric(1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%n == 1 .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 4 + successful_updates .and. &
               state%updates_since_fresh == 2 + successful_updates, &
               'real order-one deletion is rejected without mutation', failures)
    call empty_state%delete_symmetric(1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real deletion rejects an unfactorized state', failures)

    write(*,'(a,i0,2(a,es12.4))') '  real deletion wp=', wp, &
      ' max residual=', max_residual, &
      ' max orthogonality=', max_orthogonality
  end subroutine test_real_delete

  ! Apply the same index sequence to a dense Hermitian shifted matrix. Complex
  ! off-diagonal phases ensure row deletion and unitary rotations are tested,
  ! rather than merely exercising the real subset of the complex kernels.
  subroutine test_complex_delete(failures)
    integer, intent(inout) :: failures
    integer, parameter :: initial_n = 5, capacity = 6, num_deletes = 4
    integer, parameter :: delete_indices(num_deletes) = [3, 1, 3, 2]
    type(qr_complex_state) :: state, empty_state
    real(wp), parameter :: shift = -0.15_wp
    complex(wp) :: h(initial_n,initial_n), s(initial_n,initial_n)
    complex(wp) :: expected(capacity,capacity), shifted_column(capacity)
    complex(wp) :: h_column(capacity), s_column(capacity)
    complex(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: residual, unitarity, lower_error, tolerance
    real(wp) :: max_residual, max_unitarity
    integer :: delete_number, i, idx, info, j, new_n, old_n
    integer :: successful_updates
    logical :: factors_accurate, inactive_storage_cleared

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    expected = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do j = 1, initial_n
      expected(j,j) = cmplx(4.0_wp + real(j,wp), 0.0_wp, kind=wp)
      s(j,j) = cmplx(1.0_wp + 0.1_wp * real(j,wp), 0.0_wp, kind=wp)
      do i = j + 1, initial_n
        expected(i,j) = cmplx(0.025_wp * real(i + j,wp), &
                              0.015_wp * real(i - 2 * j,wp), kind=wp)
        expected(j,i) = conjg(expected(i,j))
        s(i,j) = cmplx(0.004_wp * real(i + j,wp), &
                       0.003_wp * real(i - j,wp), kind=wp)
        s(j,i) = conjg(s(i,j))
      end do
    end do
    h = expected(1:initial_n,1:initial_n) + &
        cmplx(shift, 0.0_wp, kind=wp) * s

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, &
               'complex deletion state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'complex deletion base factorizes', failures)
    state%structural_updates = 6
    state%updates_since_fresh = 3
    successful_updates = 0

    do i = 1, initial_n
      shifted_column(i) = cmplx(0.02_wp * real(i - 2,wp), &
                                0.018_wp * real(3 - i,wp), kind=wp)
      s_column(i) = cmplx(0.003_wp * real(i - 1,wp), &
                          0.002_wp * real(2 * i - 5,wp), kind=wp)
    end do
    shifted_column(capacity) = cmplx(8.0_wp, 0.0_wp, kind=wp)
    s_column(capacity) = cmplx(1.4_wp, 0.0_wp, kind=wp)
    h_column = shifted_column + &
               cmplx(shift, 0.0_wp, kind=wp) * s_column
    expected(1:capacity,capacity) = shifted_column
    expected(capacity,1:initial_n) = conjg(shifted_column(1:initial_n))
    call state%append_symmetric(h_column, s_column, info)
    call check(info == QR_SUCCESS .and. state%n == capacity, &
               'complex append/delete round trip appends successfully', failures)
    successful_updates = successful_updates + 1
    old_n = state%n
    call state%delete_symmetric(capacity, info)
    call remove_complex_index(expected, old_n, capacity)
    successful_updates = successful_updates + 1
    call complex_factor_errors(state, expected, initial_n, residual, &
                               unitarity, lower_error)
    call check(info == QR_SUCCESS .and. state%n == initial_n .and. &
               residual <= tolerance .and. unitarity <= tolerance .and. &
               lower_error <= tolerance, &
               'complex append/delete round trip restores the original matrix', &
               failures)

    q_before = state%q
    r_before = state%r
    call state%delete_symmetric(0, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex deletion rejects index zero', failures)
    call state%delete_symmetric(state%n + 1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 6 + successful_updates .and. &
               state%updates_since_fresh == 3 + successful_updates, &
               'invalid complex deletion preserves the complete state', &
               failures)

    factors_accurate = .true.
    inactive_storage_cleared = .true.
    max_residual = residual
    max_unitarity = unitarity
    do delete_number = 1, num_deletes
      old_n = state%n
      idx = delete_indices(delete_number)
      new_n = old_n - 1
      call state%delete_symmetric(idx, info)
      call remove_complex_index(expected, old_n, idx)
      successful_updates = successful_updates + 1
      call complex_factor_errors(state, expected, new_n, residual, &
                                 unitarity, lower_error)
      factors_accurate = factors_accurate .and. info == QR_SUCCESS .and. &
        state%n == new_n .and. state%valid .and. &
        residual <= tolerance .and. unitarity <= tolerance .and. &
        lower_error <= tolerance .and. &
        state%structural_updates == 6 + successful_updates .and. &
        state%updates_since_fresh == 3 + successful_updates
      inactive_storage_cleared = inactive_storage_cleared .and. &
        all(abs(state%q(old_n,:)) <= 0.0_wp) .and. &
        all(abs(state%q(:,old_n)) <= 0.0_wp) .and. &
        all(abs(state%r(old_n,:)) <= 0.0_wp) .and. &
        all(abs(state%r(:,old_n)) <= 0.0_wp)
      max_residual = max(max_residual, residual)
      max_unitarity = max(max_unitarity, unitarity)
    end do
    call check(factors_accurate, &
               'every complex deletion preserves accurate factors and counters', &
               failures)
    call check(inactive_storage_cleared, &
               'complex deletion clears inactive factor rows and columns', &
               failures)
    call check(state%n == 1 .and. state%valid .and. &
               abs(state%shift - shift) <= 0.0_wp, &
               'complex deletion retains a valid order-one factorization', &
               failures)

    q_before = state%q
    r_before = state%r
    call state%delete_symmetric(1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%n == 1 .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 6 + successful_updates .and. &
               state%updates_since_fresh == 3 + successful_updates, &
               'complex order-one deletion is rejected without mutation', &
               failures)
    call empty_state%delete_symmetric(1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex deletion rejects an unfactorized state', failures)

    write(*,'(a,i0,2(a,es12.4))') '  complex deletion wp=', wp, &
      ' max residual=', max_residual, ' max unitarity=', max_unitarity
  end subroutine test_complex_delete

  ! Remove one row and column from a real reference matrix in place. The
  ! retained principal submatrix is packed into the leading (n-1)-square block.
  subroutine remove_real_index(matrix, n, idx)
    real(wp), intent(inout) :: matrix(:,:)
    integer, intent(in) :: n, idx
    integer :: i, j

    do j = idx, n - 1
      do i = 1, n
        matrix(i,j) = matrix(i,j + 1)
      end do
    end do
    do i = idx, n - 1
      do j = 1, n - 1
        matrix(i,j) = matrix(i + 1,j)
      end do
    end do
    matrix(n,:) = 0.0_wp
    matrix(:,n) = 0.0_wp
  end subroutine remove_real_index

  ! Complex counterpart of remove_real_index. Copying complete retained rows
  ! and columns preserves the independently constructed Hermitian reference.
  subroutine remove_complex_index(matrix, n, idx)
    complex(wp), intent(inout) :: matrix(:,:)
    integer, intent(in) :: n, idx
    integer :: i, j

    do j = idx, n - 1
      do i = 1, n
        matrix(i,j) = matrix(i,j + 1)
      end do
    end do
    do i = idx, n - 1
      do j = 1, n - 1
        matrix(i,j) = matrix(i + 1,j)
      end do
    end do
    matrix(n,:) = cmplx(0.0_wp, 0.0_wp, kind=wp)
    matrix(:,n) = cmplx(0.0_wp, 0.0_wp, kind=wp)
  end subroutine remove_complex_index

  ! Compute sign-independent quality measures for an active real QR state.
  subroutine real_factor_errors(state, matrix, n, residual, orthogonality, &
                                lower_error)
    type(qr_real_state), intent(in) :: state
    real(wp), intent(in) :: matrix(:,:)
    integer, intent(in) :: n
    real(wp), intent(out) :: residual, orthogonality, lower_error
    real(wp) :: identity(n,n)
    integer :: i, j

    identity = 0.0_wp
    do i = 1, n
      identity(i,i) = 1.0_wp
    end do
    residual = real_frobenius(matrix(1:n,1:n) - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(real_frobenius(matrix(1:n,1:n)), tiny(1.0_wp))
    orthogonality = real_frobenius(identity - &
      matmul(transpose(state%q(1:n,1:n)), state%q(1:n,1:n))) / &
      sqrt(real(n,wp))
    lower_error = 0.0_wp
    do j = 1, n - 1
      do i = j + 1, n
        lower_error = max(lower_error, abs(state%r(i,j)))
      end do
    end do
  end subroutine real_factor_errors

  ! Compute phase-independent reconstruction, unitarity, and triangularity
  ! measures for an active complex QR state.
  subroutine complex_factor_errors(state, matrix, n, residual, unitarity, &
                                   lower_error)
    type(qr_complex_state), intent(in) :: state
    complex(wp), intent(in) :: matrix(:,:)
    integer, intent(in) :: n
    real(wp), intent(out) :: residual, unitarity, lower_error
    complex(wp) :: identity(n,n)
    integer :: i, j

    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do i = 1, n
      identity(i,i) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do
    residual = complex_frobenius(matrix(1:n,1:n) - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(complex_frobenius(matrix(1:n,1:n)), tiny(1.0_wp))
    unitarity = complex_frobenius(identity - &
      matmul(conjg(transpose(state%q(1:n,1:n))), &
             state%q(1:n,1:n))) / sqrt(real(n,wp))
    lower_error = 0.0_wp
    do j = 1, n - 1
      do i = j + 1, n
        lower_error = max(lower_error, abs(state%r(i,j)))
      end do
    end do
  end subroutine complex_factor_errors

end program test_delete

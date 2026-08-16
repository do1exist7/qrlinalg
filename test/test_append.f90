! Analytical tests for symmetric and Hermitian order-increasing QR updates.
!
! The tests begin with order-two fresh factorizations and append three rows and
! columns to reach capacity. An independent shifted matrix is expanded after
! every operation and compared with Q*R; no expected factor signs or phases are
! assumed.
program test_append
  use qrlinalg
  use test_support, only: check, finish_suite, real_frobenius, &
                          complex_frobenius
  implicit none

  integer :: failures

  failures = 0
  call test_real_append(failures)
  call test_complex_append(failures)
  call finish_suite('qrlinalg symmetric append', failures)

contains

  ! Append three real symmetric rows and columns. Each new physical H and S
  ! column is formed from a prescribed shifted column, so the reference matrix
  ! follows directly from H-shift*S without invoking an independent QR method.
  subroutine test_real_append(failures)
    integer, intent(inout) :: failures
    integer, parameter :: initial_n = 2, capacity = 5, num_appends = 3
    type(qr_real_state) :: state, empty_state
    real(wp), parameter :: shift = 0.3_wp
    real(wp) :: h(initial_n,initial_n), s(initial_n,initial_n)
    real(wp) :: expected(capacity,capacity), identity(capacity,capacity)
    real(wp) :: shifted_column(capacity), h_column(capacity)
    real(wp) :: s_column(capacity), h_before(capacity), s_before(capacity)
    real(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: residual, orthogonality, lower_error, tolerance
    real(wp) :: max_residual, max_orthogonality, max_lower_error
    integer :: append_number, i, info, j, new_n, old_n
    logical :: caller_vectors_preserved, factors_accurate

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    expected = 0.0_wp
    expected(1:initial_n,1:initial_n) = &
      reshape([4.0_wp, 0.5_wp, 0.5_wp, 3.0_wp], [initial_n,initial_n])
    s = reshape([1.5_wp, 0.1_wp, 0.1_wp, 1.2_wp], &
                [initial_n,initial_n])
    h = expected(1:initial_n,1:initial_n) + shift * s
    identity = 0.0_wp
    do i = 1, capacity
      identity(i,i) = 1.0_wp
    end do

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'real append state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'real append base factorizes', failures)
    state%structural_updates = 3
    state%updates_since_fresh = 1

    ! A wrong extent must be rejected before qrinc changes either factor.
    h_column = 0.0_wp
    s_column = 0.0_wp
    q_before = state%q
    r_before = state%r
    call state%append_symmetric(h_column(1:initial_n), &
                                s_column(1:initial_n), info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%n == initial_n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp), &
               'real append rejects wrong extents without changing state', &
               failures)

    caller_vectors_preserved = .true.
    factors_accurate = .true.
    max_residual = 0.0_wp
    max_orthogonality = 0.0_wp
    max_lower_error = 0.0_wp
    do append_number = 1, num_appends
      old_n = state%n
      new_n = old_n + 1
      do i = 1, old_n
        shifted_column(i) = 0.04_wp * &
          real(mod(append_number * (2 * i + 1), 9) - 4, wp)
        s_column(i) = 0.01_wp * &
          real(mod((append_number + 1) * (i + 2), 7) - 3, wp)
      end do
      shifted_column(new_n) = 4.0_wp + 0.5_wp * real(new_n, wp)
      s_column(new_n) = 1.0_wp + 0.1_wp * real(new_n, wp)
      h_column(1:new_n) = shifted_column(1:new_n) + &
                          shift * s_column(1:new_n)
      h_before(1:new_n) = h_column(1:new_n)
      s_before(1:new_n) = s_column(1:new_n)

      expected(1:new_n,new_n) = shifted_column(1:new_n)
      expected(new_n,1:old_n) = shifted_column(1:old_n)
      call state%append_symmetric(h_column(1:new_n), &
                                  s_column(1:new_n), info)
      caller_vectors_preserved = caller_vectors_preserved .and. &
        all(abs(h_column(1:new_n) - h_before(1:new_n)) <= 0.0_wp) .and. &
        all(abs(s_column(1:new_n) - s_before(1:new_n)) <= 0.0_wp)

      residual = real_frobenius(expected(1:new_n,1:new_n) - &
        matmul(state%q(1:new_n,1:new_n), &
               state%r(1:new_n,1:new_n))) / &
        max(real_frobenius(expected(1:new_n,1:new_n)), tiny(1.0_wp))
      orthogonality = real_frobenius(identity(1:new_n,1:new_n) - &
        matmul(transpose(state%q(1:new_n,1:new_n)), &
               state%q(1:new_n,1:new_n))) / sqrt(real(new_n,wp))
      lower_error = 0.0_wp
      do j = 1, new_n - 1
        do i = j + 1, new_n
          lower_error = max(lower_error, abs(state%r(i,j)))
        end do
      end do
      factors_accurate = factors_accurate .and. info == QR_SUCCESS .and. &
        state%n == new_n .and. state%valid .and. &
        residual <= tolerance .and. orthogonality <= tolerance .and. &
        lower_error <= tolerance
      max_residual = max(max_residual, residual)
      max_orthogonality = max(max_orthogonality, orthogonality)
      max_lower_error = max(max_lower_error, lower_error)
    end do

    call check(factors_accurate, &
               'every real append preserves an accurate QR factorization', &
               failures)
    call check(caller_vectors_preserved, &
               'every real append preserves caller columns', failures)
    call check(state%n == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp .and. &
               state%structural_updates == 3 + num_appends .and. &
               state%updates_since_fresh == 1 + num_appends, &
               'real appends commit order, shift, and counters', failures)

    ! Capacity rejection occurs before either numerical insertion and leaves
    ! even the inactive portions of the factor arrays unchanged.
    q_before = state%q
    r_before = state%r
    call state%append_symmetric(h_column, s_column, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%n == capacity .and. &
               state%structural_updates == 3 + num_appends, &
               'real append at capacity preserves the complete state', failures)
    call empty_state%append_symmetric(h_column(1:3), s_column(1:3), info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'real append rejects an unfactorized state', failures)

    write(*,'(a,i0,3(a,es12.4))') '  real append wp=', wp, &
      ' max residual=', max_residual, &
      ' max orthogonality=', max_orthogonality, &
      ' max lower(R)=', max_lower_error
  end subroutine test_real_append

  ! Hermitian counterpart. Off-diagonal input columns contain changing complex
  ! values; qrinr must form their conjugate rows while retaining a real
  ! diagonal and a unitary explicit Q.
  subroutine test_complex_append(failures)
    integer, intent(inout) :: failures
    integer, parameter :: initial_n = 2, capacity = 5, num_appends = 3
    type(qr_complex_state) :: state
    real(wp), parameter :: shift = -0.2_wp
    complex(wp) :: h(initial_n,initial_n), s(initial_n,initial_n)
    complex(wp) :: expected(capacity,capacity), identity(capacity,capacity)
    complex(wp) :: shifted_column(capacity), h_column(capacity)
    complex(wp) :: s_column(capacity), h_before(capacity), s_before(capacity)
    complex(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: residual, unitarity, lower_error, hermitian_error, tolerance
    real(wp) :: max_residual, max_unitarity, max_lower_error
    integer :: append_number, i, info, j, new_n, old_n
    logical :: caller_vectors_preserved, factors_accurate

    tolerance = 100000.0_wp * epsilon(1.0_wp)
    expected = cmplx(0.0_wp, 0.0_wp, kind=wp)
    expected(1,1) = cmplx(4.0_wp, 0.0_wp, kind=wp)
    expected(1,2) = cmplx(0.5_wp, -0.3_wp, kind=wp)
    expected(2,1) = conjg(expected(1,2))
    expected(2,2) = cmplx(3.0_wp, 0.0_wp, kind=wp)
    s(1,1) = cmplx(1.5_wp, 0.0_wp, kind=wp)
    s(1,2) = cmplx(0.1_wp, 0.05_wp, kind=wp)
    s(2,1) = conjg(s(1,2))
    s(2,2) = cmplx(1.2_wp, 0.0_wp, kind=wp)
    h = expected(1:initial_n,1:initial_n) + &
        cmplx(shift, 0.0_wp, kind=wp) * s
    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do i = 1, capacity
      identity(i,i) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do

    call state%initialize(capacity, info)
    call check(info == QR_SUCCESS, 'complex append state initializes', failures)
    call state%factorize_fresh(h, s, shift, info)
    call check(info == QR_SUCCESS, 'complex append base factorizes', failures)
    state%structural_updates = 5
    state%updates_since_fresh = 2

    ! H and S diagonals are independently Hermitian requirements. Check each
    ! rejection, including complete factor and counter preservation.
    h_column = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s_column = cmplx(0.0_wp, 0.0_wp, kind=wp)
    h_column(3) = cmplx(2.0_wp, 0.01_wp, kind=wp)
    s_column(3) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    q_before = state%q
    r_before = state%r
    call state%append_symmetric(h_column(1:3), s_column(1:3), info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex append rejects a non-real H diagonal', failures)
    h_column(3) = cmplx(2.0_wp, 0.0_wp, kind=wp)
    s_column(3) = cmplx(1.0_wp, -0.01_wp, kind=wp)
    call state%append_symmetric(h_column(1:3), s_column(1:3), info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%n == initial_n .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%structural_updates == 5 .and. &
               state%updates_since_fresh == 2, &
               'invalid complex diagonals preserve the complete state', failures)

    caller_vectors_preserved = .true.
    factors_accurate = .true.
    max_residual = 0.0_wp
    max_unitarity = 0.0_wp
    max_lower_error = 0.0_wp
    do append_number = 1, num_appends
      old_n = state%n
      new_n = old_n + 1
      do i = 1, old_n
        shifted_column(i) = cmplx( &
          0.03_wp * real(mod(append_number * (i + 2), 9) - 4, wp), &
          0.02_wp * real(mod((append_number + 1) * (2 * i + 1), 7) - 3, wp), &
          kind=wp)
        s_column(i) = cmplx( &
          0.008_wp * real(mod((append_number + 2) * (i + 1), 7) - 3, wp), &
          0.006_wp * real(mod((append_number + 3) * (i + 2), 5) - 2, wp), &
          kind=wp)
      end do
      shifted_column(new_n) = &
        cmplx(4.0_wp + 0.5_wp * real(new_n,wp), 0.0_wp, kind=wp)
      s_column(new_n) = &
        cmplx(1.0_wp + 0.1_wp * real(new_n,wp), 0.0_wp, kind=wp)
      h_column(1:new_n) = shifted_column(1:new_n) + &
        cmplx(shift, 0.0_wp, kind=wp) * s_column(1:new_n)
      if (append_number == 1) then
        ! Roundoff-sized imaginary diagonal parts are accepted but discarded
        ! by the library, yielding an exactly Hermitian represented matrix.
        h_column(new_n) = h_column(new_n) + &
          cmplx(0.0_wp, 10.0_wp * epsilon(1.0_wp), kind=wp)
        s_column(new_n) = s_column(new_n) - &
          cmplx(0.0_wp, 12.0_wp * epsilon(1.0_wp), kind=wp)
      end if
      h_before(1:new_n) = h_column(1:new_n)
      s_before(1:new_n) = s_column(1:new_n)

      expected(1:new_n,new_n) = shifted_column(1:new_n)
      expected(new_n,1:old_n) = conjg(shifted_column(1:old_n))
      call state%append_symmetric(h_column(1:new_n), &
                                  s_column(1:new_n), info)
      caller_vectors_preserved = caller_vectors_preserved .and. &
        all(abs(h_column(1:new_n) - h_before(1:new_n)) <= 0.0_wp) .and. &
        all(abs(s_column(1:new_n) - s_before(1:new_n)) <= 0.0_wp)

      residual = complex_frobenius(expected(1:new_n,1:new_n) - &
        matmul(state%q(1:new_n,1:new_n), &
               state%r(1:new_n,1:new_n))) / &
        max(complex_frobenius(expected(1:new_n,1:new_n)), tiny(1.0_wp))
      unitarity = complex_frobenius(identity(1:new_n,1:new_n) - &
        matmul(conjg(transpose(state%q(1:new_n,1:new_n))), &
               state%q(1:new_n,1:new_n))) / sqrt(real(new_n,wp))
      lower_error = 0.0_wp
      do j = 1, new_n - 1
        do i = j + 1, new_n
          lower_error = max(lower_error, abs(state%r(i,j)))
        end do
      end do
      hermitian_error = maxval(abs(expected(1:new_n,1:new_n) - &
        conjg(transpose(expected(1:new_n,1:new_n)))))
      factors_accurate = factors_accurate .and. info == QR_SUCCESS .and. &
        state%n == new_n .and. state%valid .and. &
        residual <= tolerance .and. unitarity <= tolerance .and. &
        lower_error <= tolerance .and. hermitian_error <= tolerance
      max_residual = max(max_residual, residual)
      max_unitarity = max(max_unitarity, unitarity)
      max_lower_error = max(max_lower_error, lower_error)
    end do

    call check(factors_accurate, &
               'every complex append preserves an accurate QR factorization', &
               failures)
    call check(caller_vectors_preserved, &
               'every complex append preserves caller columns', failures)
    call check(state%n == capacity .and. &
               abs(state%shift - shift) <= 0.0_wp .and. &
               state%structural_updates == 5 + num_appends .and. &
               state%updates_since_fresh == 2 + num_appends, &
               'complex appends commit order, shift, and counters', failures)

    q_before = state%q
    r_before = state%r
    call state%append_symmetric(h_column, s_column, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%n == capacity .and. &
               state%structural_updates == 5 + num_appends, &
               'complex append at capacity preserves the complete state', &
               failures)

    write(*,'(a,i0,3(a,es12.4))') '  complex append wp=', wp, &
      ' max residual=', max_residual, ' max unitarity=', max_unitarity, &
      ' max lower(R)=', max_lower_error
  end subroutine test_complex_append

end program test_append

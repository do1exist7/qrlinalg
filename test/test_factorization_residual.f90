! Tests for capacity-sized matrix inputs and the public factorization residual.
program test_factorization_residual
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, &
                                                ieee_quiet_nan, ieee_value
  use qrlinalg
  use test_support, only: check, finish_suite
  implicit none

  integer :: failures

  failures = 0
  call test_real_leading_dimension(failures)
  call test_complex_leading_dimension(failures)
  call test_real_residual_lifecycle(failures)
  call test_complex_residual_lifecycle(failures)
  call test_residual_validation(failures)
  call finish_suite('qrlinalg factorization residual', failures)

contains

  subroutine test_real_leading_dimension(failures)
  !Verify that real factorization and solve use the active principal block and
  !the physical leading dimension of capacity-sized caller matrices. Quiet NaN
  !values outside the active lower triangle make any unintended numerical use
  !of upper or inactive storage visible in the returned results.
    integer, intent(inout) :: failures
    integer, parameter :: active_n = 3, capacity = 6
    type(qr_real_state) :: full_state, packed_state
    real(wp), parameter :: shift = 0.8_wp
    real(wp) :: full_h(capacity,capacity), full_s(capacity,capacity)
    real(wp) :: packed_h(active_n,active_n), packed_s(active_n,active_n)
    real(wp) :: initial(active_n), full_x(active_n), packed_x(active_n)
    real(wp) :: full_lambda, packed_lambda, full_accuracy, packed_accuracy
    real(wp) :: absolute_residual, relative_residual, nan_value, tolerance
    integer :: full_iterations, info, packed_iterations, i, j

    nan_value = ieee_value(0.0_wp, ieee_quiet_nan)
    full_h = nan_value
    full_s = nan_value
    packed_h = nan_value
    packed_s = nan_value
    do j = 1, active_n
      do i = j, active_n
        full_h(i,j) = 0.0_wp
        full_s(i,j) = 0.0_wp
        packed_h(i,j) = 0.0_wp
        packed_s(i,j) = 0.0_wp
      end do
      full_h(j,j) = real(2 * j - 1, wp)
      full_s(j,j) = 1.0_wp
      packed_h(j,j) = full_h(j,j)
      packed_s(j,j) = full_s(j,j)
    end do
    initial = [1.0_wp, 0.5_wp, -0.25_wp]
    tolerance = 10000.0_wp * epsilon(1.0_wp)

    call full_state%initialize(capacity, info)
    call full_state%factorize_fresh(full_h, full_s, shift, info, &
                                    active_order=active_n)
    call check(info == QR_SUCCESS .and. full_state%order() == active_n, &
               'real active-order factorization accepts capacity storage', &
               failures)
    call full_state%factorization_residual(full_h, full_s, initial, &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'real residual ignores inactive and upper NaN storage', failures)

    call packed_state%initialize(active_n, info)
    call packed_state%factorize_fresh(packed_h, packed_s, shift, info)
    call full_state%solve(full_s, initial, full_x, full_lambda, tolerance, 100, &
                          1, full_accuracy, full_iterations, info)
    call check(info == QR_SUCCESS, &
               'real solve accepts capacity-sized overlap storage', failures)
    call packed_state%solve(packed_s, initial, packed_x, packed_lambda, &
                            tolerance, 100, 1, packed_accuracy, &
                            packed_iterations, info)
    call check(info == QR_SUCCESS .and. &
               abs(full_lambda - packed_lambda) <= tolerance .and. &
               maxval(abs(abs(full_x) - abs(packed_x))) <= tolerance, &
               'real capacity and packed solves agree', failures)
    call check(ieee_is_nan(full_h(capacity,1)) .and. &
               ieee_is_nan(full_s(capacity,capacity)) .and. &
               ieee_is_finite(full_lambda), &
               'real inactive sentinels remain untouched', failures)
  end subroutine test_real_leading_dimension

  subroutine test_complex_leading_dimension(failures)
  !Apply the leading-dimension and inactive-storage checks to the Hermitian
  !interface. Complex NaNs in unused storage also detect accidental use of the
  !upper triangle or rows below the active principal block.
    integer, intent(inout) :: failures
    integer, parameter :: active_n = 3, capacity = 6
    type(qr_complex_state) :: full_state, packed_state
    real(wp), parameter :: shift = 0.8_wp
    complex(wp) :: full_h(capacity,capacity), full_s(capacity,capacity)
    complex(wp) :: packed_h(active_n,active_n), packed_s(active_n,active_n)
    complex(wp) :: initial(active_n), full_x(active_n), packed_x(active_n)
    complex(wp) :: nan_value
    real(wp) :: full_lambda, packed_lambda, full_accuracy, packed_accuracy
    real(wp) :: absolute_residual, relative_residual, real_nan, tolerance
    integer :: full_iterations, info, packed_iterations, i, j

    real_nan = ieee_value(0.0_wp, ieee_quiet_nan)
    nan_value = cmplx(real_nan, real_nan, kind=wp)
    full_h = nan_value
    full_s = nan_value
    packed_h = nan_value
    packed_s = nan_value
    do j = 1, active_n
      do i = j, active_n
        full_h(i,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        full_s(i,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        packed_h(i,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
        packed_s(i,j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      end do
      full_h(j,j) = cmplx(real(2 * j - 1, wp), 0.0_wp, kind=wp)
      full_s(j,j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
      packed_h(j,j) = full_h(j,j)
      packed_s(j,j) = full_s(j,j)
    end do
    initial = [cmplx(1.0_wp, 0.25_wp, kind=wp), &
               cmplx(0.5_wp, -0.1_wp, kind=wp), &
               cmplx(-0.25_wp, 0.2_wp, kind=wp)]
    tolerance = 10000.0_wp * epsilon(1.0_wp)

    call full_state%initialize(capacity, info)
    call full_state%factorize_fresh(full_h, full_s, shift, info, &
                                    active_order=active_n)
    call check(info == QR_SUCCESS .and. full_state%order() == active_n, &
               'complex active-order factorization accepts capacity storage', &
               failures)
    call full_state%factorization_residual(full_h, full_s, initial, &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'complex residual ignores inactive and upper NaN storage', &
               failures)

    call packed_state%initialize(active_n, info)
    call packed_state%factorize_fresh(packed_h, packed_s, shift, info)
    call full_state%solve(full_s, initial, full_x, full_lambda, tolerance, 100, &
                          1, full_accuracy, full_iterations, info)
    call check(info == QR_SUCCESS, &
               'complex solve accepts capacity-sized overlap storage', failures)
    call packed_state%solve(packed_s, initial, packed_x, packed_lambda, &
                            tolerance, 100, 1, packed_accuracy, &
                            packed_iterations, info)
    call check(info == QR_SUCCESS .and. &
               abs(full_lambda - packed_lambda) <= tolerance .and. &
               maxval(abs(abs(full_x) - abs(packed_x))) <= tolerance, &
               'complex capacity and packed solves agree', failures)
    call check(ieee_is_nan(real(full_h(capacity,1), wp)) .and. &
               ieee_is_nan(aimag(full_s(capacity,capacity))) .and. &
               ieee_is_finite(full_lambda), &
               'complex inactive sentinels remain untouched', failures)
  end subroutine test_complex_leading_dimension

  subroutine test_real_residual_lifecycle(failures)
  !Measure a real factorization after fresh construction, one replacement, one
  !hundred further replacements, append, arbitrary deletion, and a deliberate
  !caller/factor mismatch. A separately rebuilt state supplies the fresh-QR
  !comparison after the repeated updates.
    integer, intent(inout) :: failures
    integer, parameter :: capacity = 7, initial_n = 4, num_updates = 100
    type(qr_real_state) :: state, fresh_state
    real(wp), parameter :: shift = 0.35_wp
    real(wp) :: h(capacity,capacity), s(capacity,capacity), v(capacity)
    real(wp) :: h_before(capacity,capacity), s_before(capacity,capacity)
    real(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    real(wp) :: v_before(capacity)
    real(wp) :: delta_h(capacity), delta_s(capacity)
    real(wp) :: h_column(capacity), s_column(capacity)
    real(wp) :: absolute_residual, relative_residual, updated_residual
    real(wp) :: fresh_absolute, fresh_residual, tolerance
    integer :: info, i, idx, j, matrix_n, update_number

    h = 0.0_wp
    s = 0.0_wp
    do j = 1, initial_n
      h(j,j) = 3.0_wp + real(j,wp)
      s(j,j) = 1.0_wp + 0.05_wp * real(j,wp)
      do i = j + 1, initial_n
        h(i,j) = 0.04_wp * real(i + j,wp)
        h(j,i) = h(i,j)
        s(i,j) = 0.004_wp * real(i - j,wp)
        s(j,i) = s(i,j)
      end do
    end do
    do i = 1, capacity
      v(i) = 0.5_wp + 0.1_wp * real(i,wp)
    end do
    tolerance = 1000000.0_wp * epsilon(1.0_wp)
    matrix_n = initial_n

    call state%initialize(capacity, info)
    call state%factorize_fresh(h, s, shift, info, active_order=matrix_n)
    h_before = h
    s_before = s
    v_before = v
    q_before = state%q
    r_before = state%r
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'fresh real factorization has roundoff-scale action residual', &
               failures)
    call check(all(abs(h - h_before) <= 0.0_wp) .and. &
               all(abs(s - s_before) <= 0.0_wp) .and. &
               all(abs(v - v_before) <= 0.0_wp) .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%is_valid() .and. state%order() == matrix_n .and. &
               abs(state%get_shift() - shift) <= 0.0_wp .and. &
               state%get_update_count() == 0 .and. &
               state%get_updates_since_fresh() == 0, &
               'real residual preserves inputs, factors, and metadata', failures)

    idx = 2
    delta_h(1:matrix_n) = [0.02_wp, 0.01_wp, -0.015_wp, 0.005_wp]
    delta_s(1:matrix_n) = [0.001_wp, 0.002_wp, -0.001_wp, 0.0005_wp]
    call state%replace_symmetric(idx, delta_h(1:matrix_n), &
                                 delta_s(1:matrix_n), info)
    call apply_real_replacement(h, s, matrix_n, idx, delta_h, delta_s)
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'real residual remains small after synchronized replacement', &
               failures)

    do update_number = 1, num_updates
      idx = 1 + mod(update_number - 1, matrix_n)
      do i = 1, matrix_n
        delta_h(i) = 0.00002_wp * &
          real(mod(update_number * (i + 1), 9) - 4, wp)
        delta_s(i) = 0.000002_wp * &
          real(mod((update_number + 2) * i, 7) - 3, wp)
      end do
      call state%replace_symmetric(idx, delta_h(1:matrix_n), &
                                   delta_s(1:matrix_n), info)
      call apply_real_replacement(h, s, matrix_n, idx, delta_h, delta_s)
    end do
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, updated_residual, info)
    call fresh_state%initialize(capacity, info)
    call fresh_state%factorize_fresh(h, s, shift, info, active_order=matrix_n)
    call fresh_state%factorization_residual(h, s, v(1:matrix_n), &
      fresh_absolute, fresh_residual, info)
    call check(updated_residual <= tolerance .and. &
               fresh_residual <= tolerance .and. &
               updated_residual <= max(1000.0_wp * fresh_residual, tolerance), &
               'many real updates remain close to a fresh factorization', &
               failures)

    matrix_n = matrix_n + 1
    do i = 1, matrix_n - 1
      h_column(i) = 0.03_wp * real(i - 2,wp)
      s_column(i) = 0.002_wp * real(i,wp)
    end do
    h_column(matrix_n) = 7.5_wp
    s_column(matrix_n) = 1.3_wp
    call state%append_symmetric(h_column(1:matrix_n), &
                                s_column(1:matrix_n), info)
    call append_real_column(h, s, matrix_n, h_column, s_column)
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'real residual remains small after append', failures)

    idx = 3
    call state%delete_symmetric(idx, info)
    call remove_real_index(h, matrix_n, idx)
    call remove_real_index(s, matrix_n, idx)
    matrix_n = matrix_n - 1
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'real residual remains small after arbitrary deletion', failures)

    h(1,1) = h(1,1) + 0.5_wp
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual > 1.0e-5_wp, &
               'real residual detects caller/factor mismatch', failures)
  end subroutine test_real_residual_lifecycle

  subroutine test_complex_residual_lifecycle(failures)
  !Exercise the corresponding Hermitian update sequence with genuinely complex
  !off-diagonal changes so conjugate row reconstruction is tested throughout.
    integer, intent(inout) :: failures
    integer, parameter :: capacity = 7, initial_n = 4, num_updates = 100
    type(qr_complex_state) :: state, fresh_state
    real(wp), parameter :: shift = -0.2_wp
    complex(wp) :: h(capacity,capacity), s(capacity,capacity), v(capacity)
    complex(wp) :: h_before(capacity,capacity), s_before(capacity,capacity)
    complex(wp) :: q_before(capacity,capacity), r_before(capacity,capacity)
    complex(wp) :: v_before(capacity)
    complex(wp) :: delta_h(capacity), delta_s(capacity)
    complex(wp) :: h_column(capacity), s_column(capacity)
    real(wp) :: absolute_residual, relative_residual, updated_residual
    real(wp) :: fresh_absolute, fresh_residual, tolerance
    integer :: info, i, idx, j, matrix_n, update_number

    h = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do j = 1, initial_n
      h(j,j) = cmplx(3.0_wp + real(j,wp), 0.0_wp, kind=wp)
      s(j,j) = cmplx(1.0_wp + 0.05_wp * real(j,wp), 0.0_wp, kind=wp)
      do i = j + 1, initial_n
        h(i,j) = cmplx(0.03_wp * real(i + j,wp), &
                        0.02_wp * real(i - j,wp), kind=wp)
        h(j,i) = conjg(h(i,j))
        s(i,j) = cmplx(0.003_wp * real(i - j,wp), &
                        0.001_wp * real(i + j,wp), kind=wp)
        s(j,i) = conjg(s(i,j))
      end do
    end do
    do i = 1, capacity
      v(i) = cmplx(0.5_wp + 0.1_wp * real(i,wp), &
                   0.08_wp * real(i - 2,wp), kind=wp)
    end do
    tolerance = 1000000.0_wp * epsilon(1.0_wp)
    matrix_n = initial_n

    call state%initialize(capacity, info)
    call state%factorize_fresh(h, s, shift, info, active_order=matrix_n)
    h_before = h
    s_before = s
    v_before = v
    q_before = state%q
    r_before = state%r
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'fresh complex factorization has roundoff-scale action residual', &
               failures)
    call check(all(abs(h - h_before) <= 0.0_wp) .and. &
               all(abs(s - s_before) <= 0.0_wp) .and. &
               all(abs(v - v_before) <= 0.0_wp) .and. &
               all(abs(state%q - q_before) <= 0.0_wp) .and. &
               all(abs(state%r - r_before) <= 0.0_wp) .and. &
               state%is_valid() .and. state%order() == matrix_n .and. &
               abs(state%get_shift() - shift) <= 0.0_wp .and. &
               state%get_update_count() == 0 .and. &
               state%get_updates_since_fresh() == 0, &
               'complex residual preserves inputs, factors, and metadata', &
               failures)

    idx = 2
    delta_h(1:matrix_n) = [ &
      cmplx(0.02_wp, 0.01_wp, kind=wp), &
      cmplx(0.01_wp, 0.0_wp, kind=wp), &
      cmplx(-0.015_wp, 0.008_wp, kind=wp), &
      cmplx(0.005_wp, -0.004_wp, kind=wp)]
    delta_s(1:matrix_n) = [ &
      cmplx(0.001_wp, -0.0005_wp, kind=wp), &
      cmplx(0.002_wp, 0.0_wp, kind=wp), &
      cmplx(-0.001_wp, 0.0003_wp, kind=wp), &
      cmplx(0.0005_wp, 0.0002_wp, kind=wp)]
    call state%replace_symmetric(idx, delta_h(1:matrix_n), &
                                 delta_s(1:matrix_n), info)
    call apply_complex_replacement(h, s, matrix_n, idx, delta_h, delta_s)
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'complex residual remains small after replacement', failures)

    do update_number = 1, num_updates
      idx = 1 + mod(update_number - 1, matrix_n)
      do i = 1, matrix_n
        delta_h(i) = cmplx( &
          0.00002_wp * real(mod(update_number * (i + 1), 9) - 4, wp), &
          0.00001_wp * real(mod(update_number + 2 * i, 7) - 3, wp), kind=wp)
        delta_s(i) = cmplx( &
          0.000002_wp * real(mod((update_number + 2) * i, 7) - 3, wp), &
          0.000001_wp * real(mod(update_number + i, 5) - 2, wp), kind=wp)
      end do
      delta_h(idx) = cmplx(real(delta_h(idx), wp), 0.0_wp, kind=wp)
      delta_s(idx) = cmplx(real(delta_s(idx), wp), 0.0_wp, kind=wp)
      call state%replace_symmetric(idx, delta_h(1:matrix_n), &
                                   delta_s(1:matrix_n), info)
      call apply_complex_replacement(h, s, matrix_n, idx, delta_h, delta_s)
    end do
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, updated_residual, info)
    call fresh_state%initialize(capacity, info)
    call fresh_state%factorize_fresh(h, s, shift, info, active_order=matrix_n)
    call fresh_state%factorization_residual(h, s, v(1:matrix_n), &
      fresh_absolute, fresh_residual, info)
    call check(updated_residual <= tolerance .and. &
               fresh_residual <= tolerance .and. &
               updated_residual <= max(1000.0_wp * fresh_residual, tolerance), &
               'many complex updates remain close to fresh factorization', &
               failures)

    matrix_n = matrix_n + 1
    do i = 1, matrix_n - 1
      h_column(i) = cmplx(0.025_wp * real(i - 2,wp), &
                           0.012_wp * real(i,wp), kind=wp)
      s_column(i) = cmplx(0.002_wp * real(i,wp), &
                           -0.001_wp * real(i - 1,wp), kind=wp)
    end do
    h_column(matrix_n) = cmplx(7.5_wp, 0.0_wp, kind=wp)
    s_column(matrix_n) = cmplx(1.3_wp, 0.0_wp, kind=wp)
    call state%append_symmetric(h_column(1:matrix_n), &
                                s_column(1:matrix_n), info)
    call append_complex_column(h, s, matrix_n, h_column, s_column)
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'complex residual remains small after append', failures)

    idx = 3
    call state%delete_symmetric(idx, info)
    call remove_complex_index(h, matrix_n, idx)
    call remove_complex_index(s, matrix_n, idx)
    matrix_n = matrix_n - 1
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual <= tolerance, &
               'complex residual remains small after arbitrary deletion', &
               failures)

    h(1,1) = h(1,1) + cmplx(0.5_wp, 0.0_wp, kind=wp)
    call state%factorization_residual(h, s, v(1:matrix_n), &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. relative_residual > 1.0e-5_wp, &
               'complex residual detects caller/factor mismatch', failures)
  end subroutine test_complex_residual_lifecycle

  subroutine test_residual_validation(failures)
  !Check deterministic outputs and symbolic statuses for invalid states,
  !undersized matrices, invalid explicit orders, and zero probe vectors. A
  !small residual after every rejected factorization proves that old valid
  !factors remain represented.
    integer, intent(inout) :: failures
    integer, parameter :: n = 2, capacity = 4
    type(qr_real_state) :: empty_state, state, zero_state
    type(qr_complex_state) :: complex_empty, complex_state, complex_zero_state
    real(wp) :: h(capacity,capacity), s(capacity,capacity), v(n), zero_v(n)
    real(wp) :: short_matrix(1,1), absolute_residual, relative_residual
    complex(wp) :: complex_h(capacity,capacity), complex_s(capacity,capacity)
    complex(wp) :: complex_v(n), complex_zero_v(n), complex_short(1,1)
    integer :: info

    h = 0.0_wp
    s = 0.0_wp
    h(1,1) = 2.0_wp
    h(2,2) = 4.0_wp
    s(1,1) = 1.0_wp
    s(2,2) = 1.0_wp
    v = [1.0_wp, -0.5_wp]
    zero_v = 0.0_wp
    short_matrix = 1.0_wp

    absolute_residual = -1.0_wp
    relative_residual = -1.0_wp
    call empty_state%factorization_residual(h, s, v, absolute_residual, &
                                            relative_residual, info)
    call check(info == QR_ERR_INVALID_STATE .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'real residual rejects invalid state with defined outputs', &
               failures)

    call state%initialize(capacity, info)
    call state%factorize_fresh(h, s, 0.5_wp, info, active_order=n)
    call state%factorization_residual(short_matrix, s, v, absolute_residual, &
                                      relative_residual, info)
    call check(info == QR_ERR_DIMENSION_MISMATCH .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'real residual rejects undersized matrix safely', failures)
    call state%factorization_residual(h, s, zero_v, absolute_residual, &
                                      relative_residual, info)
    call check(info == QR_ERR_ZERO_INITIAL_VECTOR .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'real residual rejects zero probe safely', failures)
    call state%factorize_fresh(h, s, 0.5_wp, info, active_order=0)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%is_valid() .and. &
               state%order() == n, &
               'real factorization rejects nonpositive active order', failures)
    call state%factorize_fresh(h, s, 0.5_wp, info, &
                               active_order=capacity + 1)
    call check(info == QR_ERR_CAPACITY_EXCEEDED .and. state%is_valid(), &
               'real factorization rejects active order above capacity', &
               failures)
    call state%factorize_fresh(short_matrix, short_matrix, 0.5_wp, info, &
                               active_order=n)
    call check(info == QR_ERR_DIMENSION_MISMATCH .and. state%is_valid(), &
               'real factorization rejects undersized active block', failures)
    call state%factorization_residual(h, s, v, absolute_residual, &
                                      relative_residual, info)
    call check(info == QR_SUCCESS .and. &
               relative_residual <= 10000.0_wp * epsilon(1.0_wp), &
               'rejected real factorizations preserve valid factors', failures)

    complex_h = cmplx(h, 0.0_wp, kind=wp)
    complex_s = cmplx(s, 0.0_wp, kind=wp)
    complex_v = cmplx(v, 0.25_wp * v, kind=wp)
    complex_zero_v = cmplx(0.0_wp, 0.0_wp, kind=wp)
    complex_short = cmplx(1.0_wp, 0.0_wp, kind=wp)
    absolute_residual = -1.0_wp
    relative_residual = -1.0_wp
    call complex_empty%factorization_residual(complex_h, complex_s, complex_v, &
      absolute_residual, relative_residual, info)
    call check(info == QR_ERR_INVALID_STATE .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'complex residual rejects invalid state with defined outputs', &
               failures)
    call complex_state%initialize(capacity, info)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.5_wp, info, &
                                       active_order=n)
    call complex_state%factorization_residual(complex_short, complex_s, &
      complex_v, absolute_residual, relative_residual, info)
    call check(info == QR_ERR_DIMENSION_MISMATCH .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'complex residual rejects undersized matrix safely', failures)
    call complex_state%factorization_residual(complex_h, complex_s, &
      complex_zero_v, absolute_residual, relative_residual, info)
    call check(info == QR_ERR_ZERO_INITIAL_VECTOR .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'complex residual rejects zero probe safely', failures)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.5_wp, info, &
                                       active_order=0)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               complex_state%is_valid() .and. complex_state%order() == n, &
               'complex factorization rejects nonpositive active order', &
               failures)
    call complex_state%factorize_fresh(complex_short, complex_short, 0.5_wp, &
                                       info, active_order=n)
    call check(info == QR_ERR_DIMENSION_MISMATCH .and. &
               complex_state%is_valid(), &
               'complex factorization rejects undersized active block', failures)
    call complex_state%factorization_residual(complex_h, complex_s, complex_v, &
      absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. &
               relative_residual <= 10000.0_wp * epsilon(1.0_wp), &
               'rejected complex factorizations preserve valid factors', &
               failures)

    !A nonzero probe against two exactly zero actions exercises the safe-zero
    !denominator directly. Rank deficiency does not prevent residual analysis.
    h = 0.0_wp
    s = 0.0_wp
    call zero_state%initialize(capacity, info)
    call zero_state%factorize_fresh(h, s, 0.5_wp, info, active_order=n)
    call zero_state%factorization_residual(h, s, v, absolute_residual, &
                                           relative_residual, info)
    call check(info == QR_SUCCESS .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'real zero actions return zero relative residual', failures)

    complex_h = cmplx(0.0_wp, 0.0_wp, kind=wp)
    complex_s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call complex_zero_state%initialize(capacity, info)
    call complex_zero_state%factorize_fresh(complex_h, complex_s, 0.5_wp, &
                                            info, active_order=n)
    call complex_zero_state%factorization_residual(complex_h, complex_s, &
      complex_v, absolute_residual, relative_residual, info)
    call check(info == QR_SUCCESS .and. &
               abs(absolute_residual) <= 0.0_wp .and. &
               abs(relative_residual) <= 0.0_wp, &
               'complex zero actions return zero relative residual', failures)
  end subroutine test_residual_validation

  subroutine apply_real_replacement(h, s, n, idx, delta_h, delta_s)
  !Apply one symmetric column change to full caller-owned reference matrices.
    real(wp), intent(inout) :: h(:,:), s(:,:)
    integer, intent(in) :: n, idx
    real(wp), intent(in) :: delta_h(:), delta_s(:)
    integer :: i

    do i = 1, n
      h(i,idx) = h(i,idx) + delta_h(i)
      s(i,idx) = s(i,idx) + delta_s(i)
      if (i /= idx) then
        h(idx,i) = h(idx,i) + delta_h(i)
        s(idx,i) = s(idx,i) + delta_s(i)
      end if
    end do
  end subroutine apply_real_replacement

  subroutine apply_complex_replacement(h, s, n, idx, delta_h, delta_s)
  !Apply one Hermitian column change and its conjugate row to the references.
    complex(wp), intent(inout) :: h(:,:), s(:,:)
    integer, intent(in) :: n, idx
    complex(wp), intent(in) :: delta_h(:), delta_s(:)
    integer :: i

    do i = 1, n
      h(i,idx) = h(i,idx) + delta_h(i)
      s(i,idx) = s(i,idx) + delta_s(i)
      if (i /= idx) then
        h(idx,i) = h(idx,i) + conjg(delta_h(i))
        s(idx,i) = s(idx,i) + conjg(delta_s(i))
      end if
    end do
  end subroutine apply_complex_replacement

  subroutine append_real_column(h, s, n, h_column, s_column)
  !Store an appended symmetric row and column in the active reference block.
    real(wp), intent(inout) :: h(:,:), s(:,:)
    integer, intent(in) :: n
    real(wp), intent(in) :: h_column(:), s_column(:)

    h(1:n,n) = h_column(1:n)
    h(n,1:n) = h_column(1:n)
    s(1:n,n) = s_column(1:n)
    s(n,1:n) = s_column(1:n)
  end subroutine append_real_column

  subroutine append_complex_column(h, s, n, h_column, s_column)
  !Store an appended Hermitian column and conjugate row in the active block.
    complex(wp), intent(inout) :: h(:,:), s(:,:)
    integer, intent(in) :: n
    complex(wp), intent(in) :: h_column(:), s_column(:)

    h(1:n,n) = h_column(1:n)
    h(n,1:n) = conjg(h_column(1:n))
    h(n,n) = cmplx(real(h_column(n), wp), 0.0_wp, kind=wp)
    s(1:n,n) = s_column(1:n)
    s(n,1:n) = conjg(s_column(1:n))
    s(n,n) = cmplx(real(s_column(n), wp), 0.0_wp, kind=wp)
  end subroutine append_complex_column

  subroutine remove_real_index(matrix, n, idx)
  !Pack the retained real principal submatrix into the leading active block.
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
  end subroutine remove_real_index

  subroutine remove_complex_index(matrix, n, idx)
  !Pack the retained complex principal submatrix into the leading active block.
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
  end subroutine remove_complex_index

end program test_factorization_residual

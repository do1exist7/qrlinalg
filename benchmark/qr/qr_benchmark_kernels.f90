module qr_benchmark_kernels
  use iso_fortran_env, only: int64, real64
  use wp_def, only: wp
  use qrlinalg, only: QR_SUCCESS, qr_complex_state, qr_real_state
  use benchmark_support
  implicit none
  private

  public :: run_qr_append
  public :: run_qr_deletion
  public :: run_qr_factorization
  public :: run_qr_full_solve
  public :: run_qr_replacement
  public :: run_qr_solve

contains

  subroutine run_qr_factorization()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call factorization_real(options)
    case ('complex')
      call factorization_complex(options)
    end select
  end subroutine run_qr_factorization

  subroutine run_qr_solve()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call solve_real(options, .false.)
    case ('complex')
      call solve_complex(options, .false.)
    end select
  end subroutine run_qr_solve

  subroutine run_qr_full_solve()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call solve_real(options, .true.)
    case ('complex')
      call solve_complex(options, .true.)
    end select
  end subroutine run_qr_full_solve

  subroutine run_qr_replacement()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call replacement_real(options)
    case ('complex')
      call replacement_complex(options)
    end select
  end subroutine run_qr_replacement

  subroutine run_qr_append()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call append_real(options)
    case ('complex')
      call append_complex(options)
    end select
  end subroutine run_qr_append

  subroutine run_qr_deletion()
    type(benchmark_options) :: options

    call read_benchmark_options(options)
    select case (trim(options%kind))
    case ('real')
      call deletion_real(options)
    case ('complex')
      call deletion_complex(options)
    end select
  end subroutine run_qr_deletion

  subroutine factorization_real(options)
    type(benchmark_options), intent(in) :: options
    type(qr_real_state) :: state
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('real', h_order)
    allocate(initial(h_order), x(h_order))
    initial = 1.0_wp
    call state%initialize(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR initialization failed')

    !The untimed call brings code and data pages into memory before sampling.
    call state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR warm-up factorization failed')
    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      call state%factorize_fresh(h, s, shift, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do

    call validate_real_state(state, h, s, initial, x, eigenvalue, &
                             relative_accuracy, iterations, info)
    call fill_result(result_data, 'factorization', 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine factorization_real

  subroutine factorization_complex(options)
    type(benchmark_options), intent(in) :: options
    type(qr_complex_state) :: state
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('complex', h_order)
    allocate(initial(h_order), x(h_order))
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call state%initialize(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR initialization failed')

    call state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR warm-up factorization failed')
    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      call state%factorize_fresh(h, s, shift, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do

    call validate_complex_state(state, h, s, initial, x, eigenvalue, &
                                relative_accuracy, iterations, info)
    call fill_result(result_data, 'factorization', 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_complex(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine factorization_complex

  subroutine solve_real(options, include_factorization)
    type(benchmark_options), intent(in) :: options
    logical, intent(in) :: include_factorization
    type(qr_real_state) :: state
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('real', h_order)
    allocate(initial(h_order), x(h_order))
    initial = 1.0_wp
    call state%initialize(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR initialization failed')

    call state%factorize_fresh(h, s, shift, info)
    if (info == QR_SUCCESS) then
      call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                       benchmark_max_iter, benchmark_norm_mode, &
                       relative_accuracy, iterations, info)
    end if
    if (info /= QR_SUCCESS) call fail_benchmark('QR solve warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      if (include_factorization) call state%factorize_fresh(h, s, shift, info)
      if (info == QR_SUCCESS) then
        call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                         benchmark_max_iter, benchmark_norm_mode, &
                         relative_accuracy, iterations, info)
      end if
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do

    if (include_factorization) then
      operation = 'full_solve'
    else
      operation = 'solve'
    end if
    call fill_result(result_data, operation, 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine solve_real

  subroutine solve_complex(options, include_factorization)
    type(benchmark_options), intent(in) :: options
    logical, intent(in) :: include_factorization
    type(qr_complex_state) :: state
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('complex', h_order)
    allocate(initial(h_order), x(h_order))
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call state%initialize(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR initialization failed')

    call state%factorize_fresh(h, s, shift, info)
    if (info == QR_SUCCESS) then
      call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                       benchmark_max_iter, benchmark_norm_mode, &
                       relative_accuracy, iterations, info)
    end if
    if (info /= QR_SUCCESS) call fail_benchmark('QR solve warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      if (include_factorization) call state%factorize_fresh(h, s, shift, info)
      if (info == QR_SUCCESS) then
        call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                         benchmark_max_iter, benchmark_norm_mode, &
                         relative_accuracy, iterations, info)
      end if
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do

    if (include_factorization) then
      operation = 'full_solve'
    else
      operation = 'solve'
    end if
    call fill_result(result_data, operation, 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_complex(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine solve_complex

  subroutine replacement_real(options)
    type(benchmark_options), intent(in) :: options
    type(qr_real_state) :: state
    type(benchmark_result) :: result_data
    real(wp), allocatable :: delta_h(:), delta_s(:), h(:,:), h_check(:,:)
    real(wp), allocatable :: initial(:), s(:,:), s_check(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, idx, info, iterations, repetition, s_order

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('real', h_order)
    idx = (h_order + 1) / 2
    allocate(delta_h(h_order), delta_s(h_order), h_check(h_order,h_order))
    allocate(initial(h_order), s_check(h_order,h_order), x(h_order))
    delta_h = 1.0e-6_wp * h(:,idx)
    delta_s = 1.0e-7_wp * s(:,idx)
    initial = 1.0_wp
    call state%initialize(h_order, info)
    if (info == QR_SUCCESS) call state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR replacement setup failed')

    !Warm both signs so the timed sequence starts from the original factors.
    call state%replace_symmetric(idx, delta_h, delta_s, info)
    if (info == QR_SUCCESS) call state%replace_symmetric(idx, -delta_h, -delta_s, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR replacement warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      call state%replace_symmetric(idx, delta_h, delta_s, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
      delta_h = -delta_h
      delta_s = -delta_s
    end do

    h_check = h
    s_check = s
    if (mod(options%repetitions, 2) == 1) then
      call apply_real_replacement(h_check, idx, -delta_h)
      call apply_real_replacement(s_check, idx, -delta_s)
    end if
    call validate_real_state(state, h_check, s_check, initial, x, eigenvalue, &
                             relative_accuracy, iterations, info)
    call fill_result(result_data, 'replacement', 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h_check, s_check, x, eigenvalue))
    call emit_result(result_data)
  end subroutine replacement_real

  subroutine replacement_complex(options)
    type(benchmark_options), intent(in) :: options
    type(qr_complex_state) :: state
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: delta_h(:), delta_s(:), h(:,:), h_check(:,:)
    complex(wp), allocatable :: initial(:), s(:,:), s_check(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, idx, info, iterations, repetition, s_order

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    shift = benchmark_shift('complex', h_order)
    idx = (h_order + 1) / 2
    allocate(delta_h(h_order), delta_s(h_order), h_check(h_order,h_order))
    allocate(initial(h_order), s_check(h_order,h_order), x(h_order))
    delta_h = 1.0e-6_wp * h(:,idx)
    delta_s = 1.0e-7_wp * s(:,idx)
    delta_h(idx) = cmplx(real(delta_h(idx), wp), 0.0_wp, kind=wp)
    delta_s(idx) = cmplx(real(delta_s(idx), wp), 0.0_wp, kind=wp)
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call state%initialize(h_order, info)
    if (info == QR_SUCCESS) call state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR replacement setup failed')

    call state%replace_symmetric(idx, delta_h, delta_s, info)
    if (info == QR_SUCCESS) call state%replace_symmetric(idx, -delta_h, -delta_s, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR replacement warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      call timer_start(clock_value)
      call state%replace_symmetric(idx, delta_h, delta_s, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
      delta_h = -delta_h
      delta_s = -delta_s
    end do

    h_check = h
    s_check = s
    if (mod(options%repetitions, 2) == 1) then
      call apply_complex_replacement(h_check, idx, -delta_h)
      call apply_complex_replacement(s_check, idx, -delta_s)
    end if
    call validate_complex_state(state, h_check, s_check, initial, x, &
                                eigenvalue, relative_accuracy, iterations, info)
    call fill_result(result_data, 'replacement', 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_complex(h_check, s_check, x, eigenvalue))
    call emit_result(result_data)
  end subroutine replacement_complex

  subroutine append_real(options)
    type(benchmark_options), intent(in) :: options
    type(qr_real_state) :: primed_state, state
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('append benchmark requires n >= 2')
    shift = benchmark_shift('real', h_order)
    allocate(initial(h_order), x(h_order))
    initial = 1.0_wp
    call primed_state%initialize(h_order, info)
    if (info == QR_SUCCESS) then
      call primed_state%factorize_fresh(h(1:h_order-1,1:h_order-1), &
                                        s(1:h_order-1,1:h_order-1), shift, info)
    end if
    if (info /= QR_SUCCESS) call fail_benchmark('QR append setup failed')

    state = primed_state
    call state%append_symmetric(h(:,h_order), s(:,h_order), info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR append warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      state = primed_state
      call timer_start(clock_value)
      call state%append_symmetric(h(:,h_order), s(:,h_order), info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do
    call validate_real_state(state, h, s, initial, x, eigenvalue, &
                             relative_accuracy, iterations, info)
    call fill_result(result_data, 'append', 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine append_real

  subroutine append_complex(options)
    type(benchmark_options), intent(in) :: options
    type(qr_complex_state) :: primed_state, state
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('append benchmark requires n >= 2')
    shift = benchmark_shift('complex', h_order)
    allocate(initial(h_order), x(h_order))
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call primed_state%initialize(h_order, info)
    if (info == QR_SUCCESS) then
      call primed_state%factorize_fresh(h(1:h_order-1,1:h_order-1), &
                                        s(1:h_order-1,1:h_order-1), shift, info)
    end if
    if (info /= QR_SUCCESS) call fail_benchmark('QR append setup failed')

    state = primed_state
    call state%append_symmetric(h(:,h_order), s(:,h_order), info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR append warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      state = primed_state
      call timer_start(clock_value)
      call state%append_symmetric(h(:,h_order), s(:,h_order), info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do
    call validate_complex_state(state, h, s, initial, x, eigenvalue, &
                                relative_accuracy, iterations, info)
    call fill_result(result_data, 'append', 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_complex(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine append_complex

  subroutine deletion_real(options)
    type(benchmark_options), intent(in) :: options
    type(qr_real_state) :: primed_state, state
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: final_n, h_order, info, iterations, repetition, s_order

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('deletion benchmark requires n >= 2')
    final_n = h_order - 1
    shift = benchmark_shift('real', h_order)
    allocate(initial(final_n), x(final_n))
    initial = 1.0_wp
    call primed_state%initialize(h_order, info)
    if (info == QR_SUCCESS) call primed_state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR deletion setup failed')

    state = primed_state
    call state%delete_symmetric(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR deletion warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      state = primed_state
      call timer_start(clock_value)
      call state%delete_symmetric(h_order, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do
    call validate_real_state(state, h(1:final_n,1:final_n), &
                             s(1:final_n,1:final_n), initial, x, eigenvalue, &
                             relative_accuracy, iterations, info)
    call fill_result(result_data, 'deletion', 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, normalized_residual_real( &
                     h(1:final_n,1:final_n), s(1:final_n,1:final_n), x, eigenvalue))
    call emit_result(result_data)
  end subroutine deletion_real

  subroutine deletion_complex(options)
    type(benchmark_options), intent(in) :: options
    type(qr_complex_state) :: primed_state, state
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), s(:,:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: final_n, h_order, info, iterations, repetition, s_order

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('deletion benchmark requires n >= 2')
    final_n = h_order - 1
    shift = benchmark_shift('complex', h_order)
    allocate(initial(final_n), x(final_n))
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call primed_state%initialize(h_order, info)
    if (info == QR_SUCCESS) call primed_state%factorize_fresh(h, s, shift, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR deletion setup failed')

    state = primed_state
    call state%delete_symmetric(h_order, info)
    if (info /= QR_SUCCESS) call fail_benchmark('QR deletion warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      state = primed_state
      call timer_start(clock_value)
      call state%delete_symmetric(h_order, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= QR_SUCCESS) exit
    end do
    call validate_complex_state(state, h(1:final_n,1:final_n), &
                                s(1:final_n,1:final_n), initial, x, eigenvalue, &
                                relative_accuracy, iterations, info)
    call fill_result(result_data, 'deletion', 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, normalized_residual_complex( &
                     h(1:final_n,1:final_n), s(1:final_n,1:final_n), x, eigenvalue))
    call emit_result(result_data)
  end subroutine deletion_complex

  subroutine validate_real_state(state, h, s, initial, x, eigenvalue, &
                                 relative_accuracy, iterations, info)
    type(qr_real_state), intent(inout) :: state
    real(wp), intent(in) :: h(:,:), s(:,:), initial(:)
    real(wp), intent(out) :: x(:), eigenvalue, relative_accuracy
    integer, intent(out) :: iterations, info

    call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                     benchmark_max_iter, benchmark_norm_mode, &
                     relative_accuracy, iterations, info)
  end subroutine validate_real_state

  subroutine validate_complex_state(state, h, s, initial, x, eigenvalue, &
                                    relative_accuracy, iterations, info)
    type(qr_complex_state), intent(inout) :: state
    complex(wp), intent(in) :: h(:,:), s(:,:), initial(:)
    complex(wp), intent(out) :: x(:)
    real(wp), intent(out) :: eigenvalue, relative_accuracy
    integer, intent(out) :: iterations, info

    call state%solve(s, initial, x, eigenvalue, benchmark_tolerance, &
                     benchmark_max_iter, benchmark_norm_mode, &
                     relative_accuracy, iterations, info)
  end subroutine validate_complex_state

  subroutine apply_real_replacement(matrix, idx, delta)
    real(wp), intent(inout) :: matrix(:,:)
    integer, intent(in) :: idx
    real(wp), intent(in) :: delta(:)
    integer :: i

    do i = 1, size(delta)
      matrix(i,idx) = matrix(i,idx) + delta(i)
      if (i /= idx) matrix(idx,i) = matrix(idx,i) + delta(i)
    end do
  end subroutine apply_real_replacement

  subroutine apply_complex_replacement(matrix, idx, delta)
    complex(wp), intent(inout) :: matrix(:,:)
    integer, intent(in) :: idx
    complex(wp), intent(in) :: delta(:)
    integer :: i

    do i = 1, size(delta)
      matrix(i,idx) = matrix(i,idx) + delta(i)
      if (i /= idx) matrix(idx,i) = matrix(idx,i) + conjg(delta(i))
    end do
  end subroutine apply_complex_replacement

  subroutine fill_result(result_data, operation, problem_kind, matrix_n, &
                         shift, repetitions, total_seconds, status, iterations, &
                         eigenvalue, relative_accuracy, residual)
    type(benchmark_result), intent(out) :: result_data
    character(len=*), intent(in) :: operation, problem_kind
    integer, intent(in) :: matrix_n, repetitions, status, iterations
    real(wp), intent(in) :: shift, eigenvalue, relative_accuracy, residual
    real(real64), intent(in) :: total_seconds

    result_data%implementation = 'qr'
    result_data%operation = operation
    result_data%kind = problem_kind
    result_data%matrix_n = matrix_n
    result_data%shift = shift
    result_data%repetitions = repetitions
    result_data%total_seconds = total_seconds
    result_data%status = status
    result_data%iterations = iterations
    result_data%eigenvalue = eigenvalue
    result_data%relative_accuracy = relative_accuracy
    result_data%residual = residual
  end subroutine fill_result

  subroutine require_matching_orders(h_order, s_order)
    integer, intent(in) :: h_order, s_order
    if (h_order /= s_order) call fail_benchmark('H and S orders differ')
  end subroutine require_matching_orders

end module qr_benchmark_kernels

module ldlt_benchmark_kernels
  use iso_fortran_env, only: int64, real64
  use wp_def, only: wp
  use linalg
  use benchmark_support
  implicit none
  private

  public :: run_ldlt_append_and_solve
  public :: run_ldlt_append_factorization
  public :: run_ldlt_factorization
  public :: run_ldlt_full_solve
  public :: run_ldlt_solve

contains

  subroutine run_ldlt_factorization()
    type(benchmark_options) :: options
    call read_benchmark_options(options)
    call configure_serial_original()
    select case (trim(options%kind))
    case ('real')
      call factorization_real(options)
    case ('complex')
      call factorization_complex(options)
    end select
  end subroutine run_ldlt_factorization

  subroutine run_ldlt_solve()
    type(benchmark_options) :: options
    call read_benchmark_options(options)
    call configure_serial_original()
    select case (trim(options%kind))
    case ('real')
      call solve_real(options, .false.)
    case ('complex')
      call solve_complex(options, .false.)
    end select
  end subroutine run_ldlt_solve

  subroutine run_ldlt_full_solve()
    type(benchmark_options) :: options
    call read_benchmark_options(options)
    call configure_serial_original()
    select case (trim(options%kind))
    case ('real')
      call solve_real(options, .true.)
    case ('complex')
      call solve_complex(options, .true.)
    end select
  end subroutine run_ldlt_full_solve

  subroutine run_ldlt_append_factorization()
    type(benchmark_options) :: options
    call read_benchmark_options(options)
    call configure_serial_original()
    select case (trim(options%kind))
    case ('real')
      call append_real(options, .false.)
    case ('complex')
      call append_complex(options, .false.)
    end select
  end subroutine run_ldlt_append_factorization

  subroutine run_ldlt_append_and_solve()
    type(benchmark_options) :: options
    call read_benchmark_options(options)
    call configure_serial_original()
    select case (trim(options%kind))
    case ('real')
      call append_real(options, .true.)
    case ('complex')
      call append_complex(options, .true.)
    end select
  end subroutine run_ldlt_append_and_solve

  subroutine configure_serial_original()
    Glob_NumOfProcs = 1
    Glob_ProcID = 0
    Glob_MPIErrCode = 0
  end subroutine configure_serial_original

  subroutine factorization_real(options)
    type(benchmark_options), intent(in) :: options
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), invd(:), matrix(:,:)
    real(wp), allocatable :: matrix_master(:,:), s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    call linalg_setparam(h_order)
    shift = benchmark_shift('real', h_order)
    allocate(matrix_master(h_order,h_order), matrix(h_order,h_order))
    allocate(invd(h_order), work(h_order), initial(h_order), x(h_order))
    matrix_master = h - shift * s
    initial = 1.0_wp

    matrix = matrix_master
    invd = 0.0_wp
    call LDLTF(1, h_order, matrix, h_order, invd, work, info)
    if (info /= 0) call fail_benchmark('LDLT warm-up factorization failed')
    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      matrix = matrix_master
      invd = 0.0_wp
      call timer_start(clock_value)
      call LDLTF(1, h_order, matrix, h_order, invd, work, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
    end do

    call validate_real_factor(matrix, invd, s, shift, initial, x, eigenvalue, &
                              relative_accuracy, iterations, info)
    call fill_result(result_data, 'factorization', 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine factorization_real

  subroutine factorization_complex(options)
    type(benchmark_options), intent(in) :: options
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), invd(:), matrix(:,:)
    complex(wp), allocatable :: matrix_master(:,:), s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    call linalg_setparam(h_order)
    shift = benchmark_shift('complex', h_order)
    allocate(matrix_master(h_order,h_order), matrix(h_order,h_order))
    allocate(invd(h_order), work(h_order), initial(h_order), x(h_order))
    matrix_master = h - cmplx(shift, 0.0_wp, kind=wp) * s
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)

    matrix = matrix_master
    invd = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call LDLHF(1, h_order, matrix, h_order, invd, work, info)
    if (info /= 0) call fail_benchmark('LDLH warm-up factorization failed')
    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      matrix = matrix_master
      invd = cmplx(0.0_wp, 0.0_wp, kind=wp)
      call timer_start(clock_value)
      call LDLHF(1, h_order, matrix, h_order, invd, work, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
    end do

    call validate_complex_factor(matrix, invd, s, shift, initial, x, eigenvalue, &
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
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), invd(:), matrix(:,:)
    real(wp), allocatable :: matrix_master(:,:), overlap(:,:), s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: factor_start, h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    call linalg_setparam(h_order)
    shift = benchmark_shift('real', h_order)
    allocate(matrix_master(h_order,h_order), matrix(h_order,h_order))
    allocate(overlap(h_order,h_order), invd(h_order), work(h_order))
    allocate(initial(h_order), x(h_order))
    matrix_master = h - shift * s
    matrix = matrix_master
    invd = 0.0_wp
    call LDLTF(1, h_order, matrix, h_order, invd, work, info)
    if (info /= 0) call fail_benchmark('LDLT solve setup failed')

    overlap = s
    initial = 1.0_wp
    call GSEPIIS(h_order + 1, h_order, matrix, h_order, invd, overlap, &
                 h_order, shift, initial, work, benchmark_tolerance, &
                 eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                 benchmark_norm_mode, iterations, info)
    if (info /= 0) call fail_benchmark('LDLT solve warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      if (include_factorization) then
        matrix = matrix_master
        invd = 0.0_wp
        factor_start = 1
      else
        factor_start = h_order + 1
      end if
      overlap = s
      initial = 1.0_wp
      call timer_start(clock_value)
      call GSEPIIS(factor_start, h_order, matrix, h_order, invd, overlap, &
                   h_order, shift, initial, work, benchmark_tolerance, &
                   eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                   benchmark_norm_mode, iterations, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
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
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), invd(:), matrix(:,:)
    complex(wp), allocatable :: matrix_master(:,:), overlap(:,:), s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: factor_start, h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    call linalg_setparam(h_order)
    shift = benchmark_shift('complex', h_order)
    allocate(matrix_master(h_order,h_order), matrix(h_order,h_order))
    allocate(overlap(h_order,h_order), invd(h_order), work(h_order))
    allocate(initial(h_order), x(h_order))
    matrix_master = h - cmplx(shift, 0.0_wp, kind=wp) * s
    matrix = matrix_master
    invd = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call LDLHF(1, h_order, matrix, h_order, invd, work, info)
    if (info /= 0) call fail_benchmark('LDLH solve setup failed')

    overlap = s
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call GHEPIIS(h_order + 1, h_order, matrix, h_order, invd, overlap, &
                 h_order, shift, initial, work, benchmark_tolerance, &
                 eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                 benchmark_norm_mode, iterations, info)
    if (info /= 0) call fail_benchmark('LDLH solve warm-up failed')

    total_seconds = 0.0_real64
    do repetition = 1, options%repetitions
      if (include_factorization) then
        matrix = matrix_master
        invd = cmplx(0.0_wp, 0.0_wp, kind=wp)
        factor_start = 1
      else
        factor_start = h_order + 1
      end if
      overlap = s
      initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
      call timer_start(clock_value)
      call GHEPIIS(factor_start, h_order, matrix, h_order, invd, overlap, &
                   h_order, shift, initial, work, benchmark_tolerance, &
                   eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                   benchmark_norm_mode, iterations, info)
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
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

  subroutine append_real(options, include_solve)
    type(benchmark_options), intent(in) :: options
    logical, intent(in) :: include_solve
    type(benchmark_result) :: result_data
    real(wp), allocatable :: h(:,:), initial(:), invd(:), invd_primed(:)
    real(wp), allocatable :: matrix(:,:), matrix_primed(:,:), overlap(:,:)
    real(wp), allocatable :: s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_real_matrix(options%h_file, h, h_order)
    call read_real_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('append benchmark requires n >= 2')
    call linalg_setparam(h_order)
    shift = benchmark_shift('real', h_order)
    allocate(matrix(h_order,h_order), matrix_primed(h_order,h_order))
    allocate(overlap(h_order,h_order), invd(h_order), invd_primed(h_order))
    allocate(work(h_order), initial(h_order), x(h_order))
    matrix_primed = h - shift * s
    invd_primed = 0.0_wp
    call LDLTF(1, h_order - 1, matrix_primed, h_order, invd_primed, work, info)
    if (info /= 0) call fail_benchmark('LDLT leading factorization failed')

    matrix = matrix_primed
    invd = invd_primed
    overlap = s
    initial = 1.0_wp
    if (include_solve) then
      call GSEPIIS(h_order, h_order, matrix, h_order, invd, overlap, h_order, &
                   shift, initial, work, benchmark_tolerance, eigenvalue, x, &
                   relative_accuracy, benchmark_max_iter, benchmark_norm_mode, &
                   iterations, info)
    else
      call LDLTF(h_order, h_order, matrix, h_order, invd, work, info)
    end if
    if (info /= 0) call fail_benchmark('LDLT append warm-up failed')

    total_seconds = 0.0_real64
    iterations = 0
    eigenvalue = 0.0_wp
    relative_accuracy = 0.0_wp
    do repetition = 1, options%repetitions
      matrix = matrix_primed
      invd = invd_primed
      overlap = s
      initial = 1.0_wp
      call timer_start(clock_value)
      if (include_solve) then
        call GSEPIIS(h_order, h_order, matrix, h_order, invd, overlap, h_order, &
                     shift, initial, work, benchmark_tolerance, eigenvalue, x, &
                     relative_accuracy, benchmark_max_iter, benchmark_norm_mode, &
                     iterations, info)
      else
        call LDLTF(h_order, h_order, matrix, h_order, invd, work, info)
      end if
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
    end do

    if (.not. include_solve .and. info == 0) then
      call validate_real_factor(matrix, invd, s, shift, initial, x, eigenvalue, &
                                relative_accuracy, iterations, info)
      operation = 'append_factorization'
    else
      operation = 'append_and_solve'
    end if
    call fill_result(result_data, operation, 'real', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_real(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine append_real

  subroutine append_complex(options, include_solve)
    type(benchmark_options), intent(in) :: options
    logical, intent(in) :: include_solve
    type(benchmark_result) :: result_data
    complex(wp), allocatable :: h(:,:), initial(:), invd(:), invd_primed(:)
    complex(wp), allocatable :: matrix(:,:), matrix_primed(:,:), overlap(:,:)
    complex(wp), allocatable :: s(:,:), work(:), x(:)
    real(wp) :: eigenvalue, relative_accuracy, shift
    real(real64) :: total_seconds
    integer(int64) :: clock_value
    integer :: h_order, info, iterations, repetition, s_order
    character(len=32) :: operation

    call read_complex_matrix(options%h_file, h, h_order)
    call read_complex_matrix(options%s_file, s, s_order)
    call require_matching_orders(h_order, s_order)
    if (h_order < 2) call fail_benchmark('append benchmark requires n >= 2')
    call linalg_setparam(h_order)
    shift = benchmark_shift('complex', h_order)
    allocate(matrix(h_order,h_order), matrix_primed(h_order,h_order))
    allocate(overlap(h_order,h_order), invd(h_order), invd_primed(h_order))
    allocate(work(h_order), initial(h_order), x(h_order))
    matrix_primed = h - cmplx(shift, 0.0_wp, kind=wp) * s
    invd_primed = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call LDLHF(1, h_order - 1, matrix_primed, h_order, invd_primed, work, info)
    if (info /= 0) call fail_benchmark('LDLH leading factorization failed')

    matrix = matrix_primed
    invd = invd_primed
    overlap = s
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    if (include_solve) then
      call GHEPIIS(h_order, h_order, matrix, h_order, invd, overlap, h_order, &
                   shift, initial, work, benchmark_tolerance, eigenvalue, x, &
                   relative_accuracy, benchmark_max_iter, benchmark_norm_mode, &
                   iterations, info)
    else
      call LDLHF(h_order, h_order, matrix, h_order, invd, work, info)
    end if
    if (info /= 0) call fail_benchmark('LDLH append warm-up failed')

    total_seconds = 0.0_real64
    iterations = 0
    eigenvalue = 0.0_wp
    relative_accuracy = 0.0_wp
    do repetition = 1, options%repetitions
      matrix = matrix_primed
      invd = invd_primed
      overlap = s
      initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
      call timer_start(clock_value)
      if (include_solve) then
        call GHEPIIS(h_order, h_order, matrix, h_order, invd, overlap, h_order, &
                     shift, initial, work, benchmark_tolerance, eigenvalue, x, &
                     relative_accuracy, benchmark_max_iter, benchmark_norm_mode, &
                     iterations, info)
      else
        call LDLHF(h_order, h_order, matrix, h_order, invd, work, info)
      end if
      total_seconds = total_seconds + timer_elapsed(clock_value)
      if (info /= 0) exit
    end do

    if (.not. include_solve .and. info == 0) then
      call validate_complex_factor(matrix, invd, s, shift, initial, x, eigenvalue, &
                                   relative_accuracy, iterations, info)
      operation = 'append_factorization'
    else
      operation = 'append_and_solve'
    end if
    call fill_result(result_data, operation, 'complex', h_order, shift, &
                     options%repetitions, total_seconds, info, iterations, &
                     eigenvalue, relative_accuracy, &
                     normalized_residual_complex(h, s, x, eigenvalue))
    call emit_result(result_data)
  end subroutine append_complex

  subroutine validate_real_factor(matrix, invd, s, shift, initial, x, &
                                  eigenvalue, relative_accuracy, iterations, info)
    real(wp), intent(inout) :: matrix(:,:), invd(:)
    real(wp), intent(in) :: s(:,:), shift
    real(wp), intent(out) :: initial(:), x(:), eigenvalue, relative_accuracy
    integer, intent(out) :: iterations, info
    real(wp), allocatable :: overlap(:,:), work(:)
    integer :: matrix_n

    matrix_n = size(invd)
    allocate(overlap(matrix_n,matrix_n), work(matrix_n))
    overlap = s
    initial = 1.0_wp
    call GSEPIIS(matrix_n + 1, matrix_n, matrix, matrix_n, invd, overlap, &
                 matrix_n, shift, initial, work, benchmark_tolerance, &
                 eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                 benchmark_norm_mode, iterations, info)
  end subroutine validate_real_factor

  subroutine validate_complex_factor(matrix, invd, s, shift, initial, x, &
                                     eigenvalue, relative_accuracy, iterations, info)
    complex(wp), intent(inout) :: matrix(:,:), invd(:)
    complex(wp), intent(in) :: s(:,:)
    real(wp), intent(in) :: shift
    complex(wp), intent(out) :: initial(:), x(:)
    real(wp), intent(out) :: eigenvalue, relative_accuracy
    integer, intent(out) :: iterations, info
    complex(wp), allocatable :: overlap(:,:), work(:)
    integer :: matrix_n

    matrix_n = size(invd)
    allocate(overlap(matrix_n,matrix_n), work(matrix_n))
    overlap = s
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call GHEPIIS(matrix_n + 1, matrix_n, matrix, matrix_n, invd, overlap, &
                 matrix_n, shift, initial, work, benchmark_tolerance, &
                 eigenvalue, x, relative_accuracy, benchmark_max_iter, &
                 benchmark_norm_mode, iterations, info)
  end subroutine validate_complex_factor

  subroutine fill_result(result_data, operation, problem_kind, matrix_n, &
                         shift, repetitions, total_seconds, status, iterations, &
                         eigenvalue, relative_accuracy, residual)
    type(benchmark_result), intent(out) :: result_data
    character(len=*), intent(in) :: operation, problem_kind
    integer, intent(in) :: matrix_n, repetitions, status, iterations
    real(wp), intent(in) :: shift, eigenvalue, relative_accuracy, residual
    real(real64), intent(in) :: total_seconds

    result_data%implementation = 'ldlt'
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

end module ldlt_benchmark_kernels

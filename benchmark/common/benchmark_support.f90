module benchmark_support
  use iso_fortran_env, only: error_unit, int64, real64
  use wp_def, only: wp
  implicit none
  private

  real(wp), parameter, public :: benchmark_tolerance = 1.0e-12_wp
  integer, parameter, public :: benchmark_max_iter = 30
  integer, parameter, public :: benchmark_norm_mode = 1
  integer, save :: perf_control_unit = -1
  logical, save :: perf_control_checked = .false.
  logical, save :: perf_control_enabled = .false.
  character(len=1024), save :: perf_ack_path = ''

  type, public :: benchmark_options
    character(len=16) :: kind = ''
    character(len=1024) :: h_file = ''
    character(len=1024) :: s_file = ''
    integer :: repetitions = 0
  end type benchmark_options

  type, public :: benchmark_result
    character(len=16) :: implementation = ''
    character(len=32) :: operation = ''
    character(len=16) :: kind = ''
    integer :: matrix_n = 0
    integer :: repetitions = 0
    real(wp) :: shift = 0.0_wp
    real(real64) :: total_seconds = 0.0_real64
    integer :: status = 0
    integer :: iterations = 0
    real(wp) :: eigenvalue = 0.0_wp
    real(wp) :: relative_accuracy = 0.0_wp
    real(wp) :: residual = 0.0_wp
  end type benchmark_result

  public :: benchmark_shift
  public :: emit_result
  public :: fail_benchmark
  public :: normalized_residual_complex
  public :: normalized_residual_real
  public :: read_benchmark_options
  public :: read_complex_matrix
  public :: read_real_matrix
  public :: timer_elapsed
  public :: timer_start

contains

  subroutine read_benchmark_options(options)
    type(benchmark_options), intent(out) :: options
    character(len=64) :: argument
    integer :: io_status

    if (command_argument_count() /= 4) then
      call fail_benchmark('expected KIND H_FILE S_FILE REPETITIONS')
    end if

    call get_command_argument(1, options%kind)
    call get_command_argument(2, options%h_file)
    call get_command_argument(3, options%s_file)
    call get_command_argument(4, argument)
    read(argument, *, iostat=io_status) options%repetitions
    if (io_status /= 0 .or. options%repetitions <= 0) then
      call fail_benchmark('REPETITIONS must be a positive integer')
    end if
    if (trim(options%kind) /= 'real' .and. &
        trim(options%kind) /= 'complex') then
      call fail_benchmark('KIND must be real or complex')
    end if
  end subroutine read_benchmark_options

  real(wp) function benchmark_shift(problem_kind, matrix_n) result(shift)
    character(len=*), intent(in) :: problem_kind
    integer, intent(in) :: matrix_n

    select case (trim(problem_kind))
    case ('real')
      shift = -7.4_wp
      if (matrix_n >= 100) shift = -7.33473_wp
      if (matrix_n >= 200) shift = -7.33491_wp
      if (matrix_n >= 500) shift = -7.33493_wp
      if (matrix_n >= 1000) shift = -7.33493_wp
      if (matrix_n >= 2000) shift = -7.33493_wp
      if (matrix_n >= 4000) shift = -7.334928_wp
    case ('complex')
      shift = -7.8_wp
      if (matrix_n >= 100) shift = -7.77461_wp
      if (matrix_n >= 200) shift = -7.780087_wp
      if (matrix_n >= 500) shift = -7.782207_wp
      if (matrix_n >= 1000) shift = -7.7830033_wp
      if (matrix_n >= 2000) shift = -7.7831307_wp
      if (matrix_n >= 4000) shift = -7.7831672_wp
    case default
      call fail_benchmark('cannot select a shift for unknown KIND')
    end select
  end function benchmark_shift

  subroutine read_real_matrix(file_name, matrix, matrix_n)
    character(len=*), intent(in) :: file_name
    real(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_n
    real(real64), allocatable :: stored_column(:)
    integer :: column_index, io_status, row_index, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) call fail_benchmark('cannot open ' // trim(file_name))
    read(unit, iostat=io_status) matrix_n
    if (io_status /= 0 .or. matrix_n <= 0) then
      call fail_benchmark('invalid matrix order in ' // trim(file_name))
    end if

    allocate(matrix(matrix_n,matrix_n), stored_column(matrix_n))
    do column_index = 1, matrix_n
      read(unit, iostat=io_status) stored_column(1:column_index)
      if (io_status /= 0) call fail_benchmark('truncated ' // trim(file_name))
      matrix(1:column_index,column_index) = &
        real(stored_column(1:column_index), wp)
    end do
    close(unit)

    do column_index = 1, matrix_n
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = matrix(row_index,column_index)
      end do
    end do
  end subroutine read_real_matrix

  subroutine read_complex_matrix(file_name, matrix, matrix_n)
    character(len=*), intent(in) :: file_name
    complex(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_n
    complex(real64), allocatable :: stored_column(:)
    integer :: column_index, io_status, row_index, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) call fail_benchmark('cannot open ' // trim(file_name))
    read(unit, iostat=io_status) matrix_n
    if (io_status /= 0 .or. matrix_n <= 0) then
      call fail_benchmark('invalid matrix order in ' // trim(file_name))
    end if

    allocate(matrix(matrix_n,matrix_n), stored_column(matrix_n))
    do column_index = 1, matrix_n
      read(unit, iostat=io_status) stored_column(1:column_index)
      if (io_status /= 0) call fail_benchmark('truncated ' // trim(file_name))
      do row_index = 1, column_index
        matrix(row_index,column_index) = &
          cmplx(real(stored_column(row_index), real64), &
                aimag(stored_column(row_index)), kind=wp)
      end do
    end do
    close(unit)

    do column_index = 1, matrix_n
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = conjg(matrix(row_index,column_index))
      end do
      matrix(column_index,column_index) = &
        cmplx(real(matrix(column_index,column_index), wp), 0.0_wp, kind=wp)
    end do
  end subroutine read_complex_matrix

  subroutine timer_start(clock_value)
    integer(int64), intent(out) :: clock_value

    call set_perf_state('enable')
    call system_clock(clock_value)
  end subroutine timer_start

  real(real64) function timer_elapsed(clock_value) result(seconds)
    integer(int64), intent(in) :: clock_value
    integer(int64) :: clock_end, clock_rate

    call system_clock(clock_end, clock_rate)
    call set_perf_state('disable')
    seconds = real(clock_end - clock_value, real64) / real(clock_rate, real64)
  end function timer_elapsed

  subroutine set_perf_state(command)
    character(len=*), intent(in) :: command
    character(len=32) :: acknowledgement
    character(len=1024) :: acknowledgement_path, control_path
    integer :: environment_status, io_status, path_length, perf_ack_unit

    if (.not. perf_control_checked) then
      perf_control_checked = .true.
      call get_environment_variable('QRLINALG_PERF_CONTROL', control_path, &
                                    length=path_length, status=environment_status)
      if (environment_status == 0 .and. path_length > 0) then
        control_path = control_path(1:path_length)
        call get_environment_variable('QRLINALG_PERF_ACK', acknowledgement_path, &
                                      length=path_length, status=environment_status)
        if (environment_status /= 0 .or. path_length <= 0) then
          call fail_benchmark('QRLINALG_PERF_ACK is missing')
        end if
        acknowledgement_path = acknowledgement_path(1:path_length)
        perf_ack_path = trim(acknowledgement_path)
        open(newunit=perf_control_unit, file=trim(control_path), status='old', &
             action='write', form='formatted', iostat=io_status)
        if (io_status /= 0) call fail_benchmark('cannot open perf control FIFO')
        perf_control_enabled = .true.
      end if
    end if

    if (.not. perf_control_enabled) return
    write(perf_control_unit,'(a)',iostat=io_status) trim(command)
    if (io_status /= 0) call fail_benchmark('cannot send perf control command')
    flush(perf_control_unit)
    open(newunit=perf_ack_unit, file=trim(perf_ack_path), status='old', &
         action='read', form='formatted', iostat=io_status)
    if (io_status /= 0) call fail_benchmark('cannot open perf acknowledgement FIFO')
    read(perf_ack_unit,'(a)',iostat=io_status) acknowledgement
    close(perf_ack_unit)
    if (io_status /= 0 .or. trim(acknowledgement) /= 'ack') then
      call fail_benchmark('perf did not acknowledge ' // trim(command))
    end if
  end subroutine set_perf_state

  real(wp) function normalized_residual_real(h, s, eigenvector, eigenvalue) &
      result(residual)
    real(wp), intent(in) :: h(:,:), s(:,:), eigenvector(:), eigenvalue
    real(wp), allocatable :: residual_vector(:)
    real(wp) :: denominator

    allocate(residual_vector(size(eigenvector)))
    residual_vector = matmul(h, eigenvector) - &
                      eigenvalue * matmul(s, eigenvector)
    denominator = (frobenius_norm_real(h) + abs(eigenvalue) * &
                   frobenius_norm_real(s)) * vector_norm_real(eigenvector)
    if (denominator > tiny(1.0_wp)) then
      residual = vector_norm_real(residual_vector) / denominator
    else
      residual = huge(1.0_wp)
    end if
  end function normalized_residual_real

  real(wp) function normalized_residual_complex(h, s, eigenvector, eigenvalue) &
      result(residual)
    complex(wp), intent(in) :: h(:,:), s(:,:), eigenvector(:)
    real(wp), intent(in) :: eigenvalue
    complex(wp), allocatable :: residual_vector(:)
    real(wp) :: denominator

    allocate(residual_vector(size(eigenvector)))
    residual_vector = matmul(h, eigenvector) - &
      cmplx(eigenvalue, 0.0_wp, kind=wp) * matmul(s, eigenvector)
    denominator = (frobenius_norm_complex(h) + abs(eigenvalue) * &
                   frobenius_norm_complex(s)) * &
                  vector_norm_complex(eigenvector)
    if (denominator > tiny(1.0_wp)) then
      residual = vector_norm_complex(residual_vector) / denominator
    else
      residual = huge(1.0_wp)
    end if
  end function normalized_residual_complex

  subroutine emit_result(result_data)
    type(benchmark_result), intent(in) :: result_data
    real(real64) :: milliseconds_per_operation, seconds_per_operation

    seconds_per_operation = result_data%total_seconds / &
                            real(result_data%repetitions, real64)
    milliseconds_per_operation = 1000.0_real64 * seconds_per_operation

    write(*,'(a)') '------------------------------------------------------------'
    write(*,'(a,a)')       'Implementation          : ', trim(result_data%implementation)
    write(*,'(a,a)')       'Operation               : ', trim(result_data%operation)
    write(*,'(a,a)')       'Problem kind            : ', trim(result_data%kind)
    write(*,'(a,i0)')      'Working precision       : ', wp
    write(*,'(a,i0)')      'Matrix order            : ', result_data%matrix_n
    write(*,'(a,es24.16)') 'Initial eigenvalue shift: ', result_data%shift
    write(*,'(a,i0)')      'Timed operations        : ', result_data%repetitions
    write(*,'(a,i0)')      'Status                  : ', result_data%status
    write(*,'(a,i0)')      'Inverse iterations      : ', result_data%iterations
    write(*,'(a,es24.16)') 'Converged eigenvalue    : ', result_data%eigenvalue
    write(*,'(a,es24.16)') 'Relative accuracy       : ', result_data%relative_accuracy
    write(*,'(a,es24.16)') 'Normalized residual     : ', result_data%residual
    write(*,'(a,f16.6,a)') 'Total timed duration    : ', &
                            result_data%total_seconds, ' s'
    write(*,'(a,f16.6,a)') 'Mean operation duration : ', &
                            milliseconds_per_operation, ' ms'
    write(*,'(a)') '------------------------------------------------------------'
    write(*,'(a,",",a,",",a,",",a,",",i0,",",i0,",",es24.16,",",i0,",",' // &
             'es24.16,",",es24.16,",",i0,",",i0,",",es24.16,",",es24.16,' // &
             '",",es24.16)') &
      'BENCH_RESULT', trim(result_data%implementation), &
      trim(result_data%operation), trim(result_data%kind), wp, &
      result_data%matrix_n, result_data%shift, result_data%repetitions, &
      result_data%total_seconds, seconds_per_operation, result_data%status, &
      result_data%iterations, result_data%eigenvalue, &
      result_data%relative_accuracy, result_data%residual
  end subroutine emit_result

  subroutine fail_benchmark(message)
    character(len=*), intent(in) :: message

    write(error_unit,'(a)') 'benchmark error: ' // trim(message)
    error stop 1
  end subroutine fail_benchmark

  real(wp) function vector_norm_real(vector) result(norm_value)
    real(wp), intent(in) :: vector(:)
    norm_value = sqrt(sum(vector * vector))
  end function vector_norm_real

  real(wp) function vector_norm_complex(vector) result(norm_value)
    complex(wp), intent(in) :: vector(:)
    norm_value = sqrt(sum(real(conjg(vector) * vector, wp)))
  end function vector_norm_complex

  real(wp) function frobenius_norm_real(matrix) result(norm_value)
    real(wp), intent(in) :: matrix(:,:)
    norm_value = sqrt(sum(matrix * matrix))
  end function frobenius_norm_real

  real(wp) function frobenius_norm_complex(matrix) result(norm_value)
    complex(wp), intent(in) :: matrix(:,:)
    norm_value = sqrt(sum(real(conjg(matrix) * matrix, wp)))
  end function frobenius_norm_complex

end module benchmark_support

! Separate timing driver for real DGEQRF and DORGQR.
!
! Usage:
!   qr_stage_benchmark N REPETITIONS_PER_SAMPLE SAMPLES
!
! Matrix restoration is outside each timed call.  The output reports median
! and minimum seconds per LAPACK invocation for the same compact factors.
program qr_stage_benchmark
  use iso_fortran_env, only: int64, real64
  use wp_def, only: wp
  implicit none

  interface
    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      import :: wp
      integer, intent(in) :: m, n, lda, lwork
      real(wp), intent(inout) :: a(lda,*), work(*)
      real(wp), intent(out) :: tau(*)
      integer, intent(out) :: info
    end subroutine dgeqrf

    subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
      import :: wp
      integer, intent(in) :: m, n, k, lda, lwork
      real(wp), intent(inout) :: a(lda,*), work(*)
      real(wp), intent(in) :: tau(*)
      integer, intent(out) :: info
    end subroutine dorgqr
  end interface

  real(wp), allocatable :: a(:,:), compact_factors(:,:), original(:,:)
  real(wp), allocatable :: tau(:), work(:)
  real(wp) :: query_work(1)
  real(real64), allocatable :: factor_seconds(:), generate_seconds(:)
  real(real64) :: elapsed_sum
  integer(int64) :: clock_begin, clock_end, clock_rate
  integer :: info, lwork, matrix_n, repetition, repetitions_per_sample
  integer :: sample_count, sample_index

  call read_integer_argument(1, matrix_n)
  call read_integer_argument(2, repetitions_per_sample)
  call read_integer_argument(3, sample_count)
  if (matrix_n < 1 .or. repetitions_per_sample < 1 .or. sample_count < 1) &
    error stop 'invalid QR stage benchmark arguments'

  allocate(a(matrix_n,matrix_n), compact_factors(matrix_n,matrix_n))
  allocate(original(matrix_n,matrix_n), tau(matrix_n))
  allocate(factor_seconds(sample_count), generate_seconds(sample_count))
  call fill_matrix(original)

  a = original
  call dgeqrf(matrix_n, matrix_n, a, matrix_n, tau, query_work, -1, info)
  if (info /= 0) error stop 'DGEQRF workspace query failed'
  lwork = max(1, int(query_work(1)))
  call dorgqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, query_work, &
              -1, info)
  if (info /= 0) error stop 'DORGQR workspace query failed'
  lwork = max(lwork, int(query_work(1)))
  allocate(work(lwork))

  a = original
  call dgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work, lwork, info)
  if (info /= 0) error stop 'DGEQRF preparation failed'
  compact_factors = a
  a = compact_factors
  call dorgqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, work, lwork, info)
  if (info /= 0) error stop 'DORGQR warm-up failed'

  call system_clock(count_rate=clock_rate)
  do sample_index = 1, sample_count
    elapsed_sum = 0.0_real64
    do repetition = 1, repetitions_per_sample
      a = original
      call system_clock(clock_begin)
      call dgeqrf(matrix_n, matrix_n, a, matrix_n, tau, work, lwork, info)
      call system_clock(clock_end)
      if (info /= 0) error stop 'timed DGEQRF failed'
      elapsed_sum = elapsed_sum + real(clock_end-clock_begin, real64) / &
                                    real(clock_rate, real64)
    end do
    factor_seconds(sample_index) = elapsed_sum / repetitions_per_sample

    elapsed_sum = 0.0_real64
    do repetition = 1, repetitions_per_sample
      a = compact_factors
      call system_clock(clock_begin)
      call dorgqr(matrix_n, matrix_n, matrix_n, a, matrix_n, tau, work, &
                  lwork, info)
      call system_clock(clock_end)
      if (info /= 0) error stop 'timed DORGQR failed'
      elapsed_sum = elapsed_sum + real(clock_end-clock_begin, real64) / &
                                    real(clock_rate, real64)
    end do
    generate_seconds(sample_index) = elapsed_sum / repetitions_per_sample
  end do
  call sort_values(factor_seconds)
  call sort_values(generate_seconds)

  write(*,'(a,i0,2(",",i0),2(",",es24.16))') &
    'QR_STAGE,DGEQRF,', matrix_n, repetitions_per_sample, sample_count, &
    factor_seconds((sample_count+1)/2), factor_seconds(1)
  write(*,'(a,i0,2(",",i0),2(",",es24.16))') &
    'QR_STAGE,DORGQR,', matrix_n, repetitions_per_sample, sample_count, &
    generate_seconds((sample_count+1)/2), generate_seconds(1)

contains

  subroutine fill_matrix(matrix)
    real(wp), intent(out) :: matrix(:,:)
    integer :: row, column

    do column = 1, size(matrix,2)
      do row = 1, size(matrix,1)
        matrix(row,column) = &
          real(mod(17*row + 31*column, 101) - 50, wp) / 53.0_wp
      end do
      matrix(column,column) = matrix(column,column) + real(size(matrix,1), wp)
    end do
  end subroutine fill_matrix

  subroutine sort_values(values)
    real(real64), intent(inout) :: values(:)
    real(real64) :: value
    integer :: first, insertion

    do first = 2, size(values)
      value = values(first)
      insertion = first - 1
      do while (insertion >= 1)
        if (values(insertion) <= value) exit
        values(insertion+1) = values(insertion)
        insertion = insertion - 1
      end do
      values(insertion+1) = value
    end do
  end subroutine sort_values

  subroutine read_integer_argument(position, value)
    integer, intent(in) :: position
    integer, intent(out) :: value
    character(len=64) :: argument
    integer :: io_status

    call get_command_argument(position, argument)
    read(argument,*,iostat=io_status) value
    if (io_status /= 0) error stop 'invalid integer benchmark argument'
  end subroutine read_integer_argument

end program qr_stage_benchmark

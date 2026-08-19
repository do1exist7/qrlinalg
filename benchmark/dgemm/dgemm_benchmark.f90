! Isolated benchmark for the bundled precision-generic DGEMM implementation.
!
! Usage:
!   dgemm_benchmark TRANSA TRANSB M N K ALPHA BETA CALLS_PER_SAMPLE SAMPLES
!
! Each reported sample contains enough repeated calls to amortize timer noise.
! The final record reports the median duration per call and a direct-reference
! error from an independently initialized single product.
program dgemm_benchmark
  use iso_fortran_env, only: int64, real64
  use wp_def, only: wp
  implicit none

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc)
      import :: wp
      character, intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      real(wp), intent(in) :: alpha, beta
      real(wp), intent(in) :: a(lda,*), b(ldb,*)
      real(wp), intent(inout) :: c(ldc,*)
    end subroutine dgemm
  end interface

  real(wp), allocatable :: a(:,:), b(:,:), c(:,:), c_initial(:,:), expected(:,:)
  real(wp) :: alpha, beta, error, product_sum
  real(real64), allocatable :: sample_seconds(:)
  integer(int64) :: clock_begin, clock_end, clock_rate
  integer :: a_columns, a_rows, b_columns, b_rows, calls_per_sample
  integer :: i, j, k, l, lda, ldb, ldc, m, n, repetition, sample_count
  integer :: sample_index
  character :: transa, transb

  call read_character_argument(1, transa)
  call read_character_argument(2, transb)
  call read_integer_argument(3, m)
  call read_integer_argument(4, n)
  call read_integer_argument(5, k)
  call read_real_argument(6, alpha)
  call read_real_argument(7, beta)
  call read_integer_argument(8, calls_per_sample)
  call read_integer_argument(9, sample_count)
  if (m < 0 .or. n < 0 .or. k < 0 .or. calls_per_sample < 1 .or. &
      sample_count < 1) error stop 'invalid benchmark dimensions or repetitions'
  if (index('NTC', transa) == 0 .or. index('NTC', transb) == 0) &
    error stop 'transpose arguments must be N, T, or C'

  if (transa == 'N') then
    a_rows = m
    a_columns = k
  else
    a_rows = k
    a_columns = m
  end if
  if (transb == 'N') then
    b_rows = k
    b_columns = n
  else
    b_rows = n
    b_columns = k
  end if
  lda = max(1, a_rows) + 3
  ldb = max(1, b_rows) + 5
  ldc = max(1, m) + 7
  allocate(a(lda,max(1,a_columns)), b(ldb,max(1,b_columns)))
  allocate(c(ldc,max(1,n)), c_initial(ldc,max(1,n)))
  allocate(expected(ldc,max(1,n)), sample_seconds(sample_count))
  call fill_matrix(a, 3)
  call fill_matrix(b, 17)
  call fill_matrix(c_initial, 29)

  ! Validate one call before timing.  beta=0 is handled without referencing the
  ! corresponding initial C element, matching the BLAS ownership contract.
  c = c_initial
  expected = c_initial
  do j = 1, n
    do i = 1, m
      product_sum = 0.0_wp
      do l = 1, k
        product_sum = product_sum + matrix_element(a, transa, i, l) * &
                                      matrix_element(b, transb, l, j)
      end do
      if (abs(beta) <= 0.0_wp) then
        expected(i,j) = alpha * product_sum
      else
        expected(i,j) = alpha * product_sum + beta * c_initial(i,j)
      end if
    end do
  end do
  call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  if (m > 0 .and. n > 0) then
    error = maxval(abs(c(1:m,1:n) - expected(1:m,1:n)))
  else
    error = 0.0_wp
  end if

  c = c_initial
  call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call system_clock(count_rate=clock_rate)
  do sample_index = 1, sample_count
    call system_clock(clock_begin)
    do repetition = 1, calls_per_sample
      call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    end do
    call system_clock(clock_end)
    sample_seconds(sample_index) = real(clock_end - clock_begin, real64) / &
                                   real(clock_rate * calls_per_sample, real64)
  end do
  call sort_values(sample_seconds)

  write(*,'(a,a,a,a,3(",",i0),2(",",es24.16),2(",",i0),3(",",es24.16))') &
    'DGEMM_BENCH,', transa, ',', transb, m, n, k, alpha, beta, &
    calls_per_sample, sample_count, &
    real(sample_seconds((sample_count+1)/2), wp), real(sample_seconds(1), wp), &
    error

contains

  subroutine fill_matrix(matrix, seed)
    real(wp), intent(out) :: matrix(:,:)
    integer, intent(in) :: seed
    integer :: row, column

    do column = 1, size(matrix,2)
      do row = 1, size(matrix,1)
        matrix(row,column) = &
          real(mod(seed + 17*row + 31*column, 101) - 50, wp) / 53.0_wp
      end do
    end do
  end subroutine fill_matrix

  function matrix_element(matrix, transpose_mode, row, column) result(value)
    real(wp), intent(in) :: matrix(:,:)
    character, intent(in) :: transpose_mode
    integer, intent(in) :: row, column
    real(wp) :: value

    if (transpose_mode == 'N') then
      value = matrix(row,column)
    else
      value = matrix(column,row)
    end if
  end function matrix_element

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

  subroutine read_character_argument(position, value)
    integer, intent(in) :: position
    character, intent(out) :: value
    character(len=32) :: argument

    call get_command_argument(position, argument)
    if (len_trim(argument) < 1) error stop 'missing benchmark argument'
    value = argument(1:1)
  end subroutine read_character_argument

  subroutine read_integer_argument(position, value)
    integer, intent(in) :: position
    integer, intent(out) :: value
    character(len=64) :: argument
    integer :: io_status

    call get_command_argument(position, argument)
    read(argument,*,iostat=io_status) value
    if (io_status /= 0) error stop 'invalid integer benchmark argument'
  end subroutine read_integer_argument

  subroutine read_real_argument(position, value)
    integer, intent(in) :: position
    real(wp), intent(out) :: value
    character(len=64) :: argument
    integer :: io_status

    call get_command_argument(position, argument)
    read(argument,*,iostat=io_status) value
    if (io_status /= 0) error stop 'invalid real benchmark argument'
  end subroutine read_real_argument

end program dgemm_benchmark

! Semantic tests for the bundled precision-generic real matrix product.
!
! The test reference is an intentionally direct element-by-element product.
! It exercises every real transpose spelling, padded leading dimensions,
! rectangular shapes, scalar special cases, and both optimized and fallback
! paths without depending on the compiler's MATMUL implementation.
program test_dgemm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
                                           ieee_value
  use qrlinalg, only: wp
  use test_support, only: check, finish_suite
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

  character, parameter :: transposes(3) = ['N', 'T', 'C']
  real(wp), parameter :: scalars(4) = [0.0_wp, 1.0_wp, -1.0_wp, 0.375_wp]
  integer :: alpha_index, beta_index, failures, k, m, n
  integer :: transa_index, transb_index

  failures = 0

  ! Orders zero through three cover quick returns, scalar products, odd tails,
  ! and dimensions smaller than both optimized micro-kernel extents.
  do transa_index = 1, size(transposes)
    do transb_index = 1, size(transposes)
      do m = 0, 3
        do n = 0, 3
          do k = 0, 3
            do alpha_index = 1, size(scalars)
              do beta_index = 1, size(scalars)
                call run_case(transposes(transa_index), &
                              transposes(transb_index), m, n, k, &
                              scalars(alpha_index), scalars(beta_index), 2, &
                              failures)
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  ! These larger rectangular cases enter the four-column TN and NT paths and
  ! the NN and TT fallbacks with nonminimal leading dimensions.
  call run_case('T', 'N', 7, 8, 9, 1.0_wp, 1.0_wp, 3, failures)
  call run_case('T', 'N', 8, 7, 5, -1.0_wp, 0.375_wp, 4, failures)
  call run_case('N', 'T', 11, 8, 7, -1.0_wp, 1.0_wp, 5, failures)
  call run_case('N', 'T', 8, 11, 6, 0.375_wp, -1.0_wp, 2, failures)
  call run_case('N', 'N', 9, 7, 11, 0.375_wp, 1.0_wp, 3, failures)
  call run_case('C', 'T', 6, 9, 8, -1.0_wp, 0.375_wp, 4, failures)

  ! These QR-shaped cases exceed the OpenMP work threshold.  In an OpenMP
  ! build they verify the threaded TN and NT paths with padded storage while
  ! the same source remains a serial regression test in the default build.
  call run_case('T', 'N', 256, 32, 256, 1.0_wp, 1.0_wp, 3, failures)
  call run_case('N', 'T', 256, 256, 32, -1.0_wp, 1.0_wp, 3, failures)

  ! A quiet NaN detects an accidental read of C in the beta=0 paths: any
  ! arithmetic involving the prior value would contaminate the finite result.
  call run_beta_zero_nan_case('T', 'N', failures)
  call run_beta_zero_nan_case('N', 'T', failures)
  call run_beta_zero_nan_case('N', 'N', failures)
  call run_beta_zero_nan_case('T', 'T', failures)

  call finish_suite('bundled DGEMM semantics', failures)

contains

  ! Compare one valid DGEMM call with a direct reference implementation.
  ! Padding is initialized to a sentinel and checked after the call, so only
  ! the active M-by-N output may change.  The error allowance scales with K,
  ! machine epsilon, and an absolute product bound rather than requiring an
  ! order-dependent bitwise result.
  subroutine run_case(transa, transb, m, n, k, alpha, beta, padding, failures)
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, padding
    real(wp), intent(in) :: alpha, beta
    integer, intent(inout) :: failures
    real(wp), allocatable :: a(:,:), b(:,:), c(:,:), c_before(:,:), expected(:,:)
    real(wp) :: absolute_bound, error, product_sum, tolerance
    real(wp), parameter :: sentinel = -9876.5_wp
    integer :: a_columns, a_rows, b_columns, b_rows, i, j, l
    integer :: lda, ldb, ldc
    character(len=160) :: message
    logical :: active_values_match, padding_is_unchanged

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

    lda = max(1, a_rows) + padding
    ldb = max(1, b_rows) + padding + 1
    ldc = max(1, m) + padding + 2
    allocate(a(lda,max(1,a_columns)+1), b(ldb,max(1,b_columns)+1))
    allocate(c(ldc,max(1,n)+1), c_before(ldc,max(1,n)+1))
    allocate(expected(ldc,max(1,n)+1))

    call fill_matrix(a, 3)
    call fill_matrix(b, 11)
    c = sentinel
    do j = 1, n
      do i = 1, m
        c(i,j) = real(mod(17*i + 13*j, 29) - 14, wp) / 9.0_wp
      end do
    end do
    c_before = c
    expected = c

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
          expected(i,j) = alpha * product_sum + beta * c_before(i,j)
        end if
      end do
    end do

    call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)

    absolute_bound = abs(alpha) * real(max(1,k), wp) * &
                     max(1.0_wp, maxval(abs(a))) * &
                     max(1.0_wp, maxval(abs(b))) + &
                     abs(beta) * max(1.0_wp, maxval(abs(c_before)))
    tolerance = 64.0_wp * epsilon(1.0_wp) * real(max(1,k), wp) * &
                max(1.0_wp, absolute_bound)
    if (m > 0 .and. n > 0) then
      error = maxval(abs(c(1:m,1:n) - expected(1:m,1:n)))
      active_values_match = error <= tolerance
    else
      error = 0.0_wp
      active_values_match = .true.
    end if
    padding_is_unchanged = .true.
    do j = 1, size(c,2)
      do i = 1, size(c,1)
        if (i > m .or. j > n) then
          padding_is_unchanged = padding_is_unchanged .and. &
                                 abs(c(i,j) - c_before(i,j)) <= 0.0_wp
        end if
      end do
    end do

    write(message,'(a,a,a,a,3(a,i0),2(a,es10.2))') &
      'DGEMM ', transa, transb, ' case', ' m=', m, ' n=', n, ' k=', k, &
      ' alpha=', alpha, ' beta=', beta
    call check(active_values_match, trim(message) // ' matches reference', failures)
    call check(padding_is_unchanged, trim(message) // ' preserves padding', &
               failures)
  end subroutine run_case

  ! Verify beta=0 without defining the active input values of C.  Quiet NaNs
  ! make the check portable across builds that do not enable invalid traps.
  subroutine run_beta_zero_nan_case(transa, transb, failures)
    character, intent(in) :: transa, transb
    integer, intent(inout) :: failures
    integer, parameter :: m = 7, n = 8, k = 5, padding = 3
    real(wp), allocatable :: a(:,:), b(:,:), c(:,:)
    integer :: a_columns, a_rows, b_columns, b_rows, lda, ldb, ldc

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
    lda = max(1, a_rows) + padding
    ldb = max(1, b_rows) + padding
    ldc = m + padding
    allocate(a(lda,a_columns), b(ldb,b_columns), c(ldc,n))
    call fill_matrix(a, 5)
    call fill_matrix(b, 19)
    c = -123.0_wp
    c(1:m,1:n) = ieee_value(0.0_wp, ieee_quiet_nan)

    call dgemm(transa, transb, m, n, k, 0.375_wp, a, lda, b, ldb, &
               0.0_wp, c, ldc)
    call check(all(ieee_is_finite(c(1:m,1:n))), &
               'beta=0 does not read C for ' // transa // transb, failures)
  end subroutine run_beta_zero_nan_case

  ! Fill every physical array element, including leading-dimension padding,
  ! with bounded deterministic values suitable for every supported real kind.
  subroutine fill_matrix(matrix, seed)
    real(wp), intent(out) :: matrix(:,:)
    integer, intent(in) :: seed
    integer :: i, j

    do j = 1, size(matrix,2)
      do i = 1, size(matrix,1)
        matrix(i,j) = real(mod(seed + 7*i + 23*j, 37) - 18, wp) / 11.0_wp
      end do
    end do
  end subroutine fill_matrix

  ! Return op(matrix)(row,column) for the real BLAS transpose conventions.
  ! Conjugate transpose is identical to transpose for real data.
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

end program test_dgemm

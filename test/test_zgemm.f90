! Semantic tests for the bundled precision-generic complex matrix product.
!
! A direct scalar reference exercises every transpose and conjugation mode,
! padded leading dimensions, rectangular optimized paths, scalar special
! cases, and micro-kernel tails without relying on MATMUL or external BLAS.
program test_zgemm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, &
                                           ieee_value
  use qrlinalg, only: wp
  use test_support, only: check, finish_suite
  implicit none

  interface
    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc)
      import :: wp
      character, intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      complex(wp), intent(in) :: alpha, beta
      complex(wp), intent(in) :: a(lda,*), b(ldb,*)
      complex(wp), intent(inout) :: c(ldc,*)
    end subroutine zgemm
  end interface

  character, parameter :: transposes(3) = ['N', 'T', 'C']
  complex(wp), parameter :: scalars(4) = [ &
    cmplx(0.0_wp, 0.0_wp, kind=wp), cmplx(1.0_wp, 0.0_wp, kind=wp), &
    cmplx(-1.0_wp, 0.0_wp, kind=wp), cmplx(0.375_wp, -0.25_wp, kind=wp)]
  integer :: alpha_index, beta_index, failures, k, m, n
  integer :: transa_index, transb_index

  failures = 0
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

  ! Rectangular dimensions enter the complete CN and NC tiles and their tails.
  call run_case('C', 'N', 7, 8, 9, scalars(2), scalars(2), 3, failures)
  call run_case('C', 'N', 8, 7, 17, scalars(4), scalars(4), 4, failures)
  call run_case('N', 'C', 11, 8, 7, scalars(3), scalars(2), 5, failures)
  call run_case('N', 'C', 8, 11, 32, scalars(4), scalars(3), 2, failures)
  call run_case('N', 'N', 9, 7, 11, scalars(4), scalars(2), 3, failures)
  call run_case('T', 'C', 6, 9, 8, scalars(3), scalars(4), 4, failures)

  call run_beta_zero_nan_case('C', 'N', failures)
  call run_beta_zero_nan_case('N', 'C', failures)
  call run_beta_zero_nan_case('N', 'N', failures)
  call run_beta_zero_nan_case('T', 'T', failures)

  call finish_suite('bundled ZGEMM semantics', failures)

contains

  subroutine run_case(transa, transb, m, n, k, alpha, beta, padding, failures)
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, padding
    complex(wp), intent(in) :: alpha, beta
    integer, intent(inout) :: failures
    complex(wp), allocatable :: a(:,:), b(:,:), c(:,:), c_before(:,:)
    complex(wp), allocatable :: expected(:,:)
    real(wp) :: absolute_bound, error, tolerance
    complex(wp) :: product_sum
    complex(wp), parameter :: sentinel = &
      cmplx(-9876.5_wp, 123.25_wp, kind=wp)
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
        c(i,j) = cmplx(real(mod(17*i+13*j,29)-14,wp)/9.0_wp, &
                       real(mod(11*i+19*j,31)-15,wp)/10.0_wp, kind=wp)
      end do
    end do
    c_before = c
    expected = c

    do j = 1, n
      do i = 1, m
        product_sum = cmplx(0.0_wp, 0.0_wp, kind=wp)
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

    call zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    absolute_bound = abs(alpha) * real(max(1,k),wp) * &
                     max(1.0_wp,maxval(abs(a))) * &
                     max(1.0_wp,maxval(abs(b))) + &
                     abs(beta) * max(1.0_wp,maxval(abs(c_before)))
    tolerance = 128.0_wp * epsilon(1.0_wp) * real(max(1,k),wp) * &
                max(1.0_wp,absolute_bound)
    if (m > 0 .and. n > 0) then
      error = maxval(abs(c(1:m,1:n)-expected(1:m,1:n)))
      active_values_match = error <= tolerance
    else
      active_values_match = .true.
    end if
    padding_is_unchanged = .true.
    do j = 1, size(c,2)
      do i = 1, size(c,1)
        if (i > m .or. j > n) then
          padding_is_unchanged = padding_is_unchanged .and. &
                                 abs(c(i,j)-c_before(i,j)) <= 0.0_wp
        end if
      end do
    end do
    write(message,'(a,a,a,a,3(a,i0))') 'ZGEMM ', transa, transb, ' case', &
      ' m=', m, ' n=', n, ' k=', k
    call check(active_values_match, trim(message)//' matches reference', failures)
    call check(padding_is_unchanged, trim(message)//' preserves padding', failures)
  end subroutine run_case

  subroutine run_beta_zero_nan_case(transa, transb, failures)
    character, intent(in) :: transa, transb
    integer, intent(inout) :: failures
    integer, parameter :: m=7, n=8, k=5, padding=3
    complex(wp), allocatable :: a(:,:), b(:,:), c(:,:)
    real(wp) :: quiet_nan
    integer :: a_columns, a_rows, b_columns, b_rows, lda, ldb, ldc

    if (transa == 'N') then
      a_rows=m; a_columns=k
    else
      a_rows=k; a_columns=m
    end if
    if (transb == 'N') then
      b_rows=k; b_columns=n
    else
      b_rows=n; b_columns=k
    end if
    lda=max(1,a_rows)+padding; ldb=max(1,b_rows)+padding; ldc=m+padding
    allocate(a(lda,a_columns), b(ldb,b_columns), c(ldc,n))
    call fill_matrix(a,5)
    call fill_matrix(b,19)
    quiet_nan = ieee_value(0.0_wp,ieee_quiet_nan)
    c = cmplx(-123.0_wp,45.0_wp,kind=wp)
    c(1:m,1:n) = cmplx(quiet_nan,quiet_nan,kind=wp)
    call zgemm(transa,transb,m,n,k,cmplx(0.375_wp,-0.25_wp,kind=wp), &
                a,lda,b,ldb,cmplx(0.0_wp,0.0_wp,kind=wp),c,ldc)
    call check(all(ieee_is_finite(real(c(1:m,1:n),wp))) .and. &
               all(ieee_is_finite(aimag(c(1:m,1:n)))), &
               'beta=0 does not read C for '//transa//transb, failures)
  end subroutine run_beta_zero_nan_case

  subroutine fill_matrix(matrix, seed)
    complex(wp), intent(out) :: matrix(:,:)
    integer, intent(in) :: seed
    integer :: i, j
    do j = 1, size(matrix,2)
      do i = 1, size(matrix,1)
        matrix(i,j) = cmplx(real(mod(seed+7*i+23*j,37)-18,wp)/11.0_wp, &
                            real(mod(seed+13*i+17*j,41)-20,wp)/12.0_wp, &
                            kind=wp)
      end do
    end do
  end subroutine fill_matrix

  function matrix_element(matrix, transpose_mode, row, column) result(value)
    complex(wp), intent(in) :: matrix(:,:)
    character, intent(in) :: transpose_mode
    integer, intent(in) :: row, column
    complex(wp) :: value
    select case (transpose_mode)
    case ('N')
      value = matrix(row,column)
    case ('C')
      value = conjg(matrix(column,row))
    case default
      value = matrix(column,row)
    end select
  end function matrix_element

end program test_zgemm

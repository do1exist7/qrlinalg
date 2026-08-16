! Differential-test driver for the QR-backed inverse-iteration library.
!
! The driver reads one symmetric or Hermitian generalized eigenproblem in the
! binary upper-triangle format used by orig/claude/main.f90. Runtime arguments
! select the scalar domain and inverse-iteration controls:
!
!   qrlinalg_driver KIND H_FILE S_FILE SHIFT TOL MAX_ITER NORM_MODE
!
! KIND is "real" or "complex". Matrix values in the files are binary64 and
! are converted to the compile-time working kind wp while reading. Output is a
! tagged, whitespace-separated format consumed by compare_outputs.py.
program qrlinalg_differential_driver
  use iso_fortran_env, only: error_unit, real64
  use qrlinalg
  implicit none

  character(len=16) :: problem_kind
  character(len=1024) :: h_file, s_file
  real(wp) :: shift, tolerance
  integer :: max_iter, norm_mode

  call read_arguments(problem_kind, h_file, s_file, shift, tolerance, &
                      max_iter, norm_mode)

  select case (trim(problem_kind))
  case ('real')
    call solve_real_case(h_file, s_file, shift, tolerance, max_iter, norm_mode)
  case ('complex')
    call solve_complex_case(h_file, s_file, shift, tolerance, max_iter, norm_mode)
  case default
    call fail('KIND must be real or complex')
  end select

contains

  ! Parse and validate the command-line interface shared with the original
  ! implementation's differential driver.
  subroutine read_arguments(problem_kind, h_file, s_file, shift, tolerance, &
                            max_iter, norm_mode)
    character(len=*), intent(out) :: problem_kind, h_file, s_file
    real(wp), intent(out) :: shift, tolerance
    integer, intent(out) :: max_iter, norm_mode
    character(len=256) :: argument
    integer :: io_status

    if (command_argument_count() /= 7) then
      call fail('expected KIND H_FILE S_FILE SHIFT TOL MAX_ITER NORM_MODE')
    end if

    call get_command_argument(1, problem_kind)
    call get_command_argument(2, h_file)
    call get_command_argument(3, s_file)
    call get_command_argument(4, argument)
    read(argument, *, iostat=io_status) shift
    if (io_status /= 0) call fail('SHIFT is not a working-precision real')
    call get_command_argument(5, argument)
    read(argument, *, iostat=io_status) tolerance
    if (io_status /= 0) call fail('TOL is not a working-precision real')
    call get_command_argument(6, argument)
    read(argument, *, iostat=io_status) max_iter
    if (io_status /= 0 .or. max_iter <= 0) &
      call fail('MAX_ITER must be a positive integer')
    call get_command_argument(7, argument)
    read(argument, *, iostat=io_status) norm_mode
    if (io_status /= 0) call fail('NORM_MODE must be an integer')
  end subroutine read_arguments

  ! Factorize and solve one real generalized symmetric problem. The residual
  ! uses the pristine input matrices and is normalized by
  !
  !   (||H||_F + |lambda|*||S||_F) * ||x||_2.
  subroutine solve_real_case(h_file, s_file, shift, tolerance, max_iter, &
                             norm_mode)
    character(len=*), intent(in) :: h_file, s_file
    real(wp), intent(in) :: shift, tolerance
    integer, intent(in) :: max_iter, norm_mode
    type(qr_real_state) :: state
    real(wp), allocatable :: h(:,:), s(:,:), initial(:), eigenvector(:)
    real(wp), allocatable :: residual_vector(:)
    real(wp) :: eigenvalue, relative_accuracy, residual
    integer :: h_order, s_order, info, num_iter

    call read_real_matrix(h_file, h, h_order)
    call read_real_matrix(s_file, s, s_order)
    if (h_order /= s_order) call fail('H and S have different orders')

    allocate(initial(h_order), eigenvector(h_order), residual_vector(h_order))
    initial = 1.0_wp
    eigenvector = 0.0_wp
    eigenvalue = 0.0_wp
    relative_accuracy = huge(1.0_wp)
    residual = huge(1.0_wp)
    num_iter = 0

    call state%initialize(h_order, info)
    if (info == QR_SUCCESS) &
      call state%factorize_fresh(h, s, shift, info)
    if (info == QR_SUCCESS) then
      call state%solve(s, initial, eigenvector, eigenvalue, tolerance, &
                       max_iter, norm_mode, relative_accuracy, num_iter, info)
    end if

    if (info == QR_SUCCESS .or. info == QR_ERR_NO_CONVERGENCE) then
      residual_vector = matmul(h, eigenvector) - &
                        eigenvalue * matmul(s, eigenvector)
      residual = vector_norm_real(residual_vector) / &
                 ((matrix_norm_real(h) + abs(eigenvalue) * &
                   matrix_norm_real(s)) * vector_norm_real(eigenvector))
    end if

    call write_real_result(h_order, info, num_iter, eigenvalue, &
                           relative_accuracy, residual, eigenvector)
  end subroutine solve_real_case

  ! Factorize and solve one complex generalized Hermitian problem. The
  ! residual definition is the complex analogue of solve_real_case.
  subroutine solve_complex_case(h_file, s_file, shift, tolerance, max_iter, &
                                norm_mode)
    character(len=*), intent(in) :: h_file, s_file
    real(wp), intent(in) :: shift, tolerance
    integer, intent(in) :: max_iter, norm_mode
    type(qr_complex_state) :: state
    complex(wp), allocatable :: h(:,:), s(:,:), initial(:), eigenvector(:)
    complex(wp), allocatable :: residual_vector(:)
    real(wp) :: eigenvalue, relative_accuracy, residual
    integer :: h_order, s_order, info, num_iter

    call read_complex_matrix(h_file, h, h_order)
    call read_complex_matrix(s_file, s, s_order)
    if (h_order /= s_order) call fail('H and S have different orders')

    allocate(initial(h_order), eigenvector(h_order), residual_vector(h_order))
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    eigenvector = cmplx(0.0_wp, 0.0_wp, kind=wp)
    eigenvalue = 0.0_wp
    relative_accuracy = huge(1.0_wp)
    residual = huge(1.0_wp)
    num_iter = 0

    call state%initialize(h_order, info)
    if (info == QR_SUCCESS) &
      call state%factorize_fresh(h, s, shift, info)
    if (info == QR_SUCCESS) then
      call state%solve(s, initial, eigenvector, eigenvalue, tolerance, &
                       max_iter, norm_mode, relative_accuracy, num_iter, info)
    end if

    if (info == QR_SUCCESS .or. info == QR_ERR_NO_CONVERGENCE) then
      residual_vector = matmul(h, eigenvector) - &
                        cmplx(eigenvalue, 0.0_wp, kind=wp) * &
                        matmul(s, eigenvector)
      residual = vector_norm_complex(residual_vector) / &
                 ((matrix_norm_complex(h) + abs(eigenvalue) * &
                   matrix_norm_complex(s)) * &
                  vector_norm_complex(eigenvector))
    end if

    call write_complex_result(h_order, info, num_iter, eigenvalue, &
                              relative_accuracy, residual, eigenvector)
  end subroutine solve_complex_case

  ! Read a symmetric matrix stored as a four-byte default integer order
  ! followed by upper-triangle binary64 columns. Conversion to wp occurs one
  ! column at a time, keeping auxiliary storage linear in the matrix order.
  subroutine read_real_matrix(file_name, matrix, matrix_order)
    character(len=*), intent(in) :: file_name
    real(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_order
    real(real64), allocatable :: column(:)
    integer :: column_index, row_index, io_status, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) call fail('cannot open real matrix: ' // trim(file_name))
    read(unit, iostat=io_status) matrix_order
    if (io_status /= 0 .or. matrix_order <= 0) &
      call fail('invalid matrix order in: ' // trim(file_name))
    allocate(matrix(matrix_order,matrix_order), column(matrix_order))
    do column_index = 1, matrix_order
      read(unit, iostat=io_status) column(1:column_index)
      if (io_status /= 0) call fail('truncated matrix file: ' // trim(file_name))
      matrix(1:column_index,column_index) = column(1:column_index)
    end do
    close(unit)

    do column_index = 1, matrix_order
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = matrix(row_index,column_index)
      end do
    end do
  end subroutine read_real_matrix

  ! Read the Hermitian equivalent of read_real_matrix. Stored binary64
  ! complex values are mirrored with conjugation to reconstruct the lower
  ! triangle used by both implementations.
  subroutine read_complex_matrix(file_name, matrix, matrix_order)
    character(len=*), intent(in) :: file_name
    complex(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_order
    complex(real64), allocatable :: column(:)
    integer :: column_index, row_index, io_status, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) &
      call fail('cannot open complex matrix: ' // trim(file_name))
    read(unit, iostat=io_status) matrix_order
    if (io_status /= 0 .or. matrix_order <= 0) &
      call fail('invalid matrix order in: ' // trim(file_name))
    allocate(matrix(matrix_order,matrix_order), column(matrix_order))
    do column_index = 1, matrix_order
      read(unit, iostat=io_status) column(1:column_index)
      if (io_status /= 0) call fail('truncated matrix file: ' // trim(file_name))
      matrix(1:column_index,column_index) = column(1:column_index)
    end do
    close(unit)

    do column_index = 1, matrix_order
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = conjg(matrix(row_index,column_index))
      end do
    end do
  end subroutine read_complex_matrix

  ! Write a complete real result using a precision-independent tagged format.
  subroutine write_real_result(matrix_order, info, num_iter, eigenvalue, &
                               relative_accuracy, residual, eigenvector)
    integer, intent(in) :: matrix_order, info, num_iter
    real(wp), intent(in) :: eigenvalue, relative_accuracy, residual
    real(wp), intent(in) :: eigenvector(:)
    integer :: index

    write(*,'(a,1x,a,3(1x,i0),4(1x,es45.36e4))') &
      'DIFF_RESULT', 'real', matrix_order, info, num_iter, eigenvalue, &
      relative_accuracy, residual, epsilon(1.0_wp)
    do index = 1, matrix_order
      write(*,'(a,1x,i0,2(1x,es45.36e4))') &
        'DIFF_VECTOR', index, eigenvector(index), 0.0_wp
    end do
  end subroutine write_real_result

  ! Write a complete complex result with separate real and imaginary columns.
  subroutine write_complex_result(matrix_order, info, num_iter, eigenvalue, &
                                  relative_accuracy, residual, eigenvector)
    integer, intent(in) :: matrix_order, info, num_iter
    real(wp), intent(in) :: eigenvalue, relative_accuracy, residual
    complex(wp), intent(in) :: eigenvector(:)
    integer :: index

    write(*,'(a,1x,a,3(1x,i0),4(1x,es45.36e4))') &
      'DIFF_RESULT', 'complex', matrix_order, info, num_iter, eigenvalue, &
      relative_accuracy, residual, epsilon(1.0_wp)
    do index = 1, matrix_order
      write(*,'(a,1x,i0,2(1x,es45.36e4))') 'DIFF_VECTOR', index, &
        real(eigenvector(index), wp), aimag(eigenvector(index))
    end do
  end subroutine write_complex_result

  pure function vector_norm_real(vector) result(norm)
    real(wp), intent(in) :: vector(:)
    real(wp) :: norm
    norm = sqrt(sum(vector * vector))
  end function vector_norm_real

  pure function vector_norm_complex(vector) result(norm)
    complex(wp), intent(in) :: vector(:)
    real(wp) :: norm
    norm = sqrt(sum(abs(vector)**2))
  end function vector_norm_complex

  pure function matrix_norm_real(matrix) result(norm)
    real(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm
    norm = sqrt(sum(matrix * matrix))
  end function matrix_norm_real

  pure function matrix_norm_complex(matrix) result(norm)
    complex(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm
    norm = sqrt(sum(abs(matrix)**2))
  end function matrix_norm_complex

  subroutine fail(message)
    character(len=*), intent(in) :: message
    write(error_unit,'(a)') 'qrlinalg differential driver: ' // trim(message)
    error stop 2
  end subroutine fail

end program qrlinalg_differential_driver

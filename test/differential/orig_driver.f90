! Differential-test driver for the original GSEPIIS and GHEPIIS routines.
!
! This program deliberately provides only the full-factorization path k=1,
! which is mathematically comparable with qrlinalg%factorize_fresh followed by
! qrlinalg%solve. The benchmark and update experiments in
! orig/claude/main.f90 are not part of the differential comparison.
program original_differential_driver
  use iso_fortran_env, only: error_unit, real64
  use linalg
  implicit none

  character(len=16) :: problem_kind
  character(len=1024) :: h_file, s_file
  real(wp) :: shift, tolerance
  integer :: max_iter, norm_mode

  Glob_NumOfProcs = 1
  Glob_ProcID = 0
  Glob_MPIErrCode = 0

  call read_arguments(problem_kind, h_file, s_file, shift, tolerance, &
                      max_iter, norm_mode)

  select case (trim(problem_kind))
  case ('real')
    call solve_real_case(h_file, s_file, shift, tolerance, max_iter, norm_mode)
  case ('complex')
    call solve_complex_case(h_file, s_file, shift, tolerance, max_iter, norm_mode)
  case default
    call abort_driver('KIND must be real or complex')
  end select

contains

  ! Parse the command line used by the serial qrlinalg driver.
  subroutine read_arguments(problem_kind, h_file, s_file, shift, tolerance, &
                            max_iter, norm_mode)
    character(len=*), intent(out) :: problem_kind, h_file, s_file
    real(wp), intent(out) :: shift, tolerance
    integer, intent(out) :: max_iter, norm_mode
    character(len=256) :: argument
    integer :: io_status

    if (command_argument_count() /= 7) then
      call abort_driver('expected KIND H_FILE S_FILE SHIFT TOL MAX_ITER NORM_MODE')
    end if

    call get_command_argument(1, problem_kind)
    call get_command_argument(2, h_file)
    call get_command_argument(3, s_file)
    call get_command_argument(4, argument)
    read(argument, *, iostat=io_status) shift
    if (io_status /= 0) call abort_driver('SHIFT is not a working-precision real')
    call get_command_argument(5, argument)
    read(argument, *, iostat=io_status) tolerance
    if (io_status /= 0) call abort_driver('TOL is not a working-precision real')
    call get_command_argument(6, argument)
    read(argument, *, iostat=io_status) max_iter
    if (io_status /= 0 .or. max_iter <= 0) &
      call abort_driver('MAX_ITER must be a positive integer')
    call get_command_argument(7, argument)
    read(argument, *, iostat=io_status) norm_mode
    if (io_status /= 0) call abort_driver('NORM_MODE must be an integer')
  end subroutine read_arguments

  ! Solve one real problem with the original full-factorization interface.
  ! Pristine H and S are retained solely for the generalized residual; GSEPIIS
  ! receives independent mutable M, B, and starting-vector arrays.
  subroutine solve_real_case(h_file, s_file, shift, tolerance, max_iter, &
                             norm_mode)
    character(len=*), intent(in) :: h_file, s_file
    real(wp), intent(in) :: shift, tolerance
    integer, intent(in) :: max_iter, norm_mode
    real(wp), allocatable :: h(:,:), s(:,:), matrix(:,:), overlap(:,:)
    real(wp), allocatable :: inverse_diagonal(:), initial(:), work(:)
    real(wp), allocatable :: eigenvector(:), residual_vector(:)
    real(wp) :: eigenvalue, relative_accuracy, residual
    integer :: h_order, s_order, error_code, num_iter

    call read_real_matrix(h_file, h, h_order)
    call read_real_matrix(s_file, s, s_order)
    if (h_order /= s_order) call abort_driver('H and S have different orders')
    call linalg_setparam(h_order)

    allocate(matrix(h_order,h_order), overlap(h_order,h_order))
    allocate(inverse_diagonal(h_order), initial(h_order), work(h_order))
    allocate(eigenvector(h_order), residual_vector(h_order))
    matrix = h - shift * s
    overlap = s
    inverse_diagonal = 0.0_wp
    initial = 1.0_wp
    work = 0.0_wp
    eigenvector = 0.0_wp

    call GSEPIIS(1, h_order, matrix, h_order, inverse_diagonal, overlap, &
                 h_order, shift, initial, work, tolerance, eigenvalue, &
                 eigenvector, relative_accuracy, max_iter, norm_mode, &
                 num_iter, error_code)

    residual = huge(1.0_wp)
    if (error_code == 0 .or. error_code == 2) then
      residual_vector = matmul(h, eigenvector) - &
                        eigenvalue * matmul(s, eigenvector)
      residual = vector_norm_real(residual_vector) / &
                 ((matrix_norm_real(h) + abs(eigenvalue) * &
                   matrix_norm_real(s)) * vector_norm_real(eigenvector))
    end if

    if (Glob_ProcID == 0) then
      call write_real_result(h_order, error_code, num_iter, eigenvalue, &
                             relative_accuracy, residual, eigenvector)
    end if
  end subroutine solve_real_case

  ! Solve the Hermitian counterpart with GHEPIIS and preserve the same
  ! residual and output conventions as the real path.
  subroutine solve_complex_case(h_file, s_file, shift, tolerance, max_iter, &
                                norm_mode)
    character(len=*), intent(in) :: h_file, s_file
    real(wp), intent(in) :: shift, tolerance
    integer, intent(in) :: max_iter, norm_mode
    complex(wp), allocatable :: h(:,:), s(:,:), matrix(:,:), overlap(:,:)
    complex(wp), allocatable :: inverse_diagonal(:), initial(:), work(:)
    complex(wp), allocatable :: eigenvector(:), residual_vector(:)
    real(wp) :: eigenvalue, relative_accuracy, residual
    integer :: h_order, s_order, error_code, num_iter

    call read_complex_matrix(h_file, h, h_order)
    call read_complex_matrix(s_file, s, s_order)
    if (h_order /= s_order) call abort_driver('H and S have different orders')
    call linalg_setparam(h_order)

    allocate(matrix(h_order,h_order), overlap(h_order,h_order))
    allocate(inverse_diagonal(h_order), initial(h_order), work(h_order))
    allocate(eigenvector(h_order), residual_vector(h_order))
    matrix = h - cmplx(shift, 0.0_wp, kind=wp) * s
    overlap = s
    inverse_diagonal = cmplx(0.0_wp, 0.0_wp, kind=wp)
    initial = cmplx(1.0_wp, 0.0_wp, kind=wp)
    work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    eigenvector = cmplx(0.0_wp, 0.0_wp, kind=wp)

    call GHEPIIS(1, h_order, matrix, h_order, inverse_diagonal, overlap, &
                 h_order, shift, initial, work, tolerance, eigenvalue, &
                 eigenvector, relative_accuracy, max_iter, norm_mode, &
                 num_iter, error_code)

    residual = huge(1.0_wp)
    if (error_code == 0 .or. error_code == 2) then
      residual_vector = matmul(h, eigenvector) - &
                        cmplx(eigenvalue, 0.0_wp, kind=wp) * &
                        matmul(s, eigenvector)
      residual = vector_norm_complex(residual_vector) / &
                 ((matrix_norm_complex(h) + abs(eigenvalue) * &
                   matrix_norm_complex(s)) * &
                  vector_norm_complex(eigenvector))
    end if

    if (Glob_ProcID == 0) then
      call write_complex_result(h_order, error_code, num_iter, eigenvalue, &
                                relative_accuracy, residual, eigenvector)
    end if
  end subroutine solve_complex_case

  ! Read a symmetric upper-triangle binary64 matrix. The file layout is
  ! identical to the existing original demonstration driver.
  subroutine read_real_matrix(file_name, matrix, matrix_order)
    character(len=*), intent(in) :: file_name
    real(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_order
    real(real64), allocatable :: column(:)
    integer :: column_index, row_index, io_status, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) &
      call abort_driver('cannot open real matrix: ' // trim(file_name))
    read(unit, iostat=io_status) matrix_order
    if (io_status /= 0 .or. matrix_order <= 0) &
      call abort_driver('invalid matrix order in: ' // trim(file_name))
    allocate(matrix(matrix_order,matrix_order), column(matrix_order))
    do column_index = 1, matrix_order
      read(unit, iostat=io_status) column(1:column_index)
      if (io_status /= 0) &
        call abort_driver('truncated matrix file: ' // trim(file_name))
      matrix(1:column_index,column_index) = column(1:column_index)
    end do
    close(unit)
    do column_index = 1, matrix_order
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = matrix(row_index,column_index)
      end do
    end do
  end subroutine read_real_matrix

  ! Read the Hermitian equivalent of read_real_matrix and reconstruct the
  ! unstored lower triangle with conjugate symmetry.
  subroutine read_complex_matrix(file_name, matrix, matrix_order)
    character(len=*), intent(in) :: file_name
    complex(wp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: matrix_order
    complex(real64), allocatable :: column(:)
    integer :: column_index, row_index, io_status, unit

    open(newunit=unit, file=trim(file_name), status='old', action='read', &
         access='stream', form='unformatted', iostat=io_status)
    if (io_status /= 0) &
      call abort_driver('cannot open complex matrix: ' // trim(file_name))
    read(unit, iostat=io_status) matrix_order
    if (io_status /= 0 .or. matrix_order <= 0) &
      call abort_driver('invalid matrix order in: ' // trim(file_name))
    allocate(matrix(matrix_order,matrix_order), column(matrix_order))
    do column_index = 1, matrix_order
      read(unit, iostat=io_status) column(1:column_index)
      if (io_status /= 0) &
        call abort_driver('truncated matrix file: ' // trim(file_name))
      matrix(1:column_index,column_index) = column(1:column_index)
    end do
    close(unit)
    do column_index = 1, matrix_order
      do row_index = 1, column_index - 1
        matrix(column_index,row_index) = conjg(matrix(row_index,column_index))
      end do
    end do
  end subroutine read_complex_matrix

  subroutine write_real_result(matrix_order, error_code, num_iter, eigenvalue, &
                               relative_accuracy, residual, eigenvector)
    integer, intent(in) :: matrix_order, error_code, num_iter
    real(wp), intent(in) :: eigenvalue, relative_accuracy, residual
    real(wp), intent(in) :: eigenvector(:)
    integer :: index

    write(*,'(a,1x,a,3(1x,i0),4(1x,es45.36e4))') &
      'DIFF_RESULT', 'real', matrix_order, error_code, num_iter, eigenvalue, &
      relative_accuracy, residual, epsilon(1.0_wp)
    do index = 1, matrix_order
      write(*,'(a,1x,i0,2(1x,es45.36e4))') &
        'DIFF_VECTOR', index, eigenvector(index), 0.0_wp
    end do
  end subroutine write_real_result

  subroutine write_complex_result(matrix_order, error_code, num_iter, &
                                  eigenvalue, relative_accuracy, residual, &
                                  eigenvector)
    integer, intent(in) :: matrix_order, error_code, num_iter
    real(wp), intent(in) :: eigenvalue, relative_accuracy, residual
    complex(wp), intent(in) :: eigenvector(:)
    integer :: index

    write(*,'(a,1x,a,3(1x,i0),4(1x,es45.36e4))') &
      'DIFF_RESULT', 'complex', matrix_order, error_code, num_iter, eigenvalue, &
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

  ! Terminate the standalone reference driver after a concise diagnostic.
  subroutine abort_driver(message)
    character(len=*), intent(in) :: message
    write(error_unit,'(a)') 'original differential driver: ' // trim(message)
    error stop 2
  end subroutine abort_driver

end program original_differential_driver

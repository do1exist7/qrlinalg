!A compact example of the complete qrlinalg state lifecycle.
program qrlinalg_example
  use qrlinalg, only: QR_SUCCESS, qr_real_state, qr_status_message, wp
  implicit none

  type(qr_real_state), allocatable :: qr
  real(wp) :: h(2,2), s(2,2), h3(3,3), s3(3,3)
  real(wp) :: delta_h(2), delta_s(2), h_column(3), s_column(3)
  real(wp) :: initial(3), x(3), lambda, rel_acc
  integer :: info, num_iter

  !Only the lower triangles of H and S are required.
  h = 0.0_wp
  h(1,1) = 1.0_wp
  h(2,1) = 0.1_wp
  h(2,2) = 3.0_wp
  s = 0.0_wp
  s(1,1) = 1.0_wp
  s(2,2) = 1.0_wp
  initial = 1.0_wp

  !initialize reserves storage up to capacity; factorize_fresh selects the
  !active order and shift.
  allocate(qr)
  call qr%initialize(capacity=3, info=info)
  call check_info('initialize', info)
  call qr%factorize_fresh(h, s, shift=0.8_wp, info=info)
  call check_info('factorize', info)
  call solve_and_print('initial factors', s, initial(1:2), x(1:2))

  !The caller owns H and S, so apply the same physical replacement locally.
  delta_h = [0.05_wp, 0.02_wp]
  delta_s = 0.0_wp
  call qr%replace_symmetric(1, delta_h, delta_s, info)
  call check_info('replace', info)
  h(1,1) = h(1,1) + delta_h(1)
  h(2,1) = h(2,1) + delta_h(2)

  !Append a row and column using the spare capacity.
  h_column = [0.01_wp, 0.03_wp, 5.0_wp]
  s_column = [0.0_wp, 0.0_wp, 1.0_wp]
  call qr%append_symmetric(h_column, s_column, info)
  call check_info('append', info)
  h3 = 0.0_wp
  s3 = 0.0_wp
  h3(1:2,1:2) = h
  s3(1:2,1:2) = s
  h3(3,:) = h_column
  s3(3,:) = s_column
  call solve_and_print('after append', s3, initial, x)

  call qr%delete_symmetric(3, info)
  call check_info('delete', info)

  !A shift change is not a cheap row/column update. Reuse the allocation and
  !perform a fresh factorization of the current physical H and S.
  call qr%factorize_fresh(h, s, shift=2.5_wp, info=info)
  call check_info('refactorize at new shift', info)
  call solve_and_print('after shift change', s, initial(1:2), x(1:2))

  !CLEAR releases the private arrays while retaining the state object itself.
  !The same object can then begin a new lifetime with a different capacity.
  call qr%clear()
  call qr%initialize(capacity=2, info=info)
  call check_info('reinitialize', info)
  call qr%factorize_fresh(h, s, shift=0.5_wp, info=info)
  call check_info('factorize new state', info)
  call solve_and_print('new state', s, initial(1:2), x(1:2))

  !CLEAR permits deterministic early release. DEALLOCATE would also release
  !all allocatable components automatically if CLEAR were omitted here.
  call qr%clear()
  deallocate(qr)

contains

  subroutine solve_and_print(label, overlap, starting_vector, eigenvector)
  !Solve one generalized eigenproblem and print its principal diagnostics.
    character(len=*), intent(in) :: label
    real(wp), intent(in) :: overlap(:,:), starting_vector(:)
    real(wp), intent(out) :: eigenvector(:)

    call qr%solve(overlap, starting_vector, eigenvector, lambda, &
                  tol=1.0e-12_wp, max_iter=30, norm_mode=1, &
                  rel_acc=rel_acc, num_iter=num_iter, info=info)
    call check_info('solve', info)
    write(*,'(/,a)') trim(label)
    write(*,'(a,es18.10)') '  eigenvalue = ', lambda
    write(*,'(a,i0,a,es10.2)') '  iterations = ', num_iter, &
                               ', relative accuracy = ', rel_acc
  end subroutine solve_and_print

  subroutine check_info(operation, status)
  !Stop this example if a qrlinalg operation reports an error.
    character(len=*), intent(in) :: operation
    integer, intent(in) :: status

    if (status /= QR_SUCCESS) then
      write(*,'(a,a,a,i0,a,a)') 'ERROR: ', trim(operation), ', info=', status, &
                                ': ', trim(qr_status_message(status))
      error stop 1
    end if
  end subroutine check_info

end program qrlinalg_example

! Public status-code value tests for qrlinalg.
!
! Existing numeric values are part of the external interface. New diagnostic
! categories must therefore be appended without renumbering earlier statuses.
program test_status_codes
  use qrlinalg
  use test_support, only: check, finish_suite
  implicit none

  character(len=128), parameter :: expected_message(0:11) = [character(len=128) :: &
    'operation completed successfully', &
    'invalid scalar control, index, or mathematical input property', &
    'state storage allocation failed', &
    'requested numerical operation is not implemented', &
    'LAPACK workspace query or QR factorization failed', &
    'shifted factorization or generated iterate is numerically unusable', &
    'inverse iteration did not converge within max_iter; approximation is available', &
    'QR state is uninitialized, invalid, or missing owned storage', &
    'caller array extents are empty or incompatible', &
    'requested active order exceeds initialized capacity', &
    'inverse-iteration starting vector is numerically zero', &
    'overlap quadratic form is non-positive or numerically zero']
  integer :: failures, info

  failures = 0
  call check(QR_SUCCESS == 0, 'QR_SUCCESS retains value zero', failures)
  call check(QR_ERR_INVALID_ARGUMENT == 1, &
             'QR_ERR_INVALID_ARGUMENT retains value one', failures)
  call check(QR_ERR_ALLOCATION == 2, &
             'QR_ERR_ALLOCATION retains value two', failures)
  call check(QR_ERR_NOT_IMPLEMENTED == 3, &
             'QR_ERR_NOT_IMPLEMENTED retains value three', failures)
  call check(QR_ERR_FACTORIZATION == 4, &
             'QR_ERR_FACTORIZATION retains value four', failures)
  call check(QR_ERR_SINGULAR == 5, &
             'QR_ERR_SINGULAR retains value five', failures)
  call check(QR_ERR_NO_CONVERGENCE == 6, &
             'QR_ERR_NO_CONVERGENCE retains value six', failures)
  call check(QR_ERR_INVALID_STATE == 7, &
             'QR_ERR_INVALID_STATE has appended value seven', failures)
  call check(QR_ERR_DIMENSION_MISMATCH == 8, &
             'QR_ERR_DIMENSION_MISMATCH has appended value eight', failures)
  call check(QR_ERR_CAPACITY_EXCEEDED == 9, &
             'QR_ERR_CAPACITY_EXCEEDED has appended value nine', failures)
  call check(QR_ERR_ZERO_INITIAL_VECTOR == 10, &
             'QR_ERR_ZERO_INITIAL_VECTOR has appended value ten', failures)
  call check(QR_ERR_NONPOSITIVE_OVERLAP == 11, &
             'QR_ERR_NONPOSITIVE_OVERLAP has appended value eleven', failures)

  do info = QR_SUCCESS, QR_ERR_NONPOSITIVE_OVERLAP
    call check(trim(qr_status_message(info)) == trim(expected_message(info)), &
               'status message agrees with its stable description', failures)
  end do
  call check(trim(qr_status_message(-1)) == &
             'unrecognized qrlinalg status code', &
             'unknown status has a stable fallback description', failures)
  call check(len(qr_status_message(QR_SUCCESS)) == 96, &
             'status message result retains its documented fixed length', failures)
  call finish_suite('qrlinalg status codes', failures)

end program test_status_codes

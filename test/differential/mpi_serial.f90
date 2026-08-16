! Minimal compile-time MPI compatibility layer for the original linalg module.
!
! Differential testing executes exactly one serial process and calls
! linalg_setparam with Glob_NumOfProcs=1 before either eigensolver. The original
! module consequently selects only its serial numerical branches. These names
! satisfy references in parallel branches that remain present in the object
! file without introducing an MPI compiler or runtime dependency.
module mpi
  implicit none

  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: MPI_INTEGER = 1
  integer, parameter :: MPI_DOUBLE_PRECISION = 2
  integer, parameter :: MPI_REAL16 = 3
  integer, parameter :: MPI_SUM = 4
  integer, parameter :: MPI_MAX = 5

contains

  ! The serial path must never require a reduction. Unlimited polymorphic,
  ! assumed-rank buffers reproduce MPI's choice-buffer interface sufficiently
  ! for the compiler to check every otherwise-unused parallel call site.
  subroutine MPI_ALLREDUCE(send_buffer, receive_buffer, count, datatype, &
                           operation, communicator, error_code)
    type(*), dimension(..), intent(in) :: send_buffer, receive_buffer
    integer, intent(in) :: count, datatype, operation, communicator
    integer, intent(out) :: error_code
    error_code = 0
  end subroutine MPI_ALLREDUCE

  ! Broadcasts are likewise unreachable when all original mode flags select
  ! the single-process implementation.
  subroutine MPI_BCAST(buffer, count, datatype, root, communicator, error_code)
    type(*), dimension(..), intent(in) :: buffer
    integer, intent(in) :: count, datatype, root, communicator
    integer, intent(out) :: error_code
    error_code = 0
  end subroutine MPI_BCAST

  subroutine MPI_BARRIER(communicator, error_code)
    integer, intent(in) :: communicator
    integer, intent(out) :: error_code
    error_code = 0
  end subroutine MPI_BARRIER

  ! Provide the timing symbol referenced by original calibration code. The
  ! calibration is disabled for one process, but a valid implementation makes
  ! the compatibility module complete and safe for diagnostic timing calls.
  real(8) function MPI_WTIME()
    integer :: clock_count, clock_rate

    call system_clock(clock_count, clock_rate)
    MPI_WTIME = real(clock_count, 8) / real(clock_rate, 8)
  end function MPI_WTIME

end module mpi

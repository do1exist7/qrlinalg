! Public metadata-query tests for qrlinalg.
!
! This program uses only type-bound public queries to verify that an external
! owner can observe the QR lifecycle without accessing private components or
! maintaining a duplicate description of the factorization state.
program test_metadata
  use iso_fortran_env, only: int64
  use qrlinalg
  use test_support, only: check, finish_suite
  implicit none

  integer :: failures

  failures = 0
  call test_real_metadata(failures)
  call test_complex_metadata(failures)
  call finish_suite('qrlinalg public metadata', failures)

contains

  ! Verify the real-state queries before initialization, after allocation and
  ! factorization, across a structural update and fresh refactorization, and
  ! after reinitialization resets the state lifetime.
  subroutine test_real_metadata(failures)
    integer, intent(inout) :: failures
    type(qr_real_state) :: state
    real(wp) :: delta_h(2), delta_s(2), h(2,2), s(2,2)
    integer :: info

    call check(.not. state%is_valid() .and. state%order() == 0 .and. &
               state%get_capacity() == 0, &
               'default real state reports no storage or factors', failures)
    call check(abs(state%get_shift()) <= tiny(1.0_wp) .and. &
               state%get_update_count() == 0_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'default real state reports default numerical metadata', &
               failures)

    call state%initialize(4, info)
    call check(info == QR_SUCCESS .and. .not. state%is_valid() .and. &
               state%order() == 0 .and. state%get_capacity() == 4, &
               'initialized real state exposes reserved capacity only', &
               failures)

    h = 0.0_wp
    h(1,1) = 1.0_wp
    h(2,1) = 0.2_wp
    h(2,2) = 3.0_wp
    s = 0.0_wp
    s(1,1) = 1.0_wp
    s(2,2) = 1.0_wp
    call state%factorize_fresh(h, s, 0.25_wp, info)
    call check(info == QR_SUCCESS .and. state%is_valid() .and. &
               state%order() == 2 .and. state%get_capacity() == 4 .and. &
               abs(state%get_shift() - 0.25_wp) <= epsilon(1.0_wp), &
               'real fresh factorization publishes order and shift', failures)

    delta_h = [0.01_wp, -0.02_wp]
    delta_s = 0.0_wp
    call state%replace_symmetric(1, delta_h, delta_s, info)
    call check(info == QR_SUCCESS .and. &
               state%get_update_count() == 1_int64 .and. &
               state%get_updates_since_fresh() == 1_int64, &
               'real replacement publishes both update counters', failures)

    h(1,1) = h(1,1) + delta_h(1)
    h(2,1) = h(2,1) + delta_h(2)
    call state%factorize_fresh(h, s, 1.25_wp, info)
    call check(info == QR_SUCCESS .and. &
               abs(state%get_shift() - 1.25_wp) <= epsilon(1.0_wp) .and. &
               state%get_update_count() == 1_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'real fresh factorization resets only the refresh counter', &
               failures)

    call state%delete_symmetric(0, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. state%is_valid() .and. &
               state%order() == 2 .and. &
               state%get_update_count() == 1_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'rejected real update preserves queried metadata', failures)

    call state%initialize(3, info)
    call check(info == QR_SUCCESS .and. .not. state%is_valid() .and. &
               state%order() == 0 .and. state%get_capacity() == 3 .and. &
               abs(state%get_shift()) <= tiny(1.0_wp) .and. &
               state%get_update_count() == 0_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'real reinitialization resets the complete queried lifetime', &
               failures)

    call state%clear()
    call check(.not. state%is_valid() .and. state%order() == 0 .and. &
               state%get_capacity() == 0 .and. &
               abs(state%get_shift()) <= tiny(1.0_wp) .and. &
               state%get_update_count() == 0_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'real clear restores the default public metadata', failures)
    call state%clear()
    call check(.not. state%is_valid() .and. state%get_capacity() == 0, &
               'real clear is idempotent for an empty state', failures)

    call state%initialize(2, info)
    call state%factorize_fresh(h, s, -0.75_wp, info)
    call check(info == QR_SUCCESS .and. state%is_valid() .and. &
               state%order() == 2 .and. state%get_capacity() == 2 .and. &
               abs(state%get_shift() + 0.75_wp) <= epsilon(1.0_wp), &
               'real state can be initialized and factorized after clear', &
               failures)
  end subroutine test_real_metadata

  ! Verify that the complex state exposes the same metadata semantics and that
  ! its counters remain independent of the real implementation.
  subroutine test_complex_metadata(failures)
    integer, intent(inout) :: failures
    type(qr_complex_state) :: state
    complex(wp) :: delta_h(2), delta_s(2), h(2,2), s(2,2)
    integer :: info

    call check(.not. state%is_valid() .and. state%order() == 0 .and. &
               state%get_capacity() == 0 .and. &
               state%get_update_count() == 0_int64, &
               'default complex state reports empty metadata', failures)

    call state%initialize(4, info)
    h = cmplx(0.0_wp, 0.0_wp, kind=wp)
    h(1,1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    h(2,1) = cmplx(0.2_wp, -0.1_wp, kind=wp)
    h(2,2) = cmplx(3.0_wp, 0.0_wp, kind=wp)
    s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s(1,1) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    s(2,2) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call state%factorize_fresh(h, s, -0.5_wp, info)
    call check(info == QR_SUCCESS .and. state%is_valid() .and. &
               state%order() == 2 .and. state%get_capacity() == 4 .and. &
               abs(state%get_shift() + 0.5_wp) <= epsilon(1.0_wp), &
               'complex fresh factorization publishes order and shift', &
               failures)

    delta_h = [cmplx(0.01_wp, 0.0_wp, kind=wp), &
               cmplx(-0.02_wp, 0.03_wp, kind=wp)]
    delta_s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    call state%replace_symmetric(1, delta_h, delta_s, info)
    call check(info == QR_SUCCESS .and. &
               state%get_update_count() == 1_int64 .and. &
               state%get_updates_since_fresh() == 1_int64, &
               'complex replacement publishes both update counters', failures)

    call state%initialize(3, info)
    call check(info == QR_SUCCESS .and. .not. state%is_valid() .and. &
               state%order() == 0 .and. state%get_capacity() == 3 .and. &
               abs(state%get_shift()) <= tiny(1.0_wp) .and. &
               state%get_update_count() == 0_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'complex reinitialization resets queried metadata', failures)

    call state%clear()
    call check(.not. state%is_valid() .and. state%order() == 0 .and. &
               state%get_capacity() == 0 .and. &
               abs(state%get_shift()) <= tiny(1.0_wp) .and. &
               state%get_update_count() == 0_int64 .and. &
               state%get_updates_since_fresh() == 0_int64, &
               'complex clear restores the default public metadata', failures)
    call state%clear()
    call check(.not. state%is_valid() .and. state%get_capacity() == 0, &
               'complex clear is idempotent for an empty state', failures)

    call state%initialize(2, info)
    call state%factorize_fresh(h, s, 0.75_wp, info)
    call check(info == QR_SUCCESS .and. state%is_valid() .and. &
               state%order() == 2 .and. state%get_capacity() == 2 .and. &
               abs(state%get_shift() - 0.75_wp) <= epsilon(1.0_wp), &
               'complex state can be initialized and factorized after clear', &
               failures)
  end subroutine test_complex_metadata

end program test_metadata

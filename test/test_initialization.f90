! Initialization and empty-state contract tests for qrlinalg.
!
! This program is linked to a QRLINALG_TESTING build so it can verify private
! storage extents without exposing those components in the production API.
program test_initialization
  use qrlinalg
  use test_support, only: check, finish_suite
  implicit none

  integer :: failures

  failures = 0
  call test_invalid_initialization(failures)
  call test_real_initialization(failures)
  call test_complex_initialization(failures)
  call test_state_lifetimes(failures)
  call finish_suite('qrlinalg initialization', failures)

contains

  ! Verify rejection of invalid initialization and use of an uninitialized
  ! state. These are recoverable library errors and must never terminate the
  ! caller or leave partially allocated state behind.
  subroutine test_invalid_initialization(failures)
    integer, intent(inout) :: failures
    type(qr_real_state) :: state
    type(qr_complex_state) :: complex_state
    real(wp) :: h(1,1), s(1,1)
    complex(wp) :: complex_h(1,1), complex_s(1,1)
    integer :: info

    h = 1.0_wp
    s = 1.0_wp
    call state%factorize_fresh(h, s, 0.0_wp, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'factorization rejects an uninitialized state', failures)

    call state%initialize(0, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'initialize rejects zero capacity', failures)
    call check(state%capacity == 0 .and. state%n == 0, &
               'failed initialization leaves zero dimensions', failures)
    call check(.not. state%valid, &
               'failed initialization leaves factors invalid', failures)
    call check(.not. allocated(state%q) .and. .not. allocated(state%r), &
               'failed initialization owns no matrix storage', failures)

    complex_h = cmplx(1.0_wp, 0.0_wp, kind=wp)
    complex_s = cmplx(1.0_wp, 0.0_wp, kind=wp)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.0_wp, info)
    call check(info == QR_ERR_INVALID_ARGUMENT, &
               'complex factorization rejects an uninitialized state', failures)
    call complex_state%initialize(-1, info)
    call check(info == QR_ERR_INVALID_ARGUMENT .and. &
               complex_state%capacity == 0 .and. complex_state%n == 0 .and. &
               .not. complex_state%valid, &
               'complex initialize rejects negative capacity safely', failures)
    call check(.not. allocated(complex_state%q) .and. &
               .not. allocated(complex_state%r), &
               'failed complex initialization owns no matrix storage', failures)
  end subroutine test_invalid_initialization

  ! Verify all storage and metadata established for an initialized real state.
  ! Workspace sizes are part of the allocation contract even though the exact
  ! LAPACK factor workspace may exceed its documented minimum.
  subroutine test_real_initialization(failures)
    integer, intent(inout) :: failures
    integer, parameter :: capacity = 5
    type(qr_real_state) :: state
    integer :: info

    call state%initialize(capacity=capacity, info=info)
    call check(info == QR_SUCCESS, 'real initialize succeeds', failures)
    call check(state%capacity == capacity .and. state%n == 0, &
               'real initialize records capacity and zero active order', failures)
    call check(.not. state%valid, &
               'real initialize does not claim valid factors', failures)
    call check(all(shape(state%q) == [capacity, capacity]) .and. &
               all(shape(state%r) == [capacity, capacity]), &
               'real Q and R use full capacity storage', failures)
    call check(size(state%tau) == capacity, &
               'real tau has capacity elements', failures)
    call check(size(state%factor_work) >= capacity, &
               'real LAPACK workspace meets its minimum size', failures)
    call check(size(state%update_work) == 4 * capacity .and. &
               size(state%solve_work) == 2 * capacity, &
               'real update and solve workspaces have planned sizes', failures)
    call check(state%structural_updates == 0 .and. &
               state%updates_since_fresh == 0, &
               'real initialize clears update counters', failures)
  end subroutine test_real_initialization

  ! Verify the complex state owns the same logical storage as the real state,
  ! together with the real workspace required by complex QR-update rotations.
  subroutine test_complex_initialization(failures)
    integer, intent(inout) :: failures
    integer, parameter :: capacity = 5
    type(qr_complex_state) :: state
    integer :: info

    call state%initialize(capacity=capacity, info=info)
    call check(info == QR_SUCCESS, 'complex initialize succeeds', failures)
    call check(state%capacity == capacity .and. state%n == 0 .and. &
               .not. state%valid, &
               'complex initialize records empty-state metadata', failures)
    call check(all(shape(state%q) == [capacity, capacity]) .and. &
               all(shape(state%r) == [capacity, capacity]), &
               'complex Q and R use full capacity storage', failures)
    call check(size(state%tau) == capacity .and. &
               size(state%factor_work) >= capacity, &
               'complex factorization arrays meet LAPACK requirements', failures)
    call check(size(state%update_work) == 4 * capacity .and. &
               size(state%solve_work) == 2 * capacity .and. &
               size(state%real_work) == capacity, &
               'complex update, solve, and real workspaces have planned sizes', &
               failures)
    call check(state%structural_updates == 0 .and. &
               state%updates_since_fresh == 0, &
               'complex initialize clears update counters', failures)
  end subroutine test_complex_initialization

  ! Verify the explicit allocate-update-deallocate-allocate lifecycle used by
  ! long-running callers that want deterministic release of state storage.
  ! Intrinsic deallocation of the derived object recursively releases all
  ! private allocatable components, after which a new object may use a new
  ! capacity and represent factors at a different shift.
  subroutine test_state_lifetimes(failures)
    integer, intent(inout) :: failures
    type(qr_real_state), allocatable :: real_state
    type(qr_complex_state), allocatable :: complex_state
    real(wp) :: delta_h(2), delta_s(2), h(2,2), s(2,2)
    complex(wp) :: complex_delta_h(2), complex_delta_s(2)
    complex(wp) :: complex_h(2,2), complex_s(2,2)
    integer :: info

    h = 0.0_wp
    h(1,1) = 1.0_wp
    h(2,1) = 0.1_wp
    h(2,2) = 3.0_wp
    s = 0.0_wp
    s(1,1) = 1.0_wp
    s(2,2) = 1.0_wp
    delta_h = [0.01_wp, 0.02_wp]
    delta_s = 0.0_wp

    allocate(real_state)
    call real_state%initialize(3, info)
    call real_state%factorize_fresh(h, s, 0.5_wp, info)
    call real_state%replace_symmetric(1, delta_h, delta_s, info)
    call check(info == QR_SUCCESS .and. &
               real_state%updates_since_fresh == 1, &
               'real allocatable state completes an update', failures)
    deallocate(real_state)
    call check(.not. allocated(real_state), &
               'real state object deallocates explicitly', failures)

    allocate(real_state)
    call real_state%initialize(2, info)
    call real_state%factorize_fresh(h, s, 2.0_wp, info)
    call check(info == QR_SUCCESS .and. real_state%capacity == 2 .and. &
               real_state%n == 2 .and. &
               abs(real_state%shift - 2.0_wp) <= epsilon(1.0_wp) .and. &
               real_state%updates_since_fresh == 0, &
               'real state reallocates with new capacity and shift', failures)
    deallocate(real_state)

    complex_h = cmplx(h, 0.0_wp, kind=wp)
    complex_s = cmplx(s, 0.0_wp, kind=wp)
    complex_delta_h = cmplx(delta_h, 0.0_wp, kind=wp)
    complex_delta_s = cmplx(delta_s, 0.0_wp, kind=wp)
    allocate(complex_state)
    call complex_state%initialize(3, info)
    call complex_state%factorize_fresh(complex_h, complex_s, 0.5_wp, info)
    call complex_state%replace_symmetric(1, complex_delta_h, &
                                         complex_delta_s, info)
    call check(info == QR_SUCCESS .and. &
               complex_state%updates_since_fresh == 1, &
               'complex allocatable state completes an update', failures)
    deallocate(complex_state)
    call check(.not. allocated(complex_state), &
               'complex state object deallocates explicitly', failures)

    allocate(complex_state)
    call complex_state%initialize(2, info)
    call complex_state%factorize_fresh(complex_h, complex_s, 2.0_wp, info)
    call check(info == QR_SUCCESS .and. complex_state%capacity == 2 .and. &
               complex_state%n == 2 .and. &
               abs(complex_state%shift - 2.0_wp) <= epsilon(1.0_wp) .and. &
               complex_state%updates_since_fresh == 0, &
               'complex state reallocates with new capacity and shift', failures)
    deallocate(complex_state)
  end subroutine test_state_lifetimes

end program test_initialization

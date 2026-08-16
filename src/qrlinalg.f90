module qrlinalg
  use iso_fortran_env, only: int64
  use wp_def, only: wp
  implicit none
  private

  integer, parameter, public :: QR_SUCCESS = 0
  integer, parameter, public :: QR_ERR_INVALID_ARGUMENT = 1
  integer, parameter, public :: QR_ERR_ALLOCATION = 2
  integer, parameter, public :: QR_ERR_NOT_IMPLEMENTED = 3

  public :: wp
  public :: qr_real_state
  public :: qr_complex_state

  type, public :: qr_real_state
    private
    integer :: n = 0
    integer :: capacity = 0
    real(wp) :: shift = 0.0_wp
    logical :: valid = .false.
    real(wp), allocatable :: q(:,:)
    real(wp), allocatable :: r(:,:)
    real(wp), allocatable :: tau(:)
    real(wp), allocatable :: factor_work(:)
    real(wp), allocatable :: update_work(:)
    real(wp), allocatable :: solve_work(:)
    integer(int64) :: structural_updates = 0_int64
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => real_initialize
    procedure :: factorize_fresh => real_factorize_fresh
    procedure :: replace_symmetric => real_replace_symmetric
    procedure :: append_symmetric => real_append_symmetric
    procedure :: delete_symmetric => real_delete_symmetric
    procedure :: solve => real_solve
  end type qr_real_state

  type, public :: qr_complex_state
    private
    integer :: n = 0
    integer :: capacity = 0
    real(wp) :: shift = 0.0_wp
    logical :: valid = .false.
    complex(wp), allocatable :: q(:,:)
    complex(wp), allocatable :: r(:,:)
    complex(wp), allocatable :: tau(:)
    complex(wp), allocatable :: factor_work(:)
    complex(wp), allocatable :: update_work(:)
    complex(wp), allocatable :: solve_work(:)
    real(wp), allocatable :: real_work(:)
    integer(int64) :: structural_updates = 0_int64
    integer(int64) :: updates_since_fresh = 0_int64
  contains
    procedure :: initialize => complex_initialize
    procedure :: factorize_fresh => complex_factorize_fresh
    procedure :: replace_symmetric => complex_replace_symmetric
    procedure :: append_symmetric => complex_append_symmetric
    procedure :: delete_symmetric => complex_delete_symmetric
    procedure :: solve => complex_solve
  end type qr_complex_state

contains

  subroutine clear_real_state(self)
    class(qr_real_state), intent(inout) :: self

    if (allocated(self%q)) deallocate(self%q)
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%tau)) deallocate(self%tau)
    if (allocated(self%factor_work)) deallocate(self%factor_work)
    if (allocated(self%update_work)) deallocate(self%update_work)
    if (allocated(self%solve_work)) deallocate(self%solve_work)
    self%n = 0
    self%capacity = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_real_state

  subroutine clear_complex_state(self)
    class(qr_complex_state), intent(inout) :: self

    if (allocated(self%q)) deallocate(self%q)
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%tau)) deallocate(self%tau)
    if (allocated(self%factor_work)) deallocate(self%factor_work)
    if (allocated(self%update_work)) deallocate(self%update_work)
    if (allocated(self%solve_work)) deallocate(self%solve_work)
    if (allocated(self%real_work)) deallocate(self%real_work)
    self%n = 0
    self%capacity = 0
    self%shift = 0.0_wp
    self%valid = .false.
    self%structural_updates = 0_int64
    self%updates_since_fresh = 0_int64
  end subroutine clear_complex_state

  subroutine real_initialize(self, max_n, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status

    call clear_real_state(self)
    if (max_n <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    allocate(self%q(max_n, max_n), self%r(max_n, max_n), &
             self%tau(max_n), self%factor_work(max_n), &
             self%update_work(4 * max_n), self%solve_work(2 * max_n), &
             stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = 0.0_wp
    self%r = 0.0_wp
    self%tau = 0.0_wp
    self%factor_work = 0.0_wp
    self%update_work = 0.0_wp
    self%solve_work = 0.0_wp
    self%capacity = max_n
    info = QR_SUCCESS
  end subroutine real_initialize

  subroutine complex_initialize(self, max_n, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status

    call clear_complex_state(self)
    if (max_n <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    allocate(self%q(max_n, max_n), self%r(max_n, max_n), &
             self%tau(max_n), self%factor_work(max_n), &
             self%update_work(4 * max_n), self%solve_work(2 * max_n), &
             self%real_work(max_n), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%r = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%tau = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%factor_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%update_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%solve_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%real_work = 0.0_wp
    self%capacity = max_n
    info = QR_SUCCESS
  end subroutine complex_initialize

  subroutine real_factorize_fresh(self, h, s, shift, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_factorize_fresh

  subroutine complex_factorize_fresh(self, h, s, shift, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_factorize_fresh

  subroutine real_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_replace_symmetric

  subroutine complex_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_replace_symmetric

  subroutine real_append_symmetric(self, h_column, s_column, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_append_symmetric

  subroutine complex_append_symmetric(self, h_column, s_column, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_append_symmetric

  subroutine real_delete_symmetric(self, idx, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_delete_symmetric

  subroutine complex_delete_symmetric(self, idx, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_delete_symmetric

  subroutine real_solve(self, s, v_initial, x, lambda, tol, max_iter, &
                        norm_mode, rel_acc, num_iter, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: s(:,:), v_initial(:)
    real(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter, info

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = 0.0_wp
    num_iter = 0
    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_solve

  subroutine complex_solve(self, s, v_initial, x, lambda, tol, max_iter, &
                           norm_mode, rel_acc, num_iter, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: s(:,:), v_initial(:)
    complex(wp), intent(out) :: x(:)
    real(wp), intent(out) :: lambda
    real(wp), intent(in) :: tol
    integer, intent(in) :: max_iter, norm_mode
    real(wp), intent(out) :: rel_acc
    integer, intent(out) :: num_iter, info

    x = cmplx(0.0_wp, 0.0_wp, kind=wp)
    lambda = 0.0_wp
    rel_acc = 0.0_wp
    num_iter = 0
    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_solve

end module qrlinalg

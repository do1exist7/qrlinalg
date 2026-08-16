!> Serial, reusable QR state for shifted generalized eigenproblems.
!!
!! The module owns explicit QR factors and all workspace needed by fresh
!! factorization, structural QR updates, and inverse iteration. The caller
!! retains ownership of the physical H and S matrices. A state never stores
!! pointers to those matrices and never permanently stores H - shift*S.
module qrlinalg
  use iso_fortran_env, only: int64
  use wp_def, only: wp
  implicit none
  private

  integer, parameter, public :: QR_SUCCESS = 0
    !! The requested operation completed successfully.
  integer, parameter, public :: QR_ERR_INVALID_ARGUMENT = 1
    !! An argument, state, dimension, or capacity is invalid.
  integer, parameter, public :: QR_ERR_ALLOCATION = 2
    !! Initialization could not allocate all required state storage.
  integer, parameter, public :: QR_ERR_NOT_IMPLEMENTED = 3
    !! The requested numerical operation is a deliberate version 0.1 stub.
  integer, parameter, public :: QR_ERR_FACTORIZATION = 4
    !! A bundled LAPACK factorization or workspace query reported an error.
  integer, parameter, public :: QR_ERR_SINGULAR = 5
    !! The stored shifted factorization is singular or numerically unusable.
  integer, parameter, public :: QR_ERR_NO_CONVERGENCE = 6
    !! Inverse iteration reached max_iter before satisfying the tolerance.

  public :: wp
  public :: qr_real_state
  public :: qr_complex_state

  ! Explicit interfaces for the traditional external BLAS/LAPACK routines in
  ! src/qrupdate. These declarations provide compile-time checking; the linker
  ! resolves the external symbols from the bundled objects. Both declarations
  ! and implementations import the same build-selected wp_def.
  interface
    ! Compute a compact real QR factorization in A and TAU.
    subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, lda, lwork
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(out) :: tau(*)
      real(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dgeqrf

    ! Generate explicit real Q from the reflectors returned by DGEQRF.
    subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, k, lda, lwork
      real(wp), intent(inout) :: a(lda, *)
      real(wp), intent(in) :: tau(*)
      real(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dorgqr

    ! Compute a compact complex QR factorization in A and TAU.
    subroutine zgeqrf(m, n, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, lda, lwork
      complex(wp), intent(inout) :: a(lda, *)
      complex(wp), intent(out) :: tau(*)
      complex(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zgeqrf

    ! Generate explicit complex Q from the reflectors returned by ZGEQRF.
    subroutine zungqr(m, n, k, a, lda, tau, work, lwork, info)
      import wp
      integer, intent(in) :: m, n, k, lda, lwork
      complex(wp), intent(inout) :: a(lda, *)
      complex(wp), intent(in) :: tau(*)
      complex(wp), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine zungqr

    ! Multiply a real symmetric matrix by a vector using one stored triangle.
    subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda, incx, incy
      real(wp), intent(in) :: alpha, beta
      real(wp), intent(in) :: a(lda, *), x(*)
      real(wp), intent(inout) :: y(*)
    end subroutine dsymv

    ! Multiply a complex Hermitian matrix by a vector using one triangle.
    subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, lda, incx, incy
      complex(wp), intent(in) :: alpha, beta
      complex(wp), intent(in) :: a(lda, *), x(*)
      complex(wp), intent(inout) :: y(*)
    end subroutine zhemv

    ! General real matrix-vector product, used for Q and R applications.
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      real(wp), intent(in) :: alpha, beta
      real(wp), intent(in) :: a(lda, *), x(*)
      real(wp), intent(inout) :: y(*)
    end subroutine dgemv

    ! General complex matrix-vector product, including conjugate transpose.
    subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      import wp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      complex(wp), intent(in) :: alpha, beta
      complex(wp), intent(in) :: a(lda, *), x(*)
      complex(wp), intent(inout) :: y(*)
    end subroutine zgemv

    ! Solve an upper-triangular real system in place.
    subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      import wp
      character(len=1), intent(in) :: uplo, trans, diag
      integer, intent(in) :: n, lda, incx
      real(wp), intent(in) :: a(lda, *)
      real(wp), intent(inout) :: x(*)
    end subroutine dtrsv

    ! Solve an upper-triangular complex system in place.
    subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
      import wp
      character(len=1), intent(in) :: uplo, trans, diag
      integer, intent(in) :: n, lda, incx
      complex(wp), intent(in) :: a(lda, *)
      complex(wp), intent(inout) :: x(*)
    end subroutine ztrsv
  end interface

  type, public :: qr_real_state
#ifndef QRLINALG_TESTING
    private
#endif
    ! State components remain private in normal builds. The dedicated test
    ! build defines QRLINALG_TESTING so numerical invariants can be checked
    ! without adding factor-copying accessors to the production API.
    ! Current active order of the valid factors. Zero means no active factors.
    integer :: n = 0
    ! Largest order that fits in the allocated matrices and work arrays.
    integer :: capacity = 0
    ! Shift represented by Q*R = H - shift*S when valid is true.
    real(wp) :: shift = 0.0_wp
    ! True only when Q, R, n, and shift describe a successful factorization.
    logical :: valid = .false.
    ! Explicit orthogonal factor, stored with leading dimension capacity.
    real(wp), allocatable :: q(:,:)
    ! Explicit upper-triangular factor, also full leading-dimension storage.
    real(wp), allocatable :: r(:,:)
    ! Householder scalar factors produced temporarily by DGEQRF.
    real(wp), allocatable :: tau(:)
    ! Shared optimally sized workspace for DGEQRF and DORGQR.
    real(wp), allocatable :: factor_work(:)
    ! Storage for destructive qrupdate vector arguments and kernel workspace.
    real(wp), allocatable :: update_work(:)
    ! Two capacity-length vectors used by inverse iteration and quotients.
    real(wp), allocatable :: solve_work(:)
    ! Count of successful append, replacement, and deletion operations.
    integer(int64) :: structural_updates = 0_int64
    ! Successful structural updates since the last fresh factorization.
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
#ifndef QRLINALG_TESTING
    private
#endif
    ! See qr_real_state for the test-only component-access policy.
    ! Complex-state metadata has the same meaning as in qr_real_state.
    integer :: n = 0
    integer :: capacity = 0
    real(wp) :: shift = 0.0_wp
    logical :: valid = .false.
    ! Explicit unitary Q and upper-triangular R.
    complex(wp), allocatable :: q(:,:)
    complex(wp), allocatable :: r(:,:)
    ! Householder scalars and optimally sized complex LAPACK workspace.
    complex(wp), allocatable :: tau(:)
    complex(wp), allocatable :: factor_work(:)
    ! Destructive qrupdate copies and two-vector inverse-iteration workspace.
    complex(wp), allocatable :: update_work(:)
    complex(wp), allocatable :: solve_work(:)
    ! Real cosines and other real workspace required by complex qrupdate.
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

  !> Release all storage owned by a real state and restore its metadata and
  !! counters to the uninitialized defaults. This private helper makes
  !! reinitialization and allocation-failure cleanup safe.
  !!
  !! Input/output parameter:
  !!   self - State to clear. It is valid to pass an already empty state.
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

  !> Release all complex and real workspace owned by a complex state and
  !! restore its metadata and counters to the uninitialized defaults.
  !!
  !! Input/output parameter:
  !!   self - State to clear. It is valid to pass an already empty state.
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

  !> Initialize a real QR state for matrices up to order `max_n`.
  !!
  !! Existing storage and counters are cleared first. The routine allocates
  !! full Q and R matrices plus all reusable work arrays, but does not create
  !! valid factors or set an active dimension. DGEQRF and DORGQR workspace
  !! queries are made for square matrices of order `max_n`; factor_work is
  !! allocated once to the larger recommended size.
  !!
  !! Input parameter:
  !!   max_n - Maximum matrix order that this state must support; must be > 0.
  !!
  !! Output parameter:
  !!   info - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, QR_ERR_ALLOCATION, or
  !!          QR_ERR_FACTORIZATION if a LAPACK workspace query fails.
  subroutine real_initialize(self, max_n, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    real(wp) :: work_query(1)

    ! Reinitialization discards old factors and resets lifetime counters.
    call clear_real_state(self)
    if (max_n <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    ! Allocate all fixed-size storage except factor_work. Its size is selected
    ! below using LAPACK's query mode rather than assuming a block size.
    allocate(self%q(max_n, max_n), self%r(max_n, max_n), &
             self%tau(max_n), self%update_work(4 * max_n), &
             self%solve_work(2 * max_n), &
             stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = 0.0_wp
    self%r = 0.0_wp
    self%tau = 0.0_wp
    self%update_work = 0.0_wp
    self%solve_work = 0.0_wp

    ! LWORK=-1 asks LAPACK to return its preferred workspace in WORK(1)
    ! without factorizing or otherwise referencing the matrix contents.
    call dgeqrf(max_n, max_n, self%q, max_n, self%tau, work_query, -1, &
                lapack_info)
    if (lapack_info /= 0) then
      call clear_real_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    factor_lwork = max(max_n, ceiling(work_query(1)))

    call dorgqr(max_n, max_n, max_n, self%q, max_n, self%tau, &
                work_query, -1, lapack_info)
    if (lapack_info /= 0) then
      call clear_real_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    generate_q_lwork = max(max_n, ceiling(work_query(1)))

    optimal_lwork = max(factor_lwork, generate_q_lwork)
    allocate(self%factor_work(optimal_lwork), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_real_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if
    self%factor_work = 0.0_wp

    ! Allocation alone does not create factors; factorize_fresh commits n,
    ! shift, and valid only after both LAPACK stages succeed.
    self%capacity = max_n
    info = QR_SUCCESS
  end subroutine real_initialize

  !> Initialize a complex QR state for matrices up to order `max_n`.
  !!
  !! In addition to complex factorization, update, and solve storage, this
  !! allocates the real workspace needed by complex plane rotations. ZGEQRF
  !! and ZUNGQR workspace queries determine the reusable factor_work size.
  !! The state remains factor-invalid with active dimension zero until a
  !! successful fresh factorization.
  !!
  !! Input and output parameters have the same contracts and status codes as
  !! real_initialize.
  subroutine complex_initialize(self, max_n, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    complex(wp) :: work_query(1)

    ! Reinitialization discards old factors and resets lifetime counters.
    call clear_complex_state(self)
    if (max_n <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    allocate(self%q(max_n, max_n), self%r(max_n, max_n), &
             self%tau(max_n), self%update_work(4 * max_n), &
             self%solve_work(2 * max_n), &
             self%real_work(max_n), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if

    self%q = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%r = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%tau = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%update_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%solve_work = cmplx(0.0_wp, 0.0_wp, kind=wp)
    self%real_work = 0.0_wp

    ! Complex LAPACK returns the integer-valued recommendation in the real
    ! part of WORK(1). Use the larger recommendation from the two stages.
    call zgeqrf(max_n, max_n, self%q, max_n, self%tau, work_query, -1, &
                lapack_info)
    if (lapack_info /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    factor_lwork = max(max_n, ceiling(real(work_query(1), wp)))

    call zungqr(max_n, max_n, max_n, self%q, max_n, self%tau, &
                work_query, -1, lapack_info)
    if (lapack_info /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_FACTORIZATION
      return
    end if
    generate_q_lwork = max(max_n, ceiling(real(work_query(1), wp)))

    optimal_lwork = max(factor_lwork, generate_q_lwork)
    allocate(self%factor_work(optimal_lwork), stat=allocation_status)
    if (allocation_status /= 0) then
      call clear_complex_state(self)
      info = QR_ERR_ALLOCATION
      return
    end if
    self%factor_work = cmplx(0.0_wp, 0.0_wp, kind=wp)

    self%capacity = max_n
    info = QR_SUCCESS
  end subroutine complex_initialize

  !> Build fresh real factors of `H - shift*S` in state-owned Q and R.
  !!
  !! The shifted matrix is formed directly in Q, factored by `DGEQRF`, copied
  !! into the upper triangle of R, and replaced in Q by the explicit orthogonal
  !! factor from `DORGQR`. H and S remain caller-owned and are not retained.
  !! On success this records the active dimension and shift, marks the factors
  !! valid, and resets `updates_since_fresh`; the lifetime structural-update
  !! count is unchanged.
  !!
  !! Input parameters:
  !!   h     - Square symmetric Hamiltonian matrix. It is never modified.
  !!   s     - Square symmetric overlap matrix of the same order as h. It is
  !!           never modified.
  !!   shift - Approximate eigenvalue represented by the new factors.
  !!
  !! Output parameter:
  !!   info  - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, or
  !!           QR_ERR_FACTORIZATION.
  !!
  !! A failed argument check preserves any existing valid factors. Once the
  !! state buffers are overwritten, however, a LAPACK failure leaves the state
  !! explicitly invalid rather than exposing partial factors.
  subroutine real_factorize_fresh(self, h, s, shift, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer :: i, j, lapack_info, matrix_n

    ! Validate everything before touching an existing factorization.
    info = QR_ERR_INVALID_ARGUMENT
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
    if (size(h, 1) /= size(h, 2)) return
    if (size(s, 1) /= size(s, 2)) return
    if (size(h, 1) /= size(s, 1)) return

    matrix_n = size(h, 1)
    if (matrix_n <= 0 .or. matrix_n > self%capacity) return

    ! From this point onward the old factorization cannot be recovered.
    self%valid = .false.
    self%n = 0

    ! Form M = H - shift*S directly in Q. The inner index is the first
    ! Fortran array index, so each column is traversed contiguously. An
    ! explicit loop also prevents creation of a full temporary M matrix.
    do j = 1, matrix_n
      do i = 1, matrix_n
        self%q(i, j) = h(i, j) - shift * s(i, j)
      end do
    end do

    ! DGEQRF overwrites Q with compact Householder QR storage: R occupies the
    ! upper triangle and reflector vectors occupy the strict lower triangle.
    call dgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    ! Preserve R before DORGQR overwrites the compact reflectors with explicit
    ! Q. Explicitly clear R below the diagonal so consumers never observe
    ! stale values left by an earlier, larger factorization.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = 0.0_wp
      end do
    end do

    ! Expand the Householder representation in Q into an ordinary dense,
    ! explicit orthogonal matrix, as required by the qrupdate backend.
    call dorgqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    ! Commit metadata only after both LAPACK operations have succeeded.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine real_factorize_fresh

  !> Build fresh complex factors of the Hermitian matrix `H - shift*S`.
  !!
  !! The shifted matrix is formed directly in Q and passed through `ZGEQRF`
  !! and `ZUNGQR`, producing explicit unitary Q and upper-triangular R. The
  !! ownership and state transitions match `real_factorize_fresh`.
  !!
  !! Input parameters:
  !!   h     - Square Hermitian Hamiltonian matrix. It is never modified.
  !!   s     - Square Hermitian overlap matrix of the same order as h. It is
  !!           never modified.
  !!   shift - Real approximate eigenvalue represented by the factors.
  !!
  !! Output parameter:
  !!   info  - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, or
  !!           QR_ERR_FACTORIZATION.
  subroutine complex_factorize_fresh(self, h, s, shift, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer :: i, j, lapack_info, matrix_n

    ! Validate dimensions and initialization before invalidating old factors.
    info = QR_ERR_INVALID_ARGUMENT
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
    if (size(h, 1) /= size(h, 2)) return
    if (size(s, 1) /= size(s, 2)) return
    if (size(h, 1) /= size(s, 1)) return

    matrix_n = size(h, 1)
    if (matrix_n <= 0 .or. matrix_n > self%capacity) return

    self%valid = .false.
    self%n = 0

    ! Form the shifted complex matrix directly in the Q buffer. The shift is
    ! explicitly converted to complex(wp) to keep all arithmetic kind-correct.
    do j = 1, matrix_n
      do i = 1, matrix_n
        self%q(i, j) = h(i, j) - &
                       cmplx(shift, 0.0_wp, kind=wp) * s(i, j)
      end do
    end do

    ! ZGEQRF returns R above the diagonal and unitary Householder reflectors
    ! below it, with their scalar coefficients in tau.
    call zgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    ! Copy R before the compact reflector storage is expanded into Q.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      end do
    end do

    ! Generate the explicit unitary factor required by complex qrupdate.
    call zungqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    ! Commit the new factorization atomically at the metadata level.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine complex_factorize_fresh

  !> Replace real symmetric row and column `idx` using physical H and S
  !! changes supplied by the caller.
  !!
  !! The future update uses `d = delta_h - shift*delta_s` and two `qr1up`
  !! operations, with the second update's `idx` component removed so the
  !! diagonal is counted once. Caller arrays will first be copied to reusable
  !! state workspace because qrupdate may modify them. Version 0.1 performs no
  !! update, changes no counters, and returns `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameters:
  !!   idx     - One-based index of the row and column being replaced.
  !!   delta_h - Physical change in H(:,idx), including its diagonal entry.
  !!   delta_s - Physical change in S(:,idx), including its diagonal entry.
  !!
  !! Output parameter:
  !!   info    - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine real_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_replace_symmetric

  !> Replace complex Hermitian row and column `idx` using physical H and S
  !! changes supplied by the caller.
  !!
  !! The future implementation uses two conjugate rank-one QR updates,
  !! preserves Hermitian symmetry, and rejects a non-real diagonal change
  !! outside a precision-scaled tolerance. Version 0.1 leaves the state and
  !! counters unchanged and returns `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameters:
  !!   idx     - One-based index of the Hermitian row and column to replace.
  !!   delta_h - Change in H(:,idx); its diagonal component must be real to
  !!             precision-scaled tolerance in the future implementation.
  !!   delta_s - Corresponding change in S(:,idx).
  !!
  !! Output parameter:
  !!   info    - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine complex_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_replace_symmetric

  !> Append a real symmetric row and column at active index `n+1`.
  !!
  !! A future implementation will combine `h_column` and `s_column` with the
  !! stored shift, then use `qrinc` for the new column and `qrinr` for its
  !! matching row without allocating. Version 0.1 leaves the dimension and
  !! counters unchanged and returns `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameters:
  !!   h_column - New H column, including the new diagonal element at n+1.
  !!   s_column - New S column, including the new diagonal element at n+1.
  !!
  !! Output parameter:
  !!   info     - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine real_append_symmetric(self, h_column, s_column, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_append_symmetric

  !> Append a complex Hermitian row and column at active index `n+1`.
  !!
  !! The future operation will preserve conjugate symmetry while applying
  !! `qrinc` followed by `qrinr` from preallocated workspace. Version 0.1
  !! leaves the state unchanged and returns `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameters:
  !!   h_column - New H column, including the real diagonal value at n+1.
  !!   s_column - New S column, including the real diagonal value at n+1.
  !!
  !! Output parameter:
  !!   info     - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine complex_append_symmetric(self, h_column, s_column, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_append_symmetric

  !> Delete real symmetric row and column `idx` from the active factors.
  !!
  !! A future implementation will apply `qrdec` to the selected column and
  !! `qrder` to the matching row, then reduce the active dimension. Version
  !! 0.1 leaves the state and counters unchanged and returns
  !! `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameter:
  !!   idx  - One-based active row and column index to delete.
  !!
  !! Output parameter:
  !!   info - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine real_delete_symmetric(self, idx, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_delete_symmetric

  !> Delete complex Hermitian row and column `idx` from the active factors.
  !!
  !! The future implementation will use `qrdec` and `qrder` while preserving
  !! the remaining unitary/triangular factors. Version 0.1 leaves the state
  !! and counters unchanged and returns `QR_ERR_NOT_IMPLEMENTED`.
  !!
  !! Input parameter:
  !!   idx  - One-based active Hermitian row and column index to delete.
  !!
  !! Output parameter:
  !!   info - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  subroutine complex_delete_symmetric(self, idx, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_delete_symmetric

  !> Perform real generalized inverse iteration using the stored QR factors.
  !!
  !! This is the GSEPIIS inverse iteration expressed through the state's QR
  !! factorization instead of the original LDL^T factorization. Every step
  !! computes `S*v`, applies `transpose(Q)`, and solves `R*x=y`, which is
  !! algebraically `(H-shift*S)*x=S*v`. The convergence estimate and stopping
  !! rules intentionally match GSEPIIS.
  !!
  !! Input parameters:
  !!   s         - Caller-owned symmetric overlap matrix.
  !!   v_initial - Nonzero starting vector; it is not modified.
  !!   tol       - Requested relative convergence tolerance.
  !!   max_iter  - Maximum allowed inverse iterations.
  !!   norm_mode - Zero selects `x^T*S*x=1`; one selects `x^T*x=1`; every
  !!               other value retains the iteration scaling in which the
  !!               largest absolute component is one.
  !!
  !! Output parameters:
  !!   x         - Computed eigenvector in the requested normalization.
  !!   lambda    - Shift plus the Rayleigh quotient of the shifted matrix.
  !!   rel_acc   - Direction-change estimate from the final iteration.
  !!   num_iter  - Number of completed inverse iterations.
  !!   info      - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, QR_ERR_SINGULAR, or
  !!               QR_ERR_NO_CONVERGENCE. The last status still returns the
  !!               best eigenpair produced within max_iter.
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
    integer :: matrix_n
    real(wp) :: coefficient, current_norm_squared, eigenvector_norm_squared
    real(wp) :: max_component, norm_of_diff, norm_of_diff_previous
    real(wp) :: overlap_norm_squared, shifted_numerator
    logical :: not_converged

    x = 0.0_wp
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    info = QR_ERR_INVALID_ARGUMENT

    ! Validate the complete state and caller dimensions before using a BLAS
    ! kernel. S is read through its lower triangle, exactly as in GSEPIIS.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) /= matrix_n .or. size(s, 2) /= matrix_n) return
    if (size(v_initial) /= matrix_n .or. size(x) /= matrix_n) return
    if (max_iter <= 0) return
    if (real_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) return

    ! Unlike LDL^T factorization, xGEQRF does not report a singular matrix.
    ! Detect an unusable diagonal before DTRSV can divide by it.
    if (real_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    ! The first half of solve_work is the current iterate v. The second half
    ! receives Q^T*S*v. Both buffers were allocated by initialize, so solve
    ! performs no allocation regardless of the number of iterations.
    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    do while (not_converged .and. num_iter < max_iter)
      call dsymv('L', matrix_n, 1.0_wp, s, matrix_n, &
                 self%solve_work(1:matrix_n), 1, 0.0_wp, x, 1)
      call dgemv('T', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
                 x, 1, 0.0_wp, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      call dtrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      ! GSEPIIS scales every iterate by its largest absolute component. This
      ! protects the repeated inverse solves from overflow without changing
      ! the represented direction.
      max_component = maxval(abs(x))
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / max_component

      ! Remove the component of x parallel to the previous iterate and use
      ! the relative Euclidean norm of the remainder as the convergence
      ! estimate. The coefficient is (x^T*v)/(v^T*v), as in GSEPIIS.
      current_norm_squared = real_norm_squared( &
                               matrix_n, self%solve_work(1:matrix_n))
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = dot_product(x, self%solve_work(1:matrix_n)) / &
                    current_norm_squared
      eigenvector_norm_squared = real_norm_squared(matrix_n, x)
      norm_of_diff = real_direction_difference( &
                       matrix_n, x, coefficient, &
                       self%solve_work(1:matrix_n), &
                       eigenvector_norm_squared)

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        ! A negative tolerance asks for the most accurate result available,
        ! but never accepts less accuracy than abs(tol). Requiring the error
        ! to turn upward forces at least one confirming iteration.
        if (norm_of_diff > norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        self%solve_work(1:matrix_n) = x
      end if
    end do

    rel_acc = norm_of_diff
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    ! Compute x^T*S*x before reusing the first workspace vector. S is assumed
    ! positive definite by the generalized symmetric eigenproblem.
    call dsymv('L', matrix_n, 1.0_wp, s, matrix_n, x, 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = dot_product(x, &
                                       self%solve_work(1:matrix_n))
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    ! The state does not retain M=H-shift*S. Reconstruct only its action on x:
    ! M*x = Q*(R*x), then add the stored shift to the Rayleigh quotient.
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%r, self%capacity, &
               x, 1, 0.0_wp, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1)
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    shifted_numerator = dot_product(x, self%solve_work(1:matrix_n))
    lambda = self%shift + shifted_numerator / overlap_norm_squared

    select case (norm_mode)
    case (0)
      x = x / sqrt(overlap_norm_squared)
    case (1)
      eigenvector_norm_squared = real_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / sqrt(eigenvector_norm_squared)
    end select
  end subroutine real_solve

  !> Perform complex generalized inverse iteration using the stored QR
  !! factors.
  !!
  !! This is the Hermitian GHEPIIS iteration, with `Q^H` and a complex
  !! triangular solve replacing the original LDL^H solve. Complex phases,
  !! the convergence estimate, and the normalization rules follow GHEPIIS.
  !!
  !! Parameters have the same roles as in real_solve, with complex Hermitian
  !! S and complex starting/output vectors. The eigenvalue, tolerance, and
  !! relative-accuracy estimate remain real(wp).
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
    integer :: matrix_n
    real(wp) :: current_norm_squared, eigenvector_norm_squared
    real(wp) :: max_component, norm_of_diff, norm_of_diff_previous
    real(wp) :: overlap_norm_squared, shifted_numerator
    complex(wp) :: coefficient
    complex(wp), parameter :: complex_one = (1.0_wp, 0.0_wp)
    complex(wp), parameter :: complex_zero = (0.0_wp, 0.0_wp)
    logical :: not_converged

    x = cmplx(0.0_wp, 0.0_wp, kind=wp)
    lambda = 0.0_wp
    rel_acc = huge(1.0_wp)
    num_iter = 0
    info = QR_ERR_INVALID_ARGUMENT

    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) /= matrix_n .or. size(s, 2) /= matrix_n) return
    if (size(v_initial) /= matrix_n .or. size(x) /= matrix_n) return
    if (max_iter <= 0) return
    if (complex_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) return

    if (complex_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    do while (not_converged .and. num_iter < max_iter)
      call zhemv('L', matrix_n, complex_one, s, matrix_n, &
                 self%solve_work(1:matrix_n), 1, complex_zero, x, 1)
      call zgemv('C', matrix_n, matrix_n, complex_one, self%q, &
                 self%capacity, x, 1, complex_zero, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      call ztrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      ! GHEPIIS uses the largest real or imaginary component, rather than the
      ! largest complex modulus, so retain that detail here.
      max_component = complex_max_abs_real_or_imag(matrix_n, x)
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(max_component, 0.0_wp, kind=wp)

      current_norm_squared = complex_norm_squared( &
                               matrix_n, self%solve_work(1:matrix_n))
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = dot_product(x, self%solve_work(1:matrix_n)) / &
                    cmplx(current_norm_squared, 0.0_wp, kind=wp)
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      norm_of_diff = complex_direction_difference( &
                       matrix_n, x, coefficient, &
                       self%solve_work(1:matrix_n), &
                       eigenvector_norm_squared)

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        if (norm_of_diff > norm_of_diff_previous .and. &
            norm_of_diff <= abs(tol)) not_converged = .false.
        norm_of_diff_previous = norm_of_diff
      end if

      num_iter = num_iter + 1
      if (not_converged .and. num_iter < max_iter) then
        self%solve_work(1:matrix_n) = x
      end if
    end do

    rel_acc = norm_of_diff
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    call zhemv('L', matrix_n, complex_one, s, matrix_n, x, 1, &
               complex_zero, self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = real( &
                               dot_product(x, &
                                 self%solve_work(1:matrix_n)), wp)
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    call zgemv('N', matrix_n, matrix_n, complex_one, self%r, &
               self%capacity, x, 1, complex_zero, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1)
    call zgemv('N', matrix_n, matrix_n, complex_one, self%q, &
               self%capacity, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1, &
               complex_zero, self%solve_work(1:matrix_n), 1)
    shifted_numerator = real( &
                            dot_product(x, &
                              self%solve_work(1:matrix_n)), wp)
    lambda = self%shift + shifted_numerator / overlap_norm_squared

    select case (norm_mode)
    case (0)
      x = x / cmplx(sqrt(overlap_norm_squared), 0.0_wp, kind=wp)
    case (1)
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(sqrt(eigenvector_norm_squared), 0.0_wp, kind=wp)
    end select
  end subroutine complex_solve

  !> Return the squared Euclidean norm of a real vector.
  !!
  !! This deliberately mirrors RDotProdItself from the pristine linalg
  !! implementation. Keeping the reduction here avoids an intermediate array
  !! and keeps the helper generic in wp.
  function real_norm_squared(n, x) result(norm_squared)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(:)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + x(i) * x(i)
    end do
  end function real_norm_squared

  !> Return the real squared Euclidean norm x^H*x of a complex vector.
  function complex_norm_squared(n, x) result(norm_squared)
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(:)
    real(wp) :: norm_squared
    integer :: i

    norm_squared = 0.0_wp
    do i = 1, n
      norm_squared = norm_squared + real(x(i), wp)**2 + aimag(x(i))**2
    end do
  end function complex_norm_squared

  !> Compute ||x-alpha*y||_2 / ||x||_2 for real inverse iteration.
  function real_direction_difference(n, x, alpha, y, x_norm_squared) &
      result(relative_difference)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(:), alpha, y(:), x_norm_squared
    real(wp) :: relative_difference
    real(wp) :: difference, difference_norm_squared
    integer :: i

    difference_norm_squared = 0.0_wp
    do i = 1, n
      difference = x(i) - alpha * y(i)
      difference_norm_squared = difference_norm_squared + &
                                difference * difference
    end do
    relative_difference = sqrt(difference_norm_squared / x_norm_squared)
  end function real_direction_difference

  !> Compute ||x-alpha*y||_2 / ||x||_2 for complex inverse iteration.
  function complex_direction_difference(n, x, alpha, y, x_norm_squared) &
      result(relative_difference)
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(:), alpha, y(:)
    real(wp), intent(in) :: x_norm_squared
    real(wp) :: relative_difference, difference_norm_squared
    complex(wp) :: difference
    integer :: i

    difference_norm_squared = 0.0_wp
    do i = 1, n
      difference = x(i) - alpha * y(i)
      difference_norm_squared = difference_norm_squared + abs(difference)**2
    end do
    relative_difference = sqrt(difference_norm_squared / x_norm_squared)
  end function complex_direction_difference

  !> Return the scale used by GHEPIIS: the largest magnitude among every
  !! real and imaginary component, not the largest complex modulus.
  function complex_max_abs_real_or_imag(n, x) result(max_component)
    integer, intent(in) :: n
    complex(wp), intent(in) :: x(:)
    real(wp) :: max_component
    integer :: i

    max_component = 0.0_wp
    do i = 1, n
      max_component = max(max_component, abs(real(x(i), wp)), &
                          abs(aimag(x(i))))
    end do
  end function complex_max_abs_real_or_imag

  !> Detect a zero or precision-scale-small diagonal in a real R factor.
  !! DGEQRF itself reports only invalid arguments, so this check supplies the
  !! singular-matrix error that the original LDL^T factorization provided.
  function real_upper_factor_is_singular(r, n) result(is_singular)
    real(wp), intent(in) :: r(:,:)
    integer, intent(in) :: n
    logical :: is_singular
    real(wp) :: factor_scale, threshold
    integer :: i

    factor_scale = maxval(abs(r(1:n,1:n)))
    threshold = max(tiny(1.0_wp), epsilon(1.0_wp) * factor_scale)
    is_singular = .false.
    do i = 1, n
      if (abs(r(i,i)) <= threshold) then
        is_singular = .true.
        return
      end if
    end do
  end function real_upper_factor_is_singular

  !> Complex counterpart of real_upper_factor_is_singular.
  function complex_upper_factor_is_singular(r, n) result(is_singular)
    complex(wp), intent(in) :: r(:,:)
    integer, intent(in) :: n
    logical :: is_singular
    real(wp) :: factor_scale, threshold
    integer :: i

    factor_scale = maxval(abs(r(1:n,1:n)))
    threshold = max(tiny(1.0_wp), epsilon(1.0_wp) * factor_scale)
    is_singular = .false.
    do i = 1, n
      if (abs(r(i,i)) <= threshold) then
        is_singular = .true.
        return
      end if
    end do
  end function complex_upper_factor_is_singular

end module qrlinalg

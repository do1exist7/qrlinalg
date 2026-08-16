!Module qrlinalg provides serial QR factorizations and inverse iteration for
!the generalized symmetric and Hermitian eigenvalue problems
!
!                 H*x = lambda*S*x .
!
!A real state represents a symmetric problem and a complex state represents a
!Hermitian problem. For a real shift sigma, both states store the factorization
!
!                 H - sigma*S = Q*R ,
!
!where Q is explicit and orthogonal or unitary and R is upper triangular. The
!explicit representation is required by the structural QR-update routines and
!also permits inverse iteration without retaining a separate shifted matrix.
!
!The caller owns H, S, and all vectors passed through the public interface. A
!state does not retain pointers to caller arrays. It owns Q, R, the Householder
!coefficients used while constructing Q, and all factorization, update, and
!solve workspace. After initialization, factorization and inverse iteration do
!not allocate memory.
!
!The working kind wp is selected when the library is compiled. The bundled
!BLAS, LAPACK, qrupdate sources, interfaces in this module, and public arrays
!all use that same kind. No MPI object or process-dependent state is stored.
module qrlinalg
  use iso_fortran_env, only: int64
  use wp_def, only: wp
  use qrupdate, only: qr1up
  implicit none
  private

  integer, parameter, public :: QR_SUCCESS = 0
    !The requested operation completed successfully.
  integer, parameter, public :: QR_ERR_INVALID_ARGUMENT = 1
    !An argument, state, dimension, or capacity is invalid.
  integer, parameter, public :: QR_ERR_ALLOCATION = 2
    !Initialization could not allocate all required state storage.
  integer, parameter, public :: QR_ERR_NOT_IMPLEMENTED = 3
    !The requested numerical operation is not implemented in this version.
  integer, parameter, public :: QR_ERR_FACTORIZATION = 4
    !A bundled LAPACK factorization or workspace query reported an error.
  integer, parameter, public :: QR_ERR_SINGULAR = 5
    !The stored shifted factorization is singular or numerically unusable.
  integer, parameter, public :: QR_ERR_NO_CONVERGENCE = 6
    !Inverse iteration reached max_iter before satisfying the tolerance.

  public :: wp
  public :: qr_real_state
  public :: qr_complex_state

  !Interfaces to the external BLAS and LAPACK routines supplied in
  !src/qrupdate. Explicit interfaces allow the compiler to verify argument
  !types and working kinds at every call site. The routine names retain their
  !traditional D/Z prefixes, but their scalar kind is the compile-time wp;
  !there is no assumption that a D routine always uses eight-byte reals.
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
    !The components are private in a normal library build. QRLINALG_TESTING is
    !defined only for the white-box test build, where direct access is needed
    !to verify factor and workspace invariants without adding copying accessors
    !to the public numerical interface.
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
    !The component-access policy and metadata definitions are identical to
    !qr_real_state. Matrix and vector storage is complex; shift and the
    !rotation workspace used by complex qrupdate kernels remain real.
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

  !Subroutine clear_real_state releases every allocation owned by a real QR
  !state and restores the state to its default, uninitialized condition. It is
  !used before reinitialization and after a partial allocation failure. Calling
  !the routine for an already empty state is valid.
  !
  !  Input/output parameter:
  !    self - The real QR state. On exit all allocatable components are
  !           unallocated; n and capacity are zero; valid is false; shift and
  !           both structural-update counters are zero.
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

  !Subroutine clear_complex_state releases every allocation owned by a complex
  !QR state and restores the state to its default, uninitialized condition. In
  !addition to the complex arrays it releases the real workspace required by
  !complex plane rotations. Calling the routine for an empty state is valid.
  !
  !  Input/output parameter:
  !    self - The complex QR state. On exit all allocatable components are
  !           unallocated and all metadata and counters have their default
  !           values.
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

  !Subroutine real_initialize prepares a real QR state for symmetric matrices
  !of order not greater than max_n. Any factorization and counters previously
  !held by the state are discarded.
  !
  !The following storage is allocated:
  !  - full max_n by max_n arrays for the explicit factors Q and R;
  !  - max_n Householder coefficients in tau;
  !  - 4*max_n elements for structural QR updates;
  !  - 2*max_n elements for inverse iteration and Rayleigh quotients;
  !  - one factorization workspace large enough for both DGEQRF and DORGQR.
  !
  !The size of factor_work is obtained by querying DGEQRF and DORGQR with
  !LWORK=-1 for a square matrix of order max_n. The larger recommended size is
  !allocated once and reused. Initialization does not form a QR factorization;
  !therefore n is zero and valid is false on successful exit.
  !
  !  Input parameter:
  !    max_n - Maximum matrix order supported by the state. It must be
  !            positive.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Existing allocations and factorization
  !            metadata are destroyed before max_n is validated.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when all storage is ready;
  !            QR_ERR_INVALID_ARGUMENT when max_n is not positive;
  !            QR_ERR_ALLOCATION when an allocation fails;
  !            QR_ERR_FACTORIZATION when a LAPACK workspace query fails.
  !
  !If initialization fails, self is returned in the empty state described for
  !clear_real_state; no partial allocation remains owned by the object.
  subroutine real_initialize(self, max_n, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    real(wp) :: work_query(1)

    !Discard all previous storage before constructing the new state.
    call clear_real_state(self)
    if (max_n <= 0) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    !Allocate every array whose extent follows directly from max_n.
    !factor_work is allocated after the two workspace queries below.
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

    !Query the workspace recommended for the compact Householder
    !factorization. LWORK=-1 performs no factorization and returns the
    !recommendation in WORK(1).
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

    !Only capacity is committed here. The active order, represented shift, and
    !validity flag are committed by a successful factorize_fresh operation.
    self%capacity = max_n
    info = QR_SUCCESS
  end subroutine real_initialize

  !Subroutine complex_initialize prepares a complex QR state for Hermitian
  !matrices of order not greater than max_n. It follows the allocation and
  !failure semantics of real_initialize, with the following differences:
  !  - Q, R, tau, factor_work, update_work, and solve_work are complex(wp);
  !  - real_work contains max_n real(wp) elements required by the complex
  !    qrupdate rotation kernels;
  !  - ZGEQRF and ZUNGQR supply the two factorization-work recommendations.
  !
  !  Input parameter:
  !    max_n - Maximum matrix order supported by the state; must be positive.
  !
  !  Input/output parameter:
  !    self  - The state to initialize. Its old factors, work arrays, metadata,
  !            and counters are discarded.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, QR_ERR_ALLOCATION, or
  !            QR_ERR_FACTORIZATION, with the meanings documented for
  !            real_initialize.
  !
  !On failure self is empty. On success capacity=max_n, n=0, and valid=false.
  subroutine complex_initialize(self, max_n, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: max_n
    integer, intent(out) :: info
    integer :: allocation_status, factor_lwork, generate_q_lwork
    integer :: lapack_info, optimal_lwork
    complex(wp) :: work_query(1)

    !Discard all previous complex and real state storage.
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

    !In a complex LAPACK workspace query the recommended integer workspace
    !length is returned in the real part of WORK(1). Query both stages and use
    !the larger recommendation.
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

  !Subroutine real_factorize_fresh constructs a complete QR factorization of
  !the shifted real symmetric matrix
  !
  !                     M = H - shift*S = Q*R .
  !
  !The matrix M is formed directly in the state-owned Q array. DGEQRF replaces
  !M by the upper-triangular factor R and a compact set of Householder vectors.
  !The upper triangle is copied to the state-owned R array before DORGQR
  !replaces the compact representation by the explicit orthogonal matrix Q.
  !The strict lower triangle of R is explicitly set to zero.
  !
  !No temporary matrix of order n is created, and neither H nor S is modified
  !or retained. The state must have been initialized with capacity at least n.
  !H and S are required to be square, of equal order, and symmetric. Symmetry
  !is a mathematical precondition and is not checked element by element.
  !
  !  Input parameters:
  !    h     - The n by n real symmetric Hamiltonian matrix H.
  !    s     - The n by n real symmetric overlap matrix S.
  !    shift - The real shift represented by the factorization.
  !
  !  Input/output parameter:
  !    self  - An initialized real QR state. On successful exit:
  !              self%n = n,
  !              self%shift = shift,
  !              self%valid = true,
  !              self%Q*self%R = H-shift*S to working precision,
  !              self%updates_since_fresh = 0.
  !            structural_updates is a lifetime counter and is not reset.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS when both LAPACK stages succeed;
  !            QR_ERR_INVALID_ARGUMENT for an uninitialized state, a
  !            nonsquare or inconsistent matrix, an empty matrix, or n greater
  !            than capacity;
  !            QR_ERR_FACTORIZATION when DGEQRF or DORGQR reports an error.
  !
  !All argument checks are completed before existing factors are overwritten.
  !An argument error therefore preserves the previous factorization. A LAPACK
  !error occurs after the state buffers have been modified; in that case n is
  !returned as zero and valid is false so that partial factors cannot be used.
  subroutine real_factorize_fresh(self, h, s, shift, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer :: i, j, lapack_info, matrix_n

    !Check initialization, matrix shapes, active order, and capacity before
    !overwriting any component of a previously valid factorization.
    info = QR_ERR_INVALID_ARGUMENT
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
    if (size(h, 1) /= size(h, 2)) return
    if (size(s, 1) /= size(s, 2)) return
    if (size(h, 1) /= size(s, 1)) return

    matrix_n = size(h, 1)
    if (matrix_n <= 0 .or. matrix_n > self%capacity) return

    !The Q buffer is about to be overwritten. Mark the factors invalid until
    !both the factorization and explicit-Q generation have completed.
    self%valid = .false.
    self%n = 0

    !Form M directly in Q. The first array index is the inner loop so each
    !column is written contiguously. The explicit loops prevent creation of a
    !full array temporary for H-shift*S.
    do j = 1, matrix_n
      do i = 1, matrix_n
        self%q(i, j) = h(i, j) - shift * s(i, j)
      end do
    end do

    !Compute the compact Householder representation. On exit from DGEQRF the
    !upper triangle contains R; the strict lower triangle and tau describe Q.
    call dgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Copy R before DORGQR destroys the compact reflector representation. Clear
    !the strict lower triangle so all entries in the active R block have a
    !defined triangular meaning, including after refactorization at smaller n.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = 0.0_wp
      end do
    end do

    !Expand the Householder vectors and tau into the explicit orthogonal matrix
    !Q required by inverse iteration and the qrupdate kernels.
    call dorgqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Commit the active factorization only after every numerical stage succeeds.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine real_factorize_fresh

  !Subroutine complex_factorize_fresh constructs a complete QR factorization
  !of the shifted complex Hermitian matrix
  !
  !                     M = H - shift*S = Q*R ,
  !
  !where Q is unitary and R is upper triangular. ZGEQRF produces a compact
  !Householder representation, R is copied from its upper triangle, and ZUNGQR
  !generates explicit Q. The routine allocates no memory and does not modify or
  !retain H and S.
  !
  !  Input parameters:
  !    h     - The n by n complex Hermitian Hamiltonian matrix H.
  !    s     - The n by n complex Hermitian overlap matrix S.
  !    shift - The real shift represented by the factorization.
  !
  !  Input/output parameter:
  !    self  - An initialized complex QR state with capacity at least n. The
  !            successful state transitions and counter rules are identical to
  !            real_factorize_fresh, with unitary Q in place of orthogonal Q.
  !
  !  Output parameter:
  !    info  - QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, or QR_ERR_FACTORIZATION,
  !            under the conditions documented for real_factorize_fresh.
  !
  !Hermitian structure is a mathematical precondition and is not checked. An
  !invalid argument preserves existing factors; a failure after ZGEQRF begins
  !leaves n=0 and valid=false.
  subroutine complex_factorize_fresh(self, h, s, shift, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h(:,:), s(:,:)
    real(wp), intent(in) :: shift
    integer, intent(out) :: info
    integer :: i, j, lapack_info, matrix_n

    !Validate all state and dimension requirements before modifying Q or R.
    info = QR_ERR_INVALID_ARGUMENT
    if (self%capacity <= 0 .or. .not. allocated(self%q)) return
    if (size(h, 1) /= size(h, 2)) return
    if (size(s, 1) /= size(s, 2)) return
    if (size(h, 1) /= size(s, 1)) return

    matrix_n = size(h, 1)
    if (matrix_n <= 0 .or. matrix_n > self%capacity) return

    self%valid = .false.
    self%n = 0

    !Form M directly in Q. The explicit conversion of shift makes the complex
    !working kind independent of compiler rules for mixed-kind expressions.
    do j = 1, matrix_n
      do i = 1, matrix_n
        self%q(i, j) = h(i, j) - &
                       cmplx(shift, 0.0_wp, kind=wp) * s(i, j)
      end do
    end do

    !Compute R and the compact unitary Householder representation of Q.
    call zgeqrf(matrix_n, matrix_n, self%q, self%capacity, self%tau, &
                self%factor_work, size(self%factor_work), lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Copy the active upper triangle to R and define its lower triangle as zero
    !before the reflector storage is replaced by explicit Q.
    do j = 1, matrix_n
      do i = 1, j
        self%r(i, j) = self%q(i, j)
      end do
      do i = j + 1, matrix_n
        self%r(i, j) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      end do
    end do

    !Generate the explicit unitary matrix Q from the reflectors and tau.
    call zungqr(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
                self%tau, self%factor_work, size(self%factor_work), &
                lapack_info)
    if (lapack_info /= 0) then
      info = QR_ERR_FACTORIZATION
      return
    end if

    !Publish the active order, shift, and validity only after ZUNGQR succeeds.
    self%n = matrix_n
    self%shift = shift
    self%valid = .true.
    self%updates_since_fresh = 0_int64
    info = QR_SUCCESS
  end subroutine complex_factorize_fresh

  !Subroutine real_replace_symmetric replaces one row and the corresponding
  !column of a real symmetric problem by updating the stored QR factors. If
  !delta_h and delta_s denote the changes in the physical H and S columns, the
  !change represented by the QR state is
  !
  !                 d = delta_h - self%shift*delta_s .
  !
  !A symmetric row-and-column replacement can be expressed as
  !
  !                 d*e_idx^T + e_idx*(d-d(idx)*e_idx)^T .
  !
  !The subtraction of d(idx)*e_idx prevents the diagonal change from being
  !applied twice. The operation applies these two rank-one terms sequentially
  !with qr1up. qr1up overwrites its u and v arguments, so both calls use
  !state-owned copies and leave delta_h and delta_s unchanged.
  !
  !  Input parameters:
  !    idx     - One-based index of the replaced row and column.
  !    delta_h - Change in H(:,idx), including the diagonal element.
  !    delta_s - Change in S(:,idx), including the diagonal element.
  !
  !  Input/output parameter:
  !    self    - A valid QR state of active order n. On success Q and R
  !              represent the updated shifted matrix, valid, n, capacity, and
  !              shift are preserved, and both structural-update counters are
  !              incremented by one.
  !
  !  Output parameter:
  !    info    - QR_SUCCESS when both rank-one updates are applied;
  !              QR_ERR_INVALID_ARGUMENT when the state is invalid, idx is
  !              outside 1:n, either change vector has length other than n, or
  !              required state workspace is unavailable.
  !
  !All validation precedes modification of Q or R, so QR_ERR_INVALID_ARGUMENT
  !preserves the complete state. The validated qr1up calls have no numerical
  !failure return. No allocation is performed.
  subroutine real_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    real(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info
    integer :: i, matrix_n

    info = QR_ERR_INVALID_ARGUMENT
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    matrix_n = self%n
    if (idx < 1 .or. idx > matrix_n) return
    if (size(delta_h) /= matrix_n .or. size(delta_s) /= matrix_n) return
    if (size(self%update_work) < 4 * matrix_n) return

    !First apply d*e_idx^T. The update vectors occupy the first two workspace
    !blocks and qr1up uses the remaining two blocks as rotation workspace.
    do i = 1, matrix_n
      self%update_work(i) = delta_h(i) - self%shift * delta_s(i)
      self%update_work(matrix_n + i) = 0.0_wp
    end do
    self%update_work(matrix_n + idx) = 1.0_wp
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:4 * matrix_n))

    !Then apply e_idx*(d-d(idx)*e_idx)^T. Reconstruct d because qr1up is
    !permitted to overwrite both vectors supplied to the first call.
    do i = 1, matrix_n
      self%update_work(i) = 0.0_wp
      self%update_work(matrix_n + i) = delta_h(i) - &
                                           self%shift * delta_s(i)
    end do
    self%update_work(idx) = 1.0_wp
    self%update_work(matrix_n + idx) = 0.0_wp
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:4 * matrix_n))

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine real_replace_symmetric

  !Subroutine complex_replace_symmetric is the Hermitian counterpart of
  !real_replace_symmetric. For d=delta_h-self%shift*delta_s, the represented
  !Hermitian change is
  !
  !                 d*e_idx^H + e_idx*(d-d(idx)*e_idx)^H .
  !
  !The diagonal of a Hermitian matrix is real. The represented diagonal change
  !is therefore rejected when its imaginary part exceeds
  !
  !       100*epsilon(1.0_wp)*max(1,maxval(abs(d))).
  !
  !An accepted roundoff-sized imaginary part is discarded before either
  !rank-one update. qr1up interprets its complex update as u*v^H, so the two
  !terms reconstruct the conjugate row without modifying caller arrays.
  !
  !  Input parameters:
  !    idx     - One-based index of the replaced row and column.
  !    delta_h - Change in H(:,idx), including its nominally real diagonal.
  !    delta_s - Change in S(:,idx), including its nominally real diagonal.
  !
  !  Input/output parameter:
  !    self    - A valid complex QR state. Successful state and counter changes
  !              are the same as for real_replace_symmetric.
  !
  !  Output parameter:
  !    info    - QR_SUCCESS when the Hermitian replacement is complete;
  !              QR_ERR_INVALID_ARGUMENT under the real-routine validation
  !              conditions or when the represented diagonal change is not
  !              real within the tolerance above.
  !
  !Every rejection occurs before Q or R is modified. The routine allocates no
  !memory and increments each counter once, rather than once per rank-one term.
  subroutine complex_replace_symmetric(self, idx, delta_h, delta_s, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    complex(wp), intent(in) :: delta_h(:), delta_s(:)
    integer, intent(out) :: info
    complex(wp) :: change
    real(wp) :: diagonal_tolerance, update_scale
    integer :: i, matrix_n

    info = QR_ERR_INVALID_ARGUMENT
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%update_work)) return
    if (.not. allocated(self%real_work)) return
    matrix_n = self%n
    if (idx < 1 .or. idx > matrix_n) return
    if (size(delta_h) /= matrix_n .or. size(delta_s) /= matrix_n) return
    if (size(self%update_work) < 3 * matrix_n) return
    if (size(self%real_work) < matrix_n) return

    !Determine the scale and validate the diagonal before either qr1up call.
    update_scale = 1.0_wp
    do i = 1, matrix_n
      change = delta_h(i) - &
               cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
      update_scale = max(update_scale, abs(change))
    end do
    change = delta_h(idx) - &
             cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(idx)
    diagonal_tolerance = 100.0_wp * epsilon(1.0_wp) * update_scale
    if (abs(aimag(change)) > diagonal_tolerance) return

    !Apply d*e_idx^H. The first three complex workspace blocks contain u, v,
    !and the qr1up work vector; real_work stores the rotation cosines.
    do i = 1, matrix_n
      self%update_work(i) = delta_h(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
      self%update_work(matrix_n + i) = &
        cmplx(0.0_wp, 0.0_wp, kind=wp)
    end do
    self%update_work(idx) = &
      cmplx(real(self%update_work(idx), wp), 0.0_wp, kind=wp)
    self%update_work(matrix_n + idx) = &
      cmplx(1.0_wp, 0.0_wp, kind=wp)
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), &
               self%real_work(1:matrix_n))

    !Apply e_idx*(d-d(idx)*e_idx)^H. Re-form d after the destructive first
    !call and set its diagonal component to exact complex zero.
    do i = 1, matrix_n
      self%update_work(i) = cmplx(0.0_wp, 0.0_wp, kind=wp)
      self%update_work(matrix_n + i) = delta_h(i) - &
        cmplx(self%shift, 0.0_wp, kind=wp) * delta_s(i)
    end do
    self%update_work(idx) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    self%update_work(matrix_n + idx) = &
      cmplx(0.0_wp, 0.0_wp, kind=wp)
    call qr1up(matrix_n, matrix_n, matrix_n, self%q, self%capacity, &
               self%r, self%capacity, self%update_work(1:matrix_n), &
               self%update_work(matrix_n + 1:2 * matrix_n), &
               self%update_work(2 * matrix_n + 1:3 * matrix_n), &
               self%real_work(1:matrix_n))

    self%structural_updates = self%structural_updates + 1_int64
    self%updates_since_fresh = self%updates_since_fresh + 1_int64
    info = QR_SUCCESS
  end subroutine complex_replace_symmetric

  !Subroutine real_append_symmetric is the reserved interface for increasing
  !the active real symmetric problem from order n to n+1. h_column and
  !s_column contain the new physical columns through the new diagonal element.
  !The column represented by the shifted QR factorization is
  !h_column-self%shift*s_column. A complete implementation applies qrinc for
  !the new column and qrinr for the corresponding symmetric row, using only
  !the workspace allocated for the state.
  !
  !  Input parameters:
  !    h_column - New H column of length self%n+1.
  !    s_column - New S column of length self%n+1.
  !
  !  Input/output parameter:
  !    self     - The QR state whose active order is to be increased.
  !
  !  Output parameter:
  !    info     - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  !
  !The present stub leaves n, the factors, and both counters unchanged.
  subroutine real_append_symmetric(self, h_column, s_column, info)
    class(qr_real_state), intent(inout) :: self
    real(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_append_symmetric

  !Subroutine complex_append_symmetric is the Hermitian counterpart of
  !real_append_symmetric. The supplied columns determine the new column and,
  !by conjugation, the new row of H-self%shift*S. Their final elements are
  !required to be real to a working-precision tolerance. The numerical
  !path applies qrinc followed by qrinr without allocating memory.
  !
  !  Input parameters:
  !    h_column - New complex H column of length self%n+1, including the real
  !               diagonal element.
  !    s_column - New complex S column of length self%n+1, including the real
  !               diagonal element.
  !
  !  Input/output parameter:
  !    self     - The complex QR state whose order is to be increased.
  !
  !  Output parameter:
  !    info     - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  !
  !The present stub leaves the complete state unchanged.
  subroutine complex_append_symmetric(self, h_column, s_column, info)
    class(qr_complex_state), intent(inout) :: self
    complex(wp), intent(in) :: h_column(:), s_column(:)
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_append_symmetric

  !Subroutine real_delete_symmetric is the reserved interface for deleting row
  !idx and column idx from an active real symmetric factorization. A complete
  !implementation removes the selected column with qrdec, removes the matching
  !row with qrder, and decreases the active order after both operations
  !succeed.
  !
  !  Input parameter:
  !    idx  - One-based row and column index in the active range 1:self%n.
  !
  !  Input/output parameter:
  !    self - The QR state from which the row and column are to be removed.
  !
  !  Output parameter:
  !    info - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  !
  !The present stub leaves n, the factors, and both counters unchanged.
  subroutine real_delete_symmetric(self, idx, info)
    class(qr_real_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine real_delete_symmetric

  !Subroutine complex_delete_symmetric is the Hermitian counterpart of
  !real_delete_symmetric. The complete operation uses complex qrdec and qrder
  !kernels and preserves an explicit unitary Q and upper-triangular R for the
  !remaining principal submatrix.
  !
  !  Input parameter:
  !    idx  - One-based row and column index in the active range 1:self%n.
  !
  !  Input/output parameter:
  !    self - The complex QR state from which the row and column are removed.
  !
  !  Output parameter:
  !    info - QR_ERR_NOT_IMPLEMENTED in version 0.1.
  !
  !The present stub leaves the complete state unchanged.
  subroutine complex_delete_symmetric(self, idx, info)
    class(qr_complex_state), intent(inout) :: self
    integer, intent(in) :: idx
    integer, intent(out) :: info

    info = QR_ERR_NOT_IMPLEMENTED
  end subroutine complex_delete_symmetric

  !Subroutine real_solve finds one eigenvalue and its eigenvector for the real
  !generalized symmetric eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration. On entry self must contain valid factors
  !
  !                        H - shift*S = Q*R .
  !
  !The desired eigenvalue should be closer to shift than every other
  !eigenvalue. The convergence rate is governed principally by the ratio of
  !the distance from shift to the desired eigenvalue and the distance from
  !shift to the next closest eigenvalue. shift must not coincide with an
  !eigenvalue, because H-shift*S is then singular.
  !
  !A nondegenerate eigenvalue is assumed. If several eigenvalues are separated
  !only at the scale of working-precision roundoff, inverse iteration may
  !return a vector in their joint invariant subspace and may reach max_iter
  !without satisfying the requested directional accuracy. In that case the
  !best vector and Rayleigh quotient obtained are still returned with
  !QR_ERR_NO_CONVERGENCE.
  !
  !For a current vector v, one iteration solves
  !
  !             (H-shift*S)*x = S*v .
  !
  !Using the stored factors, the operation is performed as
  !
  !             w = S*v,
  !             y = Q^T*w,
  !             R*x = y.
  !
  !The new vector is divided by max(abs(x)). Its change of direction relative
  !to v is measured by
  !
  !             alpha   = (x^T*v)/(v^T*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0, iteration stops when rel_acc<=tol. For tol<=0, iteration
  !continues until rel_acc begins to increase and the current value is not
  !larger than abs(tol). The latter rule requests the most accurate attainable
  !direction subject to the floor abs(tol), and necessarily performs at least
  !one additional comparison iteration.
  !
  !After iteration, the eigenvalue is obtained from the shifted Rayleigh
  !quotient
  !
  ! lambda = shift + (x^T*(H-shift*S)*x)/(x^T*S*x).
  !
  !The state does not store H-shift*S separately. Its action on x is evaluated
  !as Q*(R*x). Only the lower triangle of S is referenced by DSYMV; the upper
  !triangle may be absent or contain unrelated values.
  !
  !  Input parameters:
  !    s         - The n by n real symmetric overlap matrix. The lower
  !                triangle, including the diagonal, must be defined. S is
  !                expected to be positive definite.
  !    v_initial - A nonzero starting approximation of length n. A vector with
  !                a substantial component in the desired eigendirection
  !                generally converges faster. The array is not modified.
  !    tol       - Directional convergence tolerance. Its sign selects the
  !                stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^T*S*x = 1;
  !                  1: x^T*x = 1;
  !                other: max(abs(x)) = 1.
  !
  !  Input/output parameter:
  !    self      - A valid real QR state of active order n. Q, R, n, shift,
  !                validity, and structural-update counters are unchanged.
  !                solve_work is overwritten and remains internal scratch.
  !
  !  Output parameters:
  !    x         - The final eigenvector approximation in the normalization
  !                selected by norm_mode. It has length n.
  !    lambda    - The Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - The directional difference from the final iteration. This
  !                is a convergence estimate, not a rigorously bounded error.
  !    num_iter  - Number of inverse iterations performed.
  !    info      - QR_SUCCESS when the stopping criterion is satisfied;
  !                QR_ERR_INVALID_ARGUMENT for an invalid state, dimension,
  !                iteration limit, starting vector, or non-positive x^T*S*x;
  !                QR_ERR_SINGULAR when R cannot be used safely or an
  !                iteration produces a numerically zero vector;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. In this
  !                case x, lambda, rel_acc, and num_iter describe the best
  !                approximation reached.
  !
  !No allocation is performed. On an error detected before the first
  !iteration, x and lambda are zero, num_iter is zero, and rel_acc is huge.
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

    !Verify the factorization state, array dimensions, iteration limit, and
    !starting-vector norm before entering a BLAS routine.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) /= matrix_n .or. size(s, 2) /= matrix_n) return
    if (size(v_initial) /= matrix_n .or. size(x) /= matrix_n) return
    if (max_iter <= 0) return
    if (real_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) return

    !DGEQRF can complete for a rank-deficient matrix. Test the diagonal of R
    !explicitly before DTRSV performs divisions during inverse iteration.
    if (real_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !Workspace layout during iteration:
    !  solve_work(1:n)       contains the current vector v;
    !  solve_work(n+1:2*n)   contains Q^T*S*v and triangular-solve scratch.
    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform inverse iterations until the direction criterion is satisfied or
    !the caller-supplied iteration limit is exhausted.
    do while (not_converged .and. num_iter < max_iter)
      !Form the right-hand side w=S*v from the lower triangle of S.
      call dsymv('L', matrix_n, 1.0_wp, s, matrix_n, &
                 self%solve_work(1:matrix_n), 1, 0.0_wp, x, 1)
      !Transform w by Q^T. The result is the right-hand side of R*x=y.
      call dgemv('T', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
                 x, 1, 0.0_wp, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      !Solve the upper-triangular system in place to obtain the new iterate.
      call dtrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      !Scale the solution so its largest absolute component is one. Scaling
      !does not change the eigendirection and limits growth in repeated solves.
      max_component = maxval(abs(x))
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / max_component

      !Compute the least-direction-change coefficient and the relative norm of
      !the component of x not parallel to the previous iterate v.
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
        !For a non-positive tolerance, accept only after the direction error
        !has passed its minimum and remains within abs(tol).
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
    !Failure to meet the stopping rule is nonfatal: the final iterate is still
    !used to calculate and normalize an eigenpair approximation.
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    !Compute x^T*S*x for the Rayleigh quotient and, for norm_mode=0, final
    !normalization. A non-positive result violates the positive-definite S
    !precondition or indicates unusable numerical data.
    call dsymv('L', matrix_n, 1.0_wp, s, matrix_n, x, 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = dot_product(x, &
                                       self%solve_work(1:matrix_n))
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    !Evaluate M*x as Q*(R*x), form x^T*M*x, and add the stored shift to the
    !quotient. This requires two matrix-vector products but no stored M matrix.
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%r, self%capacity, &
               x, 1, 0.0_wp, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1)
    call dgemv('N', matrix_n, matrix_n, 1.0_wp, self%q, self%capacity, &
               self%solve_work(matrix_n + 1:2 * matrix_n), 1, 0.0_wp, &
               self%solve_work(1:matrix_n), 1)
    shifted_numerator = dot_product(x, self%solve_work(1:matrix_n))
    lambda = self%shift + shifted_numerator / overlap_norm_squared

    !Apply the normalization requested by the caller. No scaling is required
    !for other norm_mode values because each inverse iterate already has unit
    !largest-component magnitude.
    select case (norm_mode)
    case (0) !Normalize so that x^T*S*x=1.
      x = x / sqrt(overlap_norm_squared)
    case (1) !Normalize so that x^T*x=1.
      eigenvector_norm_squared = real_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / sqrt(eigenvector_norm_squared)
    end select
  end subroutine real_solve

  !Subroutine complex_solve finds one eigenvalue and its eigenvector for the
  !complex generalized Hermitian eigenvalue problem
  !
  !                        H*x = lambda*S*x
  !
  !by shifted inverse iteration. On entry self must contain valid factors
  !
  !                        H - shift*S = Q*R ,
  !
  !where Q is unitary and R is upper triangular. The desired eigenvalue must be
  !closer to shift than the remaining eigenvalues for ordinary nondegenerate
  !inverse iteration to select it. shift must not equal an eigenvalue.
  !Pathologically close or degenerate eigenvalues may yield a vector in their
  !joint invariant subspace and QR_ERR_NO_CONVERGENCE rather than a uniquely
  !determined eigenvector.
  !
  !For a current vector v, the iteration solves
  !
  !             (H-shift*S)*x = S*v
  !
  !through the three operations
  !
  !             w = S*v,
  !             y = Q^H*w,
  !             R*x = y.
  !
  !The vector is divided by the largest magnitude among all of its real and
  !imaginary components. This is not the same as division by max(abs(x)); it
  !preserves the component-scaling convention of the Hermitian inverse
  !iteration interface. Directional convergence is estimated by
  !
  !             alpha   = (v^H*x)/(v^H*v),
  !             rel_acc = ||x-alpha*v||_2/||x||_2 .
  !
  !For tol>0 the first rel_acc<=tol terminates iteration. For tol<=0, the
  !routine waits until rel_acc begins to increase and is no greater than
  !abs(tol), thereby seeking the smallest attainable direction change subject
  !to the requested floor.
  !
  !The final real eigenvalue approximation is the shifted Hermitian Rayleigh
  !quotient
  !
  ! lambda = shift + real(x^H*(H-shift*S)*x)/(x^H*S*x).
  !
  !The product (H-shift*S)*x is evaluated as Q*(R*x). Only the lower triangle
  !of S is referenced by ZHEMV; diagonal elements of S are assumed real and S
  !is expected to be positive definite.
  !
  !  Input parameters:
  !    s         - The n by n complex Hermitian overlap matrix. Its lower
  !                triangle and real diagonal must be defined.
  !    v_initial - A nonzero complex starting approximation of length n. It is
  !                not modified.
  !    tol       - Real directional convergence tolerance. Its sign selects
  !                the stopping rule described above.
  !    max_iter  - Maximum number of inverse iterations; must be positive.
  !    norm_mode - Required normalization of x on exit:
  !                  0: x^H*S*x = 1;
  !                  1: x^H*x = 1;
  !                other: the largest magnitude among every real and imaginary
  !                       component of x is one.
  !
  !  Input/output parameter:
  !    self      - A valid complex QR state of active order n. Numerical
  !                factors and public metadata are unchanged; solve_work is
  !                overwritten as private scratch storage.
  !
  !  Output parameters:
  !    x         - Final complex eigenvector approximation in the requested
  !                normalization.
  !    lambda    - Real Rayleigh-quotient eigenvalue approximation.
  !    rel_acc   - Real direction-change estimate from the final iteration.
  !    num_iter  - Number of inverse iterations performed.
  !    info      - QR_SUCCESS when convergence is detected;
  !                QR_ERR_INVALID_ARGUMENT for an invalid state, dimensions,
  !                iteration limit, starting vector, or non-positive x^H*S*x;
  !                QR_ERR_SINGULAR for an unusable R or zero iterate;
  !                QR_ERR_NO_CONVERGENCE when max_iter is reached. The latter
  !                status still returns the final eigenpair approximation.
  !
  !No allocation is performed. Before-iteration errors return zero x and
  !lambda, zero num_iter, and huge rel_acc.
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

    !Verify all state, dimension, iteration-limit, and starting-vector
    !requirements before using the stored factors or calling BLAS.
    if (.not. self%valid .or. self%n <= 0) return
    if (.not. allocated(self%q) .or. .not. allocated(self%r)) return
    if (.not. allocated(self%solve_work)) return
    matrix_n = self%n
    if (size(s, 1) /= matrix_n .or. size(s, 2) /= matrix_n) return
    if (size(v_initial) /= matrix_n .or. size(x) /= matrix_n) return
    if (max_iter <= 0) return
    if (complex_norm_squared(matrix_n, v_initial) <= tiny(1.0_wp)) return

    !ZGEQRF does not use a positive INFO value to report rank deficiency. The
    !diagonal of R is therefore tested explicitly before ZTRSV is called.
    if (complex_upper_factor_is_singular(self%r, matrix_n)) then
      info = QR_ERR_SINGULAR
      return
    end if

    !Workspace layout during iteration:
    !  solve_work(1:n)       contains the current vector v;
    !  solve_work(n+1:2*n)   contains Q^H*S*v and solve scratch.
    self%solve_work(1:matrix_n) = v_initial
    norm_of_diff_previous = huge(1.0_wp)
    not_converged = .true.

    !Perform Hermitian inverse iterations until convergence or max_iter.
    do while (not_converged .and. num_iter < max_iter)
      !Form w=S*v from the stored lower triangle of the Hermitian matrix S.
      call zhemv('L', matrix_n, complex_one, s, matrix_n, &
                 self%solve_work(1:matrix_n), 1, complex_zero, x, 1)
      !Apply Q^H to w to obtain the right-hand side of R*x=y.
      call zgemv('C', matrix_n, matrix_n, complex_one, self%q, &
                 self%capacity, x, 1, complex_zero, &
                 self%solve_work(matrix_n + 1:2 * matrix_n), 1)
      x = self%solve_work(matrix_n + 1:2 * matrix_n)
      !Solve the complex upper-triangular system for the new iterate.
      call ztrsv('U', 'N', 'N', matrix_n, self%r, self%capacity, x, 1)

      !Scale with the largest real or imaginary component. This bounds both
      !parts of every element without changing the complex eigendirection.
      max_component = complex_max_abs_real_or_imag(matrix_n, x)
      if (max_component <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(max_component, 0.0_wp, kind=wp)

      !Compute alpha=(v^H*x)/(v^H*v) and the prescribed relative difference
      !between x and alpha*v. Conjugating the previous iterate, rather than
      !the new iterate, makes this comparison invariant under the arbitrary
      !complex phase of an eigenvector.
      current_norm_squared = complex_norm_squared( &
                               matrix_n, self%solve_work(1:matrix_n))
      if (current_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      coefficient = dot_product(self%solve_work(1:matrix_n), x) / &
                    cmplx(current_norm_squared, 0.0_wp, kind=wp)
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      norm_of_diff = complex_direction_difference( &
                       matrix_n, x, coefficient, &
                       self%solve_work(1:matrix_n), &
                       eigenvector_norm_squared)

      if (tol > 0.0_wp) then
        if (norm_of_diff <= tol) not_converged = .false.
      else
        !A non-positive tolerance accepts a result only after the direction
        !error turns upward while remaining within abs(tol).
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
    !A nonconverged final iterate remains a valid approximation and is carried
    !through the Rayleigh-quotient and normalization calculations below.
    if (not_converged) then
      info = QR_ERR_NO_CONVERGENCE
    else
      info = QR_SUCCESS
    end if

    !Compute the Hermitian quadratic form x^H*S*x. Its real value is used both
    !in the Rayleigh quotient and, for norm_mode=0, in final normalization.
    call zhemv('L', matrix_n, complex_one, s, matrix_n, x, 1, &
               complex_zero, self%solve_work(1:matrix_n), 1)
    overlap_norm_squared = real( &
                               dot_product(x, &
                                 self%solve_work(1:matrix_n)), wp)
    if (overlap_norm_squared <= tiny(1.0_wp)) then
      info = QR_ERR_INVALID_ARGUMENT
      return
    end if

    !Evaluate M*x=Q*(R*x), form real(x^H*M*x), and add the stored shift. The
    !imaginary roundoff part of the Hermitian quadratic form is discarded.
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

    !Apply S or Euclidean normalization when requested. For every other mode,
    !retain the real/imaginary component scaling established in the iteration.
    select case (norm_mode)
    case (0) !Normalize so that x^H*S*x=1.
      x = x / cmplx(sqrt(overlap_norm_squared), 0.0_wp, kind=wp)
    case (1) !Normalize so that x^H*x=1.
      eigenvector_norm_squared = complex_norm_squared(matrix_n, x)
      if (eigenvector_norm_squared <= tiny(1.0_wp)) then
        info = QR_ERR_SINGULAR
        return
      end if
      x = x / cmplx(sqrt(eigenvector_norm_squared), 0.0_wp, kind=wp)
    end select
  end subroutine complex_solve

  !Function real_norm_squared computes the real Euclidean inner product
  !
  !                         x^T*x = sum(x(i)^2)
  !
  !over the first n elements of x. The result is the squared norm; the square
  !root is not taken. The reduction is performed in wp and does not allocate an
  !array temporary.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Real vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - x^T*x for x(1:n).
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

  !Function complex_norm_squared computes the Hermitian Euclidean inner
  !product
  !
  !              x^H*x = sum(real(x(i))^2+imag(x(i))^2)
  !
  !over the first n elements of a complex vector. The mathematically real
  !quantity is accumulated directly in real(wp), avoiding a complex reduction
  !and discarding no computed imaginary part.
  !
  !  Input parameters:
  !    n - Number of vector elements included in the reduction.
  !    x - Complex vector containing at least n elements.
  !
  !  Result:
  !    norm_squared - The real value x^H*x for x(1:n).
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

  !Function real_direction_difference computes the scale-independent change of
  !direction between two real vectors,
  !
  !                 ||x-alpha*y||_2 / ||x||_2 .
  !
  !In inverse iteration alpha=(x^T*y)/(y^T*y), so alpha*y is the component of
  !x parallel to the previous iterate y. The result measures only the remaining
  !directional change and is unaffected by real rescaling of an eigenvector.
  !
  !  Input parameters:
  !    n              - Number of elements included in the calculation.
  !    x              - New real iterate, containing at least n elements.
  !    alpha          - Scalar projection coefficient multiplying y.
  !    y              - Previous real iterate, containing at least n elements.
  !    x_norm_squared - Precomputed positive value x^T*x.
  !
  !  Result:
  !    relative_difference - Relative Euclidean norm shown above.
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

  !Function complex_direction_difference computes the convergence measure
  !used to compare two successive complex inverse iterates,
  !
  !                 ||x-alpha*y||_2 / ||x||_2 .
  !
  !For Hermitian inverse iteration alpha=(y^H*x)/(y^H*y). Each squared
  !magnitude abs(x(i)-alpha*y(i))^2 is accumulated in real(wp).
  !
  !  Input parameters:
  !    n              - Number of elements included in the calculation.
  !    x              - New complex iterate, containing at least n elements.
  !    alpha          - Complex projection coefficient multiplying y.
  !    y              - Previous complex iterate, containing at least n
  !                     elements.
  !    x_norm_squared - Precomputed positive real value x^H*x.
  !
  !  Result:
  !    relative_difference - Relative Euclidean norm shown above.
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

  !Function complex_max_abs_real_or_imag returns
  !
  !       max_i( max(abs(real(x(i))),abs(imag(x(i)))) )
  !
  !for the first n elements of x. This component norm is used to scale complex
  !inverse iterates. It differs from maxval(abs(x)), which uses the Euclidean
  !modulus of each complex element.
  !
  !  Input parameters:
  !    n - Number of vector elements to inspect.
  !    x - Complex vector containing at least n elements.
  !
  !  Result:
  !    max_component - Largest magnitude of any real or imaginary component.
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

  !Function real_upper_factor_is_singular tests whether a real
  !upper-triangular factor can be used safely by an unguarded triangular solve.
  !The active factor scale is
  !
  !                     scale = max(abs(R(1:n,1:n)))
  !
  !and a diagonal element is considered unusable when
  !
  !        abs(R(i,i)) <= max(tiny(1.0_wp),epsilon(1.0_wp)*scale).
  !
  !This test detects exact rank deficiency as well as a diagonal that is lost
  !at the relative resolution of the stored factor. xGEQRF does not report
  !rank deficiency through INFO, so the test is required before DTRSV.
  !
  !  Input parameters:
  !    r - Real array containing the active upper-triangular factor.
  !    n - Active order of the factor.
  !
  !  Result:
  !    is_singular - True when at least one active diagonal element satisfies
  !                  the threshold above; false otherwise.
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

  !Function complex_upper_factor_is_singular tests whether a complex
  !upper-triangular factor can be used safely by ZTRSV. The factor scale is the
  !largest complex modulus in R(1:n,1:n), and the diagonal threshold is
  !
  !        max(tiny(1.0_wp),epsilon(1.0_wp)*factor_scale).
  !
  !  Input parameters:
  !    r - Complex array containing the active upper-triangular factor.
  !    n - Active order of the factor.
  !
  !  Result:
  !    is_singular - True when a diagonal modulus is not greater than the
  !                  threshold; false otherwise.
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

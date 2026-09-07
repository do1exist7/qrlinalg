# ECGPACK integration handoff

This document is the implementation contract for integrating qrlinalg into
ECGPACK. It is intended to be sufficient for another agent to perform the
integration without relying on repository history or prior conversation.

## Scope

qrlinalg provides serial state objects for real symmetric and complex
Hermitian generalized eigenproblems

```text
H*x = lambda*S*x.
```

For one fixed real shift `sigma`, a state stores explicit factors

```text
H - sigma*S = Q*R.
```

The state can run shifted inverse iteration and can update the factors after a
symmetric/Hermitian replacement, append, or principal deletion. It contains
no MPI state. Optional OpenMP is confined to sufficiently large internal GEMM
calls used by fresh QR factorization.

The integration must preserve ECGPACK ownership of the physical `H` and `S`
matrices. qrlinalg does not store them and does not retain pointers to caller
arrays.

## Recommended integration architecture

Prefer source-level integration using ECGPACK's existing `wp_def` module as
the single precision authority:

1. Import `src/qrlinalg.f90` and the complete `src/qrupdate/` directory.
2. Do not compile qrlinalg's `src/wp_def_*.f90` files if ECGPACK already
   supplies a compatible module named `wp_def` with public integer parameter
   `wp`.
3. Ensure ECGPACK's `wp` is exactly the selected kind `8`, `10`, or `16`.
4. Compile `BLAS.f` with `QRLINALG_WP` set to that same value.
5. Select exactly one provider for the traditional BLAS/LAPACK symbols used
   by qrlinalg. Do not link two objects defining `DGEMM`, `ZGEMM`, `DGEQRF`,
   `ZGEQRF`, and related routines.

Building qrlinalg as a separate static archive is also possible, but both
projects use the module name `wp_def` and may contain traditional Netlib symbol
names. Source-level integration makes these ownership conflicts explicit and
is less dependent on module search paths and archive link order.

Never replace the bundled precision-generic BLAS/LAPACK with an ordinary
binary64 system library for `wp=10` or `wp=16`. Despite their D/Z names, the
bundled routines use `real(wp)` and `complex(wp)`.

The imported qrupdate sources are GPL-3.0-or-later and the bundled Netlib
sources retain their notices. Preserve `src/qrupdate/COPYING`, upstream
metadata, and source notices, and confirm that ECGPACK's distribution terms
are compatible before redistributing an integrated tree.

## Compilation order and flags

After the canonical `wp_def` module, compile in this order:

```text
wp_def
  -> src/qrupdate/BLAS.f
  -> src/qrupdate/LAPACK.f
  -> src/qrupdate/qrupdate_linalg.f90
  -> src/qrupdate/qrupdate_error.f90
  -> src/qrupdate/qrupdate_real.f90
  -> src/qrupdate/qrupdate_complex.f90
  -> src/qrupdate/qrupdate.f90
  -> src/qrlinalg.f90
  -> ECGPACK callers
```

`BLAS.f` requires fixed-form preprocessing and a matching macro:

```text
gfortran:   -cpp -DQRLINALG_WP=<8|10|16> -ffixed-line-length-none
ifort/ifx:  -fpp -DQRLINALG_WP=<8|10|16> -extend-source
nvfortran:  -Mpreprocess -DQRLINALG_WP=<8|10|16> -Mextend
```

`qrlinalg.f90` must also be preprocessed so `QRLINALG_TESTING` remains absent
from production builds. Never define `QRLINALG_TESTING` in ECGPACK; it exposes
private components solely for qrlinalg's white-box tests.

The current standalone Make interface is:

```sh
make PREC=8 CONFIG=release
make PREC=10 CONFIG=release
make PREC=16 CONFIG=release
make OPENMP=1 PREC=8 CONFIG=release
```

An OpenMP build must also link the final ECGPACK executable with the compiler's
OpenMP flag/runtime. Serial and OpenMP objects must not be mixed.

## Public module surface

Import only the required symbols:

```fortran
use qrlinalg, only: wp, qr_real_state, qr_complex_state, &
                    QR_SUCCESS, QR_ERR_INVALID_ARGUMENT, &
                    QR_ERR_ALLOCATION, QR_ERR_NOT_IMPLEMENTED, &
                    QR_ERR_FACTORIZATION, QR_ERR_SINGULAR, &
                    QR_ERR_NO_CONVERGENCE, QR_ERR_INVALID_STATE, &
                    QR_ERR_DIMENSION_MISMATCH, &
                    QR_ERR_CAPACITY_EXCEEDED, &
                    QR_ERR_ZERO_INITIAL_VECTOR, &
                    QR_ERR_NONPOSITIVE_OVERLAP, qr_status_message
```

The public status values are:

| Symbol | Value | Meaning |
|---|---:|---|
| `QR_SUCCESS` | 0 | Operation completed successfully |
| `QR_ERR_INVALID_ARGUMENT` | 1 | Invalid scalar control, index, or mathematical input property |
| `QR_ERR_ALLOCATION` | 2 | Initialization could not allocate all storage |
| `QR_ERR_NOT_IMPLEMENTED` | 3 | Reserved for future API expansion |
| `QR_ERR_FACTORIZATION` | 4 | LAPACK workspace query or factorization stage failed |
| `QR_ERR_SINGULAR` | 5 | Shifted factorization or iteration is numerically unusable |
| `QR_ERR_NO_CONVERGENCE` | 6 | Iteration limit reached; final approximation is still returned |
| `QR_ERR_INVALID_STATE` | 7 | State is uninitialized, invalid, or missing required owned storage |
| `QR_ERR_DIMENSION_MISMATCH` | 8 | Caller array extents are empty or incompatible |
| `QR_ERR_CAPACITY_EXCEEDED` | 9 | Requested active order exceeds initialized capacity |
| `QR_ERR_ZERO_INITIAL_VECTOR` | 10 | Starting vector is numerically zero in the working precision |
| `QR_ERR_NONPOSITIVE_OVERLAP` | 11 | Final overlap norm is non-positive or numerically zero |

Values `0:6` retain their existing assignments. ECGPACK should compare against
symbols rather than integer literals and should include a default branch for
future status values.

`qr_status_message(info)` is a pure, allocation-free module function that
returns a fixed-length, blank-padded description for diagnostics. For example:

```fortran
if (info /= QR_SUCCESS) then
  write(*,'(a)') trim(qr_status_message(info))
end if
```

Unknown integers have a stable fallback description. ECGPACK should use the
integer symbols for control flow and treat the returned text as descriptive;
the library does not decide whether a status is fatal or retryable.

Two independent concrete state types are provided:

```fortran
type(qr_real_state)    :: real_state
type(qr_complex_state) :: complex_state
```

Both expose the same type-bound method names:

```fortran
call state%initialize(capacity, info)
call state%clear()
call state%factorize_fresh(H, S, shift, info)
call state%factorize_fresh(H, S, shift, info, active_order=n)
call state%replace_symmetric(idx, delta_h, delta_s, info)
call state%append_symmetric(h_column, s_column, info)
call state%delete_symmetric(idx, info)
call state%solve(S, v_initial, x, lambda, tol, max_iter, norm_mode, &
                 rel_acc, num_iter, info)
call state%factorization_residual(H, S, probe, absolute_residual, &
                                  relative_residual, info)

valid = state%is_valid()
n = state%order()
capacity = state%get_capacity()
shift = state%get_shift()
total_updates = state%get_update_count()
updates_since_fresh = state%get_updates_since_fresh()
```

There is no public base type or runtime real/complex dispatch. ECGPACK should
select the concrete type at the same point where it currently selects
`GSEPIIS` or `GHEPIIS`.

State components are private. Production callers cannot inspect `Q`, `R`, or
workspace, but the type-bound queries above expose lifecycle metadata without
allowing mutation. ECGPACK should use them as the authoritative description of
the QR state instead of duplicating order, capacity, shift, validity, or update
counters. `get_shift()` is meaningful as a represented shift only while
`is_valid()` is true. Do not use test preprocessing to bypass the private
factor-storage boundary. If factor access is genuinely required, add a
separately reviewed public API.

Avoid intrinsic assignment between populated states: allocatable component
assignment would deep-copy both dense factors and all workspaces. Pass the
owning state by reference and keep one independent state per concurrent solve.

## State lifecycle

The intended lifetime is:

```text
initialize maximum capacity
  -> factorize physical H and S at one shift
  -> solve and/or perform structural updates
  -> factorize_fresh again when the shift changes or drift policy fires
  -> clear when deterministic early release is required
  -> automatic scope cleanup or DEALLOCATE of an allocatable state object
```

`initialize(capacity, info)` allocates full `capacity`-by-`capacity` Q and R,
all factorization/update workspaces, and two solve vectors. It does not select
an active order or shift. Choose capacity as the largest basis size that the
state must reach without reallocation, not merely the initial order.

Calling `initialize` again releases all existing state before allocating the
new capacity. An allocation or workspace-query failure therefore leaves an
empty state; it does not preserve the previous factorization.

`clear()` releases every state-owned allocation and resets order, capacity,
shift, validity, and both update counters. It is valid for an already empty
state and returns no status. A cleared object remains a valid Fortran object
and may be initialized again with a different capacity.

Without its optional argument, `factorize_fresh` chooses the active order from
the equal square shapes of H and S. ECGPACK may instead pass its capacity-sized
physical arrays with `active_order=n`; only their leading active principal
blocks are then used. The routine reuses initialized storage, resets the
state-internal updates-since-fresh counter, and does not allocate. A shift
change is the dense change
`-delta_shift*S`; it always requires `factorize_fresh`. Do not call a
row/column update to represent a shift change.

Intrinsic deallocation of an allocatable state recursively releases all
private allocatable components. A non-allocatable local state is cleaned up
automatically when its scope ends. `clear()` is therefore required only when a
caller wants to release a long-lived state's storage before its own lifetime
ends. There is no custom final procedure.

## Matrix storage contract

For `factorize_fresh(H,S,...)`, `solve(S,...)`, and
`factorization_residual(H,S,...)`, only the lower triangle of the leading
active block, including the diagonal, is defined and read. The upper triangle
and capacity storage outside that block may be uninitialized. Complex inputs
are interpreted as Hermitian and their opposite triangle is produced by
conjugation or by the Hermitian BLAS operation.

The complex diagonals of physical H and S must be real. Roundoff-sized
imaginary parts are accepted and discarded; larger values return
`QR_ERR_INVALID_ARGUMENT` before existing factors are modified.

ECGPACK may pass full capacity-sized allocatable arrays directly. qrlinalg
obtains the physical leading dimension from each array descriptor and uses only
the first `n` rows and columns. This avoids constructing packed
`H(1:n,1:n)` and `S(1:n,1:n)` temporaries. The matrix dummies are contiguous;
a compiler may create a correct temporary for a legacy noncontiguous section,
but full allocatable matrices retain their physical leading dimensions.

The factorization input must be physical H, not an already shifted matrix:

```fortran
call state%factorize_fresh(H, S, shift, info)
! Capacity-sized storage with active order n:
call state%factorize_fresh(H, S, shift, info, active_order=n)
```

Do not pass `M=H-shift*S` as H; qrlinalg forms the shift internally and doing
so would subtract the shift twice.

## API details

### `initialize`

```fortran
call state%initialize(capacity, info)
```

- `capacity` must be positive.
- Success creates an initialized but unfactorized state.
- Normal factorization, solve, and update operations allocate no storage.
- Possible results: success, invalid argument, allocation failure, or
  factorization error from the LAPACK workspace queries.

### `factorize_fresh`

```fortran
call state%factorize_fresh(H, S, shift, info)
```

- Without `active_order`, H and S must be square, equal-sized, nonempty, and no
  larger than capacity.
- With `active_order=n`, n must be positive and no larger than capacity, and H
  and S must each have at least n rows and columns. Their physical extents may
  otherwise differ.
- Only their lower triangles are read; neither matrix is modified or retained.
- On success, the state represents `H-shift*S=Q*R` and is ready to solve or
  update.
- Invalid state, dimension mismatch, capacity exceeded, and invalid
  mathematical input are reported separately and preserve existing valid
  factors.
- A LAPACK failure happens after factor storage is overwritten and leaves the
  state invalid and without an active order.

### `solve`

```fortran
call state%solve(S, v_initial, x, lambda, tol, max_iter, norm_mode, &
                 rel_acc, num_iter, info)
```

Inputs:

- `S(:,:)`: at least n rows and columns; the leading active lower triangle is
  the same overlap matrix represented by the factors and is expected positive
  definite.
- `v_initial(n)`: nonzero starting vector with a component in the desired
  eigendirection. It is not modified.
- `tol`: direction-change tolerance.
- `max_iter`: positive iteration limit.
- `norm_mode`: final eigenvector normalization.

Outputs:

- `x(n)`: eigenvector approximation.
- `lambda`: real shifted-Rayleigh-quotient eigenvalue.
- `rel_acc`: direction-change estimate, not a rigorous error bound.
- `num_iter`: inverse iterations performed.
- `info`: status.

### `factorization_residual`

```fortran
call state%factorization_residual(H, S, v, absolute_residual, &
                                  relative_residual, info)
```

This operation evaluates the discrepancy between `(H-shift*S)*v` and
`Q*(R*v)` without exposing the private factors. H and S may be capacity-sized;
only their leading active lower triangles are read. The nonzero probe vector v
must have exactly `state%order()` elements. Both outputs are initialized to
zero on validation errors, and the call does not change factors, metadata, or
update counters. It allocates no memory and never refactorizes automatically.

ECGPACK should choose probe vectors and thresholds appropriate to its refresh
policy. A small result covers only the tested action direction and is not a
full matrix-norm guarantee. A typical policy can combine this diagnostic with
`get_updates_since_fresh()` and explicitly call `factorize_fresh` when needed.

Normalization modes are:

| Mode | Real | Complex |
|---:|---|---|
| 0 | `x^T*S*x=1` | `x^H*S*x=1` |
| 1 | `x^T*x=1` | `x^H*x=1` |
| other | `max(abs(x))=1` | largest real or imaginary component magnitude is one |

For `tol>0`, iteration stops at `rel_acc<=tol`. For `tol<=0`, it continues
until the estimate begins to increase while remaining no larger than
`abs(tol)`. This negative-tolerance mode forces at least one comparison
iteration and can be materially more expensive.

The desired eigenvalue must be closer to `shift` than competing eigenvalues.
An exact eigenvalue shift makes the factors singular. Degenerate or tightly
clustered eigenvalues may return an invariant-subspace vector or
`QR_ERR_NO_CONVERGENCE`.

`QR_ERR_NO_CONVERGENCE` is a usable-output status: x, lambda, rel_acc, and
num_iter describe the final approximation. Errors detected before iteration
return zero x and lambda, zero iterations, and huge rel_acc. The state itself
is unchanged by every solve result.
`QR_ERR_ZERO_INITIAL_VECTOR` is detected before iteration, whereas
`QR_ERR_NONPOSITIVE_OVERLAP` is detected after an iterate is formed and can
therefore return `num_iter>0`.

### `replace_symmetric`

```fortran
call state%replace_symmetric(idx, delta_h, delta_s, info)
```

- `idx` is one-based and must lie in `1:n`.
- Each delta vector has length n and represents the change in the conceptual
  full column `H(:,idx)` or `S(:,idx)`, including the diagonal.
- The state applies `d=delta_h-shift*delta_s` to both the column and its
  symmetric/Hermitian row without double-counting the diagonal.
- Complex represented diagonal change must be real within a precision-scaled
  tolerance.
- Invalid state, vector extent, and index or complex-diagonal errors return
  `QR_ERR_INVALID_STATE`, `QR_ERR_DIMENSION_MISMATCH`, and
  `QR_ERR_INVALID_ARGUMENT`, respectively. Every rejected call leaves the
  complete state unchanged.

When ECGPACK stores only the complex lower triangle, construct a conceptual
column delta as follows:

```text
i >= idx: delta_h(i) = change in stored H(i,idx)
i <  idx: delta_h(i) = conjugate(change in stored H(idx,i))
```

Apply the same rule to S. After qrlinalg reports success, update ECGPACK's
physical H and S consistently. Updating caller matrices only after success
keeps the physical matrices synchronized with rejected-state semantics.
Preallocate the required physical matrix capacity before invoking qrlinalg so
the caller-side commit cannot subsequently fail for lack of storage.

### `append_symmetric`

```fortran
call state%append_symmetric(h_column, s_column, info)
```

- The current state must have valid order n and spare capacity.
- Both columns have length `n+1`.
- Elements `1:n` are the conceptual new columns `H(1:n,n+1)` and
  `S(1:n,n+1)`; element `n+1` is the new diagonal.
- Invalid state, exhausted capacity, and incorrect column extent return
  `QR_ERR_INVALID_STATE`, `QR_ERR_CAPACITY_EXCEEDED`, and
  `QR_ERR_DIMENSION_MISMATCH`, respectively. A non-real complex diagonal
  returns `QR_ERR_INVALID_ARGUMENT`.
- Complex H and S diagonals are validated separately as real.
- Success increases the active order by one; rejection preserves all state.

For lower-triangle-only complex storage, write the successful append back as:

```fortran
do i = 1, n
  H(n+1,i) = conjg(h_column(i))
  S(n+1,i) = conjg(s_column(i))
end do
H(n+1,n+1) = cmplx(real(h_column(n+1),wp),0.0_wp,kind=wp)
S(n+1,n+1) = cmplx(real(s_column(n+1),wp),0.0_wp,kind=wp)
```

The conjugation is essential: the API vector is a conceptual column, while
the stored bottom row belongs to the lower triangle.

### `delete_symmetric`

```fortran
call state%delete_symmetric(idx, info)
```

- Removes the principal row and column at one-based index `idx`.
- Success decreases the active order by one.
- Deletion from order one is rejected because qrlinalg has no valid order-zero
  factorization.
- Invalid indices or order-one deletion return `QR_ERR_INVALID_ARGUMENT`; an
  invalid state returns `QR_ERR_INVALID_STATE`. Both leave the complete state
  unchanged.

After success, delete the same principal row and column from ECGPACK's H and S
storage. The new active order is available from `state%order()`.

## Mapping from `GSEPIIS` and `GHEPIIS`

The old interface combines factor extension and solve:

```fortran
call GSEPIIS(k,n,M,nM,invD,B,nB,apprlambda,v,w,tol, &
             lambda,x,relacc,maxiter,specifnorm,numiter,errorcode)
```

The qrlinalg mapping is:

| Original concept | qrlinalg concept |
|---|---|
| `apprlambda` | `shift` passed to `factorize_fresh` |
| Physical A | H passed separately and retained by ECGPACK |
| B | S passed to factorization and solve |
| Destructive shifted/factorized M | Private Q/R state |
| `invD`, `w` | Removed; workspace is private |
| Input/output `v` | Read-only `v_initial` |
| `x`, `lambda`, `RelAcc`, `NumIter` | Corresponding solve outputs |
| `SpecifNorm` | `norm_mode` with the same values |
| old error 0 | `QR_SUCCESS` |
| old error 1 | normally `QR_ERR_SINGULAR` |
| old error 2 | `QR_ERR_NO_CONVERGENCE` |

For a fresh full solve, replace the old combined call with:

```fortran
call state%factorize_fresh(H, S, apprlambda, info)
if (info == QR_SUCCESS) then
  call state%solve(S, v, x, lambda, tol, maxiter, specifnorm, &
                   relacc, numiter, info)
end if
```

The original caller often precomputes `M=A-apprlambda*B` and permits GSEPIIS
or GHEPIIS to destroy it. Do not carry that behavior into the QR backend.
Retain the physical H and S so structural updates and later shift changes can
be represented correctly.

The original `k` parameter means that an LDLT/LDLH leading submatrix is already
factorized and the routine should extend it. qrlinalg has no `k` argument:

- For a fresh order-n solve, factorize the complete n-by-n physical matrices.
- To grow a valid QR state from order n-1 to n, call `append_symmetric` with
  the physical nth column.
- qrlinalg does not represent an order-zero factorization, so the first active
  matrix must be created with `factorize_fresh` at order at least one.

Iteration counts and rel_acc may differ slightly from the original solver.
Compare real eigenvectors after sign alignment and complex eigenvectors after
phase alignment. The complex direction coefficient intentionally uses
`v^H*x`, which is phase invariant for QR-backed iteration.

## Minimal owning context

ECGPACK can keep qrlinalg state in its solver context without duplicating the
metadata owned by that state. A representative real integration is:

```fortran
module ecgpack_qr_backend
  use qrlinalg, only: wp, qr_real_state, QR_SUCCESS
  implicit none

  type :: ecgpack_real_qr_context
    type(qr_real_state) :: factors
  end type ecgpack_real_qr_context

contains

  subroutine prepare_context(context, H, S, shift, capacity, info)
    type(ecgpack_real_qr_context), intent(inout) :: context
    real(wp), intent(in) :: H(:,:), S(:,:), shift
    integer, intent(in) :: capacity
    integer, intent(out) :: info

    call context%factors%initialize(capacity, info)
    if (info /= QR_SUCCESS) return
    call context%factors%factorize_fresh(H, S, shift, info)
  end subroutine prepare_context
end module ecgpack_qr_backend
```

Production code should separate initialization from refactorization so an
already allocated context can reuse its capacity when only the shift changes.
The example above combines them only to show state ownership. Subsequent code
can test `context%factors%is_valid()` and read the remaining metadata directly
from the state.

## Update drift and refactorization policy

Structural updates accumulate floating-point drift even when mathematically
reversed. The order-eight stress test remained usable after 20,000 public
updates at every precision, but binary64 factor/orthogonality defects reached
roughly `1e-13` to `2e-12` depending on operation type.

Do not encode 20,000 or any other qrlinalg-wide correctness limit. ECGPACK
should maintain an external updates-since-fresh counter and choose a policy
using its actual dimensions, conditioning, shift proximity, and accuracy
requirements. A fresh factorization resets that external counter. Consider
triggering a fresh factorization when either:

- the external update count reaches an empirically calibrated threshold;
- an eigenpair residual becomes unacceptable;
- solve convergence degrades unexpectedly; or
- the shift changes.

## MPI and OpenMP behavior

qrlinalg has no communicator, collective, rank test, or MPI allocation. For an
MPI build, either keep one independent QR state per rank that owns a local
problem or invoke the state only on the rank that owns the matrices. The
integration layer, not qrlinalg, decides how results are communicated.

Mutable states are not safe for simultaneous calls. Separate OpenMP tasks or
threads require separate states. `OPENMP=1` only enables internal GEMM
worksharing during sufficiently large fresh factorizations; it does not make a
state concurrently callable.

If ECGPACK already parallelizes over MPI ranks, eigenproblems, or outer OpenMP
regions, default to the serial qrlinalg build until oversubscription has been
measured. On the optimization host, four physical threads helped order-1000
factorization, while small matrices and eight-thread SMT did not reliably
benefit.

## Memory and performance expectations

Initialization allocates dense Q and R at full capacity. For scalar byte size
`b` and capacity `c`, factor storage alone is approximately:

```text
real state:    2*b*c^2 bytes
complex state: 4*b*c^2 bytes
```

This excludes caller-owned H/S and vector workspaces. Do not initialize every
MPI rank or every inactive solver context at the global maximum capacity
without accounting for aggregate memory.

Fresh factorization is cubic and stores explicit Q. Replacement is implemented
as two rank-one QR updates; append and deletion use qrupdate row/column
operations. A shift change remains cubic because it requires fresh factors.

## Integration acceptance checklist

Before enabling the QR backend in ECGPACK:

1. Build the integrated source at `wp=8`, `wp=10`, and `wp=16` with exactly one
   `wp_def`, BLAS, and LAPACK provider.
2. Run the qrlinalg semantic suite in serial; if OpenMP is enabled, run it with
   multiple physical threads as well.
3. Verify fresh real and complex solves on analytical noncommuting H/S cases.
4. Pass poisoned or undefined upper triangles and confirm results depend only
   on the lower triangles.
5. Compare against GSEPIIS/GHEPIIS using identical physical H, S, shift,
   starting vector, tolerance, iteration limit, and normalization.
6. Align real vectors by sign and complex vectors by phase before comparison;
   compare generalized residuals, not raw vector components alone.
7. Exercise the real and complex lifecycle: initialize, fresh factorization,
   solve, replacement, append, deletion, fresh factorization at a new shift,
   and cleanup.
8. Confirm every successful structural state change is mirrored exactly once
   in ECGPACK's physical H/S storage and in the queried qrlinalg metadata.
9. Confirm rejected operations leave ECGPACK's matrices synchronized with the
   unchanged qrlinalg state and its queried metadata.
10. Test an iteration-limit result and retain its usable eigenpair when
    `info==QR_ERR_NO_CONVERGENCE`.
11. Test a singular shift and ensure ECGPACK handles `QR_ERR_SINGULAR` without
    attempting to use the returned zero output.
12. Measure hidden array-section copies if ECGPACK matrices use leading
    dimensions larger than the active order.
13. For MPI/OpenMP configurations, verify state ownership and thread/rank
    placement without oversubscription.

## Current limitations that must remain visible

- Only real symmetric and complex Hermitian problems are supported.
- S is expected positive definite for inverse iteration and overlap
  normalization.
- A shift equal or extremely close to an eigenvalue can be singular.
- The desired eigenvalue must be the one closest to the shift for ordinary
  convergence.
- Degenerate and tightly clustered eigenvalues need not yield a unique vector.
- There is no order-zero state, cheap shift update, batch solve, factor getter,
  serialization format, MPI backend, or concurrent access to one state.
- OpenMP accelerates selected fresh-factorization GEMM paths only; structural
  qrupdate operations are not threaded.

The standalone lifecycle example is in `example/example.f90`, and the
data-driven original-solver comparison is under `test/differential/`.

# qrlinalg

`qrlinalg` is the serial QR-state layer intended for ECGPACK generalized
symmetric and Hermitian eigenproblems. Version 0.1.0 implements initialization
and fresh real/complex QR factorization plus generalized inverse iteration.
Symmetric/Hermitian row-and-column replacement, end-appending, and principal-
submatrix deletion update the stored factors in place.

The project vendors the generic-precision `qrupdate-ng` sources under
`src/qrupdate/`, copied from `linalg/src/qrupdate`. That directory records the
upstream commit and retains the QR-update GPL license and the Netlib notices.

## Build

Make selects one working precision at compile time and uses the matching
`wp_def` module for the complete library:

```sh
make                         # release, wp=8, gfortran
make PREC=10 CONFIG=debug
make PREC=16 COMPILER=ifx
make check                   # debug builds at wp=8, 10, and 16
make release OPENMP=1        # opt-in threaded Level-3 QR updates
```

Supported compilers are `gfortran`, `ifort`, `ifx`, and `nvfortran`. Serial
artifacts are written to `build/<configuration>-wp<kind>/`; `OPENMP=1` uses a
separate `-omp` build directory and activates coarse threading in sufficiently
large GEMM calls made by blocked fresh QR. The default remains serial. OpenMP
does not parallelize qrupdate operations and does not make one mutable state
safe for simultaneous calls. The library never uses MPI or a conventional
fixed-double system BLAS/LAPACK.

An fpm 0.13-or-newer build defaults to `wp=8`:

```sh
fpm build
```

For `wp=10` or `wp=16`, change the single `QRLINALG_WP` macro in `fpm.toml`
before building. This is a compile-time choice; one library contains one `wp`.

## Benchmarks

Build and run the independent QR or pristine-LDLT timing suites directly:

```sh
./benchmark/run_qr.sh --sizes 100,200,500,1000 --repetitions 5
./benchmark/run_ldlt.sh --sizes 100,200,500,1000 --repetitions 5
```

Each operation has a separate executable and prints verbose convergence and
timing metrics plus one machine-readable CSV record. These scripts do not use
the Makefile. Perf collection, visualization, and comparison limitations are
documented in `benchmark/README.md`.

## Tests

Run the complete permanent test suite with:

```sh
make check
```

`make test` is an equivalent target. Both commands create checked debug builds
and run every test at all three supported working kinds: `wp=8`, `wp=10`, and
`wp=16`. A nonzero exit status means at least one assertion failed.

During development, run the complete suite at one working kind with:

```sh
make check-one PREC=8
make check-one PREC=10
make check-one PREC=16
```

Pass `COMPILER=ifx`, `COMPILER=ifort`, or `COMPILER=nvfortran` to either form
to exercise another supported compiler. After `make check-one PREC=8`, an
individual executable can be rerun directly, for example:

```sh
./build/test-wp8/test_replacement
```

The seven independently reported executables cover:

- `test_initialization`: invalid initialization, state metadata, and every
  state-owned workspace extent;
- `test_factorization`: real and complex analytical QR reconstruction,
  orthogonality/unitarity, triangularity, caller ownership, and rejected-call
  state preservation;
- `test_replacement`: 100 successive real symmetric updates and 100 successive
  complex Hermitian updates, checking analytical reconstruction and factor
  quality after every update, caller ownership, counters, and rejected updates;
- `test_append`: repeated real symmetric and complex Hermitian order increases
  through full state capacity, with analytical reconstruction after every
  append and checks of diagonal, ownership, counter, and rejection contracts;
- `test_delete`: real and complex append/delete round trips followed by middle,
  first, and last principal-submatrix deletions down to order one, checking
  reconstruction, factor quality, counters, inactive storage, and rejection
  contracts after every operation;
- `test_inverse_iteration`: noncommuting generalized eigenproblems with known
  eigensystems, all normalization modes, both stopping rules, nonconvergence
  with a usable approximation, and recoverable error paths.
- `test_inverse_iteration_failures`: exact degenerate eigenspaces, equidistant
  shifts, missing target components, slowly converging eigenvalue clusters,
  precision-scale singular shifts, indefinite overlap matrices, and
  roundoff-scale starting vectors in both real and complex arithmetic.

The numerical assertions use precision-scaled tolerances. QR factors are
compared through reconstruction, triangularity, and orthogonality/unitarity,
not against one arbitrary choice of column signs or complex phases. Analytical
eigenvectors are compared in a sign- or phase-insensitive manner. Test inputs
are deterministic, so a failure can be reproduced without recording a random
seed.

The failure-regime suite intentionally expects different statuses for
different mathematical limitations. Degeneracy or a start confined to one
invariant subspace can produce a valid eigenpair and `QR_SUCCESS`. Oscillation
and insufficient separation return `QR_ERR_NO_CONVERGENCE` with a finite
approximation. Singular shifted factors and violations of the positive-
definite-overlap or nonzero-start preconditions return `QR_ERR_SINGULAR` or
`QR_ERR_INVALID_ARGUMENT`.

Shared assertions and reference norms live in `test/test_support.f90`. The test
build exposes private state components with `QRLINALG_TESTING` solely so these
invariants can be inspected without enlarging the public API. Normal Make and
fpm builds retain private components and contain no test-only code.

For an optimized regression run at one precision, use:

```sh
make CONFIG=release PREC=8 BUILD_DIR=build/test-release-wp8 \
  QRLINALG_TEST_FLAGS=-DQRLINALG_TESTING test-one
```

### Differential tests against the original solver

The optional differential harness runs identical datasets through qrlinalg and
the original `GSEPIIS`/`GHEPIIS` implementation under `orig/claude`. It is kept
separate from `make check` because the datasets and original reference sources
are not part of the standalone library:

```sh
make differential-build PREC=8
make compare-orig PREC=8 DATA_DIR=/path/to/data
```

`DATA_DIR` must contain `cases.csv`. Relative matrix paths are resolved from
the directory containing the manifest. Required columns are:

```text
case_id,kind,h_file,s_file,shift,tol,max_iter,norm_mode
real_100,real,H_real.dat,S_real.dat,-7.33473,1e-12,30,1
complex_100,complex,H_complex.dat,S_complex.dat,-7.77461,1e-12,30,1
```

Optional `eigen_rtol`, `vector_rtol`, and `residual_tol` columns override the
precision-scaled defaults for an individual case. Blank optional fields use
the defaults. `case_id` may contain letters, digits, dots, underscores, and
hyphens.

Each matrix file contains a four-byte default integer order followed by the
upper triangle, column by column, as binary64 real or complex values. The
drivers convert those values to the selected `wp` and reconstruct the lower
triangle by symmetry or Hermitian symmetry. H and S must have the same order.

The comparison requires compatible solver statuses, checks eigenvalues and
generalized residuals, and compares Euclidean-normalized eigenvectors after
real sign or complex phase alignment. Iteration counts and relative-accuracy
estimates are reported but are not required to be identical. Raw tagged output
and diagnostics from both implementations are retained under
`build/differential-<configuration>-wp<kind>/results/`.

Both drivers are serial. The original `linalg.f90` is compiled unchanged
against a test-only compatibility module and is configured to select its
existing single-process numerical branches. The differential harness therefore
does not require an MPI compiler, launcher, or runtime.

## State and ownership

The module exposes two independent concrete types, `qr_real_state` and
`qr_complex_state`. There is no inheritance and no parameterized derived type.
Each state privately owns:

- explicit full `Q` and `R` arrays sized to its initialized capacity;
- factorization (`tau` and work), QR-update, and solve workspaces;
- active dimension, capacity, represented shift, and factor-validity state;
- lifetime structural-update and updates-since-fresh counters.

The complex state also owns real workspace required by complex rotations.
Initialization queries both LAPACK factorization stages at the maximum basis
size and allocates the larger recommended workspace. All storage is therefore
allocated before factorization, update, or solve loops begin.

ECGPACK remains the owner of `H` and `S`. A QR state neither copies nor retains
pointers to them, and it never permanently stores `M = H - shift*S`. A fresh
factorization reads only the lower triangles of `H` and `S`, reconstructs the
full symmetric or Hermitian `M` directly in the `Q` buffer, extracts `R`, and
generates explicit `Q`. The upper triangles of caller matrices are ignored.

Approximate factor storage, for real scalar size `b` and capacity `c`, is
`2*b*c^2` bytes for a real state and `4*b*c^2` bytes for a complex state.
Caller-owned `H` and `S` bring the application total to roughly four real or
four complex full matrices.

## Numerical path

Fresh factors represent

```text
H - shift*S = Q*R.
```

Inverse iteration computes `w = S*v`, then `y = Q**T*w` for real data or
`y = Q**H*w` for complex data, followed by the triangular solve `R*x = y`.
The final shifted Rayleigh quotient reuses the factors without storing M:

```text
M*x = Q*(R*x)
lambda = shift + (x**H*M*x)/(x**H*S*x).
```

A shift change is the dense update `-delta_shift*S`, so it invalidates the
existing factors and requires a fresh factorization.

For replacement at index `i`, the physical change is
`d = delta_h - shift*delta_s`. The implementation applies two `qr1up`
operations without double-counting the diagonal:

```text
d*e_i**T + e_i*(d-d(i)*e_i)**T       real
d*e_i**H + e_i*(d-d(i)*e_i)**H       complex
```

Complex replacement preserves Hermitian symmetry and requires a real diagonal
change within a precision-scaled tolerance. Append inserts the shifted column
at `n+1` with `qrinc`, then inserts its symmetric or conjugate-transposed row
with `qrinr`. Complex H and S diagonal inputs must each be real within a
precision-scaled tolerance. Deletion removes the selected column with `qrdec`
and then the corresponding row with `qrder`. Because QR-update routines modify
vector arguments, caller inputs are copied into state-owned workspace.

## API status

`initialize(capacity, info)` reserves storage for matrices up to `capacity`,
leaves the active order at zero until `factorize_fresh`, and returns
`QR_SUCCESS`, `QR_ERR_INVALID_ARGUMENT`, or
`QR_ERR_ALLOCATION`. It may return `QR_ERR_FACTORIZATION` if either bundled
LAPACK workspace query rejects the requested configuration.

`factorize_fresh(H, S, shift, info)` forms `H-shift*S` directly in the Q
buffer from the lower triangles of H and S, calls the bundled `xGEQRF`,
extracts upper-triangular R, and calls the bundled `xORGQR/xUNGQR` to generate
explicit Q. Complex diagonal inputs must be real within a precision-scaled
tolerance. It allocates nothing and does not retain H or S. It returns
`QR_SUCCESS`, `QR_ERR_INVALID_ARGUMENT`, or `QR_ERR_FACTORIZATION`.

`solve(S, v_initial, ...)` implements the mathematical iteration and stopping
rules of the pristine `GSEPIIS`/`GHEPIIS` routines through the stored QR
factors. A positive tolerance stops at the requested direction-change estimate;
a negative tolerance continues until that estimate begins to worsen, while
still requiring accuracy `abs(tol)`. Normalization mode 0 produces unit S norm,
mode 1 produces unit Euclidean norm, and any other value retains unit-largest-
component scaling. The routine allocates nothing and leaves S and the initial
vector unchanged.

Successful solves return `QR_SUCCESS`. Invalid dimensions, state, iteration
limit, starting vector, or non-positive S norm return
`QR_ERR_INVALID_ARGUMENT`; an unusable triangular factor returns
`QR_ERR_SINGULAR`. `QR_ERR_NO_CONVERGENCE` still returns the best eigenpair
obtained within `max_iter`.

Both state types expose the same numerical method names:

```fortran
call qr%factorize_fresh(H, S, shift, info)
call qr%replace_symmetric(idx, delta_h, delta_s, info)
call qr%append_symmetric(h_column, s_column, info)
call qr%delete_symmetric(idx, info)
call qr%solve(S, v_initial, x, lambda, tol, max_iter, norm_mode, &
              rel_acc, num_iter, info)
```

`replace_symmetric` applies the physical column changes `delta_h` and `delta_s`
at the existing shift without retaining either vector. It preserves the active
order and increments both update counters once after the complete symmetric or
Hermitian operation. A complex diagonal change must be real within a
precision-scaled tolerance. Invalid calls preserve the complete state.

`append_symmetric` accepts new H and S columns of length `n+1`, expands the
stored factors without fresh factorization, and increments both update counters
once. It rejects an unfactorized or full-capacity state, incorrect extents, and
non-real complex diagonal inputs without changing the state or caller arrays.

`delete_symmetric` removes the principal row and column at a one-based active
index and decreases the order by one. It rejects invalid indices,
unfactorized states, and deletion from order one without changing the state.
The order-one restriction keeps every valid factorization nonempty.

Recoverable errors never use `error stop`, and rejected updates do not
increment counters.

The mutable states are not thread-safe. Threads or tasks must use independent
states. They contain no MPI communicator or branch and should not be replicated
accidentally on every rank; a future distributed implementation belongs in a
separate state or backend.

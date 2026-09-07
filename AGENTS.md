# Repository instructions

## Project purpose

`qrlinalg` is a small, independent, serial Fortran library for QR-backed
linear algebra in generalized symmetric and Hermitian eigenvalue problems:

```text
H*x = lambda*S*x.
```

For a fixed real shift `sigma`, a state stores

```text
H - sigma*S = Q*R,
```

with explicit orthogonal/unitary `Q` and upper-triangular `R`. The explicit
factors support shifted inverse iteration and row/column QR updates without
recomputing a full factorization.

The library is deliberately independent of MPI. Do not introduce MPI modules,
communicators, collectives, compiler wrappers, conditional MPI branches, or
`mpirun` into this repository.

## Current implementation status

The following operations are implemented for both `qr_real_state` and
`qr_complex_state`:

- `initialize(capacity, info)`;
- `clear()`;
- `factorize_fresh(H, S, shift, info[, active_order])`;
- `replace_symmetric(idx, delta_h, delta_s, info)`;
- `append_symmetric(h_column, s_column, info)`;
- `delete_symmetric(idx, info)`;
- `solve(S, v_initial, x, lambda, tol, max_iter, norm_mode, rel_acc, num_iter,
  info)`.
- `factorization_residual(H, S, v, absolute_residual, relative_residual,
  info)`.

Both state types also expose the same pure read-only metadata queries:

- `is_valid()`;
- `order()`;
- `get_capacity()`;
- `get_shift()`;
- `get_update_count()`;
- `get_updates_since_fresh()`.

`get_shift()` describes the represented shift only when `is_valid()` is true.
The update-query results have kind `int64`; fresh factorization resets only
the updates-since-fresh count, while initialization resets both counters.

Deletion from an order-one factorization is invalid because the library does
not represent a valid order-zero factorization. Rejected structural operations
leave factors, dimensions, validity, shifts, and counters unchanged.
`clear()` is idempotent, releases all state-owned allocations, resets all
metadata and counters, and returns no status. A cleared object must be
initialized before it can be factorized again.

## Repository layout and source ownership

- `src/qrlinalg.f90` is the library implementation and public API.
- `src/wp_def_{8,10,16}.f90` select the compile-time working kind `wp`.
- `src/qrupdate/` is the bundled, generic-precision BLAS/LAPACK and qrupdate
  implementation. It was imported from the neighboring linalg project.
- `test/test_support.f90` contains assertions and precision-generic reference
  norms shared by the independent test executables.
- `test/test_initialization.f90`, `test/test_factorization.f90`,
  `test/test_replacement.f90`, `test/test_append.f90`,
  `test/test_delete.f90`, and
  `test/test_inverse_iteration.f90` contain white-box state-contract and
  analytical numerical tests grouped by behavior.
- `test/test_metadata.f90` verifies the public real and complex metadata
  queries across the complete state lifecycle without accessing components.
- `test/test_status_codes.f90` fixes the public numeric status assignments and
  prevents existing status values from being renumbered.
- `test/test_factorization_residual.f90` verifies capacity-sized active blocks,
  physical leading dimensions, inactive sentinels, drift measurement across
  structural updates, fresh-factor comparisons, and validation behavior.
- `test/test_inverse_iteration_failures.f90` contains analytical tests for
  degenerate, clustered, oscillatory, singular, and precondition-violating
  inverse-iteration regimes.
- `test/differential/` contains the optional data-driven comparison harness
  for qrlinalg and the original `GSEPIIS`/`GHEPIIS` implementation.
- `benchmark/` contains standalone operation-specific qrlinalg and
  pristine-LDLT timing drivers, shell build/run/perf entry points, Python
  visualization tools, shared dataset input support, and benchmark guidance.
- `orig/` contains reference implementations, including the pristine
  `GSEPIIS` and `GHEPIIS` algorithms. Treat all files under `orig/` as
  read-only reference material.
- `README.md` describes the public behavior and build interface.
- `Makefile` is the primary build and test entry point.
- `fpm.toml` provides a secondary fpm library build.
- `build/` contains generated objects, module files, archives, and test
  executables. Never commit its contents.

Do not reformat or broadly edit `src/qrupdate/BLAS.f`,
`src/qrupdate/LAPACK.f`, or the qrupdate sources as part of unrelated work.
Preserve their licensing files and notices. Change bundled numerical sources
only when a task specifically requires it and a focused regression test
demonstrates the need.

## Precision-generic requirements

Working precision is selected at compile time and is always one of `8`, `10`,
or `16`. A single library build contains one value of `wp`.

All new numerical declarations and constants must be kind-correct:

```fortran
real(wp)
complex(wp)
0.0_wp
1.0_wp
cmplx(real_part, imag_part, kind=wp)
epsilon(1.0_wp)
tiny(1.0_wp)
```

Do not introduce unqualified real literals into working-precision arithmetic.
Do not use `double precision`, `real(8)`, or fixed-double complex types in the
library implementation.

The bundled routines retain traditional Netlib names such as `DGEQRF`,
`DORGQR`, `ZGEQRF`, `ZUNGQR`, `DSYMV`, and `ZHEMV`, but their implementation
uses the selected `wp`. Their D/Z prefixes do not imply fixed eight-byte
storage in this project.

Never replace the bundled BLAS/LAPACK with a conventional system
double-precision library for `wp=10` or `wp=16`. Any new external BLAS/LAPACK
call from `qrlinalg.f90` must have an explicit interface importing `wp`.

## State and ownership invariants

The concrete real and complex state types are independent; do not introduce
inheritance or a parameterized derived type merely to remove duplicated code.
Their public method names and behavior should remain parallel.

Each initialized state owns:

- full `capacity` by `capacity` arrays for explicit `Q` and `R`;
- Householder coefficients `tau`;
- reusable factorization workspace;
- reusable structural-update and residual-action workspace;
- two capacity-length solve vectors;
- active order, capacity, shift, validity, and update counters.

The complex state additionally owns the real workspace needed by complex
qrupdate rotations.

Caller ownership rules are strict:

- the caller owns `H`, `S`, starting vectors, and physical update vectors;
- the state must not retain pointers to caller arrays;
- the state must not permanently store `H`, `S`, or
  `M = H - shift*S`;
- successful factorization and solve operations must not modify caller input
  arrays;
- the mutable state is not thread-safe; independent threads/tasks require
  independent states.

Initialization allocates all normal operating storage. Fresh factorization,
inverse iteration, and structural updates must not allocate within their
numerical loops.

## Fresh QR factorization

Fresh factorization represents

```text
M = H - shift*S = Q*R.
```

The optional final `active_order=n` argument permits H and S to retain larger
physical extents. When present, n must be positive, fit within state capacity,
and fit within both extents of each matrix; only the leading n-by-n blocks are
used. When absent, the original equal-square shape validation applies. Matrix
dummies are contiguous and BLAS leading dimensions come from their descriptors
so full capacity-sized allocatable arrays require no packed active temporary.

Form `M` directly in the state-owned `Q` buffer. Do not create an additional
full shifted-matrix temporary. Only the lower triangles of caller-owned `H`
and `S`, including their diagonals, are defined and may be referenced. Generate
the opposite triangle of `M` by symmetry or conjugate symmetry. Traverse each
stored lower-triangular column with the first index in the inner loop so input
elements are contiguous.

For complex inputs, validate the imaginary parts of the H and S diagonals
separately against precision-scaled lower-triangle norms before factor storage
is overwritten. Discard accepted roundoff-sized imaginary diagonal parts so
the LAPACK input is exactly Hermitian.

The required factorization stages are:

```text
xGEQRF(Q) -> compact Householder Q and upper-triangular R
copy upper triangle to state R and clear the strict lower triangle
xORGQR/xUNGQR(Q) -> explicit orthogonal/unitary Q
```

Use the workspace allocated by `initialize`. Initialization must query both
factorization stages with `LWORK=-1` at `capacity` and allocate the larger
recommendation, never a hard-coded assumed block size.

All shape, initialization, and capacity checks must happen before existing
factors are overwritten. An invalid argument preserves an existing valid
factorization. Once factor storage has been overwritten, a LAPACK failure must
leave `valid=.false.` and `n=0`; partial factors must never be exposed as valid.

A successful fresh factorization:

- sets `n` and `shift`;
- sets `valid=.true.`;
- resets `updates_since_fresh` to zero;
- preserves the lifetime `structural_updates` counter.

## Inverse-iteration mathematics

The solve methods preserve the mathematical behavior of the reference
`GSEPIIS` and `GHEPIIS` inverse iterations. For the current vector `v`, solve

```text
(H - shift*S)*x = S*v
```

through stored QR factors:

```text
w = S*v
y = Q**T*w       real
y = Q**H*w       complex
R*x = y
```

Use bundled `DSYMV`/`ZHEMV` for `S*v`; only the lower triangle of `S` is
referenced. S may have extents larger than the active order, and its physical
first extent must be passed as the BLAS leading dimension. Use `DGEMV`/`ZGEMV`
for applying `Q` and `DTRSV`/`ZTRSV` for the triangular solve.

## Factorization-drift residual

Both state types expose `factorization_residual`, which compares the physical
action `(H-shift*S)*v` with the factor action `Q*(R*v)` for one nonzero probe
vector. It returns the Euclidean difference norm and that norm divided by the
sum of both action norms plus `tiny(1.0_wp)`. H and S may be capacity-sized;
only their leading active lower triangles are read with descriptor-derived
leading dimensions.

The operation reuses three vector slices of state-owned update workspace and
must not allocate, expose factors, change metadata or counters, or trigger
automatic factorization. It initializes both residual outputs to zero on every
error path. The caller owns probe selection, tolerance selection, and refresh
policy; one action residual is not a full matrix-norm bound.

Real iterates are scaled by `max(abs(x))`. Complex iterates are scaled by the
largest magnitude among every real and imaginary component, not by the
largest complex modulus.

The direction-change estimate is:

```text
alpha   = (x**T*v)/(v**T*v)       real
alpha   = (v**H*x)/(v**H*v)       complex
rel_acc = ||x-alpha*v||_2 / ||x||_2
```

The complex coefficient must conjugate the previous iterate `v`, not the new
iterate `x`. This makes the direction comparison invariant under arbitrary
complex phase and prevents a roundoff-level phase change from causing false
nonconvergence in the QR-backed solve. The literal `x**H*v` expression in the
original GHEPIIS helper is not suitable after replacing its LDLH solve with QR.

Stopping behavior must remain compatible with the reference algorithms:

- `tol > 0`: stop when `rel_acc <= tol`;
- `tol <= 0`: stop only after `rel_acc` begins to increase and the current
  value is at most `abs(tol)`; this forces at least one comparison iteration;
- reaching `max_iter` returns `QR_ERR_NO_CONVERGENCE` but still returns and
  normalizes the final eigenpair approximation.

Compute the eigenvalue from the shifted Rayleigh quotient:

```text
M*x    = Q*(R*x)
lambda = shift + (x**T*M*x)/(x**T*S*x)       real
lambda = shift + real(x**H*M*x)/(x**H*S*x)  complex
```

Do not reconstruct or retain the full shifted matrix merely to calculate the
quotient.

Normalization modes are part of the API contract:

- `norm_mode == 0`: unit overlap norm, `x**T*S*x=1` or `x**H*S*x=1`;
- `norm_mode == 1`: unit Euclidean norm;
- every other value: retain the iteration's largest-component scaling.

The desired eigenvalue must be closer to `shift` than competing eigenvalues.
A shift equal to an eigenvalue makes the factors singular. Degenerate or
roundoff-scale-separated eigenvalues may produce a vector in an invariant
subspace and a nonconvergence status; do not silently change the shift or the
algorithm to conceal that behavior.

## Structural-update mathematics

Changing the shift adds the dense matrix `-delta_shift*S`. A shift change is
not a cheap row/column update and requires fresh factorization.

For replacement at active index `i`, physical changes `delta_h` and `delta_s`
produce

```text
d = delta_h - shift*delta_s.
```

The symmetric/Hermitian matrix change is represented without double-counting
the diagonal by:

```text
d*e_i**T + e_i*(d-d(i)*e_i)**T       real
d*e_i**H + e_i*(d-d(i)*e_i)**H       complex
```

The implementation uses two rank-one `qr1up` operations. Caller vectors are
copied to state-owned workspace first because qrupdate routines may modify
vector arguments.

Append at `n+1` uses `qrinc` followed by `qrinr`. Deletion uses `qrdec` followed
by `qrder`, rejects deletion from order one, and clears newly inactive factor
storage. Complex operations must preserve Hermitian symmetry and must reject a
diagonal change with an imaginary part larger than a precision-scaled
tolerance.

Increment `structural_updates` and `updates_since_fresh` only after a complete,
successful structural operation. Never count a rejected or partially failed
operation.

## Error and state-transition rules

Public status codes are defined in `qrlinalg.f90`:

- `QR_SUCCESS`;
- `QR_ERR_INVALID_ARGUMENT`;
- `QR_ERR_ALLOCATION`;
- `QR_ERR_NOT_IMPLEMENTED`;
- `QR_ERR_FACTORIZATION`;
- `QR_ERR_SINGULAR`;
- `QR_ERR_NO_CONVERGENCE`;
- `QR_ERR_INVALID_STATE`;
- `QR_ERR_DIMENSION_MISMATCH`;
- `QR_ERR_CAPACITY_EXCEEDED`;
- `QR_ERR_ZERO_INITIAL_VECTOR`;
- `QR_ERR_NONPOSITIVE_OVERLAP`.

The numeric values `0:6` are retained for compatibility and new status values
are appended. Use `QR_ERR_INVALID_ARGUMENT` only for invalid scalar controls,
indices, and mathematical input properties. Distinguish invalid state, caller
array extents, exhausted capacity, a numerically zero initial vector, and a
non-positive overlap quadratic form with their dedicated codes.

The pure module function `qr_status_message(info)` returns the stable
human-readable description of a known status and a defined fallback for an
unknown integer. It performs no allocation or I/O. Caller control flow must use
the integer symbols; classifying conditions as fatal or retryable remains
caller policy.

Recoverable library errors return through `info`. Do not use `error stop`,
`stop`, or process termination in library routines.

Validate complete public arguments before destructive state changes whenever
possible. Document whether each failure preserves old factors or invalidates
the state. Keep output values finite and defined on recoverable errors.

`xGEQRF` does not report numerical rank deficiency through positive `INFO`.
Before triangular solves, test the active diagonal of `R` against the existing
precision-scaled threshold. Do not call an unguarded triangular solve with an
unusable diagonal.

## Fortran implementation style

- Use `implicit none` and explicit `intent` on every dummy argument.
- Keep the module default `private`; expose only intentional API symbols.
- Keep real and complex implementations adjacent and behaviorally parallel.
- Place private helper subroutines and functions at the end of
  `src/qrlinalg.f90`, after the public type-bound method implementations.
- Prefer clear loops and explicit workspace slices over hidden allocation or
  large array temporaries in performance-sensitive paths.
- Preserve Fortran column-major access order in nested matrix loops.
- Use descriptive local names such as `matrix_n`, `norm_of_diff`, and
  `overlap_norm_squared`; avoid unexplained one-letter temporaries except for
  standard mathematical matrices/vectors.
- Keep lines within normal free-form Fortran limits where practical.
- Do not add MPI, global mutable configuration, or hidden I/O to the library.

### Comment style

Comments are part of the library interface and must follow the explanatory
style of `orig/claude/linalg.f90`. They are written for any future reader, not
as a conversation with the current user or as a diary of the implementation.

Each nontrivial subroutine or function should have a standalone header that
contains, as applicable:

The `subroutine` or `function` statement must appear before its explanatory
header. Place the header immediately after the complete routine statement and
before the dummy-argument and local declarations.

1. The routine's purpose and the mathematical problem it solves.
2. The defining equations and important algorithmic assumptions.
3. A complete `Input parameters`, `Input/output parameter`, `Output
   parameters`, and/or `Result` section.
4. Array dimensions, triangle conventions, normalization modes, and ownership.
5. State changes, workspace use, allocation behavior, and counter behavior.
6. Every possible status code and whether useful outputs remain available.
7. Numerical limitations such as singular shifts, degeneracy, or
   precision-scaled thresholds.

Inline comments should explain the mathematical or storage reason for an
operation: for example, why `R` must be copied before explicit `Q` is generated,
why the lower triangle is cleared, or how workspace is partitioned. Do not
merely restate the next assignment.

Avoid conversational or historical phrasing such as:

- "this is done for you";
- "we now do";
- "unlike the old code";
- "this deliberately mirrors";
- descriptions of what an agent changed during the current task.

Historical references are appropriate only when they define a required
compatibility contract, such as preserving `GSEPIIS`/`GHEPIIS` stopping
semantics. Do not be sparse with comments in numerical routines, but keep every
comment technically precise and useful outside the context of the current
conversation.

## Build rules

The primary build interface is Make:

```sh
make                              # release, gfortran, wp=8
make release PREC=8
make debug PREC=10
make release PREC=16 COMPILER=ifx
make check
make compare-orig PREC=8 DATA_DIR=/path/to/data
./benchmark/run_qr.sh --precision 8 --repetitions 5
./benchmark/run_ldlt.sh --precision 8 --repetitions 5
```

Supported compiler selectors are `gfortran`, `ifort`, `ifx`, and `nvfortran`.
Generated files belong under `build/<configuration>-wp<kind>/` or the
dedicated test build directories.

The dependency order is:

```text
selected wp_def
  -> bundled BLAS/LAPACK
  -> qrupdate_linalg, qrupdate_error, qrupdate_real, qrupdate_complex
  -> qrupdate
  -> qrlinalg
  -> tests
```

`fpm build` is a secondary check when fpm is installed. `fpm.toml` defaults to
`QRLINALG_WP=8`; selecting `wp=10` or `wp=16` requires changing that single
macro before the fpm build. Do not claim an fpm test was run when fpm is absent.

## Testing requirements

`make check` and `make test` run every test executable in checked builds at all
supported working kinds:

```text
wp=8
wp=10
wp=16
```

Run `make check` after any numerical, interface, workspace, preprocessing, or
comment-layout change in `src/qrlinalg.f90`. Also run `make release` after
changes that may behave differently under optimization. For performance-
sensitive or aliasing changes, run an optimized test executable as well.

The permanent tests must continue to cover:

- valid and invalid initialization;
- idempotent public cleanup, complete storage release, and reinitialization
  after cleanup;
- workspace extents and initial metadata;
- real and complex fresh QR factorization;
- lower-triangle-only H and S input, including poisoned unused upper entries;
- rejection of non-real complex H and S diagonals before state mutation;
- `Q*R` reconstruction;
- orthogonality/unitarity and triangularity;
- preservation of caller matrices and vectors;
- state validity, shifts, dimensions, and counters;
- repeated real and complex replacement and append updates, with analytical
  reconstruction checked after every operation;
- real and complex inverse iteration on analytical generalized eigenproblems;
- noncommuting analytical `H` and `S`, not only diagonal or simultaneously
  diagonalizable trivial cases;
- all three eigenvector normalization modes;
- positive and negative tolerance behavior;
- iteration-limit nonconvergence with a usable returned approximation;
- singular factors, zero starting vectors, and invalid state errors;
- stable public status values and distinct state, dimension, capacity,
  zero-start, and non-positive-overlap failures, including exact status
  descriptions and the unknown-code fallback;
- exact degeneracy, equidistant shifts, missing target components, slowly
  converging clusters, indefinite overlap matrices, and roundoff-scale starts.

Scale numerical tolerances with `epsilon(1.0_wp)`. Compare eigenvectors in a
sign- or phase-invariant manner. Do not compare a valid QR factorization to one
fixed choice of column signs/phases. Prefer reconstruction residuals,
orthogonality/unitarity, triangularity, generalized eigenpair residuals, and
known analytical eigenvalues.

`QRLINALG_TESTING` may expose otherwise private state components only in the
dedicated white-box test build. Normal Make and fpm library builds must retain
private components and contain no alternate numerical implementation.

The optional `make compare-orig` target is not part of `make check`. It requires
`orig/claude`, Python 3, and a data manifest. Compile both drivers at the same
`wp` and pass identical shift, tolerance, iteration-limit, normalization, and
starting-vector conventions. Preserve the original solver sources; use the
dedicated driver and serial MPI compatibility layer in `test/differential/`
for machine-readable output. The compatibility layer must remain test-only and
must force the original linalg module onto its existing serial branches.

Differential comparison must account for mathematical nonuniqueness. Align
real eigenvectors by sign and complex eigenvectors by phase, compare normalized
vectors, and use precision-scaled tolerances for eigenvalues and residuals.
Iteration counts and rough relative-accuracy estimates are diagnostic because
roundoff can move the stopping decision by an iteration. Store generated
driver output only below `build/`.

## Editing and version-control workflow

Before editing:

- inspect `git status` and preserve unrelated user changes;
- read the relevant implementation, tests, and corresponding routine in
  `orig/` when compatibility matters;
- identify whether a change affects real, complex, or both paths;
- identify every working precision and error path affected.

While editing:

- make focused changes with `apply_patch`;
- do not modify generated build products;
- do not weaken tests to accommodate an unexplained numerical difference;
- keep README behavior descriptions synchronized with public API changes;
- preserve the independent repository and do not import MPI-era global state.

After editing:

- run `git diff --check`;
- run the smallest focused test while iterating;
- run `make check` before handoff;
- run `make release` for production-build verification;
- inspect `git status` and ensure only intended source/documentation changes
  remain;
- report unavailable toolchains honestly.

Create commits only when the user requests them. Keep commits focused and use
messages describing the completed behavior, not the editing process. Never
commit `build/` artifacts.

## Definition of done

A change is complete only when:

- the mathematical and ownership invariants above still hold;
- real and complex paths remain consistent where both apply;
- all arithmetic is generic in `wp`;
- recoverable failures leave documented, safe state and outputs;
- numerical routines allocate no unexpected hot-path storage;
- comments are standalone, extensive, and technically precise;
- analytical tests pass for `wp=8`, `wp=10`, and `wp=16`;
- the production release build succeeds;
- documentation matches the implemented API;
- the working tree contains no accidental or generated files.

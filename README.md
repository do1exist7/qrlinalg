# qrlinalg

`qrlinalg` is the serial QR-state layer intended for ECGPACK generalized
symmetric and Hermitian eigenproblems. Version 0.1 is an architectural
skeleton: initialization allocates reusable storage, while every numerical
operation returns `QR_ERR_NOT_IMPLEMENTED` and leaves finite, safe outputs.

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
```

Supported compilers are `gfortran`, `ifort`, `ifx`, and `nvfortran`. Artifacts
are written to `build/<configuration>-wp<kind>/`. The build is serial and links
the bundled generic-precision BLAS/LAPACK; it never uses MPI or a conventional
fixed-double system BLAS/LAPACK.

An fpm 0.13-or-newer build defaults to `wp=8`:

```sh
fpm build
```

For `wp=10` or `wp=16`, change the single `QRLINALG_WP` macro in `fpm.toml`
before building. This is a compile-time choice; one library contains one `wp`.

## State and ownership

The module exposes two independent concrete types, `qr_real_state` and
`qr_complex_state`. There is no inheritance and no parameterized derived type.
Each state privately owns:

- explicit full `Q` and `R` arrays sized to its initialized capacity;
- factorization (`tau` and work), QR-update, and solve workspaces;
- active dimension, capacity, represented shift, and factor-validity state;
- lifetime structural-update and updates-since-fresh counters.

The complex state also owns real workspace required by complex rotations.
Initialization allocates for the expected maximum basis size so later update
and solve loops need not allocate.

ECGPACK remains the owner of `H` and `S`. A QR state neither copies nor retains
pointers to them, and it never permanently stores `M = H - shift*S`. A future
fresh factorization will form `M` directly in the `Q` buffer, extract `R`, and
generate explicit `Q`.

Approximate factor storage, for real scalar size `b` and capacity `c`, is
`2*b*c^2` bytes for a real state and `4*b*c^2` bytes for a complex state.
Caller-owned `H` and `S` bring the application total to roughly four real or
four complex full matrices.

## Intended numerical path

Fresh factors represent

```text
H - shift*S = Q*R.
```

Inverse iteration will compute `w = S*v`, then `y = Q**T*w` for real data or
`y = Q**H*w` for complex data, followed by the triangular solve `R*x = y`.
The final shifted Rayleigh quotient can reuse the factors:

```text
M*x = Q*(R*x)
lambda = shift + (x**H*M*x)/(x**H*S*x).
```

A shift change is the dense update `-delta_shift*S`, so it invalidates the
existing factors and requires a fresh factorization.

For replacement at index `i`, the physical change is
`d = delta_h - shift*delta_s`. The future implementation will apply two
`qr1up` operations without double-counting the diagonal:

```text
d*e_i**T + e_i*(d-d(i)*e_i)**T       real
d*e_i**H + e_i*(d-d(i)*e_i)**H       complex
```

Complex replacement must preserve Hermitian symmetry and require a real
diagonal change within a precision-scaled tolerance. Append means insertion at
`n+1` and will use `qrinc` then `qrinr`; deletion will use `qrdec` then `qrder`.
Because QR-update routines modify vector arguments, caller inputs will first be
copied into state-owned workspace.

## API status

`initialize(max_n, info)` is implemented. It validates positive capacity,
allocates and clears all storage, and returns `QR_SUCCESS`,
`QR_ERR_INVALID_ARGUMENT`, or `QR_ERR_ALLOCATION`.

Both state types expose the same numerical method names:

```fortran
call qr%factorize_fresh(H, S, shift, info)
call qr%replace_symmetric(idx, delta_h, delta_s, info)
call qr%append_symmetric(h_column, s_column, info)
call qr%delete_symmetric(idx, info)
call qr%solve(S, v_initial, x, lambda, tol, max_iter, norm_mode, &
              rel_acc, num_iter, info)
```

In v0.1 these methods return `QR_ERR_NOT_IMPLEMENTED`. `solve` also sets `x`,
`lambda`, `rel_acc`, and `num_iter` to zero. Recoverable errors never use
`error stop`, and failed update stubs do not increment counters.

The mutable states are not thread-safe. Threads or tasks must use independent
states. They contain no MPI communicator or branch and should not be replicated
accidentally on every rank; a future distributed implementation belongs in a
separate state or backend.

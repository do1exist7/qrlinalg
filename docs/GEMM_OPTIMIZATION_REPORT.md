# DGEMM and ZGEMM optimization report

## Purpose

This report describes the retained optimization strategies in the bundled
precision-generic `DGEMM` and `ZGEMM` implementations, the QR workloads that
motivated them, and the alternatives that were measured and rejected. It
covers the serial work in commits `dde28f1`, `648b0f0`, and `07f883d`, plus
the current opt-in OpenMP changes.

The bundled routines use traditional D/Z names, but their element type is the
compile-time working kind `wp`, which is always 8, 10, or 16. A conventional
binary64 BLAS therefore cannot replace them for the extended kinds.

## Workload that drove the changes

The blocked QR implementation reaches GEMM through `xLARFB` in two expensive
phases:

1. GEQRF applies a block of Householder reflectors to the trailing matrix.
2. ORGQR/UNGQR applies the stored reflectors while constructing explicit Q.

At order 1000 with block size 32, the highest-weight calls are approximately:

| QR phase | Real GEMM | Complex GEMM | M | N | K |
|---|---|---|---:|---:|---:|
| GEQRF trailing update | `DGEMM(T,N)` | `ZGEMM(C,N)` | 968 | 32 | 968 |
| Explicit Q generation | `DGEMM(N,T)` | `ZGEMM(N,C)` | 968 | 968 | 32 |

The optimization effort consequently targeted TN/NT for real arithmetic and
CN/NC for complex arithmetic. NN, TT, and other transpose combinations retain
their reference-style loops because they were not material in the measured QR
profile.

## Precision isolation

`BLAS.f` requires the preprocessor macro `QRLINALG_WP` to be 8, 10, or 16.
Make, fpm, differential, and benchmark builds pass the value corresponding to
the selected `wp_def` module.

The preprocessor selects one hot-loop implementation at compile time. This is
important for both performance and maintenance:

- a wp=8 tile is not compiled into wp=10 or wp=16 code;
- changing a kernel for one precision cannot silently perturb another;
- the compiler does not have to eliminate a runtime test of the Fortran kind;
- register and cache choices can reflect the actual representation used by
  the production compiler.

Kind 10 must not be interpreted as a ten-byte storage promise. On the measured
gfortran/x86 build it uses extended-real x87 arithmetic, which imposes very
different register constraints from binary64 SIMD.

## DGEMM strategies

### `wp=8`, transpose/no-transpose: 2-by-4 register tile

For

```text
C = alpha*A**T*B + beta*C,
```

the wp=8 kernel computes two rows by four columns of C at a time. Eight
independent accumulators are held for the complete K reduction.

For each K position the loop loads two adjacent values from A and four values
from B, then reuses them across eight products. Both A and B are traversed
contiguously down Fortran columns. Compared with one scalar dot product at a
time, the tile provides:

- more independent accumulation chains for instruction-level parallelism;
- reuse of each A value across four output columns;
- reuse of each B value across two output rows;
- one final read/write decision for each C element instead of touching C in
  the reduction loop.

Odd rows and output columns not divisible by four use scalar tails. The tails
preserve arbitrary valid dimensions and leading dimensions.

### `wp=10`, transpose/no-transpose: 1-by-4 register tile

The first attempt reused the wp=8 2-by-4 tile. It regressed at wp=10 because
eight extended-real accumulators plus current operands exceeded the practical
x87 register capacity and caused spills.

The retained wp=10 kernel instead computes one row by four columns. Four
accumulators reuse one A value across four B columns and fit the measured x87
register constraints. On the laptop's order-1000 dominant TN shape, the
retained form improved the isolated kernel from 45.24 to 26.60 ms and helped
reduce the combined GEQRF plus DORGQR stages from 2.363 to 1.462 seconds.

### `wp=16`, transpose/no-transpose: scalar reference ordering

The quadruple-precision TN path retains scalar dot products. No wider tile was
kept without evidence that the compiler's software or scalar lowering could
support it efficiently. This is a deliberate precision-specific fallback, not
an assumption that the wp=8 tile is portable to wider arithmetic.

### `wp=10`, no-transpose/transpose: 4-by-1 accumulation

For

```text
C = alpha*A*B**T + beta*C,
```

the wp=10 kernel holds four adjacent rows of one output column in accumulators
for the complete K reduction. Each B value is reused across four contiguous A
values. C is read at most once and written once for each output element.

The dominant NT shape improved from 51.86 to 26.69 ms on the measured host.
Scalar row tails use the same K order and beta rules.

### `wp=8` and `wp=16`, no-transpose/transpose: four columns plus row blocking

These paths first scale or clear C once, then update four adjacent output
columns together. For each K position, one contiguous A vector is reused with
four B values.

Rows are processed in fixed blocks of `ROW_BLOCK=256`. The intent is to keep
four C segments and one A segment in a private cache while a K panel is
consumed. For an element size E and output-column tile NR, the approximate hot
set is:

```text
ROW_BLOCK * (NR + 1) * E.
```

The value 256 is a conservative default measured on the laptop, not a
portable optimum. It should be retested after changing hardware, compiler,
element representation, or output tile width.

## ZGEMM strategies

### `wp=8`, conjugate-transpose/no-transpose: 1-by-4 tile

For

```text
C = alpha*A**H*B + beta*C,
```

the wp=8 kernel computes one row by four output columns. At each K position it
loads and conjugates one A value once, then reuses it with four B values. Four
independent complex accumulators remove the single long dependency chain and
keep C out of the reduction loop.

The highest-weight CN shape was measured near 968-by-32-by-968. Its paired
isolated timing improved from roughly 51-53 ms to 33-34 ms, about 1.55x. The
complete complex fresh-QR benchmark also improved while the real path remained
unchanged.

Columns not divisible by four use the scalar reference loop. wp=10 and wp=16
also retain the scalar CN arithmetic because the binary64 complex tile was not
assumed to be appropriate for their wider real components.

### No retained complex N,C microkernel

A 4-by-1 `ZGEMM(N,C)` candidate was tested because N,C dominates explicit-Q
construction. Its isolated median was 20.78 ms versus 20.09 ms for the
reference implementation. The candidate was rejected, so the current N,C
arithmetic loop must not be described as serially optimized.

N,C is nevertheless eligible for the current coarse OpenMP worksharing because
its output columns are independent. This changes ownership across columns but
does not change the per-column arithmetic kernel.

## Alpha and beta handling

The optimized paths preserve the BLAS interface and its important ownership
rules:

- quick returns and argument validation remain in the parent routine;
- `beta=0` does not read the previous value of C;
- `beta=1` avoids a redundant multiplication;
- arbitrary alpha and beta values remain supported;
- rectangular matrices, zero sizes, padded leading dimensions, and scalar
  tile remainders remain valid;
- real `C` and complex `C` are written only within their declared M-by-N
  output regions.

The microtiles preserve each individual K reduction order. OpenMP assigns
whole output columns or column tiles to one thread, so it does not split a dot
product into a parallel reduction or require atomics.

## OpenMP strategy

OpenMP is opt-in. Default builds treat the directives as comments and remain
serial. `OPENMP=1` selects compiler-appropriate flags and separate `-omp` build
directories so serial and threaded objects cannot be mixed accidentally.

Only the QR-dominant paths are threaded:

| Routine | Arithmetic path | Parallel unit |
|---|---|---|
| DGEMM | T,N | one output column or four-column tile |
| DGEMM | N,T | one output column or four-column tile |
| ZGEMM | C,N | one output column or four-column tile |
| ZGEMM | N,C | one output column |

Static scheduling is used because output tiles have similar work and own
disjoint memory. There is no shared accumulation, packing buffer, allocation,
or synchronization within a tile.

Parallel regions are guarded by two conditions:

1. at least eight output columns/tasks are available; and
2. the approximate multiply work is at least 2,000,000.

The threshold was chosen from actual QR calls. Order-200 calls peak below
about 0.9 million multiplies and should remain serial, while order-500 and
order-1000 trailing updates are large enough to amortize a parallel region.
There is no runtime autotuning or hardware detection.

On four physical cores at wp=8 and order 1000, fresh QR improved from 320.230
to 203.091 ms real and from 1158.649 to 537.331 ms complex. Eight-thread SMT
results were unstable on the four-core host, so four physical threads are the
recommended setting there.

## What was intentionally not changed

- The blocked QR algorithm, its block size, and explicit Q construction were
  not changed by the GEMM work.
- NN, TT, and unprofiled transpose combinations were not specialized.
- No packing buffers, compiler intrinsics, assembly, external binary64 BLAS,
  runtime CPU detection, or hidden allocation was introduced.
- TRMM was not threaded because it remained about 2% of the four-thread QR
  profile.
- qrupdate replacement, append, and deletion kernels were not threaded.
- The rejected wp=10 2-by-4 and complex N,C 4-by-1 candidates were not left
  behind as inactive alternatives.

After OpenMP, GEMM remains the largest useful kernel, but OpenMP wait/barrier
cost and serial panel/DGEMV work are visible. The profile does not justify
another broad BLAS rewrite or TRMM optimization.

## Correctness and measurement support

`test/test_dgemm.f90` covers all N/T/C combinations, sizes zero through three,
rectangular optimized and fallback paths, padded leading dimensions, diverse
alpha/beta values, output padding, and NaN-poisoned C for `beta=0`.

`test/test_zgemm.f90` provides the corresponding complex coverage, including
conjugate transpose and complex alpha/beta values. Both tests include
above-threshold QR-shaped cases so OpenMP builds exercise the parallel TN/NT
and CN/NC paths rather than validating only serialized small calls.

The standalone tools under `benchmark/dgemm/` provide:

- isolated DGEMM and ZGEMM shape timings;
- separate GEQRF and ORGQR/UNGQR stage timings;
- a wp=8/gfortran call-shape profiler;
- opt-in OpenMP builds separate from the public operation benchmarks.

Serial and OpenMP checked builds pass at wp=8, wp=10, and wp=16. The detailed
end-to-end timings, counters, thread scaling, Claude comparison, and pristine
BLAS/LAPACK comparison are recorded in the
[benchmark report](../benchmark/report.md).

## Retuning on another machine

The current tiles, row block, and OpenMP threshold are defaults supported by
one machine's measurements. Retune them in this order:

1. rebuild with the production compiler and native target flags;
2. record element sizes, cache topology, physical-core/SMT topology, and
   affinity;
3. regenerate the QR GEMM call-shape histogram;
4. measure unchanged isolated kernels and complete QR first;
5. tune one precision and transpose path at a time;
6. retain a candidate only if it improves complete QR beyond run-to-run noise
   and passes all-precision semantic tests;
7. remeasure the OpenMP threshold only after the serial kernel is fixed.

Do not copy the wp=8 tile into extended precision, assume kind 10 occupies ten
bytes, or treat `ROW_BLOCK=256` and 2,000,000 multiplies as universal hardware
constants.

## Source map

- `src/qrupdate/BLAS.f`: retained DGEMM/ZGEMM kernels and OpenMP directives.
- `test/test_dgemm.f90`, `test/test_zgemm.f90`: semantic regression coverage.
- `benchmark/dgemm/`: isolated kernels, stage timings, and shape profiling.
- [`benchmark/report.md`](../benchmark/report.md): current complete benchmark
  results.
- [`DGEMM_HPC_HANDOFF.md`](DGEMM_HPC_HANDOFF.md): detailed measurement and
  machine-retuning protocol.

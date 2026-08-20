# Bundled DGEMM optimization: HPC handoff

For a consolidated description of the retained DGEMM, ZGEMM, and OpenMP
strategies, see the [GEMM optimization report](GEMM_OPTIMIZATION_REPORT.md).
This handoff focuses on measured history and the protocol for further
machine-specific work.

## Saved baseline

Work from branch `optimize/qr-hotpaths`. Commit `dde28f1` (`Optimize bundled
DGEMM hot paths`) is the recoverable performance baseline. It contains the
optimized `src/qrupdate/BLAS.f`, semantic tests in `test/test_dgemm.f90`, and
isolated DGEMM, QR-stage, and call-shape tools under `benchmark/dgemm/`.

Rebuild on the HPC node. The build uses `-O3 -march=native` with gfortran, so
laptop objects and executables are neither saved nor portable. OpenBLAS is only
a `wp=8` performance ceiling; it cannot replace the generic bundled BLAS for
`wp=10` or `wp=16`.

## Established measurements

These WSL2 laptop medians show that the first optimization is valuable, but
must not be used to select block sizes for different hardware.

| Measurement | Original | `dde28f1` | Result |
| --- | ---: | ---: | ---: |
| Real QR with explicit Q, n=1000, wp=8 | 1121.6 ms | 389.8 ms | 2.88x |
| TN, M=968 N=32 K=968, wp=8 | - | 10.43 ms | OpenBLAS 1.915 ms; 5.44x gap |
| NT, M=968 N=968 K=32, wp=8 | - | 5.228 ms | OpenBLAS 1.698 ms; 3.08x gap |

Single-threaded OpenBLAS completed the full QR in 137.1 ms. It may use
different floating-point transformations, so this is headroom rather than an
equal implementation comparison.

At `wp=10`, enabling the eight-accumulator 2-by-4 TN kernel regressed
performance because extended-real scalar arithmetic exhausted the x87 register
stack. A follow-up uses a 1-by-4 TN tile and a 4-by-1 NT tile, each with four
accumulators. On the same pinned WSL2 laptop with gfortran 15.2, the retained
follow-up produced these medians:

| wp=10 measurement | Previous | Follow-up | Speedup |
| --- | ---: | ---: | ---: |
| TN, M=968 N=32 K=968 | 45.24 ms | 26.60 ms | 1.70x |
| NT, M=968 N=968 K=32 | 51.86 ms | 26.69 ms | 1.94x |
| DGEQRF stage, n=1000 | 1.170 s | 0.721 s | 1.62x |
| DORGQR stage, n=1000 | 1.193 s | 0.741 s | 1.61x |

The complete two-stage QR median improved from 2.363 s to 1.462 s (1.62x).
Earlier timings in a noisier session were slower in absolute terms, so only
contemporaneous pinned pairs were used for the table.

A subsequent `wp=8` complex optimization targets only `ZGEMM(C,N)`.  Runtime
shape profiling confirmed the highest-weight call as M=968, N=32, K=968 with
LDA=LDB=LDC=1000, alpha=1, and beta=1.  A 1-by-4 output tile reuses one
conjugated A value across four B columns and maintains four independent complex
accumulators.  Two order-reversed pinned comparisons measured 51.08-53.42 ms
before and 33.16-34.35 ms after, approximately 1.55x faster.

The attempted 4-by-1 `ZGEMM(N,C)` tile did not survive the retention gate: its
20.78 ms median was slightly slower than the 20.09 ms reference, so the
original N,C implementation remains.  Paired public fresh-QR medians were:

| wp=8 fresh factorization | Before complex C,N | After | Speedup |
| --- | ---: | ---: | ---: |
| Real n=500 | 50.06 ms | 47.17 ms | 1.06x (noise-sensitive) |
| Real n=1000 | 372.99 ms | 375.52 ms | 0.99x (unchanged path) |
| Complex n=500 | 192.84 ms | 165.92 ms | 1.16x |
| Complex n=1000 | 2.035 s | 1.884 s | 1.08x |

Internally gated counters at complex n=1000 fell from 6.061 to 4.814 billion
cycles and from 9.156 to 8.347 billion instructions; IPC increased from 1.51
to 1.73.  ZGEMM remained dominant but fell from 93.97% to 92.34% of sampled
cycles, corresponding to about a 22% reduction in its absolute cycle count.
Wall-clock timings on this VM varied materially with host load, so counters,
paired medians, minima, and validation must be considered together.

Validation passed `make check` at `wp=8`, `wp=10`, and `wp=16`, plus a release
build. The DGEMM test covers N/T/C, dimensions zero through three, rectangular
fast and fallback paths, padded leading dimensions, alpha/beta 0, 1, -1, and
0.375, output padding, and NaN-poisoned C when beta=0.
The parallel ZGEMM test covers all transpose/conjugation pairs, complex
alpha/beta values, optimized-path tiles and tails, padded leading dimensions,
and the same beta=0 ownership rule.

## Current algorithms

`BLAS.f` selects precision-specific hot kernels with the `QRLINALG_WP`
preprocessor macro. Each build therefore presents only one TN and one NT
implementation to the Fortran compiler; it does not rely on optimization of a
Fortran condition involving the `wp` parameter. Missing or unsupported macro
values stop compilation. Make, fpm, differential, and benchmark builds pass
the same macro value used to select `wp_def`.

The `wp=8` TN (`A**T*B`) fast path computes a 2-by-4 C microtile with eight
accumulators. Its K loop reads `A(:,i:i+1)` and `B(:,j:j+3)` contiguously,
reuses two A and four B values, and writes each C element once. At `wp=10`, a
1-by-4 tile uses four accumulators, reuses one A value across four outputs, and
fits the current operands without x87 spills. Odd M and N remainders are
scalar. `wp=16` retains the original scalar dot-product path. There is no
packing, explicit SIMD, intrinsic, OpenMP, or allocation.

At `wp=10`, NT (`A*B**T`) computes four adjacent rows of one output column.
The four C values remain in accumulators for the complete K reduction, while
each B value is reused across four contiguous A values. C is read at most once
and written once instead of being updated from memory for every K. Scalar row
remainders preserve the same K order and beta handling. Other kinds first
scale or zero C, group four output columns, and tile M with `ROW_BLOCK=256`.

NN and TT remain essentially the Netlib reference paths, as do scalar tails.
This is intentional because the measured QR workload was dominated by TN and
NT. Re-profile before specializing anything else.

For complex arithmetic, only the `wp=8` C,N path is specialized.  It computes
four adjacent output columns for one row with four independent accumulators;
conjugation is applied once when A is loaded.  Column remainders and every
other precision or transpose combination retain the reference loops.  The N,C
candidate was measured and rejected, so it must not be described as optimized.

## Facts versus unresolved hypotheses

Known:

- DGEQRF/DORGQR reach DGEMM through blocked LAPACK/DLARFB calls.
- The laptop's dominant shapes were TN 968x32x968 and NT 968x968x32.
- The current serial implementation materially improves end-to-end QR.
- DGEMM creates no workspace or hidden packing buffers.
- Both inputs in the current TN inner reduction are unit stride.
- NT streams A and four C columns contiguously.

Not yet measured on the target hardware:

- the actual DGEMM shape distribution and LAPACK block size;
- which loops auto-vectorize, vector width, and whether FMA is generated;
- spills, dependency-chain cost, and address-generation overhead;
- reliable L1/LLC miss rates, bandwidth, instructions, cycles, and IPC;
- whether TN is compute/front-end/latency/cache limited;
- target-appropriate row and K-panel block sizes;
- whether DGEMM remains the top bottleneck after the next improvement.

Do not report these hypotheses as findings until target-node evidence exists.

## HPC baseline protocol

Record commit, node/CPU, cache hierarchy, NUMA/SMT topology, compiler version,
full compile/link flags, scheduler allocation, affinity, frequency policy, `wp`,
and element byte size. Use a dedicated node where possible and pin one serial
process to a physical core, avoiding its SMT sibling.

```sh
git switch optimize/qr-hotpaths
git rev-parse HEAD
make check
make release PREC=8
./benchmark/dgemm/build.sh --precision 8 --compiler gfortran

# gfortran/wp=8-only call-shape diagnostic
build/benchmarks/wp8/dgemm/bin/qr_stage_profile 1000 1 1 \
  2>build/benchmarks/wp8/dgemm/call-shapes.csv
```

For isolated timings use the driver's two warmups, at least seven samples, and
enough calls for each sample to last 0.5-1.0 seconds. Report median, minimum,
and spread. Allocation/initialization is outside the timed region. Measure QR
at n=256, 512, 1000, 1500, and 2000 where practical, plus the highest-weight
DGEMM shapes. Repeat `wp=10` using its production compiler.

Collect comparable `perf stat` (or vendor tool) counters for bundled DGEMM and
one-thread OpenBLAS with identical inputs and affinity: cycles, instructions,
branches/misses, L1 loads/misses, LLC loads/misses, elapsed time, and available
floating-point/memory-bandwidth events. Force every relevant OpenBLAS thread
variable to one. Label unsupported or multiplexed events instead of drawing
conclusions from them.

Generate a separate compiler vectorization report without changing production
semantics (`-fopt-info-vec-optimized` and `-fopt-info-vec-missed` for
gfortran). Inspect enough assembly to classify scalar/packed arithmetic, FMA,
spills, remainder overhead, and inner-loop branches/address updates. Save
summaries and raw command metadata in a tracked results Markdown file; never
commit generated objects or `build/`.

## Cache model

For NT with row block `MB`, output-column microtile `NR`, and element bytes
`E`, approximate the primary hot set as:

```text
MB * (NR + 1) * E
```

This accounts for NR C segments plus one A segment, before overhead. Aim below
the full private L1 capacity (roughly 50-75% is a reasonable starting range),
then benchmark a small geometric candidate set. Derive separate candidates for
`wp=8` and `wp=10`; do not assume kind 10 occupies ten bytes.

TN is constrained more by SIMD width, accumulator registers, and reduction
latency than by this NT model. Any packed panel must fit its intended cache and
its copy must be amortized over sufficient output reuse. Since current TN A
and B accesses are already contiguous, packing is a hypothesis: its possible
benefits are alignment, simpler addressing, and SIMD-friendly layout, not
repairing a strided load.

## Retuning after moving to another machine

Treat the current tile shapes and `ROW_BLOCK=256` as defaults for the laptop,
not portable constants. Retune only after rebuilding with the production
compiler and `-march=native` (or the target compiler's equivalent). Do not copy
objects, vectorization conclusions, or timing thresholds from the previous
machine.

Tune in this order:

1. Build and validate the unchanged branch for each required `wp`. Pin one
   process to one physical core, keep frequency policy and SMT placement fixed,
   and record `storage_size(0.0_wp)/8`; kind 10 must not be assumed to occupy
   ten bytes.
2. Regenerate the QR call-shape histogram. Use its highest-weight TN and NT
   shapes rather than assuming that 968x32x968 and 968x968x32 remain dominant.
   A change to LAPACK `NB` or the compiler can change these shapes.
3. Record unchanged isolated-kernel and complete-QR medians. Use at least two
   warmups and seven timed samples, with enough calls for 0.5-1.0 seconds per
   sample. A typical pinned invocation is:

   ```sh
   taskset -c PHYSICAL_CORE \
     build/benchmarks/wp8/dgemm/bin/dgemm_benchmark \
     T N M N K 1 0 CALLS_PER_SAMPLE 9

   taskset -c PHYSICAL_CORE \
     build/benchmarks/wp8/dgemm/bin/zgemm_benchmark \
     C N M N K 1 0 1 0 CALLS_PER_SAMPLE 9
   ```

   Run the corresponding NT case and repeat with the production builds for
   `wp=10` and `wp=16`. Use OpenBLAS only as a one-thread `wp=8` ceiling.
4. Tune the TN register tile first. For `wp=8`, start with a bounded set such
   as `MR={1,2,4}` and `NR={2,4,8}`, rejecting combinations whose accumulator
   count causes spills. For `wp=10` with gfortran/x87, start with at most four
   accumulators, for example 1x2, 1x4, and 2x2; the observed 2x4 regression is
   the reason for this limit. Treat `wp=16` independently because its lowering
   may be scalar or library-call based. Inspect generated assembly instead of
   inferring SIMD width from kind size. Apply the same bounded process to the
   complex C,N tile; do not assume its current 1x4 shape follows the real TN
   optimum on another compiler or register file.
5. Tune NT microtile shape, then `ROW_BLOCK`. Derive the first `MB` estimate
   from the cache model above using about 50-75% of private L1. Benchmark the
   nearest practical values and one value on either side; `64, 128, 256, 512`
   is a reasonable initial geometric sweep, not a required set. Retest if NR
   changes because NR changes the hot-set estimate.
6. Keep LAPACK `NB=32` while tuning BLAS so only one layer changes at a time.
   After retaining BLAS parameters, optionally compare `NB={16,32,48,64}` with
   complete DGEQRF, DORGQR, and public fresh-factorization timings. `NB` is a
   LAPACK/QR parameter rather than a DGEMM parameter; changing it changes the
   DGEMM workload and requires a new call-shape profile and workspace check.

Retain a parameter only when its median improvement is larger than run-to-run
noise, survives a second pinned run, and improves the complete QR workload.
Record minimum and spread as well as median. Reject a faster isolated kernel if
it regresses important smaller shapes, changes benchmark validation, or merely
moves cost into another phase. After each retained candidate run the focused
DGEMM tests, `make check`, and a release build at `wp=8`, `wp=10`, and `wp=16`.

Keep machine-specific choices compile-time and precision-specific. A `wp=8`
tile must not be selected in `wp=10` or `wp=16` through a runtime kind branch.
If several machines must remain supported, prefer documented build-time
presets or conservative defaults over hidden runtime CPU detection. For every
retained preset, record CPU model, compiler/version, flags, cache sizes, tile
shape, `ROW_BLOCK`, LAPACK `NB`, dominant shapes, isolated medians, and full QR
medians in this document or a linked tracked results file.

## Ranked next work

1. **Measure TN vectorization and counters, then run one `wp=8` TN experiment.**
   This has the largest measured gap and highest likely benefit. If assembly is
   scalar/address-heavy, test a bounded packed panel with a SIMD-friendly
   microkernel. If it is already vectorized but dependency-limited, test only a
   larger register tile. Difficulty is medium; rounding risk is low-to-medium
   because reduction order may change. Preserve the extended-kind fallback.

2. **Tune/specialize NT for small K at `wp=8`.** The retained `wp=10` 4-by-1
   accumulator tile resolves the repeated-C-traffic bottleneck locally. For
   `wp=8`, test NR values matched to the target vector/register file and MB
   values derived from the cache model. K=32 may not amortize general packing.

3. **Re-profile the full QR after every retained change.** Only specialize NN,
   TT, or another kernel if its measured contribution becomes material. This
   measurement is cheap and prevents optimizing a former bottleneck.

4. **The first OpenMP experiment is complete.** The retained opt-in build
   distributes independent C columns or four-column tiles only in the
   QR-dominant TN/NT and CN/NC GEMM paths. A two-million-multiply threshold
   keeps small calls serial. On the four-core laptop, order-1000 wp=8 fresh QR
   improved by 1.58x real and 2.16x complex relative to the optimized serial
   build. TRMM remained about 2%, so it was not threaded. See
   `benchmark/report.md` for all-mode comparisons, scaling, counters, and
   limitations. Re-measure the threshold and affinity after moving machines.

## Correctness and retention gate

Every candidate must preserve N/T/C, arbitrary valid leading dimensions,
rectangular and zero-size calls, all alpha/beta cases, and no C read when
beta=0. Run the focused DGEMM test during iteration, then `make check` for all
kinds and a release build. Run complete QR and QR-update/residual tests. Use
K- and epsilon-scaled error bounds rather than bitwise equality. Revert any
candidate that fails correctness, does not improve complete QR reproducibly,
or seriously regresses small matrices.

## Smallest decisive next experiment

First verify both new `wp=10` tiles on the target node with its production
compiler, actual shape histogram, counters, and complete QR timings. Then
change only TN at `wp=8`, choosing between a packed/SIMD-friendly panel and a
larger unpacked register tile according to target evidence. Compare against
commit `dde28f1` and the current branch in the isolated benchmark and complete
QR. This preserves the proven `wp=8` result while separating target-specific
work from the retained extended-real optimization.

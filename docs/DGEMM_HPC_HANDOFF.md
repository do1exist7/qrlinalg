# Bundled DGEMM optimization: HPC handoff

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

At `wp=10`, enabling the fused TN kernel regressed performance. The committed
dispatch therefore uses it only when `wp == 8`; extended kinds retain scalar
TN. The tiled NT path remains enabled for all kinds. Laptop n=1000 timing
improved from about 2.3915 s to 2.2551 s (1.06x). Isolated final timings were
about 51.05 ms for TN and 54.47 ms for NT. NT row blocks 128, 256, and 512 gave
only a noisy two-percent advantage to 512, so the shared value remained 256.

Validation passed `make check` at `wp=8`, `wp=10`, and `wp=16`, plus a release
build. The DGEMM test covers N/T/C, dimensions zero through three, rectangular
fast and fallback paths, padded leading dimensions, alpha/beta 0, 1, -1, and
0.375, output padding, and NaN-poisoned C when beta=0.

## Current algorithms

The `wp=8` TN (`A**T*B`) fast path computes a 2-by-4 C microtile with eight
accumulators. Its K loop reads `A(:,i:i+1)` and `B(:,j:j+3)` contiguously,
reuses two A and four B values, and writes each C element once. Odd M and N
remainders are scalar. There is no packing, explicit SIMD, intrinsic, OpenMP,
or allocation. `wp=10` and `wp=16` use the original scalar dot-product path.

NT (`A*B**T`) first scales/zeros C, groups four output columns, and tiles M
with `ROW_BLOCK=256`. For every K it loads four B scalars and walks contiguously
through one A column and four C columns. It updates C once per K. At `wp=8`,
the principal hot set is roughly 10 KiB. A compiler may store `real(kind=10)`
in 16 bytes, doubling that footprint; measure `storage_size(1.0_wp)/8`.

NN and TT remain essentially the Netlib reference paths, as do scalar tails.
This is intentional because the measured QR workload was dominated by TN and
NT. Re-profile before specializing anything else.

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

## Ranked next work

1. **Measure TN vectorization and counters, then run one `wp=8` TN experiment.**
   This has the largest measured gap and highest likely benefit. If assembly is
   scalar/address-heavy, test a bounded packed panel with a SIMD-friendly
   microkernel. If it is already vectorized but dependency-limited, test only a
   larger register tile. Difficulty is medium; rounding risk is low-to-medium
   because reduction order may change. Preserve the extended-kind fallback.

2. **Tune/specialize NT for small K.** Test NR values matched to the target
   vector/register file and MB values derived from the cache model. K=32 may
   not amortize general packing. Likely benefit is medium, difficulty low-to-
   medium, and numerical risk low if K order stays unchanged.

3. **Re-profile the full QR after every retained change.** Only specialize NN,
   TT, or another kernel if its measured contribution becomes material. This
   measurement is cheap and prevents optimizing a former bottleneck.

4. **Consider OpenMP only after the serial kernels stabilize.** Parallelize
   independent C tiles above a workload threshold, preserve a non-OpenMP build,
   and check for caller oversubscription. This is a build/policy change and
   should not be mixed with the next serial experiment.

## Correctness and retention gate

Every candidate must preserve N/T/C, arbitrary valid leading dimensions,
rectangular and zero-size calls, all alpha/beta cases, and no C read when
beta=0. Run the focused DGEMM test during iteration, then `make check` for all
kinds and a release build. Run complete QR and QR-update/residual tests. Use
K- and epsilon-scaled error bounds rather than bitwise equality. Revert any
candidate that fails correctness, does not improve complete QR reproducibly,
or seriously regresses small matrices.

## Smallest decisive next experiment

Do not modify production DGEMM before collecting a target-node shape histogram,
vectorization report, assembly classification, and counters for TN
968x32x968 and NT 968x968x32 (or the cluster's actual dominant equivalents).
Then change only TN at `wp=8`, choosing between a packed/SIMD-friendly panel and
a larger unpacked register tile according to that evidence. Compare it against
commit `dde28f1` in the isolated TN benchmark and complete QR. This preserves
the proven 2.88x result while decisively testing the leading bottleneck.

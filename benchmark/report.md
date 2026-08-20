# QR optimization and reference comparison

Generated on 2026-08-20 from branch `optimize/qr-hotpaths` at commit
`07f883d`, with the uncommitted OpenMP changes described below.

## Result

Coarse OpenMP parallelism in the QR-dominant GEMM output-column loops is worth
retaining. At order 1000, the four-physical-core build reduced fresh QR from
320.230 to 203.091 ms for real arithmetic (1.58x) and from 1158.649 to
537.331 ms for complex arithmetic (2.16x). Against QR built with the
pre-optimization BLAS, the corresponding end-to-end speedups are 4.16x and
2.84x.

The change does not alter the QR algorithm, block size, explicit dense Q, or
qrupdate kernels. The original Claude LDLT/LDLH implementation remains about
5.2-5.4x faster than four-core QR for fresh factorization at order 1000. That
is not an equal-arithmetic kernel comparison: QR forms an explicit Q, whereas
LDLT/LDLH does not. QR replacement and deletion remain valuable because they
avoid a complete refactorization.

## Measurement configuration

- Host: WSL2 Linux 6.6.87.2, Intel Core i7-8565U, four physical cores and two
  hardware threads per core.
- Cache reported by the host: 32 KiB L1 data, 256 KiB L2, 8 MiB L3.
- Compiler: GNU Fortran 15.2.0, `-O3 -march=native`; OpenMP adds `-fopenmp`.
- Precision for the full comparison: `wp=8`.
- Serial processes were pinned to CPU 0. The retained threaded comparison
  used physical CPUs 0, 2, 4, and 6 with `OMP_PROC_BIND=true`.
- Identical datasets, shifts, starting vectors, tolerances, and validation
  were used for every variant.
- Repetitions per executable were 20, 10, 5, and 3 at orders 100, 200, 500,
  and 1000 respectively. Each value is the executable's mean timed duration.
- Input, allocation, warm-up, state preparation, destructive-input copies,
  and post-operation validation are outside the timed region.

The four variants are optimized serial QR; optimized four-core OpenMP QR; the
same QR linked with `BLAS.f` from pre-optimization commit `701c016`; and the
pristine Claude LDLT/LDLH implementation. `LAPACK.f` has no diff between
`701c016` and the current commit, so the third case is the requested pristine
BLAS/LAPACK baseline without changing the QR algorithm.

## Complete timings

All values are milliseconds per operation.

### Optimized serial QR

| Kind | n | Factorization | Solve | Full solve | Replacement | Append | Deletion |
|---|---:|---:|---:|---:|---:|---:|---:|
| real | 100 | 0.6385 | 0.0640 | 0.6695 | 0.0871 | 0.0325 | 0.0276 |
| real | 200 | 3.4439 | 0.3172 | 3.8922 | 0.4082 | 0.1364 | 0.0807 |
| real | 500 | 41.9219 | 2.9863 | 43.4608 | 2.5018 | 0.8924 | 0.7960 |
| real | 1000 | 320.2297 | 11.0308 | 322.1794 | 10.2080 | 3.8760 | 3.2418 |
| complex | 100 | 1.1708 | 0.1934 | 1.5344 | 0.2077 | 0.1009 | 0.0541 |
| complex | 200 | 8.8908 | 0.8271 | 9.7557 | 0.6540 | 0.2251 | 0.1675 |
| complex | 500 | 131.0781 | 5.9086 | 140.1172 | 4.2225 | 1.6399 | 1.7610 |
| complex | 1000 | 1158.6493 | 24.9758 | 1181.6876 | 18.0916 | 7.2164 | 6.3743 |

### Optimized OpenMP QR, four physical cores

| Kind | n | Factorization | Solve | Full solve | Replacement | Append | Deletion |
|---|---:|---:|---:|---:|---:|---:|---:|
| real | 100 | 0.8652 | 0.0763 | 0.7158 | 0.1084 | 0.0341 | 0.0364 |
| real | 200 | 3.2010 | 0.3401 | 3.6938 | 0.3420 | 0.1204 | 0.0827 |
| real | 500 | 31.1663 | 2.7617 | 32.5081 | 2.4605 | 0.8953 | 0.6492 |
| real | 1000 | 203.0914 | 12.6849 | 212.2671 | 10.3649 | 5.1629 | 3.6344 |
| complex | 100 | 1.4682 | 0.2267 | 1.4610 | 0.1606 | 0.0567 | 0.0473 |
| complex | 200 | 8.9657 | 0.7070 | 9.9047 | 0.6658 | 0.2930 | 0.1532 |
| complex | 500 | 80.5133 | 6.3826 | 87.9513 | 4.1762 | 1.8889 | 1.6647 |
| complex | 1000 | 537.3311 | 26.1655 | 564.4145 | 18.9191 | 8.4315 | 6.6924 |

Only factorization and full solve are intended to benefit: the other QR
operations do not enter a GEMM above the threshold. Their small differences
are run-to-run noise, code-layout effects, or OpenMP-runtime overhead, not
parallel qrupdate speedups. The default serial build remains preferable for
small matrices and update-heavy workloads.

### QR with pristine BLAS/LAPACK

| Kind | n | Factorization | Solve | Full solve | Replacement | Append | Deletion |
|---|---:|---:|---:|---:|---:|---:|---:|
| real | 100 | 0.7100 | 0.0707 | 0.7574 | 0.0957 | 0.0329 | 0.0267 |
| real | 200 | 5.4315 | 0.2762 | 5.8229 | 0.3240 | 0.1053 | 0.0808 |
| real | 500 | 96.8691 | 2.1086 | 97.8644 | 1.9419 | 0.8545 | 0.6650 |
| real | 1000 | 845.7597 | 10.5260 | 871.8703 | 9.8672 | 4.3868 | 3.3822 |
| complex | 100 | 1.3127 | 0.2021 | 1.5733 | 0.1742 | 0.0594 | 0.0409 |
| complex | 200 | 10.2311 | 0.7071 | 10.8800 | 0.5449 | 0.3148 | 0.1559 |
| complex | 500 | 161.3057 | 5.8283 | 156.3947 | 4.4025 | 1.5950 | 1.5936 |
| complex | 1000 | 1523.4973 | 32.0920 | 1486.6888 | 18.5849 | 7.6723 | 6.3962 |

The pristine BLAS changes fresh factorization materially but has little
systematic effect on qrupdate operations. Isolated means can invert slightly,
as in complex order-500 full solve versus factorization, because these short
runs are exploratory means rather than paired confidence intervals.

### Original Claude LDLT/LDLH

| Kind | n | Factorization | Solve | Full solve | Append factorization | Append and solve |
|---|---:|---:|---:|---:|---:|---:|
| real | 100 | 0.1037 | 0.0408 | 0.1428 | 0.0029 | 0.0403 |
| real | 200 | 0.6270 | 0.1959 | 0.9606 | 0.0138 | 0.2303 |
| real | 500 | 5.8281 | 1.3639 | 7.5874 | 0.1317 | 1.4244 |
| real | 1000 | 39.2104 | 6.9885 | 44.5864 | 0.7076 | 7.2073 |
| complex | 100 | 0.1740 | 0.0900 | 0.2289 | 0.0052 | 0.0853 |
| complex | 200 | 1.0586 | 0.3348 | 1.4092 | 0.0244 | 0.3658 |
| complex | 500 | 12.1911 | 2.5054 | 15.1007 | 0.3372 | 3.3336 |
| complex | 1000 | 100.2980 | 12.1554 | 101.8781 | 1.0369 | 13.1541 |

Claude's implementation has no native replacement or deletion operation.
When a physical row/column changes, its applicable baseline is a complete
factorization rebuild, and that distinction must remain explicit.

## End-to-end comparisons

### Fresh factorization speedup

| Kind | n | Serial / pristine | OpenMP4 / serial | OpenMP4 / pristine |
|---|---:|---:|---:|---:|
| real | 100 | 1.11x | 0.74x | 0.82x |
| real | 200 | 1.58x | 1.08x | 1.70x |
| real | 500 | 2.31x | 1.35x | 3.11x |
| real | 1000 | 2.64x | 1.58x | 4.16x |
| complex | 100 | 1.12x | 0.80x | 0.89x |
| complex | 200 | 1.15x | 0.99x | 1.14x |
| complex | 500 | 1.23x | 1.63x | 2.00x |
| complex | 1000 | 1.31x | 2.16x | 2.84x |

Values below 1.0 indicate a regression. The two-million-multiply threshold
keeps the important order-100 and order-200 GEMMs serial. Small-size changes
therefore demonstrate why OpenMP is opt-in; they are not evidence for lowering
the threshold.

At order 1000, Claude LDLT/LDLH is 5.18x faster than four-core real QR and
5.36x faster than four-core complex QR for fresh factorization. For full
factorization plus solve the gaps are 4.76x and 5.54x. Conversely, comparing
QR's native update to the explicitly labelled Claude rebuild baseline:

| Kind | QR operation | OpenMP QR (ms) | Claude rebuild (ms) | QR advantage |
|---|---|---:|---:|---:|
| real | Replacement | 10.3649 | 39.2104 | 3.78x |
| real | Deletion | 3.6344 | 39.2104 | 10.79x |
| complex | Replacement | 18.9191 | 100.2980 | 5.30x |
| complex | Deletion | 6.6924 | 100.2980 | 14.99x |

## OpenMP implementation and threshold

The retained directives distribute independent output columns or four-column
tiles with static scheduling. Each thread owns disjoint columns of C for the
complete K reduction, so no reduction, atomic operation, temporary matrix, or
floating-point reassociation across threads is required. The targeted paths
are real `DGEMM(T,N)` and `DGEMM(N,T)`, and complex `ZGEMM(C,N)` and
`ZGEMM(N,C)`. These are the calls made by blocked `xLARFB` in both GEQRF and
explicit-Q generation.

Parallel work is enabled only when there are at least eight output tasks and
the approximate multiply count is at least 2,000,000. The comparison was
chosen from actual QR calls: order 200 peaks below about 0.9 million
multiplies, while order 500 and 1000 contain much larger trailing updates.
There is no runtime autotuning or hardware detection.

| Phase | Real call | Complex call | M | N | K |
|---|---|---|---:|---:|---:|
| GEQRF trailing update | `DGEMM(T,N)` | `ZGEMM(C,N)` | 968 | 32 | 968 |
| Explicit Q | `DGEMM(N,T)` | `ZGEMM(N,C)` | 968 | 968 | 32 |

## Thread scaling

A separate pinned order-1000 scaling run produced:

| Kind | Threads | Time (ms) | Speedup vs OpenMP-1 | Efficiency |
|---|---:|---:|---:|---:|
| real | 1 | 324.997 | 1.00x | 100% |
| real | 2 | 245.580 | 1.32x | 66% |
| real | 4 | 194.790 | 1.67x | 42% |
| real | 8 SMT | 193.810 | 1.68x | 21% |
| complex | 1 | 1194.250 | 1.00x | 100% |
| complex | 2 | 742.760 | 1.61x | 80% |
| complex | 4 | 638.870 | 1.87x | 47% |
| complex | 8 SMT | 563.060 | 2.12x | 27% |

Eight threads use SMT, not eight physical cores. The SMT result was not
robust: another stage-level run regressed from 193 to 325 ms real and from 706
to 774 ms complex. Four physical threads are therefore the defensible setting
on this host. Frequency and host-load variation in WSL2 also explains the
absolute-time difference between this scaling run and the complete comparison.

Representative extended-precision stage measurements at order 500 were:

| Precision/kind | Serial stages (ms) | OpenMP-1 (ms) | OpenMP-2 (ms) | OpenMP-4 (ms) | OpenMP-4 speedup |
|---|---:|---:|---:|---:|---:|
| wp=10 real | 209.86 | 220.51 | 145.47 | 123.76 | 1.78x |
| wp=10 complex | 750.82 | 794.00 | 492.61 | 409.77 | 1.94x |
| wp=16 real | 7545 | 7831 | 5712 | 4682 | 1.67x |

The wp=16 complex order-500 timing was stopped after exceeding 30 seconds per
calibration run. Its threaded CN/NC paths are covered by the all-precision
correctness suite, but no performance claim is made for that expensive case.

## Hardware counters and post-thread profile

The following `perf stat` data use order 1000. Counts are normalized per
factorization; IPC and miss rates are aggregate task values that include all
OpenMP workers.

| Kind/build | Time (ms) | Instructions | Cycles | IPC | Cache misses | Cache miss rate | Branch miss rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| real OpenMP-1 | 304.169 | 3.025e9 | 1.453e9 | 2.08 | 22.95e6 | 20.03% | 1.13% |
| real OpenMP-4 | 176.493 | 3.333e9 | 2.606e9 | 1.28 | 11.39e6 | 13.02% | 0.60% |
| complex OpenMP-1 | 1152.561 | 10.473e9 | 6.108e9 | 1.71 | 91.29e6 | 17.16% | 0.42% |
| complex OpenMP-4 | 525.580 | 11.035e9 | 8.577e9 | 1.29 | 46.17e6 | 8.72% | 0.31% |

Aggregate cycles increase and IPC falls because worker cycles and OpenMP
synchronization are summed across threads; those values must not be read as
single-core execution efficiency. Wall time is the scaling metric. The
roughly halved cache-miss count is consistent with the static, disjoint-column
ownership improving concurrent locality rather than creating shared writes.

| Kind | Serial top symbols | Four-thread top symbols | Classification |
|---|---|---|---|
| real | DGEMM 85.44%, DGEMV 9.56%, DTRMM 1.48% | outlined DGEMM loops 60.57%, DGEMV 5.55%, DTRMM 1.80%; unresolved libgomp wait/synchronization samples about 28% | still BLAS-3 dominated, now limited partly by parallel overhead and serial panel work |
| complex | ZGEMM 90.24%, ZGEMV 3.99%, ZTRMM 1.99% | outlined ZGEMM loops 73.16%, ZGEMV 2.40%, ZTRMM 1.97%; unresolved libgomp wait/synchronization samples about 18% | still BLAS-3 dominated, with stronger useful scaling than real |

TRMM remains about 2% after threading. Parallelizing it cannot materially
change end-to-end time and would add another threshold and parallel region, so
the optimization stops at GEMM. The next bottleneck is the combination of
remaining GEMM work, OpenMP synchronization, and serial panel/DGEMV work, not
a newly dominant scalar kernel.

## Correctness and rejected extensions

Every full-comparison run returned status 0. Inverse iteration took four to
eight iterations and the largest observed normalized residual was about
`9.1e-18`. The serial and OpenMP debug suites passed for `wp=8`, `wp=10`, and
`wp=16`, including new above-threshold real TN/NT and complex CN/NC GEMM cases
with padded leading dimensions. Optimized release builds succeeded for serial
wp=8 and OpenMP wp=8/10/16.

The following extensions were deliberately rejected or deferred:

- no threading of small GEMM calls, because order 100/200 do not amortize it;
- no TRMM threading, because its measured share is about 2%;
- no qrupdate threading, because it is outside the measured fresh-QR
  bottleneck and one mutable state remains non-thread-safe;
- no eight-thread default on a four-core CPU, because SMT scaling was unstable;
- no runtime autotuning, CPU detection, QR redesign, block-size change, or
  further BLAS microkernel work.

## Files changed

- `src/qrupdate/BLAS.f`: guarded OpenMP worksharing on the four QR-dominant
  GEMM paths.
- `Makefile`: opt-in `OPENMP=1` flags and separate `-omp` build directories.
- `benchmark/build.sh`, `benchmark/run_suite.sh`, `benchmark/run_perf.sh`, and
  `benchmark/dgemm/build.sh`: opt-in threaded builds and result locations.
- `test/test_dgemm.f90` and `test/test_zgemm.f90`: above-threshold QR-shaped
  semantic tests.
- `README.md`, `benchmark/README.md`, and `docs/DGEMM_HPC_HANDOFF.md`: usage,
  affinity, threshold, and measurement notes.

## Limitations

- These are short exploratory means on a WSL2 laptop, not confidence
  intervals from an isolated HPC node. Re-run the documented protocol after
  moving machines.
- The complete four-way comparison is wp=8 so all implementations use the
  same binary64 arithmetic and input data. wp=10 and wp=16 are validated and
  represented by bounded stage measurements, not an exhaustive all-mode run.
- QR fresh factorization includes shifted-matrix formation and explicit Q;
  Claude factorization receives a prepared shifted matrix and produces no Q.
- The OpenMP build changes only fresh-factorization GEMMs. It is not expected
  to accelerate solves or structural updates.
- Dataset files store binary64 values even when a wider working kind is used.

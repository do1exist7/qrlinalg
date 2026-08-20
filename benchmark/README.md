# Standalone benchmarks

The benchmark subsystem builds and runs independently of the project
Makefile. QR and the pristine serial LDLT/LDLH implementation use separate
executables and result files. Every operation prints a readable report and a
final `BENCH_RESULT` CSV record.

All matrices come from `data/data_<n>/`. The inverse-iteration shifts,
tolerance, iteration limit, normalization, and all-ones starting vector match
`orig/claude/main.f90`. File input, allocation, state preparation, copies of
destructive LDLT inputs, warm-up, and residual validation are not timed.

## Build and run

```sh
./benchmark/build.sh all --precision 8
./benchmark/run_qr.sh --sizes 100,200,500,1000 --repetitions 5
./benchmark/run_ldlt.sh --sizes 100,200,500,1000 --repetitions 5
./benchmark/run_qr.sh --kind complex --operation replacement --repetitions 100
./benchmark/run_ldlt.sh --kind real --operation append_and_solve

# Build/run the opt-in OpenMP QR variant in a separate qr-openmp directory.
OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./benchmark/run_qr.sh --openmp --operation factorization \
  --sizes 500,1000 --repetitions 5
```

The compiler defaults to `gfortran`; `--compiler` also accepts `ifort`, `ifx`,
and `nvfortran`. Use `--repetitions auto` for the size- and precision-dependent
scaling from the original driver, and `--no-build` to reuse executables.
`--openmp` is available for QR builds and runs only; it keeps serial and
threaded benchmark objects separate. Set `OMP_NUM_THREADS` explicitly and use
`OMP_PROC_BIND`/`OMP_PLACES` appropriate to the measured machine.

QR operations are `factorization`, `solve`, `full_solve`, `replacement`,
`append`, and `deletion`. LDLT operations are `factorization`, `solve`,
`full_solve`, `append_factorization`, and `append_and_solve`.

`solve` reuses complete factors. `full_solve` includes factorization and
inverse iteration. `append_factorization` measures only `LDLTF(n,n)` or
`LDLHF(n,n)` after priming the leading block. `append_and_solve` preserves the
original `GSEPIIS`/`GHEPIIS` path with `k=n`.

QR append and deletion copy a prepared state outside the timed region and
repeat the same final-index operation. Replacement alternates a deterministic
data-derived perturbation and its negative. Every sample therefore has the
same active order, and no inverse structural operation is timed.

LDLT has no native replacement or deletion operation. Comparisons may use a
full LDLT factorization as a rebuild baseline, but must label it explicitly.

## Visualization

```sh
./benchmark/plot_timings.py \
  build/benchmarks/wp8/results/qr.csv \
  build/benchmarks/wp8/results/ldlt.csv \
  --output build/benchmarks/wp8/plots

./benchmark/generate_report.py \
  --qr build/benchmarks/wp8/results/qr.csv \
  --ldlt build/benchmarks/wp8/results/ldlt.csv \
  --plots build/benchmarks/wp8/plots \
  --perf-root build/benchmarks/wp8/perf \
  --output build/benchmarks/wp8/report.md
```

The script creates independent timing plots and prints all eigenvalues,
iteration counts, residuals, and statuses. It requires Matplotlib.

The tracked [`report.md`](report.md) records the 2026-08-20 comparison of the
optimized serial and OpenMP QR builds, QR with the pre-optimization bundled
BLAS/LAPACK, and the original Claude LDLT/LDLH implementation across every
operation supported by the benchmark drivers.

## Perf

```sh
./benchmark/run_perf.sh \
  --implementation qr --operation replacement --kind real --size 1000 \
  --kernel-repetitions 100 --perf-repetitions 5

OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores \
  ./benchmark/run_perf.sh --openmp --implementation qr \
  --operation factorization --kind real --size 1000 \
  --kernel-repetitions 5 --perf-repetitions 3

./benchmark/plot_perf.py \
  build/benchmarks/wp8/perf/qr/replacement/real_1000.csv \
  --output build/benchmarks/wp8/plots/qr-replacement-perf.png
```

The perf script starts counters disabled and uses perf control FIFOs to enable
them at the same internal boundaries used for wall-clock timing. Dataset input,
state copies, warm-up, and validation are therefore excluded from both timing
and counters. The adjacent `.result.txt` retains the verbose numerical result.
The plotting and report scripts read that result to normalize aggregate perf
counters to one numerical operation.

## DGEMM and QR-stage benchmarks

The focused DGEMM build provides an isolated matrix-product driver, separate
`DGEQRF` and `DORGQR` timings, and an optional wp=8/gfortran link-time call-shape
profiler.  Each timing driver performs untimed warm-up and reports the median
of independently timed samples.  Choose the calls per sample so every sample
lasts at least 0.5 seconds on the machine being measured.

```sh
./benchmark/dgemm/build.sh --precision 8
./benchmark/dgemm/build.sh --precision 8 --openmp

# Dominant order-1000 DLARFB shapes measured in the real QR pipeline.
taskset -c 0 build/benchmarks/wp8/dgemm/bin/dgemm_benchmark \
  T N 968 32 968 1 1 100 7
taskset -c 0 build/benchmarks/wp8/dgemm/bin/dgemm_benchmark \
  N T 968 968 32 -1 1 100 7

# Separate QR stages; copies that restore the input are not timed.
taskset -c 0 build/benchmarks/wp8/dgemm/bin/qr_stage_benchmark 1000 3 7

# Explicitly enabled diagnostic; aggregated shapes are written to stderr.
build/benchmarks/wp8/dgemm/bin/qr_stage_profile 1000 1 1 \
  2>build/benchmarks/wp8/dgemm/call-shapes.csv
```

The profile wrapper is never linked into the library or normal benchmarks. Its
C interface is intentionally restricted to the gfortran ABI at `wp=8`; the
numerical timing drivers themselves support every compiler and precision
accepted by the main benchmark suite.

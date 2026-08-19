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
```

The compiler defaults to `gfortran`; `--compiler` also accepts `ifort`, `ifx`,
and `nvfortran`. Use `--repetitions auto` for the size- and precision-dependent
scaling from the original driver, and `--no-build` to reuse executables.

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

## Perf

```sh
./benchmark/run_perf.sh \
  --implementation qr --operation replacement --kind real --size 1000 \
  --kernel-repetitions 100 --perf-repetitions 5

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

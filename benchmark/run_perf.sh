#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
implementation=''
operation=''
kind=real
size=1000
precision=8
compiler=gfortran
kernel_repetitions=100
perf_repetitions=5
output=''
rebuild=1

usage() {
  cat <<'EOF'
Usage: benchmark/run_perf.sh --implementation qr|ldlt --operation NAME [options]
  --kind real|complex
  --size N
  --precision 8|10|16
  --compiler NAME
  --kernel-repetitions N
  --perf-repetitions N
  --output FILE
  --no-build
EOF
}

while (($#)); do
  case "$1" in
    --implementation) implementation=$2; shift 2 ;;
    --operation) operation=$2; shift 2 ;;
    --kind) kind=$2; shift 2 ;;
    --size) size=$2; shift 2 ;;
    --precision) precision=$2; shift 2 ;;
    --compiler) compiler=$2; shift 2 ;;
    --kernel-repetitions) kernel_repetitions=$2; shift 2 ;;
    --perf-repetitions) perf_repetitions=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    --no-build) rebuild=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ "$implementation" == qr || "$implementation" == ldlt ]] || { usage >&2; exit 2; }
[[ -n "$operation" ]] || { usage >&2; exit 2; }
[[ "$kind" == real || "$kind" == complex ]] || { usage >&2; exit 2; }
[[ "$size" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid size" >&2; exit 2; }
[[ "$kernel_repetitions" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid kernel repetitions" >&2; exit 2; }
[[ "$perf_repetitions" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid perf repetitions" >&2; exit 2; }
command -v perf >/dev/null || { echo "perf is not installed" >&2; exit 1; }

if ((rebuild)); then
  "$repository_dir/benchmark/build.sh" "$implementation" \
    --precision "$precision" --compiler "$compiler"
fi

executable="$repository_dir/build/benchmarks/wp${precision}/${implementation}/bin/${implementation}_${operation}"
h_file="$repository_dir/data/data_${size}/H_${kind}.dat"
s_file="$repository_dir/data/data_${size}/S_${kind}.dat"
[[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
[[ -f "$h_file" && -f "$s_file" ]] || { echo "Missing dataset in data/data_${size}" >&2; exit 1; }

if [[ -z "$output" ]]; then
  output="$repository_dir/build/benchmarks/wp${precision}/perf/${implementation}/${operation}/${kind}_${size}.csv"
fi
mkdir -p "$(dirname "$output")"
result_output="${output%.csv}.result.txt"

echo "Profiling implementation=$implementation operation=$operation kind=$kind n=$size"
echo "Kernel repetitions=$kernel_repetitions, perf samples=$perf_repetitions"
echo "Counters are enabled only across the driver's internal timed operation."

control_dir=$(mktemp -d)
control_fifo="$control_dir/control.fifo"
ack_fifo="$control_dir/ack.fifo"
mkfifo "$control_fifo" "$ack_fifo"
# Holding both FIFOs open prevents either perf or the child process from
# blocking while they independently open their respective endpoints.
exec 8<>"$control_fifo"
exec 9<>"$ack_fifo"
export QRLINALG_PERF_CONTROL="$control_fifo"
export QRLINALG_PERF_ACK="$ack_fifo"
cleanup_control() {
  exec 8>&-
  exec 9>&-
  rm -f "$control_fifo" "$ack_fifo"
  rmdir "$control_dir"
}
trap cleanup_control EXIT

perf stat -x, -r "$perf_repetitions" --delay=-1 \
  --control "fifo:$control_fifo,$ack_fifo" \
  -e task-clock,cycles,instructions,branches,branch-misses,cache-references,cache-misses,page-faults,context-switches \
  -o "$output" \
  "$executable" "$kind" "$h_file" "$s_file" "$kernel_repetitions" | tee "$result_output"

unset QRLINALG_PERF_CONTROL QRLINALG_PERF_ACK

rg -q '^BENCH_RESULT,' "$result_output" || { echo "No BENCH_RESULT was produced" >&2; exit 1; }
status=$(awk -F, '/^BENCH_RESULT,/{gsub(/ /,"",$11); value=$11} END{print value}' "$result_output")
[[ "$status" == 0 ]] || { echo "Benchmark returned status $status" >&2; exit 1; }
metrics_valid=$(awk -F, '/^BENCH_RESULT,/{
  iterations=$12+0; residual=$15+0
  invalid=($13 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $15 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/)
  value=(!invalid && iterations >= 1 && residual >= 0 && residual <= 1.0e-8)
} END{print value}' "$result_output")
[[ "$metrics_valid" == 1 ]] || { echo "Invalid convergence metrics" >&2; exit 1; }
echo "Perf counters written to $output"
echo "Verbose numerical result written to $result_output"

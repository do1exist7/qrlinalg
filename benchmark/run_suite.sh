#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
implementation=${1:?internal error: implementation is required}
shift
precision=8
compiler=gfortran
sizes=100,200,500,1000
kinds=real,complex
operations=all
repetitions=5
output=''
rebuild=1

usage() {
  cat <<EOF
Usage: benchmark/run_${implementation}.sh [options]
  --precision 8|10|16
  --compiler NAME
  --sizes 100,200,500,1000
  --kind real|complex|both
  --operation NAME|all
  --repetitions INTEGER|auto
  --output FILE
  --no-build
EOF
}

while (($#)); do
  case "$1" in
    --precision) precision=$2; shift 2 ;;
    --compiler) compiler=$2; shift 2 ;;
    --sizes) sizes=$2; shift 2 ;;
    --kind)
      case "$2" in real|complex) kinds=$2 ;; both) kinds=real,complex ;; *) usage >&2; exit 2 ;; esac
      shift 2
      ;;
    --operation) operations=$2; shift 2 ;;
    --repetitions) repetitions=$2; shift 2 ;;
    --output) output=$2; shift 2 ;;
    --no-build) rebuild=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$precision" in 8|10|16) ;; *) echo "Precision must be 8, 10, or 16" >&2; exit 2 ;; esac
if [[ "$repetitions" != auto && ! "$repetitions" =~ ^[1-9][0-9]*$ ]]; then
  echo "Repetitions must be a positive integer or auto" >&2
  exit 2
fi

if [[ "$implementation" == qr ]]; then
  available=(factorization solve full_solve replacement append deletion)
else
  available=(factorization solve full_solve append_factorization append_and_solve)
fi
if [[ "$operations" == all ]]; then
  selected_operations=("${available[@]}")
else
  IFS=',' read -r -a selected_operations <<< "$operations"
  for operation in "${selected_operations[@]}"; do
    [[ " ${available[*]} " == *" $operation "* ]] || {
      echo "Unsupported $implementation operation: $operation" >&2
      echo "Available: ${available[*]}" >&2
      exit 2
    }
  done
fi

if ((rebuild)); then
  "$repository_dir/benchmark/build.sh" "$implementation" \
    --precision "$precision" --compiler "$compiler"
fi

if [[ -z "$output" ]]; then
  output="$repository_dir/build/benchmarks/wp${precision}/results/${implementation}.csv"
fi
mkdir -p "$(dirname "$output")"
header='record,implementation,operation,kind,precision,n,shift,repetitions,total_seconds,seconds_per_operation,status,iterations,eigenvalue,relative_accuracy,residual'
printf '%s\n' "$header" > "$output"

automatic_repetitions() {
  local operation=$1 n=$2 base full update
  case "$precision" in 8) base=48 ;; 10) base=16 ;; 16) base=4 ;; esac
  full=$((base * 1000000000 / (n * n * n)))
  ((full < 1)) && full=1
  update=$((full * n / 40))
  ((update < 1)) && update=1
  case "$operation" in factorization|full_solve) echo "$full" ;; *) echo "$update" ;; esac
}

IFS=',' read -r -a selected_sizes <<< "$sizes"
IFS=',' read -r -a selected_kinds <<< "$kinds"
for n in "${selected_sizes[@]}"; do
  [[ "$n" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid size: $n" >&2; exit 2; }
  for kind in "${selected_kinds[@]}"; do
    h_file="$repository_dir/data/data_${n}/H_${kind}.dat"
    s_file="$repository_dir/data/data_${n}/S_${kind}.dat"
    [[ -f "$h_file" && -f "$s_file" ]] || {
      echo "Missing dataset for $kind order $n below data/" >&2
      exit 1
    }
    for operation in "${selected_operations[@]}"; do
      run_repetitions=$repetitions
      if [[ "$run_repetitions" == auto ]]; then
        run_repetitions=$(automatic_repetitions "$operation" "$n")
      fi
      executable="$repository_dir/build/benchmarks/wp${precision}/${implementation}/bin/${implementation}_${operation}"
      [[ -x "$executable" ]] || { echo "Missing executable: $executable" >&2; exit 1; }
      echo
      echo "[$implementation] $operation, kind=$kind, n=$n, repetitions=$run_repetitions"
      echo "  H: $h_file"
      echo "  S: $s_file"
      run_output=$("$executable" "$kind" "$h_file" "$s_file" "$run_repetitions")
      printf '%s\n' "$run_output"
      result_line=$(printf '%s\n' "$run_output" | awk '/^BENCH_RESULT,/{line=$0} END{print line}')
      [[ -n "$result_line" ]] || { echo "Executable produced no BENCH_RESULT" >&2; exit 1; }
      status=$(printf '%s\n' "$result_line" | awk -F, '{gsub(/ /,"",$11); print $11}')
      [[ "$status" == 0 ]] || { echo "Benchmark returned status $status" >&2; exit 1; }
      metrics_valid=$(printf '%s\n' "$result_line" | awk -F, '
        {
          iterations=$12+0; residual=$15+0
          invalid=($13 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/ || $15 ~ /[Nn][Aa][Nn]|[Ii][Nn][Ff]/)
          if (!invalid && iterations >= 1 && residual >= 0 && residual <= 1.0e-8) print 1
          else print 0
        }')
      [[ "$metrics_valid" == 1 ]] || {
        echo "Benchmark produced invalid convergence metrics" >&2
        exit 1
      }
      printf '%s\n' "$result_line" >> "$output"
    done
  done
done

echo
echo "Results written to $output"

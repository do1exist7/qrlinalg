#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
precision=8
compiler=gfortran
openmp=0

usage() {
  echo "Usage: $0 [--precision 8|10|16] [--compiler NAME] [--openmp]"
}

while (($#)); do
  case "$1" in
    --precision) precision=$2; shift 2 ;;
    --compiler) compiler=$2; shift 2 ;;
    --openmp) openmp=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$precision" in 8|10|16) ;; *) echo "Precision must be 8, 10, or 16" >&2; exit 2 ;; esac
command -v "$compiler" >/dev/null || { echo "Compiler not found: $compiler" >&2; exit 1; }

case "$compiler" in
  gfortran)
    optimization_flags=(-O3 -march=native)
    fixed_flags=(-ffixed-line-length-none)
    preprocess_flags=(-cpp "-DQRLINALG_WP=$precision")
    module_flag=-J
    openmp_flag=-fopenmp
    ;;
  ifort)
    optimization_flags=(-O3 -ip -fp-model precise)
    fixed_flags=(-extend-source)
    preprocess_flags=(-fpp "-DQRLINALG_WP=$precision")
    module_flag=-module
    openmp_flag=-qopenmp
    ;;
  ifx)
    optimization_flags=(-O3 -fp-model precise)
    fixed_flags=(-extend-source)
    preprocess_flags=(-fpp "-DQRLINALG_WP=$precision")
    module_flag=-module
    openmp_flag=-qopenmp
    ;;
  nvfortran)
    optimization_flags=(-O3 -tp=native)
    fixed_flags=(-Mextend)
    preprocess_flags=(-Mpreprocess "-DQRLINALG_WP=$precision")
    module_flag=-module
    openmp_flag=-mp
    ;;
  *) echo "Unsupported compiler: $compiler" >&2; exit 2 ;;
esac

variant_suffix=''
if ((openmp)); then
  optimization_flags+=("$openmp_flag")
  variant_suffix=-openmp
fi

output_dir="$repository_dir/build/benchmarks/wp${precision}/dgemm$variant_suffix"
module_dir="$output_dir/mod"
object_dir="$output_dir/obj"
binary_dir="$output_dir/bin"
mkdir -p "$module_dir" "$object_dir" "$binary_dir"
module_flags=("$module_flag" "$module_dir" -I "$module_dir")

"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "${preprocess_flags[@]}" -c "$repository_dir/src/wp_def_${precision}.f90" \
  -o "$object_dir/wp_def.o"
"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "${fixed_flags[@]}" "${preprocess_flags[@]}" \
  -c "$repository_dir/src/qrupdate/BLAS.f" \
  -o "$object_dir/blas.o"
"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "${fixed_flags[@]}" -c "$repository_dir/src/qrupdate/LAPACK.f" \
  -o "$object_dir/lapack.o"

"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "$repository_dir/benchmark/dgemm/dgemm_benchmark.f90" \
  "$object_dir/wp_def.o" "$object_dir/blas.o" \
  -o "$binary_dir/dgemm_benchmark"
"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "$repository_dir/benchmark/dgemm/zgemm_benchmark.f90" \
  "$object_dir/wp_def.o" "$object_dir/blas.o" \
  -o "$binary_dir/zgemm_benchmark"
"$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
  "$repository_dir/benchmark/dgemm/qr_stage_benchmark.f90" \
  "$object_dir/wp_def.o" "$object_dir/blas.o" "$object_dir/lapack.o" \
  -o "$binary_dir/qr_stage_benchmark"

# The call-shape wrapper follows the gfortran wp=8 ABI and remains outside the
# library.  It is linked only into this explicitly named diagnostic executable.
if [[ "$compiler" == gfortran && "$precision" == 8 ]] && command -v gcc >/dev/null; then
  gcc -O2 -Wall -Wextra -c \
    "$repository_dir/benchmark/dgemm/dgemm_profile_wrapper.c" \
    -o "$object_dir/dgemm_profile_wrapper.o"
  "$compiler" "${optimization_flags[@]}" "${module_flags[@]}" \
    "$repository_dir/benchmark/dgemm/qr_stage_benchmark.f90" \
    "$object_dir/wp_def.o" "$object_dir/blas.o" "$object_dir/lapack.o" \
    "$object_dir/dgemm_profile_wrapper.o" -Wl,--wrap=dgemm_ \
    -o "$binary_dir/qr_stage_profile"
fi

echo "Built DGEMM benchmarks in $binary_dir"

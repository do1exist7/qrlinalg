#!/usr/bin/env bash
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
implementation=${1:-all}
shift || true
precision=8
compiler=gfortran
clean=0

usage() {
  echo "Usage: $0 [qr|ldlt|all] [--precision 8|10|16] [--compiler NAME] [--clean]"
}

while (($#)); do
  case "$1" in
    --precision) precision=$2; shift 2 ;;
    --compiler) compiler=$2; shift 2 ;;
    --clean) clean=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$implementation" in qr|ldlt|all) ;; *) usage >&2; exit 2 ;; esac
case "$precision" in 8|10|16) ;; *) echo "Precision must be 8, 10, or 16" >&2; exit 2 ;; esac
command -v "$compiler" >/dev/null || { echo "Compiler not found: $compiler" >&2; exit 1; }

build_root="$repository_dir/build/benchmarks/wp${precision}"
if ((clean)); then
  case "$implementation" in
    qr|ldlt) rm -rf "$build_root/$implementation" ;;
    all) rm -rf "$build_root/qr" "$build_root/ldlt" ;;
  esac
fi

case "$compiler" in
  gfortran)
    opt_flags=(-O3 -g -march=native)
    warning_flags=(-Wall -Wextra -Wno-unused-dummy-argument)
    fixed_flags=(-ffixed-line-length-none)
    preprocess_flags=(-cpp "-DQRLINALG_WP=$precision")
    module_output_flag=J
    ;;
  ifort)
    opt_flags=(-O3 -ip -fp-model precise)
    warning_flags=(-warn all)
    fixed_flags=(-extend-source)
    preprocess_flags=(-fpp "-DQRLINALG_WP=$precision")
    module_output_flag=module
    ;;
  ifx)
    opt_flags=(-O3 -fp-model precise)
    warning_flags=(-warn all)
    fixed_flags=(-extend-source)
    preprocess_flags=(-fpp "-DQRLINALG_WP=$precision")
    module_output_flag=module
    ;;
  nvfortran)
    opt_flags=(-O3 -tp=native)
    warning_flags=(-Minform=warn)
    fixed_flags=(-Mextend)
    preprocess_flags=(-Mpreprocess "-DQRLINALG_WP=$precision")
    module_output_flag=module
    ;;
  *) echo "Unsupported compiler: $compiler" >&2; exit 2 ;;
esac

compile_qr() {
  local out="$build_root/qr" mod="$build_root/qr/mod" obj="$build_root/qr/obj"
  mkdir -p "$out/bin" "$mod" "$obj"
  local mf=("-$module_output_flag" "$mod" -I "$mod")
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${preprocess_flags[@]}" -c \
    "$repository_dir/src/wp_def_${precision}.f90" -o "$obj/wp_def.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${fixed_flags[@]}" -c \
    "$repository_dir/src/qrupdate/BLAS.f" -o "$obj/blas.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${fixed_flags[@]}" -c \
    "$repository_dir/src/qrupdate/LAPACK.f" -o "$obj/lapack.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/src/qrupdate/qrupdate_linalg.f90" -o "$obj/qrupdate_linalg.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/src/qrupdate/qrupdate_error.f90" -o "$obj/qrupdate_error.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/src/qrupdate/qrupdate_real.f90" -o "$obj/qrupdate_real.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/src/qrupdate/qrupdate_complex.f90" -o "$obj/qrupdate_complex.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/src/qrupdate/qrupdate.f90" -o "$obj/qrupdate.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" \
    "${preprocess_flags[@]}" -c "$repository_dir/src/qrlinalg.f90" -o "$obj/qrlinalg.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" -c \
    "$repository_dir/benchmark/common/benchmark_support.f90" -o "$obj/benchmark_support.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" -c \
    "$repository_dir/benchmark/qr/qr_benchmark_kernels.f90" -o "$obj/qr_kernels.o"

  local common_objects=("$obj/wp_def.o" "$obj/blas.o" "$obj/lapack.o" \
    "$obj/qrupdate_linalg.o" "$obj/qrupdate_error.o" "$obj/qrupdate_real.o" \
    "$obj/qrupdate_complex.o" "$obj/qrupdate.o" "$obj/qrlinalg.o" \
    "$obj/benchmark_support.o" "$obj/qr_kernels.o")
  local operation
  for operation in factorization solve full_solve replacement append deletion; do
    "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" \
      "$repository_dir/benchmark/qr/${operation}.f90" "${common_objects[@]}" \
      -o "$out/bin/qr_${operation}"
  done
  echo "Built QR benchmarks in $out/bin"
}

compile_ldlt() {
  local out="$build_root/ldlt" mod="$build_root/ldlt/mod" obj="$build_root/ldlt/obj"
  mkdir -p "$out/bin" "$mod" "$obj"
  local mf=("-$module_output_flag" "$mod" -I "$mod")
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/test/differential/mpi_serial.f90" -o "$obj/mpi_serial.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/orig/claude/wp_def_${precision}.f90" -o "$obj/wp_def.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/orig/claude/globvars.f90" -o "$obj/globvars.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${fixed_flags[@]}" -c \
    "$repository_dir/src/qrupdate/BLAS.f" -o "$obj/blas.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" -c \
    "$repository_dir/orig/claude/linalg.f90" -o "$obj/linalg.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" -c \
    "$repository_dir/benchmark/common/benchmark_support.f90" -o "$obj/benchmark_support.o"
  "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" -c \
    "$repository_dir/benchmark/ldlt/ldlt_benchmark_kernels.f90" -o "$obj/ldlt_kernels.o"

  local common_objects=("$obj/mpi_serial.o" "$obj/wp_def.o" "$obj/globvars.o" \
    "$obj/blas.o" "$obj/linalg.o" "$obj/benchmark_support.o" "$obj/ldlt_kernels.o")
  local operation
  for operation in factorization solve full_solve append_factorization append_and_solve; do
    "$compiler" "${opt_flags[@]}" "${mf[@]}" "${warning_flags[@]}" \
      "$repository_dir/benchmark/ldlt/${operation}.f90" "${common_objects[@]}" \
      -o "$out/bin/ldlt_${operation}"
  done
  echo "Built LDLT benchmarks in $out/bin"
}

case "$implementation" in
  qr) compile_qr ;;
  ldlt) compile_ldlt ;;
  all) compile_qr; compile_ldlt ;;
esac

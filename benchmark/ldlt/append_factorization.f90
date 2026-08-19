program ldlt_append_factorization_benchmark
  use ldlt_benchmark_kernels, only: run_ldlt_append_factorization
  implicit none
  call run_ldlt_append_factorization()
end program ldlt_append_factorization_benchmark

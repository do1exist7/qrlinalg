program ldlt_factorization_benchmark
  use ldlt_benchmark_kernels, only: run_ldlt_factorization
  implicit none
  call run_ldlt_factorization()
end program ldlt_factorization_benchmark

program ldlt_full_solve_benchmark
  use ldlt_benchmark_kernels, only: run_ldlt_full_solve
  implicit none
  call run_ldlt_full_solve()
end program ldlt_full_solve_benchmark

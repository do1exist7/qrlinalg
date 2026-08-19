program ldlt_solve_benchmark
  use ldlt_benchmark_kernels, only: run_ldlt_solve
  implicit none
  call run_ldlt_solve()
end program ldlt_solve_benchmark

program ldlt_append_and_solve_benchmark
  use ldlt_benchmark_kernels, only: run_ldlt_append_and_solve
  implicit none
  call run_ldlt_append_and_solve()
end program ldlt_append_and_solve_benchmark

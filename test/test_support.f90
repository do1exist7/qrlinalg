! Shared reporting and numerical utilities for the qrlinalg test programs.
!
! Keeping these small facilities in one module lets each numerical area use
! an independent executable without duplicating assertion or norm logic.
module test_support
  use qrlinalg, only: wp
  implicit none
  private

  public :: check
  public :: finish_suite
  public :: real_frobenius
  public :: complex_frobenius

contains

  ! Record a failed condition without aborting immediately, allowing one test
  ! program to report every violated invariant in its suite.
  subroutine check(condition, message, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write(*,'(a)') '  failed: ' // message
    end if
  end subroutine check

  ! Print the result of one test executable and return a failing process status
  ! to Make when any assertion in that executable failed.
  subroutine finish_suite(suite_name, failures)
    character(len=*), intent(in) :: suite_name
    integer, intent(in) :: failures

    if (failures /= 0) then
      write(*,'(a,a,a,i0,a,i0)') 'FAIL: ', trim(suite_name), ', wp=', wp, &
        ', failed checks=', failures
      error stop 1
    end if
    write(*,'(a,a,a,i0)') 'PASS: ', trim(suite_name), ', wp=', wp
  end subroutine finish_suite

  ! Frobenius norm of a real matrix, evaluated entirely in the selected
  ! working kind so the same reference calculation is valid for every build.
  function real_frobenius(a) result(norm)
    real(wp), intent(in) :: a(:,:)
    real(wp) :: norm

    norm = sqrt(sum(a * a))
  end function real_frobenius

  ! Frobenius norm of a complex matrix. ABS produces real(wp), avoiding any
  ! fixed-double dependency in the analytical reference calculation.
  function complex_frobenius(a) result(norm)
    complex(wp), intent(in) :: a(:,:)
    real(wp) :: norm

    norm = sqrt(sum(abs(a)**2))
  end function complex_frobenius

end module test_support

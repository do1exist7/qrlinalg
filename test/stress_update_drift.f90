! Deterministic stress measurements for numerical drift in structural QR updates.
!
! Every operation pair returns the represented matrix to its initial value.
! Consequently, differences between the untouched fresh factors and the
! repeatedly updated factors measure accumulated update error rather than a
! changing matrix or eigenproblem.
program stress_update_drift
  use qrlinalg
  implicit none

  integer, parameter :: default_cycles = 10000
  integer :: max_cycles, parse_status
  character(len=64) :: argument

  max_cycles = default_cycles
  if (command_argument_count() >= 1) then
    call get_command_argument(1, argument)
    read(argument, *, iostat=parse_status) max_cycles
    if (parse_status /= 0 .or. max_cycles <= 0) then
      write(*,'(a)') 'stress_update_drift: cycle count must be a positive integer'
      error stop 1
    end if
  end if

  write(*,'(a)') 'scenario,arithmetic,wp,cycles,public_updates,' // &
    'update_scale,factor_residual,orthogonality,forward_operator,' // &
    'residual_operator,max_backward,fresh_factor_residual,' // &
    'fresh_orthogonality,fresh_forward_operator,' // &
    'fresh_residual_operator,fresh_max_backward,backward_degradation'

  call stress_real_replacement('replacement-small', &
                               sqrt(epsilon(1.0_wp)), max_cycles)
  call stress_real_replacement('replacement-moderate', 0.01_wp, max_cycles)
  call stress_complex_replacement('replacement-small', &
                                  sqrt(epsilon(1.0_wp)), max_cycles)
  call stress_complex_replacement('replacement-moderate', 0.01_wp, max_cycles)
  call stress_real_append_delete(max_cycles)
  call stress_complex_append_delete(max_cycles)

contains

  subroutine stress_real_replacement(scenario, update_scale, max_cycles)
  !Measure drift caused by reversible real symmetric replacements. Each cycle
  !applies one deterministic row-and-column change and its exact negative, so
  !the physical shifted matrix at every reporting checkpoint is the original
  !well-conditioned matrix.
  !
  !  Input parameters:
  !    scenario     - CSV label distinguishing the update-magnitude regime.
  !    update_scale - Scale of the deterministic replacement vector.
  !    max_cycles   - Number of positive/negative replacement pairs.
  !
  !The routine allocates state storage only during initialization. Both public
  !replacement calls in a cycle are counted by updates_since_fresh. Any public
  !operation failure terminates this standalone diagnostic executable.
    character(len=*), intent(in) :: scenario
    real(wp), intent(in) :: update_scale
    integer, intent(in) :: max_cycles
    integer, parameter :: n = 8
    type(qr_real_state) :: updated_state, fresh_state
    real(wp), parameter :: shift = 0.375_wp
    real(wp) :: h(n,n), s(n,n), shifted_matrix(n,n)
    real(wp) :: delta_h(n), delta_s(n)
    integer :: cycle, i, idx, info

    call build_real_problem(shifted_matrix, s, h, shift)
    call updated_state%initialize(n, info)
    call require_success(info, scenario, 'real updated initialization')
    call fresh_state%initialize(n, info)
    call require_success(info, scenario, 'real fresh initialization')
    call updated_state%factorize_fresh(h, s, shift, info)
    call require_success(info, scenario, 'real updated factorization')
    call fresh_state%factorize_fresh(h, s, shift, info)
    call require_success(info, scenario, 'real fresh factorization')

    call report_real(scenario, 0, update_scale, shifted_matrix, &
                     updated_state, fresh_state)
    delta_s = 0.0_wp
    do cycle = 1, max_cycles
      idx = mod(cycle - 1, n) + 1
      do i = 1, n
        delta_h(i) = update_scale * &
          real(mod(17 * cycle + 11 * i, 11) - 5, wp) / 5.0_wp
      end do

      call updated_state%replace_symmetric(idx, delta_h, delta_s, info)
      call require_success(info, scenario, 'positive real replacement')
      call updated_state%replace_symmetric(idx, -delta_h, delta_s, info)
      call require_success(info, scenario, 'negative real replacement')

      if (is_checkpoint(cycle, max_cycles)) then
        call report_real(scenario, cycle, update_scale, shifted_matrix, &
                         updated_state, fresh_state)
      end if
    end do
  end subroutine stress_real_replacement

  subroutine stress_complex_replacement(scenario, update_scale, max_cycles)
  !Measure drift caused by reversible complex Hermitian replacements. The
  !diagonal component of every update is exactly real, and applying the
  !negative vector restores the same physical Hermitian matrix after each
  !cycle.
  !
  !  Input parameters:
  !    scenario     - CSV label distinguishing the update-magnitude regime.
  !    update_scale - Scale of the deterministic Hermitian update vector.
  !    max_cycles   - Number of positive/negative replacement pairs.
  !
  !The fresh state remains untouched and supplies the working-precision
  !baseline for every checkpoint. No allocation occurs inside the update loop.
    character(len=*), intent(in) :: scenario
    real(wp), intent(in) :: update_scale
    integer, intent(in) :: max_cycles
    integer, parameter :: n = 8
    type(qr_complex_state) :: updated_state, fresh_state
    real(wp), parameter :: shift = 0.375_wp
    complex(wp) :: h(n,n), s(n,n), shifted_matrix(n,n)
    complex(wp) :: delta_h(n), delta_s(n)
    real(wp) :: real_part, imaginary_part
    integer :: cycle, i, idx, info

    call build_complex_problem(shifted_matrix, s, h, shift)
    call updated_state%initialize(n, info)
    call require_success(info, scenario, 'complex updated initialization')
    call fresh_state%initialize(n, info)
    call require_success(info, scenario, 'complex fresh initialization')
    call updated_state%factorize_fresh(h, s, shift, info)
    call require_success(info, scenario, 'complex updated factorization')
    call fresh_state%factorize_fresh(h, s, shift, info)
    call require_success(info, scenario, 'complex fresh factorization')

    call report_complex(scenario, 0, update_scale, shifted_matrix, &
                        updated_state, fresh_state)
    delta_s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do cycle = 1, max_cycles
      idx = mod(cycle - 1, n) + 1
      do i = 1, n
        real_part = update_scale * &
          real(mod(13 * cycle + 7 * i, 11) - 5, wp) / 5.0_wp
        imaginary_part = update_scale * &
          real(mod(19 * cycle + 3 * i, 9) - 4, wp) / 5.0_wp
        delta_h(i) = cmplx(real_part, imaginary_part, kind=wp)
      end do
      delta_h(idx) = cmplx(real(delta_h(idx), wp), 0.0_wp, kind=wp)

      call updated_state%replace_symmetric(idx, delta_h, delta_s, info)
      call require_success(info, scenario, 'positive complex replacement')
      call updated_state%replace_symmetric(idx, -delta_h, delta_s, info)
      call require_success(info, scenario, 'negative complex replacement')

      if (is_checkpoint(cycle, max_cycles)) then
        call report_complex(scenario, cycle, update_scale, shifted_matrix, &
                            updated_state, fresh_state)
      end if
    end do
  end subroutine stress_complex_replacement

  subroutine stress_real_append_delete(max_cycles)
  !Measure drift from repeatedly appending and deleting the same final real
  !row and column. The active order is n at every checkpoint and the represented
  !matrix is exactly the initial leading principal matrix.
  !
  !  Input parameter:
  !    max_cycles - Number of append/delete round trips.
  !
  !The appended overlap column extends the identity, while the shifted column
  !is chosen to keep the intermediate matrix strictly diagonally dominant.
    integer, intent(in) :: max_cycles
    integer, parameter :: n = 8, capacity = n + 1
    character(len=*), parameter :: operation_pair = 'append-delete'
    type(qr_real_state) :: updated_state, fresh_base_state
    type(qr_real_state) :: fresh_appended_state
    real(wp), parameter :: shift = 0.375_wp
    real(wp) :: h(n,n), s(n,n), shifted_matrix(n,n)
    real(wp) :: appended_h(capacity,capacity)
    real(wp) :: appended_s(capacity,capacity)
    real(wp) :: appended_shifted_matrix(capacity,capacity)
    real(wp) :: h_column(capacity), s_column(capacity)
    integer :: cycle, i, info

    call build_real_problem(shifted_matrix, s, h, shift)
    call updated_state%initialize(capacity, info)
    call require_success(info, operation_pair, 'real updated initialization')
    call fresh_base_state%initialize(capacity, info)
    call require_success(info, operation_pair, 'real fresh initialization')
    call fresh_appended_state%initialize(capacity, info)
    call require_success(info, operation_pair, &
                         'real appended-fresh initialization')
    call updated_state%factorize_fresh(h, s, shift, info)
    call require_success(info, operation_pair, 'real updated factorization')
    call fresh_base_state%factorize_fresh(h, s, shift, info)
    call require_success(info, operation_pair, 'real fresh factorization')

    s_column = 0.0_wp
    s_column(capacity) = 1.0_wp
    do i = 1, n
      h_column(i) = 0.02_wp * real(mod(7 * i, 9) - 4, wp)
    end do
    h_column(capacity) = 2.8_wp + shift

    appended_shifted_matrix = 0.0_wp
    appended_shifted_matrix(1:n,1:n) = shifted_matrix
    do i = 1, n
      appended_shifted_matrix(i,capacity) = h_column(i) - &
                                             shift * s_column(i)
      appended_shifted_matrix(capacity,i) = &
        appended_shifted_matrix(i,capacity)
    end do
    appended_shifted_matrix(capacity,capacity) = h_column(capacity) - &
                                                  shift * s_column(capacity)
    appended_s = 0.0_wp
    do i = 1, capacity
      appended_s(i,i) = 1.0_wp
    end do
    appended_h = appended_shifted_matrix + shift * appended_s
    call fresh_appended_state%factorize_fresh(appended_h, appended_s, &
                                               shift, info)
    call require_success(info, operation_pair, &
                         'real appended-fresh factorization')

    call report_real('delete', 0, 0.0_wp, shifted_matrix, updated_state, &
                     fresh_base_state)
    do cycle = 1, max_cycles
      call updated_state%append_symmetric(h_column, s_column, info)
      call require_success(info, operation_pair, 'real append')
      if (is_checkpoint(cycle, max_cycles)) then
        call report_real('append', cycle, 0.0_wp, &
                         appended_shifted_matrix, updated_state, &
                         fresh_appended_state)
      end if

      call updated_state%delete_symmetric(capacity, info)
      call require_success(info, operation_pair, 'real delete')

      if (is_checkpoint(cycle, max_cycles)) then
        call report_real('delete', cycle, 0.0_wp, shifted_matrix, &
                         updated_state, fresh_base_state)
      end if
    end do
  end subroutine stress_real_append_delete

  subroutine stress_complex_append_delete(max_cycles)
  !Measure drift from repeatedly appending and deleting the same final complex
  !Hermitian row and column. The operation pair restores the initial matrix and
  !order before every measurement.
  !
  !  Input parameter:
  !    max_cycles - Number of append/delete round trips.
  !
  !The physical H and S diagonal entries are exactly real, satisfying the
  !Hermitian append contract independently rather than through cancellation.
    integer, intent(in) :: max_cycles
    integer, parameter :: n = 8, capacity = n + 1
    character(len=*), parameter :: operation_pair = 'append-delete'
    type(qr_complex_state) :: updated_state, fresh_base_state
    type(qr_complex_state) :: fresh_appended_state
    real(wp), parameter :: shift = 0.375_wp
    complex(wp) :: h(n,n), s(n,n), shifted_matrix(n,n)
    complex(wp) :: appended_h(capacity,capacity)
    complex(wp) :: appended_s(capacity,capacity)
    complex(wp) :: appended_shifted_matrix(capacity,capacity)
    complex(wp) :: h_column(capacity), s_column(capacity)
    real(wp) :: real_part, imaginary_part
    integer :: cycle, i, info

    call build_complex_problem(shifted_matrix, s, h, shift)
    call updated_state%initialize(capacity, info)
    call require_success(info, operation_pair, 'complex updated initialization')
    call fresh_base_state%initialize(capacity, info)
    call require_success(info, operation_pair, 'complex fresh initialization')
    call fresh_appended_state%initialize(capacity, info)
    call require_success(info, operation_pair, &
                         'complex appended-fresh initialization')
    call updated_state%factorize_fresh(h, s, shift, info)
    call require_success(info, operation_pair, 'complex updated factorization')
    call fresh_base_state%factorize_fresh(h, s, shift, info)
    call require_success(info, operation_pair, 'complex fresh factorization')

    s_column = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s_column(capacity) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    do i = 1, n
      real_part = 0.015_wp * real(mod(7 * i, 9) - 4, wp)
      imaginary_part = 0.01_wp * real(mod(5 * i, 7) - 3, wp)
      h_column(i) = cmplx(real_part, imaginary_part, kind=wp)
    end do
    h_column(capacity) = cmplx(2.8_wp + shift, 0.0_wp, kind=wp)

    appended_shifted_matrix = cmplx(0.0_wp, 0.0_wp, kind=wp)
    appended_shifted_matrix(1:n,1:n) = shifted_matrix
    do i = 1, n
      appended_shifted_matrix(i,capacity) = h_column(i) - &
        cmplx(shift, 0.0_wp, kind=wp) * s_column(i)
      appended_shifted_matrix(capacity,i) = &
        conjg(appended_shifted_matrix(i,capacity))
    end do
    appended_shifted_matrix(capacity,capacity) = h_column(capacity) - &
      cmplx(shift, 0.0_wp, kind=wp) * s_column(capacity)
    appended_s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do i = 1, capacity
      appended_s(i,i) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do
    appended_h = appended_shifted_matrix + &
                 cmplx(shift, 0.0_wp, kind=wp) * appended_s
    call fresh_appended_state%factorize_fresh(appended_h, appended_s, &
                                               shift, info)
    call require_success(info, operation_pair, &
                         'complex appended-fresh factorization')

    call report_complex('delete', 0, 0.0_wp, shifted_matrix, &
                        updated_state, fresh_base_state)
    do cycle = 1, max_cycles
      call updated_state%append_symmetric(h_column, s_column, info)
      call require_success(info, operation_pair, 'complex append')
      if (is_checkpoint(cycle, max_cycles)) then
        call report_complex('append', cycle, 0.0_wp, &
                            appended_shifted_matrix, updated_state, &
                            fresh_appended_state)
      end if

      call updated_state%delete_symmetric(capacity, info)
      call require_success(info, operation_pair, 'complex delete')

      if (is_checkpoint(cycle, max_cycles)) then
        call report_complex('delete', cycle, 0.0_wp, shifted_matrix, &
                            updated_state, fresh_base_state)
      end if
    end do
  end subroutine stress_complex_append_delete

  subroutine build_real_problem(shifted_matrix, s, h, shift)
  !Construct a deterministic, dense, strictly diagonally dominant real shifted
  !matrix with S=I and H=shifted_matrix+shift*S. Strict diagonal dominance
  !keeps the primary experiment away from conditioning-driven solve error.
  !
  !  Input parameter:
  !    shift - Fixed shift represented by the factorization.
  !
  !  Output parameters:
  !    shifted_matrix - Full symmetric matrix used for external measurements.
  !    s              - Identity overlap matrix.
  !    h              - Physical symmetric H matrix.
    real(wp), intent(out) :: shifted_matrix(:,:), s(:,:), h(:,:)
    real(wp), intent(in) :: shift
    integer :: i, j, n

    n = size(shifted_matrix, 1)
    shifted_matrix = 0.0_wp
    s = 0.0_wp
    do j = 1, n
      shifted_matrix(j,j) = 2.0_wp + 0.05_wp * real(j, wp)
      s(j,j) = 1.0_wp
      do i = j + 1, n
        shifted_matrix(i,j) = 0.015_wp * &
          real(mod(13 * i + 7 * j, 11) - 5, wp)
        shifted_matrix(j,i) = shifted_matrix(i,j)
      end do
    end do
    h = shifted_matrix + shift * s
  end subroutine build_real_problem

  subroutine build_complex_problem(shifted_matrix, s, h, shift)
  !Construct the Hermitian counterpart of build_real_problem. The diagonal is
  !real and dominant, and deterministic complex off-diagonal entries exercise
  !the conjugate-transpose path without introducing an ill-conditioned matrix.
  !
  !  Input parameter:
  !    shift - Fixed real shift represented by the factorization.
  !
  !  Output parameters:
  !    shifted_matrix - Full Hermitian matrix for external measurements.
  !    s              - Complex identity overlap matrix.
  !    h              - Physical Hermitian H matrix.
    complex(wp), intent(out) :: shifted_matrix(:,:), s(:,:), h(:,:)
    real(wp), intent(in) :: shift
    real(wp) :: real_part, imaginary_part
    integer :: i, j, n

    n = size(shifted_matrix, 1)
    shifted_matrix = cmplx(0.0_wp, 0.0_wp, kind=wp)
    s = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do j = 1, n
      shifted_matrix(j,j) = cmplx(2.0_wp + 0.05_wp * real(j, wp), &
                                  0.0_wp, kind=wp)
      s(j,j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
      do i = j + 1, n
        real_part = 0.012_wp * real(mod(11 * i + 5 * j, 11) - 5, wp)
        imaginary_part = 0.008_wp * &
          real(mod(7 * i + 13 * j, 9) - 4, wp)
        shifted_matrix(i,j) = cmplx(real_part, imaginary_part, kind=wp)
        shifted_matrix(j,i) = conjg(shifted_matrix(i,j))
      end do
    end do
    h = shifted_matrix + cmplx(shift, 0.0_wp, kind=wp) * s
  end subroutine build_complex_problem

  subroutine report_real(scenario, cycles, update_scale, shifted_matrix, &
                         updated_state, fresh_state)
  !Measure and print one real update-drift checkpoint. The untouched fresh
  !state supplies a numerical floor at the same working precision and for the
  !same matrix, while updates_since_fresh reports the number of public
  !structural operations rather than positive/negative operation pairs.
    character(len=*), intent(in) :: scenario
    integer, intent(in) :: cycles
    real(wp), intent(in) :: update_scale, shifted_matrix(:,:)
    type(qr_real_state), intent(in) :: updated_state, fresh_state
    real(wp) :: factor_residual, orthogonality, forward_operator
    real(wp) :: residual_operator, max_backward
    real(wp) :: fresh_factor_residual, fresh_orthogonality
    real(wp) :: fresh_forward_operator, fresh_residual_operator
    real(wp) :: fresh_max_backward, backward_degradation, error_floor

    call measure_real(updated_state, shifted_matrix, factor_residual, &
                      orthogonality, forward_operator, residual_operator, &
                      max_backward)
    call measure_real(fresh_state, shifted_matrix, fresh_factor_residual, &
                      fresh_orthogonality, fresh_forward_operator, &
                      fresh_residual_operator, fresh_max_backward)
    error_floor = real(size(shifted_matrix, 1), wp) * epsilon(1.0_wp)
    backward_degradation = max(max_backward, error_floor) / &
                           max(fresh_max_backward, error_floor)

    write(*,'(a,",",a,",",i0,",",i0,",",i0,12(",",es16.8))') &
      trim(scenario), 'real', wp, cycles, updated_state%updates_since_fresh, &
      update_scale, factor_residual, orthogonality, forward_operator, &
      residual_operator, max_backward, fresh_factor_residual, &
      fresh_orthogonality, fresh_forward_operator, fresh_residual_operator, &
      fresh_max_backward, backward_degradation
  end subroutine report_real

  subroutine report_complex(scenario, cycles, update_scale, shifted_matrix, &
                            updated_state, fresh_state)
  !Measure and print one complex update-drift checkpoint using the same metric
  !definitions and fresh-factor baseline as the real path.
    character(len=*), intent(in) :: scenario
    integer, intent(in) :: cycles
    real(wp), intent(in) :: update_scale
    complex(wp), intent(in) :: shifted_matrix(:,:)
    type(qr_complex_state), intent(in) :: updated_state, fresh_state
    real(wp) :: factor_residual, unitarity, forward_operator
    real(wp) :: residual_operator, max_backward
    real(wp) :: fresh_factor_residual, fresh_unitarity
    real(wp) :: fresh_forward_operator, fresh_residual_operator
    real(wp) :: fresh_max_backward, backward_degradation, error_floor

    call measure_complex(updated_state, shifted_matrix, factor_residual, &
                         unitarity, forward_operator, residual_operator, &
                         max_backward)
    call measure_complex(fresh_state, shifted_matrix, fresh_factor_residual, &
                         fresh_unitarity, fresh_forward_operator, &
                         fresh_residual_operator, fresh_max_backward)
    error_floor = real(size(shifted_matrix, 1), wp) * epsilon(1.0_wp)
    backward_degradation = max(max_backward, error_floor) / &
                           max(fresh_max_backward, error_floor)

    write(*,'(a,",",a,",",i0,",",i0,",",i0,12(",",es16.8))') &
      trim(scenario), 'complex', wp, cycles, &
      updated_state%updates_since_fresh, update_scale, factor_residual, &
      unitarity, forward_operator, residual_operator, max_backward, &
      fresh_factor_residual, fresh_unitarity, fresh_forward_operator, &
      fresh_residual_operator, fresh_max_backward, backward_degradation
  end subroutine report_complex

  subroutine measure_real(state, shifted_matrix, factor_residual, &
                          orthogonality, forward_operator, residual_operator, &
                          max_backward)
  !Evaluate factor quality and the complete direct-solve operator for one real
  !state. If B=R^{-1}Q^T, then B*M-I measures forward error for exact systems
  !b=M*x, while M*B-I measures residual propagation over all right-hand sides.
  !The maximum normwise backward error is additionally evaluated for the n
  !known solutions given by the coordinate vectors.
  !
  !  Input parameters:
  !    state          - Valid state whose active order equals the matrix order.
  !    shifted_matrix - Authoritative full symmetric M matrix.
  !
  !  Output parameters:
  !    factor_residual  - ||M-Q*R||_F/||M||_F.
  !    orthogonality    - ||I-Q^T*Q||_F/sqrt(n).
  !    forward_operator - ||B*M-I||_F/sqrt(n).
  !    residual_operator- ||M*B-I||_F/sqrt(n).
  !    max_backward     - Largest coordinate-system backward error.
    type(qr_real_state), intent(in) :: state
    real(wp), intent(in) :: shifted_matrix(:,:)
    real(wp), intent(out) :: factor_residual, orthogonality
    real(wp), intent(out) :: forward_operator, residual_operator
    real(wp), intent(out) :: max_backward
    integer :: i, j, n
    logical :: usable
    real(wp) :: inverse_action(size(shifted_matrix,1), &
                               size(shifted_matrix,2))
    real(wp) :: forward_map(size(shifted_matrix,1), &
                            size(shifted_matrix,2))
    real(wp) :: residual_map(size(shifted_matrix,1), &
                             size(shifted_matrix,2))
    real(wp) :: identity(size(shifted_matrix,1),size(shifted_matrix,2))
    real(wp) :: rhs(size(shifted_matrix,1)), solution(size(shifted_matrix,1))
    real(wp) :: residual(size(shifted_matrix,1)), denominator, matrix_norm

    n = size(shifted_matrix, 1)
    identity = 0.0_wp
    do i = 1, n
      identity(i,i) = 1.0_wp
    end do

    do j = 1, n
      rhs = 0.0_wp
      rhs(j) = 1.0_wp
      call apply_real_factors(state, rhs, inverse_action(:,j), usable)
      if (.not. usable) then
        factor_residual = huge(1.0_wp)
        orthogonality = huge(1.0_wp)
        forward_operator = huge(1.0_wp)
        residual_operator = huge(1.0_wp)
        max_backward = huge(1.0_wp)
        return
      end if
    end do

    factor_residual = real_frobenius(shifted_matrix - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(real_frobenius(shifted_matrix), tiny(1.0_wp))
    orthogonality = real_frobenius(identity - &
      matmul(transpose(state%q(1:n,1:n)), state%q(1:n,1:n))) / &
      sqrt(real(n, wp))

    forward_map = matmul(inverse_action, shifted_matrix)
    residual_map = matmul(shifted_matrix, inverse_action)
    forward_operator = real_frobenius(forward_map - identity) / &
                       sqrt(real(n, wp))
    residual_operator = real_frobenius(residual_map - identity) / &
                        sqrt(real(n, wp))

    matrix_norm = real_matrix_infinity_norm(shifted_matrix)
    max_backward = 0.0_wp
    do j = 1, n
      rhs = shifted_matrix(:,j)
      solution = forward_map(:,j)
      residual = rhs - matmul(shifted_matrix, solution)
      denominator = matrix_norm * maxval(abs(solution)) + maxval(abs(rhs))
      max_backward = max(max_backward, maxval(abs(residual)) / &
                         max(denominator, tiny(1.0_wp)))
    end do
  end subroutine measure_real

  subroutine measure_complex(state, shifted_matrix, factor_residual, &
                             unitarity, forward_operator, residual_operator, &
                             max_backward)
  !Evaluate the complex Hermitian factor and solve-operator metrics. The
  !effective inverse is B=R^{-1}Q^H, and every norm is accumulated in real(wp)
  !so results remain comparable across all supported working kinds.
    type(qr_complex_state), intent(in) :: state
    complex(wp), intent(in) :: shifted_matrix(:,:)
    real(wp), intent(out) :: factor_residual, unitarity
    real(wp), intent(out) :: forward_operator, residual_operator
    real(wp), intent(out) :: max_backward
    integer :: i, j, n
    logical :: usable
    complex(wp) :: inverse_action(size(shifted_matrix,1), &
                                  size(shifted_matrix,2))
    complex(wp) :: forward_map(size(shifted_matrix,1), &
                               size(shifted_matrix,2))
    complex(wp) :: residual_map(size(shifted_matrix,1), &
                                size(shifted_matrix,2))
    complex(wp) :: identity(size(shifted_matrix,1),size(shifted_matrix,2))
    complex(wp) :: rhs(size(shifted_matrix,1))
    complex(wp) :: solution(size(shifted_matrix,1))
    complex(wp) :: residual(size(shifted_matrix,1))
    real(wp) :: denominator, matrix_norm

    n = size(shifted_matrix, 1)
    identity = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do i = 1, n
      identity(i,i) = cmplx(1.0_wp, 0.0_wp, kind=wp)
    end do

    do j = 1, n
      rhs = cmplx(0.0_wp, 0.0_wp, kind=wp)
      rhs(j) = cmplx(1.0_wp, 0.0_wp, kind=wp)
      call apply_complex_factors(state, rhs, inverse_action(:,j), usable)
      if (.not. usable) then
        factor_residual = huge(1.0_wp)
        unitarity = huge(1.0_wp)
        forward_operator = huge(1.0_wp)
        residual_operator = huge(1.0_wp)
        max_backward = huge(1.0_wp)
        return
      end if
    end do

    factor_residual = complex_frobenius(shifted_matrix - &
      matmul(state%q(1:n,1:n), state%r(1:n,1:n))) / &
      max(complex_frobenius(shifted_matrix), tiny(1.0_wp))
    unitarity = complex_frobenius(identity - &
      matmul(conjg(transpose(state%q(1:n,1:n))), &
             state%q(1:n,1:n))) / sqrt(real(n, wp))

    forward_map = matmul(inverse_action, shifted_matrix)
    residual_map = matmul(shifted_matrix, inverse_action)
    forward_operator = complex_frobenius(forward_map - identity) / &
                       sqrt(real(n, wp))
    residual_operator = complex_frobenius(residual_map - identity) / &
                        sqrt(real(n, wp))

    matrix_norm = complex_matrix_infinity_norm(shifted_matrix)
    max_backward = 0.0_wp
    do j = 1, n
      rhs = shifted_matrix(:,j)
      solution = forward_map(:,j)
      residual = rhs - matmul(shifted_matrix, solution)
      denominator = matrix_norm * maxval(abs(solution)) + maxval(abs(rhs))
      max_backward = max(max_backward, maxval(abs(residual)) / &
                         max(denominator, tiny(1.0_wp)))
    end do
  end subroutine measure_complex

  subroutine apply_real_factors(state, rhs, solution, usable)
  !Apply B=R^{-1}Q^T to one real right-hand side using explicit loops. This
  !test-only implementation isolates the stored-factor solve from inverse
  !iteration, normalization, and stopping behavior.
  !
  !  Input parameters:
  !    state - Valid real QR state.
  !    rhs   - Right-hand side of active length n.
  !
  !  Output parameters:
  !    solution - Computed direct solution when usable is true.
  !    usable   - False when an active R diagonal cannot be divided safely.
    type(qr_real_state), intent(in) :: state
    real(wp), intent(in) :: rhs(:)
    real(wp), intent(out) :: solution(:)
    logical, intent(out) :: usable
    real(wp) :: transformed(size(rhs)), diagonal_threshold, factor_scale
    integer :: i, j, n

    n = size(rhs)
    transformed = 0.0_wp
    do i = 1, n
      do j = 1, n
        transformed(i) = transformed(i) + state%q(j,i) * rhs(j)
      end do
    end do

    factor_scale = maxval(abs(state%r(1:n,1:n)))
    diagonal_threshold = max(tiny(1.0_wp), &
                             epsilon(1.0_wp) * factor_scale)
    solution = transformed
    usable = .true.
    do i = n, 1, -1
      if (abs(state%r(i,i)) <= diagonal_threshold) then
        usable = .false.
        solution = 0.0_wp
        return
      end if
      do j = i + 1, n
        solution(i) = solution(i) - state%r(i,j) * solution(j)
      end do
      solution(i) = solution(i) / state%r(i,i)
    end do
  end subroutine apply_real_factors

  subroutine apply_complex_factors(state, rhs, solution, usable)
  !Apply B=R^{-1}Q^H to one complex right-hand side. The conjugation of Q is
  !essential because the state represents a unitary rather than orthogonal
  !factorization.
    type(qr_complex_state), intent(in) :: state
    complex(wp), intent(in) :: rhs(:)
    complex(wp), intent(out) :: solution(:)
    logical, intent(out) :: usable
    complex(wp) :: transformed(size(rhs))
    real(wp) :: diagonal_threshold, factor_scale
    integer :: i, j, n

    n = size(rhs)
    transformed = cmplx(0.0_wp, 0.0_wp, kind=wp)
    do i = 1, n
      do j = 1, n
        transformed(i) = transformed(i) + conjg(state%q(j,i)) * rhs(j)
      end do
    end do

    factor_scale = maxval(abs(state%r(1:n,1:n)))
    diagonal_threshold = max(tiny(1.0_wp), &
                             epsilon(1.0_wp) * factor_scale)
    solution = transformed
    usable = .true.
    do i = n, 1, -1
      if (abs(state%r(i,i)) <= diagonal_threshold) then
        usable = .false.
        solution = cmplx(0.0_wp, 0.0_wp, kind=wp)
        return
      end if
      do j = i + 1, n
        solution(i) = solution(i) - state%r(i,j) * solution(j)
      end do
      solution(i) = solution(i) / state%r(i,i)
    end do
  end subroutine apply_complex_factors

  function real_frobenius(matrix) result(norm)
  !Return the Frobenius norm of a real matrix in the selected working kind.
    real(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm

    norm = sqrt(sum(matrix * matrix))
  end function real_frobenius

  function complex_frobenius(matrix) result(norm)
  !Return the Frobenius norm of a complex matrix in the selected working kind.
    complex(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm

    norm = sqrt(sum(abs(matrix)**2))
  end function complex_frobenius

  function real_matrix_infinity_norm(matrix) result(norm)
  !Return max_i sum_j abs(matrix(i,j)) for a real matrix.
    real(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm
    integer :: i

    norm = 0.0_wp
    do i = 1, size(matrix, 1)
      norm = max(norm, sum(abs(matrix(i,:))))
    end do
  end function real_matrix_infinity_norm

  function complex_matrix_infinity_norm(matrix) result(norm)
  !Return max_i sum_j abs(matrix(i,j)) for a complex matrix.
    complex(wp), intent(in) :: matrix(:,:)
    real(wp) :: norm
    integer :: i

    norm = 0.0_wp
    do i = 1, size(matrix, 1)
      norm = max(norm, sum(abs(matrix(i,:))))
    end do
  end function complex_matrix_infinity_norm

  function is_checkpoint(cycle, max_cycles) result(report_now)
  !Select powers of two and the final cycle for logarithmically spaced output.
  !This keeps long runs readable while retaining the onset and growth of drift.
    integer, intent(in) :: cycle, max_cycles
    logical :: report_now
    integer :: reduced_cycle

    reduced_cycle = cycle
    do while (reduced_cycle > 1 .and. mod(reduced_cycle, 2) == 0)
      reduced_cycle = reduced_cycle / 2
    end do
    report_now = reduced_cycle == 1 .or. cycle == max_cycles
  end function is_checkpoint

  subroutine require_success(info, scenario, operation)
  !Terminate the standalone stress diagnostic when a structural operation
  !fails, because subsequent numerical rows would no longer describe the
  !requested reversible sequence.
    integer, intent(in) :: info
    character(len=*), intent(in) :: scenario, operation

    if (info /= QR_SUCCESS) then
      write(*,'(3(a,1x),i0)') 'stress_update_drift:', trim(scenario), &
        trim(operation), info
      error stop 1
    end if
  end subroutine require_success

end program stress_update_drift

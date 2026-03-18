! Test program for gjr_module.
!
! 1. Numerical gradient check: compares analytical gradient from gjr_obj
!    with central finite differences at the starting point.
! 2. Parameter recovery: simulates GJR-GARCH data with known parameters,
!    fits by BFGS, and reports estimated vs true values.

program xgjr
    use kind_mod,    only: dp
    use gjr_module,  only: gjr_simulate, gjr_set_data, gjr_obj, &
                            gjr_transform, gjr_inv_transform
    use bfgs_module, only: bfgs_minimize
    implicit none

    integer,  parameter :: nobs      = 2000
    integer,  parameter :: seed_val  = 42
    integer,  parameter :: max_iter  = 500
    real(dp), parameter :: gtol      = 1.0e-7_dp

    ! true parameters for simulation
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.05_dp
    real(dp), parameter :: true_gamma = 0.10_dp  ! asymmetry / leverage
    real(dp), parameter :: true_beta  = 0.85_dp

    real(dp) :: y(nobs)
    real(dp) :: p(4), p_pert(4), f, f_plus, f_minus, g(4), g_fd(4)
    real(dp) :: omega, alpha, gamma, beta, uncond_var, rel_err
    real(dp), parameter :: h_fd = 1.0e-5_dp  ! finite-difference step
    integer  :: i, iter
    logical  :: converged

    ! ---- simulate ----
    call gjr_simulate(true_omega, true_alpha, true_gamma, true_beta, &
                      nobs, seed_val, y)
    call gjr_set_data(y, nobs)

    ! ---- gradient check at a point away from the optimum ----
    call gjr_inv_transform(0.05_dp, 0.07_dp, 0.08_dp, 0.82_dp, p)
    call gjr_obj(p, 4, f, g)

    do i = 1, 4
        p_pert    = p
        p_pert(i) = p(i) + h_fd
        call gjr_obj(p_pert, 4, f_plus, g_fd)
        p_pert    = p
        p_pert(i) = p(i) - h_fd
        call gjr_obj(p_pert, 4, f_minus, g_fd)
        g_fd(i) = (f_plus - f_minus) / (2.0_dp * h_fd)
    end do

    print '(A)', ""
    print '(A)', " Gradient check (central finite differences, h=1e-5):"
    print '(A)', "   i   analytical     fin-diff       rel.error"
    do i = 1, 4
        rel_err = abs(g(i) - g_fd(i)) / (abs(g_fd(i)) + 1.0e-15_dp)
        print '(I4,2ES15.6,ES13.3)', i, g(i), g_fd(i), rel_err
    end do

    ! ---- fit ----
    call gjr_inv_transform(0.05_dp, 0.07_dp, 0.08_dp, 0.82_dp, p)
    call bfgs_minimize(gjr_obj, p, 4, max_iter, gtol, f, iter, converged)
    call gjr_transform(p, omega, alpha, gamma, beta)
    uncond_var = omega / (1.0_dp - alpha - 0.5_dp*gamma - beta)

    print '(A)', ""
    print '(A)', " Parameter recovery:"
    print '(A,I0)', "   nobs = ", nobs
    print '(A)', ""
    print '(A)', "              omega      alpha      gamma       beta  uncond_vol"
    print '(A,5F11.4)', "   True      ", true_omega, true_alpha, true_gamma, true_beta, &
        sqrt(true_omega / (1.0_dp - true_alpha - 0.5_dp*true_gamma - true_beta))
    print '(A,5F11.4)', "   Estimated ", omega, alpha, gamma, beta, sqrt(uncond_var)
    print '(A)', ""
    print '(A,I0,A,L1)', "   Iterations: ", iter, "   Converged: ", converged
    print '(A,ES12.4)',  "   NLL/n:      ", f

end program xgjr

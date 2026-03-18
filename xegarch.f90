! Test program for egarch_module.
!
! 1. Numerical gradient check: compares analytical gradient from egarch_obj
!    with central finite differences at the starting point.
! 2. Parameter recovery: simulates EGARCH data with known parameters,
!    fits by BFGS, and reports estimated vs true values.

program xegarch
    use kind_mod,      only: dp
    use egarch_module, only: egarch_simulate, egarch_set_data, egarch_obj, &
                              egarch_transform, egarch_inv_transform
    use bfgs_module,   only: bfgs_minimize
    implicit none

    integer,  parameter :: nobs      = 2000
    integer,  parameter :: seed_val  = 42
    integer,  parameter :: max_iter  = 500
    real(dp), parameter :: gtol      = 1.0e-7_dp

    ! true parameters for simulation
    ! log(h_1) = omega/(1-beta) = -0.03/0.03 = -1.0  =>  h_1 = exp(-1) ~ 0.37
    real(dp), parameter :: true_omega = -0.03_dp
    real(dp), parameter :: true_alpha =  0.10_dp
    real(dp), parameter :: true_gamma = -0.08_dp  ! negative -> leverage effect
    real(dp), parameter :: true_beta  =  0.97_dp

    real(dp) :: y(nobs)
    real(dp) :: p(4), p_pert(4), f, f_plus, f_minus, g(4), g_fd(4)
    real(dp) :: omega, alpha, gamma, beta, rel_err
    real(dp), parameter :: h_fd = 1.0e-5_dp
    integer  :: i, iter
    logical  :: converged

    ! ---- simulate ----
    call egarch_simulate(true_omega, true_alpha, true_gamma, true_beta, &
                         nobs, seed_val, y)
    call egarch_set_data(y, nobs)

    ! ---- gradient check at a point away from the optimum ----
    call egarch_inv_transform(-0.05_dp, 0.08_dp, -0.05_dp, 0.90_dp, p)
    call egarch_obj(p, 4, f, g)

    do i = 1, 4
        p_pert    = p
        p_pert(i) = p(i) + h_fd
        call egarch_obj(p_pert, 4, f_plus, g_fd)
        p_pert    = p
        p_pert(i) = p(i) - h_fd
        call egarch_obj(p_pert, 4, f_minus, g_fd)
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
    call egarch_inv_transform(-0.05_dp, 0.08_dp, -0.05_dp, 0.90_dp, p)
    call bfgs_minimize(egarch_obj, p, 4, max_iter, gtol, f, iter, converged)
    call egarch_transform(p, omega, alpha, gamma, beta)

    print '(A)', ""
    print '(A)', " Parameter recovery:"
    print '(A,I0)', "   nobs = ", nobs
    print '(A)', ""
    print '(A)', "              omega      alpha      gamma       beta   exp_lh_unc"
    print '(A,5F11.4)', "   True      ", true_omega, true_alpha, true_gamma, true_beta, &
        exp(true_omega / (1.0_dp - true_beta))
    print '(A,5F11.4)', "   Estimated ", omega, alpha, gamma, beta, &
        exp(omega / (1.0_dp - beta))
    print '(A)', ""
    print '(A,I0,A,L1)', "   Iterations: ", iter, "   Converged: ", converged
    print '(A,ES12.4)',  "   NLL/n:      ", f

end program xegarch

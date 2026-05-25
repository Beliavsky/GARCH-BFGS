! Test program for nagarch_mod.
!
! 1. Numerical gradient check: compares analytical gradient from nagarch_obj
!    with central finite differences at the starting point.
! 2. Parameter recovery: simulates NAGARCH data with known parameters,
!    fits by BFGS, and reports estimated vs true values.

program xnagarch
    use kind_mod,       only: dp
    use nagarch_mod, only: nagarch_simulate, nagarch_set_data, nagarch_obj, &
                               nagarch_transform, nagarch_inv_transform
    use bfgs_mod,    only: bfgs_minimize
    implicit none

    integer,  parameter :: nobs      = 2000
    integer,  parameter :: seed_val  = 42
    integer,  parameter :: max_iter  = 500
    real(dp), parameter :: gtol      = 1.0e-7_dp

    ! true parameters for simulation
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.05_dp
    real(dp), parameter :: true_beta  = 0.90_dp
    real(dp), parameter :: true_theta = 0.50_dp  ! positive -> leverage effect

    real(dp) :: y(nobs)
    real(dp) :: p(4), p_pert(4), f, f_plus, f_minus, g(4), g_fd(4)
    real(dp) :: omega, alpha, beta, theta, uncond_vol, rel_err
    real(dp), parameter :: h_fd = 1.0e-5_dp  ! finite-difference step
    integer  :: i, iter
    logical  :: converged

    ! ---- simulate ----
    call nagarch_simulate(true_omega, true_alpha, true_beta, true_theta, &
                           nobs, seed_val, y)
    call nagarch_set_data(y, nobs)

    ! ---- gradient check at a point away from the optimum ----
    call nagarch_inv_transform(0.05_dp, 0.08_dp, 0.85_dp, 0.30_dp, p)
    call nagarch_obj(p, 4, f, g)

    do i = 1, 4
        p_pert    = p
        p_pert(i) = p(i) + h_fd
        call nagarch_obj(p_pert, 4, f_plus, g_fd)
        p_pert    = p
        p_pert(i) = p(i) - h_fd
        call nagarch_obj(p_pert, 4, f_minus, g_fd)
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
    call nagarch_inv_transform(0.05_dp, 0.08_dp, 0.85_dp, 0.30_dp, p)
    call bfgs_minimize(nagarch_obj, p, 4, max_iter, gtol, f, iter, converged)
    call nagarch_transform(p, omega, alpha, beta, theta)
    uncond_vol = sqrt(omega / (1.0_dp - alpha*(1.0_dp + theta**2) - beta))

    print '(A)', ""
    print '(A)', " Parameter recovery:"
    print '(A,I0)', "   nobs = ", nobs
    print '(A)', ""
    print '(A)', "              omega      alpha       beta      theta  uncond_vol"
    print '(A,5F11.4)', "   True      ", true_omega, true_alpha, true_beta, true_theta, &
        sqrt(true_omega / (1.0_dp - true_alpha*(1.0_dp + true_theta**2) - true_beta))
    print '(A,5F11.4)', "   Estimated ", omega, alpha, beta, theta, uncond_vol
    print '(A)', ""
    print '(A,I0,A,L1)', "   Iterations: ", iter, "   Converged: ", converged
    print '(A,ES12.4)',  "   NLL/n:      ", f

end program xnagarch

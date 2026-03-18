! Simulate a GARCH(1,1) process and estimate parameters by maximum likelihood
! using BFGS with analytical gradients.

program garch_main
    use kind_mod,    only: dp
    use garch_module, only: garch_simulate, garch_set_data, garch_inv_transform, &
                            garch_obj, garch_transform
    use bfgs_module,  only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_obs     = 2000
    integer,  parameter :: seed_val  = 42
    integer,  parameter :: max_iter  = 500
    real(dp), parameter :: gtol      = 1.0e-7_dp

    ! ---- true GARCH(1,1) parameters ----
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.10_dp
    real(dp), parameter :: true_beta  = 0.85_dp

    real(dp) :: y(n_obs)
    real(dp) :: p(3), omega, alpha, beta
    real(dp) :: f_opt, f_check, g_final(3)
    logical  :: converged
    integer  :: n_iter

    real :: t0, t1, t_sim, t_opt, t_diag, t_total

    ! ---- simulate ----
    call cpu_time(t0)
    call garch_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
    call garch_set_data(y, n_obs)
    call cpu_time(t1)
    t_sim = t1 - t0

    ! ---- starting values (slightly different from truth) ----
    call garch_inv_transform(0.05_dp, 0.10_dp, 0.80_dp, p)

    ! ---- optimise ----
    call cpu_time(t0)
    call bfgs_minimize(garch_obj, p, 3, max_iter, gtol, f_opt, n_iter, converged)
    call cpu_time(t1)
    t_opt = t1 - t0

    call garch_transform(p, omega, alpha, beta)

    ! ---- recompute at solution for gradient diagnostics ----
    call cpu_time(t0)
    call garch_obj(p, 3, f_check, g_final)
    call cpu_time(t1)
    t_diag = t1 - t0

    t_total = t_sim + t_opt + t_diag

    ! ---- results ----
    print '(A)', ""
    print '(A)', " GARCH(1,1) — MLE via BFGS with analytical gradients"
    print '(A,I0)', " Observations  : ", n_obs
    print '(A)', ""
    print '(A)', "  Parameter      True    Estimated"
    print '(A)', "  ---------   -------   ---------"
    print '(A,F9.4,F12.4)', "  omega    ", true_omega, omega
    print '(A,F9.4,F12.4)', "  alpha    ", true_alpha, alpha
    print '(A,F9.4,F12.4)', "  beta     ", true_beta,  beta
    print '(A,F9.4,F12.4)', "  alpha+b  ", true_alpha + true_beta, alpha + beta
    print '(A)', ""
    print '(A,F12.2)',  " Log-likelihood : ", -f_check
    print '(A,I0)',     " Iterations     : ", n_iter
    if (converged) then
        print '(A)', " Converged      : yes"
    else
        print '(A)', " Converged      : no (hit iteration limit)"
    end if
    print '(A,ES12.4)', " Gradient norm  : ", norm2(g_final)
    print '(A)', ""

    ! ---- timing summary ----
    print '(A)', "  Section       seconds     fraction"
    print '(A)', "  -------       -------     --------"
    print '(A,F9.4,F12.4)', "  simulate  ", t_sim,  t_sim  / t_total
    print '(A,F9.4,F12.4)', "  optimise  ", t_opt,  t_opt  / t_total
    print '(A,F9.4,F12.4)', "  diag      ", t_diag, t_diag / t_total
    print '(A,F9.4)',        "  total     ", t_total
    print '(A)', ""

end program garch_main

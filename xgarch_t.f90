! Simulate GARCH(1,1) data for each entry in dof_list, then fit two models:
!   (1) GARCH-Normal  (3 parameters: omega, alpha, beta)
!   (2) GARCH-t       (4 parameters: omega, alpha, beta, nu)
! dof_list entry 0.0 means simulate with Normal noise; otherwise Student-t(nu).

program xgarch_t
    use kind_mod,       only: dp
    use garch_module,   only: garch_simulate, garch_set_data, garch_inv_transform, &
                               garch_obj, garch_transform
    use garch_t_module, only: garch_t_simulate, garch_t_set_data, &
                               garch_t_inv_transform, garch_t_obj, garch_t_transform
    use bfgs_module,    only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_obs    = 20000
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true GARCH(1,1) parameters (same for all simulations) ----
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.10_dp
    real(dp), parameter :: true_beta  = 0.85_dp

    ! ---- degrees of freedom to loop over; 0.0 => Normal noise ----
    real(dp), parameter :: dof_list(*) = [ 0.0_dp, 5.0_dp, 8.0_dp, 15.0_dp ]

    real(dp) :: y(n_obs)
    real(dp) :: p_n(3), p_t(4)
    real(dp) :: omega_n, alpha_n, beta_n
    real(dp) :: omega_t, alpha_t, beta_t, nu_t
    real(dp) :: f_n, f_t, g_n(3), g_t(4)
    logical  :: conv_n, conv_t
    integer  :: iter_n, iter_t
    integer  :: ck0, ck1, ck_rate
    real(dp) :: t_normal, t_student
    real(dp) :: y_std, uncvol_true, uncvol_n, uncvol_t
    real(dp)          :: sim_nu
    integer           :: i_dof
    character(len=11) :: nu_t_str

    call system_clock(count_rate=ck_rate)

    uncvol_true = sqrt(true_omega / (1.0_dp - true_alpha - true_beta))

    print '(A)', ""
    print '(A)', " GARCH(1,1): Normal vs Student-t innovations"
    print '(A,I0)', " Observations : ", n_obs

    do i_dof = 1, size(dof_list)
        sim_nu = dof_list(i_dof)

        ! ---- simulate ----
        if (sim_nu == 0.0_dp) then
            call garch_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
        else
            call garch_t_simulate(true_omega, true_alpha, true_beta, sim_nu, &
                                  n_obs, seed_val, y)
        end if
        call garch_set_data(y, n_obs)
        call garch_t_set_data(y, n_obs)

        ! ---- fit GARCH-Normal ----
        call garch_inv_transform(0.05_dp, 0.10_dp, 0.80_dp, p_n)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_obj, p_n, 3, max_iter, gtol, f_n, iter_n, conv_n)
        call system_clock(count=ck1)
        t_normal = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_transform(p_n, omega_n, alpha_n, beta_n)
        call garch_obj(p_n, 3, f_n, g_n)

        ! ---- fit GARCH-t, warm-started from normal estimates ----
        call garch_t_inv_transform(omega_n, alpha_n, beta_n, 8.0_dp, p_t)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_t_obj, p_t, 4, max_iter, gtol, f_t, iter_t, conv_t)
        call system_clock(count=ck1)
        t_student = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_t_transform(p_t, omega_t, alpha_t, beta_t, nu_t)
        call garch_t_obj(p_t, 4, f_t, g_t)

        ! ---- derived quantities ----
        y_std    = sqrt(sum((y - sum(y)/n_obs)**2) / (n_obs - 1))
        uncvol_n = sqrt(omega_n / (1.0_dp - alpha_n - beta_n))
        uncvol_t = sqrt(omega_t / (1.0_dp - alpha_t - beta_t))

        ! ---- results ----
        print '(A)', ""
        if (sim_nu == 0.0_dp) then
            print '(A)', " Simulated noise: Normal"
        else
            print '(A,F6.1)', " Simulated noise: Student-t, nu =", sim_nu
        end if
        print '(A,F10.4)', " Simulated data std dev: ", y_std
        print '(A)', ""

        ! parameter table
        print '(A)', "               True      Normal   Student-t"
        print '(A)', "  ---------  --------   --------   --------"
        print '(A,F10.4,F11.4,F11.4)', "  omega    ", true_omega, omega_n, omega_t
        print '(A,F10.4,F11.4,F11.4)', "  alpha    ", true_alpha, alpha_n, alpha_t
        print '(A,F10.4,F11.4,F11.4)', "  beta     ", true_beta,  beta_n,  beta_t
        if (nu_t > 99.0_dp) then
            nu_t_str = "      large"
        else
            write(nu_t_str, '(F11.4)') nu_t
        end if
        if (sim_nu == 0.0_dp) then
            print '(A,A10,A11,A)', "  nu       ", "       ---", "        ---", nu_t_str
        else
            print '(A,F10.4,A11,A)', "  nu       ", sim_nu, "        ---", nu_t_str
        end if
        print '(A,F10.4,F11.4,F11.4)', "  uncond vol", uncvol_true, uncvol_n, uncvol_t
        print '(A)', ""

        ! estimation summary
        print '(A)', "                     Normal    Student-t"
        print '(A)', "  -------------   ---------   ---------"
        print '(A,F11.4,F12.4)',   "  -loglik/n  ", f_n,        f_t
        print '(A,I11,I12)',       "  iterations ", iter_n,     iter_t
        print '(A,ES11.3,ES12.3)', "  grad norm  ", norm2(g_n), norm2(g_t)
        print '(A,F11.6,F12.6)',   "  time (s)   ", t_normal,   t_student
        print '(A)', ""
        if (.not. conv_n) print '(A)', "  Note: Normal fit did not converge"
        if (.not. conv_t) print '(A)', "  Note: t fit did not converge"

    end do

end program xgarch_t

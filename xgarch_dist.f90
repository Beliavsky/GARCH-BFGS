! Simulate GARCH(1,1) data from several noise distributions, then fit seven models:
!   (1) GARCH-Normal    (3 parameters: omega, alpha, beta)
!   (2) GARCH-t         (4 parameters: omega, alpha, beta, nu)
!   (3) GARCH-Sech      (3 parameters: omega, alpha, beta)
!   (4) GARCH-GED       (4 parameters: omega, alpha, beta, nu)
!   (5) GARCH-Laplace   (3 parameters: omega, alpha, beta)
!   (6) GARCH-Logistic  (3 parameters: omega, alpha, beta)
!   (7) GARCH-NIG       (4 parameters: omega, alpha, beta, alp)

program xgarch_dist
    use kind_mod,             only: dp
    use garch_module,         only: garch_simulate, garch_set_data, garch_inv_transform, &
                                     garch_obj, garch_transform
    use garch_t_module,       only: garch_t_simulate, garch_t_set_data, &
                                     garch_t_inv_transform, garch_t_obj, garch_t_transform
    use garch_sech_module,    only: garch_sech_simulate, garch_sech_set_data, garch_sech_obj
    use garch_ged_module,     only: garch_ged_simulate, garch_ged_set_data, &
                                     garch_ged_obj, garch_ged_transform, garch_ged_inv_transform
    use garch_laplace_module,  only: garch_laplace_simulate, garch_laplace_set_data, &
                                     garch_laplace_obj
    use garch_logistic_module, only: garch_logistic_simulate, garch_logistic_set_data, &
                                     garch_logistic_obj
    use garch_nig_module,      only: garch_nig_simulate, garch_nig_set_data, &
                                     garch_nig_obj, garch_nig_transform, garch_nig_inv_transform
    use bfgs_module,          only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_obs    = 10000
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true GARCH(1,1) parameters (same for all simulations) ----
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.10_dp
    real(dp), parameter :: true_beta  = 0.85_dp

    ! ---- simulation distribution codes ----
    integer, parameter :: dist_normal   = 0
    integer, parameter :: dist_t        = 1
    integer, parameter :: dist_sech     = 2
    integer, parameter :: dist_ged      = 3
    integer, parameter :: dist_laplace  = 4
    integer, parameter :: dist_logistic = 5
    integer, parameter :: dist_nig      = 6

    ! ---- fitted distributions (order determines table columns) ----
    integer,  parameter :: n_fit = 7
    character(len=9), parameter :: fit_names(n_fit) = &
        ["Normal   ", "Student-t", "Sech     ", "GED      ", "Laplace  ", "Logistic ", "NIG      "]

    ! ---- simulation cases: type + shape parameter (ignored for Normal/Sech/Laplace/Logistic) ----
    integer,  parameter :: sim_type(*) = [ dist_sech,     dist_normal,  &
                                            dist_t,        dist_t,       &
                                            dist_ged,      dist_ged,     &
                                            dist_laplace,  dist_logistic,&
                                            dist_nig,      dist_nig      ]
    real(dp), parameter :: sim_par(*)  = [ 0.0_dp, 0.0_dp,  &
                                            5.0_dp, 15.0_dp, &
                                            1.0_dp, 1.5_dp,  &
                                            0.0_dp, 0.0_dp,  &
                                            2.0_dp, 5.0_dp   ]

    real(dp) :: y(n_obs)

    ! fit results
    real(dp) :: p_n(3), p_t(4), p_s(3), p_g(4), p_l(3), p_lo(3), p_ni(4)
    real(dp) :: omega_n,  alpha_n,  beta_n
    real(dp) :: omega_t,  alpha_t,  beta_t,  nu_t
    real(dp) :: omega_s,  alpha_s,  beta_s
    real(dp) :: omega_g,  alpha_g,  beta_g,  nu_g
    real(dp) :: omega_l,  alpha_l,  beta_l
    real(dp) :: omega_lo, alpha_lo, beta_lo
    real(dp) :: omega_ni, alpha_ni, beta_ni, alp_ni
    real(dp) :: f_n, f_t, f_s, f_g, f_l, f_lo, f_ni
    real(dp) :: g_n(3), g_t(4), g_s(3), g_g(4), g_l(3), g_lo(3), g_ni(4)
    logical  :: conv_n, conv_t, conv_s, conv_g, conv_l, conv_lo, conv_ni
    integer  :: iter_n, iter_t, iter_s, iter_g, iter_l, iter_lo, iter_ni
    integer  :: ck0, ck1, ck_rate
    real(dp) :: t_n, t_t, t_s, t_g, t_l, t_lo, t_ni

    ! derived
    real(dp) :: y_std
    real(dp) :: uncvol_true, uncvol_n, uncvol_t, uncvol_s, uncvol_g, &
                uncvol_l, uncvol_lo, uncvol_ni
    real(dp) :: aic_n, aic_t, aic_s, aic_g, aic_l, aic_lo, aic_ni
    real(dp) :: bic_n, bic_t, bic_s, bic_g, bic_l, bic_lo, bic_ni
    real(dp) :: lnn
    real(dp) :: par
    integer  :: i_sim
    integer  :: rank_n,  rank_t,  rank_s,  rank_g,  rank_l,  rank_lo,  rank_ni
    integer  :: arank_n, arank_t, arank_s, arank_g, arank_l, arank_lo, arank_ni
    integer  :: brank_n, brank_t, brank_s, brank_g, brank_l, brank_lo, brank_ni
    character(len=11) :: nu_t_str

    call system_clock(count_rate=ck_rate)

    uncvol_true = sqrt(true_omega / (1.0_dp - true_alpha - true_beta))
    lnn         = log(real(n_obs, dp))

    print '(A)', ""
    print '(A)', " GARCH(1,1): Normal / Student-t / Sech / GED / Laplace / Logistic / NIG"
    print '(A,I0)', " Observations : ", n_obs

    do i_sim = 1, size(sim_type)
        par = sim_par(i_sim)

        ! ---- simulate ----
        select case (sim_type(i_sim))
        case (dist_normal)
            call garch_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
        case (dist_t)
            call garch_t_simulate(true_omega, true_alpha, true_beta, par, n_obs, seed_val, y)
        case (dist_sech)
            call garch_sech_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
        case (dist_ged)
            call garch_ged_simulate(true_omega, true_alpha, true_beta, par, n_obs, seed_val, y)
        case (dist_laplace)
            call garch_laplace_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
        case (dist_logistic)
            call garch_logistic_simulate(true_omega, true_alpha, true_beta, n_obs, seed_val, y)
        case (dist_nig)
            call garch_nig_simulate(true_omega, true_alpha, true_beta, par, n_obs, seed_val, y)
        end select
        call garch_set_data(y, n_obs)
        call garch_t_set_data(y, n_obs)
        call garch_sech_set_data(y, n_obs)
        call garch_ged_set_data(y, n_obs)
        call garch_laplace_set_data(y, n_obs)
        call garch_logistic_set_data(y, n_obs)
        call garch_nig_set_data(y, n_obs)

        ! ---- fit GARCH-Normal ----
        call garch_inv_transform(0.05_dp, 0.10_dp, 0.80_dp, p_n)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_obj, p_n, 3, max_iter, gtol, f_n, iter_n, conv_n)
        call system_clock(count=ck1)
        t_n = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_transform(p_n, omega_n, alpha_n, beta_n)
        call garch_obj(p_n, 3, f_n, g_n)

        ! ---- fit GARCH-t, warm-started from Normal estimates ----
        call garch_t_inv_transform(omega_n, alpha_n, beta_n, 8.0_dp, p_t)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_t_obj, p_t, 4, max_iter, gtol, f_t, iter_t, conv_t)
        call system_clock(count=ck1)
        t_t = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_t_transform(p_t, omega_t, alpha_t, beta_t, nu_t)
        call garch_t_obj(p_t, 4, f_t, g_t)

        ! ---- fit GARCH-Sech, warm-started from Normal estimates ----
        p_s = p_n
        call system_clock(count=ck0)
        call bfgs_minimize(garch_sech_obj, p_s, 3, max_iter, gtol, f_s, iter_s, conv_s)
        call system_clock(count=ck1)
        t_s = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_transform(p_s, omega_s, alpha_s, beta_s)
        call garch_sech_obj(p_s, 3, f_s, g_s)

        ! ---- fit GARCH-GED, warm-started from Normal estimates + nu=2 ----
        call garch_ged_inv_transform(omega_n, alpha_n, beta_n, 2.0_dp, p_g)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_ged_obj, p_g, 4, max_iter, gtol, f_g, iter_g, conv_g)
        call system_clock(count=ck1)
        t_g = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_ged_transform(p_g, omega_g, alpha_g, beta_g, nu_g)
        call garch_ged_obj(p_g, 4, f_g, g_g)

        ! ---- fit GARCH-Laplace, warm-started from Normal estimates ----
        p_l = p_n
        call system_clock(count=ck0)
        call bfgs_minimize(garch_laplace_obj, p_l, 3, max_iter, gtol, f_l, iter_l, conv_l)
        call system_clock(count=ck1)
        t_l = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_transform(p_l, omega_l, alpha_l, beta_l)
        call garch_laplace_obj(p_l, 3, f_l, g_l)

        ! ---- fit GARCH-Logistic, warm-started from Normal estimates ----
        p_lo = p_n
        call system_clock(count=ck0)
        call bfgs_minimize(garch_logistic_obj, p_lo, 3, max_iter, gtol, f_lo, iter_lo, conv_lo)
        call system_clock(count=ck1)
        t_lo = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_transform(p_lo, omega_lo, alpha_lo, beta_lo)
        call garch_logistic_obj(p_lo, 3, f_lo, g_lo)

        ! ---- fit GARCH-NIG, warm-started from Normal estimates + alp=3 ----
        call garch_nig_inv_transform(omega_n, alpha_n, beta_n, 3.0_dp, p_ni)   ! alp=3 in (0.1,20)
        call system_clock(count=ck0)
        call bfgs_minimize(garch_nig_obj, p_ni, 4, max_iter, gtol, f_ni, iter_ni, conv_ni)
        call system_clock(count=ck1)
        t_ni = real(ck1 - ck0, dp) / real(ck_rate, dp)
        call garch_nig_transform(p_ni, omega_ni, alpha_ni, beta_ni, alp_ni)
        call garch_nig_obj(p_ni, 4, f_ni, g_ni)

        ! ---- derived quantities ----
        y_std      = sqrt(sum((y - sum(y)/n_obs)**2) / (n_obs - 1))
        uncvol_n   = sqrt(omega_n  / (1.0_dp - alpha_n  - beta_n))
        uncvol_t   = sqrt(omega_t  / (1.0_dp - alpha_t  - beta_t))
        uncvol_s   = sqrt(omega_s  / (1.0_dp - alpha_s  - beta_s))
        uncvol_g   = sqrt(omega_g  / (1.0_dp - alpha_g  - beta_g))
        uncvol_l   = sqrt(omega_l  / (1.0_dp - alpha_l  - beta_l))
        uncvol_lo  = sqrt(omega_lo / (1.0_dp - alpha_lo - beta_lo))
        uncvol_ni  = sqrt(omega_ni / (1.0_dp - alpha_ni - beta_ni))

        rank_n  = 1 + count([f_t, f_s, f_g, f_l, f_lo, f_ni] < f_n)
        rank_t  = 1 + count([f_n, f_s, f_g, f_l, f_lo, f_ni] < f_t)
        rank_s  = 1 + count([f_n, f_t, f_g, f_l, f_lo, f_ni] < f_s)
        rank_g  = 1 + count([f_n, f_t, f_s, f_l, f_lo, f_ni] < f_g)
        rank_l  = 1 + count([f_n, f_t, f_s, f_g, f_lo, f_ni] < f_l)
        rank_lo = 1 + count([f_n, f_t, f_s, f_g, f_l,  f_ni] < f_lo)
        rank_ni = 1 + count([f_n, f_t, f_s, f_g, f_l,  f_lo] < f_ni)

        aic_n  = 2*3 + 2*n_obs*f_n
        aic_t  = 2*4 + 2*n_obs*f_t
        aic_s  = 2*3 + 2*n_obs*f_s
        aic_g  = 2*4 + 2*n_obs*f_g
        aic_l  = 2*3 + 2*n_obs*f_l
        aic_lo = 2*3 + 2*n_obs*f_lo
        aic_ni = 2*4 + 2*n_obs*f_ni
        arank_n  = 1 + count([aic_t, aic_s, aic_g, aic_l, aic_lo, aic_ni] < aic_n)
        arank_t  = 1 + count([aic_n, aic_s, aic_g, aic_l, aic_lo, aic_ni] < aic_t)
        arank_s  = 1 + count([aic_n, aic_t, aic_g, aic_l, aic_lo, aic_ni] < aic_s)
        arank_g  = 1 + count([aic_n, aic_t, aic_s, aic_l, aic_lo, aic_ni] < aic_g)
        arank_l  = 1 + count([aic_n, aic_t, aic_s, aic_g, aic_lo, aic_ni] < aic_l)
        arank_lo = 1 + count([aic_n, aic_t, aic_s, aic_g, aic_l,  aic_ni] < aic_lo)
        arank_ni = 1 + count([aic_n, aic_t, aic_s, aic_g, aic_l,  aic_lo] < aic_ni)

        bic_n  = 3*lnn + 2*n_obs*f_n
        bic_t  = 4*lnn + 2*n_obs*f_t
        bic_s  = 3*lnn + 2*n_obs*f_s
        bic_g  = 4*lnn + 2*n_obs*f_g
        bic_l  = 3*lnn + 2*n_obs*f_l
        bic_lo = 3*lnn + 2*n_obs*f_lo
        bic_ni = 4*lnn + 2*n_obs*f_ni
        brank_n  = 1 + count([bic_t, bic_s, bic_g, bic_l, bic_lo, bic_ni] < bic_n)
        brank_t  = 1 + count([bic_n, bic_s, bic_g, bic_l, bic_lo, bic_ni] < bic_t)
        brank_s  = 1 + count([bic_n, bic_t, bic_g, bic_l, bic_lo, bic_ni] < bic_s)
        brank_g  = 1 + count([bic_n, bic_t, bic_s, bic_l, bic_lo, bic_ni] < bic_g)
        brank_l  = 1 + count([bic_n, bic_t, bic_s, bic_g, bic_lo, bic_ni] < bic_l)
        brank_lo = 1 + count([bic_n, bic_t, bic_s, bic_g, bic_l,  bic_ni] < bic_lo)
        brank_ni = 1 + count([bic_n, bic_t, bic_s, bic_g, bic_l,  bic_lo] < bic_ni)

        if (nu_t > 99.0_dp) then
            nu_t_str = "      large"
        else
            write(nu_t_str, '(F11.4)') nu_t
        end if

        ! ---- results ----
        print '(A)', ""
        select case (sim_type(i_sim))
        case (dist_normal)
            print '(A)', " Simulated noise: Normal"
        case (dist_t)
            print '(A,F6.1)', " Simulated noise: Student-t, nu =", par
        case (dist_sech)
            print '(A)', " Simulated noise: Sech"
        case (dist_ged)
            print '(A,F6.1)', " Simulated noise: GED, nu =", par
        case (dist_laplace)
            print '(A)', " Simulated noise: Laplace"
        case (dist_logistic)
            print '(A)', " Simulated noise: Logistic"
        case (dist_nig)
            print '(A,F6.1)', " Simulated noise: NIG, alp =", par
        end select
        print '(A,F10.4)', " Simulated data std dev: ", y_std
        print '(A)', ""

        ! parameter table: label(11)+True(10)+7*F11.4
        print '(A)', "                 True     Normal  Student-t       Sech        GED    Laplace   Logistic        NIG"
        print '(A)', "  ---------  --------   --------   --------   --------   --------   --------   --------   --------"
        print '(A,F10.4,7F11.4)', "  omega    ", true_omega, omega_n, omega_t, omega_s, omega_g, omega_l, omega_lo, omega_ni
        print '(A,F10.4,7F11.4)', "  alpha    ", true_alpha, alpha_n, alpha_t, alpha_s, alpha_g, alpha_l, alpha_lo, alpha_ni
        print '(A,F10.4,7F11.4)', "  beta     ", true_beta,  beta_n,  beta_t,  beta_s,  beta_g,  beta_l,  beta_lo,  beta_ni
        if (sim_type(i_sim) == dist_t) then
            print '(A,F10.4,A11,A,5A11)', &
                "  t nu     ", par, "        ---", nu_t_str, &
                "        ---", "        ---", "        ---", "        ---", "        ---"
        else
            print '(A,A10,A11,A,5A11)', &
                "  t nu     ", "       ---", "        ---", nu_t_str, &
                "        ---", "        ---", "        ---", "        ---", "        ---"
        end if
        if (sim_type(i_sim) == dist_ged) then
            print '(A,F10.4,3A11,F11.4,3A11)', &
                "  GED nu   ", par, "        ---", "        ---", "        ---", nu_g, &
                "        ---", "        ---", "        ---"
        else
            print '(A,A10,3A11,F11.4,3A11)', &
                "  GED nu   ", "       ---", "        ---", "        ---", "        ---", nu_g, &
                "        ---", "        ---", "        ---"
        end if
        if (sim_type(i_sim) == dist_nig) then
            print '(A,F10.4,6A11,F11.4)', &
                "  NIG alp  ", par, "        ---", "        ---", "        ---", "        ---", "        ---", "        ---", alp_ni
        else
            print '(A,A10,6A11,F11.4)', &
                "  NIG alp  ", "       ---", "        ---", "        ---", "        ---", "        ---", "        ---", "        ---", alp_ni
        end if
        print '(A,F10.4,7F11.4)', "  uncond vol", uncvol_true, uncvol_n, uncvol_t, uncvol_s, uncvol_g, uncvol_l, uncvol_lo, uncvol_ni
        print '(A)', ""

        ! estimation summary: label(13)+Normal(11)+6*F12
        print '(A)', "                  Normal   Student-t        Sech         GED     Laplace    Logistic         NIG"
        print '(A)', "  -------------   -------   ---------   ---------   ---------   ---------   ---------   ---------"
        print '(A,I11,6I12,2X,A)', &
            "  rank       ", rank_n,  rank_t,  rank_s,  rank_g,  rank_l,  rank_lo,  rank_ni, &
            fit_names(minloc([rank_n,  rank_t,  rank_s,  rank_g,  rank_l,  rank_lo,  rank_ni],  1))
        print '(A,I11,6I12,2X,A)', &
            "  AIC rank   ", arank_n, arank_t, arank_s, arank_g, arank_l, arank_lo, arank_ni, &
            fit_names(minloc([arank_n, arank_t, arank_s, arank_g, arank_l, arank_lo, arank_ni], 1))
        print '(A,I11,6I12,2X,A)', &
            "  BIC rank   ", brank_n, brank_t, brank_s, brank_g, brank_l, brank_lo, brank_ni, &
            fit_names(minloc([brank_n, brank_t, brank_s, brank_g, brank_l, brank_lo, brank_ni], 1))
        print '(A,F11.4,6F12.4)', &
            "  -loglik/n  ", f_n,   f_t,   f_s,   f_g,   f_l,   f_lo,   f_ni
        print '(A,F11.1,6F12.1)', &
            "  AIC        ", aic_n,  aic_t,  aic_s,  aic_g,  aic_l,  aic_lo,  aic_ni
        print '(A,F11.1,6F12.1)', &
            "  BIC        ", bic_n,  bic_t,  bic_s,  bic_g,  bic_l,  bic_lo,  bic_ni
        print '(A,I11,6I12)', &
            "  iterations ", iter_n, iter_t, iter_s, iter_g, iter_l, iter_lo, iter_ni
        print '(A,ES11.3,6ES12.3)', &
            "  grad norm  ", norm2(g_n), norm2(g_t), norm2(g_s), norm2(g_g), norm2(g_l), norm2(g_lo), norm2(g_ni)
        print '(A,F11.6,6F12.6)', &
            "  time (s)   ", t_n,   t_t,   t_s,   t_g,   t_l,   t_lo,   t_ni
        print '(A)', ""
        if (.not. conv_n)  print '(A)', "  Note: Normal fit did not converge"
        if (.not. conv_t)  print '(A)', "  Note: t fit did not converge"
        if (.not. conv_s)  print '(A)', "  Note: Sech fit did not converge"
        if (.not. conv_g)  print '(A)', "  Note: GED fit did not converge"
        if (.not. conv_l)  print '(A)', "  Note: Laplace fit did not converge"
        if (.not. conv_lo) print '(A)', "  Note: Logistic fit did not converge"
        if (.not. conv_ni) print '(A)', "  Note: NIG fit did not converge"

    end do

end program xgarch_dist

! GARCH(1,1) fitted to multiple noise distributions with a user-selectable distribution set.
! For each enabled simulation distribution, nsim replications are run.
! Per-replication results are printed followed by a summary over all replications.
! Configuration, derived types, and print subroutines are in garch_choose_mod.

program xgarch_choose_dist
    use kind_mod,              only: dp
    use garch_choose_mod,      only: garch_config_t, fit_result_t, acc_t, &
                                      n_fit, np_fit, sim_type, sim_par, &
                                      i_normal, i_t, i_sech, i_ged, i_laplace, i_logistic, i_nig, &
                                      dist_normal, dist_t, dist_sech, dist_ged, &
                                      dist_laplace, dist_logistic, dist_nig, &
                                      print_results, print_summary, print_selection_matrices
    use garch_module,          only: garch_simulate, garch_set_data, &
                                      garch_inv_transform, garch_obj, garch_transform
    use garch_t_module,        only: garch_t_simulate, garch_t_set_data, &
                                      garch_t_inv_transform, garch_t_obj, garch_t_transform
    use garch_sech_module,     only: garch_sech_simulate, garch_sech_set_data, garch_sech_obj
    use garch_ged_module,      only: garch_ged_simulate, garch_ged_set_data, &
                                      garch_ged_obj, garch_ged_transform, garch_ged_inv_transform
    use garch_laplace_module,  only: garch_laplace_simulate, garch_laplace_set_data, &
                                      garch_laplace_obj
    use garch_logistic_module, only: garch_logistic_simulate, garch_logistic_set_data, &
                                      garch_logistic_obj
    use garch_nig_module,      only: garch_nig_simulate, garch_nig_set_data, &
                                      garch_nig_obj, garch_nig_transform, garch_nig_inv_transform
    use bfgs_module,           only: bfgs_minimize
    implicit none

    ! ---- which distributions to simulate from (indexed 0:6 matching dist_* constants) ----
    logical :: sim_on(0:6) = .false.

    ! ---- configuration: override any fields before the simulation loop ----
    type(garch_config_t) :: cfg
    type(fit_result_t) :: res(n_fit)
    type(acc_t)        :: acc
    type(acc_t)        :: acc_all(size(sim_type))         ! one entry per active simulation case
    integer  :: active_sim_types(size(sim_type))          ! dist_* code for each active case
    real(dp) :: active_sim_pars(size(sim_type))           ! shape parameter for each active case
    integer  :: n_active                                  ! number of active simulation cases run
    real(dp), allocatable :: y(:)
    real(dp) :: p(4), f, g(4)
    real(dp) :: omega_ws, alpha_ws, beta_ws
    real(dp) :: y_std, par, lnn
    integer  :: ck0, ck1, ck_rate, ck_start, ck_end
    integer  :: i_sim, i_rep, i, j, seed_rep

    call system_clock(count_rate=ck_rate)
    call system_clock(count=ck_start)
    n_active = 0

    ! ---- select distributions to simulate from ----
    sim_on([dist_normal, dist_t]) = .true.
    cfg%nobs = 1000 ! # of observations in each simulation
    cfg%nsim =    2 ! # of simulations
    cfg%fit_on(i_nig) = .false.
    allocate(y(cfg%nobs))
    lnn = log(real(cfg%nobs, dp))

    print '(A)', ""
    print '(A)', " GARCH(1,1) with selectable noise distributions"
    print '(A,I0)', " Observations : ", cfg%nobs
    print '(A,I0)', " Replications : ", cfg%nsim

    do i_sim = 1, size(sim_type)
        if (.not. sim_on(sim_type(i_sim))) cycle
        par = sim_par(i_sim)

        ! zero accumulators for this simulation case
        acc%omega      = 0.0_dp;  acc%alpha      = 0.0_dp;  acc%beta       = 0.0_dp
        acc%extra_par  = 0.0_dp;  acc%uncond_vol = 0.0_dp;  acc%nll_n      = 0.0_dp
        acc%time_s     = 0.0_dp;  acc%iterations = 0
        acc%rank1_ll   = 0;       acc%rank1_aic  = 0;       acc%rank1_bic  = 0
        acc%n_conv     = 0

        do i_rep = 1, cfg%nsim
            seed_rep = cfg%seed_val + i_rep - 1

            ! ---- simulate ----
            select case (sim_type(i_sim))
            case (dist_normal)
                call garch_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, cfg%nobs, seed_rep, y)
            case (dist_t)
                call garch_t_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, par, cfg%nobs, seed_rep, y)
            case (dist_sech)
                call garch_sech_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, cfg%nobs, seed_rep, y)
            case (dist_ged)
                call garch_ged_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, par, cfg%nobs, seed_rep, y)
            case (dist_laplace)
                call garch_laplace_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, cfg%nobs, seed_rep, y)
            case (dist_logistic)
                call garch_logistic_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, cfg%nobs, seed_rep, y)
            case (dist_nig)
                call garch_nig_simulate(cfg%true_omega, cfg%true_alpha, cfg%true_beta, par, cfg%nobs, seed_rep, y)
            end select

            call garch_set_data(y, cfg%nobs)
            call garch_t_set_data(y, cfg%nobs)
            call garch_sech_set_data(y, cfg%nobs)
            call garch_ged_set_data(y, cfg%nobs)
            call garch_laplace_set_data(y, cfg%nobs)
            call garch_logistic_set_data(y, cfg%nobs)
            call garch_nig_set_data(y, cfg%nobs)

            y_std = sqrt(sum((y - sum(y)/cfg%nobs)**2) / (cfg%nobs - 1))

            do i = 1, n_fit
                res(i)%fitted = .false.
            end do

            omega_ws = 0.05_dp;  alpha_ws = 0.10_dp;  beta_ws = 0.80_dp

            ! ---- fit Normal ----
            if (cfg%fit_on(i_normal)) then
                call garch_inv_transform(0.05_dp, 0.10_dp, 0.80_dp, p(1:3))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_obj, p(1:3), 3, cfg%max_iter, cfg%gtol, f, &
                                    res(i_normal)%iterations, res(i_normal)%converged)
                call system_clock(count=ck1)
                call garch_transform(p(1:3), res(i_normal)%omega, res(i_normal)%alpha, res(i_normal)%beta)
                call garch_obj(p(1:3), 3, f, g(1:3))
                res(i_normal)%nll_n     = f
                res(i_normal)%grad_norm = norm2(g(1:3))
                res(i_normal)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_normal)%fitted    = .true.
                omega_ws = res(i_normal)%omega
                alpha_ws = res(i_normal)%alpha
                beta_ws  = res(i_normal)%beta
            end if

            ! ---- fit Student-t ----
            if (cfg%fit_on(i_t)) then
                call garch_t_inv_transform(omega_ws, alpha_ws, beta_ws, 8.0_dp, p(1:4))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_t_obj, p(1:4), 4, cfg%max_iter, cfg%gtol, f, &
                                    res(i_t)%iterations, res(i_t)%converged)
                call system_clock(count=ck1)
                call garch_t_transform(p(1:4), res(i_t)%omega, res(i_t)%alpha, res(i_t)%beta, &
                                       res(i_t)%extra_par)
                call garch_t_obj(p(1:4), 4, f, g(1:4))
                res(i_t)%nll_n     = f
                res(i_t)%grad_norm = norm2(g(1:4))
                res(i_t)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_t)%fitted    = .true.
            end if

            ! ---- fit Sech ----
            if (cfg%fit_on(i_sech)) then
                call garch_inv_transform(omega_ws, alpha_ws, beta_ws, p(1:3))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_sech_obj, p(1:3), 3, cfg%max_iter, cfg%gtol, f, &
                                    res(i_sech)%iterations, res(i_sech)%converged)
                call system_clock(count=ck1)
                call garch_transform(p(1:3), res(i_sech)%omega, res(i_sech)%alpha, res(i_sech)%beta)
                call garch_sech_obj(p(1:3), 3, f, g(1:3))
                res(i_sech)%nll_n     = f
                res(i_sech)%grad_norm = norm2(g(1:3))
                res(i_sech)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_sech)%fitted    = .true.
            end if

            ! ---- fit GED ----
            if (cfg%fit_on(i_ged)) then
                call garch_ged_inv_transform(omega_ws, alpha_ws, beta_ws, 2.0_dp, p(1:4))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_ged_obj, p(1:4), 4, cfg%max_iter, cfg%gtol, f, &
                                    res(i_ged)%iterations, res(i_ged)%converged)
                call system_clock(count=ck1)
                call garch_ged_transform(p(1:4), res(i_ged)%omega, res(i_ged)%alpha, res(i_ged)%beta, &
                                         res(i_ged)%extra_par)
                call garch_ged_obj(p(1:4), 4, f, g(1:4))
                res(i_ged)%nll_n     = f
                res(i_ged)%grad_norm = norm2(g(1:4))
                res(i_ged)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_ged)%fitted    = .true.
            end if

            ! ---- fit Laplace ----
            if (cfg%fit_on(i_laplace)) then
                call garch_inv_transform(omega_ws, alpha_ws, beta_ws, p(1:3))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_laplace_obj, p(1:3), 3, cfg%max_iter, cfg%gtol, f, &
                                    res(i_laplace)%iterations, res(i_laplace)%converged)
                call system_clock(count=ck1)
                call garch_transform(p(1:3), res(i_laplace)%omega, res(i_laplace)%alpha, res(i_laplace)%beta)
                call garch_laplace_obj(p(1:3), 3, f, g(1:3))
                res(i_laplace)%nll_n     = f
                res(i_laplace)%grad_norm = norm2(g(1:3))
                res(i_laplace)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_laplace)%fitted    = .true.
            end if

            ! ---- fit Logistic ----
            if (cfg%fit_on(i_logistic)) then
                call garch_inv_transform(omega_ws, alpha_ws, beta_ws, p(1:3))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_logistic_obj, p(1:3), 3, cfg%max_iter, cfg%gtol, f, &
                                    res(i_logistic)%iterations, res(i_logistic)%converged)
                call system_clock(count=ck1)
                call garch_transform(p(1:3), res(i_logistic)%omega, res(i_logistic)%alpha, res(i_logistic)%beta)
                call garch_logistic_obj(p(1:3), 3, f, g(1:3))
                res(i_logistic)%nll_n     = f
                res(i_logistic)%grad_norm = norm2(g(1:3))
                res(i_logistic)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_logistic)%fitted    = .true.
            end if

            ! ---- fit NIG ----
            if (cfg%fit_on(i_nig)) then
                call garch_nig_inv_transform(omega_ws, alpha_ws, beta_ws, 3.0_dp, p(1:4))
                call system_clock(count=ck0)
                call bfgs_minimize(garch_nig_obj, p(1:4), 4, cfg%max_iter, cfg%gtol, f, &
                                    res(i_nig)%iterations, res(i_nig)%converged)
                call system_clock(count=ck1)
                call garch_nig_transform(p(1:4), res(i_nig)%omega, res(i_nig)%alpha, res(i_nig)%beta, &
                                         res(i_nig)%extra_par)
                call garch_nig_obj(p(1:4), 4, f, g(1:4))
                res(i_nig)%nll_n     = f
                res(i_nig)%grad_norm = norm2(g(1:4))
                res(i_nig)%time_s    = real(ck1 - ck0, dp) / real(ck_rate, dp)
                res(i_nig)%fitted    = .true.
            end if

            ! ---- uncond vol and ranks ----
            do i = 1, n_fit
                if (.not. res(i)%fitted) cycle
                res(i)%uncond_vol = sqrt(res(i)%omega / (1.0_dp - res(i)%alpha - res(i)%beta))
            end do
            do i = 1, n_fit
                if (.not. res(i)%fitted) cycle
                res(i)%rank_ll  = 1;  res(i)%rank_aic = 1;  res(i)%rank_bic = 1
                do j = 1, n_fit
                    if (.not. res(j)%fitted) cycle
                    if (res(j)%nll_n < res(i)%nll_n) res(i)%rank_ll = res(i)%rank_ll + 1
                    if (res(j)%nll_n + (np_fit(j) - np_fit(i)) * 2.0_dp / (2*cfg%nobs) &
                        < res(i)%nll_n) res(i)%rank_aic = res(i)%rank_aic + 1
                    if (res(j)%nll_n + (np_fit(j) - np_fit(i)) * lnn / (2*cfg%nobs) &
                        < res(i)%nll_n) res(i)%rank_bic = res(i)%rank_bic + 1
                end do
            end do

            call print_results(res, cfg, sim_type(i_sim), par, y_std, i_rep)

            ! ---- accumulate ----
            do i = 1, n_fit
                if (.not. res(i)%fitted) cycle
                acc%omega(i)      = acc%omega(i)      + res(i)%omega
                acc%alpha(i)      = acc%alpha(i)      + res(i)%alpha
                acc%beta(i)       = acc%beta(i)       + res(i)%beta
                acc%extra_par(i)  = acc%extra_par(i)  + res(i)%extra_par
                acc%uncond_vol(i) = acc%uncond_vol(i) + res(i)%uncond_vol
                acc%nll_n(i)      = acc%nll_n(i)      + res(i)%nll_n
                acc%time_s(i)     = acc%time_s(i)     + res(i)%time_s
                acc%iterations(i) = acc%iterations(i) + res(i)%iterations
                if (res(i)%rank_ll  == 1) acc%rank1_ll(i)  = acc%rank1_ll(i)  + 1
                if (res(i)%rank_aic == 1) acc%rank1_aic(i) = acc%rank1_aic(i) + 1
                if (res(i)%rank_bic == 1) acc%rank1_bic(i) = acc%rank1_bic(i) + 1
                if (res(i)%converged)     acc%n_conv(i)    = acc%n_conv(i)    + 1
            end do

        end do  ! i_rep

        call print_summary(acc, cfg, sim_type(i_sim), par)

        n_active = n_active + 1
        acc_all(n_active)          = acc
        active_sim_types(n_active) = sim_type(i_sim)
        active_sim_pars(n_active)  = par

    end do  ! i_sim

    call print_selection_matrices(acc_all(1:n_active), n_active, cfg, &
                                   active_sim_types(1:n_active), active_sim_pars(1:n_active))

    call system_clock(count=ck_end)
    print '(A,F8.3,A)', " Total time: ", real(ck_end - ck_start, dp) / real(ck_rate, dp), " s"

end program xgarch_choose_dist

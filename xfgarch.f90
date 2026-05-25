! Fit family GARCH models from fgarch.f90 to the ETF return panel.
! The default is a NAGARCH-like twist model; pass --full to estimate the full
! Hentschel FGARCH-Normal model with lambda and nu free.

program xfgarch
    use kind_mod,   only: dp
    use csv_mod,    only: read_price_csv, print_price_sample_info
    use stats_mod,  only: mean, sd
    use garch_types_mod, only: garch_params_t
    use fgarch_mod, only: fg_dist_normal, fgarch_set_data, fgarch_set_dist, &
                          fgarch_np, fgarch_obj, fgarch_transform, fgarch_inv_transform, &
                          fgarch_vol_ann, fgarch_skew_kurt
    use garch_fit_mod, only: fit_fgarch_twist
    use bfgs_mod,   only: bfgs_minimize
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    integer,  parameter :: n_dist = 1
    integer,  parameter :: n_start = 1

    integer,           allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp),          allocatable :: prices(:,:), ret(:)
    integer,           parameter :: dist_ids(n_dist) = [fg_dist_normal]
    character(len=9),   parameter :: dist_names(n_dist) = ["Normal   "]

    integer :: nprices, ncols, nobs, icol, idist
    real(dp) :: ret_mean, ret_std
    real(dp) :: fopt, omega, alpha, beta, lam, nu, theta, twist, nu_t
    real(dp) :: vol_ann, logl, aic, bic, skew, ekurt
    integer  :: niter, np
    integer  :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s
    logical  :: converged
    logical  :: fit_full
    character(len=32) :: arg, model_name
    type(garch_params_t) :: params

    call system_clock(clock_start, clock_rate)
    fit_full = .false.
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        fit_full = trim(arg) == "--full" .or. trim(arg) == "full"
    end if

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))

    call print_price_sample_info(prices_file, dates, ncols)
    if (fit_full) then
        write(*,'(A)') "Mode: full FGARCH-Normal"
        write(*,'(A)') "Asset     Model       Dist          omega   alpha    beta     lam      nu   theta   twist    nu_t  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt"
        write(*,'(A)') repeat("-", 170)
    else
        write(*,'(A)') "Mode: NAGARCH-twist FGARCH-Normal; c=0 nests NAGARCH"
        write(*,'(A)') "Asset          omega   alpha    beta   theta   twist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt"
        write(*,'(A)') repeat("-", 130)
    end if

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        ret_std = sd(ret)
        call fgarch_set_data(ret, nobs)

        do idist = 1, n_dist
            if (fit_full) then
                call fit_fgarch(dist_ids(idist), ret_std, fopt, omega, alpha, beta, lam, nu, theta, twist, nu_t, &
                                vol_ann, skew, ekurt, niter, converged)
                np = 7
                model_name = "FGARCH"
            else
                call fit_fgarch_twist(ret, ret_std, max_iter, gtol, fopt, params, vol_ann, skew, ekurt, &
                                      niter, converged)
                omega = params%omega
                alpha = params%alpha
                beta = params%beta
                theta = params%theta
                twist = params%twist
                lam = 2.0_dp
                nu = 2.0_dp
                nu_t = 0.0_dp
                np = 5
                model_name = "FGTWIST"
            end if
            logl = -real(nobs, dp) * fopt
            aic  = 2.0_dp * real(np, dp) - 2.0_dp * logl
            bic  = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl

            if (fit_full) then
                write(*,'(A9,1X,A10,1X,A9,ES12.3,6F8.4,F8.2,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3)') &
                    trim(col_names(icol)), trim(model_name), dist_names(idist), omega, alpha, beta, lam, nu, &
                    theta, twist, nu_t, vol_ann, logl, aic, bic, niter, converged, skew, ekurt
            else
                write(*,'(A9,ES12.3,4F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3)') &
                    trim(col_names(icol)), omega, alpha, beta, theta, twist, &
                    vol_ann, logl, aic, bic, niter, converged, skew, ekurt
            end if
        end do
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    subroutine fit_fgarch(idist, ret_std, f_best, omega_best, alpha_best, beta_best, lam_best, &
                          nu_best, b_best, c_best, nu_t_best, vol_ann_best, skew_best, ekurt_best, &
                          niter_best, converged_best)
        integer,  intent(in)  :: idist
        real(dp), intent(in)  :: ret_std
        real(dp), intent(out) :: f_best, omega_best, alpha_best, beta_best, lam_best, nu_best
        real(dp), intent(out) :: b_best, c_best, nu_t_best, vol_ann_best, skew_best, ekurt_best
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best

        real(dp) :: start_alpha(n_start), start_beta(n_start), start_lam(n_start)
        real(dp) :: start_nu(n_start), start_b(n_start), start_c(n_start)
        real(dp) :: omega0, f_try, p0_full(8)
        real(dp), allocatable :: p(:), p0(:), p_best(:)
        integer :: istart, np, niter_try
        logical :: converged_try

        start_alpha = [0.05_dp]
        start_beta  = [0.90_dp]
        start_lam   = [2.00_dp]
        start_nu    = [2.00_dp]
        start_b     = [0.00_dp]
        start_c     = [0.00_dp]

        call fgarch_set_dist(idist)
        np = fgarch_np()
        allocate(p(np), p0(np), p_best(np))

        f_best = huge(1.0_dp)
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = max(1.0e-8_dp, ret_std * (1.0_dp - start_beta(istart) - 0.8_dp*start_alpha(istart)))
            call fgarch_inv_transform(omega0, start_alpha(istart), start_beta(istart), &
                                      start_lam(istart), start_nu(istart), start_b(istart), &
                                      start_c(istart), 8.0_dp, p0_full)
            p0 = p0_full(1:np)
            p = p0
            call bfgs_minimize(fgarch_obj, p, np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        call fgarch_transform(p_best, np, omega_best, alpha_best, beta_best, lam_best, &
                              nu_best, b_best, c_best, nu_t_best)
        call fgarch_vol_ann(p_best, np, vol_ann_best)
        call fgarch_skew_kurt(p_best, np, skew_best, ekurt_best)

        deallocate(p, p0, p_best)
    end subroutine fit_fgarch

end program xfgarch

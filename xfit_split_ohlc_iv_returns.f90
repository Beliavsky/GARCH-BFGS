! Fit daily GARCH-family variance models using a split CO/OC return likelihood.
! The news impact curve uses close-close returns, while co_frac allocates daily
! conditional variance between close-to-open and open-to-close returns.

program xfit_split_ohlc_iv_returns
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, sqrt2
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv, read_ohlc_csv, print_price_sample_info
    use stats_mod,      only: mean
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t
    use bfgs_mod,       only: bfgs_minimize
    use vol_forecast_compare_mod, only: print_implied_vol_correlations
    implicit none

    character(len=*), parameter :: prices_file = "prices_ohlc.csv"
    character(len=*), parameter :: implied_vol_file = "vix_spy.csv"
    integer, parameter :: symbol_len = 16
    integer, parameter :: model_len = 16
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 120
    real(dp), parameter :: gtol = 1.0e-6_dp
    logical,  parameter :: fit_all_assets = .true.
    logical,  parameter :: flip_log_return_sign = .true.
    logical,  parameter :: standardize_log_return_by_prior_iv = .true.
    character(len=model_len), parameter :: models(*) = [character(len=model_len) :: &
        "SPLIT_SYMM", "SPLIT_NAGARCH", "SPLIT_GJR", "SPLIT_FGTWIST"]
    character(len=symbol_len), parameter :: fit_assets(*) = [character(len=symbol_len) :: "SPY"]
    character(len=symbol_len), parameter :: asset_iv_assets(*) = [character(len=symbol_len) :: "SPY"]
    character(len=symbol_len), parameter :: asset_iv_indices(*) = [character(len=symbol_len) :: "VIX"]
    character(len=model_len), parameter :: extra_series_names(*) = [character(len=model_len) :: "LOG_RETURN"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:), iv_dates(:)
    character(len=32), allocatable :: col_names(:), iv_col_names(:)
    real(dp), allocatable :: prices(:,:), open_prices(:,:), high_prices(:,:), low_prices(:,:), iv_values(:,:)
    real(dp), allocatable :: ret(:), ret_co(:), ret_oc(:), vol_forecast(:), vol_forecasts(:,:,:)
    real(dp), allocatable :: extra_series(:,:,:)
    real(dp) :: extra_signs(size(extra_series_names))
    character(len=model_len) :: fit_model_names(n_model)
    logical :: fit_model_have(n_model)
    integer :: nprices, ncols, nobs, icol, imod, clock_start, clock_end, clock_rate
    integer :: niter, np
    real(dp) :: fopt, logl_fit, logl_cc, aic, bic, persist, vol_ann, co_frac, skew, ekurt, elapsed_s
    type(garch_params_t) :: params
    logical :: converged

    character(len=model_len) :: active_model
    real(dp), allocatable :: active_cc(:), active_co(:), active_oc(:)

    call system_clock(clock_start, clock_rate)
    call nagarch_set_news_impact(.false.)
    fit_model_names = ""
    fit_model_have = .false.

    if (fit_all_assets) then
        call read_ohlc_csv(prices_file, dates, col_names, open_prices, prices, high_prices=high_prices, low_prices=low_prices)
    else
        call read_ohlc_csv(prices_file, dates, col_names, open_prices, prices, fit_assets, high_prices, low_prices)
    end if
    call read_price_csv(implied_vol_file, iv_dates, iv_col_names, iv_values)

    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), ret_co(nobs), ret_oc(nobs), vol_forecast(nobs), vol_forecasts(nobs,ncols,n_model))
    allocate(active_cc(nobs), active_co(nobs), active_oc(nobs))
    allocate(extra_series(nprices,ncols,size(extra_series_names)))
    vol_forecasts = 0.0_dp
    extra_series(:,:,1) = prices
    extra_signs = 1.0_dp
    if (flip_log_return_sign) extra_signs(1) = -1.0_dp

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A,A)') "Implied vol file: ", trim(implied_vol_file)
    call print_fit_assets()
    write(*,'(A)') "Split likelihood: h_co = co_frac*h_daily, h_oc = (1-co_frac)*h_daily; news impact uses CC returns."
    write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist co_frac persist  vol_ann%   logL_fit     logL_cc         AIC         BIC  iter conv    skew   ekurt"
    write(*,'(A)') repeat("-", 176)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret = ret - mean(ret)
        ret_co = log(open_prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_co = ret_co - mean(ret_co)
        ret_oc = log(prices(2:nprices,icol) / open_prices(2:nprices,icol))
        ret_oc = ret_oc - mean(ret_oc)

        active_cc = ret
        active_co = ret_co
        active_oc = ret_oc

        do imod = 1, n_model
            active_model = uppercase(trim(models(imod)))
            call fit_split_model(active_model, fopt, params, co_frac, persist, niter, converged, np)
            logl_fit = -real(nobs, dp) * fopt
            logl_cc = cc_loglik(ret, params, persist, active_model)
            aic = 2.0_dp * real(np, dp) - 2.0_dp * logl_fit
            bic = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl_fit
            call daily_vol_forecast(ret, params, persist, active_model, vol_forecast)
            vol_ann = sqrt(sum(vol_forecast**2) / real(nobs, dp))
            call cc_skew_kurt(ret, params, persist, active_model, skew, ekurt)
            vol_forecasts(:,icol,imod) = vol_forecast
            fit_model_names(imod) = active_model
            fit_model_have(imod) = .true.

            write(*,'(A16,1X,A9,ES12.3,7F8.4,F10.2,4F12.2,I6,1X,L1,2F9.3)') &
                trim(active_model), trim(col_names(icol)), params%omega, params%alpha, params%gamma, &
                params%beta, params%theta, params%twist, co_frac, persist, vol_ann, logl_fit, logl_cc, &
                aic, bic, niter, converged, skew, ekurt
        end do
    end do

    call print_implied_vol_correlations(implied_vol_file, col_names, fit_model_names, fit_model_have, &
                                        dates(2:nprices), vol_forecasts, iv_dates, iv_col_names, iv_values, &
                                        asset_iv_assets, asset_iv_indices)
    call print_implied_vol_correlations(implied_vol_file, col_names, fit_model_names, fit_model_have, &
                                        dates(2:nprices), vol_forecasts, iv_dates, iv_col_names, iv_values, &
                                        asset_iv_assets, asset_iv_indices, "log_diff", [1, 5, 22], &
                                        extra_series_names, dates, extra_series, extra_signs, &
                                        standardize_log_return_by_prior_iv)

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    subroutine fit_split_model(model_name, f_best, params_best, co_frac_best, persist_best, &
                               niter_best, converged_best, np_out)
        character(len=*), intent(in) :: model_name
        real(dp), intent(out) :: f_best, co_frac_best, persist_best
        type(garch_params_t), intent(out) :: params_best
        integer, intent(out) :: niter_best, np_out
        logical, intent(out) :: converged_best
        integer, parameter :: n_start = 4, max_np = 6
        real(dp) :: p(max_np), p0(max_np), p_best(max_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        select case (trim(model_name))
        case ("SPLIT_SYMM")
            np_out = 4
        case ("SPLIT_NAGARCH", "SPLIT_GJR")
            np_out = 5
        case ("SPLIT_FGTWIST")
            np_out = 6
        case default
            print '(A,A)', "Unknown split model: ", trim(model_name)
            error stop
        end select

        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.
        do istart = 1, n_start
            call split_start(model_name, istart, p0(1:np_out))
            p(1:np_out) = p0(1:np_out)
            call bfgs_minimize(split_obj, p, np_out, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best(1:np_out) = p(1:np_out)
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        call split_transform(model_name, p_best(1:np_out), params_best, co_frac_best, persist_best)
    end subroutine fit_split_model

    subroutine split_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: step, fp, fm
        real(dp) :: pp(np), pm(np)
        integer :: i

        f = split_value(p, np)
        do i = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(i)))
            pp = p
            pm = p
            pp(i) = pp(i) + step
            pm(i) = pm(i) - step
            fp = split_value(pp, np)
            fm = split_value(pm, np)
            g(i) = (fp - fm) / (2.0_dp*step)
        end do
    end subroutine split_obj

    real(dp) function split_value(p, np)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        type(garch_params_t) :: params
        real(dp) :: co_frac, persist, h, h_co, h_oc
        integer :: t

        call split_transform(active_model, p, params, co_frac, persist)
        h = initial_variance(active_cc, params, persist)
        split_value = 0.0_dp
        do t = 1, size(active_cc)
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                split_value = huge(1.0_dp)
                return
            end if
            h = max(h, 1.0e-12_dp)
            h_co = max(co_frac*h, 1.0e-12_dp)
            h_oc = max((1.0_dp - co_frac)*h, 1.0e-12_dp)
            split_value = split_value + log_sqrt_2pi + 0.5_dp*(log(h_co) + active_co(t)**2/h_co) + &
                          log_sqrt_2pi + 0.5_dp*(log(h_oc) + active_oc(t)**2/h_oc)
            h = next_variance(active_cc(t), h, params, active_model)
        end do
        split_value = split_value / real(size(active_cc), dp)
    end function split_value

    subroutine split_start(model_name, istart, p)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: istart
        real(dp), intent(out) :: p(:)
        real(dp), parameter :: omega0(4) = [1.0e-6_dp, 1.0e-6_dp, 5.0e-6_dp, 5.0e-6_dp]
        real(dp), parameter :: a0(4) = [0.04_dp, 0.06_dp, 0.08_dp, 0.10_dp]
        real(dp), parameter :: b0(4) = [0.90_dp, 0.88_dp, 0.85_dp, 0.80_dp]
        real(dp), parameter :: th0(4) = [-0.50_dp, 0.00_dp, 0.50_dp, 1.00_dp]
        real(dp), parameter :: tw0(4) = [0.00_dp, 0.00_dp, -0.20_dp, 0.20_dp]
        type(garch_params_t) :: params
        real(dp) :: co_frac

        params = garch_params_t()
        params%omega = omega0(istart)
        params%alpha = a0(istart)
        params%beta = b0(istart)
        co_frac = 0.25_dp
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            call split_inv_transform(model_name, params, co_frac, p)
        case ("SPLIT_NAGARCH")
            params%theta = th0(istart)
            call split_inv_transform(model_name, params, co_frac, p)
        case ("SPLIT_GJR")
            params%alpha = 0.03_dp + 0.01_dp*real(istart - 1, dp)
            params%gamma = 0.08_dp + 0.03_dp*real(istart - 1, dp)
            params%beta = 0.88_dp - 0.02_dp*real(istart - 1, dp)
            call split_inv_transform(model_name, params, co_frac, p)
        case ("SPLIT_FGTWIST")
            params%theta = th0(istart)
            params%twist = tw0(istart)
            call split_inv_transform(model_name, params, co_frac, p)
        end select
    end subroutine split_start

    subroutine split_transform(model_name, p, params, co_frac, persist)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: p(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: co_frac, persist
        real(dp) :: e2, e3, e4, s, aa, moment

        params = garch_params_t()
        params%omega = exp(p(1))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            call simplex2(p(2), p(3), params%alpha, params%beta)
            co_frac = logistic(p(4))
            persist = params%alpha + params%beta
        case ("SPLIT_NAGARCH")
            params%theta = p(4)
            moment = 1.0_dp + params%theta**2
            call simplex2(p(2), p(3), aa, params%beta)
            params%alpha = aa / moment
            co_frac = logistic(p(5))
            persist = aa + params%beta
        case ("SPLIT_GJR")
            e2 = exp(p(2))
            e3 = exp(p(3))
            e4 = exp(p(4))
            s = 1.0_dp + e2 + e3 + e4
            params%alpha = e2 / s
            params%gamma = 2.0_dp * e3 / s
            params%beta = e4 / s
            co_frac = logistic(p(5))
            persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case ("SPLIT_FGTWIST")
            params%theta = p(4)
            params%twist = p(5)
            moment = fgarch_twist_moment_local(params%theta, params%twist)
            call simplex2(p(2), p(3), aa, params%beta)
            params%alpha = aa / moment
            co_frac = logistic(p(6))
            persist = aa + params%beta
        end select
    end subroutine split_transform

    subroutine split_inv_transform(model_name, params, co_frac, p)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: co_frac
        real(dp), intent(out) :: p(:)
        real(dp) :: aa, bg, slack, moment

        p(1) = log(max(params%omega, 1.0e-12_dp))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            slack = max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp)
            p(2) = log(max(params%alpha, 1.0e-12_dp) / slack)
            p(3) = log(max(params%beta, 1.0e-12_dp) / slack)
            p(4) = logit(co_frac)
        case ("SPLIT_NAGARCH")
            moment = 1.0_dp + params%theta**2
            aa = params%alpha * moment
            slack = max(1.0_dp - aa - params%beta, 1.0e-8_dp)
            p(2) = log(max(aa, 1.0e-12_dp) / slack)
            p(3) = log(max(params%beta, 1.0e-12_dp) / slack)
            p(4) = params%theta
            p(5) = logit(co_frac)
        case ("SPLIT_GJR")
            bg = 0.5_dp * params%gamma
            slack = max(1.0_dp - params%alpha - bg - params%beta, 1.0e-8_dp)
            p(2) = log(max(params%alpha, 1.0e-12_dp) / slack)
            p(3) = log(max(bg, 1.0e-12_dp) / slack)
            p(4) = log(max(params%beta, 1.0e-12_dp) / slack)
            p(5) = logit(co_frac)
        case ("SPLIT_FGTWIST")
            moment = fgarch_twist_moment_local(params%theta, params%twist)
            aa = params%alpha * moment
            slack = max(1.0_dp - aa - params%beta, 1.0e-8_dp)
            p(2) = log(max(aa, 1.0e-12_dp) / slack)
            p(3) = log(max(params%beta, 1.0e-12_dp) / slack)
            p(4) = params%theta
            p(5) = params%twist
            p(6) = logit(co_frac)
        end select
    end subroutine split_inv_transform

    subroutine simplex2(pa, pb, a, b)
        real(dp), intent(in) :: pa, pb
        real(dp), intent(out) :: a, b
        real(dp) :: ea, eb, s

        ea = exp(pa)
        eb = exp(pb)
        s = 1.0_dp + ea + eb
        a = ea / s
        b = eb / s
    end subroutine simplex2

    real(dp) function logistic(x)
        real(dp), intent(in) :: x
        logistic = 1.0_dp / (1.0_dp + exp(-x))
        logistic = min(max(logistic, 1.0e-5_dp), 1.0_dp - 1.0e-5_dp)
    end function logistic

    real(dp) function logit(x)
        real(dp), intent(in) :: x
        real(dp) :: xx

        xx = min(max(x, 1.0e-5_dp), 1.0_dp - 1.0e-5_dp)
        logit = log(xx / (1.0_dp - xx))
    end function logit

    real(dp) function initial_variance(y, params, persist)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params

        initial_variance = max(params%omega / max(1.0_dp - persist, 1.0e-8_dp), &
            sum(y**2) / real(size(y), dp), 1.0e-12_dp)
    end function initial_variance

    real(dp) function next_variance(y, h, params, model_name)
        real(dp), intent(in) :: y, h
        type(garch_params_t), intent(in) :: params
        character(len=*), intent(in) :: model_name
        real(dp) :: sqrth, r, ind, z, q

        sqrth = sqrt(max(h, 1.0e-12_dp))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            next_variance = params%omega + params%alpha*y**2 + params%beta*h
        case ("SPLIT_NAGARCH")
            r = y - params%theta*sqrth
            next_variance = params%omega + params%alpha*r**2 + params%beta*h
        case ("SPLIT_GJR")
            ind = merge(1.0_dp, 0.0_dp, y < 0.0_dp)
            next_variance = params%omega + (params%alpha + params%gamma*ind)*y**2 + params%beta*h
        case ("SPLIT_FGTWIST")
            z = y / sqrth
            q = abs(z - params%theta) - params%twist*(z - params%theta)
            next_variance = params%omega + params%alpha*h*q**2 + params%beta*h
        case default
            next_variance = h
        end select
        next_variance = max(next_variance, 1.0e-12_dp)
    end function next_variance

    real(dp) function cc_loglik(y, params, persist, model_name)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params
        character(len=*), intent(in) :: model_name
        real(dp) :: h
        integer :: t

        h = initial_variance(y, params, persist)
        cc_loglik = 0.0_dp
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            cc_loglik = cc_loglik - (log_sqrt_2pi + 0.5_dp*(log(h) + y(t)**2/h))
            h = next_variance(y(t), h, params, model_name)
        end do
    end function cc_loglik

    subroutine daily_vol_forecast(y, params, persist, model_name, vol)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params
        character(len=*), intent(in) :: model_name
        real(dp), intent(out) :: vol(:)
        real(dp) :: h
        integer :: t

        h = initial_variance(y, params, persist)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            vol(t) = sqrt(trading_days*h) * 100.0_dp
            h = next_variance(y(t), h, params, model_name)
        end do
    end subroutine daily_vol_forecast

    subroutine cc_skew_kurt(y, params, persist, model_name, skew, ekurt)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params
        character(len=*), intent(in) :: model_name
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: z(:)
        real(dp) :: h
        integer :: t

        allocate(z(size(y)))
        h = initial_variance(y, params, persist)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            z(t) = y(t) / sqrt(h)
            h = next_variance(y(t), h, params, model_name)
        end do
        call moments(z, skew, ekurt)
        deallocate(z)
    end subroutine cc_skew_kurt

    subroutine moments(z, skew, ekurt)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: avg, m2, m3, m4, dz, rn
        integer :: t

        rn = real(size(z), dp)
        avg = sum(z) / rn
        m2 = 0.0_dp
        m3 = 0.0_dp
        m4 = 0.0_dp
        do t = 1, size(z)
            dz = z(t) - avg
            m2 = m2 + dz**2
            m3 = m3 + dz**3
            m4 = m4 + dz**4
        end do
        m2 = m2 / rn
        if (m2 <= 0.0_dp) then
            skew = 0.0_dp
            ekurt = 0.0_dp
        else
            skew = (m3 / rn) / m2**1.5_dp
            ekurt = (m4 / rn) / m2**2 - 3.0_dp
        end if
    end subroutine moments

    real(dp) function fgarch_twist_moment_local(theta, twist)
        real(dp), intent(in) :: theta, twist
        real(dp) :: Phi, ph, A, B, one_minus, one_plus

        Phi = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
        ph  = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
        A = (1.0_dp + theta**2) * (1.0_dp - Phi) - theta * ph
        B = (1.0_dp + theta**2) * Phi + theta * ph
        one_minus = 1.0_dp - twist
        one_plus  = 1.0_dp + twist
        fgarch_twist_moment_local = max(one_minus**2 * A + one_plus**2 * B, 1.0e-8_dp)
    end function fgarch_twist_moment_local

    subroutine print_fit_assets()
        integer :: i

        if (fit_all_assets) then
            write(*,'(A)') "Fitting split GARCH models for assets: all price-file assets"
        else
            write(*,'(A)', advance='no') "Fitting split GARCH models for assets:"
            do i = 1, size(fit_assets)
                write(*,'(1X,A)', advance='no') trim(fit_assets(i))
            end do
            write(*,*)
        end if
    end subroutine print_fit_assets

end program xfit_split_ohlc_iv_returns

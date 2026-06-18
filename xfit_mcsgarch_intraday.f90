! Fit rugarch-style multiplicative component sGARCH to intraday OHLCV prices.
!
! The fitted returns are regular-session close-to-close intraday log returns.
! Overnight returns are excluded from the intraday return vector, but enter the
! exogenous daily variance proxy used for that trading day:
!   DailyVar_t = realized_var_{t-1, intraday} + overnight_return_t^2

module fit_mcsgarch_intraday_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use date_mod, only: date_label, yyyymmdd, seconds_per_minute, seconds_per_hour
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, intraday_bin_ids, &
                               default_bar_minutes, default_session_start_seconds
    use garch_mcsgarch_mod, only: mcsgarch_params_t, mcsgarch_fit_result_t, mcsgarch_remaining_var, &
                                  fit_mcsgarch, fit_mcsgarch_nagarch, fit_mcsgarch_gjr, &
                                  fit_mcsgarch_t, fit_mcsgarch_nagarch_t, fit_mcsgarch_gjr_t, &
                                  fit_mcsgarch_fs_skewt, fit_mcsgarch_nagarch_fs_skewt, fit_mcsgarch_gjr_fs_skewt
    use stats_mod, only: mean
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    use path_utils_mod, only: basename_without_extension
    implicit none
    private

    character(len=*), parameter :: file_pattern = "c:\python\intraday_prices\*.csv"
    real(dp), parameter :: min_daily_var = 1.0e-12_dp
    real(dp), parameter :: trading_days_per_year = 252.0_dp
    logical, parameter :: print_diurnal_curve = .false.
    logical, parameter :: write_diurnal_curve_csv = .true.
    logical, parameter :: write_forecast_history_csv = .true.
    logical, parameter :: smooth_diurnal_curve = .true.
    logical, parameter :: print_squared_noise_acf = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    integer, parameter :: squared_noise_acf_lags(7) = [1, 2, 3, 6, 12, 39, 78]
    character(len=12), parameter :: model_names(*) = [character(len=12) :: &
        "MCSGARCH", "MCSNAGARCH", "MCSGJRGARCH"]
    character(len=8), parameter :: dist_names(*) = [character(len=8) :: &
        "NORMAL", "T", "FS_SKEWT"]

    type :: comparison_row_t
        character(len=36) :: model = ""
        character(len=8) :: dist = ""
        integer :: k = 0
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: nu = 0.0_dp
        real(dp) :: xi = 0.0_dp
        real(dp) :: fit_sec = 0.0_dp
        real(dp) :: skew = 0.0_dp
        real(dp) :: ex_kurt = 0.0_dp
        integer :: aic_rank = 0
        integer :: bic_rank = 0
    end type comparison_row_t

    public :: run_fit_mcsgarch_intraday

contains

    ! Fit all configured MCS-GARCH models to each input file.
    subroutine run_fit_mcsgarch_intraday()
        character(len=MAX_PATH_LEN), allocatable :: filenames(:)
        integer :: i

        call collect_input_filenames(filenames, &
            file_pattern=file_pattern, &
            default_filenames=[character(len=MAX_PATH_LEN) :: &
                "c:\python\intraday_prices\spy_5min_databento.csv"])
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call fit_one_file(trim(filenames(i)))
        end do
        deallocate(filenames)
    end subroutine run_fit_mcsgarch_intraday

    ! Read intraday prices, fit the MCS-GARCH model, and write diagnostic outputs.
    subroutine fit_one_file(filename)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        real(dp), allocatable :: returns(:), daily_var(:), diurnal_var(:), q(:)
        real(dp), allocatable :: diurnal_fit(:,:), q_fit(:,:)
        real(dp), allocatable :: fit_sec(:)
        integer, allocatable :: bin_id(:), return_dates(:)
        type(mcsgarch_fit_result_t), allocatable :: fit(:)
        type(mcsgarch_fit_result_t) :: best_fit
        integer :: max_iter, nobs, nfit, imodel, idist, ifit
        real(dp) :: gtol, t_start, t0, t1, t_end, read_sec, elapsed_sec
        character(len=512) :: stem

        call cpu_time(t_start)
        stem = basename_without_extension(filename)
        max_iter = 120
        gtol = 1.0e-5_dp

        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), bars)
        call cpu_time(t1)
        read_sec = t1 - t0
        call filter_intraday_session(bars, regular_bars)
        call build_mcsgarch_inputs(regular_bars, returns, daily_var, bin_id, return_dates)
        nobs = size(returns)
        nfit = size(model_names) * size(dist_names)
        allocate(diurnal_var(nobs), q(nobs), diurnal_fit(nobs, nfit), q_fit(nobs, nfit), fit(nfit), fit_sec(nfit))

        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                call cpu_time(t0)
                call fit_configured_model(trim(model_names(imodel)), trim(dist_names(idist)), &
                                          returns, daily_var, bin_id, max_iter, gtol, &
                                          fit(ifit), diurnal_fit(:, ifit), q_fit(:, ifit))
                call cpu_time(t1)
                fit_sec(ifit) = t1 - t0
            end do
        end do
        call select_best_fit_path(fit, diurnal_fit, q_fit, diurnal_var, q, best_fit)
        call cpu_time(t_end)
        elapsed_sec = t_end - t_start
        call print_fit_summary(trim(filename), regular_bars%nobs(), returns, daily_var, bin_id, return_dates, &
                               fit, fit_sec, read_sec, elapsed_sec)
        call print_model_comparison(returns, daily_var, bin_id, fit, fit_sec, diurnal_fit, q_fit)
        if (print_squared_noise_acf) call print_squared_noise_acf_table(returns, daily_var, bin_id, diurnal_fit, q_fit)
        if (print_diurnal_curve) call print_diurnal_variance_curve(bin_id, diurnal_var, &
                                                                    smooth_diurnal_curve, &
                                                                    diurnal_smooth_half_width)
        if (write_diurnal_curve_csv) call write_diurnal_variance_curve_csv(trim(stem) // "_mcsgarch_diurnal.csv", &
                                                                           bin_id, diurnal_var)
        if (write_forecast_history_csv) call write_forecast_history_csv_file(trim(stem) // "_mcsgarch_forecasts.csv", &
                                                                             returns, daily_var, bin_id, &
                                                                             return_dates, diurnal_var, q, &
                                                                             best_fit%params, best_fit%persist)

        deallocate(returns, daily_var, bin_id, return_dates, diurnal_var, q, diurnal_fit, q_fit, fit, fit_sec)
    end subroutine fit_one_file

    ! Fit one configured dynamic model and innovation distribution.
    subroutine fit_configured_model(model_name, dist_name, returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: returns(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: fit
        real(dp), intent(out) :: diurnal_var(:), q(:)

        select case (trim(model_name))
        case ("MCSGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                    smooth_diurnal=smooth_diurnal_curve, &
                                    smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                           smooth_diurnal=smooth_diurnal_curve, &
                                           smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                  smooth_diurnal=smooth_diurnal_curve, &
                                  smooth_half_width=diurnal_smooth_half_width)
            end if
        case ("MCSNAGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_nagarch_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                            smooth_diurnal=smooth_diurnal_curve, &
                                            smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_nagarch_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                                   smooth_diurnal=smooth_diurnal_curve, &
                                                   smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch_nagarch(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                          smooth_diurnal=smooth_diurnal_curve, &
                                          smooth_half_width=diurnal_smooth_half_width)
            end if
        case ("MCSGJRGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_gjr_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                        smooth_diurnal=smooth_diurnal_curve, &
                                        smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_gjr_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                               smooth_diurnal=smooth_diurnal_curve, &
                                               smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch_gjr(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                      smooth_diurnal=smooth_diurnal_curve, &
                                      smooth_half_width=diurnal_smooth_half_width)
            end if
        case default
            error stop "fit_configured_model: unsupported model"
        end select
    end subroutine fit_configured_model

    ! Convert regular-session bars into within-session returns and daily variance proxies.
    subroutine build_mcsgarch_inputs(bars, returns, daily_var, bin_id, return_dates)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: returns(:), daily_var(:)
        integer, allocatable, intent(out) :: bin_id(:), return_dates(:)
        integer, allocatable :: bar_bins(:), day_index(:), day_dates(:), day_first(:), day_last(:)
        real(dp), allocatable :: intraday_rv(:), overnight_sq(:), daily_proxy(:)
        real(dp), allocatable :: returns_all(:), daily_var_all(:)
        integer, allocatable :: bin_all(:), dates_all(:)
        integer :: n, ndays, i, k, d
        real(dp) :: fallback_daily_var

        n = bars%nobs()
        if (n < 3) error stop "build_mcsgarch_inputs: not enough regular-session bars"
        call intraday_bin_ids(bars, bar_bins)
        call map_days(bars, day_index, day_dates, day_first, day_last)
        ndays = size(day_dates)
        allocate(intraday_rv(ndays), overnight_sq(ndays), daily_proxy(ndays))
        intraday_rv = 0.0_dp
        overnight_sq = 0.0_dp

        do i = 2, n
            if (day_index(i) == day_index(i - 1)) then
                intraday_rv(day_index(i)) = intraday_rv(day_index(i)) + log(bars%close(i) / bars%close(i - 1))**2
            end if
        end do
        do d = 2, ndays
            overnight_sq(d) = log(bars%open(day_first(d)) / bars%close(day_last(d - 1)))**2
        end do

        fallback_daily_var = max(mean(intraday_rv(max(1, min(2, ndays)):ndays)), min_daily_var)
        daily_proxy(1) = fallback_daily_var
        do d = 2, ndays
            daily_proxy(d) = max(intraday_rv(d - 1) + overnight_sq(d), min_daily_var)
        end do

        allocate(returns_all(n - 1), daily_var_all(n - 1), bin_all(n - 1), dates_all(n - 1))
        k = 0
        do i = 2, n
            if (day_index(i) /= day_index(i - 1)) cycle
            k = k + 1
            returns_all(k) = log(bars%close(i) / bars%close(i - 1))
            daily_var_all(k) = daily_proxy(day_index(i))
            bin_all(k) = bar_bins(i)
            dates_all(k) = day_dates(day_index(i))
        end do
        if (k < 3) error stop "build_mcsgarch_inputs: not enough intraday returns"
        allocate(returns(k), daily_var(k), bin_id(k), return_dates(k))
        returns = returns_all(1:k)
        daily_var = daily_var_all(1:k)
        bin_id = bin_all(1:k)
        return_dates = dates_all(1:k)

        deallocate(bar_bins, day_index, day_dates, day_first, day_last, intraday_rv, overnight_sq, daily_proxy)
        deallocate(returns_all, daily_var_all, bin_all, dates_all)
    end subroutine build_mcsgarch_inputs

    ! Map each bar to a trading day and record the first and last bar of each day.
    subroutine map_days(bars, day_index, day_dates, day_first, day_last)
        type(ohlcv_series_t), intent(in) :: bars
        integer, allocatable, intent(out) :: day_index(:), day_dates(:), day_first(:), day_last(:)
        integer :: n, i, ndays, current_date

        n = bars%nobs()
        allocate(day_index(n), day_dates(n), day_first(n), day_last(n))
        ndays = 0
        current_date = -1
        do i = 1, n
            if (yyyymmdd(bars%timestamp(i)%date) /= current_date) then
                ndays = ndays + 1
                current_date = yyyymmdd(bars%timestamp(i)%date)
                day_dates(ndays) = current_date
                day_first(ndays) = i
                if (ndays > 1) day_last(ndays - 1) = i - 1
            end if
            day_index(i) = ndays
        end do
        day_last(ndays) = n
        day_dates = day_dates(1:ndays)
        day_first = day_first(1:ndays)
        day_last = day_last(1:ndays)
    end subroutine map_days

    ! Select the fitted variance path from the configured fit with the highest log likelihood.
    subroutine select_best_fit_path(fit, diurnal_fit, q_fit, diurnal_best, q_best, best_fit)
        type(mcsgarch_fit_result_t), intent(in) :: fit(:)
        real(dp), intent(in) :: diurnal_fit(:,:), q_fit(:,:)
        real(dp), intent(out) :: diurnal_best(:), q_best(:)
        type(mcsgarch_fit_result_t), intent(out) :: best_fit
        integer :: best_index

        best_index = maxloc(fit%loglik, dim=1)
        diurnal_best = diurnal_fit(:, best_index)
        q_best = q_fit(:, best_index)
        best_fit = fit(best_index)
    end subroutine select_best_fit_path

    ! Print the fitted MCS-GARCH parameter rows and data/timing summary.
    subroutine print_fit_summary(filename, n_regular_bars, returns, daily_var, bin_id, return_dates, fit, fit_sec, &
                                 read_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_regular_bars
        real(dp), intent(in) :: returns(:), daily_var(:)
        integer, intent(in) :: bin_id(:), return_dates(:)
        type(mcsgarch_fit_result_t), intent(in) :: fit(:)
        real(dp), intent(in) :: fit_sec(:)
        real(dp), intent(in) :: read_sec, elapsed_sec
        integer :: ifit, imodel, idist

        print '(A)', "MCS-GARCH intraday fits"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0)', "Regular-session bars: ", n_regular_bars
        print '(A,I0,A,A,A,A)', "Intraday returns used: ", size(returns), " from ", &
              date_label(return_dates(1)), " to ", date_label(return_dates(size(return_dates)))
        print '(A)', "Returns are within-session close-to-close log returns; overnight returns enter DailyVar only."
        print '(A,I0,A,I0)', "Intraday bin range: ", minval(bin_id), " to ", maxval(bin_id)
        print '(A,ES12.4,A,ES12.4)', "Mean DailyVar: ", mean(daily_var), "  mean return variance proxy: ", mean(returns**2)
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------"
        print '(A12,1X,A8,1X,A12,7(1X,A10),1X,A8,1X,A5,1X,A13,1X,A9)', &
              "Model", "Dist", "omega", "alpha", "gamma", "beta", "theta", "nu", "xi", "persist", "iter", "conv", &
              "logL", "fit_sec"
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------"
        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                call print_fit_row(trim(model_names(imodel)), trim(dist_names(idist)), fit(ifit), fit_sec(ifit))
            end do
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------"
        print '(A,F10.3)', "Data read seconds: ", read_sec
        print '(A,F10.3)', "Elapsed seconds:   ", elapsed_sec
    end subroutine print_fit_summary

    ! Print one fitted dynamic-model parameter row.
    subroutine print_fit_row(model_name, dist_name, fit, fit_sec)
        character(len=*), intent(in) :: model_name, dist_name
        type(mcsgarch_fit_result_t), intent(in) :: fit
        real(dp), intent(in) :: fit_sec

        print '(A12,1X,A8,1X,ES12.4,7(1X,F10.4),1X,I8,5X,L1,1X,F13.3,1X,F9.3)', &
              model_name, dist_name, fit%params%omega, fit%params%alpha, fit%params%gamma, &
              fit%params%beta, fit%params%theta, fit%nu, fit%xi, fit%persist, fit%niter, fit%converged, fit%loglik, &
              fit_sec
    end subroutine print_fit_row

    ! Compare nested intraday volatility models with IC values and residual moments.
    subroutine print_model_comparison(returns, daily_var, bin_id, fit, fit_sec, diurnal_fit, q_fit)
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_fit(:,:), q_fit(:,:)
        integer, intent(in) :: bin_id(:)
        type(mcsgarch_fit_result_t), intent(in) :: fit(:)
        real(dp), intent(in) :: fit_sec(:)
        real(dp), allocatable :: h(:), unit_scale(:)
        type(comparison_row_t), allocatable :: rows(:)
        real(dp) :: sigma2
        integer :: nobs, k, nbins_used, ifit, imodel, idist, irow, aic_best, bic_best

        nobs = size(returns)
        allocate(h(nobs), unit_scale(nobs), rows(3 + size(fit)))
        unit_scale = 1.0_dp

        sigma2 = max(sum(returns**2) / real(nobs, dp), min_daily_var)
        h = sigma2
        k = 1
        irow = 1
        call fill_model_comparison_row(rows(irow), "constant intraday volatility", "NORMAL", k, returns, h)

        call bin_scaled_variance_path(returns, unit_scale, bin_id, h, nbins_used)
        k = nbins_used
        irow = irow + 1
        call fill_model_comparison_row(rows(irow), "diurnal multiplier only", "NORMAL", k, returns, h)

        call bin_scaled_variance_path(returns, daily_var, bin_id, h, nbins_used)
        k = nbins_used
        irow = irow + 1
        call fill_model_comparison_row(rows(irow), "daily volatility x diurnal", "NORMAL", k, returns, h)

        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                irow = irow + 1
                h = max(daily_var * diurnal_fit(:, ifit) * q_fit(:, ifit), min_daily_var)
                k = count_populated_bins(bin_id) + model_param_count(trim(model_names(imodel))) + &
                    dist_param_count(trim(dist_names(idist)))
                call fill_model_comparison_row(rows(irow), model_table_label(trim(model_names(imodel))), &
                                               trim(dist_names(idist)), k, returns, h, &
                                               loglik_override=fit(ifit)%loglik, nu=fit(ifit)%nu, xi=fit(ifit)%xi, &
                                               fit_sec=fit_sec(ifit))
            end do
        end do

        call rank_model_comparison_rows(rows)
        aic_best = minloc(rows%aic, dim=1)
        bic_best = minloc(rows%bic, dim=1)

        print '(A)', ""
        print '(A)', "Intraday volatility model comparison"
        print '(A)', "------------------------------------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A8,1X,A4,1X,A14,1X,A14,1X,A14,1X,A7,1X,A7,1X,A9,1X,A8,1X,A8,1X,A8,1X,A9)', &
              "Model", "Dist", "k", "logL", "AIC", "BIC", "nu", "xi", "skew", "ex_kurt", "AIC_rank", "BIC_rank", &
              "fit_sec"
        print '(A)', "------------------------------------------------------------------------------------------------------------------------------------------"
        do irow = 1, size(rows)
            call print_model_comparison_row(rows(irow))
        end do
        print '(A)', "------------------------------------------------------------------------------------------------------------------------------------------"
        print '(A,A,A,A)', "AIC selects: ", trim(rows(aic_best)%model), " ", trim(rows(aic_best)%dist)
        print '(A,A,A,A)', "BIC selects: ", trim(rows(bic_best)%model), " ", trim(rows(bic_best)%dist)

        deallocate(h, unit_scale, rows)
    end subroutine print_model_comparison

    ! Assign AIC and BIC ranks, with rank 1 being lowest IC.
    subroutine rank_model_comparison_rows(rows)
        type(comparison_row_t), intent(inout) :: rows(:)
        integer :: i

        do i = 1, size(rows)
            rows(i)%aic_rank = 1 + count(rows%aic < rows(i)%aic)
            rows(i)%bic_rank = 1 + count(rows%bic < rows(i)%bic)
        end do
    end subroutine rank_model_comparison_rows

    ! Return the dynamic-parameter count for an MCS model name.
    integer function model_param_count(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_param_count = 3
        case ("MCSNAGARCH", "MCSGJRGARCH")
            model_param_count = 4
        case default
            error stop "model_param_count: unsupported model"
        end select
    end function model_param_count

    ! Return the innovation-distribution parameter count.
    integer function dist_param_count(dist_name)
        character(len=*), intent(in) :: dist_name

        select case (trim(dist_name))
        case ("NORMAL")
            dist_param_count = 0
        case ("T")
            dist_param_count = 1
        case ("FS_SKEWT")
            dist_param_count = 2
        case default
            error stop "dist_param_count: unsupported distribution"
        end select
    end function dist_param_count

    ! Return display label used in the comparison table.
    character(len=36) function model_table_label(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_table_label = "MCS-GARCH"
        case ("MCSNAGARCH")
            model_table_label = "MCS-NAGARCH"
        case ("MCSGJRGARCH")
            model_table_label = "MCS-GJRGARCH"
        case default
            model_table_label = trim(model_name)
        end select
    end function model_table_label

    ! Fill one model-comparison row from a variance path and standardized residual moments.
    subroutine fill_model_comparison_row(row, model_name, dist_name, k, returns, h, loglik_override, nu, xi, fit_sec)
        type(comparison_row_t), intent(out) :: row
        character(len=*), intent(in) :: model_name, dist_name
        integer, intent(in) :: k
        real(dp), intent(in) :: returns(:), h(:)
        real(dp), intent(in), optional :: loglik_override, nu, xi, fit_sec
        integer :: nobs

        nobs = size(returns)
        row%model = model_name
        row%dist = dist_name
        row%k = k
        if (present(loglik_override)) then
            row%loglik = loglik_override
        else
            row%loglik = gaussian_loglik(returns, h)
        end if
        call standardized_residual_moments(returns, h, row%skew, row%ex_kurt)
        row%nu = 0.0_dp
        if (present(nu)) row%nu = nu
        row%xi = 0.0_dp
        if (present(xi)) row%xi = xi
        row%fit_sec = 0.0_dp
        if (present(fit_sec)) row%fit_sec = fit_sec

        row%aic = -2.0_dp*row%loglik + 2.0_dp*real(k, dp)
        row%bic = -2.0_dp*row%loglik + real(k, dp)*log(real(nobs, dp))
    end subroutine fill_model_comparison_row

    ! Print one ranked model-comparison row.
    subroutine print_model_comparison_row(row)
        type(comparison_row_t), intent(in) :: row

        print '(A32,1X,A8,1X,I4,1X,F14.3,1X,F14.3,1X,F14.3,1X,F7.3,1X,F7.3,1X,F9.4,1X,F8.4,1X,I8,1X,I8,1X,F9.3)', &
              row%model, row%dist, row%k, row%loglik, row%aic, row%bic, row%nu, row%xi, row%skew, row%ex_kurt, &
              row%aic_rank, row%bic_rank, row%fit_sec
    end subroutine print_model_comparison_row

    ! Print autocorrelations of squared standardized intraday residuals for baselines and fitted models.
    subroutine print_squared_noise_acf_table(returns, daily_var, bin_id, diurnal_fit, q_fit)
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_fit(:,:), q_fit(:,:)
        integer, intent(in) :: bin_id(:)
        real(dp), allocatable :: h(:), unit_scale(:)
        real(dp) :: sigma2
        integer :: ifit, imodel, idist, nbins_used, nobs

        nobs = size(returns)
        allocate(h(nobs), unit_scale(nobs))
        unit_scale = 1.0_dp
        print '(A)', ""
        print '(A)', "Squared standardized intraday noise autocorrelations"
        print '(A)', "------------------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A8,7(1X,A8))', "Model", "Dist", "acf_1", "acf_2", "acf_3", "acf_6", "acf_12", &
              "acf_39", "acf_78"
        print '(A)', "------------------------------------------------------------------------------------------------------------------------"

        sigma2 = max(sum(returns**2) / real(nobs, dp), min_daily_var)
        h = sigma2
        call print_squared_noise_acf_row("constant intraday volatility", "NORMAL", returns, h)

        call bin_scaled_variance_path(returns, unit_scale, bin_id, h, nbins_used)
        call print_squared_noise_acf_row("diurnal multiplier only", "NORMAL", returns, h)

        call bin_scaled_variance_path(returns, daily_var, bin_id, h, nbins_used)
        call print_squared_noise_acf_row("daily volatility x diurnal", "NORMAL", returns, h)

        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                h = max(daily_var * diurnal_fit(:, ifit) * q_fit(:, ifit), min_daily_var)
                call print_squared_noise_acf_row(model_table_label(trim(model_names(imodel))), trim(dist_names(idist)), &
                                                 returns, h)
            end do
        end do
        print '(A)', "------------------------------------------------------------------------------------------------------------------------"
        deallocate(h, unit_scale)
    end subroutine print_squared_noise_acf_table

    ! Print one squared standardized residual autocorrelation row.
    subroutine print_squared_noise_acf_row(model_name, dist_name, returns, h)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: returns(:), h(:)
        real(dp) :: acf(size(squared_noise_acf_lags))

        call squared_standardized_residual_acf(returns, h, squared_noise_acf_lags, acf)
        print '(A32,1X,A8,7(1X,F8.4))', model_name, dist_name, acf
    end subroutine print_squared_noise_acf_row

    ! Compute autocorrelations of squared standardized residuals at requested lags.
    subroutine squared_standardized_residual_acf(returns, h, lags, acf)
        real(dp), intent(in) :: returns(:), h(:)
        integer, intent(in) :: lags(:)
        real(dp), intent(out) :: acf(:)
        real(dp), allocatable :: x(:)
        real(dp) :: xbar, denom, numer
        integer :: i, ilag, lag, n

        if (size(returns) /= size(h)) error stop "squared_standardized_residual_acf: array sizes differ"
        if (size(lags) /= size(acf)) error stop "squared_standardized_residual_acf: lag/acf sizes differ"
        n = size(returns)
        allocate(x(n))
        do i = 1, n
            x(i) = returns(i)**2 / max(h(i), min_daily_var)
        end do
        xbar = sum(x) / real(n, dp)
        denom = sum((x - xbar)**2)

        do ilag = 1, size(lags)
            lag = lags(ilag)
            if (lag < 1 .or. lag >= n .or. denom <= min_daily_var) then
                acf(ilag) = 0.0_dp
            else
                numer = 0.0_dp
                do i = lag + 1, n
                    numer = numer + (x(i) - xbar) * (x(i - lag) - xbar)
                end do
                acf(ilag) = numer / denom
            end if
        end do
        deallocate(x)
    end subroutine squared_standardized_residual_acf

    ! Estimate a per-bin variance multiplier around a supplied variance scale.
    subroutine bin_scaled_variance_path(returns, scale, bin_id, h, nbins_used)
        real(dp), intent(in) :: returns(:), scale(:)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(out) :: h(:)
        integer, intent(out) :: nbins_used
        real(dp), allocatable :: bin_sum(:)
        integer, allocatable :: bin_count(:)
        integer :: b, i, nbins
        real(dp) :: fallback

        if (size(returns) /= size(scale) .or. size(returns) /= size(bin_id) .or. &
            size(returns) /= size(h)) then
            error stop "bin_scaled_variance_path: array sizes differ"
        end if
        nbins = maxval(bin_id)
        allocate(bin_sum(nbins), bin_count(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do i = 1, size(returns)
            b = bin_id(i)
            bin_sum(b) = bin_sum(b) + returns(i)**2 / max(scale(i), min_daily_var)
            bin_count(b) = bin_count(b) + 1
        end do

        nbins_used = count(bin_count > 0)
        fallback = max(sum(returns**2) / real(size(returns), dp), min_daily_var)
        do b = 1, nbins
            if (bin_count(b) > 0) then
                bin_sum(b) = max(bin_sum(b) / real(bin_count(b), dp), min_daily_var)
            else
                bin_sum(b) = fallback
            end if
        end do
        do i = 1, size(returns)
            h(i) = max(scale(i) * bin_sum(bin_id(i)), min_daily_var)
        end do
        deallocate(bin_sum, bin_count)
    end subroutine bin_scaled_variance_path

    ! Return the zero-mean Gaussian log likelihood for a variance path.
    real(dp) function gaussian_loglik(returns, h)
        real(dp), intent(in) :: returns(:), h(:)
        integer :: i

        if (size(returns) /= size(h)) error stop "gaussian_loglik: array sizes differ"
        gaussian_loglik = 0.0_dp
        do i = 1, size(returns)
            gaussian_loglik = gaussian_loglik - log_sqrt_2pi - 0.5_dp*log(max(h(i), min_daily_var)) - &
                              0.5_dp*returns(i)**2 / max(h(i), min_daily_var)
        end do
    end function gaussian_loglik

    ! Compute skewness and excess kurtosis of standardized residuals.
    subroutine standardized_residual_moments(returns, h, skew, ex_kurt)
        real(dp), intent(in) :: returns(:), h(:)
        real(dp), intent(out) :: skew, ex_kurt
        real(dp) :: z, m1, m2, m3, m4, centered
        integer :: i, n

        if (size(returns) /= size(h)) error stop "standardized_residual_moments: array sizes differ"
        n = size(returns)
        m1 = 0.0_dp
        do i = 1, n
            m1 = m1 + returns(i) / sqrt(max(h(i), min_daily_var))
        end do
        m1 = m1 / real(n, dp)

        m2 = 0.0_dp
        m3 = 0.0_dp
        m4 = 0.0_dp
        do i = 1, n
            z = returns(i) / sqrt(max(h(i), min_daily_var))
            centered = z - m1
            m2 = m2 + centered**2
            m3 = m3 + centered**3
            m4 = m4 + centered**4
        end do
        m2 = max(m2 / real(n, dp), min_daily_var)
        m3 = m3 / real(n, dp)
        m4 = m4 / real(n, dp)
        skew = m3 / m2**1.5_dp
        ex_kurt = m4 / m2**2 - 3.0_dp
    end subroutine standardized_residual_moments

    ! Count the number of intraday bins represented in the return sample.
    integer function count_populated_bins(bin_id)
        integer, intent(in) :: bin_id(:)
        integer, allocatable :: bin_count(:)
        integer :: i

        allocate(bin_count(maxval(bin_id)))
        bin_count = 0
        do i = 1, size(bin_id)
            bin_count(bin_id(i)) = bin_count(bin_id(i)) + 1
        end do
        count_populated_bins = count(bin_count > 0)
        deallocate(bin_count)
    end function count_populated_bins

    ! Print the estimated intraday seasonal variance curve by time-of-day bin.
    subroutine print_diurnal_variance_curve(bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(in) :: diurnal_var(:)
        logical, intent(in) :: smooth_diurnal
        integer, intent(in) :: smooth_half_width
        real(dp), allocatable :: bin_sum(:)
        integer, allocatable :: bin_count(:)
        integer :: b, i, nbins

        nbins = maxval(bin_id)
        allocate(bin_sum(nbins), bin_count(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do i = 1, size(bin_id)
            b = bin_id(i)
            bin_sum(b) = bin_sum(b) + diurnal_var(i)
            bin_count(b) = bin_count(b) + 1
        end do

        print '(A)', ""
        if (smooth_diurnal) then
            print '(A,I0,A)', "Estimated intraday diurnal volatility curve (log-variance triangular smooth, half-width ", &
                  smooth_half_width, ")"
        else
            print '(A)', "Estimated intraday diurnal volatility curve (raw)"
        end if
        print '(A)', "--------------------------------------------"
        print '(A6,1X,A8,1X,A8,1X,A10,1X,A10)', "bin", "time", "count", "var_mult", "vol_mult"
        print '(A)', "--------------------------------------------"
        do b = 1, nbins
            if (bin_count(b) < 1) cycle
            bin_sum(b) = bin_sum(b) / real(bin_count(b), dp)
            print '(I6,1X,A8,1X,I8,1X,F10.4,1X,F10.4)', b, &
                  time_of_day_label(bin_end_seconds(b)), bin_count(b), bin_sum(b), sqrt(bin_sum(b))
        end do
        print '(A)', "--------------------------------------------"

        deallocate(bin_sum, bin_count)
    end subroutine print_diurnal_variance_curve

    ! Write the intraday seasonal variance curve to a CSV file.
    subroutine write_diurnal_variance_curve_csv(filename, bin_id, diurnal_var)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: bin_id(:)
        real(dp), intent(in) :: diurnal_var(:)
        real(dp), allocatable :: bin_sum(:)
        integer, allocatable :: bin_count(:)
        integer :: b, i, nbins, unit

        nbins = maxval(bin_id)
        allocate(bin_sum(nbins), bin_count(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do i = 1, size(bin_id)
            b = bin_id(i)
            bin_sum(b) = bin_sum(b) + diurnal_var(i)
            bin_count(b) = bin_count(b) + 1
        end do

        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit,'(A)') "bin,time,count,var_mult,vol_mult"
        do b = 1, nbins
            if (bin_count(b) < 1) cycle
            bin_sum(b) = bin_sum(b) / real(bin_count(b), dp)
            write(unit,'(I0,A,A,A,I0,A,F12.6,A,F12.6)') b, ",", time_of_day_label(bin_end_seconds(b)), ",", &
                bin_count(b), ",", bin_sum(b), ",", sqrt(bin_sum(b))
        end do
        close(unit)
        print '(A,A)', "Diurnal curve CSV: ", trim(filename)

        deallocate(bin_sum, bin_count)
    end subroutine write_diurnal_variance_curve_csv

    ! Write intraday return, MCS-GARCH forecasts, and remaining-session variance forecasts to CSV.
    !
    ! Columns added beyond the per-bar fit: remaining_var_garch (GARCH multi-step
    ! remaining-session variance), remaining_var_diurnal (diurnal-only baseline),
    ! and realized_remaining_var (sum of squared returns from the next bar to session end).
    subroutine write_forecast_history_csv_file(filename, returns, daily_var, bin_id, return_dates, &
                                               diurnal_var, q, params, persist)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), q(:), persist
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: bin_id(:), return_dates(:)
        integer :: unit, i, b, nobs, intraday_returns_per_day, nbins
        real(dp) :: intraday_var, daily_vol_ann_pct, intraday_vol_ann_pct
        real(dp) :: rem_garch, rem_diurnal, realized_rem
        real(dp), allocatable :: bin_diurnal(:), bin_tail_sum(:), realized_remaining(:)

        nobs = size(returns)
        if (size(daily_var) /= nobs .or. size(diurnal_var) /= nobs .or. size(q) /= nobs .or. &
            size(bin_id) /= nobs .or. size(return_dates) /= nobs) then
            error stop "write_forecast_history_csv_file: array sizes differ"
        end if

        intraday_returns_per_day = maxval(bin_id) - minval(bin_id) + 1
        nbins = maxval(bin_id)

        ! Build diurnal curve and cumulative tail sums indexed by bin number.
        allocate(bin_diurnal(nbins), bin_tail_sum(nbins))
        bin_diurnal = 0.0_dp
        do i = 1, nobs
            bin_diurnal(bin_id(i)) = diurnal_var(i)
        end do
        bin_tail_sum(nbins) = 0.0_dp
        do b = nbins - 1, 1, -1
            bin_tail_sum(b) = bin_tail_sum(b + 1) + bin_diurnal(b + 1)
        end do

        ! Cumulative sum of squared returns from each bar to the end of its session day.
        allocate(realized_remaining(nobs))
        realized_remaining(nobs) = 0.0_dp
        do i = nobs - 1, 1, -1
            if (return_dates(i + 1) == return_dates(i)) then
                realized_remaining(i) = realized_remaining(i + 1) + returns(i + 1)**2
            else
                realized_remaining(i) = 0.0_dp
            end if
        end do

        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit,'(A)') "date,time,bin,return_pct,daily_var,daily_vol_ann_pct,diurnal_var,q," // &
            "intraday_var,intraday_vol_ann_pct,remaining_var_garch,remaining_var_diurnal,realized_remaining_var"
        do i = 1, nobs
            b = bin_id(i)
            intraday_var = max(daily_var(i) * diurnal_var(i) * q(i), min_daily_var)
            daily_vol_ann_pct = 100.0_dp * sqrt(max(daily_var(i), min_daily_var) * trading_days_per_year)
            intraday_vol_ann_pct = 100.0_dp * sqrt(intraday_var * real(intraday_returns_per_day, dp) * &
                                                   trading_days_per_year)
            rem_garch = mcsgarch_remaining_var(params, persist, q(i), daily_var(i), bin_diurnal(b + 1:nbins))
            rem_diurnal = max(daily_var(i) * bin_tail_sum(b), 0.0_dp)
            realized_rem = realized_remaining(i)
            write(unit,'(A,",",A,",",I0,2(",",ES16.8),",",F12.6,2(",",ES16.8),",",ES16.8,",",F12.6,3(",",ES16.8))') &
                date_label(return_dates(i)), time_of_day_label(bin_end_seconds(bin_id(i))), b, &
                100.0_dp*returns(i), daily_var(i), daily_vol_ann_pct, diurnal_var(i), q(i), &
                intraday_var, intraday_vol_ann_pct, rem_garch, rem_diurnal, realized_rem
        end do
        close(unit)
        deallocate(bin_diurnal, bin_tail_sum, realized_remaining)
        print '(A,A)', "Forecast history CSV: ", trim(filename)
    end subroutine write_forecast_history_csv_file

    ! Convert an intraday bin number to its ending second of the trading day.
    pure integer function bin_end_seconds(bin)
        integer, intent(in) :: bin

        bin_end_seconds = default_session_start_seconds + &
                          (bin - 1) * default_bar_minutes * seconds_per_minute
    end function bin_end_seconds

    ! Format seconds since midnight as an HH:MM:SS label.
    pure function time_of_day_label(seconds) result(label)
        integer, intent(in) :: seconds
        character(len=8) :: label

        write(label,'(I2.2,A,I2.2,A,I2.2)') seconds / seconds_per_hour, ":", &
            mod(seconds, seconds_per_hour) / seconds_per_minute, ":", mod(seconds, seconds_per_minute)
    end function time_of_day_label

end module fit_mcsgarch_intraday_mod

program xfit_mcsgarch_intraday
    use date_mod, only: print_program_header
    use fit_mcsgarch_intraday_mod, only: run_fit_mcsgarch_intraday
    implicit none
    call print_program_header("xfit_mcsgarch_intraday.f90")
    call run_fit_mcsgarch_intraday()
end program xfit_mcsgarch_intraday

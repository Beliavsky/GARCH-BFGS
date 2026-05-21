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
    use garch_mcsgarch_mod, only: mcsgarch_fit_result_t, fit_mcsgarch, fit_mcsgarch_nagarch
    use stats_mod, only: mean
    implicit none
    private

    real(dp), parameter :: min_daily_var = 1.0e-12_dp
    real(dp), parameter :: trading_days_per_year = 252.0_dp
    logical, parameter :: print_diurnal_curve = .true.
    logical, parameter :: write_diurnal_curve_csv = .true.
    logical, parameter :: write_forecast_history_csv = .true.
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    character(len=*), parameter :: diurnal_curve_csv_file = "mcsgarch_diurnal_curve.csv"
    character(len=*), parameter :: forecast_history_csv_file = "mcsgarch_intraday_forecasts.csv"

    public :: run_fit_mcsgarch_intraday

contains

    ! Read intraday prices, fit the MCS-GARCH model, and write diagnostic outputs.
    subroutine run_fit_mcsgarch_intraday()
        character(len=256) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        real(dp), allocatable :: returns(:), daily_var(:), diurnal_var(:), q(:)
        real(dp), allocatable :: diurnal_sym(:), q_sym(:), diurnal_nag(:), q_nag(:)
        integer, allocatable :: bin_id(:), return_dates(:)
        type(mcsgarch_fit_result_t) :: fit_sym, fit_nag
        integer :: nargs, max_iter, nobs
        real(dp) :: gtol, t_start, t0, t1, t_end, read_sec, elapsed_sec

        call cpu_time(t_start)
        filename = "c:\python\intraday_prices\spy_5min_databento.csv"
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, filename)
        max_iter = 120
        gtol = 1.0e-5_dp

        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), bars)
        call cpu_time(t1)
        read_sec = t1 - t0
        call filter_intraday_session(bars, regular_bars)
        call build_mcsgarch_inputs(regular_bars, returns, daily_var, bin_id, return_dates)
        nobs = size(returns)
        allocate(diurnal_var(nobs), q(nobs), diurnal_sym(nobs), q_sym(nobs), diurnal_nag(nobs), q_nag(nobs))

        call fit_mcsgarch(returns, daily_var, bin_id, max_iter, gtol, fit_sym, diurnal_sym, q_sym, &
                          smooth_diurnal=smooth_diurnal_curve, &
                          smooth_half_width=diurnal_smooth_half_width)
        call fit_mcsgarch_nagarch(returns, daily_var, bin_id, max_iter, gtol, fit_nag, diurnal_nag, q_nag, &
                                  smooth_diurnal=smooth_diurnal_curve, &
                                  smooth_half_width=diurnal_smooth_half_width)
        call select_best_fit_path(fit_sym, diurnal_sym, q_sym, fit_nag, diurnal_nag, q_nag, diurnal_var, q)
        call cpu_time(t_end)
        elapsed_sec = t_end - t_start
        call print_fit_summary(trim(filename), regular_bars%nobs(), returns, daily_var, bin_id, return_dates, &
                               fit_sym, fit_nag, read_sec, elapsed_sec)
        call print_model_comparison(returns, daily_var, bin_id, fit_sym, diurnal_sym, q_sym, fit_nag, diurnal_nag, q_nag)
        if (print_diurnal_curve) call print_diurnal_variance_curve(bin_id, diurnal_var, &
                                                                    smooth_diurnal_curve, &
                                                                    diurnal_smooth_half_width)
        if (write_diurnal_curve_csv) call write_diurnal_variance_curve_csv(diurnal_curve_csv_file, &
                                                                           bin_id, diurnal_var)
        if (write_forecast_history_csv) call write_forecast_history_csv_file(forecast_history_csv_file, &
                                                                             returns, daily_var, bin_id, &
                                                                             return_dates, diurnal_var, q)

        deallocate(returns, daily_var, bin_id, return_dates, diurnal_var, q, diurnal_sym, q_sym, diurnal_nag, q_nag)
    end subroutine run_fit_mcsgarch_intraday

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

    ! Select the fitted variance path from the model with the higher log likelihood.
    subroutine select_best_fit_path(fit_sym, diurnal_sym, q_sym, fit_nag, diurnal_nag, q_nag, diurnal_best, q_best)
        type(mcsgarch_fit_result_t), intent(in) :: fit_sym, fit_nag
        real(dp), intent(in) :: diurnal_sym(:), q_sym(:), diurnal_nag(:), q_nag(:)
        real(dp), intent(out) :: diurnal_best(:), q_best(:)

        if (fit_nag%loglik > fit_sym%loglik) then
            diurnal_best = diurnal_nag
            q_best = q_nag
        else
            diurnal_best = diurnal_sym
            q_best = q_sym
        end if
    end subroutine select_best_fit_path

    ! Print the fitted MCS-GARCH parameter rows and data/timing summary.
    subroutine print_fit_summary(filename, n_regular_bars, returns, daily_var, bin_id, return_dates, fit, &
                                 fit_nag, read_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_regular_bars
        real(dp), intent(in) :: returns(:), daily_var(:)
        integer, intent(in) :: bin_id(:), return_dates(:)
        type(mcsgarch_fit_result_t), intent(in) :: fit, fit_nag
        real(dp), intent(in) :: read_sec, elapsed_sec

        print '(A)', "MCS-GARCH intraday fits"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0)', "Regular-session bars: ", n_regular_bars
        print '(A,I0,A,A,A,A)', "Intraday returns used: ", size(returns), " from ", &
              date_label(return_dates(1)), " to ", date_label(return_dates(size(return_dates)))
        print '(A)', "Returns are within-session close-to-close log returns; overnight returns enter DailyVar only."
        print '(A,I0,A,I0)', "Intraday bin range: ", minval(bin_id), " to ", maxval(bin_id)
        print '(A,ES12.4,A,ES12.4)', "Mean DailyVar: ", mean(daily_var), "  mean return variance proxy: ", mean(returns**2)
        print '(A)', "-----------------------------------------------------------------------------------------"
        print '(A,10X,A,8X,A,8X,A,8X,A,7X,A,6X,A,8X,A,8X,A)', &
              "Model", "omega", "alpha", "beta", "theta", "persist", "iter", "conv", "logL"
        print '(A)', "-----------------------------------------------------------------------------------------"
        print '(A,1X,ES12.4,4(1X,F10.4),1X,I8,5X,L1,1X,F13.3)', "MCSGARCH", &
              fit%params%omega, fit%params%alpha, fit%params%beta, fit%params%theta, fit%persist, &
              fit%niter, fit%converged, fit%loglik
        print '(A,1X,ES12.4,4(1X,F10.4),1X,I8,5X,L1,1X,F13.3)', "MCSNAGARCH", &
              fit_nag%params%omega, fit_nag%params%alpha, fit_nag%params%beta, fit_nag%params%theta, &
              fit_nag%persist, fit_nag%niter, fit_nag%converged, fit_nag%loglik
        print '(A)', "-----------------------------------------------------------------------------------------"
        print '(A,F10.3)', "Data read seconds: ", read_sec
        print '(A,F10.3)', "Elapsed seconds:   ", elapsed_sec
    end subroutine print_fit_summary

    ! Compare nested intraday volatility models with Gaussian logL, AIC, and BIC.
    subroutine print_model_comparison(returns, daily_var, bin_id, fit_sym, diurnal_sym, q_sym, &
                                      fit_nag, diurnal_nag, q_nag)
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_sym(:), q_sym(:), diurnal_nag(:), q_nag(:)
        integer, intent(in) :: bin_id(:)
        type(mcsgarch_fit_result_t), intent(in) :: fit_sym, fit_nag
        real(dp), allocatable :: h(:), unit_scale(:)
        real(dp) :: loglik, aic, bic, sigma2
        integer :: nobs, k, nbins_used

        nobs = size(returns)
        allocate(h(nobs), unit_scale(nobs))
        unit_scale = 1.0_dp

        print '(A)', ""
        print '(A)', "Intraday volatility model comparison"
        print '(A)', "----------------------------------------------------------------------"
        print '(A36,1X,A8,1X,A14,1X,A14,1X,A14)', "Model", "k", "logL", "AIC", "BIC"
        print '(A)', "----------------------------------------------------------------------"

        sigma2 = max(sum(returns**2) / real(nobs, dp), min_daily_var)
        h = sigma2
        k = 1
        call print_model_comparison_row("constant intraday volatility", k, gaussian_loglik(returns, h), nobs)

        call bin_scaled_variance_path(returns, unit_scale, bin_id, h, nbins_used)
        k = nbins_used
        call print_model_comparison_row("diurnal multiplier only", k, gaussian_loglik(returns, h), nobs)

        call bin_scaled_variance_path(returns, daily_var, bin_id, h, nbins_used)
        k = nbins_used
        call print_model_comparison_row("daily volatility x diurnal", k, gaussian_loglik(returns, h), nobs)

        h = max(daily_var * diurnal_sym * q_sym, min_daily_var)
        k = count_populated_bins(bin_id) + 3
        loglik = fit_sym%loglik
        aic = -2.0_dp*loglik + 2.0_dp*real(k, dp)
        bic = -2.0_dp*loglik + real(k, dp)*log(real(nobs, dp))
        print '(A36,1X,I8,1X,F14.3,1X,F14.3,1X,F14.3)', "MCS-GARCH", k, loglik, aic, bic

        h = max(daily_var * diurnal_nag * q_nag, min_daily_var)
        k = count_populated_bins(bin_id) + 4
        loglik = fit_nag%loglik
        aic = -2.0_dp*loglik + 2.0_dp*real(k, dp)
        bic = -2.0_dp*loglik + real(k, dp)*log(real(nobs, dp))
        print '(A36,1X,I8,1X,F14.3,1X,F14.3,1X,F14.3)', "MCS-NAGARCH", k, loglik, aic, bic

        print '(A)', "----------------------------------------------------------------------"
        deallocate(h, unit_scale)
    end subroutine print_model_comparison

    ! Print one row of the intraday volatility model comparison table.
    subroutine print_model_comparison_row(model_name, k, loglik, nobs)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: k, nobs
        real(dp), intent(in) :: loglik
        real(dp) :: aic, bic

        aic = -2.0_dp*loglik + 2.0_dp*real(k, dp)
        bic = -2.0_dp*loglik + real(k, dp)*log(real(nobs, dp))
        print '(A36,1X,I8,1X,F14.3,1X,F14.3,1X,F14.3)', model_name, k, loglik, aic, bic
    end subroutine print_model_comparison_row

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

    ! Write intraday return, daily-volatility, and MCS-GARCH forecast history to CSV.
    subroutine write_forecast_history_csv_file(filename, returns, daily_var, bin_id, return_dates, diurnal_var, q)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), q(:)
        integer, intent(in) :: bin_id(:), return_dates(:)
        integer :: unit, i, intraday_returns_per_day
        real(dp) :: intraday_var, daily_vol_ann_pct, intraday_vol_ann_pct

        if (size(returns) /= size(daily_var) .or. size(returns) /= size(diurnal_var) .or. &
            size(returns) /= size(q) .or. size(returns) /= size(bin_id) .or. &
            size(returns) /= size(return_dates)) then
            error stop "write_forecast_history_csv_file: array sizes differ"
        end if

        intraday_returns_per_day = maxval(bin_id) - minval(bin_id) + 1
        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit,'(A)') "date,time,bin,return_pct,daily_var,daily_vol_ann_pct,diurnal_var,q,intraday_var,intraday_vol_ann_pct"
        do i = 1, size(returns)
            intraday_var = max(daily_var(i) * diurnal_var(i) * q(i), min_daily_var)
            daily_vol_ann_pct = 100.0_dp * sqrt(max(daily_var(i), min_daily_var) * trading_days_per_year)
            intraday_vol_ann_pct = 100.0_dp * sqrt(intraday_var * real(intraday_returns_per_day, dp) * &
                                                   trading_days_per_year)
            write(unit,'(A,",",A,",",I0,",",ES16.8,",",ES16.8,",",F12.6,",",ES16.8,",",ES16.8,",",ES16.8,",",F12.6)') &
                date_label(return_dates(i)), time_of_day_label(bin_end_seconds(bin_id(i))), bin_id(i), &
                100.0_dp*returns(i), daily_var(i), daily_vol_ann_pct, diurnal_var(i), q(i), &
                intraday_var, intraday_vol_ann_pct
        end do
        close(unit)
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
    use fit_mcsgarch_intraday_mod, only: run_fit_mcsgarch_intraday
    implicit none
    call run_fit_mcsgarch_intraday()
end program xfit_mcsgarch_intraday

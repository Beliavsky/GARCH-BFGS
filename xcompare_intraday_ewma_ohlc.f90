! Compare simple 5-minute volatility EWMA baselines using close-close and OHLC proxies.

module compare_intraday_ewma_ohlc_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use date_mod, only: date_label, yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, intraday_bin_ids
    use intraday_vol_baseline_mod, only: fit_lag1_diurnal_baseline, fit_lag1_diurnal_intraday_ewma_baseline, &
                                         parkinson_variance_proxy, garman_klass_variance_proxy, &
                                         intraday_ewma_multiplier_from_proxy, intraday_variance_forecast
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    implicit none
    private

    character(len=*), parameter :: file_pattern = "c:\python\intraday_prices\*.csv"
    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: intraday_lambda = 0.94_dp
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2

    type :: ewma_row_t
        character(len=24) :: model = ""
        integer :: k = 0
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
    end type ewma_row_t

    public :: run_compare_intraday_ewma_ohlc

contains

    ! Compare EWMA update proxies for each input file.
    subroutine run_compare_intraday_ewma_ohlc()
        character(len=MAX_PATH_LEN), allocatable :: filenames(:)
        integer :: i

        call collect_input_filenames(filenames, &
            file_pattern=file_pattern, &
            default_filenames=[character(len=MAX_PATH_LEN) :: &
                "c:\python\intraday_prices\spy_5min_databento.csv"])
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call compare_one_file(trim(filenames(i)))
        end do
        deallocate(filenames)
    end subroutine run_compare_intraday_ewma_ohlc

    ! Read intraday OHLCV data and compare EWMA update proxies by Gaussian likelihood.
    subroutine compare_one_file(filename)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        real(dp), allocatable :: returns(:), open_px(:), high_px(:), low_px(:), close_px(:)
        real(dp), allocatable :: daily_var(:), diurnal_var(:), q(:), h(:), proxy(:)
        integer, allocatable :: bin_id(:), date_id(:)
        type(ewma_row_t) :: rows(4)
        integer :: nobs
        real(dp) :: t0, t1, read_sec, elapsed_sec

        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), bars)
        call cpu_time(t1)
        read_sec = t1 - t0
        call filter_intraday_session(bars, regular_bars)
        call build_return_ohlc_inputs(regular_bars, returns, open_px, high_px, low_px, close_px, bin_id, date_id)
        nobs = size(returns)
        allocate(daily_var(nobs), diurnal_var(nobs), q(nobs), h(nobs), proxy(nobs))

        call fit_lag1_diurnal_baseline(returns, date_id, bin_id, daily_var, diurnal_var, h, &
                                       smooth_diurnal=smooth_diurnal_curve, &
                                       smooth_half_width=diurnal_smooth_half_width)
        call fill_row(rows(1), "no intraday update", 0, returns, h)

        call fit_lag1_diurnal_intraday_ewma_baseline(returns, date_id, bin_id, intraday_lambda, &
                                                     daily_var, diurnal_var, q, h, &
                                                     smooth_diurnal=smooth_diurnal_curve, &
                                                     smooth_half_width=diurnal_smooth_half_width)
        call fill_row(rows(2), "EWMA close-close", 1, returns, h)

        call parkinson_variance_proxy(high_px, low_px, proxy)
        call intraday_ewma_multiplier_from_proxy(proxy, date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(rows(3), "EWMA Parkinson", 1, returns, h)

        call garman_klass_variance_proxy(open_px, high_px, low_px, close_px, proxy)
        call intraday_ewma_multiplier_from_proxy(proxy, date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(rows(4), "EWMA Garman-Klass", 1, returns, h)

        call cpu_time(t1)
        elapsed_sec = t1 - t0
        call print_summary(trim(filename), regular_bars%nobs(), returns, date_id, bin_id, rows, read_sec, elapsed_sec)

        deallocate(returns, open_px, high_px, low_px, close_px, daily_var, diurnal_var, q, h, proxy, bin_id, date_id)
    end subroutine compare_one_file

    ! Convert regular-session bars to aligned returns and OHLC bars.
    subroutine build_return_ohlc_inputs(bars, returns, open_px, high_px, low_px, close_px, bin_id, date_id)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: returns(:), open_px(:), high_px(:), low_px(:), close_px(:)
        integer, allocatable, intent(out) :: bin_id(:), date_id(:)
        integer, allocatable :: bar_bins(:), ret_bin(:), ret_date(:)
        real(dp), allocatable :: ret_all(:), open_all(:), high_all(:), low_all(:), close_all(:)
        integer :: n, i, k, this_date, prev_date

        n = bars%nobs()
        if (n < 3) error stop "build_return_ohlc_inputs: not enough regular-session bars"
        call intraday_bin_ids(bars, bar_bins)
        allocate(ret_all(n - 1), open_all(n - 1), high_all(n - 1), low_all(n - 1), close_all(n - 1))
        allocate(ret_bin(n - 1), ret_date(n - 1))
        k = 0
        do i = 2, n
            this_date = yyyymmdd(bars%timestamp(i)%date)
            prev_date = yyyymmdd(bars%timestamp(i - 1)%date)
            if (this_date /= prev_date) cycle
            k = k + 1
            ret_all(k) = log(bars%close(i) / bars%close(i - 1))
            open_all(k) = bars%open(i)
            high_all(k) = bars%high(i)
            low_all(k) = bars%low(i)
            close_all(k) = bars%close(i)
            ret_bin(k) = bar_bins(i)
            ret_date(k) = this_date
        end do
        if (k < 3) error stop "build_return_ohlc_inputs: not enough intraday returns"
        allocate(returns(k), open_px(k), high_px(k), low_px(k), close_px(k), bin_id(k), date_id(k))
        returns = ret_all(1:k)
        open_px = open_all(1:k)
        high_px = high_all(1:k)
        low_px = low_all(1:k)
        close_px = close_all(1:k)
        bin_id = ret_bin(1:k)
        date_id = ret_date(1:k)
        deallocate(bar_bins, ret_all, open_all, high_all, low_all, close_all, ret_bin, ret_date)
    end subroutine build_return_ohlc_inputs

    ! Fill one likelihood comparison row.
    subroutine fill_row(row, model, k, returns, h)
        type(ewma_row_t), intent(out) :: row
        character(len=*), intent(in) :: model
        integer, intent(in) :: k
        real(dp), intent(in) :: returns(:), h(:)

        row%model = model
        row%k = k
        row%loglik = gaussian_loglik(returns, h)
        row%aic = -2.0_dp*row%loglik + 2.0_dp*real(k, dp)
        row%bic = -2.0_dp*row%loglik + log(real(size(returns), dp))*real(k, dp)
    end subroutine fill_row

    ! Gaussian zero-mean log likelihood for returns and variance forecasts.
    real(dp) function gaussian_loglik(returns, h)
        real(dp), intent(in) :: returns(:), h(:)
        integer :: i

        if (size(returns) /= size(h)) error stop "gaussian_loglik: array sizes differ"
        gaussian_loglik = 0.0_dp
        do i = 1, size(returns)
            gaussian_loglik = gaussian_loglik - log_sqrt_2pi - 0.5_dp*log(max(h(i), min_var)) - &
                              0.5_dp*returns(i)**2 / max(h(i), min_var)
        end do
    end function gaussian_loglik

    ! Print the EWMA proxy comparison table.
    subroutine print_summary(filename, n_regular_bars, returns, date_id, bin_id, rows, read_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_regular_bars, date_id(:), bin_id(:)
        real(dp), intent(in) :: returns(:), read_sec, elapsed_sec
        type(ewma_row_t), intent(in) :: rows(:)
        integer :: i, aic_best, bic_best

        aic_best = minloc(rows%aic, dim=1)
        bic_best = minloc(rows%bic, dim=1)
        print '(A)', "Intraday EWMA OHLC volatility proxy comparison"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0)', "Regular-session bars: ", n_regular_bars
        print '(A,I0,A,A,A,A)', "Intraday returns: ", size(returns), " from ", date_label(date_id(1)), &
              " to ", date_label(date_id(size(date_id)))
        print '(A,I0)', "Intraday bins: ", maxval(bin_id)
        print '(A,F8.4)', "Intraday EWMA lambda: ", intraday_lambda
        print '(A)', "--------------------------------------------------------------------------"
        print '(A24,1X,A5,1X,A14,1X,A14,1X,A14)', "Model", "k", "logL", "AIC", "BIC"
        print '(A)', "--------------------------------------------------------------------------"
        do i = 1, size(rows)
            print '(A24,1X,I5,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), rows(i)%k, &
                  rows(i)%loglik, rows(i)%aic, rows(i)%bic
        end do
        print '(A)', "--------------------------------------------------------------------------"
        print '(A,A)', "AIC selects: ", trim(rows(aic_best)%model)
        print '(A,A)', "BIC selects: ", trim(rows(bic_best)%model)
        print '(A,F10.3)', "Data read seconds: ", read_sec
        print '(A,F10.3)', "Elapsed seconds:   ", elapsed_sec
    end subroutine print_summary

end module compare_intraday_ewma_ohlc_mod

program xcompare_intraday_ewma_ohlc
    use date_mod, only: print_program_header
    use compare_intraday_ewma_ohlc_mod, only: run_compare_intraday_ewma_ohlc
    implicit none

    call print_program_header("xcompare_intraday_ewma_ohlc.f90")
    call run_compare_intraday_ewma_ohlc()
end program xcompare_intraday_ewma_ohlc

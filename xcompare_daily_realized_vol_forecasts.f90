! Compare daily volatility forecasts built from intraday realized measures.
!
! The target is daily close-to-close return likelihood.  Realized measures are
! computed from regular-session intraday OHLC bars and used as lagged predictors.

module compare_daily_realized_vol_forecasts_mod
    use kind_mod, only: dp
    use date_mod, only: print_program_header, date_label
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session
    use csv_mod, only: read_price_csv
    use intraday_realized_measures_mod, only: daily_realized_panel_t, build_daily_realized_panel, &
                                              select_realized_measure
    use implied_vol_utils_mod, only: adjust_implied_vol_weekday_levels
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod, only: fit_symm_garch, fit_nagarch, symm_garch_persist, nagarch_persist
    use garch_forecast_mod, only: symm_garch_variance_path, nagarch_variance_path, riskmetrics_ewma_variance_path
    use realized_vol_forecast_mod, only: ewma_affine_result_t, har_variance_result_t, log_har_variance_result_t, &
                                         harq_variance_result_t, harj_variance_result_t, harqj_variance_result_t, &
                                         har_negret_variance_result_t, &
                                         semivar_har_variance_result_t, &
                                         midas_variance_result_t, heavy_variance_result_t, &
                                         fit_ewma_affine_variance, fit_har_variance, fit_log_har_variance, &
                                         fit_harq_variance, fit_harj_variance, fit_harqj_variance, &
                                         fit_har_negret_variance, fit_semivar_har_variance, &
                                         fit_midas_variance, fit_heavy_variance, &
                                         gaussian_variance_loglik, qlike_loss
    use realized_garch_mod, only: realized_garch_result_t, fit_realized_garch
    use stats_mod, only: mean, correlation, demean_first
    use time_series_compare_mod, only: aligned_index_pairs
    use rank_mod, only: rank_desc
    use program_utils_mod, only: read_integer_arg, elapsed_since
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: riskmetrics_lambda = 0.94_dp
    real(dp), parameter :: gtol = 1.0e-5_dp
    integer, parameter :: default_ntest_days = 250
    integer, parameter :: max_iter = 300
    integer, parameter :: midas_lag_count = 22
    real(dp), parameter :: trading_days_per_year = 252.0_dp
    logical, parameter :: adjust_implied_vol_day_of_week = .true.

    character(len=16), parameter :: measure_names(*) = [character(len=16) :: &
        "RV", "BPV", "RSV_NEG", "RSV_POS", "PARKINSON", "GARMAN_KLASS"]

    type :: score_row_t
        character(len=16) :: measure = ""
        character(len=16) :: model = ""
        integer :: k = 0
        real(dp) :: lambda = 0.0_dp
        real(dp) :: a = 0.0_dp
        real(dp) :: b = 0.0_dp
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: qlike = 0.0_dp
        real(dp) :: vol_ann = 0.0_dp
        integer :: niter = 0
        integer :: aic_rank = 0
        integer :: bic_rank = 0
        logical :: converged = .false.
        real(dp), allocatable :: h_test(:)
    end type score_row_t

    public :: run_compare_daily_realized_vol_forecasts

contains

    subroutine run_compare_daily_realized_vol_forecasts()
        character(len=256) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_realized_panel_t) :: daily
        character(len=32), allocatable :: implied_vol_names(:)
        integer, allocatable :: implied_vol_dates(:)
        real(dp), allocatable :: implied_vol_values(:, :)
        real(dp), allocatable :: ret_cc(:), h(:), x(:), jump(:)
        type(score_row_t), allocatable :: rows(:)
        type(garch_params_t) :: garch_params, nag_params
        integer :: nargs, ntest_days, nobs, train_nobs, test_first, irow, imeas
        integer :: niter_garch, niter_nag
        logical :: conv_garch, conv_nag
        real(dp) :: f_garch, f_nag, t0, t1, read_sec, fit_sec, elapsed_sec

        call cpu_time(t0)
        filename = "c:\python\intraday_prices\spy_5min_databento.csv"
        ntest_days = default_ntest_days
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, filename)
        if (nargs >= 2) call read_integer_arg(2, ntest_days)

        call print_program_header("xcompare_daily_realized_vol_forecasts.f90")
        call cpu_time(t1)
        call read_intraday_prices_csv(filename, bars)
        call filter_intraday_session(bars, regular_bars)
        call read_price_csv("vix_spy.csv", implied_vol_dates, implied_vol_names, implied_vol_values, &
                            selected_cols=[character(len=3) :: "VIX"])
        read_sec = elapsed_since(t1)

        call build_daily_realized_panel(regular_bars, daily)
        nobs = size(daily%date) - 1
        if (nobs <= ntest_days + 30) error stop "xcompare_daily_realized_vol_forecasts: not enough daily observations"
        train_nobs = nobs - ntest_days
        test_first = train_nobs + 1
        allocate(ret_cc(nobs), h(nobs), x(nobs), jump(nobs), rows(7 + 11*size(measure_names)))
        ret_cc = log(daily%close(2:nobs + 1) / daily%close(1:nobs))
        jump = max(daily%rv(1:nobs) - daily%bpv(1:nobs), 0.0_dp)
        call demean_first(ret_cc, train_nobs)

        call cpu_time(t1)
        call fit_symm_garch(ret_cc(1:train_nobs), max_iter, gtol, f_garch, garch_params, niter_garch, conv_garch)
        call fit_nagarch(ret_cc(1:train_nobs), max_iter, gtol, f_nag, nag_params, niter_nag, conv_nag)
        call symm_garch_variance_path(ret_cc, garch_params, h)
        call fill_row(rows(1), "CC", "GARCH", 3, 0.0_dp, 0.0_dp, 0.0_dp, ret_cc(test_first:nobs), &
                      h(test_first:nobs), niter_garch, conv_garch)
        call riskmetrics_ewma_variance_path(ret_cc, riskmetrics_lambda, train_nobs, h)
        call fill_row(rows(2), "CC", "EWMA_0P94", 0, riskmetrics_lambda, 0.0_dp, 0.0_dp, ret_cc(test_first:nobs), &
                      h(test_first:nobs), 0, .true.)
        call nagarch_variance_path(ret_cc, nag_params, h)
        call fill_row(rows(3), "CC", "NAGARCH", 4, 0.0_dp, 0.0_dp, 0.0_dp, ret_cc(test_first:nobs), &
                      h(test_first:nobs), niter_nag, conv_nag)

        irow = 3
        do imeas = 1, size(measure_names)
            call select_realized_measure(daily, measure_names(imeas), x)
            irow = irow + 1
            call fit_and_score_measure(rows(irow), measure_names(imeas), "EWMA_0P94", ret_cc, x, &
                                       train_nobs, test_first, .false., h)
            irow = irow + 1
            call fit_and_score_measure(rows(irow), measure_names(imeas), "EWMA_FIT", ret_cc, x, &
                                       train_nobs, test_first, .true., h)
            irow = irow + 1
            call fit_and_score_har(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_log_har(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_harq(rows(irow), measure_names(imeas), ret_cc, x, daily%rq(1:nobs), train_nobs, &
                                    test_first, h)
            irow = irow + 1
            call fit_and_score_harj(rows(irow), measure_names(imeas), ret_cc, x, jump, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_harqj(rows(irow), measure_names(imeas), ret_cc, x, daily%rq(1:nobs), jump, &
                                     train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_har_negret(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_midas(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_heavy(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_realgarch(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
        end do
        irow = irow + 1
        call fit_and_score_semivar_har(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_harcj(rows(irow), daily, ret_cc, jump, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_char(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_charq(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        fit_sec = elapsed_since(t1)
        elapsed_sec = elapsed_since(t0)

        call rank_rows(rows)
        call print_summary(filename, daily, train_nobs, ntest_days, rows, read_sec, fit_sec, elapsed_sec)
        call print_implied_vol_correlation_table(daily%date(test_first + 1:nobs + 1), rows, implied_vol_dates, &
                                                 implied_vol_values(:, 1), "VIX", adjust_implied_vol_day_of_week)
        deallocate(ret_cc, h, x, jump, rows)
    end subroutine run_compare_daily_realized_vol_forecasts

    subroutine fit_and_score_measure(row, measure, model, ret_cc, x, train_nobs, test_first, fit_lambda, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure, model
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        logical, intent(in) :: fit_lambda
        real(dp), intent(out) :: h(:)
        type(ewma_affine_result_t) :: result
        integer :: nobs, k

        nobs = size(ret_cc)
        call fit_ewma_affine_variance(ret_cc, x, train_nobs, fit_lambda, riskmetrics_lambda, max_iter, gtol, result, h)
        k = merge(3, 2, fit_lambda)
        call fill_row(row, measure, model, k, result%lambda, result%a, result%b, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_measure

    subroutine fit_and_score_har(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(har_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_har_variance(ret_cc, x, train_nobs, result, h)
        call fill_row(row, measure, "HAR_RV", 4, 0.0_dp, result%coef(1), 0.0_dp, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_har

    subroutine fit_and_score_log_har(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(log_har_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_log_har_variance(ret_cc, x, train_nobs, result, h)
        call fill_row(row, measure, "LOG_HAR", 4, 0.0_dp, result%coef(1), result%coef(2), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_log_har

    subroutine fit_and_score_harq(row, measure, ret_cc, x, rq, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:), rq(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harq_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harq_variance(ret_cc, x, rq, train_nobs, result, h)
        call fill_row(row, measure, "HARQ", 5, 0.0_dp, result%coef(1), result%coef(5), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_harq

    subroutine fit_and_score_harj(row, measure, ret_cc, x, jump, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:), jump(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harj_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harj_variance(ret_cc, x, jump, train_nobs, result, h)
        call fill_row(row, measure, "HARJ", 7, 0.0_dp, result%coef(1), result%coef(5), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_harj

    subroutine fit_and_score_harqj(row, measure, ret_cc, x, rq, jump, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:), rq(:), jump(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harqj_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harqj_variance(ret_cc, x, rq, jump, train_nobs, result, h)
        call fill_row(row, measure, "HARQJ", 8, 0.0_dp, result%coef(1), result%coef(8), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_harqj

    subroutine fit_and_score_harcj(row, daily, ret_cc, jump, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        type(daily_realized_panel_t), intent(in) :: daily
        real(dp), intent(in) :: ret_cc(:), jump(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harj_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harj_variance(ret_cc, daily%bpv(1:nobs), jump, train_nobs, result, h)
        call fill_row(row, "BPV_JUMP", "HARCJ", 7, 0.0_dp, result%coef(1), result%coef(5), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_harcj

    subroutine fit_and_score_char(row, daily, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        type(daily_realized_panel_t), intent(in) :: daily
        real(dp), intent(in) :: ret_cc(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(har_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_har_variance(ret_cc, daily%bpv(1:nobs), train_nobs, result, h)
        call fill_row(row, "BPV", "CHAR", 4, 0.0_dp, result%coef(1), 0.0_dp, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_char

    subroutine fit_and_score_charq(row, daily, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        type(daily_realized_panel_t), intent(in) :: daily
        real(dp), intent(in) :: ret_cc(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harq_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harq_variance(ret_cc, daily%bpv(1:nobs), daily%rq(1:nobs), train_nobs, result, h)
        call fill_row(row, "BPV", "CHARQ", 5, 0.0_dp, result%coef(1), result%coef(5), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_charq

    subroutine fit_and_score_har_negret(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(har_negret_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_har_negret_variance(ret_cc, x, train_nobs, result, h)
        call fill_row(row, measure, "HAR_NEGRET", 5, 0.0_dp, result%coef(1), result%coef(5), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_har_negret

    subroutine fit_and_score_midas(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(midas_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_midas_variance(ret_cc, x, train_nobs, midas_lag_count, result, h)
        call fill_row(row, measure, "MIDAS_RV", 4, 0.0_dp, result%a, result%b, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_midas

    subroutine fit_and_score_heavy(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(heavy_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_heavy_variance(ret_cc, x, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, measure, "HEAVY", 6, 0.0_dp, result%alpha, result%beta, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_heavy

    subroutine fit_and_score_semivar_har(row, daily, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        type(daily_realized_panel_t), intent(in) :: daily
        real(dp), intent(in) :: ret_cc(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(semivar_har_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_semivar_har_variance(ret_cc, daily%rsv_pos(1:nobs), daily%rsv_neg(1:nobs), train_nobs, result, h)
        call fill_row(row, "RSV_POS_NEG", "SEMIVAR_HAR", 7, 0.0_dp, result%coef(1), 0.0_dp, &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_semivar_har

    subroutine fit_and_score_realgarch(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(realized_garch_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_realized_garch(ret_cc, x, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, measure, "REALGARCH", 8, 0.0_dp, result%params%omega, result%persist, &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
    end subroutine fit_and_score_realgarch

    subroutine fill_row(row, measure, model, k, lambda, a, b, y, h, niter, converged)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure, model
        integer, intent(in) :: k, niter
        real(dp), intent(in) :: lambda, a, b, y(:), h(:)
        logical, intent(in) :: converged
        integer :: n

        n = size(y)
        row%measure = measure
        row%model = model
        row%k = k
        row%lambda = lambda
        row%a = a
        row%b = b
        row%loglik = gaussian_variance_loglik(y, h)
        row%aic = 2.0_dp*real(k, dp) - 2.0_dp*row%loglik
        row%bic = log(real(n, dp))*real(k, dp) - 2.0_dp*row%loglik
        row%qlike = qlike_loss(y, h)
        row%vol_ann = 100.0_dp*sqrt(trading_days_per_year*mean(h))
        row%niter = niter
        row%converged = converged
        allocate(row%h_test(n))
        row%h_test = h
    end subroutine fill_row

    subroutine print_implied_vol_correlation_table(forecast_dates, rows, implied_vol_dates, implied_vol_close, &
                                                   implied_vol_label, adjust_weekday)
        integer, intent(in) :: forecast_dates(:), implied_vol_dates(:)
        type(score_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: implied_vol_close(:)
        character(len=*), intent(in) :: implied_vol_label
        logical, intent(in) :: adjust_weekday
        real(dp), allocatable :: model_vol(:), implied_vol(:), dlog_model(:), dlog_implied_vol(:)
        real(dp), allocatable :: implied_vol_used(:)
        real(dp), allocatable :: corr_level(:), corr_dlog(:)
        integer, allocatable :: corr_rank(:), dlog_rank(:), forecast_idx(:), implied_vol_idx(:)
        integer :: i, j, nmatch, ndlog

        allocate(implied_vol_used(size(implied_vol_close)))
        if (adjust_weekday) then
            call adjust_implied_vol_weekday_levels(implied_vol_dates, implied_vol_close, implied_vol_used)
        else
            implied_vol_used = implied_vol_close
        end if
        call aligned_index_pairs(forecast_dates, implied_vol_dates, forecast_idx, implied_vol_idx)
        nmatch = size(forecast_idx)
        if (nmatch < 2) then
            print '(A,A,/)', trim(implied_vol_label), " correlation table skipped: fewer than 2 aligned observations."
            return
        end if

        allocate(model_vol(nmatch), implied_vol(nmatch))
        allocate(corr_level(size(rows)), corr_dlog(size(rows)), corr_rank(size(rows)), dlog_rank(size(rows)))
        corr_level = -huge(1.0_dp)
        corr_dlog = -huge(1.0_dp)
        do i = 1, size(rows)
            do j = 1, nmatch
                if (implied_vol_used(implied_vol_idx(j)) > 0.0_dp .and. rows(i)%h_test(forecast_idx(j)) > 0.0_dp) then
                    model_vol(j) = 100.0_dp*sqrt(trading_days_per_year*rows(i)%h_test(forecast_idx(j)))
                    implied_vol(j) = implied_vol_used(implied_vol_idx(j))
                else
                    model_vol(j) = 0.0_dp
                    implied_vol(j) = 0.0_dp
                end if
            end do
            corr_level(i) = correlation(model_vol, implied_vol)
            ndlog = nmatch - 1
            if (ndlog >= 2) then
                allocate(dlog_model(ndlog), dlog_implied_vol(ndlog))
                dlog_model = log(max(model_vol(2:nmatch), min_var) / max(model_vol(1:nmatch-1), min_var))
                dlog_implied_vol = log(max(implied_vol(2:nmatch), min_var) / max(implied_vol(1:nmatch-1), min_var))
                corr_dlog(i) = correlation(dlog_model, dlog_implied_vol)
                deallocate(dlog_model, dlog_implied_vol)
            end if
        end do
        call rank_desc(corr_level, corr_rank)
        call rank_desc(corr_dlog, dlog_rank)

        print '(A)', ""
        if (adjust_weekday) then
            print '(A,A,A)', "Forecast correlation with weekday-adjusted ", trim(implied_vol_label), " close"
        else
            print '(A,A,A)', "Forecast correlation with ", trim(implied_vol_label), " close"
        end if
        print '(A,A,A,A,A,I0)', "Date range: ", date_label(forecast_dates(forecast_idx(1))), &
              " to ", date_label(forecast_dates(forecast_idx(nmatch))), "  observations: ", nmatch
        print '(A)', "------------------------------------------------------------------------------------------------------"
        print '(A16,1X,A16,1X,A10,1X,A10,1X,A9,1X,A9)', &
              "Measure", "Model", "corr_vol", "corr_dlog", "corr_rank", "dlog_rank"
        print '(A)', "------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            print '(A16,1X,A16,1X,F10.4,1X,F10.4,1X,I9,1X,I9)', rows(i)%measure, rows(i)%model, &
                  corr_level(i), corr_dlog(i), corr_rank(i), dlog_rank(i)
        end do
        print '(A)', "------------------------------------------------------------------------------------------------------"
    end subroutine print_implied_vol_correlation_table

    subroutine rank_rows(rows)
        type(score_row_t), intent(inout) :: rows(:)
        integer :: i

        do i = 1, size(rows)
            rows(i)%aic_rank = 1 + count(rows%aic < rows(i)%aic)
            rows(i)%bic_rank = 1 + count(rows%bic < rows(i)%bic)
        end do
    end subroutine rank_rows

    subroutine print_summary(filename, daily, train_nobs, ntest_days, rows, read_sec, fit_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        type(daily_realized_panel_t), intent(in) :: daily
        integer, intent(in) :: train_nobs, ntest_days
        type(score_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: read_sec, fit_sec, elapsed_sec
        integer :: i, best

        best = maxloc(rows%loglik, dim=1)
        print '(A)', "Daily CC volatility forecasts from intraday realized measures"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0,A,A,A,A)', "Daily bars: ", size(daily%date), " from ", date_label(daily%date(1)), " to ", &
              date_label(daily%date(size(daily%date)))
        print '(A,I0,A,I0,A,A)', "Training observations: ", train_nobs, "  test observations: ", ntest_days, &
              "  test starts ", date_label(daily%date(train_nobs + 2))
        print '(A,I0)', "MIDAS lag count: ", midas_lag_count
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------------------------------"
        print '(A16,1X,A16,1X,A3,1X,A8,1X,A10,1X,A10,1X,A12,1X,A12,1X,A12,1X,A8,1X,A8,1X,A9,1X,A9,1X,A5,1X,A4)', &
              "Measure", "Model", "k", "lambda", "a", "b", "logL", "AIC", "BIC", "AIC_rank", "BIC_rank", &
              "QLIKE", "vol_ann%", "iter", "conv"
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            print '(A16,1X,A16,1X,I3,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,F12.3,1X,F12.3,1X,F12.3,1X,I8,1X,I8,1X,F9.5,1X,F9.3,1X,I5,4X,L1)', &
                  rows(i)%measure, rows(i)%model, rows(i)%k, rows(i)%lambda, rows(i)%a, rows(i)%b, rows(i)%loglik, &
                  rows(i)%aic, rows(i)%bic, rows(i)%aic_rank, rows(i)%bic_rank, rows(i)%qlike, rows(i)%vol_ann, &
                  rows(i)%niter, rows(i)%converged
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------------------------------------------------------"
        print '(A,A,1X,A)', "Best logL: ", trim(rows(best)%measure), trim(rows(best)%model)
        print '(A,F8.3,2X,A,F8.3,2X,A,F8.3)', "read_sec=", read_sec, "fit_sec=", fit_sec, "elapsed_sec=", elapsed_sec
    end subroutine print_summary

end module compare_daily_realized_vol_forecasts_mod

program xcompare_daily_realized_vol_forecasts
    use compare_daily_realized_vol_forecasts_mod, only: run_compare_daily_realized_vol_forecasts
    implicit none

    call run_compare_daily_realized_vol_forecasts()
end program xcompare_daily_realized_vol_forecasts

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
    use garch_fit_mod, only: fit_ewma, fit_symm_garch, fit_qgarch, fit_csgarch, fit_nagarch, fit_gjr, &
                             fit_gjr_signed, fit_egarch, fit_fgarch_twist, fit_aparch, fit_harch, fit_tgarch, &
                             fit_avgarch, fit_riskmetrics2006, fit_midas_hyperbolic, fit_midas_hyperbolic_asym
    use garch_fit_mod, only: ewma_persist, symm_garch_persist, qgarch_persist, csgarch_persist, &
                             nagarch_persist, gjr_persist, egarch_persist, fgarch_twist_persist, &
                             aparch_persist, harch_persist, tgarch_persist, avgarch_persist, &
                             riskmetrics2006_persist, midas_hyperbolic_persist, midas_hyperbolic_asym_persist
    use garch_forecast_mod, only: symm_garch_variance_path, nagarch_variance_path, gjr_variance_path, &
                                  egarch_variance_path, model_vol_forecast
    use realized_vol_forecast_mod, only: ewma_affine_result_t, affine_variance_result_t, affine2_variance_result_t, &
                                         har_variance_result_t, harx_variance_result_t, harx_lev_variance_result_t, &
                                         log_har_variance_result_t, &
                                         sqrt_har_variance_result_t, harq_variance_result_t, harj_variance_result_t, &
                                         harqj_variance_result_t, har_negret_variance_result_t, &
                                         har_lev_variance_result_t, &
                                         semivar_har_variance_result_t, &
                                         midas_variance_result_t, heavy_variance_result_t, &
                                         fit_ewma_affine_variance, fit_affine_variance, fit_affine2_variance, &
                                         fit_har_variance, fit_harx_variance, fit_harx_lev_variance, &
                                         fit_log_har_variance, fit_sqrt_har_variance, fit_harq_variance, fit_harj_variance, &
                                         fit_harqj_variance, fit_har_negret_variance, fit_har_lev_variance, &
                                         fit_semivar_har_variance, &
                                         fit_midas_variance, fit_heavy_variance, &
                                         gaussian_variance_loglik, qlike_loss
    use realized_garch_mod, only: realized_garch_result_t, fit_realized_garch, &
                                  realized_egarch_result_t, fit_realized_egarch
    use msgarch_mod, only: msgarch_result_t, fit_msgarch
    use msnagarch_mod, only: msnagarch_result_t, fit_msnagarch
    use stats_mod, only: mean, sd, correlation, demean_first, nonnegative_least_squares
    use regression_diagnostics_mod, only: active_nnls_tstats, predictor_correlations, rank_nonzero_by_tstat, &
                                          print_response_predictor_corr_matrix
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
    logical, parameter :: fit_implied_vol_only = .true.
    logical, parameter :: fit_harx_implied_vol = .true.
    logical, parameter :: fit_harx_implied_vol_leverage = .true.
    logical, parameter :: fit_iv_augmented_models = .false.
    logical, parameter :: fit_implied_vol_nnls = .false.
    logical, parameter :: print_implied_vol_nnls_corr_matrix = .true.
    logical, parameter :: implied_vol_nnls_intercept = .true.
    integer, parameter :: implied_vol_nnls_max_iter = 1000
    real(dp), parameter :: implied_vol_nnls_tol = 1.0e-10_dp

    character(len=16), parameter :: measure_names(*) = [character(len=16) :: &
        "RV", "BPV", "RSV_NEG", "RSV_POS", "PARKINSON", "GARMAN_KLASS"]
    character(len=16), parameter :: cc_garch_models(*) = [character(len=16) :: &
        "GARCH", "NAGARCH", "GJR", "EGARCH", "EWMA_FIT", "RM2006", &
        "QGARCH", "CSGARCH", "GJR_SIGNED", "APARCH", "HARCH", &
        "TGARCH", "AVGARCH", "FGTWIST", "MIDASHYP", "MIDASHYP_ASYM"]

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
        real(dp), allocatable :: h_all(:)
    end type score_row_t

    public :: run_compare_daily_realized_vol_forecasts

contains

    subroutine run_compare_daily_realized_vol_forecasts()
        character(len=256) :: intraday_prices_file, implied_vol_file
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_realized_panel_t) :: daily
        character(len=32), allocatable :: implied_vol_names(:)
        integer, allocatable :: implied_vol_dates(:)
        real(dp), allocatable :: implied_vol_values(:, :)
        real(dp), allocatable :: ret_cc(:), h(:), x(:), jump(:), implied_var(:)
        type(score_row_t), allocatable :: rows(:)
        integer :: nargs, ntest_days, nobs, train_nobs, test_first, irow, imeas, imodel
        real(dp) :: t0, t1, read_sec, fit_sec, elapsed_sec

        call cpu_time(t0)
        intraday_prices_file = "c:\python\databento\data_1min\spy_1min_databento.csv" ! "c:\python\intraday_prices\spy_5min_databento.csv"
        implied_vol_file = "vix_spy.csv"
        ntest_days = default_ntest_days
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, intraday_prices_file)
        if (nargs >= 2) call read_integer_arg(2, ntest_days)

        call print_program_header("xcompare_daily_realized_vol_forecasts.f90")
        call cpu_time(t1)
        call read_intraday_prices_csv(intraday_prices_file, bars)
        call filter_intraday_session(bars, regular_bars)
        call read_price_csv(implied_vol_file, implied_vol_dates, implied_vol_names, implied_vol_values, &
                            selected_cols=[character(len=3) :: "VIX"])
        read_sec = elapsed_since(t1)

        call build_daily_realized_panel(regular_bars, daily)
        nobs = size(daily%date) - 1
        if (nobs <= ntest_days + 30) error stop "xcompare_daily_realized_vol_forecasts: not enough daily observations"
        train_nobs = nobs - ntest_days
        test_first = train_nobs + 1
        allocate(ret_cc(nobs), h(nobs), x(nobs), jump(nobs), implied_var(nobs), &
                 rows(2*(size(cc_garch_models) + merge(1, 0, fit_implied_vol_only) + &
                      (14 + merge(1, 0, fit_harx_implied_vol) + &
                       merge(1, 0, fit_harx_implied_vol_leverage))*size(measure_names))))
        ret_cc = log(daily%close(2:nobs + 1) / daily%close(1:nobs))
        jump = max(daily%rv(1:nobs) - daily%bpv(1:nobs), 0.0_dp)
        call build_lagged_implied_variance(daily%date, implied_vol_dates, implied_vol_values(:, 1), &
                                           adjust_implied_vol_day_of_week, implied_var)
        call demean_first(ret_cc, train_nobs)

        call cpu_time(t1)
        irow = 0
        do imodel = 1, size(cc_garch_models)
            irow = irow + 1
            call fit_and_score_cc_garch(rows(irow), trim(cc_garch_models(imodel)), ret_cc, train_nobs, test_first, h)
        end do
        if (fit_implied_vol_only) then
            irow = irow + 1
            call fit_and_score_iv_only(rows(irow), ret_cc, implied_var, train_nobs, test_first, h)
        end if
        irow = irow + 1
        call fit_and_score_msgarch(rows(irow), ret_cc, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_msnagarch(rows(irow), ret_cc, train_nobs, test_first, h)
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
            if (fit_harx_implied_vol) then
                irow = irow + 1
                call fit_and_score_harx_iv(rows(irow), measure_names(imeas), ret_cc, x, implied_var, &
                                           train_nobs, test_first, h)
            end if
            if (fit_harx_implied_vol_leverage) then
                irow = irow + 1
                call fit_and_score_harx_iv_lev(rows(irow), measure_names(imeas), ret_cc, x, implied_var, &
                                               train_nobs, test_first, h)
            end if
            irow = irow + 1
            call fit_and_score_log_har(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_sqrt_har(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
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
            call fit_and_score_har_lev(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_midas(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_heavy(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_realgarch(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
            irow = irow + 1
            call fit_and_score_realegarch(rows(irow), measure_names(imeas), ret_cc, x, train_nobs, test_first, h)
        end do
        irow = irow + 1
        call fit_and_score_semivar_har(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_harcj(rows(irow), daily, ret_cc, jump, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_char(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        irow = irow + 1
        call fit_and_score_charq(rows(irow), daily, ret_cc, train_nobs, test_first, h)
        if (fit_iv_augmented_models) call append_iv_augmented_rows(rows, irow, ret_cc, implied_var, train_nobs, test_first, h)
        fit_sec = elapsed_since(t1)
        elapsed_sec = elapsed_since(t0)

        call rank_rows(rows(1:irow))
        call print_summary(intraday_prices_file, implied_vol_file, daily, train_nobs, ntest_days, rows(1:irow), read_sec, fit_sec, elapsed_sec)
        call print_implied_vol_correlation_table(daily%date(test_first + 1:nobs + 1), rows(1:irow), implied_vol_dates, &
                                                 implied_vol_values(:, 1), "VIX", adjust_implied_vol_day_of_week)
        deallocate(ret_cc, h, x, jump, implied_var, rows)
    end subroutine run_compare_daily_realized_vol_forecasts

    subroutine fit_and_score_cc_garch(row, model, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: model
        real(dp), intent(in) :: ret_cc(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(garch_params_t) :: params
        character(len=16) :: path_model, row_model
        real(dp) :: fopt, ret_std, persist, vol_ann, skew, ekurt
        integer :: nobs, niter
        logical :: converged

        nobs = size(ret_cc)
        params = garch_params_t()
        row_model = cc_garch_row_model(model)
        path_model = cc_garch_path_model(model)
        select case (trim(row_model))
        case ("EWMA_FIT")
            call fit_ewma(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = ewma_persist(params)
        case ("GARCH")
            call fit_symm_garch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = symm_garch_persist(params)
        case ("QGARCH")
            call fit_qgarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = qgarch_persist(params)
        case ("CSGARCH")
            call fit_csgarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = csgarch_persist(params)
        case ("NAGARCH")
            call fit_nagarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = nagarch_persist(params)
        case ("GJR")
            call fit_gjr(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = gjr_persist(params)
        case ("GJR_SIGNED")
            call fit_gjr_signed(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = gjr_persist(params)
        case ("EGARCH")
            call fit_egarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = egarch_persist(params)
        case ("APARCH")
            call fit_aparch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = aparch_persist(params)
        case ("HARCH")
            call fit_harch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = harch_persist(params)
        case ("TGARCH")
            call fit_tgarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = tgarch_persist(params)
        case ("AVGARCH")
            call fit_avgarch(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = avgarch_persist(params)
        case ("FGTWIST")
            ret_std = sd(ret_cc(1:train_nobs))
            call fit_fgarch_twist(ret_cc(1:train_nobs), ret_std, max_iter, gtol, fopt, params, vol_ann, skew, ekurt, &
                                  niter, converged)
            persist = fgarch_twist_persist(params)
        case ("RM2006")
            call fit_riskmetrics2006(ret_cc(1:train_nobs), fopt, params, niter, converged)
            persist = riskmetrics2006_persist()
        case ("MIDASHYP")
            call fit_midas_hyperbolic(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = midas_hyperbolic_persist(params)
        case ("MIDASHYP_ASYM")
            call fit_midas_hyperbolic_asym(ret_cc(1:train_nobs), max_iter, gtol, fopt, params, niter, converged)
            persist = midas_hyperbolic_asym_persist(params)
        case default
            error stop "fit_and_score_cc_garch: unsupported model"
        end select

        call cc_garch_variance_path(path_model, ret_cc, params, persist, h)
        call fill_row(row, "CC", row_model, cc_garch_param_count(row_model), 0.0_dp, 0.0_dp, 0.0_dp, &
                      ret_cc(test_first:nobs), h(test_first:nobs), niter, converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_cc_garch

    subroutine cc_garch_variance_path(model, ret_cc, params, persist, h)
        character(len=*), intent(in) :: model
        real(dp), intent(in) :: ret_cc(:), persist
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp), allocatable :: vol(:)

        select case (trim(model))
        case ("SYMM_GARCH")
            call symm_garch_variance_path(ret_cc, params, h)
        case ("NAGARCH")
            call nagarch_variance_path(ret_cc, params, h)
        case ("GJR_GARCH")
            call gjr_variance_path(ret_cc, params, h)
        case ("EGARCH")
            call egarch_variance_path(ret_cc, params, h)
        case default
            allocate(vol(size(ret_cc)))
            call model_vol_forecast(trim(model), ret_cc, params, persist, trading_days_per_year, vol)
            h = max((vol / 100.0_dp)**2 / trading_days_per_year, min_var)
            deallocate(vol)
        end select
    end subroutine cc_garch_variance_path

    pure character(len=16) function cc_garch_row_model(model)
        character(len=*), intent(in) :: model

        select case (trim(model))
        case ("SYMM_GARCH")
            cc_garch_row_model = "GARCH"
        case ("GJR_GARCH")
            cc_garch_row_model = "GJR"
        case default
            cc_garch_row_model = trim(model)
        end select
    end function cc_garch_row_model

    pure character(len=16) function cc_garch_path_model(model)
        character(len=*), intent(in) :: model

        select case (trim(model))
        case ("GARCH")
            cc_garch_path_model = "SYMM_GARCH"
        case ("GJR", "GJR_SIGNED")
            cc_garch_path_model = "GJR_GARCH"
        case ("EWMA_FIT")
            cc_garch_path_model = "EWMA"
        case default
            cc_garch_path_model = trim(model)
        end select
    end function cc_garch_path_model

    pure integer function cc_garch_param_count(model)
        character(len=*), intent(in) :: model

        select case (trim(model))
        case ("EWMA_FIT")
            cc_garch_param_count = 1
        case ("GARCH")
            cc_garch_param_count = 3
        case ("RM2006")
            cc_garch_param_count = 0
        case ("QGARCH", "NAGARCH", "GJR", "GJR_SIGNED", "EGARCH", "HARCH", "TGARCH", "MIDASHYP_ASYM")
            cc_garch_param_count = 4
        case ("CSGARCH", "APARCH", "AVGARCH", "FGTWIST")
            cc_garch_param_count = 5
        case ("MIDASHYP")
            cc_garch_param_count = 3
        case default
            cc_garch_param_count = 0
        end select
    end function cc_garch_param_count

    subroutine fit_and_score_msnagarch(row, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        real(dp), intent(in)  :: ret_cc(:)
        integer,  intent(in)  :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(msnagarch_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_msnagarch(ret_cc, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, "CC", "MSNAGARCH", 11, 0.0_dp, result%persist(1), result%persist(2), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_msnagarch

    subroutine fit_and_score_msgarch(row, ret_cc, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        real(dp), intent(in)  :: ret_cc(:)
        integer,  intent(in)  :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(msgarch_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_msgarch(ret_cc, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, "CC", "MSGARCH", 9, 0.0_dp, result%persist(1), result%persist(2), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_msgarch

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
        call store_full_path(row, h)
    end subroutine fit_and_score_measure

    subroutine fit_and_score_iv_only(row, ret_cc, implied_var, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        real(dp), intent(in) :: ret_cc(:), implied_var(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(affine_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_affine_variance(ret_cc, implied_var, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, "IV", "IV_ONLY", 2, 0.0_dp, result%a, result%b, ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_iv_only

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
        call store_full_path(row, h)
    end subroutine fit_and_score_har

    subroutine fit_and_score_harx_iv(row, measure, ret_cc, x, implied_var, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:), implied_var(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harx_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harx_variance(ret_cc, x, implied_var, train_nobs, result, h)
        call fill_row(row, measure, "HARX_IV", 5, 0.0_dp, result%coef(1), result%coef(5), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_harx_iv

    subroutine fit_and_score_harx_iv_lev(row, measure, ret_cc, x, implied_var, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:), implied_var(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(harx_lev_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_harx_lev_variance(ret_cc, x, implied_var, train_nobs, result, h)
        call fill_row(row, measure, "HARX_IV_LEV", 8, 0.0_dp, result%coef(1), result%coef(5), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_harx_iv_lev

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
        call store_full_path(row, h)
    end subroutine fit_and_score_log_har

    subroutine fit_and_score_sqrt_har(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(sqrt_har_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_sqrt_har_variance(ret_cc, x, train_nobs, result, h)
        call fill_row(row, measure, "SQRT_HAR", 4, 0.0_dp, result%coef(1), result%coef(2), &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_sqrt_har

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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
    end subroutine fit_and_score_har_negret

    subroutine fit_and_score_har_lev(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(har_lev_variance_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_har_lev_variance(ret_cc, x, train_nobs, result, h)
        call fill_row(row, measure, "HAR_LEV", 7, 0.0_dp, result%coef(1), result%coef(5), ret_cc(test_first:nobs), &
                      h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_har_lev

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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
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
        call store_full_path(row, h)
    end subroutine fit_and_score_realgarch

    subroutine fit_and_score_realegarch(row, measure, ret_cc, x, train_nobs, test_first, h)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: measure
        real(dp), intent(in) :: ret_cc(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(realized_egarch_result_t) :: result
        integer :: nobs

        nobs = size(ret_cc)
        call fit_realized_egarch(ret_cc, x, train_nobs, max_iter, gtol, result, h)
        call fill_row(row, measure, "REALEGARCH", 9, 0.0_dp, result%params%omega, result%persist, &
                      ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
        call store_full_path(row, h)
    end subroutine fit_and_score_realegarch

    subroutine append_iv_augmented_rows(rows, irow, ret_cc, implied_var, train_nobs, test_first, h)
        type(score_row_t), intent(inout) :: rows(:)
        integer, intent(inout) :: irow
        real(dp), intent(in) :: ret_cc(:), implied_var(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(affine2_variance_result_t) :: result
        character(len=16) :: model_name
        integer :: base_n, i, nobs

        nobs = size(ret_cc)
        base_n = irow
        do i = 1, base_n
            if (.not. allocated(rows(i)%h_all)) cycle
            if (index(rows(i)%model, "IV") > 0 .or. trim(rows(i)%measure) == "IV") cycle
            if (irow >= size(rows)) error stop "append_iv_augmented_rows: rows array too small"
            call fit_affine2_variance(ret_cc, rows(i)%h_all, implied_var, train_nobs, max_iter, gtol, result, h)
            model_name = iv_augmented_model_name(rows(i)%model)
            irow = irow + 1
            call fill_row(rows(irow), rows(i)%measure, model_name, rows(i)%k + 3, 0.0_dp, result%a, result%b2, &
                          ret_cc(test_first:nobs), h(test_first:nobs), result%niter, result%converged)
            call store_full_path(rows(irow), h)
        end do
    end subroutine append_iv_augmented_rows

    pure character(len=16) function iv_augmented_model_name(model)
        character(len=*), intent(in) :: model

        iv_augmented_model_name = adjustl(trim(model) // "_IV")
    end function iv_augmented_model_name

    subroutine build_lagged_implied_variance(daily_dates, implied_dates, implied_close, adjust_weekday, implied_var)
        integer, intent(in) :: daily_dates(:), implied_dates(:)
        real(dp), intent(in) :: implied_close(:)
        logical, intent(in) :: adjust_weekday
        real(dp), intent(out) :: implied_var(:)
        real(dp), allocatable :: implied_used(:)
        real(dp) :: fallback, vol
        integer :: i, j

        if (size(daily_dates) < size(implied_var)) error stop "build_lagged_implied_variance: too few daily dates"
        if (size(implied_dates) /= size(implied_close)) error stop "build_lagged_implied_variance: implied size mismatch"
        allocate(implied_used(size(implied_close)))
        if (adjust_weekday) then
            call adjust_implied_vol_weekday_levels(implied_dates, implied_close, implied_used)
        else
            implied_used = implied_close
        end if
        if (count(implied_used > 0.0_dp) > 0) then
            fallback = sum(implied_used, mask=implied_used > 0.0_dp) / real(count(implied_used > 0.0_dp), dp)
        else
            fallback = 20.0_dp
        end if

        if (size(implied_dates) < 1) then
            implied_var = max((fallback / 100.0_dp)**2 / trading_days_per_year, min_var)
            deallocate(implied_used)
            return
        end if
        j = 1
        do i = 1, size(implied_var)
            do while (j < size(implied_dates))
                if (implied_dates(j + 1) > daily_dates(i)) exit
                j = j + 1
            end do
            vol = fallback
            if (implied_dates(j) <= daily_dates(i)) then
                if (implied_used(j) > 0.0_dp) vol = implied_used(j)
            end if
            implied_var(i) = max((vol / 100.0_dp)**2 / trading_days_per_year, min_var)
        end do
        deallocate(implied_used)
    end subroutine build_lagged_implied_variance

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

    subroutine store_full_path(row, h)
        type(score_row_t), intent(inout) :: row
        real(dp), intent(in) :: h(:)

        if (allocated(row%h_all)) deallocate(row%h_all)
        allocate(row%h_all(size(h)))
        row%h_all = h
    end subroutine store_full_path

    subroutine print_implied_vol_correlation_table(forecast_dates, rows, implied_vol_dates, implied_vol_close, &
                                                   implied_vol_label, adjust_weekday)
        integer, intent(in) :: forecast_dates(:), implied_vol_dates(:)
        type(score_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: implied_vol_close(:)
        character(len=*), intent(in) :: implied_vol_label
        logical, intent(in) :: adjust_weekday
        real(dp), allocatable :: model_vol(:), implied_vol(:), dlog_model(:), dlog_implied_vol(:)
        real(dp), allocatable :: model_vol_mat(:, :)
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

        allocate(model_vol(nmatch), implied_vol(nmatch), model_vol_mat(nmatch, size(rows)))
        allocate(corr_level(size(rows)), corr_dlog(size(rows)), corr_rank(size(rows)), dlog_rank(size(rows)))
        corr_level = -huge(1.0_dp)
        corr_dlog = -huge(1.0_dp)
        do j = 1, nmatch
            implied_vol(j) = max(implied_vol_used(implied_vol_idx(j)), 0.0_dp)
        end do
        do i = 1, size(rows)
            do j = 1, nmatch
                if (implied_vol_used(implied_vol_idx(j)) > 0.0_dp .and. rows(i)%h_test(forecast_idx(j)) > 0.0_dp) then
                    model_vol(j) = 100.0_dp*sqrt(trading_days_per_year*rows(i)%h_test(forecast_idx(j)))
                else
                    model_vol(j) = 0.0_dp
                end if
            end do
            model_vol_mat(:, i) = model_vol
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
        if (fit_implied_vol_nnls) then
            call print_implied_vol_nnls_table(rows, model_vol_mat, implied_vol, implied_vol_label)
        end if
    end subroutine print_implied_vol_correlation_table

    subroutine print_implied_vol_nnls_table(rows, model_vol_mat, implied_vol, implied_vol_label)
        type(score_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: model_vol_mat(:, :), implied_vol(:)
        character(len=*), intent(in) :: implied_vol_label
        real(dp), allocatable :: beta(:), fitted(:), xreg(:, :), tstat(:), pred_corr(:)
        real(dp), allocatable :: retained_x(:, :)
        character(len=32), allocatable :: retained_labels(:)
        real(dp) :: sse, rmse, corr_fit
        integer, allocatable :: order(:)
        integer :: i, ii, k, ncoef, nmodel, niter, nnz, offset

        nmodel = size(rows)
        ncoef = nmodel + merge(1, 0, implied_vol_nnls_intercept)
        allocate(beta(ncoef), fitted(size(implied_vol)), xreg(size(implied_vol), ncoef), tstat(ncoef), pred_corr(ncoef))
        if (implied_vol_nnls_intercept) then
            xreg(:, 1) = 1.0_dp
            xreg(:, 2:ncoef) = model_vol_mat
        else
            xreg = model_vol_mat
        end if

        call nonnegative_least_squares(xreg, implied_vol, beta, sse, implied_vol_nnls_max_iter, &
                                       implied_vol_nnls_tol, niter)
        fitted = matmul(xreg, beta)
        rmse = sqrt(sse / real(size(implied_vol), dp))
        corr_fit = correlation(fitted, implied_vol)
        nnz = count(beta > 1.0e-8_dp)
        call active_nnls_tstats(xreg, implied_vol, beta, sse, tstat)
        call predictor_correlations(xreg, implied_vol, pred_corr)
        offset = merge(1, 0, implied_vol_nnls_intercept)

        print '(A)', ""
        print '(A,A,A)', "Non-negative regression of ", trim(implied_vol_label), " close on model forecasts"
        print '(A,I0,A,I0,A,F8.4,A,F8.4)', "models: ", nmodel, "  nonzero: ", nnz, "  corr: ", corr_fit, &
              "  rmse: ", rmse
        print '(A,I0,A,ES10.3)', "iterations: ", niter, "  sse: ", sse
        print '(A)', "---------------------------------------------------------------------------------------"
        print '(A16,1X,A16,1X,A12,1X,A10,1X,A10)', "Measure", "Model", "coef", "t_stat", "corr"
        print '(A)', "---------------------------------------------------------------------------------------"
        if (implied_vol_nnls_intercept .and. beta(1) > 1.0e-8_dp) then
            print '(A16,1X,A16,1X,F12.6,1X,F10.3,1X,F10.4)', "-", "INTERCEPT", beta(1), tstat(1), pred_corr(1)
        end if
        allocate(order(nmodel))
        call rank_nonzero_by_tstat(beta(offset + 1:offset + nmodel), tstat(offset + 1:offset + nmodel), order)
        do ii = 1, nmodel
            i = order(ii)
            if (i < 1) cycle
            if (beta(i + offset) > 1.0e-8_dp) then
                print '(A16,1X,A16,1X,F12.6,1X,F10.3,1X,F10.4)', rows(i)%measure, rows(i)%model, &
                      beta(i + offset), tstat(i + offset), pred_corr(i + offset)
            end if
        end do
        print '(A)', "---------------------------------------------------------------------------------------"
        if (print_implied_vol_nnls_corr_matrix) then
            k = count(beta(offset + 1:offset + nmodel) > 1.0e-8_dp)
            if (k > 0) then
                allocate(retained_x(size(implied_vol), k), retained_labels(k))
                k = 0
                do ii = 1, nmodel
                    i = order(ii)
                    if (i < 1) cycle
                    if (beta(i + offset) <= 1.0e-8_dp) cycle
                    k = k + 1
                    retained_x(:, k) = model_vol_mat(:, i)
                    retained_labels(k) = predictor_matrix_label(rows(i))
                end do
                call print_response_predictor_corr_matrix("Correlation matrix of response and retained NNLS predictors", &
                                                          implied_vol_label, retained_labels, implied_vol, retained_x)
                deallocate(retained_x, retained_labels)
            end if
        end if
        deallocate(beta, fitted, xreg, tstat, pred_corr, order)
    end subroutine print_implied_vol_nnls_table

    pure character(len=32) function predictor_matrix_label(row)
        type(score_row_t), intent(in) :: row

        predictor_matrix_label = adjustl(trim(row%measure) // ":" // trim(row%model))
    end function predictor_matrix_label

    subroutine rank_rows(rows)
        type(score_row_t), intent(inout) :: rows(:)
        integer :: i

        do i = 1, size(rows)
            rows(i)%aic_rank = 1 + count(rows%aic < rows(i)%aic)
            rows(i)%bic_rank = 1 + count(rows%bic < rows(i)%bic)
        end do
    end subroutine rank_rows

    subroutine print_summary(intraday_prices_file, implied_vol_file, daily, train_nobs, ntest_days, rows, read_sec, fit_sec, elapsed_sec)
        character(len=*), intent(in) :: intraday_prices_file, implied_vol_file
        type(daily_realized_panel_t), intent(in) :: daily
        integer, intent(in) :: train_nobs, ntest_days
        type(score_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: read_sec, fit_sec, elapsed_sec
        integer :: i, best

        best = maxloc(rows%loglik, dim=1)
        print '(A)', "Daily CC volatility forecasts from intraday realized measures"
        print '(A,A)', "Intraday prices file: ", trim(intraday_prices_file)
        print '(A,A)', "Implied vol file:     ", trim(implied_vol_file)
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

! Fit close-close and OHLC GARCH-family Normal models to ETF returns and
! compare model volatility forecasts with mapped implied-volatility indices.
! Edit models(:), fit_assets(:), asset_iv_assets(:), and asset_iv_indices(:) to configure.

program xfit_garch_ohlc_iv_returns
    use kind_mod,       only: dp
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv, read_ohlc_csv, print_price_sample_info
    use stats_mod,      only: mean, sd
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use vol_forecast_compare_mod, only: print_implied_vol_correlations
    use garch_forecast_mod, only: model_vol_forecast, vol_forecast_stats_t, summarize_vol_forecast, &
                                  print_vol_forecast_table
    use model_selection_mod, only: print_model_selection_counts, print_model_fit_times
    use garch_fit_mod,  only: fit_symm_garch, fit_qgarch, fit_figarch, fit_nagarch, fit_rgarch, fit_regarch1, fit_regarch2, &
                              fit_rgarch_meas, fit_carr_park, fit_nagarch_range, fit_fgarch_twist_range, &
                              fit_gjr, fit_fgarch_twist, fit_ewma, &
                              fit_aewma_nag, fit_aewma_twist, &
                              garch_skew_kurt, qgarch_skew_kurt, figarch_skew_kurt, nagarch_skew_kurt, gjr_skew_kurt, &
                              rgarch_skew_kurt, regarch1_skew_kurt, regarch2_skew_kurt, &
                              rgarch_meas_skew_kurt, carr_park_skew_kurt, &
                              nagarch_range_skew_kurt, fgarch_twist_range_skew_kurt, &
                              ewma_skew_kurt, aewma_nag_skew_kurt, aewma_twist_skew_kurt, &
                              symm_garch_persist, qgarch_persist, figarch_persist, nagarch_persist, gjr_persist, &
                              rgarch_persist, regarch1_persist, regarch2_persist, rgarch_meas_persist, &
                              carr_park_persist, fgarch_twist_persist, ewma_persist, &
                              aewma_nag_persist, aewma_twist_persist, qgarch_mean_variance, figarch_variance
    implicit none

    character(len=*), parameter :: prices_file = "prices_ohlc.csv"
    character(len=*), parameter :: implied_vol_file = "vix_spy.csv"
    integer, parameter :: symbol_len = 16
    integer, parameter :: model_len = 16
    real(dp), parameter :: trading_days = 252.0_dp
    real(dp), parameter :: co_oc_corr = 0.0_dp
    integer,  parameter :: max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp
    logical,  parameter :: print_vol_forecast_stats = .false.
    logical,  parameter :: fit_all_assets = .true.
    logical,  parameter :: flip_log_return_sign = .true.
    logical,  parameter :: standardize_log_return_by_prior_iv = .true.
    character(len=model_len), parameter :: models(*) = [character(len=model_len) :: &
        "EWMA", "AEWMA_NAG", "AEWMA_TWIST", "SYMM_GARCH", "QGARCH", "FIGARCH", "NAGARCH", "GJR_GARCH", "FGTWIST", "REGARCH1", &
        "REGARCH2", "CARR_PARK", "CARR_GK", "OHLC_SYMM", "OHLC_NAGARCH", "OHLC_RGARCH", "OHLC_NAG_RANGE", &
        "OHLC_FGTW_RANGE", "OHLC_GJR", "OHLC_FIGARCH", "OHLC_FGTWIST"] ! , "RGARCH_MEAS" 
    character(len=symbol_len), parameter :: fit_assets(*) = [character(len=symbol_len) :: "SPY"]
    character(len=symbol_len), parameter :: asset_iv_assets(*) = [character(len=symbol_len) :: "SPY"]
    character(len=symbol_len), parameter :: asset_iv_indices(*) = [character(len=symbol_len) :: "VIX"]
    character(len=model_len), parameter :: extra_series_names(*) = [character(len=model_len) :: "LOG_RETURN"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:), iv_dates(:)
    character(len=32), allocatable :: col_names(:), iv_col_names(:)
    real(dp), allocatable :: prices(:,:), open_prices(:,:), high_prices(:,:), low_prices(:,:), iv_values(:,:)
    real(dp), allocatable :: ret(:), ret_co(:), ret_oc(:), range_var(:), gk_var(:), log_range(:), variance_tmp(:)
    type(garch_fit_result_t) :: rows(n_model)
    integer :: row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    logical :: row_comparable(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: conv_count(n_model), conv_total(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    integer :: fit_count(n_model), fit_clock_start, fit_clock_end
    character(len=256) :: aic_symbols(n_model)
    character(len=256) :: failed_symbols(n_model)
    character(len=model_len) :: fit_model_names(n_model)
    character(len=model_len), allocatable :: vf_model(:)
    integer :: nprices, ncols, nobs, icol, imod, ifit, jfit, nfit
    integer :: niter, np, clock_start, clock_end, clock_rate
    real(dp) :: ret_mean, ret_std, fopt
    real(dp) :: persist, h_unc, vol_ann, logl, aic, bic, skew, ekurt, elapsed_s
    real(dp) :: fit_seconds(n_model), sec_per_asset_arr(n_model)
    real(dp), allocatable :: vol_forecast(:), vol_forecasts(:,:,:)
    real(dp), allocatable :: extra_series(:,:,:)
    real(dp) :: extra_signs(size(extra_series_names))
    type(vol_forecast_stats_t), allocatable :: vf_stats(:,:)
    type(garch_params_t) :: params
    logical :: converged
    logical :: fit_model_have(n_model)
    logical, allocatable :: vf_have(:,:)
    character(len=model_len) :: model

    call system_clock(clock_start, clock_rate)
    call nagarch_set_news_impact(.false.)
    aic_wins = 0
    bic_wins = 0
    aic_rank_sum = 0
    bic_rank_sum = 0
    rank_count = 0
    param_count = 0
    fit_count = 0
    conv_count = 0
    conv_total = 0
    fit_seconds = 0.0_dp
    aic_symbols = ""
    failed_symbols = ""
    fit_model_have = .false.
    fit_model_names = ""

    if (fit_all_assets) then
        call read_ohlc_csv(prices_file, dates, col_names, open_prices, prices, high_prices=high_prices, low_prices=low_prices)
    else
        call read_ohlc_csv(prices_file, dates, col_names, open_prices, prices, fit_assets, high_prices, low_prices)
    end if
    call read_price_csv(implied_vol_file, iv_dates, iv_col_names, iv_values)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), ret_co(nobs), ret_oc(nobs), range_var(nobs), gk_var(nobs), log_range(nobs), variance_tmp(nobs), vol_forecast(nobs), &
             vol_forecasts(nobs,ncols,n_model))
    allocate(extra_series(nprices,ncols,size(extra_series_names)))
    vol_forecasts = 0.0_dp
    extra_series(:,:,1) = prices
    extra_signs = 1.0_dp
    if (flip_log_return_sign) extra_signs(1) = -1.0_dp
    if (print_vol_forecast_stats) then
        allocate(vf_stats(ncols,n_model), vf_have(ncols,n_model))
        allocate(vf_model(n_model))
        vf_have = .false.
        vf_model = ""
    end if

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A,A)') "Implied vol file: ", trim(implied_vol_file)
    write(*,'(A,F7.3)') "OHLC CO/OC variance forecast correlation: ", co_oc_corr
    write(*,'(A)') "OHLC model logL/AIC/BIC use separate CO and OC likelihoods; compare those criteria to CC rows with care."
    write(*,'(A)') "REGARCH1/REGARCH2 logL/AIC/BIC use the log-range likelihood and are excluded from AIC/BIC ranks and wins."
    write(*,'(A)') "RGARCH_MEAS logL/AIC/BIC use joint return/log-range likelihood and are excluded from AIC/BIC ranks and wins."
    write(*,'(A)') "CARR_PARK/CARR_GK logL/AIC/BIC use range-estimator likelihoods and are excluded from AIC/BIC ranks and wins."
    call print_fit_assets()
    write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist   scale  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
    write(*,'(A)') repeat("-", 180)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        ret_std = sd(ret)
        ret_co = log(open_prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_co = ret_co - mean(ret_co)
        ret_oc = log(prices(2:nprices,icol) / open_prices(2:nprices,icol))
        ret_oc = ret_oc - mean(ret_oc)
        range_var = (log(high_prices(2:nprices,icol) / low_prices(2:nprices,icol)))**2 / (4.0_dp*log(2.0_dp))
        gk_var = max(0.5_dp*(log(high_prices(2:nprices,icol) / low_prices(2:nprices,icol)))**2 - &
                     (2.0_dp*log(2.0_dp) - 1.0_dp)* &
                     (log(prices(2:nprices,icol) / open_prices(2:nprices,icol)))**2, 1.0e-12_dp)
        log_range = log(max(log(high_prices(2:nprices,icol) / low_prices(2:nprices,icol)), 1.0e-12_dp))
        nfit = 0

        do imod = 1, size(models)
            model = uppercase(trim(models(imod)))
            params = garch_params_t()
            call system_clock(fit_clock_start)

            select case (trim(model))
            case ("EWMA", "IGARCH", "IGARCH_EWMA")
                call fit_ewma(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = ewma_persist(params)
                h_unc = sum(ret**2) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 1
                call ewma_skew_kurt(ret, params, skew, ekurt)
                model = "EWMA"
            case ("AEWMA_NAG", "ASYM_EWMA_NAG")
                call fit_aewma_nag(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = aewma_nag_persist(params)
                h_unc = sum(ret**2) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 3
                call aewma_nag_skew_kurt(ret, params, skew, ekurt)
                model = "AEWMA_NAG"
            case ("AEWMA_TWIST", "ASYM_EWMA_TWIST")
                call fit_aewma_twist(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = aewma_twist_persist(params)
                h_unc = sum(ret**2) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call aewma_twist_skew_kurt(ret, params, skew, ekurt)
                model = "AEWMA_TWIST"
            case ("SYMM_GARCH", "GARCH", "SYMM")
                call fit_symm_garch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = symm_garch_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 3
                call garch_skew_kurt(ret, params, skew, ekurt)
                model = "SYMM_GARCH"
            case ("QGARCH", "QUADRATIC_GARCH")
                call fit_qgarch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = qgarch_persist(params)
                h_unc = qgarch_mean_variance(params)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call qgarch_skew_kurt(ret, params, skew, ekurt)
                model = "QGARCH"
            case ("FIGARCH")
                call fit_figarch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = figarch_persist(params)
                call figarch_variance(ret, params, variance_tmp)
                h_unc = sum(variance_tmp) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call figarch_skew_kurt(ret, params, skew, ekurt)
            case ("NAGARCH")
                call fit_nagarch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = nagarch_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call nagarch_skew_kurt(ret, params, skew, ekurt)
            case ("GJR_GARCH", "GJR")
                call fit_gjr(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = gjr_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call gjr_skew_kurt(ret, params, skew, ekurt)
                model = "GJR_GARCH"
            case ("FGTWIST", "FGARCH_TWIST")
                call fit_fgarch_twist(ret, ret_std, max_iter, gtol, fopt, params, &
                                      vol_ann, skew, ekurt, niter, converged)
                persist = fgarch_twist_persist(params)
                np = 5
                model = "FGTWIST"
            case ("REGARCH1", "RANGE_EGARCH")
                call fit_regarch1(ret, log_range, max_iter, gtol, fopt, params, niter, converged)
                persist = regarch1_persist(params)
                vol_ann = regarch1_vol_ann(ret, log_range, params)
                np = 4
                call regarch1_skew_kurt(ret, log_range, params, skew, ekurt)
                model = "REGARCH1"
            case ("REGARCH2", "RANGE_EGARCH2")
                call fit_regarch2(ret, log_range, max_iter, gtol, fopt, params, niter, converged)
                persist = regarch2_persist(params)
                vol_ann = regarch2_vol_ann(ret, log_range, params)
                np = 7
                call regarch2_skew_kurt(ret, log_range, params, skew, ekurt)
                model = "REGARCH2"
            case ("RGARCH_MEAS", "REALGARCH_RANGE")
                call fit_rgarch_meas(ret, log_range, max_iter, gtol, fopt, params, niter, converged)
                persist = rgarch_meas_persist(params)
                vol_ann = rgarch_meas_vol_ann(ret, log_range, params)
                np = 9
                call rgarch_meas_skew_kurt(ret, log_range, params, skew, ekurt)
                model = "RGARCH_MEAS"
            case ("CARR_PARK", "GARCH_PARK_R", "CARR")
                call fit_carr_park(range_var, max_iter, gtol, fopt, params, niter, converged)
                persist = carr_park_persist(params)
                vol_ann = carr_park_vol_ann(range_var, params)
                np = 3
                call carr_park_skew_kurt(ret, range_var, params, skew, ekurt)
                model = "CARR_PARK"
            case ("CARR_GK", "GARCH_GK_R")
                call fit_carr_park(gk_var, max_iter, gtol, fopt, params, niter, converged)
                persist = carr_park_persist(params)
                vol_ann = carr_park_vol_ann(gk_var, params)
                np = 3
                call carr_park_skew_kurt(ret, gk_var, params, skew, ekurt)
                model = "CARR_GK"
            case ("OHLC_SYMM", "OHLC_SYMM_GARCH")
                call fit_ohlc_model("SYMM_GARCH", ret_co, ret_oc, fopt, params, vol_ann, persist, &
                                    skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_SYMM"
            case ("OHLC_NAGARCH")
                call fit_ohlc_model("NAGARCH", ret_co, ret_oc, fopt, params, vol_ann, persist, &
                                    skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_NAGARCH"
            case ("OHLC_RGARCH", "OHLC_RANGE_GARCH")
                call fit_ohlc_rgarch_model(ret_co, ret_oc, range_var, fopt, params, vol_ann, persist, &
                                           skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_RGARCH"
            case ("OHLC_NAG_RANGE", "OHLC_NAGARCH_RANGE")
                call fit_ohlc_nag_range_model(ret_co, ret_oc, range_var, fopt, params, vol_ann, persist, &
                                              skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_NAG_RANGE"
            case ("OHLC_FGTW_RANGE", "OHLC_FGTWIST_RANGE", "OHLC_FG_RANGE")
                call fit_ohlc_fgarch_twist_range_model(ret_co, ret_oc, range_var, fopt, params, vol_ann, persist, &
                                                       skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_FGTW_RANGE"
            case ("OHLC_GJR", "OHLC_GJR_GARCH")
                call fit_ohlc_model("GJR_GARCH", ret_co, ret_oc, fopt, params, vol_ann, persist, &
                                    skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_GJR"
            case ("OHLC_FIGARCH")
                call fit_ohlc_model("FIGARCH", ret_co, ret_oc, fopt, params, vol_ann, persist, &
                                    skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_FIGARCH"
            case ("OHLC_FGTWIST")
                call fit_ohlc_model("FGTWIST", ret_co, ret_oc, fopt, params, vol_ann, persist, &
                                    skew, ekurt, niter, converged, np, vol_forecast)
                model = "OHLC_FGTWIST"
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select

            call system_clock(fit_clock_end)
            fit_seconds(imod) = fit_seconds(imod) + real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)
            fit_count(imod) = fit_count(imod) + 1
            conv_total(imod) = conv_total(imod) + 1
            if (converged) then
                conv_count(imod) = conv_count(imod) + 1
            else if (len_trim(failed_symbols(imod)) == 0) then
                failed_symbols(imod) = trim(col_names(icol))
            else
                failed_symbols(imod) = trim(failed_symbols(imod)) // " " // trim(col_names(icol))
            end if

            logl = -real(nobs, dp) * fopt
            aic = 2.0_dp * real(np, dp) - 2.0_dp * logl
            bic = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl
            if (trim(model) == "REGARCH1") then
                call regarch1_vol_forecast(ret, log_range, params, vol_forecast)
            else if (trim(model) == "REGARCH2") then
                call regarch2_vol_forecast(ret, log_range, params, vol_forecast)
            else if (trim(model) == "RGARCH_MEAS") then
                call rgarch_meas_vol_forecast(ret, log_range, params, vol_forecast)
            else if (trim(model) == "CARR_PARK") then
                call carr_park_vol_forecast(range_var, params, vol_forecast)
            else if (trim(model) == "CARR_GK") then
                call carr_park_vol_forecast(gk_var, params, vol_forecast)
            else if (index(trim(model), "OHLC_") /= 1) then
                call model_vol_forecast(trim(model), ret, params, persist, trading_days, vol_forecast)
            end if
            vol_forecasts(:,icol,imod) = vol_forecast
            fit_model_names(imod) = model
            fit_model_have(imod) = .true.
            if (print_vol_forecast_stats) then
                call summarize_vol_forecast(vol_forecast, dates(2:nprices), vf_stats(icol,imod))
                vf_model(imod) = model
                vf_have(icol,imod) = .true.
            end if

            nfit = nfit + 1
            rows(nfit) = garch_fit_result_t()
            rows(nfit)%model = model
            row_model_idx(nfit) = imod
            rows(nfit)%params = params
            rows(nfit)%persist = persist
            rows(nfit)%vol_ann = vol_ann
            rows(nfit)%logl = logl
            rows(nfit)%aic = aic
            rows(nfit)%bic = bic
            rows(nfit)%nparam = np
            rows(nfit)%niter = niter
            rows(nfit)%converged = converged
            rows(nfit)%skew = skew
            rows(nfit)%ekurt = ekurt
            row_comparable(nfit) = trim(model) /= "REGARCH1" .and. trim(model) /= "REGARCH2" .and. &
                                   trim(model) /= "RGARCH_MEAS" .and. trim(model) /= "CARR_PARK" .and. &
                                   trim(model) /= "CARR_GK"
            param_count(imod) = np
        end do

        do ifit = 1, nfit
            row_aic_rank(ifit) = 1
            row_bic_rank(ifit) = 1
            if (.not. row_comparable(ifit)) then
                row_aic_rank(ifit) = 0
                row_bic_rank(ifit) = 0
                cycle
            end if
            do jfit = 1, nfit
                if (.not. row_comparable(jfit)) cycle
                if (rows(jfit)%aic < rows(ifit)%aic) row_aic_rank(ifit) = row_aic_rank(ifit) + 1
                if (rows(jfit)%bic < rows(ifit)%bic) row_bic_rank(ifit) = row_bic_rank(ifit) + 1
            end do
            if (row_aic_rank(ifit) == 1) then
                aic_wins(row_model_idx(ifit)) = aic_wins(row_model_idx(ifit)) + 1
                if (len_trim(aic_symbols(row_model_idx(ifit))) == 0) then
                    aic_symbols(row_model_idx(ifit)) = trim(col_names(icol))
                else
                    aic_symbols(row_model_idx(ifit)) = trim(aic_symbols(row_model_idx(ifit))) // " " // trim(col_names(icol))
                end if
            end if
            if (row_bic_rank(ifit) == 1) bic_wins(row_model_idx(ifit)) = bic_wins(row_model_idx(ifit)) + 1
            aic_rank_sum(row_model_idx(ifit)) = aic_rank_sum(row_model_idx(ifit)) + row_aic_rank(ifit)
            bic_rank_sum(row_model_idx(ifit)) = bic_rank_sum(row_model_idx(ifit)) + row_bic_rank(ifit)
            rank_count(row_model_idx(ifit)) = rank_count(row_model_idx(ifit)) + 1
        end do

        do ifit = 1, nfit
            write(*,'(A16,1X,A9,ES12.3,7F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F10.3,2I9)') &
                trim(rows(ifit)%model), trim(col_names(icol)), rows(ifit)%params%omega, rows(ifit)%params%alpha, &
                rows(ifit)%params%gamma, rows(ifit)%params%beta, rows(ifit)%params%theta, &
                rows(ifit)%params%twist, rows(ifit)%params%scale, rows(ifit)%persist, &
                rows(ifit)%vol_ann, rows(ifit)%logl, rows(ifit)%aic, rows(ifit)%bic, rows(ifit)%niter, &
                rows(ifit)%converged, rows(ifit)%skew, rows(ifit)%ekurt, row_aic_rank(ifit), row_bic_rank(ifit)
        end do
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    do imod = 1, size(models)
        sec_per_asset_arr(imod) = sec_per_asset(imod)
    end do
    call print_model_selection_counts(models, param_count, aic_wins, bic_wins, aic_symbols, &
                                      aic_rank_sum, bic_rank_sum, rank_count, sec_per_asset_arr)
    call print_model_fit_times(models, param_count, fit_seconds)
    call print_convergence_failures()
    if (print_vol_forecast_stats) call print_vol_forecast_table(col_names, vf_model, vf_stats, vf_have)
    call print_implied_vol_correlations(implied_vol_file, col_names, fit_model_names, fit_model_have, &
                                        dates(2:nprices), vol_forecasts, iv_dates, iv_col_names, iv_values, &
                                        asset_iv_assets, asset_iv_indices)
    call print_implied_vol_correlations(implied_vol_file, col_names, fit_model_names, fit_model_have, &
                                        dates(2:nprices), vol_forecasts, iv_dates, iv_col_names, iv_values, &
                                        asset_iv_assets, asset_iv_indices, "log_diff", [1, 5, 22], &
                                        extra_series_names, dates, extra_series, extra_signs, &
                                        standardize_log_return_by_prior_iv)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    real(dp) function sec_per_asset(imod_in)
        integer, intent(in) :: imod_in

        if (fit_count(imod_in) > 0) then
            sec_per_asset = fit_seconds(imod_in) / real(fit_count(imod_in), dp)
        else
            sec_per_asset = 0.0_dp
        end if
    end function sec_per_asset

    subroutine print_convergence_failures()
        integer :: i
        logical :: any_failed

        any_failed = .false.
        do i = 1, size(models)
            if (conv_total(i) > 0 .and. conv_count(i) < conv_total(i)) any_failed = .true.
        end do
        if (.not. any_failed) return

        write(*,'(/,A)') "Convergence failures:"
        write(*,'(A)') "Model              conv  total  fit_sec  failed_symbols"
        write(*,'(A)') repeat("-", 68)
        do i = 1, size(models)
            if (conv_total(i) == 0) cycle
            if (conv_count(i) == conv_total(i)) cycle
            write(*,'(A16,2I7,F9.3,2X,A)') trim(uppercase(trim(models(i)))), conv_count(i), conv_total(i), &
                fit_seconds(i), &
                trim(failed_symbols(i))
        end do
    end subroutine print_convergence_failures

    subroutine fit_ohlc_rgarch_model(y_co, y_oc, x_range, f_combo, params_out, vol_ann_out, persist_out, &
                                     skew_out, ekurt_out, niter_out, converged_out, np_out, vol_out)
        real(dp), intent(in) :: y_co(:), y_oc(:), x_range(:)
        real(dp), intent(out) :: f_combo, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out, np_out
        logical, intent(out) :: converged_out
        real(dp), intent(out) :: vol_out(:)
        type(garch_params_t) :: p_co, p_oc
        real(dp), allocatable :: vol_co(:), vol_oc(:)
        real(dp) :: f_co, f_oc, persist_co, persist_oc, vol_ann_co, vol_ann_oc
        real(dp) :: skew_co, skew_oc, ekurt_co, ekurt_oc
        integer :: niter_co, niter_oc, t
        logical :: conv_co, conv_oc

        allocate(vol_co(size(y_co)), vol_oc(size(y_oc)))
        call fit_one_rgarch(y_co, x_range, f_co, p_co, vol_ann_co, persist_co, &
                            skew_co, ekurt_co, niter_co, conv_co)
        call fit_one_rgarch(y_oc, x_range, f_oc, p_oc, vol_ann_oc, persist_oc, &
                            skew_oc, ekurt_oc, niter_oc, conv_oc)
        call rgarch_vol_forecast(y_co, x_range, p_co, vol_co)
        call rgarch_vol_forecast(y_oc, x_range, p_oc, vol_oc)
        do t = 1, size(vol_out)
            vol_out(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do

        f_combo = f_co + f_oc
        params_out = p_co
        persist_out = 0.5_dp * (persist_co + persist_oc)
        vol_ann_out = sqrt(sum(vol_out**2) / real(size(vol_out), dp))
        skew_out = 0.5_dp * (skew_co + skew_oc)
        ekurt_out = 0.5_dp * (ekurt_co + ekurt_oc)
        niter_out = niter_co + niter_oc
        converged_out = conv_co .and. conv_oc
        np_out = 6
        deallocate(vol_co, vol_oc)
    end subroutine fit_ohlc_rgarch_model

    subroutine fit_one_rgarch(y, x_range, f_out, params_out, vol_ann_out, persist_out, &
                              skew_out, ekurt_out, niter_out, converged_out)
        real(dp), intent(in) :: y(:), x_range(:)
        real(dp), intent(out) :: f_out, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out
        logical, intent(out) :: converged_out
        real(dp) :: h0

        call fit_rgarch(y, x_range, 2*max_iter, gtol, f_out, params_out, niter_out, converged_out)
        persist_out = rgarch_persist(params_out)
        h0 = rgarch_h0(y, x_range, params_out)
        vol_ann_out = sqrt(trading_days * h0) * 100.0_dp
        call rgarch_skew_kurt(y, x_range, params_out, skew_out, ekurt_out)
    end subroutine fit_one_rgarch

    subroutine fit_ohlc_nag_range_model(y_co, y_oc, x_range, f_combo, params_out, vol_ann_out, persist_out, &
                                        skew_out, ekurt_out, niter_out, converged_out, np_out, vol_out)
        real(dp), intent(in) :: y_co(:), y_oc(:), x_range(:)
        real(dp), intent(out) :: f_combo, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out, np_out
        logical, intent(out) :: converged_out
        real(dp), intent(out) :: vol_out(:)
        type(garch_params_t) :: p_co, p_oc
        real(dp), allocatable :: vol_co(:), vol_oc(:)
        real(dp) :: f_co, f_oc, persist_co, persist_oc, vol_ann_co, vol_ann_oc
        real(dp) :: skew_co, skew_oc, ekurt_co, ekurt_oc
        integer :: niter_co, niter_oc, t
        logical :: conv_co, conv_oc

        allocate(vol_co(size(y_co)), vol_oc(size(y_oc)))
        call fit_one_nag_range(y_co, x_range, f_co, p_co, vol_ann_co, persist_co, &
                               skew_co, ekurt_co, niter_co, conv_co)
        call fit_one_nag_range(y_oc, x_range, f_oc, p_oc, vol_ann_oc, persist_oc, &
                               skew_oc, ekurt_oc, niter_oc, conv_oc)
        call nagarch_range_vol_forecast(y_co, x_range, p_co, persist_co, vol_co)
        call nagarch_range_vol_forecast(y_oc, x_range, p_oc, persist_oc, vol_oc)
        do t = 1, size(vol_out)
            vol_out(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do

        f_combo = f_co + f_oc
        params_out = p_co
        persist_out = 0.5_dp * (persist_co + persist_oc)
        vol_ann_out = sqrt(sum(vol_out**2) / real(size(vol_out), dp))
        skew_out = 0.5_dp * (skew_co + skew_oc)
        ekurt_out = 0.5_dp * (ekurt_co + ekurt_oc)
        niter_out = niter_co + niter_oc
        converged_out = conv_co .and. conv_oc
        np_out = 10
        deallocate(vol_co, vol_oc)
    end subroutine fit_ohlc_nag_range_model

    subroutine fit_one_nag_range(y, x_range, f_out, params_out, vol_ann_out, persist_out, &
                                 skew_out, ekurt_out, niter_out, converged_out)
        real(dp), intent(in) :: y(:), x_range(:)
        real(dp), intent(out) :: f_out, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out
        logical, intent(out) :: converged_out
        real(dp) :: h0

        call fit_nagarch_range(y, x_range, 2*max_iter, gtol, f_out, params_out, niter_out, converged_out)
        persist_out = nagarch_persist(params_out)
        h0 = nagarch_range_h0(y, x_range, params_out, persist_out)
        vol_ann_out = sqrt(trading_days * h0) * 100.0_dp
        call nagarch_range_skew_kurt(y, x_range, params_out, skew_out, ekurt_out)
    end subroutine fit_one_nag_range

    subroutine fit_ohlc_fgarch_twist_range_model(y_co, y_oc, x_range, f_combo, params_out, vol_ann_out, persist_out, &
                                                 skew_out, ekurt_out, niter_out, converged_out, np_out, vol_out)
        real(dp), intent(in) :: y_co(:), y_oc(:), x_range(:)
        real(dp), intent(out) :: f_combo, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out, np_out
        logical, intent(out) :: converged_out
        real(dp), intent(out) :: vol_out(:)
        type(garch_params_t) :: p_co, p_oc
        real(dp), allocatable :: vol_co(:), vol_oc(:)
        real(dp) :: f_co, f_oc, persist_co, persist_oc, vol_ann_co, vol_ann_oc
        real(dp) :: skew_co, skew_oc, ekurt_co, ekurt_oc
        integer :: niter_co, niter_oc, t
        logical :: conv_co, conv_oc

        allocate(vol_co(size(y_co)), vol_oc(size(y_oc)))
        call fit_one_fgarch_twist_range(y_co, x_range, f_co, p_co, vol_ann_co, persist_co, &
                                        skew_co, ekurt_co, niter_co, conv_co)
        call fit_one_fgarch_twist_range(y_oc, x_range, f_oc, p_oc, vol_ann_oc, persist_oc, &
                                        skew_oc, ekurt_oc, niter_oc, conv_oc)
        call fgarch_twist_range_vol_forecast(y_co, x_range, p_co, persist_co, vol_co)
        call fgarch_twist_range_vol_forecast(y_oc, x_range, p_oc, persist_oc, vol_oc)
        do t = 1, size(vol_out)
            vol_out(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do

        f_combo = f_co + f_oc
        params_out = p_co
        persist_out = 0.5_dp * (persist_co + persist_oc)
        vol_ann_out = sqrt(sum(vol_out**2) / real(size(vol_out), dp))
        skew_out = 0.5_dp * (skew_co + skew_oc)
        ekurt_out = 0.5_dp * (ekurt_co + ekurt_oc)
        niter_out = niter_co + niter_oc
        converged_out = conv_co .and. conv_oc
        np_out = 12
        deallocate(vol_co, vol_oc)
    end subroutine fit_ohlc_fgarch_twist_range_model

    subroutine fit_one_fgarch_twist_range(y, x_range, f_out, params_out, vol_ann_out, persist_out, &
                                          skew_out, ekurt_out, niter_out, converged_out)
        real(dp), intent(in) :: y(:), x_range(:)
        real(dp), intent(out) :: f_out, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out
        logical, intent(out) :: converged_out
        real(dp) :: h0

        call fit_fgarch_twist_range(y, x_range, 2*max_iter, gtol, f_out, params_out, niter_out, converged_out)
        persist_out = fgarch_twist_persist(params_out)
        h0 = fgarch_twist_range_h0(y, x_range, params_out, persist_out)
        vol_ann_out = sqrt(trading_days * h0) * 100.0_dp
        call fgarch_twist_range_skew_kurt(y, x_range, params_out, skew_out, ekurt_out)
    end subroutine fit_one_fgarch_twist_range

    subroutine fit_ohlc_model(component_model, y_co, y_oc, f_combo, params_out, vol_ann_out, persist_out, &
                              skew_out, ekurt_out, niter_out, converged_out, np_out, vol_out)
        character(len=*), intent(in) :: component_model
        real(dp), intent(in) :: y_co(:), y_oc(:)
        real(dp), intent(out) :: f_combo, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out, np_out
        logical, intent(out) :: converged_out
        real(dp), intent(out) :: vol_out(:)
        type(garch_params_t) :: p_co, p_oc
        real(dp), allocatable :: vol_co(:), vol_oc(:)
        real(dp) :: f_co, f_oc, persist_co, persist_oc, vol_ann_co, vol_ann_oc
        real(dp) :: skew_co, skew_oc, ekurt_co, ekurt_oc
        integer :: niter_co, niter_oc, np_component, t
        logical :: conv_co, conv_oc

        allocate(vol_co(size(y_co)), vol_oc(size(y_oc)))
        call fit_one_model(component_model, y_co, f_co, p_co, vol_ann_co, persist_co, &
                           skew_co, ekurt_co, niter_co, conv_co, np_component)
        call fit_one_model(component_model, y_oc, f_oc, p_oc, vol_ann_oc, persist_oc, &
                           skew_oc, ekurt_oc, niter_oc, conv_oc, np_component)
        call model_vol_forecast(component_model, y_co, p_co, persist_co, trading_days, vol_co)
        call model_vol_forecast(component_model, y_oc, p_oc, persist_oc, trading_days, vol_oc)
        do t = 1, size(vol_out)
            vol_out(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do

        f_combo = f_co + f_oc
        params_out = p_co
        persist_out = 0.5_dp * (persist_co + persist_oc)
        vol_ann_out = sqrt(sum(vol_out**2) / real(size(vol_out), dp))
        skew_out = 0.5_dp * (skew_co + skew_oc)
        ekurt_out = 0.5_dp * (ekurt_co + ekurt_oc)
        niter_out = niter_co + niter_oc
        converged_out = conv_co .and. conv_oc
        np_out = 2 * np_component
        deallocate(vol_co, vol_oc)
    end subroutine fit_ohlc_model

    subroutine fit_one_model(model_name, y, f_out, params_out, vol_ann_out, persist_out, &
                             skew_out, ekurt_out, niter_out, converged_out, np_out)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: f_out, vol_ann_out, persist_out, skew_out, ekurt_out
        type(garch_params_t), intent(out) :: params_out
        integer, intent(out) :: niter_out, np_out
        logical, intent(out) :: converged_out
        real(dp) :: h_unc, y_sd
        real(dp), allocatable :: variance_local(:)

        select case (trim(model_name))
        case ("SYMM_GARCH")
            call fit_symm_garch(y, max_iter, gtol, f_out, params_out, niter_out, converged_out)
            persist_out = symm_garch_persist(params_out)
            h_unc = params_out%omega / max(1.0_dp - persist_out, 1.0e-8_dp)
            vol_ann_out = sqrt(trading_days * h_unc) * 100.0_dp
            np_out = 3
            call garch_skew_kurt(y, params_out, skew_out, ekurt_out)
        case ("QGARCH")
            call fit_qgarch(y, max_iter, gtol, f_out, params_out, niter_out, converged_out)
            persist_out = qgarch_persist(params_out)
            h_unc = qgarch_mean_variance(params_out)
            vol_ann_out = sqrt(trading_days * h_unc) * 100.0_dp
            np_out = 4
            call qgarch_skew_kurt(y, params_out, skew_out, ekurt_out)
        case ("FIGARCH")
            call fit_figarch(y, max_iter, gtol, f_out, params_out, niter_out, converged_out)
            persist_out = figarch_persist(params_out)
            allocate(variance_local(size(y)))
            call figarch_variance(y, params_out, variance_local)
            h_unc = sum(variance_local) / real(size(y), dp)
            deallocate(variance_local)
            vol_ann_out = sqrt(trading_days * h_unc) * 100.0_dp
            np_out = 4
            call figarch_skew_kurt(y, params_out, skew_out, ekurt_out)
        case ("NAGARCH")
            call fit_nagarch(y, max_iter, gtol, f_out, params_out, niter_out, converged_out)
            persist_out = nagarch_persist(params_out)
            h_unc = params_out%omega / max(1.0_dp - persist_out, 1.0e-8_dp)
            vol_ann_out = sqrt(trading_days * h_unc) * 100.0_dp
            np_out = 4
            call nagarch_skew_kurt(y, params_out, skew_out, ekurt_out)
        case ("GJR_GARCH")
            call fit_gjr(y, max_iter, gtol, f_out, params_out, niter_out, converged_out)
            persist_out = gjr_persist(params_out)
            h_unc = params_out%omega / max(1.0_dp - persist_out, 1.0e-8_dp)
            vol_ann_out = sqrt(trading_days * h_unc) * 100.0_dp
            np_out = 4
            call gjr_skew_kurt(y, params_out, skew_out, ekurt_out)
        case ("FGTWIST")
            y_sd = sd(y)
            call fit_fgarch_twist(y, y_sd, max_iter, gtol, f_out, params_out, &
                                  vol_ann_out, skew_out, ekurt_out, niter_out, converged_out)
            persist_out = fgarch_twist_persist(params_out)
            np_out = 5
        case default
            f_out = huge(1.0_dp)
            params_out = garch_params_t()
            vol_ann_out = 0.0_dp
            persist_out = 0.0_dp
            skew_out = 0.0_dp
            ekurt_out = 0.0_dp
            niter_out = 0
            converged_out = .false.
            np_out = 0
        end select
    end subroutine fit_one_model

    subroutine nagarch_range_vol_forecast(y, x_range, params, persist, vol)
        real(dp), intent(in) :: y(:), x_range(:), persist
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h, sqrth, r
        integer :: t

        h = nagarch_range_h0(y, x_range, params, persist)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            vol(t) = sqrt(trading_days * h) * 100.0_dp
            r = y(t) - params%theta*sqrth
            h = params%omega + params%alpha*r**2 + params%gamma*x_range(t) + params%beta*h
        end do
    end subroutine nagarch_range_vol_forecast

    subroutine fgarch_twist_range_vol_forecast(y, x_range, params, persist, vol)
        real(dp), intent(in) :: y(:), x_range(:), persist
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h, sqrth, z, q
        integer :: t

        h = fgarch_twist_range_h0(y, x_range, params, persist)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            vol(t) = sqrt(trading_days * h) * 100.0_dp
            z = y(t) / sqrth
            q = abs(z - params%theta) - params%twist*(z - params%theta)
            h = params%omega + params%alpha*h*q**2 + params%gamma*x_range(t) + params%beta*h
        end do
    end subroutine fgarch_twist_range_vol_forecast

    subroutine rgarch_vol_forecast(y, x_range, params, vol)
        real(dp), intent(in) :: y(:), x_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h
        integer :: t

        h = rgarch_h0(y, x_range, params)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            vol(t) = sqrt(trading_days * h) * 100.0_dp
            h = params%omega + params%alpha*x_range(t) + params%beta*h
        end do
    end subroutine rgarch_vol_forecast

    subroutine carr_park_vol_forecast(x_range, params, vol)
        real(dp), intent(in) :: x_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h
        integer :: t

        h = carr_park_h0(x_range, params)
        do t = 1, size(x_range)
            h = max(h, 1.0e-12_dp)
            vol(t) = sqrt(trading_days * h) * 100.0_dp
            h = params%omega + params%alpha*x_range(t) + params%beta*h
        end do
    end subroutine carr_park_vol_forecast

    subroutine regarch1_vol_forecast(y, x_log_range, params, vol)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: lv, sig, z, x_range
        integer :: t

        lv = regarch1_h0(x_log_range)
        do t = 1, size(y)
            sig = exp(lv)
            vol(t) = sqrt(trading_days) * sig * 100.0_dp
            z = y(t) / max(sig, 1.0e-8_dp)
            x_range = (x_log_range(t) - 0.43_dp - lv) / 0.29_dp
            lv = params%omega + params%beta*lv + params%alpha*x_range + params%gamma*z
        end do
    end subroutine regarch1_vol_forecast

    subroutine regarch2_vol_forecast(y, x_log_range, params, vol)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: lh, lq, sig, z, x_range, lh_next, lq_next
        integer :: t

        lh = regarch1_h0(x_log_range)
        lq = lh
        do t = 1, size(y)
            sig = exp(lh)
            vol(t) = sqrt(trading_days) * sig * 100.0_dp
            z = y(t) / max(sig, 1.0e-8_dp)
            x_range = (x_log_range(t) - 0.43_dp - lh) / 0.29_dp
            lh_next = params%beta*lh + (1.0_dp - params%beta)*lq + params%alpha*x_range + params%gamma*z
            lq_next = (1.0_dp - params%scale)*params%omega + params%scale*lq + &
                      params%theta*x_range + params%twist*z
            lh = lh_next
            lq = lq_next
        end do
    end subroutine regarch2_vol_forecast

    subroutine rgarch_meas_vol_forecast(y, x_log_range, params, vol)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: lv, sig, z, meas_hat, u
        integer :: t

        lv = log(max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)) / 2.0_dp
        do t = 1, size(y)
            sig = exp(lv)
            vol(t) = sqrt(trading_days) * sig * 100.0_dp
            z = y(t) / max(sig, 1.0e-8_dp)
            meas_hat = params%theta + params%twist*lv + params%extra1*z + params%extra2*(z**2 - 1.0_dp)
            u = (x_log_range(t) - meas_hat) / max(params%scale, 1.0e-8_dp)
            lv = params%omega + params%beta*lv + params%alpha*u + params%gamma*z
        end do
    end subroutine rgarch_meas_vol_forecast

    real(dp) function regarch1_vol_ann(y, x_log_range, params)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: vol(:)

        allocate(vol(size(x_log_range)))
        call regarch1_vol_forecast(y, x_log_range, params, vol)
        regarch1_vol_ann = sqrt(sum(vol**2) / real(size(vol), dp))
        deallocate(vol)
    end function regarch1_vol_ann

    real(dp) function regarch2_vol_ann(y, x_log_range, params)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: vol(:)

        allocate(vol(size(x_log_range)))
        call regarch2_vol_forecast(y, x_log_range, params, vol)
        regarch2_vol_ann = sqrt(sum(vol**2) / real(size(vol), dp))
        deallocate(vol)
    end function regarch2_vol_ann

    real(dp) function rgarch_meas_vol_ann(y, x_log_range, params)
        real(dp), intent(in) :: y(:), x_log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: vol(:)

        allocate(vol(size(x_log_range)))
        call rgarch_meas_vol_forecast(y, x_log_range, params, vol)
        rgarch_meas_vol_ann = sqrt(sum(vol**2) / real(size(vol), dp))
        deallocate(vol)
    end function rgarch_meas_vol_ann

    real(dp) function carr_park_vol_ann(x_range, params)
        real(dp), intent(in) :: x_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: vol(:)

        allocate(vol(size(x_range)))
        call carr_park_vol_forecast(x_range, params, vol)
        carr_park_vol_ann = sqrt(sum(vol**2) / real(size(vol), dp))
        deallocate(vol)
    end function carr_park_vol_ann

    real(dp) function regarch1_h0(x_log_range)
        real(dp), intent(in) :: x_log_range(:)

        regarch1_h0 = sum(x_log_range - 0.43_dp) / real(size(x_log_range), dp)
    end function regarch1_h0

    real(dp) function rgarch_h0(y, x_range, params)
        real(dp), intent(in) :: y(:), x_range(:)
        type(garch_params_t), intent(in) :: params

        rgarch_h0 = max((params%omega + params%alpha*sum(x_range)/real(size(x_range), dp)) / &
            max(1.0_dp - params%beta, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function rgarch_h0

    real(dp) function carr_park_h0(x_range, params)
        real(dp), intent(in) :: x_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp) :: xbar

        xbar = sum(max(x_range, 1.0e-12_dp)) / real(size(x_range), dp)
        carr_park_h0 = max((params%omega + params%alpha*xbar) / &
            max(1.0_dp - params%beta, 1.0e-8_dp), xbar, 1.0e-12_dp)
    end function carr_park_h0

    real(dp) function nagarch_range_h0(y, x_range, params, persist)
        real(dp), intent(in) :: y(:), x_range(:), persist
        type(garch_params_t), intent(in) :: params

        nagarch_range_h0 = max((params%omega + params%gamma*sum(x_range)/real(size(x_range), dp)) / &
            max(1.0_dp - persist, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function nagarch_range_h0

    real(dp) function fgarch_twist_range_h0(y, x_range, params, persist)
        real(dp), intent(in) :: y(:), x_range(:), persist
        type(garch_params_t), intent(in) :: params

        fgarch_twist_range_h0 = max((params%omega + params%gamma*sum(x_range)/real(size(x_range), dp)) / &
            max(1.0_dp - persist, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function fgarch_twist_range_h0

    subroutine print_fit_assets()
        integer :: i

        if (fit_all_assets) then
            write(*,'(A)') "Fitting GARCH models for assets: all price-file assets"
        else
            write(*,'(A)', advance='no') "Fitting GARCH models for assets:"
            do i = 1, size(fit_assets)
                write(*,'(1X,A)', advance='no') trim(fit_assets(i))
            end do
            write(*,*)
        end if
    end subroutine print_fit_assets

end program xfit_garch_ohlc_iv_returns

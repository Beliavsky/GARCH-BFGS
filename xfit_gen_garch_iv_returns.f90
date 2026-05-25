! Fit selected GARCH-family Normal models to ETF demeaned log returns and
! compare model volatility forecasts with mapped implied-volatility indices.
! Edit models(:), fit_assets(:), asset_iv_assets(:), and asset_iv_indices(:) to configure.

program xfit_gen_garch_iv_returns
    use kind_mod,       only: dp
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv, print_price_sample_info
    use stats_mod,      only: mean, sd
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use vol_forecast_compare_mod, only: print_implied_vol_correlations
    use garch_forecast_mod, only: model_vol_forecast, vol_forecast_stats_t, summarize_vol_forecast, &
                                  print_vol_forecast_table
    use model_selection_mod, only: print_model_selection_counts
    use garch_fit_mod,  only: fit_symm_garch, fit_qgarch, fit_figarch, fit_fi_nagarch, &
                              fit_nagarch, fit_gjr, fit_fgarch_twist, fit_aparch, fit_harch, &
                              fit_riskmetrics2006, fit_ewma, &
                              fit_midas_hyperbolic, fit_midas_hyperbolic_asym, &
                              fit_aewma_nag, fit_aewma_twist, &
                              garch_skew_kurt, qgarch_skew_kurt, figarch_skew_kurt, fi_nagarch_skew_kurt, &
                              nagarch_skew_kurt, gjr_skew_kurt, aparch_skew_kurt, &
                              harch_skew_kurt, riskmetrics2006_skew_kurt, riskmetrics2006_variance, &
                              figarch_variance, fi_nagarch_variance, &
                              midas_hyperbolic_skew_kurt, midas_hyperbolic_asym_skew_kurt, &
                              ewma_skew_kurt, aewma_nag_skew_kurt, aewma_twist_skew_kurt, &
                              symm_garch_persist, qgarch_persist, figarch_persist, fi_nagarch_persist, &
                              nagarch_persist, gjr_persist, &
                              fgarch_twist_persist, ewma_persist, &
                              aewma_nag_persist, aewma_twist_persist, aparch_persist, harch_persist, &
                              riskmetrics2006_persist, midas_hyperbolic_persist, midas_hyperbolic_asym_persist, &
                              aparch_mean_variance, qgarch_mean_variance
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    character(len=*), parameter :: implied_vol_file = "vix_spy.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    logical,  parameter :: print_vol_forecast_stats = .false.
    logical,  parameter :: fit_all_assets = .true.
    logical,  parameter :: flip_log_return_sign = .true.
    logical,  parameter :: standardize_log_return_by_prior_iv = .true.
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "EWMA", "AEWMA_NAG", "AEWMA_TWIST", "SYMM_GARCH", "QGARCH", "FIGARCH", "FI_NAGARCH", "NAGARCH", &
        "GJR_GARCH", "FGTWIST", "APARCH", "HARCH", "RM2006", "MIDASHYP", "MIDASHYP_ASYM"]
    character(len=16), parameter :: fit_assets(*) = [character(len=16) :: "SPY"]
    character(len=16), parameter :: asset_iv_assets(*) = [character(len=16) :: "SPY"]
    character(len=16), parameter :: asset_iv_indices(*) = [character(len=16) :: "VIX"]
    character(len=16), parameter :: extra_series_names(*) = [character(len=16) :: "LOG_RETURN"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:), iv_dates(:)
    character(len=32), allocatable :: col_names(:), iv_col_names(:)
    real(dp), allocatable :: prices(:,:), iv_values(:,:), ret(:)
    type(garch_fit_result_t) :: rows(n_model)
    integer :: row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    character(len=256) :: aic_symbols(n_model)
    character(len=16) :: fit_model_names(n_model)
    character(len=16), allocatable :: vf_model(:)
    integer :: nprices, ncols, nobs, icol, imod, ifit, jfit, nfit
    integer :: niter, np, clock_start, clock_end, clock_rate
    real(dp) :: ret_mean, ret_std, fopt
    real(dp) :: persist, h_unc, vol_ann, logl, aic, bic, skew, ekurt, elapsed_s
    real(dp), allocatable :: vol_forecast(:), vol_forecasts(:,:,:), variance_tmp(:)
    real(dp), allocatable :: extra_series(:,:,:)
    real(dp) :: extra_signs(size(extra_series_names))
    type(vol_forecast_stats_t), allocatable :: vf_stats(:,:)
    type(garch_params_t) :: params
    logical :: converged
    logical :: fit_model_have(n_model)
    logical, allocatable :: vf_have(:,:)
    character(len=16) :: model

    call system_clock(clock_start, clock_rate)
    call nagarch_set_news_impact(.false.)
    aic_wins = 0
    bic_wins = 0
    aic_rank_sum = 0
    bic_rank_sum = 0
    rank_count = 0
    param_count = 0
    aic_symbols = ""
    fit_model_have = .false.
    fit_model_names = ""

    if (fit_all_assets) then
        call read_price_csv(prices_file, dates, col_names, prices)
    else
        call read_price_csv(prices_file, dates, col_names, prices, fit_assets)
    end if
    call read_price_csv(implied_vol_file, iv_dates, iv_col_names, iv_values)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), vol_forecast(nobs), vol_forecasts(nobs,ncols,n_model))
    allocate(variance_tmp(nobs))
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
    call print_fit_assets()
    write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist   scale  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
    write(*,'(A)') repeat("-", 180)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        ret_std = sd(ret)
        nfit = 0

        do imod = 1, size(models)
            model = uppercase(trim(models(imod)))
            params = garch_params_t()

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
            case ("FI_NAGARCH", "FINAGARCH")
                call fit_fi_nagarch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = fi_nagarch_persist(params)
                call fi_nagarch_variance(ret, params, variance_tmp)
                h_unc = sum(variance_tmp) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 5
                call fi_nagarch_skew_kurt(ret, params, skew, ekurt)
                model = "FI_NAGARCH"
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
            case ("APARCH")
                call fit_aparch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = aparch_persist(params)
                h_unc = aparch_mean_variance(ret, params)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 5
                call aparch_skew_kurt(ret, params, skew, ekurt)
            case ("HARCH")
                call fit_harch(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = harch_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call harch_skew_kurt(ret, params, skew, ekurt)
            case ("RM2006", "RISKMETRICS2006")
                call fit_riskmetrics2006(ret, fopt, params, niter, converged)
                persist = riskmetrics2006_persist()
                call riskmetrics2006_variance(ret, variance_tmp)
                h_unc = sum(variance_tmp) / real(nobs, dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 0
                call riskmetrics2006_skew_kurt(ret, skew, ekurt)
                model = "RM2006"
            case ("MIDASHYP", "MIDAS_HYPERBOLIC")
                call fit_midas_hyperbolic(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = midas_hyperbolic_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 3
                call midas_hyperbolic_skew_kurt(ret, params, skew, ekurt)
                model = "MIDASHYP"
            case ("MIDASHYP_ASYM", "MIDAS_HYPERBOLIC_ASYM")
                call fit_midas_hyperbolic_asym(ret, max_iter, gtol, fopt, params, niter, converged)
                persist = midas_hyperbolic_asym_persist(params)
                h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call midas_hyperbolic_asym_skew_kurt(ret, params, skew, ekurt)
                model = "MIDASHYP_ASYM"
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select

            logl = -real(nobs, dp) * fopt
            aic = 2.0_dp * real(np, dp) - 2.0_dp * logl
            bic = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl
            call model_vol_forecast(trim(model), ret, params, persist, trading_days, vol_forecast)
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
            param_count(imod) = np
        end do

        do ifit = 1, nfit
            row_aic_rank(ifit) = 1
            row_bic_rank(ifit) = 1
            do jfit = 1, nfit
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
            write(*,'(A16,1X,A9,ES12.3,7F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3,2I9)') &
                trim(rows(ifit)%model), trim(col_names(icol)), rows(ifit)%params%omega, rows(ifit)%params%alpha, &
                rows(ifit)%params%gamma, rows(ifit)%params%beta, rows(ifit)%params%theta, &
                rows(ifit)%params%twist, rows(ifit)%params%scale, rows(ifit)%persist, &
                rows(ifit)%vol_ann, rows(ifit)%logl, rows(ifit)%aic, rows(ifit)%bic, rows(ifit)%niter, &
                rows(ifit)%converged, rows(ifit)%skew, rows(ifit)%ekurt, row_aic_rank(ifit), row_bic_rank(ifit)
        end do
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    call print_model_selection_counts(models, param_count, aic_wins, bic_wins, aic_symbols, &
                                      aic_rank_sum, bic_rank_sum, rank_count)
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

end program xfit_gen_garch_iv_returns

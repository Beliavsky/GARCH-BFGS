! Fit selected GARCH-family Normal models to ETF demeaned log returns.
! Edit the models(:) character array to choose which model rows are produced.

program xfit_gen_garch_returns
    use date_mod, only: print_program_header
    use kind_mod,       only: dp
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv, print_price_sample_info
    use stats_mod,      only: mean, sd
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use garch_forecast_mod, only: model_vol_forecast, vol_forecast_stats_t, summarize_vol_forecast, &
                                  print_vol_forecast_table, finalize_return_garch_fit
    use model_selection_mod, only: print_model_selection_counts, print_model_fit_times
    use garch_fit_mod,  only: fit_symm_garch, fit_symm_garch_pq, fit_qgarch, fit_figarch, fit_fi_nagarch, &
                              fit_figarch_t, fit_fi_nagarch_t, fit_garch_m, &
                              fit_nagarch, fit_nagarch_pq, fit_gjr, fit_fgarch_twist, fit_aparch, fit_harch, fit_tgarch, &
                              fit_avgarch, &
                              fit_csgarch, fit_riskmetrics2006, &
                              fit_midas_hyperbolic, fit_midas_hyperbolic_asym
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    integer,  parameter :: max_p = 1
    integer,  parameter :: max_q = 1
    real(dp), parameter :: gtol = 1.0e-7_dp
    logical,  parameter :: print_vol_forecast_stats = .false.
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "GARCH_M", "SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2", &
        "QGARCH", "FIGARCH", "FIGARCH_T", "FI_NAGARCH", "FI_NAGARCH_T", &
        "NAGARCH", "NAGARCH_2_1", "NAGARCH_1_2", "NAGARCH_2_2", &
        "GJR_GARCH", "CSGARCH", "FGTWIST", "APARCH", "HARCH", "TGARCH", "AVGARCH", &
        "RM2006", "MIDASHYP", "MIDASHYP_ASYM"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:)
    type(garch_fit_result_t) :: rows(n_model)
    integer :: row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    character(len=256) :: aic_symbols(n_model)
    character(len=16), allocatable :: vf_model(:)
    integer :: nprices, ncols, nobs, icol, imod, ifit, jfit, nfit
    integer :: niter, clock_start, clock_end, clock_rate, fit_clock_start, fit_clock_end
    real(dp) :: ret_mean, ret_std, fopt
    real(dp) :: vol_ann, skew, ekurt, elapsed_s
    real(dp) :: fit_seconds(n_model)
    real(dp), allocatable :: vol_forecast(:)
    type(vol_forecast_stats_t), allocatable :: vf_stats(:,:)
    type(garch_params_t) :: params
    logical :: converged
    logical, allocatable :: vf_have(:,:)
    character(len=16) :: model

    call print_program_header("xfit_gen_garch_returns.f90")
    call system_clock(clock_start, clock_rate)
    print*,"max_iter:", max_iter ! debug
    call nagarch_set_news_impact(.false.)
    aic_wins = 0
    bic_wins = 0
    aic_rank_sum = 0
    bic_rank_sum = 0
    rank_count = 0
    param_count = 0
    fit_seconds = 0.0_dp
    aic_symbols = ""

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))
    if (print_vol_forecast_stats) then
        allocate(vol_forecast(nobs))
        allocate(vf_stats(ncols,n_model), vf_have(ncols,n_model))
        allocate(vf_model(n_model))
        vf_have = .false.
        vf_model = ""
    end if

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
    write(*,'(A)') repeat("-", 172)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        ret_std = sd(ret)
        nfit = 0

        do imod = 1, size(models)
            model = uppercase(trim(models(imod)))
            params = garch_params_t()
            call system_clock(fit_clock_start)

            select case (trim(model))
            case ("SYMM_GARCH", "GARCH", "SYMM")
                call fit_symm_garch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "SYMM_GARCH"
            case ("GARCH_M")
                call fit_garch_m(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH_2_1")
                if (max_p < 2 .or. max_q < 1) cycle
                call fit_symm_garch_pq(ret, 2, 1, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH_1_2")
                if (max_p < 1 .or. max_q < 2) cycle
                call fit_symm_garch_pq(ret, 1, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH_2_2")
                if (max_p < 2 .or. max_q < 2) cycle
                call fit_symm_garch_pq(ret, 2, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("QGARCH", "QUADRATIC_GARCH")
                call fit_qgarch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "QGARCH"
            case ("FIGARCH")
                call fit_figarch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("FIGARCH_T")
                call fit_figarch_t(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("FI_NAGARCH", "FINAGARCH")
                call fit_fi_nagarch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "FI_NAGARCH"
            case ("FI_NAGARCH_T")
                call fit_fi_nagarch_t(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "FI_NAGARCH_T"
            case ("NAGARCH")
                call fit_nagarch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("NAGARCH_2_1")
                if (max_p < 2 .or. max_q < 1) cycle
                call fit_nagarch_pq(ret, 2, 1, max_iter, gtol, fopt, params, niter, converged)
            case ("NAGARCH_1_2")
                if (max_p < 1 .or. max_q < 2) cycle
                call fit_nagarch_pq(ret, 1, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("NAGARCH_2_2")
                if (max_p < 2 .or. max_q < 2) cycle
                call fit_nagarch_pq(ret, 2, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("GJR_GARCH", "GJR")
                call fit_gjr(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "GJR_GARCH"
            case ("CSGARCH", "COMPONENT_GARCH")
                call fit_csgarch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "CSGARCH"
            case ("FGTWIST", "FGARCH_TWIST")
                call fit_fgarch_twist(ret, ret_std, max_iter, gtol, fopt, params, &
                                      vol_ann, skew, ekurt, niter, converged)
                model = "FGTWIST"
            case ("APARCH")
                call fit_aparch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("HARCH")
                call fit_harch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("TGARCH")
                call fit_tgarch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("AVGARCH")
                call fit_avgarch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("RM2006", "RISKMETRICS2006")
                call fit_riskmetrics2006(ret, fopt, params, niter, converged)
                model = "RM2006"
            case ("MIDASHYP", "MIDAS_HYPERBOLIC")
                call fit_midas_hyperbolic(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "MIDASHYP"
            case ("MIDASHYP_ASYM", "MIDAS_HYPERBOLIC_ASYM")
                call fit_midas_hyperbolic_asym(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "MIDASHYP_ASYM"
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select
            call system_clock(fit_clock_end)
            fit_seconds(imod) = fit_seconds(imod) + real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)

            nfit = nfit + 1
            row_model_idx(nfit) = imod
            if (trim(model) == "FGTWIST") then
                call finalize_return_garch_fit(model, ret, params, fopt, niter, converged, trading_days, rows(nfit), &
                                               vol_ann, skew, ekurt)
            else
                call finalize_return_garch_fit(model, ret, params, fopt, niter, converged, trading_days, rows(nfit))
            end if
            if (print_vol_forecast_stats) then
                call model_vol_forecast(trim(model), ret, rows(nfit)%params, rows(nfit)%persist, trading_days, vol_forecast)
                call summarize_vol_forecast(vol_forecast, dates(2:nprices), vf_stats(icol,imod))
                vf_model(imod) = model
                vf_have(icol,imod) = .true.
            end if

            param_count(imod) = rows(nfit)%nparam
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
            write(*,'(A16,1X,A9,ES12.3,6F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3,2I9)') &
                trim(rows(ifit)%model), trim(col_names(icol)), rows(ifit)%params%omega, rows(ifit)%params%alpha, &
                rows(ifit)%params%gamma, rows(ifit)%params%beta, rows(ifit)%params%theta, &
                rows(ifit)%params%twist, rows(ifit)%persist, &
                rows(ifit)%vol_ann, rows(ifit)%logl, rows(ifit)%aic, rows(ifit)%bic, rows(ifit)%niter, &
                rows(ifit)%converged, rows(ifit)%skew, rows(ifit)%ekurt, row_aic_rank(ifit), row_bic_rank(ifit)
        end do
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    call print_model_selection_counts(models, param_count, aic_wins, bic_wins, aic_symbols, &
                                      aic_rank_sum, bic_rank_sum, rank_count)
    call print_model_fit_times(models, param_count, fit_seconds)
    if (print_vol_forecast_stats) call print_vol_forecast_table(col_names, vf_model, vf_stats, vf_have)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

end program xfit_gen_garch_returns

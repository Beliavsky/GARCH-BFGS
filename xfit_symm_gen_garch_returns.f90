! Compare symmetric-news-impact GARCH-type models on demeaned ETF log returns.
! Edit models(:), max_p, and max_q to choose which rows are produced.

program xfit_symm_gen_garch_returns
    use kind_mod,        only: dp
    use strings_mod,     only: uppercase
    use csv_mod,         only: read_price_csv, print_price_sample_info
    use stats_mod,       only: mean
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use garch_forecast_mod, only: finalize_return_garch_fit
    use model_selection_mod, only: print_model_selection_counts, print_model_fit_times
    use garch_fit_mod,   only: fit_ewma, fit_symm_garch, fit_symm_garch_pq, fit_figarch, &
                               fit_csgarch, fit_harch, fit_riskmetrics2006, &
                               fit_midas_hyperbolic
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    integer,  parameter :: max_p = 2
    integer,  parameter :: max_q = 2
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2", &
        "FIGARCH", "CSGARCH", "HARCH", "RM2006", "MIDASHYP", "EWMA"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:)
    type(garch_fit_result_t) :: rows(n_model)
    integer :: row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    character(len=256) :: aic_symbols(n_model)
    integer :: nprices, ncols, nobs, icol, imod, ifit, jfit, nfit
    integer :: niter, clock_start, clock_end, clock_rate, fit_clock_start, fit_clock_end
    real(dp) :: ret_mean, fopt, elapsed_s
    real(dp) :: fit_seconds(n_model)
    type(garch_params_t) :: params
    logical :: converged
    character(len=16) :: model

    call system_clock(clock_start, clock_rate)
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

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A,I0,A,I0,A)') "Higher-order symmetric GARCH enabled through max_p=", max_p, ", max_q=", max_q, "."
    write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
    write(*,'(A)') repeat("-", 172)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        nfit = 0

        do imod = 1, size(models)
            model = uppercase(trim(models(imod)))
            params = garch_params_t()
            call system_clock(fit_clock_start)

            select case (trim(model))
            case ("EWMA")
                call fit_ewma(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH", "GARCH", "SYMM")
                call fit_symm_garch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "SYMM_GARCH"
            case ("SYMM_GARCH_2_1")
                if (max_p < 2 .or. max_q < 1) cycle
                call fit_symm_garch_pq(ret, 2, 1, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH_1_2")
                if (max_p < 1 .or. max_q < 2) cycle
                call fit_symm_garch_pq(ret, 1, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("SYMM_GARCH_2_2")
                if (max_p < 2 .or. max_q < 2) cycle
                call fit_symm_garch_pq(ret, 2, 2, max_iter, gtol, fopt, params, niter, converged)
            case ("FIGARCH")
                call fit_figarch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("CSGARCH", "COMPONENT_GARCH")
                call fit_csgarch(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "CSGARCH"
            case ("HARCH")
                call fit_harch(ret, max_iter, gtol, fopt, params, niter, converged)
            case ("RM2006", "RISKMETRICS2006")
                call fit_riskmetrics2006(ret, fopt, params, niter, converged)
                model = "RM2006"
            case ("MIDASHYP", "MIDAS_HYPERBOLIC")
                call fit_midas_hyperbolic(ret, max_iter, gtol, fopt, params, niter, converged)
                model = "MIDASHYP"
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select
            call system_clock(fit_clock_end)
            fit_seconds(imod) = fit_seconds(imod) + real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)

            nfit = nfit + 1
            row_model_idx(nfit) = imod
            call finalize_return_garch_fit(model, ret, params, fopt, niter, converged, trading_days, rows(nfit))
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
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

end program xfit_symm_gen_garch_returns

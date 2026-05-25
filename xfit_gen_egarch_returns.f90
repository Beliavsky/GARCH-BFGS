! Fit selected GARCH/EGARCH-family Normal models to ETF demeaned log returns.
! Edit the models(:) character array to choose which model rows are produced.

program xfit_gen_egarch_returns
    use kind_mod,       only: dp
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv, print_price_sample_info
    use stats_mod,      only: mean, sd
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use model_selection_mod, only: print_model_selection_counts
    use garch_fit_mod,  only: fit_symm_garch, fit_qgarch, fit_figarch, fit_nagarch, fit_gjr, fit_egarch, fit_fgarch_twist, &
                              figarch_variance, &
                              garch_skew_kurt, qgarch_skew_kurt, figarch_skew_kurt, nagarch_skew_kurt, gjr_skew_kurt, egarch_skew_kurt, &
                              symm_garch_persist, qgarch_persist, figarch_persist, nagarch_persist, gjr_persist, egarch_persist, &
                              fgarch_twist_persist, qgarch_mean_variance
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "QGARCH", "FIGARCH", "NAGARCH", "GJR_GARCH", "EGARCH", "FGTWIST"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:), variance_tmp(:)
    type(garch_fit_result_t) :: rows(n_model)
    integer :: row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    character(len=256) :: aic_symbols(n_model)
    integer :: nprices, ncols, nobs, icol, imod, ifit, jfit, nfit
    integer :: niter, np, clock_start, clock_end, clock_rate
    real(dp) :: ret_mean, ret_std, fopt
    real(dp) :: persist, h_unc, vol_ann, logl, aic, bic, skew, ekurt, elapsed_s
    type(garch_params_t) :: params
    logical :: converged
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

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), variance_tmp(nobs))

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

            select case (trim(model))
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
            case ("EGARCH")
                call fit_egarch(ret, max_iter, 1.0e-6_dp, fopt, params, niter, converged)
                persist = egarch_persist(params)
                h_unc = exp(params%omega / (1.0_dp - params%beta))
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
                np = 4
                call egarch_skew_kurt(ret, params, skew, ekurt)
            case ("FGTWIST", "FGARCH_TWIST")
                call fit_fgarch_twist(ret, ret_std, max_iter, gtol, fopt, params, &
                                      vol_ann, skew, ekurt, niter, converged)
                persist = fgarch_twist_persist(params)
                np = 5
                model = "FGTWIST"
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select

            logl = -real(nobs, dp) * fopt
            aic = 2.0_dp * real(np, dp) - 2.0_dp * logl
            bic = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl

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
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

end program xfit_gen_egarch_returns

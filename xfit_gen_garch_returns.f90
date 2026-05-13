! Fit selected GARCH-family Normal models to ETF demeaned log returns.
! Edit the models(:) character array to choose which model rows are produced.

program xfit_gen_garch_returns
    use kind_mod,       only: dp
    use strings_mod,    only: uppercase
    use csv_mod,        only: read_price_csv
    use stats_mod,      only: mean, sd
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,  only: fit_symm_garch, fit_nagarch, fit_gjr, fit_fgarch_twist, &
                              garch_skew_kurt, nagarch_skew_kurt, gjr_skew_kurt, &
                              symm_garch_persist, nagarch_persist, gjr_persist, fgarch_twist_persist
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "NAGARCH", "GJR_GARCH", "FGTWIST"]
    integer, parameter :: n_model = size(models)

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:)
    character(len=16) :: row_model(n_model)
    type(garch_params_t) :: row_params(n_model)
    real(dp) :: row_persist(n_model), row_vol_ann(n_model)
    real(dp) :: row_logl(n_model), row_aic(n_model), row_bic(n_model), row_skew(n_model), row_ekurt(n_model)
    integer :: row_iter(n_model), row_np(n_model), row_aic_rank(n_model), row_bic_rank(n_model), row_model_idx(n_model)
    logical :: row_conv(n_model)
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
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
    param_count = 0
    aic_symbols = ""

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))

    write(*,'(A,I0,A,I0,A)') "Using ", nobs, " demeaned log returns for ", ncols, " assets"
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
            case default
                write(*,'(A,A)') "Skipping unknown model: ", trim(models(imod))
                cycle
            end select

            logl = -real(nobs, dp) * fopt
            aic = 2.0_dp * real(np, dp) - 2.0_dp * logl
            bic = log(real(nobs, dp)) * real(np, dp) - 2.0_dp * logl

            nfit = nfit + 1
            row_model(nfit) = model
            row_model_idx(nfit) = imod
            row_params(nfit) = params
            row_persist(nfit) = persist
            row_vol_ann(nfit) = vol_ann
            row_logl(nfit) = logl
            row_aic(nfit) = aic
            row_bic(nfit) = bic
            row_np(nfit) = np
            row_iter(nfit) = niter
            row_conv(nfit) = converged
            row_skew(nfit) = skew
            row_ekurt(nfit) = ekurt
            param_count(imod) = np
        end do

        do ifit = 1, nfit
            row_aic_rank(ifit) = 1
            row_bic_rank(ifit) = 1
            do jfit = 1, nfit
                if (row_aic(jfit) < row_aic(ifit)) row_aic_rank(ifit) = row_aic_rank(ifit) + 1
                if (row_bic(jfit) < row_bic(ifit)) row_bic_rank(ifit) = row_bic_rank(ifit) + 1
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
        end do

        do ifit = 1, nfit
            write(*,'(A16,1X,A9,ES12.3,6F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3,2I9)') &
                trim(row_model(ifit)), trim(col_names(icol)), row_params(ifit)%omega, row_params(ifit)%alpha, &
                row_params(ifit)%gamma, row_params(ifit)%beta, row_params(ifit)%theta, &
                row_params(ifit)%twist, row_persist(ifit), &
                row_vol_ann(ifit), row_logl(ifit), row_aic(ifit), row_bic(ifit), row_iter(ifit), &
                row_conv(ifit), row_skew(ifit), row_ekurt(ifit), row_aic_rank(ifit), row_bic_rank(ifit)
        end do
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A)') "Model selection counts:"
    write(*,'(A)') "Model            #param  AIC_wins  BIC_wins  AIC_symbols"
    write(*,'(A)') repeat("-", 72)
    do imod = 1, size(models)
        write(*,'(A16,I8,2I10,2X,A)') trim(uppercase(trim(models(imod)))), param_count(imod), &
            aic_wins(imod), bic_wins(imod), trim(aic_symbols(imod))
    end do
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

end program xfit_gen_garch_returns

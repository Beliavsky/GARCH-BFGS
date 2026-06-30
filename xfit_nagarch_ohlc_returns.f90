! Fit close-close NAGARCH and OHLC NAGARCH components to adjusted OHLC prices.
! CO is log(Open_t / Close_{t-1}); OC is log(Close_t / Open_t).

program xfit_nagarch_ohlc_returns
    use date_mod, only: print_program_header
    use kind_mod,       only: dp
    use csv_mod,        only: read_ohlc_csv, print_price_sample_info
    use stats_mod,      only: mean
    use nagarch_mod,    only: nagarch_set_news_impact
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,  only: fit_nagarch, fit_rgarch, fit_nagarch_range, nagarch_skew_kurt, &
                              rgarch_skew_kurt, nagarch_range_skew_kurt, nagarch_persist, rgarch_persist
    implicit none

    character(len=*), parameter :: prices_file = "prices_ohlc.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    real(dp), parameter :: co_oc_corr = 0.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    logical,  parameter :: fit_all_assets = .true.
    character(len=16), parameter :: fit_assets(*) = [character(len=16) :: "SPY"]

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: asset_names(:)
    real(dp), allocatable :: open_prices(:,:), high_prices(:,:), low_prices(:,:), close_prices(:,:)
    real(dp), allocatable :: ret_cc(:), ret_co(:), ret_oc(:), range_var(:), vol_co(:), vol_oc(:), vol_ohlc(:)
    type(garch_params_t), allocatable :: p_co_all(:), p_oc_all(:), p_co_rg_all(:), p_oc_rg_all(:)
    type(garch_params_t), allocatable :: p_co_rng_all(:), p_oc_rng_all(:)
    real(dp), allocatable :: persist_co_all(:), persist_oc_all(:), persist_co_rg_all(:), persist_oc_rg_all(:)
    real(dp), allocatable :: persist_co_rng_all(:), persist_oc_rng_all(:)
    integer :: nprices, nassets, nobs, iasset
    integer :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s

    type(garch_params_t) :: p_cc, p_cc_rg, p_co, p_oc, p_co_rg, p_oc_rg, p_co_rng, p_oc_rng
    real(dp) :: f_cc, f_cc_rg, f_co, f_oc, f_co_rg, f_oc_rg, f_co_rng, f_oc_rng
    real(dp) :: persist_cc, persist_cc_rg, persist_co, persist_oc, persist_co_rg, persist_oc_rg, persist_co_rng, persist_oc_rng
    real(dp) :: vol_ann_cc, vol_ann_cc_rg, vol_ann_co, vol_ann_oc, vol_ann_ohlc, vol_ann_co_rg, vol_ann_oc_rg
    real(dp) :: vol_ann_ohlc_rg, vol_ann_co_rng, vol_ann_oc_rng, vol_ann_ohlc_rng
    real(dp) :: logl_cc, logl_cc_rg, logl_co, logl_oc, logl_ohlc, logl_co_rg, logl_oc_rg, logl_ohlc_rg
    real(dp) :: logl_co_rng, logl_oc_rng, logl_ohlc_rng
    real(dp) :: aic_cc, bic_cc, aic_cc_rg, bic_cc_rg, aic_co, bic_co, aic_oc, bic_oc, aic_ohlc, bic_ohlc
    real(dp) :: aic_co_rg, bic_co_rg, aic_oc_rg, bic_oc_rg, aic_ohlc_rg, bic_ohlc_rg, aic_co_rng, bic_co_rng
    real(dp) :: aic_oc_rng, bic_oc_rng, aic_ohlc_rng, bic_ohlc_rng
    real(dp) :: skew_cc, ekurt_cc, skew_cc_rg, ekurt_cc_rg, skew_co, ekurt_co, skew_oc, ekurt_oc
    real(dp) :: skew_co_rg, ekurt_co_rg, skew_oc_rg, ekurt_oc_rg, skew_co_rng, ekurt_co_rng, skew_oc_rng, ekurt_oc_rng
    integer :: niter_cc, niter_cc_rg, niter_co, niter_oc, niter_co_rg, niter_oc_rg, niter_co_rng, niter_oc_rng, t
    logical :: conv_cc, conv_cc_rg, conv_co, conv_oc, conv_co_rg, conv_oc_rg, conv_co_rng, conv_oc_rng

    call print_program_header("xfit_nagarch_ohlc_returns.f90")
    call system_clock(clock_start, clock_rate)
    call nagarch_set_news_impact(.false.)
    if (fit_all_assets) then
        call read_ohlc_csv(prices_file, dates, asset_names, open_prices, close_prices, &
                           high_prices=high_prices, low_prices=low_prices)
    else
        call read_ohlc_csv(prices_file, dates, asset_names, open_prices, close_prices, fit_assets, high_prices, low_prices)
    end if

    nprices = size(close_prices, 1)
    nassets = size(close_prices, 2)
    nobs = nprices - 1
    allocate(ret_cc(nobs), ret_co(nobs), ret_oc(nobs), range_var(nobs), vol_co(nobs), vol_oc(nobs), vol_ohlc(nobs))
    allocate(p_co_all(nassets), p_oc_all(nassets), p_co_rg_all(nassets), p_oc_rg_all(nassets))
    allocate(p_co_rng_all(nassets), p_oc_rng_all(nassets))
    allocate(persist_co_all(nassets), persist_oc_all(nassets), persist_co_rg_all(nassets), persist_oc_rg_all(nassets))
    allocate(persist_co_rng_all(nassets), persist_oc_rng_all(nassets))

    call print_price_sample_info(prices_file, dates, nassets)
    write(*,'(A,F7.3)') "OHLC CO/OC variance forecast correlation: ", co_oc_corr
    call print_parameter_descriptions()
    write(*,'(A)') "Asset       Piece          omega    alpha    gamma     beta    theta  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
    write(*,'(A)') repeat("-", 163)

    do iasset = 1, nassets
        ret_cc = log(close_prices(2:nprices,iasset) / close_prices(1:nprices-1,iasset))
        ret_cc = ret_cc - mean(ret_cc)
        ret_co = log(open_prices(2:nprices,iasset) / close_prices(1:nprices-1,iasset))
        ret_co = ret_co - mean(ret_co)
        ret_oc = log(close_prices(2:nprices,iasset) / open_prices(2:nprices,iasset))
        ret_oc = ret_oc - mean(ret_oc)
        range_var = (log(high_prices(2:nprices,iasset) / low_prices(2:nprices,iasset)))**2 / (4.0_dp*log(2.0_dp))

        call fit_nagarch_piece(ret_cc, p_cc, f_cc, persist_cc, vol_ann_cc, logl_cc, aic_cc, bic_cc, &
                               niter_cc, conv_cc, skew_cc, ekurt_cc)
        call fit_rgarch_piece(ret_cc, range_var, p_cc_rg, f_cc_rg, persist_cc_rg, vol_ann_cc_rg, &
                              logl_cc_rg, aic_cc_rg, bic_cc_rg, niter_cc_rg, conv_cc_rg, skew_cc_rg, ekurt_cc_rg)
        call fit_nagarch_piece(ret_co, p_co, f_co, persist_co, vol_ann_co, logl_co, aic_co, bic_co, &
                               niter_co, conv_co, skew_co, ekurt_co)
        call fit_nagarch_piece(ret_oc, p_oc, f_oc, persist_oc, vol_ann_oc, logl_oc, aic_oc, bic_oc, &
                               niter_oc, conv_oc, skew_oc, ekurt_oc)
        call fit_rgarch_piece(ret_co, range_var, p_co_rg, f_co_rg, persist_co_rg, vol_ann_co_rg, &
                              logl_co_rg, aic_co_rg, bic_co_rg, niter_co_rg, conv_co_rg, skew_co_rg, ekurt_co_rg)
        call fit_rgarch_piece(ret_oc, range_var, p_oc_rg, f_oc_rg, persist_oc_rg, vol_ann_oc_rg, &
                              logl_oc_rg, aic_oc_rg, bic_oc_rg, niter_oc_rg, conv_oc_rg, skew_oc_rg, ekurt_oc_rg)
        call fit_nagarch_range_piece(ret_co, range_var, p_co_rng, f_co_rng, persist_co_rng, vol_ann_co_rng, &
                                     logl_co_rng, aic_co_rng, bic_co_rng, niter_co_rng, conv_co_rng, &
                                     skew_co_rng, ekurt_co_rng)
        call fit_nagarch_range_piece(ret_oc, range_var, p_oc_rng, f_oc_rng, persist_oc_rng, vol_ann_oc_rng, &
                                     logl_oc_rng, aic_oc_rng, bic_oc_rng, niter_oc_rng, conv_oc_rng, &
                                     skew_oc_rng, ekurt_oc_rng)
        p_co_all(iasset) = p_co
        p_oc_all(iasset) = p_oc
        p_co_rg_all(iasset) = p_co_rg
        p_oc_rg_all(iasset) = p_oc_rg
        p_co_rng_all(iasset) = p_co_rng
        p_oc_rng_all(iasset) = p_oc_rng
        persist_co_all(iasset) = persist_co
        persist_oc_all(iasset) = persist_oc
        persist_co_rg_all(iasset) = persist_co_rg
        persist_oc_rg_all(iasset) = persist_oc_rg
        persist_co_rng_all(iasset) = persist_co_rng
        persist_oc_rng_all(iasset) = persist_oc_rng

        call nagarch_vol_forecast(ret_co, p_co, persist_co, vol_co)
        call nagarch_vol_forecast(ret_oc, p_oc, persist_oc, vol_oc)
        do t = 1, nobs
            vol_ohlc(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do
        vol_ann_ohlc = sqrt(sum(vol_ohlc**2) / real(nobs, dp))
        logl_ohlc = logl_co + logl_oc
        aic_ohlc = 2.0_dp*8.0_dp - 2.0_dp*logl_ohlc
        bic_ohlc = log(real(nobs, dp))*8.0_dp - 2.0_dp*logl_ohlc
        call rgarch_vol_forecast(ret_co, range_var, p_co_rg, vol_co)
        call rgarch_vol_forecast(ret_oc, range_var, p_oc_rg, vol_oc)
        do t = 1, nobs
            vol_ohlc(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do
        vol_ann_ohlc_rg = sqrt(sum(vol_ohlc**2) / real(nobs, dp))
        logl_ohlc_rg = logl_co_rg + logl_oc_rg
        aic_ohlc_rg = 2.0_dp*6.0_dp - 2.0_dp*logl_ohlc_rg
        bic_ohlc_rg = log(real(nobs, dp))*6.0_dp - 2.0_dp*logl_ohlc_rg

        call print_piece(asset_names(iasset), "CC_NAGARCH", p_cc, persist_cc, vol_ann_cc, logl_cc, &
                         aic_cc, bic_cc, niter_cc, conv_cc, skew_cc, ekurt_cc, &
                         rank_in_values(aic_cc, [aic_cc, aic_cc_rg]), rank_in_values(bic_cc, [bic_cc, bic_cc_rg]))
        call print_piece(asset_names(iasset), "CC_RGARCH", p_cc_rg, persist_cc_rg, vol_ann_cc_rg, logl_cc_rg, &
                         aic_cc_rg, bic_cc_rg, niter_cc_rg, conv_cc_rg, skew_cc_rg, ekurt_cc_rg, &
                         rank_in_values(aic_cc_rg, [aic_cc, aic_cc_rg]), rank_in_values(bic_cc_rg, [bic_cc, bic_cc_rg]))
        call print_piece(asset_names(iasset), "CO_NAGARCH", p_co, persist_co, vol_ann_co, logl_co, &
                         aic_co, bic_co, niter_co, conv_co, skew_co, ekurt_co, &
                         rank_in_values(aic_co, [aic_co, aic_co_rg, aic_co_rng]), &
                         rank_in_values(bic_co, [bic_co, bic_co_rg, bic_co_rng]))
        call print_piece(asset_names(iasset), "OC_NAGARCH", p_oc, persist_oc, vol_ann_oc, logl_oc, &
                         aic_oc, bic_oc, niter_oc, conv_oc, skew_oc, ekurt_oc, &
                         rank_in_values(aic_oc, [aic_oc, aic_oc_rg, aic_oc_rng]), &
                         rank_in_values(bic_oc, [bic_oc, bic_oc_rg, bic_oc_rng]))
        call print_piece(asset_names(iasset), "OHLC_TOTAL", garch_params_t(), 0.5_dp*(persist_co + persist_oc), &
                         vol_ann_ohlc, logl_ohlc, aic_ohlc, bic_ohlc, niter_co + niter_oc, &
                         conv_co .and. conv_oc, 0.5_dp*(skew_co + skew_oc), 0.5_dp*(ekurt_co + ekurt_oc), &
                         rank_in_values(aic_ohlc, [aic_ohlc, aic_ohlc_rg, aic_ohlc_rng]), &
                         rank_in_values(bic_ohlc, [bic_ohlc, bic_ohlc_rg, bic_ohlc_rng]))
        call print_piece(asset_names(iasset), "CO_RGARCH", p_co_rg, persist_co_rg, vol_ann_co_rg, logl_co_rg, &
                         aic_co_rg, bic_co_rg, niter_co_rg, conv_co_rg, skew_co_rg, ekurt_co_rg, &
                         rank_in_values(aic_co_rg, [aic_co, aic_co_rg, aic_co_rng]), &
                         rank_in_values(bic_co_rg, [bic_co, bic_co_rg, bic_co_rng]))
        call print_piece(asset_names(iasset), "OC_RGARCH", p_oc_rg, persist_oc_rg, vol_ann_oc_rg, logl_oc_rg, &
                         aic_oc_rg, bic_oc_rg, niter_oc_rg, conv_oc_rg, skew_oc_rg, ekurt_oc_rg, &
                         rank_in_values(aic_oc_rg, [aic_oc, aic_oc_rg, aic_oc_rng]), &
                         rank_in_values(bic_oc_rg, [bic_oc, bic_oc_rg, bic_oc_rng]))
        call print_piece(asset_names(iasset), "OHLC_RGARCH", garch_params_t(), 0.5_dp*(persist_co_rg + persist_oc_rg), &
                         vol_ann_ohlc_rg, logl_ohlc_rg, aic_ohlc_rg, bic_ohlc_rg, niter_co_rg + niter_oc_rg, &
                         conv_co_rg .and. conv_oc_rg, 0.5_dp*(skew_co_rg + skew_oc_rg), &
                         0.5_dp*(ekurt_co_rg + ekurt_oc_rg), &
                         rank_in_values(aic_ohlc_rg, [aic_ohlc, aic_ohlc_rg, aic_ohlc_rng]), &
                         rank_in_values(bic_ohlc_rg, [bic_ohlc, bic_ohlc_rg, bic_ohlc_rng]))
        call nagarch_range_vol_forecast(ret_co, range_var, p_co_rng, persist_co_rng, vol_co)
        call nagarch_range_vol_forecast(ret_oc, range_var, p_oc_rng, persist_oc_rng, vol_oc)
        do t = 1, nobs
            vol_ohlc(t) = sqrt(max(vol_co(t)**2 + vol_oc(t)**2 + 2.0_dp*co_oc_corr*vol_co(t)*vol_oc(t), 0.0_dp))
        end do
        vol_ann_ohlc_rng = sqrt(sum(vol_ohlc**2) / real(nobs, dp))
        logl_ohlc_rng = logl_co_rng + logl_oc_rng
        aic_ohlc_rng = 2.0_dp*10.0_dp - 2.0_dp*logl_ohlc_rng
        bic_ohlc_rng = log(real(nobs, dp))*10.0_dp - 2.0_dp*logl_ohlc_rng
        call print_piece(asset_names(iasset), "CO_NAG_RANGE", p_co_rng, persist_co_rng, vol_ann_co_rng, &
                         logl_co_rng, aic_co_rng, bic_co_rng, niter_co_rng, conv_co_rng, skew_co_rng, ekurt_co_rng, &
                         rank_in_values(aic_co_rng, [aic_co, aic_co_rg, aic_co_rng]), &
                         rank_in_values(bic_co_rng, [bic_co, bic_co_rg, bic_co_rng]))
        call print_piece(asset_names(iasset), "OC_NAG_RANGE", p_oc_rng, persist_oc_rng, vol_ann_oc_rng, &
                         logl_oc_rng, aic_oc_rng, bic_oc_rng, niter_oc_rng, conv_oc_rng, skew_oc_rng, ekurt_oc_rng, &
                         rank_in_values(aic_oc_rng, [aic_oc, aic_oc_rg, aic_oc_rng]), &
                         rank_in_values(bic_oc_rng, [bic_oc, bic_oc_rg, bic_oc_rng]))
        call print_piece(asset_names(iasset), "OHLC_RANGE", garch_params_t(), &
                         0.5_dp*(persist_co_rng + persist_oc_rng), vol_ann_ohlc_rng, logl_ohlc_rng, &
                         aic_ohlc_rng, bic_ohlc_rng, niter_co_rng + niter_oc_rng, conv_co_rng .and. conv_oc_rng, &
                         0.5_dp*(skew_co_rng + skew_oc_rng), 0.5_dp*(ekurt_co_rng + ekurt_oc_rng), &
                         rank_in_values(aic_ohlc_rng, [aic_ohlc, aic_ohlc_rg, aic_ohlc_rng]), &
                         rank_in_values(bic_ohlc_rng, [bic_ohlc, bic_ohlc_rg, bic_ohlc_rng]))
    end do

    call print_ohlc_component_parameter_table()

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    subroutine print_parameter_descriptions()
        write(*,'(/,A)') "Parameter descriptions:"
        write(*,'(A)') "  omega: variance intercept, the baseline variance contribution each period."
        write(*,'(A)') "  alpha: response to squared return news, or to Parkinson high-low variance in RGARCH rows."
        write(*,'(A)') "  gamma: response to the high-low range variance term; zero for models without range input."
        write(*,'(A)') "  beta : persistence of the previous conditional variance forecast."
        write(*,'(A,/)') "  theta: NAGARCH leverage shift in the news term, where nonzero values make positive and negative returns affect variance differently."
    end subroutine print_parameter_descriptions

    subroutine fit_nagarch_piece(y, params, fopt, persist, vol_ann, logl, aic, bic, niter, converged, skew, ekurt)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: fopt, persist, vol_ann, logl, aic, bic, skew, ekurt
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        real(dp) :: h_unc

        call fit_nagarch(y, max_iter, gtol, fopt, params, niter, converged)
        persist = nagarch_persist(params)
        h_unc = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
        vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        logl = -real(size(y), dp) * fopt
        aic = 2.0_dp*4.0_dp - 2.0_dp*logl
        bic = log(real(size(y), dp))*4.0_dp - 2.0_dp*logl
        call nagarch_skew_kurt(y, params, skew, ekurt)
    end subroutine fit_nagarch_piece

    subroutine fit_nagarch_range_piece(y, range_var, params, fopt, persist, vol_ann, logl, aic, bic, &
                                       niter, converged, skew, ekurt)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: fopt, persist, vol_ann, logl, aic, bic, skew, ekurt
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        real(dp) :: h0

        call fit_nagarch_range(y, range_var, 2*max_iter, gtol, fopt, params, niter, converged)
        persist = nagarch_persist(params)
        h0 = nagarch_range_h0(y, range_var, params, persist)
        vol_ann = sqrt(trading_days * h0) * 100.0_dp
        logl = -real(size(y), dp) * fopt
        aic = 2.0_dp*5.0_dp - 2.0_dp*logl
        bic = log(real(size(y), dp))*5.0_dp - 2.0_dp*logl
        call nagarch_range_skew_kurt(y, range_var, params, skew, ekurt)
    end subroutine fit_nagarch_range_piece

    subroutine fit_rgarch_piece(y, range_var, params, fopt, persist, vol_ann, logl, aic, bic, &
                                niter, converged, skew, ekurt)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: fopt, persist, vol_ann, logl, aic, bic, skew, ekurt
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        real(dp) :: h0

        call fit_rgarch(y, range_var, 2*max_iter, gtol, fopt, params, niter, converged)
        persist = rgarch_persist(params)
        h0 = rgarch_h0(y, range_var, params)
        vol_ann = sqrt(trading_days * h0) * 100.0_dp
        logl = -real(size(y), dp) * fopt
        aic = 2.0_dp*3.0_dp - 2.0_dp*logl
        bic = log(real(size(y), dp))*3.0_dp - 2.0_dp*logl
        call rgarch_skew_kurt(y, range_var, params, skew, ekurt)
    end subroutine fit_rgarch_piece

    subroutine nagarch_vol_forecast(y, params, persist, vol)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h, sqrth, r
        integer :: t

        h = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            vol(t) = sqrt(trading_days*h) * 100.0_dp
            r = y(t) - params%theta*sqrth
            h = params%omega + params%alpha*r**2 + params%beta*h
        end do
    end subroutine nagarch_vol_forecast

    subroutine rgarch_vol_forecast(y, range_var, params, vol)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h
        integer :: t

        h = rgarch_h0(y, range_var, params)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            vol(t) = sqrt(trading_days*h) * 100.0_dp
            h = params%omega + params%alpha*range_var(t) + params%beta*h
        end do
    end subroutine rgarch_vol_forecast

    subroutine nagarch_range_vol_forecast(y, range_var, params, persist, vol)
        real(dp), intent(in) :: y(:), range_var(:), persist
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h, sqrth, r
        integer :: t

        h = nagarch_range_h0(y, range_var, params, persist)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            vol(t) = sqrt(trading_days*h) * 100.0_dp
            r = y(t) - params%theta*sqrth
            h = params%omega + params%alpha*r**2 + params%gamma*range_var(t) + params%beta*h
        end do
    end subroutine nagarch_range_vol_forecast

    real(dp) function nagarch_range_h0(y, range_var, params, persist)
        real(dp), intent(in) :: y(:), range_var(:), persist
        type(garch_params_t), intent(in) :: params

        nagarch_range_h0 = max((params%omega + params%gamma*sum(range_var)/real(size(range_var), dp)) / &
            max(1.0_dp - persist, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function nagarch_range_h0

    real(dp) function rgarch_h0(y, range_var, params)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params

        rgarch_h0 = max((params%omega + params%alpha*sum(range_var)/real(size(range_var), dp)) / &
            max(1.0_dp - params%beta, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function rgarch_h0

    integer function rank_in_values(value, values)
        real(dp), intent(in) :: value, values(:)
        integer :: i

        rank_in_values = 1
        do i = 1, size(values)
            if (values(i) < value) rank_in_values = rank_in_values + 1
        end do
    end function rank_in_values

    subroutine print_piece(asset, piece, params, persist, vol_ann, logl, aic, bic, niter, converged, skew, ekurt, &
                           aic_rank, bic_rank)
        character(len=*), intent(in) :: asset, piece
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: persist, vol_ann, logl, aic, bic, skew, ekurt
        integer, intent(in) :: niter, aic_rank, bic_rank
        logical, intent(in) :: converged

        write(*,'(A10,1X,A12,ES12.3,5F9.4,F12.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3,2I9)') &
            trim(asset), trim(piece), params%omega, params%alpha, params%gamma, params%beta, params%theta, &
            persist, vol_ann, logl, aic, bic, niter, converged, skew, ekurt, aic_rank, bic_rank
    end subroutine print_piece

    subroutine print_ohlc_component_parameter_table()
        integer :: i

        write(*,'(/,A)') "OHLC component parameters:"
        write(*,'(A)') "Asset       Model          Component      omega    alpha    gamma     beta    theta  persist"
        write(*,'(A)') repeat("-", 92)
        do i = 1, nassets
            call print_component_params(asset_names(i), "OHLC_TOTAL", "CO", p_co_all(i), persist_co_all(i))
            call print_component_params(asset_names(i), "OHLC_TOTAL", "OC", p_oc_all(i), persist_oc_all(i))
            call print_component_params(asset_names(i), "OHLC_RGARCH", "CO", p_co_rg_all(i), persist_co_rg_all(i))
            call print_component_params(asset_names(i), "OHLC_RGARCH", "OC", p_oc_rg_all(i), persist_oc_rg_all(i))
            call print_component_params(asset_names(i), "OHLC_RANGE", "CO", p_co_rng_all(i), persist_co_rng_all(i))
            call print_component_params(asset_names(i), "OHLC_RANGE", "OC", p_oc_rng_all(i), persist_oc_rng_all(i))
        end do
    end subroutine print_ohlc_component_parameter_table

    subroutine print_component_params(asset, model, component, params, persist)
        character(len=*), intent(in) :: asset, model, component
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: persist

        write(*,'(A10,1X,A12,1X,A9,ES12.3,5F9.4)') trim(asset), trim(model), trim(component), &
            params%omega, params%alpha, params%gamma, params%beta, params%theta, persist
    end subroutine print_component_params

end program xfit_nagarch_ohlc_returns

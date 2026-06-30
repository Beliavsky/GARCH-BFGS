! Driver module for comparing NAGARCH and NAGARCH-X on daily CC returns.

module nagarchx_compare_mod
    use kind_mod,       only: dp
    use date_mod,       only: date_label
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file, filter_intraday_session
    use intraday_realized_measures_mod, only: daily_realized_panel_t, build_daily_realized_panel
    use garch_types_mod,  only: garch_params_t
    use garch_fit_mod,    only: fit_nagarch, nagarch_persist
    use garch_forecast_mod, only: nagarch_variance_path
    use stats_mod,        only: mean, demean_first, gaussian_loglik, correlation
    use program_utils_mod,  only: elapsed_since
    use path_utils_mod,     only: asset_label
    use nagarchx_mod
    implicit none
    private

    real(dp), parameter :: min_var       = 1.0e-12_dp
    integer,  parameter :: default_ntest = 250
    integer,  parameter :: max_iter      = 500
    real(dp), parameter :: gtol          = 1.0e-5_dp
    real(dp), parameter :: trading_days  = 252.0_dp

    type, public :: nagarchx_result_t
        character(len=8) :: asset     = ""
        real(dp) :: logL_te1  = 0.0_dp, logL_te2  = 0.0_dp
        real(dp) :: qlike1    = 0.0_dp, qlike2    = 0.0_dp
        real(dp) :: qlike_rv1 = 0.0_dp, qlike_rv2 = 0.0_dp
        real(dp) :: corr_rv1  = 0.0_dp, corr_rv2  = 0.0_dp
    end type nagarchx_result_t

    public :: run_nagarchx_compare

contains

    subroutine run_nagarchx_compare(filenames)
        character(len=*), intent(in) :: filenames(:)
        type(nagarchx_result_t), allocatable :: results(:)
        integer :: i

        allocate(results(size(filenames)))
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call compare_one_file(trim(filenames(i)), results(i))
        end do
        if (size(filenames) > 1) call print_asset_summary(results)
        deallocate(results)
    end subroutine run_nagarchx_compare

    subroutine compare_one_file(filename, result)
        character(len=*),        intent(in)  :: filename
        type(nagarchx_result_t), intent(out) :: result
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_realized_panel_t) :: daily
        real(dp), allocatable :: ret_cc(:), rv(:), h_nag(:), h_nagx(:)
        type(garch_params_t)    :: p_nag
        type(nagarchx_params_t) :: p_nagx, p_warm
        integer  :: nobs, ntrain, test_first, te_last
        integer  :: niter_nag, niter_nagx
        logical  :: conv_nag, conv_nagx
        real(dp) :: fopt_nag, fopt_nagx, persist_nag, persist_nagx
        real(dp) :: t0, t1, read_sec, fit_sec, elapsed_sec
        real(dp) :: mu_y2, mu_rv

        call cpu_time(t0)
        call cpu_time(t1)
        call read_intraday_prices_file(filename, bars)
        call filter_intraday_session(bars, regular_bars)
        read_sec = elapsed_since(t1)

        call build_daily_realized_panel(regular_bars, daily)
        nobs = size(daily%date) - 1
        if (nobs < default_ntest + 30) error stop "xfit_nagarchx_returns: not enough daily observations"
        ntrain     = nobs - default_ntest
        test_first = ntrain + 1
        te_last    = nobs

        allocate(ret_cc(nobs), rv(nobs), h_nag(nobs), h_nagx(nobs))
        ret_cc = log(daily%close(2:nobs+1) / daily%close(1:nobs))
        rv     = daily%rv(1:nobs)
        call demean_first(ret_cc, ntrain)

        print '(A,A)',          "Input file:   ", filename
        print '(A,I0,A,A,A,A)', "Daily obs:    ", nobs, &
              " from ", date_label(daily%date(1)), " to ", date_label(daily%date(nobs+1))
        print '(A,I0,A,A,A,A)', "Training:     ", ntrain, &
              " through ", date_label(daily%date(ntrain+1)), &
              ";  test starts ", date_label(daily%date(ntrain+2))
        print '(A,I0)',          "Test obs:     ", default_ntest
        print '(A,ES10.3,A,ES10.3)', "Mean(y^2):    ", mean(ret_cc(1:ntrain)**2), &
              "  Mean(RV): ", mean(rv(1:ntrain))

        mu_y2 = max(mean(ret_cc(1:ntrain)**2), min_var)
        mu_rv  = max(mean(rv(1:ntrain)),         min_var)

        call cpu_time(t1)
        call fit_nagarch(ret_cc(1:ntrain), max_iter, gtol, fopt_nag, p_nag, niter_nag, conv_nag)
        call nagarch_variance_path(ret_cc, p_nag, h_nag)
        persist_nag = nagarch_persist(p_nag)

        p_warm = nagarchx_params_t(omega = p_nag%omega, &
                                   alpha = p_nag%alpha, &
                                   beta  = p_nag%beta,  &
                                   theta = p_nag%theta, &
                                   delta = 0.05_dp * mu_y2 / mu_rv)
        call fit_nagarchx(ret_cc, rv, ntrain, max_iter, gtol, fopt_nagx, p_nagx, &
                          niter_nagx, conv_nagx, warm=p_warm)
        call nagarchx_variance_path(ret_cc, rv, p_nagx, h_nagx)
        persist_nagx = nagarchx_persist(p_nagx)

        fit_sec     = elapsed_since(t1)
        elapsed_sec = elapsed_since(t0)

        print '(/,A)', "NAGARCH parameters:"
        print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,F7.4,A,I0,A,L1)', &
              "omega=", p_nag%omega,  "  alpha=", p_nag%alpha, &
              "  beta=",  p_nag%beta,  "  theta=", p_nag%theta, &
              "  persist=", persist_nag, "  iter=", niter_nag, "  conv=", conv_nag
        print '(A)', "NAGARCH-X parameters:"
        print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,ES10.3,A,F7.4,A,I0,A,L1)', &
              "omega=", p_nagx%omega, "  alpha=", p_nagx%alpha, &
              "  beta=",  p_nagx%beta,  "  theta=", p_nagx%theta, &
              "  delta=", p_nagx%delta, "  persist=", persist_nagx, &
              "  iter=", niter_nagx, "  conv=", conv_nagx

        ! Compute all test-set metrics for this asset.
        result%asset     = asset_label(filename)
        result%logL_te1  = gaussian_loglik(ret_cc(test_first:te_last), h_nag(test_first:te_last))
        result%logL_te2  = gaussian_loglik(ret_cc(test_first:te_last), h_nagx(test_first:te_last))
        result%qlike1    = mean(log(max(h_nag (test_first:te_last), min_var)) + &
                                ret_cc(test_first:te_last)**2 / max(h_nag (test_first:te_last), min_var))
        result%qlike2    = mean(log(max(h_nagx(test_first:te_last), min_var)) + &
                                ret_cc(test_first:te_last)**2 / max(h_nagx(test_first:te_last), min_var))
        result%qlike_rv1 = mean(log(max(h_nag (test_first:te_last), min_var)) + &
                                rv(test_first:te_last) / max(h_nag (test_first:te_last), min_var))
        result%qlike_rv2 = mean(log(max(h_nagx(test_first:te_last), min_var)) + &
                                rv(test_first:te_last) / max(h_nagx(test_first:te_last), min_var))
        result%corr_rv1  = correlation(h_nag (test_first:te_last), rv(test_first:te_last))
        result%corr_rv2  = correlation(h_nagx(test_first:te_last), rv(test_first:te_last))

        call print_comparison(ret_cc, h_nag, h_nagx, ntrain, default_ntest, &
                              niter_nag, niter_nagx, conv_nag, conv_nagx, result)

        print '(/,A,F7.3,A,F7.3,A,F7.3)', &
              "read_sec=", read_sec, "  fit_sec=", fit_sec, "  elapsed_sec=", elapsed_sec

        deallocate(ret_cc, rv, h_nag, h_nagx)
    end subroutine compare_one_file

    subroutine print_comparison(y, h1, h2, ntrain, ntest, niter1, niter2, conv1, conv2, result)
        real(dp),                intent(in) :: y(:), h1(:), h2(:)
        integer,                 intent(in) :: ntrain, ntest, niter1, niter2
        logical,                 intent(in) :: conv1, conv2
        type(nagarchx_result_t), intent(in) :: result
        real(dp) :: logL_tr1, logL_tr2, aic1, bic1, aic2, bic2, vol1, vol2, logn
        integer  :: test_first, te_last

        test_first = ntrain + 1
        te_last    = ntrain + ntest
        logL_tr1 = gaussian_loglik(y(1:ntrain), h1(1:ntrain))
        logL_tr2 = gaussian_loglik(y(1:ntrain), h2(1:ntrain))
        logn = log(real(ntrain, dp))
        aic1 = 2.0_dp*4 - 2.0_dp*logL_tr1
        aic2 = 2.0_dp*5 - 2.0_dp*logL_tr2
        bic1 = logn*4   - 2.0_dp*logL_tr1
        bic2 = logn*5   - 2.0_dp*logL_tr2
        vol1 = 100.0_dp * sqrt(trading_days * mean(h1(test_first:te_last)))
        vol2 = 100.0_dp * sqrt(trading_days * mean(h2(test_first:te_last)))

        print '(/,A)', "Daily CC return forecast comparison (AIC/BIC from training logL)"
        print '(A,I0,A,I0)', "Training obs: ", ntrain, "  test obs: ", ntest
        print '(A)', repeat("-", 108)
        print '(A10,1X,A3,2(3X,A12),2(3X,A11),3X,A9,3X,A8,3X,A4,1X,A4)', &
              "Model", "k", "logL_train", "logL_test", "AIC", "BIC", "QLIKE", "vol_ann%", "iter", "conv"
        print '(A)', repeat("-", 108)
        print '(A10,1X,I3,2(3X,F12.3),2(3X,F11.3),3X,F9.5,3X,F8.3,3X,I4,1X,L1)', &
              "NAGARCH",   4, logL_tr1, result%logL_te1, aic1, bic1, result%qlike1, vol1, niter1, conv1
        print '(A10,1X,I3,2(3X,F12.3),2(3X,F11.3),3X,F9.5,3X,F8.3,3X,I4,1X,L1)', &
              "NAGARCH-X", 5, logL_tr2, result%logL_te2, aic2, bic2, result%qlike2, vol2, niter2, conv2
        print '(A)', repeat("-", 108)
        if (result%logL_te1 >= result%logL_te2) then
            print '(A)', "Best test logL: NAGARCH"
        else
            print '(A)', "Best test logL: NAGARCH-X"
        end if

        print '(/,A)', "RV prediction (test set, h_t vs RV_t)"
        print '(A)', repeat("-", 46)
        print '(A10,3X,A9,3X,A10)', "Model", "QLIKE_RV", "Corr(h,RV)"
        print '(A)', repeat("-", 46)
        print '(A10,3X,F9.5,3X,F10.6)', "NAGARCH",   result%qlike_rv1, result%corr_rv1
        print '(A10,3X,F9.5,3X,F10.6)', "NAGARCH-X", result%qlike_rv2, result%corr_rv2
        print '(A)', repeat("-", 46)
        if (result%qlike_rv1 <= result%qlike_rv2) then
            print '(A)', "Best QLIKE_RV: NAGARCH"
        else
            print '(A)', "Best QLIKE_RV: NAGARCH-X"
        end if
    end subroutine print_comparison

    subroutine print_asset_summary(results)
        type(nagarchx_result_t), intent(in) :: results(:)
        integer :: i, n
        integer :: win_logL, win_qlike, win_qlike_rv, win_corr_rv

        n            = size(results)
        win_logL     = count(results%logL_te2  > results%logL_te1)
        win_qlike    = count(results%qlike2    < results%qlike1)
        win_qlike_rv = count(results%qlike_rv2 < results%qlike_rv1)
        win_corr_rv  = count(results%corr_rv2  > results%corr_rv1)

        print '(/,/,A,I0,A)', "Summary across ", n, " assets"

        print '(A)', "CC return forecasting metrics (test set)"
        print '(A)', repeat("-", 58)
        print '(A8,2X,A11,2X,A12,2X,A9,2X,A10)', &
              "Asset", "logL_te_NAG", "logL_te_NAGX", "QLIKE_NAG", "QLIKE_NAGX"
        print '(A)', repeat("-", 58)
        do i = 1, n
            print '(A8,2X,F11.3,2X,F12.3,2X,F9.5,2X,F10.5)', &
                  results(i)%asset, results(i)%logL_te1, results(i)%logL_te2, &
                  results(i)%qlike1, results(i)%qlike2
        end do
        print '(A)', repeat("-", 58)
        write(*, '(A,I0,A,I0)') "NAGARCH-X wins logL_test: ", win_logL,  " / ", n
        write(*, '(A,I0,A,I0)') "NAGARCH-X wins QLIKE:     ", win_qlike, " / ", n

        print '(/,A)', "RV prediction metrics (test set, h_t vs RV_t)"
        print '(A)', repeat("-", 64)
        print '(A8,2X,A12,2X,A13,2X,A11,2X,A12)', &
              "Asset", "QLIKE_RV_NAG", "QLIKE_RV_NAGX", "Corr_RV_NAG", "Corr_RV_NAGX"
        print '(A)', repeat("-", 64)
        do i = 1, n
            print '(A8,2X,F12.5,2X,F13.5,2X,F11.6,2X,F12.6)', &
                  results(i)%asset, results(i)%qlike_rv1, results(i)%qlike_rv2, &
                  results(i)%corr_rv1, results(i)%corr_rv2
        end do
        print '(A)', repeat("-", 64)
        write(*, '(A,I0,A,I0)') "NAGARCH-X wins QLIKE_RV:  ", win_qlike_rv, " / ", n
        write(*, '(A,I0,A,I0)') "NAGARCH-X wins Corr_RV:   ", win_corr_rv,  " / ", n
    end subroutine print_asset_summary

end module nagarchx_compare_mod

! Fit standard symmetric GARCH(1,1)-Normal to the ETF return panel.
! The program reads spy_efa_eem_tlt_lqd.csv, uses demeaned log returns, and
! reports fitted parameters, likelihood diagnostics, and residual moments.

program xfit_symm_garch_returns
    use kind_mod,  only: dp
    use csv_mod,   only: read_price_csv, print_price_sample_info
    use stats_mod, only: mean
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod, only: fit_symm_garch, garch_skew_kurt, symm_garch_persist
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp

    integer,           allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp),          allocatable :: prices(:,:), ret(:)

    integer :: nprices, ncols, nobs, icol
    real(dp) :: ret_mean, fopt
    real(dp) :: persist, h_unc, vol_ann, logl, aic, bic, skew, ekurt
    type(garch_params_t) :: params
    integer :: niter
    integer :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s
    logical :: converged

    call system_clock(clock_start, clock_rate)
    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A)') "Asset          omega    alpha     beta  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt"
    write(*,'(A)') repeat("-", 118)

    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean

        call fit_symm_garch(ret, max_iter, gtol, fopt, params, niter, converged)
        persist = symm_garch_persist(params)
        h_unc   = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
        vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        logl    = -real(nobs, dp) * fopt
        aic     = 2.0_dp * 3.0_dp - 2.0_dp * logl
        bic     = log(real(nobs, dp)) * 3.0_dp - 2.0_dp * logl
        call garch_skew_kurt(ret, params, skew, ekurt)

        write(*,'(A10,ES12.3,4F9.4,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3)') &
            trim(col_names(icol)), params%omega, params%alpha, params%beta, persist, vol_ann, &
            logl, aic, bic, niter, converged, skew, ekurt
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

end program xfit_symm_garch_returns

! Fit GARCH, GARCH-M, NAGARCH, and NAGARCH-M models to excess returns.
!
! The -M models include a conditional mean: exret_t = mu + lam*sqrt(h_t) + eps_t
! where lam is the price of volatility risk and mu absorbs any residual constant.
!
! The risk-free rate rf can be supplied as:
!   (1) nothing          -> named constant default_ann_rf (4% annual)
!   (2) a numeric scalar -> that annual rate (fraction, e.g. 0.05 for 5%)
!   (3) a CSV filename   -> two-column file: YYYYMMDD integer, annual rf fraction
!                          dates are matched to return dates (nearest-previous rule)
!
! Usage: xfit_garch_m_returns [prices_csv] [rf_arg] [-simple|-log]
!   prices_csv : default spy_efa_eem_tlt_lqd.csv
!   rf_arg     : numeric constant or CSV filename (default: 4% annual)
!   -simple    : use simple returns (P_t/P_{t-1} - 1) instead of log returns
!   -log       : use log returns (default)

program xfit_garch_m_returns
    use kind_mod,        only: dp
    use date_mod,        only: print_program_header
    use stats_mod,       only: mean, variance
    use csv_mod,         only: read_price_csv, print_price_sample_info, read_rf_csv, nearest_previous_rf
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,   only: fit_symm_garch, fit_nagarch, fit_garch_m, fit_nagarch_m, &
                               nagarch_persist, symm_garch_persist
    implicit none

    ! ── named constants ──────────────────────────────────────────────────────
    character(len=*), parameter :: default_prices_file  = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: default_ann_rf   = 0.04_dp   ! 4% annual risk-free rate
    logical,  parameter :: default_log_returns = .true.  ! use log returns by default
    real(dp), parameter :: trading_days     = 252.0_dp
    integer,  parameter :: max_iter        = 500
    real(dp), parameter :: gtol            = 1.0e-7_dp

    character(len=10), parameter :: model_names(4) = &
        [character(len=10) :: "SYMM_GARCH", "GARCH_M", "NAGARCH", "NAGARCH_M"]
    integer, parameter :: nmodel = 4
    integer, parameter :: nparams(nmodel) = [3, 5, 4, 6]

    ! ── local variables ──────────────────────────────────────────────────────
    integer,  allocatable :: dates(:), rf_dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), rf_series(:)
    real(dp), allocatable :: ret(:), exret(:), rf_daily(:)
    type(garch_params_t) :: params(nmodel)
    real(dp) :: fopt(nmodel), persist, h_unc, vol_ann, logL, aic, bic
    integer  :: niter(nmodel)
    logical  :: converged(nmodel)

    character(len=256) :: prices_file, rf_arg, arg
    real(dp) :: ann_rf, rf_try
    integer  :: nprices, ncols, nobs, icol, imod, ios, narg, iarg, iret
    logical  :: rf_is_file, log_returns
    character(len=256) :: hdr

    call print_program_header("xfit_garch_m_returns.f90")

    ! ── parse command-line arguments ─────────────────────────────────────────
    narg        = command_argument_count()
    prices_file = default_prices_file
    rf_arg      = ""
    log_returns = default_log_returns

    do iarg = 1, narg
        call get_command_argument(iarg, arg)
        select case (trim(arg))
        case ("-simple")
            log_returns = .false.
        case ("-log")
            log_returns = .true.
        case default
            if (iarg == 1) then
                prices_file = arg
            else if (iarg == 2 .or. len_trim(rf_arg) == 0) then
                rf_arg = arg
            end if
        end select
    end do

    ! determine rf source
    rf_is_file = .false.
    ann_rf     = default_ann_rf
    if (len_trim(rf_arg) > 0) then
        read(rf_arg, *, iostat=ios) rf_try
        if (ios == 0) then
            ann_rf = rf_try                      ! numeric constant supplied
        else
            rf_is_file = .true.                  ! treat as filename
        end if
    end if

    ! ── load prices ──────────────────────────────────────────────────────────
    call read_price_csv(trim(prices_file), dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), exret(nobs), rf_daily(nobs))

    ! ── load or construct daily rf vector ────────────────────────────────────
    if (rf_is_file) then
        call read_rf_csv(trim(rf_arg), rf_dates, rf_series)
        ! align rf_series to return dates (dates(2:nprices)) using nearest-previous rule
        do iret = 1, nobs
            rf_daily(iret) = nearest_previous_rf(dates(iret + 1), rf_dates, rf_series) / trading_days
        end do
        write(*, '(A,A)') "Risk-free rate: time series from ", trim(rf_arg)
    else
        rf_daily = ann_rf / trading_days
        if (len_trim(rf_arg) == 0) then
            write(*, '(A,F5.2,A)') "Risk-free rate: default constant ", ann_rf*100.0_dp, "% annual"
        else
            write(*, '(A,F5.2,A)') "Risk-free rate: supplied constant ", ann_rf*100.0_dp, "% annual"
        end if
    end if

    call print_price_sample_info(trim(prices_file), dates, ncols)

    ! ── column header ────────────────────────────────────────────────────────
    hdr = "Model       Asset     omega      alpha       beta     theta    lambda        mu" // &
          "   persist  vol_ann%        logL         AIC         BIC  iter conv"
    write(*, '(A)') trim(hdr)
    write(*, '(A)') repeat("-", len_trim(hdr))

    write(*, '(A)') "Return type: " // merge("log   ", "simple", log_returns)

    ! ── fit models for each asset ─────────────────────────────────────────────
    do icol = 1, ncols
        if (log_returns) then
            ret = log(prices(2:nprices, icol) / prices(1:nprices-1, icol))
        else
            ret = prices(2:nprices, icol) / prices(1:nprices-1, icol) - 1.0_dp
        end if
        exret = ret - rf_daily
        exret = exret - mean(exret)      ! demean excess returns for all models

        call fit_symm_garch(exret, max_iter, gtol, fopt(1), params(1), niter(1), converged(1))
        call fit_garch_m   (exret, max_iter, gtol, fopt(2), params(2), niter(2), converged(2))
        call fit_nagarch   (exret, max_iter, gtol, fopt(3), params(3), niter(3), converged(3))
        call fit_nagarch_m (exret, max_iter, gtol, fopt(4), params(4), niter(4), converged(4))

        do imod = 1, nmodel
            select case (imod)
            case (1, 2)
                persist = symm_garch_persist(params(imod))
            case (3, 4)
                persist = nagarch_persist(params(imod))
            end select
            h_unc   = params(imod)%omega / max(1.0_dp - persist, 1.0e-8_dp)
            vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
            logL    = -real(nobs, dp) * fopt(imod)
            aic     = 2.0_dp*real(nparams(imod), dp) - 2.0_dp*logL
            bic     = log(real(nobs, dp))*real(nparams(imod), dp) - 2.0_dp*logL

            write(*, '(A10,1X,A9,ES11.3,5F9.4,F10.2,F10.2,2F12.2,F12.2,I5,1X,L1)') &
                trim(model_names(imod)), trim(col_names(icol)), &
                params(imod)%omega, params(imod)%alpha, params(imod)%beta, &
                params(imod)%theta, params(imod)%gamma, params(imod)%twist, &
                persist, vol_ann, logL, aic, bic, niter(imod), converged(imod)
        end do
        write(*, '(A)') ""
    end do

end program xfit_garch_m_returns

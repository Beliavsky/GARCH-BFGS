! Fit AR(0..p_max) models to excess returns by OLS and select order by AIC/BIC.
!
! Mean equation: r_t - rf_t = mu + phi_1*(r_{t-1}-rf_{t-1}) + ... + phi_p*(r_{t-p}-rf_{t-p}) + eps_t
! Estimation: exact MLE = OLS (Gaussian residuals assumed).
! Information criteria: AIC = 2k - 2*logL,  BIC = k*log(neff) - 2*logL,  k = p+1.
!
! Risk-free rate rf:
!   (1) nothing -> named constant default_ann_rf (4% annual)
!   (2) -rf numeric -> that annual rate as a fraction (e.g. 0.05 for 5%)
!   (3) -rf filename -> two-column CSV: YYYYMMDD integer, annual rf fraction
!
! Usage: xfit_ar_returns [prices_csv] [-p pmax] [-rf value_or_file] [-simple|-log]
!   prices_csv : default spy_efa_eem_tlt_lqd.csv
!   -p pmax    : maximum AR order to consider (default 10)
!   -rf ...    : risk-free rate constant or time-series file
!   -simple    : use simple returns (P_t/P_{t-1} - 1)
!   -log       : use log returns log(P_t/P_{t-1}) (default)

program xfit_ar_returns
    use kind_mod,       only: dp
    use date_mod,       only: print_program_header
    use stats_mod,      only: mean, fit_ar_ols
    use csv_mod,        only: read_price_csv, print_price_sample_info, read_rf_csv, nearest_previous_rf
    implicit none

    ! ── named constants ──────────────────────────────────────────────────────
    character(len=*), parameter :: default_prices_file  = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: default_ann_rf     = 0.04_dp   ! 4% annual risk-free rate
    logical,  parameter :: default_log_returns = .true.   ! log returns by default
    integer,  parameter :: default_p_max       = 10       ! maximum AR order

    real(dp), parameter :: trading_days = 252.0_dp

    ! ── local variables ──────────────────────────────────────────────────────
    integer,  allocatable :: dates(:), rf_dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), rf_series(:)
    real(dp), allocatable :: ret(:), exret(:), rf_daily(:)
    real(dp), allocatable :: betas(:,:)          ! betas(k, p+1) = AR(p) coefficients
    real(dp), allocatable :: aic_v(:), bic_v(:), logl_v(:), sigma_v(:)

    character(len=256) :: prices_file, rf_arg, arg
    real(dp) :: ann_rf, rf_try, sigma2, logl_val
    integer  :: nprices, ncols, nobs, icol, p, p_max, ios, narg, iarg, iret
    integer  :: p_aic, p_bic
    logical  :: rf_is_file, log_returns
    character(len=2) :: sel

    call print_program_header("xfit_ar_returns.f90")

    ! ── defaults ─────────────────────────────────────────────────────────────
    prices_file = default_prices_file
    rf_arg      = ""
    log_returns = default_log_returns
    p_max       = default_p_max
    ann_rf      = default_ann_rf
    rf_is_file  = .false.

    ! ── parse command-line arguments ─────────────────────────────────────────
    narg = command_argument_count()
    iarg = 0
    do while (iarg < narg)
        iarg = iarg + 1
        call get_command_argument(iarg, arg)
        select case (trim(arg))
        case ("-simple");  log_returns = .false.
        case ("-log");     log_returns = .true.
        case ("-p")
            iarg = iarg + 1
            call get_command_argument(iarg, arg)
            read(arg, *, iostat=ios) p_max
            if (ios /= 0 .or. p_max < 0) error stop "invalid p_max"
        case ("-rf")
            iarg = iarg + 1
            call get_command_argument(iarg, rf_arg)
        case default
            if (iarg == 1) prices_file = arg
        end select
    end do

    ! ── resolve risk-free rate ────────────────────────────────────────────────
    if (len_trim(rf_arg) > 0) then
        read(rf_arg, *, iostat=ios) rf_try
        if (ios == 0) then
            ann_rf = rf_try
        else
            rf_is_file = .true.
        end if
    end if

    ! ── load prices ───────────────────────────────────────────────────────────
    call read_price_csv(trim(prices_file), dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), exret(nobs), rf_daily(nobs))
    allocate(betas(p_max+1, p_max+1))
    allocate(aic_v(p_max+1), bic_v(p_max+1), logl_v(p_max+1), sigma_v(p_max+1))

    ! ── build daily rf vector ─────────────────────────────────────────────────
    if (rf_is_file) then
        call read_rf_csv(trim(rf_arg), rf_dates, rf_series)
        do iret = 1, nobs
            rf_daily(iret) = nearest_previous_rf(dates(iret+1), rf_dates, rf_series) / trading_days
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

    write(*, '(A)') "Return type: " // merge("log   ", "simple", log_returns)
    write(*, '(A,I0)') "Maximum AR order: ", p_max
    call print_price_sample_info(trim(prices_file), dates, ncols)

    ! ── fit and select for each asset ─────────────────────────────────────────
    do icol = 1, ncols
        if (log_returns) then
            ret = log(prices(2:nprices, icol) / prices(1:nprices-1, icol))
        else
            ret = prices(2:nprices, icol) / prices(1:nprices-1, icol) - 1.0_dp
        end if
        exret = ret - rf_daily

        do p = 0, p_max
            call fit_ar_ols(exret, nobs, p, betas(1:p+1, p+1), sigma2, logl_val)
            sigma_v(p+1) = sqrt(sigma2) * 100.0_dp
            logl_v(p+1)  = logl_val
            aic_v(p+1)   = 2.0_dp*real(p+1, dp) - 2.0_dp*logl_val
            bic_v(p+1)   = log(real(nobs-p, dp))*real(p+1, dp) - 2.0_dp*logl_val
        end do

        p_aic = minloc(aic_v(1:p_max+1), 1) - 1
        p_bic = minloc(bic_v(1:p_max+1), 1) - 1

        write(*, '(/,A,A,A,I0,A)') &
            "Asset: ", trim(col_names(icol)), "  (n=", nobs, " excess returns)"
        write(*, '(2X,A)') "p    neff   sigma%      logL       AIC       BIC  sel"
        write(*, '(2X,A)') repeat("-", 55)
        do p = 0, p_max
            sel = "  "
            if (p == p_aic) sel(1:1) = "A"
            if (p == p_bic) sel(2:2) = "B"
            write(*, '(2X,I2,I8,F9.4,F11.2,2F10.2,2X,A2)') &
                p, nobs-p, sigma_v(p+1), logl_v(p+1), aic_v(p+1), bic_v(p+1), sel
        end do
        call print_selected("AIC", p_aic, betas)
        call print_selected("BIC", p_bic, betas)
    end do

contains

    subroutine print_selected(ic_name, p_sel, b)
        character(len=*), intent(in) :: ic_name
        integer,          intent(in) :: p_sel
        real(dp),         intent(in) :: b(:,:)
        integer :: i
        write(*, '(4X,A,A,I0,A,F9.5)', advance='no') &
            ic_name, " selected p=", p_sel, "  mu=", b(1, p_sel+1)
        do i = 1, p_sel
            write(*, '(A,I0,A,F8.4)', advance='no') "  phi_", i, "=", b(i+1, p_sel+1)
        end do
        write(*, *)
    end subroutine print_selected

end program xfit_ar_returns

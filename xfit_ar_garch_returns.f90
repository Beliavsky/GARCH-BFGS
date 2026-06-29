! Fit AR(p)-GARCH(1,1), AR(p)-NAGARCH(1,1), and AR(p)-GJR-GARCH(1,1) to excess returns.
!
! Mean equation: r_t = mu + phi_1*r_{t-1} + ... + phi_p*r_{t-p} + eps_t
! GARCH:   h_t = omega + alpha*eps_{t-1}^2 + beta*h_{t-1}
! NAGARCH: h_t = omega + alpha*(eps_{t-1} - theta*sqrt(h_{t-1}))^2 + beta*h_{t-1}
! GJR:     h_t = omega + (alpha + gamma*I(eps_{t-1}<0))*eps_{t-1}^2 + beta*h_{t-1}
!
! Noise distribution selected by -noise NORMAL (default) or -noise T.
! For Student-t fits, degrees of freedom nu are stored in params%extra1.
!
! Risk-free rate rf:
!   (1) nothing -> named constant default_ann_rf (4% annual)
!   (2) -rf numeric -> that annual rate as a fraction (e.g. 0.05 for 5%)
!   (3) -rf filename -> two-column CSV: YYYYMMDD integer, annual rf fraction
!
! Usage: xfit_ar_garch_returns [prices_csv] [-p ar_order] [-rf value_or_file]
!                               [-simple|-log] [-noise NORMAL|T]

program xfit_ar_garch_returns
    use kind_mod,        only: dp
    use date_mod,        only: print_program_header
    use stats_mod,       only: mean
    use csv_mod,         only: read_price_csv, print_price_sample_info, read_rf_csv, nearest_previous_rf
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,   only: fit_ar_garch, fit_ar_nagarch, fit_ar_gjr, &
                               fit_ar_garch_t, fit_ar_nagarch_t, fit_ar_gjr_t, &
                               symm_garch_persist, nagarch_persist, gjr_persist
    implicit none

    ! ── named constants ──────────────────────────────────────────────────────
    character(len=*), parameter :: default_prices_file  = "spy_efa_eem_tlt_lqd.csv"
    character(len=*), parameter :: default_noise        = "NORMAL"
    real(dp), parameter :: default_ann_rf     = 0.04_dp   ! 4% annual risk-free rate
    logical,  parameter :: default_log_returns = .true.   ! log returns by default
    integer,  parameter :: default_ar_p        = 3        ! AR order in mean equation
    real(dp), parameter :: trading_days        = 252.0_dp
    integer,  parameter :: max_iter            = 500
    real(dp), parameter :: gtol                = 1.0e-7_dp

    character(len=12), parameter :: model_names(3) = &
        [character(len=12) :: "AR_GARCH", "AR_NAGARCH", "AR_GJR"]
    integer, parameter :: nmodel = 3

    ! ── local variables ──────────────────────────────────────────────────────
    integer,  allocatable :: dates(:), rf_dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), rf_series(:)
    real(dp), allocatable :: ret(:), exret(:), rf_daily(:)
    type(garch_params_t) :: params(nmodel)
    real(dp) :: fopt(nmodel), persist, h_unc, vol_ann, logL, aic, bic, nu_val
    integer  :: niter(nmodel), nparam(nmodel)
    logical  :: converged(nmodel)

    character(len=256) :: prices_file, rf_arg, arg
    character(len=8)   :: noise
    real(dp) :: ann_rf, rf_try
    integer  :: nprices, ncols, nobs, icol, imod, ios, narg, iarg, iret
    integer  :: ar_p, i
    logical  :: rf_is_file, log_returns
    character(len=256) :: hdr

    call print_program_header("xfit_ar_garch_returns.f90")

    ! ── defaults ─────────────────────────────────────────────────────────────
    prices_file = default_prices_file
    noise       = default_noise
    rf_arg      = ""
    log_returns = default_log_returns
    ar_p        = default_ar_p
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
            read(arg, *, iostat=ios) ar_p
            if (ios /= 0 .or. ar_p < 0) error stop "invalid AR order"
        case ("-rf")
            iarg = iarg + 1
            call get_command_argument(iarg, rf_arg)
        case ("-noise")
            iarg = iarg + 1
            call get_command_argument(iarg, arg)
            noise = trim(arg)
            if (noise /= "NORMAL" .and. noise /= "T") error stop "noise must be NORMAL or T"
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
    write(*, '(A,I0)') "AR order in mean equation: ", ar_p
    write(*, '(A,A)') "Noise distribution: ", trim(noise)
    call print_price_sample_info(trim(prices_file), dates, ncols)

    ! ── nparam for each model ─────────────────────────────────────────────────
    ! GARCH/N=ar_p+4, GARCH/T=ar_p+5; NAGARCH/N=ar_p+5, NAGARCH/T=ar_p+6
    ! GJR/N=ar_p+5, GJR/T=ar_p+6
    if (trim(noise) == "T") then
        nparam(1) = ar_p + 5
        nparam(2) = ar_p + 6
        nparam(3) = ar_p + 6
    else
        nparam(1) = ar_p + 4
        nparam(2) = ar_p + 5
        nparam(3) = ar_p + 5
    end if

    ! ── header ───────────────────────────────────────────────────────────────
    hdr = "Model          Asset     omega      alpha     gamma       beta     theta        mu" // &
          "   persist  vol_ann%        nu        logL         AIC         BIC  iter conv"
    write(*, '(A)') trim(hdr)
    write(*, '(A)') repeat("-", len_trim(hdr))

    ! ── fit for each asset ────────────────────────────────────────────────────
    do icol = 1, ncols
        if (log_returns) then
            ret = log(prices(2:nprices, icol) / prices(1:nprices-1, icol))
        else
            ret = prices(2:nprices, icol) / prices(1:nprices-1, icol) - 1.0_dp
        end if
        exret = ret - rf_daily

        if (trim(noise) == "T") then
            call fit_ar_garch_t  (exret, ar_p, max_iter, gtol, fopt(1), params(1), niter(1), converged(1))
            call fit_ar_nagarch_t(exret, ar_p, max_iter, gtol, fopt(2), params(2), niter(2), converged(2))
            call fit_ar_gjr_t    (exret, ar_p, max_iter, gtol, fopt(3), params(3), niter(3), converged(3))
        else
            call fit_ar_garch  (exret, ar_p, max_iter, gtol, fopt(1), params(1), niter(1), converged(1))
            call fit_ar_nagarch(exret, ar_p, max_iter, gtol, fopt(2), params(2), niter(2), converged(2))
            call fit_ar_gjr    (exret, ar_p, max_iter, gtol, fopt(3), params(3), niter(3), converged(3))
        end if

        do imod = 1, nmodel
            select case (imod)
            case (1);  persist = symm_garch_persist(params(imod))
            case (2);  persist = nagarch_persist(params(imod))
            case (3);  persist = gjr_persist(params(imod))
            end select
            h_unc   = params(imod)%omega / max(1.0_dp - persist, 1.0e-8_dp)
            vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
            logL    = -real(nobs, dp) * fopt(imod)
            aic     = 2.0_dp*real(nparam(imod), dp) - 2.0_dp*logL
            bic     = log(real(nobs, dp))*real(nparam(imod), dp) - 2.0_dp*logL
            nu_val  = params(imod)%extra1   ! 0 for Normal, dof for Student-t

            write(*, '(A12,2X,A9,ES11.3,5F9.4,2F10.2,F9.2,3F12.2,I5,1X,L1)') &
                trim(model_names(imod)), trim(col_names(icol)), &
                params(imod)%omega, params(imod)%alpha, params(imod)%gamma, &
                params(imod)%beta, params(imod)%theta, params(imod)%twist, &
                persist, vol_ann, nu_val, logL, aic, bic, niter(imod), converged(imod)
        end do

        ! print AR coefficients for each model
        if (ar_p > 0) then
            do imod = 1, nmodel
                write(*, '(4X,A,A)', advance='no') trim(model_names(imod)), " AR coefs: mu="
                write(*, '(F8.5)', advance='no') params(imod)%twist
                if (allocated(params(imod)%ar_coefs)) then
                    do i = 1, ar_p
                        write(*, '(A,I0,A,F8.4)', advance='no') &
                            "  phi_", i, "=", params(imod)%ar_coefs(i)
                    end do
                end if
                write(*, *)
            end do
        end if
        write(*, *)
    end do

end program xfit_ar_garch_returns

! Fit time-varying parameter AR(p) model (TVP-AR) to excess returns via Kalman filter MLE.
!
! State space model:
!   Observation: r_t = Z_t*x_t + eps_t,  Z_t = (1, r_{t-1}, ..., r_{t-p}),  eps_t ~ N(0,sigma^2)
!   State:       x_t = c + F*x_{t-1} + eta_t,  F = diag(rho_0..rho_p),  eta_t ~ N(0,Q)
!
! x_t = (mu_t, phi_{1,t}, ..., phi_{p,t})' — p+1 time-varying coefficients
! Q   = diag(q_0, ..., q_p),  c_i = (1 - rho_i)*x_bar_i  (mean-reverting)
!
! Parameters (3*(p+1)+1 total, raw unconstrained):
!   x_bar(0:p)  : long-run state means (unconstrained)
!   rho_raw(0:p): persistence, rho_i = tanh(rho_raw_i) in (-1,1)
!   log_q(0:p)  : log state noise variances, q_i = exp(log_q_i)
!   log_sigma2  : log observation noise variance, sigma^2 = exp(log_sigma2)
!
! Warm-started from OLS AR(p) fit; gradient via central finite differences.
!
! Usage: xfit_tvp_ar_returns [prices_csv] [-p ar_order] [-rf value_or_file] [-simple|-log]

module tvp_ar_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    implicit none
    private

    public :: tvp_ar_set_data, tvp_ar_nll, tvp_ar_nll_f, tvp_ar_filter

    real(dp), allocatable, save :: tvp_obs(:)
    integer,  save :: tvp_nobs = 0
    integer,  save :: tvp_p    = 0

contains

    subroutine tvp_ar_set_data(y, p)
        ! Store time series and AR order for use by tvp_ar_obj and tvp_ar_filter.
        real(dp), intent(in) :: y(:)  ! return series (length nobs)
        integer,  intent(in) :: p     ! AR order
        tvp_obs  = y
        tvp_nobs = size(y)
        tvp_p    = p
    end subroutine tvp_ar_set_data

    subroutine tvp_ar_nll(p_raw, nll)
        ! Total Kalman-filter negative log-likelihood for TVP-AR(p).
        real(dp), intent(in)  :: p_raw(:)
        real(dp), intent(out) :: nll
        integer  :: m, t, i, j
        real(dp) :: x_bar_v(tvp_p+1), rho(tvp_p+1), q(tvp_p+1), sigma2
        real(dp) :: x_filt(tvp_p+1), x_pred(tvp_p+1)
        real(dp) :: P_filt(tvp_p+1, tvp_p+1), P_pred(tvp_p+1, tvp_p+1)
        real(dp) :: Z(tvp_p+1), Pz(tvp_p+1)
        real(dp) :: v, S
        real(dp), parameter :: S_floor = 1.0e-14_dp

        m       = tvp_p + 1
        x_bar_v = p_raw(1:m)
        rho     = tanh(p_raw(m+1:2*m))
        q       = exp(p_raw(2*m+1:3*m))
        sigma2  = exp(p_raw(3*m+1))

        ! Stationary initial state
        x_filt   = x_bar_v
        P_filt   = 0.0_dp
        do i = 1, m
            P_filt(i,i) = q(i) / max(1.0_dp - rho(i)**2, 1.0e-8_dp)
        end do

        nll = 0.0_dp
        do t = m, tvp_nobs    ! first observation at t = p+1
            ! --- predict ---
            do i = 1, m
                x_pred(i) = (1.0_dp - rho(i))*x_bar_v(i) + rho(i)*x_filt(i)
            end do
            do i = 1, m
                do j = 1, m
                    P_pred(i,j) = rho(i)*rho(j)*P_filt(i,j)
                end do
                P_pred(i,i) = P_pred(i,i) + q(i)
            end do

            ! --- observation vector Z_t = (1, r_{t-1}, ..., r_{t-p}) ---
            Z(1) = 1.0_dp
            do i = 1, tvp_p
                Z(i+1) = tvp_obs(t-i)
            end do

            ! --- scalar innovation v = r_t - Z_t * x_pred ---
            v = tvp_obs(t)
            do i = 1, m
                v = v - Z(i)*x_pred(i)
            end do

            ! --- Pz = P_pred * Z ---
            do i = 1, m
                Pz(i) = 0.0_dp
                do j = 1, m
                    Pz(i) = Pz(i) + P_pred(i,j)*Z(j)
                end do
            end do

            ! --- innovation variance S = Z' * Pz + sigma2 ---
            S = sigma2
            do i = 1, m
                S = S + Z(i)*Pz(i)
            end do
            S = max(S, S_floor)

            nll = nll + log_sqrt_2pi + 0.5_dp*(log(S) + v**2/S)

            ! --- update: gain K = Pz/S ---
            do i = 1, m
                x_filt(i) = x_pred(i) + Pz(i)/S * v
            end do
            do i = 1, m
                do j = 1, m
                    P_filt(i,j) = P_pred(i,j) - Pz(i)*Pz(j)/S
                end do
            end do
        end do
    end subroutine tvp_ar_nll

    subroutine tvp_ar_nll_f(p_raw, np, f)
        ! Wrapper matching obj_f_iface for bfgs_minimize_fd.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p_raw(np)
        real(dp), intent(out) :: f
        call tvp_ar_nll(p_raw(1:np), f)
    end subroutine tvp_ar_nll_f

    subroutine tvp_ar_filter(p_raw, x_filt_all, neff)
        ! Full filtered state path after optimization (duplicate pass for reporting).
        real(dp), intent(in)  :: p_raw(:)
        real(dp), allocatable, intent(out) :: x_filt_all(:,:)  ! (neff, p+1)
        integer,  intent(out) :: neff
        integer  :: m, t, i, j, it
        real(dp) :: x_bar_v(tvp_p+1), rho(tvp_p+1), q(tvp_p+1), sigma2
        real(dp) :: x_filt(tvp_p+1), x_pred(tvp_p+1)
        real(dp) :: P_filt(tvp_p+1, tvp_p+1), P_pred(tvp_p+1, tvp_p+1)
        real(dp) :: Z(tvp_p+1), Pz(tvp_p+1)
        real(dp) :: v, S
        real(dp), parameter :: S_floor = 1.0e-14_dp

        m    = tvp_p + 1
        neff = tvp_nobs - tvp_p
        allocate(x_filt_all(neff, m))

        x_bar_v = p_raw(1:m)
        rho     = tanh(p_raw(m+1:2*m))
        q       = exp(p_raw(2*m+1:3*m))
        sigma2  = exp(p_raw(3*m+1))

        x_filt = x_bar_v
        P_filt = 0.0_dp
        do i = 1, m
            P_filt(i,i) = q(i) / max(1.0_dp - rho(i)**2, 1.0e-8_dp)
        end do

        it = 0
        do t = m, tvp_nobs
            do i = 1, m
                x_pred(i) = (1.0_dp - rho(i))*x_bar_v(i) + rho(i)*x_filt(i)
            end do
            do i = 1, m
                do j = 1, m
                    P_pred(i,j) = rho(i)*rho(j)*P_filt(i,j)
                end do
                P_pred(i,i) = P_pred(i,i) + q(i)
            end do
            Z(1) = 1.0_dp
            do i = 1, tvp_p
                Z(i+1) = tvp_obs(t-i)
            end do
            v = tvp_obs(t)
            do i = 1, m
                v = v - Z(i)*x_pred(i)
            end do
            do i = 1, m
                Pz(i) = 0.0_dp
                do j = 1, m
                    Pz(i) = Pz(i) + P_pred(i,j)*Z(j)
                end do
            end do
            S = sigma2
            do i = 1, m
                S = S + Z(i)*Pz(i)
            end do
            S = max(S, S_floor)
            do i = 1, m
                x_filt(i) = x_pred(i) + Pz(i)/S * v
            end do
            do i = 1, m
                do j = 1, m
                    P_filt(i,j) = P_pred(i,j) - Pz(i)*Pz(j)/S
                end do
            end do
            it = it + 1
            x_filt_all(it, :) = x_filt
        end do
    end subroutine tvp_ar_filter

end module tvp_ar_mod

! ─────────────────────────────────────────────────────────────────────────────

program xfit_tvp_ar_returns
    use kind_mod,       only: dp
    use date_mod,       only: print_program_header
    use stats_mod,      only: fit_ar_ols
    use csv_mod,        only: read_price_csv, print_price_sample_info, read_rf_csv, nearest_previous_rf
    use bfgs_mod,       only: bfgs_minimize_fd
    use tvp_ar_mod,     only: tvp_ar_set_data, tvp_ar_nll_f, tvp_ar_filter
    implicit none

    ! ── named constants ──────────────────────────────────────────────────────
    character(len=*), parameter :: default_prices_file  = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: default_ann_rf     = 0.04_dp
    logical,  parameter :: default_log_returns = .true.
    integer,  parameter :: default_ar_p        = 3
    real(dp), parameter :: trading_days        = 252.0_dp
    integer,  parameter :: max_iter            = 2000
    real(dp), parameter :: gtol_fd             = 1.0e-4_dp   ! loose: matches numerical gradient noise

    ! ── local variables ──────────────────────────────────────────────────────
    integer,  allocatable :: dates(:), rf_dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), rf_series(:)
    real(dp), allocatable :: ret(:), exret(:), rf_daily(:)
    real(dp), allocatable :: x_filt_all(:,:)
    real(dp), allocatable :: p0(:), p_best(:), beta_ols(:)
    real(dp), allocatable :: x_bar_v(:), rho_v(:), q_v(:), mean_filt(:)

    character(len=256) :: prices_file, rf_arg, arg
    real(dp) :: ann_rf, rf_try, sigma2_ols, logl_ols, nll_best, logL, aic, bic, sigma_obs
    integer  :: nprices, ncols, nobs, neff, icol, ios, narg, iarg, iret
    integer  :: ar_p, m, np, niter_best, i
    logical  :: rf_is_file, log_returns, converged_best

    call print_program_header("xfit_tvp_ar_returns.f90")

    ! ── defaults ─────────────────────────────────────────────────────────────
    prices_file = default_prices_file
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

    m  = ar_p + 1       ! number of state components (intercept + p AR lags)
    np = 3*m + 1        ! total parameters

    allocate(p0(np), p_best(np), beta_ols(m))
    allocate(x_bar_v(m), rho_v(m), q_v(m), mean_filt(m))

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
    write(*, '(A,I0)') "AR order: ", ar_p
    write(*, '(A,I0)') "Parameters: ", np
    call print_price_sample_info(trim(prices_file), dates, ncols)

    ! ── fit for each asset ────────────────────────────────────────────────────
    do icol = 1, ncols
        if (log_returns) then
            ret = log(prices(2:nprices, icol) / prices(1:nprices-1, icol))
        else
            ret = prices(2:nprices, icol) / prices(1:nprices-1, icol) - 1.0_dp
        end if
        exret = ret - rf_daily
        neff  = nobs - ar_p

        ! --- warm start from OLS AR(p): moderate process noise, high persistence ---
        call fit_ar_ols(exret, nobs, ar_p, beta_ols, sigma2_ols, logl_ols)
        p0(1:m)       = beta_ols                              ! x_bar from OLS
        p0(m+1:2*m)   = atanh(0.95_dp)                       ! rho_raw -> rho=0.95
        p0(2*m+1:3*m) = log(max(sigma2_ols * 0.01_dp, 1.0e-16_dp))   ! q = 1% of obs var
        p0(3*m+1)     = log(max(sigma2_ols, 1.0e-16_dp))     ! log(sigma^2)

        call tvp_ar_set_data(exret, ar_p)
        p_best = p0
        call bfgs_minimize_fd(tvp_ar_nll_f, p_best, np, max_iter, gtol_fd, nll_best, niter_best, converged_best)

        ! --- decode optimized parameters ---
        x_bar_v  = p_best(1:m)
        rho_v    = tanh(p_best(m+1:2*m))
        q_v      = exp(p_best(2*m+1:3*m))
        sigma_obs = sqrt(exp(p_best(3*m+1)))
        logL     = -nll_best
        aic      = 2.0_dp*real(np, dp) - 2.0_dp*logL
        bic      = log(real(neff, dp))*real(np, dp) - 2.0_dp*logL

        ! --- filtered path for mean/final coefficients ---
        call tvp_ar_filter(p_best, x_filt_all, neff)
        do i = 1, m
            mean_filt(i) = sum(x_filt_all(:, i)) / real(neff, dp)
        end do

        ! --- output ---
        write(*, '(/,A,A,A,I0,A,I0,A)') &
            "Asset: ", trim(col_names(icol)), "  (nobs=", nobs, "  neff=", neff, ")"
        write(*, '(A,F12.2,A,F12.2,A,F12.2,A,I5,A,L1)') &
            "TVP-AR(" // trim(itoa(ar_p)) // ")  logL=", logL, &
            "  AIC=", aic, "  BIC=", bic, "  iter=", niter_best, "  conv=", converged_best

        write(*, '(/,2X,A)') "Component    x_bar          rho        sigma_q"
        write(*, '(2X,A)') repeat("-", 55)
        write(*, '(2X,A12,3ES14.5)') "mu", x_bar_v(1), rho_v(1), sqrt(q_v(1))
        do i = 1, ar_p
            write(*, '(2X,A6,I0,A5,3ES14.5)') "phi_", i, "     ", x_bar_v(i+1), rho_v(i+1), sqrt(q_v(i+1))
        end do
        write(*, '(2X,A,ES14.5)') "sigma_obs:  ", sigma_obs

        write(*, '(/,4X,A)', advance='no') "Mean filtered:  mu="
        write(*, '(ES11.4)', advance='no') mean_filt(1)
        do i = 1, ar_p
            write(*, '(A,I0,A,ES11.4)', advance='no') "  phi_", i, "=", mean_filt(i+1)
        end do
        write(*, *)

        write(*, '(4X,A)', advance='no') "Final filtered: mu="
        write(*, '(ES11.4)', advance='no') x_filt_all(neff, 1)
        do i = 1, ar_p
            write(*, '(A,I0,A,ES11.4)', advance='no') "  phi_", i, "=", x_filt_all(neff, i+1)
        end do
        write(*, *)

        deallocate(x_filt_all)
    end do

contains

    pure function itoa(n) result(s)
        integer, intent(in) :: n
        character(len=12) :: s
        write(s, '(I0)') n
        s = adjustl(s)
    end function itoa

end program xfit_tvp_ar_returns

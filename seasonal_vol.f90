! Calendar-day annual seasonal volatility tests.
!
! The core test regresses log((r_t - mean(r))^2 + eps) on an intercept and
! annual Fourier terms evaluated from the calendar day-of-year.  It reports the
! usual nested-model F statistic for the joint null that all seasonal sine/cosine
! coefficients are zero.  This first version assumes calendar-day observations;
! business-day returns can reuse the Fourier basis later with trading-calendar
! controls added separately.

module seasonal_vol_mod
    use kind_mod, only: dp
    use stats_mod, only: mean, variance
    use linalg_mod, only: chol_factor, chol_solve_vec
    implicit none
    private

    public :: seasonal_vol_result_t, test_annual_vol_seasonality_calendar

    type :: seasonal_vol_result_t
        integer :: nobs = 0
        integer :: nharm = 0
        integer :: df_num = 0
        integer :: df_den = 0
        real(dp) :: rss_null = 0.0_dp
        real(dp) :: rss_full = 0.0_dp
        real(dp) :: r2 = 0.0_dp
        real(dp) :: delta_aic = 0.0_dp
        real(dp) :: delta_bic = 0.0_dp
        real(dp) :: f_stat = 0.0_dp
        real(dp) :: p_value = 1.0_dp
        real(dp) :: seasonal_vol_avg = 0.0_dp
        real(dp) :: seasonal_vol_sd = 0.0_dp
        real(dp) :: seasonal_vol_min = 0.0_dp
        real(dp) :: seasonal_vol_max = 0.0_dp
        character(len=5) :: seasonal_vol_min_date = ""
        character(len=5) :: seasonal_vol_max_date = ""
        logical :: ok = .false.
    end type seasonal_vol_result_t

contains

    subroutine test_annual_vol_seasonality_calendar(dates, returns, nharm, result, annualization)
        integer, intent(in) :: dates(:)
        real(dp), intent(in) :: returns(:)
        integer, intent(in) :: nharm
        type(seasonal_vol_result_t), intent(out) :: result
        real(dp), optional, intent(in) :: annualization
        real(dp), allocatable :: y(:), x(:,:), beta(:)
        real(dp) :: ret_mean, eps, rss_null, rss_full, ann
        integer :: n, k, q, t

        n = size(returns)
        k = 1 + 2*nharm
        q = 2*nharm
        result%nobs = n
        result%nharm = nharm
        result%df_num = q
        result%df_den = n - k

        if (size(dates) /= n .or. nharm < 1 .or. n <= k) return

        allocate(y(n), x(n,k), beta(k))
        ret_mean = mean(returns)
        eps = max(1.0e-12_dp*variance(returns), 1.0e-300_dp)
        do t = 1, n
            y(t) = log((returns(t) - ret_mean)**2 + eps)
        end do
        call annual_fourier_design(dates, nharm, x)

        rss_null = sum((y - mean(y))**2)
        call ols_fit(x, y, beta, rss_full, result%ok)
        if (.not. result%ok .or. rss_full <= 0.0_dp) return

        result%rss_null = rss_null
        result%rss_full = rss_full
        result%r2 = max(0.0_dp, 1.0_dp - rss_full / max(rss_null, 1.0e-300_dp))
        result%delta_aic = real(n, dp)*log(rss_full / rss_null) + 2.0_dp*real(q, dp)
        result%delta_bic = real(n, dp)*log(rss_full / rss_null) + log(real(n, dp))*real(q, dp)
        result%f_stat = ((rss_null - rss_full) / real(q, dp)) / (rss_full / real(n - k, dp))
        result%f_stat = max(result%f_stat, 0.0_dp)
        result%p_value = f_survival(result%f_stat, q, n - k)
        ann = 365.0_dp
        if (present(annualization)) ann = annualization
        call seasonal_vol_profile(beta, nharm, ann, variance(returns), result)
    end subroutine test_annual_vol_seasonality_calendar

    subroutine seasonal_vol_profile(beta, nharm, annualization, target_variance, result)
        real(dp), intent(in) :: beta(:)
        real(dp), intent(in) :: annualization, target_variance
        integer, intent(in) :: nharm
        type(seasonal_vol_result_t), intent(inout) :: result
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: raw_var(365), vols(365), phase, log_var, scale, vmean
        integer :: doy, h, col, imin, imax

        do doy = 1, 365
            phase = 2.0_dp*pi*real(doy - 1, dp) / 365.0_dp
            log_var = beta(1)
            do h = 1, nharm
                col = 2*h
                log_var = log_var + beta(col)*cos(real(h, dp)*phase) + beta(col + 1)*sin(real(h, dp)*phase)
            end do
            raw_var(doy) = exp(log_var)
        end do

        scale = max(target_variance, 0.0_dp) / max(sum(raw_var) / real(size(raw_var), dp), 1.0e-300_dp)
        vols = sqrt(max(scale*raw_var, 0.0_dp) * annualization) * 100.0_dp

        vmean = sum(vols) / real(size(vols), dp)
        result%seasonal_vol_avg = vmean
        result%seasonal_vol_sd = sqrt(sum((vols - vmean)**2) / real(size(vols) - 1, dp))
        imin = minloc(vols, dim=1)
        imax = maxloc(vols, dim=1)
        result%seasonal_vol_min = vols(imin)
        result%seasonal_vol_max = vols(imax)
        result%seasonal_vol_min_date = doy_label(imin)
        result%seasonal_vol_max_date = doy_label(imax)
    end subroutine seasonal_vol_profile

    subroutine annual_fourier_design(dates, nharm, x)
        integer, intent(in) :: dates(:)
        integer, intent(in) :: nharm
        real(dp), intent(out) :: x(:,:)
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: phase, year_len
        integer :: t, h, col, doy, year

        x(:,1) = 1.0_dp
        do t = 1, size(dates)
            call split_yyyymmdd(dates(t), year=year)
            doy = day_of_year(dates(t))
            year_len = merge(366.0_dp, 365.0_dp, is_leap_year(year))
            phase = 2.0_dp*pi*real(doy - 1, dp) / year_len
            do h = 1, nharm
                col = 2*h
                x(t,col) = cos(real(h, dp)*phase)
                x(t,col + 1) = sin(real(h, dp)*phase)
            end do
        end do
    end subroutine annual_fourier_design

    subroutine ols_fit(x, y, beta, rss, ok)
        real(dp), intent(in) :: x(:,:), y(:)
        real(dp), intent(out) :: beta(:), rss
        logical, intent(out) :: ok
        real(dp), allocatable :: xtx(:,:), xty(:), L(:,:), resid(:)
        integer :: n, k

        n = size(x, 1)
        k = size(x, 2)
        allocate(xtx(k,k), xty(k), L(k,k), resid(n))
        xtx = matmul(transpose(x), x)
        xty = matmul(transpose(x), y)
        call chol_factor(xtx, k, L, ok)
        if (.not. ok) then
            rss = huge(1.0_dp)
            beta = 0.0_dp
            return
        end if
        call chol_solve_vec(L, k, xty, beta)
        resid = y - matmul(x, beta)
        rss = sum(resid**2)
    end subroutine ols_fit

    integer function day_of_year(yyyymmdd)
        integer, intent(in) :: yyyymmdd
        integer, parameter :: month_days_common(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
        integer :: year, month, day

        call split_yyyymmdd(yyyymmdd, year, month, day)
        day_of_year = day
        if (month > 1) day_of_year = day_of_year + sum(month_days_common(1:month-1))
        if (month > 2 .and. is_leap_year(year)) day_of_year = day_of_year + 1
    end function day_of_year

    character(len=5) function doy_label(doy)
        integer, intent(in) :: doy
        integer, parameter :: month_days_common(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
        integer :: month, day, remaining

        remaining = doy
        month = 1
        do while (month < 12 .and. remaining > month_days_common(month))
            remaining = remaining - month_days_common(month)
            month = month + 1
        end do
        day = remaining
        write(doy_label, '(I2.2,A1,I2.2)') month, "-", day
    end function doy_label

    subroutine split_yyyymmdd(yyyymmdd, year, month, day)
        integer, intent(in) :: yyyymmdd
        integer, optional, intent(out) :: year, month, day
        integer :: yy, mm, dd

        yy = yyyymmdd / 10000
        mm = mod(yyyymmdd / 100, 100)
        dd = mod(yyyymmdd, 100)
        if (present(year)) year = yy
        if (present(month)) month = mm
        if (present(day)) day = dd
    end subroutine split_yyyymmdd

    logical function is_leap_year(year)
        integer, intent(in) :: year

        is_leap_year = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. mod(year, 400) == 0
    end function is_leap_year

    real(dp) function f_survival(f, df1, df2)
        real(dp), intent(in) :: f
        integer, intent(in) :: df1, df2
        real(dp) :: x

        if (f <= 0.0_dp) then
            f_survival = 1.0_dp
            return
        end if
        x = real(df2, dp) / (real(df2, dp) + real(df1, dp)*f)
        f_survival = regularized_beta(x, 0.5_dp*real(df2, dp), 0.5_dp*real(df1, dp))
        f_survival = min(max(f_survival, 0.0_dp), 1.0_dp)
    end function f_survival

    real(dp) function regularized_beta(x, a, b)
        real(dp), intent(in) :: x, a, b
        real(dp) :: bt

        if (x <= 0.0_dp) then
            regularized_beta = 0.0_dp
        else if (x >= 1.0_dp) then
            regularized_beta = 1.0_dp
        else
            bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp - x))
            if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
                regularized_beta = bt*beta_contfrac(x, a, b) / a
            else
                regularized_beta = 1.0_dp - bt*beta_contfrac(1.0_dp - x, b, a) / b
            end if
        end if
    end function regularized_beta

    real(dp) function beta_contfrac(x, a, b)
        real(dp), intent(in) :: x, a, b
        integer, parameter :: max_iter = 200
        real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
        real(dp) :: qab, qap, qam, c, d, h, aa, del
        integer :: m, m2

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp/d
        h = d
        do m = 1, max_iter
            m2 = 2*m
            aa = real(m, dp)*(b - real(m, dp))*x / ((qam + real(m2, dp))*(a + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            h = h*d*c
            aa = -(a + real(m, dp))*(qab + real(m, dp))*x / ((a + real(m2, dp))*(qap + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        beta_contfrac = h
    end function beta_contfrac

end module seasonal_vol_mod

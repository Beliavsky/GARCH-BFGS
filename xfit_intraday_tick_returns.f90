! Fit simple integer-valued models to 1-second close-to-close tick changes.
!
! Prices are converted to half-penny ticks:
!   x_t = round((close_t - close_{t-1}) / 0.005)
! Consecutive observations from different calendar dates are skipped.

program xfit_intraday_tick_returns
    use kind_mod, only: dp
    use date_mod, only: yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv
    implicit none

    real(dp), parameter :: tick_size = 0.005_dp
    real(dp), parameter :: min_prob = 1.0e-300_dp
    real(dp), parameter :: sqrt2 = 1.41421356237309504880_dp
    integer, parameter :: coord_iter = 3
    integer, parameter :: golden_iter = 28
    character(len=256) :: filename
    type(ohlcv_series_t) :: series
    integer, allocatable :: ticks(:), values(:), counts(:)
    integer :: nobs, nticks, nval
    integer :: opt_mode = 0
    real(dp) :: opt_sigma = 1.0_dp, opt_nu = 5.0_dp, opt_p0 = 0.0_dp

    filename = "c:\python\databento\spy_1s_databento.csv"
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    call read_intraday_prices_csv(trim(filename), series)
    call close_tick_changes(series, ticks, nticks)
    call tick_frequencies(ticks, values, counts, nval)
    nobs = sum(counts)

    print '(A)', "Intraday integer tick-return model comparison"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,F8.4)', "Tick size: ", tick_size
    print '(A,I0)', "Close-to-close tick changes: ", nobs
    print '(A,I0)', "Unique tick changes: ", nval
    print '(A)', ""
    call print_model_table(values, counts, nobs)

contains

    ! Convert within-day close-to-close changes to integer tick changes.
    subroutine close_tick_changes(series, ticks, nticks)
        type(ohlcv_series_t), intent(in) :: series
        integer, allocatable, intent(out) :: ticks(:)
        integer, intent(out) :: nticks
        integer :: i

        allocate(ticks(max(series%nobs() - 1, 0)))
        nticks = 0
        do i = 2, series%nobs()
            if (yyyymmdd(series%timestamp(i)%date) /= yyyymmdd(series%timestamp(i - 1)%date)) cycle
            nticks = nticks + 1
            ticks(nticks) = nint((series%close(i) - series%close(i - 1)) / tick_size)
        end do
        if (nticks < 1) error stop "close_tick_changes: no within-day close-to-close changes"
    end subroutine close_tick_changes

    ! Count frequencies of integer tick changes.
    subroutine tick_frequencies(ticks, values, counts, nval)
        integer, intent(in) :: ticks(:)
        integer, allocatable, intent(out) :: values(:), counts(:)
        integer, intent(out) :: nval
        integer :: i, j

        allocate(values(size(ticks)), counts(size(ticks)))
        nval = 0
        do i = 1, size(ticks)
            if (nval > 0) then
                j = find_value(values(1:nval), ticks(i))
            else
                j = 0
            end if
            if (j == 0) then
                nval = nval + 1
                values(nval) = ticks(i)
                counts(nval) = 1
            else
                counts(j) = counts(j) + 1
            end if
        end do
        call sort_values(values(1:nval), counts(1:nval))
    end subroutine tick_frequencies

    ! Print fitted model comparison for tick-change frequencies.
    subroutine print_model_table(values, counts, nobs)
        integer, intent(in) :: values(:), counts(:), nobs
        real(dp) :: sigma_n, sigma_t, nu_t, p0, scale_lap
        real(dp) :: ll

        print '(A)', "-------------------------------------------------------------------------------------------"
        print '(A28,1X,A5,1X,A10,1X,A10,1X,A10,1X,A14,1X,A14,1X,A14)', &
              "Model", "k", "sigma", "nu", "p_zero", "logL", "AIC", "BIC"
        print '(A)', "-------------------------------------------------------------------------------------------"

        sigma_n = fit_sigma_rounded_normal()
        ll = rounded_normal_loglik(values, counts, sigma_n)
        call print_row("rounded normal", 1, sigma_n, 0.0_dp, 0.0_dp, ll, nobs)

        call fit_zero_inflated_rounded_normal(sigma_n, p0)
        ll = zi_rounded_normal_loglik(values, counts, sigma_n, p0)
        call print_row("ZI rounded normal", 2, sigma_n, 0.0_dp, p0, ll, nobs)

        call fit_rounded_t(sigma_t, nu_t)
        ll = rounded_t_loglik(values, counts, sigma_t, nu_t)
        call print_row("rounded Student-t", 2, sigma_t, nu_t, 0.0_dp, ll, nobs)

        call fit_zero_inflated_rounded_t(sigma_t, nu_t, p0)
        ll = zi_rounded_t_loglik(values, counts, sigma_t, nu_t, p0)
        call print_row("ZI rounded Student-t", 3, sigma_t, nu_t, p0, ll, nobs)

        scale_lap = fit_discrete_laplace()
        ll = discrete_laplace_loglik(values, counts, scale_lap)
        call print_row("discrete Laplace", 1, scale_lap, 0.0_dp, 0.0_dp, ll, nobs)

        print '(A)', "-------------------------------------------------------------------------------------------"
        print '(A)', "sigma is in half-penny tick units. ZI means zero-inflated."
    end subroutine print_model_table

    ! Print one model comparison row.
    subroutine print_row(name, k, sigma, nu, p_zero, loglik, nobs)
        character(len=*), intent(in) :: name
        integer, intent(in) :: k, nobs
        real(dp), intent(in) :: sigma, nu, p_zero, loglik
        real(dp) :: aic, bic

        aic = -2.0_dp*loglik + 2.0_dp*real(k, dp)
        bic = -2.0_dp*loglik + log(real(nobs, dp))*real(k, dp)
        if (nu > 0.0_dp) then
            print '(A28,1X,I5,1X,F10.4,1X,F10.3,1X,F10.4,1X,F14.3,1X,F14.3,1X,F14.3)', &
                  trim(name), k, sigma, nu, p_zero, loglik, aic, bic
        else
            print '(A28,1X,I5,1X,F10.4,1X,A10,1X,F10.4,1X,F14.3,1X,F14.3,1X,F14.3)', &
                  trim(name), k, sigma, "-", p_zero, loglik, aic, bic
        end if
    end subroutine print_row

    ! Fit sigma for rounded normal likelihood.
    real(dp) function fit_sigma_rounded_normal() result(sigma)
        opt_mode = 1
        sigma = exp(maximize_1d(-4.0_dp, 5.0_dp))
    end function fit_sigma_rounded_normal

    ! Fit sigma and zero mass for rounded normal by coordinate search.
    subroutine fit_zero_inflated_rounded_normal(sigma, p0)
        real(dp), intent(out) :: sigma, p0
        integer :: iter

        sigma = fit_sigma_rounded_normal()
        p0 = 0.1_dp
        do iter = 1, coord_iter
            opt_mode = 2
            opt_p0 = p0
            sigma = exp(maximize_1d(log(max(sigma, 1.0e-6_dp)) - 3.0_dp, log(sigma) + 3.0_dp))
            opt_mode = 3
            opt_sigma = sigma
            p0 = logistic(maximize_1d(-12.0_dp, 12.0_dp))
        end do
    end subroutine fit_zero_inflated_rounded_normal

    ! Fit sigma and nu for rounded standardised Student-t.
    subroutine fit_rounded_t(sigma, nu)
        real(dp), intent(out) :: sigma, nu
        integer :: iter

        sigma = fit_sigma_rounded_normal()
        nu = 5.0_dp
        do iter = 1, coord_iter
            opt_mode = 4
            opt_nu = nu
            sigma = exp(maximize_1d(log(max(sigma, 1.0e-6_dp)) - 3.0_dp, log(sigma) + 3.0_dp))
            opt_mode = 5
            opt_sigma = sigma
            nu = nu_from_q(maximize_1d(-12.0_dp, 12.0_dp))
        end do
    end subroutine fit_rounded_t

    ! Fit sigma, nu, and zero mass for rounded standardised Student-t.
    subroutine fit_zero_inflated_rounded_t(sigma, nu, p0)
        real(dp), intent(out) :: sigma, nu, p0
        integer :: iter

        call fit_rounded_t(sigma, nu)
        p0 = 0.1_dp
        do iter = 1, coord_iter
            opt_mode = 6
            opt_nu = nu
            opt_p0 = p0
            sigma = exp(maximize_1d(log(max(sigma, 1.0e-6_dp)) - 3.0_dp, log(sigma) + 3.0_dp))
            opt_mode = 7
            opt_sigma = sigma
            opt_p0 = p0
            nu = nu_from_q(maximize_1d(-12.0_dp, 12.0_dp))
            opt_mode = 8
            opt_sigma = sigma
            opt_nu = nu
            p0 = logistic(maximize_1d(-12.0_dp, 12.0_dp))
        end do
    end subroutine fit_zero_inflated_rounded_t

    ! Fit symmetric discrete Laplace scale parameter.
    real(dp) function fit_discrete_laplace() result(scale)
        opt_mode = 9
        scale = exp(maximize_1d(-6.0_dp, 5.0_dp))
    end function fit_discrete_laplace

    ! Rounded normal log likelihood for frequency-counted ticks.
    real(dp) function rounded_normal_loglik(values, counts, sigma) result(ll)
        integer, intent(in) :: values(:), counts(:)
        real(dp), intent(in) :: sigma
        integer :: i
        real(dp) :: p

        ll = 0.0_dp
        do i = 1, size(values)
            p = normal_cdf((real(values(i), dp) + 0.5_dp)/sigma) - &
                normal_cdf((real(values(i), dp) - 0.5_dp)/sigma)
            ll = ll + real(counts(i), dp)*log(max(p, min_prob))
        end do
    end function rounded_normal_loglik

    ! Zero-inflated rounded normal log likelihood.
    real(dp) function zi_rounded_normal_loglik(values, counts, sigma, p0) result(ll)
        integer, intent(in) :: values(:), counts(:)
        real(dp), intent(in) :: sigma, p0
        integer :: i
        real(dp) :: p

        ll = 0.0_dp
        do i = 1, size(values)
            p = normal_cdf((real(values(i), dp) + 0.5_dp)/sigma) - &
                normal_cdf((real(values(i), dp) - 0.5_dp)/sigma)
            if (values(i) == 0) p = p0 + (1.0_dp - p0)*p
            if (values(i) /= 0) p = (1.0_dp - p0)*p
            ll = ll + real(counts(i), dp)*log(max(p, min_prob))
        end do
    end function zi_rounded_normal_loglik

    ! Rounded standardised Student-t log likelihood for frequency-counted ticks.
    real(dp) function rounded_t_loglik(values, counts, sigma, nu) result(ll)
        integer, intent(in) :: values(:), counts(:)
        real(dp), intent(in) :: sigma, nu
        integer :: i
        real(dp) :: p

        ll = 0.0_dp
        do i = 1, size(values)
            p = std_t_cdf((real(values(i), dp) + 0.5_dp)/sigma, nu) - &
                std_t_cdf((real(values(i), dp) - 0.5_dp)/sigma, nu)
            ll = ll + real(counts(i), dp)*log(max(p, min_prob))
        end do
    end function rounded_t_loglik

    ! Zero-inflated rounded Student-t log likelihood.
    real(dp) function zi_rounded_t_loglik(values, counts, sigma, nu, p0) result(ll)
        integer, intent(in) :: values(:), counts(:)
        real(dp), intent(in) :: sigma, nu, p0
        integer :: i
        real(dp) :: p

        ll = 0.0_dp
        do i = 1, size(values)
            p = std_t_cdf((real(values(i), dp) + 0.5_dp)/sigma, nu) - &
                std_t_cdf((real(values(i), dp) - 0.5_dp)/sigma, nu)
            if (values(i) == 0) p = p0 + (1.0_dp - p0)*p
            if (values(i) /= 0) p = (1.0_dp - p0)*p
            ll = ll + real(counts(i), dp)*log(max(p, min_prob))
        end do
    end function zi_rounded_t_loglik

    ! Symmetric discrete Laplace log likelihood with alpha=exp(-1/scale).
    real(dp) function discrete_laplace_loglik(values, counts, scale) result(ll)
        integer, intent(in) :: values(:), counts(:)
        real(dp), intent(in) :: scale
        real(dp) :: alpha, logc
        integer :: i

        alpha = exp(-1.0_dp / max(scale, 1.0e-8_dp))
        logc = log(1.0_dp - alpha) - log(1.0_dp + alpha)
        ll = 0.0_dp
        do i = 1, size(values)
            ll = ll + real(counts(i), dp)*(logc + abs(values(i))*log(alpha))
        end do
    end function discrete_laplace_loglik

    ! Maximize a scalar function on [a,b] with golden-section search.
    real(dp) function maximize_1d(a, b) result(xmax)
        real(dp), intent(in) :: a, b
        real(dp), parameter :: inv_phi = 0.61803398874989484820_dp
        real(dp) :: lo, hi, c, d, fc, fd
        integer :: iter

        lo = a
        hi = b
        c = hi - inv_phi*(hi - lo)
        d = lo + inv_phi*(hi - lo)
        fc = objective_1d(c)
        fd = objective_1d(d)
        do iter = 1, golden_iter
            if (fc < fd) then
                lo = c
                c = d
                fc = fd
                d = lo + inv_phi*(hi - lo)
                fd = objective_1d(d)
            else
                hi = d
                d = c
                fd = fc
                c = hi - inv_phi*(hi - lo)
                fc = objective_1d(c)
            end if
        end do
        xmax = 0.5_dp*(lo + hi)
    end function maximize_1d

    ! Objective selected by opt_mode for scalar coordinate optimization.
    real(dp) function objective_1d(q) result(ll)
        real(dp), intent(in) :: q

        select case (opt_mode)
        case (1)
            ll = rounded_normal_loglik(values, counts, exp(q))
        case (2)
            ll = zi_rounded_normal_loglik(values, counts, exp(q), opt_p0)
        case (3)
            ll = zi_rounded_normal_loglik(values, counts, opt_sigma, logistic(q))
        case (4)
            ll = rounded_t_loglik(values, counts, exp(q), opt_nu)
        case (5)
            ll = rounded_t_loglik(values, counts, opt_sigma, nu_from_q(q))
        case (6)
            ll = zi_rounded_t_loglik(values, counts, exp(q), opt_nu, opt_p0)
        case (7)
            ll = zi_rounded_t_loglik(values, counts, opt_sigma, nu_from_q(q), opt_p0)
        case (8)
            ll = zi_rounded_t_loglik(values, counts, opt_sigma, opt_nu, logistic(q))
        case (9)
            ll = discrete_laplace_loglik(values, counts, exp(q))
        case default
            error stop "objective_1d: opt_mode not set"
        end select
    end function objective_1d

    ! Standard normal CDF.
    pure real(dp) function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        p = 0.5_dp*(1.0_dp + erf(x / sqrt2))
    end function normal_cdf

    ! Unit-variance standardised Student-t CDF.
    real(dp) function std_t_cdf(x, nu) result(p)
        real(dp), intent(in) :: x, nu
        real(dp) :: scale

        scale = sqrt((nu - 2.0_dp) / nu)
        p = raw_t_cdf(x / scale, nu)
    end function std_t_cdf

    ! Ordinary Student-t CDF using the regularized incomplete beta function.
    real(dp) function raw_t_cdf(x, nu) result(p)
        real(dp), intent(in) :: x, nu
        real(dp) :: z, ib

        if (x == 0.0_dp) then
            p = 0.5_dp
            return
        end if
        z = nu / (nu + x*x)
        ib = regularized_beta(0.5_dp*nu, 0.5_dp, z)
        if (x > 0.0_dp) then
            p = 1.0_dp - 0.5_dp*ib
        else
            p = 0.5_dp*ib
        end if
    end function raw_t_cdf

    ! Regularized incomplete beta I_x(a,b).
    real(dp) function regularized_beta(a, b, x) result(bt)
        real(dp), intent(in) :: a, b, x
        real(dp) :: front

        if (x <= 0.0_dp) then
            bt = 0.0_dp
            return
        else if (x >= 1.0_dp) then
            bt = 1.0_dp
            return
        end if
        front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp - x))
        if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
            bt = front * beta_cf(a, b, x) / a
        else
            bt = 1.0_dp - front * beta_cf(b, a, 1.0_dp - x) / b
        end if
    end function regularized_beta

    ! Continued fraction for incomplete beta.
    real(dp) function beta_cf(a, b, x) result(cf)
        real(dp), intent(in) :: a, b, x
        integer, parameter :: max_iter = 200
        real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
        integer :: m, m2
        real(dp) :: aa, c, d, del, h, qab, qam, qap

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp / d
        h = d
        do m = 1, max_iter
            m2 = 2*m
            aa = real(m, dp)*(b - real(m, dp))*x / ((qam + real(m2, dp))*(a + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            h = h*d*c
            aa = -(a + real(m, dp))*(qab + real(m, dp))*x / ((a + real(m2, dp))*(qap + real(m2, dp)))
            d = 1.0_dp + aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        cf = h
    end function beta_cf

    ! Map real line to nu in (2.05, 100).
    pure real(dp) function nu_from_q(q) result(nu)
        real(dp), intent(in) :: q
        nu = 2.05_dp + 97.95_dp / (1.0_dp + exp(-q))
    end function nu_from_q

    ! Logistic transform.
    pure real(dp) function logistic(q) result(p)
        real(dp), intent(in) :: q
        p = 1.0_dp / (1.0_dp + exp(-q))
    end function logistic

    ! Return index of value in x, or zero if absent.
    pure integer function find_value(x, value) result(idx)
        integer, intent(in) :: x(:), value
        integer :: i

        idx = 0
        do i = 1, size(x)
            if (x(i) == value) then
                idx = i
                return
            end if
        end do
    end function find_value

    ! Sort integer values ascending with counts carried along.
    subroutine sort_values(x, counts)
        integer, intent(inout) :: x(:), counts(:)
        integer :: i, j, xv, cv

        do i = 2, size(x)
            xv = x(i)
            cv = counts(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= xv) exit
                x(j + 1) = x(j)
                counts(j + 1) = counts(j)
                j = j - 1
            end do
            x(j + 1) = xv
            counts(j + 1) = cv
        end do
    end subroutine sort_values
end program xfit_intraday_tick_returns

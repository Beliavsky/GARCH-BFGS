! Simple baseline models for intraday volatility forecasts.
!
! The baseline forecast has the form
!     Var(r_{d,j}) = daily_var_d * diurnal_var_j
! where daily_var_d is a daily variance forecast known before the session and
! diurnal_var_j is a normalized deterministic time-of-day variance multiplier.

module intraday_vol_baseline_mod
    use kind_mod, only: dp
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    public :: daily_variance_lag1
    public :: daily_variance_ewma
    public :: estimate_diurnal_variance_baseline
    public :: intraday_ewma_multiplier
    public :: parkinson_variance_proxy
    public :: garman_klass_variance_proxy
    public :: intraday_ewma_multiplier_from_proxy
    public :: intraday_variance_forecast
    public :: intraday_ewma_variance_forecast
    public :: fit_lag1_diurnal_baseline
    public :: fit_ewma_diurnal_baseline
    public :: fit_lag1_diurnal_intraday_ewma_baseline
    public :: fit_ewma_diurnal_intraday_ewma_baseline
    public :: remaining_session_var_diurnal

contains

    ! Forecast each intraday observation's daily variance using the prior day's realized variance.
    subroutine daily_variance_lag1(y, date_id, daily_var)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: daily_var(:)
        real(dp), allocatable :: day_rv(:), day_forecast(:)
        integer, allocatable :: obs_day(:)
        integer :: ndays, i, d
        real(dp) :: fallback

        call daily_realized_variance(y, date_id, day_rv, obs_day, ndays)
        allocate(day_forecast(ndays))
        fallback = max(sum(day_rv) / real(ndays, dp), min_var)
        day_forecast(1) = fallback
        do d = 2, ndays
            day_forecast(d) = max(day_rv(d - 1), min_var)
        end do
        do i = 1, size(y)
            daily_var(i) = day_forecast(obs_day(i))
        end do
        deallocate(day_rv, obs_day, day_forecast)
    end subroutine daily_variance_lag1

    ! Forecast each intraday observation's daily variance using an EWMA of prior daily realized variance.
    subroutine daily_variance_ewma(y, date_id, lambda, daily_var)
        real(dp), intent(in) :: y(:), lambda
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: daily_var(:)
        real(dp), allocatable :: day_rv(:), day_forecast(:)
        integer, allocatable :: obs_day(:)
        integer :: ndays, i, d
        real(dp) :: lam, fallback

        call daily_realized_variance(y, date_id, day_rv, obs_day, ndays)
        allocate(day_forecast(ndays))
        lam = min(max(lambda, 0.0_dp), 0.9999_dp)
        fallback = max(sum(day_rv) / real(ndays, dp), min_var)
        day_forecast(1) = fallback
        do d = 2, ndays
            day_forecast(d) = max(lam*day_forecast(d - 1) + (1.0_dp - lam)*day_rv(d - 1), min_var)
        end do
        do i = 1, size(y)
            daily_var(i) = day_forecast(obs_day(i))
        end do
        deallocate(day_rv, obs_day, day_forecast)
    end subroutine daily_variance_ewma

    ! Estimate normalized deterministic time-of-day variance multipliers.
    subroutine estimate_diurnal_variance_baseline(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(out) :: diurnal_var(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width
        real(dp), allocatable :: bin_sum(:), bin_curve(:)
        integer, allocatable :: bin_count(:)
        integer :: i, b, nbins, half_width
        logical :: do_smooth
        real(dp) :: fallback, norm

        call check_same_size(size(y), size(daily_var), size(bin_id), size(diurnal_var), &
                             "estimate_diurnal_variance_baseline")
        if (minval(bin_id) < 1) error stop "estimate_diurnal_variance_baseline: bin_id must be positive"
        nbins = maxval(bin_id)
        allocate(bin_sum(nbins), bin_count(nbins), bin_curve(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do i = 1, size(y)
            b = bin_id(i)
            bin_sum(b) = bin_sum(b) + y(i)**2 / max(daily_var(i), min_var)
            bin_count(b) = bin_count(b) + 1
        end do
        fallback = max(sum(bin_sum) / real(max(sum(bin_count), 1), dp), min_var)
        do b = 1, nbins
            if (bin_count(b) > 0) then
                bin_curve(b) = max(bin_sum(b) / real(bin_count(b), dp), min_var)
            else
                bin_curve(b) = fallback
            end if
        end do
        do_smooth = .false.
        if (present(smooth_diurnal)) do_smooth = smooth_diurnal
        half_width = 2
        if (present(smooth_half_width)) half_width = max(0, smooth_half_width)
        if (do_smooth) call smooth_log_curve(bin_curve, bin_count, half_width)
        do i = 1, size(y)
            diurnal_var(i) = max(bin_curve(bin_id(i)), min_var)
        end do
        norm = max(sum(diurnal_var) / real(size(diurnal_var), dp), min_var)
        diurnal_var = max(diurnal_var / norm, min_var)
        deallocate(bin_sum, bin_count, bin_curve)
    end subroutine estimate_diurnal_variance_baseline

    ! Estimate a within-day EWMA multiplier that resets to q=1 at each day's first observation.
    subroutine intraday_ewma_multiplier(y, date_id, daily_var, diurnal_var, lambda, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:), lambda
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: q(:)
        integer :: i
        real(dp) :: lam, z2

        if (size(y) /= size(date_id) .or. size(y) /= size(daily_var) .or. &
            size(y) /= size(diurnal_var) .or. size(q) /= size(y)) then
            error stop "intraday_ewma_multiplier: array sizes differ"
        end if
        lam = min(max(lambda, 0.0_dp), 0.9999_dp)
        q(1) = 1.0_dp
        do i = 2, size(y)
            if (date_id(i) /= date_id(i - 1)) then
                q(i) = 1.0_dp
            else
                z2 = y(i - 1)**2 / max(daily_var(i - 1) * diurnal_var(i - 1), min_var)
                q(i) = max(lam*q(i - 1) + (1.0_dp - lam)*z2, min_var)
            end if
        end do
    end subroutine intraday_ewma_multiplier

    ! Compute Parkinson high-low variance proxies for each OHLC bar.
    subroutine parkinson_variance_proxy(high_price, low_price, proxy)
        real(dp), intent(in) :: high_price(:), low_price(:)
        real(dp), intent(out) :: proxy(:)
        integer :: i

        if (size(high_price) /= size(low_price) .or. size(proxy) /= size(high_price)) then
            error stop "parkinson_variance_proxy: array sizes differ"
        end if
        do i = 1, size(proxy)
            proxy(i) = max(log(max(high_price(i), min_var) / max(low_price(i), min_var))**2 / &
                           (4.0_dp * log(2.0_dp)), min_var)
        end do
    end subroutine parkinson_variance_proxy

    ! Compute Garman-Klass OHLC variance proxies for each OHLC bar.
    subroutine garman_klass_variance_proxy(open_price, high_price, low_price, close_price, proxy)
        real(dp), intent(in) :: open_price(:), high_price(:), low_price(:), close_price(:)
        real(dp), intent(out) :: proxy(:)
        integer :: i
        real(dp) :: hl, co

        if (size(open_price) /= size(high_price) .or. size(open_price) /= size(low_price) .or. &
            size(open_price) /= size(close_price) .or. size(proxy) /= size(open_price)) then
            error stop "garman_klass_variance_proxy: array sizes differ"
        end if
        do i = 1, size(proxy)
            hl = log(max(high_price(i), min_var) / max(low_price(i), min_var))
            co = log(max(close_price(i), min_var) / max(open_price(i), min_var))
            proxy(i) = max(0.5_dp*hl**2 - (2.0_dp*log(2.0_dp) - 1.0_dp)*co**2, min_var)
        end do
    end subroutine garman_klass_variance_proxy

    ! Estimate a within-day EWMA multiplier from a supplied variance proxy.
    subroutine intraday_ewma_multiplier_from_proxy(proxy, date_id, daily_var, diurnal_var, lambda, q)
        real(dp), intent(in) :: proxy(:), daily_var(:), diurnal_var(:), lambda
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: q(:)
        integer :: i
        real(dp) :: lam, z2

        if (size(proxy) /= size(date_id) .or. size(proxy) /= size(daily_var) .or. &
            size(proxy) /= size(diurnal_var) .or. size(q) /= size(proxy)) then
            error stop "intraday_ewma_multiplier_from_proxy: array sizes differ"
        end if
        lam = min(max(lambda, 0.0_dp), 0.9999_dp)
        q(1) = 1.0_dp
        do i = 2, size(proxy)
            if (date_id(i) /= date_id(i - 1)) then
                q(i) = 1.0_dp
            else
                z2 = proxy(i - 1) / max(daily_var(i - 1) * diurnal_var(i - 1), min_var)
                q(i) = max(lam*q(i - 1) + (1.0_dp - lam)*z2, min_var)
            end if
        end do
    end subroutine intraday_ewma_multiplier_from_proxy

    ! Combine daily variance and diurnal multipliers into intraday variance forecasts.
    subroutine intraday_variance_forecast(daily_var, diurnal_var, h)
        real(dp), intent(in) :: daily_var(:), diurnal_var(:)
        real(dp), intent(out) :: h(:)

        if (size(daily_var) /= size(diurnal_var) .or. size(h) /= size(daily_var)) then
            error stop "intraday_variance_forecast: array sizes differ"
        end if
        h = max(daily_var * diurnal_var, min_var)
    end subroutine intraday_variance_forecast

    ! Combine daily variance, diurnal multipliers, and an intraday EWMA multiplier into h_t.
    subroutine intraday_ewma_variance_forecast(y, date_id, daily_var, diurnal_var, lambda, q, h)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:), lambda
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: q(:), h(:)

        call intraday_ewma_multiplier(y, date_id, daily_var, diurnal_var, lambda, q)
        if (size(h) /= size(q)) error stop "intraday_ewma_variance_forecast: array sizes differ"
        h = max(daily_var * diurnal_var * q, min_var)
    end subroutine intraday_ewma_variance_forecast

    ! Fit the lag-1 daily RV plus diurnal multiplier baseline and return h_t.
    subroutine fit_lag1_diurnal_baseline(y, date_id, bin_id, daily_var, diurnal_var, h, smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: date_id(:), bin_id(:)
        real(dp), intent(out) :: daily_var(:), diurnal_var(:), h(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call daily_variance_lag1(y, date_id, daily_var)
        call estimate_diurnal_variance_baseline(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        call intraday_variance_forecast(daily_var, diurnal_var, h)
    end subroutine fit_lag1_diurnal_baseline

    ! Fit the EWMA daily RV plus diurnal multiplier baseline and return h_t.
    subroutine fit_ewma_diurnal_baseline(y, date_id, bin_id, lambda, daily_var, diurnal_var, h, &
                                         smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), lambda
        integer, intent(in) :: date_id(:), bin_id(:)
        real(dp), intent(out) :: daily_var(:), diurnal_var(:), h(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call daily_variance_ewma(y, date_id, lambda, daily_var)
        call estimate_diurnal_variance_baseline(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        call intraday_variance_forecast(daily_var, diurnal_var, h)
    end subroutine fit_ewma_diurnal_baseline

    ! Fit lag-1 daily RV, diurnal multipliers, and a within-day EWMA multiplier.
    subroutine fit_lag1_diurnal_intraday_ewma_baseline(y, date_id, bin_id, intraday_lambda, daily_var, diurnal_var, q, h, &
                                                       smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), intraday_lambda
        integer, intent(in) :: date_id(:), bin_id(:)
        real(dp), intent(out) :: daily_var(:), diurnal_var(:), q(:), h(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call daily_variance_lag1(y, date_id, daily_var)
        call estimate_diurnal_variance_baseline(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        call intraday_ewma_variance_forecast(y, date_id, daily_var, diurnal_var, intraday_lambda, q, h)
    end subroutine fit_lag1_diurnal_intraday_ewma_baseline

    ! Fit EWMA daily RV, diurnal multipliers, and a within-day EWMA multiplier.
    subroutine fit_ewma_diurnal_intraday_ewma_baseline(y, date_id, bin_id, daily_lambda, intraday_lambda, &
                                                       daily_var, diurnal_var, q, h, smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_lambda, intraday_lambda
        integer, intent(in) :: date_id(:), bin_id(:)
        real(dp), intent(out) :: daily_var(:), diurnal_var(:), q(:), h(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call daily_variance_ewma(y, date_id, daily_lambda, daily_var)
        call estimate_diurnal_variance_baseline(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        call intraday_ewma_variance_forecast(y, date_id, daily_var, diurnal_var, intraday_lambda, q, h)
    end subroutine fit_ewma_diurnal_intraday_ewma_baseline

    ! Aggregate squared intraday returns by consecutive date identifiers.
    subroutine daily_realized_variance(y, date_id, day_rv, obs_day, ndays)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: date_id(:)
        real(dp), allocatable, intent(out) :: day_rv(:)
        integer, allocatable, intent(out) :: obs_day(:)
        integer, intent(out) :: ndays
        integer :: i, current_date

        if (size(y) /= size(date_id)) error stop "daily_realized_variance: array sizes differ"
        if (size(y) < 1) error stop "daily_realized_variance: empty input"
        allocate(day_rv(size(y)), obs_day(size(y)))
        day_rv = 0.0_dp
        ndays = 0
        current_date = huge(1)
        do i = 1, size(y)
            if (date_id(i) /= current_date) then
                ndays = ndays + 1
                current_date = date_id(i)
            end if
            obs_day(i) = ndays
            day_rv(ndays) = day_rv(ndays) + y(i)**2
        end do
        day_rv = day_rv(1:ndays)
    end subroutine daily_realized_variance

    ! Smooth a positive bin curve using a triangular weighted average in log variance.
    subroutine smooth_log_curve(bin_curve, bin_count, half_width)
        real(dp), intent(inout) :: bin_curve(:)
        integer, intent(in) :: bin_count(:), half_width
        real(dp), allocatable :: smoothed(:)
        integer :: b, j, lo, hi
        real(dp) :: w, wsum, xsum

        if (half_width < 1) return
        allocate(smoothed(size(bin_curve)))
        do b = 1, size(bin_curve)
            lo = max(1, b - half_width)
            hi = min(size(bin_curve), b + half_width)
            wsum = 0.0_dp
            xsum = 0.0_dp
            do j = lo, hi
                if (bin_count(j) < 1) cycle
                w = real(half_width + 1 - abs(j - b), dp)
                wsum = wsum + w
                xsum = xsum + w*log(max(bin_curve(j), min_var))
            end do
            if (wsum > 0.0_dp) then
                smoothed(b) = exp(xsum / wsum)
            else
                smoothed(b) = bin_curve(b)
            end if
        end do
        bin_curve = max(smoothed, min_var)
        deallocate(smoothed)
    end subroutine smooth_log_curve

    ! Validate the common one-dimensional input/output array sizes.
    subroutine check_same_size(n1, n2, n3, n4, caller)
        integer, intent(in) :: n1, n2, n3, n4
        character(len=*), intent(in) :: caller

        if (n1 /= n2 .or. n1 /= n3 .or. n1 /= n4) then
            error stop caller // ": array sizes differ"
        end if
    end subroutine check_same_size

    ! Expected remaining-session variance at each bar using only the diurnal multiplier.
    !
    ! At bar i with bin b, the forecast is daily_var(i) * sum_{j=b+1}^{nbins} diurnal_curve(j).
    ! This is the GARCH-free baseline: no intraday dynamics, purely time-of-day weighting.
    subroutine remaining_session_var_diurnal(daily_var, diurnal_var, bin_id, h_remaining)
        real(dp), intent(in) :: daily_var(:), diurnal_var(:)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(out) :: h_remaining(:)
        real(dp), allocatable :: bin_curve(:), bin_tail_sum(:)
        integer :: i, b, nbins

        if (size(daily_var) /= size(diurnal_var) .or. size(daily_var) /= size(bin_id) .or. &
            size(daily_var) /= size(h_remaining)) then
            error stop "remaining_session_var_diurnal: array sizes differ"
        end if
        if (minval(bin_id) < 1) error stop "remaining_session_var_diurnal: bin_id must be positive"
        nbins = maxval(bin_id)
        allocate(bin_curve(nbins), bin_tail_sum(nbins))
        bin_curve = 0.0_dp
        do i = 1, size(daily_var)
            bin_curve(bin_id(i)) = diurnal_var(i)
        end do
        bin_tail_sum(nbins) = 0.0_dp
        do b = nbins - 1, 1, -1
            bin_tail_sum(b) = bin_tail_sum(b + 1) + bin_curve(b + 1)
        end do
        do i = 1, size(daily_var)
            h_remaining(i) = max(daily_var(i) * bin_tail_sum(bin_id(i)), 0.0_dp)
        end do
        deallocate(bin_curve, bin_tail_sum)
    end subroutine remaining_session_var_diurnal

end module intraday_vol_baseline_mod

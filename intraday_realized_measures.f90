! Daily realized measures computed from intraday OHLC bars.

module intraday_realized_measures_mod
    use kind_mod, only: dp
    use date_mod, only: yyyymmdd
    use market_data_mod, only: ohlcv_series_t
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    type, public :: daily_realized_panel_t
        integer, allocatable :: date(:)
        real(dp), allocatable :: open(:)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: close(:)
        real(dp), allocatable :: rv(:)
        real(dp), allocatable :: bpv(:)
        real(dp), allocatable :: rsv_neg(:)
        real(dp), allocatable :: rsv_pos(:)
        real(dp), allocatable :: park(:)
        real(dp), allocatable :: gk(:)
        real(dp), allocatable :: rq(:)
    end type daily_realized_panel_t

    public :: build_daily_realized_panel, select_realized_measure, daily_log_rv
    public :: daily_log_rv_and_rsv_neg, daily_rv_rsv_and_returns

contains

    subroutine build_daily_realized_panel(bars, daily)
        ! Aggregate intraday OHLC bars to daily OHLC and realized variance measures.
        type(ohlcv_series_t), intent(in) :: bars
        type(daily_realized_panel_t), intent(out) :: daily
        integer, allocatable :: day_index(:), first(:), last(:)
        integer :: n, i, d, ndays, current_date, nday_bars
        real(dp) :: r, rprev, hl, co

        n = bars%nobs()
        allocate(day_index(n), first(n), last(n))
        ndays = 0
        current_date = -1
        do i = 1, n
            if (yyyymmdd(bars%timestamp(i)%date) /= current_date) then
                ndays = ndays + 1
                current_date = yyyymmdd(bars%timestamp(i)%date)
                first(ndays) = i
                if (ndays > 1) last(ndays - 1) = i - 1
            end if
            day_index(i) = ndays
        end do
        last(ndays) = n
        first = first(1:ndays)
        last = last(1:ndays)

        allocate(daily%date(ndays), daily%open(ndays), daily%high(ndays), daily%low(ndays), daily%close(ndays))
        allocate(daily%rv(ndays), daily%bpv(ndays), daily%rsv_neg(ndays), daily%rsv_pos(ndays), daily%park(ndays), &
                 daily%gk(ndays), daily%rq(ndays))
        daily%rv = 0.0_dp
        daily%bpv = 0.0_dp
        daily%rsv_neg = 0.0_dp
        daily%rsv_pos = 0.0_dp
        daily%park = 0.0_dp
        daily%gk = 0.0_dp
        daily%rq = 0.0_dp
        do d = 1, ndays
            nday_bars = last(d) - first(d) + 1
            daily%date(d) = yyyymmdd(bars%timestamp(first(d))%date)
            daily%open(d) = bars%open(first(d))
            daily%high(d) = maxval(bars%high(first(d):last(d)))
            daily%low(d) = minval(bars%low(first(d):last(d)))
            daily%close(d) = bars%close(last(d))
            rprev = 0.0_dp
            do i = first(d), last(d)
                r = log(bars%close(i) / bars%open(i))
                daily%rv(d) = daily%rv(d) + r**2
                daily%rq(d) = daily%rq(d) + r**4
                if (i > first(d)) daily%bpv(d) = daily%bpv(d) + 0.5_dp*acos(-1.0_dp)*abs(r)*abs(rprev)
                rprev = r
                if (r < 0.0_dp) then
                    daily%rsv_neg(d) = daily%rsv_neg(d) + r**2
                else
                    daily%rsv_pos(d) = daily%rsv_pos(d) + r**2
                end if
                hl = log(bars%high(i) / bars%low(i))
                co = log(bars%close(i) / bars%open(i))
                daily%park(d) = daily%park(d) + hl**2 / (4.0_dp*log(2.0_dp))
                daily%gk(d) = daily%gk(d) + max(0.5_dp*hl**2 - (2.0_dp*log(2.0_dp) - 1.0_dp)*co**2, min_var)
            end do
            daily%rq(d) = max(real(nday_bars, dp)*daily%rq(d) / 3.0_dp, min_var)
        end do
        deallocate(day_index, first, last)
    end subroutine build_daily_realized_panel

    subroutine select_realized_measure(daily, name, x)
        ! Copy a named realized measure from a daily panel.
        type(daily_realized_panel_t), intent(in) :: daily
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: x(:)

        select case (trim(name))
        case ("RV")
            x = daily%rv(1:size(x))
        case ("BPV")
            x = daily%bpv(1:size(x))
        case ("RSV_NEG")
            x = daily%rsv_neg(1:size(x))
        case ("RSV_POS")
            x = daily%rsv_pos(1:size(x))
        case ("PARKINSON")
            x = daily%park(1:size(x))
        case ("GARMAN_KLASS")
            x = daily%gk(1:size(x))
        case ("RQ")
            x = daily%rq(1:size(x))
        case default
            error stop "select_realized_measure: unknown measure"
        end select
        x = max(x, min_var)
    end subroutine select_realized_measure

    ! Compute daily log realized variance from close-to-close intraday bar returns.
    ! rv_day = sum of squared log(close_i / close_{i-1}) for consecutive same-day bars.
    subroutine daily_log_rv(bars, log_rv, dates, ndays)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: log_rv(:)
        integer,  allocatable, intent(out) :: dates(:)
        integer,  intent(out) :: ndays

        integer :: i, n, cur_date, prev_date
        real(dp) :: r, rv_today
        integer,  allocatable :: tmp_dates(:)
        real(dp), allocatable :: tmp_rv(:)

        n = bars%nobs()
        allocate(tmp_dates(n), tmp_rv(n))
        ndays    = 0
        prev_date = -1
        rv_today  = 0.0_dp

        do i = 1, n
            cur_date = yyyymmdd(bars%timestamp(i)%date)
            if (cur_date /= prev_date) then
                if (prev_date > 0 .and. rv_today > 0.0_dp) then
                    ndays = ndays + 1
                    tmp_dates(ndays) = prev_date
                    tmp_rv(ndays)    = log(rv_today)
                end if
                rv_today  = 0.0_dp
                prev_date = cur_date
            else
                r        = log(bars%close(i) / bars%close(i-1))
                rv_today = rv_today + r*r
            end if
        end do
        if (prev_date > 0 .and. rv_today > 0.0_dp) then
            ndays = ndays + 1
            tmp_dates(ndays) = prev_date
            tmp_rv(ndays)    = log(rv_today)
        end if

        allocate(log_rv(ndays), dates(ndays))
        log_rv = tmp_rv(1:ndays)
        dates  = tmp_dates(1:ndays)
        deallocate(tmp_dates, tmp_rv)
    end subroutine daily_log_rv

    subroutine daily_log_rv_and_rsv_neg(bars, log_rv, log_rsv_neg, dates, ndays)
        ! Compute log(RV) and log(RSV^-) from intraday bars in a single pass.
        ! RSV^- = sum of squared intraday CC returns on down-ticks (negative semi-variance).
        type(ohlcv_series_t), intent(in)   :: bars
        real(dp), allocatable, intent(out) :: log_rv(:)       ! daily log realized variance
        real(dp), allocatable, intent(out) :: log_rsv_neg(:)  ! daily log negative semi-variance
        integer,  allocatable, intent(out) :: dates(:)
        integer,              intent(out)  :: ndays

        integer :: i, n, cur_date, prev_date
        real(dp) :: r, rv_today, rsv_neg_today
        integer,  allocatable :: tmp_dates(:)
        real(dp), allocatable :: tmp_rv(:), tmp_neg(:)
        real(dp), parameter :: min_rsv = 1.0e-12_dp

        n = bars%nobs()
        allocate(tmp_dates(n), tmp_rv(n), tmp_neg(n))
        ndays         = 0
        prev_date     = -1
        rv_today      = 0.0_dp
        rsv_neg_today = 0.0_dp

        do i = 1, n
            cur_date = yyyymmdd(bars%timestamp(i)%date)
            if (cur_date /= prev_date) then
                if (prev_date > 0 .and. rv_today > 0.0_dp) then
                    ndays             = ndays + 1
                    tmp_dates(ndays)  = prev_date
                    tmp_rv(ndays)     = log(rv_today)
                    tmp_neg(ndays)    = log(max(rsv_neg_today, min_rsv))
                end if
                rv_today      = 0.0_dp
                rsv_neg_today = 0.0_dp
                prev_date     = cur_date
            else
                r        = log(bars%close(i) / bars%close(i-1))
                rv_today = rv_today + r*r
                if (r < 0.0_dp) rsv_neg_today = rsv_neg_today + r*r
            end if
        end do
        if (prev_date > 0 .and. rv_today > 0.0_dp) then
            ndays            = ndays + 1
            tmp_dates(ndays) = prev_date
            tmp_rv(ndays)    = log(rv_today)
            tmp_neg(ndays)   = log(max(rsv_neg_today, min_rsv))
        end if

        allocate(log_rv(ndays), log_rsv_neg(ndays), dates(ndays))
        log_rv      = tmp_rv(1:ndays)
        log_rsv_neg = tmp_neg(1:ndays)
        dates       = tmp_dates(1:ndays)
        deallocate(tmp_dates, tmp_rv, tmp_neg)
    end subroutine daily_log_rv_and_rsv_neg

    subroutine daily_rv_rsv_and_returns(bars, log_rv, log_rsv_neg, ret_cc, dates, ndays)
        ! Compute log(RV), log(RSV^-), and daily close-to-close log-returns in one pass.
        type(ohlcv_series_t), intent(in)   :: bars
        real(dp), allocatable, intent(out) :: log_rv(:)       ! daily log realized variance
        real(dp), allocatable, intent(out) :: log_rsv_neg(:)  ! daily log negative semi-variance
        real(dp), allocatable, intent(out) :: ret_cc(:)       ! daily close-to-close log-return
        integer,  allocatable, intent(out) :: dates(:)
        integer,              intent(out)  :: ndays

        integer :: i, n, cur_date, prev_date
        real(dp) :: r, rv_today, rsv_neg_today, last_close, prev_day_close
        integer,  allocatable :: tmp_dates(:)
        real(dp), allocatable :: tmp_rv(:), tmp_neg(:), tmp_ret(:)
        real(dp), parameter :: min_rsv = 1.0e-12_dp

        n = bars%nobs()
        allocate(tmp_dates(n), tmp_rv(n), tmp_neg(n), tmp_ret(n))
        ndays          = 0
        prev_date      = -1
        rv_today       = 0.0_dp
        rsv_neg_today  = 0.0_dp
        last_close     = 0.0_dp
        prev_day_close = 0.0_dp

        do i = 1, n
            cur_date = yyyymmdd(bars%timestamp(i)%date)
            if (cur_date /= prev_date) then
                if (prev_date > 0 .and. rv_today > 0.0_dp) then
                    ndays             = ndays + 1
                    tmp_dates(ndays)  = prev_date
                    tmp_rv(ndays)     = log(rv_today)
                    tmp_neg(ndays)    = log(max(rsv_neg_today, min_rsv))
                    ! last_close holds bars%close(i-1) = close of last bar of prev day
                    if (prev_day_close > 0.0_dp) then
                        tmp_ret(ndays) = log(last_close / prev_day_close)
                    else
                        tmp_ret(ndays) = 0.0_dp
                    end if
                    prev_day_close = last_close
                end if
                rv_today      = 0.0_dp
                rsv_neg_today = 0.0_dp
                prev_date     = cur_date
            else
                r        = log(bars%close(i) / bars%close(i-1))
                rv_today = rv_today + r*r
                if (r < 0.0_dp) rsv_neg_today = rsv_neg_today + r*r
            end if
            last_close = bars%close(i)
        end do
        if (prev_date > 0 .and. rv_today > 0.0_dp) then
            ndays             = ndays + 1
            tmp_dates(ndays)  = prev_date
            tmp_rv(ndays)     = log(rv_today)
            tmp_neg(ndays)    = log(max(rsv_neg_today, min_rsv))
            if (prev_day_close > 0.0_dp) then
                tmp_ret(ndays) = log(last_close / prev_day_close)
            else
                tmp_ret(ndays) = 0.0_dp
            end if
        end if

        allocate(log_rv(ndays), log_rsv_neg(ndays), ret_cc(ndays), dates(ndays))
        log_rv      = tmp_rv(1:ndays)
        log_rsv_neg = tmp_neg(1:ndays)
        ret_cc      = tmp_ret(1:ndays)
        dates       = tmp_dates(1:ndays)
        deallocate(tmp_dates, tmp_rv, tmp_neg, tmp_ret)
    end subroutine daily_rv_rsv_and_returns

end module intraday_realized_measures_mod

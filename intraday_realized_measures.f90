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

    public :: build_daily_realized_panel, select_realized_measure

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

end module intraday_realized_measures_mod

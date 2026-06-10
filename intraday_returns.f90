! Helpers for building and aligning intraday return series.

module intraday_returns_mod
    use iso_fortran_env, only: int64
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, resample_ohlcv_series, default_session_start_seconds, &
                               default_session_end_seconds
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: trading_days_per_year_default = 252.0_dp
    integer, parameter :: seconds_per_day = 86400

    type, public :: return_series_t
        integer(int64), allocatable :: key(:)
        real(dp), allocatable :: ret(:)
    end type return_series_t

    public :: returns_at_frequency
    public :: close_open_returns
    public :: aligned_return_matrix
    public :: timestamp_key
    public :: realized_vol_ann_pct
    public :: frequency_label
    public :: seconds_per_day

contains

    ! Build log(close/open) returns at the requested aggregation scale.
    subroutine returns_at_frequency(series, source_seconds, target_seconds, returns, session_start_seconds)
        type(ohlcv_series_t), intent(in) :: series
        integer, intent(in) :: source_seconds, target_seconds
        type(return_series_t), intent(out) :: returns
        integer, intent(in), optional :: session_start_seconds
        type(ohlcv_series_t) :: sampled
        integer :: start_sec

        start_sec = default_session_start_seconds
        if (present(session_start_seconds)) start_sec = session_start_seconds
        if (target_seconds == source_seconds) then
            call close_open_returns(series, returns)
        else
            call resample_ohlcv_series(series, sampled, source_seconds, target_seconds, start_sec)
            call close_open_returns(sampled, returns)
        end if
    end subroutine returns_at_frequency

    ! Compute log(close/open) returns and timestamp keys for OHLC bars.
    subroutine close_open_returns(series, returns)
        type(ohlcv_series_t), intent(in) :: series
        type(return_series_t), intent(out) :: returns
        integer :: i, n

        n = series%nobs()
        allocate(returns%key(n), returns%ret(n))
        do i = 1, n
            returns%key(i) = timestamp_key(series, i)
            returns%ret(i) = log(max(series%close(i), min_var) / max(series%open(i), min_var))
        end do
    end subroutine close_open_returns

    ! Align return series by timestamp key using the first asset as the reference.
    subroutine aligned_return_matrix(returns, x)
        type(return_series_t), intent(in) :: returns(:)
        real(dp), allocatable, intent(out) :: x(:, :)
        real(dp), allocatable :: tmp(:, :)
        integer, allocatable :: pos(:)
        integer :: nasset, iref, iasset, k, nmatch

        nasset = size(returns)
        allocate(tmp(size(returns(1)%key), nasset), pos(nasset))
        pos = 1
        nmatch = 0
        do iref = 1, size(returns(1)%key)
            tmp(nmatch + 1, 1) = returns(1)%ret(iref)
            do iasset = 2, nasset
                k = advance_to_key(returns(iasset)%key, pos(iasset), returns(1)%key(iref))
                if (k == 0) exit
                tmp(nmatch + 1, iasset) = returns(iasset)%ret(k)
            end do
            if (iasset > nasset) nmatch = nmatch + 1
        end do
        allocate(x(nmatch, nasset))
        if (nmatch > 0) x = tmp(1:nmatch, :)
        deallocate(tmp, pos)
    end subroutine aligned_return_matrix

    ! Timestamp key for alignment.
    integer(int64) function timestamp_key(series, i) result(key)
        type(ohlcv_series_t), intent(in) :: series
        integer, intent(in) :: i

        key = int(series%timestamp(i)%date%to_yyyymmdd(), int64)*100000_int64 + &
              int(series%timestamp(i)%seconds_since_midnight(), int64)
    end function timestamp_key

    ! Advance a sorted key pointer to target and return its index, or zero.
    integer function advance_to_key(keys, pos, target) result(idx)
        integer(int64), intent(in) :: keys(:), target
        integer, intent(inout) :: pos

        do while (pos <= size(keys))
            if (keys(pos) >= target) exit
            pos = pos + 1
        end do
        idx = 0
        if (pos <= size(keys)) then
            if (keys(pos) == target) idx = pos
        end if
    end function advance_to_key

    ! Annualized realized volatility percent from returns at a fixed frequency.
    pure real(dp) function realized_vol_ann_pct(ret, target_seconds, regular_session, &
                                                trading_days_per_year) result(vol)
        real(dp), intent(in) :: ret(:)
        integer, intent(in) :: target_seconds
        logical, intent(in), optional :: regular_session
        real(dp), intent(in), optional :: trading_days_per_year
        real(dp) :: periods_per_year, tdays
        logical :: is_regular

        if (target_seconds < 1) error stop "realized_vol_ann_pct: target_seconds must be positive"
        tdays = trading_days_per_year_default
        if (present(trading_days_per_year)) tdays = trading_days_per_year
        is_regular = .false.
        if (present(regular_session)) is_regular = regular_session
        if (target_seconds >= seconds_per_day) then
            periods_per_year = tdays
        else if (is_regular) then
            periods_per_year = tdays * real(default_session_end_seconds - default_session_start_seconds, dp) / &
                               real(target_seconds, dp)
        else
            periods_per_year = tdays * real(seconds_per_day, dp) / real(target_seconds, dp)
        end if
        vol = 100.0_dp * sqrt(max(sum(ret**2) / real(size(ret), dp), min_var) * periods_per_year)
    end function realized_vol_ann_pct

    ! Label a frequency in minutes, hours, or days.
    pure character(len=16) function frequency_label(seconds) result(label)
        integer, intent(in) :: seconds

        label = ""
        if (seconds == seconds_per_day) then
            label = "daily"
        else if (mod(seconds, 3600) == 0) then
            write(label, '(I0,A)') seconds / 3600, "h"
        else
            write(label, '(I0,A)') seconds / 60, "m"
        end if
    end function frequency_label

end module intraday_returns_mod

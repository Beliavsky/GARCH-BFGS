! Compute autocorrelations of intraday high-low price ranges
! for several sampling frequencies.

program xacf_price_ranges
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, &
                               resample_ohlcv_series
    implicit none

    integer, parameter :: max_lag = 10
    integer, parameter :: source_seconds = 1
    integer, parameter :: freq_seconds(*) = [1, 10, 30, 60, 150, 300]
    character(len=256) :: filename
    type(ohlcv_series_t) :: raw_series, regular_series, sampled_series
    real(dp), allocatable :: price_range(:)
    real(dp) :: acf(max_lag)
    integer :: ifreq, nrange

    filename = "c:\python\databento\spy_1s_databento.csv"
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    call read_intraday_prices_csv(trim(filename), raw_series)
    call filter_intraday_session(raw_series, regular_series)

    print '(A)', "Autocorrelations of high-low price ranges"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,I0)', "Regular-session 1s rows: ", regular_series%nobs()
    print '(A)', "Ranges are computed after regular-session filtering and frequency resampling."
    print '(A)', ""
    print '(A)', "------------------------------------------------------------------------------------------------------------"
    print '(A8,1X,A10,10(1X,A8))', "freq_sec", "nrange", "acf_1", "acf_2", "acf_3", "acf_4", "acf_5", &
          "acf_6", "acf_7", "acf_8", "acf_9", "acf_10"
    print '(A)', "------------------------------------------------------------------------------------------------------------"

    do ifreq = 1, size(freq_seconds)
        if (freq_seconds(ifreq) == source_seconds) then
            call high_low_ranges(regular_series, price_range, nrange)
        else
            call resample_ohlcv_series(regular_series, sampled_series, source_seconds, freq_seconds(ifreq))
            call high_low_ranges(sampled_series, price_range, nrange)
        end if
        call autocorr(price_range(1:nrange), acf)
        print '(I8,1X,I10,10(1X,F8.4))', freq_seconds(ifreq), nrange, acf
        if (allocated(price_range)) deallocate(price_range)
    end do
    print '(A)', "------------------------------------------------------------------------------------------------------------"

contains

    ! Compute high-low price ranges for each bar.
    subroutine high_low_ranges(series, price_range, nrange)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), allocatable, intent(out) :: price_range(:)
        integer, intent(out) :: nrange

        nrange = series%nobs()
        if (nrange < max_lag + 2) error stop "high_low_ranges: not enough observations"
        allocate(price_range(nrange))
        price_range = max(series%high - series%low, 0.0_dp)
    end subroutine high_low_ranges

    ! Compute autocorrelations at lags 1:size(acf).
    subroutine autocorr(x, acf)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: acf(:)
        real(dp) :: mu, denom, num
        integer :: lag, i, n

        n = size(x)
        mu = sum(x) / real(n, dp)
        denom = sum((x - mu)**2)
        if (denom <= 0.0_dp) then
            acf = 0.0_dp
            return
        end if
        do lag = 1, size(acf)
            num = 0.0_dp
            do i = lag + 1, n
                num = num + (x(i) - mu)*(x(i - lag) - mu)
            end do
            acf(lag) = num / denom
        end do
    end subroutine autocorr
end program xacf_price_ranges

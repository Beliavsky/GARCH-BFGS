! Compute autocorrelations of signed close-to-close price changes
! for several intraday sampling frequencies.

program xacf_price_changes
    use kind_mod, only: dp
    use date_mod, only: yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, &
                               resample_ohlcv_series
    implicit none

    integer, parameter :: max_lag = 10
    integer, parameter :: source_seconds = 1
    integer, parameter :: freq_seconds(*) = [1, 10, 30, 60, 150, 300]
    character(len=256) :: filename
    type(ohlcv_series_t) :: raw_series, regular_series, sampled_series
    real(dp), allocatable :: change(:)
    real(dp) :: acf(max_lag)
    integer :: ifreq, nchange

    filename = "c:\python\databento\spy_1s_databento.csv"
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    call read_intraday_prices_csv(trim(filename), raw_series)
    call filter_intraday_session(raw_series, regular_series)

    print '(A)', "Autocorrelations of signed close-to-close price changes"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,I0)', "Regular-session 1s rows: ", regular_series%nobs()
    print '(A)', "Overnight and cross-date changes are excluded."
    print '(A)', ""
    print '(A)', "------------------------------------------------------------------------------------------------------------"
    print '(A8,1X,A10,10(1X,A8))', "freq_sec", "nchange", "acf_1", "acf_2", "acf_3", "acf_4", "acf_5", &
          "acf_6", "acf_7", "acf_8", "acf_9", "acf_10"
    print '(A)', "------------------------------------------------------------------------------------------------------------"

    do ifreq = 1, size(freq_seconds)
        if (freq_seconds(ifreq) == source_seconds) then
            call close_changes(regular_series, change, nchange)
        else
            call resample_ohlcv_series(regular_series, sampled_series, source_seconds, freq_seconds(ifreq))
            call close_changes(sampled_series, change, nchange)
        end if
        call autocorr(change(1:nchange), acf)
        print '(I8,1X,I10,10(1X,F8.4))', freq_seconds(ifreq), nchange, acf
        if (allocated(change)) deallocate(change)
    end do
    print '(A)', "------------------------------------------------------------------------------------------------------------"

contains

    ! Compute within-day signed close-to-close price changes.
    subroutine close_changes(series, change, nchange)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), allocatable, intent(out) :: change(:)
        integer, intent(out) :: nchange
        integer :: i

        allocate(change(max(series%nobs() - 1, 0)))
        nchange = 0
        do i = 2, series%nobs()
            if (yyyymmdd(series%timestamp(i)%date) /= yyyymmdd(series%timestamp(i - 1)%date)) cycle
            nchange = nchange + 1
            change(nchange) = series%close(i) - series%close(i - 1)
        end do
        if (nchange < max_lag + 2) error stop "close_changes: not enough observations"
    end subroutine close_changes

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
end program xacf_price_changes

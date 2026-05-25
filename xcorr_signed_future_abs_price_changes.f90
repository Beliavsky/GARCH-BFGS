! Correlate signed close-to-close price changes with future absolute changes
! for several intraday sampling frequencies.

program xcorr_signed_future_abs_price_changes
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
    real(dp), allocatable :: signed_change(:)
    integer, allocatable :: change_date(:)
    real(dp) :: corr(max_lag)
    integer :: ifreq, nchange

    filename = "c:\python\databento\spy_1s_databento.csv"
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    call read_intraday_prices_csv(trim(filename), raw_series)
    call filter_intraday_session(raw_series, regular_series)

    print '(A)', "Correlation of signed price changes with future absolute price changes"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,I0)', "Regular-session 1s rows: ", regular_series%nobs()
    print '(A)', "Pairs require both changes to be within the same trading date."
    print '(A)', ""
    print '(A)', "------------------------------------------------------------------------------------------------------------"
    print '(A8,1X,A10,10(1X,A8))', "freq_sec", "nchange", "lead_1", "lead_2", "lead_3", "lead_4", "lead_5", &
          "lead_6", "lead_7", "lead_8", "lead_9", "lead_10"
    print '(A)', "------------------------------------------------------------------------------------------------------------"

    do ifreq = 1, size(freq_seconds)
        if (freq_seconds(ifreq) == source_seconds) then
            call close_signed_changes(regular_series, signed_change, change_date, nchange)
        else
            call resample_ohlcv_series(regular_series, sampled_series, source_seconds, freq_seconds(ifreq))
            call close_signed_changes(sampled_series, signed_change, change_date, nchange)
        end if
        call signed_future_abs_corr(signed_change(1:nchange), change_date(1:nchange), corr)
        print '(I8,1X,I10,10(1X,F8.4))', freq_seconds(ifreq), nchange, corr
        if (allocated(signed_change)) deallocate(signed_change, change_date)
    end do
    print '(A)', "------------------------------------------------------------------------------------------------------------"

contains

    ! Compute within-day signed close-to-close price changes and their dates.
    subroutine close_signed_changes(series, signed_change, change_date, nchange)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), allocatable, intent(out) :: signed_change(:)
        integer, allocatable, intent(out) :: change_date(:)
        integer, intent(out) :: nchange
        integer :: i, date_i

        allocate(signed_change(max(series%nobs() - 1, 0)), change_date(max(series%nobs() - 1, 0)))
        nchange = 0
        do i = 2, series%nobs()
            date_i = yyyymmdd(series%timestamp(i)%date)
            if (date_i /= yyyymmdd(series%timestamp(i - 1)%date)) cycle
            nchange = nchange + 1
            signed_change(nchange) = series%close(i) - series%close(i - 1)
            change_date(nchange) = date_i
        end do
        if (nchange < max_lag + 2) error stop "close_signed_changes: not enough observations"
    end subroutine close_signed_changes

    ! Correlate x_t with abs(x_{t+lag}), excluding cross-date pairs.
    subroutine signed_future_abs_corr(x, change_date, corr)
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: change_date(:)
        real(dp), intent(out) :: corr(:)
        real(dp), allocatable :: lhs(:), rhs(:)
        integer :: lag, i, n, k

        n = size(x)
        if (size(change_date) /= n) error stop "signed_future_abs_corr: array sizes differ"
        allocate(lhs(n), rhs(n))
        do lag = 1, size(corr)
            k = 0
            do i = 1, n - lag
                if (change_date(i) /= change_date(i + lag)) cycle
                k = k + 1
                lhs(k) = x(i)
                rhs(k) = abs(x(i + lag))
            end do
            if (k < 3) then
                corr(lag) = 0.0_dp
            else
                corr(lag) = correlation(lhs(1:k), rhs(1:k))
            end if
        end do
        deallocate(lhs, rhs)
    end subroutine signed_future_abs_corr

    ! Pearson correlation.
    real(dp) function correlation(x, y) result(rho)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: xmean, ymean, xden, yden, num

        if (size(x) /= size(y)) error stop "correlation: array sizes differ"
        xmean = sum(x) / real(size(x), dp)
        ymean = sum(y) / real(size(y), dp)
        num = sum((x - xmean)*(y - ymean))
        xden = sum((x - xmean)**2)
        yden = sum((y - ymean)**2)
        if (xden <= 0.0_dp .or. yden <= 0.0_dp) then
            rho = 0.0_dp
        else
            rho = num / sqrt(xden*yden)
        end if
    end function correlation
end program xcorr_signed_future_abs_price_changes

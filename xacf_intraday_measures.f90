! Compute ACFs for signed price changes, absolute price changes,
! and high-low ranges at several intraday sampling frequencies.

program xacf_intraday_measures
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, ohlcv_tick_series_t, read_intraday_prices_csv, &
                               read_ohlcv_tick_stream, convert_tick_to_ohlcv_series, &
                               filter_intraday_session, resample_ohlcv_series, &
                               intraday_close_changes, intraday_high_low_ranges
    use path_utils_mod, only: has_extension
    use stats_mod, only: autocorr, print_acf_table
    implicit none

    integer, parameter :: max_lag = 10
    integer, parameter :: source_seconds = 1
    integer, parameter :: freq_seconds(*) = [1, 10, 30, 60, 150, 300]
    character(len=*), parameter :: default_csv_file = "c:\python\databento\spy_1s_databento.csv"
    character(len=*), parameter :: default_bin_file = "spy_1s_databento.bin"
    character(len=256) :: filename
    type(ohlcv_series_t) :: raw_series, regular_series, sampled_series
    type(ohlcv_tick_series_t) :: tick_series
    real(dp), allocatable :: signed_change(:), abs_change(:), price_range(:)
    real(dp), allocatable :: acf_signed(:, :), acf_abs(:, :), acf_range(:, :)
    integer :: ifreq, nchange, nrange
    integer :: nchange_by_freq(size(freq_seconds)), nrange_by_freq(size(freq_seconds))
    real(dp) :: t0, t1, read_sec, series_sec, acf_sec, elapsed_sec
    real(dp) :: t_section0, t_section1
    logical :: default_bin_exists

    inquire(file=default_bin_file, exist=default_bin_exists)
    if (default_bin_exists) then
        filename = default_bin_file
    else
        filename = default_csv_file
    end if
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    allocate(acf_signed(size(freq_seconds), max_lag), acf_abs(size(freq_seconds), max_lag), &
             acf_range(size(freq_seconds), max_lag))
    nchange_by_freq = 0
    nrange_by_freq = 0
    read_sec = 0.0_dp
    series_sec = 0.0_dp
    acf_sec = 0.0_dp

    call cpu_time(t0)
    call cpu_time(t_section0)
    if (has_extension(filename, ".bin")) then
        call read_ohlcv_tick_stream(filename, tick_series)
        call convert_tick_to_ohlcv_series(tick_series, raw_series)
    else
        call read_intraday_prices_csv(filename, raw_series)
    end if
    call filter_intraday_session(raw_series, regular_series)
    call cpu_time(t_section1)
    read_sec = t_section1 - t_section0

    do ifreq = 1, size(freq_seconds)
        call cpu_time(t_section0)
        if (freq_seconds(ifreq) == source_seconds) then
            call intraday_close_changes(regular_series, signed_change, abs_change, nchange)
            call intraday_high_low_ranges(regular_series, price_range, nrange)
        else
            call resample_ohlcv_series(regular_series, sampled_series, source_seconds, freq_seconds(ifreq))
            call intraday_close_changes(sampled_series, signed_change, abs_change, nchange)
            call intraday_high_low_ranges(sampled_series, price_range, nrange)
        end if
        if (nchange < max_lag + 2) error stop "xacf_intraday_measures: not enough close changes"
        if (nrange < max_lag + 2) error stop "xacf_intraday_measures: not enough ranges"
        call cpu_time(t_section1)
        series_sec = series_sec + (t_section1 - t_section0)

        call cpu_time(t_section0)
        call autocorr(signed_change(1:nchange), acf_signed(ifreq, :))
        call autocorr(abs_change(1:nchange), acf_abs(ifreq, :))
        call autocorr(price_range(1:nrange), acf_range(ifreq, :))
        call cpu_time(t_section1)
        acf_sec = acf_sec + (t_section1 - t_section0)

        nchange_by_freq(ifreq) = nchange
        nrange_by_freq(ifreq) = nrange
        if (allocated(signed_change)) deallocate(signed_change, abs_change, price_range)
    end do

    call cpu_time(t1)
    elapsed_sec = t1 - t0

    print '(A)', "Intraday ACF diagnostics"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,I0)', "Regular-session 1s rows: ", regular_series%nobs()
    print '(A,/)', "Overnight and cross-date close-to-close changes are excluded."

    call print_acf_table("Signed close-to-close price changes", "freq_sec", "nchange", &
                         freq_seconds, nchange_by_freq, acf_signed)
    call print_acf_table("Absolute close-to-close price changes", "freq_sec", "nchange", &
                         freq_seconds, nchange_by_freq, acf_abs)
    call print_acf_table("High-low price ranges", "freq_sec", "nrange", &
                         freq_seconds, nrange_by_freq, acf_range)

    print "(a)", "Times (sec)"
    print "(a20,f10.3)", "read_filter", read_sec, "build_series", series_sec, &
        "compute_acf", acf_sec, "elapsed", elapsed_sec

end program xacf_intraday_measures

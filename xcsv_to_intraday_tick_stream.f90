! Convert an intraday OHLCV CSV file to compact unformatted stream tick data.
program xcsv_to_intraday_tick_stream
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, ohlcv_tick_series_t, read_intraday_prices_csv, &
        convert_ohlcv_to_tick_series, write_ohlcv_tick_stream, default_price_tick_size
    use path_utils_mod, only: basename_with_extension
    implicit none

    character(len=512) :: input_file, output_file, arg
    type(ohlcv_series_t) :: text_series
    type(ohlcv_tick_series_t) :: tick_series
    real(dp) :: tick_size, t_start, t0, t1, read_text_sec, convert_sec, write_sec, elapsed_sec
    integer :: max_obs

    input_file = "c:\python\databento\spy_1s_databento.csv"
    tick_size = default_price_tick_size
    max_obs = 0
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) tick_size
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *) max_obs
    end if
    if (tick_size <= 0.0_dp) error stop "xcsv_to_intraday_tick_stream: tick_size must be positive"
    if (max_obs < 0) error stop "xcsv_to_intraday_tick_stream: max_obs must be nonnegative"
    output_file = basename_with_extension(input_file, ".bin")

    call cpu_time(t_start)

    call cpu_time(t0)
    if (max_obs > 0) then
        call read_intraday_prices_csv(input_file, text_series, max_obs=max_obs)
    else
        call read_intraday_prices_csv(input_file, text_series)
    end if
    call cpu_time(t1)
    read_text_sec = t1 - t0

    call cpu_time(t0)
    call convert_ohlcv_to_tick_series(text_series, tick_series, tick_size)
    call cpu_time(t1)
    convert_sec = t1 - t0

    call cpu_time(t0)
    call write_ohlcv_tick_stream(output_file, tick_series)
    call cpu_time(t1)
    write_sec = t1 - t0

    call cpu_time(t1)
    elapsed_sec = t1 - t_start

    print '(A)', "Intraday OHLCV CSV to tick stream"
    print '(A,A)', "Input CSV:       ", trim(input_file)
    print '(A,A)', "Output BIN:      ", trim(output_file)
    print '(A,I0)', "Observations:    ", tick_series%nobs()
    print '(A,F10.6)', "Tick size:       ", tick_series%tick_size
    if (max_obs > 0) print '(A,I0)', "Max observations:", max_obs
    print *
    print '(A)', "Times (sec)"
    print '(A24,1X,F10.3)', "read_text", read_text_sec
    print '(A24,1X,F10.3)', "convert", convert_sec
    print '(A24,1X,F10.3)', "write_stream", write_sec
    print '(A24,1X,F10.3)', "elapsed", elapsed_sec
end program xcsv_to_intraday_tick_stream

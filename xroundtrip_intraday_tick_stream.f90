! Round-trip test for compact intraday OHLCV tick data saved as an unformatted stream.

program xroundtrip_intraday_tick_stream
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, ohlcv_tick_series_t, read_intraday_prices_csv, &
        convert_ohlcv_to_tick_series, write_ohlcv_tick_stream, read_ohlcv_tick_stream, &
        equal_ohlcv_tick_series
    implicit none

    character(len=512) :: input_file, stream_file
    type(ohlcv_series_t) :: text_series
    type(ohlcv_tick_series_t) :: tick_series, tick_series_read
    real(dp), parameter :: tick_size = 0.005_dp
    real(dp) :: t_start, t0, t1, read_text_sec, convert_sec, write_sec, read_stream_sec
    real(dp) :: compare_sec, elapsed_sec
    logical :: same

    input_file = "c:\python\databento\spy_1s_databento.csv"
    stream_file = "intraday_tick_stream_test.bin"
    if (command_argument_count() >= 1) call get_command_argument(1, input_file)
    if (command_argument_count() >= 2) call get_command_argument(2, stream_file)

    call cpu_time(t_start)

    call cpu_time(t0)
    call read_intraday_prices_csv(trim(input_file), text_series)
    call cpu_time(t1)
    read_text_sec = t1 - t0

    call cpu_time(t0)
    call convert_ohlcv_to_tick_series(text_series, tick_series, tick_size)
    call cpu_time(t1)
    convert_sec = t1 - t0

    call cpu_time(t0)
    call write_ohlcv_tick_stream(trim(stream_file), tick_series)
    call cpu_time(t1)
    write_sec = t1 - t0

    call cpu_time(t0)
    call read_ohlcv_tick_stream(trim(stream_file), tick_series_read)
    call cpu_time(t1)
    read_stream_sec = t1 - t0

    call cpu_time(t0)
    same = equal_ohlcv_tick_series(tick_series, tick_series_read)
    call cpu_time(t1)
    compare_sec = t1 - t0

    call cpu_time(t1)
    elapsed_sec = t1 - t_start

    print '(A)', "Intraday OHLCV tick stream round-trip test"
    print '(A,A)', "Input text file:  ", trim(input_file)
    print '(A,A)', "Stream file:      ", trim(stream_file)
    print '(A,I0)', "Observations:     ", tick_series%nobs()
    print '(A,F10.6)', "Tick size:        ", tick_series%tick_size
    print '(A,L1)', "Round trip equal: ", same
    print *
    print '(A)', "Timing"
    print '(A)', "--------------------------------"
    print '(A24,1X,F10.3)', "read_text_sec", read_text_sec
    print '(A24,1X,F10.3)', "convert_sec", convert_sec
    print '(A24,1X,F10.3)', "write_stream_sec", write_sec
    print '(A24,1X,F10.3)', "read_stream_sec", read_stream_sec
    print '(A24,1X,F10.3)', "compare_sec", compare_sec
    print '(A24,1X,F10.3)', "elapsed_sec", elapsed_sec
    if (.not. same) error stop "xroundtrip_intraday_tick_stream: round-trip mismatch"
end program xroundtrip_intraday_tick_stream

! Convert an intraday OHLCV CSV file to compact unformatted stream tick data.
program xcsv_to_intraday_tick_stream
    use date_mod, only: print_program_header
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, ohlcv_tick_series_t, read_intraday_prices_csv_auto, &
        convert_ohlcv_to_tick_series, write_ohlcv_tick_stream, write_ohlcv_dp_stream, default_price_tick_size
    use path_utils_mod, only: basename_with_extension, basename_without_extension, dirname, &
        csv_files_in_dir, resolve_filename
    use glob_mod, only: glob, MAX_PATH_LEN
    use strings_mod, only: uppercase
    implicit none

    character(len=512), allocatable :: input_files(:)
    character(len=512) :: arg, data_dir, out_dir
    character(len=8) :: write_format
    logical :: same_dir
    real(dp) :: tick_size
    integer :: iarg, nfile, ifile
    integer :: max_obs

    call print_program_header("xcsv_to_intraday_tick_stream.f90")
    call parse_arguments(input_files, data_dir, out_dir)
    do ifile = 1, size(input_files)
        if (ifile > 1) print '(A)', ""
        call convert_one_file(input_files(ifile), out_dir)
    end do
    deallocate(input_files)

contains

    ! Convert one CSV file to basename.bin; directory is out_dir or same as input.
    subroutine convert_one_file(input_file, out_dir)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(in) :: out_dir
        character(len=512) :: output_file, eff_out_dir
        character(len=25) :: dt
        type(ohlcv_series_t) :: text_series
        type(ohlcv_tick_series_t) :: tick_series
        real(dp) :: t_start, t0, t1, read_text_sec, convert_sec, write_sec, elapsed_sec
        integer :: n_obs

        if (same_dir) then
            eff_out_dir = dirname(input_file)
        else
            eff_out_dir = out_dir
        end if
        output_file = resolve_filename(basename_with_extension(input_file, ".bin"), eff_out_dir)

        call cpu_time(t_start)

        call cpu_time(t0)
        if (max_obs > 0) then
            call read_intraday_prices_csv_auto(input_file, text_series, max_obs=max_obs)
        else
            call read_intraday_prices_csv_auto(input_file, text_series)
        end if
        call cpu_time(t1)
        read_text_sec = t1 - t0

        call cpu_time(t0)
        if (write_format == "dp") then
            call write_ohlcv_dp_stream(output_file, text_series)
            convert_sec = 0.0_dp
        else
            call convert_ohlcv_to_tick_series(text_series, tick_series, tick_size)
            call cpu_time(t1)
            convert_sec = t1 - t0
            call cpu_time(t0)
            call write_ohlcv_tick_stream(output_file, tick_series)
        end if
        call cpu_time(t1)
        write_sec = t1 - t0

        call cpu_time(t1)
        elapsed_sec = t1 - t_start

        if (write_format == "dp") then
            print '(A)', "Intraday OHLCV CSV to double-precision stream"
        else
            print '(A)', "Intraday OHLCV CSV to tick stream"
        end if
        print '(A,A)', "Input CSV:       ", trim(input_file)
        print '(A,A)', "Output BIN:      ", trim(output_file)
        n_obs = text_series%nobs()
        print '(A,I0)', "Observations:    ", n_obs
        if (write_format /= "dp") print '(A,F10.6)', "Tick size:       ", tick_series%tick_size
        if (max_obs > 0) print '(A,I0)', "Max observations:", max_obs
        print *
        print '(A)', "                   Datetime      Open      High       Low     Close       Volume"
        dt = text_series%timestamp(1)%to_str()
        print '("First ",A19,4(1X,F9.4),1X,I12)', dt(1:19), &
            text_series%open(1), text_series%high(1), text_series%low(1), &
            text_series%close(1), nint(text_series%volume(1))
        dt = text_series%timestamp(n_obs)%to_str()
        print '(" Last ",A19,4(1X,F9.4),1X,I12)', dt(1:19), &
            text_series%open(n_obs), text_series%high(n_obs), text_series%low(n_obs), &
            text_series%close(n_obs), nint(text_series%volume(n_obs))
        print *
        print '(A)', "Times (sec)"
        print '(A24,1X,F10.3)', "read_text", read_text_sec
        if (write_format /= "dp") print '(A24,1X,F10.3)', "convert", convert_sec
        print '(A24,1X,F10.3)', "write_stream", write_sec
        print '(A24,1X,F10.3)', "elapsed", elapsed_sec
    end subroutine convert_one_file

    subroutine parse_arguments(input_files, data_dir, out_dir)
        character(len=512), allocatable, intent(out) :: input_files(:)
        character(len=512), intent(out) :: data_dir, out_dir
        character(len=512), allocatable :: raw(:), expanded(:)
        character(len=MAX_PATH_LEN), allocatable :: glob_matches(:)
        integer :: nargs, i

        tick_size = default_price_tick_size
        max_obs = 0
        write_format = "tick"
        same_dir = .false.
        data_dir = ""
        out_dir = ""
        nargs = command_argument_count()
        allocate(raw(max(nargs, 1)))
        nfile = 0
        do iarg = 1, nargs
            call get_command_argument(iarg, arg)
            if (arg == "--same-dir") then
                same_dir = .true.
            else if (index(arg, "dir=") == 1) then
                data_dir = arg(5:)
            else if (index(arg, "out_dir=") == 1) then
                out_dir = arg(9:)
            else if (index(arg, "tick_size=") == 1) then
                read(arg(11:), *) tick_size
            else if (index(arg, "max_obs=") == 1) then
                read(arg(9:), *) max_obs
            else if (index(arg, "format=") == 1) then
                write_format = arg(8:)
            else
                nfile = nfile + 1
                raw(nfile) = arg
            end if
        end do
        if (nfile == 0) then
            if (len_trim(data_dir) > 0) then
                call csv_files_in_dir(data_dir, expanded)
                if (size(expanded) < 1) error stop "xcsv_to_intraday_tick_stream: no CSV files found in dir"
                call keep_price_csv_files(expanded, input_files)
                deallocate(expanded, raw)
                if (size(input_files) < 1) error stop "xcsv_to_intraday_tick_stream: no price CSV files found in dir"
                if (tick_size <= 0.0_dp) error stop "xcsv_to_intraday_tick_stream: tick_size must be positive"
                if (max_obs < 0) error stop "xcsv_to_intraday_tick_stream: max_obs must be nonnegative"
                return
            else
                nfile = 1
                raw(1) = "c:\python\databento\spy_1s_databento.csv"
            end if
        end if
        if (tick_size <= 0.0_dp) error stop "xcsv_to_intraday_tick_stream: tick_size must be positive"
        if (max_obs < 0) error stop "xcsv_to_intraday_tick_stream: max_obs must be nonnegative"
        ! Expand any glob patterns in the positional arguments.
        allocate(expanded(0))
        do i = 1, nfile
            if (scan(raw(i), "*?") > 0) then
                call glob(trim(raw(i)), glob_matches)
                if (size(glob_matches) > 0) then
                    expanded = [character(len=512) :: expanded, glob_matches]
                else
                    expanded = [expanded, raw(i:i)]
                end if
            else
                expanded = [expanded, resolve_filename(raw(i), data_dir)]
            end if
        end do
        if (size(expanded) < 1) error stop "xcsv_to_intraday_tick_stream: no files matched"
        call keep_price_csv_files(expanded, input_files)
        if (size(input_files) < 1) error stop "xcsv_to_intraday_tick_stream: no price CSV files selected"
        deallocate(expanded, raw)
    end subroutine parse_arguments

    ! Keep likely price CSV files, skipping known summary/report files.
    subroutine keep_price_csv_files(files_in, files_out)
        character(len=512), intent(in) :: files_in(:)
        character(len=512), allocatable, intent(out) :: files_out(:)
        character(len=512), allocatable :: tmp(:)
        character(len=512) :: base
        integer :: i, n

        allocate(tmp(size(files_in)))
        n = 0
        do i = 1, size(files_in)
            base = uppercase(basename_without_extension(files_in(i)))
            if (index(base, "SUMMARY") > 0) cycle
            n = n + 1
            tmp(n) = files_in(i)
        end do
        allocate(files_out(n))
        if (n > 0) files_out = tmp(1:n)
        deallocate(tmp)
    end subroutine keep_price_csv_files
end program xcsv_to_intraday_tick_stream

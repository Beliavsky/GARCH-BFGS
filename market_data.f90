! Market data containers and readers for intraday OHLCV bars.

module market_data_mod
    use iso_fortran_env, only: int64
    use kind_mod, only: dp
    use date_mod, only: date_t, date_time_t, date_from_iso, datetime_from_iso, seconds_per_minute, seconds_per_hour, yyyymmdd
    use path_utils_mod, only: has_extension
    use strings_mod, only: split_string, uppercase
    implicit none
    private

    integer, parameter :: default_bar_minutes = 5
    integer, parameter :: default_session_start_seconds = 9*seconds_per_hour + 30*seconds_per_minute ! market opens at 9:30 AM
    integer, parameter :: default_session_end_seconds = 16*seconds_per_hour ! market closes at 4:00 PM
    real(dp), parameter :: default_price_tick_size = 0.001_dp
    character(len=*), parameter :: ohlcv_tick_stream_magic = "OHLCVTICKSTREAM"
    integer, parameter :: ohlcv_tick_stream_version = 1
    integer, parameter :: ohlcv_dp_stream_version = 2

    type, public :: ohlcv_series_t
        type(date_time_t), allocatable :: timestamp(:)
        real(dp), allocatable :: open(:)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: close(:)
        real(dp), allocatable :: volume(:)
    contains
        procedure :: nobs
    end type ohlcv_series_t

    type, public :: ohlcv_tick_series_t
        real(dp) :: tick_size = default_price_tick_size
        type(date_t), allocatable :: date(:)
        integer, allocatable :: seconds_after_midnight(:)
        integer(int64), allocatable :: open_tick(:)
        integer(int64), allocatable :: high_tick(:)
        integer(int64), allocatable :: low_tick(:)
        integer(int64), allocatable :: close_tick(:)
        integer(int64), allocatable :: volume(:)
    contains
        procedure :: nobs => nobs_tick
    end type ohlcv_tick_series_t

    public :: read_intraday_prices_csv
    public :: read_intraday_prices_csv_auto
    public :: read_intraday_prices_file
    public :: read_continuous_intraday_prices_csv
    public :: intraday_bin_ids
    public :: infer_bar_interval_seconds
    public :: resample_ohlcv_series
    public :: check_prices_on_penny_grid
    public :: filter_intraday_session
    public :: intraday_close_changes
    public :: intraday_high_low_ranges
    public :: convert_ohlcv_to_tick_series
    public :: convert_tick_to_ohlcv_series
    public :: write_ohlcv_tick_stream
    public :: read_ohlcv_tick_stream
    public :: write_ohlcv_dp_stream
    public :: read_ohlcv_dp_stream
    public :: peek_ohlcv_stream_version
    public :: ohlcv_tick_stream_version
    public :: ohlcv_dp_stream_version
    public :: equal_ohlcv_tick_series
    public :: default_bar_minutes
    public :: default_session_start_seconds, default_session_end_seconds
    public :: default_price_tick_size

contains

    subroutine read_intraday_prices_csv(filename, series, max_obs)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series
        integer, intent(in), optional :: max_obs
        integer :: io, unit, nrows, i, data_limit
        integer :: idt, iopen, ihigh, ilow, iclose, ivolume
        character(len=4096) :: line
        character(:), allocatable :: tokens(:), header(:)

        data_limit = huge(1)
        if (present(max_obs)) then
            if (max_obs < 1) error stop "read_intraday_prices_csv: max_obs must be positive"
            data_limit = max_obs
        end if

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "market_data_mod: cannot open ", trim(filename)
            error stop
        end if

        read(unit, '(A)', iostat=io) line
        if (io /= 0) error stop "read_intraday_prices_csv: empty file"
        call split_string(line, ",", header)
        idt = find_col_any(header, [character(len=16) :: "Datetime", "ts_event"])
        iopen = find_col(header, "Open")
        ihigh = find_col(header, "High")
        ilow = find_col(header, "Low")
        iclose = find_col(header, "Close")
        ivolume = find_col(header, "Volume")
        if (idt == 0 .or. iopen == 0 .or. ihigh == 0 .or. ilow == 0 .or. iclose == 0 .or. ivolume == 0) then
            error stop "read_intraday_prices_csv: expected Datetime, Open, High, Low, Close, Volume columns"
        end if

        nrows = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            nrows = nrows + 1
            if (nrows >= data_limit) exit
        end do
        if (nrows < 1) error stop "read_intraday_prices_csv: no data rows"

        allocate(series%timestamp(nrows), series%open(nrows), series%high(nrows), &
                 series%low(nrows), series%close(nrows), series%volume(nrows))

        rewind(unit)
        read(unit, '(A)')
        do i = 1, nrows
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            call split_string(line, ",", tokens)
            if (size(tokens) < max(idt, iopen, ihigh, ilow, iclose, ivolume)) then
                error stop "read_intraday_prices_csv: short data row"
            end if
            series%timestamp(i) = datetime_from_iso(tokens(idt))
            read(tokens(iopen), *) series%open(i)
            read(tokens(ihigh), *) series%high(i)
            read(tokens(ilow), *) series%low(i)
            read(tokens(iclose), *) series%close(i)
            read(tokens(ivolume), *) series%volume(i)
        end do
        close(unit)
    end subroutine read_intraday_prices_csv

    ! Read either headered OHLCV CSV or headerless continuous futures CSV.
    subroutine read_intraday_prices_csv_auto(filename, series, max_obs)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series
        integer, intent(in), optional :: max_obs
        integer :: io, unit
        character(len=4096) :: line
        character(:), allocatable :: tokens(:)

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "market_data_mod: cannot open ", trim(filename)
            error stop
        end if
        read(unit, '(A)', iostat=io) line
        close(unit)
        if (io /= 0) error stop "read_intraday_prices_csv_auto: empty file"
        call split_string(line, ",", tokens)
        if (size(tokens) >= 7 .and. yyyymmdd(date_from_iso(tokens(1))) > 0) then
            call read_continuous_intraday_prices_csv(filename, series, max_obs)
        else
            call read_intraday_prices_csv(filename, series, max_obs)
        end if
    end subroutine read_intraday_prices_csv_auto

    ! Read intraday OHLCV from CSV or compact unformatted stream based on suffix.
    ! Binary files are auto-detected as v1 (tick) or v2 (dp) by the version field.
    subroutine read_intraday_prices_file(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series

        if (has_extension(filename, ".bin")) then
            if (peek_ohlcv_stream_version(filename) == ohlcv_dp_stream_version) then
                call read_ohlcv_dp_stream(filename, series)
            else
                block
                    type(ohlcv_tick_series_t) :: tick_series
                    call read_ohlcv_tick_stream(filename, tick_series)
                    call convert_tick_to_ohlcv_series(tick_series, series)
                end block
            end if
        else
            call read_intraday_prices_csv_auto(filename, series)
        end if
    end subroutine read_intraday_prices_file


    ! Read headerless continuous futures bars: date,time,open,high,low,close,volume,symbol.
    subroutine read_continuous_intraday_prices_csv(filename, series, max_obs)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series
        integer, intent(in), optional :: max_obs
        integer :: io, unit, nrows, i, data_limit
        character(len=4096) :: line
        character(:), allocatable :: tokens(:)

        data_limit = huge(1)
        if (present(max_obs)) then
            if (max_obs < 1) error stop "read_continuous_intraday_prices_csv: max_obs must be positive"
            data_limit = max_obs
        end if

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "market_data_mod: cannot open ", trim(filename)
            error stop
        end if

        nrows = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0) exit
            if (trim(line) == "") cycle
            nrows = nrows + 1
            if (nrows >= data_limit) exit
        end do
        if (nrows < 1) error stop "read_continuous_intraday_prices_csv: no data rows"

        allocate(series%timestamp(nrows), series%open(nrows), series%high(nrows), &
                 series%low(nrows), series%close(nrows), series%volume(nrows))

        rewind(unit)
        i = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0) exit
            if (trim(line) == "") cycle
            i = i + 1
            call split_string(line, ",", tokens)
            if (size(tokens) < 7) error stop "read_continuous_intraday_prices_csv: short data row"
            series%timestamp(i)%date = date_from_iso(tokens(1))
            if (yyyymmdd(series%timestamp(i)%date) == 0) then
                error stop "read_continuous_intraday_prices_csv: invalid date"
            end if
            call set_time_from_hhmmss(tokens(2), series%timestamp(i))
            read(tokens(3), *) series%open(i)
            read(tokens(4), *) series%high(i)
            read(tokens(5), *) series%low(i)
            read(tokens(6), *) series%close(i)
            read(tokens(7), *) series%volume(i)
            if (i >= nrows) exit
        end do
        close(unit)
    end subroutine read_continuous_intraday_prices_csv

    pure subroutine intraday_bin_ids(series, bin_id, session_start_seconds, interval_seconds)
        type(ohlcv_series_t), intent(in) :: series
        integer, allocatable, intent(out) :: bin_id(:)
        integer, intent(in), optional :: session_start_seconds, interval_seconds
        integer :: i, start_sec, step_sec, sec

        start_sec = default_session_start_seconds
        step_sec = default_bar_minutes * seconds_per_minute
        if (present(session_start_seconds)) start_sec = session_start_seconds
        if (present(interval_seconds)) step_sec = interval_seconds
        if (step_sec < 1) error stop "intraday_bin_ids: interval_seconds must be positive"
        allocate(bin_id(series%nobs()))
        do i = 1, size(bin_id)
            sec = series%timestamp(i)%seconds_since_midnight()
            bin_id(i) = (sec - start_sec) / step_sec + 1
        end do
    end subroutine intraday_bin_ids

    ! Infer the bar interval as the minimum positive same-day gap across all timestamps.
    ! Using the minimum (not the first gap found) correctly handles missing bars,
    ! which create gaps that are multiples of the true interval.
    pure integer function infer_bar_interval_seconds(series) result(interval_seconds)
        type(ohlcv_series_t), intent(in) :: series
        integer :: i, s0, s1

        interval_seconds = default_bar_minutes * seconds_per_minute
        do i = 2, series%nobs()
            if (yyyymmdd(series%timestamp(i)%date) == yyyymmdd(series%timestamp(i - 1)%date)) then
                s0 = series%timestamp(i - 1)%seconds_since_midnight()
                s1 = series%timestamp(i)%seconds_since_midnight()
                if (s1 > s0) interval_seconds = min(interval_seconds, s1 - s0)
            end if
        end do
    end function infer_bar_interval_seconds

    ! Aggregate OHLCV bars from source_seconds to a lower frequency target_seconds.
    subroutine resample_ohlcv_series(series, resampled, source_seconds, target_seconds, session_start_seconds)
        type(ohlcv_series_t), intent(in) :: series
        type(ohlcv_series_t), intent(out) :: resampled
        integer, intent(in) :: source_seconds, target_seconds
        integer, intent(in), optional :: session_start_seconds
        type(date_time_t), allocatable :: timestamp_tmp(:)
        real(dp), allocatable :: open_tmp(:), high_tmp(:), low_tmp(:), close_tmp(:), volume_tmp(:)
        integer :: start_sec, i, k, sec, bin_id, current_bin
        integer :: bucket_start_sec
        type(date_time_t) :: bucket_timestamp

        if (source_seconds < 1) error stop "resample_ohlcv_series: source_seconds must be positive"
        if (target_seconds < 1) error stop "resample_ohlcv_series: target_seconds must be positive"
        if (target_seconds < source_seconds) then
            error stop "resample_ohlcv_series: target_seconds must be >= source_seconds"
        end if
        if (mod(target_seconds, source_seconds) /= 0) then
            error stop "resample_ohlcv_series: target_seconds must be an integer multiple of source_seconds"
        end if
        if (series%nobs() < 1) error stop "resample_ohlcv_series: empty input series"

        start_sec = default_session_start_seconds
        if (present(session_start_seconds)) start_sec = session_start_seconds

        allocate(timestamp_tmp(series%nobs()), open_tmp(series%nobs()), high_tmp(series%nobs()), &
                 low_tmp(series%nobs()), close_tmp(series%nobs()), volume_tmp(series%nobs()))

        k = 0
        current_bin = huge(1)
        do i = 1, series%nobs()
            sec = series%timestamp(i)%seconds_since_midnight()
            bin_id = (sec - start_sec) / target_seconds + 1
            if (sec < start_sec) bin_id = bin_id - 1
            if (k == 0) then
                k = k + 1
                current_bin = bin_id
                bucket_start_sec = start_sec + (bin_id - 1)*target_seconds
                bucket_timestamp = series%timestamp(i)
                call set_time_from_seconds(bucket_timestamp, bucket_start_sec)
                timestamp_tmp(k) = bucket_timestamp
                open_tmp(k) = series%open(i)
                high_tmp(k) = series%high(i)
                low_tmp(k) = series%low(i)
                close_tmp(k) = series%close(i)
                volume_tmp(k) = series%volume(i)
            else if (yyyymmdd(series%timestamp(i)%date) /= yyyymmdd(timestamp_tmp(k)%date) .or. &
                     bin_id /= current_bin) then
                k = k + 1
                current_bin = bin_id
                bucket_start_sec = start_sec + (bin_id - 1)*target_seconds
                bucket_timestamp = series%timestamp(i)
                call set_time_from_seconds(bucket_timestamp, bucket_start_sec)
                timestamp_tmp(k) = bucket_timestamp
                open_tmp(k) = series%open(i)
                high_tmp(k) = series%high(i)
                low_tmp(k) = series%low(i)
                close_tmp(k) = series%close(i)
                volume_tmp(k) = series%volume(i)
            else
                high_tmp(k) = max(high_tmp(k), series%high(i))
                low_tmp(k) = min(low_tmp(k), series%low(i))
                close_tmp(k) = series%close(i)
                volume_tmp(k) = volume_tmp(k) + series%volume(i)
            end if
        end do

        allocate(resampled%timestamp(k), resampled%open(k), resampled%high(k), &
                 resampled%low(k), resampled%close(k), resampled%volume(k))
        resampled%timestamp = timestamp_tmp(1:k)
        resampled%open = open_tmp(1:k)
        resampled%high = high_tmp(1:k)
        resampled%low = low_tmp(1:k)
        resampled%close = close_tmp(1:k)
        resampled%volume = volume_tmp(1:k)
        deallocate(timestamp_tmp, open_tmp, high_tmp, low_tmp, close_tmp, volume_tmp)
    end subroutine resample_ohlcv_series

    ! Check whether all OHLC prices are integer pennies within an optional tolerance.
    subroutine check_prices_on_penny_grid(series, ok, n_bad, max_abs_error, tolerance)
        type(ohlcv_series_t), intent(in) :: series
        logical, intent(out) :: ok(4)
        integer, intent(out) :: n_bad(4)
        real(dp), intent(out) :: max_abs_error(4)
        real(dp), intent(in), optional :: tolerance
        real(dp) :: tol
        integer :: i

        tol = 1.0e-8_dp
        if (present(tolerance)) tol = tolerance
        n_bad = 0
        max_abs_error = 0.0_dp
        do i = 1, series%nobs()
            call update_penny_grid_error(series%open(i), tol, n_bad(1), max_abs_error(1))
            call update_penny_grid_error(series%high(i), tol, n_bad(2), max_abs_error(2))
            call update_penny_grid_error(series%low(i), tol, n_bad(3), max_abs_error(3))
            call update_penny_grid_error(series%close(i), tol, n_bad(4), max_abs_error(4))
        end do
        ok = n_bad == 0
    end subroutine check_prices_on_penny_grid

    subroutine filter_intraday_session(series, filtered, session_start_seconds, session_end_seconds, include_end)
        type(ohlcv_series_t), intent(in) :: series
        type(ohlcv_series_t), intent(out) :: filtered
        integer, intent(in), optional :: session_start_seconds, session_end_seconds
        logical, intent(in), optional :: include_end
        integer :: start_sec, end_sec, i, j, n_keep, sec
        logical :: include_end_

        start_sec = default_session_start_seconds
        end_sec = default_session_end_seconds
        include_end_ = .false.
        if (present(session_start_seconds)) start_sec = session_start_seconds
        if (present(session_end_seconds)) end_sec = session_end_seconds
        if (present(include_end)) include_end_ = include_end
        if (end_sec < start_sec) error stop "filter_intraday_session: end before start"

        n_keep = 0
        do i = 1, series%nobs()
            sec = series%timestamp(i)%seconds_since_midnight()
            if (in_session(sec, start_sec, end_sec, include_end_)) n_keep = n_keep + 1
        end do
        allocate(filtered%timestamp(n_keep), filtered%open(n_keep), filtered%high(n_keep), &
                 filtered%low(n_keep), filtered%close(n_keep), filtered%volume(n_keep))
        j = 0
        do i = 1, series%nobs()
            sec = series%timestamp(i)%seconds_since_midnight()
            if (.not. in_session(sec, start_sec, end_sec, include_end_)) cycle
            j = j + 1
            filtered%timestamp(j) = series%timestamp(i)
            filtered%open(j) = series%open(i)
            filtered%high(j) = series%high(i)
            filtered%low(j) = series%low(i)
            filtered%close(j) = series%close(i)
            filtered%volume(j) = series%volume(i)
        end do
    end subroutine filter_intraday_session

    ! Compute within-day signed and absolute close-to-close price changes.
    subroutine intraday_close_changes(series, signed_change, abs_change, nchange)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), allocatable, intent(out) :: signed_change(:), abs_change(:)
        integer, intent(out) :: nchange
        integer :: i

        allocate(signed_change(max(series%nobs() - 1, 0)), abs_change(max(series%nobs() - 1, 0)))
        nchange = 0
        do i = 2, series%nobs()
            if (yyyymmdd(series%timestamp(i)%date) /= yyyymmdd(series%timestamp(i - 1)%date)) cycle
            nchange = nchange + 1
            signed_change(nchange) = series%close(i) - series%close(i - 1)
            abs_change(nchange) = abs(signed_change(nchange))
        end do
    end subroutine intraday_close_changes

    ! Compute the high-low price range for each OHLCV bar.
    subroutine intraday_high_low_ranges(series, price_range, nrange)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), allocatable, intent(out) :: price_range(:)
        integer, intent(out) :: nrange

        nrange = series%nobs()
        allocate(price_range(nrange))
        price_range = max(series%high - series%low, 0.0_dp)
    end subroutine intraday_high_low_ranges

    ! Convert real-valued OHLCV bars to an integer tick-price representation.
    subroutine convert_ohlcv_to_tick_series(series, tick_series, tick_size)
        type(ohlcv_series_t), intent(in) :: series
        type(ohlcv_tick_series_t), intent(out) :: tick_series
        real(dp), intent(in), optional :: tick_size
        real(dp) :: scale
        integer :: i, n

        tick_series%tick_size = default_price_tick_size
        if (present(tick_size)) tick_series%tick_size = tick_size
        if (tick_series%tick_size <= 0.0_dp) error stop "convert_ohlcv_to_tick_series: tick_size must be positive"

        n = series%nobs()
        scale = 1.0_dp / tick_series%tick_size
        allocate(tick_series%date(n), tick_series%seconds_after_midnight(n), &
                 tick_series%open_tick(n), tick_series%high_tick(n), tick_series%low_tick(n), &
                 tick_series%close_tick(n), tick_series%volume(n))
        do i = 1, n
            tick_series%date(i) = series%timestamp(i)%date
            tick_series%seconds_after_midnight(i) = series%timestamp(i)%seconds_since_midnight()
            tick_series%open_tick(i) = nint(scale * series%open(i), int64)
            tick_series%high_tick(i) = nint(scale * series%high(i), int64)
            tick_series%low_tick(i) = nint(scale * series%low(i), int64)
            tick_series%close_tick(i) = nint(scale * series%close(i), int64)
            tick_series%volume(i) = nint(series%volume(i), int64)
        end do
    end subroutine convert_ohlcv_to_tick_series

    ! Convert an integer tick OHLCV series to the real-valued OHLCV representation.
    subroutine convert_tick_to_ohlcv_series(tick_series, series)
        type(ohlcv_tick_series_t), intent(in) :: tick_series
        type(ohlcv_series_t), intent(out) :: series
        integer :: i, n, sec

        if (tick_series%tick_size <= 0.0_dp) error stop "convert_tick_to_ohlcv_series: tick_size must be positive"
        n = tick_series%nobs()
        allocate(series%timestamp(n), series%open(n), series%high(n), &
                 series%low(n), series%close(n), series%volume(n))
        do i = 1, n
            sec = tick_series%seconds_after_midnight(i)
            if (sec < 0 .or. sec >= 24*seconds_per_hour) then
                error stop "convert_tick_to_ohlcv_series: seconds_after_midnight outside current date"
            end if
            series%timestamp(i)%date = tick_series%date(i)
            series%timestamp(i)%hour = sec / seconds_per_hour
            series%timestamp(i)%minute = mod(sec, seconds_per_hour) / seconds_per_minute
            series%timestamp(i)%second = mod(sec, seconds_per_minute)
            series%open(i) = tick_series%tick_size * real(tick_series%open_tick(i), dp)
            series%high(i) = tick_series%tick_size * real(tick_series%high_tick(i), dp)
            series%low(i) = tick_series%tick_size * real(tick_series%low_tick(i), dp)
            series%close(i) = tick_series%tick_size * real(tick_series%close_tick(i), dp)
            series%volume(i) = real(tick_series%volume(i), dp)
        end do
    end subroutine convert_tick_to_ohlcv_series

    ! Save an integer tick OHLCV series using an explicit unformatted stream layout.
    subroutine write_ohlcv_tick_stream(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_tick_series_t), intent(in) :: series
        integer :: io, unit, n, i
        integer, allocatable :: year(:), month(:), day(:)
        character(len=16) :: magic

        n = series%nobs()
        if (n < 1) error stop "write_ohlcv_tick_stream: empty input series"
        magic = ohlcv_tick_stream_magic
        allocate(year(n), month(n), day(n))
        do i = 1, n
            year(i) = series%date(i)%year
            month(i) = series%date(i)%month
            day(i) = series%date(i)%day
        end do

        open(newunit=unit, file=filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=io)
        if (io /= 0) then
            print '(A,A)', "write_ohlcv_tick_stream: cannot open ", trim(filename)
            error stop
        end if

        write(unit) magic
        write(unit) ohlcv_tick_stream_version
        write(unit) n
        write(unit) series%tick_size
        write(unit) year
        write(unit) month
        write(unit) day
        write(unit) series%seconds_after_midnight
        write(unit) series%open_tick
        write(unit) series%high_tick
        write(unit) series%low_tick
        write(unit) series%close_tick
        write(unit) series%volume
        close(unit)
        deallocate(year, month, day)
    end subroutine write_ohlcv_tick_stream

    ! Read an integer tick OHLCV series written by write_ohlcv_tick_stream.
    subroutine read_ohlcv_tick_stream(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_tick_series_t), intent(out) :: series
        integer :: io, unit, version, n, i
        integer, allocatable :: year(:), month(:), day(:)
        character(len=16) :: magic

        open(newunit=unit, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "read_ohlcv_tick_stream: cannot open ", trim(filename)
            error stop
        end if

        read(unit, iostat=io) magic
        if (io /= 0 .or. trim(magic) /= ohlcv_tick_stream_magic) then
            error stop "read_ohlcv_tick_stream: invalid file magic"
        end if
        read(unit) version
        if (version /= ohlcv_tick_stream_version) then
            error stop "read_ohlcv_tick_stream: unsupported file version"
        end if
        read(unit) n
        if (n < 1) error stop "read_ohlcv_tick_stream: invalid observation count"
        read(unit) series%tick_size

        allocate(year(n), month(n), day(n), series%date(n), series%seconds_after_midnight(n), &
                 series%open_tick(n), series%high_tick(n), series%low_tick(n), &
                 series%close_tick(n), series%volume(n))
        read(unit) year
        read(unit) month
        read(unit) day
        read(unit) series%seconds_after_midnight
        read(unit) series%open_tick
        read(unit) series%high_tick
        read(unit) series%low_tick
        read(unit) series%close_tick
        read(unit) series%volume
        close(unit)

        do i = 1, n
            series%date(i)%year = year(i)
            series%date(i)%month = month(i)
            series%date(i)%day = day(i)
        end do
        deallocate(year, month, day)
    end subroutine read_ohlcv_tick_stream

    ! Return true if two integer tick OHLCV series have identical metadata and columns.
    pure logical function equal_ohlcv_tick_series(a, b)
        type(ohlcv_tick_series_t), intent(in) :: a, b
        integer :: i, n

        equal_ohlcv_tick_series = .false.
        n = a%nobs()
        if (n /= b%nobs()) return
        if (a%tick_size /= b%tick_size) return
        if (allocated(a%date) .neqv. allocated(b%date)) return
        if (allocated(a%seconds_after_midnight) .neqv. allocated(b%seconds_after_midnight)) return
        if (allocated(a%open_tick) .neqv. allocated(b%open_tick)) return
        if (allocated(a%high_tick) .neqv. allocated(b%high_tick)) return
        if (allocated(a%low_tick) .neqv. allocated(b%low_tick)) return
        if (allocated(a%close_tick) .neqv. allocated(b%close_tick)) return
        if (allocated(a%volume) .neqv. allocated(b%volume)) return
        if (n == 0) then
            equal_ohlcv_tick_series = .true.
            return
        end if
        if (any(a%seconds_after_midnight /= b%seconds_after_midnight)) return
        if (any(a%open_tick /= b%open_tick)) return
        if (any(a%high_tick /= b%high_tick)) return
        if (any(a%low_tick /= b%low_tick)) return
        if (any(a%close_tick /= b%close_tick)) return
        if (any(a%volume /= b%volume)) return
        do i = 1, n
            if (a%date(i)%year /= b%date(i)%year) return
            if (a%date(i)%month /= b%date(i)%month) return
            if (a%date(i)%day /= b%date(i)%day) return
        end do
        equal_ohlcv_tick_series = .true.
    end function equal_ohlcv_tick_series

    pure elemental integer function nobs(this)
        class(ohlcv_series_t), intent(in) :: this

        if (allocated(this%close)) then
            nobs = size(this%close)
        else
            nobs = 0
        end if
    end function nobs

    pure elemental integer function nobs_tick(this)
        class(ohlcv_tick_series_t), intent(in) :: this

        if (allocated(this%close_tick)) then
            nobs_tick = size(this%close_tick)
        else
            nobs_tick = 0
        end if
    end function nobs_tick

    pure logical function in_session(sec, start_sec, end_sec, include_end)
        integer, intent(in) :: sec, start_sec, end_sec
        logical, intent(in) :: include_end

        if (include_end) then
            in_session = sec >= start_sec .and. sec <= end_sec
        else
            in_session = sec >= start_sec .and. sec < end_sec
        end if
    end function in_session

    pure subroutine update_penny_grid_error(price, tolerance, n_bad, max_abs_error)
        real(dp), intent(in) :: price, tolerance
        integer, intent(inout) :: n_bad
        real(dp), intent(inout) :: max_abs_error
        real(dp) :: err

        err = abs(100.0_dp*price - anint(100.0_dp*price))
        max_abs_error = max(max_abs_error, err)
        if (err > tolerance) n_bad = n_bad + 1
    end subroutine update_penny_grid_error

    pure subroutine set_time_from_seconds(timestamp, seconds)
        type(date_time_t), intent(inout) :: timestamp
        integer, intent(in) :: seconds

        if (seconds < 0 .or. seconds >= 24*seconds_per_hour) then
            error stop "set_time_from_seconds: seconds outside current date"
        end if
        timestamp%hour = seconds / seconds_per_hour
        timestamp%minute = mod(seconds, seconds_per_hour) / seconds_per_minute
        timestamp%second = mod(seconds, seconds_per_minute)
    end subroutine set_time_from_seconds

    subroutine set_time_from_hhmmss(s, timestamp)
        character(len=*), intent(in) :: s
        type(date_time_t), intent(inout) :: timestamp
        character(len=len(s)) :: t
        integer :: hh, mm, ss

        t = adjustl(s)
        if (len_trim(t) < 5) error stop "set_time_from_hhmmss: time too short"
        read(t(1:2), *) hh
        read(t(4:5), *) mm
        ss = 0
        if (len_trim(t) >= 8) read(t(7:8), *) ss
        if (hh < 0 .or. hh > 23 .or. mm < 0 .or. mm > 59 .or. ss < 0 .or. ss > 59) then
            error stop "set_time_from_hhmmss: invalid time"
        end if
        timestamp%hour = hh
        timestamp%minute = mm
        timestamp%second = ss
    end subroutine set_time_from_hhmmss


    pure integer function find_col(names, target)
        character(len=*), intent(in) :: names(:), target

        find_col = findloc(uppercase(names), uppercase(target), dim=1)
    end function find_col

    pure integer function find_col_any(names, targets)
        character(len=*), intent(in) :: names(:), targets(:)
        integer :: i

        find_col_any = 0
        do i = 1, size(targets)
            find_col_any = find_col(names, targets(i))
            if (find_col_any > 0) return
        end do
    end function find_col_any

    ! Read only the version integer from a binary OHLCV stream without loading data.
    integer function peek_ohlcv_stream_version(filename)
        character(len=*), intent(in) :: filename
        integer :: io, unit, version
        character(len=16) :: magic

        open(newunit=unit, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=io)
        if (io /= 0) error stop "peek_ohlcv_stream_version: cannot open file"
        read(unit, iostat=io) magic
        if (io /= 0 .or. trim(magic) /= ohlcv_tick_stream_magic) &
            error stop "peek_ohlcv_stream_version: invalid file magic"
        read(unit) version
        close(unit)
        peek_ohlcv_stream_version = version
    end function peek_ohlcv_stream_version

    ! Save an ohlcv_series_t directly as real(dp) (version 2), no tick_size needed.
    subroutine write_ohlcv_dp_stream(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(in) :: series
        integer :: io, unit, n, i
        integer, allocatable :: year(:), month(:), day(:), seconds(:)
        character(len=16) :: magic

        n = series%nobs()
        if (n < 1) error stop "write_ohlcv_dp_stream: empty input series"
        magic = ohlcv_tick_stream_magic
        allocate(year(n), month(n), day(n), seconds(n))
        do i = 1, n
            year(i) = series%timestamp(i)%date%year
            month(i) = series%timestamp(i)%date%month
            day(i) = series%timestamp(i)%date%day
            seconds(i) = series%timestamp(i)%seconds_since_midnight()
        end do

        open(newunit=unit, file=filename, access='stream', form='unformatted', &
             status='replace', action='write', iostat=io)
        if (io /= 0) then
            print '(A,A)', "write_ohlcv_dp_stream: cannot open ", trim(filename)
            error stop
        end if

        write(unit) magic
        write(unit) ohlcv_dp_stream_version
        write(unit) n
        write(unit) year
        write(unit) month
        write(unit) day
        write(unit) seconds
        write(unit) series%open
        write(unit) series%high
        write(unit) series%low
        write(unit) series%close
        write(unit) series%volume
        close(unit)
        deallocate(year, month, day, seconds)
    end subroutine write_ohlcv_dp_stream

    ! Read a double-precision OHLCV series written by write_ohlcv_dp_stream.
    subroutine read_ohlcv_dp_stream(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series
        integer :: io, unit, version, n, i
        integer, allocatable :: year(:), month(:), day(:), seconds(:)
        character(len=16) :: magic

        open(newunit=unit, file=filename, access='stream', form='unformatted', &
             status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "read_ohlcv_dp_stream: cannot open ", trim(filename)
            error stop
        end if

        read(unit, iostat=io) magic
        if (io /= 0 .or. trim(magic) /= ohlcv_tick_stream_magic) &
            error stop "read_ohlcv_dp_stream: invalid file magic"
        read(unit) version
        if (version /= ohlcv_dp_stream_version) &
            error stop "read_ohlcv_dp_stream: unsupported file version"
        read(unit) n
        if (n < 1) error stop "read_ohlcv_dp_stream: invalid observation count"

        allocate(year(n), month(n), day(n), seconds(n))
        allocate(series%timestamp(n), series%open(n), series%high(n), &
                 series%low(n), series%close(n), series%volume(n))
        read(unit) year
        read(unit) month
        read(unit) day
        read(unit) seconds
        read(unit) series%open
        read(unit) series%high
        read(unit) series%low
        read(unit) series%close
        read(unit) series%volume
        close(unit)

        do i = 1, n
            series%timestamp(i)%date%year = year(i)
            series%timestamp(i)%date%month = month(i)
            series%timestamp(i)%date%day = day(i)
            call set_time_from_seconds(series%timestamp(i), seconds(i))
        end do
        deallocate(year, month, day, seconds)
    end subroutine read_ohlcv_dp_stream

end module market_data_mod

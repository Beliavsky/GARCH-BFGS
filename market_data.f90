! Market data containers and readers for intraday OHLCV bars.

module market_data_mod
    use kind_mod, only: dp
    use date_mod, only: date_time_t, datetime_from_iso, seconds_per_minute, seconds_per_hour
    use strings_mod, only: split_string, uppercase
    implicit none
    private

    integer, parameter :: default_bar_minutes = 5
    integer, parameter :: default_session_start_seconds = 9*seconds_per_hour + 30*seconds_per_minute ! market opens at 9:30 AM
    integer, parameter :: default_session_end_seconds = 16*seconds_per_hour ! market closes at 4:00 PM

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

    public :: read_intraday_prices_csv
    public :: intraday_bin_ids
    public :: filter_intraday_session
    public :: default_bar_minutes
    public :: default_session_start_seconds, default_session_end_seconds

contains

    subroutine read_intraday_prices_csv(filename, series)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t), intent(out) :: series
        integer :: io, unit, nrows, i
        integer :: idt, iopen, ihigh, ilow, iclose, ivolume
        character(len=4096) :: line
        character(:), allocatable :: tokens(:), header(:)

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "market_data_mod: cannot open ", trim(filename)
            error stop
        end if

        read(unit, '(A)', iostat=io) line
        if (io /= 0) error stop "read_intraday_prices_csv: empty file"
        call split_string(line, ",", header)
        idt = find_col(header, "Datetime")
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

    pure elemental integer function nobs(this)
        class(ohlcv_series_t), intent(in) :: this

        if (allocated(this%close)) then
            nobs = size(this%close)
        else
            nobs = 0
        end if
    end function nobs

    pure logical function in_session(sec, start_sec, end_sec, include_end)
        integer, intent(in) :: sec, start_sec, end_sec
        logical, intent(in) :: include_end

        if (include_end) then
            in_session = sec >= start_sec .and. sec <= end_sec
        else
            in_session = sec >= start_sec .and. sec < end_sec
        end if
    end function in_session

    pure integer function find_col(names, target)
        character(len=*), intent(in) :: names(:), target

        find_col = findloc(uppercase(names), uppercase(target), dim=1)
    end function find_col

end module market_data_mod

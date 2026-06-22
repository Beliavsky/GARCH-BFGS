! CSV reading utilities for financial price data.
!
! split_string is copied verbatim from util.f90 in the public-domain
! Beliavsky/DataFrame library (github.com/Beliavsky/DataFrame).
!
! read_price_csv reads a CSV whose first column is a date in YYYY-MM-DD
! format and whose remaining columns are real-valued prices or levels.
! Dates are stored as integers in YYYYMMDD format (dashes removed).
!
! log_returns(prices) returns log(p_t / p_{t-1}).

module csv_mod
    use kind_mod, only: dp
    use strings_mod, only: split_string, uppercase
    implicit none
    private
    public :: read_price_csv, read_ohlc_csv, read_numeric_csv, log_returns, print_price_sample_info, date_label

contains

    ! ------------------------------------------------------------------
    ! Convert a date string in YYYY-MM-DD format to an integer YYYYMMDD.
    ! ------------------------------------------------------------------
    integer function date_to_int(s)
        character(len=*), intent(in) :: s
        character(len=8) :: s8
        s8 = s(1:4) // s(6:7) // s(9:10)
        read(s8, *) date_to_int
    end function date_to_int

    ! ------------------------------------------------------------------
    ! Read a CSV file with:
    !   - a header row (first field ignored, remaining fields are col names)
    !   - a date column (YYYY-MM-DD) stored as integer YYYYMMDD in dates(:)
    !   - ncols real-valued columns stored in values(:, 1:ncols)
    ! col_names is allocated to size ncols with entries trimmed from the header.
    ! ------------------------------------------------------------------
    subroutine read_price_csv(filename, dates, col_names, values, selected_cols, max_col)
        character(len=*),              intent(in)  :: filename
        integer,          allocatable, intent(out) :: dates(:)
        character(len=32), allocatable, intent(out) :: col_names(:)
        real(dp),         allocatable, intent(out) :: values(:,:)
        character(len=*), optional, intent(in) :: selected_cols(:)
        integer, optional, intent(in) :: max_col

        integer :: io, unit, i, j, nrows, ncols, nfile_cols, nrequested
        integer, allocatable :: selected_idx(:)
        character(len=32), allocatable :: file_col_names(:)
        character(len=1024) :: line
        character(:), allocatable :: tokens(:)

        ! 1) Open file.
        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "csv_mod: cannot open ", trim(filename)
            error stop
        end if

        ! 2) Read header: first token is date label (ignored), rest are column names.
        read(unit, '(A)') line
        call split_string(line, ",", tokens)
        nfile_cols = size(tokens) - 1
        if (nfile_cols < 1) then
            print '(A)', "csv_mod: no data columns in header"
            error stop
        end if
        allocate(file_col_names(nfile_cols))
        do i = 1, nfile_cols
            file_col_names(i) = adjustl(tokens(i+1))
        end do

        if (present(selected_cols)) then
            nrequested = size(selected_cols)
            if (nrequested < 1) then
                print '(A)', "csv_mod: selected_cols is empty"
                error stop
            end if
            ncols = nrequested
            if (present(max_col)) then
                if (max_col < 1) then
                    print '(A)', "csv_mod: max_col must be positive"
                    error stop
                end if
                ncols = min(ncols, max_col)
            end if
            allocate(selected_idx(ncols), col_names(ncols))
            do i = 1, ncols
                selected_idx(i) = find_col(file_col_names, selected_cols(i))
                if (selected_idx(i) == 0) then
                    print '(A,A,A,A)', "csv_mod: selected column ", trim(selected_cols(i)), &
                        " not found in ", trim(filename)
                    error stop
                end if
                col_names(i) = adjustl(file_col_names(selected_idx(i)))
            end do
        else
            ncols = nfile_cols
            if (present(max_col)) then
                if (max_col < 1) then
                    print '(A)', "csv_mod: max_col must be positive"
                    error stop
                end if
                ncols = min(ncols, max_col)
            end if
            allocate(selected_idx(ncols), col_names(ncols))
            do i = 1, ncols
                selected_idx(i) = i
                col_names(i) = file_col_names(i)
            end do
        end if

        ! 3) Count data rows (skip blank lines at end of file).
        nrows = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            nrows = nrows + 1
        end do
        if (nrows == 0) then
            print '(A)', "csv_mod: no data rows found"
            error stop
        end if

        ! 4) Rewind, skip header, allocate, read data.
        rewind(unit)
        read(unit, '(A)')   ! skip header

        allocate(dates(nrows), values(nrows, ncols))

        do i = 1, nrows
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            call split_string(line, ",", tokens)
            dates(i) = date_to_int(tokens(1))
            do j = 1, ncols
                read(tokens(selected_idx(j)+1), *) values(i, j)
            end do
        end do

        close(unit)
    end subroutine read_price_csv

    subroutine read_ohlc_csv(filename, dates, asset_names, open_prices, close_prices, selected_assets, high_prices, low_prices)
        character(len=*),               intent(in)  :: filename
        integer,           allocatable, intent(out) :: dates(:)
        character(len=32), allocatable, intent(out) :: asset_names(:)
        real(dp),          allocatable, intent(out) :: open_prices(:,:), close_prices(:,:)
        character(len=*), optional, intent(in) :: selected_assets(:)
        real(dp), optional, allocatable, intent(out) :: high_prices(:,:), low_prices(:,:)

        integer :: io, unit, i, j, nrows, nassets, nfile_cols, iopen, ihigh, ilow, iclose, iadj
        integer, allocatable :: open_idx(:), high_idx(:), low_idx(:), close_idx(:), adj_idx(:)
        real(dp) :: open_raw, high_raw, low_raw, close_raw, adj_close, adj_factor
        character(len=4096) :: line
        character(len=32), allocatable :: file_symbols(:), file_fields(:)
        character(:), allocatable :: symbol_tokens(:), field_tokens(:), tokens(:)

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "csv_mod: cannot open ", trim(filename)
            error stop
        end if

        read(unit, '(A)') line
        call split_string(line, ",", symbol_tokens)
        read(unit, '(A)') line
        call split_string(line, ",", field_tokens)
        read(unit, '(A)') line

        nfile_cols = size(symbol_tokens) - 1
        if (nfile_cols < 1 .or. size(field_tokens) - 1 /= nfile_cols) then
            print '(A)', "csv_mod: invalid OHLC header"
            error stop
        end if
        allocate(file_symbols(nfile_cols), file_fields(nfile_cols))
        do i = 1, nfile_cols
            file_symbols(i) = adjustl(symbol_tokens(i+1))
            file_fields(i) = adjustl(field_tokens(i+1))
        end do

        if (present(selected_assets)) then
            nassets = size(selected_assets)
            allocate(asset_names(nassets))
            do i = 1, nassets
                asset_names(i) = adjustl(selected_assets(i))
            end do
        else
            nassets = count_unique_symbols(file_symbols)
            allocate(asset_names(nassets))
            call unique_symbols(file_symbols, asset_names)
        end if

        allocate(open_idx(nassets), high_idx(nassets), low_idx(nassets), close_idx(nassets), adj_idx(nassets))
        do i = 1, nassets
            iopen  = find_ohlc_col(file_symbols, file_fields, asset_names(i), "Open")
            ihigh  = find_ohlc_col(file_symbols, file_fields, asset_names(i), "High")
            ilow   = find_ohlc_col(file_symbols, file_fields, asset_names(i), "Low")
            iclose = find_ohlc_col(file_symbols, file_fields, asset_names(i), "Close")
            iadj   = find_ohlc_col(file_symbols, file_fields, asset_names(i), "Adj Close")
            if (iopen == 0 .or. ihigh == 0 .or. ilow == 0 .or. iclose == 0 .or. iadj == 0) then
                print '(A,A,A)', "csv_mod: incomplete OHLC fields for ", trim(asset_names(i)), " in ", trim(filename)
                error stop
            end if
            open_idx(i) = iopen
            high_idx(i) = ihigh
            low_idx(i) = ilow
            close_idx(i) = iclose
            adj_idx(i) = iadj
        end do

        nrows = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            nrows = nrows + 1
        end do
        if (nrows == 0) then
            print '(A)', "csv_mod: no OHLC data rows found"
            error stop
        end if

        rewind(unit)
        read(unit, '(A)')
        read(unit, '(A)')
        read(unit, '(A)')
        allocate(dates(nrows), open_prices(nrows,nassets), close_prices(nrows,nassets))
        if (present(high_prices)) allocate(high_prices(nrows,nassets))
        if (present(low_prices)) allocate(low_prices(nrows,nassets))
        do i = 1, nrows
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            call split_string(line, ",", tokens)
            dates(i) = date_to_int(tokens(1))
            do j = 1, nassets
                read(tokens(open_idx(j)+1), *) open_raw
                read(tokens(high_idx(j)+1), *) high_raw
                read(tokens(low_idx(j)+1), *) low_raw
                read(tokens(close_idx(j)+1), *) close_raw
                read(tokens(adj_idx(j)+1), *) adj_close
                adj_factor = adj_close / close_raw
                open_prices(i,j) = open_raw * adj_factor
                if (present(high_prices)) high_prices(i,j) = high_raw * adj_factor
                if (present(low_prices)) low_prices(i,j) = low_raw * adj_factor
                close_prices(i,j) = adj_close
            end do
        end do
        close(unit)
    end subroutine read_ohlc_csv

    subroutine read_numeric_csv(filename, row_index, col_names, values, has_header, has_index_col, &
                                index_min, index_max, max_col, max_rows)
        ! Read a numeric CSV, optionally using the first column as an integer row index.
        character(len=*),              intent(in)  :: filename
        integer,          allocatable, intent(out) :: row_index(:)
        character(len=32), allocatable, intent(out) :: col_names(:)
        real(dp),         allocatable, intent(out) :: values(:,:)
        logical, optional, intent(in) :: has_header, has_index_col
        integer, optional, intent(in) :: index_min, index_max, max_col, max_rows

        logical :: header, index_col, keep
        integer :: io, unit, irow, iout, j, nrows, ncols, ntokens, first_data_col, idx
        integer :: min_idx, max_idx, row_limit
        character(len=4096) :: line
        character(:), allocatable :: tokens(:)

        header = .true.
        index_col = .false.
        if (present(has_header)) header = has_header
        if (present(has_index_col)) index_col = has_index_col

        min_idx = -huge(min_idx)
        max_idx = huge(max_idx)
        if (present(index_min)) min_idx = index_min
        if (present(index_max)) max_idx = index_max
        row_limit = huge(row_limit)
        if (present(max_rows)) then
            if (max_rows < 1) then
                print '(A)', "csv_mod: max_rows must be positive"
                error stop
            end if
            row_limit = max_rows
        end if

        open(newunit=unit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print '(A,A)', "csv_mod: cannot open ", trim(filename)
            error stop
        end if

        read(unit, '(A)', iostat=io) line
        if (io /= 0 .or. trim(line) == "") then
            print '(A)', "csv_mod: empty numeric CSV"
            error stop
        end if
        call split_string(line, ",", tokens)
        ntokens = size(tokens)
        first_data_col = 1
        if (index_col) first_data_col = 2
        ncols = ntokens - first_data_col + 1
        if (present(max_col)) then
            if (max_col < 1) then
                print '(A)', "csv_mod: max_col must be positive"
                error stop
            end if
            ncols = min(ncols, max_col)
        end if
        if (ncols < 1) then
            print '(A)', "csv_mod: no numeric columns found"
            error stop
        end if

        allocate(col_names(ncols))
        if (header) then
            do j = 1, ncols
                col_names(j) = adjustl(tokens(first_data_col + j - 1))
            end do
        else
            do j = 1, ncols
                write(col_names(j), '(A,I0)') "x", j
            end do
            rewind(unit)
        end if

        nrows = 0
        irow = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            irow = irow + 1
            call split_string(line, ",", tokens)
            if (size(tokens) < first_data_col + ncols - 1) cycle
            if (index_col) then
                read(tokens(1), *) idx
            else
                idx = irow
            end if
            keep = idx >= min_idx .and. idx <= max_idx
            if (keep) then
                nrows = nrows + 1
                if (nrows >= row_limit) exit
            end if
        end do
        if (nrows == 0) then
            print '(A)', "csv_mod: no numeric rows selected"
            error stop
        end if

        allocate(row_index(nrows), values(nrows, ncols))
        rewind(unit)
        if (header) read(unit, '(A)')

        irow = 0
        iout = 0
        do
            read(unit, '(A)', iostat=io) line
            if (io /= 0 .or. trim(line) == "") exit
            irow = irow + 1
            call split_string(line, ",", tokens)
            if (size(tokens) < first_data_col + ncols - 1) cycle
            if (index_col) then
                read(tokens(1), *) idx
            else
                idx = irow
            end if
            keep = idx >= min_idx .and. idx <= max_idx
            if (.not. keep) cycle
            iout = iout + 1
            row_index(iout) = idx
            do j = 1, ncols
                read(tokens(first_data_col + j - 1), *) values(iout, j)
            end do
            if (iout >= nrows) exit
        end do

        close(unit)
    end subroutine read_numeric_csv

    ! ------------------------------------------------------------------
    ! Compute log returns: r(i) = log(prices(i+1) / prices(i)).
    ! Returns an array of length n-1.
    ! ------------------------------------------------------------------
    function log_returns(prices) result(r)
        real(dp), intent(in) :: prices(:)
        real(dp), allocatable :: r(:)
        integer :: n
        n = size(prices)
        if (n < 2) then
            allocate(r(0))
            return
        end if
        r = log(prices(2:n) / prices(1:n-1))
    end function log_returns

    subroutine print_price_sample_info(filename, dates, n_assets, n_returns, return_desc)
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: dates(:)
        integer,          intent(in) :: n_assets
        integer, optional, intent(in) :: n_returns
        character(len=*), optional, intent(in) :: return_desc

        integer :: nret, first_return_idx, last_return_idx
        character(len=64) :: desc

        nret = size(dates) - 1
        if (present(n_returns)) nret = n_returns

        desc = "demeaned log returns"
        if (present(return_desc)) desc = return_desc

        first_return_idx = size(dates) - nret + 1
        last_return_idx  = size(dates)

        write(*,'(A,A)') "Prices file: ", trim(filename)
        write(*,'(A,I0,1X,A,A,I0,A,A,A,A)') "Using ", nret, trim(desc), " for ", n_assets, &
            " assets from ", date_label(dates(first_return_idx)), " to ", date_label(dates(last_return_idx))
    end subroutine print_price_sample_info

    integer function find_col(names, target)
        character(len=*), intent(in) :: names(:)
        character(len=*), intent(in) :: target

        find_col = findloc(uppercase(names), uppercase(target), dim=1)
    end function find_col

    integer function find_ohlc_col(symbol_tokens, field_tokens, symbol, field)
        character(len=*), intent(in) :: symbol_tokens(:), field_tokens(:)
        character(len=*), intent(in) :: symbol, field
        integer :: i

        find_ohlc_col = 0
        do i = 1, size(symbol_tokens)
            if (trim(uppercase(symbol_tokens(i))) == trim(uppercase(symbol)) .and. &
                trim(uppercase(field_tokens(i))) == trim(uppercase(field))) then
                find_ohlc_col = i
                return
            end if
        end do
    end function find_ohlc_col

    integer function count_unique_symbols(symbols)
        character(len=*), intent(in) :: symbols(:)
        integer :: i, j
        logical :: seen

        count_unique_symbols = 0
        do i = 1, size(symbols)
            if (len_trim(symbols(i)) == 0) cycle
            seen = .false.
            do j = 1, i - 1
                if (trim(uppercase(symbols(j))) == trim(uppercase(symbols(i)))) seen = .true.
            end do
            if (.not. seen) count_unique_symbols = count_unique_symbols + 1
        end do
    end function count_unique_symbols

    subroutine unique_symbols(symbols, out)
        character(len=*), intent(in) :: symbols(:)
        character(len=*), intent(out) :: out(:)
        integer :: i, j, n
        logical :: seen

        n = 0
        do i = 1, size(symbols)
            if (len_trim(symbols(i)) == 0) cycle
            seen = .false.
            do j = 1, n
                if (trim(uppercase(out(j))) == trim(uppercase(symbols(i)))) seen = .true.
            end do
            if (.not. seen) then
                n = n + 1
                out(n) = adjustl(symbols(i))
            end if
        end do
    end subroutine unique_symbols

    function date_label(yyyymmdd) result(label)
        integer, intent(in) :: yyyymmdd
        character(len=10) :: label

        write(label,'(I4.4,A,I2.2,A,I2.2)') yyyymmdd / 10000, "-", &
            mod(yyyymmdd / 100, 100), "-", mod(yyyymmdd, 100)
    end function date_label

end module csv_mod

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
    implicit none
    private
    public :: read_price_csv, log_returns

contains

    ! ------------------------------------------------------------------
    ! Copied from util.f90 in github.com/Beliavsky/DataFrame (public domain).
    ! Splits str at each occurrence of the single-character delimiter delim
    ! and returns the pieces in the allocatable array tokens.
    ! ------------------------------------------------------------------
    subroutine split_string(str, delim, tokens)
        character(len=*), intent(in)           :: str
        character(len=*), intent(in)           :: delim
        character(:), allocatable, intent(out) :: tokens(:)
        integer :: start, pos, i, count, n

        n = len_trim(str)
        if (n == 0) then
            allocate(character(len=0) :: tokens(1))
            tokens(1) = ""
            return
        end if

        ! First pass: count tokens.
        count = 0
        start = 1
        do
            pos = index(str(start:), delim)
            if (pos == 0) then
                count = count + 1
                exit
            else
                count = count + 1
                start = start + pos
            end if
        end do

        ! Allocate tokens; each token gets the full length of the input.
        allocate(character(len=n) :: tokens(count))

        ! Second pass: extract tokens.
        start = 1
        i = 1
        do
            pos = index(str(start:), delim)
            if (pos == 0) then
                tokens(i) = adjustl(str(start:))
                exit
            else
                tokens(i) = adjustl(str(start:start+pos-2))
                start = start + pos
                i = i + 1
            end if
        end do
    end subroutine split_string

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
    subroutine read_price_csv(filename, dates, col_names, values)
        character(len=*),              intent(in)  :: filename
        integer,          allocatable, intent(out) :: dates(:)
        character(len=32), allocatable, intent(out) :: col_names(:)
        real(dp),         allocatable, intent(out) :: values(:,:)

        integer :: io, unit, i, j, nrows, ncols
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
        ncols = size(tokens) - 1
        if (ncols < 1) then
            print '(A)', "csv_mod: no data columns in header"
            error stop
        end if
        allocate(col_names(ncols))
        do i = 1, ncols
            col_names(i) = adjustl(tokens(i+1))
        end do

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
                read(tokens(j+1), *) values(i, j)
            end do
        end do

        close(unit)
    end subroutine read_price_csv

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

end module csv_mod

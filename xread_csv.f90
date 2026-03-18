program xread_csv
    ! Test program for csv_mod: read vix_spy.csv, report dimensions and
    ! first/last few rows, then compute and print log-return statistics.
    use kind_mod,  only: dp
    use csv_mod,   only: read_price_csv, log_returns
    implicit none

    integer,           allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp),          allocatable :: values(:,:), r(:)
    integer :: nrows, ncols, i, j

    call read_price_csv("vix_spy.csv", dates, col_names, values)

    nrows = size(dates)
    ncols = size(col_names)

    write(*, '(A,I0,A,I0,A)') "Read ", nrows, " rows x ", ncols, " columns"
    write(*, '(A)', advance='no') "Columns:"
    do i = 1, ncols
        write(*, '(2X,A)', advance='no') trim(col_names(i))
    end do
    write(*,*)

    ! Print first 3 and last 3 rows.
    write(*, '(/,A)') "First 3 rows:"
    write(*, '(A10,*(2X,A10))') "Date", (trim(col_names(j)), j=1,ncols)
    do i = 1, min(3, nrows)
        write(*, '(I8,*(2X,F10.4))') dates(i), (values(i,j), j=1,ncols)
    end do

    write(*, '(/,A)') "Last 3 rows:"
    write(*, '(A10,*(2X,A10))') "Date", (trim(col_names(j)), j=1,ncols)
    do i = max(1, nrows-2), nrows
        write(*, '(I8,*(2X,F10.4))') dates(i), (values(i,j), j=1,ncols)
    end do

    ! Log-return statistics for each column.
    write(*, '(/,A)') "Log-return statistics (mean, std dev, min, max):"
    write(*, '(A12,4(2X,A12))') "Series", "mean", "std", "min", "max"
    do j = 1, ncols
        r = log_returns(values(:, j))
        write(*, '(A12,4(2X,F12.6))') trim(col_names(j)), &
            sum(r)/size(r), &
            sqrt(sum((r - sum(r)/size(r))**2) / (size(r)-1)), &
            minval(r), maxval(r)
    end do

end program xread_csv

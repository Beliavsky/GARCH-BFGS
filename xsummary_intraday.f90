! Print summary statistics for intraday OHLCV files (CSV or binary).
program xsummary_intraday
    use date_mod, only: print_program_header
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file
    use intraday_summary_mod, only: file_summary_t, scale_ret, &
        compute_intraday_summary, print_intraday_summary, print_summary_table, &
        print_time_gap_table, print_time_gaps
    use glob_mod, only: glob, MAX_PATH_LEN
    implicit none

    character(len=MAX_PATH_LEN), allocatable :: input_files(:)
    type(file_summary_t), allocatable        :: summaries(:)
    type(ohlcv_series_t) :: series
    real(dp) :: t0, t1
    integer :: ifile, nfiles

    call print_program_header("xsummary_intraday.f90")
    call parse_arguments(input_files)
    print '(A,F6.1)', "Return scale:    ", scale_ret
    nfiles = size(input_files)
    allocate(summaries(nfiles))
    do ifile = 1, nfiles
        if (ifile > 1) print '(A)', ""
        call cpu_time(t0)
        call read_intraday_prices_file(trim(input_files(ifile)), series)
        call cpu_time(t1)
        call compute_intraday_summary(trim(input_files(ifile)), series, summaries(ifile))
        summaries(ifile)%read_sec = t1 - t0
        call print_intraday_summary(summaries(ifile))
        if (print_time_gaps) call print_time_gap_table(series%timestamp)
    end do
    deallocate(input_files)

    if (nfiles > 1) then
        print '(A)', ""
        call print_summary_table(summaries)
    end if
    deallocate(summaries)

contains

    subroutine parse_arguments(files)
        character(len=MAX_PATH_LEN), allocatable, intent(out) :: files(:)
        character(len=MAX_PATH_LEN) :: arg
        character(len=MAX_PATH_LEN), allocatable :: glob_matches(:), expanded(:)
        integer :: nargs, iarg

        nargs = command_argument_count()
        if (nargs == 0) then
            files = [character(len=MAX_PATH_LEN) :: "c:\python\databento\spy_1s_databento.csv"]
            return
        end if
        allocate(expanded(0))
        do iarg = 1, nargs
            call get_command_argument(iarg, arg)
            if (scan(arg, "*?") > 0) then
                call glob(trim(arg), glob_matches)
                if (size(glob_matches) > 0) then
                    expanded = [character(len=MAX_PATH_LEN) :: expanded, glob_matches]
                else
                    expanded = [character(len=MAX_PATH_LEN) :: expanded, arg]
                end if
            else
                expanded = [character(len=MAX_PATH_LEN) :: expanded, arg]
            end if
        end do
        if (size(expanded) < 1) error stop "xsummary_intraday: no files matched"
        files = expanded
    end subroutine parse_arguments

end program xsummary_intraday

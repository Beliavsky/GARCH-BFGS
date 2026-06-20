! Intraday OHLCV file summary: computation and formatted output.

module intraday_summary_mod
    use kind_mod, only: dp
    use date_mod, only: date_time_t, operator(-)
    use market_data_mod, only: ohlcv_series_t, infer_bar_interval_seconds, &
        peek_ohlcv_stream_version, ohlcv_tick_stream_version, ohlcv_dp_stream_version
    use path_utils_mod, only: has_extension, basename
    use stats_mod, only: sd, mean
    use glob_mod, only: MAX_PATH_LEN
    implicit none
    private

    integer, parameter  :: trading_days_per_year = 252
    real(dp), parameter :: scale_ret = 100.0_dp
    integer, parameter  :: k_extremes = 5
    logical, parameter  :: print_time_gaps = .true.
    integer, parameter  :: max_gap_rows = 10

    type, public :: bar_t
        character(len=19) :: dt     = ""
        real(dp)          :: open   = 0.0_dp
        real(dp)          :: high   = 0.0_dp
        real(dp)          :: low    = 0.0_dp
        real(dp)          :: close  = 0.0_dp
        real(dp)          :: volume = 0.0_dp
    end type bar_t

    type, public :: ret_bar_t
        character(len=19) :: dt0     = ""
        character(len=19) :: dt1     = ""
        integer           :: gap_sec = 0
        real(dp)          :: ret     = 0.0_dp
        real(dp)          :: volume  = 0.0_dp
    end type ret_bar_t

    type, public :: file_summary_t
        character(len=MAX_PATH_LEN) :: filename    = ""
        character(len=20)           :: format_str  = ""
        integer                     :: n_obs       = 0
        integer                     :: n_days      = 0
        integer                     :: bar_sec     = 0
        real(dp)                    :: close_min   = 0.0_dp
        real(dp)                    :: close_max   = 0.0_dp
        real(dp)                    :: close_ratio = 0.0_dp
        real(dp)                    :: ann_ret     = 0.0_dp
        real(dp)                    :: ann_vol     = 0.0_dp
        real(dp)                    :: min_ret     = 0.0_dp
        real(dp)                    :: max_ret     = 0.0_dp
        character(len=10)           :: date_first  = ""
        character(len=10)           :: date_last   = ""
        real(dp)                    :: read_sec    = 0.0_dp
        type(bar_t)                 :: first_bar
        type(bar_t)                 :: last_bar
        integer                     :: n_extremes  = 0
        type(ret_bar_t)             :: top_ret(k_extremes)
        type(ret_bar_t)             :: bot_ret(k_extremes)
    end type file_summary_t

    public :: scale_ret, k_extremes, print_time_gaps, max_gap_rows
    public :: compute_intraday_summary
    public :: print_intraday_summary
    public :: print_summary_table
    public :: print_time_gap_table

contains

    subroutine compute_intraday_summary(filename, series, summary)
        character(len=*), intent(in)      :: filename
        type(ohlcv_series_t), intent(in)  :: series
        type(file_summary_t), intent(out) :: summary
        real(dp), allocatable :: returns(:), tmp_vals(:)
        integer, allocatable  :: bar_idx(:), ord(:)
        real(dp) :: bars_per_day, rtmp
        integer  :: n_obs, n_days, n_ret, ver, i, n_k, j, jj, itmp

        n_obs = series%nobs()
        n_days = 1
        do i = 2, n_obs
            associate(d => series%timestamp(i)%date, d1 => series%timestamp(i-1)%date)
                if (d%year /= d1%year .or. d%month /= d1%month .or. d%day /= d1%day) &
                    n_days = n_days + 1
            end associate
        end do

        allocate(returns(max(n_obs - 1, 1)), bar_idx(max(n_obs - 1, 1)))
        n_ret = 0
        do i = 2, n_obs
            associate(d => series%timestamp(i)%date, d1 => series%timestamp(i-1)%date)
                if (d%year == d1%year .and. d%month == d1%month .and. d%day == d1%day) then
                    n_ret = n_ret + 1
                    returns(n_ret) = scale_ret * log(series%close(i) / series%close(i-1))
                    bar_idx(n_ret) = i
                end if
            end associate
        end do
        bars_per_day = real(n_ret, dp) / real(n_days, dp)
        if (n_ret > 1) then
            summary%ann_ret = mean(returns(1:n_ret)) * bars_per_day * real(trading_days_per_year, dp)
            summary%ann_vol = sd(returns(1:n_ret)) * sqrt(bars_per_day * real(trading_days_per_year, dp))
            summary%min_ret = minval(returns(1:n_ret))
            summary%max_ret = maxval(returns(1:n_ret))
        else
            summary%ann_ret = 0.0_dp
            summary%ann_vol = 0.0_dp
            summary%min_ret = 0.0_dp
            summary%max_ret = 0.0_dp
        end if

        n_k = min(k_extremes, n_ret)
        summary%n_extremes = n_k
        if (n_k >= 1) then
            allocate(ord(n_ret), tmp_vals(n_ret))
            do j = 1, n_ret
                ord(j) = j
                tmp_vals(j) = returns(j)
            end do
            do j = 1, n_k
                jj = j
                do i = j + 1, n_ret
                    if (tmp_vals(i) < tmp_vals(jj)) jj = i
                end do
                rtmp = tmp_vals(j); tmp_vals(j) = tmp_vals(jj); tmp_vals(jj) = rtmp
                itmp = ord(j);      ord(j)      = ord(jj);      ord(jj)      = itmp
            end do
            do j = 1, n_k
                call fill_ret_bar(series, bar_idx(ord(j)), tmp_vals(j), summary%bot_ret(j))
            end do
            do j = 1, n_ret
                ord(j) = j
                tmp_vals(j) = returns(j)
            end do
            do j = 1, n_k
                jj = j
                do i = j + 1, n_ret
                    if (tmp_vals(i) > tmp_vals(jj)) jj = i
                end do
                rtmp = tmp_vals(j); tmp_vals(j) = tmp_vals(jj); tmp_vals(jj) = rtmp
                itmp = ord(j);      ord(j)      = ord(jj);      ord(jj)      = itmp
            end do
            do j = 1, n_k
                call fill_ret_bar(series, bar_idx(ord(j)), tmp_vals(j), summary%top_ret(j))
            end do
            deallocate(ord, tmp_vals)
        end if
        deallocate(bar_idx, returns)

        if (has_extension(filename, ".bin")) then
            ver = peek_ohlcv_stream_version(filename)
            if (ver == ohlcv_dp_stream_version) then
                summary%format_str = "binary dp (v2)"
            else if (ver == ohlcv_tick_stream_version) then
                summary%format_str = "binary tick (v1)"
            else
                write(summary%format_str, '(A,I0,A)') "binary (v", ver, ")"
            end if
        else
            summary%format_str = "csv"
        end if

        summary%filename    = filename
        summary%n_obs       = n_obs
        summary%n_days      = n_days
        summary%bar_sec     = infer_bar_interval_seconds(series)
        summary%close_min   = minval(series%close)
        summary%close_max   = maxval(series%close)
        summary%close_ratio = summary%close_max / summary%close_min
        write(summary%date_first, '(I4.4,"-",I2.2,"-",I2.2)') &
            series%timestamp(1)%date%year, series%timestamp(1)%date%month, &
            series%timestamp(1)%date%day
        write(summary%date_last, '(I4.4,"-",I2.2,"-",I2.2)') &
            series%timestamp(n_obs)%date%year, series%timestamp(n_obs)%date%month, &
            series%timestamp(n_obs)%date%day
        call fill_bar(series, 1, summary%first_bar)
        call fill_bar(series, n_obs, summary%last_bar)
    end subroutine compute_intraday_summary

    subroutine fill_bar(series, idx, bar)
        type(ohlcv_series_t), intent(in) :: series
        integer, intent(in)              :: idx
        type(bar_t), intent(out)         :: bar

        write(bar%dt, '(I4.4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
            series%timestamp(idx)%date%year, series%timestamp(idx)%date%month, &
            series%timestamp(idx)%date%day,  series%timestamp(idx)%hour, &
            series%timestamp(idx)%minute,    series%timestamp(idx)%second
        bar%open   = series%open(idx)
        bar%high   = series%high(idx)
        bar%low    = series%low(idx)
        bar%close  = series%close(idx)
        bar%volume = series%volume(idx)
    end subroutine fill_bar

    subroutine fill_ret_bar(series, idx, ret_val, rb)
        type(ohlcv_series_t), intent(in) :: series
        integer, intent(in)              :: idx
        real(dp), intent(in)             :: ret_val
        type(ret_bar_t), intent(out)     :: rb

        write(rb%dt0, '(I4.4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
            series%timestamp(idx-1)%date%year, series%timestamp(idx-1)%date%month, &
            series%timestamp(idx-1)%date%day,  series%timestamp(idx-1)%hour, &
            series%timestamp(idx-1)%minute,    series%timestamp(idx-1)%second
        write(rb%dt1, '(I4.4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2,":",I2.2)') &
            series%timestamp(idx)%date%year, series%timestamp(idx)%date%month, &
            series%timestamp(idx)%date%day,  series%timestamp(idx)%hour, &
            series%timestamp(idx)%minute,    series%timestamp(idx)%second
        rb%gap_sec = (series%timestamp(idx)%hour   - series%timestamp(idx-1)%hour)   * 3600 + &
                    (series%timestamp(idx)%minute - series%timestamp(idx-1)%minute) * 60   + &
                    (series%timestamp(idx)%second - series%timestamp(idx-1)%second)
        rb%ret    = ret_val
        rb%volume = series%volume(idx)
    end subroutine fill_ret_bar

    subroutine print_intraday_summary(s)
        type(file_summary_t), intent(in) :: s
        integer :: j

        print '(A,A)',  "File:            ", trim(s%filename)
        print '(A,A)',  "Format:          ", trim(s%format_str)
        print '(A,I0,A)', "Bar interval:    ", s%bar_sec, " s"
        print '(A,I0)', "Observations:    ", s%n_obs
        print '(A,I0)', "Days:            ", s%n_days
        print '(A,F10.1)', "Obs/day:         ", real(s%n_obs, dp) / real(s%n_days, dp)
        print *
        print '(A)', "                 Datetime      Open      High       Low     Close       Volume"
        print '("First ",A19,4(1X,F9.4),1X,I12)', s%first_bar%dt, &
            s%first_bar%open, s%first_bar%high, s%first_bar%low, s%first_bar%close, &
            nint(s%first_bar%volume)
        print '(" Last ",A19,4(1X,F9.4),1X,I12)', s%last_bar%dt, &
            s%last_bar%open, s%last_bar%high, s%last_bar%low, s%last_bar%close, &
            nint(s%last_bar%volume)
        print *
        print '(A,F9.4,A,F9.4,A,F7.3)', "Close range:     ", s%close_min, &
            " to ", s%close_max, "  ratio:", s%close_ratio
        print '(A,F8.4)', "Ann ret (intra): ", s%ann_ret
        print '(A,F8.4)', "Ann vol (intra): ", s%ann_vol
        print '(A,F8.4,A,F8.4)', "Ret range:       ", s%min_ret, " to ", s%max_ret
        if (k_extremes >= 1 .and. s%n_extremes >= 1) then
            print *
            print '(A)', "    From                To                 Gap(s)      Return      Volume"
            print '(A,I0,A)', "Top ", s%n_extremes, " returns:"
            do j = 1, s%n_extremes
                print '(4X,A19,1X,A19,I6,F12.4,I12)', s%top_ret(j)%dt0, s%top_ret(j)%dt1, &
                    s%top_ret(j)%gap_sec, s%top_ret(j)%ret, nint(s%top_ret(j)%volume)
            end do
            print '(A,I0,A)', "Bottom ", s%n_extremes, " returns:"
            do j = 1, s%n_extremes
                print '(4X,A19,1X,A19,I6,F12.4,I12)', s%bot_ret(j)%dt0, s%bot_ret(j)%dt1, &
                    s%bot_ret(j)%gap_sec, s%bot_ret(j)%ret, nint(s%bot_ret(j)%volume)
            end do
        end if
        print '(A,F10.3)', "Read seconds:    ", s%read_sec
    end subroutine print_intraday_summary

    ! Print a frequency table of consecutive-observation time gaps.
    ! Only timestamps are used; no price/volume data needed.
    subroutine print_time_gap_table(timestamps)
        type(date_time_t), intent(in) :: timestamps(:)
        integer, allocatable :: gaps(:), ugap(:), ucnt(:)
        real(dp) :: pct, cum_pct
        integer  :: n, n_gaps, n_uniq, i, j, max_gap_val, max_gap_cnt

        n = size(timestamps)
        if (n < 2) return
        n_gaps = n - 1
        allocate(gaps(n_gaps))
        do i = 2, n
            gaps(i-1) = (timestamps(i)%date - timestamps(i-1)%date) * 86400 + &
                        timestamps(i)%seconds_since_midnight() - &
                        timestamps(i-1)%seconds_since_midnight()
        end do
        call qsort_int(gaps, 1, n_gaps)

        ! Build unique (gap_val, count) pairs from sorted gaps array.
        allocate(ugap(n_gaps), ucnt(n_gaps))
        n_uniq = 0
        i = 1
        do while (i <= n_gaps)
            n_uniq = n_uniq + 1
            j = i
            do while (j < n_gaps)
                if (gaps(j+1) /= gaps(i)) exit
                j = j + 1
            end do
            ugap(n_uniq) = gaps(i)
            ucnt(n_uniq) = j - i + 1
            i = j + 1
        end do
        deallocate(gaps)

        max_gap_val = ugap(n_uniq)   ! last in ascending-sorted list
        max_gap_cnt = ucnt(n_uniq)

        ! Sort by count descending; carry gap values along.
        call qsort_int_pair_desc(ucnt, ugap, 1, n_uniq)

        print '(A)', "Gap frequency (seconds between observations):"
        print '(A8,A12,A10,A10)', "  Gap(s)", "       Count", "  Fraction", "  Cum.Freq"
        cum_pct = 0.0_dp
        do i = 1, min(n_uniq, max_gap_rows)
            pct     = 100.0_dp * real(ucnt(i), dp) / real(n_gaps, dp)
            cum_pct = cum_pct + pct
            print '(I8,I12,F9.2,A1,F9.2,A1)', ugap(i), ucnt(i), pct, "%", cum_pct, "%"
        end do
        pct = 100.0_dp * real(max_gap_cnt, dp) / real(n_gaps, dp)
        print '(A,I0,A,I0,A,F6.2,A)', &
            "Largest gap: ", max_gap_val, " s  count: ", max_gap_cnt, "  fraction: ", pct, "%"
        deallocate(ugap, ucnt)
    end subroutine print_time_gap_table

    subroutine print_summary_table(sums)
        type(file_summary_t), intent(in) :: sums(:)
        character(len=512) :: base
        integer :: i, col

        col = 8
        do i = 1, size(sums)
            col = max(col, len_trim(basename(sums(i)%filename)))
        end do

        print '(A,I0,A)', "Summary across ", size(sums), " files:"
        print '(A)', ""
        print '(A)', repeat("-", col + 105)
        write(*, '(A)', advance='no') repeat(" ", col - 8) // "Filename"
        print '(2X,A8,A7,A7,A12,A12,A9,A9,A7,A8,A8,A8,A8)', &
            "     Obs", "   Days", "Bar(s)", "   First    ", "    Last    ", &
            "   ClMin", "   ClMax", "  Ratio", " AnnRet", " AnnVol", " MinRet", " MaxRet"
        print '(A)', repeat("-", col + 105)
        do i = 1, size(sums)
            base = basename(sums(i)%filename)
            write(*, '(A)', advance='no') trim(base) // repeat(" ", col - len_trim(base) + 2)
            print '(I8,I7,I7,2X,A10,2X,A10,2(1X,F8.3),F7.3,4(F8.4))', &
                sums(i)%n_obs, sums(i)%n_days, sums(i)%bar_sec, &
                sums(i)%date_first, sums(i)%date_last, &
                sums(i)%close_min, sums(i)%close_max, &
                sums(i)%close_ratio, sums(i)%ann_ret, sums(i)%ann_vol, &
                sums(i)%min_ret, sums(i)%max_ret
        end do
        print '(A)', repeat("-", col + 105)
    end subroutine print_summary_table

    recursive subroutine qsort_int_pair_desc(key, val, lo, hi)
        integer, intent(inout) :: key(:), val(:)
        integer, intent(in)    :: lo, hi
        integer :: pivot, ii, jj, tk, tv
        if (lo >= hi) return
        pivot = key((lo + hi) / 2)
        ii = lo; jj = hi
        do while (ii <= jj)
            do while (key(ii) > pivot); ii = ii + 1; end do
            do while (key(jj) < pivot); jj = jj - 1; end do
            if (ii <= jj) then
                tk = key(ii); key(ii) = key(jj); key(jj) = tk
                tv = val(ii); val(ii) = val(jj); val(jj) = tv
                ii = ii + 1; jj = jj - 1
            end if
        end do
        if (lo < jj) call qsort_int_pair_desc(key, val, lo, jj)
        if (ii < hi) call qsort_int_pair_desc(key, val, ii, hi)
    end subroutine qsort_int_pair_desc

    recursive subroutine qsort_int(a, lo, hi)
        integer, intent(inout) :: a(:)
        integer, intent(in)    :: lo, hi
        integer :: pivot, ii, jj, tmp
        if (lo >= hi) return
        pivot = a((lo + hi) / 2)
        ii = lo; jj = hi
        do while (ii <= jj)
            do while (a(ii) < pivot); ii = ii + 1; end do
            do while (a(jj) > pivot); jj = jj - 1; end do
            if (ii <= jj) then
                tmp = a(ii); a(ii) = a(jj); a(jj) = tmp
                ii = ii + 1; jj = jj - 1
            end if
        end do
        if (lo < jj) call qsort_int(a, lo, jj)
        if (ii < hi) call qsort_int(a, ii, hi)
    end subroutine qsort_int

end module intraday_summary_mod

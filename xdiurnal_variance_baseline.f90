! Estimate a deterministic intraday diurnal variance curve from OHLCV prices.

module xdiurnal_variance_baseline_mod
    use kind_mod, only: dp
    use date_mod, only: print_program_header, yyyymmdd, seconds_per_hour, seconds_per_minute, date_label
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file, filter_intraday_session, &
                               intraday_bin_ids, infer_bar_interval_seconds, default_session_start_seconds
    use intraday_vol_baseline_mod, only: daily_variance_lag1, estimate_diurnal_variance_baseline
    implicit none
    private

    character(len=*), parameter :: default_input_file = "c:\python\intraday_prices\spy_5min_databento.csv"
    character(len=*), parameter :: default_output_file = "diurnal_variance_baseline.csv"
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2

    public :: run_diurnal_variance_baseline

contains

    ! Read intraday prices and estimate a normalized time-of-day variance multiplier.
    subroutine run_diurnal_variance_baseline()
        character(len=512) :: input_file, output_file
        type(ohlcv_series_t) :: raw_bars, regular_bars
        real(dp), allocatable :: returns(:), daily_var(:), diurnal_var(:)
        integer, allocatable :: bin_id(:), date_id(:)
        integer :: nobs, bar_seconds
        real(dp) :: t0, t1, tread, tcompute

        call print_program_header("xdiurnal_variance_baseline.f90")
        call parse_arguments(input_file, output_file)

        call cpu_time(t0)
        call read_intraday_prices_file(input_file, raw_bars)
        call filter_intraday_session(raw_bars, regular_bars)
        call cpu_time(t1)
        tread = t1 - t0

        call cpu_time(t0)
        call build_return_inputs(regular_bars, returns, bin_id, date_id)
        nobs = size(returns)
        allocate(daily_var(nobs), diurnal_var(nobs))
        call daily_variance_lag1(returns, date_id, daily_var)
        call estimate_diurnal_variance_baseline(returns, daily_var, bin_id, diurnal_var, &
                                                smooth_diurnal=smooth_diurnal_curve, &
                                                smooth_half_width=diurnal_smooth_half_width)
        call cpu_time(t1)
        tcompute = t1 - t0

        bar_seconds = infer_bar_interval_seconds(regular_bars)
        call print_summary(input_file, output_file, raw_bars%nobs(), regular_bars%nobs(), returns, date_id, &
                           bin_id, daily_var, diurnal_var, bar_seconds, tread, tcompute)
        call write_diurnal_csv(output_file, bin_id, diurnal_var, bar_seconds)

        deallocate(returns, daily_var, diurnal_var, bin_id, date_id)
    end subroutine run_diurnal_variance_baseline

    ! Parse optional input and output filenames.
    subroutine parse_arguments(input_file, output_file)
        character(len=*), intent(out) :: input_file, output_file
        integer :: nargs

        input_file = default_input_file
        output_file = default_output_file
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, input_file)
        if (nargs >= 2) call get_command_argument(2, output_file)
    end subroutine parse_arguments

    ! Convert regular-session bars to within-day log close-to-close returns and bin identifiers.
    subroutine build_return_inputs(bars, returns, bin_id, date_id)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: returns(:)
        integer, allocatable, intent(out) :: bin_id(:), date_id(:)
        integer, allocatable :: bar_bins(:), ret_bin(:), ret_date(:)
        real(dp), allocatable :: ret_all(:)
        integer :: n, i, k, this_date, prev_date, bar_seconds

        n = bars%nobs()
        if (n < 3) error stop "build_return_inputs: not enough regular-session bars"
        bar_seconds = infer_bar_interval_seconds(bars)
        call intraday_bin_ids(bars, bar_bins, interval_seconds=bar_seconds)
        allocate(ret_all(n - 1), ret_bin(n - 1), ret_date(n - 1))
        k = 0
        do i = 2, n
            this_date = yyyymmdd(bars%timestamp(i)%date)
            prev_date = yyyymmdd(bars%timestamp(i - 1)%date)
            if (this_date /= prev_date) cycle
            k = k + 1
            ret_all(k) = log(bars%close(i) / bars%close(i - 1))
            ret_bin(k) = bar_bins(i)
            ret_date(k) = this_date
        end do
        if (k < 3) error stop "build_return_inputs: not enough intraday returns"
        allocate(returns(k), bin_id(k), date_id(k))
        returns = ret_all(1:k)
        bin_id = ret_bin(1:k)
        date_id = ret_date(1:k)
        deallocate(bar_bins, ret_all, ret_bin, ret_date)
    end subroutine build_return_inputs

    ! Print input details and the estimated diurnal variance curve by bin.
    subroutine print_summary(input_file, output_file, nraw, nregular, returns, date_id, bin_id, daily_var, &
                             diurnal_var, bar_seconds, read_sec, compute_sec)
        character(len=*), intent(in) :: input_file, output_file
        integer, intent(in) :: nraw, nregular, date_id(:), bin_id(:), bar_seconds
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), read_sec, compute_sec
        integer :: nbins, b, n
        real(dp) :: rv, avg_daily_var

        nbins = maxval(bin_id)
        n = size(returns)
        rv = sum(returns**2)
        avg_daily_var = sum(daily_var) / real(n, dp)

        print '(A,A)', "Input file: ", trim(input_file)
        print '(A,A)', "Output CSV: ", trim(output_file)
        print '(A,I0)', "Rows read: ", nraw
        print '(A,I0)', "Regular-session bars: ", nregular
        print '(A,I0,A,A,A,A)', "Intraday returns: ", n, " from ", date_label(date_id(1)), &
              " to ", date_label(date_id(n))
        print '(A,I0)', "Bar seconds: ", bar_seconds
        print '(A,I0)', "Intraday bins: ", nbins
        print '(A,I0)', "Populated return bins: ", count_populated_bins(bin_id)
        print '(A,L1,A,I0)', "Smoothed: ", smooth_diurnal_curve, "  half_width: ", diurnal_smooth_half_width
        print '(A,ES12.4)', "Realized variance sum: ", rv
        print '(A,ES12.4)', "Mean lag-1 daily variance forecast: ", avg_daily_var
        print '(A,F10.3)', "Read seconds: ", read_sec
        print '(A,F10.3)', "Compute seconds: ", compute_sec
        print '(A)', ""
        print '(A)', "Estimated diurnal variance multiplier"
        print '(A)', "------------------------------------------------"
        print '(A8,1X,A8,1X,A10,1X,A14)', "bin", "time", "nobs", "diurnal_var"
        print '(A)', "------------------------------------------------"
        do b = 1, nbins
            call print_bin_row(b, bin_id, diurnal_var, bar_seconds)
        end do
        print '(A)', "------------------------------------------------"
    end subroutine print_summary

    ! Print one row of the bin-level diurnal curve.
    subroutine print_bin_row(b, bin_id, diurnal_var, bar_seconds)
        integer, intent(in) :: b, bin_id(:), bar_seconds
        real(dp), intent(in) :: diurnal_var(:)
        integer :: count_b
        real(dp) :: avg_b

        count_b = count(bin_id == b)
        if (count_b < 1) return
        if (count_b > 0) then
            avg_b = sum(diurnal_var, mask=bin_id == b) / real(count_b, dp)
        else
            avg_b = 0.0_dp
        end if
        print '(I8,1X,A8,1X,I10,1X,F14.6)', b, time_label(default_session_start_seconds + (b - 1)*bar_seconds), &
              count_b, avg_b
    end subroutine print_bin_row

    ! Write the bin-level diurnal curve to CSV.
    subroutine write_diurnal_csv(filename, bin_id, diurnal_var, bar_seconds)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: bin_id(:), bar_seconds
        real(dp), intent(in) :: diurnal_var(:)
        integer :: unit, io, b, nbins, count_b
        real(dp) :: avg_b

        nbins = maxval(bin_id)
        open(newunit=unit, file=filename, status='replace', action='write', iostat=io)
        if (io /= 0) then
            print '(A,A)', "write_diurnal_csv: cannot open ", trim(filename)
            error stop
        end if
        write(unit, '(A)') "bin,time,nobs,diurnal_var"
        do b = 1, nbins
            count_b = count(bin_id == b)
            if (count_b < 1) cycle
            if (count_b > 0) then
                avg_b = sum(diurnal_var, mask=bin_id == b) / real(count_b, dp)
            else
                avg_b = 0.0_dp
            end if
            write(unit, '(I0,A,A,A,I0,A,ES18.10)') b, ",", time_label(default_session_start_seconds + (b - 1)*bar_seconds), &
                                                    ",", count_b, ",", avg_b
        end do
        close(unit)
    end subroutine write_diurnal_csv

    ! Count bins with at least one return observation.
    integer function count_populated_bins(bin_id)
        integer, intent(in) :: bin_id(:)
        integer, allocatable :: bin_count(:)
        integer :: i

        allocate(bin_count(maxval(bin_id)))
        bin_count = 0
        do i = 1, size(bin_id)
            bin_count(bin_id(i)) = bin_count(bin_id(i)) + 1
        end do
        count_populated_bins = count(bin_count > 0)
        deallocate(bin_count)
    end function count_populated_bins

    ! Format seconds after midnight as HH:MM.
    pure function time_label(seconds) result(label)
        integer, intent(in) :: seconds
        character(len=8) :: label
        integer :: hh, mm

        hh = modulo(seconds, 24*seconds_per_hour) / seconds_per_hour
        mm = modulo(seconds, seconds_per_hour) / seconds_per_minute
        write(label, '(I2.2,A,I2.2)') hh, ":", mm
    end function time_label

end module xdiurnal_variance_baseline_mod

program xdiurnal_variance_baseline
    use xdiurnal_variance_baseline_mod, only: run_diurnal_variance_baseline
    implicit none

    call run_diurnal_variance_baseline()
end program xdiurnal_variance_baseline

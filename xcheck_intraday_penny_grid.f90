! Check whether intraday OHLC prices lie on the integer-penny grid.

program xcheck_intraday_penny_grid
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, check_prices_on_penny_grid
    implicit none

    character(len=256) :: filename
    character(len=5), parameter :: field_name(4) = [character(len=5) :: "open", "high", "low", "close"]
    integer, parameter :: k_bad_examples = 10
    real(dp), parameter :: max_abs_change_table = 0.10_dp
    type(ohlcv_series_t) :: series
    logical :: ok(4)
    integer :: n_bad(4)
    real(dp) :: max_abs_error(4)

    filename = "c:\python\databento\spy_1s_databento.csv"
    if (command_argument_count() >= 1) call get_command_argument(1, filename)

    call read_intraday_prices_csv(trim(filename), series)
    call check_prices_on_penny_grid(series, ok, n_bad, max_abs_error)

    print '(A)', "Intraday penny-grid check"
    print '(A,A)', "Input file: ", trim(filename)
    print '(A,I0)', "Rows read: ", series%nobs()
    print '(A)', ""
    call print_grid_summary("Penny-grid check", field_name, ok, n_bad, max_abs_error, series%nobs(), "max_err_pennies")
    print '(A)', ""
    call print_bad_examples(series, k_bad_examples, 100.0_dp, "penny", "err_o", "err_h", "err_l", "err_c")
    print '(A)', ""
    call check_prices_on_grid(series, 200.0_dp, ok, n_bad, max_abs_error)
    call print_grid_summary("Half-penny-grid check", field_name, ok, n_bad, max_abs_error, series%nobs(), "max_err_halfp")
    print '(A)', ""
    call print_bad_examples(series, k_bad_examples, 200.0_dp, "half-penny", "err_o", "err_h", "err_l", "err_c")
    print '(A)', ""
    call print_close_change_frequency(series, 200.0_dp, "0.5 penny")
    print '(A)', ""
    call print_close_change_frequency(series, 100.0_dp, "1 penny")

contains

    ! Print one OHLC grid-check summary table.
    subroutine print_grid_summary(title, field_name, ok, n_bad, max_abs_error, nobs, error_label)
        character(len=*), intent(in) :: title, field_name(4), error_label
        logical, intent(in) :: ok(4)
        integer, intent(in) :: n_bad(4), nobs
        real(dp), intent(in) :: max_abs_error(4)
        real(dp) :: frac_bad
        integer :: i

        print '(A)', trim(title)
        print '(A)', "----------------------------------------------------------"
        print '(A8,1X,A6,1X,A12,1X,A10,1X,A16)', "field", "ok", "n_bad", "frac_bad", trim(error_label)
        print '(A)', "----------------------------------------------------------"
        do i = 1, 4
            frac_bad = real(n_bad(i), dp) / real(max(nobs, 1), dp)
            print '(A8,1X,L6,1X,I12,1X,F10.6,1X,ES16.6)', field_name(i), ok(i), n_bad(i), frac_bad, max_abs_error(i)
        end do
        print '(A)', "----------------------------------------------------------"
    end subroutine print_grid_summary

    ! Check whether all OHLC prices lie on a grid measured by units_per_dollar.
    subroutine check_prices_on_grid(series, units_per_dollar, ok, n_bad, max_abs_error)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), intent(in) :: units_per_dollar
        logical, intent(out) :: ok(4)
        integer, intent(out) :: n_bad(4)
        real(dp), intent(out) :: max_abs_error(4)
        real(dp) :: err(4)
        integer :: i

        n_bad = 0
        max_abs_error = 0.0_dp
        do i = 1, series%nobs()
            err = grid_errors(series%open(i), series%high(i), series%low(i), series%close(i), units_per_dollar)
            where (err > 1.0e-8_dp) n_bad = n_bad + 1
            max_abs_error = max(max_abs_error, err)
        end do
        ok = n_bad == 0
    end subroutine check_prices_on_grid

    ! Print the first k observations where at least one OHLC field is off the requested grid.
    subroutine print_bad_examples(series, k, units_per_dollar, grid_name, err1, err2, err3, err4)
        type(ohlcv_series_t), intent(in) :: series
        integer, intent(in) :: k
        real(dp), intent(in) :: units_per_dollar
        character(len=*), intent(in) :: grid_name, err1, err2, err3, err4
        real(dp) :: err(4)
        integer :: i, first_bad, nprint

        first_bad = 0
        do i = 1, series%nobs()
            err = grid_errors(series%open(i), series%high(i), series%low(i), series%close(i), units_per_dollar)
            if (any(err > 1.0e-8_dp)) then
                first_bad = i
                exit
            end if
        end do

        if (first_bad == 0) then
            print '(A,A,A)', "No OHLC observations with non-", trim(grid_name), " increments."
            return
        end if

        print '(A,I0,A,A,A)', "First ", k, " OHLC observations with non-", trim(grid_name), " increments"
        print '(A)', "---------------------------------------------------------------------------------------------------------------"
        print '(A6,1X,A19,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10,1X,A10)', &
              "row", "timestamp", "open", "high", "low", "close", err1, err2, err3, err4
        print '(A)', "---------------------------------------------------------------------------------------------------------------"
        nprint = 0
        do i = first_bad, series%nobs()
            err = grid_errors(series%open(i), series%high(i), series%low(i), series%close(i), units_per_dollar)
            if (any(err > 1.0e-8_dp)) then
                nprint = nprint + 1
                print '(I6,1X,A19,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.4,1X,F10.6,1X,F10.6,1X,F10.6,1X,F10.6)', &
                      i, series%timestamp(i)%to_str(), series%open(i), series%high(i), series%low(i), series%close(i), &
                      err(1), err(2), err(3), err(4)
                if (nprint >= k) exit
            end if
        end do
        print '(A)', "---------------------------------------------------------------------------------------------------------------"
    end subroutine print_bad_examples

    ! Print a frequency table of close-to-close price changes in grid units.
    subroutine print_close_change_frequency(series, units_per_dollar, grid_name)
        type(ohlcv_series_t), intent(in) :: series
        real(dp), intent(in) :: units_per_dollar
        character(len=*), intent(in) :: grid_name
        integer, allocatable :: count_ticks(:)
        integer :: i, n, tick, max_tick, offset, lower_count, upper_count
        real(dp) :: frac

        n = max(series%nobs() - 1, 0)
        if (n < 1) then
            print '(A,A,A)', "Close-to-close price change frequencies at ", trim(grid_name), " increments: no returns"
            return
        end if
        max_tick = nint(units_per_dollar * max_abs_change_table)
        allocate(count_ticks(-max_tick:max_tick))
        count_ticks = 0
        lower_count = 0
        upper_count = 0
        do i = 1, n
            tick = nint(units_per_dollar * (series%close(i + 1) - series%close(i)))
            if (tick < -max_tick) then
                lower_count = lower_count + 1
            else if (tick > max_tick) then
                upper_count = upper_count + 1
            else
                count_ticks(tick) = count_ticks(tick) + 1
            end if
        end do

        print '(A,A,A)', "Close-to-close price change frequencies at ", trim(grid_name), " increments"
        print '(A)', "-----------------------------------------------"
        print '(A12,1X,A12,1X,A12)', "change", "count", "frac"
        print '(A)', "-----------------------------------------------"
        frac = real(lower_count, dp) / real(n, dp)
        print '(A12,1X,I12,1X,F12.6)', "< range", lower_count, frac
        do offset = -max_tick, max_tick
            frac = real(count_ticks(offset), dp) / real(n, dp)
            print '(F12.4,1X,I12,1X,F12.6)', real(offset, dp) / units_per_dollar, count_ticks(offset), frac
        end do
        frac = real(upper_count, dp) / real(n, dp)
        print '(A12,1X,I12,1X,F12.6)', "> range", upper_count, frac
        print '(A)', "-----------------------------------------------"
        deallocate(count_ticks)
    end subroutine print_close_change_frequency

    ! Return absolute distance from a price grid in grid units.
    pure function grid_errors(open_price, high_price, low_price, close_price, units_per_dollar) result(err)
        real(dp), intent(in) :: open_price, high_price, low_price, close_price, units_per_dollar
        real(dp) :: err(4)

        err(1) = grid_error(open_price, units_per_dollar)
        err(2) = grid_error(high_price, units_per_dollar)
        err(3) = grid_error(low_price, units_per_dollar)
        err(4) = grid_error(close_price, units_per_dollar)
    end function grid_errors

    ! Return absolute distance from a price grid in grid units.
    pure real(dp) function grid_error(price, units_per_dollar) result(err)
        real(dp), intent(in) :: price, units_per_dollar

        err = abs(units_per_dollar*price - anint(units_per_dollar*price))
    end function grid_error
end program xcheck_intraday_penny_grid

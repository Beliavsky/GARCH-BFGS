! Smoke-test reader for Databento/Yahoo-style intraday OHLCV CSV files.

module xread_intraday_prices_mod
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, intraday_bin_ids, &
                               filter_intraday_session
    implicit none
    private
    public :: run_xread_intraday_prices

contains

    subroutine run_xread_intraday_prices()
        character(len=256) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        integer, allocatable :: bin_id(:)
        integer :: n, n_regular, nargs

        filename = "c:\python\intraday_prices\spy_5min_databento.csv"
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, filename)

        call read_intraday_prices_csv(trim(filename), bars)
        n = bars%nobs()
        call filter_intraday_session(bars, regular_bars)
        call intraday_bin_ids(regular_bars, bin_id)
        n_regular = regular_bars%nobs()

        print '(A,A)', "Intraday file: ", trim(filename)
        print '(A,I0)', "Rows read: ", n
        print '(A,A,1X,A)', "First/last timestamp: ", trim(bars%timestamp(1)%to_str()), trim(bars%timestamp(n)%to_str())
        print '(A,2F12.4)', "First open/close: ", bars%open(1), bars%close(1)
        print '(A,2F12.4)', "Last  open/close: ", bars%open(n), bars%close(n)
        print '(A,I0)', "Regular-session rows: ", n_regular
        print '(A,A,1X,A)', "Regular first/last timestamp: ", &
              trim(regular_bars%timestamp(1)%to_str()), trim(regular_bars%timestamp(n_regular)%to_str())
        print '(A,2I8)', "First/last regular bin_id: ", bin_id(1), bin_id(n_regular)
        print '(A,F14.0)', "Regular-session volume: ", sum(regular_bars%volume)

        deallocate(bin_id)
    end subroutine run_xread_intraday_prices

end module xread_intraday_prices_mod

program xread_intraday_prices
    use xread_intraday_prices_mod, only: run_xread_intraday_prices
    implicit none
    call run_xread_intraday_prices()
end program xread_intraday_prices

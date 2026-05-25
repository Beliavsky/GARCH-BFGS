! Test annual seasonal volatility in calendar-day return series.
!
! Usage:
!   xseasonal_vol_calendar.exe [csv_file] [nharm] [max_col] [PRICE|RETURNS]
!
! The CSV is expected to have a YYYY-MM-DD date column followed by one or more
! numeric columns.  In PRICE mode, the program first computes log returns from
! adjacent prices.  In RETURNS mode, the columns are used directly as returns.
! The test uses annual Fourier terms with nharm harmonics in a regression of
! log((r_t - mean(r))^2 + eps) on the seasonal terms.

program xseasonal_vol_calendar
    use kind_mod, only: dp
    use csv_mod, only: read_price_csv, date_label
    use seasonal_vol_mod, only: seasonal_vol_result_t, test_annual_vol_seasonality_calendar
    implicit none

    character(len=*), parameter :: default_csv_file = "spy_efa_eem_tlt_lqd.csv"
    character(len=*), parameter :: default_input_kind = "PRICE"
    integer, parameter :: default_nharm = 4
    integer, parameter :: default_max_col = 20

    integer, allocatable :: dates(:), test_dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: values(:,:), returns(:,:)
    type(seasonal_vol_result_t) :: result
    real(dp) :: annualization
    character(len=256) :: csv_file, arg, input_kind
    integer :: nharm, max_col, ncols, nobs, icol, ios

    csv_file = default_csv_file
    input_kind = default_input_kind
    nharm = default_nharm
    max_col = default_max_col

    if (command_argument_count() >= 1) call get_command_argument(1, csv_file)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *, iostat=ios) nharm
        if (ios /= 0 .or. nharm < 1) then
            write(*,'(A)') "Usage: xseasonal_vol_calendar.exe [returns_csv] [nharm>=1] [max_col>=1]"
            error stop "xseasonal_vol_calendar: invalid nharm"
        end if
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        read(arg, *, iostat=ios) max_col
        if (ios /= 0 .or. max_col < 1) then
            write(*,'(A)') "Usage: xseasonal_vol_calendar.exe [csv_file] [nharm>=1] [max_col>=1] [PRICE|RETURNS]"
            error stop "xseasonal_vol_calendar: invalid max_col"
        end if
    end if
    if (command_argument_count() >= 4) call get_command_argument(4, input_kind)

    call read_price_csv(trim(csv_file), dates, col_names, values, max_col=max_col)
    ncols = size(values, 2)
    if (trim(input_kind) == "PRICE" .or. trim(input_kind) == "price") then
        nobs = size(values, 1) - 1
        allocate(returns(nobs,ncols), test_dates(nobs))
        returns = log(values(2:size(values,1),:) / values(1:size(values,1)-1,:))
        test_dates = dates(2:size(dates))
        annualization = 252.0_dp
    else if (trim(input_kind) == "RETURNS" .or. trim(input_kind) == "returns") then
        nobs = size(values, 1)
        allocate(returns(nobs,ncols), test_dates(nobs))
        returns = values
        test_dates = dates
        annualization = 365.0_dp
    else
        write(*,'(A)') "Usage: xseasonal_vol_calendar.exe [csv_file] [nharm>=1] [max_col>=1] [PRICE|RETURNS]"
        error stop "xseasonal_vol_calendar: input kind must be PRICE or RETURNS"
    end if

    write(*,'(A,A)') "Input file: ", trim(csv_file)
    write(*,'(A,A)') "Input kind: ", trim(input_kind)
    write(*,'(A,I0,A,I0,A,A,A,A)') "Using ", nobs, " return observations for ", ncols, &
        " series from ", date_label(test_dates(1)), " to ", date_label(test_dates(size(test_dates)))

    write(*,'(A,I0,A,I0,A,F0.0)') "Annual Fourier harmonics: ", nharm, "  columns tested: ", ncols, &
        "  vol annualization: ", annualization
    write(*,'(/,A)') "Calendar-day annual seasonal volatility test"
    write(*,'(A)') repeat("-", 190)
    write(*,'(A16,1X,A8,1X,A8,1X,A8,1X,A10,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A7,1X,A12,1X,A7,1X,A8)') &
        "Series", "nobs", "df_num", "df_den", "F_stat", "p_value", "R2", "dAIC", "dBIC", &
        "vol_avg%", "vol_sd%", "vol_min%", "min_day", "vol_max%", "max_day", "ok"
    write(*,'(A)') repeat("-", 190)
    do icol = 1, ncols
        call test_annual_vol_seasonality_calendar(test_dates, returns(:,icol), nharm, result, annualization)
        write(*,'(A16,1X,I8,1X,I8,1X,I8,1X,F10.4,1X,ES12.4,1X,F12.6,1X,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,A7,1X,F12.4,1X,A7,1X,L8)') &
            trim(col_names(icol)), result%nobs, result%df_num, result%df_den, &
            result%f_stat, result%p_value, result%r2, result%delta_aic, result%delta_bic, &
            result%seasonal_vol_avg, result%seasonal_vol_sd, result%seasonal_vol_min, &
            result%seasonal_vol_min_date, result%seasonal_vol_max, result%seasonal_vol_max_date, result%ok
    end do
    write(*,'(A)') repeat("-", 190)

end program xseasonal_vol_calendar

! Fit the Rough Fractional Stochastic Volatility (RFSV) model to intraday data.
!
! Model (Gatheral, Jaisson, Rosenbaum 2018):
!   log(RV_t) = nu + sigma * W^H(t)
! where W^H is standard fractional Brownian motion with Hurst index H.
!
! Estimation:
!   nu    = mean(log RV)
!   H, sigma from log-MSD regression:
!     log E[|log(RV_{t+k}) - log(RV_t)|^2] = log(sigma^2) + 2H * log(k)
!
! Forecasting:
!   Direct AR(ar_p) regression of log(RV_{t+h}) on log(RV_t), ..., log(RV_{t-p+1}).
!   Log-HAR (Corsi) benchmark with daily / weekly / monthly predictors.
!   OOS evaluation on the last n_oos trading days.
!
! Input: intraday OHLCV CSV (any bar frequency; daily RV computed from within-day
!        log close-to-close returns).
!
! Usage: xfit_rfsv_intraday.exe [file1 [file2 ...]]
!   Files may be glob patterns.  With no arguments file_pattern below is used.

program xfit_rfsv_intraday
    use kind_mod,                       only: dp
    use date_mod,                       only: print_program_header, date_label
    use strings_mod,                    only: int_to_str
    use input_files_mod,                only: collect_input_filenames, MAX_PATH_LEN
    use market_data_mod,                only: ohlcv_series_t, read_intraday_prices_file, &
                                              filter_intraday_session
    use intraday_realized_measures_mod, only: daily_log_rv_and_rsv_neg
    use rough_sv_mod,                   only: estimate_rfsv
    use rfsv_intraday_mod,              only: print_log_msd_table, eval_oos_one_horizon
    use path_utils_mod,                 only: asset_label
    use stats_mod,                      only: mean
    implicit none

    character(len=*), parameter :: file_pattern = &
        "c:\python\databento\data_1min\*.bin" ! "c:\python\intraday_prices\spy_5min_databento.csv"
    integer, parameter :: ar_p        = 22            ! AR order for forecasting
    integer, parameter :: n_oos       = 252           ! OOS evaluation window (1 trading year)
    integer, parameter :: msd_lags(5) = [1, 5, 10, 22, 44]  ! lags for log-MSD regression
    integer, parameter :: fc_horz(3)  = [1, 5, 22]          ! forecast horizons (days)
    logical, parameter :: use_rfsv    = .true.        ! include true RFSV fGn forecasts

    character(len=MAX_PATH_LEN), allocatable :: filenames(:)
    character(len=8), allocatable :: tickers(:)
    type(ohlcv_series_t) :: raw_series, session_series
    real(dp), allocatable :: log_rv(:), log_rsv_neg(:)
    real(dp), allocatable :: rmse_all(:,:,:), qlike_all(:,:,:)
    integer,  allocatable :: rv_dates(:)
    integer :: ndays, ntrain, ifile, ih, iacc, nfiles, i
    real(dp) :: h_hat, sigma_hat, nu_hat, implied_vol
    real(dp) :: t_start, t_end

    call print_program_header("xfit_rfsv_intraday.f90")
    call cpu_time(t_start)
    call collect_input_filenames(filenames, &
        file_pattern=file_pattern, &
        default_filenames=[character(len=MAX_PATH_LEN) :: file_pattern])
    nfiles = size(filenames)
    allocate(tickers(nfiles), rmse_all(nfiles, size(fc_horz), 4), qlike_all(nfiles, size(fc_horz), 4))
    iacc = 0

    do ifile = 1, nfiles

        ! ---- load data ----
        call read_intraday_prices_file(trim(filenames(ifile)), raw_series)
        call filter_intraday_session(raw_series, session_series)
        call daily_log_rv_and_rsv_neg(session_series, log_rv, log_rsv_neg, rv_dates, ndays)

        ntrain = ndays - n_oos
        if (ntrain < ar_p + maxval(fc_horz) + 5) then
            write(*,'(A,A)') "Skipping (too few days): ", trim(filenames(ifile))
            cycle
        end if
        iacc = iacc + 1
        tickers(iacc) = asset_label(filenames(ifile))

        write(*,'(A,A)') "File  : ", trim(filenames(ifile))
        write(*,'(A,I0,A,A,A,A,A)') "Days  : ", ndays, &
            "  [", date_label(rv_dates(1)), " - ", date_label(rv_dates(ndays)), "]"
        write(*,'(A,I0,A,I0,A)') "Train : ", ntrain, " days,  OOS: ", n_oos, " days"
        write(*,*)

        ! ---- estimate RFSV parameters on training data ----
        call estimate_rfsv(log_rv(1:ntrain), msd_lags, h_hat, sigma_hat, nu_hat)

        implied_vol = sqrt(252.0_dp * exp(nu_hat + 0.5_dp * sigma_hat**2)) * 100.0_dp
        write(*,'(A)') "RFSV model:  log(RV_t) = nu + sigma * fBM^H(t)"
        write(*,'(A)') repeat("-", 44)
        write(*,'(A,F8.4)') "  H     (roughness)     =", h_hat
        write(*,'(A,F8.4)') "  sigma (vol of log-RV) =", sigma_hat
        write(*,'(A,F8.4)') "  nu    (mean log-RV)   =", nu_hat
        write(*,'(A,F7.2,A)') "  implied annual vol    =", implied_vol, "  %"
        write(*,'(A)') repeat("-", 44)
        write(*,*)

        ! ---- log-MSD diagnostic table ----
        call print_log_msd_table(log_rv(1:ntrain), msd_lags, h_hat, sigma_hat)

        ! ---- OOS forecast comparison: AR(p) vs log-HAR vs RFSV ----
        write(*,'(A)') "OOS forecast evaluation (last " // &
            trim(int_to_str(n_oos)) // " trading days)"
        write(*,*)
        if (use_rfsv) then
            write(*,'(A6,4(2(1X,A10)))') "Horiz", &
                "RMSE_AR", "QLIKE_AR", "RMSE_HAR", "QLIKE_HAR", &
                "RMSE_HARA", "QLIKE_HARA", "RMSE_RFSV", "QLIKE_RFSV"
            write(*,'(A)') repeat("-", 94)
        else
            write(*,'(A6,3(2(1X,A10)))') "Horiz", &
                "RMSE_AR", "QLIKE_AR", "RMSE_HAR", "QLIKE_HAR", "RMSE_HARA", "QLIKE_HARA"
            write(*,'(A)') repeat("-", 72)
        end if
        do ih = 1, size(fc_horz)
            call eval_oos_one_horizon(log_rv, log_rsv_neg, ndays, ntrain, ar_p, fc_horz(ih), use_rfsv, h_hat, &
                rmse_all(iacc,ih,1), qlike_all(iacc,ih,1), &
                rmse_all(iacc,ih,2), qlike_all(iacc,ih,2), &
                rmse_all(iacc,ih,3), qlike_all(iacc,ih,3), &
                rmse_all(iacc,ih,4), qlike_all(iacc,ih,4))
        end do
        if (use_rfsv) then
            write(*,"(A)") repeat("-", 94)
        else
            write(*,"(A)") repeat("-", 72)
        end if
        write(*,"(/,A)") "RMSE and QLIKE are for log(RV).  QLIKE = exp(e) - e - 1, e = actual - fcst."
        write(*,*)

    end do

    if (iacc > 0) then
        write(*,*)
        write(*,'(A,I0,A)') "=== Cross-asset OOS summary (", iacc, " assets) ==="
        do ih = 1, size(fc_horz)
            write(*,*)
            write(*,'(A,I0)') "Horizon: ", fc_horz(ih)
            if (use_rfsv) then
                write(*,'(A8,4(2(1X,A10)))') "Asset   ", &
                    "RMSE_AR", "QLIKE_AR", "RMSE_HAR", "QLIKE_HAR", &
                    "RMSE_HARA", "QLIKE_HARA", "RMSE_RFSV", "QLIKE_RFSV"
                write(*,'(A)') repeat("-", 96)
                do i = 1, iacc
                    write(*,'(A8,4(2(1X,F10.4)))') tickers(i), &
                        rmse_all(i,ih,1), qlike_all(i,ih,1), &
                        rmse_all(i,ih,2), qlike_all(i,ih,2), &
                        rmse_all(i,ih,3), qlike_all(i,ih,3), &
                        rmse_all(i,ih,4), qlike_all(i,ih,4)
                end do
                write(*,'(A)') repeat("-", 96)
                write(*,'(A8,4(2(1X,F10.4)))') "Mean    ", &
                    mean(rmse_all(1:iacc,ih,1)),  mean(qlike_all(1:iacc,ih,1)), &
                    mean(rmse_all(1:iacc,ih,2)),  mean(qlike_all(1:iacc,ih,2)), &
                    mean(rmse_all(1:iacc,ih,3)),  mean(qlike_all(1:iacc,ih,3)), &
                    mean(rmse_all(1:iacc,ih,4)),  mean(qlike_all(1:iacc,ih,4))
            else
                write(*,'(A8,3(2(1X,A10)))') "Asset   ", &
                    "RMSE_AR", "QLIKE_AR", "RMSE_HAR", "QLIKE_HAR", "RMSE_HARA", "QLIKE_HARA"
                write(*,'(A)') repeat("-", 74)
                do i = 1, iacc
                    write(*,'(A8,3(2(1X,F10.4)))') tickers(i), &
                        rmse_all(i,ih,1), qlike_all(i,ih,1), &
                        rmse_all(i,ih,2), qlike_all(i,ih,2), &
                        rmse_all(i,ih,3), qlike_all(i,ih,3)
                end do
                write(*,'(A)') repeat("-", 74)
                write(*,'(A8,3(2(1X,F10.4)))') "Mean    ", &
                    mean(rmse_all(1:iacc,ih,1)),  mean(qlike_all(1:iacc,ih,1)), &
                    mean(rmse_all(1:iacc,ih,2)),  mean(qlike_all(1:iacc,ih,2)), &
                    mean(rmse_all(1:iacc,ih,3)),  mean(qlike_all(1:iacc,ih,3))
            end if
        end do
    end if

    deallocate(tickers, rmse_all, qlike_all, filenames)
    call cpu_time(t_end)
    write(*,'(A,F10.3)') "Overall elapsed seconds:", t_end - t_start

end program xfit_rfsv_intraday

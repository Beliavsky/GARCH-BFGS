! Compare realized-variance models (AR, HAR, HAR-A, RFSV) with user-specified GARCH models.
!
! Models:
!   AR(p)   : direct OLS regression of log(RV_{t+h}) on log(RV_t), ..., log(RV_{t-p+1})
!   HAR     : direct OLS log-HAR (daily + monthly log-RV predictors)
!   HAR-A   : HAR augmented with log(RSV^-) to capture leverage
!   RFSV    : fractional Brownian motion Yule-Walker forecast
!   GARCH   : any model listed in garch_names; fit on daily returns, forecast variance
!             converted to log-scale with a training-sample level adjustment
!
! GARCH log-variance forecast: log(E[h_{t+k} | F_t]) + adj
!   where h is the conditional variance of daily log-returns and adj corrects the level.
!
! Output: per-asset transposed table (rows = models, columns = horizons × metrics)
!         followed by a cross-asset mean summary.
!
! Usage: xcompare_rv_vol_models.exe [file1 [file2 ...]]   (glob patterns accepted)

program xcompare_rv_vol_models
    use kind_mod,                       only: dp
    use date_mod,                       only: print_program_header, date_label
    use strings_mod,                    only: int_to_str
    use input_files_mod,                only: collect_input_filenames, MAX_PATH_LEN
    use market_data_mod,                only: ohlcv_series_t, read_intraday_prices_file, &
                                              filter_intraday_session
    use intraday_realized_measures_mod, only: daily_rv_rsv_and_returns
    use rough_sv_mod,                   only: estimate_rfsv
    use rfsv_intraday_mod,              only: print_log_msd_table, eval_oos_one_horizon
    use path_utils_mod,                 only: asset_label
    use stats_mod,                      only: mean
    use garch_rv_compare_mod,           only: garch_rv_model_t, fit_garch_rv_model, &
                                              eval_garch_oos_one_horizon
    implicit none

    ! ---- user-specified GARCH models to include in the comparison ----
    character(len=10), parameter :: garch_names(*) = [character(len=10) :: "NAGARCH", "GJR", "REALGARCH"]

    character(len=*), parameter :: file_pattern = &
        "c:\python\databento\data_1min\*.bin"
    integer, parameter :: ar_p        = 22
    integer, parameter :: n_oos       = 252
    integer, parameter :: msd_lags(5) = [1, 5, 10, 22, 44]
    integer, parameter :: fc_horz(3)  = [1, 5, 22]

    integer, parameter :: n_garch     = size(garch_names)
    integer, parameter :: n_base      = 4          ! AR, HAR, HAR-A, RFSV
    integer, parameter :: n_models    = n_base + n_garch
    integer, parameter :: nhorizons   = size(fc_horz)

    character(len=10) :: model_names(n_models)
    character(len=MAX_PATH_LEN), allocatable :: filenames(:)
    character(len=8), allocatable :: tickers(:)
    type(ohlcv_series_t) :: raw_series, session_series
    real(dp), allocatable :: log_rv(:), log_rsv_neg(:), ret_cc(:)
    integer,  allocatable :: rv_dates(:)
    real(dp), allocatable :: rmse_all(:,:,:), qlike_all(:,:,:)
    type(garch_rv_model_t) :: garch_models(n_garch)
    integer :: ndays, ntrain, ifile, ih, iacc, nfiles, im, ig
    real(dp) :: h_hat, sigma_hat, nu_hat, implied_vol
    real(dp) :: t_start, t_end

    call print_program_header("xcompare_rv_vol_models.f90")
    call cpu_time(t_start)

    ! Build model name list: base models followed by GARCH models
    model_names(1) = "AR"
    model_names(2) = "HAR"
    model_names(3) = "HAR-A"
    model_names(4) = "RFSV"
    do ig = 1, n_garch
        model_names(n_base + ig) = garch_names(ig)
    end do

    call collect_input_filenames(filenames, &
        file_pattern=file_pattern, &
        default_filenames=[character(len=MAX_PATH_LEN) :: file_pattern])
    nfiles = size(filenames)
    allocate(tickers(nfiles), &
             rmse_all(nfiles, nhorizons, n_models), &
             qlike_all(nfiles, nhorizons, n_models))
    iacc = 0

    do ifile = 1, nfiles

        ! ---- load intraday data ----
        call read_intraday_prices_file(trim(filenames(ifile)), raw_series)
        call filter_intraday_session(raw_series, session_series)
        call daily_rv_rsv_and_returns(session_series, log_rv, log_rsv_neg, ret_cc, rv_dates, ndays)

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

        ! ---- RFSV parameters ----
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
        call print_log_msd_table(log_rv(1:ntrain), msd_lags, h_hat, sigma_hat)

        ! ---- fit GARCH models on training data ----
        if (n_garch > 0) then
            write(*,'(A,I0,A)') "GARCH models (fitted on ", ntrain, " training days):"
            do ig = 1, n_garch
                call fit_garch_rv_model(ret_cc, log_rv, ntrain, garch_names(ig), garch_models(ig))
                write(*,'(2X,A10,A,F6.4,A,ES9.2,A,F7.4,A,F7.4)') garch_names(ig), &
                    "  persist=", garch_models(ig)%persist, &
                    "  omega=",   garch_models(ig)%omega, &
                    "  b=",       garch_models(ig)%ols_b, &
                    "  a=",       garch_models(ig)%ols_a
            end do
            write(*,*)
        end if

        ! ---- OOS evaluation for each horizon ----
        do ih = 1, nhorizons
            call eval_oos_one_horizon(log_rv, log_rsv_neg, ndays, ntrain, ar_p, fc_horz(ih), &
                .true., h_hat, &
                rmse_all(iacc,ih,1), qlike_all(iacc,ih,1), &
                rmse_all(iacc,ih,2), qlike_all(iacc,ih,2), &
                rmse_all(iacc,ih,3), qlike_all(iacc,ih,3), &
                rmse_all(iacc,ih,4), qlike_all(iacc,ih,4), &
                print_row=.false.)
            do ig = 1, n_garch
                call eval_garch_oos_one_horizon(log_rv, ndays, ntrain, fc_horz(ih), &
                    garch_models(ig), &
                    rmse_all(iacc,ih,n_base+ig), qlike_all(iacc,ih,n_base+ig))
            end do
        end do

        ! ---- print transposed OOS table ----
        write(*,'(A,I0,A)') "OOS forecast evaluation (last " // &
            trim(int_to_str(n_oos)) // " trading days)"
        write(*,*)
        call print_rv_table(model_names, rmse_all(iacc,:,:), qlike_all(iacc,:,:), &
                            n_models, nhorizons, fc_horz)
        write(*,'(A)') "RMSE and QLIKE are for log(RV).  QLIKE = exp(e) - e - 1, e = actual - fcst."
        write(*,*)

    end do

    ! ---- cross-asset mean summary ----
    if (iacc > 0) then
        write(*,*)
        write(*,'(A,I0,A)') "=== Cross-asset OOS mean (", iacc, " assets) ==="
        write(*,*)
        ! Per-asset rows grouped by asset
        do im = 1, iacc
            write(*,'(A,A)') "Asset: ", trim(tickers(im))
            call print_rv_table(model_names, rmse_all(im,:,:), qlike_all(im,:,:), &
                                n_models, nhorizons, fc_horz)
        end do
        ! Overall mean
        block
            real(dp) :: rmse_mean(nhorizons, n_models), qlike_mean(nhorizons, n_models)
            integer  :: ih2, im2
            do ih2 = 1, nhorizons
                do im2 = 1, n_models
                    rmse_mean(ih2, im2)  = mean(rmse_all(1:iacc, ih2, im2))
                    qlike_mean(ih2, im2) = mean(qlike_all(1:iacc, ih2, im2))
                end do
            end do
            write(*,'(A,I0,A)') "Mean (", iacc, " assets):"
            call print_rv_table(model_names, rmse_mean, qlike_mean, n_models, nhorizons, fc_horz)
        end block
    end if

    deallocate(tickers, rmse_all, qlike_all, filenames)
    call cpu_time(t_end)
    write(*,'(A,F10.3)') "Overall elapsed seconds:", t_end - t_start

contains

    subroutine print_rv_table(names, rmse, qlike, nm, nh, horizons)
        ! Print a transposed OOS table: rows = models, columns = RMSE/QLIKE per horizon.
        character(len=10), intent(in) :: names(:)
        real(dp),          intent(in) :: rmse(:,:)     ! (nh, nm)
        real(dp),          intent(in) :: qlike(:,:)    ! (nh, nm)
        integer,           intent(in) :: nm, nh
        integer,           intent(in) :: horizons(:)
        character(len=8) :: h1, h2, h3
        integer :: im, ih

        write(h1,'(A,I0)') "RMSE_",  horizons(1)
        write(h2,'(A,I0)') "RMSE_",  horizons(2)
        write(h3,'(A,I0)') "RMSE_",  horizons(3)
        write(*,'(A10,3(2(1X,A8)))') "Model     ", &
            adjustr(h1), "QLIKE_" // trim(int_to_str(horizons(1))), &
            adjustr(h2), "QLIKE_" // trim(int_to_str(horizons(2))), &
            adjustr(h3), "QLIKE_" // trim(int_to_str(horizons(3)))
        write(*,'(A)') repeat("-", 64)
        do im = 1, nm
            write(*,'(A10,3(2(1X,F8.4)))') names(im), &
                (rmse(ih,im), qlike(ih,im), ih=1,nh)
        end do
        write(*,'(A)') repeat("-", 64)
        write(*,*)
    end subroutine print_rv_table

end program xcompare_rv_vol_models

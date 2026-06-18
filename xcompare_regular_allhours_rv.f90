! Compare regular-session and all-hours realized volatility as predictors of
! future close-to-close volatility.
!
! The target is the Gaussian likelihood of daily close-to-close returns built
! from the all-hours intraday file.  Predictors are prior-day realized
! variance from regular hours, all hours, outside regular hours, and a
! two-predictor regular+outside model.

module compare_regular_allhours_rv_mod
    use kind_mod, only: dp
    use date_mod, only: date_label, print_program_header
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv_auto, filter_intraday_session
    use intraday_realized_measures_mod, only: daily_realized_panel_t, build_daily_realized_panel
    use realized_vol_forecast_mod, only: ewma_affine_result_t, har_variance_result_t, affine2_variance_result_t, &
                                         fit_ewma_affine_variance, fit_har_variance, fit_affine2_variance, &
                                         gaussian_variance_loglik, qlike_loss
    use stats_mod, only: demean_first
    use program_utils_mod, only: elapsed_since
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: riskmetrics_lambda = 0.94_dp
    real(dp), parameter :: gtol = 1.0e-5_dp
    integer, parameter :: default_ntest_days = 250
    integer, parameter :: max_iter = 300

    type :: rv_compare_row_t
        character(len=16) :: predictor = ""
        character(len=16) :: model = ""
        integer :: k = 0
        real(dp) :: lambda = 0.0_dp
        real(dp) :: a = 0.0_dp
        real(dp) :: b1 = 0.0_dp
        real(dp) :: b2 = 0.0_dp
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: qlike = 0.0_dp
        integer :: niter = 0
        integer :: aic_rank = 0
        integer :: bic_rank = 0
        logical :: converged = .false.
    end type rv_compare_row_t

    public :: run_compare_regular_allhours_rv

contains

    ! Fit the comparison to command-line files or a default ES/JY/TY batch.
    subroutine run_compare_regular_allhours_rv()
        character(len=MAX_PATH_LEN), allocatable :: filenames(:)
        real(dp) :: t0, elapsed_sec
        integer :: i

        call print_program_header("xcompare_regular_allhours_rv.f90")
        call collect_input_filenames(filenames, &
            default_filenames=[character(len=MAX_PATH_LEN) :: &
                "c:\python\intraday_prices\continuous\ES.csv", &
                "c:\python\intraday_prices\continuous\JY.csv", &
                "c:\python\intraday_prices\continuous\TY.csv"])
        call cpu_time(t0)
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call compare_one_file(trim(filenames(i)))
        end do
        elapsed_sec = elapsed_since(t0)
        print '(A,F10.3)', "Total elapsed seconds: ", elapsed_sec
        deallocate(filenames)
    end subroutine run_compare_regular_allhours_rv

    ! Read one intraday file, build daily RV panels, and print forecast scores.
    subroutine compare_one_file(filename)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_realized_panel_t) :: all_daily, regular_daily
        real(dp), allocatable :: ret_cc(:), regular_rv(:), all_rv(:), outside_rv(:), h(:)
        integer, allocatable :: forecast_dates(:)
        type(rv_compare_row_t), allocatable :: rows(:)
        integer :: nobs, train_nobs, test_first, irow
        real(dp) :: t0, read_sec, build_sec, fit_sec

        call cpu_time(t0)
        call read_intraday_prices_csv_auto(filename, bars)
        call filter_intraday_session(bars, regular_bars)
        read_sec = elapsed_since(t0)

        call cpu_time(t0)
        call build_daily_realized_panel(bars, all_daily)
        call build_daily_realized_panel(regular_bars, regular_daily)
        call build_aligned_targets(all_daily, regular_daily, ret_cc, regular_rv, all_rv, outside_rv, forecast_dates)
        build_sec = elapsed_since(t0)

        nobs = size(ret_cc)
        if (nobs <= default_ntest_days + 30) error stop "compare_one_file: not enough daily observations"
        train_nobs = nobs - default_ntest_days
        test_first = train_nobs + 1
        call demean_first(ret_cc, train_nobs)

        allocate(h(nobs), rows(10))
        call cpu_time(t0)
        irow = 0
        call score_ewma(rows, irow, "REGULAR_RV", ret_cc, regular_rv, train_nobs, test_first, .false., h)
        call score_ewma(rows, irow, "REGULAR_RV", ret_cc, regular_rv, train_nobs, test_first, .true., h)
        call score_har(rows, irow, "REGULAR_RV", ret_cc, regular_rv, train_nobs, test_first, h)
        call score_ewma(rows, irow, "ALL_HOURS_RV", ret_cc, all_rv, train_nobs, test_first, .false., h)
        call score_ewma(rows, irow, "ALL_HOURS_RV", ret_cc, all_rv, train_nobs, test_first, .true., h)
        call score_har(rows, irow, "ALL_HOURS_RV", ret_cc, all_rv, train_nobs, test_first, h)
        call score_ewma(rows, irow, "OUTSIDE_RV", ret_cc, outside_rv, train_nobs, test_first, .false., h)
        call score_ewma(rows, irow, "OUTSIDE_RV", ret_cc, outside_rv, train_nobs, test_first, .true., h)
        call score_har(rows, irow, "OUTSIDE_RV", ret_cc, outside_rv, train_nobs, test_first, h)
        call score_regular_outside_affine(rows, irow, ret_cc, regular_rv, outside_rv, train_nobs, test_first, h)
        fit_sec = elapsed_since(t0)

        call rank_rows(rows(1:irow))
        call print_summary(filename, bars%nobs(), regular_bars%nobs(), forecast_dates, train_nobs, &
                           rows(1:irow), read_sec, build_sec, fit_sec)
        deallocate(ret_cc, regular_rv, all_rv, outside_rv, forecast_dates, h, rows)
    end subroutine compare_one_file

    ! Align all-hours close-to-close returns with same-date regular/all-hours RV.
    subroutine build_aligned_targets(all_daily, regular_daily, ret_cc, regular_rv, all_rv, outside_rv, forecast_dates)
        type(daily_realized_panel_t), intent(in) :: all_daily, regular_daily
        real(dp), allocatable, intent(out) :: ret_cc(:), regular_rv(:), all_rv(:), outside_rv(:)
        integer, allocatable, intent(out) :: forecast_dates(:)
        integer :: i, j, k, nmax

        nmax = max(size(all_daily%date) - 1, 0)
        allocate(ret_cc(nmax), regular_rv(nmax), all_rv(nmax), outside_rv(nmax), forecast_dates(nmax))
        k = 0
        do i = 1, size(all_daily%date) - 1
            j = find_date(regular_daily%date, all_daily%date(i))
            if (j < 1) cycle
            k = k + 1
            ret_cc(k) = log(all_daily%close(i + 1) / all_daily%close(i))
            regular_rv(k) = max(regular_daily%rv(j), min_var)
            all_rv(k) = max(all_daily%rv(i), min_var)
            outside_rv(k) = max(all_rv(k) - regular_rv(k), min_var)
            forecast_dates(k) = all_daily%date(i + 1)
        end do
        ret_cc = ret_cc(1:k)
        regular_rv = regular_rv(1:k)
        all_rv = all_rv(1:k)
        outside_rv = outside_rv(1:k)
        forecast_dates = forecast_dates(1:k)
    end subroutine build_aligned_targets

    ! Fit and score an affine EWMA realized-variance forecast.
    subroutine score_ewma(rows, irow, predictor, y, x, train_nobs, test_first, fit_lambda, h)
        type(rv_compare_row_t), intent(inout) :: rows(:)
        integer, intent(inout) :: irow
        character(len=*), intent(in) :: predictor
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        logical, intent(in) :: fit_lambda
        real(dp), intent(out) :: h(:)
        type(ewma_affine_result_t) :: fit

        irow = irow + 1
        call fit_ewma_affine_variance(y, x, train_nobs, fit_lambda, riskmetrics_lambda, max_iter, gtol, fit, h)
        rows(irow)%predictor = predictor
        rows(irow)%model = merge("EWMA_FIT       ", "EWMA_0P94      ", fit_lambda)
        rows(irow)%k = merge(3, 2, fit_lambda)
        rows(irow)%lambda = fit%lambda
        rows(irow)%a = fit%a
        rows(irow)%b1 = fit%b
        rows(irow)%niter = fit%niter
        rows(irow)%converged = fit%converged
        call score_test_period(rows(irow), y, h, test_first)
    end subroutine score_ewma

    ! Fit and score a HAR realized-variance forecast.
    subroutine score_har(rows, irow, predictor, y, x, train_nobs, test_first, h)
        type(rv_compare_row_t), intent(inout) :: rows(:)
        integer, intent(inout) :: irow
        character(len=*), intent(in) :: predictor
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(har_variance_result_t) :: fit

        irow = irow + 1
        call fit_har_variance(y, x, train_nobs, fit, h)
        rows(irow)%predictor = predictor
        rows(irow)%model = "HAR"
        rows(irow)%k = 4
        rows(irow)%a = fit%coef(1)
        rows(irow)%b1 = fit%coef(2)
        rows(irow)%b2 = fit%coef(3)
        rows(irow)%niter = fit%niter
        rows(irow)%converged = fit%converged
        call score_test_period(rows(irow), y, h, test_first)
    end subroutine score_har

    ! Fit and score h_t = a + b1*regular_RV_t + b2*outside_RV_t.
    subroutine score_regular_outside_affine(rows, irow, y, regular_rv, outside_rv, train_nobs, test_first, h)
        type(rv_compare_row_t), intent(inout) :: rows(:)
        integer, intent(inout) :: irow
        real(dp), intent(in) :: y(:), regular_rv(:), outside_rv(:)
        integer, intent(in) :: train_nobs, test_first
        real(dp), intent(out) :: h(:)
        type(affine2_variance_result_t) :: fit

        irow = irow + 1
        call fit_affine2_variance(y, regular_rv, outside_rv, train_nobs, max_iter, gtol, fit, h)
        rows(irow)%predictor = "REG_PLUS_OUT"
        rows(irow)%model = "AFFINE2"
        rows(irow)%k = 3
        rows(irow)%a = fit%a
        rows(irow)%b1 = fit%b1
        rows(irow)%b2 = fit%b2
        rows(irow)%niter = fit%niter
        rows(irow)%converged = fit%converged
        call score_test_period(rows(irow), y, h, test_first)
    end subroutine score_regular_outside_affine

    ! Score a fitted variance path on the out-of-sample test period.
    subroutine score_test_period(row, y, h, test_first)
        type(rv_compare_row_t), intent(inout) :: row
        real(dp), intent(in) :: y(:), h(:)
        integer, intent(in) :: test_first
        integer :: ntest

        ntest = size(y) - test_first + 1
        row%loglik = gaussian_variance_loglik(y(test_first:), h(test_first:))
        row%aic = -2.0_dp*row%loglik + 2.0_dp*real(row%k, dp)
        row%bic = -2.0_dp*row%loglik + log(real(ntest, dp))*real(row%k, dp)
        row%qlike = qlike_loss(y(test_first:), h(test_first:))
    end subroutine score_test_period

    ! Print one file's forecast comparison.
    subroutine print_summary(filename, n_all_bars, n_regular_bars, forecast_dates, train_nobs, rows, &
                             read_sec, build_sec, fit_sec)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_all_bars, n_regular_bars, forecast_dates(:), train_nobs
        type(rv_compare_row_t), intent(in) :: rows(:)
        real(dp), intent(in) :: read_sec, build_sec, fit_sec
        integer :: i, aic_best, bic_best

        aic_best = minloc(rows%aic, dim=1)
        bic_best = minloc(rows%bic, dim=1)
        print '(A)', "Regular-hours versus all-hours realized volatility forecast comparison"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0,A,I0)', "All-hours bars: ", n_all_bars, "  regular-session bars: ", n_regular_bars
        print '(A,I0,A,I0,A,A,A,A)', "Daily observations: ", size(forecast_dates), "  train: ", train_nobs, &
              "  test dates: ", date_label(forecast_dates(train_nobs + 1)), " to ", &
              date_label(forecast_dates(size(forecast_dates)))
        print '(A)', "------------------------------------------------------------------------------------------------------------"
        print '(A16,1X,A16,1X,A3,1X,A8,3(1X,A10),1X,A12,1X,A12,1X,A12,1X,A8,1X,A4,1X,A8,1X,A8)', &
              "Predictor", "Model", "k", "lambda", "a", "b1", "b2", "logL_test", "AIC_test", "BIC_test", &
              "QLIKE", "conv", "AIC_rank", "BIC_rank"
        print '(A)', "------------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            print '(A16,1X,A16,1X,I3,1X,F8.4,3(1X,ES10.3),1X,F12.3,1X,F12.3,1X,F12.3,1X,F8.4,4X,L1,1X,I8,1X,I8)', &
                  rows(i)%predictor, rows(i)%model, rows(i)%k, rows(i)%lambda, rows(i)%a, rows(i)%b1, rows(i)%b2, &
                  rows(i)%loglik, rows(i)%aic, rows(i)%bic, rows(i)%qlike, rows(i)%converged, &
                  rows(i)%aic_rank, rows(i)%bic_rank
        end do
        print '(A)', "------------------------------------------------------------------------------------------------------------"
        print '(A,A,1X,A)', "AIC selects: ", trim(rows(aic_best)%predictor), trim(rows(aic_best)%model)
        print '(A,A,1X,A)', "BIC selects: ", trim(rows(bic_best)%predictor), trim(rows(bic_best)%model)
        print '(A,F9.3,A,F9.3,A,F9.3)', "seconds read/build/fit: ", read_sec, " ", build_sec, " ", fit_sec
    end subroutine print_summary

    ! Rank rows by ascending AIC and BIC.
    subroutine rank_rows(rows)
        type(rv_compare_row_t), intent(inout) :: rows(:)
        integer :: i, j

        do i = 1, size(rows)
            rows(i)%aic_rank = 1
            rows(i)%bic_rank = 1
            do j = 1, size(rows)
                if (rows(j)%aic < rows(i)%aic) rows(i)%aic_rank = rows(i)%aic_rank + 1
                if (rows(j)%bic < rows(i)%bic) rows(i)%bic_rank = rows(i)%bic_rank + 1
            end do
        end do
    end subroutine rank_rows

    ! Return the location of a date in a sorted date vector, or zero.
    pure integer function find_date(dates, target)
        integer, intent(in) :: dates(:), target

        find_date = findloc(dates, target, dim=1)
    end function find_date

end module compare_regular_allhours_rv_mod

program xcompare_regular_allhours_rv
    use compare_regular_allhours_rv_mod, only: run_compare_regular_allhours_rv
    implicit none

    call run_compare_regular_allhours_rv()
end program xcompare_regular_allhours_rv

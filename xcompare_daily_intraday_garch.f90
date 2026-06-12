! Compare daily OHLC NAGARCH and intraday MCS-NAGARCH forecasts on daily CC returns.
!
! The program reads one intraday OHLCV file, aggregates it to daily OHLC bars,
! fits daily NAGARCH benchmarks and an intraday MCS-NAGARCH model on the training
! sample, then scores one-day-ahead variance forecasts on test close-to-close
! returns using a Gaussian likelihood.

module compare_daily_intraday_garch_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use date_mod, only: print_program_header, yyyymmdd, date_label
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, intraday_bin_ids
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod, only: fit_nagarch, nagarch_persist
    use garch_mcsgarch_mod, only: mcsgarch_params_t, mcsgarch_fit_result_t, fit_mcsgarch_nagarch
    use stats_mod, only: mean, demean_first
    use program_utils_mod, only: read_integer_arg, elapsed_since
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: gtol = 1.0e-5_dp
    integer, parameter :: default_ntest_days = 250
    integer, parameter :: max_iter_daily = 500
    integer, parameter :: max_iter_mcs = 120
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2

    type :: daily_ohlc_t
        integer, allocatable :: date(:)
        real(dp), allocatable :: open(:)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: close(:)
        real(dp), allocatable :: volume(:)
    end type daily_ohlc_t

    type :: score_row_t
        character(len=32) :: model = ""
        integer :: k = 0
        real(dp) :: scale = 1.0_dp
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: qlike = 0.0_dp
        real(dp) :: vol_ann = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type score_row_t

    public :: run_compare_daily_intraday_garch

contains

    subroutine run_compare_daily_intraday_garch()
        character(len=256) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_ohlc_t) :: daily
        real(dp), allocatable :: ret_cc(:), ret_co(:), ret_oc(:)
        real(dp), allocatable :: h_cc(:), h_co(:), h_oc(:), h_ohlc(:)
        real(dp), allocatable :: intra_ret(:), intra_scale(:), diurnal(:), q(:), diurnal_bin(:)
        real(dp), allocatable :: train_ret(:), train_scale(:)
        real(dp), allocatable :: mcs_daily_var(:), mcs_train_var(:)
        integer, allocatable :: intra_bin(:), intra_day(:), train_bin(:), train_day(:)
        integer, allocatable :: day_index(:), day_first(:), day_last(:), bar_bins(:)
        real(dp), allocatable :: intraday_rv(:)
        type(garch_params_t) :: p_cc, p_co, p_oc
        type(mcsgarch_fit_result_t) :: mcs_fit
        type(score_row_t) :: rows(4)
        integer :: nargs, ntest_days, ndays, nobs, train_end_day, train_nobs, test_first_obs, test_nobs
        integer :: niter_cc, niter_co, niter_oc, ntrain_intra
        logical :: conv_cc, conv_co, conv_oc
        real(dp) :: f_cc, f_co, f_oc, mcs_scale
        real(dp) :: t0, t1, read_sec, fit_sec, elapsed_sec

        call cpu_time(t0)
        filename = "c:\python\intraday_prices\spy_5min_databento.csv"
        ntest_days = default_ntest_days
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, filename)
        if (nargs >= 2) call read_integer_arg(2, ntest_days)

        call print_program_header("xcompare_daily_intraday_garch.f90")
        call cpu_time(t1)
        call read_intraday_prices_csv(filename, bars)
        call filter_intraday_session(bars, regular_bars)
        read_sec = elapsed_since(t1)

        call aggregate_daily_ohlc(regular_bars, daily, day_index, day_first, day_last)
        ndays = size(daily%date)
        if (ndays < ntest_days + 30) error stop "xcompare_daily_intraday_garch: not enough trading days"
        train_end_day = ndays - ntest_days
        nobs = ndays - 1
        train_nobs = train_end_day - 1
        test_first_obs = train_nobs + 1
        test_nobs = nobs - train_nobs

        allocate(ret_cc(nobs), ret_co(nobs), ret_oc(nobs), h_cc(nobs), h_co(nobs), h_oc(nobs), h_ohlc(nobs))
        call daily_returns(daily, ret_cc, ret_co, ret_oc)
        call demean_first(ret_cc, train_nobs)
        call demean_first(ret_co, train_nobs)
        call demean_first(ret_oc, train_nobs)

        call cpu_time(t1)
        call fit_nagarch(ret_cc(1:train_nobs), max_iter_daily, gtol, f_cc, p_cc, niter_cc, conv_cc)
        call fit_nagarch(ret_co(1:train_nobs), max_iter_daily, gtol, f_co, p_co, niter_co, conv_co)
        call fit_nagarch(ret_oc(1:train_nobs), max_iter_daily, gtol, f_oc, p_oc, niter_oc, conv_oc)
        call nagarch_variance_path(ret_cc, p_cc, h_cc)
        call nagarch_variance_path(ret_co, p_co, h_co)
        call nagarch_variance_path(ret_oc, p_oc, h_oc)
        h_ohlc = h_co + h_oc

        call intraday_bin_ids(regular_bars, bar_bins)
        call intraday_arrays(regular_bars, day_index, bar_bins, intraday_rv, intra_ret, intra_scale, intra_bin, intra_day)
        call select_intraday_train(intra_ret, intra_scale, intra_bin, intra_day, train_end_day, &
                                   train_ret, train_scale, train_bin, train_day)
        ntrain_intra = size(train_ret)
        allocate(diurnal(ntrain_intra), q(ntrain_intra))
        call fit_mcsgarch_nagarch(train_ret, train_scale, train_bin, max_iter_mcs, gtol, mcs_fit, diurnal, q, &
                                  smooth_diurnal=smooth_diurnal_curve, smooth_half_width=diurnal_smooth_half_width)
        call diurnal_curve_from_fit(train_bin, diurnal, diurnal_bin)
        call mcs_filtered_daily_variance(train_scale, train_day, diurnal, q, train_end_day, mcs_train_var)
        call mcs_daily_forecast_path(regular_bars, day_first, day_last, bar_bins, intraday_rv, &
                                     train_end_day, mcs_fit%params, diurnal_bin, q(size(q)), &
                                     train_ret(size(train_ret)), train_scale(size(train_scale)), diurnal(size(diurnal)), &
                                     mcs_daily_var)
        mcs_scale = sum(ret_cc(1:train_nobs)**2) / max(sum(mcs_train_var(2:train_end_day)), min_var)
        fit_sec = elapsed_since(t1)

        call fill_score_row(rows(1), "daily CC NAGARCH", 4, 1.0_dp, ret_cc(test_first_obs:nobs), &
                            h_cc(test_first_obs:nobs), niter_cc, conv_cc)
        call fill_score_row(rows(2), "daily OHLC NAGARCH", 8, 1.0_dp, ret_cc(test_first_obs:nobs), &
                            h_ohlc(test_first_obs:nobs), niter_co + niter_oc, conv_co .and. conv_oc)
        call fill_score_row(rows(3), "MCS-NAGARCH intraday", 4 + size(diurnal_bin), 1.0_dp, &
                            ret_cc(test_first_obs:nobs), mcs_daily_var(train_end_day + 1:ndays), &
                            mcs_fit%niter, mcs_fit%converged)
        call fill_score_row(rows(4), "MCS-NAGARCH scaled CC", 5 + size(diurnal_bin), mcs_scale, &
                            ret_cc(test_first_obs:nobs), mcs_scale*mcs_daily_var(train_end_day + 1:ndays), &
                            mcs_fit%niter, mcs_fit%converged)

        elapsed_sec = elapsed_since(t0)
        call print_summary(filename, daily, train_end_day, ntest_days, ntrain_intra, mcs_fit, read_sec, fit_sec, elapsed_sec)
        call print_score_table(rows, test_nobs)

        deallocate(ret_cc, ret_co, ret_oc, h_cc, h_co, h_oc, h_ohlc)
        deallocate(intra_ret, intra_scale, intra_bin, intra_day, train_ret, train_scale, train_bin, train_day)
        deallocate(diurnal, q, diurnal_bin, mcs_daily_var, mcs_train_var, intraday_rv, day_index, day_first, day_last, bar_bins)
    end subroutine run_compare_daily_intraday_garch

    subroutine aggregate_daily_ohlc(bars, daily, day_index, day_first, day_last)
        type(ohlcv_series_t), intent(in) :: bars
        type(daily_ohlc_t), intent(out) :: daily
        integer, allocatable, intent(out) :: day_index(:), day_first(:), day_last(:)
        integer :: n, i, ndays, d, current_date

        n = bars%nobs()
        if (n < 2) error stop "aggregate_daily_ohlc: not enough bars"
        allocate(day_index(n), day_first(n), day_last(n))
        ndays = 0
        current_date = -1
        do i = 1, n
            if (yyyymmdd(bars%timestamp(i)%date) /= current_date) then
                ndays = ndays + 1
                current_date = yyyymmdd(bars%timestamp(i)%date)
                day_first(ndays) = i
                if (ndays > 1) day_last(ndays - 1) = i - 1
            end if
            day_index(i) = ndays
        end do
        day_last(ndays) = n
        day_first = day_first(1:ndays)
        day_last = day_last(1:ndays)

        allocate(daily%date(ndays), daily%open(ndays), daily%high(ndays), daily%low(ndays), &
                 daily%close(ndays), daily%volume(ndays))
        do d = 1, ndays
            daily%date(d) = yyyymmdd(bars%timestamp(day_first(d))%date)
            daily%open(d) = bars%open(day_first(d))
            daily%high(d) = maxval(bars%high(day_first(d):day_last(d)))
            daily%low(d) = minval(bars%low(day_first(d):day_last(d)))
            daily%close(d) = bars%close(day_last(d))
            daily%volume(d) = sum(bars%volume(day_first(d):day_last(d)))
        end do
    end subroutine aggregate_daily_ohlc

    subroutine daily_returns(daily, ret_cc, ret_co, ret_oc)
        type(daily_ohlc_t), intent(in) :: daily
        real(dp), intent(out) :: ret_cc(:), ret_co(:), ret_oc(:)
        integer :: i

        do i = 2, size(daily%date)
            ret_cc(i - 1) = log(daily%close(i) / daily%close(i - 1))
            ret_co(i - 1) = log(daily%open(i) / daily%close(i - 1))
            ret_oc(i - 1) = log(daily%close(i) / daily%open(i))
        end do
    end subroutine daily_returns

    subroutine intraday_arrays(bars, day_index, bar_bins, intraday_rv, returns, daily_scale, bin_id, ret_day)
        type(ohlcv_series_t), intent(in) :: bars
        integer, intent(in) :: day_index(:), bar_bins(:)
        real(dp), allocatable, intent(out) :: intraday_rv(:), returns(:), daily_scale(:)
        integer, allocatable, intent(out) :: bin_id(:), ret_day(:)
        integer :: n, ndays, i, d
        real(dp) :: fallback

        n = bars%nobs()
        ndays = maxval(day_index)
        allocate(intraday_rv(ndays), returns(n), daily_scale(n), bin_id(n), ret_day(n))
        intraday_rv = 0.0_dp
        do i = 1, n
            d = day_index(i)
            returns(i) = log(bars%close(i) / bars%open(i))
            intraday_rv(d) = intraday_rv(d) + returns(i)**2
            bin_id(i) = bar_bins(i)
            ret_day(i) = d
        end do
        fallback = max(mean(intraday_rv), min_var)
        do i = 1, n
            d = day_index(i)
            if (d == 1) then
                daily_scale(i) = fallback
            else
                daily_scale(i) = max(intraday_rv(d - 1), min_var)
            end if
        end do
    end subroutine intraday_arrays

    subroutine select_intraday_train(returns, daily_scale, bin_id, ret_day, train_end_day, &
                                     train_ret, train_scale, train_bin, train_day)
        real(dp), intent(in) :: returns(:), daily_scale(:)
        integer, intent(in) :: bin_id(:), ret_day(:), train_end_day
        real(dp), allocatable, intent(out) :: train_ret(:), train_scale(:)
        integer, allocatable, intent(out) :: train_bin(:), train_day(:)
        integer :: n, i, k

        n = count(ret_day >= 2 .and. ret_day <= train_end_day)
        if (n < 20) error stop "select_intraday_train: not enough training intraday returns"
        allocate(train_ret(n), train_scale(n), train_bin(n), train_day(n))
        k = 0
        do i = 1, size(returns)
            if (ret_day(i) < 2 .or. ret_day(i) > train_end_day) cycle
            k = k + 1
            train_ret(k) = returns(i)
            train_scale(k) = daily_scale(i)
            train_bin(k) = bin_id(i)
            train_day(k) = ret_day(i)
        end do
    end subroutine select_intraday_train

    subroutine diurnal_curve_from_fit(bin_id, diurnal, diurnal_bin)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(in) :: diurnal(:)
        real(dp), allocatable, intent(out) :: diurnal_bin(:)
        integer, allocatable :: bin_count(:)
        integer :: i, nbins

        nbins = maxval(bin_id)
        allocate(diurnal_bin(nbins), bin_count(nbins))
        diurnal_bin = 0.0_dp
        bin_count = 0
        do i = 1, size(bin_id)
            diurnal_bin(bin_id(i)) = diurnal_bin(bin_id(i)) + diurnal(i)
            bin_count(bin_id(i)) = bin_count(bin_id(i)) + 1
        end do
        do i = 1, nbins
            if (bin_count(i) > 0) then
                diurnal_bin(i) = diurnal_bin(i) / real(bin_count(i), dp)
            else
                diurnal_bin(i) = 1.0_dp
            end if
        end do
        deallocate(bin_count)
    end subroutine diurnal_curve_from_fit

    subroutine mcs_filtered_daily_variance(scale, ret_day, diurnal, q, train_end_day, daily_var)
        real(dp), intent(in) :: scale(:), diurnal(:), q(:)
        integer, intent(in) :: ret_day(:), train_end_day
        real(dp), allocatable, intent(out) :: daily_var(:)
        integer :: i

        allocate(daily_var(train_end_day))
        daily_var = 0.0_dp
        do i = 1, size(scale)
            daily_var(ret_day(i)) = daily_var(ret_day(i)) + scale(i)*diurnal(i)*q(i)
        end do
        daily_var = max(daily_var, min_var)
    end subroutine mcs_filtered_daily_variance

    subroutine mcs_daily_forecast_path(bars, day_first, day_last, bar_bins, intraday_rv, train_end_day, &
                                       params, diurnal_bin, q_last_train, ret_last_train, scale_last_train, &
                                       diurnal_last_train, daily_var)
        type(ohlcv_series_t), intent(in) :: bars
        integer, intent(in) :: day_first(:), day_last(:), bar_bins(:), train_end_day
        real(dp), intent(in) :: intraday_rv(:), diurnal_bin(:), q_last_train, ret_last_train, scale_last_train
        real(dp), intent(in) :: diurnal_last_train
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), allocatable, intent(out) :: daily_var(:)
        real(dp) :: qprev, eprev, qnow, base, persistence
        integer :: ndays, d, i, b

        ndays = size(intraday_rv)
        allocate(daily_var(ndays))
        daily_var = 0.0_dp
        qprev = max(q_last_train, min_var)
        eprev = ret_last_train / sqrt(max(scale_last_train*diurnal_last_train, min_var))
        persistence = params%alpha*(1.0_dp + params%theta**2) + params%beta

        do d = train_end_day + 1, ndays
            base = max(intraday_rv(d - 1), min_var)
            qnow = qprev
            do i = day_first(d), day_last(d)
                b = min(max(bar_bins(i), 1), size(diurnal_bin))
                if (i == day_first(d)) then
                    qnow = nagarch_q_next(params, qprev, eprev)
                else
                    qnow = max(params%omega + persistence*qnow, min_var)
                end if
                daily_var(d) = daily_var(d) + base*diurnal_bin(b)*qnow
            end do
            daily_var(d) = max(daily_var(d), min_var)
            do i = day_first(d), day_last(d)
                b = min(max(bar_bins(i), 1), size(diurnal_bin))
                qnow = nagarch_q_next(params, qprev, eprev)
                eprev = log(bars%close(i) / bars%open(i)) / sqrt(max(base*diurnal_bin(b), min_var))
                qprev = qnow
            end do
        end do
    end subroutine mcs_daily_forecast_path

    subroutine nagarch_variance_path(y, params, h)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: persist, sqrth
        integer :: i

        persist = nagarch_persist(params)
        h(1) = max(params%omega / max(1.0_dp - persist, 1.0e-8_dp), mean(y**2))
        do i = 2, size(y)
            sqrth = sqrt(max(h(i - 1), min_var))
            h(i) = max(params%omega + params%alpha*(y(i - 1) - params%theta*sqrth)**2 + params%beta*h(i - 1), min_var)
        end do
    end subroutine nagarch_variance_path

    real(dp) function nagarch_q_next(params, qprev, eprev)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(in) :: qprev, eprev

        nagarch_q_next = max(params%omega + params%alpha*(eprev - params%theta*sqrt(max(qprev, min_var)))**2 + &
                             params%beta*qprev, min_var)
    end function nagarch_q_next

    subroutine fill_score_row(row, model, k, scale, y, h, niter, converged)
        type(score_row_t), intent(out) :: row
        character(len=*), intent(in) :: model
        integer, intent(in) :: k, niter
        real(dp), intent(in) :: scale, y(:), h(:)
        logical, intent(in) :: converged
        integer :: n

        n = size(y)
        row%model = model
        row%k = k
        row%scale = scale
        row%loglik = gaussian_loglik(y, h)
        row%aic = 2.0_dp*real(k, dp) - 2.0_dp*row%loglik
        row%bic = log(real(n, dp))*real(k, dp) - 2.0_dp*row%loglik
        row%qlike = mean(log(max(h, min_var)) + y**2 / max(h, min_var))
        row%vol_ann = 100.0_dp*sqrt(252.0_dp*mean(h))
        row%niter = niter
        row%converged = converged
    end subroutine fill_score_row

    real(dp) function gaussian_loglik(y, h)
        real(dp), intent(in) :: y(:), h(:)
        integer :: i

        gaussian_loglik = 0.0_dp
        do i = 1, size(y)
            gaussian_loglik = gaussian_loglik - log_sqrt_2pi - 0.5_dp*log(max(h(i), min_var)) - &
                              0.5_dp*y(i)**2 / max(h(i), min_var)
        end do
    end function gaussian_loglik

    subroutine print_summary(filename, daily, train_end_day, ntest_days, ntrain_intra, mcs_fit, read_sec, fit_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        type(daily_ohlc_t), intent(in) :: daily
        integer, intent(in) :: train_end_day, ntest_days, ntrain_intra
        type(mcsgarch_fit_result_t), intent(in) :: mcs_fit
        real(dp), intent(in) :: read_sec, fit_sec, elapsed_sec

        print '(A)', "Daily OHLC NAGARCH vs intraday MCS-NAGARCH forecast comparison"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0,A,A,A,A)', "Daily bars: ", size(daily%date), " from ", date_label(daily%date(1)), " to ", &
              date_label(daily%date(size(daily%date)))
        print '(A,I0,A,A,A,A)', "Training days: ", train_end_day, " through ", date_label(daily%date(train_end_day)), &
              "; test starts ", date_label(daily%date(train_end_day + 1))
        print '(A,I0,A,I0)', "Test days: ", ntest_days, "  training intraday bars: ", ntrain_intra
        print '(A)', "Daily OHLC forecasts use CO and OC NAGARCH variance sums; MCS forecasts use prior-day intraday RV."
        print '(A,F8.4,1X,A,F8.4,1X,A,F8.4,1X,A,F8.4,1X,A,F8.4,1X,A,I0,1X,A,L1)', &
              "MCS-NAGARCH omega=", mcs_fit%params%omega, "alpha=", mcs_fit%params%alpha, &
              "beta=", mcs_fit%params%beta, "theta=", mcs_fit%params%theta, "persist=", mcs_fit%persist, &
              "iter=", mcs_fit%niter, "conv=", mcs_fit%converged
        print '(A,F8.3,2X,A,F8.3,2X,A,F8.3)', "read_sec=", read_sec, "fit_sec=", fit_sec, "elapsed_sec=", elapsed_sec
    end subroutine print_summary

    subroutine print_score_table(rows, ntest)
        type(score_row_t), intent(in) :: rows(:)
        integer, intent(in) :: ntest
        integer :: i, best

        best = maxloc(rows%loglik, dim=1)
        print '(A)', ""
        print '(A)', "Daily close-to-close return forecast likelihood comparison"
        print '(A,I0)', "Test observations: ", ntest
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A5,1X,A8,1X,A13,1X,A13,1X,A13,1X,A10,1X,A9,1X,A6,1X,A4)', &
              "Model", "k", "scale", "logL", "AIC", "BIC", "QLIKE", "vol_ann%", "iter", "conv"
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            print '(A32,1X,I5,1X,F8.4,1X,F13.3,1X,F13.3,1X,F13.3,1X,F10.5,1X,F9.3,1X,I6,4X,L1)', &
                  rows(i)%model, rows(i)%k, rows(i)%scale, rows(i)%loglik, rows(i)%aic, rows(i)%bic, &
                  rows(i)%qlike, rows(i)%vol_ann, rows(i)%niter, rows(i)%converged
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        print '(A,A)', "Best logL: ", trim(rows(best)%model)
    end subroutine print_score_table

end module compare_daily_intraday_garch_mod

program xcompare_daily_intraday_garch
    use compare_daily_intraday_garch_mod, only: run_compare_daily_intraday_garch
    implicit none

    call run_compare_daily_intraday_garch()
end program xcompare_daily_intraday_garch

! Compare target-frequency intraday volatility forecasts using lagged EWMA
! realized-variance proxies computed at the target frequency or higher.

module compare_intraday_ewma_freq_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use date_mod, only: date_label, yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, &
                               default_session_start_seconds
    use intraday_vol_baseline_mod, only: fit_lag1_diurnal_baseline, intraday_ewma_multiplier_from_proxy, &
                                         intraday_variance_forecast
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: intraday_lambda = 0.94_dp
    real(dp), parameter :: scale_log_search_half_width = 8.0_dp
    real(dp), parameter :: t_nu_min = 2.05_dp
    real(dp), parameter :: t_nu_max = 100.0_dp
    logical, parameter :: smooth_diurnal_curve = .true.
    logical, parameter :: compare_student_t_noise = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    integer, parameter :: target_seconds(*) = [23400, 11700, 7800, 3900, 1800, 900, 600, 300, 150, 60]
    integer, parameter :: predictor_seconds(*) = [23400, 11700, 7800, 3900, 1800, 900, 600, 300, 150, 60, 30, 10, 1]
    integer, parameter :: proxy_cc = 1
    integer, parameter :: proxy_parkinson = 2
    integer, parameter :: proxy_garman_klass = 3
    character(len=16), parameter :: proxy_names(*) = [character(len=16) :: &
        "CC_RV", "PARKINSON", "GARMAN_KLASS"]

    type :: bar_series_t
        integer, allocatable :: date_id(:)
        integer, allocatable :: bin_id(:)
        integer, allocatable :: end_seconds(:)
        real(dp), allocatable :: open(:)
        real(dp), allocatable :: high(:)
        real(dp), allocatable :: low(:)
        real(dp), allocatable :: close(:)
    end type bar_series_t

    type :: target_returns_t
        real(dp), allocatable :: returns(:)
        integer, allocatable :: date_id(:)
        integer, allocatable :: bin_id(:)
    end type target_returns_t

    type :: ewma_row_t
        character(len=32) :: model = ""
        character(len=16) :: proxy = ""
        integer :: predictor_sec = 0
        integer :: k = 0
        real(dp) :: proxy_scale = 0.0_dp
        real(dp) :: loglik = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: t_nu = 0.0_dp
        real(dp) :: t_loglik = 0.0_dp
        real(dp) :: t_aic = 0.0_dp
        real(dp) :: t_bic = 0.0_dp
        real(dp) :: proxy_mean = 0.0_dp
        real(dp) :: proxy_to_target_var = 0.0_dp
        real(dp) :: proxy_max = 0.0_dp
        real(dp) :: proxy_nonzero_pct = 0.0_dp
    end type ewma_row_t

    public :: run_compare_intraday_ewma_freq

contains

    ! Read data and compare predictor frequencies for every configured target frequency.
    subroutine run_compare_intraday_ewma_freq()
        character(len=256) :: filename
        type(ohlcv_series_t) :: raw_bars, regular_bars
        integer :: nargs, itarget
        real(dp) :: t0, t1, read_sec, elapsed_sec

        filename = "c:\python\databento\spy_1s_databento.csv"
        nargs = command_argument_count()
        if (nargs >= 1) call get_command_argument(1, filename)

        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), raw_bars)
        call filter_intraday_session(raw_bars, regular_bars)
        call cpu_time(t1)
        read_sec = t1 - t0

        print '(A)', "Intraday EWMA frequency comparison"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0)', "Regular-session source bars: ", regular_bars%nobs()
        print '(A,F8.4)', "Intraday EWMA lambda: ", intraday_lambda
        print '(A)', ""

        do itarget = 1, size(target_seconds)
            call compare_one_target(regular_bars, target_seconds(itarget))
        end do

        call cpu_time(t1)
        elapsed_sec = t1 - t0
        print '(A,F10.3)', "Data read seconds: ", read_sec
        print '(A,F10.3)', "Elapsed seconds:   ", elapsed_sec
    end subroutine run_compare_intraday_ewma_freq

    ! Compare all valid predictor frequencies for one target frequency.
    subroutine compare_one_target(source_bars, target_sec)
        type(ohlcv_series_t), intent(in) :: source_bars
        integer, intent(in) :: target_sec
        type(bar_series_t) :: target_bars, pred_bars
        type(target_returns_t) :: target
        type(ewma_row_t), allocatable :: rows(:)
        real(dp), allocatable :: daily_var(:), diurnal_var(:), q(:), h(:), proxy(:), asym_proxy_mat(:,:), &
                                 proxy_mat(:,:), combined_proxy(:), combined_weights(:), neg_mask(:)
        integer :: ipred, iproxy, nvalid, ncombo, icombo, irow, aic_best, bic_best, nret

        call aggregate_bars(source_bars, target_sec, target_bars)
        nret = count_target_returns(target_bars)
        if (nret < 3) then
            print '(A,I0,A,I0,A)', "Target interval: ", target_sec, " seconds skipped (", nret, " within-day returns)"
            return
        end if
        call build_target_returns(target_bars, target)
        allocate(daily_var(size(target%returns)), diurnal_var(size(target%returns)), q(size(target%returns)), &
                 h(size(target%returns)), proxy(size(target%returns)))

        nvalid = 1
        do ipred = 1, size(predictor_seconds)
            if (valid_predictor(target_sec, predictor_seconds(ipred))) nvalid = nvalid + 2*size(proxy_names)
        end do
        ncombo = (nvalid - 1) / 2
        allocate(rows(nvalid + merge(2, 0, ncombo > 1)))
        if (ncombo < 1) error stop "compare_one_target: no valid proxy predictors"
        allocate(proxy_mat(size(target%returns), ncombo), combined_proxy(size(target%returns)), &
                 combined_weights(2*ncombo), neg_mask(size(target%returns)))
        neg_mask = negative_return_mask(target%returns)

        call fit_lag1_diurnal_baseline(target%returns, target%date_id, target%bin_id, daily_var, diurnal_var, h, &
                                       smooth_diurnal=smooth_diurnal_curve, &
                                       smooth_half_width=diurnal_smooth_half_width)
        call fill_row(rows(1), "no intraday update", "-", 0, 0, target%returns, h)

        irow = 1
        icombo = 0
        do ipred = 1, size(predictor_seconds)
            if (.not. valid_predictor(target_sec, predictor_seconds(ipred))) cycle
            call aggregate_bars(source_bars, predictor_seconds(ipred), pred_bars)
            do iproxy = 1, size(proxy_names)
                irow = irow + 1
                call predictor_proxy_for_target(pred_bars, target_sec, target, iproxy, proxy)
                icombo = icombo + 1
                proxy_mat(:, icombo) = proxy
                call fit_scaled_ewma_row(rows(irow), trim(proxy_names(iproxy)), predictor_seconds(ipred), target, &
                                         daily_var, diurnal_var, proxy, q, h)
                irow = irow + 1
                allocate(asym_proxy_mat(size(target%returns), 2))
                asym_proxy_mat(:, 1) = proxy
                asym_proxy_mat(:, 2) = proxy * negative_return_mask(target%returns)
                call fit_asymmetric_ewma_row(rows(irow), trim(proxy_names(iproxy)), predictor_seconds(ipred), target, &
                                             daily_var, diurnal_var, asym_proxy_mat, combined_weights(1:2), &
                                             combined_proxy, q, h)
                deallocate(asym_proxy_mat)
            end do
        end do
        if (ncombo > 1) then
            irow = irow + 1
            call fit_combined_ewma_row(rows(irow), target, daily_var, diurnal_var, proxy_mat, combined_weights(1:ncombo), &
                                       combined_proxy, q, h)
            irow = irow + 1
            allocate(asym_proxy_mat(size(target%returns), 2*ncombo))
            asym_proxy_mat(:, 1:ncombo) = proxy_mat
            asym_proxy_mat(:, ncombo+1:2*ncombo) = proxy_mat * spread(neg_mask, dim=2, ncopies=ncombo)
            call fit_asymmetric_combined_ewma_row(rows(irow), target, daily_var, diurnal_var, asym_proxy_mat, &
                                                  combined_weights, combined_proxy, q, h)
            deallocate(asym_proxy_mat)
        end if

        aic_best = minloc(rows%aic, dim=1)
        bic_best = minloc(rows%bic, dim=1)
        call print_target_table(target_sec, target, rows, aic_best, bic_best)

        deallocate(proxy_mat, combined_proxy, combined_weights, neg_mask)
        deallocate(daily_var, diurnal_var, q, h, proxy, rows)
    end subroutine compare_one_target

    ! Aggregate source OHLCV bars into interval-second OHLC bars by trading date and bucket.
    subroutine aggregate_bars(source_bars, interval_sec, bars)
        type(ohlcv_series_t), intent(in) :: source_bars
        integer, intent(in) :: interval_sec
        type(bar_series_t), intent(out) :: bars
        integer, allocatable :: date_tmp(:), bin_tmp(:), end_tmp(:)
        real(dp), allocatable :: open_tmp(:), high_tmp(:), low_tmp(:), close_tmp(:)
        integer :: i, k, this_date, this_bin, sec

        if (interval_sec < 1) error stop "aggregate_bars: interval_sec must be positive"
        allocate(date_tmp(source_bars%nobs()), bin_tmp(source_bars%nobs()), end_tmp(source_bars%nobs()))
        allocate(open_tmp(source_bars%nobs()), high_tmp(source_bars%nobs()), low_tmp(source_bars%nobs()), &
                 close_tmp(source_bars%nobs()))
        k = 0
        do i = 1, source_bars%nobs()
            this_date = yyyymmdd(source_bars%timestamp(i)%date)
            sec = source_bars%timestamp(i)%seconds_since_midnight()
            this_bin = (sec - default_session_start_seconds) / interval_sec + 1
            if (this_bin < 1) cycle
            if (k < 1) then
                k = k + 1
                date_tmp(k) = this_date
                bin_tmp(k) = this_bin
                end_tmp(k) = default_session_start_seconds + this_bin*interval_sec
                open_tmp(k) = source_bars%open(i)
                high_tmp(k) = source_bars%high(i)
                low_tmp(k) = source_bars%low(i)
                close_tmp(k) = source_bars%close(i)
            else if (date_tmp(k) /= this_date .or. bin_tmp(k) /= this_bin) then
                k = k + 1
                date_tmp(k) = this_date
                bin_tmp(k) = this_bin
                end_tmp(k) = default_session_start_seconds + this_bin*interval_sec
                open_tmp(k) = source_bars%open(i)
                high_tmp(k) = source_bars%high(i)
                low_tmp(k) = source_bars%low(i)
                close_tmp(k) = source_bars%close(i)
            else
                high_tmp(k) = max(high_tmp(k), source_bars%high(i))
                low_tmp(k) = min(low_tmp(k), source_bars%low(i))
                close_tmp(k) = source_bars%close(i)
            end if
        end do
        if (k < 3) error stop "aggregate_bars: not enough aggregated bars"
        allocate(bars%date_id(k), bars%bin_id(k), bars%end_seconds(k), bars%open(k), bars%high(k), bars%low(k), bars%close(k))
        bars%date_id = date_tmp(1:k)
        bars%bin_id = bin_tmp(1:k)
        bars%end_seconds = end_tmp(1:k)
        bars%open = open_tmp(1:k)
        bars%high = high_tmp(1:k)
        bars%low = low_tmp(1:k)
        bars%close = close_tmp(1:k)
        deallocate(date_tmp, bin_tmp, end_tmp, open_tmp, high_tmp, low_tmp, close_tmp)
    end subroutine aggregate_bars

    ! Count close-to-close returns available within each trading day.
    integer function count_target_returns(bars) result(nret)
        type(bar_series_t), intent(in) :: bars
        integer :: i

        nret = 0
        do i = 2, size(bars%close)
            if (bars%date_id(i) == bars%date_id(i - 1)) nret = nret + 1
        end do
    end function count_target_returns

    ! Build close-to-close returns from target bars within each trading day.
    subroutine build_target_returns(bars, target)
        type(bar_series_t), intent(in) :: bars
        type(target_returns_t), intent(out) :: target
        real(dp), allocatable :: ret_tmp(:)
        integer, allocatable :: date_tmp(:), bin_tmp(:)
        integer :: i, k

        allocate(ret_tmp(size(bars%close) - 1), date_tmp(size(bars%close) - 1), bin_tmp(size(bars%close) - 1))
        k = 0
        do i = 2, size(bars%close)
            if (bars%date_id(i) /= bars%date_id(i - 1)) cycle
            k = k + 1
            ret_tmp(k) = log(bars%close(i) / bars%close(i - 1))
            date_tmp(k) = bars%date_id(i)
            bin_tmp(k) = bars%bin_id(i)
        end do
        if (k < 3) error stop "build_target_returns: not enough target returns"
        allocate(target%returns(k), target%date_id(k), target%bin_id(k))
        target%returns = ret_tmp(1:k)
        target%date_id = date_tmp(1:k)
        target%bin_id = bin_tmp(1:k)
        deallocate(ret_tmp, date_tmp, bin_tmp)
    end subroutine build_target_returns

    ! Sum predictor-frequency variance proxies into the matching target-frequency bucket.
    subroutine predictor_proxy_for_target(pred_bars, target_sec, target, proxy_kind, proxy)
        type(bar_series_t), intent(in) :: pred_bars
        integer, intent(in) :: target_sec, proxy_kind
        type(target_returns_t), intent(in) :: target
        real(dp), intent(out) :: proxy(:)
        integer :: i, idx, match_idx, target_bin
        real(dp) :: r, proxy_i

        if (size(proxy) /= size(target%returns)) error stop "predictor_proxy_for_target: array sizes differ"
        proxy = 0.0_dp
        idx = 1
        do i = 2, size(pred_bars%close)
            if (pred_bars%date_id(i) /= pred_bars%date_id(i - 1)) cycle
            r = log(pred_bars%close(i) / pred_bars%close(i - 1))
            select case (proxy_kind)
            case (proxy_cc)
                proxy_i = r**2
            case (proxy_parkinson)
                proxy_i = log(max(pred_bars%high(i), min_var) / max(pred_bars%low(i), min_var))**2 / &
                          (4.0_dp * log(2.0_dp))
            case (proxy_garman_klass)
                proxy_i = garman_klass_proxy_one(pred_bars%open(i), pred_bars%high(i), &
                                                 pred_bars%low(i), pred_bars%close(i))
            case default
                error stop "predictor_proxy_for_target: unsupported proxy kind"
            end select
            target_bin = (pred_bars%end_seconds(i) - default_session_start_seconds - 1) / target_sec + 1
            call advance_target_index(target, pred_bars%date_id(i), target_bin, idx, match_idx)
            if (match_idx > 0) proxy(match_idx) = proxy(match_idx) + proxy_i
        end do
        proxy = max(proxy, min_var)
    end subroutine predictor_proxy_for_target

    ! Compute one Garman-Klass OHLC variance proxy.
    real(dp) function garman_klass_proxy_one(open_price, high_price, low_price, close_price)
        real(dp), intent(in) :: open_price, high_price, low_price, close_price
        real(dp) :: hl, co

        hl = log(max(high_price, min_var) / max(low_price, min_var))
        co = log(max(close_price, min_var) / max(open_price, min_var))
        garman_klass_proxy_one = max(0.5_dp*hl**2 - (2.0_dp*log(2.0_dp) - 1.0_dp)*co**2, min_var)
    end function garman_klass_proxy_one

    ! Advance a monotone target index to the supplied date/bin key.
    subroutine advance_target_index(target, date_id, bin_id, idx, match_idx)
        type(target_returns_t), intent(in) :: target
        integer, intent(in) :: date_id, bin_id
        integer, intent(inout) :: idx
        integer, intent(out) :: match_idx

        if (idx < 1) idx = 1
        do while (idx <= size(target%returns))
            if (target%date_id(idx) > date_id) exit
            if (target%date_id(idx) == date_id .and. target%bin_id(idx) >= bin_id) exit
            idx = idx + 1
        end do
        match_idx = 0
        if (idx <= size(target%returns)) then
            if (target%date_id(idx) == date_id .and. target%bin_id(idx) == bin_id) match_idx = idx
        end if
    end subroutine advance_target_index

    ! Return whether predictor_sec is valid for target_sec.
    logical function valid_predictor(target_sec, predictor_sec)
        integer, intent(in) :: target_sec, predictor_sec

        valid_predictor = predictor_sec <= target_sec .and. mod(target_sec, predictor_sec) == 0
    end function valid_predictor

    ! Return 1 for negative returns and 0 otherwise.
    function negative_return_mask(returns) result(mask)
        real(dp), intent(in) :: returns(:)
        real(dp), allocatable :: mask(:)

        allocate(mask(size(returns)))
        mask = merge(1.0_dp, 0.0_dp, returns < 0.0_dp)
    end function negative_return_mask

    ! Fit a single-proxy asymmetric EWMA row with separate non-negative downside weight.
    subroutine fit_asymmetric_ewma_row(row, proxy_name, predictor_sec, target, daily_var, diurnal_var, proxy_mat, &
                                       weights, combined_proxy, q, h)
        type(ewma_row_t), intent(out) :: row
        character(len=*), intent(in) :: proxy_name
        integer, intent(in) :: predictor_sec
        type(target_returns_t), intent(in) :: target
        real(dp), intent(in) :: daily_var(:), diurnal_var(:), proxy_mat(:, :)
        real(dp), intent(out) :: weights(:), combined_proxy(:), q(:), h(:)

        call fit_nonnegative_proxy_weights(target%returns, target%date_id, daily_var, diurnal_var, proxy_mat, &
                                           weights, combined_proxy, q, h)
        call intraday_ewma_multiplier_from_proxy(combined_proxy, target%date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(row, "EWMA asym", proxy_name, predictor_sec, 3, target%returns, h, combined_proxy)
        row%proxy_scale = row%proxy_to_target_var
    end subroutine fit_asymmetric_ewma_row

    ! Fit a non-negative weighted combination of all valid realized-variance proxies.
    subroutine fit_combined_ewma_row(row, target, daily_var, diurnal_var, proxy_mat, weights, combined_proxy, q, h)
        type(ewma_row_t), intent(out) :: row
        type(target_returns_t), intent(in) :: target
        real(dp), intent(in) :: daily_var(:), diurnal_var(:), proxy_mat(:, :)
        real(dp), intent(out) :: weights(:), combined_proxy(:), q(:), h(:)

        call fit_nonnegative_proxy_weights(target%returns, target%date_id, daily_var, diurnal_var, proxy_mat, &
                                           weights, combined_proxy, q, h)
        call intraday_ewma_multiplier_from_proxy(combined_proxy, target%date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(row, "EWMA combined", "NONNEG", -1, size(weights), target%returns, h, combined_proxy)
        row%proxy_scale = row%proxy_to_target_var
    end subroutine fit_combined_ewma_row

    ! Fit a combined model with baseline and downside proxy columns.
    subroutine fit_asymmetric_combined_ewma_row(row, target, daily_var, diurnal_var, proxy_mat, weights, &
                                                combined_proxy, q, h)
        type(ewma_row_t), intent(out) :: row
        type(target_returns_t), intent(in) :: target
        real(dp), intent(in) :: daily_var(:), diurnal_var(:), proxy_mat(:, :)
        real(dp), intent(out) :: weights(:), combined_proxy(:), q(:), h(:)

        call fit_nonnegative_proxy_weights(target%returns, target%date_id, daily_var, diurnal_var, proxy_mat, &
                                           weights, combined_proxy, q, h)
        call intraday_ewma_multiplier_from_proxy(combined_proxy, target%date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(row, "EWMA combined asym", "NONNEG_ASYM", -1, size(weights), target%returns, h, combined_proxy)
        row%proxy_scale = row%proxy_to_target_var
    end subroutine fit_asymmetric_combined_ewma_row

    ! Fit non-negative proxy weights by coordinate search on the Gaussian likelihood.
    subroutine fit_nonnegative_proxy_weights(returns, date_id, daily_var, diurnal_var, proxy_mat, weights, &
                                             combined_proxy, q, h)
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), proxy_mat(:, :)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: weights(:), combined_proxy(:), q(:), h(:)
        real(dp), allocatable :: upper(:), proxy_mean(:)
        real(dp) :: target_var
        integer :: j, iter, nproxy

        nproxy = size(proxy_mat, dim=2)
        if (size(weights) /= nproxy) error stop "fit_nonnegative_proxy_weights: weight size differs"
        if (size(combined_proxy) /= size(returns)) error stop "fit_nonnegative_proxy_weights: proxy size differs"
        allocate(upper(nproxy), proxy_mean(nproxy))

        target_var = max(sum(returns**2) / real(size(returns), dp), min_var)
        do j = 1, nproxy
            proxy_mean(j) = max(sum(proxy_mat(:, j)) / real(size(returns), dp), min_var)
            upper(j) = 20.0_dp * target_var / proxy_mean(j)
            weights(j) = target_var / (real(nproxy, dp) * proxy_mean(j))
        end do
        call weighted_proxy(proxy_mat, weights, combined_proxy)

        do iter = 1, 6
            do j = 1, nproxy
                call optimize_one_nonnegative_weight(j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                                     upper(j), weights, combined_proxy, q, h)
            end do
        end do
        combined_proxy = max(combined_proxy, min_var)
        deallocate(upper, proxy_mean)
    end subroutine fit_nonnegative_proxy_weights

    ! Optimize one non-negative coordinate while holding the other weights fixed.
    subroutine optimize_one_nonnegative_weight(j, returns, date_id, daily_var, diurnal_var, proxy_mat, upper, &
                                               weights, combined_proxy, q, h)
        integer, intent(in) :: j
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), proxy_mat(:, :), upper
        integer, intent(in) :: date_id(:)
        real(dp), intent(inout) :: weights(:), combined_proxy(:)
        real(dp), intent(out) :: q(:), h(:)
        real(dp), parameter :: inv_phi = 0.61803398874989484820_dp
        real(dp) :: a, b, c, d, fc, fd, f0, best_w, best_f
        integer :: iter

        combined_proxy = max(combined_proxy - weights(j)*proxy_mat(:, j), min_var)
        a = 0.0_dp
        b = max(upper, min_var)
        c = b - inv_phi*(b - a)
        d = a + inv_phi*(b - a)
        fc = combined_proxy_weight_loglik(c, j, returns, date_id, daily_var, diurnal_var, proxy_mat, combined_proxy, q, h)
        fd = combined_proxy_weight_loglik(d, j, returns, date_id, daily_var, diurnal_var, proxy_mat, combined_proxy, q, h)
        do iter = 1, 24
            if (fc < fd) then
                a = c
                c = d
                fc = fd
                d = a + inv_phi*(b - a)
                fd = combined_proxy_weight_loglik(d, j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                                  combined_proxy, q, h)
            else
                b = d
                d = c
                fd = fc
                c = b - inv_phi*(b - a)
                fc = combined_proxy_weight_loglik(c, j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                                  combined_proxy, q, h)
            end if
        end do
        best_w = 0.5_dp*(a + b)
        best_f = combined_proxy_weight_loglik(best_w, j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                             combined_proxy, q, h)
        f0 = combined_proxy_weight_loglik(0.0_dp, j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                          combined_proxy, q, h)
        if (f0 > best_f) best_w = 0.0_dp
        weights(j) = best_w
        combined_proxy = max(combined_proxy + weights(j)*proxy_mat(:, j), min_var)
    end subroutine optimize_one_nonnegative_weight

    ! Return likelihood after adding one weighted proxy column to a base proxy.
    real(dp) function combined_proxy_weight_loglik(weight, j, returns, date_id, daily_var, diurnal_var, proxy_mat, &
                                                   base_proxy, q, h)
        real(dp), intent(in) :: weight, returns(:), daily_var(:), diurnal_var(:), proxy_mat(:, :), base_proxy(:)
        integer, intent(in) :: j, date_id(:)
        real(dp), intent(out) :: q(:), h(:)
        real(dp), allocatable :: trial_proxy(:)

        allocate(trial_proxy(size(base_proxy)))
        trial_proxy = max(base_proxy + weight*proxy_mat(:, j), min_var)
        call intraday_ewma_multiplier_from_proxy(trial_proxy, date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        combined_proxy_weight_loglik = gaussian_loglik(returns, h)
        deallocate(trial_proxy)
    end function combined_proxy_weight_loglik

    ! Compute a weighted sum of proxy columns.
    subroutine weighted_proxy(proxy_mat, weights, proxy)
        real(dp), intent(in) :: proxy_mat(:, :), weights(:)
        real(dp), intent(out) :: proxy(:)
        integer :: j

        if (size(weights) /= size(proxy_mat, dim=2)) error stop "weighted_proxy: weight size differs"
        if (size(proxy) /= size(proxy_mat, dim=1)) error stop "weighted_proxy: proxy size differs"
        proxy = 0.0_dp
        do j = 1, size(weights)
            proxy = proxy + weights(j)*proxy_mat(:, j)
        end do
        proxy = max(proxy, min_var)
    end subroutine weighted_proxy

    ! Fit a positive scale for a realized-variance proxy and fill one EWMA row.
    subroutine fit_scaled_ewma_row(row, proxy_name, predictor_sec, target, daily_var, diurnal_var, proxy, q, h)
        type(ewma_row_t), intent(out) :: row
        character(len=*), intent(in) :: proxy_name
        integer, intent(in) :: predictor_sec
        type(target_returns_t), intent(in) :: target
        real(dp), intent(in) :: daily_var(:), diurnal_var(:), proxy(:)
        real(dp), intent(out) :: q(:), h(:)
        real(dp) :: scale

        scale = fitted_proxy_scale(target%returns, target%date_id, daily_var, diurnal_var, proxy, q, h)
        call intraday_ewma_multiplier_from_proxy(scale*proxy, target%date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        call fill_row(row, "EWMA", proxy_name, predictor_sec, 2, target%returns, h, proxy)
        row%proxy_scale = scale
    end subroutine fit_scaled_ewma_row

    ! Estimate the positive variance-proxy scale that maximizes Gaussian likelihood.
    real(dp) function fitted_proxy_scale(returns, date_id, daily_var, diurnal_var, proxy, q, h) result(scale)
        real(dp), intent(in) :: returns(:), daily_var(:), diurnal_var(:), proxy(:)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: q(:), h(:)
        real(dp), parameter :: inv_phi = 0.61803398874989484820_dp
        real(dp) :: target_var, proxy_mean, center, a, b, c, d, fc, fd
        integer :: iter

        target_var = max(sum(returns**2) / real(size(returns), dp), min_var)
        proxy_mean = max(sum(proxy) / real(size(proxy), dp), min_var)
        center = log(target_var / proxy_mean)
        a = center - scale_log_search_half_width
        b = center + scale_log_search_half_width
        c = b - inv_phi*(b - a)
        d = a + inv_phi*(b - a)
        fc = scaled_proxy_loglik(c, returns, date_id, daily_var, diurnal_var, proxy, q, h)
        fd = scaled_proxy_loglik(d, returns, date_id, daily_var, diurnal_var, proxy, q, h)
        do iter = 1, 48
            if (fc < fd) then
                a = c
                c = d
                fc = fd
                d = a + inv_phi*(b - a)
                fd = scaled_proxy_loglik(d, returns, date_id, daily_var, diurnal_var, proxy, q, h)
            else
                b = d
                d = c
                fd = fc
                c = b - inv_phi*(b - a)
                fc = scaled_proxy_loglik(c, returns, date_id, daily_var, diurnal_var, proxy, q, h)
            end if
        end do
        scale = exp(0.5_dp*(a + b))
    end function fitted_proxy_scale

    ! Return likelihood for one log-scale value.
    real(dp) function scaled_proxy_loglik(log_scale, returns, date_id, daily_var, diurnal_var, proxy, q, h)
        real(dp), intent(in) :: log_scale, returns(:), daily_var(:), diurnal_var(:), proxy(:)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: q(:), h(:)

        call intraday_ewma_multiplier_from_proxy(exp(log_scale)*proxy, date_id, daily_var, diurnal_var, intraday_lambda, q)
        call intraday_variance_forecast(daily_var, diurnal_var*q, h)
        scaled_proxy_loglik = gaussian_loglik(returns, h)
    end function scaled_proxy_loglik

    ! Fill one likelihood comparison row.
    subroutine fill_row(row, model, proxy_name, predictor_sec, k, returns, h, proxy_values)
        type(ewma_row_t), intent(out) :: row
        character(len=*), intent(in) :: model, proxy_name
        integer, intent(in) :: predictor_sec, k
        real(dp), intent(in) :: returns(:), h(:)
        real(dp), optional, intent(in) :: proxy_values(:)
        real(dp) :: target_var

        row%model = model
        row%proxy = proxy_name
        row%predictor_sec = predictor_sec
        row%k = k
        if (k == 0) row%proxy_scale = 0.0_dp
        row%loglik = gaussian_loglik(returns, h)
        row%aic = -2.0_dp*row%loglik + 2.0_dp*real(k, dp)
        row%bic = -2.0_dp*row%loglik + log(real(size(returns), dp))*real(k, dp)
        if (compare_student_t_noise) then
            row%t_nu = fitted_t_nu(returns, h)
            row%t_loglik = student_t_loglik(returns, h, row%t_nu)
            row%t_aic = -2.0_dp*row%t_loglik + 2.0_dp*real(k + 1, dp)
            row%t_bic = -2.0_dp*row%t_loglik + log(real(size(returns), dp))*real(k + 1, dp)
        end if
        if (present(proxy_values)) then
            if (size(proxy_values) /= size(returns)) error stop "fill_row: proxy size differs"
            target_var = max(sum(returns**2) / real(size(returns), dp), min_var)
            row%proxy_mean = sum(proxy_values) / real(size(proxy_values), dp)
            row%proxy_to_target_var = row%proxy_mean / target_var
            row%proxy_max = maxval(proxy_values)
            row%proxy_nonzero_pct = 100.0_dp * real(count(proxy_values > min_var), dp) / real(size(proxy_values), dp)
        end if
    end subroutine fill_row

    ! Gaussian zero-mean log likelihood for returns and variance forecasts.
    real(dp) function gaussian_loglik(returns, h)
        real(dp), intent(in) :: returns(:), h(:)
        integer :: i

        if (size(returns) /= size(h)) error stop "gaussian_loglik: array sizes differ"
        gaussian_loglik = 0.0_dp
        do i = 1, size(returns)
            gaussian_loglik = gaussian_loglik - log_sqrt_2pi - 0.5_dp*log(max(h(i), min_var)) - &
                              0.5_dp*returns(i)**2 / max(h(i), min_var)
        end do
    end function gaussian_loglik

    ! Fit Student-t degrees of freedom for fixed variance forecasts.
    real(dp) function fitted_t_nu(returns, h) result(nu)
        real(dp), intent(in) :: returns(:), h(:)
        real(dp), parameter :: inv_phi = 0.61803398874989484820_dp
        real(dp) :: a, b, c, d, fc, fd
        integer :: iter

        a = -12.0_dp
        b = 12.0_dp
        c = b - inv_phi*(b - a)
        d = a + inv_phi*(b - a)
        fc = student_t_loglik_from_q(c, returns, h)
        fd = student_t_loglik_from_q(d, returns, h)
        do iter = 1, 48
            if (fc < fd) then
                a = c
                c = d
                fc = fd
                d = a + inv_phi*(b - a)
                fd = student_t_loglik_from_q(d, returns, h)
            else
                b = d
                d = c
                fd = fc
                c = b - inv_phi*(b - a)
                fc = student_t_loglik_from_q(c, returns, h)
            end if
        end do
        nu = t_nu_from_unconstrained(0.5_dp*(a + b))
    end function fitted_t_nu

    ! Student-t log likelihood for an unconstrained degrees-of-freedom value.
    real(dp) function student_t_loglik_from_q(qnu, returns, h)
        real(dp), intent(in) :: qnu, returns(:), h(:)

        student_t_loglik_from_q = student_t_loglik(returns, h, t_nu_from_unconstrained(qnu))
    end function student_t_loglik_from_q

    ! Map real line to the configured Student-t degrees-of-freedom interval.
    real(dp) function t_nu_from_unconstrained(qnu) result(nu)
        real(dp), intent(in) :: qnu

        nu = t_nu_min + (t_nu_max - t_nu_min) / (1.0_dp + exp(-qnu))
    end function t_nu_from_unconstrained

    ! Student-t zero-mean log likelihood for unit-variance t innovations.
    real(dp) function student_t_loglik(returns, h, nu)
        real(dp), intent(in) :: returns(:), h(:), nu
        real(dp) :: c, hi
        integer :: i

        if (size(returns) /= size(h)) error stop "student_t_loglik: array sizes differ"
        c = log_gamma(0.5_dp*(nu + 1.0_dp)) - log_gamma(0.5_dp*nu) - 0.5_dp*log(acos(-1.0_dp)*(nu - 2.0_dp))
        student_t_loglik = 0.0_dp
        do i = 1, size(returns)
            hi = max(h(i), min_var)
            student_t_loglik = student_t_loglik + c - 0.5_dp*log(hi) - &
                               0.5_dp*(nu + 1.0_dp)*log(1.0_dp + returns(i)**2 / ((nu - 2.0_dp)*hi))
        end do
    end function student_t_loglik

    ! Print one target-frequency comparison table.
    subroutine print_target_table(target_sec, target, rows, aic_best, bic_best)
        integer, intent(in) :: target_sec, aic_best, bic_best
        type(target_returns_t), intent(in) :: target
        type(ewma_row_t), intent(in) :: rows(:)
        integer :: i

        print '(A,I0,A)', "Target interval: ", target_sec, " seconds"
        print '(A,I0,A,A,A,A)', "Target returns: ", size(target%returns), " from ", date_label(target%date_id(1)), &
              " to ", date_label(target%date_id(size(target%date_id)))
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A16,1X,A10,1X,A5,1X,A12,1X,A14,1X,A14,1X,A14)', &
              "Model", "Proxy", "pred_sec", "k", "scale", "logL", "AIC", "BIC"
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            if (rows(i)%predictor_sec > 0) then
                print '(A32,1X,A16,1X,I10,1X,I5,1X,ES12.4,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), rows(i)%predictor_sec, rows(i)%k, rows(i)%proxy_scale, &
                      rows(i)%loglik, rows(i)%aic, rows(i)%bic
            else if (rows(i)%predictor_sec < 0) then
                print '(A32,1X,A16,1X,A10,1X,I5,1X,ES12.4,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), "multi", rows(i)%k, rows(i)%proxy_scale, &
                      rows(i)%loglik, rows(i)%aic, rows(i)%bic
            else
                print '(A32,1X,A16,1X,A10,1X,I5,1X,A12,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), "-", rows(i)%k, "-", rows(i)%loglik, rows(i)%aic, rows(i)%bic
            end if
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        call print_selected("AIC", rows(aic_best))
        call print_selected("BIC", rows(bic_best))
        print '(A)', ""
        if (compare_student_t_noise) call print_student_t_table(rows)
        call print_proxy_diagnostics(target_sec, rows)
    end subroutine print_target_table

    ! Print a selected model row label.
    subroutine print_selected(criterion, row)
        character(len=*), intent(in) :: criterion
        type(ewma_row_t), intent(in) :: row

        if (row%predictor_sec > 0) then
            print '(A,A,A,A,A,A,I0)', criterion, " selects: ", trim(row%model), " ", trim(row%proxy), &
                  " pred_sec=", row%predictor_sec
        else if (row%predictor_sec < 0) then
            print '(A,A,A,A,A)', criterion, " selects: ", trim(row%model), " ", trim(row%proxy)
        else
            print '(A,A,A)', criterion, " selects: ", trim(row%model)
        end if
    end subroutine print_selected

    ! Print Student-t likelihood scores for the same variance forecasts.
    subroutine print_student_t_table(rows)
        type(ewma_row_t), intent(in) :: rows(:)
        integer :: i, aic_best, bic_best

        aic_best = minloc(rows%t_aic, dim=1)
        bic_best = minloc(rows%t_bic, dim=1)
        print '(A)', "Student-t noise scores for the same variance forecasts"
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A16,1X,A10,1X,A5,1X,A9,1X,A14,1X,A14,1X,A14)', &
              "Model", "Proxy", "pred_sec", "k", "nu", "logL", "AIC", "BIC"
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            if (rows(i)%predictor_sec > 0) then
                print '(A32,1X,A16,1X,I10,1X,I5,1X,F9.3,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), rows(i)%predictor_sec, rows(i)%k + 1, rows(i)%t_nu, &
                      rows(i)%t_loglik, rows(i)%t_aic, rows(i)%t_bic
            else if (rows(i)%predictor_sec < 0) then
                print '(A32,1X,A16,1X,A10,1X,I5,1X,F9.3,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), "multi", rows(i)%k + 1, rows(i)%t_nu, &
                      rows(i)%t_loglik, rows(i)%t_aic, rows(i)%t_bic
            else
                print '(A32,1X,A16,1X,A10,1X,I5,1X,F9.3,1X,F14.3,1X,F14.3,1X,F14.3)', trim(rows(i)%model), &
                      trim(rows(i)%proxy), "-", rows(i)%k + 1, rows(i)%t_nu, rows(i)%t_loglik, rows(i)%t_aic, rows(i)%t_bic
            end if
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------------------"
        call print_t_selected("Student-t AIC", rows(aic_best))
        call print_t_selected("Student-t BIC", rows(bic_best))
        print '(A)', ""
    end subroutine print_student_t_table

    ! Print a selected Student-t model row label.
    subroutine print_t_selected(criterion, row)
        character(len=*), intent(in) :: criterion
        type(ewma_row_t), intent(in) :: row

        if (row%predictor_sec > 0) then
            print '(A,A,A,A,A,A,I0,A,F6.2)', criterion, " selects: ", trim(row%model), " ", trim(row%proxy), &
                  " pred_sec=", row%predictor_sec, " nu=", row%t_nu
        else if (row%predictor_sec < 0) then
            print '(A,A,A,A,A,A,F6.2)', criterion, " selects: ", trim(row%model), " ", trim(row%proxy), " nu=", row%t_nu
        else
            print '(A,A,A,A,F6.2)', criterion, " selects: ", trim(row%model), " nu=", row%t_nu
        end if
    end subroutine print_t_selected

    ! Print proxy scale diagnostics for EWMA rows.
    subroutine print_proxy_diagnostics(target_sec, rows)
        integer, intent(in) :: target_sec
        type(ewma_row_t), intent(in) :: rows(:)
        integer :: i

        print '(A,I0,A)', "Proxy diagnostics for target interval: ", target_sec, " seconds"
        print '(A)', "------------------------------------------------------------------------------------------------"
        print '(A16,1X,A10,1X,A14,1X,A16,1X,A14,1X,A12)', &
              "Proxy", "pred_sec", "mean_proxy", "mean/targetvar", "max_proxy", "nonzero_pct"
        print '(A)', "------------------------------------------------------------------------------------------------"
        do i = 1, size(rows)
            if (rows(i)%predictor_sec <= 0) cycle
            print '(A16,1X,I10,1X,ES14.5,1X,F16.5,1X,ES14.5,1X,F12.3)', &
                  trim(rows(i)%proxy), rows(i)%predictor_sec, rows(i)%proxy_mean, &
                  rows(i)%proxy_to_target_var, rows(i)%proxy_max, rows(i)%proxy_nonzero_pct
        end do
        print '(A)', "------------------------------------------------------------------------------------------------"
        print '(A)', ""
    end subroutine print_proxy_diagnostics

end module compare_intraday_ewma_freq_mod

program xcompare_intraday_ewma_freq
    use compare_intraday_ewma_freq_mod, only: run_compare_intraday_ewma_freq
    implicit none

    call run_compare_intraday_ewma_freq()
end program xcompare_intraday_ewma_freq

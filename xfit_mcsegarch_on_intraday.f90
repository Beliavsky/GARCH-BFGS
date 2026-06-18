! Fit joint overnight and intraday multiplicative-component EGARCH volatility models.
!
! The fitted likelihood includes close-to-open overnight returns and
! regular-session intraday close-to-close returns.  Overnight returns have an
! estimated variance weight, while intraday returns use a diurnal curve and
! EGARCH intraday dynamics.

module fit_mcsegarch_on_intraday_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi, pi
    use date_mod, only: date_label, yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, intraday_bin_ids
    use garch_mcsgarch_mod, only: estimate_diurnal_variance
    use distributions_mod, only: pdf_fs_skewt
    use bfgs_mod, only: bfgs_minimize
    use stats_mod, only: mean, variance
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    implicit none
    private

    character(len=*), parameter :: file_pattern = "c:\python\intraday_prices\*.csv"
    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: beta_abs_max = 0.999_dp
    real(dp), parameter :: egarch_abs_z_mean = sqrt(2.0_dp / pi)
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t = 2
    integer, parameter :: dist_fs_skewt = 3
    integer, parameter :: lambda_zero = 1
    integer, parameter :: lambda_est = 2
    integer, parameter :: lambda_one = 3
    character(len=12), parameter :: model_names(*) = [character(len=12) :: &
        "MCSEGARCH"]
    character(len=8), parameter :: dist_names(*) = [character(len=8) :: &
        "NORMAL", "T", "FS_SKEWT"]
    character(len=6), parameter :: lambda_names(*) = [character(len=6) :: &
        "ON0", "ONEST", "ON1"]

    type :: mcsegarch_params_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.06_dp
        real(dp) :: gamma = 0.0_dp
        real(dp) :: beta = 0.90_dp
    end type mcsegarch_params_t

    type :: mcsegarch_fit_result_t
        type(mcsegarch_params_t) :: params
        real(dp) :: nll = huge(1.0_dp)
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        real(dp) :: nu = 0.0_dp
        real(dp) :: xi = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type mcsegarch_fit_result_t

    type :: on_intraday_fit_t
        type(mcsegarch_fit_result_t) :: fit
        real(dp) :: overnight_weight = 1.0_dp
        real(dp) :: lambda_on = 1.0_dp
        real(dp) :: theta_on = 0.0_dp
        real(dp) :: gamma_on = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: fit_sec = 0.0_dp
    end type on_intraday_fit_t

    real(dp), allocatable, save :: obj_intra_ret(:), obj_intra_base(:), obj_intra_on_ret(:), obj_diurnal_var(:)
    real(dp), allocatable, save :: obj_on_ret(:), obj_on_scale(:)
    integer, save :: obj_dist = dist_normal
    integer, save :: obj_lambda_mode = lambda_one

    public :: run_fit_mcsegarch_on_intraday

contains

    ! Fit all configured models to each input file obtained from command-line, glob, or default.
    subroutine run_fit_mcsegarch_on_intraday()
        character(len=MAX_PATH_LEN), allocatable :: filenames(:)
        integer :: i

        call collect_input_filenames(filenames, &
            file_pattern=file_pattern, &
            default_filenames=[character(len=MAX_PATH_LEN) :: &
                "c:\python\intraday_prices\spy_5min_databento.csv"])
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call fit_one_file(trim(filenames(i)))
        end do
        deallocate(filenames)
    end subroutine run_fit_mcsegarch_on_intraday

    ! Read OHLCV prices, build overnight/intraday returns, and fit all configured models.
    subroutine fit_one_file(filename)
        character(len=*), intent(in) :: filename
        type(ohlcv_series_t) :: bars, regular_bars
        real(dp), allocatable :: intra_ret(:), intra_base(:), intra_on_ret(:), on_ret(:), on_scale(:)
        real(dp), allocatable :: diurnal_fit(:,:), q_fit(:,:)
        integer, allocatable :: bin_id(:), intra_dates(:), on_dates(:)
        type(on_intraday_fit_t), allocatable :: results(:)
        integer :: nfit, ifit, imodel, idist, ilambda, max_iter
        real(dp) :: gtol, t0, t1, read_sec, elapsed_sec

        max_iter = 120
        gtol = 1.0e-5_dp

        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), bars)
        call cpu_time(t1)
        read_sec = t1 - t0
        call filter_intraday_session(bars, regular_bars)
        call build_on_intraday_inputs(regular_bars, intra_ret, intra_base, intra_on_ret, bin_id, intra_dates, &
                                      on_ret, on_scale, on_dates)

        nfit = size(model_names) * size(dist_names) * size(lambda_names)
        allocate(results(nfit), diurnal_fit(size(intra_ret), nfit), q_fit(size(intra_ret), nfit))

        call cpu_time(t0)
        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                do ilambda = 1, size(lambda_names)
                    ifit = ifit + 1
                    call fit_one_joint_model(trim(model_names(imodel)), trim(dist_names(idist)), ilambda, &
                                             intra_ret, intra_base, intra_on_ret, bin_id, on_ret, on_scale, &
                                             max_iter, gtol, results(ifit), diurnal_fit(:, ifit), q_fit(:, ifit))
                end do
            end do
        end do
        call cpu_time(t1)
        elapsed_sec = read_sec + (t1 - t0)

        call print_summary(trim(filename), regular_bars%nobs(), intra_ret, intra_base, intra_on_ret, bin_id, intra_dates, &
                           on_ret, on_scale, on_dates, results, read_sec, elapsed_sec)

        deallocate(intra_ret, intra_base, intra_on_ret, bin_id, intra_dates, on_ret, on_scale, on_dates)
        deallocate(results, diurnal_fit, q_fit)
    end subroutine fit_one_file

    ! Fit one dynamic model and innovation distribution to the joint likelihood.
    subroutine fit_one_joint_model(model_name, dist_name, lambda_mode, intra_ret, intra_base, intra_on_ret, &
                                   bin_id, on_ret, on_scale, max_iter, gtol, result, diurnal_var, q)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: intra_ret(:), intra_base(:), intra_on_ret(:), on_ret(:), on_scale(:), gtol
        integer, intent(in) :: lambda_mode, bin_id(:), max_iter
        type(on_intraday_fit_t), intent(out) :: result
        real(dp), intent(out) :: diurnal_var(:), q(:)
        real(dp), allocatable :: p(:), p_start(:), grad(:), initial_intra_scale(:), fitted_intra_scale(:)
        real(dp) :: t0, t1
        integer :: dist, np, niter, nobs
        logical :: converged

        if (trim(model_name) /= "MCSEGARCH") error stop "fit_one_joint_model: unsupported model"
        dist = dist_id(dist_name)
        np = model_param_count() + 1 + lambda_param_count(lambda_mode) + dist_param_count(dist)
        allocate(p(np), p_start(np), grad(np), initial_intra_scale(size(intra_ret)), fitted_intra_scale(size(intra_ret)))
        call intraday_scale_from_overnight(intra_base, intra_on_ret, initial_lambda_on(lambda_mode), &
                                           initial_gamma_on(lambda_mode), initial_intra_scale)
        call estimate_diurnal_variance(intra_ret, initial_intra_scale, bin_id, diurnal_var, &
                                       smooth_diurnal_curve, diurnal_smooth_half_width)
        call set_objective_data(intra_ret, intra_base, intra_on_ret, diurnal_var, on_ret, on_scale, &
                                dist, lambda_mode)
        call pack_joint_params(initial_params(intra_ret, initial_intra_scale, diurnal_var), &
                               initial_overnight_weight(on_ret, on_scale), initial_lambda_on(lambda_mode), &
                               initial_gamma_on(lambda_mode), dist, lambda_mode, p)
        p_start = p

        call cpu_time(t0)
        call bfgs_minimize(joint_obj, p, np, max_iter, gtol, result%fit%nll, niter, converged)
        call cpu_time(t1)
        result%fit_sec = t1 - t0
        call unpack_joint_params(p, dist, result%fit%params, result%overnight_weight, &
                                 result%lambda_on, result%gamma_on, result%fit%nu, result%fit%xi, lambda_mode)
        result%theta_on = 0.0_dp
        call joint_obj(p, np, result%fit%nll, grad)
        if (ieee_bad(result%fit%nll) .or. result%fit%nll >= huge(1.0_dp) / 1000.0_dp) then
            p = p_start
            call unpack_joint_params(p, dist, result%fit%params, result%overnight_weight, &
                                     result%lambda_on, result%gamma_on, result%fit%nu, result%fit%xi, lambda_mode)
            call joint_obj(p, np, result%fit%nll, grad)
            converged = .false.
        end if
        call intraday_scale_from_overnight(intra_base, intra_on_ret, result%lambda_on, result%gamma_on, &
                                           fitted_intra_scale)
        call mcsegarch_filter(intra_ret, fitted_intra_scale, diurnal_var, result%fit%params, q)

        nobs = size(intra_ret) + size(on_ret)
        result%fit%loglik = -result%fit%nll * real(nobs, dp)
        result%fit%persist = result%fit%params%beta
        result%fit%niter = niter
        result%fit%converged = converged
        result%aic = -2.0_dp*result%fit%loglik + 2.0_dp*real(model_param_count() + 1 + &
                     lambda_param_count(lambda_mode) + &
                     dist_param_count(dist) + count_populated_bins(bin_id), dp)
        result%bic = -2.0_dp*result%fit%loglik + log(real(nobs, dp)) * real(model_param_count() + 1 + &
                     lambda_param_count(lambda_mode) + &
                     dist_param_count(dist) + count_populated_bins(bin_id), dp)
        deallocate(p, p_start, grad, initial_intra_scale, fitted_intra_scale)
    end subroutine fit_one_joint_model

    ! Store data used by the BFGS objective.
    subroutine set_objective_data(intra_ret, intra_base, intra_on_ret, diurnal_var, on_ret, on_scale, &
                                  dist, lambda_mode)
        real(dp), intent(in) :: intra_ret(:), intra_base(:), intra_on_ret(:), diurnal_var(:), on_ret(:), on_scale(:)
        integer, intent(in) :: dist, lambda_mode

        if (allocated(obj_intra_ret)) deallocate(obj_intra_ret, obj_intra_base, obj_intra_on_ret, obj_diurnal_var)
        if (allocated(obj_on_ret)) deallocate(obj_on_ret, obj_on_scale)
        allocate(obj_intra_ret(size(intra_ret)), obj_intra_base(size(intra_base)), obj_intra_on_ret(size(intra_on_ret)), &
                 obj_diurnal_var(size(diurnal_var)), obj_on_ret(size(on_ret)), obj_on_scale(size(on_scale)))
        obj_intra_ret = intra_ret
        obj_intra_base = max(intra_base, min_var)
        obj_intra_on_ret = intra_on_ret
        obj_diurnal_var = max(diurnal_var, min_var)
        obj_on_ret = on_ret
        obj_on_scale = max(on_scale, min_var)
        obj_dist = dist
        obj_lambda_mode = lambda_mode
    end subroutine set_objective_data

    ! BFGS objective and finite-difference gradient for the joint likelihood.
    subroutine joint_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = joint_nll_from_p(p)
        if (ieee_bad(f) .or. f >= huge(1.0_dp) / 1000.0_dp) then
            f = huge(1.0_dp) / 1000.0_dp
            g = 0.0_dp
            return
        end if
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = joint_nll_from_p(pp)
            fm = joint_nll_from_p(pm)
            if (ieee_bad(fp) .or. ieee_bad(fm) .or. fp >= huge(1.0_dp) / 1000.0_dp .or. &
                fm >= huge(1.0_dp) / 1000.0_dp) then
                g(j) = 0.0_dp
            else
                g(j) = (fp - fm) / (2.0_dp * step)
            end if
        end do
        deallocate(pp, pm)
    end subroutine joint_obj

    ! Return average negative log likelihood for packed parameters.
    real(dp) function joint_nll_from_p(p)
        real(dp), intent(in) :: p(:)
        type(mcsegarch_params_t) :: params
        real(dp), allocatable :: q(:), intra_scale(:)
        real(dp) :: on_weight, lambda_on, gamma_on, nu, xi, loss, h
        integer :: i

        call unpack_joint_params(p, obj_dist, params, on_weight, lambda_on, gamma_on, nu, xi, obj_lambda_mode)
        if (.not. params_valid(params) .or. on_weight <= 0.0_dp .or. &
            lambda_on < 0.0_dp .or. gamma_on < 0.0_dp) then
            joint_nll_from_p = huge(1.0_dp) / 10.0_dp
            return
        end if
        allocate(q(size(obj_intra_ret)), intra_scale(size(obj_intra_ret)))
        call intraday_scale_from_overnight(obj_intra_base, obj_intra_on_ret, lambda_on, gamma_on, intra_scale)
        call mcsegarch_filter(obj_intra_ret, intra_scale, obj_diurnal_var, params, q)
        loss = 0.0_dp
        do i = 1, size(obj_on_ret)
            h = max(obj_on_scale(i) * on_weight, min_var)
            loss = loss + innovation_nll(obj_on_ret(i), h, obj_dist, nu, xi)
        end do
        do i = 1, size(obj_intra_ret)
            h = max(intra_scale(i) * obj_diurnal_var(i) * q(i), min_var)
            loss = loss + innovation_nll(obj_intra_ret(i), h, obj_dist, nu, xi)
        end do
        if (ieee_bad(loss)) then
            joint_nll_from_p = huge(1.0_dp) / 10.0_dp
            deallocate(q, intra_scale)
            return
        end if
        joint_nll_from_p = loss / real(size(obj_on_ret) + size(obj_intra_ret), dp)
        deallocate(q, intra_scale)
    end function joint_nll_from_p

    ! Convert regular-session bars into overnight and intraday model inputs.
    subroutine build_on_intraday_inputs(bars, intra_ret, intra_base, intra_on_ret, bin_id, intra_dates, &
                                        on_ret, on_scale, on_dates)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: intra_ret(:), intra_base(:), intra_on_ret(:), on_ret(:), on_scale(:)
        integer, allocatable, intent(out) :: bin_id(:), intra_dates(:), on_dates(:)
        integer, allocatable :: bar_bins(:), day_index(:), day_dates(:), day_first(:), day_last(:)
        real(dp), allocatable :: intraday_rv(:), overnight_ret(:), intra_day_base(:), on_day_scale(:)
        real(dp), allocatable :: intra_all(:), intra_base_all(:), intra_on_ret_all(:), on_all(:), on_scale_all(:)
        integer, allocatable :: bin_all(:), intra_dates_all(:), on_dates_all(:)
        integer :: n, ndays, i, k, d, kon
        real(dp) :: fallback

        n = bars%nobs()
        if (n < 3) error stop "build_on_intraday_inputs: not enough regular-session bars"
        call intraday_bin_ids(bars, bar_bins)
        call map_days(bars, day_index, day_dates, day_first, day_last)
        ndays = size(day_dates)
        allocate(intraday_rv(ndays), overnight_ret(ndays), intra_day_base(ndays), on_day_scale(ndays))
        intraday_rv = 0.0_dp
        overnight_ret = 0.0_dp
        do i = 2, n
            if (day_index(i) == day_index(i - 1)) then
                intraday_rv(day_index(i)) = intraday_rv(day_index(i)) + log(bars%close(i) / bars%close(i - 1))**2
            end if
        end do
        do d = 2, ndays
            overnight_ret(d) = log(bars%open(day_first(d)) / bars%close(day_last(d - 1)))
        end do
        fallback = max(mean(intraday_rv(max(1, min(2, ndays)):ndays)), min_var)
        on_day_scale(1) = fallback
        intra_day_base(1) = fallback
        do d = 2, ndays
            on_day_scale(d) = max(intraday_rv(d - 1), min_var)
            intra_day_base(d) = max(intraday_rv(d - 1), min_var)
        end do

        allocate(intra_all(n - 1), intra_base_all(n - 1), intra_on_ret_all(n - 1), bin_all(n - 1), intra_dates_all(n - 1))
        allocate(on_all(max(1, ndays - 1)), on_scale_all(max(1, ndays - 1)), on_dates_all(max(1, ndays - 1)))
        k = 0
        do i = 2, n
            if (day_index(i) /= day_index(i - 1)) cycle
            k = k + 1
            intra_all(k) = log(bars%close(i) / bars%close(i - 1))
            intra_base_all(k) = intra_day_base(day_index(i))
            intra_on_ret_all(k) = overnight_ret(day_index(i))
            bin_all(k) = bar_bins(i)
            intra_dates_all(k) = day_dates(day_index(i))
        end do
        kon = 0
        do d = 2, ndays
            kon = kon + 1
            on_all(kon) = log(bars%open(day_first(d)) / bars%close(day_last(d - 1)))
            on_scale_all(kon) = on_day_scale(d)
            on_dates_all(kon) = day_dates(d)
        end do
        if (k < 3 .or. kon < 2) error stop "build_on_intraday_inputs: not enough returns"
        allocate(intra_ret(k), intra_base(k), intra_on_ret(k), bin_id(k), intra_dates(k), on_ret(kon), on_scale(kon), &
                 on_dates(kon))
        intra_ret = intra_all(1:k)
        intra_base = intra_base_all(1:k)
        intra_on_ret = intra_on_ret_all(1:k)
        bin_id = bin_all(1:k)
        intra_dates = intra_dates_all(1:k)
        on_ret = on_all(1:kon)
        on_scale = on_scale_all(1:kon)
        on_dates = on_dates_all(1:kon)

        deallocate(bar_bins, day_index, day_dates, day_first, day_last, intraday_rv, overnight_ret)
        deallocate(intra_day_base, on_day_scale, intra_all, intra_base_all, intra_on_ret_all, bin_all, intra_dates_all)
        deallocate(on_all, on_scale_all, on_dates_all)
    end subroutine build_on_intraday_inputs

    ! Map each bar to a trading day and record the first and last bar of each day.
    subroutine map_days(bars, day_index, day_dates, day_first, day_last)
        type(ohlcv_series_t), intent(in) :: bars
        integer, allocatable, intent(out) :: day_index(:), day_dates(:), day_first(:), day_last(:)
        integer :: n, i, ndays, current_date

        n = bars%nobs()
        allocate(day_index(n), day_dates(n), day_first(n), day_last(n))
        ndays = 0
        current_date = -1
        do i = 1, n
            if (yyyymmdd(bars%timestamp(i)%date) /= current_date) then
                ndays = ndays + 1
                current_date = yyyymmdd(bars%timestamp(i)%date)
                day_dates(ndays) = current_date
                day_first(ndays) = i
                if (ndays > 1) day_last(ndays - 1) = i - 1
            end if
            day_index(i) = ndays
        end do
        day_last(ndays) = n
        day_dates = day_dates(1:ndays)
        day_first = day_first(1:ndays)
        day_last = day_last(1:ndays)
    end subroutine map_days

    ! Print data summary, fitted parameters, and IC ranking table.
    subroutine print_summary(filename, n_bars, intra_ret, intra_base, intra_on_ret, bin_id, intra_dates, on_ret, on_scale, &
                             on_dates, results, read_sec, elapsed_sec)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_bars
        real(dp), intent(in) :: intra_ret(:), intra_base(:), intra_on_ret(:), on_ret(:), on_scale(:)
        integer, intent(in) :: bin_id(:), intra_dates(:), on_dates(:)
        type(on_intraday_fit_t), intent(in) :: results(:)
        real(dp), intent(in) :: read_sec, elapsed_sec
        integer :: ifit, imodel, idist, ilambda, aic_best, bic_best

        print '(A)', "Joint overnight/intraday MCS-EGARCH fits"
        print '(A,A)', "Input file: ", trim(filename)
        print '(A,I0)', "Regular-session bars: ", n_bars
        print '(A,I0,A,A,A,A)', "Intraday returns: ", size(intra_ret), " from ", &
              date_label(intra_dates(1)), " to ", date_label(intra_dates(size(intra_dates)))
        print '(A,I0,A,A,A,A)', "Overnight returns: ", size(on_ret), " from ", &
              date_label(on_dates(1)), " to ", date_label(on_dates(size(on_dates)))
        print '(A,I0)', "Intraday bins: ", count_populated_bins(bin_id)
        print '(A,ES12.4,A,ES12.4,A,ES12.4)', "Mean intraday base scale: ", mean(intra_base), &
              "  mean overnight shock: ", mean(intra_on_ret**2), &
              "  mean overnight scale: ", mean(on_scale)
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------"
        print '(A12,1X,A8,1X,A6,1X,A12,9(1X,A10),1X,A8,1X,A5,1X,A13,1X,A9)', &
              "Model", "Dist", "Lambda", "omega", "alpha", "gamma", "beta", "on_wt", "lambda_on", "gamma_on", "nu", &
              "xi", "persist", "iter", "conv", "logL", "fit_sec"
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------"
        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                do ilambda = 1, size(lambda_names)
                    ifit = ifit + 1
                    call print_fit_row(trim(model_names(imodel)), trim(dist_names(idist)), &
                                       trim(lambda_names(ilambda)), results(ifit))
                end do
            end do
        end do
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------"
        aic_best = minloc(results%aic, dim=1)
        bic_best = minloc(results%bic, dim=1)
        call print_ic_table(results)
        call print_lambda_impact_table(results)
        print '(A,A,A,A)', "AIC selects: ", model_table_label(trim(model_names(model_index(aic_best)))), &
              " ", trim(dist_names(dist_index(aic_best))) // " " // trim(lambda_names(lambda_index(aic_best)))
        print '(A,A,A,A)', "BIC selects: ", model_table_label(trim(model_names(model_index(bic_best)))), &
              " ", trim(dist_names(dist_index(bic_best))) // " " // trim(lambda_names(lambda_index(bic_best)))
        print '(A,F10.3)', "Data read seconds: ", read_sec
        print '(A,F10.3)', "Elapsed seconds:   ", elapsed_sec
    end subroutine print_summary

    ! Print one fitted parameter row.
    subroutine print_fit_row(model_name, dist_name, lambda_name, result)
        character(len=*), intent(in) :: model_name, dist_name, lambda_name
        type(on_intraday_fit_t), intent(in) :: result

        print '(A12,1X,A8,1X,A6,1X,ES12.4,9(1X,F10.4),1X,I8,5X,L1,1X,F13.3,1X,F9.3)', &
              model_name, dist_name, lambda_name, result%fit%params%omega, result%fit%params%alpha, &
              result%fit%params%gamma, result%fit%params%beta, result%overnight_weight, result%lambda_on, &
              result%gamma_on, result%fit%nu, result%fit%xi, result%fit%persist, &
              result%fit%niter, result%fit%converged, result%fit%loglik, result%fit_sec
    end subroutine print_fit_row

    ! Print the information criteria comparison table.
    subroutine print_ic_table(results)
        type(on_intraday_fit_t), intent(in) :: results(:)
        integer :: i

        print '(A)', ""
        print '(A)', "Joint overnight/intraday MCS-EGARCH model comparison"
        print '(A)', "----------------------------------------------------------------------------------------------"
        print '(A18,1X,A8,1X,A6,1X,A10,1X,A10,1X,A14,1X,A14,1X,A14,1X,A8,1X,A8,1X,A8)', &
              "Model", "Dist", "Lambda", "lambda_on", "gamma_on", "logL", "AIC", "BIC", "conv", "iter", &
              "fit_sec"
        print '(A)', "----------------------------------------------------------------------------------------------"
        do i = 1, size(results)
            print '(A18,1X,A8,1X,A6,1X,F10.4,1X,F10.4,1X,F14.3,1X,F14.3,1X,F14.3,3X,L1,1X,I8,1X,F8.3)', &
                  model_table_label(trim(model_names(model_index(i)))), trim(dist_names(dist_index(i))), &
                  trim(lambda_names(lambda_index(i))), results(i)%lambda_on, results(i)%gamma_on, &
                  results(i)%fit%loglik, results(i)%aic, results(i)%bic, results(i)%fit%converged, &
                  results(i)%fit%niter, results(i)%fit_sec
        end do
        print '(A)', "----------------------------------------------------------------------------------------------"
    end subroutine print_ic_table

    ! Print likelihood-ratio diagnostics for estimated overnight impact versus lambda_on = 0.
    subroutine print_lambda_impact_table(results)
        type(on_intraday_fit_t), intent(in) :: results(:)
        integer :: imodel, idist, idx0, idxest, df
        real(dp) :: lr_stat, p_value, delta_aic, delta_bic

        print '(A)', ""
        print '(A)', "Overnight impact on next-day intraday volatility"
        print '(A)', "----------------------------------------------------------------------------------------------------"
        print '(A18,1X,A8,1X,A10,1X,A10,1X,A4,1X,A12,1X,A12,1X,A12,1X,A12)', &
              "Model", "Dist", "lambda_on", "gamma_on", "df", "LR_vs_ON0", "p_value", &
              "dAIC_ONEST", "dBIC_ONEST"
        print '(A)', "----------------------------------------------------------------------------------------------------"
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                idx0 = flat_index(imodel, idist, lambda_zero)
                idxest = flat_index(imodel, idist, lambda_est)
                df = lambda_param_count(lambda_est)
                lr_stat = max(2.0_dp * (results(idxest)%fit%loglik - results(idx0)%fit%loglik), 0.0_dp)
                if (df == 2) then
                    p_value = exp(-0.5_dp * lr_stat)
                else
                    p_value = erfc(sqrt(0.5_dp * lr_stat))
                end if
                delta_aic = results(idxest)%aic - results(idx0)%aic
                delta_bic = results(idxest)%bic - results(idx0)%bic
                print '(A18,1X,A8,1X,F10.4,1X,F10.4,1X,I4,1X,F12.3,1X,ES12.4,1X,F12.3,1X,F12.3)', &
                      model_table_label(trim(model_names(imodel))), trim(dist_names(idist)), &
                      results(idxest)%lambda_on, results(idxest)%gamma_on, df, &
                      lr_stat, p_value, delta_aic, delta_bic
            end do
        end do
        print '(A)', "----------------------------------------------------------------------------------------------------"
    end subroutine print_lambda_impact_table

    ! Return the flattened result index for model, distribution, and lambda mode indices.
    integer function flat_index(imodel, idist, ilambda)
        integer, intent(in) :: imodel, idist, ilambda

        flat_index = ((idist - 1) * size(model_names) + imodel - 1) * size(lambda_names) + ilambda
    end function flat_index

    ! Return the configured model index for a flattened fit index.
    integer function model_index(ifit)
        integer, intent(in) :: ifit

        model_index = mod((ifit - 1) / size(lambda_names), size(model_names)) + 1
    end function model_index

    ! Return the configured distribution index for a flattened fit index.
    integer function dist_index(ifit)
        integer, intent(in) :: ifit

        dist_index = (ifit - 1) / (size(model_names) * size(lambda_names)) + 1
    end function dist_index

    ! Return the configured overnight-impact index for a flattened fit index.
    integer function lambda_index(ifit)
        integer, intent(in) :: ifit

        lambda_index = mod(ifit - 1, size(lambda_names)) + 1
    end function lambda_index

    ! Filter the intraday EGARCH multiplicative component q_t.
    subroutine mcsegarch_filter(y, daily_var, diurnal_var, params, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsegarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)
        real(dp) :: log_q, hscale, e, z
        integer :: t

        if (size(y) /= size(daily_var) .or. size(y) /= size(diurnal_var) .or. size(q) /= size(y)) then
            error stop "mcsegarch_filter: array sizes differ"
        end if
        log_q = params%omega / max(1.0_dp - params%beta, 1.0e-6_dp)
        do t = 1, size(y)
            log_q = max(min(log_q, 50.0_dp), -50.0_dp)
            q(t) = max(exp(log_q), min_var)
            hscale = max(daily_var(t) * diurnal_var(t), min_var)
            e = y(t) / sqrt(hscale)
            z = e / sqrt(q(t))
            log_q = params%omega + params%beta * log_q + &
                    params%alpha * (abs(z) - egarch_abs_z_mean) + params%gamma * z
        end do
    end subroutine mcsegarch_filter

    ! Build the intraday daily scale using an overnight EGARCH-style log impact.
    subroutine intraday_scale_from_overnight(base_scale, overnight_return, lambda_on, gamma_on, scale)
        real(dp), intent(in) :: base_scale(:), overnight_return(:), lambda_on, gamma_on
        real(dp), intent(out) :: scale(:)
        integer :: i
        real(dp) :: z, log_impact

        if (size(base_scale) /= size(overnight_return) .or. size(scale) /= size(base_scale)) then
            error stop "intraday_scale_from_overnight: array sizes differ"
        end if
        do i = 1, size(base_scale)
            z = overnight_return(i) / sqrt(max(base_scale(i), min_var))
            log_impact = lambda_on * (abs(z) - egarch_abs_z_mean) + gamma_on * z
            scale(i) = max(base_scale(i) * exp(max(min(log_impact, 50.0_dp), -50.0_dp)), min_var)
        end do
    end subroutine intraday_scale_from_overnight

    ! Pack constrained EGARCH, overnight, and distribution parameters.
    subroutine pack_joint_params(params, overnight_weight, lambda_on, gamma_on, dist, lambda_mode, p)
        type(mcsegarch_params_t), intent(in) :: params
        real(dp), intent(in) :: overnight_weight, lambda_on, gamma_on
        integer, intent(in) :: dist, lambda_mode
        real(dp), intent(out) :: p(:)
        integer :: j

        p(1) = params%omega
        p(2) = params%alpha
        p(3) = params%gamma
        p(4) = 0.5_dp * log((1.0_dp + params%beta) / (1.0_dp - params%beta))
        j = model_param_count() + 1
        p(j) = log(max(overnight_weight, min_var))
        if (lambda_mode == lambda_est) then
            j = j + 1
            p(j) = log(max(lambda_on, min_var))
            j = j + 1
            p(j) = gamma_on
        end if
        if (dist == dist_t) p(j + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
        if (dist == dist_fs_skewt) then
            p(j + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
            p(j + 2) = 0.0_dp
        end if
    end subroutine pack_joint_params

    ! Unpack optimizer variables into constrained EGARCH parameters.
    subroutine unpack_joint_params(p, dist, params, overnight_weight, lambda_on, gamma_on, nu, xi, lambda_mode)
        real(dp), intent(in) :: p(:)
        integer, intent(in) :: dist, lambda_mode
        type(mcsegarch_params_t), intent(out) :: params
        real(dp), intent(out) :: overnight_weight, lambda_on, gamma_on, nu, xi
        integer :: j

        params%omega = max(min(p(1), 50.0_dp), -50.0_dp)
        params%alpha = max(min(p(2), 50.0_dp), -50.0_dp)
        params%gamma = max(min(p(3), 50.0_dp), -50.0_dp)
        params%beta = beta_abs_max * tanh(p(4))
        j = model_param_count() + 1
        overnight_weight = exp(max(min(p(j), 50.0_dp), -50.0_dp))
        select case (lambda_mode)
        case (lambda_zero)
            lambda_on = 0.0_dp
            gamma_on = 0.0_dp
        case (lambda_est)
            j = j + 1
            lambda_on = exp(max(min(p(j), 50.0_dp), -50.0_dp))
            j = j + 1
            gamma_on = max(min(p(j), 50.0_dp), -50.0_dp)
        case (lambda_one)
            lambda_on = 1.0_dp
            gamma_on = 0.0_dp
        case default
            error stop "unpack_joint_params: unsupported lambda mode"
        end select
        if (dist == dist_t .or. dist == dist_fs_skewt) then
            nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(j + 1)))
        else
            nu = 0.0_dp
        end if
        if (dist == dist_fs_skewt) then
            xi = exp(max(min(p(j + 2), 20.0_dp), -20.0_dp))
        else
            xi = 0.0_dp
        end if
    end subroutine unpack_joint_params

    ! Negative log likelihood contribution for one zero-mean innovation.
    real(dp) function innovation_nll(y, h, dist, nu, xi)
        real(dp), intent(in) :: y, h, nu, xi
        integer, intent(in) :: dist
        real(dp) :: hh, z, pdf

        hh = max(h, min_var)
        if (dist == dist_t) then
            innovation_nll = log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu + 1.0_dp)) + &
                             0.5_dp*log(pi*(nu - 2.0_dp)) + 0.5_dp*log(hh) + &
                             0.5_dp*(nu + 1.0_dp)*log(1.0_dp + y**2 / ((nu - 2.0_dp)*hh))
        else if (dist == dist_fs_skewt) then
            z = y / sqrt(hh)
            pdf = max(pdf_fs_skewt(z, nu, xi), min_pdf)
            innovation_nll = 0.5_dp*log(hh) - log(pdf)
        else
            innovation_nll = log_sqrt_2pi + 0.5_dp*log(hh) + 0.5_dp*y**2 / hh
        end if
    end function innovation_nll

    ! Return true for NaN or effectively infinite values.
    logical function ieee_bad(x)
        real(dp), intent(in) :: x

        ieee_bad = (x /= x) .or. abs(x) > huge(1.0_dp) / 100.0_dp
    end function ieee_bad

    ! Build initial dynamic parameters from scaled intraday residual variance.
    type(mcsegarch_params_t) function initial_params(intra_ret, intra_scale, diurnal_var)
        real(dp), intent(in) :: intra_ret(:), intra_scale(:), diurnal_var(:)
        real(dp) :: vnorm

        vnorm = max(variance(intra_ret / sqrt(max(intra_scale * diurnal_var, min_var))), min_var)
        initial_params = mcsegarch_params_t(omega=log(vnorm) * (1.0_dp - 0.90_dp), alpha=0.06_dp, &
                                            gamma=-0.05_dp, beta=0.90_dp)
    end function initial_params

    ! Estimate the starting overnight variance weight by method of moments.
    real(dp) function initial_overnight_weight(on_ret, on_scale)
        real(dp), intent(in) :: on_ret(:), on_scale(:)

        initial_overnight_weight = max(mean(on_ret**2 / max(on_scale, min_var)), min_var)
    end function initial_overnight_weight

    ! Return whether parameters satisfy positivity and persistence constraints.
    logical function params_valid(params)
        type(mcsegarch_params_t), intent(in) :: params

        params_valid = abs(params%beta) < 1.0_dp
    end function params_valid

    ! Return distribution id for a configured distribution name.
    integer function dist_id(dist_name)
        character(len=*), intent(in) :: dist_name

        select case (trim(dist_name))
        case ("NORMAL")
            dist_id = dist_normal
        case ("T")
            dist_id = dist_t
        case ("FS_SKEWT")
            dist_id = dist_fs_skewt
        case default
            error stop "dist_id: unsupported distribution"
        end select
    end function dist_id

    ! Return the dynamic-parameter count for MCS-EGARCH.
    integer function model_param_count()

        model_param_count = 4
    end function model_param_count

    ! Return the innovation-distribution parameter count.
    integer function dist_param_count(dist)
        integer, intent(in) :: dist

        select case (dist)
        case (dist_normal)
            dist_param_count = 0
        case (dist_t)
            dist_param_count = 1
        case (dist_fs_skewt)
            dist_param_count = 2
        case default
            error stop "dist_param_count: unsupported distribution"
        end select
    end function dist_param_count

    ! Return the number of estimated parameters for the overnight intraday-impact mode.
    integer function lambda_param_count(lambda_mode)
        integer, intent(in) :: lambda_mode

        if (lambda_mode == lambda_est) then
            lambda_param_count = 2
        else
            lambda_param_count = 0
        end if
    end function lambda_param_count

    ! Return the starting or fixed overnight impact coefficient.
    real(dp) function initial_lambda_on(lambda_mode)
        integer, intent(in) :: lambda_mode

        select case (lambda_mode)
        case (lambda_zero)
            initial_lambda_on = 0.0_dp
        case (lambda_est)
            initial_lambda_on = 0.10_dp
        case (lambda_one)
            initial_lambda_on = 1.0_dp
        case default
            error stop "initial_lambda_on: unsupported lambda mode"
        end select
    end function initial_lambda_on

    ! Return the starting or fixed overnight GJR impact coefficient.
    real(dp) function initial_gamma_on(lambda_mode)
        integer, intent(in) :: lambda_mode

        if (lambda_mode < 1) error stop "initial_gamma_on: unsupported lambda mode"
        initial_gamma_on = 0.0_dp
    end function initial_gamma_on

    ! Return display label used in output tables.
    character(len=36) function model_table_label(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSEGARCH")
            model_table_label = "MCS-EGARCH"
        case default
            model_table_label = trim(model_name)
        end select
    end function model_table_label

    ! Count the number of intraday bins represented in the sample.
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

end module fit_mcsegarch_on_intraday_mod

program xfit_mcsegarch_on_intraday
    use date_mod, only: print_program_header
    use fit_mcsegarch_on_intraday_mod, only: run_fit_mcsegarch_on_intraday
    implicit none
    call print_program_header("xfit_mcsegarch_on_intraday.f90")
    call run_fit_mcsegarch_on_intraday()
end program xfit_mcsegarch_on_intraday
